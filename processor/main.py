"""
atac-regulome: per-label chromatin accessibility of the Multiome Gene Regulation workflow.

Reads atac-prep's <sample_id>.h5ad files (port `atac`), labels.csv (port `labels`, e.g. cell
types from sc-annotate) and optionally the RNA side (port `cells`), and writes under
OUTPUT_DIR/<LABEL>/

  peaks.bed, peaks/<label>.bed   the merged peak set; each label's own peaks
  marker_peaks.csv               peaks more accessible in a label than in the others
  motifs.csv                     JASPAR motif enrichment of each label's marker peaks
  links.csv                      peak-to-gene links (only with RNA linked)
  *.png, summary.json, provenance.json

See processor/regulome.py for the method.
"""
from __future__ import annotations

import json
import logging
import math
import os
import platform
import re
import shutil
import sys
import tempfile
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("atac-regulome")
for noisy in ("matplotlib", "PIL", "fontTools", "numba", "h5py", "anndata"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

warnings.filterwarnings("ignore", category=FutureWarning)   # scikit-learn deprecations inside SnapATAC2
APP_NAME = "atac-regulome"
REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULTS = {
    "LABEL_COLUMN": "cell_type",
    "LABEL": "",
    "MIN_CELLS": "50",
    "PEAK_QVALUE": "0.05",
    "MAX_CELLS_PER_TEST": "500",
    "MARKER_FDR": "0.01",
    "LINK_WINDOW_KB": "250",
    "LINK_MIN_COR": "0.3",
    "PSEUDOCELL_SIZE": "50",
    "SEED": "42",
    "THREADS": "0",
}
PARAMETERS = tuple(DEFAULTS)
CITATIONS = [   # each checked on Crossref and PubMed (README "Citations")
    {"doi": "10.1038/s41592-023-02139-9", "pmid": "38191932", "text": "Zhang K et al. Nat Methods 2024"},
    {"doi": "10.1186/gb-2008-9-9-r137", "pmid": "18798982", "text": "Zhang Y et al. Genome Biol 2008"},
    {"doi": "10.1093/nar/gkaf1209", "pmid": "41325984", "text": "Ovek Baydar D et al. Nucleic Acids Res 2026"},
    {"doi": "10.1093/nar/gkac1071", "pmid": "36420896", "text": "Frankish A et al. Nucleic Acids Res 2023"},
    {"doi": "10.1038/s41576-023-00586-w", "pmid": "37002403", "text": "Heumos L et al. Nat Rev Genet 2023"},
]


class ParameterError(ValueError):
    pass


def read_parameters(env=os.environ) -> dict:
    raw = {k: (env.get(k) if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p: dict = {"LABEL_COLUMN": raw["LABEL_COLUMN"].strip()}
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", p["LABEL_COLUMN"]):
        raise ParameterError(f"LABEL_COLUMN must be a column name, got {p['LABEL_COLUMN']!r}")
    p["LABEL"] = raw["LABEL"].strip() or p["LABEL_COLUMN"].replace("_", "-")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", p["LABEL"]):
        raise ParameterError(f"LABEL must be a folder name (letters, digits, . _ -), got {p['LABEL']!r}")
    try:
        for k in ("MIN_CELLS", "MAX_CELLS_PER_TEST", "LINK_WINDOW_KB", "PSEUDOCELL_SIZE", "SEED", "THREADS"):
            p[k] = int(raw[k])
        for k in ("PEAK_QVALUE", "MARKER_FDR", "LINK_MIN_COR"):
            p[k] = float(raw[k])
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if p["MIN_CELLS"] < 1 or p["MAX_CELLS_PER_TEST"] < 10 or p["PSEUDOCELL_SIZE"] < 2 or p["THREADS"] < 0:
        raise ParameterError("need MIN_CELLS >= 1, MAX_CELLS_PER_TEST >= 10, PSEUDOCELL_SIZE >= 2, THREADS >= 0")
    if not 1 <= p["LINK_WINDOW_KB"] <= 1000:
        raise ParameterError("LINK_WINDOW_KB must be 1-1000")
    for k in ("PEAK_QVALUE", "MARKER_FDR"):
        if not 0 < p[k] < 1:
            raise ParameterError(f"{k} must be between 0 and 1")
    if not -1 <= p["LINK_MIN_COR"] <= 1:
        raise ParameterError("LINK_MIN_COR must be between -1 and 1")
    return p


def write_json(path: Path, data) -> None:
    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if hasattr(o, "item"):
            o = o.item()
        if isinstance(o, float) and not math.isfinite(o):
            return None
        return str(o) if isinstance(o, Path) else o

    with open(path, "w") as fh:
        json.dump(clean(data), fh, indent=2)
        fh.write("\n")


def tool_versions() -> dict:
    from importlib.metadata import PackageNotFoundError, version

    out = {"python": platform.python_version()}
    for name in ("snapatac2", "MACS3", "anndata", "numpy", "polars", "pandas", "scipy", "scikit-learn", "pyfaidx",
                 "igraph", "matplotlib"):
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            out[name] = None
    return out


def peak_rss_mb() -> float | None:
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)
    except (ImportError, OSError):
        return None


def cgroup_peak_mb() -> float | None:
    try:
        return round(int(open("/sys/fs/cgroup/memory.peak").read()) / 2**20, 1)
    except (OSError, ValueError):
        return None


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def genome_of(files: list[Path]) -> str:
    import snapatac2 as snap

    names = set()
    for f in files:
        a = snap.read(f, backed="r")
        try:
            names.add(str((a.uns["pennsieve"] or {}).get("genome", "")))
        finally:
            a.close()
    if len(names) != 1 or "" in names:
        raise ValueError(f"the atac files must name one genome in uns['pennsieve']['genome'] (found {sorted(names)})")
    return names.pop()


def main() -> int:
    started = datetime.now(timezone.utc)
    t0 = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("bad parameter: %s", exc)
        return 2
    threads = params["THREADS"] or int(os.environ.get("OMP_NUM_THREADS") or 2)

    import numpy as np
    import pandas as pd
    import snapatac2 as snap

    from processor import figures as F
    from processor import inputs as I
    from processor import regulome as R

    timings: dict[str, float] = {}
    try:
        found = I.discover(input_dir)
        if not found.atac:
            raise I.InputError("no atac-prep .h5ad (with obsm['fragment_paired']) found")
        if found.labels is None:
            raise I.InputError("no labels.csv found")
        labels_of = I.read_labels(found.inp, found.labels, params["LABEL_COLUMN"])
        genome = genome_of(found.atac)
    except (I.InputError, OSError, ValueError) as exc:
        log.error("cannot use the inputs: %s", exc)
        return 1
    gdir = Path(os.environ.get("GENOMES_DIR", "/opt/genomes")) / genome
    motif_file = Path(os.environ.get("MOTIFS_FILE", "/opt/motifs/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt"))
    chrom_sizes = {c: int(n) for c, n in (line.split("\t") for line in (gdir / "chrom_sizes.tsv").read_text().splitlines())}
    rel = lambda p: str(p.relative_to(input_dir))  # noqa: E731
    log.info("input: %d atac file(s), labels %s (column %s), %d RNA file(s), genome %s, threads %d", len(found.atac),
             rel(found.labels), params["LABEL_COLUMN"], len(found.rna), genome, threads)

    out = output_dir / params["LABEL"]
    out.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="regulome-", dir="/tmp"))
    ds = None
    try:
        t = time.monotonic()
        ds = snap.AnnDataSet(adatas=[(f.stem, snap.read(f, backed="r")) for f in found.atac], filename=work / "data.h5ads")
        cells = np.asarray(ds.obs_names)
        lab = pd.Series(cells).map(labels_of).fillna("").to_numpy()
        counts = pd.Series(lab[lab != ""]).value_counts()
        kept = [str(x) for x in counts.index if counts[x] >= params["MIN_CELLS"]]
        small = {str(x): int(counts[x]) for x in counts.index if counts[x] < params["MIN_CELLS"]}
        log.info("%d cells, %d labelled (%d labels), %d labels with >= %d cells: %s; left out: %s", len(cells),
                 int((lab != "").sum()), len(counts), len(kept), params["MIN_CELLS"],
                 ", ".join(f"{k} {counts[k]}" for k in kept), ", ".join(f"{k} {v}" for k, v in small.items()) or "none")
        if not kept:
            raise ValueError(f"no label has {params['MIN_CELLS']} or more cells (labels matched {int((lab != '').sum())} "
                             "of the atac cells; obs_id must equal the atac obs_names <sample_id>:<barcode>)")
        n_called, merged = R.call_peaks(ds, list(lab), kept, params, chrom_sizes, threads)
        timings["peaks"] = round(time.monotonic() - t, 1)
        peak_names = merged["Peaks"].to_list()
        log.info("MACS3: %s peaks; merged set %d peaks (%.1fs)", ", ".join(f"{k} {v}" for k, v in n_called.items()),
                 len(peak_names), timings["peaks"])

        t = time.monotonic()
        pm = snap.pp.make_peak_matrix(ds, use_rep=peak_names)
        sel = np.flatnonzero(np.isin(lab, kept))
        pm = pm[sel].copy()
        pm.obs["label"] = pd.Categorical(lab[sel])
        network = R.peak_network(peak_names, gdir / "annotation.gff3.gz", params["LINK_WINDOW_KB"])
        names = R.gene_names(gdir / "annotation.gff3.gz")
        nearest = R.nearest_genes(network, names)
        markers = R.marker_peaks(pm, kept, merged, params, params["SEED"])
        markers["nearest_gene"] = markers["peak"].map(lambda p: nearest.get(p, ("", None))[0])
        markers["distance_to_tss"] = markers["peak"].map(lambda p: nearest.get(p, ("", None))[1])
        markers = pd.concat([markers[["label", "peak"]], R.split_peak(markers["peak"]),
                             markers.drop(columns=["label", "peak"])], axis=1)
        timings["markers"] = round(time.monotonic() - t, 1)
        log.info("marker peaks (FDR < %g): %s (%.1fs)", params["MARKER_FDR"],
                 ", ".join(f"{k} {int((markers['label'] == k).sum())}" for k in kept), timings["markers"])

        t = time.monotonic()
        mot = R.motifs(markers, motif_file, gdir, work)
        timings["motifs"] = round(time.monotonic() - t, 1)

        t = time.monotonic()
        link_info: dict = {"rna_linked": bool(found.rna)}
        lk = pd.DataFrame()
        if found.rna:
            lk, info = R.links(pm, found.rna, network, names, params, params["SEED"])
            link_info.update(info)
            called = {p: ";".join(t_ for t_ in kept if v) for p, *vals in merged.select(["Peaks"] + kept).iter_rows()
                      for v in [any(vals)] if False}
            per_peak = merged.select(["Peaks"] + kept).to_pandas().set_index("Peaks")
            lk = pd.concat([lk[["peak"]], R.split_peak(lk["peak"]), lk.drop(columns=["peak"])], axis=1)
            lk["called_in"] = lk["peak"].map(lambda p: ";".join(c for c in kept if per_peak.at[p, c]))
            del called
            lk.to_csv(out / "links.csv", index=False)
        timings["links"] = round(time.monotonic() - t, 1)

        # files
        bed = R.split_peak(peak_names).assign(name=peak_names)
        bed.to_csv(out / "peaks.bed", sep="\t", header=False, index=False)
        (out / "peaks").mkdir(exist_ok=True)
        files = {}
        for k in kept:
            fn = f"{R.safe_name(k)}.bed"
            bed[merged[k].to_numpy()].to_csv(out / "peaks" / fn, sep="\t", header=False, index=False)
            files[k] = f"peaks/{fn}"
        markers.to_csv(out / "marker_peaks.csv", index=False)
        mot.to_csv(out / "motifs.csv", index=False)
        stats = pd.DataFrame({"label": kept, "n_cells": [int(counts[k]) for k in kept],
                              "n_peaks_called": [n_called.get(k, 0) for k in kept],
                              "n_peaks_in_merged_set": [int(merged[k].sum()) for k in kept],
                              "n_marker_peaks": [int((markers["label"] == k).sum()) for k in kept]})
        pseudo = snap.tl.aggregate_X(pm, groupby="label", normalize="RPM")
        pb = pd.DataFrame(np.asarray(pseudo.X), index=list(pseudo.obs_names), columns=list(pseudo.var_names))
        figs = [F.peaks_per_label(stats, out / "peaks_per_label.png"),
                F.marker_heatmap(pb, markers, out / "marker_peaks_heatmap.png"),
                F.motif_heatmap(mot, out / "motifs_heatmap.png", R.MOTIF_FDR),
                F.links_plot(lk, out / "links.png") if found.rna else None]
        top_motifs = {k: g[(g["adjusted_p_value"] < R.MOTIF_FDR) & (g["log2_fold_change"] > 0)]["motif_name"].head(5).tolist()
                      for k, g in mot.groupby("label", sort=False)}
        summary = {
            "app": APP_NAME, "version": app_version(), "label_column": params["LABEL_COLUMN"], "label": params["LABEL"],
            "genome": genome, "n_cells": int(len(cells)), "n_cells_labelled": int((lab != "").sum()),
            "n_cells_used": int(len(sel)), "labels_left_out": small, "n_peaks": len(peak_names),
            "labels": [{**r, "peaks_bed": files[r["label"]], "top_motifs": top_motifs.get(r["label"], []),
                        "n_links": int(lk["called_in"].str.contains(re.escape(r["label"])).sum()) if len(lk) else 0}
                       for r in stats.to_dict("records")],
            "n_marker_peaks": len(markers), "n_motif_tests": len(mot),
            "n_enriched_motifs": int(((mot["adjusted_p_value"] < R.MOTIF_FDR) & (mot["log2_fold_change"] > 0)).sum()) if len(mot) else 0,
            "links": {**link_info, "n_links": len(lk)}, "figures": [f.name for f in figs if f is not None],
        }
        write_json(out / "summary.json", summary)
    except Exception as exc:  # noqa: BLE001
        if isinstance(exc, ValueError):
            log.error("FAILED: %s", exc)
        else:
            log.exception("FAILED: %s", exc)
        return 1
    finally:
        if ds is not None:
            ds.close()
        shutil.rmtree(work, ignore_errors=True)

    write_json(out / "provenance.json", {
        "app": APP_NAME, "version": app_version(), "tools": tool_versions(), "parameters": params,
        "inputs": [rel(f) for f in found.atac] + [rel(found.labels)] + [rel(f) for f in found.rna],
        "started": started.isoformat(timespec="seconds"), "finished": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "seconds": round(time.monotonic() - t0, 1), "timings": timings, "threads": threads, "genome": genome,
        "motifs": motif_file.name, "samples": [{"sample_id": f.stem, "input": rel(f), "status": "ok"} for f in found.atac],
        "peak_rss_mb": peak_rss_mb(), "cgroup_peak_mb": cgroup_peak_mb(), "citations": CITATIONS,
    })
    log.info("%s: %d peaks, %d marker peaks, %d motif tests, %d links; peak RSS %s MB, %.1fs", out.name, summary["n_peaks"],
             summary["n_marker_peaks"], summary["n_motif_tests"], summary["links"]["n_links"], peak_rss_mb(),
             time.monotonic() - t0)
    return 0
