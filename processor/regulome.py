"""
The analysis, with SnapATAC2 (Zhang et al. 2024) throughout:

1. The samples' atac-prep files joined in a backed AnnDataSet (nothing copied); cells get their
   label from labels.csv by obs_names; labels with fewer than MIN_CELLS cells are left out.
2. Pseudobulk peaks per label: snap.tl.macs3 (MACS3 on the label's pooled Tn5 insertions,
   q < PEAK_QVALUE), chrM peaks dropped, snap.tl.merge_peaks: summits ± 250 bp, overlaps
   resolved by significance, giving one non-overlapping 501 bp peak set.
3. Peak matrix: snap.pp.make_peak_matrix over the labelled cells.
4. Marker peaks per label: snap.tl.diff_test (logistic-regression likelihood-ratio test on
   binarised accessibility, log total counts as covariate) of up to MAX_CELLS_PER_TEST cells of the
   label against an equal-share background from every other label, over the peaks called in the
   label, positive direction; BH-adjusted p < MARKER_FDR. This is the SnapATAC2 tutorial's design.
5. Motifs: snap.tl.motif_enrichment of every label's marker peaks against the union of all
   labels' marker peaks (hypergeometric), JASPAR CORE vertebrates PFMs, sequences from the
   genome FASTA baked into the image.
6. Peak-to-gene links (RNA linked): cells matched by obs_names; pseudo-cells by
   snap.tl.aggregate_cells on the RNA embedding; snap.tl.init_network_from_annotation (peaks
   within LINK_WINDOW_KB of a protein-coding gene's TSS) and, as snap.tl.add_cor_scores does, the
   Spearman correlation over pseudo-cells of RPM accessibility and CPM expression, for peaks and
   genes detected in >= 20 % of the pseudo-cells (computed here:
   add_cor_scores misaligns columns in 2.10.0); a t-approximation p-value per link,
   BH-adjusted; links with rho >= LINK_MIN_COR and FDR < 0.05 are kept.
"""
from __future__ import annotations

import gzip
import logging
import re
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

log = logging.getLogger("atac-regulome")
LINK_FDR = 0.05
MOTIF_FDR = 0.05
MIN_MARKERS_FOR_MOTIFS = 10
MIN_DETECTED_FRACTION = 0.2
GC_BINS = 20
BACKGROUND_FACTOR = 3


def bh(p: np.ndarray) -> np.ndarray:
    from scipy.stats import false_discovery_control

    return false_discovery_control(p, method="bh") if len(p) else np.array([])


def split_peak(names) -> pd.DataFrame:
    parts = pd.Series(list(names), dtype=str).str.extract(r"^(.+):(\d+)-(\d+)$")
    return pd.DataFrame({"chrom": parts[0], "start": parts[1].astype(np.int64), "end": parts[2].astype(np.int64)})


def safe_name(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("_") or "label"


def gene_names(annotation: Path) -> dict[str, str]:
    """Upper-case gene name -> the annotation's spelling (SnapATAC2's network upper-cases them)."""
    out = {}
    with gzip.open(annotation, "rt") as fh:
        for line in fh:
            if "\tgene\t" in line:
                m = re.search(r"gene_name=([^;\n]+)", line)
                if m:
                    out.setdefault(m.group(1).upper(), m.group(1))
    return out


UNUSED = "__not_analysed__"


def call_peaks(ds, groups: list[str], labels: list[str], params: dict, chrom_sizes: dict, threads: int,
               blacklist: Path | None):
    """MACS3 per label. Cells outside `labels` form one extra group whose peaks are discarded:
    snap.tl.macs3's `selections` returns nothing in SnapATAC2 2.10.0 (it compares label names with
    the internal group indices)."""
    import polars as pl
    import snapatac2 as snap

    groups = [g if g in labels else UNUSED for g in groups]
    peaks = snap.tl.macs3(ds, groupby=groups, qvalue=params["PEAK_QVALUE"], inplace=False, n_jobs=threads,
                          tempdir=Path("/tmp"), blacklist=blacklist)
    tempfile.tempdir = None   # with n_jobs=1 snap.tl.macs3 leaves it pointing at its deleted temporary folder
    peaks = {k: v.filter(pl.col("chrom") != "chrM") for k, v in peaks.items() if k in labels}
    merged = snap.tl.merge_peaks(peaks, chrom_sizes)
    return {k: peaks[k].height for k in labels if k in peaks}, merged


def peak_network(peaks: list[str], annotation: Path, window_kb: int):
    import snapatac2 as snap

    return snap.tl.init_network_from_annotation(peaks, annotation, upstream=window_kb * 1000,
                                                downstream=window_kb * 1000, id_type="gene_name")


def nearest_genes(network, names: dict[str, str]) -> dict[str, tuple[str, int]]:
    out: dict[str, tuple[str, int]] = {}
    for u, v, e in network.weighted_edge_list():
        peak, gene = network[u].id, names.get(network[v].id, network[v].id)
        if peak not in out or abs(e.distance) < abs(out[peak][1]):
            out[peak] = (gene, int(e.distance))
    return out


def lrt_pvalues(cov: np.ndarray, z: np.ndarray, Y: sp.csc_matrix, chunk: int = 4096, iters: int = 100) -> np.ndarray:
    """
    The p-values of snap.tl.diff_test's likelihood-ratio test, for every column of Y at once: per column a
    logistic regression of the binarised accessibility on [1, cov] (null) and on [1, cov, z] (alternative),
    unpenalised, chi-square with 1 df on twice the log-likelihood gain. SnapATAC2 2.10.0 fits the two models
    per peak with scikit-learn, one after another (~40 peaks/s on a Fargate vCPU); here Newton's method runs
    on a block of columns at a time with matrix products, giving the same maximum-likelihood fits.
    Columns whose fit does not come out finite are fitted by SnapATAC2's own function.
    """
    from scipy.stats import chi2

    n = len(cov)
    cov = cov - cov.mean()     # centred: the same likelihoods, a better-conditioned Newton step
    designs = (np.c_[np.ones(n), cov], np.c_[np.ones(n), cov, z])
    out = np.empty(Y.shape[1])
    for s in range(0, Y.shape[1], chunk):
        y = (Y[:, s:s + chunk].toarray() > 0).astype(np.float64)
        ll = []
        for D in designs:
            q = D.shape[1]
            B = np.zeros((q, y.shape[1]))
            pairs = [(i, j, D[:, i] * D[:, j]) for i in range(q) for j in range(i, q)]
            for _ in range(iters):
                p = 1 / (1 + np.exp(-np.clip(D @ B, -30, 30)))
                w = p * (1 - p)
                H = np.empty((y.shape[1], q, q))
                for i, j, dd in pairs:
                    H[:, i, j] = H[:, j, i] = dd @ w
                step = np.linalg.solve(H + 1e-9 * np.eye(q), (D.T @ (y - p)).T[..., None])[..., 0].T
                # damped along the Newton direction: a separable column walks off instead of jumping
                B += step * np.minimum(1, 5 / np.maximum(np.abs(step).max(axis=0), 1e-300))
                if np.abs(step).max() < 1e-8:
                    break
            p = np.clip(1 / (1 + np.exp(-np.clip(D @ B, -30, 30))), 1e-15, 1 - 1e-15)
            ll.append((y * np.log(p) + (1 - y) * np.log1p(-p)).sum(axis=0))
        out[s:s + chunk] = chi2.sf(np.maximum(2 * (ll[1] - ll[0]), 0), 1)
    bad = np.flatnonzero(~np.isfinite(out))
    if len(bad):
        from snapatac2.tools._diff import _likelihood_ratio_test_many

        log.warning("%d column(s) refitted one by one", len(bad))
        out[bad] = _likelihood_ratio_test_many(cov[:, None], z[:, None], Y[:, bad].copy())
    return out


def diff_test(X: sp.spmatrix, var_names, fg_rows: np.ndarray, bg_rows: np.ndarray, features: np.ndarray,
              min_pct: float = 0.05, min_log_fc: float = 0.25) -> pd.DataFrame:
    """
    snap.tl.diff_test(direction="positive") with its defaults, over the rows fg_rows (group 1) and bg_rows
    (group 2) of the cell-by-peak matrix X: a peak is tested when it is open in >= min_pct of either group
    and log2 of its RPM (count + 1, over all peaks) in group 1 over group 2 is >= min_log_fc; the covariate is
    log1p of each cell's total count; BH-adjusted. Returns peak, log2_fold_change, p_value, adjusted_p_value.
    """
    M = sp.csc_matrix(X[np.r_[fg_rows, bg_rows]])
    n1 = len(fg_rows)
    rpm = lambda m: (lambda x: x / (x.sum() / 1e6))(np.ravel(m.sum(axis=0)) + 1)  # noqa: E731
    log_fc = np.log2(rpm(M[:n1]) / rpm(M[n1:]))
    open1 = np.ravel((M[:n1] != 0).sum(axis=0))
    open2 = np.ravel((M[n1:] != 0).sum(axis=0))
    keep = ((open1 / n1 >= min_pct) | (open2 / (M.shape[0] - n1) >= min_pct)) & (log_fc >= min_log_fc)
    features = np.asarray(features)[keep[features]]
    cols = ["peak", "log2_fold_change", "p_value", "adjusted_p_value"]
    log.info("%d of the peaks pass min_pct %g / min_log_fc %g; testing them", len(features), min_pct, min_log_fc)
    if not len(features):
        return pd.DataFrame(columns=cols)
    cov = np.log1p(np.ravel(M.sum(axis=1)))
    z = np.r_[np.zeros(n1), np.ones(M.shape[0] - n1)]
    p = lrt_pvalues(cov, z, M[:, features])
    df = pd.DataFrame({"peak": np.asarray(var_names)[features], "log2_fold_change": log_fc[features], "p_value": p,
                       "adjusted_p_value": bh(p)})
    return df.sort_values("adjusted_p_value", kind="stable").reset_index(drop=True)


def marker_peaks(pm, labels: list[str], merged, params: dict, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    lab = np.asarray(pm.obs["label"].astype(str))
    X = sp.csr_matrix(pm.X)
    cap = params["MAX_CELLS_PER_TEST"]
    out = []
    for t in labels:
        others = [o for o in labels if o != t]
        if not others:
            break
        fg = np.flatnonzero(lab == t)
        fg = rng.choice(fg, size=min(cap, len(fg)), replace=False)
        share = max(1, -(-cap // len(others)))
        bg = np.concatenate([rng.choice(np.flatnonzero(lab == o), size=min(share, int((lab == o).sum())), replace=False)
                             for o in others])
        feats = np.flatnonzero(merged[t].to_numpy())
        # a peak open in every tested cell cannot be tested (the logistic fit needs both outcomes)
        sub = sp.csc_matrix(X[np.r_[fg, bg]][:, feats])
        feats = feats[np.diff(sub.indptr) < len(fg) + len(bg)]
        if not len(feats):
            continue
        df = diff_test(X, pm.var_names, fg, bg, feats)
        df = df[df["adjusted_p_value"] < params["MARKER_FDR"]]
        df.insert(0, "label", t)
        df["n_cells"], df["n_background"] = len(fg), len(bg)
        out.append(df)
    cols = ["label", "peak", "log2_fold_change", "p_value", "adjusted_p_value", "n_cells", "n_background"]
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=cols)


def plain_fasta(genome_dir: Path, work: Path) -> Path:
    """genome.fa.gz of the image, decompressed into the work folder with its prebuilt index (pyfaidx
    reads a compressed chromosome from its start for every region, far too slow for thousands)."""
    import shutil

    fa = work / "genome.fa"
    with gzip.open(genome_dir / "genome.fa.gz", "rb") as fin, open(fa, "wb") as fout:
        shutil.copyfileobj(fin, fout, 1 << 24)
    shutil.copyfile(genome_dir / "genome.fa.fai", work / "genome.fa.fai")    # newer than genome.fa: reused
    return fa


def gc_content(fasta: Path, peaks: list[str]) -> np.ndarray:
    from pyfaidx import Fasta

    loc = split_peak(peaks)
    gc = np.empty(len(peaks))
    with Fasta(str(fasta), one_based_attributes=False, sequence_always_upper=True) as g:
        for i, (c, a, b) in enumerate(zip(loc["chrom"], loc["start"], loc["end"])):
            seq = g[c][int(a):int(b)].seq
            gc[i] = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
    return gc


def gc_matched_background(peaks: list[str], gc: np.ndarray, markers: set[str], fg: list[str], seed: int) -> list[str]:
    """Peaks that are no label's marker, drawn so their GC distribution (20 bins) matches `fg`'s:
    BACKGROUND_FACTOR times as many per bin as `fg` has (a larger background steadies the binomial
    test's expected rate), fewer when the bin runs out."""
    rng = np.random.default_rng(seed)
    bins = np.minimum((gc * GC_BINS).astype(int), GC_BINS - 1)
    pool = np.array([p not in markers for p in peaks])
    bin_of = dict(zip(peaks, bins))
    want = BACKGROUND_FACTOR * np.bincount([bin_of[p] for p in fg], minlength=GC_BINS)
    out = []
    for b in range(GC_BINS):
        cand = np.flatnonzero(pool & (bins == b))
        if want[b] and len(cand):
            out += [peaks[i] for i in rng.choice(cand, size=min(int(want[b]), len(cand)), replace=False)]
    return out


def motifs(markers: pd.DataFrame, peaks: list[str], motif_file: Path, genome_dir: Path, work: Path,
           seed: int) -> tuple[pd.DataFrame, dict]:
    import snapatac2 as snap

    regions = {t: list(g["peak"]) for t, g in markers.groupby("label", sort=False)
               if len(g) >= MIN_MARKERS_FOR_MOTIFS}
    cols = ["label", "motif_id", "motif_name", "log2_fold_change", "p_value", "adjusted_p_value"]
    if not regions:
        log.warning("motif enrichment skipped: no label has >= %d marker peaks", MIN_MARKERS_FOR_MOTIFS)
        return pd.DataFrame(columns=cols), {}
    fasta = plain_fasta(genome_dir, work)
    gc = gc_content(fasta, peaks)
    gc_of = dict(zip(peaks, gc))
    union = {p for r in regions.values() for p in r}
    pfms = snap.read_motifs(str(motif_file))
    for m in pfms:                       # MEME "MOTIF MA0004.1 Arnt": SnapATAC2 keeps both words in id
        m.id, _, name = m.id.partition(" ")
        m.name = name or m.id
    # each label against its own GC-matched draw from the peaks that are no label's marker
    backgrounds = {t: gc_matched_background(peaks, gc, union, fg, seed) for t, fg in regions.items()}
    info = {t: {"marker_gc": round(float(np.mean([gc_of[p] for p in fg])), 3), "background_peaks": len(backgrounds[t]),
                "background_gc": round(float(np.mean([gc_of[p] for p in backgrounds[t]])), 3) if backgrounds[t] else None}
            for t, fg in regions.items()}
    out = [pd.DataFrame({"label": t, **d})[cols]
           for t, d in motif_enrichment(pfms, regions, backgrounds, fasta).items()]
    return pd.concat(out, ignore_index=True).sort_values(["label", "p_value"], kind="stable"), info


def motif_enrichment(pfms, regions: dict[str, list[str]], backgrounds: dict[str, list[str]], fasta: Path) -> dict[str, dict]:
    """
    snap.tl.motif_enrichment(method="binomial") of each label's regions against that label's own background,
    with every distinct sequence scanned once: SnapATAC2 2.10.0 takes one background per call, so a call per
    label scanned the overlapping backgrounds again for each label. Per motif and label, the same counts,
    fold change and binomial p-value as SnapATAC2, BH-adjusted per label.
    """
    from pyfaidx import Fasta
    from scipy.stats import binom
    from snapatac2._utils import fetch_seq
    from snapatac2.tools._diff import _p_adjust_bh

    every = list(dict.fromkeys(p for t in regions for p in (*regions[t], *backgrounds[t])))
    index = {p: i for i, p in enumerate(every)}
    log.info("scanning %d motifs over %d distinct sequences (%d label + background sequences in all)", len(pfms),
             len(every), sum(len(regions[t]) + len(backgrounds[t]) for t in regions))
    genome = Fasta(str(fasta), one_based_attributes=False)
    sequences = [fetch_seq(genome, p) for p in every]
    fg_idx = {t: np.array([index[p] for p in regions[t]]) for t in regions}
    bg_idx = {t: np.array([index[p] for p in backgrounds[t]], dtype=int) for t in regions}
    counts = {t: ([], []) for t in regions}
    for m in pfms:
        bound = np.asarray(m.with_nucl_prob().exists(sequences), dtype=bool)
        for t in regions:
            counts[t][0].append(int(bound[fg_idx[t]].sum()))
            counts[t][1].append(int(bound[bg_idx[t]].sum()))
    out = {}
    for t in regions:
        n_fg, n_bg = np.array(counts[t][0]), np.array(counts[t][1])
        N_fg, N_bg = len(regions[t]), len(backgrounds[t])
        with np.errstate(divide="ignore", invalid="ignore"):
            fc = np.where(n_fg == 0, np.where(n_bg == 0, 0.0, -np.inf),
                          np.where(n_bg > 0, np.log2((n_fg / N_fg) / (n_bg / max(N_bg, 1))), np.inf))
            rate = n_bg / N_bg if N_bg else np.full(len(n_bg), np.nan)
            pval = np.where(fc >= 0, binom.sf(n_fg - 1, N_fg, rate), binom.cdf(n_fg, N_fg, rate))
        pval = np.clip(pval, 1e-300, 1)
        out[t] = {"motif_id": [m.id for m in pfms], "motif_name": [m.name for m in pfms], "log2_fold_change": fc,
                  "p_value": pval, "adjusted_p_value": _p_adjust_bh(list(pval))}
    return out


def rna_matrix(rna_files: list[Path], cells: list[str]) -> tuple[list[str], sp.csr_matrix, np.ndarray, list[str]]:
    """The cells of `cells` (in that order) found in the RNA files: their counts (layers['counts'],
    else X) over the genes every file has, and an embedding: sc-integrate's obsm X_integrated (one
    space for all samples) when every file has it, else a 30-component truncated SVD of log CPM over
    all matched cells (per-sample X_pca are separate spaces and cannot be stacked)."""
    import anndata as ad

    want = set(cells)
    names, mats, embs, gene_lists = [], [], [], []
    for f in rna_files:
        a = ad.read_h5ad(f)
        hit = a.obs_names.isin(want)
        if not hit.any():
            continue
        a = a[hit]
        names += list(a.obs_names)
        mats.append((list(a.var_names), sp.csr_matrix(a.layers["counts"] if "counts" in a.layers else a.X)))
        embs.append(np.asarray(a.obsm["X_integrated"]) if "X_integrated" in a.obsm else None)
        gene_lists.append(list(a.var_names))
    if not names:
        return [], sp.csr_matrix((0, 0)), np.zeros((0, 0)), []
    common = set.intersection(*(set(g) for g in gene_lists))
    genes = [g for g in gene_lists[0] if g in common]
    cols = []
    for g, M in mats:
        index = {x: i for i, x in enumerate(g)}
        cols.append(M[:, [index[x] for x in genes]])
    X = sp.vstack(cols).tocsr()
    order = {n: i for i, n in enumerate(names)}
    matched = [c for c in cells if c in order]
    sel = [order[c] for c in matched]
    X = X[sel]
    if all(e is not None for e in embs) and len({e.shape[1] for e in embs}) == 1:
        emb = np.vstack(embs)[sel]
    else:
        from sklearn.decomposition import TruncatedSVD

        L = sp.csr_matrix(X.multiply(1e4 / np.maximum(np.asarray(X.sum(axis=1)).ravel(), 1)[:, None]))
        L.data = np.log1p(L.data)
        emb = TruncatedSVD(n_components=min(30, min(L.shape) - 1), random_state=0).fit_transform(L)
    return matched, X, emb, genes


def links(pm, rna_files: list[Path], network, names: dict[str, str], params: dict, seed: int) -> tuple[pd.DataFrame, dict]:
    import snapatac2 as snap
    from scipy.stats import t as tdist

    cols = ["peak", "gene", "distance_to_tss", "correlation", "p_value", "adjusted_p_value"]
    matched, X, emb, genes = rna_matrix(rna_files, list(pm.obs_names))
    info = {"n_matched_cells": len(matched)}
    if len(matched) < 2 * params["PSEUDOCELL_SIZE"]:
        info["skipped"] = (f"only {len(matched)} nuclei are in both the ATAC and the RNA data; links need at least "
                           f"{2 * params['PSEUDOCELL_SIZE']} (twice PSEUDOCELL_SIZE). Check that the RNA cells are the same "
                           "nuclei, with the same barcodes and sample ids")
        log.warning("peak-to-gene links skipped: %s", info["skipped"])
        return pd.DataFrame(columns=cols), info
    pseudo = np.asarray(snap.tl.aggregate_cells(np.ascontiguousarray(emb, dtype=np.float64), inplace=False,
                                                min_cluster_size=params["PSEUDOCELL_SIZE"], random_state=seed), dtype=str)
    ok = pseudo != "-1"
    groups = sorted(set(pseudo[ok]))
    n = len(groups)
    info.update(n_pseudocells=n, n_genes_rna=len(genes))
    if n < 5:
        info["skipped"] = (f"the {len(matched):,} nuclei in both the ATAC and the RNA data formed {n} pseudo-cells of at least "
                           f"{params['PSEUDOCELL_SIZE']} nuclei; links need at least 5. More nuclei, or a smaller "
                           "PSEUDOCELL_SIZE, gives more pseudo-cells")
        log.warning("peak-to-gene links skipped: %s", info["skipped"])
        return pd.DataFrame(columns=cols), info
    gidx = {g: i for i, g in enumerate(groups)}
    G = sp.csr_matrix((np.ones(ok.sum()), ([gidx[p] for p in pseudo[ok]], np.flatnonzero(ok))), shape=(n, len(matched)))
    rpm = lambda M: sp.csr_matrix(M.multiply(1e6 / np.maximum(np.asarray(M.sum(axis=1)).ravel(), 1)[:, None]))  # noqa: E731
    acc = rpm(G @ sp.csr_matrix(pm[matched].X)).toarray()
    expr = rpm(G @ X).toarray()
    # Spearman over pseudo-cells for every peak-gene edge of the network (SnapATAC2's add_cor_scores does the
    # same, but in 2.10.0 it misaligns columns when every peak of the matrix is in the network)
    from scipy.stats import rankdata

    def ranked(M):
        R_ = rankdata(M, axis=0)
        R_ -= R_.mean(axis=0)
        norm = np.linalg.norm(R_, axis=0)
        return R_ / np.where(norm > 0, norm, np.nan)

    # a peak or gene seen in only a few pseudo-cells gives ranks dominated by ties (ρ near 1 from one shared
    # pseudo-cell): both must be non-zero in at least MIN_DETECTED_FRACTION of the pseudo-cells (and >= 5)
    need = max(5, int(np.ceil(MIN_DETECTED_FRACTION * n)))
    acc_ok, expr_ok = (acc > 0).sum(axis=0) >= need, (expr > 0).sum(axis=0) >= need
    acc_r, expr_r = ranked(acc), ranked(expr)
    pcol = {p: i for i, p in enumerate(pm.var_names)}
    gcol = {}
    for i, g in enumerate(genes):
        gcol.setdefault(g.upper(), i)
    edges = [(network[u].id, network[v].id, int(e.distance)) for u, v, e in network.weighted_edge_list()
             if network[v].id in gcol and acc_ok[pcol[network[u].id]] and expr_ok[gcol[network[v].id]]]
    pi = np.array([pcol[p] for p, _, _ in edges], dtype=int)
    gi = np.array([gcol[g] for _, g, _ in edges], dtype=int)
    cor = np.einsum("ij,ij->j", acc_r[:, pi], expr_r[:, gi]) if len(edges) else np.array([])
    rows = [(p, names.get(g, g), d, float(c)) for (p, g, d), c in zip(edges, cor) if np.isfinite(c)]
    df = pd.DataFrame(rows, columns=["peak", "gene", "distance_to_tss", "correlation"])
    r = df["correlation"].clip(-0.999999, 0.999999).to_numpy()
    df["p_value"] = 2 * tdist.sf(np.abs(r * np.sqrt((n - 2) / (1 - r ** 2))), n - 2)
    df["adjusted_p_value"] = bh(df["p_value"].to_numpy())
    info["n_tested"] = len(df)
    df = df[(df["correlation"] >= params["LINK_MIN_COR"]) & (df["adjusted_p_value"] < LINK_FDR)]
    return df.sort_values(["adjusted_p_value", "peak"]).reset_index(drop=True), info
