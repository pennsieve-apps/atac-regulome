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
   Spearman correlation over pseudo-cells of RPM accessibility and CPM expression (computed here:
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


def bh(p: np.ndarray) -> np.ndarray:
    from snapatac2.tools._diff import _p_adjust_bh

    return _p_adjust_bh(p) if len(p) else np.array([])


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


def call_peaks(ds, groups: list[str], labels: list[str], params: dict, chrom_sizes: dict, threads: int):
    """MACS3 per label. Cells outside `labels` form one extra group whose peaks are discarded:
    snap.tl.macs3's `selections` returns nothing in SnapATAC2 2.10.0 (it compares label names with
    the internal group indices)."""
    import polars as pl
    import snapatac2 as snap

    groups = [g if g in labels else UNUSED for g in groups]
    peaks = snap.tl.macs3(ds, groupby=groups, qvalue=params["PEAK_QVALUE"], inplace=False, n_jobs=threads,
                          tempdir=Path("/tmp"))
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


def marker_peaks(pm, labels: list[str], merged, params: dict, seed: int) -> pd.DataFrame:
    import snapatac2 as snap

    rng = np.random.default_rng(seed)
    lab = np.asarray(pm.obs["label"].astype(str))
    names = np.asarray(pm.obs_names)
    cap = params["MAX_CELLS_PER_TEST"]
    out = []
    for t in labels:
        others = [o for o in labels if o != t]
        if not others:
            break
        fg = names[lab == t]
        fg = rng.choice(fg, size=min(cap, len(fg)), replace=False)
        share = max(1, -(-cap // len(others)))
        bg = np.concatenate([rng.choice(names[lab == o], size=min(share, int((lab == o).sum())), replace=False)
                             for o in others])
        feats = np.flatnonzero(merged[t].to_numpy())
        # a peak open in every tested cell cannot be tested (the logistic fit needs both outcomes)
        rows = np.flatnonzero(np.isin(names, np.concatenate([fg, bg])))
        sub = sp.csc_matrix(pm.X[rows][:, feats])
        feats = feats[np.diff(sub.indptr) < len(rows)]
        if not len(feats):
            continue
        res = snap.tl.diff_test(pm, cell_group1=list(fg), cell_group2=list(bg), features=list(feats),
                                direction="positive")
        if res.height == 0:
            continue
        df = res.to_pandas().rename(columns={"feature name": "peak", "log2(fold_change)": "log2_fold_change",
                                             "p-value": "p_value", "adjusted p-value": "adjusted_p_value"})
        df = df[df["adjusted_p_value"] < params["MARKER_FDR"]]
        df.insert(0, "label", t)
        df["n_cells"], df["n_background"] = len(fg), len(bg)
        out.append(df)
    cols = ["label", "peak", "log2_fold_change", "p_value", "adjusted_p_value", "n_cells", "n_background"]
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=cols)


def motifs(markers: pd.DataFrame, motif_file: Path, fasta: Path) -> pd.DataFrame:
    import snapatac2 as snap

    regions = {t: list(g["peak"]) for t, g in markers.groupby("label", sort=False)
               if len(g) >= MIN_MARKERS_FOR_MOTIFS}
    cols = ["label", "motif_id", "motif_name", "log2_fold_change", "p_value", "adjusted_p_value"]
    if len(regions) < 2:
        log.warning("motif enrichment skipped: needs >= 2 labels with >= %d marker peaks (have %d)",
                    MIN_MARKERS_FOR_MOTIFS, len(regions))
        return pd.DataFrame(columns=cols)
    pfms = snap.read_motifs(str(motif_file))
    for m in pfms:                       # MEME "MOTIF MA0004.1 Arnt": SnapATAC2 keeps both words in id
        m.id, _, name = m.id.partition(" ")
        m.name = name or m.id
    res = snap.tl.motif_enrichment(pfms, regions, fasta)
    out = []
    for t, df in res.items():
        d = df.to_pandas().rename(columns={"id": "motif_id", "name": "motif_name", "log2(fold change)": "log2_fold_change",
                                           "p-value": "p_value", "adjusted p-value": "adjusted_p_value"})
        d.insert(0, "label", t)
        out.append(d[cols])
    return pd.concat(out, ignore_index=True).sort_values(["label", "p_value"], kind="stable")


def rna_matrix(rna_files: list[Path], cells: list[str]) -> tuple[list[str], sp.csr_matrix, np.ndarray, list[str]]:
    """The cells of `cells` (in that order) found in the RNA files: their counts (layers['counts'],
    else X) over the genes every file has, and an embedding (obsm X_integrated, else X_pca, else a
    30-component truncated SVD of log CPM)."""
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
        rep = next((k for k in ("X_integrated", "X_pca") if k in a.obsm), None)
        embs.append(np.asarray(a.obsm[rep]) if rep else None)
        gene_lists.append(list(a.var_names))
    if not names:
        return [], sp.csr_matrix((0, 0)), np.zeros((0, 0)), []
    common = set.intersection(*(set(g) for g in gene_lists))
    genes = [g for g in gene_lists[0] if g in common]
    X = sp.vstack([M[:, [{x: i for i, x in enumerate(g)}[x] for x in genes]] for g, M in mats]).tocsr()
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
        log.warning("peak-to-gene links skipped: %d cells matched between ATAC and RNA", len(matched))
        return pd.DataFrame(columns=cols), info
    pseudo = np.asarray(snap.tl.aggregate_cells(np.ascontiguousarray(emb, dtype=np.float64), inplace=False,
                                                min_cluster_size=params["PSEUDOCELL_SIZE"], random_state=seed), dtype=str)
    ok = pseudo != "-1"
    groups = sorted(set(pseudo[ok]))
    n = len(groups)
    info.update(n_pseudocells=n, n_genes_rna=len(genes))
    if n < 5:
        log.warning("peak-to-gene links skipped: %d pseudo-cells (need >= 5)", n)
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

    acc_r, expr_r = ranked(acc), ranked(expr)
    pcol = {p: i for i, p in enumerate(pm.var_names)}
    gcol = {}
    for i, g in enumerate(genes):
        gcol.setdefault(g.upper(), i)
    edges = [(network[u].id, network[v].id, int(e.distance)) for u, v, e in network.weighted_edge_list()
             if network[v].id in gcol]
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
