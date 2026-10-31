"""Figures (PNG, 150 dpi) for the report."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

DPI = 150


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def peaks_per_label(stats: pd.DataFrame, path: Path) -> Path:
    """Peaks MACS3 called per label and how many of them are marker peaks."""
    plt = _plt()
    f, a = plt.subplots(figsize=(max(5, 0.5 * len(stats) + 2), 4))
    x = np.arange(len(stats))
    a.bar(x - 0.2, stats["n_peaks_called"], 0.4, color="#95a5a6", label="peaks called")
    a.bar(x + 0.2, stats["n_marker_peaks"], 0.4, color="#2c3e50", label="marker peaks")
    a.set_xticks(x, stats["label"], rotation=45, ha="right", fontsize=8)
    a.set_ylabel("peaks")
    a.legend(fontsize=8)
    a.set_title("Peaks per label")
    f.tight_layout()
    f.savefig(path, dpi=DPI)
    plt.close(f)
    return path


def marker_heatmap(pseudobulk: pd.DataFrame, markers: pd.DataFrame, path: Path, top: int = 20) -> Path | None:
    """Accessibility (log2 RPM + 1, z-scored per peak) of each label's top marker peaks across labels."""
    if markers.empty:
        return None
    plt = _plt()
    labels = list(pseudobulk.index)
    top_of = markers.sort_values("adjusted_p_value").groupby("label").head(top)
    peaks = list(dict.fromkeys(p for t in labels for p in top_of.loc[top_of["label"] == t, "peak"]))   # grouped by label
    m = np.log2(pseudobulk[peaks].to_numpy(dtype=float) + 1)
    z = (m - m.mean(axis=0)) / np.where(m.std(axis=0) > 0, m.std(axis=0), 1)
    f, a = plt.subplots(figsize=(8, max(3, 0.35 * len(labels) + 1.5)))
    im = a.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5, interpolation="nearest")
    a.set_yticks(range(len(labels)), labels, fontsize=8)
    a.set_xticks([])
    a.set_xlabel(f"top {top} marker peaks per label ({len(peaks)} peaks)")
    f.colorbar(im, ax=a, label="z-score")
    a.set_title("Marker peak accessibility")
    f.tight_layout()
    f.savefig(path, dpi=DPI)
    plt.close(f)
    return path


def motif_heatmap(motifs: pd.DataFrame, path: Path, fdr: float, top: int = 5) -> Path | None:
    """-log10 adjusted p of each label's top enriched motifs, across labels."""
    up = motifs[(motifs["adjusted_p_value"] < fdr) & (motifs["log2_fold_change"] > 0)]
    if up.empty:
        return None
    plt = _plt()
    names = up.sort_values("p_value").groupby("label", sort=False).head(top)["motif_name"].drop_duplicates().tolist()
    mat = (motifs[motifs["motif_name"].isin(names)].assign(score=lambda d: np.where(
        d["log2_fold_change"] > 0, -np.log10(d["adjusted_p_value"].clip(lower=1e-300)), 0.0))
        .pivot_table(index="label", columns="motif_name", values="score", aggfunc="max")[names])
    f, a = plt.subplots(figsize=(max(5, 0.35 * len(names) + 2), max(3, 0.35 * len(mat) + 1.5)))
    im = a.imshow(mat.to_numpy(), aspect="auto", cmap="viridis", interpolation="nearest")
    a.set_yticks(range(len(mat)), mat.index, fontsize=8)
    a.set_xticks(range(len(names)), names, rotation=90, fontsize=7)
    f.colorbar(im, ax=a, label="-log10 adjusted p (enriched)")
    a.set_title("Motif enrichment in marker peaks (JASPAR)")
    f.tight_layout()
    f.savefig(path, dpi=DPI)
    plt.close(f)
    return path


def links_plot(links: pd.DataFrame, path: Path) -> Path | None:
    """Kept peak-to-gene links: correlation against distance to the gene's TSS."""
    if links.empty:
        return None
    plt = _plt()
    f, a = plt.subplots(figsize=(6, 4))
    a.scatter(links["distance_to_tss"] / 1000, links["correlation"], s=6, alpha=0.5, c="#2c3e50", lw=0)
    a.set_xlabel("distance to TSS (kb)")
    a.set_ylabel("Spearman correlation (pseudo-cells)")
    a.set_title(f"Peak-to-gene links ({len(links):,})")
    f.tight_layout()
    f.savefig(path, dpi=DPI)
    plt.close(f)
    return path
