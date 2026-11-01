"""The batched marker test and the single-scan motif enrichment against SnapATAC2's own functions."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from processor import regulome as R
from tests.make_fixture import make_fixture


def binary_matrix(n1: int, n2: int, k: int, seed: int) -> tuple[np.ndarray, sp.csr_matrix]:
    """Cells x peaks counts with a depth effect, a group effect on a third of the peaks and a few separable peaks."""
    rng = np.random.default_rng(seed)
    z = np.r_[np.zeros(n1), np.ones(n2)]
    depth = rng.lognormal(8, 0.6, n1 + n2)
    base = rng.uniform(0.02, 0.4, k)
    eff = np.where(rng.random(k) < 0.3, rng.normal(0, 1.5, k), 0)
    logit = np.log(base / (1 - base))[None, :] + eff[None, :] * z[:, None] + 0.5 * (np.log(depth) - 8)[:, None]
    Y = (rng.random((n1 + n2, k)) < 1 / (1 + np.exp(-logit))).astype(np.float64)
    Y[:n1, :3] = 0                       # open only in group 2: perfectly separable
    Y[n1:, :3] = 1
    Y[:, 3] = 0
    Y[:5, 3] = 1                         # open in five cells
    counts = Y * rng.integers(1, 4, Y.shape)
    return z, sp.csr_matrix(counts)


def test_lrt_matches_snapatac2():
    """SnapATAC2's per-peak scikit-learn fits stop at L-BFGS's default tolerance; the batched Newton fits
    are the exact maximum-likelihood ones: equal to scikit-learn run to a tight tolerance, and close to
    SnapATAC2's p-values."""
    from scipy.stats import chi2
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import log_loss
    from snapatac2.tools._diff import _likelihood_ratio_test_many

    z, M = binary_matrix(120, 150, 300, 1)
    M = sp.csc_matrix(M)
    cov = np.log1p(np.ravel(M.sum(axis=1)))
    new = R.lrt_pvalues(cov, z, M)
    assert np.isfinite(new).all()
    for i in range(0, 300, 15):
        y = (M[:, i].toarray().ravel() > 0).astype(int)
        ll = [-log_loss(y, LogisticRegression(C=np.inf, tol=1e-12, max_iter=100_000).fit(X, y).predict_proba(X),
                        normalize=False) for X in (cov[:, None], np.c_[cov, z])]
        exact = chi2.sf(max(2 * (ll[1] - ll[0]), 0), 1)
        if exact < 1e-30:                # separable: both far beyond any threshold
            assert new[i] < 1e-25
        else:
            assert new[i] == pytest.approx(exact, rel=1e-3)
    ref = np.array(_likelihood_ratio_test_many(cov[:, None], z[:, None], M.copy()))
    assert np.abs(np.log10(np.maximum(ref, 1e-300)) - np.log10(np.maximum(new, 1e-300))).max() < 0.1


def test_diff_test_matches_snapatac2():
    import anndata as ad
    import snapatac2 as snap

    z, M = binary_matrix(200, 260, 400, 2)
    names = [f"chr1:{i * 1000}-{i * 1000 + 501}" for i in range(M.shape[1])]
    a = ad.AnnData(X=M, obs=pd.DataFrame(index=[f"c{i}" for i in range(M.shape[0])]), var=pd.DataFrame(index=names))
    fg, bg = np.flatnonzero(z == 1), np.flatnonzero(z == 0)       # group 1 = the open-in-group-2 peaks' side
    feats = np.arange(0, M.shape[1], 2)
    ref = snap.tl.diff_test(a, cell_group1=list(a.obs_names[fg]), cell_group2=list(a.obs_names[bg]),
                            features=list(feats), direction="positive").to_pandas()
    new = R.diff_test(M, a.var_names, fg, bg, feats)
    assert set(new["peak"]) == set(ref["feature name"])
    m = new.merge(ref, left_on="peak", right_on="feature name")
    assert np.allclose(m["log2_fold_change"], m["log2(fold_change)"])
    lp = lambda x: -np.log10(np.maximum(x.to_numpy(dtype=float), 1e-300))  # noqa: E731
    assert np.abs(lp(m["p_value"]) - lp(m["p-value"])).max() < 0.1
    assert np.abs(lp(m["adjusted_p_value"]) - lp(m["adjusted p-value"])).max() < 0.1


@pytest.fixture(scope="module")
def genome(tmp_path_factory):
    trees = make_fixture(tmp_path_factory.mktemp("fx"), "flat")
    genome_dir = next(trees["genomes"].rglob("genome.fa.gz")).parent
    return R.plain_fasta(genome_dir, tmp_path_factory.mktemp("work")), trees["motifs"]


def test_motif_enrichment_matches_snapatac2(genome):
    import snapatac2 as snap

    fasta, motif_file = genome
    pfms = snap.read_motifs(str(motif_file))
    rng = np.random.default_rng(3)
    peaks = [f"chr{c}:{s}-{s + 501}" for c, s in zip(rng.integers(1, 3, 600), rng.integers(1000, 1_900_000, 600))]
    peaks = list(dict.fromkeys(peaks))
    regions = {"A": peaks[:80], "B": peaks[60:150]}                                    # overlapping foregrounds
    backgrounds = {"A": peaks[200:450], "B": peaks[300:600]}                           # overlapping backgrounds
    new = R.motif_enrichment(pfms, regions, backgrounds, fasta)
    for t in regions:
        ref = snap.tl.motif_enrichment(pfms, {t: regions[t]}, str(fasta), background=backgrounds[t],
                                       method="binomial")[t].to_pandas()
        assert list(new[t]["motif_id"]) == list(ref["id"])
        assert np.array_equal(np.asarray(new[t]["log2_fold_change"]), ref["log2(fold change)"].to_numpy())
        assert np.allclose(new[t]["p_value"], ref["p-value"], rtol=1e-12)
        assert np.allclose(new[t]["adjusted_p_value"], ref["adjusted p-value"], rtol=1e-12)
