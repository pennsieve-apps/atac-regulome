"""
A tiny synthetic input in the shapes atac-regulome reads, plus a tiny genome and motif file:

  python -m tests.make_fixture data/fixture [--layout linked|flat|both]

writes data/fixture/genomes/GRCh38/ (annotation.gff3.gz, chrom_sizes.tsv, genome.fa.bgz: two
short random chromosomes with 60 genes; read through GENOMES_DIR), data/fixture/motifs.meme (three
motifs; read through MOTIFS_FILE) and data/fixture/<layout>/ with

  atac/<sample_id>.h5ad   two samples as atac-prep writes them (SnapATAC2 import of synthetic
                          fragments, obs_names <sample_id>:<barcode>, uns["pennsieve"]["genome"])
  labels/labels.csv       sample_id, obs_id, cell_type (as sc-annotate writes it)
  cells/<sample_id>.h5ad  the RNA side (counts in layers["counts"], log-normalized X)
  (flat: the same files at the root)

Three cell types (TypeA, TypeB, TypeC; 100 cells per sample each) and a 10-cell RareType. Each
main type has 40 sites only it opens; TypeA's sites carry the MOTIF_A consensus, and the gene
nearest each type's sites (LINKED[type]) is expressed in that type only.
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

CHROMS = {"chr1": 3_000_000, "chr2": 2_000_000, "chrM": 16_569}
TYPES = ("TypeA", "TypeB", "TypeC")
N_PER_TYPE, N_RARE, N_SITES, N_GENES = 100, 10, 40, 60
SAMPLES = ("s1", "s2")
MOTIF_A = "TTGACCGTCAAGGTCA"      # planted in TypeA's sites (motif MA9001.1 FIXA)
MEME = """MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

"""


def pwm(consensus: str, motif_id: str, name: str) -> str:
    rows = []
    for ch in consensus:
        p = {b: 0.02 for b in "ACGT"}
        p[ch] = 0.94
        rows.append(" " + "  ".join(f"{p[b]:.6f}" for b in "ACGT"))
    return (f"MOTIF {motif_id} {name}\nletter-probability matrix: alength= 4 w= {len(consensus)} nsites= 20 E= 0\n"
            + "\n".join(rows) + "\n\n")


def layout_genome(rng):
    """Genes and type-specific sites: each type's sites cluster 5-30 kb from its linked genes."""
    genes = []
    for i in range(N_GENES):
        chrom = "chr1" if i % 3 else "chr2"
        genes.append((chrom, int(rng.integers(50_000, CHROMS[chrom] - 50_000)), "+" if i % 2 else "-", f"GENE{i:02d}"))
    genes.sort()
    linked = {t: [genes[3 * k + j] for k in range(8)] for j, t in enumerate(TYPES)}   # 8 genes per type
    sites = {}
    for t in TYPES:
        s = []
        for k in range(N_SITES):
            chrom, tss, _, _ = linked[t][k % 8]
            s.append((chrom, tss + int(rng.choice([-1, 1])) * int(rng.integers(5_000, 30_000))))
        sites[t] = s
    return genes, linked, sites


def write_genome(root: Path, genes, sites, rng) -> Path:
    from Bio import bgzf

    d = root / "genomes" / "GRCh38"
    d.mkdir(parents=True, exist_ok=True)
    (d / "chrom_sizes.tsv").write_text("".join(f"{c}\t{n}\n" for c, n in CHROMS.items()))
    with gzip.open(d / "annotation.gff3.gz", "wt") as fh:
        fh.write("##gff-version 3\n")
        for i, (chrom, tss, strand, name) in enumerate(genes):
            start, end = (tss, tss + 5000) if strand == "+" else (tss - 5000, tss)
            gid, tid = f"ENSG{i:011d}.1", f"ENST{i:011d}.1"
            common = f"gene_id={gid};gene_type=protein_coding;gene_name={name}"
            fh.write(f"{chrom}\tHAVANA\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid};{common}\n")
            fh.write(f"{chrom}\tHAVANA\ttranscript\t{start}\t{end}\t.\t{strand}\t.\tID={tid};Parent={gid};{common};"
                     f"transcript_id={tid};transcript_type=protein_coding\n")
            fh.write(f"{chrom}\tHAVANA\texon\t{start}\t{end}\t.\t{strand}\t.\tID=exon:{tid}:1;Parent={tid};{common};"
                     f"transcript_id={tid};transcript_type=protein_coding\n")
    seqs = {c: rng.choice(np.array(list("ACGT")), n) for c, n in CHROMS.items()}
    for chrom, centre in sites["TypeA"]:
        for off in (-60, 0, 60):
            p = centre + off - len(MOTIF_A) // 2
            seqs[chrom][p:p + len(MOTIF_A)] = list(MOTIF_A)
    with bgzf.BgzfWriter(str(d / "genome.fa.bgz"), "wb") as fh:
        for c, s in seqs.items():
            fh.write(f">{c}\n".encode())
            text = "".join(s)
            for i in range(0, len(text), 60):
                fh.write((text[i:i + 60] + "\n").encode())
    from pyfaidx import Fasta

    Fasta(str(d / "genome.fa.bgz")).close()      # .fai and .gzi, as the image builds them
    (root / "motifs.meme").write_text(MEME + pwm(MOTIF_A, "MA9001.1", "FIXA") + pwm("CCGCGGACGTACGCGG", "MA9002.1", "FIXB")
                                      + pwm("AATTTAAATTTAAATT", "MA9003.1", "FIXC"))
    return d.parent


def cell_fragments(bc: str, n: int, own_sites, genes, rng) -> pd.DataFrame:
    """30 % at the cell type's own sites, 20 % at promoters, 3 % chrM, the rest uniform."""
    kind = rng.random(n)
    length = np.where(rng.random(n) < 0.7, rng.integers(40, 150, n), rng.integers(180, 260, n))
    s = rng.integers(len(own_sites), size=n) if own_sites else np.zeros(n, int)
    g = rng.integers(len(genes), size=n)
    site_chrom = np.array([x[0] for x in own_sites] or ["chr1"])[s]
    site_pos = np.array([x[1] for x in own_sites] or [0])[s]
    gene_chrom, gene_pos = np.array([x[0] for x in genes])[g], np.array([x[1] for x in genes])[g]
    uni_chrom = np.where(rng.random(n) < 0.6, "chr1", "chr2")
    uni_pos = (rng.random(n) * (np.where(uni_chrom == "chr1", CHROMS["chr1"], CHROMS["chr2"]) - 500)).astype(int)
    site = (kind < 0.30) & bool(own_sites)
    prom = (kind >= 0.30) & (kind < 0.50)
    mito = (kind >= 0.50) & (kind < 0.53)
    chrom = np.where(site, site_chrom, np.where(prom, gene_chrom, np.where(mito, "chrM", uni_chrom)))
    centre = np.where(site, site_pos, np.where(prom, gene_pos, np.where(mito, rng.integers(500, 16_000, n), uni_pos)))
    start = np.maximum(0, centre + rng.normal(0, 100, n).astype(int) - length // 2)
    return pd.DataFrame({"chrom": chrom, "start": start, "end": start + length, "barcode": bc, "count": 1})


def barcodes(n: int, rng) -> list[str]:
    seen: set[str] = set()
    while len(seen) < n:
        seen.add("".join(rng.choice(list("ACGT"), 16)) + "-1")
    return sorted(seen)


def write_atac(path: Path, frags: pd.DataFrame, sid: str, genome_dir: Path, tmp: Path) -> None:
    """The sample as atac-prep writes it: SnapATAC2 import, TSSE, 500 bp tiles, prefixed obs_names."""
    import polars as pl
    import snapatac2 as snap

    order = {c: i for i, c in enumerate(CHROMS)}
    frags = frags.assign(o=frags["chrom"].map(order)).sort_values(["o", "start", "end"]).drop(columns="o")
    fpath = tmp / f"{sid}_fragments.tsv.gz"
    with gzip.open(fpath, "wt") as fh:
        frags.to_csv(fh, sep="\t", header=False, index=False)
    a = snap.pp.import_fragments(fpath, chrom_sizes=CHROMS, file=path, sorted_by_barcode=False, min_num_fragments=0,
                                 n_jobs=1)
    snap.metrics.tsse(a, genome_dir / "GRCh38" / "annotation.gff3.gz")
    snap.pp.add_tile_matrix(a, bin_size=500, n_jobs=1)
    a.obs_names = [f"{sid}:{b}" for b in a.obs_names]
    a.obs["sample_id"] = pl.Series([sid] * a.n_obs, dtype=pl.Categorical)
    a.uns["pennsieve"] = {"genome": "GRCh38", "sample_id": sid, "cells_from": "rna", "organism": "human"}
    a.close()


def write_rna(path: Path, sid: str, cells: list[tuple[str, str]], genes, linked, rng) -> None:
    import anndata as ad
    import scipy.sparse as sp

    names = [g[3] for g in genes]
    lam = np.full((len(cells), len(genes)), 1.0)
    for i, (_, t) in enumerate(cells):
        for g in linked.get(t, []):
            lam[i, names.index(g[3])] = 12.0
    counts = sp.csr_matrix(rng.poisson(lam).astype(np.float32))
    tot = np.asarray(counts.sum(axis=1)).ravel()
    X = sp.csr_matrix(counts.multiply(1e4 / np.maximum(tot, 1)[:, None]))
    X.data = np.log1p(X.data)
    a = ad.AnnData(X=X.astype(np.float32), obs=pd.DataFrame(index=[f"{sid}:{b}" for b, _ in cells]),
                   var=pd.DataFrame(index=names))
    a.layers["counts"] = counts
    a.obs["sample_id"] = pd.Categorical([sid] * len(cells))
    a.uns["pennsieve"] = {"organism": "human"}
    a.write_h5ad(path)


TRUTH: dict = {}


def make_fixture(out: Path, layout: str = "both", seed: int = 0) -> dict[str, Path]:
    rng = np.random.default_rng(seed)
    genes, linked, sites = layout_genome(rng)
    genome_dir = write_genome(out, genes, sites, rng)
    TRUTH.update(linked={t: [g[3] for g in linked[t]] for t in TYPES}, sites=sites)
    tmp = out / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    built: dict[str, dict] = {}
    labels = []
    for sid in SAMPLES:
        bcs = barcodes(N_PER_TYPE * len(TYPES) + N_RARE, rng)
        cells = [(b, TYPES[i // N_PER_TYPE] if i < N_PER_TYPE * len(TYPES) else "RareType") for i, b in enumerate(bcs)]
        frags = pd.concat([cell_fragments(b, int(rng.integers(1500, 2000)), sites.get(t, []), genes, rng)
                           for b, t in cells], ignore_index=True)
        atac = tmp / f"{sid}.h5ad"
        write_atac(atac, frags, sid, genome_dir, tmp)
        rna = tmp / f"{sid}.rna.h5ad"
        write_rna(rna, sid, cells, genes, linked, rng)
        built[sid] = {"atac": atac, "rna": rna}
        labels += [{"sample_id": sid, "obs_id": f"{sid}:{b}", "cell_type": t, "cell_type_score": 0.9} for b, t in cells]
    trees = {"genomes": genome_dir, "motifs": out / "motifs.meme"}
    for lay in (("linked", "flat") if layout == "both" else (layout,)):
        root = out / lay
        sub = (lambda p: root / p) if lay == "linked" else (lambda p: root)
        for d in ("atac", "labels", "cells"):
            sub(d).mkdir(parents=True, exist_ok=True)
        for sid, f in built.items():
            (sub("atac") / f"{sid}.h5ad").write_bytes(f["atac"].read_bytes())
            (sub("cells") / (f"{sid}.h5ad" if lay == "linked" else f"{sid}_rna.h5ad")).write_bytes(f["rna"].read_bytes())
        pd.DataFrame(labels).to_csv(sub("labels") / "labels.csv", index=False)
        trees[lay] = root
    return trees


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path)
    ap.add_argument("--layout", choices=("linked", "flat", "both"), default="both")
    args = ap.parse_args()
    for name, path in make_fixture(args.out, args.layout).items():
        print(f"{name}: {path}")
