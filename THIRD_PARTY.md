# Third-party software and data

This app is Apache-2.0 (see `LICENSE`). The image contains the following third-party components,
installed from PyPI at the versions pinned in `requirements.txt` (MACS3 built from its pinned
source release), and reference data downloaded at build time; their licences apply to the image.

**GPL notice:** the image contains python-igraph 1.0.0 (GPL-2.0-or-later), which SnapATAC2's
`aggregate_cells` uses for the Leiden clustering behind the pseudo-cells. Our code is Apache-2.0 and
calls igraph through SnapATAC2; whoever redistributes the image must meet the GPL's terms for
igraph (its source is at https://github.com/igraph/python-igraph, tag 1.0.0).

| Component | What | Licence | Source |
|---|---|---|---|
| SnapATAC2 2.10.0 | `AnnDataSet`, `tl.macs3`, `tl.merge_peaks`, `pp.make_peak_matrix`, `tl.diff_test`, `tl.motif_enrichment`, `read_motifs`, `tl.init_network_from_annotation`, `tl.aggregate_cells`, `tl.aggregate_X` | MIT | https://github.com/scverse/SnapATAC2 |
| MACS3 3.0.5 | Peak calling (through SnapATAC2) | BSD-3-Clause | https://github.com/macs3-project/MACS |
| python-igraph 1.0.0, texttable 1.7.2 | Leiden clustering for pseudo-cells | GPL-2.0-or-later; MIT | https://igraph.org/python |
| pyfaidx 0.8.2 | Random access to the genome FASTA (decompressed to `/tmp` for the motif step) | BSD-3-Clause | PyPI |
| hmmlearn 0.3.3 | Dependency of MACS3 | BSD-3-Clause | PyPI |
| multiprocess 0.70.19, dill 0.4.1 | SnapATAC2 runs MACS3 in worker processes | BSD-3-Clause | PyPI |
| anndata 0.13.4, h5py 3.16.0 | RNA `.h5ad`, peak matrix | BSD-3-Clause | https://github.com/scverse/anndata |
| polars 1.44.2, pyarrow 26.0.0 | SnapATAC2's data frames | MIT; Apache-2.0 | PyPI |
| NumPy 2.5.4, SciPy 1.18.1, pandas 2.3.3, scikit-learn 1.9.1 | Arrays, statistics, logistic regression (diff_test), truncated SVD | BSD-3-Clause | PyPI |
| Matplotlib 3.11.2 and its dependencies | PNG figures | Matplotlib licence (PSF-based) and BSD/MIT | https://matplotlib.org |
| rustworkx 0.18.1, pooch, scverse-misc, typeguard, natsort, tqdm and the rest of `requirements.txt` | Dependencies of the above | Apache-2.0 / BSD / MIT / MPL-2.0 | PyPI |
| python:3.12.16-slim-bookworm | Base image | PSF; Debian packages under their own free licences | https://hub.docker.com/_/python |

Build stage only (not in the final image): gcc 12.2.0, libc6-dev, zlib1g-dev, setuptools 84.2.1,
Cython 3.3.0 (to build the MACS3 wheel).

## Data in the image

| File | What | Licence / terms | Source, checksum |
|---|---|---|---|
| `/opt/motifs/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt` | JASPAR 2026 CORE vertebrates, non-redundant, 1,019 PFMs | CC BY 4.0 (stated on jaspar.elixir.no, checked 2026-10-31) | jaspar.elixir.no/download/data/2026/CORE/, sha256 ae4f77a1… |
| `/opt/genomes/GRCh38/genome.fa.gz` (+ `genome.fa.fai`) | GENCODE GRCh38 primary assembly (release 41), chr1-22, X, Y, M, recompressed | GENCODE / Genome Reference Consortium data, freely available (EMBL-EBI terms of use) | …/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz, sha256 4fac949d… (SnapATAC2's pin) |
| `/opt/genomes/GRCm39/genome.fa.gz` (+ `genome.fa.fai`) | GENCODE GRCm39 primary assembly (release M30), chr1-19, X, Y, M | as above | …/Gencode_mouse/release_M30/GRCm39.primary_assembly.genome.fa.gz, sha256 3b923c06… |
| `/opt/genomes/*/annotation.gff3.gz` | GENCODE basic annotation (v41 / vM30) | as above | sha256 b82a655b… / 6f433e26… |
| `/opt/genomes/*/chrom_sizes.tsv` | Primary chromosome lengths | SnapATAC2's tables (MIT) | written at build time |

## Not bundled

Demo data: Yang L, Gereau RW, *A single-cell multi-omic atlas of human sympathetic ganglia*,
Pennsieve Discover dataset 465 v1, doi:10.26275/oqxx-si9h (CC BY 4.0).
