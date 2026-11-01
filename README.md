# atac-regulome

Gene regulation per cell type from single-nucleus ATAC, the second step of the **Multiome Gene
Regulation** workflow:

- **B. Multiome Gene Regulation:** `atac-prep → atac-regulome (+ labels from sc-annotate) → sc-report`

Genes are switched on and off by stretches of DNA (promoters, enhancers) that are open in the cell
types that use them. Knowing which regions are open in, say, satellite glia but not in neurons
points to the regulatory elements, and the transcription factors binding them, that make each cell
type what it is, and to where disease-associated DNA variants may act. This app takes the ATAC half
of a multiome dataset (from atac-prep) and cell-type labels for the same nuclei (from sc-annotate,
or any labels), and for every cell type:

1. pools its nuclei and finds the open regions (**peaks**) with MACS3;
2. merges all cell types' peaks into one peak set and counts every nucleus in it;
3. finds **marker peaks**, regions more open in this cell type than in the others;
4. tests which transcription-factor **motifs** (JASPAR 2026) are over-represented in those marker
   peaks, hinting at the factors that drive the cell type;
5. when the RNA half is linked, **links peaks to genes**: a peak whose openness rises and falls
   with a nearby gene's expression across groups of similar nuclei is a candidate regulator of it.

For the computational analyst: SnapATAC2 end to end. A backed `AnnDataSet` over the samples (no
copy); `tl.macs3` per label (pseudobulk of the label's Tn5 insertions, q < 0.05; ENCODE exclusion
list v2 regions removed for GRCh38; chrM peaks dropped);
`tl.merge_peaks` (summit ± 250 bp, overlaps resolved by significance); `pp.make_peak_matrix`;
`tl.diff_test`'s test (logistic-regression likelihood-ratio test, one label against an equal-share
background of the others, over the peaks called in the label, BH FDR < 0.01) as in SnapATAC2's
differential-accessibility tutorial, fitted for all peaks at once (see Runtime); `tl.motif_enrichment`'s test
(binomial, each label's marker peaks against its own GC-matched background drawn from the merged peak
set; every distinct sequence scanned once for all labels) with the JASPAR 2026 CORE vertebrates
non-redundant PFMs and sequences from the GENCODE genome built into the image;
`tl.init_network_from_annotation` (peaks within 250 kb of protein-coding TSSs) and Spearman
correlation across `tl.aggregate_cells` pseudo-cells of RPM accessibility and CPM expression.

## Method details

- **Labels.** `labels.csv` is joined on `obs_id` = the ATAC `obs_names` (`<sample_id>:<barcode>`,
  the same names on both sides of a multiome run). Unlabelled nuclei and labels with fewer than
  `MIN_CELLS` nuclei are left out (listed in `summary.json`). sc-annotate's `labels.csv` works as is;
  the demo used ATAC clusters (see Results).
- **Peaks.** MACS3 runs on each label's pooled insertions (SnapATAC2 shifts −100 bp and extends 200 bp,
  `nolambda` off, summits called, fold enrichment ≥ 1). Labels are pooled across samples; peaks are
  not required to replicate per sample (ArchR's pseudo-replicate scheme is not used).
  `peaks/<label>.bed` holds the merged-set peaks that overlap the label's own calls.
- **Marker peaks.** Per label at most `MAX_CELLS_PER_TEST` nuclei (seeded sample) against as many
  background nuclei drawn equally from the other labels (the SnapATAC2 tutorial draws 30 per other
  label; equal shares keep large labels from dominating). Tested peaks: those called in the label,
  detected in ≥ 5 % of either group, log2 fold change ≥ 0.25 (SnapATAC2's pre-filter), not open in
  every tested nucleus. `nearest_gene` / `distance_to_tss`: the closest protein-coding gene within
  `LINK_WINDOW_KB` (empty when none).
- **Motifs.** Labels with ≥ 10 marker peaks are tested, each against its own background drawn from
  the merged peaks that are no label's marker peak, matched to that label's marker GC content
  (20 bins, three times as many per bin as the markers have): marker peaks of some cell types are GC-rich
  promoters and CpG islands more often than average, and a GC-blind background would make every
  GC-rich motif (SP/KLF, E2F, ZNF) look enriched. `summary.json → motif_background` gives, per
  label, the background size and the mean GC of the markers and of the background. A motif open in every cell type (e.g. CTCF) shows as
  enriched in all labels against this background.
  `motifs.csv` has every motif × label; the summary and figure count `adjusted_p_value < 0.05`
  with positive log2 fold change.
- **Links** (RNA linked): nuclei present in both modalities; pseudo-cells from
  `snap.tl.aggregate_cells` (iterative Leiden, ≥ `PSEUDOCELL_SIZE` nuclei each) on the RNA
  embedding (`obsm["X_integrated"]` from sc-integrate when every RNA file has it; else a
  30-component truncated SVD of log CPM over all matched nuclei, because per-sample `X_pca` from
  sc-qc are separate spaces); per pseudo-cell, peak RPM and gene CPM (from `layers["counts"]`, else `X`);
  Spearman ρ for every peak–gene pair of the network whose peak and gene are both detected in
  ≥ 20 % (and ≥ 5) of the pseudo-cells (ties otherwise make ρ meaningless); p from the t approximation with n − 2
  degrees of freedom, BH-adjusted over all pairs. Kept: ρ ≥ `LINK_MIN_COR` and FDR < 0.05
  (positive links only, as ArchR's peak-to-gene links). These defaults are looser than ArchR's
  (ρ > 0.45, FDR < 1e-4, over ~500 metacells of 100 nuclei): with tens of pseudo-cells, ArchR's cut
  keeps almost nothing; raise `LINK_MIN_COR` on large datasets. SnapATAC2 2.10.0's `add_cor_scores`
  misaligns columns when every peak of the matrix is in the network (it then takes the matrix in
  storage order but pairs it with the network's node order), so the same correlation is computed
  here directly.
- **Two SnapATAC2 2.10.0 workarounds**, kept small and commented in `processor/regulome.py`:
  `tl.macs3(selections=…)` returns nothing (it compares label names with internal indices), so
  unlabelled nuclei form one extra group whose peaks are discarded; with one worker `tl.macs3`
  leaves Python's temp directory pointing at a deleted folder, which is reset.

Deviations from the best-practices guide (Heumos et al. 2023): it recommends peak calling per cell
type with MACS2 and pseudo-replicates, chromVAR for motif activity and ArchR-style peak-to-gene
links. This app calls peaks per type with MACS3 (MACS2's successor) without replicate filtering,
tests motif over-representation in marker peaks instead of per-cell chromVAR deviations (one
table per cell type is what the report shows), and builds links on SnapATAC2's pseudo-cells.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `atac` | input (subfolder) | `<sample_id>.h5ad` from atac-prep | All samples must name the same genome in `uns["pennsieve"]["genome"]`. |
| `labels` | input (subfolder) | `labels.csv` (`obs_id` + `LABEL_COLUMN`) | From sc-annotate (`cell_type`) or sc-integrate (`cluster`). |
| `cells` | input (subfolder, optional) | RNA `<sample_id>.h5ad` (sc-qc / sc-integrate / sc-annotate) | Only for peak-to-gene links. |
| `regulome` | output | `<LABEL>/…` | Below. Unlinked, everything may sit side by side at the input root (ATAC files are recognised by their fragments). |

## Outputs (`OUTPUT_DIR/<LABEL>/`)

| File | Content |
|---|---|
| `peaks.bed` | The merged peak set: `chrom, start, end, name` (`chrom:start-end`), 501 bp each, GRCh38/GRCm39 coordinates (0-based starts). |
| `peaks/<label>.bed` | The merged peaks called in that label (file names: label with characters outside `A-Za-z0-9_.-` replaced by `_`; mapping in `summary.json`). |
| `marker_peaks.csv` | `label, peak, chrom, start, end, log2_fold_change, p_value, adjusted_p_value, n_cells, n_background, nearest_gene, distance_to_tss` (FDR < `MARKER_FDR`). |
| `motifs.csv` | `label, motif_id, motif_name, log2_fold_change, p_value, adjusted_p_value` (JASPAR id and TF name; every motif × label tested). |
| `links.csv` | RNA linked only: `peak, chrom, start, end, gene, distance_to_tss, correlation, p_value, adjusted_p_value, called_in` (labels whose peak calls contain the peak). |
| `peaks_per_label.png` | Peaks called and marker peaks per label. |
| `marker_peaks_heatmap.png` | z-scored pseudobulk accessibility of each label's top 20 marker peaks across labels. |
| `motifs_heatmap.png` | −log10 adjusted p of each label's top 5 enriched motifs across labels. |
| `links.png` | Kept links: correlation against distance to TSS. |
| `summary.json` | Cells used and left out, peaks, per label: cells, peaks, marker peaks, top motifs, links, BED file; link statistics. |
| `provenance.json` | Version, tools, parameters, inputs, timings, peak memory, citations. |

## Parameters

| Parameter | Default | What it does |
|---|---|---|
| `LABEL_COLUMN` | `cell_type` | Column of `labels.csv` to group by. |
| `LABEL` | (empty) | Output folder; empty = `LABEL_COLUMN` with `_` → `-` (`cell-type`). |
| `MIN_CELLS` | `50` | Labels with fewer nuclei are left out. |
| `PEAK_QVALUE` | `0.05` | MACS3 q-value (SnapATAC2's default). |
| `MAX_CELLS_PER_TEST` | `500` | Nuclei per label (and background) in a marker test. |
| `MARKER_FDR` | `0.01` | Marker peak FDR (SnapATAC2 tutorial). |
| `LINK_WINDOW_KB` | `250` | Peak–TSS distance for links and nearest genes (SnapATAC2's default, as ArchR). |
| `LINK_MIN_COR` | `0.3` | Minimum Spearman ρ of a kept link. |
| `PSEUDOCELL_SIZE` | `50` | Minimum nuclei per pseudo-cell (SnapATAC2's default). |
| `SEED` | `42` | Marker-test sampling, pseudo-cell clustering. |
| `THREADS` | `0` | Threads and MACS3 worker processes (each worker needs ~1 GB, see Runtime); 0 = CPU quota. |

## Results on the demo data

**Labels.** Discover 465's `meta.data.csv` has the authors' cell types, but they do not line up with
its cell names: per-cell counts match the matrix (nCount_RNA r = 0.9997), while the "sample id"
column is filled in blocks and the cell types show no marker-gene specificity in the RNA (C1QA,
CD163, MPZ, PECAM1, DCN are detected equally in every "type"). So the demo labels are SnapATAC2
Leiden clusters of the ATAC data (`scripts/demo_labels.py`, `LABEL_COLUMN=cluster`); RNA markers
name them: C01, C02 SOX10⁺ glia (satellite glia / non-myelinating Schwann cells), C03 macrophages
(CD163, F13A1, MRC1, PTPRC), C04 fibroblasts and vascular cells (DCN, PDGFRB; PECAM1, VWF), C05
myelinating Schwann cells (MPZ, PRX), C06 T cells (PTPRC, CD3E).

**Input.** atac-prep's output for the two-library cut of Discover 465, restricted to chr1–chr3
(so the run takes minutes under emulation; see Runtime for the whole genome): 2,987 nuclei
(1,609 + 1,378; fewer than on the whole genome because the per-nucleus fragment counts on three
chromosomes fall below 1,000 for some), the RNA side linked. Outputs in
`workflow-single-cell/data/work/sympathetic/atac-regulome/regulome/cluster/`.

| Label | Nuclei | MACS3 peaks | Marker peaks (FDR < 0.01) | Top JASPAR motifs in marker peaks (GC-matched background) |
|---|---|---|---|---|
| C01 glia | 1,016 | 24,578 | 879 | SOX10, NFIA, NFIX, NR3C2, ZNF281 |
| C02 glia | 905 | 32,152 | 995 | SOX10, ZNF281, PGR, ZNF148, NFIA |
| C03 macrophages | 400 | 28,344 | 2,117 | SPIB, Spi1 (PU.1), IKZF2, EHF, ELF3 |
| C04 fibroblast / vascular | 331 | 32,272 | 1,770 | BHLHE22, NEUROD2, OLIG2, TCF12, TWIST2 (E-box bHLH) |
| C05 myelinating Schwann | 276 | 19,240 | 107 | HNF1B, SRF, E2F7, GMEB2, HES6 |
| C06 T cells | 55 | 5,948 | 793 | HOXA11, HOXC12, HOXD11, TFCP2, GABPA |

Merged set (ENCODE exclusion list removed): 78,877 peaks; 863 motif × label enrichments at
FDR < 0.05. The well-populated labels get the expected lineage factors (SOX10 and NFI in glia,
PU.1/SPI and ETS in macrophages, E-box bHLH incl. TWIST2 in mesenchymal cells); the two small ones
are noisy (C05 has only 107 marker peaks, C06 only 55 nuclei), which is what to expect from labels
that small. The backgrounds match each label's GC content (e.g. C05 markers 0.597, background
0.594). Peak-to-gene links: 26 pseudo-cells
from 2,983 matched nuclei, 150,030 peak–gene pairs tested (both detected in ≥ 20 % of
pseudo-cells), 28 kept (ρ ≥ 0.3, FDR < 0.05) — few, as expected with 26 pseudo-cells; a full
dataset (tens of thousands of nuclei, hundreds of pseudo-cells) has far more power. Before the
detection filter, genes seen in a single pseudo-cell gave spurious ρ ≈ 1 links; the filter removed
them.

The same nuclei on the **whole genome** (run stopped in the motif step before the FASTA fix; peaks
and markers completed): 156,721 merged peaks; marker peaks C01 2,639, C02 3,119, C03 7,488,
C04 5,700, C05 830, C06 2,716.

## Runtime

`computeTypes: [standard]`, 8 vCPU / 32 GiB, 12 h timeout. **The 32 GiB and 12 h are estimates**,
extrapolated from the chr1–3 cut below to a 40,000-nucleus dataset; a platform run on full data will
measure them (a whole-genome run cannot finish under emulation in the 8 GB of Docker here). Measured in the linux/amd64 image under
emulation on an M-series Mac (Docker `--user 1000`, `--memory 6g --cpus 4`, `THREADS=2`, `/tmp`
on the container's disk):

| Input | Wall time | Steps | Peak RSS | cgroup peak |
|---|---|---|---|---|
| chr1–3 cut: 2,987 nuclei, 79 k peaks, RNA linked | 669 s | peaks 48 s, markers 120 s, motifs 494 s (decompressing GRCh38, GC of all peaks, 6 labels × markers + 3× GC-matched background × 1,019 motifs), links 3 s | 1.45 GB | 3.44 GB (incl. page cache of the 3.1 GB decompressed genome) |
| whole-genome cut: 4,001 nuclei, 157 k peaks | peaks 226 s; marker tests ≈ 55 min (container clock) | — | — | killed at 6 GB with `THREADS=4` during MACS3 |

What drives memory and time on full data (estimates, extrapolated from these runs):
- **MACS3** runs one worker process per label in parallel (up to `THREADS`); measured natively on
  the whole-genome cut: the parent 2.8 GB while exporting insertions, each worker ~0.9 GB for
  ~500 nuclei. Workers scale with a label's fragments, so a 20,000-nucleus label with ~10,000
  fragments each is ~20× this cut's largest label: several GB per worker. `THREADS` caps how many
  run at once.
- **Peak matrix**: in memory, ~12 bytes per non-zero (CSR); 40,000 nuclei × ~5,000 open peaks each
  ≈ 2.4 GB. It is counted for every nucleus of the samples and then subset to the labelled ones
  (counting only those would need a copy of the fragments).
- **Marker tests.** SnapATAC2 2.10.0's `tl.diff_test` fits two scikit-learn logistic regressions per
  peak, one after another: on the platform (v0.1.0, 4 vCPU) 54,000 peaks of 3 labels took 22 min. Since
  v0.1.1 the same two models (null: intercept + log1p total count; alternative: + group) are fitted by
  Newton's method for a block of 4,096 peaks at a time with matrix products (`regulome.lrt_pvalues`):
  11× faster on one core (8,000 peaks: 23.3 s → 2.1 s), more with BLAS threads. They are the exact
  maximum-likelihood fits; SnapATAC2's L-BFGS fits stop at scikit-learn's default tolerance, so its
  p-values differ slightly (within 0.1 on the log10 scale; tests/test_equivalence.py checks both).
  The peak filter (`min_pct` 0.05, `min_log_fc` 0.25, RPM with pseudo-count 1) and BH are SnapATAC2's.
- **Motifs**: the genome is decompressed to `/tmp` (3.1 GB GRCh38, 2.7 GB GRCm39); scanning is
  ~14 µs per sequence × motif natively on one core (measured: 1.4 s for 100 motifs × 1,000 sequences)
  and SnapATAC2 spreads it over the CPUs. Each label's markers are compared with three times as many
  GC-matched background peaks; `tl.motif_enrichment` takes one background per call, so v0.1.0 scanned
  the overlapping backgrounds once per label (on the platform, 3 labels: 114,600 sequences, 29.5 min
  at 4 vCPU). Since v0.1.1 every distinct sequence is scanned once for all labels
  (`regulome.motif_enrichment`, the same counts and binomial test, checked against SnapATAC2).
- **Disk** in `/tmp`: the decompressed genome plus MACS3's exported insertions of the labelled nuclei.

32 GiB and 8 vCPU leave room for ~10 labels of a 40,000-nucleus multiome dataset at `THREADS=4`;
lower `THREADS` before raising memory if MACS3 runs out.

Image: about 2.8 GB of uncompressed layers (`docker history`: genomes 1.74 GB, Python packages
861 MB, base ~150 MB); Docker Desktop's `docker images` shows 4.79 GB because its containerd store
counts the compressed and the unpacked content together.

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python       # pytest on the synthetic fixture (no Docker)
make run PYTHON=.venv/bin/python        # build the image, run it over data/fixture/linked -> data/output
make validate PYTHON=.venv/bin/python   # app.yml against the manifest schema, CITATION.cff
```

The fixture (`tests/make_fixture.py`, `sample-input/README.md`) plants type-specific sites, a motif
in one type's sites and type-specific genes next to them; the tests check that the marker peaks are
those sites, that the planted motif is that type's top motif and that the links recover the planted
genes (≥ 90 % precision), in both input layouts, plus parameter and label errors.

## Next: pain-GWAS SNP overlap (prepared, not built)

The per-label BED files are the input such an analysis needs; nothing else has to change here.

1. **Coordinates.** `peaks.bed` and `peaks/<label>.bed` are GRCh38 (GENCODE primary assembly,
   `chr`-prefixed, 0-based starts). GWAS summary statistics on GRCh37/hg19 must be lifted over first
   (UCSC liftOver or CrossMap with the hg19ToHg38 chain).
2. **Simple overlap** (lead SNPs and their LD proxies, e.g. r² ≥ 0.8 in the matching 1000 Genomes
   population): `bedtools intersect -u -a snps.bed -b peaks/<label>.bed` per label, enrichment
   against the merged set or matched random regions (`bedtools shuffle -incl peaks.bed`), or GREGOR.
   Report per label the fraction of credible SNPs in its peaks; `marker_peaks.csv` and `links.csv`
   then name the nearby and linked genes of every hit.
3. **Heritability enrichment** (the stronger test): stratified LD score regression with each
   label's peaks as an annotation, conditioned on the baseline model and on the merged peak set
   (`peaks.bed`) so that a cell type wins only for its specific regions — the LDSC-SEG design of
   Finucane et al. (2018). Peaks are usually extended by ±500 bp (`bedtools slop`) to capture flanks.
4. **Which GWAS.** Chronic pain phenotypes (multisite chronic pain, back, knee, neck/shoulder,
   headache, osteoarthritis pain) from UK Biobank-scale studies or the consortia's own; the
   choice and their licences belong to that future step, and so do the references.
5. **Where.** A small light step (`atac-gwas`, Lambda-capable for simple overlap; S-LDSC needs the
   LD scores and baseline model baked in, ~10 GB) taking `regulome` + a GWAS port, or a section of
   sc-report.

## Limitations

- **Accessibility is not activity.** An open region in a cell type may be poised, repressed by
  other means, or regulate a gene outside the window. Links are correlations across pseudo-cells of
  the same dataset, not proof of regulation; they are strongest when the labels span very different
  cell types and weakest within one type.
- **Motif enrichment is not binding.** Families of factors share motifs (e.g. SOX, FOX, AP-1); a
  hit names the family more reliably than the member. Expression of the factor (RNA) narrows it.
- **Rare labels.** A label with few nuclei has few fragments: fewer peaks, few or no marker peaks
  (the demo's 55-nucleus T-cell cluster, C06, has a fifth of the peaks of the others). `MIN_CELLS` removes the worst cases; merge sub-types
  for a first look.
- **Pooled samples.** Peaks and markers pool all samples of a label; sample or donor effects are
  not modelled, and conditions are not compared (that is st-compare's job for RNA).
- **Blacklist for human only.** Peaks in the ENCODE exclusion list v2 (Amemiya et al. 2019) are
  removed for GRCh38; no such list exists for GRCm39, so mouse peaks in artefact-prone regions
  (satellite repeats, rDNA) are kept (logged as a warning). chrM peaks are dropped for both.
- **Genomes.** GRCh38 and GRCm39 only (as atac-prep).

## Security

Inputs are read through `processor/safe_inputs.py`, vendored verbatim from
`workflow-single-cell/shared/` (do not edit the copy): only files that resolve inside `INPUT_DIR` or
the platform's links are used, `labels.csv` is size-capped, `.h5ad` files with external HDF5 links
or storage are refused. The container runs as uid 1000 with a read-only root and writes only to
`OUTPUT_DIR` and `/tmp`; genomes and motifs are read from the image.

## Licences

Apache-2.0 for this app (`LICENSE`). The image bundles **python-igraph (GPL-2.0-or-later)**, used
through SnapATAC2 for the pseudo-cell clustering; see `THIRD_PARTY.md` for it and every other
component and data file (JASPAR 2026: CC BY 4.0; MACS3: BSD-3-Clause; SnapATAC2: MIT).

## Credits

- **SnapATAC2** (Kai Zhang, Bing Ren and colleagues; scverse): the whole analysis — peak calling
  wrapper, merging, peak matrix, differential tests, motif scanning, gene network and pseudo-cells.
- **MACS / MACS3** (Yong Zhang, Tao Liu, X. Shirley Liu and colleagues; the MACS3 project): peak
  calling on each label's pseudobulk.
- **JASPAR** (the JASPAR consortium; 2026 release): the transcription-factor motifs.
- **ENCODE Blacklist** (Haley Amemiya, Anshul Kundaje, Alan Boyle; ENCODE DCC): the GRCh38
  exclusion list removed from the peaks.
- **GENCODE** (the GENCODE consortium): genome sequence and gene annotation for TSSs and links.
- **anndata** (scverse): reading the RNA side.
- **Single-cell best practices** (Heumos et al.): the reference for the analysis choices above.
- **LDSC-SEG** (Hilary Finucane, Alkes Price and colleagues): the design proposed for the GWAS step.
- **Demo data:** Lite Yang and Robert W. Gereau, *A single-cell multi-omic atlas of human
  sympathetic ganglia*, SPARC / Pennsieve Discover 465.

## Citations

> Zhang K, Zemke NR, Armand EJ, Ren B. **A fast, scalable and versatile tool for analysis of
> single-cell omics data.** Nat Methods 21(2):217-227 (2024). https://doi.org/10.1038/s41592-023-02139-9,
> PMID 38191932.

> Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M,
> Li W, Liu XS. **Model-based Analysis of ChIP-Seq (MACS).** Genome Biol 9(9):R137 (2008).
> https://doi.org/10.1186/gb-2008-9-9-r137, PMID 18798982.

> Ovek Baydar D, Rauluseviciute I, Aronsen DR, Blanc-Mathieu R, et al. (34 authors). **JASPAR 2026:
> expansion of transcription factor binding profiles and integration of deep learning models.**
> Nucleic Acids Res 54(D1):D184-D193 (2026). https://doi.org/10.1093/nar/gkaf1209, PMID 41325984.

> Amemiya HM, Kundaje A, Boyle AP. **The ENCODE Blacklist: Identification of Problematic Regions of
> the Genome.** Sci Rep 9(1):9354 (2019). https://doi.org/10.1038/s41598-019-45839-z, PMID 31249361.

> Frankish A, Carbonell-Sala S, Diekhans M, et al. (62 authors). **GENCODE: reference annotation for
> the human and mouse genomes in 2023.** Nucleic Acids Res 51(D1):D942-D949 (2023).
> https://doi.org/10.1093/nar/gkac1071, PMID 36420896.

> Heumos L, Schaar AC, Lance C, et al. (43 authors incl. the Single-cell Best Practices
> Consortium). **Best practices for single-cell analysis across modalities.** Nat Rev Genet
> 24(8):550-572 (2023). https://doi.org/10.1038/s41576-023-00586-w, PMID 37002403.

> Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. **anndata: Access and store annotated data
> matrices.** J Open Source Softw 9(101):4371 (2024). https://doi.org/10.21105/joss.04371 (not in
> PubMed; checked on Crossref).

> Finucane HK, Reshef YA, Anttila V, Slowikowski K, Gusev A, et al. (22 authors). **Heritability
> enrichment of specifically expressed genes identifies disease-relevant tissues and cell types.**
> Nat Genet 50(4):621-629 (2018). https://doi.org/10.1038/s41588-018-0081-4, PMID 29632380.

Demo data:

> Yang L, Gereau R. **A single-cell multi-omic atlas of human sympathetic ganglia** (Version 1)
> [Data set]. SPARC Portal (2025). https://doi.org/10.26275/oqxx-si9h, CC BY 4.0 (checked on
> DataCite and the Pennsieve Discover API, dataset 465). Probably described by: Yang L, Dourson AJ,
> Tao R, et al. (18 authors), Gereau RW. **A multi-omic atlas of human autonomic and sensory ganglia
> implicates cell types in peripheral neuropathies.** bioRxiv (2025).
> https://doi.org/10.1101/2025.09.18.677119, PMID 41000833 (the dataset record does not link it).

> Yang L, Meriau P, Cavalli V, Gereau RW. **Nuclear extraction using opti-prep gradient
> centrifugation for 10X Genomics Multiome Assay v1.** protocols.io (2024).
> https://doi.org/10.17504/protocols.io.5jyl8q4y8l2w/v1 (the dataset's protocol; checked on Crossref).

JASPAR 2026 is the current release (the 2024 release, doi:10.1093/nar/gkad1059, PMID 37962376, was
the one first planned). Every DOI was checked on Crossref (DataCite for the
dataset) and every PMID on PubMed; full author lists are in `CITATION.cff`, generated from the
Crossref records. JASPAR's CC BY 4.0 licence was read on jaspar.elixir.no.
