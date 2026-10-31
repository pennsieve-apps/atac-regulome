# sample-input

atac-regulome reads atac-prep's output (`atac/<sample_id>.h5ad`), a `labels.csv` (`obs_id` =
`<sample_id>:<barcode>` plus the label column, as sc-annotate writes it) and optionally the RNA
side (`cells/`). Contract: `workflow-single-cell/docs/conventions.md`.

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/{genomes,motifs.meme,linked,flat}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked/ -> data/output/cell-type/
```

Two samples, three cell types of 100 cells each per sample and a 10-cell `RareType`. Each main type
opens 40 sites of its own; TypeA's carry a planted motif (`FIXA` in `motifs.meme`), and the 8 genes
next to each type's sites are expressed only in that type. The tests check that the marker peaks are
those sites, that FIXA is TypeA's top motif and that the links recover the planted genes.

## Real data

The demo cut (`atac-prep/scripts/cut_sympathetic.py`, Discover 465) gives the ATAC side
(`work/sympathetic/atac-prep/atac/`) and the RNA side (`work/sympathetic/cut/cells/`). The authors'
cell types in Discover 465's `meta.data.csv` do not line up with its cell names (the RNA shows no
marker specificity for them), so the demo labels come from the ATAC data itself:
`scripts/demo_labels.py` (SnapATAC2 spectral embedding + Leiden, `cluster` column, the shape
sc-integrate writes) → `work/sympathetic/labels-atac-clusters/labels.csv`.

```sh
python scripts/demo_labels.py .../atac-prep/atac .../labels-atac-clusters/labels.csv
docker run --rm --platform linux/amd64 --memory 6g --network none --read-only --tmpfs /tmp \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output -e LABEL_COLUMN=cluster -e THREADS=2 \
  -v .../atac-prep/atac:/data/input/atac:ro -v .../labels-atac-clusters:/data/input/labels:ro \
  -v .../cut/cells:/data/input/cells:ro -v .../work/sympathetic/atac-regulome/regulome:/data/output atac-regulome:dev
```
