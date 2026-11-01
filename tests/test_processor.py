"""The synthetic fixture through the processor: both input layouts, with and without RNA, parameter checks."""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pytest

from processor import main as app
from tests.make_fixture import BLACKLISTED, TRUTH, TYPES, make_fixture

FILES = {"peaks.bed", "marker_peaks.csv", "motifs.csv", "summary.json", "provenance.json", "peaks_per_label.png",
         "marker_peaks_heatmap.png", "motifs_heatmap.png"}


@pytest.fixture(scope="session")
def fx(tmp_path_factory) -> dict[str, Path]:
    return make_fixture(tmp_path_factory.mktemp("fixture"), "both")


def run(input_dir: Path, output_dir: Path, fx, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("GENOMES_DIR", str(fx["genomes"]))
    monkeypatch.setenv("MOTIFS_FILE", str(fx["motifs"]))
    monkeypatch.setenv("MPLCONFIGDIR", "/tmp/matplotlib")
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    for name, value in {"THREADS": 1, "MIN_CELLS": 20, "PSEUDOCELL_SIZE": 10, **params}.items():
        monkeypatch.setenv(name, str(value))
    return app.main()


def rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_regulome(fx, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(fx[layout], out, fx, monkeypatch) == 0
    d = out / "cell-type"
    assert FILES | {"links.csv", "links.png"} <= {p.name for p in d.iterdir()}
    assert sorted(p.name for p in (d / "peaks").iterdir()) == sorted(f"{t}.bed" for t in (*TYPES, "RareType"))
    peaks = (d / "peaks.bed").read_text().splitlines()
    assert peaks and all(len(line.split("\t")) == 4 for line in peaks)
    markers = rows(d / "marker_peaks.csv")
    for t in TYPES:
        mine = [m for m in markers if m["label"] == t]
        assert len(mine) >= 30          # 40 planted sites per type
        # each marker peak lies at one of the type's planted sites
        sites = TRUTH["sites"][t]
        assert all(any(m["chrom"] == c and int(m["start"]) - 300 <= p <= int(m["end"]) + 300 for c, p in sites)
                   for m in mine)
    for chrom, centre in TRUTH["sites"]["TypeA"][:BLACKLISTED]:      # blacklisted sites have no peak
        assert not any(line.split("\t")[0] == chrom and int(line.split("\t")[1]) - 1000 < centre < int(line.split("\t")[2]) + 1000
                       for line in peaks)
    motifs = rows(d / "motifs.csv")
    top = {t: min((m for m in motifs if m["label"] == t), key=lambda m: float(m["p_value"])) for t in TYPES}
    assert top["TypeA"]["motif_name"] == "FIXA" and float(top["TypeA"]["adjusted_p_value"]) < 1e-6
    links = rows(d / "links.csv")
    for t in TYPES:   # the type's marker peaks link to the genes expressed in that type only
        mine = {m["peak"] for m in markers if m["label"] == t}
        genes = [lk["gene"] for lk in links if lk["peak"] in mine]
        assert sum(g in TRUTH["linked"][t] for g in genes) >= 0.9 * len(genes)
        assert len(set(genes) & set(TRUTH["linked"][t])) >= 6      # of the 8 planted
    summary = json.loads((d / "summary.json").read_text())
    assert summary["n_cells_used"] == 620 and summary["links"]["n_links"] == len(links)
    assert summary["notes"] == [] and "skipped" not in summary["links"]
    assert next(x for x in summary["labels"] if x["label"] == "TypeA")["top_motifs"][0] == "FIXA"
    prov = json.loads((d / "provenance.json").read_text())
    assert {"app", "version", "tools", "parameters", "inputs", "started", "finished", "samples", "peak_rss_mb",
            "citations"} <= set(prov)
    assert all("doi" in c and "text" in c for c in prov["citations"])


def test_without_rna_and_min_cells(fx, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fx["linked"] / "atac", src / "atac")
    shutil.copytree(fx["linked"] / "labels", src / "labels")
    out = tmp_path / "out"
    assert run(src, out, fx, monkeypatch, MIN_CELLS=50, LABEL="types") == 0
    d = out / "types"
    assert not (d / "links.csv").exists()
    summary = json.loads((d / "summary.json").read_text())
    assert summary["labels_left_out"] == {"RareType": 20} and summary["links"]["rna_linked"] is False


def test_links_skipped_says_why(fx, tmp_path, monkeypatch):
    """Too few pseudo-cells: no links, and summary.json says why (sc-report prints the notes)."""
    out = tmp_path / "out"
    assert run(fx["linked"], out, fx, monkeypatch, PSEUDOCELL_SIZE=200) == 0
    summary = json.loads((out / "cell-type" / "summary.json").read_text())
    assert summary["links"]["n_links"] == 0 and summary["links"]["n_pseudocells"] < 5
    assert summary["notes"] == [f"Peak-to-gene links were not computed: {summary['links']['skipped']}."]
    assert "PSEUDOCELL_SIZE" in summary["notes"][0]


def test_labels_not_matching(fx, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fx["linked"] / "atac", src / "atac")
    (src / "labels").mkdir()
    (src / "labels" / "labels.csv").write_text("sample_id,obs_id,cell_type\ns1,AAAA-1,TypeA\n")
    assert run(src, tmp_path / "out", fx, monkeypatch) == 1


@pytest.mark.parametrize("params", [{"LABEL": "../x"}, {"PEAK_QVALUE": "2"}, {"LINK_MIN_COR": "abc"},
                                    {"MIN_CELLS": "0"}, {"LABEL_COLUMN": "a b"}])
def test_bad_parameters(fx, tmp_path, monkeypatch, params):
    assert run(fx["linked"], tmp_path / "out", fx, monkeypatch, **params) == 2


def test_platform_link_layout(fx, tmp_path, monkeypatch):
    """As the platform stages a run: real port folders whose entries are symlinks into the upstream
    apps' output folders (siblings of INPUT_DIR)."""
    up = tmp_path / "upstream"
    shutil.copytree(fx["linked"], up)
    src = tmp_path / "in"
    for port in ("atac", "labels", "cells"):
        (src / port).mkdir(parents=True)
        for f in (up / port).iterdir():
            (src / port / f.name).symlink_to(f)
    out = tmp_path / "out"
    assert run(src, out, fx, monkeypatch) == 0
    prov = json.loads((out / "cell-type" / "provenance.json").read_text())
    assert "labels/labels.csv" in prov["inputs"] and len(prov["samples"]) == 2
    assert len(rows(out / "cell-type" / "links.csv")) > 0


def test_labels_that_contain_each_other(fx, tmp_path, monkeypatch):
    """n_links per label counts exact label names, not substrings ("T" in "Fibroblast_2")."""
    src = tmp_path / "in"
    shutil.copytree(fx["linked"], src)
    lab = (src / "labels" / "labels.csv").read_text()
    for old, new in (("TypeA", "T"), ("TypeB", "Fibroblast"), ("TypeC", "Fibroblast_2")):
        lab = lab.replace(f",{old},", f",{new},")
    (src / "labels" / "labels.csv").write_text(lab)
    out = tmp_path / "out"
    assert run(src, out, fx, monkeypatch) == 0
    d = out / "cell-type"
    links = rows(d / "links.csv")
    summary = json.loads((d / "summary.json").read_text())
    for x in summary["labels"]:
        assert x["n_links"] == sum(x["label"] in lk["called_in"].split(";") for lk in links)
    assert sorted(p.name for p in (d / "peaks").iterdir()) == sorted(
        ["T.bed", "Fibroblast.bed", "Fibroblast_2.bed", "RareType.bed"])


def test_external_hdf5_link_refused(fx, tmp_path, monkeypatch):
    import h5py

    src = tmp_path / "in"
    shutil.copytree(fx["linked"], src)
    with h5py.File(src / "atac" / "s1.h5ad", "a") as h:
        h["obsm"]["leak"] = h5py.ExternalLink("/etc/hosts", "/")
    assert run(src, tmp_path / "out", fx, monkeypatch) == 1


def test_missing_label_column_is_a_parameter_error(fx, tmp_path, monkeypatch):
    assert run(fx["linked"], tmp_path / "out", fx, monkeypatch, LABEL_COLUMN="no_such_column") == 2
