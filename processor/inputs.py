"""
Finding the inputs under INPUT_DIR (symlinks followed, nothing outside INPUT_DIR and its staged
entries): atac-prep's <sample_id>.h5ad (port `atac`), labels.csv (port `labels`) and the RNA
side's .h5ad (port `cells`). Linked, each port is its own folder; unlinked, everything sits side
by side and the ATAC files are told from the RNA ones by their content (SnapATAC2 fragments).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from processor.safe_inputs import InputDir

MAX_LABELS_BYTES = 2 * 1024 ** 3


class InputError(ValueError):
    pass


def is_atac(path: Path) -> bool:
    import h5py

    try:
        with h5py.File(path, "r") as h:
            return "obsm" in h and "fragment_paired" in h["obsm"]
    except OSError:
        return False


@dataclass
class Found:
    inp: InputDir
    atac: list[Path]
    labels: Path | None
    rna: list[Path]


def discover(root: Path) -> Found:
    inp = InputDir(root)
    files = inp.walk()
    port = lambda name, f: not (root / name).is_dir() or f.is_relative_to(root / name)  # noqa: E731
    h5ads = [f for f in files if f.suffix.lower() == ".h5ad"]
    atac = [f for f in h5ads if port("atac", f) and is_atac(f)]
    rna = [f for f in h5ads if f not in atac and port("cells", f) and not is_atac(f)]
    labels = [f for f in files if f.name == "labels.csv" and port("labels", f)]
    if len(labels) > 1:
        raise InputError("several labels.csv: " + ", ".join(str(p.relative_to(root)) for p in labels))
    for f in atac + rna:
        inp.check_h5(f)               # no external links or external storage reaching outside the input
    return Found(inp=inp, atac=atac, labels=labels[0] if labels else None, rna=rna)


def read_labels(inp: InputDir, path: Path, column: str) -> pd.Series:
    """labels.csv -> label per obs_id (empty labels dropped), size-capped."""
    df = inp.read_csv(path, MAX_LABELS_BYTES, dtype=str, keep_default_na=False)
    for c in ("obs_id", column):
        if c not in df.columns:
            raise InputError(f"{path.name} has no {c!r} column (columns: {', '.join(df.columns)})")
    df = df[["obs_id", column]]
    df = df[df[column].str.strip() != ""]
    if df["obs_id"].duplicated().any():
        raise InputError(f"{path.name}: obs_id values repeat")
    return df.set_index("obs_id")[column].str.strip()
