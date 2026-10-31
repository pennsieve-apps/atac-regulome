"""
Finding the inputs under INPUT_DIR (symlinks followed, nothing outside INPUT_DIR and its staged
entries): atac-prep's <sample_id>.h5ad (port `atac`), labels.csv (port `labels`) and the RNA
side's .h5ad (port `cells`). Linked, each port is its own folder; unlinked, everything sits side
by side and the ATAC files are told from the RNA ones by their content (SnapATAC2 fragments).
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

MAX_LABELS_BYTES = 2 * 1024 ** 3


class InputError(ValueError):
    pass


def walk(root: Path) -> list[Path]:
    allowed = [Path(os.path.realpath(root))] + [Path(os.path.realpath(e)) for e in root.iterdir() if e.is_symlink()]
    found, seen = [], set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen:
            dirnames[:] = []
            continue
        seen.add(real)
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        found.extend(Path(dirpath) / f for f in sorted(filenames) if not f.startswith("."))
    return [p for p in found if p.is_file() and any(Path(os.path.realpath(p)).is_relative_to(a) for a in allowed)]


def is_atac(path: Path) -> bool:
    import h5py

    try:
        with h5py.File(path, "r") as h:
            return "obsm" in h and "fragment_paired" in h["obsm"]
    except OSError:
        return False


@dataclass
class Found:
    atac: list[Path]
    labels: Path | None
    rna: list[Path]


def discover(root: Path) -> Found:
    files = walk(root)
    port = lambda name, f: not (root / name).is_dir() or f.is_relative_to(root / name)  # noqa: E731
    h5ads = [f for f in files if f.suffix.lower() == ".h5ad"]
    atac = [f for f in h5ads if port("atac", f) and is_atac(f)]
    rna = [f for f in h5ads if f not in atac and port("cells", f) and not is_atac(f)]
    labels = [f for f in files if f.name == "labels.csv" and port("labels", f)]
    if len(labels) > 1:
        raise InputError("several labels.csv: " + ", ".join(str(p.relative_to(root)) for p in labels))
    return Found(atac=atac, labels=labels[0] if labels else None, rna=rna)


def read_labels(path: Path, column: str) -> pd.Series:
    """labels.csv -> label per obs_id (empty labels dropped)."""
    if path.stat().st_size > MAX_LABELS_BYTES:
        raise InputError(f"{path.name} is larger than 2 GB")
    head = pd.read_csv(path, nrows=0).columns
    for c in ("obs_id", column):
        if c not in head:
            raise InputError(f"{path.name} has no {c!r} column (columns: {', '.join(head)})")
    df = pd.read_csv(path, usecols=["obs_id", column], dtype=str, keep_default_na=False)
    df = df[df[column].str.strip() != ""]
    if df["obs_id"].duplicated().any():
        raise InputError(f"{path.name}: obs_id values repeat")
    return df.set_index("obs_id")[column].str.strip()
