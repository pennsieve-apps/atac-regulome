"""
Demo only: labels for the Discover 465 cut from the ATAC data itself, because the authors'
meta.data.csv cell types do not line up with its cell names (see README "Results"). SnapATAC2's
standard clustering of atac-prep's tile matrices (select_features, spectral, knn, leiden), written
as labels.csv with a `cluster` column, the shape sc-integrate writes.

  python scripts/demo_labels.py <atac-prep atac/ folder> <out labels.csv>
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import snapatac2 as snap

if __name__ == "__main__":
    atac, out = Path(sys.argv[1]), Path(sys.argv[2])
    with tempfile.TemporaryDirectory() as tmp:
        ds = snap.AnnDataSet(adatas=[(f.stem, snap.read(f, backed="r")) for f in sorted(atac.glob("*.h5ad"))],
                             filename=Path(tmp) / "d.h5ads")
        snap.pp.select_features(ds, n_features=250000)
        snap.tl.spectral(ds, random_state=0)
        snap.pp.knn(ds, random_state=0)
        snap.tl.leiden(ds, resolution=0.5, random_state=0)
        names = np.asarray(ds.obs_names)
        cl = np.asarray(ds.obs["leiden"])
        ds.close()
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"sample_id": [n.split(":")[0] for n in names], "obs_id": names,
                  "cluster": [f"C{int(c) + 1:02d}" for c in cl]}).to_csv(out, index=False)
    print(pd.Series(cl).value_counts().to_dict())
