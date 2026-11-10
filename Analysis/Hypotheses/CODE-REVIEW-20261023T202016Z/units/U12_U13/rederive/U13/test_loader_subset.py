"""Quick end-to-end test of u13_spatial.py on a small spot subset (login-node sized).

Reads each h5ad in backed mode, keeps the first N spots of every sample, and runs
main() with a tiny null (B=3) so every step after the loader is exercised once.
Usage: python test_loader_subset.py <dataset> <out_dir> [spots_per_sample] [family]
"""
import sys
from pathlib import Path

import anndata
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import u13_spatial as m  # noqa: E402

dataset, out = sys.argv[1], sys.argv[2]
per = int(sys.argv[3]) if len(sys.argv) > 3 else 300
_read = anndata.read_h5ad


def read_subset(path, *a, **k):
    b = _read(path, backed="r")
    sid = b.obs["sample_id"].astype(str).to_numpy()
    keep = np.sort(np.concatenate([np.flatnonzero(sid == s)[:per] for s in np.unique(sid)]))
    sub = b[keep].to_memory()
    b.file.close()
    print(f"subset {dataset}: {sub.n_obs} spots x {sub.n_vars} genes from {len(np.unique(sid))} samples", flush=True)
    return sub


anndata.read_h5ad = read_subset
m.sha256_file = lambda p: "skipped-in-test"
sys.argv = ["u13_spatial.py", "--dataset", dataset, "--seed", str(m.PRODUCER_SEED[dataset]),
            "--order", "registry", "--workers", "1", "--n-null", "3", "--n-sens", "2",
            "--out", out] + ([] if (len(sys.argv) > 4 and sys.argv[4] == "family") else ["--skip-family"])
sys.exit(m.main())
