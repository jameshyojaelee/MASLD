"""Small matrix checks for the DANB count guard without loading an atlas."""
import ast
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.sparse import csr_matrix, issparse

root = Path(__file__).resolve().parents[5]
path = root / "Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py"
tree = ast.parse(path.read_text())
function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "ensure_raw_layer")
namespace = dict(np=np, issparse=issparse, ad=SimpleNamespace(AnnData=object))
exec(compile(ast.Module(body=[function], type_ignores=[]), str(path), "exec"), namespace)
guard = namespace["ensure_raw_layer"]
raw = np.array([[0, 2], [1, 3]], dtype=float)
for matrix in (raw, csr_matrix(raw)):
    adata = SimpleNamespace(X=matrix, layers={}, obs={})
    guard(adata)
    assert "counts" in adata.layers
    assert np.array_equal(adata.obs["n_counts"], np.array([2, 4]))
for matrix in (np.log1p(raw), np.array([[0, -1]]), np.array([[0, np.nan]]), csr_matrix(np.log1p(raw))):
    adata = SimpleNamespace(X=matrix, layers={}, obs={})
    try:
        guard(adata)
    except ValueError:
        assert "counts" not in adata.layers
    else:
        raise AssertionError("Non-count X was accepted")
adata = SimpleNamespace(X=np.log1p(raw), layers={"counts": csr_matrix(raw)}, obs={})
guard(adata)
assert np.array_equal(adata.obs["n_counts"], np.array([2, 4]))
adata = SimpleNamespace(X=raw, layers={"counts": np.log1p(raw)}, obs={})
try:
    guard(adata)
except ValueError:
    pass
else:
    raise AssertionError("A normalized counts layer was accepted")
print("PASS: dense/sparse integer counts accepted; normalized, fractional, negative and nonfinite inputs rejected before assigning counts")
