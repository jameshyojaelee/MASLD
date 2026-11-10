"""Small exact-null/invariant checks; synthetic donors, no biological raw data."""
import importlib.util
import itertools
from pathlib import Path
import numpy as np

P = Path(__file__).resolve().parents[3] / "X1-20261028T103410Z/scripts/x1_selection_aware.py"
spec = importlib.util.spec_from_file_location("selection", P)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
rng = np.random.default_rng(20261107)
genes = np.array(["HPX", "KNG1", "FBP1", "TEST"])
raw = rng.uniform(.05, 2, (6, 3, 4))
raw[:, 0, :3] = 20  # a finite, invariant hepatocyte spillover reference
n_int = np.array([[100, 40, 0], [100, 40, 0], [100, 40, 40],
                  [100, 0, 40], [100, 0, 40], [100, 40, 40]])
a = dict(genes=genes, types=np.array(["Hepatocyte", "A", "B"]),
         donors=np.array(list("abcdef")), m_int=raw.copy(), m_all=raw * 1.1,
         n_all=np.full((6, 3), 100), n_int=n_int)
engine = module.SelectionReplay(a, minimum_arm=2)
masks = np.array([np.isin(np.arange(6), c) for c in itertools.combinations(range(6), 3)])
indices, use = engine.choices(masks)
assert use[:, 0].all() and np.any(use[:, 1]) and not np.all(use[:, 1])
bulk = engine.statistics("TEST", masks, indices, use)
for i, mask in enumerate(masks):
    one = engine.statistics("TEST", mask[None, :])
    for key in bulk:
        np.testing.assert_allclose(bulk[key][i], one[key][0])
    # Under the fixed exchangeable null, every possible observed assignment is
    # equally likely. Exact tail p values must be superuniform, even with selection.
values = bulk["within"]
pvalues = np.array([np.mean(np.abs(values) >= abs(v) - 1e-12) for v in values])
for alpha in np.arange(1, 21) / 20:
    assert np.mean(pvalues <= alpha) <= alpha + 1e-12
assert np.isfinite(values).all()
np.testing.assert_array_equal(module.bounded_difference(np.array([0., 1., 0.]),
                                                       np.array([0., 0., 1.])), [0., 1., -1.])
np.testing.assert_array_equal(module.holm([.01, .02, 1]), [.03, .04, 1])
# Resampling one donor twice must count both donor draws in the eligibility rule.
bi = np.array([[0, 0, 1, 3, 3, 4]])
bm = np.array([[False, False, False, True, True, True]])
_, buse = engine.choices(bm, bi)
assert not buse[0, 1] and not buse[0, 2]
print("PASS: vectorized replay, label-dependent cellsets, exact-null calibration, zeros, full Holm family, bootstrap draw eligibility")
