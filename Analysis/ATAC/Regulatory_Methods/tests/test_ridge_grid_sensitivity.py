"""Nested donor selection and original-readout identity invariants."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd
from scipy import sparse

root = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ridge = load("ridge_sensitivity", root / "scripts/ridge_grid_sensitivity.py")
core = load("original_core_fixture", root / "scripts/network_compare.py")


class RidgeInvariants(unittest.TestCase):
    def tearDown(self):
        ridge.CORE = None
        ridge.CONTEXT = None

    def test_tuning_averages_donors_equally(self):
        # Donor a has three validation rows; donor b has one. Their mean losses
        # are 1 and 9, so their equal-donor mean is 5, not the row mean 3.
        donors = np.array(["a", "a", "a", "b", "c", "d"])
        train = np.arange(6)
        inner = np.array([0, 0, 0, 0, 1, 1])
        y = np.zeros(6)
        design = np.array([[1.], [1.], [1.], [3.], [2.], [2.]])
        ridge.CORE = SimpleNamespace(weighted_ridge=lambda x, y, d, a, prediction: prediction[:, 0])
        losses = ridge.tune({0: design, 1: design}, y, donors, train, inner, [1., 10.])
        self.assertEqual(losses, {1.: 4.5, 10.: 4.5})

    def test_selection_is_relearned_inside_every_inner_training_fold(self):
        donors = np.repeat(np.array(["a", "b", "c", "d", "held"]), 2)
        values = np.array([[1., 2., 3.], [2., 1., 4.], [4., 2., 1.], [2., 4., 1.],
                           [3., 1., 2.], [1., 3., 2.], [2., 2., 3.], [4., 1., 2.],
                           [1., 2., 4.], [3., 2., 1.]])
        x = sparse.csr_matrix(values)
        train, test = np.arange(8), np.arange(8, 10)
        inner = np.array([0, 0, 1, 1, 2, 2, 0, 0])
        metadata = pd.DataFrame(dict(donor_id=donors))
        selections = []

        def select(method, candidates, prom, peaks, gene, train_x, train_y, train_donors, seed, max_links, graph):
            selections.append((method, tuple(train_donors), graph, train_x.shape[0]))
            return np.array([1]), np.array([.5])

        ridge.CORE = SimpleNamespace(**{name: getattr(core, name) for name in
            ["log_columns", "donor_weights", "stable_seed", "design_matrix", "weighted_ridge", "save_predictions"]},
            selected_indices=select)
        ridge.CONTEXT = dict(rna_lookup={"gene": 0}, geometry=SimpleNamespace(query=lambda *args: (np.array([0]), np.array([1, 2]))),
            window=250000, promoter_bp=2000, x=x, feature_x=x.tocsc(), train=train, test=test,
            inner=inner, donors=donors, rna_normal=sparse.csc_matrix(np.arange(1., 11.)[:, None]),
            fold=0, seed=42, peaks=pd.DataFrame(dict(peak_id=["prom", "dist1", "dist2"])),
            max_links=10, training_x=x[train].tocsc(),
            inner_training_x={k: x[train[inner != k]].tocsc() for k in range(3)},
            cicero={"outer": "outer", **{f"inner{k}": f"inner{k}" for k in range(3)}},
            metadata=metadata, lineage="hepatocyte")
        _, _, rows, links, tuning = ridge.fit_gene((1, {"gene_id": "gene"}))
        self.assertEqual(len(selections), 12)
        for method in ridge.METHODS:
            stages = [item for item in selections if item[0] == method]
            self.assertEqual(set(stages[0][1]), set(donors[train]))
            self.assertEqual(stages[0][2], "outer")
            for k, stage in enumerate(stages[1:]):
                self.assertEqual(set(stage[1]), set(donors[train[inner != k]]))
                self.assertTrue(set(stage[1]).isdisjoint(donors[train[inner == k]]))
                self.assertNotIn("held", stage[1])
                self.assertEqual(stage[2], f"inner{k}")
        self.assertEqual(len(rows), 8)  # One held donor, four methods, two grids.
        self.assertEqual(len(links), 3)
        self.assertEqual(len(tuning), 8)
        self.assertEqual({row["grid"] for row in rows}, {"wide", "original_replay"})
        self.assertTrue(all(row["donor_id"] == "held" for row in rows))

    def test_anchor_rejects_coverage_alpha_and_error_changes(self):
        row = dict(outer_fold=0, gene_id="g", endpoint="rna", method="distance", donor_id="d",
                   n_metacells=2, input_representation="binary", target_unit="log1p_RNA_counts_per_10000",
                   ridge_alpha=100., observed=1., predicted=2., squared_error=1., standardized_squared_error=.5,
                   grid="original_replay")
        expected = pd.DataFrame([row]).drop(columns="grid")
        self.assertEqual(ridge.verify_anchor([row], expected), (0., 1))
        for key, changed in [("donor_id", "other"), ("ridge_alpha", 10.), ("squared_error", 2.)]:
            actual = {**row, key: changed}
            with self.subTest(field=key), self.assertRaises(ValueError):
                ridge.verify_anchor([actual], expected)
        with self.assertRaises(ValueError):
            ridge.verify_anchor([], expected)
        with self.assertRaises(ValueError):
            ridge.verify_anchor([row, row], expected)


if __name__ == "__main__":
    unittest.main()
