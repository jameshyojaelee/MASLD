#!/usr/bin/env python3
"""Unit checks for the task-native inductive scVI/scANVI lane."""

from __future__ import annotations

from pathlib import Path
import unittest

import numpy as np
from scipy import sparse

from masld_bench.adapters import scvi_scanvi_inductive as lane


ROOT = Path(__file__).parents[2]
ACQUIRE = ROOT / "slurm/acquire_scvi_tools_1_4_2_inductive_runtime.sbatch"
PROBE = ROOT / "scripts/probe_scvi_scanvi_1_4_2_inductive.py"


class SCVISCANVIInductiveTests(unittest.TestCase):
    def test_raw_counts_reject_noninteger_values(self) -> None:
        lane.validate_integer_like_counts(
            sparse.csr_matrix(np.asarray([[0.0, 1.0], [2.0, 3.0]]))
        )
        with self.assertRaises(lane.InductiveSCVIError):
            lane.validate_integer_like_counts(
                sparse.csr_matrix(np.asarray([[0.0, 1.25]]))
            )

    def test_donor_class_weights_equalize_each_observed_stratum(self) -> None:
        donors = ["d1", "d1", "d2", "d2", "d2", "d3"]
        labels = ["immune", "immune", "immune", "immune", "immune", "hepatocyte"]
        weights = lane.donor_class_weights(donors, labels)
        totals = {}
        for donor, label, weight in zip(donors, labels, weights, strict=True):
            totals[(donor, label)] = totals.get((donor, label), 0.0) + weight
        np.testing.assert_allclose(list(totals.values()), list(totals.values())[0])
        self.assertAlmostEqual(float(weights.mean()), 1.0)

    def test_donor_class_weights_accept_numpy_arrays(self) -> None:
        """The production caller passes arrays, not lists.

        The list path above never exercised the length guard the way the fold
        loop does, so a guard written for lists raised on every real call.
        """
        donors = np.asarray(["d1", "d1", "d2"], dtype=str)
        labels = np.asarray(["immune", "immune", "hepatocyte"], dtype=str)
        weights = lane.donor_class_weights(donors, labels)
        self.assertEqual(len(weights), 3)
        self.assertAlmostEqual(float(weights.mean()), 1.0)
        with self.assertRaises(lane.InductiveSCVIError):
            lane.donor_class_weights(np.asarray([], dtype=str), np.asarray([], dtype=str))

    def test_weighted_knn_is_normalized_and_uses_frozen_roster(self) -> None:
        reference = np.asarray(
            [[class_index * 10.0, offset] for class_index in range(5) for offset in (0.0, 0.1, 0.2)],
            dtype=np.float32,
        )
        labels = [label for label in lane.ROSTER for _ in range(3)]
        query = np.asarray([[class_index * 10.0, 0.05] for class_index in range(5)])
        values = lane.weighted_knn_probabilities(
            reference,
            labels,
            np.ones(len(reference)),
            query,
            neighbors=3,
            n_jobs=1,
        )
        self.assertEqual(values.shape, (5, 5))
        np.testing.assert_allclose(values.sum(axis=1), 1.0)
        np.testing.assert_array_equal(np.argmax(values, axis=1), np.arange(5))

    def test_runtime_identity_and_query_firewall_are_explicit(self) -> None:
        text = ACQUIRE.read_text(encoding="utf-8")
        header = "\n".join(
            line for line in text.splitlines() if line.startswith("#SBATCH")
        )
        self.assertIn("--partition=io", header)
        self.assertIn("--account=nslab", header)
        self.assertIn("--qos=nslab", header)
        self.assertIn("--cpus-per-task=2", header)
        self.assertIn("--mem=64G", header)
        self.assertNotIn("--gres=gpu", header)
        self.assertNotIn("--array", header)
        self.assertIn("scvi-tools==1.4.2", text)
        self.assertIn(
            "c453b75b0fa0d222bf0eb3a531607d17f5ea0dbe5b848ab12fd8cc24599884e5",
            text,
        )
        probe = PROBE.read_text(encoding="utf-8")
        self.assertIn("direct_latent_inference", probe)
        self.assertIn("direct_scanvi_probabilities", probe)
        self.assertNotIn("load_query_data", probe)
        self.assertNotIn("prepare_query_anndata", probe)


if __name__ == "__main__":
    unittest.main()
