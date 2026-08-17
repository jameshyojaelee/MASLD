from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from masld_cl.control_evaluation import (
    _neighbor_jaccard_loss,
    _reference_f1_change_ci,
    _shift_and_standard_error,
)


class TestControlEvaluation(unittest.TestCase):
    def test_control_only_metrics_use_donors_and_are_reproducible(self):
        rows = []
        latent = []
        rng = np.random.default_rng(17)
        for role, count in (("reference", 6), ("control", 6)):
            for donor_index in range(count):
                donor = f"{role}_{donor_index}"
                study = f"{role}_study_{donor_index % 2}"
                for cell in range(6):
                    label = "A" if cell < 3 else "B"
                    rows.append({
                        "cell_id": f"{donor}_{cell}", "donor_id": donor, "dataset": study,
                        "audit_cell_type": label, "predicted_cell_type": label,
                        "strict_reference": role == "reference",
                        "primary_query": role == "control", "query_control": role == "control",
                    })
                    center = 0.0 if role == "reference" else 0.5
                    latent.append([center + donor_index * 0.05, cell * 0.02] + rng.normal(0, .001, 2).tolist())
        cells = pd.DataFrame(rows)
        latent = np.asarray(latent)
        shift_a = _shift_and_standard_error(cells, latent, 200, 31)
        shift_b = _shift_and_standard_error(cells, latent, 200, 31)
        self.assertEqual(shift_a, shift_b)
        self.assertGreater(shift_a[0], 0)

    def test_reference_f1_and_neighborhood_identity(self):
        cells = pd.DataFrame({
            "cell_id": [f"c{i}" for i in range(40)],
            "donor_id": [f"d{i // 10}" for i in range(40)],
            "dataset": [f"s{(i // 10) % 2}" for i in range(40)],
            "audit_cell_type": ["A" if i % 2 else "B" for i in range(40)],
            "predicted_cell_type": ["A" if i % 2 else "B" for i in range(40)],
        })
        result = _reference_f1_change_ci(cells, cells.copy(), 200, 17)
        self.assertEqual(result["estimate"], 0)
        latent = np.random.default_rng(17).normal(size=(40, 5))
        loss = _neighbor_jaccard_loss(latent, latent.copy(), cells, k=5, maximum_cells=40, seed=17)
        self.assertAlmostEqual(loss, 0.0)
