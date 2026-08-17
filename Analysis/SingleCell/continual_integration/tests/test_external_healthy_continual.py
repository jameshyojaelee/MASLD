import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from masld_cl.config import load_config
from masld_cl.external_healthy_continual import (
    _external_healthy_centroids,
    _v22_gate_decision,
    load_external_healthy_continual_policy,
)
from masld_cl.external_healthy_anchor import (
    _proper_rigid_alignment,
    load_external_healthy_anchor_policy,
)
from masld_cl.external_healthy_direct_freeze import (
    load_external_healthy_direct_freeze_policy,
)


class TestExternalHealthyContinual(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(__file__).resolve().parents[1]
        cls.config = load_config(cls.root / "config_v1.json")
        _, cls.policy, decision = load_external_healthy_continual_policy(
            cls.config,
            cls.root / "reference" / "external_healthy_continual_policy_v22.json",
        )
        assert decision["selected_roster"] == "common_strict7"
        _, cls.anchor_policy, parent = load_external_healthy_anchor_policy(
            cls.config,
            cls.root / "reference" / "external_healthy_anchor_policy_v23.json",
        )
        assert parent["decision"] == "reject_paper_setting_pilot"
        _, cls.direct_policy, parent = load_external_healthy_direct_freeze_policy(
            cls.config,
            cls.root / "reference" / "external_healthy_direct_freeze_policy_v24.json",
        )
        assert parent["decision"] == "reject_paper_setting_pilot"

    def test_every_prespecified_gate_is_required(self):
        arguments = dict(
            reference_f1={"ci_low": 0.0},
            jaccard_loss=0.01,
            query_macro_change=0.0,
            query_lineage_changes={"a": 0.0, "b": 0.0},
            external_macro_change=0.0,
            external_lineage_changes={"a": 0.0, "b": 0.0},
            primary_alignment={"improvement": 0.01, "ci_low": 0.001},
            external_alignment={"improvement": 0.11, "ci_low": 0.001},
        )
        gates, decision = _v22_gate_decision(self.policy, **arguments)
        self.assertEqual(len(gates), 8)
        self.assertTrue(all(item["pass"] for item in gates))
        self.assertEqual(decision, "pass_paper_setting_pilot")
        arguments["external_alignment"] = {"improvement": 0.099, "ci_low": 0.05}
        gates, decision = _v22_gate_decision(self.policy, **arguments)
        self.assertFalse(gates[-1]["pass"])
        self.assertEqual(decision, "reject_paper_setting_pilot")

    def test_external_alignment_uses_seven_plus_nineteen_donors(self):
        cells = pd.DataFrame({
            "cell_id": (
                [f"ref-cell-{i}" for i in range(7)]
                + [f"HLiCA|cell-{i}" for i in range(19)]
            ),
            "donor_id": (
                [f"ref-{i}" for i in range(7)]
                + [f"HLiCA|donor-{i}" for i in range(19)]
            ),
            "dataset": ["reference"] * 7 + ["external"] * 19,
            "strict_reference": [True] * 7 + [False] * 19,
        })
        bundle = ({}, np.arange(52, dtype=float).reshape(26, 2), cells)
        _, values, _, reference, external = _external_healthy_centroids(bundle)
        self.assertEqual(values.shape, (26, 2))
        self.assertEqual(int(reference.sum()), 7)
        self.assertEqual(int(external.sum()), 19)

    def test_proper_rigid_alignment_recovers_rotation_and_translation(self):
        rng = np.random.default_rng(17)
        source = rng.normal(size=(100, 3))
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(q) < 0:
            q[:, -1] *= -1
        target = source @ q + np.asarray([2.0, -1.0, 0.5])
        source_mean, target_mean, rotation, _, error = _proper_rigid_alignment(
            source, target
        )
        observed = (source - source_mean) @ rotation + target_mean
        self.assertLess(error, 1e-12)
        self.assertGreater(np.linalg.det(rotation), 0)
        self.assertTrue(np.allclose(observed, target, atol=1e-10))


if __name__ == "__main__":
    unittest.main()
