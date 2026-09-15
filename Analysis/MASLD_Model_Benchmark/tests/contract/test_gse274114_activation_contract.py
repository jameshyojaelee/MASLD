from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import unittest

from masld_bench.registry import load_dataset_manifest, load_task_spec


ROOT = Path(__file__).resolve().parents[2]


class GSE274114ActivationContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.dataset = load_dataset_manifest(
            ROOT / "config/datasets/gse274114_mash_hbv.toml"
        )
        cls.task = load_task_spec(
            ROOT / "config/evaluation/gse274114_within_instrument_etiology_task.toml"
        )

    def test_participant_and_source_topology(self) -> None:
        self.assertEqual(self.dataset.biological_unit, "participant")
        self.assertEqual(self.dataset.expected_biological_units, 39)
        self.assertEqual(
            set(self.dataset.accession), {"GSE274114", "PRJNA1144975", "SRP524555"}
        )
        self.assertTrue(self.dataset.admission_blocking)

    def test_task_prohibits_instrument_confounded_claims(self) -> None:
        parameters = self.task.evaluator_parameters
        self.assertFalse(parameters["four_class_performance_allowed"])
        self.assertFalse(parameters["mash_vs_non_mash_performance_allowed"])
        self.assertFalse(parameters["participant_metadata_inputs_available"])
        self.assertEqual(parameters["par_y_x_two_row_targets"], 44)
        self.assertEqual(parameters["unmapped_source_rows"], 1230)
        self.assertIn("No four-class", self.task.claim_gate)

    def test_hbv_only_is_not_registered_as_masld_negative(self) -> None:
        gates = " ".join(self.task.admission_gates)
        self.assertIn("never a MASLD-negative", gates)
        self.assertIn("CTRL versus ENEG", gates)
        self.assertIn("NASH versus ENEG_NASH", gates)

    def test_promotion_gate_hash_is_exact(self) -> None:
        path = ROOT / "config/evaluation/gse274114_within_instrument_etiology_promotion_gate.json"
        self.assertEqual(self.task.promotion_gate_config_sha256, sha256(path.read_bytes()).hexdigest())


if __name__ == "__main__":
    unittest.main()
