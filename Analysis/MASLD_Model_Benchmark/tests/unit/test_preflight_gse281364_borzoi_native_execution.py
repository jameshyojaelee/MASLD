from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.preflight_gse281364_borzoi_native_execution import (
    BorzoiExecutionError,
    EXPECTED_GATE_ORDER,
    observe_local_port,
    validate_checkpoints,
    validate_execution_gates,
    validate_native_parity,
    validate_runtime,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "config/gse281364_borzoi_native_execution.json"


class BorzoiNativeExecutionPreflightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

    def test_current_contract_is_fail_closed(self) -> None:
        self.assertEqual(
            tuple(self.config["required_execution_artifacts"]), EXPECTED_GATE_ORDER
        )
        self.assertEqual(
            {
                value["artifacts_sha256"]
                for value in self.config["required_execution_artifacts"].values()
            },
            {"UNRESOLVED"},
        )
        self.assertEqual(
            self.config["status"],
            "blocked_missing_terms_weights_runtime_reference_parity_native_parity",
        )
        self.assertFalse(
            self.config["official_release_contract"][
                "westminster_required_for_zero_shot_inference"
            ]
        )
        self.assertEqual(
            self.config["official_release_contract"]["baskerville_revision"],
            "UNRESOLVED",
        )
        self.assertFalse(
            self.config["policy_transition"]["current_model_forward_authorized"]
        )

    def test_unbound_gate_census_cannot_be_executable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gates, blockers = validate_execution_gates(root, self.config)
        self.assertEqual(tuple(gates), EXPECTED_GATE_ORDER)
        self.assertTrue(all(not gate["passed"] for gate in gates.values()))
        self.assertEqual(
            {value for value in blockers if value.startswith("artifact_hash_unresolved")},
            {f"artifact_hash_unresolved:{name}" for name in EXPECTED_GATE_ORDER},
        )
        self.assertEqual(
            {value for value in blockers if value.startswith("artifact_missing")},
            {f"artifact_missing:{name}" for name in EXPECTED_GATE_ORDER},
        )

    def test_registered_local_port_is_never_native(self) -> None:
        observation = observe_local_port(self.config)
        self.assertTrue(observation["observed"])
        self.assertFalse(observation["eligible_as_native_runtime"])
        self.assertFalse(observation["modules_imported"])
        self.assertFalse(observation["checkpoint_bytes_opened"])
        self.assertEqual(
            observation["identities"]["tensorflow_distribution_count"], "0"
        )
        self.assertEqual(
            observation["identities"]["baskerville_distribution_count"], "0"
        )

    def test_pytorch_port_cannot_satisfy_native_runtime_receipt(self) -> None:
        release = self.config["official_release_contract"]
        receipt = {
            "status": "pass_offline_native_inference_runtime",
            "borzoi_revision": release["borzoi_revision"],
            "python_version": "3.11.9",
            "tensorflow_version": "0.0.0",
            "h5py_version": "0.0.0",
            "numpy_version": "2.4.6",
            "westminster_required": False,
            "network_access": False,
            "baskerville_revision": "0" * 40,
            "python_executable": {"path": "bin/python", "sha256": "1" * 64},
            "prediction_entrypoint": {
                "path": "drivers/predict.py",
                "sha256": "2" * 64,
            },
        }
        with self.assertRaisesRegex(BorzoiExecutionError, "python_version"):
            validate_runtime(receipt, self.config)

    def test_checkpoint_gate_requires_all_exact_members(self) -> None:
        members = []
        for configured in self.config["official_checkpoint_members"]:
            members.append(
                {
                    **configured,
                    "path": f"weights/f{configured['replicate']}/model0_best.h5",
                    "sha256": str(configured["replicate"] + 1) * 64,
                }
            )
        receipt = {
            "status": "pass_strict_native_checkpoint_bundle",
            "external_or_soft_hdf5_links_present": False,
            "strict_tensor_inventory_passed": True,
            "members": members,
        }
        result = validate_checkpoints(receipt, self.config)
        self.assertEqual([value["replicate"] for value in result["members"]], [0, 1, 2, 3])
        receipt["members"][2]["generation"] = "changed"
        with self.assertRaisesRegex(BorzoiExecutionError, "generation"):
            validate_checkpoints(receipt, self.config)

    def test_numeric_parity_requires_every_native_path(self) -> None:
        receipt = {
            "status": "pass_all_member_native_numeric_parity",
            "strict_restore_passed": True,
            "forward_passed": True,
            "reverse_complement_passed": True,
            "strand_pair_restore_passed": True,
            "inverse_transform_passed": True,
            "zero_shift_passed": True,
            "determinism_passed": True,
            "members_passed": [0, 1, 2, 3],
            "human_tracks": 7611,
            "raw_output_bins": 16384,
            "released_prediction_bins": 16352,
            "central_training_loss_bins": 6144,
            "numeric_tolerance": {"passed": True, "rtol": 1e-5, "atol": 1e-6},
        }
        self.assertEqual(
            validate_native_parity(receipt, self.config)["members_passed"],
            [0, 1, 2, 3],
        )
        receipt["strand_pair_restore_passed"] = False
        with self.assertRaisesRegex(BorzoiExecutionError, "strand_pair_restore_passed"):
            validate_native_parity(receipt, self.config)


if __name__ == "__main__":
    unittest.main()
