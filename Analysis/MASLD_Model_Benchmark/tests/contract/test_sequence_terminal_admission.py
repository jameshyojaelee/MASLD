from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tomllib
import unittest

from masld_bench.contracts import ModelManifest


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config"
AUTHORITY_PATH = CONFIG / "evaluation/sequence_terminal_admission.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_toml(path: Path) -> dict:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def model_entry(path: Path, model_id: str) -> tuple[dict, ModelManifest]:
    family = read_toml(path)
    entry = next(model for model in family["models"] if model["model_id"] == model_id)
    return entry, ModelManifest.from_family_entry(family, entry)


class SequenceTerminalAdmissionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.authority = read_json(AUTHORITY_PATH)

    def test_terminal_bundle_hashes_and_all_eight_receipts_are_exact(self) -> None:
        expected = {
            "bpnet": (
                "14af86ac22e67bc29f8a7b85ee3a049b63f875ba5423210cc8e539d29b9fe9e8",
                "predictions/summary.json",
            ),
            "sequence_cnn_control": (
                "b379da6299a7fd9b42c77da7da0150403a6268683cc8ecd1a35dd3ad983a33a6",
                "predictions/fixed_ccre/summary.json",
            ),
        }
        for model_id, (expected_hash, prediction_summary) in expected.items():
            binding = self.authority["models"][model_id]
            bundle = ROOT / binding["bundle_path"]
            self.assertEqual(digest(bundle / "ARTIFACTS.json"), expected_hash)
            self.assertEqual(binding["bundle_artifacts_sha256"], expected_hash)
            self.assertTrue((bundle / "COMPLETE").is_file())
            summary = read_json(bundle / "bundle_summary.json")
            self.assertEqual(summary["terminal_status"], "passed")
            self.assertEqual(summary["logical_tasks"], 8)
            self.assertEqual(summary["status_counts"], {"passed": 8})
            receipts = sorted((bundle / "task_receipts").glob("*.json"))
            self.assertEqual(len(receipts), 8)
            for receipt_path in receipts:
                receipt = read_json(receipt_path)
                self.assertEqual(receipt["status"], "passed")
                self.assertEqual(receipt["model_id"], model_id)
                artifact = Path(receipt["expected_artifact"])
                self.assertTrue((artifact / "COMPLETE").is_file())
                self.assertEqual(
                    digest(artifact / "ARTIFACTS.json"),
                    receipt["artifacts_sha256"],
                )
                self.assertTrue((artifact / prediction_summary).is_file())

    def test_frozen_subordinate_receipts_prove_required_false_firewall(self) -> None:
        required_false = self.authority["required_false_firewall"]
        self.assertTrue(required_false)
        self.assertTrue(all(value is False for value in required_false.values()))
        for model_id, binding in self.authority["models"].items():
            bundle = ROOT / binding["bundle_path"]
            for receipt_path in sorted((bundle / "task_receipts").glob("*.json")):
                artifact = Path(read_json(receipt_path)["expected_artifact"])
                metadata = read_json(artifact / "ARTIFACTS.json")["metadata"]
                absent = {
                    field for field in required_false if field not in metadata
                }
                self.assertEqual(
                    absent,
                    set(binding["root_metadata_missing_false_fields"]),
                    f"root firewall gap changed for {artifact}",
                )
                if model_id == "bpnet":
                    training = read_json(artifact / "model/training_receipt.json")
                    prediction = read_json(artifact / "predictions/summary.json")
                    validation = read_json(artifact / "validation/summary.json")
                    self.assertIs(training["evaluator_outcomes_exposed"], False)
                    self.assertIs(training["observed_held_donor_atac_exposed"], False)
                    self.assertIs(prediction["benchmark_metrics_calculated"], False)
                    self.assertIs(prediction["observed_atac_input_exposed"], False)
                    self.assertIs(validation["evaluator_outcomes_exposed"], False)
                    self.assertIs(
                        validation["held_donor_outcomes_used_for_training"], False
                    )
                else:
                    training = read_json(
                        artifact / "model/sequence_cnn_control.training_receipt.json"
                    )
                    prediction = read_json(
                        artifact / "predictions/fixed_ccre/summary.json"
                    )
                    contract = read_json(
                        artifact
                        / "predictions/fixed_ccre/prediction_contract.json"
                    )
                    self.assertIs(training["held_donor_atac_exposed"], False)
                    self.assertIs(training["test_outcomes_used"], False)
                    self.assertIs(prediction["benchmark_metrics_calculated"], False)
                    self.assertIs(prediction["observed_atac_input_exposed"], False)
                    self.assertIs(contract["observed_atac_input_exposed"], False)

    def test_registry_admission_is_narrow_and_uncertainty_is_preserved(self) -> None:
        bpnet_entry, bpnet = model_entry(
            CONFIG / "models/regulatory_local.toml", "bpnet"
        )
        self.assertIs(bpnet.admission_blocking, False)
        self.assertTrue(bpnet.execution.ready)
        self.assertEqual(bpnet_entry["exposure_status"], "target_label_unexposed")
        self.assertNotEqual(bpnet_entry["exposure_status"], "clean_declared")
        self.assertEqual(bpnet.execution.executable_tasks, ("rna_conditioned_atac",))
        for artifact in bpnet_entry["execution"]["evidence_artifacts"]:
            path = CONFIG / artifact["path"]
            self.assertEqual(digest(path), artifact["sha256"])
            self.assertEqual(path.stat().st_size, artifact["size_bytes"])
        environment = bpnet_entry["execution"]["environment_artifact"]
        environment_path = Path(environment["path"])
        self.assertEqual(digest(environment_path), environment["sha256"])
        self.assertEqual(environment_path.stat().st_size, environment["size_bytes"])
        for authority in (
            "code_license",
            "weights_license",
            "derivative_weights_license",
        ):
            artifact = next(
                item
                for item in bpnet_entry["execution"]["evidence_artifacts"]
                if item["role"] == f"model_authority:bpnet:{authority}"
            )
            license_authority = read_json(CONFIG / artifact["path"])
            self.assertEqual(license_authority["model_id"], "bpnet")
            self.assertEqual(license_authority["authority"], authority)
            self.assertIs(license_authority["use_allowed"], True)
        cnn_entry, cnn = model_entry(
            CONFIG / "models/mandatory_baselines.toml", "sequence_cnn_control"
        )
        self.assertIs(cnn.admission_blocking, True)
        self.assertEqual(cnn_entry["license_status"], "UNRESOLVED")
        self.assertEqual(cnn_entry["exposure_status"], "unknown")
        self.assertIn(
            self.authority["models"]["sequence_cnn_control"][
                "bundle_artifacts_sha256"
            ],
            cnn_entry["checkpoint_revision"],
        )
        capability = read_toml(
            CONFIG / "evaluation/rna_conditioned_atac_capabilities.toml"
        )
        self.assertIn("bpnet", capability["profile_fixture_passed_models"])
        self.assertNotIn(
            "sequence_cnn_control", capability["profile_fixture_passed_models"]
        )

    def test_no_terminal_artifact_is_promoted_beyond_development(self) -> None:
        firewall = self.authority["required_false_firewall"]
        for field in (
            "external_evaluation_complete",
            "sealed_evaluation_complete",
            "champion_eligible",
        ):
            self.assertIs(firewall[field], False)
        for binding in self.authority["models"].values():
            self.assertIs(binding["variant_output_admitted"], False)
            self.assertIs(binding["complete_five_fold_five_seed_finalist"], False)
            self.assertIs(binding["exposure_upgrade_to_clean_declared"], False)
        self.assertIs(
            self.authority["models"]["bpnet"][
                "development_profile_prediction_admitted"
            ],
            True,
        )
        self.assertIs(
            self.authority["models"]["sequence_cnn_control"][
                "development_profile_prediction_admitted"
            ],
            True,
        )
        self.assertIs(
            self.authority["models"]["sequence_cnn_control"]["comparison_only"],
            True,
        )
        self.assertIs(
            self.authority["models"]["sequence_cnn_control"]["release_eligible"],
            False,
        )


if __name__ == "__main__":
    unittest.main()
