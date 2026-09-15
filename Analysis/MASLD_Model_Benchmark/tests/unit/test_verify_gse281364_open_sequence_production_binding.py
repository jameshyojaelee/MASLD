from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

from masld_bench.artifacts import freeze_tree
from scripts import verify_gse281364_open_sequence_production_binding as binding


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class ProductionBindingTests(unittest.TestCase):
    def build_fixture(self, temporary: str) -> tuple[Path, Path, str]:
        root = Path(temporary) / "benchmark"
        (root / "executions").mkdir(parents=True)
        for index, relative in enumerate(binding.REQUIRED_SOURCES):
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"validated-source-{index}\n", encoding="utf-8")
        validation = root / "executions/model-check-224-12345"
        (validation / "validation").mkdir(parents=True)
        source_rows = [
            f"{digest(root / relative)}  {relative}"
            for relative in binding.REQUIRED_SOURCES
        ]
        (validation / "source.sha256").write_text(
            "\n".join(source_rows) + "\n", encoding="utf-8"
        )
        wrapper_sha256 = digest(root / binding.PRODUCTION_WRAPPER)
        config_sha256 = digest(root / binding.CONFIG)
        source_sha256 = digest(validation / "source.sha256")
        receipt = {
            "schema_version": "masld-bench-gse281364-open-sequence-validation-v2",
            "status": "passed_bound_production_wrapper",
            "campaign_scope": "partial_common_head_campaign",
            "mandatory_baselines_complete": False,
            "mandatory_baselines_missing": list(binding.MISSING_BASELINES),
            "shortlist_blocked": True,
            "finalist_claim_blocked": True,
            "complementarity_blocked": True,
            "conditional_trigger_blocked": True,
            "outcomes_read": False,
            "model_fit": False,
            "predictions_generated": False,
            "production_wrapper_submitted": False,
            "production_wrapper_sha256": wrapper_sha256,
            "config_sha256": config_sha256,
            "source_manifest_sha256": source_sha256,
        }
        (validation / "validation/receipt.json").write_text(
            json.dumps(receipt, sort_keys=True) + "\n", encoding="utf-8"
        )
        artifacts_sha256 = freeze_tree(
            validation,
            {
                "artifact_class": binding.VALIDATION_ARTIFACT_CLASS,
                **{key: value for key, value in receipt.items() if key != "schema_version"},
            },
        )
        return root, validation, artifacts_sha256

    def test_exact_frozen_binding_passes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, validation, expected = self.build_fixture(temporary)
            receipt = binding.verify_binding(
                root=root,
                validation_root=validation,
                validation_sha256=expected,
            )
            self.assertEqual(receipt["status"], "pass_exact_frozen_validation_binding")
            self.assertEqual(receipt["validated_source_count"], 9)
            self.assertFalse(receipt["mandatory_baselines_complete"])
            self.assertTrue(receipt["shortlist_blocked"])

    def test_wrong_validation_sha_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, validation, _ = self.build_fixture(temporary)
            with self.assertRaises(binding.BindingError):
                binding.verify_binding(
                    root=root,
                    validation_root=validation,
                    validation_sha256="0" * 64,
                )

    def test_source_drift_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, validation, expected = self.build_fixture(temporary)
            drifted = root / binding.CONFIG
            drifted.write_text("drifted\n", encoding="utf-8")
            with self.assertRaisesRegex(binding.BindingError, "drifted"):
                binding.verify_binding(
                    root=root,
                    validation_root=validation,
                    validation_sha256=expected,
                )

    def test_noncanonical_validation_location_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root, validation, expected = self.build_fixture(temporary)
            moved = root / "executions/not-v2-validation"
            validation.rename(moved)
            with self.assertRaises(binding.BindingError):
                binding.verify_binding(
                    root=root,
                    validation_root=moved,
                    validation_sha256=expected,
                )


if __name__ == "__main__":
    unittest.main()
