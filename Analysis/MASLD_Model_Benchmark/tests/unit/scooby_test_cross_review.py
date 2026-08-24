from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "scooby_cross_review.py"
SPEC = importlib.util.spec_from_file_location("scooby_cross_review_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
reviewed = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(reviewed)


class ScoobyCrossReviewTests(unittest.TestCase):
    def test_project_authorities_pass_with_fail_closed_dispositions(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            receipt = reviewed.review(ROOT, Path(directory) / "receipt")
        self.assertEqual(receipt["status"], "pass")
        self.assertEqual(receipt["rna_conditioned_atac_models"], [])
        self.assertEqual(
            receipt["observed_query_atac_models"],
            ["scooby_epicardioids", "scooby_neurips"],
        )
        self.assertEqual(receipt["checkpoint_forward_allowed"], [])
        self.assertEqual(
            receipt["production_allowed"]["scooby_neurips"],
            "metadata_and_pointer_only",
        )
        self.assertTrue(
            receipt["models"]["scooby_epicardioids"][
                "observed_atac_in_checkpoint_context"
            ]
        )
        self.assertFalse(
            receipt["models"]["scooby_onek1k"][
                "observed_atac_in_checkpoint_context"
            ]
        )

    def test_observed_atac_topology_tamper_fails_closed(self) -> None:
        model_id = "scooby_epicardioids"
        bundle = ROOT / "config" / "artifacts" / "models" / model_id
        checkpoint = reviewed._load_json(bundle / "checkpoints.json")
        crosswalk = reviewed._load_json(bundle / "development_crosswalk.json")
        exposure = reviewed._load_json(bundle / "exposure_audit.json")
        checkpoint = copy.deepcopy(checkpoint)
        checkpoint["input_contract"]["cell_context"] = "50-dimensional RNA latent"
        checkpoint["input_contract"]["topology"] = "same-cell RNA"
        checkpoint["input_contract"]["warning"] = "new context"
        registry = reviewed._model_blocks(
            ROOT / "config" / "models" / "regulatory_sequence.toml"
        )[model_id]
        capability = reviewed._model_blocks(
            ROOT / "config" / "evaluation" / "variant_to_regulation_capabilities.toml"
        )[model_id]
        with self.assertRaisesRegex(
            reviewed.ScoobyCrossReviewError, "observed-ATAC topology"
        ):
            reviewed._review_model(
                model_id, checkpoint, crosswalk, exposure, registry, capability
            )

    def test_outcome_firewall_tamper_fails_closed(self) -> None:
        firewall = {
            "review_inputs": "authority_metadata_only",
            "project_data_read": False,
            "sealed_features_read": False,
            "sealed_labels_or_outcomes_read": False,
            "gse289173_eqtl_or_ieqtl_read": False,
            "gse296875_outcomes_read": False,
            "frozen_117_hotspot_programs_read": False,
            "checkpoint_bytes_read": False,
            "observed_query_atac_results_must_be_separate": True,
            "held_atac_forbidden_in_rna_conditioned_track": True,
            "prediction_hash_committed_before_sealed_label_join": True,
            "sealed_truth": "present",
        }
        with self.assertRaisesRegex(
            reviewed.ScoobyCrossReviewError, "forbidden key"
        ):
            reviewed._validate_outcome_firewall(firewall)


if __name__ == "__main__":
    unittest.main()
