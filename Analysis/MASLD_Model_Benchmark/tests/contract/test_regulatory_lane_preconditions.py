from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class RegulatoryLanePreconditionTests(unittest.TestCase):
    def test_bpnet_chrombpnet_comparison_binds_frozen_primary_views(self) -> None:
        path = CONFIG / "evaluation/bpnet_chrombpnet_native_head_comparison.json"
        contract = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(contract["status"], "prepared_not_evaluated")
        manifest = ROOT / contract["prediction_manifest"]["path"]
        self.assertEqual(digest(manifest), contract["prediction_manifest"]["sha256"])
        for model_id in ("bpnet", "chrombpnet"):
            model = contract["models"][model_id]
            root = ROOT / model["artifact_path"]
            self.assertEqual(
                digest(root / "ARTIFACTS.json"), model["artifacts_sha256"]
            )
            self.assertIs(model["held_donor_atac_used_for_inference"], False)
            self.assertIs(model["benchmark_metrics_calculated"], False)
        chrombpnet = contract["models"]["chrombpnet"]
        admission = ROOT / chrombpnet["admission_path"]
        self.assertEqual(
            digest(admission / "ARTIFACTS.json"),
            chrombpnet["admission_artifacts_sha256"],
        )
        self.assertEqual(
            chrombpnet["prediction_subdir"], "predictions/nobias_model"
        )
        self.assertEqual(
            contract["evaluation_contract"]["authorized_role_before_finalist_lock"],
            "valid",
        )
        self.assertEqual(
            contract["evaluation_contract"]["test_role_visibility"],
            "withheld_until_finalist_selection_lock",
        )
        self.assertIs(contract["claim_boundary"]["champion_eligible"], False)

    def test_corgi_mapper_contract_is_fold_safe_and_missingness_explicit(self) -> None:
        with (CONFIG / "adaptation/corgi.toml").open("rb") as handle:
            contract = tomllib.load(handle)
        trans = contract["trans_input"]
        self.assertEqual(trans["context_width"], 2891)
        self.assertEqual(trans["fit_unit"], "outer_training_donor_by_cell_state_pseudobulk")
        self.assertIn("never silently set to zero", trans["missing_gene_rule"])
        self.assertIs(trans["sealed_selection_forbidden"], True)
        self.assertIs(
            trans["native_numeric_parity_required_before_native_lane_claim"], True
        )
        self.assertEqual(
            set(trans["arm_contracts"]),
            {
                "released_rank_quantile_from_counts",
                "length_adjusted_tpm_then_released_rank_quantile",
                "learned_count_to_rank_context_encoder",
            },
        )
        learned = trans["arm_contracts"]["learned_count_to_rank_context_encoder"]
        self.assertIs(learned["query_fit_forbidden"], True)
        self.assertIn("no sealed outcomes", learned["loss_scope"])
        self.assertIs(trans["controls"]["held_context_as_normalization_reference"], False)
        self.assertIs(trans["promotion"]["test_role_used_for_selection"], False)
        self.assertIs(trans["promotion"]["sealed_data_used_for_selection"], False)

    def test_context_backbones_remain_fail_closed(self) -> None:
        path = CONFIG / "evaluation/context_backbone_feasibility.json"
        contract = json.loads(path.read_text(encoding="utf-8"))
        self.assertIs(contract["build_authorized"], False)
        self.assertEqual(contract["trigger_result"], "NOT_EVALUATED")
        gates = contract["candidate_gates"]
        self.assertIs(gates["borzoi_ensemble"]["open_backbone_selection_eligible"], False)
        self.assertIs(gates["corgi_regular"]["open_backbone_selection_eligible"], False)
        self.assertIs(gates["corgi_plus"]["comparative_claim_eligible"], False)
        self.assertIs(gates["corgi_plus"]["conditional_component_eligible"], False)
        authorities = {
            "borzoi_ensemble": (
                "artifacts/models/borzoi/checkpoints.json",
                "artifacts/models/borzoi/exposure_audit.json",
            ),
            "corgi_regular": (
                "artifacts/models/corgi/checkpoints.json",
                "artifacts/models/corgi/exposure_audit.json",
            ),
        }
        for model_id, (checkpoint_path, exposure_path) in authorities.items():
            self.assertEqual(
                digest(CONFIG / checkpoint_path),
                gates[model_id]["checkpoint_authority_sha256"],
            )
            self.assertEqual(
                digest(CONFIG / exposure_path),
                gates[model_id]["exposure_authority_sha256"],
            )
        self.assertIs(contract["evaluator_outcomes_exposed"], False)
        self.assertIs(contract["model_training_authorized"], False)


if __name__ == "__main__":
    unittest.main()
