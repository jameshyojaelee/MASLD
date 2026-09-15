from __future__ import annotations

import unittest

from scripts.freeze_corgi_regular_crossfit_head_contract import (
    CorgiHeadContractError,
    validate_contract,
    validate_fold_plan,
)


def _folds() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for evaluation in range(5):
        fit_valid = [fold for fold in range(5) if fold != evaluation]
        rows.append(
            {
                "evaluation_fold": evaluation,
                "evaluation_base_outer_fold": (evaluation - 1) % 5,
                "fit_valid_folds": fit_valid,
                "fit_base_outer_folds": [(fold - 1) % 5 for fold in fit_valid],
            }
        )
    return rows


def _contract() -> dict[str, object]:
    return {
        "schema_version": "masld-bench-corgi-crossfit-profile-head-contract-v1",
        "status": "frozen_development_preflight_only",
        "model_id": "corgi_regular",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "lineage_id": "hepatocyte",
        "adaptation_rung": "new_assay_head",
        "biological_unit": "donor",
        "genomic_block_unit": "chromosome",
        "same_nucleus_topology_required": True,
        "gse296875_masld_diagnosis_claimed": False,
        "histology_or_disease_endpoint_used": False,
        "test_or_sealed_features_used": False,
        "test_or_sealed_outcomes_used": False,
        "champion_or_external_transfer_claim_allowed": False,
        "global_census_revision": False,
        "historical_audit_hashes_modified": False,
        "base_representation": {
            "score_name": "strand_tta_orientation_mean_softplus_regional_sum",
            "primary_context_arm": "actual_released_rank_masked",
            "context_controls": [
                "training_lineage_mean_released_rank",
                "nearest_training_released_rank",
                "shuffled_valid_released_rank",
            ],
            "context_mapper_selected_in_this_contract": False,
            "model_weights_updated": False,
            "base_predictions_are_out_of_fold": True,
        },
        "profile_head": {
            "identity": "one_parameter_geometric_profile_interpolation",
            "alpha_grid": [0.0, 0.25, 0.5, 0.75, 1.0],
            "baseline_alpha": 0.0,
            "raw_corgi_alpha": 1.0,
            "epsilon_total_mass": 1.0e-6,
            "fit_loss": "multinomial_deviance_per_insertion",
            "fit_aggregation": "equal_weight_donor_by_chromosome_units",
            "selection_rule": "one_standard_error_from_minimum_mean_fit_loss_then_smallest_alpha",
            "randomness": "none",
            "fit_context_arm": "actual_released_rank_masked_only",
            "apply_same_selected_alpha_to_context_controls": True,
            "context_arm_specific_fitting_forbidden": True,
            "held_fold_outcome_during_fit_forbidden": True,
        },
        "evaluation": {
            "primary_endpoint": "conditional_multinomial_profile_deviance_skill_over_training_only_mean",
            "secondary_endpoint": "regional_count_spearman",
            "macro_averaging_unit": "donor_by_chromosome",
            "promotion_gate_evaluated": False,
            "development_diagnostic_only": True,
            "full_family_ranking_allowed": False,
        },
        "crossfit_folds": _folds(),
        "terminal_boundaries": {
            "native_numeric_parity_established": False,
            "open_champion_eligible": False,
            "external_validation_present": False,
            "full_donor_lineage_screen_complete": False,
        },
    }


class CorgiCrossfitHeadContractTests(unittest.TestCase):
    def test_valid_fold_plan_is_crossed_and_exhaustive(self) -> None:
        rows = validate_fold_plan(_folds())
        self.assertEqual([row["evaluation_fold"] for row in rows], list(range(5)))

    def test_fold_plan_rejects_held_donor_and_locus_in_fit(self) -> None:
        rows = _folds()
        rows[2]["fit_valid_folds"] = [0, 1, 2, 4]
        with self.assertRaises(CorgiHeadContractError):
            validate_fold_plan(rows)

    def test_contract_rejects_masld_diagnosis_claim(self) -> None:
        value = _contract()
        value["gse296875_masld_diagnosis_claimed"] = True
        with self.assertRaises(CorgiHeadContractError):
            validate_contract(value)

    def test_contract_rejects_context_specific_head_fits(self) -> None:
        value = _contract()
        value["profile_head"]["context_arm_specific_fitting_forbidden"] = False
        with self.assertRaises(CorgiHeadContractError):
            validate_contract(value)

    def test_contract_rejects_sealed_outcome_use(self) -> None:
        value = _contract()
        value["test_or_sealed_outcomes_used"] = True
        with self.assertRaises(CorgiHeadContractError):
            validate_contract(value)

    def test_contract_accepts_single_deterministic_head(self) -> None:
        rows = validate_contract(_contract())
        self.assertEqual(len(rows), 5)


if __name__ == "__main__":
    unittest.main()
