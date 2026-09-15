#!/usr/bin/env python3
"""Audit the frozen negative Regular Corgi new-assay-head disposition."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "masld-bench-corgi-new-assay-head-development-disposition-v1"
NEGATIVE_STATUS = (
    "NEGATIVE_DEVELOPMENT_SMOKE_NO_INCREMENTAL_CORGI_PROFILE_CONTRIBUTION"
)
EXPECTED_CONTEXT_ARMS = {
    "actual_length_adjusted_tpm_rank_masked",
    "actual_released_rank_masked",
    "nearest_training_released_rank",
    "shuffled_valid_released_rank",
    "training_lineage_mean_released_rank",
}


class CorgiNewHeadDispositionError(ValueError):
    """Raised when the frozen result or one of its bound authorities drifts."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise CorgiNewHeadDispositionError(f"JSON authority is not an object: {path}")
    return value


def _resolve_member(root: Path, relative: str) -> Path:
    member = (root / relative).resolve(strict=True)
    try:
        member.relative_to(root)
    except ValueError as error:
        raise CorgiNewHeadDispositionError("bound authority escapes benchmark root") from error
    if not member.is_file() or member.is_symlink():
        raise CorgiNewHeadDispositionError(f"bound authority is not a regular file: {relative}")
    return member


def _assert_false(mapping: dict[str, Any], fields: tuple[str, ...]) -> None:
    for field in fields:
        if mapping.get(field) is not False:
            raise CorgiNewHeadDispositionError(f"claim boundary is not fail-closed: {field}")


def _audit_fit_grid(path: Path) -> None:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 25:
        raise CorgiNewHeadDispositionError("fit loss grid must contain 5 folds x 5 alphas")
    for fold in range(5):
        fold_rows = [row for row in rows if int(row["evaluation_fold"]) == fold]
        if len(fold_rows) != 5:
            raise CorgiNewHeadDispositionError(f"fit loss grid differs for fold {fold}")
        minimum = min(fold_rows, key=lambda row: float(row["mean_fit_loss"]))
        selected = [row for row in fold_rows if int(row["selected"]) == 1]
        if len(selected) != 1 or float(selected[0]["alpha"]) != 0.0:
            raise CorgiNewHeadDispositionError(f"selected alpha differs for fold {fold}")
        if float(minimum["alpha"]) != 0.0:
            raise CorgiNewHeadDispositionError(f"alpha zero is not fit-loss minimum in fold {fold}")


def audit_disposition(root: Path, authority_path: Path) -> dict[str, Any]:
    root = root.resolve(strict=True)
    authority = _load_json(authority_path.resolve(strict=True))
    if authority.get("schema_version") != SCHEMA_VERSION:
        raise CorgiNewHeadDispositionError("disposition schema differs")
    if authority.get("status") != NEGATIVE_STATUS:
        raise CorgiNewHeadDispositionError("negative development status differs")
    if (
        authority.get("model_id") != "corgi_regular"
        or authority.get("task_id") != "rna_conditioned_atac"
        or authority.get("dataset_id") != "gse296875"
        or authority.get("lineage_id") != "hepatocyte"
        or authority.get("adaptation_rung") != "new_assay_head"
    ):
        raise CorgiNewHeadDispositionError("result identity differs")

    hashes = authority.get("frozen_authority_hashes")
    if not isinstance(hashes, dict) or len(hashes) != 11:
        raise CorgiNewHeadDispositionError("frozen authority hash family differs")
    verified: dict[str, str] = {}
    paths: dict[str, Path] = {}
    for relative, expected in sorted(hashes.items()):
        if not isinstance(relative, str) or not isinstance(expected, str):
            raise CorgiNewHeadDispositionError("invalid frozen authority hash record")
        member = _resolve_member(root, relative)
        observed = _digest(member)
        if observed != expected:
            raise CorgiNewHeadDispositionError(f"frozen authority drifted: {relative}")
        verified[relative] = observed
        paths[relative] = member

    fit = _load_json(paths["executions/model-fit-257-21096578/fit/fit_receipt.json"])
    evaluation = _load_json(
        paths["executions/data-eval-258-21096614/evaluation/evaluation.json"]
    )
    _audit_fit_grid(paths["executions/model-fit-257-21096578/fit/fit_loss_grid.tsv"])

    fit_folds = fit.get("folds", [])
    selected = [fold.get("selected_alpha") for fold in fit_folds]
    if len(fit_folds) != 5 or selected != [0.0] * 5:
        raise CorgiNewHeadDispositionError("training-only selected alphas differ")
    if any(
        fold.get("held_fold_atac_signal_values_used_during_fit") is not False
        or fold.get("held_donor_overlap_with_fit") != 0
        or fold.get("held_genomic_fold_overlap_with_fit") != 0
        or fold.get("same_alpha_applied_to_all_context_arms") is not True
        for fold in fit_folds
    ):
        raise CorgiNewHeadDispositionError("training-only crossfit firewall differs")
    if fit.get("benchmark_metrics_computed") is not False:
        raise CorgiNewHeadDispositionError("fit stage computed benchmark metrics")

    independent = authority.get("independent_evaluation", {})
    summaries = evaluation.get("context_arm_summaries", {})
    if set(summaries) != EXPECTED_CONTEXT_ARMS:
        raise CorgiNewHeadDispositionError("context-arm family differs")
    if evaluation.get("all_five_heads_selected_training_mean_alpha_zero") is not True:
        raise CorgiNewHeadDispositionError("independent evaluator alpha-zero result differs")
    if evaluation.get("development_disposition") != independent.get(
        "development_disposition"
    ):
        raise CorgiNewHeadDispositionError("independent negative disposition differs")
    reported = float(evaluation.get("primary_macro_relative_deviance_skill"))
    if reported != float(independent.get("reported_primary_macro_relative_deviance_skill")):
        raise CorgiNewHeadDispositionError("reported numeric residual differs")
    tolerance = float(independent.get("numeric_zero_tolerance"))
    if abs(reported) > tolerance:
        raise CorgiNewHeadDispositionError("numeric residual exceeds frozen zero tolerance")
    if independent.get("reported_numeric_residual_interpretation") != (
        "floating_point_implementation_residual_only_not_biological_gain"
    ):
        raise CorgiNewHeadDispositionError("numeric residual interpretation differs")
    metric_tuples = {
        (
            summary.get("donors"),
            summary.get("units"),
            summary.get("macro_head_deviance_per_insertion"),
            summary.get("macro_regional_spearman"),
            summary.get("macro_relative_deviance_skill"),
        )
        for summary in summaries.values()
    }
    if metric_tuples != {
        (
            independent.get("donors"),
            independent.get("donor_by_chromosome_units_per_context_arm"),
            0.19253147396724268,
            0.6757702964447784,
            reported,
        )
    }:
        raise CorgiNewHeadDispositionError("alpha-zero context arms are not identical")

    boundaries = authority.get("claim_boundaries", {})
    _assert_false(
        boundaries,
        (
            "gse296875_masld_diagnosis_claimed",
            "histology_or_disease_endpoint_used",
            "test_or_sealed_features_or_outcomes_used",
            "promotion_gate_evaluated",
            "full_family_ranking_performed",
            "champion_or_external_transfer_claim_allowed",
            "global_model_census_modified",
            "global_promotion_gate_modified",
            "negative_result_may_be_recast_as_gain",
        ),
    )
    next_rung = authority.get("next_rung", {})
    if (
        next_rung.get("name") != "film_plus_head"
        or next_rung.get("prospective_contract_preparation_allowed") is not True
        or next_rung.get("training_or_prediction_authorized_by_this_disposition") is not False
        or next_rung.get("required_split_identity")
        != "same_five_donor_by_chromosome_crossfit_folds"
        or next_rung.get("required_selection") != "training_only"
        or next_rung.get("required_context_ablations")
        != "identical_to_new_assay_head_rung"
        or next_rung.get("independent_evaluator_required") is not True
    ):
        raise CorgiNewHeadDispositionError("next-rung boundary differs")

    return {
        "schema_version": "masld-bench-corgi-new-head-disposition-audit-receipt-v1",
        "status": "pass_frozen_negative_development_result",
        "authority_sha256": _digest(authority_path),
        "verified_authority_hashes": verified,
        "selected_alpha_by_fold": selected,
        "donors": independent["donors"],
        "donor_by_chromosome_units_per_context_arm": independent[
            "donor_by_chromosome_units_per_context_arm"
        ],
        "numeric_residual_is_biological_gain": False,
        "incremental_corgi_profile_contribution": False,
        "film_plus_head_execution_authorized": False,
        "sealed_or_test_data_opened_by_audit": False,
    }


def main() -> None:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument(
        "--authority",
        type=Path,
        default=default_root
        / "config/artifacts/models/corgi/new_assay_head_development_disposition_20260825.json",
    )
    arguments = parser.parse_args()
    print(
        json.dumps(
            audit_disposition(arguments.root, arguments.authority),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
