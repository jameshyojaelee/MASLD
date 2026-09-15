#!/usr/bin/env python3
"""Independently score frozen cross-fitted Corgi profile-head predictions."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from scripts.evaluate_corgi_regional_smoke import (
    CONTEXT_ARMS,
    PREDICTION_FIELDS,
    _artifact_member,
    _artifact_metadata,
    _bigwig_sums,
    _load_outcome_paths,
    multinomial_deviance_per_insertion,
    read_tsv,
    spearman_or_zero,
    verify_complete_marker,
    verify_selected_member,
)


SCHEMA_VERSION = "masld-bench-corgi-crossfit-profile-head-evaluation-v1"
FIT_SCHEMA_VERSION = "masld-bench-corgi-crossfit-profile-head-fit-v1"
PRIMARY_CONTEXT_ARM = "actual_released_rank_masked"
FOLDS = tuple(range(5))


class CorgiHeadEvaluationError(RuntimeError):
    """Raised when a fitted head, held outcome, or evaluation boundary drifts."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def _public_donor_hash(donor_id: str) -> str:
    return sha256(f"gse296875-corgi-head-eval-v1\0{donor_id}".encode()).hexdigest()


def relative_deviance_skill(baseline: float, model: float) -> float:
    if not math.isfinite(baseline) or not math.isfinite(model) or baseline < 0 or model < 0:
        raise CorgiHeadEvaluationError("deviance values differ")
    return (baseline - model) / baseline if baseline > 0 else 0.0


def validate_fit_receipt(value: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    if (
        value.get("schema_version") != FIT_SCHEMA_VERSION
        or value.get("status")
        != "pass_training_only_crossfit_head_fit_and_prediction_export"
        or value.get("model_id") != "corgi_regular"
        or value.get("task_id") != "rna_conditioned_atac"
        or value.get("dataset_id") != "gse296875"
        or value.get("lineage_id") != "hepatocyte"
        or value.get("adaptation_rung") != "new_assay_head"
        or value.get("development_atac_outcomes_read_for_training_only_fit") is not True
        or value.get("held_fold_atac_signal_values_used_during_each_head_fit") is not False
        or value.get("histology_or_disease_labels_read") is not False
        or value.get("test_or_sealed_features_or_outcomes_read") is not False
        or value.get("benchmark_metrics_computed") is not False
        or value.get("evaluation_performed") is not False
        or value.get("promotion_gate_evaluated") is not False
        or value.get("champion_or_external_claim_allowed") is not False
        or value.get("global_census_modified") is not False
    ):
        raise CorgiHeadEvaluationError("fit receipt boundary differs")
    folds = value.get("folds")
    if not isinstance(folds, list) or len(folds) != len(FOLDS):
        raise CorgiHeadEvaluationError("fit receipt fold count differs")
    for fold in folds:
        evaluation = int(fold.get("evaluation_fold", -1))
        expected_fit = [item for item in FOLDS if item != evaluation]
        if (
            evaluation not in FOLDS
            or int(fold.get("evaluation_base_outer_fold", -1)) != (evaluation - 1) % 5
            or [int(item) for item in fold.get("fit_valid_folds", [])] != expected_fit
            or int(fold.get("held_donor_overlap_with_fit", -1)) != 0
            or int(fold.get("held_genomic_fold_overlap_with_fit", -1)) != 0
            or fold.get("same_alpha_applied_to_all_context_arms") is not True
            or fold.get("held_fold_atac_signal_values_used_during_fit") is not False
            or fold.get("held_fold_atac_signal_values_used_during_export") is not False
            or fold.get("evaluation_metrics_computed") is not False
            or fold.get("test_or_sealed_features_or_outcomes_read") is not False
        ):
            raise CorgiHeadEvaluationError("fitted fold firewall differs")
    if {int(fold["evaluation_fold"]) for fold in folds} != set(FOLDS):
        raise CorgiHeadEvaluationError("fitted evaluation folds are incomplete")
    return sorted(folds, key=lambda fold: int(fold["evaluation_fold"]))


def _validate_tables(
    records: Sequence[Mapping[str, str]],
    windows: Sequence[Mapping[str, str]],
    evaluation_fold: int,
) -> None:
    if not records or not windows:
        raise CorgiHeadEvaluationError("prediction tables are empty")
    identities: set[tuple[str, str]] = set()
    arms_by_donor: dict[str, set[str]] = {}
    for index, row in enumerate(records):
        if (
            int(row["prediction_index"]) != index
            or row["lineage_id"] != "hepatocyte"
            or int(row["outer_fold"]) != evaluation_fold
            or int(row["donor_fold"]) != evaluation_fold
            or row["outcome_role"] != "valid"
            or row["context_arm"] not in CONTEXT_ARMS
        ):
            raise CorgiHeadEvaluationError("prediction record firewall differs")
        identity = (row["donor_id"], row["context_arm"])
        if identity in identities:
            raise CorgiHeadEvaluationError("prediction record is duplicated")
        identities.add(identity)
        arms_by_donor.setdefault(row["donor_id"], set()).add(row["context_arm"])
    if any(arms != set(CONTEXT_ARMS) for arms in arms_by_donor.values()):
        raise CorgiHeadEvaluationError("prediction context arms are incomplete")
    window_ids: set[tuple[str, str]] = set()
    for row in windows:
        identity = (row["window_id"], row["selection_hash"])
        if (
            int(row["genomic_fold"]) != evaluation_fold
            or identity in window_ids
            or int(row["output_end"]) - int(row["output_start"]) != 1_000
        ):
            raise CorgiHeadEvaluationError("held genomic window differs")
        window_ids.add(identity)


def evaluate(
    *, fit_root: Path, donor_bigwigs: Path, training_bigwigs: Path, output: Path
) -> Mapping[str, Any]:
    if output.exists():
        raise CorgiHeadEvaluationError("refusing to overwrite profile-head evaluation")
    for root in (fit_root, donor_bigwigs, training_bigwigs):
        verify_complete_marker(root)
    fit_metadata = _artifact_metadata(fit_root)
    donor_metadata = _artifact_metadata(donor_bigwigs)
    training_metadata = _artifact_metadata(training_bigwigs)
    if (
        fit_metadata.get("artifact_class") != "corgi_regular_crossfit_profile_head_fit"
        or fit_metadata.get("benchmark_metrics_computed") is not False
        or fit_metadata.get("evaluation_performed") is not False
        or fit_metadata.get("held_fold_atac_signal_values_used_during_each_head_fit") is not False
        or fit_metadata.get("test_or_sealed_features_or_outcomes_read") is not False
        or donor_metadata.get("artifact_class") != "gse296875_deduplicated_tn5_bigwigs"
        or donor_metadata.get("biological_unit") != "donor"
        or training_metadata.get("artifact_class") != "gse296875_training_fold_pseudobulk"
        or training_metadata.get("biological_outer_unit") != "donor"
    ):
        raise CorgiHeadEvaluationError("evaluation source boundary differs")
    fit_member = _artifact_member(fit_root, "fit/fit_receipt.json")
    fit_receipt_path = verify_selected_member(
        fit_root, "fit/fit_receipt.json", str(fit_member["sha256"])
    )
    fit_receipt = json.loads(fit_receipt_path.read_text(encoding="utf-8"))
    fold_receipts = validate_fit_receipt(fit_receipt)
    donor_paths, training_paths = _load_outcome_paths(donor_bigwigs, training_bigwigs)
    unit_rows: list[dict[str, object]] = []
    fold_rows: list[dict[str, object]] = []
    for fold_receipt in fold_receipts:
        evaluation_fold = int(fold_receipt["evaluation_fold"])
        evaluation_base = int(fold_receipt["evaluation_base_outer_fold"])
        fold_root = fit_root / f"fit/fold{evaluation_fold}"
        for filename in (
            "receipt.json",
            "prediction_records.tsv",
            "window_records.tsv",
            "profile_predictions.npz",
        ):
            relative = f"fit/fold{evaluation_fold}/{filename}"
            member = _artifact_member(fit_root, relative)
            verify_selected_member(fit_root, relative, str(member["sha256"]))
        record_fields, records = read_tsv(fold_root / "prediction_records.tsv")
        _, windows = read_tsv(fold_root / "window_records.tsv")
        if record_fields != PREDICTION_FIELDS:
            raise CorgiHeadEvaluationError("prediction record schema differs")
        _validate_tables(records, windows, evaluation_fold)
        with np.load(fold_root / "profile_predictions.npz", allow_pickle=False) as archive:
            if archive.files != ["profile_probability"]:
                raise CorgiHeadEvaluationError("profile prediction archive differs")
            profiles = np.asarray(archive["profile_probability"], dtype=np.float64)
        if (
            profiles.shape != (len(records), len(windows))
            or np.any(~np.isfinite(profiles))
            or np.any(profiles <= 0)
        ):
            raise CorgiHeadEvaluationError("profile prediction matrix differs")
        baseline = np.asarray(_bigwig_sums(training_paths[evaluation_base], windows), dtype=np.float64)
        by_contig: dict[str, list[int]] = {}
        for index, window in enumerate(windows):
            by_contig.setdefault(window["contig"], []).append(index)
        observed_cache: dict[str, np.ndarray] = {}
        for row_index, record in enumerate(records):
            donor = record["donor_id"]
            if donor not in observed_cache:
                path = donor_paths.get((donor, evaluation_fold))
                if path is None:
                    raise CorgiHeadEvaluationError("held-evaluation donor ATAC is unavailable")
                observed_cache[donor] = np.asarray(
                    _bigwig_sums(path, windows), dtype=np.float64
                )
            observed = observed_cache[donor]
            for contig, indices in sorted(by_contig.items()):
                if len(indices) < 3 or float(observed[indices].sum()) <= 0:
                    continue
                probability = profiles[row_index, indices]
                if not math.isclose(float(probability.sum()), 1.0, rel_tol=1.0e-6, abs_tol=1.0e-6):
                    raise CorgiHeadEvaluationError("profile probability is not block-normalized")
                baseline_deviance = multinomial_deviance_per_insertion(
                    observed[indices], baseline[indices]
                )
                head_deviance = multinomial_deviance_per_insertion(
                    observed[indices], probability
                )
                unit_rows.append(
                    {
                        "evaluation_fold": evaluation_fold,
                        "donor_hash": _public_donor_hash(donor),
                        "lineage_id": "hepatocyte",
                        "genomic_block": contig,
                        "context_arm": record["context_arm"],
                        "selected_alpha": fold_receipt["selected_alpha"],
                        "windows": len(indices),
                        "observed_tn5_insertions": format(float(observed[indices].sum()), ".17g"),
                        "training_mean_deviance_per_insertion": format(baseline_deviance, ".17g"),
                        "head_deviance_per_insertion": format(head_deviance, ".17g"),
                        "relative_deviance_skill": format(
                            relative_deviance_skill(baseline_deviance, head_deviance), ".17g"
                        ),
                        "regional_spearman": format(
                            spearman_or_zero(observed[indices], probability), ".17g"
                        ),
                    }
                )
        fold_rows.append(
            {
                "evaluation_fold": evaluation_fold,
                "selected_alpha": float(fold_receipt["selected_alpha"]),
                "held_donors": int(fold_receipt["held_donors"]),
                "windows": len(windows),
                "genomic_blocks": len(by_contig),
            }
        )
    if not unit_rows:
        raise CorgiHeadEvaluationError("no held donor by block units were scoreable")
    output.mkdir(parents=True, exist_ok=False)
    fields = (
        "evaluation_fold",
        "donor_hash",
        "lineage_id",
        "genomic_block",
        "context_arm",
        "selected_alpha",
        "windows",
        "observed_tn5_insertions",
        "training_mean_deviance_per_insertion",
        "head_deviance_per_insertion",
        "relative_deviance_skill",
        "regional_spearman",
    )
    _write_tsv(output / "unit_metrics.tsv", fields, unit_rows)
    summaries: dict[str, dict[str, float | int]] = {}
    for arm in CONTEXT_ARMS:
        selected = [row for row in unit_rows if row["context_arm"] == arm]
        if not selected:
            raise CorgiHeadEvaluationError("context arm has no held units")
        summaries[arm] = {
            "units": len(selected),
            "donors": len({row["donor_hash"] for row in selected}),
            "macro_relative_deviance_skill": float(
                np.mean([float(row["relative_deviance_skill"]) for row in selected])
            ),
            "macro_regional_spearman": float(
                np.mean([float(row["regional_spearman"]) for row in selected])
            ),
            "macro_head_deviance_per_insertion": float(
                np.mean([float(row["head_deviance_per_insertion"]) for row in selected])
            ),
        }
    primary = summaries[PRIMARY_CONTEXT_ARM]
    all_zero_alpha = all(float(row["selected_alpha"]) == 0.0 for row in fold_rows)
    result = {
        "schema_version": SCHEMA_VERSION,
        "status": "pass_independent_development_smoke_evaluation",
        "model_id": "corgi_regular",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "lineage_id": "hepatocyte",
        "adaptation_rung": "new_assay_head",
        "fit_artifacts_sha256": digest(fit_root / "ARTIFACTS.json"),
        "donor_bigwigs_artifacts_sha256": digest(donor_bigwigs / "ARTIFACTS.json"),
        "training_bigwigs_artifacts_sha256": digest(training_bigwigs / "ARTIFACTS.json"),
        "folds": fold_rows,
        "context_arm_summaries": summaries,
        "all_five_heads_selected_training_mean_alpha_zero": all_zero_alpha,
        "primary_macro_relative_deviance_skill": float(
            primary["macro_relative_deviance_skill"]
        ),
        "development_disposition": (
            "selected_head_is_exact_training_mean_baseline_no_incremental_corgi_profile_contribution"
            if all_zero_alpha
            and abs(float(primary["macro_relative_deviance_skill"])) <= 1.0e-6
            else "nonzero_head_requires_bounded_diagnostic_interpretation"
        ),
        "biological_unit": "donor",
        "genomic_block_unit": "chromosome",
        "macro_averaging_unit": "donor_by_chromosome",
        "development_outcomes_read": True,
        "histology_or_disease_labels_read": False,
        "test_or_sealed_features_or_outcomes_read": False,
        "promotion_gate_evaluated": False,
        "full_family_ranking_performed": False,
        "champion_or_external_claim_allowed": False,
        "scope": "existing_128_tile_per_fold_hepatocyte_development_smoke_only",
    }
    (output / "evaluation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fit-root", type=Path, required=True)
    parser.add_argument("--donor-bigwigs", type=Path, required=True)
    parser.add_argument("--training-bigwigs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(
        fit_root=args.fit_root.resolve(strict=True),
        donor_bigwigs=args.donor_bigwigs.resolve(strict=True),
        training_bigwigs=args.training_bigwigs.resolve(strict=True),
        output=args.output.resolve(strict=False),
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
