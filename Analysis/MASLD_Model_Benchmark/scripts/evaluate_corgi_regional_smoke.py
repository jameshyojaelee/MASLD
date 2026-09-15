#!/usr/bin/env python3
"""Evaluate frozen Corgi regional predictions on held-validation donor ATAC."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCHEMA_VERSION = "masld-bench-corgi-regional-smoke-evaluation-v1"
SCORE_NAME = "strand_tta_orientation_mean_softplus_regional_sum"
CONTEXT_ARMS = (
    "actual_released_rank_masked",
    "actual_length_adjusted_tpm_rank_masked",
    "training_lineage_mean_released_rank",
    "nearest_training_released_rank",
    "shuffled_valid_released_rank",
)
PREDICTION_FIELDS = (
    "prediction_index",
    "unit_index",
    "donor_id",
    "lineage_id",
    "outer_fold",
    "donor_fold",
    "context_arm",
    "context_source_unit",
    "outcome_role",
)
WINDOW_REQUIRED_FIELDS = (
    "contig",
    "output_start",
    "output_end",
    "window_id",
    "genomic_fold",
    "selection_hash",
)
BIGWIG_FIELDS = (
    "donor_id",
    "outer_fold",
    "lineage_id",
    "analysis_role",
    "nuclei",
    "unique_fragments",
    "read_support",
    "tn5_insertions",
    "nonzero_positions",
    "max_pending_positions",
    "path",
    "size_bytes",
    "sha256",
)
TRAINING_FIELDS = (
    "donor_test_fold",
    "donor_valid_fold",
    "donor_train_folds",
    "lineage_id",
    "donors",
    "nuclei",
    "unique_fragments",
    "read_support",
    "tn5_insertions",
    "nonzero_positions",
    "max_pending_positions",
    "fragment_path",
    "fragment_size_bytes",
    "fragment_sha256",
    "bigwig_path",
    "bigwig_size_bytes",
    "bigwig_sha256",
)


class CorgiEvaluationError(RuntimeError):
    """Raised when the development outcome separation or alignment differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise CorgiEvaluationError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def write_tsv(
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


def _artifact_metadata(root: Path) -> Mapping[str, Any]:
    manifest = json.loads((root / "ARTIFACTS.json").read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "masld-bench-artifacts-v1":
        raise CorgiEvaluationError(f"artifact manifest schema differs: {root}")
    metadata = manifest.get("metadata")
    if not isinstance(metadata, dict):
        raise CorgiEvaluationError(f"artifact metadata differs: {root}")
    return metadata


def verify_complete_marker(root: Path) -> None:
    manifest = json.loads((root / "ARTIFACTS.json").read_text(encoding="utf-8"))
    complete_path = root / "COMPLETE"
    if complete_path.is_symlink() or not complete_path.is_file():
        raise CorgiEvaluationError(f"frozen completion marker differs: {root}")
    complete = json.loads(complete_path.read_text(encoding="utf-8"))
    if (
        complete.get("schema_version") != "masld-bench-complete-v1"
        or complete.get("manifest_sha256") != digest(root / "ARTIFACTS.json")
        or complete.get("artifact_count") != len(manifest.get("artifacts", []))
    ):
        raise CorgiEvaluationError(f"frozen completion receipt differs: {root}")


def _artifact_member(root: Path, relative: str) -> Mapping[str, Any]:
    manifest = json.loads((root / "ARTIFACTS.json").read_text(encoding="utf-8"))
    matches = [item for item in manifest.get("artifacts", []) if item.get("path") == relative]
    if len(matches) != 1:
        raise CorgiEvaluationError(f"artifact member is not uniquely frozen: {relative}")
    return matches[0]


def verify_selected_member(root: Path, relative: str, expected_sha256: str) -> Path:
    member = _artifact_member(root, relative)
    path = (root / relative).resolve(strict=True)
    try:
        path.relative_to(root.resolve(strict=True))
    except ValueError as error:
        raise CorgiEvaluationError("artifact member escapes frozen root") from error
    if member.get("sha256") != expected_sha256 or digest(path) != expected_sha256:
        raise CorgiEvaluationError(f"selected frozen artifact differs: {relative}")
    if member.get("size_bytes") != path.stat().st_size:
        raise CorgiEvaluationError(f"selected artifact size differs: {relative}")
    return path


def smooth_probability(values: Any, epsilon: float = 1.0e-6) -> Any:
    import numpy as np

    vector = np.asarray(values, dtype=np.float64)
    if vector.ndim != 1 or vector.size < 3 or np.any(~np.isfinite(vector)) or np.any(vector < 0):
        raise CorgiEvaluationError("score vector differs")
    total = float(vector.sum())
    probability = (
        np.full(vector.shape, 1.0 / vector.size, dtype=np.float64)
        if total <= 0
        else vector / total
    )
    return (1.0 - epsilon) * probability + epsilon / vector.size


def multinomial_deviance_per_insertion(observed: Any, predicted: Any) -> float:
    import numpy as np

    truth = np.asarray(observed, dtype=np.float64)
    if truth.ndim != 1 or np.any(~np.isfinite(truth)) or np.any(truth < 0):
        raise CorgiEvaluationError("observed vector differs")
    total = float(truth.sum())
    if total <= 0:
        raise CorgiEvaluationError("observed vector has no insertions")
    probability = smooth_probability(predicted)
    if probability.shape != truth.shape:
        raise CorgiEvaluationError("observed and predicted vectors differ")
    positive = truth > 0
    expected = total * probability[positive]
    return float(2.0 * np.sum(truth[positive] * np.log(truth[positive] / expected)) / total)


def rank_average(values: Any) -> Any:
    import numpy as np

    vector = np.asarray(values, dtype=np.float64)
    order = np.argsort(vector, kind="mergesort")
    ranks = np.empty(vector.size, dtype=np.float64)
    start = 0
    while start < vector.size:
        end = start + 1
        while end < vector.size and vector[order[end]] == vector[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + end - 1) / 2.0 + 1.0
        start = end
    return ranks


def spearman_or_zero(observed: Any, predicted: Any) -> float:
    import numpy as np

    truth = np.asarray(observed, dtype=np.float64)
    score = np.asarray(predicted, dtype=np.float64)
    if truth.shape != score.shape or truth.ndim != 1 or truth.size < 3:
        raise CorgiEvaluationError("Spearman vectors differ")
    left, right = rank_average(truth), rank_average(score)
    if np.all(left == left[0]) or np.all(right == right[0]):
        return 0.0
    value = float(np.corrcoef(left, right)[0, 1])
    return value if math.isfinite(value) else 0.0


def validate_prediction_tables(
    records: Sequence[Mapping[str, str]],
    windows: Sequence[Mapping[str, str]],
    receipt: Mapping[str, Any],
) -> None:
    outer_fold = int(receipt["mapper_outer_fold"])
    valid_fold = (outer_fold + 1) % 5
    if (
        receipt.get("schema_version") != "masld-bench-corgi-outcome-aligned-prediction-v3"
        or receipt.get("status") != "pass_outcome_free_prediction"
        or receipt.get("outcome_role") != "valid"
        or receipt.get("lineage_id") != "hepatocyte"
        or receipt.get("valid_donor_fold") != valid_fold
        or receipt.get("valid_genomic_fold") != valid_fold
        or receipt.get("held_ATAC_or_other_outcomes_used") is not False
        or receipt.get("test_or_sealed_features_or_labels_read") is not False
        or receipt.get("model_fitted_or_adapted") is not False
        or receipt.get("context_arms") != list(CONTEXT_ARMS)
        or SCORE_NAME not in receipt.get("regional_score_names", [])
    ):
        raise CorgiEvaluationError("Corgi prediction receipt differs")
    if not records or len(records) != int(receipt["valid_donors"]) * len(CONTEXT_ARMS):
        raise CorgiEvaluationError("prediction record count differs")
    indices = [int(row["prediction_index"]) for row in records]
    if indices != list(range(len(records))):
        raise CorgiEvaluationError("prediction indices differ")
    identities: set[tuple[str, str]] = set()
    arms_by_donor: dict[str, set[str]] = {}
    for row in records:
        if (
            row["lineage_id"] != "hepatocyte"
            or int(row["outer_fold"]) != valid_fold
            or int(row["donor_fold"]) != valid_fold
            or row["outcome_role"] != "valid"
            or row["context_arm"] not in CONTEXT_ARMS
        ):
            raise CorgiEvaluationError("prediction record firewall differs")
        identity = (row["donor_id"], row["context_arm"])
        if identity in identities:
            raise CorgiEvaluationError("prediction record is duplicated")
        identities.add(identity)
        arms_by_donor.setdefault(row["donor_id"], set()).add(row["context_arm"])
    if any(arms != set(CONTEXT_ARMS) for arms in arms_by_donor.values()):
        raise CorgiEvaluationError("context arms are incomplete")
    if not windows or len(windows) != int(receipt["scoreable_windows"]):
        raise CorgiEvaluationError("window record count differs")
    window_ids: set[tuple[str, str]] = set()
    for row in windows:
        if int(row["genomic_fold"]) != valid_fold:
            raise CorgiEvaluationError("held genomic fold differs")
        identity = (row["window_id"], row["selection_hash"])
        if identity in window_ids or int(row["output_end"]) - int(row["output_start"]) != 1_000:
            raise CorgiEvaluationError("window identity or width differs")
        window_ids.add(identity)


def _bigwig_sums(path: Path, windows: Sequence[Mapping[str, str]]) -> list[float]:
    import pyBigWig

    values: list[float] = []
    with pyBigWig.open(str(path)) as handle:
        chroms = handle.chroms()
        for row in windows:
            contig = row["contig"]
            start, end = int(row["output_start"]), int(row["output_end"])
            if contig not in chroms or start < 0 or end > int(chroms[contig]) or end <= start:
                raise CorgiEvaluationError("outcome interval is outside the BigWig reference")
            value = handle.stats(contig, start, end, type="sum", exact=True)[0]
            values.append(0.0 if value is None else float(value))
    if any(not math.isfinite(value) or value < 0 for value in values):
        raise CorgiEvaluationError("outcome BigWig contains invalid signal")
    return values


def _public_donor_hash(donor_id: str) -> str:
    return sha256(f"gse296875-corgi-regional-smoke-v1\0{donor_id}".encode()).hexdigest()


def _load_outcome_paths(
    donor_bigwigs: Path,
    training_bigwigs: Path,
) -> tuple[dict[tuple[str, int], Path], dict[int, Path]]:
    for root, relative in (
        (donor_bigwigs, "bigwig_manifest.tsv"),
        (training_bigwigs, "training_pseudobulk_manifest.tsv"),
    ):
        member = _artifact_member(root, relative)
        verify_selected_member(root, relative, str(member["sha256"]))
    donor_fields, donor_rows = read_tsv(donor_bigwigs / "bigwig_manifest.tsv")
    training_fields, training_rows = read_tsv(training_bigwigs / "training_pseudobulk_manifest.tsv")
    if donor_fields != BIGWIG_FIELDS or training_fields != TRAINING_FIELDS:
        raise CorgiEvaluationError("outcome manifest schema differs")
    donors: dict[tuple[str, int], Path] = {}
    for row in donor_rows:
        if row["lineage_id"] != "hepatocyte":
            continue
        key = (row["donor_id"], int(row["outer_fold"]))
        if key in donors:
            raise CorgiEvaluationError("donor ATAC outcome is duplicated")
        donors[key] = verify_selected_member(donor_bigwigs, row["path"], row["sha256"])
    training: dict[int, Path] = {}
    for row in training_rows:
        if row["lineage_id"] != "hepatocyte":
            continue
        outer_fold, valid_fold = int(row["donor_test_fold"]), int(row["donor_valid_fold"])
        if valid_fold != (outer_fold + 1) % 5 or outer_fold in training:
            raise CorgiEvaluationError("training-only ATAC fold contract differs")
        training[outer_fold] = verify_selected_member(
            training_bigwigs, row["bigwig_path"], row["bigwig_sha256"]
        )
    if len(training) != 5:
        raise CorgiEvaluationError("training-only ATAC baselines are incomplete")
    return donors, training


def evaluate(
    prediction_root: Path,
    donor_bigwigs: Path,
    training_bigwigs: Path,
    output: Path,
) -> Mapping[str, Any]:
    import numpy as np

    prediction_metadata = _artifact_metadata(prediction_root)
    donor_metadata = _artifact_metadata(donor_bigwigs)
    training_metadata = _artifact_metadata(training_bigwigs)
    for root in (prediction_root, donor_bigwigs, training_bigwigs):
        verify_complete_marker(root)
    contract_member = _artifact_member(training_bigwigs, "contract.json")
    verify_selected_member(training_bigwigs, "contract.json", str(contract_member["sha256"]))
    training_contract = json.loads(
        (training_bigwigs / "contract.json").read_text(encoding="utf-8")
    )
    if (
        prediction_metadata.get("artifact_class") != "Corgi_outcome_aligned_tile_smoke_bundle_v3"
        or prediction_metadata.get("outcomes_read") is not False
        or prediction_metadata.get("logical_tasks") != 5
        or donor_metadata.get("artifact_class") != "gse296875_deduplicated_tn5_bigwigs"
        or donor_metadata.get("dataset_id") != "gse296875"
        or donor_metadata.get("biological_unit") != "donor"
        or training_metadata.get("artifact_class") != "gse296875_training_fold_pseudobulk"
        or training_metadata.get("dataset_id") != "gse296875"
        or training_metadata.get("biological_outer_unit") != "donor"
        or training_contract.get("schema_version") != "masld-bench-training-fold-pseudobulk-v1"
        or training_contract.get("donor_valid_fold_rule") != "test_plus_1_mod_5"
        or training_contract.get("held_donor_test_used") is not False
        or training_contract.get("held_donor_validation_used") is not False
    ):
        raise CorgiEvaluationError("evaluation source firewall differs")
    donor_paths, training_paths = _load_outcome_paths(donor_bigwigs, training_bigwigs)
    unit_rows: list[dict[str, object]] = []
    fold_receipts: list[dict[str, object]] = []
    for outer_fold in range(5):
        fold_root = prediction_root / f"fold{outer_fold}"
        for filename in (
            "receipt.json",
            "prediction_records.tsv",
            "window_records.tsv",
            "regional_predictions.npz",
        ):
            relative = f"fold{outer_fold}/{filename}"
            member = _artifact_member(prediction_root, relative)
            verify_selected_member(prediction_root, relative, str(member["sha256"]))
        receipt = json.loads((fold_root / "receipt.json").read_text(encoding="utf-8"))
        record_fields, records = read_tsv(fold_root / "prediction_records.tsv")
        window_fields, windows = read_tsv(fold_root / "window_records.tsv")
        if record_fields != PREDICTION_FIELDS or any(
            field not in window_fields for field in WINDOW_REQUIRED_FIELDS
        ):
            raise CorgiEvaluationError("prediction table schema differs")
        validate_prediction_tables(records, windows, receipt)
        archive = np.load(fold_root / "regional_predictions.npz", allow_pickle=False)
        if set(archive.files) != set(receipt["regional_score_names"]):
            raise CorgiEvaluationError("regional prediction archive keys differ")
        scores = np.asarray(archive[SCORE_NAME], dtype=np.float64)
        if scores.shape != (len(records), len(windows)) or np.any(~np.isfinite(scores)) or np.any(scores < 0):
            raise CorgiEvaluationError("regional prediction matrix differs")
        valid_fold = (outer_fold + 1) % 5
        baseline = np.asarray(_bigwig_sums(training_paths[outer_fold], windows), dtype=np.float64)
        indices_by_contig: dict[str, list[int]] = {}
        for index, window in enumerate(windows):
            indices_by_contig.setdefault(window["contig"], []).append(index)
        observed_cache: dict[str, Any] = {}
        for row_index, record in enumerate(records):
            donor = record["donor_id"]
            if donor not in observed_cache:
                path = donor_paths.get((donor, valid_fold))
                if path is None:
                    raise CorgiEvaluationError("held-validation donor ATAC is unavailable")
                observed_cache[donor] = np.asarray(_bigwig_sums(path, windows), dtype=np.float64)
            observed = observed_cache[donor]
            for contig, indices in sorted(indices_by_contig.items()):
                if len(indices) < 3 or float(observed[indices].sum()) <= 0:
                    continue
                baseline_deviance = multinomial_deviance_per_insertion(
                    observed[indices], baseline[indices]
                )
                model_deviance = multinomial_deviance_per_insertion(
                    observed[indices], scores[row_index, indices]
                )
                skill = (
                    (baseline_deviance - model_deviance) / baseline_deviance
                    if baseline_deviance > 0
                    else 0.0
                )
                unit_rows.append(
                    {
                        "outer_fold": outer_fold,
                        "valid_donor_fold": valid_fold,
                        "donor_hash": _public_donor_hash(donor),
                        "lineage_id": "hepatocyte",
                        "genomic_block": contig,
                        "context_arm": record["context_arm"],
                        "windows": len(indices),
                        "observed_tn5_insertions": format(float(observed[indices].sum()), ".17g"),
                        "training_mean_deviance_per_insertion": format(baseline_deviance, ".17g"),
                        "model_deviance_per_insertion": format(model_deviance, ".17g"),
                        "relative_deviance_skill": format(skill, ".17g"),
                        "regional_spearman": format(
                            spearman_or_zero(observed[indices], scores[row_index, indices]),
                            ".17g",
                        ),
                    }
                )
        fold_receipts.append(
            {
                "outer_fold": outer_fold,
                "valid_donor_fold": valid_fold,
                "donors": len(observed_cache),
                "windows": len(windows),
                "genomic_blocks": len(indices_by_contig),
            }
        )
    if not unit_rows:
        raise CorgiEvaluationError("no donor by genomic-block outcome units were scoreable")
    metric_fields = (
        "outer_fold",
        "valid_donor_fold",
        "donor_hash",
        "lineage_id",
        "genomic_block",
        "context_arm",
        "windows",
        "observed_tn5_insertions",
        "training_mean_deviance_per_insertion",
        "model_deviance_per_insertion",
        "relative_deviance_skill",
        "regional_spearman",
    )
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "unit_metrics.tsv", metric_fields, unit_rows)
    summaries: dict[str, dict[str, float | int]] = {}
    for arm in CONTEXT_ARMS:
        selected = [row for row in unit_rows if row["context_arm"] == arm]
        if not selected:
            raise CorgiEvaluationError("context arm has no scoreable units")
        summaries[arm] = {
            "units": len(selected),
            "donors": len({row["donor_hash"] for row in selected}),
            "macro_relative_deviance_skill": float(
                np.mean([float(row["relative_deviance_skill"]) for row in selected])
            ),
            "macro_regional_spearman": float(
                np.mean([float(row["regional_spearman"]) for row in selected])
            ),
            "macro_model_deviance_per_insertion": float(
                np.mean([float(row["model_deviance_per_insertion"]) for row in selected])
            ),
        }
    actual = summaries[CONTEXT_ARMS[0]]
    contrasts = {
        arm: {
            "actual_minus_comparator_relative_deviance_skill": float(
                actual["macro_relative_deviance_skill"]
                - summaries[arm]["macro_relative_deviance_skill"]
            ),
            "actual_minus_comparator_regional_spearman": float(
                actual["macro_regional_spearman"] - summaries[arm]["macro_regional_spearman"]
            ),
        }
        for arm in CONTEXT_ARMS[1:]
    }
    result = {
        "schema_version": SCHEMA_VERSION,
        "status": "pass_development_smoke_evaluation",
        "dataset_id": "gse296875",
        "model_id": "corgi_regular",
        "lineage_id": "hepatocyte",
        "prediction_artifacts_sha256": digest(prediction_root / "ARTIFACTS.json"),
        "donor_bigwig_artifacts_sha256": digest(donor_bigwigs / "ARTIFACTS.json"),
        "training_bigwig_artifacts_sha256": digest(training_bigwigs / "ARTIFACTS.json"),
        "score_name": SCORE_NAME,
        "folds": fold_receipts,
        "context_arm_summaries": summaries,
        "actual_context_contrasts": contrasts,
        "biological_unit": "donor",
        "genomic_block_unit": "chromosome",
        "macro_averaging_unit": "donor_by_chromosome",
        "development_outcomes_read": True,
        "held_validation_atac_read": True,
        "held_test_or_sealed_outcomes_read": False,
        "histology_or_disease_labels_read": False,
        "champion_claim_allowed": False,
        "promotion_gate_evaluated": False,
        "scope": "128_tile_per_fold_regional_development_smoke_only",
    }
    (output / "evaluation.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--donor-bigwigs", required=True, type=Path)
    parser.add_argument("--training-bigwigs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = evaluate(
        args.predictions.resolve(strict=True),
        args.donor_bigwigs.resolve(strict=True),
        args.training_bigwigs.resolve(strict=True),
        args.output.resolve(strict=False),
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
