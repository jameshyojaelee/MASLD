#!/usr/bin/env python3
"""Independently score sequence/accessibility models on held donor profiles."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


EPSILON_MIXTURE = 1.0e-6
DEVELOPMENT_EVALUATION_ROLE = "valid"
LEGACY_PREDICTION_MANIFEST_FIELDS = (
    "model_id",
    "root",
    "artifacts_sha256",
    "prediction_subdir",
)
PREDICTION_MANIFEST_FIELDS = (
    "candidate_id",
    "root_model_id",
    "root",
    "artifacts_sha256",
    "prediction_subdir",
    "seed",
)
COUNT_REQUIRED_FIELDS = (
    "window_id",
    "selection_hash",
    "contig",
    "output_start",
    "output_end",
    "genomic_fold",
    "role",
    "ccre_class",
)


class SequenceEvaluationError(RuntimeError):
    """Raised when prediction/outcome separation or alignment is violated."""


def validate_evaluation_role(value: str) -> str:
    """Fail closed until a separate finalist selection record authorizes test evaluation."""
    if value != DEVELOPMENT_EVALUATION_ROLE:
        raise SequenceEvaluationError(
            "only the validation role is authorized before finalist selection lock"
        )
    return value


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _decode(values: Any) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SequenceEvaluationError(f"TSV has no header: {path}")
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


def smooth_distribution(values: Any, epsilon: float = EPSILON_MIXTURE) -> Any:
    import numpy as np

    vector = np.asarray(values, dtype=np.float64)
    if vector.size == 0 or np.any(~np.isfinite(vector)) or np.any(vector < 0):
        raise SequenceEvaluationError("distribution contains invalid values")
    total = float(vector.sum())
    if total <= 0:
        probability = np.full(vector.shape, 1.0 / vector.size, dtype=np.float64)
    else:
        probability = vector / total
    return (1.0 - epsilon) * probability + epsilon / vector.size


def multinomial_deviance_per_insertion(observed: Any, predicted: Any) -> float:
    import numpy as np

    truth = np.asarray(observed, dtype=np.float64).reshape(-1)
    score = np.asarray(predicted, dtype=np.float64).reshape(-1)
    if truth.shape != score.shape or np.any(~np.isfinite(truth)) or np.any(truth < 0):
        raise SequenceEvaluationError("multinomial inputs differ")
    total = float(truth.sum())
    if total <= 0:
        raise SequenceEvaluationError("multinomial truth has zero mass")
    probability = smooth_distribution(score)
    positive = truth > 0
    expected = total * probability[positive]
    return float(2.0 * np.sum(truth[positive] * np.log(truth[positive] / expected)) / total)


def within_window_deviance_per_insertion(observed: Any, predicted: Any) -> float:
    import numpy as np

    truth = np.asarray(observed, dtype=np.float64)
    score = np.asarray(predicted, dtype=np.float64)
    if truth.ndim != 2 or truth.shape != score.shape:
        raise SequenceEvaluationError("within-window profile shapes differ")
    row_totals = truth.sum(axis=1)
    eligible = row_totals > 0
    if not np.any(eligible):
        raise SequenceEvaluationError("within-window truth has zero mass")
    truth = truth[eligible]
    score = score[eligible]
    row_totals = row_totals[eligible]
    score_totals = score.sum(axis=1)
    probability = np.empty_like(score, dtype=np.float64)
    positive_score = score_totals > 0
    probability[positive_score] = score[positive_score] / score_totals[positive_score, None]
    probability[~positive_score] = 1.0 / score.shape[1]
    probability = (1.0 - EPSILON_MIXTURE) * probability + (
        EPSILON_MIXTURE / score.shape[1]
    )
    positive = truth > 0
    expected = row_totals[:, None] * probability
    total_deviance = 2.0 * np.sum(
        truth[positive] * np.log(truth[positive] / expected[positive])
    )
    return float(total_deviance / row_totals.sum())


def spearman_or_zero(observed: Any, predicted: Any) -> float:
    import numpy as np
    from scipy.stats import spearmanr

    truth = np.asarray(observed, dtype=np.float64)
    score = np.asarray(predicted, dtype=np.float64)
    if truth.shape != score.shape or truth.size < 3:
        raise SequenceEvaluationError("Spearman inputs differ")
    if np.all(truth == truth[0]) or np.all(score == score[0]):
        return 0.0
    value = float(spearmanr(truth, score).statistic)
    return 0.0 if not math.isfinite(value) else value


def fisher_mean(values: Sequence[float]) -> float:
    import numpy as np

    if not values:
        raise SequenceEvaluationError("Fisher mean has no values")
    clipped = np.clip(np.asarray(values, dtype=np.float64), -0.999999, 0.999999)
    return float(np.tanh(np.mean(np.arctanh(clipped))))


def _donor_hash(namespace: str, donor_id: str) -> str:
    return sha256(f"{namespace}\0donor\0{donor_id}".encode()).hexdigest()


def _block_hash(namespace: str, contig: str) -> str:
    return sha256(f"{namespace}\0block\0{contig}".encode()).hexdigest()


def validate_evaluator_input_metadata(
    *,
    root: Path,
    expected_artifacts_sha256: str,
    artifact_class: str,
    dataset_id: str,
    split_id: str,
    lineage_id: str,
) -> dict[str, Any]:
    artifacts_path = root / "ARTIFACTS.json"
    if sha256_file(artifacts_path) != expected_artifacts_sha256:
        raise SequenceEvaluationError(f"evaluator input ARTIFACTS differs: {root}")
    metadata = json.loads(artifacts_path.read_text(encoding="utf-8")).get(
        "metadata", {}
    )
    if (
        metadata.get("artifact_class") != artifact_class
        or metadata.get("dataset_id") != dataset_id
        or metadata.get("split_id") != split_id
        or metadata.get("lineage_id") != lineage_id
        or metadata.get("biological_unit") != "donor"
        or metadata.get("status") != "passed"
        or metadata.get("benchmark_metrics_calculated") is not False
    ):
        raise SequenceEvaluationError(f"evaluator input identity differs: {root}")
    if artifact_class == "chrombpnet_development_base_profiles":
        if (
            metadata.get("observed_atac_outcomes") is not True
            or metadata.get("model_inference_input_eligible") is not False
        ):
            raise SequenceEvaluationError("observed profile firewall differs")
    elif (
        metadata.get("model_id") != "training_pseudobulk_mean"
        or metadata.get("held_donor_outcomes_exposed") is not False
    ):
        raise SequenceEvaluationError("training baseline firewall differs")
    return metadata


def load_prediction_manifest(path: Path) -> list[dict[str, str]]:
    fields, rows = read_tsv(path)
    if fields == LEGACY_PREDICTION_MANIFEST_FIELDS:
        rows = [
            {
                "candidate_id": row["model_id"],
                "root_model_id": row["model_id"],
                "root": row["root"],
                "artifacts_sha256": row["artifacts_sha256"],
                "prediction_subdir": row["prediction_subdir"],
                "seed": "20260824",
            }
            for row in rows
        ]
    elif fields != PREDICTION_MANIFEST_FIELDS:
        raise SequenceEvaluationError("prediction manifest schema/count differs")
    if len(rows) < 1:
        raise SequenceEvaluationError("prediction manifest has no candidates")
    if len({row["candidate_id"] for row in rows}) != len(rows):
        raise SequenceEvaluationError("prediction candidate IDs are duplicated")
    for row in rows:
        if not row["candidate_id"] or not row["root_model_id"]:
            raise SequenceEvaluationError("prediction candidate identity is empty")
        try:
            seed = int(row["seed"])
        except ValueError as error:
            raise SequenceEvaluationError("prediction seed is invalid") from error
        if seed < 0:
            raise SequenceEvaluationError("prediction seed is invalid")
    return rows


def bind_prediction_model_id(
    summary: Mapping[str, Any], expected_model_id: str
) -> str:
    """Bind an optional subview ID through the already verified root metadata."""
    declared = summary.get("model_id")
    if declared is not None and declared != expected_model_id:
        raise SequenceEvaluationError("prediction summary model ID differs")
    return expected_model_id


def _verify_prediction_root(
    row: Mapping[str, str],
    *,
    dataset_id: str,
    split_id: str,
    lineage_id: str,
) -> tuple[Path, Path]:
    root = Path(row["root"]).resolve(strict=True)
    if sha256_file(root / "ARTIFACTS.json") != row["artifacts_sha256"]:
        raise SequenceEvaluationError(
            f"prediction ARTIFACTS differs: {row['candidate_id']}"
        )
    manifest = json.loads((root / "ARTIFACTS.json").read_text(encoding="utf-8"))
    metadata = manifest.get("metadata", {})
    if (
        metadata.get("model_id") != row["root_model_id"]
        or metadata.get("dataset_id") != dataset_id
        or metadata.get("split_id") != split_id
        or metadata.get("lineage_id") != lineage_id
        or metadata.get("seed") != int(row["seed"])
        or metadata.get("status") != "passed"
        or metadata.get("benchmark_metrics_calculated") is True
        or metadata.get("evaluator_outcomes_exposed") is True
        or metadata.get("held_donor_atac_exposed") is True
        or metadata.get("test_outcomes_used") is True
    ):
        raise SequenceEvaluationError(
            f"prediction outcome firewall differs: {row['candidate_id']}"
        )
    subdir = (root / row["prediction_subdir"]).resolve(strict=True)
    if root not in subdir.parents:
        raise SequenceEvaluationError("prediction subdirectory escapes root")
    return root, subdir


def load_prediction(
    row: Mapping[str, str],
    *,
    dataset_id: str,
    split_id: str,
    lineage_id: str,
) -> dict[str, Any]:
    import h5py
    import numpy as np

    _root, subdir = _verify_prediction_root(
        row,
        dataset_id=dataset_id,
        split_id=split_id,
        lineage_id=lineage_id,
    )
    summary = json.loads((subdir / "summary.json").read_text(encoding="utf-8"))
    bind_prediction_model_id(summary, row["root_model_id"])
    role_counts = summary.get("role_counts")
    if not isinstance(role_counts, dict) or set(role_counts) != {"test", "valid"}:
        raise SequenceEvaluationError(
            f"prediction role census differs: {row['candidate_id']}"
        )
    if any(not isinstance(value, int) or value <= 0 for value in role_counts.values()):
        raise SequenceEvaluationError(
            f"prediction role census differs: {row['candidate_id']}"
        )
    windows = sum(role_counts.values())
    output_length = summary.get("output_length")
    if (
        summary.get("status") != "pass"
        or summary.get("windows") != windows
        or not isinstance(output_length, int)
        or output_length <= 0
        or summary.get("seed") != int(row["seed"])
        or summary.get("observed_atac_input_exposed") is not False
        or summary.get("benchmark_metrics_calculated") is not False
    ):
        raise SequenceEvaluationError(
            f"prediction summary differs: {row['candidate_id']}"
        )
    counts_path = subdir / str(summary["regional_counts_path"])
    profile_path = subdir / str(summary["profile_probabilities_path"])
    fields, count_rows = read_tsv(counts_path)
    if (
        tuple(fields[:8]) != COUNT_REQUIRED_FIELDS
        or "strand_averaged_mass" not in fields
        or len(count_rows) != windows
    ):
        raise SequenceEvaluationError(
            f"prediction count table differs: {row['candidate_id']}"
        )
    counts: dict[tuple[str, str], dict[str, Any]] = {}
    for value in count_rows:
        key = (value["window_id"], value["selection_hash"])
        if key in counts:
            raise SequenceEvaluationError("prediction window is duplicated")
        mass = float(value["strand_averaged_mass"])
        if not math.isfinite(mass) or mass < 0:
            raise SequenceEvaluationError("predicted regional mass is invalid")
        counts[key] = {**value, "mass": mass}
    with h5py.File(profile_path, "r") as handle:
        ids = _decode(handle["window_id"][:])
        hashes = _decode(handle["selection_hash"][:])
        probabilities = np.asarray(handle["profile_probability"][:], dtype=np.float64)
    if probabilities.shape != (windows, output_length) or len(set(zip(ids, hashes))) != windows:
        raise SequenceEvaluationError("prediction profile tensor differs")
    if (
        np.any(~np.isfinite(probabilities))
        or np.any(probabilities < 0)
        or not np.allclose(probabilities.sum(axis=1), 1.0, rtol=0, atol=2.0e-6)
    ):
        raise SequenceEvaluationError("prediction profiles are not probabilities")
    profile_by_key = {
        key: index for index, key in enumerate(zip(ids, hashes, strict=True))
    }
    if set(profile_by_key) != set(counts):
        raise SequenceEvaluationError("prediction profile/count keys differ")
    return {
        "model_id": row["candidate_id"],
        "root_model_id": row["root_model_id"],
        "seed": int(row["seed"]),
        "role_counts": role_counts,
        "zero_mass_windows": {
            role: sum(
                value["mass"] == 0 and value["role"] == role
                for value in counts.values()
            )
            for role in ("valid", "test")
        },
        "counts": counts,
        "profile_by_key": profile_by_key,
        "profiles": probabilities,
    }


def evaluate(
    *,
    prediction_manifest: Path,
    observed_profiles: Path,
    observed_artifacts_sha256: str,
    training_baseline: Path,
    baseline_artifacts_sha256: str,
    dataset_id: str,
    split_id: str,
    lineage_id: str,
    evaluation_role: str,
    output: Path,
) -> dict[str, Any]:
    import h5py
    import numpy as np

    if output.exists():
        raise SequenceEvaluationError(f"output exists: {output}")
    role = validate_evaluation_role(evaluation_role)
    validate_evaluator_input_metadata(
        root=observed_profiles,
        expected_artifacts_sha256=observed_artifacts_sha256,
        artifact_class="chrombpnet_development_base_profiles",
        dataset_id=dataset_id,
        split_id=split_id,
        lineage_id=lineage_id,
    )
    validate_evaluator_input_metadata(
        root=training_baseline,
        expected_artifacts_sha256=baseline_artifacts_sha256,
        artifact_class="sequence_training_mean_baseline",
        dataset_id=dataset_id,
        split_id=split_id,
        lineage_id=lineage_id,
    )
    rows = load_prediction_manifest(prediction_manifest)
    predictions = [
        load_prediction(
            row,
            dataset_id=dataset_id,
            split_id=split_id,
            lineage_id=lineage_id,
        )
        for row in rows
    ]
    model_ids = ["uniform_global", "training_pseudobulk_mean"] + [
        item["model_id"] for item in predictions
    ]
    namespace = sha256(
        (observed_artifacts_sha256 + baseline_artifacts_sha256 + "\0" + "\0".join(model_ids)).encode()
    ).hexdigest()
    unit_rows: list[dict[str, object]] = []
    with h5py.File(
        observed_profiles / "profiles" / "profiles.h5", "r"
    ) as observed, h5py.File(
        training_baseline / "baseline" / "training_mean_baseline.h5", "r"
    ) as baseline:
        if bool(observed.attrs.get("model_inference_input_eligible")):
            raise SequenceEvaluationError("observed profiles are model-input eligible")
        truth_group = observed[role]
        truth_ids = _decode(truth_group["window_id"][:])
        truth_hashes = _decode(truth_group["selection_hash"][:])
        donors = _decode(truth_group["donor_id"][:])
        truth_counts = truth_group["counts"]
        base_group = baseline[role]
        base_ids = _decode(base_group["window_id"][:])
        base_hashes = _decode(base_group["selection_hash"][:])
        contigs = _decode(base_group["contig"][:])
        if (
            truth_ids != base_ids
            or truth_hashes != base_hashes
            or len(contigs) != len(truth_ids)
            or len(set(donors)) != len(donors)
            or truth_counts.shape[0] != len(donors)
            or truth_counts.shape[1] != len(truth_ids)
            or truth_counts.shape[2] != 1_000
        ):
            raise SequenceEvaluationError(f"{role} baseline/outcome axes differ")
        if baseline[role]["pseudobulk_counts"].shape != truth_counts.shape[1:]:
            raise SequenceEvaluationError(f"{role} baseline profile shape differs")
        keys = list(zip(truth_ids, truth_hashes, strict=True))
        baseline_profile = np.asarray(base_group["pseudobulk_counts"][:], dtype=np.float64)
        baseline_mass = baseline_profile.sum(axis=1)
        model_intensity: dict[str, np.ndarray] = {
            "uniform_global": np.ones(baseline_profile.shape, dtype=np.float64),
            "training_pseudobulk_mean": baseline_profile,
        }
        model_mass: dict[str, np.ndarray] = {
            "uniform_global": np.ones(len(keys), dtype=np.float64),
            "training_pseudobulk_mean": baseline_mass,
        }
        for prediction in predictions:
            if prediction["role_counts"][role] != len(keys):
                raise SequenceEvaluationError("prediction/outcome role census differs")
            if prediction["profiles"].shape[1] != truth_counts.shape[2]:
                raise SequenceEvaluationError("prediction/outcome profile width differs")
            offsets = [prediction["profile_by_key"].get(key, -1) for key in keys]
            if any(offset < 0 for offset in offsets):
                raise SequenceEvaluationError("prediction lacks an outcome window")
            rows_for_role = [prediction["counts"][key] for key in keys]
            if any(row["role"] != role for row in rows_for_role):
                raise SequenceEvaluationError("prediction role differs")
            mass = np.asarray([row["mass"] for row in rows_for_role], dtype=np.float64)
            probability = prediction["profiles"][np.asarray(offsets, dtype=np.int64)]
            model_intensity[prediction["model_id"]] = probability * mass[:, None]
            model_mass[prediction["model_id"]] = mass
        block_indices = {
            contig: np.asarray(
                [index for index, value in enumerate(contigs) if value == contig],
                dtype=np.int64,
            )
            for contig in sorted(set(contigs))
        }
        for donor_offset, donor_id in enumerate(donors):
            donor_truth = truth_counts[donor_offset]
            for contig, indexes in block_indices.items():
                observed_block = np.asarray(donor_truth[indexes, :], dtype=np.float64)
                observed_mass = observed_block.sum(axis=1)
                if observed_mass.sum() <= 0:
                    raise SequenceEvaluationError("held donor/block has zero outcome mass")
                for model_id in model_ids:
                    prediction_block = model_intensity[model_id][indexes]
                    unit_rows.append(
                        {
                            "role": role,
                            "model_id": model_id,
                            "donor_hash": _donor_hash(namespace, donor_id),
                            "block_hash": _block_hash(namespace, contig),
                            "observed_insertions": int(observed_mass.sum()),
                            "windows": len(indexes),
                            "block_deviance_per_insertion": format(
                                multinomial_deviance_per_insertion(
                                    observed_block, prediction_block
                                ),
                                ".17g",
                            ),
                            "within_window_deviance_per_insertion": format(
                                within_window_deviance_per_insertion(
                                    observed_block, prediction_block
                                ),
                                ".17g",
                            ),
                            "regional_count_spearman": format(
                                spearman_or_zero(
                                    observed_mass, model_mass[model_id][indexes]
                                ),
                                ".17g",
                            ),
                        }
                    )
    output.mkdir(mode=0o750)
    fields = (
        "role",
        "model_id",
        "donor_hash",
        "block_hash",
        "observed_insertions",
        "windows",
        "block_deviance_per_insertion",
        "within_window_deviance_per_insertion",
        "regional_count_spearman",
    )
    write_tsv(
        output / "donor_block_metrics.tsv",
        fields,
        sorted(
            unit_rows,
            key=lambda row: (
                str(row["role"]),
                str(row["model_id"]),
                str(row["donor_hash"]),
                str(row["block_hash"]),
            ),
        ),
    )
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in unit_rows:
        grouped[(str(row["role"]), str(row["model_id"]))].append(row)
    summaries: dict[str, dict[str, Any]] = defaultdict(dict)
    for (role, model_id), values in sorted(grouped.items()):
        summaries[role][model_id] = {
            "donor_block_units": len(values),
            "mean_block_deviance_per_insertion": sum(
                float(row["block_deviance_per_insertion"]) for row in values
            )
            / len(values),
            "mean_within_window_deviance_per_insertion": sum(
                float(row["within_window_deviance_per_insertion"]) for row in values
            )
            / len(values),
            "fisher_z_mean_regional_count_spearman": fisher_mean(
                [float(row["regional_count_spearman"]) for row in values]
            ),
        }
    strongest_baseline = min(
        ("uniform_global", "training_pseudobulk_mean"),
        key=lambda model_id: (
            summaries["valid"][model_id]["mean_block_deviance_per_insertion"],
            model_id,
        ),
    )
    valid_ranking = sorted(
        (item["model_id"] for item in predictions),
        key=lambda model_id: (
            summaries["valid"][model_id]["mean_block_deviance_per_insertion"],
            model_id,
        ),
    )
    relative: dict[str, dict[str, float]] = {}
    for model_id in valid_ranking:
        relative[model_id] = {}
        baseline_value = summaries[role][strongest_baseline][
            "mean_block_deviance_per_insertion"
        ]
        candidate = summaries[role][model_id]["mean_block_deviance_per_insertion"]
        relative[model_id][role] = (baseline_value - candidate) / baseline_value
    result = {
        "schema_version": "masld-bench-sequence-accessibility-evaluation-v3",
        "status": "pass",
        "task_scope": "sequence_native_regulatory_development_screen",
        "dataset_id": dataset_id,
        "split_id": split_id,
        "lineage_id": lineage_id,
        "evaluation_role": role,
        "biological_unit": "donor",
        "genomic_unit": "held_contig_block",
        "models": model_ids,
        "summaries": summaries,
        "strongest_training_only_baseline_selected_on_valid": strongest_baseline,
        "candidate_valid_ranking": valid_ranking,
        "relative_block_deviance_reduction": relative,
        "smoothing": {
            "method": "identical_uniform_probability_mixture",
            "epsilon": EPSILON_MIXTURE,
        },
        "donor_balanced": True,
        "cells_used_as_independent_replicates": False,
        "observed_profiles_exposed_to_models": False,
        "candidate_seeds": {
            item["model_id"]: item["seed"] for item in predictions
        },
        "prediction_diagnostics": {
            item["model_id"]: {
                "zero_mass_windows": item["zero_mass_windows"],
            }
            for item in predictions
        },
        "seeds_used_as_biological_replicates": False,
        "one_split_screen": True,
        "test_outcomes_read": False,
        "champion_claim_allowed": False,
        "promotion_gate_evaluated": False,
        "external_evaluation": False,
        "observed_artifacts_sha256": observed_artifacts_sha256,
        "baseline_artifacts_sha256": baseline_artifacts_sha256,
        "prediction_artifacts": {
            row["candidate_id"]: row["artifacts_sha256"] for row in rows
        },
        "donor_hash_namespace": namespace,
    }
    (output / "evaluation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prediction-manifest", type=Path, required=True)
    parser.add_argument("--observed-profiles", type=Path, required=True)
    parser.add_argument("--observed-artifacts-sha256", required=True)
    parser.add_argument("--training-baseline", type=Path, required=True)
    parser.add_argument("--baseline-artifacts-sha256", required=True)
    parser.add_argument("--dataset-id", default="gse296875")
    parser.add_argument("--split-id", default="donor0_genomic0")
    parser.add_argument("--lineage-id", default="hepatocyte")
    parser.add_argument("--evaluation-role", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate(
        prediction_manifest=args.prediction_manifest,
        observed_profiles=args.observed_profiles,
        observed_artifacts_sha256=args.observed_artifacts_sha256,
        training_baseline=args.training_baseline,
        baseline_artifacts_sha256=args.baseline_artifacts_sha256,
        dataset_id=args.dataset_id,
        split_id=args.split_id,
        lineage_id=args.lineage_id,
        evaluation_role=args.evaluation_role,
        output=args.output,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
