#!/usr/bin/env python3
"""Independently evaluate frozen donor-held MultiVI smoke predictions."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


CANDIDATE_MODEL_ID = "multivi"
PAIRING_TOPOLOGY = "same_nucleus"
MODEL_IDS = (CANDIDATE_MODEL_ID, "training_lineage_mean", "training_global_mean")
LINEAGES = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)
NAMESPACE = "gse296875:rna_conditioned_atac:donor_outer:smoke:v1"
PREDICTION_FIELDS = (
    "row_hash",
    "donor_hash",
    "block_hash",
    "stratum",
    "predicted",
)


class MultiVIEvaluationError(ValueError):
    """Raised when frozen MultiVI development evaluation differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def join_hash(namespace: str, kind: str, identifier: str) -> str:
    return sha256(f"{namespace}\0{kind}\0{identifier}".encode()).hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise MultiVIEvaluationError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
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


def decode(values: Any) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def read_csr(group: Any) -> Any:
    import numpy as np
    from scipy import sparse

    return sparse.csr_matrix(
        (
            np.asarray(group["data"][:]),
            np.asarray(group["indices"][:], dtype=np.int64),
            np.asarray(group["indptr"][:], dtype=np.int64),
        ),
        shape=tuple(int(value) for value in group["shape"][:]),
    )


def normalize_profile(values: Any, pseudocount: float = 1.0e-8) -> Any:
    import numpy as np

    profile = np.asarray(values, dtype=np.float64).ravel()
    if (
        profile.size == 0
        or np.any(~np.isfinite(profile))
        or np.any(profile < 0)
        or profile.sum() <= 0
    ):
        raise MultiVIEvaluationError("profile is invalid")
    profile += pseudocount
    profile /= profile.sum()
    return profile


def deviance_per_insertion(observed: Any, predicted: Any) -> float:
    import numpy as np

    truth = np.asarray(observed, dtype=np.float64)
    score = np.asarray(predicted, dtype=np.float64)
    if truth.shape != score.shape or truth.size == 0 or np.any(truth < 0):
        raise MultiVIEvaluationError("multinomial axes or truth differ")
    total = float(truth.sum())
    score_total = float(score.sum())
    if (
        total <= 0
        or score_total <= 0
        or np.any(~np.isfinite(score))
        or np.any(score <= 0)
    ):
        raise MultiVIEvaluationError("multinomial mass differs")
    probability = score / score_total
    positive = truth > 0
    expected = total * probability[positive]
    return float(
        2.0 * np.sum(truth[positive] * np.log(truth[positive] / expected)) / total
    )


def average_precision(observed: Any, scores: Any) -> float:
    from sklearn.metrics import average_precision_score

    labels = observed > 0
    if not labels.any():
        raise MultiVIEvaluationError("AUPRC has no positive region")
    return float(average_precision_score(labels, scores))


def spearman_or_zero(observed: Any, predicted: Any) -> float:
    import numpy as np
    from scipy.stats import spearmanr

    if np.all(observed == observed[0]) or np.all(predicted == predicted[0]):
        return 0.0
    result = spearmanr(observed, predicted)
    value = float(getattr(result, "statistic", result[0]))
    return value if math.isfinite(value) else 0.0


def load_bundle(root: Path, fold: int) -> tuple[dict[str, Any], dict[str, dict[str, str]]]:
    path = root / f"fold_{fold}/prediction_bundle.json"
    bundle = json.loads(path.read_text())
    if (
        bundle.get("schema_version") != "masld-bench-prediction-bundle-v1"
        or bundle.get("task_id") != "rna_conditioned_atac"
        or bundle.get("model_id") != CANDIDATE_MODEL_ID
        or bundle.get("dataset_ids") != ["gse296875"]
        or bundle.get("split_id") != "donor_outer"
        or bundle.get("unit_id_namespace") != NAMESPACE
        or bundle.get("biological_unit") != "donor"
        or bundle.get("metadata", {}).get("held_out_fold") != fold
        or bundle.get("metadata", {}).get("held_atac_input_exposed") is not False
        or bundle.get("metadata", {}).get("observed_atac_exported") is not False
        or bundle.get("metadata", {}).get("champion_claim_allowed") is not False
    ):
        raise MultiVIEvaluationError("MultiVI prediction bundle contract differs")
    metadata = bundle["metadata"]
    if CANDIDATE_MODEL_ID == "multivi" and (
        metadata.get("query_atac_state") != "structurally_missing"
        or metadata.get("query_atac_tensor_nonzero_values") != 0
    ):
        raise MultiVIEvaluationError("MultiVI query missingness contract differs")
    if CANDIDATE_MODEL_ID == "peakvi_training_context" and (
        metadata.get("held_atac_outcome_available_to_adapter") is not False
        or metadata.get("query_rna_used") is not False
        or metadata.get("rna_conditioned_atac_eligible") is not False
        or metadata.get("sealed_inference_eligible") is not False
    ):
        raise MultiVIEvaluationError("PeakVI context-transfer firewall differs")
    record = bundle["standardized_table"]
    table = path.parent / record["path"]
    if (
        table.parent != path.parent
        or table.is_symlink()
        or table.stat().st_size != record["size_bytes"]
        or sha256_file(table) != record["sha256"]
    ):
        raise MultiVIEvaluationError("prediction table changed")
    fields, rows = read_tsv(table)
    if fields != PREDICTION_FIELDS or len(rows) != bundle["n_predictions"]:
        raise MultiVIEvaluationError("prediction rows differ")
    index = {row["row_hash"]: row for row in rows}
    if len(index) != len(rows):
        raise MultiVIEvaluationError("prediction row hash is duplicated")
    return bundle, index


def evaluate(
    prediction_root: Path,
    prepared: Path,
    source_h5: Path,
    source_h5_sha256: str,
    output: Path,
) -> dict[str, Any]:
    import h5py
    import numpy as np

    if output.exists():
        raise MultiVIEvaluationError("evaluation output exists")
    if source_h5.is_symlink() or sha256_file(source_h5) != source_h5_sha256:
        raise MultiVIEvaluationError("source HDF5 differs")
    fold_data: dict[int, tuple[dict[str, Any], dict[str, dict[str, str]]]] = {
        fold: load_bundle(prediction_root, fold) for fold in range(5)
    }
    selected_rosters: list[list[dict[str, str]]] = []
    train_rosters: dict[int, list[dict[str, str]]] = {}
    query_rosters: dict[int, list[dict[str, str]]] = {}
    for fold in range(5):
        root = prepared / f"fold_{fold}"
        peak_fields, peaks = read_tsv(root / "selected_peaks.tsv")
        train_fields, train = read_tsv(root / "training_rows.tsv")
        query_fields, query = read_tsv(root / "query_rows.tsv")
        if (
            peak_fields
            != (
                "selected_index",
                "source_index",
                "peak_id",
                "chromosome",
                "bed_start",
                "bed_end",
            )
            or train_fields
            != ("cell_id", "donor_id", "lineage", "rna_state", "atac_state")
            or query_fields != train_fields
            or len(peaks) != 10_000
            or any(row["atac_state"] != "observed" for row in train)
            or any(row["atac_state"] != "structurally_missing" for row in query)
        ):
            raise MultiVIEvaluationError("prepared fold contract differs")
        selected_rosters.append(peaks)
        train_rosters[fold] = train
        query_rosters[fold] = query
    identity = [row["peak_id"] for row in selected_rosters[0]]
    if any([row["peak_id"] for row in roster] != identity for roster in selected_rosters[1:]):
        raise MultiVIEvaluationError("fold peak rosters differ")
    peaks = selected_rosters[0]
    selected = np.asarray([int(row["source_index"]) for row in peaks], dtype=np.int64)
    with h5py.File(source_h5, "r") as handle:
        cell_ids = decode(handle["obs/cell_id"][:])
        donors = decode(handle["obs/donor_id"][:])
        lineages = decode(handle["obs/broad_label"][:])
        atac = read_csr(handle["atac/counts_csr"])
        source_peak_ids = decode(handle["atac/peak_id"][:])
    if [source_peak_ids[index] for index in selected] != identity:
        raise MultiVIEvaluationError("selected peak source join differs")
    cell_index = {cell_id: index for index, cell_id in enumerate(cell_ids)}
    if len(cell_index) != len(cell_ids):
        raise MultiVIEvaluationError("source cell ID is duplicated")
    block_positions = {
        chromosome: np.asarray(
            [index for index, peak in enumerate(peaks) if peak["chromosome"] == chromosome],
            dtype=np.int64,
        )
        for chromosome in sorted({row["chromosome"] for row in peaks})
    }
    metric_rows: list[dict[str, Any]] = []
    secondary_rows: list[dict[str, Any]] = []
    grouped_deviance: dict[tuple[str, str], list[float]] = defaultdict(list)
    expected_prediction_rows = 0
    for fold in range(5):
        bundle, prediction_index = fold_data[fold]
        train_indices = np.asarray(
            [cell_index[row["cell_id"]] for row in train_rosters[fold]], dtype=np.int64
        )
        query_indices = np.asarray(
            [cell_index[row["cell_id"]] for row in query_rosters[fold]], dtype=np.int64
        )
        if (
            {donors[index] for index in train_indices}
            & {donors[index] for index in query_indices}
            or any(
                donors[index] != row["donor_id"] or lineages[index] != row["lineage"]
                for index, row in zip(query_indices, query_rosters[fold], strict=True)
            )
        ):
            raise MultiVIEvaluationError("source donor or lineage join differs")
        train_atac = atac[train_indices][:, selected]
        global_mean = normalize_profile(np.asarray(train_atac.sum(axis=0)).ravel())
        lineage_mean: dict[str, Any] = {}
        for lineage in LINEAGES:
            positions = [
                offset
                for offset, row in enumerate(train_rosters[fold])
                if row["lineage"] == lineage
            ]
            if not positions:
                raise MultiVIEvaluationError("training lineage is empty")
            lineage_mean[lineage] = normalize_profile(
                np.asarray(train_atac[positions].sum(axis=0)).ravel()
            )
        groups: dict[tuple[str, str], list[int]] = {}
        for index in query_indices:
            groups.setdefault((donors[index], lineages[index]), []).append(index)
        for (donor, lineage), positions in sorted(groups.items()):
            observed = np.asarray(atac[positions][:, selected].sum(axis=0)).ravel().astype(float)
            if observed.sum() <= 0:
                raise MultiVIEvaluationError("held donor-lineage ATAC is empty")
            candidate = np.empty(len(peaks), dtype=np.float64)
            for peak_index, peak in enumerate(peaks):
                row_hash = join_hash(
                    NAMESPACE, "row", f"{donor}\0{lineage}\0{peak['peak_id']}"
                )
                prediction = prediction_index.pop(row_hash, None)
                if (
                    prediction is None
                    or prediction["donor_hash"] != join_hash(NAMESPACE, "unit", donor)
                    or prediction["block_hash"]
                    != join_hash(NAMESPACE, "block", peak["chromosome"])
                    or prediction["stratum"] != lineage
                ):
                    raise MultiVIEvaluationError("prediction identity join differs")
                candidate[peak_index] = float(prediction["predicted"])
            expected_prediction_rows += len(peaks)
            vectors = {
                CANDIDATE_MODEL_ID: normalize_profile(candidate),
                "training_lineage_mean": lineage_mean[lineage],
                "training_global_mean": global_mean,
            }
            donor_hash = join_hash(NAMESPACE, "unit", donor)
            for model_id, vector in vectors.items():
                secondary_rows.append(
                    {
                        "model_id": model_id,
                        "donor_hash": donor_hash,
                        "stratum": lineage,
                        "peak_auprc": format(average_precision(observed, vector), ".17g"),
                        "profile_spearman": format(spearman_or_zero(observed, vector), ".17g"),
                    }
                )
                for chromosome, peak_positions in block_positions.items():
                    block_truth = observed[peak_positions]
                    if block_truth.sum() <= 0:
                        continue
                    value = deviance_per_insertion(block_truth, vector[peak_positions])
                    grouped_deviance[(model_id, lineage)].append(value)
                    metric_rows.append(
                        {
                            "model_id": model_id,
                            "donor_hash": donor_hash,
                            "block_hash": join_hash(NAMESPACE, "block", chromosome),
                            "stratum": lineage,
                            "observed_insertions": int(block_truth.sum()),
                            "regions": len(peak_positions),
                            "deviance_per_insertion": format(value, ".17g"),
                        }
                    )
        if prediction_index:
            raise MultiVIEvaluationError("prediction bundle contains unexpected rows")
        if bundle["n_predictions"] <= 0:
            raise MultiVIEvaluationError("prediction bundle is empty")
    model_mean = {
        model_id: sum(
            float(row["deviance_per_insertion"])
            for row in metric_rows
            if row["model_id"] == model_id
        )
        / sum(row["model_id"] == model_id for row in metric_rows)
        for model_id in MODEL_IDS
    }
    strongest = min(
        ("training_lineage_mean", "training_global_mean"),
        key=lambda model_id: (model_mean[model_id], model_id),
    )
    relative = (
        model_mean[strongest] - model_mean[CANDIDATE_MODEL_ID]
    ) / model_mean[strongest]
    lineage_gain: dict[str, float] = {}
    for lineage in LINEAGES:
        baseline = grouped_deviance[(strongest, lineage)]
        candidate = grouped_deviance[(CANDIDATE_MODEL_ID, lineage)]
        if not baseline or len(baseline) != len(candidate):
            raise MultiVIEvaluationError("lineage deviance universe differs")
        baseline_mean = sum(baseline) / len(baseline)
        candidate_mean = sum(candidate) / len(candidate)
        lineage_gain[lineage] = (baseline_mean - candidate_mean) / baseline_mean
    output.mkdir(parents=True, mode=0o750)
    write_tsv(
        output / "donor_lineage_block_metrics.tsv",
        (
            "model_id",
            "donor_hash",
            "block_hash",
            "stratum",
            "observed_insertions",
            "regions",
            "deviance_per_insertion",
        ),
        sorted(
            metric_rows,
            key=lambda row: (
                row["model_id"], row["donor_hash"], row["stratum"], row["block_hash"]
            ),
        ),
    )
    write_tsv(
        output / "donor_lineage_secondary_metrics.tsv",
        ("model_id", "donor_hash", "stratum", "peak_auprc", "profile_spearman"),
        sorted(
            secondary_rows,
            key=lambda row: (row["model_id"], row["donor_hash"], row["stratum"]),
        ),
    )
    result = {
        "schema_version": f"masld-bench-{CANDIDATE_MODEL_ID}-smoke-evaluation-v1",
        "status": "pass",
        "model_id": CANDIDATE_MODEL_ID,
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "evaluation_role": "development_smoke",
        "pairing_topology": PAIRING_TOPOLOGY,
        "biological_unit": "donor",
        "n_donors": len(set(donors)),
        "n_lineages": len(LINEAGES),
        "n_regions": len(peaks),
        "n_genomic_blocks": len(block_positions),
        "prediction_rows": expected_prediction_rows,
        "mean_deviance_per_insertion": model_mean,
        "strongest_training_only_baseline": strongest,
        "candidate_relative_deviance_reduction": relative,
        "lineage_relative_deviance_reduction": lineage_gain,
        "development_gate": {
            "minimum_relative_deviance_reduction": 0.05,
            "minimum_improved_lineages": 4,
            "maximum_allowed_lineage_worsening": -0.02,
            "overall_threshold_passed": relative >= 0.05,
            "improved_lineages": sum(value > 0 for value in lineage_gain.values()),
            "no_lineage_worse_than_threshold": min(lineage_gain.values()) >= -0.02,
        },
        "prediction_frozen_before_outcomes": True,
        "held_atac_exposed_to_model": False,
        "held_atac_read_by_evaluator_only": True,
        "test_outcomes_read": False,
        "external_evaluation": False,
        "smoke_only": True,
        "one_seed_five_donor_folds": True,
        "promotion_gate_finalized": False,
        "champion_claim_allowed": False,
    }
    (output / "evaluation.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction-root", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--source-h5", type=Path, required=True)
    parser.add_argument("--source-h5-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    evaluate(
        arguments.prediction_root,
        arguments.prepared,
        arguments.source_h5,
        arguments.source_h5_sha256,
        arguments.output,
    )


if __name__ == "__main__":
    main()
