#!/usr/bin/env python
"""Independent development evaluator for RNA-conditioned ATAC smoke runs."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


class RNAATACEvaluationError(RuntimeError):
    """Raised when an OOF prediction/outcome join does not meet its requirements."""


TASK_ID = "rna_conditioned_atac"
DATASET_ID = "gse296875"
VIEW_ID = "gse296875_rna_atac_smoke_1000_v1"
MODEL_IDS = ("assay_native_pseudobulk", "mean_track", "nearest_context")
LINEAGES = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)
PREDICTION_FIELDS = (
    "row_hash",
    "donor_hash",
    "block_hash",
    "stratum",
    "predicted",
)


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _canonical_hash(value: Any) -> str:
    return sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        handle.write(_canonical_json(value))
        handle.write("\n")


def _write_tsv(
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
        for row in rows:
            writer.writerow(dict(row))


def _read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RNAATACEvaluationError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def _decode(values: Any) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def fold_index(unit_id: str, *, seed: int, outer_folds: int) -> int:
    digest = sha256(f"{seed}\0{unit_id}".encode()).digest()
    return int.from_bytes(digest[:8], "big") % outer_folds


def join_hash(namespace: str, kind: str, identifier: str) -> str:
    return sha256(f"{namespace}\0{kind}\0{identifier}".encode()).hexdigest()


def deterministic_peak_indices(peak_ids: Sequence[str], n_peaks: int) -> tuple[int, ...]:
    ranked = sorted(
        range(len(peak_ids)),
        key=lambda index: (
            sha256(f"masld-rna-atac-smoke-v1\0{peak_ids[index]}".encode()).digest(),
            peak_ids[index],
            index,
        ),
    )
    return tuple(sorted(ranked[:n_peaks]))


def _read_csr(group: Any) -> Any:
    import numpy as np
    from scipy import sparse

    return sparse.csr_matrix(
        (
            np.asarray(group["data"][:]),
            np.asarray(group["indices"][:], dtype=np.int64),
            np.asarray(group["indptr"][:], dtype=np.int64),
        ),
        shape=tuple(int(item) for item in group["shape"][:]),
    )


def _average_precision(observed: Any, scores: Any) -> float:
    import numpy as np

    labels = np.asarray(observed, dtype=np.int8)
    values = np.asarray(scores, dtype=np.float64)
    positives = int(labels.sum())
    if positives == 0:
        raise RNAATACEvaluationError("peak AUPRC is undefined without positives")
    order = np.argsort(-values, kind="stable")
    labels = labels[order]
    values = values[order]
    true_positive = 0
    false_positive = 0
    previous_recall = 0.0
    area = 0.0
    index = 0
    while index < len(labels):
        stop = index + 1
        while stop < len(labels) and values[stop] == values[index]:
            stop += 1
        group = labels[index:stop]
        true_positive += int(group.sum())
        false_positive += len(group) - int(group.sum())
        recall = true_positive / positives
        precision = true_positive / (true_positive + false_positive)
        area += (recall - previous_recall) * precision
        previous_recall = recall
        index = stop
    return area


def _multinomial_deviance(observed: Any, predicted: Any) -> float:
    import numpy as np

    truth = np.asarray(observed, dtype=np.float64)
    score = np.asarray(predicted, dtype=np.float64)
    if truth.size == 0 or truth.shape != score.shape or np.any(truth < 0) or np.any(score < 0):
        raise RNAATACEvaluationError("invalid multinomial profile")
    truth_total = float(truth.sum())
    score_total = float(score.sum())
    if truth_total <= 0 or score_total <= 0:
        raise RNAATACEvaluationError("multinomial profiles must have positive totals")
    expected = truth_total * score / score_total
    positive = truth > 0
    if np.any(expected[positive] <= 0):
        return math.inf
    return float(2.0 * np.sum(truth[positive] * np.log(truth[positive] / expected[positive])))


def _load_bundle(path: Path, *, model_id: str, fold: int) -> list[dict[str, str]]:
    try:
        bundle = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RNAATACEvaluationError(f"invalid prediction bundle: {path}") from error
    exact = {
        "schema_version": "masld-bench-prediction-bundle-v1",
        "task_id": TASK_ID,
        "model_id": model_id,
        "dataset_ids": [DATASET_ID],
        "split_id": "donor_outer",
        "row_id_field": "row_hash",
        "unit_id_field": "donor_hash",
        "biological_unit": "donor",
        "format_version": "tsv-v1",
        "missing_state": "observed",
    }
    for field, expected in exact.items():
        if bundle.get(field) != expected:
            raise RNAATACEvaluationError(f"prediction bundle {field} differs")
    metadata = bundle.get("metadata")
    if (
        not isinstance(metadata, Mapping)
        or metadata.get("held_out_fold") != fold
        or metadata.get("held_atac_input_exposed") is not False
        or metadata.get("observed_atac_exported") is not False
    ):
        raise RNAATACEvaluationError("prediction bundle held-ATAC firewall differs")
    record = bundle.get("standardized_table")
    if not isinstance(record, Mapping) or record.get("role") != f"standardized_prediction_table:{TASK_ID}":
        raise RNAATACEvaluationError("standardized prediction artifact differs")
    table = path.parent / str(record.get("path", ""))
    if (
        table.parent != path.parent
        or table.is_symlink()
        or not table.is_file()
        or table.stat().st_size != record.get("size_bytes")
        or _sha256_file(table) != record.get("sha256")
    ):
        raise RNAATACEvaluationError("standardized prediction artifact changed")
    fields, rows = _read_tsv(table)
    if fields != PREDICTION_FIELDS or len(rows) != bundle.get("n_predictions"):
        raise RNAATACEvaluationError("prediction table schema or count differs")
    if bundle.get("table_schema_sha256") != _canonical_hash(
        {"format": "tsv", "fields": list(PREDICTION_FIELDS)}
    ):
        raise RNAATACEvaluationError("prediction table schema hash differs")
    return rows


def _discover_runs(candidate: Path, execution_root: Path) -> dict[str, dict[int, Path]]:
    try:
        plan = json.loads((candidate / "plan.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RNAATACEvaluationError("campaign plan is invalid") from error
    result = {model_id: {} for model_id in MODEL_IDS}
    runs = plan.get("runs")
    if not isinstance(runs, list):
        raise RNAATACEvaluationError("campaign plan has no runs")
    for run in runs:
        if not isinstance(run, Mapping) or run.get("task_id") != TASK_ID:
            continue
        model_id = str(run.get("model_id", ""))
        if model_id not in result:
            continue
        fold = run.get("fold")
        if isinstance(fold, bool) or not isinstance(fold, int) or fold in result[model_id]:
            raise RNAATACEvaluationError("campaign has invalid duplicate fold")
        bundle = (
            execution_root
            / "runs"
            / str(run["run_id"])
            / "attempt-001"
            / "adapter_actions"
            / "003-predict"
            / "prediction_bundle.json"
        )
        if not bundle.is_file():
            raise RNAATACEvaluationError(f"prediction bundle is missing: {bundle}")
        result[model_id][fold] = bundle
    if any(set(folds) != set(range(5)) for folds in result.values()):
        raise RNAATACEvaluationError("each baseline requires all five OOF folds")
    return result


def evaluate(
    *,
    candidate: Path,
    execution_root: Path,
    h5_path: Path,
    h5_sha256: str,
    output: Path,
) -> dict[str, Any]:
    import h5py
    import numpy as np
    from scipy.stats import spearmanr

    if output.exists():
        raise RNAATACEvaluationError(f"evaluation output already exists: {output}")
    if h5_path.is_symlink() or not h5_path.is_file() or _sha256_file(h5_path) != h5_sha256:
        raise RNAATACEvaluationError("GSE296875 smoke HDF5 changed")
    run_paths = _discover_runs(candidate, execution_root)
    all_predictions: dict[str, dict[str, dict[str, str]]] = {}
    for model_id, folds in run_paths.items():
        index: dict[str, dict[str, str]] = {}
        for fold, bundle_path in sorted(folds.items()):
            for row in _load_bundle(bundle_path, model_id=model_id, fold=fold):
                if row["row_hash"] in index:
                    raise RNAATACEvaluationError("OOF prediction row appears twice")
                index[row["row_hash"]] = row
        all_predictions[model_id] = index
    with h5py.File(h5_path, "r") as handle:
        donors = _decode(handle["obs/donor_id"][:])
        lineages = _decode(handle["obs/broad_label"][:])
        peak_ids = _decode(handle["atac/peak_id"][:])
        chromosomes = _decode(handle["atac/chromosome"][:])
        atac = _read_csr(handle["atac/counts_csr"])
    first_bundle = next(iter(next(iter(run_paths.values())).values()))
    first_payload = json.loads(first_bundle.read_text(encoding="utf-8"))
    metadata = first_payload["metadata"]
    n_peaks = int(metadata["selected_peak_count"])
    namespace = str(first_payload["unit_id_namespace"])
    selected = deterministic_peak_indices(peak_ids, n_peaks)
    groups: dict[tuple[str, str], list[int]] = {}
    for index, (donor, lineage) in enumerate(zip(donors, lineages, strict=True)):
        groups.setdefault((donor, lineage), []).append(index)
    expected_row_ids: set[str] = set()
    profile_metrics: list[dict[str, Any]] = []
    outcome_by_row: dict[str, dict[str, Any]] = {}
    deviance_by_model: dict[str, list[float]] = {model_id: [] for model_id in MODEL_IDS}
    auprc_by_model: dict[str, list[float]] = {model_id: [] for model_id in MODEL_IDS}
    spearman_by_model: dict[str, list[float]] = {model_id: [] for model_id in MODEL_IDS}
    lineage_deviance: dict[str, dict[str, list[float]]] = {
        model_id: {lineage: [] for lineage in LINEAGES} for model_id in MODEL_IDS
    }
    selected_chromosomes = [chromosomes[index] for index in selected]
    block_indices = {
        chromosome: np.asarray(
            [index for index, value in enumerate(selected_chromosomes) if value == chromosome],
            dtype=np.int64,
        )
        for chromosome in sorted(set(selected_chromosomes))
    }
    if len(block_indices) < 2:
        raise RNAATACEvaluationError("smoke peak universe has fewer than two genomic blocks")
    for (donor, lineage), cell_indices in sorted(groups.items()):
        observed = np.asarray(atac[cell_indices][:, selected].sum(axis=0)).ravel().astype(float)
        if observed.sum() <= 0:
            raise RNAATACEvaluationError("observed donor-lineage profile is empty")
        model_vectors: dict[str, np.ndarray] = {}
        for model_id in MODEL_IDS:
            values = np.empty(n_peaks, dtype=np.float64)
            for peak_position, source_index in enumerate(selected):
                row_hash = join_hash(
                    namespace,
                    "row",
                    f"{donor}\0{lineage}\0{peak_ids[source_index]}",
                )
                expected_row_ids.add(row_hash)
                prediction = all_predictions[model_id].get(row_hash)
                expected_fold = fold_index(donor, seed=20260821, outer_folds=5)
                if prediction is None:
                    raise RNAATACEvaluationError("expected OOF prediction is missing")
                if (
                    prediction["donor_hash"] != join_hash(namespace, "unit", donor)
                    or prediction["block_hash"] != join_hash(namespace, "block", chromosomes[source_index])
                    or prediction["stratum"] != lineage
                ):
                    raise RNAATACEvaluationError("prediction join identity differs")
                try:
                    value = float(prediction["predicted"])
                except ValueError as error:
                    raise RNAATACEvaluationError("prediction is not numeric") from error
                if not math.isfinite(value) or value <= 0:
                    raise RNAATACEvaluationError("profile prediction must be finite and positive")
                values[peak_position] = value
                outcome_by_row.setdefault(
                    row_hash,
                    {
                        "row_hash": row_hash,
                        "donor_hash": prediction["donor_hash"],
                        "block_hash": prediction["block_hash"],
                        "stratum": lineage,
                        "observed": format(float(observed[peak_position]), ".17g"),
                        "fold": expected_fold,
                    },
                )
            model_vectors[model_id] = values
            auprc = _average_precision(observed > 0, values)
            correlation = float(spearmanr(observed, values).statistic)
            if not math.isfinite(correlation):
                correlation = 0.0
            auprc_by_model[model_id].append(auprc)
            spearman_by_model[model_id].append(correlation)
            for chromosome, positions in block_indices.items():
                block_observed = observed[positions]
                if block_observed.sum() <= 0:
                    continue
                deviance = _multinomial_deviance(block_observed, values[positions])
                if not math.isfinite(deviance):
                    raise RNAATACEvaluationError("profile deviance is not finite")
                deviance_by_model[model_id].append(deviance)
                lineage_deviance[model_id][lineage].append(deviance)
                profile_metrics.append(
                    {
                        "model_id": model_id,
                        "donor_hash": join_hash(namespace, "unit", donor),
                        "block_hash": join_hash(namespace, "block", chromosome),
                        "stratum": lineage,
                        "deviance": format(deviance, ".17g"),
                        "peak_auprc": format(auprc, ".17g"),
                        "profile_spearman": format(correlation, ".17g"),
                    }
                )
    for model_id, rows in all_predictions.items():
        if set(rows) != expected_row_ids:
            raise RNAATACEvaluationError(f"{model_id} OOF row set differs from outcomes")
    baseline_candidates = ("mean_track", "assay_native_pseudobulk")
    strongest_baseline = min(
        baseline_candidates,
        key=lambda model_id: (sum(deviance_by_model[model_id]), model_id),
    )
    candidate_model = "nearest_context"
    endpoint_rows = []
    for row_hash in sorted(outcome_by_row):
        outcome = outcome_by_row[row_hash]
        endpoint_rows.append(
            {
                "row_hash": row_hash,
                "donor_hash": outcome["donor_hash"],
                "block_hash": outcome["block_hash"],
                "stratum": outcome["stratum"],
                "observed": outcome["observed"],
                "candidate": all_predictions[candidate_model][row_hash]["predicted"],
                "baseline": all_predictions[strongest_baseline][row_hash]["predicted"],
            }
        )
    output.mkdir(parents=True, exist_ok=False)
    profile_path = output / "profile_metrics.tsv"
    endpoint_path = output / "nearest_context_vs_strongest_baseline_endpoint.tsv"
    _write_tsv(
        profile_path,
        ("model_id", "donor_hash", "block_hash", "stratum", "deviance", "peak_auprc", "profile_spearman"),
        sorted(profile_metrics, key=lambda row: (row["model_id"], row["donor_hash"], row["block_hash"], row["stratum"])),
    )
    _write_tsv(
        endpoint_path,
        ("row_hash", "donor_hash", "block_hash", "stratum", "observed", "candidate", "baseline"),
        endpoint_rows,
    )
    summary_models = {}
    for model_id in MODEL_IDS:
        summary_models[model_id] = {
            "total_multinomial_deviance": sum(deviance_by_model[model_id]),
            "mean_peak_auprc": sum(auprc_by_model[model_id]) / len(auprc_by_model[model_id]),
            "mean_profile_spearman": sum(spearman_by_model[model_id]) / len(spearman_by_model[model_id]),
            "lineage_mean_deviance": {
                lineage: sum(values) / len(values)
                for lineage, values in lineage_deviance[model_id].items()
                if values
            },
        }
    baseline_deviance = sum(deviance_by_model[strongest_baseline])
    candidate_deviance = sum(deviance_by_model[candidate_model])
    result = {
        "schema_version": "masld-bench-rna-atac-development-evaluation-v1",
        "task_id": TASK_ID,
        "dataset_id": DATASET_ID,
        "dataset_view_id": VIEW_ID,
        "input_h5_sha256": h5_sha256,
        "unit_of_inference": "donor",
        "genomic_resampling_unit": "chromosome_block_smoke_only",
        "n_donors": len(set(donors)),
        "n_donor_lineage_profiles": len(groups),
        "n_selected_peaks": n_peaks,
        "n_genomic_blocks": len(block_indices),
        "n_oof_prediction_rows_per_model": len(expected_row_ids),
        "models": summary_models,
        "strongest_simple_baseline": strongest_baseline,
        "nearest_context_relative_deviance_reduction": (
            baseline_deviance - candidate_deviance
        ) / baseline_deviance,
        "smoke_only": True,
        "champion_claim_allowed": False,
        "cells_used_as_independent_replicates": False,
        "held_atac_exposed_to_fit_or_predict": False,
        "outcomes_joined_only_by_independent_evaluator": True,
        "artifacts": {
            "profile_metrics": {
                "path": profile_path.name,
                "sha256": _sha256_file(profile_path),
                "size_bytes": profile_path.stat().st_size,
            },
            "endpoint_table": {
                "path": endpoint_path.name,
                "sha256": _sha256_file(endpoint_path),
                "size_bytes": endpoint_path.stat().st_size,
            },
        },
    }
    _write_json(output / "evaluation.json", result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--execution-root", required=True, type=Path)
    parser.add_argument("--h5", required=True, type=Path)
    parser.add_argument("--h5-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args(argv)
    result = evaluate(
        candidate=arguments.candidate,
        execution_root=arguments.execution_root,
        h5_path=arguments.h5,
        h5_sha256=arguments.h5_sha256,
        output=arguments.output,
    )
    print(_canonical_json(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
