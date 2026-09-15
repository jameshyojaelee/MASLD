#!/usr/bin/env python
"""Five-seed evaluator for the frozen v5 Cobolt RNA-to-ATAC smoke campaign.

This evaluator treats seeds as computational stability checks, never as
biological replicates.  It verifies every read-only prediction attempt before
opening the outcome HDF5, then scores each seed and the arithmetic five-seed
ensemble on donor by chromosome-block units.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from masld_bench.campaign import verify_run_execution_attempt

from . import rna_atac_cobolt_development as cobolt_v1
from . import rna_atac_development as base
from . import rna_atac_scpair_development as scpair_v1


CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v5"
MODEL_IDS = (
    "assay_native_pseudobulk",
    "cobolt",
    "mean_track",
    "nearest_context",
    "shrunken_pseudobulk",
    "shuffled_context",
    "trans_only",
)
TASK_NATIVE_BASELINES = (
    "assay_native_pseudobulk",
    "mean_track",
    "shrunken_pseudobulk",
)
SEEDS = (1103, 2107, 3109, 4111, 5113)
FOLDS = (0, 1, 2, 3, 4)
EXPECTED_RUNS = len(MODEL_IDS) * len(SEEDS) * len(FOLDS)
ENSEMBLE_LABEL = "five_seed_arithmetic_mean"


class RNAATACCoboltFiveSeedEvaluationError(RuntimeError):
    """Raised when the v5 grid or outcome-separated evaluation differs."""


@dataclass(frozen=True)
class RunCell:
    model_id: str
    seed: int
    fold: int
    run_id: str
    attempt: Path
    prediction_bundle: Path


Verifier = Callable[..., Mapping[str, Any]]


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RNAATACCoboltFiveSeedEvaluationError(f"invalid {label}: {path}") from error
    if not isinstance(value, dict):
        raise RNAATACCoboltFiveSeedEvaluationError(f"{label} is not an object")
    return value


def relative_deviance_reduction(reference: float, candidate: float) -> float:
    try:
        return cobolt_v1.relative_deviance_reduction(reference, candidate)
    except cobolt_v1.RNAATACCoboltEvaluationError as error:
        raise RNAATACCoboltFiveSeedEvaluationError(str(error)) from error


def discover_verified_run_grid(
    *,
    candidate: Path,
    execution_root: Path,
    candidate_sha256: str,
    verifier: Verifier = verify_run_execution_attempt,
) -> dict[str, dict[int, dict[int, RunCell]]]:
    """Verify the exact 7 x 5 x 5 grid before exposing prediction paths."""

    candidate = candidate.resolve(strict=True)
    execution_root = execution_root.resolve(strict=True)
    manifest = candidate / "ARTIFACTS.json"
    if base._sha256_file(manifest) != candidate_sha256:
        raise RNAATACCoboltFiveSeedEvaluationError(
            "candidate artifact hash differs from the frozen v5 campaign"
        )
    plan = _read_json(candidate / "plan.json", "campaign plan")
    campaign = plan.get("campaign")
    if not isinstance(campaign, Mapping) or campaign.get("campaign_id") != CAMPAIGN_ID:
        raise RNAATACCoboltFiveSeedEvaluationError("campaign identity is not v5")
    runs = plan.get("runs")
    if not isinstance(runs, list) or len(runs) != EXPECTED_RUNS:
        raise RNAATACCoboltFiveSeedEvaluationError(
            f"v5 requires exactly {EXPECTED_RUNS} planned runs"
        )

    planned: dict[tuple[str, int, int], tuple[str, Path]] = {}
    for run in runs:
        if not isinstance(run, Mapping):
            raise RNAATACCoboltFiveSeedEvaluationError("campaign run is not an object")
        model_id = run.get("model_id")
        seed = run.get("seed")
        fold = run.get("fold")
        run_id = run.get("run_id")
        if (
            run.get("task_id") != base.TASK_ID
            or model_id not in MODEL_IDS
            or isinstance(seed, bool)
            or seed not in SEEDS
            or isinstance(fold, bool)
            or fold not in FOLDS
            or not isinstance(run_id, str)
            or len(run_id) != 64
            or run.get("action") != ["prepare", "fit", "predict"]
        ):
            raise RNAATACCoboltFiveSeedEvaluationError(
                "campaign run differs from the v5 model/seed/fold contract"
            )
        key = (str(model_id), int(seed), int(fold))
        if key in planned:
            raise RNAATACCoboltFiveSeedEvaluationError(
                f"duplicate v5 model/seed/fold cell: {key}"
            )
        attempt = execution_root / "runs" / run_id / "attempt-001"
        planned[key] = (run_id, attempt)

    expected = {
        (model_id, seed, fold)
        for model_id in MODEL_IDS
        for seed in SEEDS
        for fold in FOLDS
    }
    if set(planned) != expected:
        missing = sorted(expected - set(planned))
        extra = sorted(set(planned) - expected)
        raise RNAATACCoboltFiveSeedEvaluationError(
            f"v5 grid is incomplete or expanded; missing={missing[:3]}, extra={extra[:3]}"
        )

    # No outcome record is parsed or joined before all 175 prediction attempts
    # pass. Attempt verification may checksum the held-back input file as an
    # opaque byte stream, but it never interprets outcome fields.
    for key in sorted(planned):
        _, attempt = planned[key]
        try:
            verifier(attempt, require_succeeded=True)
        except Exception as error:
            raise RNAATACCoboltFiveSeedEvaluationError(
                f"immutable prediction attempt failed verification for {key}: {attempt}"
            ) from error

    grid: dict[str, dict[int, dict[int, RunCell]]] = {
        model_id: {seed: {} for seed in SEEDS} for model_id in MODEL_IDS
    }
    for (model_id, seed, fold), (run_id, attempt) in sorted(planned.items()):
        bundle = attempt / "adapter_actions" / "003-predict" / "prediction_bundle.json"
        if bundle.is_symlink() or not bundle.is_file():
            raise RNAATACCoboltFiveSeedEvaluationError(
                f"verified attempt has no regular prediction bundle: {bundle}"
            )
        grid[model_id][seed][fold] = RunCell(
            model_id=model_id,
            seed=seed,
            fold=fold,
            run_id=run_id,
            attempt=attempt,
            prediction_bundle=bundle,
        )
    return grid


def _prediction_inventory(
    grid: Mapping[str, Mapping[int, Mapping[int, RunCell]]],
) -> tuple[list[dict[str, Any]], int, str, dict[str, list[dict[str, Any]]]]:
    records: list[dict[str, Any]] = []
    selected_peak_counts: set[int] = set()
    namespaces: set[str] = set()
    cobolt_fit_receipts: dict[str, list[dict[str, Any]]] = {}
    regularized_fit_selections: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for model_id in MODEL_IDS:
        for seed in SEEDS:
            for fold in FOLDS:
                cell = grid[model_id][seed][fold]
                bundle = _read_json(cell.prediction_bundle, "prediction bundle")
                metadata = bundle.get("metadata")
                table = bundle.get("standardized_table")
                namespace = bundle.get("unit_id_namespace")
                if (
                    not isinstance(metadata, Mapping)
                    or metadata.get("held_out_fold") != fold
                    or metadata.get("held_atac_input_exposed") is not False
                    or metadata.get("observed_atac_exported") is not False
                    or isinstance(metadata.get("selected_peak_count"), bool)
                    or not isinstance(metadata.get("selected_peak_count"), int)
                    or not isinstance(namespace, str)
                    or not namespace
                    or not isinstance(table, Mapping)
                    or table.get("role")
                    != f"standardized_prediction_table:{base.TASK_ID}"
                ):
                    raise RNAATACCoboltFiveSeedEvaluationError(
                        f"prediction metadata differs for {(model_id, seed, fold)}"
                    )
                table_path = cell.prediction_bundle.parent / str(table.get("path", ""))
                if (
                    table_path.parent != cell.prediction_bundle.parent
                    or table_path.is_symlink()
                    or not table_path.is_file()
                    or table_path.stat().st_size != table.get("size_bytes")
                    or base._sha256_file(table_path) != table.get("sha256")
                ):
                    raise RNAATACCoboltFiveSeedEvaluationError(
                        "standardized prediction table differs after attempt verification"
                    )
                selected_peak_counts.add(int(metadata["selected_peak_count"]))
                namespaces.add(namespace)
                records.append(
                    {
                        "model_id": model_id,
                        "seed": seed,
                        "fold": fold,
                        "run_id": cell.run_id,
                        "attempt_path": cell.attempt.as_posix(),
                        "prediction_bundle_path": cell.prediction_bundle.as_posix(),
                        "prediction_bundle_sha256": base._sha256_file(
                            cell.prediction_bundle
                        ),
                        "prediction_table_path": table_path.as_posix(),
                        "prediction_table_sha256": str(table["sha256"]),
                        "n_predictions": bundle.get("n_predictions"),
                    }
                )
            seed_paths = {
                model_id: {
                    fold: grid[model_id][seed][fold].prediction_bundle
                    for fold in FOLDS
                }
            }
            if model_id == cobolt_v1.CANDIDATE_MODEL:
                cobolt_fit_receipts[str(seed)] = cobolt_v1._cobolt_fit_receipts(
                    seed_paths
                )
        if model_id == MODEL_IDS[-1]:
            # The selection helper needs all three regularized baseline paths.
            for seed in SEEDS:
                seed_paths = {
                    baseline: {
                        fold: grid[baseline][seed][fold].prediction_bundle
                        for fold in FOLDS
                    }
                    for baseline in (
                        "shrunken_pseudobulk",
                        "shuffled_context",
                        "trans_only",
                    )
                }
                regularized_fit_selections[str(seed)] = (
                    scpair_v1._regularized_fit_selections(seed_paths)
                )
    if len(selected_peak_counts) != 1 or len(namespaces) != 1:
        raise RNAATACCoboltFiveSeedEvaluationError(
            "prediction bundles disagree on peak count or row namespace"
        )
    receipts = {
        "cobolt_fit_receipts_by_seed_and_outer_fold": cobolt_fit_receipts,
        "regularized_fit_selections_by_seed_and_outer_fold": regularized_fit_selections,
    }
    return records, next(iter(selected_peak_counts)), next(iter(namespaces)), receipts


def _load_outcomes(
    *, h5_path: Path, h5_sha256: str, selected_peak_count: int, namespace: str
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    import h5py
    import numpy as np

    if (
        h5_path.is_symlink()
        or not h5_path.is_file()
        or base._sha256_file(h5_path) != h5_sha256
    ):
        raise RNAATACCoboltFiveSeedEvaluationError("GSE296875 smoke HDF5 changed")
    with h5py.File(h5_path, "r") as handle:
        donors = base._decode(handle["obs/donor_id"][:])
        lineages = base._decode(handle["obs/broad_label"][:])
        peak_ids = base._decode(handle["atac/peak_id"][:])
        chromosomes = base._decode(handle["atac/chromosome"][:])
        atac = base._read_csr(handle["atac/counts_csr"])
    selected = base.deterministic_peak_indices(peak_ids, selected_peak_count)
    groups: dict[tuple[str, str], list[int]] = {}
    for index, (donor, lineage) in enumerate(zip(donors, lineages, strict=True)):
        groups.setdefault((donor, lineage), []).append(index)
    outcomes: list[dict[str, Any]] = []
    for (donor, lineage), cell_indices in sorted(groups.items()):
        observed = np.asarray(atac[cell_indices][:, selected].sum(axis=0)).ravel()
        if float(observed.sum()) <= 0:
            raise RNAATACCoboltFiveSeedEvaluationError(
                "observed donor-lineage profile is empty"
            )
        donor_hash = base.join_hash(namespace, "unit", donor)
        fold = base.fold_index(donor, seed=20260821, outer_folds=5)
        for peak_position, source_index in enumerate(selected):
            outcomes.append(
                {
                    "row_hash": base.join_hash(
                        namespace,
                        "row",
                        f"{donor}\0{lineage}\0{peak_ids[source_index]}",
                    ),
                    "donor_hash": donor_hash,
                    "block_hash": base.join_hash(
                        namespace, "block", chromosomes[source_index]
                    ),
                    "stratum": lineage,
                    "observed": float(observed[peak_position]),
                    "fold": fold,
                }
            )
    row_hashes = [str(row["row_hash"]) for row in outcomes]
    if len(set(row_hashes)) != len(row_hashes):
        raise RNAATACCoboltFiveSeedEvaluationError("outcome row hashes are not unique")
    genomic_blocks = {str(row["block_hash"]) for row in outcomes}
    if len(genomic_blocks) < 2:
        raise RNAATACCoboltFiveSeedEvaluationError(
            "smoke outcome has fewer than two chromosome blocks"
        )
    return outcomes, {
        "n_donors": len(set(donors)),
        "n_donor_lineage_profiles": len(groups),
        "n_selected_peaks": selected_peak_count,
        "n_genomic_blocks": len(genomic_blocks),
        "n_oof_prediction_rows_per_model_seed": len(outcomes),
    }


def _load_seed_values(
    *, cells: Mapping[int, RunCell], outcomes: Sequence[Mapping[str, Any]]
) -> Any:
    import numpy as np

    positions = {str(row["row_hash"]): index for index, row in enumerate(outcomes)}
    values = np.empty(len(outcomes), dtype=np.float64)
    seen = np.zeros(len(outcomes), dtype=np.bool_)
    model_id = next(iter(cells.values())).model_id
    for fold in FOLDS:
        cell = cells[fold]
        for row in base._load_bundle(
            cell.prediction_bundle, model_id=model_id, fold=fold
        ):
            position = positions.get(row["row_hash"])
            if position is None or seen[position]:
                raise RNAATACCoboltFiveSeedEvaluationError(
                    "prediction row is unexpected or duplicated"
                )
            outcome = outcomes[position]
            if (
                outcome["fold"] != fold
                or row["donor_hash"] != outcome["donor_hash"]
                or row["block_hash"] != outcome["block_hash"]
                or row["stratum"] != outcome["stratum"]
            ):
                raise RNAATACCoboltFiveSeedEvaluationError(
                    "prediction donor, fold, block, or lineage differs"
                )
            try:
                predicted = float(row["predicted"])
            except (TypeError, ValueError) as error:
                raise RNAATACCoboltFiveSeedEvaluationError(
                    "prediction is not numeric"
                ) from error
            if not math.isfinite(predicted) or predicted <= 0:
                raise RNAATACCoboltFiveSeedEvaluationError(
                    "prediction must be finite and positive"
                )
            values[position] = predicted
            seen[position] = True
    if not bool(seen.all()):
        raise RNAATACCoboltFiveSeedEvaluationError(
            "prediction seed does not cover the exact OOF outcome universe"
        )
    return values


def ensemble_mean(seed_values: Sequence[Sequence[float]]) -> list[float]:
    """Pure helper used by tests; production scoring uses the same arithmetic mean."""

    if len(seed_values) != len(SEEDS) or not seed_values:
        raise RNAATACCoboltFiveSeedEvaluationError("ensemble requires exactly five seeds")
    width = len(seed_values[0])
    if width == 0 or any(len(values) != width for values in seed_values):
        raise RNAATACCoboltFiveSeedEvaluationError(
            "ensemble seed prediction widths differ"
        )
    result = []
    for index in range(width):
        values = [float(seed[index]) for seed in seed_values]
        if any(not math.isfinite(value) or value <= 0 for value in values):
            raise RNAATACCoboltFiveSeedEvaluationError(
                "ensemble predictions must be finite and positive"
            )
        result.append(sum(values) / len(values))
    return result


def _score_values(
    *, outcomes: Sequence[Mapping[str, Any]], values: Any, model_id: str
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    import numpy as np
    from scipy.stats import spearmanr

    if len(values) != len(outcomes):
        raise RNAATACCoboltFiveSeedEvaluationError("score vector width differs")
    groups: dict[tuple[str, str], list[int]] = {}
    for index, row in enumerate(outcomes):
        groups.setdefault((str(row["donor_hash"]), str(row["stratum"])), []).append(
            index
        )
    block_deviances: list[float] = []
    auprcs: list[float] = []
    correlations: list[float] = []
    lineage_deviances: dict[str, list[float]] = {
        lineage: [] for lineage in base.LINEAGES
    }
    profile_rows: list[dict[str, Any]] = []
    for (donor_hash, stratum), indices in sorted(groups.items()):
        positions = np.asarray(indices, dtype=np.int64)
        observed = np.asarray(
            [float(outcomes[index]["observed"]) for index in indices],
            dtype=np.float64,
        )
        predicted = np.asarray(values[positions], dtype=np.float64)
        auprc = base._average_precision(observed > 0, predicted)
        correlation = float(spearmanr(observed, predicted).statistic)
        if not math.isfinite(correlation):
            correlation = 0.0
        auprcs.append(auprc)
        correlations.append(correlation)
        blocks: dict[str, list[int]] = {}
        for local_index, global_index in enumerate(indices):
            blocks.setdefault(str(outcomes[global_index]["block_hash"]), []).append(
                local_index
            )
        for block_hash, local_indices in sorted(blocks.items()):
            local = np.asarray(local_indices, dtype=np.int64)
            block_observed = observed[local]
            if float(block_observed.sum()) <= 0:
                continue
            deviance = base._multinomial_deviance(
                block_observed, predicted[local]
            )
            if not math.isfinite(deviance):
                raise RNAATACCoboltFiveSeedEvaluationError(
                    "profile deviance is not finite"
                )
            block_deviances.append(deviance)
            lineage_deviances[stratum].append(deviance)
            profile_rows.append(
                {
                    "model_id": model_id,
                    "donor_hash": donor_hash,
                    "block_hash": block_hash,
                    "stratum": stratum,
                    "deviance": format(deviance, ".17g"),
                    "peak_auprc": format(auprc, ".17g"),
                    "profile_spearman": format(correlation, ".17g"),
                }
            )
    if not block_deviances or not auprcs:
        raise RNAATACCoboltFiveSeedEvaluationError("no donor-block profiles scored")
    return {
        "total_multinomial_deviance": sum(block_deviances),
        "mean_peak_auprc": sum(auprcs) / len(auprcs),
        "mean_profile_spearman": sum(correlations) / len(correlations),
        "lineage_mean_deviance": {
            lineage: sum(items) / len(items)
            for lineage, items in lineage_deviances.items()
            if items
        },
    }, profile_rows


def select_strongest_task_native_baseline(
    ensemble_summaries: Mapping[str, Mapping[str, Any]],
) -> str:
    try:
        values = {
            model_id: float(ensemble_summaries[model_id]["total_multinomial_deviance"])
            for model_id in TASK_NATIVE_BASELINES
        }
    except (KeyError, TypeError, ValueError) as error:
        raise RNAATACCoboltFiveSeedEvaluationError(
            "frozen task-native baseline summaries are incomplete"
        ) from error
    if any(not math.isfinite(value) or value <= 0 for value in values.values()):
        raise RNAATACCoboltFiveSeedEvaluationError(
            "frozen task-native baseline deviance is invalid"
        )
    return min(TASK_NATIVE_BASELINES, key=lambda item: (values[item], item))


def comparison_summary(
    *,
    per_seed: Mapping[str, Mapping[int, Mapping[str, Any]]],
    ensembles: Mapping[str, Mapping[str, Any]],
    strongest_baseline: str,
) -> dict[str, Any]:
    ensemble_baseline = float(
        ensembles[strongest_baseline]["total_multinomial_deviance"]
    )
    ensemble_gains: dict[str, float] = {}
    seed_gains: dict[str, dict[str, float]] = {}
    positive_seed_counts: dict[str, int] = {}
    lineage_gains: dict[str, dict[str, float]] = {}
    for model_id in MODEL_IDS:
        ensemble_gains[model_id] = relative_deviance_reduction(
            ensemble_baseline,
            float(ensembles[model_id]["total_multinomial_deviance"]),
        )
        seed_gains[model_id] = {}
        for seed in SEEDS:
            seed_gains[model_id][str(seed)] = relative_deviance_reduction(
                float(
                    per_seed[strongest_baseline][seed][
                        "total_multinomial_deviance"
                    ]
                ),
                float(per_seed[model_id][seed]["total_multinomial_deviance"]),
            )
        positive_seed_counts[model_id] = sum(
            value > 0 for value in seed_gains[model_id].values()
        )
        lineage_gains[model_id] = {
            lineage: relative_deviance_reduction(
                float(
                    ensembles[strongest_baseline]["lineage_mean_deviance"][
                        lineage
                    ]
                ),
                float(ensembles[model_id]["lineage_mean_deviance"][lineage]),
            )
            for lineage in base.LINEAGES
        }
    return {
        "strongest_task_native_baseline": strongest_baseline,
        "ensemble_relative_deviance_reduction_vs_strongest_task_native": ensemble_gains,
        "per_seed_relative_deviance_reduction_vs_strongest_task_native": seed_gains,
        "positive_gain_seed_count": positive_seed_counts,
        "ensemble_lineage_relative_deviance_reduction_vs_strongest_task_native": lineage_gains,
    }


def _artifact_record(path: Path, *, relative_to: Path) -> dict[str, Any]:
    return {
        "path": path.relative_to(relative_to).as_posix(),
        "sha256": base._sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def evaluate(
    *,
    candidate: Path,
    execution_root: Path,
    candidate_sha256: str,
    h5_path: Path,
    h5_sha256: str,
    output: Path,
) -> dict[str, Any]:
    import numpy as np

    if output.exists():
        raise RNAATACCoboltFiveSeedEvaluationError(
            f"evaluation output already exists: {output}"
        )
    grid = discover_verified_run_grid(
        candidate=candidate,
        execution_root=execution_root,
        candidate_sha256=candidate_sha256,
    )
    inventory, n_peaks, namespace, fit_receipts = _prediction_inventory(grid)
    # This is the first outcome access. Every attempt and prediction file is
    # already verified above.
    outcomes, dimensions = _load_outcomes(
        h5_path=h5_path,
        h5_sha256=h5_sha256,
        selected_peak_count=n_peaks,
        namespace=namespace,
    )

    output.mkdir(parents=True, exist_ok=False)
    inventory_path = output / "per_seed_prediction_artifacts.tsv"
    base._write_tsv(
        inventory_path,
        (
            "model_id",
            "seed",
            "fold",
            "run_id",
            "attempt_path",
            "prediction_bundle_path",
            "prediction_bundle_sha256",
            "prediction_table_path",
            "prediction_table_sha256",
            "n_predictions",
        ),
        inventory,
    )
    ensemble_prediction_path = output / "five_seed_ensemble_predictions.tsv"
    ensemble_handle = ensemble_prediction_path.open(
        "x", encoding="utf-8", newline=""
    )
    ensemble_writer = csv.DictWriter(
        ensemble_handle,
        fieldnames=(
            "model_id",
            "row_hash",
            "donor_hash",
            "block_hash",
            "stratum",
            "predicted",
        ),
        delimiter="\t",
        lineterminator="\n",
    )
    ensemble_writer.writeheader()

    per_seed_summaries: dict[str, dict[int, dict[str, Any]]] = {
        model_id: {} for model_id in MODEL_IDS
    }
    ensemble_summaries: dict[str, dict[str, Any]] = {}
    per_seed_profile_rows: list[dict[str, Any]] = []
    ensemble_profile_rows: list[dict[str, Any]] = []
    try:
        for model_id in MODEL_IDS:
            ensemble_values = np.zeros(len(outcomes), dtype=np.float64)
            for seed in SEEDS:
                values = _load_seed_values(
                    cells=grid[model_id][seed], outcomes=outcomes
                )
                ensemble_values += values / len(SEEDS)
                summary, profiles = _score_values(
                    outcomes=outcomes,
                    values=values,
                    model_id=model_id,
                )
                per_seed_summaries[model_id][seed] = summary
                per_seed_profile_rows.extend(
                    {**row, "seed": seed} for row in profiles
                )
            ensemble_summary, profiles = _score_values(
                outcomes=outcomes,
                values=ensemble_values,
                model_id=model_id,
            )
            ensemble_summaries[model_id] = ensemble_summary
            ensemble_profile_rows.extend(profiles)
            for outcome, predicted in zip(outcomes, ensemble_values, strict=True):
                ensemble_writer.writerow(
                    {
                        "model_id": model_id,
                        "row_hash": outcome["row_hash"],
                        "donor_hash": outcome["donor_hash"],
                        "block_hash": outcome["block_hash"],
                        "stratum": outcome["stratum"],
                        "predicted": format(float(predicted), ".17g"),
                    }
                )
    finally:
        ensemble_handle.close()

    strongest = select_strongest_task_native_baseline(ensemble_summaries)
    comparisons = comparison_summary(
        per_seed=per_seed_summaries,
        ensembles=ensemble_summaries,
        strongest_baseline=strongest,
    )
    per_seed_metrics_path = output / "per_seed_metrics.tsv"
    base._write_tsv(
        per_seed_metrics_path,
        (
            "model_id",
            "seed",
            "total_multinomial_deviance",
            "mean_peak_auprc",
            "mean_profile_spearman",
            "relative_deviance_reduction_vs_strongest_task_native",
        ),
        (
            {
                "model_id": model_id,
                "seed": seed,
                "total_multinomial_deviance": format(
                    float(
                        per_seed_summaries[model_id][seed][
                            "total_multinomial_deviance"
                        ]
                    ),
                    ".17g",
                ),
                "mean_peak_auprc": format(
                    float(per_seed_summaries[model_id][seed]["mean_peak_auprc"]),
                    ".17g",
                ),
                "mean_profile_spearman": format(
                    float(
                        per_seed_summaries[model_id][seed][
                            "mean_profile_spearman"
                        ]
                    ),
                    ".17g",
                ),
                "relative_deviance_reduction_vs_strongest_task_native": format(
                    comparisons[
                        "per_seed_relative_deviance_reduction_vs_strongest_task_native"
                    ][model_id][str(seed)],
                    ".17g",
                ),
            }
            for model_id in MODEL_IDS
            for seed in SEEDS
        ),
    )
    per_seed_profile_path = output / "per_seed_profile_metrics.tsv"
    base._write_tsv(
        per_seed_profile_path,
        (
            "model_id",
            "seed",
            "donor_hash",
            "block_hash",
            "stratum",
            "deviance",
            "peak_auprc",
            "profile_spearman",
        ),
        sorted(
            per_seed_profile_rows,
            key=lambda row: (
                row["model_id"],
                row["seed"],
                row["donor_hash"],
                row["block_hash"],
                row["stratum"],
            ),
        ),
    )
    ensemble_profile_path = output / "five_seed_ensemble_profile_metrics.tsv"
    base._write_tsv(
        ensemble_profile_path,
        scpair_v1.PROFILE_FIELDS,
        sorted(
            ensemble_profile_rows,
            key=lambda row: (
                row["model_id"],
                row["donor_hash"],
                row["block_hash"],
                row["stratum"],
            ),
        ),
    )

    block_index = {
        (
            str(row["model_id"]),
            str(row["donor_hash"]),
            str(row["block_hash"]),
            str(row["stratum"]),
        ): float(row["deviance"])
        for row in ensemble_profile_rows
    }
    baseline_keys = {
        (donor, block, stratum): value
        for (model, donor, block, stratum), value in block_index.items()
        if model == strongest
    }
    paired_path = output / "ensemble_vs_strongest_task_native_donor_block.tsv"
    paired_rows = []
    for model_id in MODEL_IDS:
        model_keys = {
            (donor, block, stratum): value
            for (model, donor, block, stratum), value in block_index.items()
            if model == model_id
        }
        if set(model_keys) != set(baseline_keys):
            raise RNAATACCoboltFiveSeedEvaluationError(
                f"donor-block universe differs for {model_id}"
            )
        for key in sorted(baseline_keys):
            baseline_value = baseline_keys[key]
            model_value = model_keys[key]
            paired_rows.append(
                {
                    "model_id": model_id,
                    "baseline_model_id": strongest,
                    "donor_hash": key[0],
                    "block_hash": key[1],
                    "stratum": key[2],
                    "model_deviance": format(model_value, ".17g"),
                    "baseline_deviance": format(baseline_value, ".17g"),
                    "relative_deviance_reduction": format(
                        relative_deviance_reduction(baseline_value, model_value),
                        ".17g",
                    ),
                }
            )
    base._write_tsv(
        paired_path,
        (
            "model_id",
            "baseline_model_id",
            "donor_hash",
            "block_hash",
            "stratum",
            "model_deviance",
            "baseline_deviance",
            "relative_deviance_reduction",
        ),
        paired_rows,
    )

    cobolt_lineage = comparisons[
        "ensemble_lineage_relative_deviance_reduction_vs_strongest_task_native"
    ]["cobolt"]
    cobolt_gain = comparisons[
        "ensemble_relative_deviance_reduction_vs_strongest_task_native"
    ]["cobolt"]
    gate = {
        "minimum_relative_deviance_reduction": 0.05,
        "minimum_improved_major_lineages": 4,
        "maximum_any_lineage_worsening": 0.02,
        "required_positive_gain_seeds": 4,
        "relative_deviance_threshold_passed": cobolt_gain >= 0.05,
        "improved_lineage_count": sum(value > 0 for value in cobolt_lineage.values()),
        "no_lineage_worse_by_more_than_two_percent": all(
            value >= -0.02 for value in cobolt_lineage.values()
        ),
        "positive_gain_seed_count": comparisons["positive_gain_seed_count"][
            "cobolt"
        ],
    }
    gate["all_development_criteria_passed"] = (
        gate["relative_deviance_threshold_passed"]
        and gate["improved_lineage_count"] >= 4
        and gate["no_lineage_worse_by_more_than_two_percent"]
        and gate["positive_gain_seed_count"] >= 4
    )
    result = {
        "schema_version": "masld-bench-rna-atac-cobolt-five-seed-development-evaluation-v1",
        "campaign_id": CAMPAIGN_ID,
        "campaign_artifact_sha256": candidate_sha256,
        "task_id": base.TASK_ID,
        "dataset_id": base.DATASET_ID,
        "dataset_view_id": base.VIEW_ID,
        "input_h5_sha256": h5_sha256,
        "evaluation_roster": list(MODEL_IDS),
        "task_native_baseline_roster": list(TASK_NATIVE_BASELINES),
        "seeds": list(SEEDS),
        "folds": list(FOLDS),
        "ensemble_method": ENSEMBLE_LABEL,
        "verified_prediction_attempt_count": EXPECTED_RUNS,
        "all_prediction_attempts_verified_before_outcome_parsing_or_join": True,
        "unit_of_inference": "donor",
        "genomic_resampling_unit": "chromosome_block_smoke_only",
        "seeds_used_as_biological_replicates": False,
        "seed_role": "computational_stability_only",
        "cells_used_as_independent_replicates": False,
        "held_atac_exposed_to_fit_or_predict": False,
        "outcomes_joined_only_by_independent_evaluator": True,
        **dimensions,
        "per_seed_models": {
            model_id: {str(seed): value for seed, value in records.items()}
            for model_id, records in per_seed_summaries.items()
        },
        "five_seed_ensemble_models": ensemble_summaries,
        **comparisons,
        "cobolt_task_development_gate": gate,
        **fit_receipts,
        "comparison_policy": (
            "The strongest task-native baseline is selected once from five-seed "
            "ensemble total multinomial deviance among the frozen assay-native "
            "pseudobulk, mean-track, and shrunken-pseudobulk roster. Every model "
            "is compared with that same baseline. Seed-specific directions are "
            "stability checks and are not inferential replicates."
        ),
        "smoke_only": True,
        "champion_claim_allowed": False,
        "conditional_model_trigger_eligible": False,
        "conditional_model_trigger_ineligibility_reason": (
            "Single-study source-wide consensus-peak smoke fixture."
        ),
        "artifacts": {
            "per_seed_prediction_artifacts": _artifact_record(
                inventory_path, relative_to=output
            ),
            "five_seed_ensemble_predictions": _artifact_record(
                ensemble_prediction_path, relative_to=output
            ),
            "per_seed_metrics": _artifact_record(
                per_seed_metrics_path, relative_to=output
            ),
            "per_seed_profile_metrics": _artifact_record(
                per_seed_profile_path, relative_to=output
            ),
            "five_seed_ensemble_profile_metrics": _artifact_record(
                ensemble_profile_path, relative_to=output
            ),
            "ensemble_vs_strongest_task_native_donor_block": _artifact_record(
                paired_path, relative_to=output
            ),
        },
    }
    base._write_json(output / "evaluation.json", result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--execution-root", required=True, type=Path)
    parser.add_argument("--candidate-sha256", required=True)
    parser.add_argument("--h5", required=True, type=Path)
    parser.add_argument("--h5-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args(argv)
    result = evaluate(
        candidate=arguments.candidate,
        execution_root=arguments.execution_root,
        candidate_sha256=arguments.candidate_sha256,
        h5_path=arguments.h5,
        h5_sha256=arguments.h5_sha256,
        output=arguments.output,
    )
    print(base._canonical_json(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
