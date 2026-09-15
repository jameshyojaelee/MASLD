#!/usr/bin/env python3
"""Evaluate the exact task-native five-seed sequence rectangle."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from random import Random
from statistics import fmean
from typing import Any, Iterable, Mapping, Sequence

from scripts import evaluate_sequence_accessibility as base


VIEWS = (
    "bpnet",
    "chrombpnet_full",
    "chrombpnet_nobias",
    "sequence_cnn_control",
    "sequence_transformer_control",
)
SEEDS = (20260824, 20260825, 20260826, 20260827, 20260828)
LINEAGES = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)
FOLDS = tuple(range(5))
UNIT_FIELDS = (
    "candidate_view",
    "seed",
    "lineage_id",
    "outer_fold",
    "donor_hash",
    "block_hash",
    "candidate_block_deviance_per_insertion",
    "uniform_block_deviance_per_insertion",
    "training_mean_block_deviance_per_insertion",
    "candidate_within_window_deviance_per_insertion",
    "candidate_regional_count_spearman",
)


class FiveSeedEvaluationError(RuntimeError):
    """Raised when exact coverage, alignment, or outcome separation differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
    path: Path,
    rows: Iterable[Mapping[str, object]],
    fields: Sequence[str],
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def validate_candidate_roster(rows: Sequence[Mapping[str, str]]) -> None:
    expected = {f"{view}__s{seed}" for view in VIEWS for seed in SEEDS}
    observed = {row["candidate_id"] for row in rows}
    if len(rows) != 25 or observed != expected:
        raise FiveSeedEvaluationError("candidate roster is not the exact five-by-five rectangle")
    if any(int(row["seed"]) not in SEEDS for row in rows):
        raise FiveSeedEvaluationError("candidate seed roster differs")


def _stable_hash(namespace: str, kind: str, identifier: str) -> str:
    return sha256(f"{namespace}\0{kind}\0{identifier}".encode()).hexdigest()


def _prediction_groups(rows: Sequence[Mapping[str, str]]) -> dict[str, list[dict[str, Any]]]:
    validate_candidate_roster(rows)
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        view = row["candidate_id"].rsplit("__s", 1)[0]
        loaded = base.load_prediction(
            row,
            dataset_id="gse296875",
            split_id=str(row["split_id"]),
            lineage_id=str(row["lineage_id"]),
        )
        grouped[view].append(loaded)
    if set(grouped) != set(VIEWS):
        raise FiveSeedEvaluationError("loaded candidate views differ")
    for view, values in grouped.items():
        values.sort(key=lambda row: row["seed"])
        if tuple(row["seed"] for row in values) != SEEDS:
            raise FiveSeedEvaluationError(f"loaded seed roster differs for {view}")
    return grouped


def _aligned_view(
    predictions: Sequence[Mapping[str, Any]],
    keys: Sequence[tuple[str, str]],
    role: str,
) -> tuple[Any, Any, dict[int, tuple[Any, Any]]]:
    import numpy as np

    seed_values: dict[int, tuple[Any, Any]] = {}
    for prediction in predictions:
        if prediction["role_counts"][role] != len(keys):
            raise FiveSeedEvaluationError("prediction role census differs")
        offsets = [prediction["profile_by_key"].get(key, -1) for key in keys]
        if any(offset < 0 for offset in offsets):
            raise FiveSeedEvaluationError("prediction lacks an evaluation window")
        count_rows = [prediction["counts"][key] for key in keys]
        if any(row["role"] != role for row in count_rows):
            raise FiveSeedEvaluationError("prediction role differs")
        mass = np.asarray([row["mass"] for row in count_rows], dtype=np.float64)
        profiles = prediction["profiles"][np.asarray(offsets, dtype=np.int64)]
        if np.any(~np.isfinite(mass)) or np.any(mass < 0):
            raise FiveSeedEvaluationError("prediction mass differs")
        seed_values[int(prediction["seed"])] = (profiles, mass)
    if tuple(sorted(seed_values)) != SEEDS:
        raise FiveSeedEvaluationError("aligned seed roster differs")
    ensemble_profile = np.mean(
        np.stack([seed_values[seed][0] for seed in SEEDS]), axis=0
    )
    ensemble_mass = np.mean(
        np.stack([seed_values[seed][1] for seed in SEEDS]), axis=0
    )
    if not np.allclose(ensemble_profile.sum(axis=1), 1.0, rtol=0, atol=2.0e-6):
        raise FiveSeedEvaluationError("ensemble profiles are not probabilities")
    return ensemble_profile, ensemble_mass, seed_values


def evaluate_group(
    *,
    matrix_root: Path,
    matrix_row: Mapping[str, str],
    namespace: str,
) -> list[dict[str, object]]:
    import h5py
    import numpy as np

    lineage = matrix_row["lineage_id"]
    split_id = matrix_row["split_id"]
    fold = int(split_id.removeprefix("donor").split("_", 1)[0])
    if lineage not in LINEAGES or fold not in FOLDS or split_id != f"donor{fold}_genomic{fold}":
        raise FiveSeedEvaluationError("matrix group identity differs")
    manifest = (matrix_root / matrix_row["prediction_manifest"]).resolve(strict=True)
    if digest(manifest) != matrix_row["prediction_manifest_sha256"]:
        raise FiveSeedEvaluationError("prediction manifest hash differs")
    prediction_rows = base.load_prediction_manifest(manifest)
    enriched_rows = [
        {**row, "lineage_id": lineage, "split_id": split_id}
        for row in prediction_rows
    ]

    # Every prediction tree and tensor is checked before an outcome path is opened.
    grouped = _prediction_groups(enriched_rows)
    observed_root = Path(matrix_row["observed_root"]).resolve(strict=True)
    baseline_root = Path(matrix_row["baseline_root"]).resolve(strict=True)
    base.validate_evaluator_input_metadata(
        root=observed_root,
        expected_artifacts_sha256=matrix_row["observed_artifacts_sha256"],
        artifact_class="chrombpnet_development_base_profiles",
        dataset_id="gse296875",
        split_id=split_id,
        lineage_id=lineage,
    )
    base.validate_evaluator_input_metadata(
        root=baseline_root,
        expected_artifacts_sha256=matrix_row["baseline_artifacts_sha256"],
        artifact_class="sequence_training_mean_baseline",
        dataset_id="gse296875",
        split_id=split_id,
        lineage_id=lineage,
    )

    unit_rows: list[dict[str, object]] = []
    role = "valid"
    with h5py.File(observed_root / "profiles/profiles.h5", "r") as observed, h5py.File(
        baseline_root / "baseline/training_mean_baseline.h5", "r"
    ) as baseline:
        if bool(observed.attrs.get("model_inference_input_eligible")):
            raise FiveSeedEvaluationError("observed profiles are inference inputs")
        truth_group = observed[role]
        truth_ids = base._decode(truth_group["window_id"][:])
        truth_hashes = base._decode(truth_group["selection_hash"][:])
        donors = base._decode(truth_group["donor_id"][:])
        truth_counts = truth_group["counts"]
        base_group = baseline[role]
        base_ids = base._decode(base_group["window_id"][:])
        base_hashes = base._decode(base_group["selection_hash"][:])
        contigs = base._decode(base_group["contig"][:])
        if (
            truth_ids != base_ids
            or truth_hashes != base_hashes
            or len(set(donors)) != len(donors)
            or truth_counts.shape[0] != len(donors)
            or truth_counts.shape[1] != len(truth_ids)
            or truth_counts.shape[2] != 1_000
            or base_group["pseudobulk_counts"].shape != truth_counts.shape[1:]
        ):
            raise FiveSeedEvaluationError("baseline and outcome axes differ")
        keys = list(zip(truth_ids, truth_hashes, strict=True))
        training_profile = np.asarray(base_group["pseudobulk_counts"][:], dtype=np.float64)
        training_mass = training_profile.sum(axis=1)
        uniform_profile = np.ones(training_profile.shape, dtype=np.float64)
        uniform_mass = np.ones(len(keys), dtype=np.float64)
        views: dict[str, dict[str, tuple[Any, Any]]] = {}
        for view in VIEWS:
            ensemble_profile, ensemble_mass, seed_values = _aligned_view(
                grouped[view], keys, role
            )
            views[view] = {
                "ensemble": (ensemble_profile, ensemble_mass),
                **{str(seed): seed_values[seed] for seed in SEEDS},
            }
        block_positions = {
            contig: np.asarray(
                [index for index, value in enumerate(contigs) if value == contig],
                dtype=np.int64,
            )
            for contig in sorted(set(contigs))
        }
        for donor_offset, donor in enumerate(donors):
            donor_truth = np.asarray(truth_counts[donor_offset], dtype=np.float64)
            for contig, indexes in block_positions.items():
                observed_block = donor_truth[indexes, :]
                observed_mass = observed_block.sum(axis=1)
                if observed_mass.sum() <= 0:
                    raise FiveSeedEvaluationError("donor-block outcome mass is zero")
                uniform_deviance = base.multinomial_deviance_per_insertion(
                    observed_block, uniform_profile[indexes]
                )
                training_deviance = base.multinomial_deviance_per_insertion(
                    observed_block, training_profile[indexes]
                )
                for view in VIEWS:
                    for seed_label, (profile, mass) in views[view].items():
                        predicted = profile[indexes] * mass[indexes, None]
                        unit_rows.append(
                            {
                                "candidate_view": view,
                                "seed": seed_label,
                                "lineage_id": lineage,
                                "outer_fold": fold,
                                "donor_hash": _stable_hash(namespace, "donor", donor),
                                "block_hash": _stable_hash(namespace, "block", contig),
                                "candidate_block_deviance_per_insertion": format(
                                    base.multinomial_deviance_per_insertion(
                                        observed_block, predicted
                                    ),
                                    ".17g",
                                ),
                                "uniform_block_deviance_per_insertion": format(
                                    uniform_deviance, ".17g"
                                ),
                                "training_mean_block_deviance_per_insertion": format(
                                    training_deviance, ".17g"
                                ),
                                "candidate_within_window_deviance_per_insertion": format(
                                    base.within_window_deviance_per_insertion(
                                        observed_block, predicted
                                    ),
                                    ".17g",
                                ),
                                "candidate_regional_count_spearman": format(
                                    base.spearman_or_zero(observed_mass, mass[indexes]),
                                    ".17g",
                                ),
                            }
                        )
    return unit_rows


def _macro_mean(rows: Sequence[Mapping[str, Any]], field: str) -> float:
    strata: dict[tuple[str, int], list[float]] = defaultdict(list)
    for row in rows:
        strata[(str(row["lineage_id"]), int(row["outer_fold"]))].append(
            float(row[field])
        )
    if set(strata) != {(lineage, fold) for lineage in LINEAGES for fold in FOLDS}:
        raise FiveSeedEvaluationError("macro stratum roster differs")
    return fmean(fmean(values) for values in strata.values())


def _macro_fisher_z(rows: Sequence[Mapping[str, Any]], field: str) -> float:
    strata: dict[tuple[str, int], list[float]] = defaultdict(list)
    for row in rows:
        value = float(row[field])
        if not math.isfinite(value) or value < -1 or value > 1:
            raise FiveSeedEvaluationError("correlation value differs")
        clipped = min(max(value, -1 + 1.0e-12), 1 - 1.0e-12)
        strata[(str(row["lineage_id"]), int(row["outer_fold"]))].append(
            math.atanh(clipped)
        )
    if set(strata) != {(lineage, fold) for lineage in LINEAGES for fold in FOLDS}:
        raise FiveSeedEvaluationError("correlation macro stratum roster differs")
    return math.tanh(fmean(fmean(values) for values in strata.values()))


def _skill_rows(
    rows: Sequence[Mapping[str, Any]], baseline_field: str
) -> list[dict[str, Any]]:
    result = []
    for row in rows:
        candidate = float(row["candidate_block_deviance_per_insertion"])
        baseline = float(row[baseline_field])
        if not math.isfinite(candidate) or not math.isfinite(baseline) or baseline <= 0:
            raise FiveSeedEvaluationError("deviance skill inputs differ")
        result.append({**row, "skill": (baseline - candidate) / baseline})
    return result


def _macro_skill(rows: Sequence[Mapping[str, Any]]) -> float:
    strata: dict[tuple[str, int], list[float]] = defaultdict(list)
    for row in rows:
        strata[(str(row["lineage_id"]), int(row["outer_fold"]))].append(
            float(row["skill"])
        )
    if set(strata) != {(lineage, fold) for lineage in LINEAGES for fold in FOLDS}:
        raise FiveSeedEvaluationError("skill macro stratum roster differs")
    return fmean(fmean(values) for values in strata.values())


def stratified_two_way_bootstrap(
    rows: Sequence[Mapping[str, Any]],
    *,
    n_resamples: int,
    seed: int,
) -> dict[str, Any]:
    if n_resamples < 100:
        raise FiveSeedEvaluationError("bootstrap requires at least 100 resamples")
    cell: dict[tuple[str, int, str, str], list[float]] = defaultdict(list)
    for row in rows:
        key = (
            str(row["lineage_id"]),
            int(row["outer_fold"]),
            str(row["donor_hash"]),
            str(row["block_hash"]),
        )
        cell[key].append(float(row["skill"]))
    lineages_by_fold: dict[int, dict[str, tuple[list[str], list[str]]]] = defaultdict(dict)
    for fold in FOLDS:
        for lineage in LINEAGES:
            keys = [key for key in cell if key[0] == lineage and key[1] == fold]
            donors = sorted({key[2] for key in keys})
            blocks = sorted({key[3] for key in keys})
            if len(donors) < 2 or len(blocks) < 2:
                raise FiveSeedEvaluationError("bootstrap independent-unit count differs")
            if any((lineage, fold, donor, block) not in cell for donor in donors for block in blocks):
                raise FiveSeedEvaluationError("bootstrap donor-block grid is incomplete")
            lineages_by_fold[fold][lineage] = (donors, blocks)
        reference = lineages_by_fold[fold][LINEAGES[0]]
        if any(lineages_by_fold[fold][lineage] != reference for lineage in LINEAGES[1:]):
            raise FiveSeedEvaluationError("lineage donor-block grids differ within fold")
    estimate = fmean(
        fmean(
            fmean(cell[(lineage, fold, donor, block)])
            for donor in lineages_by_fold[fold][lineage][0]
            for block in lineages_by_fold[fold][lineage][1]
        )
        for fold in FOLDS
        for lineage in LINEAGES
    )
    rng = Random(seed)
    draws = []
    for _ in range(n_resamples):
        stratum_values = []
        for fold in FOLDS:
            donors, blocks = lineages_by_fold[fold][LINEAGES[0]]
            sampled_donors = [donors[rng.randrange(len(donors))] for _ in donors]
            sampled_blocks = [blocks[rng.randrange(len(blocks))] for _ in blocks]
            for lineage in LINEAGES:
                stratum_values.append(
                    fmean(
                        fmean(cell[(lineage, fold, donor, block)])
                        for donor in sampled_donors
                        for block in sampled_blocks
                    )
                )
        draws.append(fmean(stratum_values))
    draws.sort()
    lower = draws[int(0.025 * (n_resamples - 1))]
    upper = draws[int(0.975 * (n_resamples - 1))]
    return {
        "estimate": estimate,
        "lower": lower,
        "upper": upper,
        "confidence_level": 0.95,
        "n_resamples": n_resamples,
        "seed": seed,
        "probability_improvement": sum(value > 0 for value in draws) / n_resamples,
        "biological_unit": "donor",
        "genomic_unit": "held_contig_block",
        "lineage_macro_weighted": True,
        "fold_macro_weighted": True,
    }


def summarize(
    unit_rows: Sequence[Mapping[str, Any]],
    *,
    n_resamples: int,
    bootstrap_seed: int,
) -> dict[str, Any]:
    ensemble_rows = [row for row in unit_rows if row["seed"] == "ensemble"]
    baseline_fields = (
        "uniform_block_deviance_per_insertion",
        "training_mean_block_deviance_per_insertion",
    )
    baseline_macro = {
        field: _macro_mean(ensemble_rows, field) for field in baseline_fields
    }
    strongest_field = min(baseline_macro, key=lambda field: (baseline_macro[field], field))
    strongest_id = {
        "uniform_block_deviance_per_insertion": "uniform_global",
        "training_mean_block_deviance_per_insertion": "training_pseudobulk_mean",
    }[strongest_field]
    models = {}
    for view in VIEWS:
        view_ensemble = [row for row in ensemble_rows if row["candidate_view"] == view]
        skills = _skill_rows(view_ensemble, strongest_field)
        seed_scores = {}
        for seed in SEEDS:
            seed_rows = [
                row
                for row in unit_rows
                if row["candidate_view"] == view and row["seed"] == str(seed)
            ]
            seed_scores[str(seed)] = _macro_skill(
                _skill_rows(seed_rows, strongest_field)
            )
        lineage_scores = {
            lineage: fmean(
                fmean(
                    row["skill"]
                    for row in skills
                    if row["lineage_id"] == lineage and row["outer_fold"] == fold
                )
                for fold in FOLDS
            )
            for lineage in LINEAGES
        }
        models[view] = {
            "five_seed_ensemble_relative_deviance_reduction": _macro_skill(skills),
            "paired_two_way_bootstrap": stratified_two_way_bootstrap(
                skills, n_resamples=n_resamples, seed=bootstrap_seed
            ),
            "seed_relative_deviance_reduction": seed_scores,
            "positive_gain_seeds": sum(value > 0 for value in seed_scores.values()),
            "seed_stability_requirement_met": sum(
                value > 0 for value in seed_scores.values()
            )
            >= 4,
            "lineage_relative_deviance_reduction": lineage_scores,
            "improved_lineages": sum(value > 0 for value in lineage_scores.values()),
            "lineages_worse_than_minus_2_percent": sum(
                value < -0.02 for value in lineage_scores.values()
            ),
            "regional_count_spearman_fisher_z_macro": _macro_fisher_z(
                view_ensemble, "candidate_regional_count_spearman"
            ),
            "within_window_deviance_per_insertion": _macro_mean(
                view_ensemble, "candidate_within_window_deviance_per_insertion"
            ),
        }
    return {
        "strongest_training_only_baseline": strongest_id,
        "training_only_baseline_macro_deviance": baseline_macro,
        "models": models,
    }


def evaluate(
    *,
    matrix_root: Path,
    completion_lock: Path,
    completion_artifacts_sha256: str,
    output: Path,
    n_resamples: int,
    bootstrap_seed: int,
) -> dict[str, Any]:
    if output.exists():
        raise FiveSeedEvaluationError("evaluation output exists")
    matrix_root = matrix_root.resolve(strict=True)
    completion_lock = completion_lock.resolve(strict=True)
    if digest(completion_lock / "ARTIFACTS.json") != completion_artifacts_sha256:
        raise FiveSeedEvaluationError("completion-lock ARTIFACTS differs")
    completion = json.loads((completion_lock / "completion/completion_lock.json").read_text())
    if (
        completion.get("status") != "ready_both_rectangles_locked"
        or completion.get("evaluation_launch_authorized") is not True
        or completion.get("partial_ranking_authorized") is not False
        or completion.get("prediction_values_read") is not False
        or completion.get("observed_outcomes_read") is not False
    ):
        raise FiveSeedEvaluationError("completion lock does not authorize evaluation")
    matrix_summary = json.loads((matrix_root / "summary.json").read_text(encoding="utf-8"))
    matrix_rows = read_tsv(matrix_root / "matrix.tsv")
    if (
        matrix_summary.get("status") != "pass"
        or matrix_summary.get("evaluation_groups") != 25
        or len(matrix_rows) != 25
        or matrix_summary.get("test_outcomes_authorized") is not False
        or matrix_summary.get("sealed_data_authorized") is not False
        or any(int(row["candidate_count"]) != 25 for row in matrix_rows)
    ):
        raise FiveSeedEvaluationError("evaluation matrix is not exact")
    expected_groups = {
        (lineage, f"donor{fold}_genomic{fold}")
        for lineage in LINEAGES
        for fold in FOLDS
    }
    if {(row["lineage_id"], row["split_id"]) for row in matrix_rows} != expected_groups:
        raise FiveSeedEvaluationError("evaluation group roster differs")
    namespace = sha256(
        f"{completion_artifacts_sha256}\0sequence-five-seed-evaluation".encode()
    ).hexdigest()
    unit_rows: list[dict[str, object]] = []
    for row in sorted(matrix_rows, key=lambda item: (item["lineage_id"], item["split_id"])):
        unit_rows.extend(
            evaluate_group(matrix_root=matrix_root, matrix_row=row, namespace=namespace)
        )
    summary = summarize(
        unit_rows, n_resamples=n_resamples, bootstrap_seed=bootstrap_seed
    )
    output.mkdir(mode=0o750)
    write_tsv(output / "donor_block_metrics.tsv", unit_rows, UNIT_FIELDS)
    result = {
        "schema_version": "masld-bench-sequence-task-native-five-seed-evaluation-v1",
        "status": "pass",
        "dataset_id": "gse296875",
        "task_id": "sequence_native_regulatory_profile_development",
        "evaluation_role": "valid",
        "fit_families": 4,
        "candidate_views": list(VIEWS),
        "fits": 500,
        "prediction_views": 625,
        "evaluation_groups": 25,
        "seed_ensemble": "arithmetic_mean_profile_probability_and_positive_regional_mass_before_scoring",
        "seeds_are_biological_replicates": False,
        "uncertainty": "paired_two_way_donor_genomic_block_bootstrap_within_fold_then_equal_lineage_fold_macro",
        "bootstrap_replicates": n_resamples,
        "bootstrap_seed": bootstrap_seed,
        "summary": summary,
        "partial_rectangle_used": False,
        "cross_task_family_ranking_performed": False,
        "test_outcomes_read": False,
        "sealed_data_read": False,
        "external_evaluation": False,
        "champion_claim_allowed": False,
        "universal_claim_allowed": False,
        "completion_artifacts_sha256": completion_artifacts_sha256,
    }
    (output / "evaluation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix-root", type=Path, required=True)
    parser.add_argument("--completion-lock", type=Path, required=True)
    parser.add_argument("--completion-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-resamples", type=int, default=10_000)
    parser.add_argument("--bootstrap-seed", type=int, default=17)
    args = parser.parse_args()
    result = evaluate(
        matrix_root=args.matrix_root,
        completion_lock=args.completion_lock,
        completion_artifacts_sha256=args.completion_artifacts_sha256,
        output=args.output,
        n_resamples=args.bootstrap_resamples,
        bootstrap_seed=args.bootstrap_seed,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
