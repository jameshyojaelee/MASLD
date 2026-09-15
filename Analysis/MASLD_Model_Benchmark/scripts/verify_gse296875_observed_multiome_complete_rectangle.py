#!/usr/bin/env python3
"""Independently rederive the complete observed-multiome rectangle."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import floor
from pathlib import Path
from random import Random
from statistics import fmean
from typing import Any, Mapping, Sequence

import h5py
import numpy as np

from masld_bench.artifacts import freeze_tree, reject_symlink_components, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-observed-multiome-complete-rectangle-independent-verification-v1"
MODELS = ("masked_modality", "observed_atac_glm", "observed_atac_only", "rna_only", "shuffled_modality")
SEEDS = (20260824, 20260825, 20260826, 20260827, 20260828)
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")


class IndependentVerificationError(ValueError):
    """Raised when any output file, join, score, or summary differs."""


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            value.update(block)
    return value.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise IndependentVerificationError("JSON object required")
    return value


def _decode(values: Sequence[Any]) -> tuple[str, ...]:
    result = tuple(value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values)
    if not result or any(not value for value in result):
        raise IndependentVerificationError("string axis differs")
    return result


def _dense_csr_rows(group: Any, rows: np.ndarray) -> np.ndarray:
    shape = tuple(int(value) for value in group["shape"][:])
    data = np.asarray(group["data"][:], dtype=np.float64)
    indices = np.asarray(group["indices"][:], dtype=np.int64)
    indptr = np.asarray(group["indptr"][:], dtype=np.int64)
    if len(shape) != 2 or indptr.shape != (shape[0] + 1,) or np.any(data < 0) or not np.isfinite(data).all():
        raise IndependentVerificationError("CSR authority differs")
    output = np.zeros((len(rows), shape[1]), dtype=np.float64)
    for target_row, source_row in enumerate(rows):
        start, end = int(indptr[source_row]), int(indptr[source_row + 1])
        output[target_row, indices[start:end]] = data[start:end]
    return output


def independent_profile_skill(observed: np.ndarray, predicted: np.ndarray, pseudocount: float = 1e-8) -> np.ndarray:
    counts = np.asarray(observed, dtype=np.float64)
    estimate = np.asarray(predicted, dtype=np.float64)
    if counts.ndim != 2 or counts.shape != estimate.shape or counts.shape[1] < 2 or np.any(counts < 0) or np.any(estimate < 0) or not np.isfinite(counts).all() or not np.isfinite(estimate).all() or pseudocount <= 0:
        raise IndependentVerificationError("profile inputs differ")
    totals = counts.sum(axis=1)
    probabilities = (estimate + pseudocount) / (estimate.sum(axis=1, keepdims=True) + pseudocount * estimate.shape[1])
    model_deviance = np.zeros(counts.shape[0], dtype=np.float64)
    null_deviance = np.zeros(counts.shape[0], dtype=np.float64)
    uniform = 1.0 / estimate.shape[1]
    for row in range(counts.shape[0]):
        positive = counts[row] > 0
        if totals[row] > 0:
            model_deviance[row] = 2.0 * np.sum(counts[row, positive] * np.log(counts[row, positive] / (totals[row] * probabilities[row, positive])))
            null_deviance[row] = 2.0 * np.sum(counts[row, positive] * np.log(counts[row, positive] / (totals[row] * uniform)))
    valid = (totals > 0) & (null_deviance > 0)
    if not np.all(valid):
        raise IndependentVerificationError("all donor-lineage units must be evaluable")
    return 1.0 - model_deviance / null_deviance


def _percentile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(float(value) for value in values)
    position = (len(ordered) - 1) * probability
    lower = floor(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def independent_two_way_interval(values: Sequence[float], donors: Sequence[str], blocks: Sequence[int], *, n_resamples: int = 10_000, seed: int = 20260825) -> dict[str, Any]:
    if not (len(values) == len(donors) == len(blocks)):
        raise IndependentVerificationError("bootstrap alignment differs")
    donor_order = list(dict.fromkeys(donors))
    block_order = list(dict.fromkeys(blocks))
    if len(donor_order) != 39 or block_order != [0, 1, 2, 3, 4]:
        raise IndependentVerificationError("bootstrap unit census differs")
    cells: dict[tuple[str, int], list[float]] = {}
    for value, donor, block in zip(values, donors, blocks, strict=True):
        cells.setdefault((donor, block), []).append(float(value))
    if set(cells) != {(donor, block) for donor in donor_order for block in block_order} or any(len(cell) != 5 for cell in cells.values()):
        raise IndependentVerificationError("complete donor-block-lineage grid required")
    differences = {key: fmean(cell) for key, cell in cells.items()}
    estimate = fmean(differences.values())
    centered = {key: value - estimate for key, value in differences.items()}
    rng = Random(seed)
    boot = []
    null = []
    for _ in range(n_resamples):
        sampled_donors = [donor_order[rng.randrange(len(donor_order))] for _ in donor_order]
        sampled_blocks = [block_order[rng.randrange(len(block_order))] for _ in block_order]
        keys = [(donor, block) for donor in sampled_donors for block in sampled_blocks]
        boot.append(fmean(differences[key] for key in keys))
        null.append(fmean(centered[key] for key in keys))
    return {
        "estimate": estimate,
        "lower": _percentile(boot, 0.025),
        "upper": _percentile(boot, 0.975),
        "confidence_level": 0.95,
        "n_resamples": n_resamples,
        "n_donors": len(donor_order),
        "n_genomic_blocks": len(block_order),
        "n_cells": len(cells),
        "n_observations": len(values),
        "seed": seed,
        "bootstrap_p_value": (sum(abs(value) >= abs(estimate) for value in null) + 1) / (n_resamples + 1),
        "probability_improvement": sum(value > 0 for value in boot) / n_resamples,
    }


def _assert_close(observed: Any, expected: Any, label: str, maximum: list[float], tolerance: float = 1e-12) -> None:
    if isinstance(expected, dict):
        if not isinstance(observed, dict) or set(observed) != set(expected):
            raise IndependentVerificationError(f"mapping roster differs: {label}")
        for key in expected:
            _assert_close(observed[key], expected[key], f"{label}.{key}", maximum, tolerance)
    elif isinstance(expected, float):
        difference = abs(float(observed) - expected)
        maximum[0] = max(maximum[0], difference)
        if difference > tolerance:
            raise IndependentVerificationError(f"numeric value differs: {label}")
    elif observed != expected:
        raise IndependentVerificationError(f"value differs: {label}")


def verify(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    if config.get("schema_version") != SCHEMA or config.get("dataset_id") != "gse296875":
        raise IndependentVerificationError("verification identity differs")
    bindings = config.get("source_bindings")
    if not isinstance(bindings, dict) or set(bindings) != {"verifier", "unit_test", "sbatch"}:
        raise IndependentVerificationError("source binding roster differs")
    for label, record in bindings.items():
        path = (root / str(record["path"])).resolve(strict=True)
        path.relative_to(root)
        if _digest(path) != record["sha256"]:
            raise IndependentVerificationError(f"source drifted: {label}")
    production_record = config.get("production_artifact")
    production = reject_symlink_components(root / str(production_record["path"]), label="production artifact").resolve(strict=True)
    production.relative_to(root)
    verify_frozen_tree(production)
    if _digest(production / "ARTIFACTS.json") != production_record["artifacts_sha256"]:
        raise IndependentVerificationError("production artifact drifted")
    manifest_record = config.get("production_manifest")
    manifest_path = (root / str(manifest_record["path"])).resolve(strict=True)
    if _digest(manifest_path) != manifest_record["sha256"]:
        raise IndependentVerificationError("production manifest drifted")
    manifest = _json(manifest_path)

    surfaces: dict[tuple[int, int, int], Path] = {}
    for bundle_record in manifest["seed_bundles"]:
        bundle = (root / bundle_record["path"]).resolve(strict=True)
        verify_frozen_tree(bundle)
        if _digest(bundle / "ARTIFACTS.json") != bundle_record["artifacts_sha256"]:
            raise IndependentVerificationError("bundle drifted")
        bundle_receipt = _json(bundle / "receipt.json")
        for child in bundle_receipt["surfaces"]:
            key = (int(child["seed"]), int(child["held_genomic_fold"]), int(child["held_donor_fold"]))
            surface = (bundle / child["path"]).resolve(strict=True)
            verify_frozen_tree(surface)
            if _digest(surface / "ARTIFACTS.json") != child["artifacts_sha256"] or key in surfaces:
                raise IndependentVerificationError("surface drifted or duplicated")
            surfaces[key] = surface
    expected_surface_keys = {(seed, genomic, donor) for seed in SEEDS for genomic in range(5) for donor in range(5)}
    if set(surfaces) != expected_surface_keys:
        raise IndependentVerificationError("surface rectangle differs")

    unit_hash: list[str] = []
    lineage: list[str] = []
    genomic_fold: list[int] = []
    ensemble_skills = {model: [] for model in MODELS}
    seed_skills = {seed: {model: [] for model in MODELS} for seed in SEEDS}
    for genomic in range(5):
        for donor_fold in range(5):
            identifiers = None
            prediction_by_seed: dict[int, dict[str, np.ndarray]] = {}
            for seed in SEEDS:
                prediction_root = surfaces[(seed, genomic, donor_fold)] / "predictions"
                verify_frozen_tree(prediction_root)
                current = _json(prediction_root / "identifiers.json")
                axes = (_decode(current["row_hash"]), _decode(current["unit_hash"]), _decode(current["lineage"]), _decode(current["target_hash"]))
                if identifiers is None:
                    identifiers = axes
                elif axes != identifiers:
                    raise IndependentVerificationError("seed prediction identifiers differ")
                prediction_by_seed[seed] = {}
                for model in MODELS:
                    with (prediction_root / f"{model}.npy").open("rb") as handle:
                        matrix = np.load(handle, allow_pickle=False)
                    if matrix.shape != (len(axes[0]), len(axes[3])) or np.any(matrix < 0) or not np.isfinite(matrix).all():
                        raise IndependentVerificationError("prediction matrix differs")
                    prediction_by_seed[seed][model] = np.asarray(matrix, dtype=np.float64)
            assert identifiers is not None
            binding = _json(surfaces[(SEEDS[0], genomic, donor_fold)] / "bindings/evaluate.json")
            outcome_record = binding["evaluator_outcome"]
            evaluator = Path(outcome_record["path"]).resolve(strict=True)
            evaluator.relative_to(root)
            verify_frozen_tree(evaluator)
            if _digest(evaluator / "ARTIFACTS.json") != outcome_record["artifacts_sha256"]:
                raise IndependentVerificationError("evaluator drifted")
            with h5py.File(evaluator / "evaluator_outcomes.h5", "r") as handle:
                all_rows = _decode(handle["rows/row_hash"][:])
                all_units = _decode(handle["rows/unit_hash"][:])
                all_lineages = _decode(handle["rows/lineage"][:])
                folds = np.asarray(handle["rows/donor_fold"][:], dtype=int)
                query = np.flatnonzero(folds == donor_fold)
                targets = _decode(handle["target_atac/target_hash"][:])
                observed = _dense_csr_rows(handle["target_atac/counts_csr"], query)
            evaluator_axes = (tuple(all_rows[index] for index in query), tuple(all_units[index] for index in query), tuple(all_lineages[index] for index in query), targets)
            if evaluator_axes != identifiers:
                raise IndependentVerificationError("prediction/evaluator join differs")
            unit_hash.extend(identifiers[1])
            lineage.extend(identifiers[2])
            genomic_fold.extend([genomic] * len(identifiers[1]))
            for model in MODELS:
                matrices = [prediction_by_seed[seed][model] for seed in SEEDS]
                ensemble = np.mean(np.stack(matrices), axis=0)
                ensemble_skills[model].extend(independent_profile_skill(observed, ensemble))
                for seed, matrix in zip(SEEDS, matrices, strict=True):
                    seed_skills[seed][model].extend(independent_profile_skill(observed, matrix))

    with np.load(production / "unit_results.npz", allow_pickle=False) as archive:
        if set(archive.files) != {"unit_hash", "lineage", "genomic_fold", "model_id", "ensemble_profile_deviance_skill", "seed_profile_deviance_skill"}:
            raise IndependentVerificationError("unit result roster differs")
        frozen = {name: np.asarray(archive[name]) for name in archive.files}
    if tuple(frozen["unit_hash"].astype(str)) != tuple(unit_hash) or tuple(frozen["lineage"].astype(str)) != tuple(lineage) or not np.array_equal(frozen["genomic_fold"], np.asarray(genomic_fold)) or tuple(frozen["model_id"].astype(str)) != MODELS:
        raise IndependentVerificationError("unit result identifiers differ")
    expected_ensemble = np.column_stack([ensemble_skills[model] for model in MODELS])
    expected_seed = np.stack([np.column_stack([seed_skills[seed][model] for model in MODELS]) for seed in SEEDS])
    maximum_array_difference = max(float(np.max(np.abs(frozen["ensemble_profile_deviance_skill"] - expected_ensemble))), float(np.max(np.abs(frozen["seed_profile_deviance_skill"] - expected_seed))))
    if maximum_array_difference > 1e-12:
        raise IndependentVerificationError("unit metric rederivation differs")

    summaries = {}
    for model in MODELS:
        ensemble_values = [float(value) for value in ensemble_skills[model]]
        seed_points = {str(seed): fmean(float(value) for value in seed_skills[seed][model]) for seed in SEEDS}
        summaries[model] = {
            "five_seed_ensemble_profile_deviance_skill": fmean(ensemble_values),
            "two_way_bootstrap_vs_uniform": independent_two_way_interval(ensemble_values, unit_hash, genomic_fold),
            "seed_profile_deviance_skill": seed_points,
            "positive_uniform_gain_seeds": sum(value > 0 for value in seed_points.values()),
            "lineage_profile_deviance_skill": {name: fmean(value for value, observed_name in zip(ensemble_values, lineage, strict=True) if observed_name == name) for name in LINEAGES},
        }
    strongest = max(MODELS, key=lambda model: (summaries[model]["five_seed_ensemble_profile_deviance_skill"], -MODELS.index(model)))
    production_receipt = _json(production / "receipt.json")
    maximum_scalar_difference = [0.0]
    _assert_close(production_receipt["models"], summaries, "models", maximum_scalar_difference)
    if production_receipt.get("strongest_development_control") != strongest or production_receipt.get("surface_seed_runs") != 125 or production_receipt.get("biological_donors") != 39 or production_receipt.get("supports_external_claim") is not False or production_receipt.get("supports_champion_claim") is not False:
        raise IndependentVerificationError("production receipt disposition differs")
    return {"schema_version": "masld-bench-observed-multiome-complete-rectangle-independent-verification-receipt-v1", "dataset_id": "gse296875", "production_artifacts_sha256": production_record["artifacts_sha256"], "surface_seed_runs_verified": 125, "prediction_matrices_rederived": 625, "biological_donors_verified": 39, "lineages_verified": list(LINEAGES), "genomic_blocks_verified": 5, "bootstrap_replicates_rederived": 10000, "maximum_absolute_unit_metric_difference": maximum_array_difference, "maximum_absolute_reported_scalar_difference": maximum_scalar_difference[0], "strongest_development_control_verified": strongest, "partial_rectangle_used": False, "independent_rederivation_passed": True, "supports_external_claim": False, "supports_champion_claim": False, "sealed_outcomes_read": False, "next_gate": "specialist_comparison_against_verified_development_control"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    output = reject_symlink_components(args.output, label="independent verification output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    receipt = verify(root, _json(config_path))
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "tests/unit/test_verify_observed_multiome_complete_rectangle.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_complete_rectangle_independent_verification", "independent_rederivation_passed": True, "supports_external_claim": False, "supports_champion_claim": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
