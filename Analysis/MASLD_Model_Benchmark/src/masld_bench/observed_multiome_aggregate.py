"""Complete-rectangle aggregation for the GSE296875 observed-multiome lane."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from statistics import fmean
from typing import Any, Mapping, Sequence

import h5py
import numpy as np

from .adapters.observed_multiome_factorized import MODEL_IDS
from .artifacts import reject_symlink_components, verify_frozen_tree
from .evaluators.stats import two_way_donor_block_bootstrap
from .observed_multiome_surface import decode, profile_deviance_skill, read_csr


SEEDS = (20260824, 20260825, 20260826, 20260827, 20260828)
GENOMIC_FOLDS = (0, 1, 2, 3, 4)
DONOR_FOLDS = (0, 1, 2, 3, 4)
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
BOOTSTRAP_SEED = 20260825
BOOTSTRAP_REPLICATES = 10_000


class ObservedMultiomeAggregateError(ValueError):
    """Raised when coverage, joins, values, or inference differ."""


@dataclass(frozen=True, slots=True)
class SurfaceRecord:
    seed: int
    held_genomic_fold: int
    held_donor_fold: int
    unit_hash: tuple[str, ...]
    lineage: tuple[str, ...]
    observed: np.ndarray
    predictions: Mapping[str, np.ndarray]


def _digest(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ObservedMultiomeAggregateError("JSON object required")
    return value


def expected_surface_keys() -> set[tuple[int, int, int]]:
    return {(seed, genomic, donor) for seed in SEEDS for genomic in GENOMIC_FOLDS for donor in DONOR_FOLDS}


def preflight_bundles(root: Path, bundle_records: Sequence[Mapping[str, Any]]) -> dict[tuple[int, int, int], Path]:
    """Verify read-only complete coverage without deserializing metric arrays."""

    if len(bundle_records) != len(SEEDS):
        raise ObservedMultiomeAggregateError("exactly five seed bundles are required")
    surfaces: dict[tuple[int, int, int], Path] = {}
    seen_seeds: set[int] = set()
    for record in bundle_records:
        bundle = reject_symlink_components(root / str(record.get("path", "")), label="seed bundle").resolve(strict=True)
        bundle.relative_to(root)
        verify_frozen_tree(bundle)
        if _digest(bundle / "ARTIFACTS.json") != record.get("artifacts_sha256"):
            raise ObservedMultiomeAggregateError("seed bundle drifted")
        receipt = _json(bundle / "receipt.json")
        seed = int(receipt.get("seed", -1))
        if seed not in SEEDS or seed in seen_seeds or receipt.get("surface_count") != 25 or receipt.get("all_surfaces_complete") is not True or receipt.get("partial_ranking_authorized") is not False:
            raise ObservedMultiomeAggregateError("seed bundle disposition differs")
        seen_seeds.add(seed)
        rows = receipt.get("surfaces")
        if not isinstance(rows, list) or len(rows) != 25:
            raise ObservedMultiomeAggregateError("seed bundle surface roster differs")
        for child in rows:
            key = (int(child["seed"]), int(child["held_genomic_fold"]), int(child["held_donor_fold"]))
            if key in surfaces:
                raise ObservedMultiomeAggregateError("surface key is duplicated")
            surface = reject_symlink_components(bundle / str(child["path"]), label="surface").resolve(strict=True)
            surface.relative_to(bundle)
            verify_frozen_tree(surface)
            if _digest(surface / "ARTIFACTS.json") != child.get("artifacts_sha256"):
                raise ObservedMultiomeAggregateError("surface drifted")
            surface_receipt = _json(surface / "receipt.json")
            if surface_receipt.get("surface") != {"held_genomic_fold": key[1], "held_donor_fold": key[2], "seed": key[0], "outer_training_donors": surface_receipt["surface"]["outer_training_donors"], "outer_training_donor_lineage_units": surface_receipt["surface"]["outer_training_donor_lineage_units"], "held_donors": surface_receipt["surface"]["held_donors"], "held_donor_lineage_units": surface_receipt["surface"]["held_donor_lineage_units"], "training_targets": 1000, "held_targets": 1000}:
                raise ObservedMultiomeAggregateError("surface receipt identity differs")
            if surface_receipt.get("prediction_committed_before_evaluator") is not True or surface_receipt.get("unit_metrics_written") is not True or surface_receipt.get("partial_ranking_authorized") is not False:
                raise ObservedMultiomeAggregateError("surface firewall disposition differs")
            surfaces[key] = surface
    if seen_seeds != set(SEEDS) or set(surfaces) != expected_surface_keys():
        raise ObservedMultiomeAggregateError("complete rectangle coverage is required")
    return surfaces


def _load_surface(surface: Path, key: tuple[int, int, int], observed_cache: dict[tuple[int, int], tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...], tuple[str, ...], np.ndarray]]) -> SurfaceRecord:
    seed, genomic, donor = key
    predictions = surface / "predictions"
    verify_frozen_tree(predictions)
    identifiers = _json(predictions / "identifiers.json")
    row_hash = tuple(str(value) for value in identifiers.get("row_hash", ()))
    unit_hash = tuple(str(value) for value in identifiers.get("unit_hash", ()))
    lineage = tuple(str(value) for value in identifiers.get("lineage", ()))
    target_hash = tuple(str(value) for value in identifiers.get("target_hash", ()))
    if not row_hash or len(row_hash) != len(unit_hash) or len(row_hash) != len(lineage) or len(target_hash) != 1000:
        raise ObservedMultiomeAggregateError("prediction identifiers differ")
    model_predictions: dict[str, np.ndarray] = {}
    for model_id in MODEL_IDS:
        with (predictions / f"{model_id}.npy").open("rb") as handle:
            values = np.load(handle, allow_pickle=False)
        if values.shape != (len(row_hash), len(target_hash)) or not np.isfinite(values).all() or np.any(values < 0):
            raise ObservedMultiomeAggregateError("prediction matrix differs")
        model_predictions[model_id] = np.asarray(values, dtype=np.float64)

    cache_key = (genomic, donor)
    if cache_key not in observed_cache:
        binding = _json(surface / "bindings/evaluate.json")
        outcome_record = binding.get("evaluator_outcome")
        if not isinstance(outcome_record, dict):
            raise ObservedMultiomeAggregateError("evaluator outcome binding is missing")
        evaluator = Path(str(outcome_record.get("path", ""))).resolve(strict=True)
        project_root = surface.parents[2]
        evaluator.relative_to(project_root)
        verify_frozen_tree(evaluator)
        if _digest(evaluator / "ARTIFACTS.json") != outcome_record.get("artifacts_sha256"):
            raise ObservedMultiomeAggregateError("evaluator outcome drifted")
        with h5py.File(evaluator / "evaluator_outcomes.h5", "r") as handle:
            all_rows = decode(handle["rows/row_hash"][:], unique=True)
            all_units = decode(handle["rows/unit_hash"][:], unique=False)
            all_lineages = decode(handle["rows/lineage"][:], unique=False)
            folds = np.asarray(handle["rows/donor_fold"][:], dtype=int)
            query = np.flatnonzero(folds == donor)
            all_targets = decode(handle["target_atac/target_hash"][:], unique=True)
            observed = read_csr(handle["target_atac/counts_csr"])[query].toarray()
        observed_cache[cache_key] = (
            tuple(all_rows[index] for index in query),
            tuple(all_units[index] for index in query),
            tuple(all_lineages[index] for index in query),
            tuple(all_targets),
            np.asarray(observed, dtype=np.float64),
        )
    cached_rows, cached_units, cached_lineages, cached_targets, observed = observed_cache[cache_key]
    if (row_hash, unit_hash, lineage, target_hash) != (cached_rows, cached_units, cached_lineages, cached_targets):
        raise ObservedMultiomeAggregateError("prediction-to-outcome join differs")
    return SurfaceRecord(seed=seed, held_genomic_fold=genomic, held_donor_fold=donor, unit_hash=unit_hash, lineage=lineage, observed=observed, predictions=model_predictions)


def load_complete_records(surfaces: Mapping[tuple[int, int, int], Path]) -> list[SurfaceRecord]:
    if set(surfaces) != expected_surface_keys():
        raise ObservedMultiomeAggregateError("value loading requires a complete rectangle")
    cache: dict[tuple[int, int], tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...], tuple[str, ...], np.ndarray]] = {}
    return [_load_surface(surfaces[key], key, cache) for key in sorted(surfaces)]


def aggregate_records(records: Sequence[SurfaceRecord], *, pseudocount: float = 1e-8, n_resamples: int = BOOTSTRAP_REPLICATES, bootstrap_seed: int = BOOTSTRAP_SEED) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    indexed = {(row.seed, row.held_genomic_fold, row.held_donor_fold): row for row in records}
    if len(indexed) != len(records) or set(indexed) != expected_surface_keys():
        raise ObservedMultiomeAggregateError("aggregation requires the exact complete rectangle")
    unit_donors: list[str] = []
    unit_lineages: list[str] = []
    unit_blocks: list[int] = []
    ensemble_skills: dict[str, list[float]] = {model_id: [] for model_id in MODEL_IDS}
    seed_skills: dict[int, dict[str, list[float]]] = {seed: {model_id: [] for model_id in MODEL_IDS} for seed in SEEDS}

    for genomic in GENOMIC_FOLDS:
        for donor_fold in DONOR_FOLDS:
            group = [indexed[(seed, genomic, donor_fold)] for seed in SEEDS]
            reference = group[0]
            for row in group[1:]:
                if row.unit_hash != reference.unit_hash or row.lineage != reference.lineage or not np.array_equal(row.observed, reference.observed):
                    raise ObservedMultiomeAggregateError("seed-specific surface joins differ")
            if set(reference.lineage) != set(LINEAGES):
                raise ObservedMultiomeAggregateError("lineage roster differs")
            unit_donors.extend(reference.unit_hash)
            unit_lineages.extend(reference.lineage)
            unit_blocks.extend([genomic] * len(reference.unit_hash))
            for model_id in MODEL_IDS:
                matrices = [np.asarray(row.predictions[model_id], dtype=np.float64) for row in group]
                if any(matrix.shape != reference.observed.shape for matrix in matrices):
                    raise ObservedMultiomeAggregateError("seed prediction shape differs")
                ensemble = np.mean(np.stack(matrices), axis=0)
                skill, _, _ = profile_deviance_skill(reference.observed, ensemble, pseudocount=pseudocount)
                if not np.isfinite(skill).all():
                    raise ObservedMultiomeAggregateError("complete unit evaluation is required")
                ensemble_skills[model_id].extend(float(value) for value in skill)
                for seed, matrix in zip(SEEDS, matrices, strict=True):
                    seed_skill, _, _ = profile_deviance_skill(reference.observed, matrix, pseudocount=pseudocount)
                    if not np.isfinite(seed_skill).all():
                        raise ObservedMultiomeAggregateError("complete seed evaluation is required")
                    seed_skills[seed][model_id].extend(float(value) for value in seed_skill)

    keys = list(zip(unit_donors, unit_lineages, unit_blocks, strict=True))
    if len(keys) != 39 * 5 * 5 or len(set(keys)) != len(keys):
        raise ObservedMultiomeAggregateError("donor-lineage-block grid differs")
    donors = sorted(set(unit_donors))
    if len(donors) != 39 or set(unit_lineages) != set(LINEAGES) or set(unit_blocks) != set(GENOMIC_FOLDS):
        raise ObservedMultiomeAggregateError("independent-unit census differs")
    for donor in donors:
        for block in GENOMIC_FOLDS:
            if {lineage for unit, lineage, genomic in keys if unit == donor and genomic == block} != set(LINEAGES):
                raise ObservedMultiomeAggregateError("donor-block lineage grid is incomplete")

    summaries: dict[str, Any] = {}
    zeros = [0.0] * len(unit_donors)
    for model_id in MODEL_IDS:
        interval = two_way_donor_block_bootstrap(
            ensemble_skills[model_id],
            zeros,
            unit_donors,
            unit_blocks,
            n_resamples=n_resamples,
            seed=bootstrap_seed,
            expected_n_donors=39,
            expected_n_genomic_blocks=5,
        )
        seed_points = {str(seed): fmean(seed_skills[seed][model_id]) for seed in SEEDS}
        lineage_points = {lineage: fmean(value for value, observed_lineage in zip(ensemble_skills[model_id], unit_lineages, strict=True) if observed_lineage == lineage) for lineage in LINEAGES}
        summaries[model_id] = {
            "five_seed_ensemble_profile_deviance_skill": fmean(ensemble_skills[model_id]),
            "two_way_bootstrap_vs_uniform": interval.to_dict(),
            "seed_profile_deviance_skill": seed_points,
            "positive_uniform_gain_seeds": sum(value > 0 for value in seed_points.values()),
            "lineage_profile_deviance_skill": lineage_points,
        }
    strongest = max(MODEL_IDS, key=lambda model_id: (summaries[model_id]["five_seed_ensemble_profile_deviance_skill"], -MODEL_IDS.index(model_id)))
    receipt = {
        "schema_version": "masld-bench-observed-multiome-factorized-complete-rectangle-aggregation-v1",
        "dataset_id": "gse296875",
        "surface_seed_runs": 125,
        "prediction_matrices": 625,
        "biological_donors": 39,
        "genomic_blocks": 5,
        "lineages": list(LINEAGES),
        "unit_rows": len(unit_donors),
        "seed_ensemble": "arithmetic_mean_prediction_across_five_fixed_seeds",
        "bootstrap_replicates": n_resamples,
        "bootstrap_seed": bootstrap_seed,
        "models": summaries,
        "strongest_development_control": strongest,
        "partial_rectangle_used": False,
        "seeds_are_biological_replicates": False,
        "supports_external_claim": False,
        "supports_champion_claim": False,
        "sealed_outcomes_read": False,
    }
    arrays = {
        "unit_hash": np.asarray(unit_donors),
        "lineage": np.asarray(unit_lineages),
        "genomic_fold": np.asarray(unit_blocks, dtype=np.int8),
        "model_id": np.asarray(MODEL_IDS),
        "ensemble_profile_deviance_skill": np.column_stack([ensemble_skills[model_id] for model_id in MODEL_IDS]),
        "seed_profile_deviance_skill": np.stack([np.column_stack([seed_skills[seed][model_id] for model_id in MODEL_IDS]) for seed in SEEDS]),
    }
    return receipt, arrays
