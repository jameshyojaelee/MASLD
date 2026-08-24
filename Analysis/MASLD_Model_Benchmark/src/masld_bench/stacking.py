"""Cross-fitted development stacking and conditional-model evidence artifacts.

This module implements the only scientific path that may produce the frozen
``development_stack_gain_bundle`` records required by the conditional-model
authorization contract.  It never consumes sealed outcomes.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import csv
import json
from math import fsum, isfinite, sqrt
from pathlib import Path
import shutil
import tempfile
from typing import Any

from .artifacts import (
    ArtifactError,
    freeze_tree,
    publish_directory_noreplace,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
    write_text_exclusive,
)
from .conditional_model import (
    STACKING_EVIDENCE_SCHEMA_VERSION,
    ComplementarityEvidence,
)
from .contracts import ArtifactRef, ContractError
from .evaluators.endpoint_power import EndpointPowerError, recompute_endpoint_bootstrap
from .hashing import canonical_sha256, sha256_file
from .tournament import (
    RNA_ATAC_TASK,
    VARIANT_TASK,
    TournamentError,
    _residual_contract_from_shortlist,
    verify_development_residual_bundle,
    verify_development_shortlist,
)


class StackingError(RuntimeError):
    """Raised when development stack evidence is incomplete or forgeable."""


DEVELOPMENT_STACK_GAIN_BUNDLE_SCHEMA_VERSION = (
    "masld-bench-development-stack-gain-bundle-v1"
)
STACK_FIT_POLICY_ID = "outer_fold_convex_two_component_task_native_mse_v1"
STACK_COMPARISON_POLICY_ID = "minimum_gain_against_each_open_component_v1"
STACK_SEEDS = (1103, 2909, 4721, 6673, 8111)
STACK_WEIGHT_GRID_INTERVALS = 1000
STACK_BOOTSTRAP_SEED = 20260821
STACK_BOOTSTRAP_RESAMPLES = 10_000
_SHA256_LENGTH = 64
_COMPONENT_FIELDS = (
    "seed",
    "row_hash",
    "unit_hash",
    "block_hash",
    "stratum",
    "outer_fold",
    "study_id",
    "observed",
    "sequence_prediction",
    "context_prediction",
)
_STACK_FIELDS = (*_COMPONENT_FIELDS, "stack_prediction")
_WEIGHT_FIELDS = (
    "seed",
    "outer_fold",
    "sequence_weight",
    "context_weight",
    "fit_row_set_sha256",
    "held_row_set_sha256",
)
_TASK_CONTRACTS: Mapping[str, Mapping[str, Any]] = {
    VARIANT_TASK: {
        "endpoint_id": "mpra_allelic_direction",
        "evaluator_id": "locus_heldout_allelic_spearman_v1",
        "gain_scale": "absolute_correlation",
        "minimum_gain": 0.02,
        "study_ids": ("gse281364",),
        "strata": (),
    },
    RNA_ATAC_TASK: {
        "endpoint_id": "rna_atac_regulatory_profile",
        "evaluator_id": "rna_atac_two_way_deviance_reduction_v1",
        "gain_scale": "relative_deviance",
        "minimum_gain": 0.05,
        "study_ids": ("gse244832", "gse296875"),
        "strata": (
            "cholangiocyte",
            "fibroblast",
            "hepatocyte",
            "macrophage",
            "t_cell",
        ),
    },
}


def _safe_path(path: str | Path, label: str, *, require_file: bool = False) -> Path:
    try:
        resolved = reject_symlink_components(Path(path), label=label).resolve(strict=True)
    except (ArtifactError, OSError) as exc:
        raise StackingError(f"invalid {label}: {exc}") from exc
    if require_file and (resolved.is_symlink() or not resolved.is_file()):
        raise StackingError(f"{label} must be a regular file")
    return resolved


def _hash(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != _SHA256_LENGTH:
        raise StackingError(f"{label} must be a lowercase SHA-256 identifier")
    if any(character not in "0123456789abcdef" for character in value):
        raise StackingError(f"{label} must be a lowercase SHA-256 identifier")
    return value


def _identifier(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or any(
        character in value for character in "\x00\t\r\n"
    ):
        raise StackingError(f"{label} must be a non-empty identifier")
    return value


def _number(value: object, label: str) -> float:
    if isinstance(value, bool):
        raise StackingError(f"{label} must be finite numeric data")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise StackingError(f"{label} must be finite numeric data") from exc
    if not isfinite(result):
        raise StackingError(f"{label} must be finite numeric data")
    return result


def _read_component_rows(path: Path, task_id: str) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != _COMPONENT_FIELDS:
                raise StackingError(
                    "component table header must be exactly "
                    + " ".join(_COMPONENT_FIELDS)
                )
            rows = [dict(row) for row in reader]
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise StackingError(f"cannot read component table: {exc}") from exc
    if not rows:
        raise StackingError("component table must not be empty")
    contract = _TASK_CONTRACTS[task_id]
    expected_studies = set(contract["study_ids"])
    expected_strata = set(contract["strata"])
    identities: set[tuple[int, str]] = set()
    metadata_by_row: dict[str, tuple[str, ...]] = {}
    seeds: set[int] = set()
    folds: set[str] = set()
    studies: set[str] = set()
    strata: set[str] = set()
    for index, row in enumerate(rows):
        try:
            seed = int(row["seed"])
        except ValueError as exc:
            raise StackingError(f"component row {index} seed must be an integer") from exc
        if str(seed) != row["seed"] or seed not in STACK_SEEDS:
            raise StackingError(f"component row {index} uses an unregistered seed")
        row_hash = _hash(row["row_hash"], f"component row {index}.row_hash")
        _hash(row["unit_hash"], f"component row {index}.unit_hash")
        _hash(row["block_hash"], f"component row {index}.block_hash")
        fold = _identifier(row["outer_fold"], f"component row {index}.outer_fold")
        study = _identifier(row["study_id"], f"component row {index}.study_id")
        if study not in expected_studies:
            raise StackingError(f"component row {index} uses disallowed study {study!r}")
        stratum = row["stratum"]
        if task_id == VARIANT_TASK:
            if stratum != "all":
                raise StackingError("MPRA proxy component rows require stratum='all'")
        elif stratum not in expected_strata:
            raise StackingError(
                f"RNA-ATAC component row {index} is outside the frozen strata roster"
            )
        observed = _number(row["observed"], f"component row {index}.observed")
        sequence = _number(
            row["sequence_prediction"],
            f"component row {index}.sequence_prediction",
        )
        context = _number(
            row["context_prediction"],
            f"component row {index}.context_prediction",
        )
        if task_id == RNA_ATAC_TASK and min(observed, sequence, context) < 0.0:
            raise StackingError("RNA-ATAC component profiles must be nonnegative")
        key = (seed, row_hash)
        if key in identities:
            raise StackingError("component seed/row_hash pairs must be unique")
        identities.add(key)
        metadata = tuple(
            row[field]
            for field in (
                "unit_hash",
                "block_hash",
                "stratum",
                "outer_fold",
                "study_id",
                "observed",
            )
        )
        previous = metadata_by_row.setdefault(row_hash, metadata)
        if previous != metadata:
            raise StackingError("component row metadata or outcomes differ across seeds")
        seeds.add(seed)
        folds.add(fold)
        studies.add(study)
        strata.add(stratum)
    if tuple(sorted(seeds)) != STACK_SEEDS:
        raise StackingError("component table must contain exactly the five fixed seeds")
    row_sets = {
        seed: {row["row_hash"] for row in rows if int(row["seed"]) == seed}
        for seed in STACK_SEEDS
    }
    if len({canonical_sha256(sorted(values)) for values in row_sets.values()}) != 1:
        raise StackingError("component row sets must be identical across seeds")
    if len(folds) < 2:
        raise StackingError("cross-fitted stacking requires at least two outer folds")
    if studies != expected_studies:
        raise StackingError("component table must cover the exact proxy study roster")
    if task_id == RNA_ATAC_TASK and strata != expected_strata:
        raise StackingError("RNA-ATAC component table lacks a frozen stratum")
    return sorted(rows, key=lambda row: (int(row["seed"]), row["row_hash"]))


def _profile_normalize(
    rows: Sequence[Mapping[str, str]], field: str
) -> dict[str, float]:
    by_profile: dict[tuple[str, str, str], list[Mapping[str, str]]] = {}
    for row in rows:
        key = (row["unit_hash"], row["block_hash"], row["stratum"])
        by_profile.setdefault(key, []).append(row)
    result: dict[str, float] = {}
    for key, profile in by_profile.items():
        total = fsum(_number(row[field], field) for row in profile)
        if total <= 0.0:
            raise StackingError(
                f"RNA-ATAC {field} profile has nonpositive mass for {key}"
            )
        for row in profile:
            result[row["row_hash"]] = _number(row[field], field) / total
    return result


def _fit_loss(
    rows: Sequence[Mapping[str, str]], task_id: str, sequence_weight: float
) -> float:
    if task_id == VARIANT_TASK:
        errors = [
            (
                _number(row["observed"], "observed")
                - sequence_weight
                * _number(row["sequence_prediction"], "sequence_prediction")
                - (1.0 - sequence_weight)
                * _number(row["context_prediction"], "context_prediction")
            )
            ** 2
            for row in rows
        ]
        return fsum(errors) / len(errors)
    observed = _profile_normalize(rows, "observed")
    sequence = _profile_normalize(rows, "sequence_prediction")
    context = _profile_normalize(rows, "context_prediction")
    by_profile: dict[tuple[str, str, str], list[str]] = {}
    for row in rows:
        key = (row["unit_hash"], row["block_hash"], row["stratum"])
        by_profile.setdefault(key, []).append(row["row_hash"])
    profile_losses = []
    for row_hashes in by_profile.values():
        errors = [
            (
                observed[row_hash]
                - sequence_weight * sequence[row_hash]
                - (1.0 - sequence_weight) * context[row_hash]
            )
            ** 2
            for row_hash in row_hashes
        ]
        profile_losses.append(fsum(errors) / len(errors))
    return fsum(profile_losses) / len(profile_losses)


def _fit_convex_weight(rows: Sequence[Mapping[str, str]], task_id: str) -> float:
    candidates = []
    for index in range(STACK_WEIGHT_GRID_INTERVALS + 1):
        weight = index / STACK_WEIGHT_GRID_INTERVALS
        candidates.append((_fit_loss(rows, task_id, weight), abs(weight - 0.5), weight))
    return min(candidates)[2]


def _tsv(fields: Sequence[str], rows: Sequence[Mapping[str, object]]) -> str:
    lines = ["\t".join(fields)]
    lines.extend("\t".join(str(row[field]) for field in fields) for row in rows)
    return "\n".join(lines) + "\n"


def _endpoint_table(
    rows: Sequence[Mapping[str, object]],
    task_id: str,
    baseline_field: str,
) -> tuple[tuple[str, ...], list[dict[str, object]]]:
    if task_id == VARIANT_TASK:
        fields = (
            "row_hash",
            "unit_hash",
            "block_hash",
            "observed",
            "candidate",
            "baseline",
        )
    else:
        fields = (
            "row_hash",
            "donor_hash",
            "block_hash",
            "stratum",
            "observed",
            "candidate",
            "baseline",
        )
    endpoint_rows = []
    for row in rows:
        common: dict[str, object] = {
            "row_hash": row["row_hash"],
            "block_hash": row["block_hash"],
            "observed": row["observed"],
            "candidate": row["stack_prediction"],
            "baseline": row[baseline_field],
        }
        if task_id == VARIANT_TASK:
            common["unit_hash"] = row["unit_hash"]
        else:
            common["donor_hash"] = row["unit_hash"]
            common["stratum"] = row["stratum"]
        endpoint_rows.append(common)
    return fields, sorted(endpoint_rows, key=lambda row: str(row["row_hash"]))


def _derive_stack(
    component_rows: Sequence[Mapping[str, str]], task_id: str
) -> dict[str, Any]:
    contract = _TASK_CONTRACTS[task_id]
    stack_rows: list[dict[str, object]] = []
    weight_rows: list[dict[str, object]] = []
    endpoint_tables: dict[str, tuple[tuple[str, ...], list[dict[str, object]]]] = {}
    endpoint_results: list[dict[str, Any]] = []
    per_seed_gains: dict[int, float] = {}
    for seed in STACK_SEEDS:
        seed_rows = [row for row in component_rows if int(row["seed"]) == seed]
        folds = sorted({row["outer_fold"] for row in seed_rows})
        seed_stack_rows: list[dict[str, object]] = []
        for fold in folds:
            held = [row for row in seed_rows if row["outer_fold"] == fold]
            held_units = {row["unit_hash"] for row in held}
            held_blocks = {row["block_hash"] for row in held}
            fit = [
                row
                for row in seed_rows
                if row["unit_hash"] not in held_units
                and row["block_hash"] not in held_blocks
            ]
            if not fit:
                raise StackingError(
                    f"outer fold {fold!r} has no donor/locus and block-disjoint fit rows"
                )
            weight = _fit_convex_weight(fit, task_id)
            if task_id == RNA_ATAC_TASK:
                sequence = _profile_normalize(held, "sequence_prediction")
                context_values = _profile_normalize(held, "context_prediction")
            else:
                sequence = {
                    row["row_hash"]: _number(
                        row["sequence_prediction"], "sequence_prediction"
                    )
                    for row in held
                }
                context_values = {
                    row["row_hash"]: _number(
                        row["context_prediction"], "context_prediction"
                    )
                    for row in held
                }
            for row in held:
                row_hash = row["row_hash"]
                normalized = dict(row)
                normalized["sequence_prediction"] = sequence[row_hash]
                normalized["context_prediction"] = context_values[row_hash]
                normalized["stack_prediction"] = (
                    weight * sequence[row_hash]
                    + (1.0 - weight) * context_values[row_hash]
                )
                seed_stack_rows.append(normalized)
            weight_rows.append(
                {
                    "seed": seed,
                    "outer_fold": fold,
                    "sequence_weight": weight,
                    "context_weight": 1.0 - weight,
                    "fit_row_set_sha256": canonical_sha256(
                        sorted(row["row_hash"] for row in fit)
                    ),
                    "held_row_set_sha256": canonical_sha256(
                        sorted(row["row_hash"] for row in held)
                    ),
                }
            )
        seed_stack_rows.sort(key=lambda row: str(row["row_hash"]))
        stack_rows.extend(seed_stack_rows)
        gains = []
        for baseline_role, baseline_field in (
            ("sequence", "sequence_prediction"),
            ("context", "context_prediction"),
        ):
            fields, table_rows = _endpoint_table(
                seed_stack_rows, task_id, baseline_field
            )
            key = f"seed-{seed}--{baseline_role}"
            endpoint_tables[key] = (fields, table_rows)
        # Evaluation happens after tables are materialized in the staging tree.
    by_row: dict[str, list[Mapping[str, object]]] = {}
    for row in stack_rows:
        by_row.setdefault(str(row["row_hash"]), []).append(row)
    ensemble_rows: list[dict[str, object]] = []
    for row_hash, rows in sorted(by_row.items()):
        if len(rows) != len(STACK_SEEDS):
            raise StackingError("stack rows are not aligned across all five seeds")
        first = rows[0]
        ensemble_rows.append(
            {
                **{field: first[field] for field in _COMPONENT_FIELDS if field != "seed"},
                "seed": "ensemble",
                "sequence_prediction": fsum(
                    _number(row["sequence_prediction"], "sequence_prediction")
                    for row in rows
                )
                / len(rows),
                "context_prediction": fsum(
                    _number(row["context_prediction"], "context_prediction")
                    for row in rows
                )
                / len(rows),
                "stack_prediction": fsum(
                    _number(row["stack_prediction"], "stack_prediction")
                    for row in rows
                )
                / len(rows),
            }
        )
    for baseline_role, baseline_field in (
        ("sequence", "sequence_prediction"),
        ("context", "context_prediction"),
    ):
        fields, table_rows = _endpoint_table(ensemble_rows, task_id, baseline_field)
        endpoint_tables[f"ensemble--{baseline_role}"] = (fields, table_rows)
    return {
        "contract": contract,
        "stack_rows": sorted(
            stack_rows, key=lambda row: (int(str(row["seed"])), str(row["row_hash"]))
        ),
        "weight_rows": sorted(
            weight_rows, key=lambda row: (int(str(row["seed"])), str(row["outer_fold"]))
        ),
        "endpoint_tables": endpoint_tables,
        "endpoint_results": endpoint_results,
        "per_seed_gains": per_seed_gains,
    }


def _evaluate_stack(derived: dict[str, Any], staging: Path, task_id: str) -> None:
    contract = derived["contract"]
    parameters = {"strata": list(contract["strata"])} if task_id == RNA_ATAC_TASK else {}
    results = []
    gains_by_seed: dict[int, list[float]] = {seed: [] for seed in STACK_SEEDS}
    ensemble_gains = []
    for key, (fields, rows) in sorted(derived["endpoint_tables"].items()):
        table_path = staging / "endpoint_tables" / f"{key}.tsv"
        table_path.parent.mkdir(exist_ok=True)
        write_text_exclusive(table_path, _tsv(fields, rows), mode=0o440)
        result = recompute_endpoint_bootstrap(
            evaluator_id=contract["evaluator_id"],
            table_path=table_path,
            n_resamples=STACK_BOOTSTRAP_RESAMPLES,
            seed=STACK_BOOTSTRAP_SEED,
            parameters=parameters,
        )
        result_payload = result.to_dict()
        label, baseline_role = key.split("--", 1)
        results.append(
            {
                "evaluation_label": label,
                "baseline_role": baseline_role,
                "table_sha256": sha256_file(table_path),
                "result": result_payload,
            }
        )
        if label == "ensemble":
            ensemble_gains.append(float(result.observed_effect))
        else:
            gains_by_seed[int(label.removeprefix("seed-"))].append(
                float(result.observed_effect)
            )
    if len(ensemble_gains) != 2 or any(len(values) != 2 for values in gains_by_seed.values()):
        raise StackingError("stack evaluation did not compare both open components")
    derived["endpoint_results"] = results
    derived["per_seed_gains"] = {
        seed: min(values) for seed, values in gains_by_seed.items()
    }
    derived["ensemble_gain"] = min(ensemble_gains)


def _candidate_contracts(
    shortlist: Mapping[str, Any], sequence_candidate_id: str, context_candidate_id: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    if sequence_candidate_id == context_candidate_id:
        raise StackingError("stack components must be distinct candidates")
    try:
        sequence = _residual_contract_from_shortlist(shortlist, sequence_candidate_id)
        context = _residual_contract_from_shortlist(shortlist, context_candidate_id)
    except TournamentError as exc:
        raise StackingError(f"stack component is not open and shortlist-selected: {exc}") from exc
    shared_fields = (
        "task_id",
        "endpoint_evaluator_id",
        "endpoint_evaluator_sha256",
        "dataset_ids",
        "dataset_registry_sha256s",
        "fold",
        "row_id_field",
        "biological_unit",
    )
    for field in shared_fields:
        if sequence[field] != context[field]:
            raise StackingError(f"stack components differ in frozen {field}")
    if sequence["source_family_id"] == context["source_family_id"]:
        raise StackingError("stack components must come from distinct model families")
    return sequence, context


def _shortlist_binding(root: Path, shortlist: Mapping[str, Any]) -> dict[str, str]:
    verify_frozen_tree(root)
    return {
        "path": root.as_posix(),
        "manifest_sha256": sha256_file(root / "ARTIFACTS.json"),
        "document_sha256": sha256_file(root / "development_shortlist.json"),
        "shortlist_id": str(shortlist["shortlist_id"]),
    }


def freeze_development_stack_gain_bundle(
    *,
    development_shortlist_dir: str | Path,
    sequence_candidate_id: str,
    context_candidate_id: str,
    component_prediction_table_path: str | Path,
    output_root: str | Path,
) -> Path:
    """Fit and freeze one two-component, outer-fold cross-fitted stack."""

    shortlist_root = _safe_path(development_shortlist_dir, "development shortlist")
    try:
        shortlist = verify_development_shortlist(shortlist_root)
    except TournamentError as exc:
        raise StackingError(f"invalid development shortlist: {exc}") from exc
    task_id = str(shortlist.get("task_id"))
    if task_id not in _TASK_CONTRACTS:
        raise StackingError("stack gain bundles support only variant and RNA-ATAC tasks")
    sequence_id = _hash(sequence_candidate_id, "sequence_candidate_id")
    context_id = _hash(context_candidate_id, "context_candidate_id")
    sequence, context = _candidate_contracts(shortlist, sequence_id, context_id)
    component_path = _safe_path(
        component_prediction_table_path, "component prediction table", require_file=True
    )
    rows = _read_component_rows(component_path, task_id)
    output = reject_symlink_components(Path(output_root), label="stack output root").resolve()
    output.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".development-stack.", dir=output))
    try:
        component_copy = staging / "component_predictions.tsv"
        write_text_exclusive(component_copy, _tsv(_COMPONENT_FIELDS, rows), mode=0o440)
        derived = _derive_stack(rows, task_id)
        stack_path = staging / "stack_predictions.tsv"
        weights_path = staging / "cross_fit_weights.tsv"
        write_text_exclusive(
            stack_path, _tsv(_STACK_FIELDS, derived["stack_rows"]), mode=0o440
        )
        write_text_exclusive(
            weights_path, _tsv(_WEIGHT_FIELDS, derived["weight_rows"]), mode=0o440
        )
        _evaluate_stack(derived, staging, task_id)
        results_path = staging / "endpoint_results.json"
        write_json_exclusive(results_path, derived["endpoint_results"], mode=0o440)
        refs = {
            "component_predictions": ArtifactRef.from_path(
                component_copy,
                relative_to=staging,
                media_type="text/tab-separated-values",
                role=f"development_stack_components:{task_id}",
            ).to_dict(),
            "stack_predictions": ArtifactRef.from_path(
                stack_path,
                relative_to=staging,
                media_type="text/tab-separated-values",
                role=f"cross_fitted_stack_predictions:{task_id}",
            ).to_dict(),
            "cross_fit_weights": ArtifactRef.from_path(
                weights_path,
                relative_to=staging,
                media_type="text/tab-separated-values",
                role=f"cross_fitted_nonnegative_weights:{task_id}",
            ).to_dict(),
            "endpoint_results": ArtifactRef.from_path(
                results_path,
                relative_to=staging,
                media_type="application/json",
                role=f"development_stack_endpoint_results:{task_id}",
            ).to_dict(),
        }
        contract = derived["contract"]
        positive_seed_count = sum(
            gain > 0.0 for gain in derived["per_seed_gains"].values()
        )
        identity = {
            "schema_version": DEVELOPMENT_STACK_GAIN_BUNDLE_SCHEMA_VERSION,
            "task_id": task_id,
            "endpoint_id": contract["endpoint_id"],
            "endpoint_evaluator_id": contract["evaluator_id"],
            "endpoint_evaluator_sha256": sha256_file(
                Path(__file__).resolve().parent / "evaluators" / "endpoint_power.py"
            ),
            "gain_scale": contract["gain_scale"],
            "minimum_gain": contract["minimum_gain"],
            "observed_gain": derived["ensemble_gain"],
            "per_seed_gains": {
                str(seed): derived["per_seed_gains"][seed] for seed in STACK_SEEDS
            },
            "qualifying_seed_count": positive_seed_count,
            "evaluated_seed_count": len(STACK_SEEDS),
            "study_ids": list(contract["study_ids"]),
            "development_shortlist_binding": _shortlist_binding(
                shortlist_root, shortlist
            ),
            "sequence_candidate": sequence,
            "context_candidate": context,
            "fit_policy_id": STACK_FIT_POLICY_ID,
            "comparison_policy_id": STACK_COMPARISON_POLICY_ID,
            "weight_grid_intervals": STACK_WEIGHT_GRID_INTERVALS,
            "bootstrap_seed": STACK_BOOTSTRAP_SEED,
            "bootstrap_resamples": STACK_BOOTSTRAP_RESAMPLES,
            "cross_fitted": True,
            "nonnegative_stack": True,
            "best_open_models_compared": True,
            "created_before_outcome_unblind": True,
            "sealed_results_used": False,
            **refs,
        }
        payload = {"gain_bundle_id": canonical_sha256(identity), **identity}
        write_json_exclusive(
            staging / "development_stack_gain_bundle.json", payload, mode=0o440
        )
        freeze_tree(
            staging,
            {
                "artifact_class": "development_stack_gain_bundle",
                "gain_bundle_id": payload["gain_bundle_id"],
            },
        )
        target = output / f"development-stack-gain--{payload['gain_bundle_id']}"
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_development_stack_gain_bundle(target)
    return target


def _load_ref(payload: Mapping[str, Any], root: Path, field: str, role: str) -> Path:
    try:
        ref = ArtifactRef.from_dict(payload.get(field))
        if ref.role != role:
            raise StackingError(f"{field} has the wrong artifact role")
        return ref.validate(root, require_relative=True)
    except (ContractError, TypeError) as exc:
        raise StackingError(f"invalid {field} artifact: {exc}") from exc


def verify_development_stack_gain_bundle(path: str | Path) -> dict[str, Any]:
    """Recompute a frozen stack, every weight, and every registered endpoint."""

    root = _safe_path(path, "development stack gain bundle")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads(
            (root / "development_stack_gain_bundle.json").read_text(encoding="utf-8")
        )
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise StackingError(f"invalid development stack gain bundle: {exc}") from exc
    if not isinstance(payload, Mapping) or payload.get("schema_version") != DEVELOPMENT_STACK_GAIN_BUNDLE_SCHEMA_VERSION:
        raise StackingError("unsupported development stack gain bundle schema")
    expected_fields = {
        "gain_bundle_id",
        "schema_version",
        "task_id",
        "endpoint_id",
        "endpoint_evaluator_id",
        "endpoint_evaluator_sha256",
        "gain_scale",
        "minimum_gain",
        "observed_gain",
        "per_seed_gains",
        "qualifying_seed_count",
        "evaluated_seed_count",
        "study_ids",
        "development_shortlist_binding",
        "sequence_candidate",
        "context_candidate",
        "fit_policy_id",
        "comparison_policy_id",
        "weight_grid_intervals",
        "bootstrap_seed",
        "bootstrap_resamples",
        "cross_fitted",
        "nonnegative_stack",
        "best_open_models_compared",
        "created_before_outcome_unblind",
        "sealed_results_used",
        "component_predictions",
        "stack_predictions",
        "cross_fit_weights",
        "endpoint_results",
    }
    if set(payload) != expected_fields:
        raise StackingError("development stack gain bundle has the wrong exact schema")
    claimed = _hash(payload.get("gain_bundle_id"), "gain_bundle_id")
    identity = dict(payload)
    identity.pop("gain_bundle_id", None)
    if canonical_sha256(identity) != claimed:
        raise StackingError("development stack gain identity does not rederive")
    if manifest.get("metadata") != {
        "artifact_class": "development_stack_gain_bundle",
        "gain_bundle_id": claimed,
    }:
        raise StackingError("development stack gain manifest metadata mismatch")
    task_id = str(payload.get("task_id"))
    if task_id not in _TASK_CONTRACTS:
        raise StackingError("development stack gain has an unsupported task")
    contract = _TASK_CONTRACTS[task_id]
    fixed = {
        "endpoint_id": contract["endpoint_id"],
        "endpoint_evaluator_id": contract["evaluator_id"],
        "endpoint_evaluator_sha256": sha256_file(
            Path(__file__).resolve().parent / "evaluators" / "endpoint_power.py"
        ),
        "gain_scale": contract["gain_scale"],
        "minimum_gain": contract["minimum_gain"],
        "study_ids": list(contract["study_ids"]),
        "fit_policy_id": STACK_FIT_POLICY_ID,
        "comparison_policy_id": STACK_COMPARISON_POLICY_ID,
        "weight_grid_intervals": STACK_WEIGHT_GRID_INTERVALS,
        "bootstrap_seed": STACK_BOOTSTRAP_SEED,
        "bootstrap_resamples": STACK_BOOTSTRAP_RESAMPLES,
        "evaluated_seed_count": len(STACK_SEEDS),
        "cross_fitted": True,
        "nonnegative_stack": True,
        "best_open_models_compared": True,
        "created_before_outcome_unblind": True,
        "sealed_results_used": False,
    }
    for field, expected in fixed.items():
        if payload.get(field) != expected:
            raise StackingError(f"development stack gain changed frozen {field}")
    binding = payload.get("development_shortlist_binding")
    if not isinstance(binding, Mapping) or set(binding) != {
        "path",
        "manifest_sha256",
        "document_sha256",
        "shortlist_id",
    }:
        raise StackingError("development stack shortlist binding has the wrong schema")
    shortlist_root = _safe_path(binding["path"], "bound development shortlist")
    shortlist = verify_development_shortlist(shortlist_root)
    if _shortlist_binding(shortlist_root, shortlist) != dict(binding):
        raise StackingError("development stack shortlist binding changed")
    sequence_raw = payload.get("sequence_candidate")
    context_raw = payload.get("context_candidate")
    if not isinstance(sequence_raw, Mapping) or not isinstance(context_raw, Mapping):
        raise StackingError("development stack component contracts must be objects")
    sequence, context = _candidate_contracts(
        shortlist,
        _hash(sequence_raw.get("candidate_id"), "sequence candidate"),
        _hash(context_raw.get("candidate_id"), "context candidate"),
    )
    if payload.get("sequence_candidate") != sequence or payload.get("context_candidate") != context:
        raise StackingError("development stack component contract changed")
    component_path = _load_ref(
        payload,
        root,
        "component_predictions",
        f"development_stack_components:{task_id}",
    )
    rows = _read_component_rows(component_path, task_id)
    derived = _derive_stack(rows, task_id)
    expected_stack = _tsv(_STACK_FIELDS, derived["stack_rows"])
    expected_weights = _tsv(_WEIGHT_FIELDS, derived["weight_rows"])
    stack_path = _load_ref(
        payload,
        root,
        "stack_predictions",
        f"cross_fitted_stack_predictions:{task_id}",
    )
    weight_path = _load_ref(
        payload,
        root,
        "cross_fit_weights",
        f"cross_fitted_nonnegative_weights:{task_id}",
    )
    if stack_path.read_text(encoding="utf-8") != expected_stack:
        raise StackingError("published cross-fitted stack predictions do not rederive")
    if weight_path.read_text(encoding="utf-8") != expected_weights:
        raise StackingError("published cross-fitted stack weights do not rederive")
    # Reuse the frozen endpoint tables as a scratch-free recomputation surface.
    contract_parameters = (
        {"strata": list(contract["strata"])} if task_id == RNA_ATAC_TASK else {}
    )
    observed_results = []
    gains_by_seed: dict[int, list[float]] = {seed: [] for seed in STACK_SEEDS}
    ensemble_gains = []
    for key, (fields, table_rows) in sorted(derived["endpoint_tables"].items()):
        table_path = root / "endpoint_tables" / f"{key}.tsv"
        if table_path.is_symlink() or table_path.read_text(encoding="utf-8") != _tsv(fields, table_rows):
            raise StackingError(f"endpoint table {key} does not rederive")
        result = recompute_endpoint_bootstrap(
            evaluator_id=contract["evaluator_id"],
            table_path=table_path,
            n_resamples=STACK_BOOTSTRAP_RESAMPLES,
            seed=STACK_BOOTSTRAP_SEED,
            parameters=contract_parameters,
        )
        label, baseline_role = key.split("--", 1)
        observed_results.append(
            {
                "evaluation_label": label,
                "baseline_role": baseline_role,
                "table_sha256": sha256_file(table_path),
                "result": result.to_dict(),
            }
        )
        if label == "ensemble":
            ensemble_gains.append(float(result.observed_effect))
        else:
            gains_by_seed[int(label.removeprefix("seed-"))].append(
                float(result.observed_effect)
            )
    results_path = _load_ref(
        payload,
        root,
        "endpoint_results",
        f"development_stack_endpoint_results:{task_id}",
    )
    stored_results = json.loads(results_path.read_text(encoding="utf-8"))
    if stored_results != observed_results:
        raise StackingError("development stack endpoint results do not rederive")
    per_seed = {str(seed): min(gains_by_seed[seed]) for seed in STACK_SEEDS}
    observed_gain = min(ensemble_gains)
    qualifying = sum(value > 0.0 for value in per_seed.values())
    if (
        payload.get("per_seed_gains") != per_seed
        or payload.get("observed_gain") != observed_gain
        or payload.get("qualifying_seed_count") != qualifying
    ):
        raise StackingError("development stack gain summary does not rederive")
    return dict(payload)


def _source_binding(
    path: Path, *, role: str, artifact_class: str, filename: str
) -> dict[str, str]:
    manifest = verify_frozen_tree(path)
    if not isinstance(manifest.get("metadata"), Mapping) or manifest["metadata"].get("artifact_class") != artifact_class:
        raise StackingError(f"{role} source has the wrong artifact class")
    document = path / filename
    if document.is_symlink() or not document.is_file():
        raise StackingError(f"{role} source lacks {filename}")
    return {
        "role": role,
        "artifact_class": artifact_class,
        "path": path.as_posix(),
        "manifest_sha256": sha256_file(path / "ARTIFACTS.json"),
        "document_filename": filename,
        "document_sha256": sha256_file(document),
    }


def _residual_correlation(
    first_root: Path,
    first: Mapping[str, Any],
    second_root: Path,
    second: Mapping[str, Any],
) -> float:
    identity_fields = (
        "task_id",
        "development_shortlist_binding",
        "endpoint_evaluator_id",
        "endpoint_evaluator_sha256",
        "dataset_ids",
        "dataset_registry_sha256s",
        "fold",
        "row_id_field",
        "unit_id_namespace",
        "biological_unit",
        "n_rows",
        "row_set_sha256",
    )
    if any(first.get(field) != second.get(field) for field in identity_fields):
        raise StackingError("residual bundles do not share one held development identity")
    try:
        first_ref = ArtifactRef.from_dict(first["residual_table"])
        second_ref = ArtifactRef.from_dict(second["residual_table"])
        first_path = first_ref.validate(first_root, require_relative=True)
        second_path = second_ref.validate(second_root, require_relative=True)
    except (ContractError, KeyError, TypeError) as exc:
        raise StackingError(f"invalid residual table binding: {exc}") from exc

    def read(path: Path) -> list[tuple[str, float]]:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != ("row_hash", "residual"):
                raise StackingError("residual table has the wrong schema")
            return [(row["row_hash"], _number(row["residual"], "residual")) for row in reader]

    first_rows = read(first_path)
    second_rows = read(second_path)
    if [row[0] for row in first_rows] != [row[0] for row in second_rows] or len(first_rows) < 3:
        raise StackingError("residual rows are not identically aligned")
    x = [row[1] for row in first_rows]
    y = [row[1] for row in second_rows]
    x_mean = fsum(x) / len(x)
    y_mean = fsum(y) / len(y)
    x_centered = [value - x_mean for value in x]
    y_centered = [value - y_mean for value in y]
    denominator = sqrt(
        fsum(value * value for value in x_centered)
        * fsum(value * value for value in y_centered)
    )
    if denominator == 0.0:
        raise StackingError("residual correlation is undefined for zero variance")
    return fsum(
        first_value * second_value
        for first_value, second_value in zip(x_centered, y_centered, strict=True)
    ) / denominator


def _derive_stacking_evidence(
    shortlist_paths: Sequence[Path],
    residual_paths: Sequence[Path],
    gain_paths: Sequence[Path],
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    if len(shortlist_paths) != 2 or len(set(shortlist_paths)) != 2:
        raise StackingError("stacking evidence requires exactly two shortlists")
    if len(residual_paths) != 2 or len(set(residual_paths)) != 2:
        raise StackingError("stacking evidence requires exactly two residual bundles")
    if len(gain_paths) != 2 or len(set(gain_paths)) != 2:
        raise StackingError("stacking evidence requires exactly two gain bundles")
    shortlists = [verify_development_shortlist(path) for path in shortlist_paths]
    shortlist_by_task = {str(item["task_id"]): (path, item) for path, item in zip(shortlist_paths, shortlists, strict=True)}
    if set(shortlist_by_task) != {VARIANT_TASK, RNA_ATAC_TASK}:
        raise StackingError("stacking shortlists must cover variant and RNA-ATAC")
    gains = [verify_development_stack_gain_bundle(path) for path in gain_paths]
    gain_by_task = {str(item["task_id"]): (path, item) for path, item in zip(gain_paths, gains, strict=True)}
    if set(gain_by_task) != set(shortlist_by_task):
        raise StackingError("stack gain bundles must cover variant and RNA-ATAC")
    for task_id, (_, gain) in gain_by_task.items():
        shortlist_path, shortlist = shortlist_by_task[task_id]
        if gain["development_shortlist_binding"] != _shortlist_binding(shortlist_path, shortlist):
            raise StackingError(f"{task_id} gain binds a different shortlist")
    residuals = [verify_development_residual_bundle(path) for path in residual_paths]
    if residuals[0]["candidate_id"] == residuals[1]["candidate_id"]:
        raise StackingError("residual bundles must compare distinct candidates")
    if residuals[0]["source_family_id"] == residuals[1]["source_family_id"]:
        raise StackingError("residual bundles must compare distinct model families")
    residual_task = str(residuals[0]["task_id"])
    if residual_task not in shortlist_by_task or residuals[1]["task_id"] != residual_task:
        raise StackingError("residual bundles must belong to one comparison task")
    shortlist_path, shortlist = shortlist_by_task[residual_task]
    expected_binding = _shortlist_binding(shortlist_path, shortlist)
    if any(item["development_shortlist_binding"] != expected_binding for item in residuals):
        raise StackingError("residual bundles bind a different shortlist")
    residual_component_ids = {
        str(residuals[0]["candidate_id"]), str(residuals[1]["candidate_id"])
    }
    residual_gain = gain_by_task[residual_task][1]
    gain_component_ids = {
        str(residual_gain["sequence_candidate"]["candidate_id"]),
        str(residual_gain["context_candidate"]["candidate_id"]),
    }
    if residual_component_ids != gain_component_ids:
        raise StackingError(
            "residual complementarity and stack gain must use the same two candidates"
        )
    correlation = _residual_correlation(
        residual_paths[0], residuals[0], residual_paths[1], residuals[1]
    )
    variant_gain = gain_by_task[VARIANT_TASK][1]
    rna_gain = gain_by_task[RNA_ATAC_TASK][1]
    derived = {
        "residual_correlation": correlation,
        "relative_deviance_gains": {
            RNA_ATAC_TASK: float(rna_gain["observed_gain"])
        },
        "absolute_correlation_or_f1_gains": {
            VARIANT_TASK: float(variant_gain["observed_gain"])
        },
        "qualifying_seed_count": min(
            int(variant_gain["qualifying_seed_count"]),
            int(rna_gain["qualifying_seed_count"]),
        ),
        "evaluated_seed_count": len(STACK_SEEDS),
        "study_count": len(
            set(variant_gain["study_ids"]) | set(rna_gain["study_ids"])
        ),
        "cross_fitted": all(item["cross_fitted"] is True for item in gains),
        "nonnegative_stack": all(item["nonnegative_stack"] is True for item in gains),
        "best_open_models_compared": all(
            item["best_open_models_compared"] is True for item in gains
        ),
    }
    bindings = []
    for path in shortlist_paths:
        bindings.append(
            _source_binding(
                path,
                role="development_shortlist",
                artifact_class="development_shortlist",
                filename="development_shortlist.json",
            )
        )
    for path in residual_paths:
        bindings.append(
            _source_binding(
                path,
                role="development_residual_bundle",
                artifact_class="development_residual_bundle",
                filename="development_residual_bundle.json",
            )
        )
    for path in gain_paths:
        bindings.append(
            _source_binding(
                path,
                role="development_gain_bundle",
                artifact_class="development_stack_gain_bundle",
                filename="development_stack_gain_bundle.json",
            )
        )
    return derived, sorted(bindings, key=lambda item: (item["role"], item["path"]))


def freeze_stacking_evidence(
    *,
    development_shortlist_paths: Sequence[str | Path],
    development_residual_bundle_paths: Sequence[str | Path],
    development_stack_gain_bundle_paths: Sequence[str | Path],
    output_root: str | Path,
) -> Path:
    """Freeze the only evidence tree that may authorize the conditional model."""

    shortlist_paths = tuple(_safe_path(path, "development shortlist") for path in development_shortlist_paths)
    residual_paths = tuple(_safe_path(path, "development residual bundle") for path in development_residual_bundle_paths)
    gain_paths = tuple(_safe_path(path, "development stack gain bundle") for path in development_stack_gain_bundle_paths)
    derived, bindings = _derive_stacking_evidence(
        shortlist_paths, residual_paths, gain_paths
    )
    identity = {
        "schema_version": STACKING_EVIDENCE_SCHEMA_VERSION,
        "derived_evidence": derived,
        "source_bindings": bindings,
        "created_before_outcome_unblind": True,
        "sealed_results_used": False,
    }
    payload = {"evidence_id": canonical_sha256(identity), **identity}
    output = reject_symlink_components(Path(output_root), label="stacking evidence output root").resolve()
    output.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".stacking-evidence.", dir=output))
    try:
        write_json_exclusive(staging / "stacking_evidence.json", payload, mode=0o440)
        freeze_tree(
            staging,
            {"artifact_class": "stacking_evidence", "evidence_id": payload["evidence_id"]},
        )
        target = output / f"stacking-evidence--{payload['evidence_id']}"
        publish_directory_noreplace(staging, target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_stacking_evidence(target)
    return target


def verify_stacking_evidence(path: str | Path) -> dict[str, Any]:
    """Recursively rederive complementarity evidence from its six sources."""

    root = _safe_path(path, "stacking evidence")
    try:
        manifest = verify_frozen_tree(root)
        payload = json.loads((root / "stacking_evidence.json").read_text(encoding="utf-8"))
    except (ArtifactError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise StackingError(f"invalid stacking evidence: {exc}") from exc
    if not isinstance(payload, Mapping) or payload.get("schema_version") != STACKING_EVIDENCE_SCHEMA_VERSION:
        raise StackingError("unsupported stacking evidence schema")
    claimed = _hash(payload.get("evidence_id"), "evidence_id")
    identity = dict(payload)
    identity.pop("evidence_id", None)
    if canonical_sha256(identity) != claimed:
        raise StackingError("stacking evidence identity does not rederive")
    if manifest.get("metadata") != {
        "artifact_class": "stacking_evidence",
        "evidence_id": claimed,
    }:
        raise StackingError("stacking evidence manifest metadata mismatch")
    bindings = payload.get("source_bindings")
    if not isinstance(bindings, list) or len(bindings) != 6:
        raise StackingError("stacking evidence must bind exactly six source artifacts")
    if bindings != sorted(bindings, key=lambda item: (item["role"], item["path"])):
        raise StackingError("stacking evidence source bindings are not canonical")
    shortlist_paths = []
    residual_paths = []
    gain_paths = []
    contracts = {
        "development_shortlist": (
            shortlist_paths,
            "development_shortlist",
            "development_shortlist.json",
        ),
        "development_residual_bundle": (
            residual_paths,
            "development_residual_bundle",
            "development_residual_bundle.json",
        ),
        "development_gain_bundle": (
            gain_paths,
            "development_stack_gain_bundle",
            "development_stack_gain_bundle.json",
        ),
    }
    for binding in bindings:
        if not isinstance(binding, Mapping) or binding.get("role") not in contracts:
            raise StackingError("stacking evidence source binding has the wrong schema")
        destination, artifact_class, filename = contracts[binding["role"]]
        source = _safe_path(binding.get("path"), f"{binding['role']} source")
        expected = _source_binding(
            source,
            role=binding["role"],
            artifact_class=artifact_class,
            filename=filename,
        )
        if dict(binding) != expected:
            raise StackingError("stacking evidence source binding changed")
        destination.append(source)
    derived, expected_bindings = _derive_stacking_evidence(
        shortlist_paths, residual_paths, gain_paths
    )
    if payload.get("derived_evidence") != derived or bindings != expected_bindings:
        raise StackingError("stacking evidence does not rederive from frozen sources")
    evidence = ComplementarityEvidence.from_stacking_evidence(root)
    evidence.require_verified_artifact()
    return dict(payload)


__all__ = [
    "DEVELOPMENT_STACK_GAIN_BUNDLE_SCHEMA_VERSION",
    "STACK_COMPARISON_POLICY_ID",
    "STACK_FIT_POLICY_ID",
    "STACK_SEEDS",
    "StackingError",
    "freeze_development_stack_gain_bundle",
    "freeze_stacking_evidence",
    "verify_development_stack_gain_bundle",
    "verify_stacking_evidence",
]
