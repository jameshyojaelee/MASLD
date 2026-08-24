"""Resource-paper and project-sealed evaluation firewalls."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from enum import StrEnum
import json
from math import fsum, isfinite, sqrt
import os
from pathlib import Path
import shutil
import tempfile
import tomllib
from typing import Any, Iterable, Mapping

from .artifacts import (
    ArtifactError,
    canonical_hash,
    freeze_tree,
    publish_directory_noreplace,
    reject_symlink_components,
    sha256_file,
    verify_frozen_tree,
    write_json_exclusive,
)
from .contracts import (
    ArtifactRef,
    ContractError,
    DatasetManifest,
    DatasetRole,
    DatasetStatus,
    MissingState,
    PredictionBundle,
    SelectionLock,
)
from .evaluators.stats import (
    EmpiricalPowerGateResult,
    empirical_bootstrap_power_gate,
)
from .hashing import HashingError, canonicalize, is_sha256


class FirewallError(RuntimeError):
    """Raised when an operation crosses a declared scientific firewall."""


class CampaignPhase(StrEnum):
    PLANNING = "planning"
    DEVELOPMENT = "development"
    SELECTION_LOCK = "selection_lock"
    SEALED_INFERENCE = "sealed_inference"
    SEALED_EVALUATION = "sealed_evaluation"
    TERMINAL_REPORTING = "terminal_reporting"


class SealedArtifact(StrEnum):
    FEATURES = "features"
    LABELS = "labels"
    OUTCOMES = "outcomes"
    PREDICTIONS = "predictions"


DEFAULT_RESOURCE_AUTHORITIES = (
    "docs/README.md",
    "docs/PAPER.md",
    "docs/STATUS.md",
    "docs/RESULTS.md",
    "docs/ROADMAP.md",
    "Analysis/SingleCell/results_gpu_v2/hotspot_modules/all_modules.tsv",
    "Analysis/SingleCell/results_gpu_v2/hotspot_modules/module_names.tsv",
    "Analysis/SingleCell/results_gpu_v2/hotspot_modules/cholangiocytes/module_genes.tsv",
    "Analysis/SingleCell/results_gpu_v2/hotspot_modules/fibroblasts/module_genes.tsv",
    "Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes/module_genes.tsv",
    "Analysis/SingleCell/results_gpu_v2/hotspot_modules/macrophages/module_genes.tsv",
    "Analysis/SingleCell/results_gpu_v2/hotspot_modules/tcells/module_genes.tsv",
    "Analysis/Multimodal_Program_Projection/candidates/gene-catalog-v2-contract-2026-08-11",
    # figures/main/INDEX.md is the current six-figure routing authority and
    # names all six child directories as the manuscript figure surface.
    "figures/main",
    # The Next.js source tree is the current portal implementation authority.
    # Generated .next/out/node_modules trees are excluded. Public portal data,
    # including every per-gene JSON, is an authority and is fully inventoried.
    "masld-atlas-v2/src",
    "masld-atlas-v2/public/data",
    "masld-atlas-v2/package.json",
    "masld-atlas-v2/package-lock.json",
    "masld-atlas-v2/next.config.ts",
    "masld-atlas-v2/tsconfig.json",
    "masld-atlas-v2/eslint.config.mjs",
    # Exact Gene Catalog/portal contract producers. Runtime logs and Python
    # caches beneath scripts/portal are not protected authorities.
    "scripts/portal/gene_catalog_v2.py",
    "scripts/portal/gene_catalog_observability.py",
    "scripts/portal/emit_gene_catalog_v2_contract.py",
    "scripts/portal/validate_gene_catalog_v2_contract.py",
    # Cas13_Library_Design/README.md names these as the two current screen
    # authorities.  Protect the June-v9 build and the later gene-membership
    # freeze; no final order-ready guide authority exists yet.
    "Cas13_Library_Design/README.md",
    "Cas13_Library_Design/data/BUILD_MANIFEST.txt",
    "Cas13_Library_Design/data/guides/library.csv",
    "Cas13_Library_Design/scripts/rebuild_cas13_library.R",
    "Cas13_Library_Design/scripts/guides/guide_overlap_fix_20260713/frozen_cas13_gene_library_20260731_v1/run_v1",
)

PREDICTION_BUNDLE_SCHEMA_VERSION = "masld-bench-prediction-bundle-v1"
PREDICTION_COMMIT_SCHEMA_VERSION = "masld-bench-prediction-commit-v1"
_VARIANT_TASK_ID = "variant_to_regulation"
_VARIANT_SECONDARY_ENDPOINT_ID = "eqtl_retrieval"
_VARIANT_SECONDARY_EVALUATOR_ID = "variant_ld_block_eqtl_retrieval_auprc_v1"
_VARIANT_SECONDARY_SCORE_TRANSFORM_ID = "identity_link_score_v1"
DEVELOPMENT_POWER_EVIDENCE_SCHEMA_VERSION = (
    "masld-bench-development-power-evidence-v1"
)
BOOTSTRAP_DISTRIBUTION_SCHEMA_VERSION = (
    "masld-bench-endpoint-bootstrap-distribution-v1"
)
POWER_DECISION_SCHEMA_VERSION = "masld-bench-power-decision-v2"
SEALED_OUTCOME_BUNDLE_SCHEMA_VERSION = "masld-bench-sealed-outcome-bundle-v1"
OUTCOME_CONSUMPTION_SCHEMA_VERSION = "masld-bench-outcome-consumption-v1"
TERMINAL_OUTCOME_POLICY = "no_reselection_recalibration_threshold_change_or_repair"
EMPIRICAL_POWER_METHOD_ID = "empirical_paired_endpoint_bootstrap_v1"
PRODUCTION_POWER_RESAMPLES = 10_000
PRODUCTION_RESAMPLING_POLICY = "production_10000_endpoint_recomputations"
UNIT_TEST_RESAMPLING_POLICY = "unit_test_fixture_endpoint_recomputations"

_PROMOTION_GATE_REGISTRY = (
    Path(__file__).resolve().parents[2]
    / "config"
    / "evaluation"
    / "promotion_gates.toml"
)
_ENDPOINT_EVALUATOR_SOURCE = Path(__file__).resolve().parent / "evaluators" / "endpoint_power.py"
_TASK_ENDPOINT_EVALUATORS: Mapping[str, str] = {
    "bulk_state_transfer": "bulk_paired_spearman_gain_v1",
    "cell_state_mapping": "cell_donor_balanced_macro_f1_v1",
    "rna_conditioned_atac": "rna_atac_two_way_deviance_reduction_v1",
    "typed_evidence_graph": "graph_ld_block_auprc_gain_v1",
    "variant_to_regulation": "variant_ld_block_fisher_z_spearman_gain_v1",
}
_PREDICTION_TABLE_FIELDS: Mapping[str, tuple[str, ...]] = {
    "bulk_state_transfer": ("row_hash", "unit_hash", "predicted"),
    "cell_state_mapping": ("row_hash", "unit_hash", "predicted_class"),
    "rna_conditioned_atac": (
        "row_hash",
        "donor_hash",
        "block_hash",
        "stratum",
        "predicted",
    ),
    "typed_evidence_graph": (
        "row_hash",
        "unit_hash",
        "block_hash",
        "predicted",
    ),
    "variant_to_regulation": (
        "row_hash",
        "unit_hash",
        "block_hash",
        "stratum",
        "predicted",
    ),
}
_OUTCOME_TABLE_FIELDS: Mapping[str, tuple[str, ...]] = {
    "bulk_state_transfer": ("row_hash", "unit_hash", "observed"),
    "cell_state_mapping": ("row_hash", "unit_hash", "observed_class"),
    "rna_conditioned_atac": (
        "row_hash",
        "donor_hash",
        "block_hash",
        "stratum",
        "observed",
    ),
    "typed_evidence_graph": (
        "row_hash",
        "unit_hash",
        "block_hash",
        "observed_binary",
    ),
    "variant_to_regulation": (
        "row_hash",
        "unit_hash",
        "block_hash",
        "stratum",
        "observed",
    ),
}


def _strict_mapping(
    value: object,
    *,
    fields: frozenset[str],
    label: str,
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise FirewallError(f"{label} must be an object")
    if any(not isinstance(key, str) for key in value):
        raise FirewallError(f"{label} keys must be strings")
    unknown = sorted(set(value).difference(fields))
    missing = sorted(fields.difference(value))
    if unknown:
        raise FirewallError(f"{label} has unknown fields: {', '.join(unknown)}")
    if missing:
        raise FirewallError(f"{label} is missing fields: {', '.join(missing)}")
    return dict(value)


def _read_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise FirewallError(f"cannot read {label}: {error}") from error
    if not isinstance(payload, Mapping):
        raise FirewallError(f"{label} must contain a JSON object")
    return dict(payload)


def _configured_path(value: str | Path, label: str) -> Path:
    try:
        return reject_symlink_components(Path(value), label=label)
    except ArtifactError as error:
        raise FirewallError(str(error)) from error


def _nonempty_string(value: object, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip() or "\x00" in value:
        raise FirewallError(f"{field_name} must be a non-empty string")
    return value


def _sha256(value: object, field_name: str) -> str:
    if not is_sha256(value):
        raise FirewallError(
            f"{field_name} must be a 64-character lowercase hexadecimal SHA-256"
        )
    return str(value)


def _canonical_mapping(value: object, field_name: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise FirewallError(f"{field_name} must be an object")
    try:
        normalized = canonicalize(value)
    except (HashingError, TypeError, ValueError) as error:
        raise FirewallError(f"{field_name} is not canonical JSON: {error}") from error
    if not isinstance(normalized, dict):
        raise FirewallError(f"{field_name} must be an object")
    return normalized


def _string_tuple(
    value: object,
    field_name: str,
    *,
    allow_empty: bool = False,
    sorted_required: bool = True,
) -> tuple[str, ...]:
    if not isinstance(value, (list, tuple)):
        raise FirewallError(f"{field_name} must be an array")
    result = tuple(
        _nonempty_string(item, f"{field_name}[{index}]")
        for index, item in enumerate(value)
    )
    if not allow_empty and not result:
        raise FirewallError(f"{field_name} must not be empty")
    if len(set(result)) != len(result):
        raise FirewallError(f"{field_name} must not contain duplicates")
    if sorted_required and list(result) != sorted(result):
        raise FirewallError(f"{field_name} must be canonically sorted")
    return result


def _sha256_tuple(
    value: object, field_name: str, *, allow_empty: bool = False
) -> tuple[str, ...]:
    strings = _string_tuple(value, field_name, allow_empty=allow_empty)
    return tuple(
        _sha256(item, f"{field_name}[{index}]")
        for index, item in enumerate(strings)
    )


def _verify_selection_dir(path: Path) -> SelectionLock:
    # Imported lazily because planner imports resource_snapshot from this module.
    from .selection import SelectionError, verify_selection_lock

    try:
        return verify_selection_lock(path)
    except SelectionError as error:
        raise FirewallError(f"invalid frozen SelectionLock: {error}") from error


def assert_dataset_access(
    *,
    dataset: DatasetManifest | Mapping[str, Any],
    phase: CampaignPhase | str,
    artifact: SealedArtifact | str,
) -> None:
    """Enforce access from a validated manifest, never a caller-supplied role."""

    try:
        manifest = (
            dataset
            if isinstance(dataset, DatasetManifest)
            else DatasetManifest.from_dict(dataset)
        )
    except (ContractError, TypeError, ValueError) as error:
        raise FirewallError(f"invalid DatasetManifest: {error}") from error
    try:
        normalized_phase = CampaignPhase(phase)
        normalized_artifact = SealedArtifact(artifact)
    except ValueError as error:
        raise FirewallError(f"invalid dataset-access request: {error}") from error

    if manifest.status in {DatasetStatus.BLOCKED, DatasetStatus.DEFERRED} or manifest.role in {
        DatasetRole.BLOCKED,
        DatasetRole.DEFERRED,
    }:
        raise FirewallError(
            f"dataset {manifest.dataset_id} is {manifest.status.value} and may not be accessed"
        )
    is_withheld = (
        manifest.role is DatasetRole.WITHHELD_SEALED
        and manifest.status is DatasetStatus.WITHHELD_SEALED
    )
    if not is_withheld:
        return
    permitted = {
        CampaignPhase.SEALED_INFERENCE: {
            SealedArtifact.FEATURES,
            SealedArtifact.PREDICTIONS,
        },
        CampaignPhase.SEALED_EVALUATION: {
            SealedArtifact.LABELS,
            SealedArtifact.OUTCOMES,
            SealedArtifact.PREDICTIONS,
        },
        CampaignPhase.TERMINAL_REPORTING: {SealedArtifact.PREDICTIONS},
    }
    if normalized_artifact not in permitted.get(normalized_phase, set()):
        raise FirewallError(
            f"{normalized_artifact.value} from {manifest.dataset_id} may not be accessed "
            f"during {normalized_phase.value}"
        )


def resource_snapshot(
    repository_root: str | Path,
    paths: Iterable[str] = DEFAULT_RESOURCE_AUTHORITIES,
) -> dict[str, Any]:
    root = _configured_path(repository_root, "Resource repository root").resolve()
    records: list[dict[str, Any]] = []
    authority_roots = sorted(set(str(item) for item in paths))
    if not authority_roots:
        raise FirewallError("Resource firewall requires at least one authority root")

    def record(path: Path) -> dict[str, Any]:
        if path.is_symlink() or not path.is_file():
            raise FirewallError(
                f"protected Resource authority member is not a regular file: {path}"
            )
        before = path.stat(follow_symlinks=False)
        digest = sha256_file(path)
        after = path.stat(follow_symlinks=False)
        if (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
        ) != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
        ):
            raise FirewallError(
                f"protected Resource authority changed while hashing: {path}"
            )
        return {
            "path": path.relative_to(root).as_posix(),
            "sha256": digest,
            "size_bytes": after.st_size,
        }

    for relative in authority_roots:
        configured = _configured_path(
            root / relative, f"protected Resource authority {relative}"
        )
        path = configured.resolve()
        try:
            path.relative_to(root)
        except ValueError as error:
            raise FirewallError(f"protected path escapes repository: {relative}") from error
        if path.is_file():
            records.append(record(path))
            continue
        if path.is_dir():
            pending = [path]
            while pending:
                directory = pending.pop()
                for member in sorted(directory.iterdir(), key=lambda item: item.name):
                    if member.is_symlink():
                        raise FirewallError(
                            f"protected Resource authority may not contain symlinks: {member}"
                        )
                    if member.name == "__pycache__" or member.suffix == ".pyc":
                        continue
                    if member.is_dir():
                        pending.append(member)
                    elif member.is_file():
                        records.append(record(member))
                    else:
                        raise FirewallError(
                            "protected Resource authority contains a non-regular member: "
                            f"{member}"
                        )
            continue
        raise FirewallError(f"protected Resource authority is missing: {relative}")
    records.sort(key=lambda item: str(item["path"]))
    paths_seen = [str(item["path"]) for item in records]
    if len(set(paths_seen)) != len(paths_seen):
        raise FirewallError("Resource authority roots overlap and duplicate protected files")
    identity = {"authority_roots": authority_roots, "files": records}
    return {
        "schema_version": "masld-resource-firewall-v2",
        "authority_roots": authority_roots,
        "files": records,
        "snapshot_sha256": canonical_hash(identity),
    }


def assert_resource_unchanged(
    baseline: Mapping[str, Any],
    *,
    repository_root: str | Path,
) -> dict[str, Any]:
    if baseline.get("schema_version") != "masld-resource-firewall-v2":
        raise FirewallError("unsupported or membership-blind Resource firewall schema")
    roots = baseline.get("authority_roots")
    files = baseline.get("files")
    if not isinstance(roots, list) or not isinstance(files, list):
        raise FirewallError("Resource firewall baseline is missing roots or files")
    identity = {"authority_roots": roots, "files": files}
    if canonical_hash(identity) != baseline.get("snapshot_sha256"):
        raise FirewallError("Resource firewall baseline hash is invalid")
    observed = resource_snapshot(repository_root, (str(item) for item in roots))
    if observed["snapshot_sha256"] != baseline.get("snapshot_sha256"):
        expected = {item["path"]: item for item in baseline.get("files", [])}
        actual = {item["path"]: item for item in observed["files"]}
        changed = sorted(
            path for path in set(expected) | set(actual) if expected.get(path) != actual.get(path)
        )
        raise FirewallError(
            "current MASLD Resource authorities changed during model work: "
            + ", ".join(changed)
        )
    return observed


def _freeze_state_document(
    *,
    state_root: Path,
    collection: str,
    identifier: str,
    filename: str,
    payload: Mapping[str, Any],
    artifact_class: str,
) -> Path:
    configured_state = _configured_path(state_root, "evaluator-state root")
    collection_root = configured_state.resolve() / collection
    collection_root.mkdir(parents=True, exist_ok=True)
    target = collection_root / identifier
    if target.exists():
        raise FirewallError(f"immutable {artifact_class} already exists: {target}")
    staging = Path(tempfile.mkdtemp(prefix=f".{artifact_class}.", dir=collection_root))
    try:
        write_json_exclusive(staging / filename, payload, mode=0o440)
        freeze_tree(
            staging,
            {"artifact_class": artifact_class, "artifact_id": identifier},
        )
        try:
            publish_directory_noreplace(staging, target)
        except ArtifactError as error:
            raise FirewallError(
                f"could not publish immutable {artifact_class}: {error}"
            ) from error
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return target


def _verify_state_manifest(
    root: Path, *, artifact_class: str, artifact_id: str
) -> None:
    try:
        manifest = verify_frozen_tree(root)
    except (ArtifactError, OSError, ValueError, KeyError, TypeError) as error:
        raise FirewallError(f"invalid frozen {artifact_class}: {error}") from error
    metadata = manifest.get("metadata")
    if not isinstance(metadata, Mapping):
        raise FirewallError(f"{artifact_class} manifest metadata is missing")
    if metadata.get("artifact_class") != artifact_class:
        raise FirewallError(f"artifact tree is not a {artifact_class}")
    if metadata.get("artifact_id") != artifact_id:
        raise FirewallError(f"{artifact_class} manifest identity mismatch")


def _task_decision(lock: SelectionLock, task_id: str) -> Mapping[str, Any]:
    matches = [
        decision
        for decision in lock.task_decisions
        if decision.get("task_id") == task_id
    ]
    if len(matches) != 1:
        raise FirewallError(
            f"SelectionLock does not contain exactly one decision for task {task_id}"
        )
    return matches[0]


def _prediction_binding(
    bundle: PredictionBundle, lock: SelectionLock
) -> tuple[Mapping[str, Any], str, str]:
    if bundle.schema_version != PREDICTION_BUNDLE_SCHEMA_VERSION:
        raise FirewallError(
            f"unsupported PredictionBundle schema_version: {bundle.schema_version}"
        )
    if bundle.n_predictions < 1:
        raise FirewallError("PredictionBundle must contain at least one prediction")
    if len({artifact.path for artifact in bundle.artifacts}) != len(bundle.artifacts):
        raise FirewallError("PredictionBundle artifact paths must be unique")
    decision = _task_decision(lock, bundle.task_id)
    selected_run_ids = {
        str(record["run_id"]) for record in decision.get("selected_runs", ())
    }
    baseline_run_ids = {
        str(record["run_id"]) for record in decision.get("baseline_runs", ())
    }
    is_selected = (
        bundle.model_id == decision.get("selected_model_id")
        and bundle.run_id in selected_run_ids
    )
    is_baseline = (
        bundle.model_id == decision.get("baseline_model_id")
        and bundle.run_id in baseline_run_ids
    )
    secondary_record: Mapping[str, Any] | None = None
    secondary = decision.get("variant_secondary_evaluation")
    if bundle.task_id == _VARIANT_TASK_ID:
        if not isinstance(secondary, Mapping):
            raise FirewallError(
                "variant PredictionBundles require a locked secondary evaluation"
            )
        raw_comparators = secondary.get("comparators")
        if not isinstance(raw_comparators, (list, tuple)):
            raise FirewallError("variant secondary comparator roster is invalid")
        matching_secondary = [
            record
            for record in raw_comparators
            if isinstance(record, Mapping)
            and record.get("run_id") == bundle.run_id
            and record.get("model_id") == bundle.model_id
        ]
        if len(matching_secondary) > 1:
            raise FirewallError("variant secondary comparator roster repeats a run/model")
        if matching_secondary:
            secondary_record = matching_secondary[0]
    elif secondary is not None:
        raise FirewallError(
            "nonvariant PredictionBundle carries a variant secondary evaluation"
        )
    is_secondary = secondary_record is not None
    if is_secondary and (is_selected or is_baseline):
        raise FirewallError("PredictionBundle has ambiguous locked model roles")
    if not is_selected and not is_baseline and not is_secondary:
        raise FirewallError(
            "PredictionBundle run/model is neither selected, baseline, nor a locked "
            f"secondary comparator for {bundle.task_id}"
        )
    model_role = (
        "selected_and_baseline"
        if is_selected and is_baseline
        else "selected"
        if is_selected
        else "baseline"
        if is_baseline
        else "secondary_comparator"
    )
    if is_secondary:
        if bundle.missing_state is not MissingState.OBSERVED:
            raise FirewallError(
                "secondary comparator PredictionBundle must contain observed scores"
            )
        expected_metadata = {
            "variant_endpoint_id": _VARIANT_SECONDARY_ENDPOINT_ID,
            "variant_evaluator_id": _VARIANT_SECONDARY_EVALUATOR_ID,
            "variant_score_transform_id": _VARIANT_SECONDARY_SCORE_TRANSFORM_ID,
            "primary_endpoint_scoring_allowed": False,
            "variant_secondary_comparator_sha256": canonical_hash(secondary_record),
        }
        if any(bundle.metadata.get(key) != value for key, value in expected_metadata.items()):
            raise FirewallError(
                "secondary comparator PredictionBundle endpoint metadata differs from "
                "its SelectionLock binding"
            )
    sealed_dataset_ids = tuple(decision.get("sealed_dataset_ids", ()))
    if not sealed_dataset_ids:
        raise FirewallError(f"task {bundle.task_id} has no locked sealed dataset")
    if tuple(sorted(bundle.dataset_ids)) != sealed_dataset_ids:
        raise FirewallError(
            f"PredictionBundle datasets do not equal the sealed binding for {bundle.task_id}"
        )
    registry_sha256s = decision.get("dataset_registry_sha256s")
    if not isinstance(registry_sha256s, Mapping):
        raise FirewallError(f"SelectionLock dataset bindings are missing for {bundle.task_id}")
    binding = {
        "task_id": bundle.task_id,
        "dataset_ids": list(sealed_dataset_ids),
        "dataset_registry_sha256s": {
            dataset_id: registry_sha256s[dataset_id]
            for dataset_id in sealed_dataset_ids
        },
    }
    return decision, canonical_hash(binding), model_role


@dataclass(frozen=True, slots=True)
class PredictionCommit:
    selection_lock_id: str
    selection_lock_manifest_sha256: str
    selection_lock_dir: str
    prediction_bundle_sha256: str
    prediction_bundle_path: str
    bundle_id: str
    run_id: str
    task_id: str
    model_id: str
    model_role: str
    dataset_ids: tuple[str, ...]
    dataset_bindings_sha256: str
    split_id: str
    row_id_field: str
    unit_id_field: str
    unit_id_namespace: str
    biological_unit: str
    table_schema_sha256: str
    source_join_key_sha256: str
    n_predictions: int
    prediction_artifacts: tuple[ArtifactRef, ...]
    standardized_table: ArtifactRef
    row_ids: ArtifactRef
    prediction_commit_sha256: str

    _FIELDS = frozenset(
        {
            "schema_version",
            "selection_lock_id",
            "selection_lock_manifest_sha256",
            "selection_lock_dir",
            "prediction_bundle_sha256",
            "prediction_bundle_path",
            "bundle_id",
            "run_id",
            "task_id",
            "model_id",
            "model_role",
            "dataset_ids",
            "dataset_bindings_sha256",
            "split_id",
            "row_id_field",
            "unit_id_field",
            "unit_id_namespace",
            "biological_unit",
            "table_schema_sha256",
            "source_join_key_sha256",
            "n_predictions",
            "prediction_artifacts",
            "standardized_table",
            "row_ids",
            "prediction_commit_sha256",
        }
    )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PredictionCommit":
        raw = _strict_mapping(value, fields=cls._FIELDS, label="PredictionCommit")
        if raw["schema_version"] != PREDICTION_COMMIT_SCHEMA_VERSION:
            raise FirewallError(
                f"unsupported PredictionCommit schema_version: {raw['schema_version']}"
            )
        artifacts = raw["prediction_artifacts"]
        if not isinstance(artifacts, (list, tuple)) or not artifacts:
            raise FirewallError("PredictionCommit.prediction_artifacts must not be empty")
        try:
            parsed_artifacts = tuple(ArtifactRef.from_dict(item) for item in artifacts)
            standardized_table = ArtifactRef.from_dict(raw["standardized_table"])
            row_ids = ArtifactRef.from_dict(raw["row_ids"])
        except ContractError as error:
            raise FirewallError(f"invalid PredictionCommit artifact: {error}") from error
        if len({artifact.path for artifact in parsed_artifacts}) != len(parsed_artifacts):
            raise FirewallError("PredictionCommit artifact paths must be unique")
        if (
            sum(artifact == standardized_table for artifact in parsed_artifacts) != 1
            or sum(artifact == row_ids for artifact in parsed_artifacts) != 1
        ):
            raise FirewallError(
                "PredictionCommit standardized table and row IDs must be bound artifacts"
            )
        selection_dir = _nonempty_string(
            raw["selection_lock_dir"], "PredictionCommit.selection_lock_dir"
        )
        bundle_path = _nonempty_string(
            raw["prediction_bundle_path"], "PredictionCommit.prediction_bundle_path"
        )
        if not Path(selection_dir).is_absolute() or not Path(bundle_path).is_absolute():
            raise FirewallError("PredictionCommit source paths must be absolute")
        model_role = _nonempty_string(
            raw["model_role"], "PredictionCommit.model_role"
        )
        if model_role not in {
            "selected",
            "baseline",
            "selected_and_baseline",
            "secondary_comparator",
        }:
            raise FirewallError("PredictionCommit.model_role is invalid")
        n_predictions = raw["n_predictions"]
        if (
            isinstance(n_predictions, bool)
            or not isinstance(n_predictions, int)
            or n_predictions < 1
        ):
            raise FirewallError("PredictionCommit.n_predictions must be a positive integer")
        return cls(
            selection_lock_id=_sha256(
                raw["selection_lock_id"], "PredictionCommit.selection_lock_id"
            ),
            selection_lock_manifest_sha256=_sha256(
                raw["selection_lock_manifest_sha256"],
                "PredictionCommit.selection_lock_manifest_sha256",
            ),
            selection_lock_dir=selection_dir,
            prediction_bundle_sha256=_sha256(
                raw["prediction_bundle_sha256"],
                "PredictionCommit.prediction_bundle_sha256",
            ),
            prediction_bundle_path=bundle_path,
            bundle_id=_nonempty_string(raw["bundle_id"], "PredictionCommit.bundle_id"),
            run_id=_sha256(raw["run_id"], "PredictionCommit.run_id"),
            task_id=_nonempty_string(raw["task_id"], "PredictionCommit.task_id"),
            model_id=_nonempty_string(raw["model_id"], "PredictionCommit.model_id"),
            model_role=model_role,
            dataset_ids=_string_tuple(
                raw["dataset_ids"], "PredictionCommit.dataset_ids"
            ),
            dataset_bindings_sha256=_sha256(
                raw["dataset_bindings_sha256"],
                "PredictionCommit.dataset_bindings_sha256",
            ),
            split_id=_nonempty_string(raw["split_id"], "PredictionCommit.split_id"),
            row_id_field=_nonempty_string(
                raw["row_id_field"], "PredictionCommit.row_id_field"
            ),
            unit_id_field=_nonempty_string(
                raw["unit_id_field"], "PredictionCommit.unit_id_field"
            ),
            unit_id_namespace=_nonempty_string(
                raw["unit_id_namespace"], "PredictionCommit.unit_id_namespace"
            ),
            biological_unit=_nonempty_string(
                raw["biological_unit"], "PredictionCommit.biological_unit"
            ),
            table_schema_sha256=_sha256(
                raw["table_schema_sha256"], "PredictionCommit.table_schema_sha256"
            ),
            source_join_key_sha256=_sha256(
                raw["source_join_key_sha256"],
                "PredictionCommit.source_join_key_sha256",
            ),
            n_predictions=n_predictions,
            prediction_artifacts=parsed_artifacts,
            standardized_table=standardized_table,
            row_ids=row_ids,
            prediction_commit_sha256=_sha256(
                raw["prediction_commit_sha256"],
                "PredictionCommit.prediction_commit_sha256",
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": PREDICTION_COMMIT_SCHEMA_VERSION,
            "selection_lock_id": self.selection_lock_id,
            "selection_lock_manifest_sha256": self.selection_lock_manifest_sha256,
            "selection_lock_dir": self.selection_lock_dir,
            "prediction_bundle_sha256": self.prediction_bundle_sha256,
            "prediction_bundle_path": self.prediction_bundle_path,
            "bundle_id": self.bundle_id,
            "run_id": self.run_id,
            "task_id": self.task_id,
            "model_id": self.model_id,
            "model_role": self.model_role,
            "dataset_ids": list(self.dataset_ids),
            "dataset_bindings_sha256": self.dataset_bindings_sha256,
            "split_id": self.split_id,
            "row_id_field": self.row_id_field,
            "unit_id_field": self.unit_id_field,
            "unit_id_namespace": self.unit_id_namespace,
            "biological_unit": self.biological_unit,
            "table_schema_sha256": self.table_schema_sha256,
            "source_join_key_sha256": self.source_join_key_sha256,
            "n_predictions": self.n_predictions,
            "prediction_artifacts": [item.to_dict() for item in self.prediction_artifacts],
            "standardized_table": self.standardized_table.to_dict(),
            "row_ids": self.row_ids.to_dict(),
            "prediction_commit_sha256": self.prediction_commit_sha256,
        }


def _prediction_commit_identity(commit: PredictionCommit) -> dict[str, Any]:
    identity = commit.to_dict()
    identity.pop("prediction_commit_sha256")
    return identity


def commit_predictions(
    *,
    selection_lock_dir: str | Path,
    prediction_bundle_path: str | Path,
    evaluator_state_dir: str | Path,
) -> PredictionCommit:
    """Commit a strict PredictionBundle and all referenced artifacts."""

    selection_dir = _configured_path(selection_lock_dir, "SelectionLock").resolve()
    lock = _verify_selection_dir(selection_dir)
    bundle_path = _configured_path(
        prediction_bundle_path, "PredictionBundle"
    ).resolve()
    try:
        bundle = PredictionBundle.load_json(bundle_path)
        bundle.validate_artifacts(bundle_path.parent)
        bundle_sha256 = sha256_file(bundle_path)
    except (ContractError, OSError, ValueError) as error:
        raise FirewallError(f"invalid PredictionBundle: {error}") from error
    _, dataset_bindings_sha256, model_role = _prediction_binding(bundle, lock)
    identity = {
        "schema_version": PREDICTION_COMMIT_SCHEMA_VERSION,
        "selection_lock_id": lock.lock_id,
        "selection_lock_manifest_sha256": sha256_file(
            selection_dir / "ARTIFACTS.json"
        ),
        "selection_lock_dir": selection_dir.as_posix(),
        "prediction_bundle_sha256": bundle_sha256,
        "prediction_bundle_path": bundle_path.as_posix(),
        "bundle_id": bundle.bundle_id,
        "run_id": bundle.run_id,
        "task_id": bundle.task_id,
        "model_id": bundle.model_id,
        "model_role": model_role,
        "dataset_ids": sorted(bundle.dataset_ids),
        "dataset_bindings_sha256": dataset_bindings_sha256,
        "split_id": bundle.split_id,
        "row_id_field": bundle.row_id_field,
        "unit_id_field": bundle.unit_id_field,
        "unit_id_namespace": bundle.unit_id_namespace,
        "biological_unit": bundle.biological_unit,
        "table_schema_sha256": bundle.table_schema_sha256,
        "source_join_key_sha256": bundle.source_join_key_sha256,
        "n_predictions": bundle.n_predictions,
        "prediction_artifacts": [artifact.to_dict() for artifact in bundle.artifacts],
        "standardized_table": bundle.standardized_table.to_dict(),
        "row_ids": bundle.row_ids.to_dict(),
    }
    payload = {
        **identity,
        "prediction_commit_sha256": canonical_hash(identity),
    }
    commit = PredictionCommit.from_dict(payload)
    commit_dir = _freeze_state_document(
        state_root=_configured_path(evaluator_state_dir, "evaluator-state root"),
        collection="prediction_commits",
        identifier=commit.prediction_commit_sha256,
        filename="prediction_commit.json",
        payload=commit.to_dict(),
        artifact_class="prediction_commit",
    )
    return verify_prediction_commit(commit_dir, reverify_sources=True)


def verify_prediction_commit(
    path: str | Path, *, reverify_sources: bool = True
) -> PredictionCommit:
    root = _configured_path(path, "PredictionCommit").resolve()
    try:
        verify_frozen_tree(root)
    except (ArtifactError, OSError, ValueError, KeyError, TypeError) as error:
        raise FirewallError(f"invalid frozen prediction_commit: {error}") from error
    payload = _read_json_object(root / "prediction_commit.json", "PredictionCommit")
    commit = PredictionCommit.from_dict(payload)
    _verify_state_manifest(
        root,
        artifact_class="prediction_commit",
        artifact_id=commit.prediction_commit_sha256,
    )
    if canonical_hash(_prediction_commit_identity(commit)) != commit.prediction_commit_sha256:
        raise FirewallError("PredictionCommit identity hash mismatch")
    if not reverify_sources:
        return commit

    try:
        selection_dir = Path(commit.selection_lock_dir)
        lock = _verify_selection_dir(selection_dir)
        if lock.lock_id != commit.selection_lock_id:
            raise FirewallError("PredictionCommit selection-lock identity mismatch")
        if (
            sha256_file(selection_dir / "ARTIFACTS.json")
            != commit.selection_lock_manifest_sha256
        ):
            raise FirewallError("PredictionCommit selection-lock manifest changed")
        bundle_path = Path(commit.prediction_bundle_path)
        if sha256_file(bundle_path) != commit.prediction_bundle_sha256:
            raise FirewallError("PredictionBundle changed after prediction commit")
        bundle = PredictionBundle.load_json(bundle_path)
        bundle.validate_artifacts(bundle_path.parent)
    except FirewallError:
        raise
    except (ContractError, OSError, ValueError) as error:
        raise FirewallError(f"PredictionBundle failed reverification: {error}") from error
    _, binding_sha256, model_role = _prediction_binding(bundle, lock)
    expected = {
        "bundle_id": bundle.bundle_id,
        "run_id": bundle.run_id,
        "task_id": bundle.task_id,
        "model_id": bundle.model_id,
        "model_role": model_role,
        "dataset_ids": tuple(sorted(bundle.dataset_ids)),
        "dataset_bindings_sha256": binding_sha256,
        "split_id": bundle.split_id,
        "row_id_field": bundle.row_id_field,
        "unit_id_field": bundle.unit_id_field,
        "unit_id_namespace": bundle.unit_id_namespace,
        "biological_unit": bundle.biological_unit,
        "table_schema_sha256": bundle.table_schema_sha256,
        "source_join_key_sha256": bundle.source_join_key_sha256,
        "n_predictions": bundle.n_predictions,
        "prediction_artifacts": tuple(bundle.artifacts),
        "standardized_table": bundle.standardized_table,
        "row_ids": bundle.row_ids,
    }
    observed = {
        "bundle_id": commit.bundle_id,
        "run_id": commit.run_id,
        "task_id": commit.task_id,
        "model_id": commit.model_id,
        "model_role": commit.model_role,
        "dataset_ids": commit.dataset_ids,
        "dataset_bindings_sha256": commit.dataset_bindings_sha256,
        "split_id": commit.split_id,
        "row_id_field": commit.row_id_field,
        "unit_id_field": commit.unit_id_field,
        "unit_id_namespace": commit.unit_id_namespace,
        "biological_unit": commit.biological_unit,
        "table_schema_sha256": commit.table_schema_sha256,
        "source_join_key_sha256": commit.source_join_key_sha256,
        "n_predictions": commit.n_predictions,
        "prediction_artifacts": commit.prediction_artifacts,
        "standardized_table": commit.standardized_table,
        "row_ids": commit.row_ids,
    }
    if observed != expected:
        raise FirewallError("PredictionCommit does not match its strict PredictionBundle")
    return commit


_EMPIRICAL_POWER_RESULT_FIELDS = frozenset(
    {
        "method_id",
        "passed",
        "n_units",
        "min_units",
        "observed_effect",
        "minimum_effect",
        "empirical_standard_error",
        "critical_value",
        "estimated_power",
        "power_target",
        "alpha",
        "effective_alpha",
        "n_primary_claims",
        "two_sided",
        "n_resamples",
        "reason",
    }
)


def _validated_empirical_power_result(
    value: EmpiricalPowerGateResult | Mapping[str, Any],
) -> dict[str, Any]:
    raw_value = value.to_dict() if isinstance(value, EmpiricalPowerGateResult) else value
    raw = _strict_mapping(
        raw_value,
        fields=_EMPIRICAL_POWER_RESULT_FIELDS,
        label="EmpiricalPowerGateResult",
    )
    if raw["method_id"] != EMPIRICAL_POWER_METHOD_ID:
        raise FirewallError(
            "EmpiricalPowerGateResult.method_id is not the production method"
        )
    for field_name in ("passed", "two_sided"):
        if not isinstance(raw[field_name], bool):
            raise FirewallError(
                f"EmpiricalPowerGateResult.{field_name} must be boolean"
            )
    for field_name in (
        "n_units",
        "min_units",
        "n_primary_claims",
        "n_resamples",
    ):
        value_at_field = raw[field_name]
        if (
            isinstance(value_at_field, bool)
            or not isinstance(value_at_field, int)
            or value_at_field < 1
        ):
            raise FirewallError(
                f"EmpiricalPowerGateResult.{field_name} must be a positive integer"
            )
    for field_name in (
        "observed_effect",
        "minimum_effect",
        "empirical_standard_error",
        "critical_value",
        "estimated_power",
        "power_target",
        "alpha",
        "effective_alpha",
    ):
        value_at_field = raw[field_name]
        if (
            isinstance(value_at_field, bool)
            or not isinstance(value_at_field, (int, float))
            or not isfinite(float(value_at_field))
        ):
            raise FirewallError(
                f"EmpiricalPowerGateResult.{field_name} must be finite numeric"
            )
    _nonempty_string(raw["reason"], "EmpiricalPowerGateResult.reason")
    return _canonical_mapping(raw, "EmpiricalPowerGateResult")


def _positive_integer_mapping(value: object, field_name: str) -> dict[str, int]:
    raw = _canonical_mapping(value, field_name)
    if not raw:
        raise FirewallError(f"{field_name} must not be empty")
    result: dict[str, int] = {}
    for key, item in raw.items():
        _nonempty_string(key, f"{field_name} key")
        if isinstance(item, bool) or not isinstance(item, int) or item < 1:
            raise FirewallError(f"{field_name}[{key!r}] must be a positive integer")
        result[key] = item
    return result


def _endpoint_evaluator_source_sha256() -> str:
    sources = {
        "endpoint_power.py": sha256_file(_ENDPOINT_EVALUATOR_SOURCE),
        "metrics.py": sha256_file(_ENDPOINT_EVALUATOR_SOURCE.with_name("metrics.py")),
    }
    return canonical_hash(sources)


def _power_method_source_sha256() -> str:
    return sha256_file(_ENDPOINT_EVALUATOR_SOURCE.with_name("stats.py"))


def _load_task_power_contract(
    lock: SelectionLock, task_id: str
) -> tuple[Mapping[str, Any], Mapping[str, Any], str]:
    try:
        registry_bytes = _PROMOTION_GATE_REGISTRY.read_bytes()
        registry = tomllib.loads(registry_bytes.decode("utf-8"))
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise FirewallError(f"cannot load promotion-gate authority: {error}") from error
    if registry.get("schema_version") != "masld-bench-promotion-gates-v1":
        raise FirewallError("unsupported promotion-gate authority schema")
    tasks = registry.get("tasks")
    if not isinstance(tasks, Mapping) or not isinstance(tasks.get(task_id), Mapping):
        raise FirewallError(f"promotion-gate authority has no task {task_id}")
    policy = dict(tasks[task_id])
    if policy.get("power_method_id") != EMPIRICAL_POWER_METHOD_ID:
        raise FirewallError(
            f"task {task_id} is not authorized for the production empirical power method"
        )
    expected_evaluator = _TASK_ENDPOINT_EVALUATORS.get(task_id)
    if expected_evaluator is None:
        raise FirewallError(
            f"task {task_id} has no registered endpoint-recomputed power evaluator"
        )

    decision = _task_decision(lock, task_id)
    if decision.get("promotion_gate") != policy.get("promotion_gate_id"):
        raise FirewallError(
            f"SelectionLock promotion gate for {task_id} differs from the authority"
        )
    thresholds = decision.get("thresholds")
    if not isinstance(thresholds, Mapping):
        raise FirewallError(f"SelectionLock thresholds are missing for {task_id}")
    promotion = thresholds.get("promotion")
    if not isinstance(promotion, Mapping):
        raise FirewallError(f"SelectionLock promotion thresholds are missing for {task_id}")
    raw_power = _strict_mapping(
        promotion.get("power"),
        fields=frozenset(
            {
                "metric_id",
                "minimum_effect",
                "alpha",
                "target_power",
                "n_primary_claims",
                "min_units",
                "two_sided",
            }
        ),
        label=f"{task_id}.promotion.power",
    )
    if raw_power["metric_id"] != policy.get("power_metric"):
        raise FirewallError(f"SelectionLock power metric changed for {task_id}")
    for field_name, expected in (
        ("minimum_effect", policy.get("power_minimum_effect")),
        ("alpha", 0.05),
        ("target_power", 0.80),
    ):
        observed = raw_power[field_name]
        if (
            isinstance(observed, bool)
            or not isinstance(observed, (int, float))
            or not isfinite(float(observed))
            or float(observed) != float(expected)
        ):
            raise FirewallError(
                f"SelectionLock power {field_name} changed for {task_id}"
            )
    for field_name in ("n_primary_claims", "min_units"):
        observed = raw_power[field_name]
        if isinstance(observed, bool) or not isinstance(observed, int) or observed < 1:
            raise FirewallError(
                f"SelectionLock power {field_name} is invalid for {task_id}"
            )
    if raw_power["min_units"] < 2 or raw_power["two_sided"] is not True:
        raise FirewallError(f"SelectionLock power contract is invalid for {task_id}")
    return policy, _canonical_mapping(raw_power, f"{task_id}.promotion.power"), sha256_file(
        _PROMOTION_GATE_REGISTRY
    )


def _locked_run_ids(
    lock: SelectionLock, task_id: str
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    decision = _task_decision(lock, task_id)
    selected = tuple(
        sorted(
            _sha256(record.get("run_id"), f"{task_id}.selected_runs.run_id")
            for record in decision.get("selected_runs", ())
            if isinstance(record, Mapping)
        )
    )
    baseline = tuple(
        sorted(
            _sha256(record.get("run_id"), f"{task_id}.baseline_runs.run_id")
            for record in decision.get("baseline_runs", ())
            if isinstance(record, Mapping)
        )
    )
    if not selected or not baseline:
        raise FirewallError(f"SelectionLock has incomplete finalist runs for {task_id}")
    if len(set(selected)) != len(selected) or len(set(baseline)) != len(baseline):
        raise FirewallError(f"SelectionLock repeats finalist runs for {task_id}")
    return selected, baseline


def _validate_resampling_policy(
    policy: object, n_resamples: object, *, allow_fixture_policy: bool
) -> tuple[str, int]:
    checked_policy = _nonempty_string(policy, "resampling_policy")
    if isinstance(n_resamples, bool) or not isinstance(n_resamples, int):
        raise FirewallError("n_resamples must be an integer")
    if checked_policy == PRODUCTION_RESAMPLING_POLICY:
        if n_resamples != PRODUCTION_POWER_RESAMPLES:
            raise FirewallError(
                "production power evidence requires exactly 10,000 endpoint recomputations"
            )
    elif checked_policy == UNIT_TEST_RESAMPLING_POLICY:
        if not allow_fixture_policy:
            raise FirewallError(
                "unit-test power evidence requires explicit fixture-policy opt-in"
            )
        if not 100 <= n_resamples < PRODUCTION_POWER_RESAMPLES:
            raise FirewallError(
                "unit-test power evidence requires 100-9,999 endpoint recomputations"
            )
    else:
        raise FirewallError(f"unsupported resampling_policy: {checked_policy}")
    return checked_policy, n_resamples


_FINALIST_ENDPOINT_RESULT_FIELDS = frozenset(
    {
        "schema_version",
        "aggregation_method",
        "evaluator_id",
        "observed_effect",
        "candidate_primary_metric",
        "baseline_primary_metric",
        "candidate_primary_bootstrap",
        "baseline_primary_bootstrap",
        "paired_difference_bootstrap",
        "candidate_primary_standard_error",
        "baseline_primary_standard_error",
        "paired_effect_standard_error",
        "n_units",
        "independent_unit_counts",
        "unit_set_sha256",
        "n_rows",
        "row_set_sha256",
        "strata",
        "bootstrap_seed",
        "n_resamples",
        "selected_seed_count",
        "development_ensemble_policy_id",
        "ensemble_alignment_sha256",
        "ensemble_rows_sha256",
        "seed_stability",
    }
)
# Development power must be computed for the statistic that is actually
# deployed: the task-native endpoint of the five-seed ensemble prediction.  A
# v1 mean-of-per-seed-endpoint bundle estimates a different quantity and is
# rejected outright.
FINALIST_ENDPOINT_DISTRIBUTION_SCHEMA_VERSION = (
    "masld-bench-finalist-endpoint-distribution-v2"
)
DEVELOPMENT_ENSEMBLE_AGGREGATION = (
    "task_native_endpoint_on_mean_of_exactly_five_locked_seed_predictions"
)
DEVELOPMENT_ENSEMBLE_POLICY_ID = (
    "mean_prediction_five_seed_development_ensemble_v1"
)


def _finite_metric(value: object, field_name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not isfinite(float(value))
    ):
        raise FirewallError(f"{field_name} must be finite numeric")
    return float(value)


def _sample_standard_error(values: list[float]) -> float:
    mean = fsum(values) / len(values)
    return sqrt(fsum((value - mean) ** 2 for value in values) / (len(values) - 1))


def _validated_finalist_endpoint_result(value: object) -> dict[str, Any]:
    raw = _strict_mapping(
        value,
        fields=_FINALIST_ENDPOINT_RESULT_FIELDS,
        label="finalist endpoint distribution",
    )
    if raw["schema_version"] != FINALIST_ENDPOINT_DISTRIBUTION_SCHEMA_VERSION:
        raise FirewallError("unsupported finalist endpoint distribution schema")
    if raw["aggregation_method"] != DEVELOPMENT_ENSEMBLE_AGGREGATION:
        raise FirewallError("unsupported finalist endpoint aggregation method")
    if raw["development_ensemble_policy_id"] != DEVELOPMENT_ENSEMBLE_POLICY_ID:
        raise FirewallError("unsupported development ensemble policy")
    ensemble_alignment_sha256 = _sha256(
        raw["ensemble_alignment_sha256"], "finalist ensemble_alignment_sha256"
    )
    ensemble_rows_sha256 = _sha256(
        raw["ensemble_rows_sha256"], "finalist ensemble_rows_sha256"
    )
    stability = raw["seed_stability"]
    if not isinstance(stability, Mapping):
        raise FirewallError("finalist endpoint seed_stability must be an object")
    if stability.get("used_for_ranking") is not False:
        raise FirewallError("seed stability may never order candidates")
    if stability.get("seeds_are_biological_replicates") is not False:
        raise FirewallError("seeds are not biological replicates")
    evaluator_id = _nonempty_string(
        raw["evaluator_id"], "finalist endpoint evaluator_id"
    )
    integers: dict[str, int] = {}
    for field_name, minimum in (
        ("n_units", 1),
        ("n_rows", 1),
        ("bootstrap_seed", 0),
        ("n_resamples", 1),
        ("selected_seed_count", 1),
    ):
        item = raw[field_name]
        if isinstance(item, bool) or not isinstance(item, int) or item < minimum:
            raise FirewallError(
                f"finalist endpoint {field_name} must be an integer >= {minimum}"
            )
        integers[field_name] = item
    if integers["selected_seed_count"] != 5:
        raise FirewallError("finalist endpoint must aggregate exactly five locked seeds")
    _validate_resampling_policy(
        PRODUCTION_RESAMPLING_POLICY,
        integers["n_resamples"],
        allow_fixture_policy=False,
    )
    metrics = {
        field_name: _finite_metric(raw[field_name], f"finalist endpoint {field_name}")
        for field_name in (
            "observed_effect",
            "candidate_primary_metric",
            "baseline_primary_metric",
            "candidate_primary_standard_error",
            "baseline_primary_standard_error",
            "paired_effect_standard_error",
        )
    }
    if any(
        metrics[field_name] < 0.0
        for field_name in (
            "candidate_primary_standard_error",
            "baseline_primary_standard_error",
            "paired_effect_standard_error",
        )
    ):
        raise FirewallError("finalist endpoint standard errors cannot be negative")
    vectors: dict[str, list[float]] = {}
    for field_name in (
        "candidate_primary_bootstrap",
        "baseline_primary_bootstrap",
        "paired_difference_bootstrap",
    ):
        vector = raw[field_name]
        if not isinstance(vector, list) or len(vector) != integers["n_resamples"]:
            raise FirewallError(
                f"finalist endpoint {field_name} must contain exactly "
                f"{integers['n_resamples']} values"
            )
        vectors[field_name] = [
            _finite_metric(item, f"finalist endpoint {field_name}[{index}]")
            for index, item in enumerate(vector)
        ]
    for index, (candidate, baseline, paired) in enumerate(
        zip(
            vectors["candidate_primary_bootstrap"],
            vectors["baseline_primary_bootstrap"],
            vectors["paired_difference_bootstrap"],
            strict=True,
        )
    ):
        if abs((candidate - baseline) - paired) > 1e-12:
            raise FirewallError(
                "finalist paired bootstrap differs from candidate minus baseline "
                f"at replicate {index}"
            )
    for vector_name, standard_error_name in (
        ("candidate_primary_bootstrap", "candidate_primary_standard_error"),
        ("baseline_primary_bootstrap", "baseline_primary_standard_error"),
        ("paired_difference_bootstrap", "paired_effect_standard_error"),
    ):
        if _sample_standard_error(vectors[vector_name]) != metrics[standard_error_name]:
            raise FirewallError(
                f"finalist endpoint {standard_error_name} does not rederive"
            )
    independent_counts = _positive_integer_mapping(
        raw["independent_unit_counts"],
        "finalist endpoint independent_unit_counts",
    )
    strata = _string_tuple(
        raw["strata"], "finalist endpoint strata", allow_empty=True
    )
    result = {
        "schema_version": FINALIST_ENDPOINT_DISTRIBUTION_SCHEMA_VERSION,
        "aggregation_method": DEVELOPMENT_ENSEMBLE_AGGREGATION,
        "development_ensemble_policy_id": DEVELOPMENT_ENSEMBLE_POLICY_ID,
        "ensemble_alignment_sha256": ensemble_alignment_sha256,
        "ensemble_rows_sha256": ensemble_rows_sha256,
        "seed_stability": _canonical_mapping(
            dict(stability), "finalist endpoint seed_stability"
        ),
        "evaluator_id": evaluator_id,
        **metrics,
        **vectors,
        "n_units": integers["n_units"],
        "independent_unit_counts": independent_counts,
        "unit_set_sha256": _sha256(
            raw["unit_set_sha256"], "finalist endpoint unit_set_sha256"
        ),
        "n_rows": integers["n_rows"],
        "row_set_sha256": _sha256(
            raw["row_set_sha256"], "finalist endpoint row_set_sha256"
        ),
        "strata": list(strata),
        "bootstrap_seed": integers["bootstrap_seed"],
        "n_resamples": integers["n_resamples"],
        "selected_seed_count": integers["selected_seed_count"],
    }
    return _canonical_mapping(result, "finalist endpoint distribution")


def _verified_finalist_metric_source(
    path: str | Path,
) -> tuple[Path, dict[str, Any], dict[str, Any]]:
    try:
        from .tournament import (
            TournamentError,
            verify_finalist_development_metric_bundle,
        )
    except (ImportError, AttributeError) as error:
        raise FirewallError(
            "finalist development metric verifier is unavailable"
        ) from error
    configured = _configured_path(path, "finalist development metric bundle")
    try:
        root = configured.resolve(strict=True)
    except OSError as error:
        raise FirewallError(
            f"finalist development metric bundle is unavailable: {error}"
        ) from error
    if not root.is_dir():
        raise FirewallError("finalist development metric bundle must be a directory")
    try:
        payload = verify_finalist_development_metric_bundle(
            root, reverify_sources=True
        )
        endpoint_ref = ArtifactRef.from_dict(payload["endpoint_result"])
        endpoint_path = endpoint_ref.validate(root, require_relative=True)
    except (TournamentError, ContractError, KeyError, OSError, ValueError) as error:
        raise FirewallError(
            f"invalid finalist development metric bundle: {error}"
        ) from error
    task_id = _nonempty_string(payload.get("task_id"), "finalist metric task_id")
    if (
        endpoint_ref.media_type != "application/json"
        or endpoint_ref.role
        != f"finalist_development_endpoint_distribution:{task_id}"
    ):
        raise FirewallError("finalist endpoint distribution artifact role is invalid")
    endpoint = _validated_finalist_endpoint_result(
        _read_json_object(endpoint_path, "finalist endpoint distribution")
    )
    if canonical_hash(endpoint) != _sha256(
        payload.get("endpoint_result_sha256"),
        "finalist metric endpoint_result_sha256",
    ):
        raise FirewallError("finalist endpoint distribution hash changed")
    summary_fields = (
        "observed_effect",
        "candidate_primary_metric",
        "baseline_primary_metric",
        "n_units",
        "unit_set_sha256",
        "n_rows",
        "row_set_sha256",
    )
    if any(payload.get(field) != endpoint[field] for field in summary_fields):
        raise FirewallError("finalist metric endpoint summary changed")
    if payload.get("endpoint_evaluator_id") != endpoint["evaluator_id"]:
        raise FirewallError("finalist metric endpoint evaluator changed")
    if payload.get("power_method_id") != EMPIRICAL_POWER_METHOD_ID:
        raise FirewallError("finalist metric power method is not registered")
    if (
        payload.get("created_before_outcome_unblind") is not True
        or payload.get("sealed_results_used") is not False
    ):
        raise FirewallError("finalist metric bundle is not pre-unblind evidence")
    return root, payload, endpoint


@dataclass(frozen=True, slots=True)
class DevelopmentPowerEvidence:
    selection_lock_id: str
    selection_lock_manifest_sha256: str
    selection_lock_dir: str
    task_id: str
    selected_run_ids: tuple[str, ...]
    baseline_run_ids: tuple[str, ...]
    promotion_gate_id: str
    promotion_gate_registry_sha256: str
    power_method_id: str
    power_method_source_sha256: str
    finalist_metric_bundle_dir: str
    finalist_metric_bundle_id: str
    finalist_metric_bundle_manifest_sha256: str
    finalist_metric_bundle_document_sha256: str
    selection_candidate_ledger_id: str
    selection_candidate_ledger_manifest_sha256: str
    scientific_receipt_bindings_sha256: str
    development_outcome_bundle_id: str
    development_outcome_manifest_sha256: str
    endpoint_evaluator_id: str
    endpoint_evaluator_sha256: str
    endpoint_implementation_sha256: str
    endpoint_parameters: Mapping[str, Any]
    endpoint_parameters_sha256: str
    bootstrap_seed: int
    n_resamples: int
    resampling_policy: str
    observed_effect: float
    candidate_primary_metric: float
    baseline_primary_metric: float
    n_units: int
    independent_unit_counts: Mapping[str, int]
    unit_set_sha256: str
    n_rows: int
    row_set_sha256: str
    strata: tuple[str, ...]
    bootstrap_distribution: ArtifactRef
    endpoint_result_sha256: str
    power_result: Mapping[str, Any]
    power_result_sha256: str
    evidence_sha256: str

    _FIELDS = frozenset(
        {
            "schema_version",
            *(
                field_name
                for field_name in __annotations__
                if not field_name.startswith("_")
            ),
        }
    )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "DevelopmentPowerEvidence":
        raw = _strict_mapping(
            value, fields=cls._FIELDS, label="DevelopmentPowerEvidence"
        )
        if raw["schema_version"] != DEVELOPMENT_POWER_EVIDENCE_SCHEMA_VERSION:
            raise FirewallError(
                "unsupported DevelopmentPowerEvidence schema_version: "
                f"{raw['schema_version']}"
            )
        paths = {
            field_name: _nonempty_string(raw[field_name], field_name)
            for field_name in (
                "selection_lock_dir",
                "finalist_metric_bundle_dir",
            )
        }
        if any(not Path(path).is_absolute() for path in paths.values()):
            raise FirewallError("DevelopmentPowerEvidence source paths must be absolute")
        try:
            distribution = ArtifactRef.from_dict(raw["bootstrap_distribution"])
        except ContractError as error:
            raise FirewallError(
                f"invalid DevelopmentPowerEvidence distribution artifact: {error}"
            ) from error
        task_id = _nonempty_string(raw["task_id"], "DevelopmentPowerEvidence.task_id")
        if (
            Path(distribution.path).is_absolute()
            or distribution.media_type != "application/json"
            or distribution.role != f"finalist_endpoint_distribution:{task_id}"
        ):
            raise FirewallError(
                "DevelopmentPowerEvidence bootstrap distribution is invalid"
            )
        integers: dict[str, int] = {}
        for field_name, minimum in (
            ("bootstrap_seed", 0),
            ("n_resamples", 1),
            ("n_units", 1),
            ("n_rows", 1),
        ):
            item = raw[field_name]
            if isinstance(item, bool) or not isinstance(item, int) or item < minimum:
                raise FirewallError(
                    f"DevelopmentPowerEvidence.{field_name} must be an integer >= {minimum}"
                )
            integers[field_name] = item
        _validate_resampling_policy(
            raw["resampling_policy"],
            integers["n_resamples"],
            allow_fixture_policy=False,
        )
        if raw["power_method_id"] != EMPIRICAL_POWER_METHOD_ID:
            raise FirewallError("DevelopmentPowerEvidence power method is invalid")
        parameters = _canonical_mapping(
            raw["endpoint_parameters"],
            "DevelopmentPowerEvidence.endpoint_parameters",
        )
        return cls(
            selection_lock_id=_sha256(raw["selection_lock_id"], "selection_lock_id"),
            selection_lock_manifest_sha256=_sha256(
                raw["selection_lock_manifest_sha256"],
                "selection_lock_manifest_sha256",
            ),
            selection_lock_dir=paths["selection_lock_dir"],
            task_id=task_id,
            selected_run_ids=_sha256_tuple(raw["selected_run_ids"], "selected_run_ids"),
            baseline_run_ids=_sha256_tuple(raw["baseline_run_ids"], "baseline_run_ids"),
            promotion_gate_id=_nonempty_string(raw["promotion_gate_id"], "promotion_gate_id"),
            promotion_gate_registry_sha256=_sha256(
                raw["promotion_gate_registry_sha256"],
                "promotion_gate_registry_sha256",
            ),
            power_method_id=EMPIRICAL_POWER_METHOD_ID,
            power_method_source_sha256=_sha256(
                raw["power_method_source_sha256"], "power_method_source_sha256"
            ),
            finalist_metric_bundle_dir=paths["finalist_metric_bundle_dir"],
            finalist_metric_bundle_id=_sha256(
                raw["finalist_metric_bundle_id"], "finalist_metric_bundle_id"
            ),
            finalist_metric_bundle_manifest_sha256=_sha256(
                raw["finalist_metric_bundle_manifest_sha256"],
                "finalist_metric_bundle_manifest_sha256",
            ),
            finalist_metric_bundle_document_sha256=_sha256(
                raw["finalist_metric_bundle_document_sha256"],
                "finalist_metric_bundle_document_sha256",
            ),
            selection_candidate_ledger_id=_sha256(
                raw["selection_candidate_ledger_id"], "selection_candidate_ledger_id"
            ),
            selection_candidate_ledger_manifest_sha256=_sha256(
                raw["selection_candidate_ledger_manifest_sha256"],
                "selection_candidate_ledger_manifest_sha256",
            ),
            scientific_receipt_bindings_sha256=_sha256(
                raw["scientific_receipt_bindings_sha256"],
                "scientific_receipt_bindings_sha256",
            ),
            development_outcome_bundle_id=_sha256(
                raw["development_outcome_bundle_id"], "development_outcome_bundle_id"
            ),
            development_outcome_manifest_sha256=_sha256(
                raw["development_outcome_manifest_sha256"],
                "development_outcome_manifest_sha256",
            ),
            endpoint_evaluator_id=_nonempty_string(
                raw["endpoint_evaluator_id"], "endpoint_evaluator_id"
            ),
            endpoint_evaluator_sha256=_sha256(
                raw["endpoint_evaluator_sha256"], "endpoint_evaluator_sha256"
            ),
            endpoint_implementation_sha256=_sha256(
                raw["endpoint_implementation_sha256"],
                "endpoint_implementation_sha256",
            ),
            endpoint_parameters=parameters,
            endpoint_parameters_sha256=_sha256(
                raw["endpoint_parameters_sha256"], "endpoint_parameters_sha256"
            ),
            bootstrap_seed=integers["bootstrap_seed"],
            n_resamples=integers["n_resamples"],
            resampling_policy=PRODUCTION_RESAMPLING_POLICY,
            observed_effect=_finite_metric(raw["observed_effect"], "observed_effect"),
            candidate_primary_metric=_finite_metric(
                raw["candidate_primary_metric"], "candidate_primary_metric"
            ),
            baseline_primary_metric=_finite_metric(
                raw["baseline_primary_metric"], "baseline_primary_metric"
            ),
            n_units=integers["n_units"],
            independent_unit_counts=_positive_integer_mapping(
                raw["independent_unit_counts"], "independent_unit_counts"
            ),
            unit_set_sha256=_sha256(raw["unit_set_sha256"], "unit_set_sha256"),
            n_rows=integers["n_rows"],
            row_set_sha256=_sha256(raw["row_set_sha256"], "row_set_sha256"),
            strata=_string_tuple(raw["strata"], "strata", allow_empty=True),
            bootstrap_distribution=distribution,
            endpoint_result_sha256=_sha256(
                raw["endpoint_result_sha256"], "endpoint_result_sha256"
            ),
            power_result=_validated_empirical_power_result(raw["power_result"]),
            power_result_sha256=_sha256(
                raw["power_result_sha256"], "power_result_sha256"
            ),
            evidence_sha256=_sha256(raw["evidence_sha256"], "evidence_sha256"),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": DEVELOPMENT_POWER_EVIDENCE_SCHEMA_VERSION,
            **{
                field_name: getattr(self, field_name)
                for field_name in self.__annotations__
                if field_name
                not in {
                    "selected_run_ids",
                    "baseline_run_ids",
                    "endpoint_parameters",
                    "independent_unit_counts",
                    "strata",
                    "bootstrap_distribution",
                    "power_result",
                }
                and not field_name.startswith("_")
            },
            "selected_run_ids": list(self.selected_run_ids),
            "baseline_run_ids": list(self.baseline_run_ids),
            "endpoint_parameters": dict(self.endpoint_parameters),
            "independent_unit_counts": dict(self.independent_unit_counts),
            "strata": list(self.strata),
            "bootstrap_distribution": self.bootstrap_distribution.to_dict(),
            "power_result": dict(self.power_result),
        }


def _development_power_evidence_identity(
    evidence: DevelopmentPowerEvidence,
) -> dict[str, Any]:
    identity = evidence.to_dict()
    identity.pop("evidence_sha256")
    return identity


def _source_evidence_fields(
    *,
    source_root: Path,
    source: Mapping[str, Any],
    endpoint: Mapping[str, Any],
) -> dict[str, Any]:
    lock_binding = _strict_mapping(
        source["selection_lock_binding"],
        fields=frozenset(
            {"path", "manifest_sha256", "document_sha256", "selection_lock_id"}
        ),
        label="finalist metric SelectionLock binding",
    )
    ledger_binding = _strict_mapping(
        source["selection_candidate_ledger_binding"],
        fields=frozenset({"path", "manifest_sha256", "document_sha256", "ledger_id"}),
        label="finalist metric ledger binding",
    )
    receipt_bindings = source.get("scientific_receipt_bindings")
    if not isinstance(receipt_bindings, list) or len(receipt_bindings) != 5:
        raise FirewallError("finalist metric bundle must bind five scientific receipts")
    parameters = _canonical_mapping(
        source["endpoint_parameters"], "finalist metric endpoint_parameters"
    )
    if canonical_hash(parameters) != _sha256(
        source["endpoint_parameters_sha256"],
        "finalist metric endpoint_parameters_sha256",
    ):
        raise FirewallError("finalist metric endpoint parameter hash changed")
    return {
        "selection_lock_id": _sha256(
            lock_binding["selection_lock_id"], "selection_lock_id"
        ),
        "selection_lock_manifest_sha256": _sha256(
            lock_binding["manifest_sha256"], "selection_lock_manifest_sha256"
        ),
        "selection_lock_dir": _nonempty_string(
            lock_binding["path"], "selection_lock_dir"
        ),
        "task_id": _nonempty_string(source["task_id"], "task_id"),
        "selected_run_ids": list(_sha256_tuple(source["selected_run_ids"], "selected_run_ids")),
        "baseline_run_ids": list(_sha256_tuple(source["baseline_run_ids"], "baseline_run_ids")),
        "finalist_metric_bundle_dir": source_root.as_posix(),
        "finalist_metric_bundle_id": _sha256(
            source["metric_bundle_id"], "finalist_metric_bundle_id"
        ),
        "finalist_metric_bundle_manifest_sha256": sha256_file(
            source_root / "ARTIFACTS.json"
        ),
        "finalist_metric_bundle_document_sha256": sha256_file(
            source_root / "finalist_development_metric_bundle.json"
        ),
        "selection_candidate_ledger_id": _sha256(
            ledger_binding["ledger_id"], "selection_candidate_ledger_id"
        ),
        "selection_candidate_ledger_manifest_sha256": _sha256(
            ledger_binding["manifest_sha256"],
            "selection_candidate_ledger_manifest_sha256",
        ),
        "scientific_receipt_bindings_sha256": canonical_hash(receipt_bindings),
        "development_outcome_bundle_id": _sha256(
            source["development_outcome_bundle_id"],
            "development_outcome_bundle_id",
        ),
        "development_outcome_manifest_sha256": _sha256(
            source["development_outcome_manifest_sha256"],
            "development_outcome_manifest_sha256",
        ),
        "endpoint_evaluator_id": _nonempty_string(
            source["endpoint_evaluator_id"], "endpoint_evaluator_id"
        ),
        "endpoint_evaluator_sha256": _sha256(
            source["endpoint_evaluator_sha256"], "endpoint_evaluator_sha256"
        ),
        "endpoint_implementation_sha256": _endpoint_evaluator_source_sha256(),
        "endpoint_parameters": parameters,
        "endpoint_parameters_sha256": canonical_hash(parameters),
        "bootstrap_seed": endpoint["bootstrap_seed"],
        "n_resamples": endpoint["n_resamples"],
        "resampling_policy": PRODUCTION_RESAMPLING_POLICY,
        "observed_effect": endpoint["observed_effect"],
        "candidate_primary_metric": endpoint["candidate_primary_metric"],
        "baseline_primary_metric": endpoint["baseline_primary_metric"],
        "n_units": endpoint["n_units"],
        "independent_unit_counts": dict(endpoint["independent_unit_counts"]),
        "unit_set_sha256": endpoint["unit_set_sha256"],
        "n_rows": endpoint["n_rows"],
        "row_set_sha256": endpoint["row_set_sha256"],
        "strata": list(endpoint["strata"]),
        "endpoint_result_sha256": canonical_hash(endpoint),
    }


def _recompute_power_result(
    *, endpoint: Mapping[str, Any], power_spec: Mapping[str, Any]
) -> dict[str, Any]:
    try:
        result = empirical_bootstrap_power_gate(
            paired_difference_bootstrap=endpoint["paired_difference_bootstrap"],
            observed_effect=float(endpoint["observed_effect"]),
            n_units=int(endpoint["n_units"]),
            minimum_effect=float(power_spec["minimum_effect"]),
            alpha=float(power_spec["alpha"]),
            power_target=float(power_spec["target_power"]),
            n_primary_claims=int(power_spec["n_primary_claims"]),
            min_units=int(power_spec["min_units"]),
            two_sided=bool(power_spec["two_sided"]),
        )
    except (KeyError, TypeError, ValueError) as error:
        raise FirewallError(f"cannot derive empirical power result: {error}") from error
    return _validated_empirical_power_result(result)


def freeze_development_power_evidence(
    *,
    finalist_development_metric_bundle_dir: str | Path,
    output_root: str | Path,
) -> Path:
    """Freeze power only from a recursively verified five-seed metric bundle."""

    source_root, source, endpoint = _verified_finalist_metric_source(
        finalist_development_metric_bundle_dir
    )
    source_fields = _source_evidence_fields(
        source_root=source_root, source=source, endpoint=endpoint
    )
    selection_dir = _configured_path(
        source_fields["selection_lock_dir"], "DevelopmentPowerEvidence SelectionLock"
    ).resolve()
    lock = _verify_selection_dir(selection_dir)
    policy, power_spec, registry_sha256 = _load_task_power_contract(
        lock, source_fields["task_id"]
    )
    selected_run_ids, baseline_run_ids = _locked_run_ids(
        lock, source_fields["task_id"]
    )
    if (
        lock.lock_id != source_fields["selection_lock_id"]
        or selected_run_ids != tuple(source_fields["selected_run_ids"])
        or baseline_run_ids != tuple(source_fields["baseline_run_ids"])
    ):
        raise FirewallError("finalist metric bundle does not bind the locked finalists")
    expected_evaluator = _TASK_ENDPOINT_EVALUATORS[source_fields["task_id"]]
    if source_fields["endpoint_evaluator_id"] != expected_evaluator:
        raise FirewallError("finalist metric bundle uses the wrong endpoint evaluator")
    power_result = _recompute_power_result(endpoint=endpoint, power_spec=power_spec)
    output = _configured_path(output_root, "development-power output root").resolve()
    collection = output / "development_power_evidence"
    collection.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".development_power_evidence.", dir=collection))
    try:
        distribution_path = staging / "bootstrap_distribution.json"
        write_json_exclusive(distribution_path, endpoint, mode=0o440)
        distribution_ref = ArtifactRef.from_path(
            distribution_path,
            relative_to=staging,
            media_type="application/json",
            role=f"finalist_endpoint_distribution:{source_fields['task_id']}",
        )
        identity = {
            "schema_version": DEVELOPMENT_POWER_EVIDENCE_SCHEMA_VERSION,
            **source_fields,
            "promotion_gate_id": policy["promotion_gate_id"],
            "promotion_gate_registry_sha256": registry_sha256,
            "power_method_id": EMPIRICAL_POWER_METHOD_ID,
            "power_method_source_sha256": _power_method_source_sha256(),
            "bootstrap_distribution": distribution_ref.to_dict(),
            "power_result": power_result,
            "power_result_sha256": canonical_hash(power_result),
        }
        evidence_sha256 = canonical_hash(identity)
        evidence = DevelopmentPowerEvidence.from_dict(
            {**identity, "evidence_sha256": evidence_sha256}
        )
        write_json_exclusive(
            staging / "development_power_evidence.json",
            evidence.to_dict(),
            mode=0o440,
        )
        freeze_tree(
            staging,
            {
                "artifact_class": "development_power_evidence",
                "artifact_id": evidence_sha256,
            },
        )
        target = collection / evidence_sha256
        try:
            publish_directory_noreplace(staging, target)
        except ArtifactError as error:
            raise FirewallError(
                f"could not publish immutable development_power_evidence: {error}"
            ) from error
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_development_power_evidence(target, reverify_sources=True)
    return target


def verify_development_power_evidence(
    path: str | Path, *, reverify_sources: bool = True
) -> DevelopmentPowerEvidence:
    configured = _configured_path(path, "development-power evidence")
    try:
        root = configured.resolve(strict=True)
        verify_frozen_tree(root)
    except (ArtifactError, OSError, ValueError, KeyError, TypeError) as error:
        raise FirewallError(
            f"invalid frozen development_power_evidence: {error}"
        ) from error
    evidence = DevelopmentPowerEvidence.from_dict(
        _read_json_object(
            root / "development_power_evidence.json", "DevelopmentPowerEvidence"
        )
    )
    _verify_state_manifest(
        root,
        artifact_class="development_power_evidence",
        artifact_id=evidence.evidence_sha256,
    )
    if canonical_hash(_development_power_evidence_identity(evidence)) != evidence.evidence_sha256:
        raise FirewallError("DevelopmentPowerEvidence identity hash mismatch")
    if canonical_hash(evidence.power_result) != evidence.power_result_sha256:
        raise FirewallError("DevelopmentPowerEvidence power-result hash mismatch")
    if canonical_hash(evidence.endpoint_parameters) != evidence.endpoint_parameters_sha256:
        raise FirewallError("DevelopmentPowerEvidence parameter hash mismatch")
    try:
        distribution_path = evidence.bootstrap_distribution.validate(
            root, require_relative=True
        )
    except ContractError as error:
        raise FirewallError(
            f"DevelopmentPowerEvidence distribution changed: {error}"
        ) from error
    endpoint = _validated_finalist_endpoint_result(
        _read_json_object(distribution_path, "finalist endpoint distribution")
    )
    if canonical_hash(endpoint) != evidence.endpoint_result_sha256:
        raise FirewallError("DevelopmentPowerEvidence endpoint-result hash mismatch")
    expected_summary = {
        "bootstrap_seed": evidence.bootstrap_seed,
        "n_resamples": evidence.n_resamples,
        "observed_effect": evidence.observed_effect,
        "candidate_primary_metric": evidence.candidate_primary_metric,
        "baseline_primary_metric": evidence.baseline_primary_metric,
        "n_units": evidence.n_units,
        "independent_unit_counts": dict(evidence.independent_unit_counts),
        "unit_set_sha256": evidence.unit_set_sha256,
        "n_rows": evidence.n_rows,
        "row_set_sha256": evidence.row_set_sha256,
        "strata": list(evidence.strata),
    }
    if any(endpoint[field] != expected for field, expected in expected_summary.items()):
        raise FirewallError("DevelopmentPowerEvidence endpoint summary changed")
    if not reverify_sources:
        return evidence
    source_root, source, source_endpoint = _verified_finalist_metric_source(
        evidence.finalist_metric_bundle_dir
    )
    source_fields = _source_evidence_fields(
        source_root=source_root, source=source, endpoint=source_endpoint
    )
    serialized_evidence = evidence.to_dict()
    observed_source_fields = {
        field_name: serialized_evidence[field_name]
        for field_name in source_fields
    }
    if observed_source_fields != source_fields or source_endpoint != endpoint:
        raise FirewallError("DevelopmentPowerEvidence finalist source changed")
    selection_dir = _configured_path(
        evidence.selection_lock_dir, "DevelopmentPowerEvidence SelectionLock"
    )
    lock = _verify_selection_dir(selection_dir)
    policy, power_spec, registry_sha256 = _load_task_power_contract(
        lock, evidence.task_id
    )
    if (
        evidence.promotion_gate_id != policy["promotion_gate_id"]
        or evidence.promotion_gate_registry_sha256 != registry_sha256
        or evidence.power_method_id != policy["power_method_id"]
        or evidence.power_method_source_sha256 != _power_method_source_sha256()
        or evidence.endpoint_implementation_sha256
        != _endpoint_evaluator_source_sha256()
    ):
        raise FirewallError("DevelopmentPowerEvidence executable authority changed")
    recomputed_power = _recompute_power_result(
        endpoint=source_endpoint, power_spec=power_spec
    )
    if recomputed_power != evidence.power_result:
        raise FirewallError("DevelopmentPowerEvidence power result changed")
    return evidence


_PREDICTION_BINDING_FIELDS = frozenset(
    {
        "run_id",
        "model_role",
        "prediction_commit_sha256",
        "prediction_bundle_sha256",
    }
)


def _parse_prediction_binding(value: object, label: str) -> dict[str, str]:
    raw = _strict_mapping(value, fields=_PREDICTION_BINDING_FIELDS, label=label)
    model_role = _nonempty_string(raw["model_role"], f"{label}.model_role")
    if model_role not in {
        "selected",
        "baseline",
        "selected_and_baseline",
        "secondary_comparator",
    }:
        raise FirewallError(f"{label}.model_role is invalid")
    return {
        "run_id": _sha256(raw["run_id"], f"{label}.run_id"),
        "model_role": model_role,
        "prediction_commit_sha256": _sha256(
            raw["prediction_commit_sha256"],
            f"{label}.prediction_commit_sha256",
        ),
        "prediction_bundle_sha256": _sha256(
            raw["prediction_bundle_sha256"],
            f"{label}.prediction_bundle_sha256",
        ),
    }


@dataclass(frozen=True, slots=True)
class PowerDecision:
    selection_lock_id: str
    selection_lock_manifest_sha256: str
    task_id: str
    dataset_ids: tuple[str, ...]
    dataset_bindings_sha256: str
    split_id: str
    row_id_field: str
    unit_id_field: str
    unit_id_namespace: str
    biological_unit: str
    table_schema_sha256: str
    source_join_key_sha256: str
    n_predictions: int
    prediction_bindings: tuple[Mapping[str, str], ...]
    prediction_set_sha256: str
    development_power_evidence_dir: str
    development_power_evidence_sha256: str
    development_power_evidence_manifest_sha256: str
    power_method_id: str
    power_result: Mapping[str, Any]
    power_result_sha256: str
    power_decision_sha256: str

    _FIELDS = frozenset(
        {
            "schema_version",
            "selection_lock_id",
            "selection_lock_manifest_sha256",
            "task_id",
            "dataset_ids",
            "dataset_bindings_sha256",
            "split_id",
            "row_id_field",
            "unit_id_field",
            "unit_id_namespace",
            "biological_unit",
            "table_schema_sha256",
            "source_join_key_sha256",
            "n_predictions",
            "prediction_bindings",
            "prediction_set_sha256",
            "development_power_evidence_dir",
            "development_power_evidence_sha256",
            "development_power_evidence_manifest_sha256",
            "power_method_id",
            "power_result",
            "power_result_sha256",
            "power_decision_sha256",
        }
    )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PowerDecision":
        raw = _strict_mapping(value, fields=cls._FIELDS, label="PowerDecision")
        if raw["schema_version"] != POWER_DECISION_SCHEMA_VERSION:
            raise FirewallError(
                f"unsupported PowerDecision schema_version: {raw['schema_version']}"
            )
        bindings_raw = raw["prediction_bindings"]
        if not isinstance(bindings_raw, (list, tuple)) or not bindings_raw:
            raise FirewallError("PowerDecision.prediction_bindings must not be empty")
        bindings = tuple(
            _parse_prediction_binding(item, f"prediction_bindings[{index}]")
            for index, item in enumerate(bindings_raw)
        )
        if [item["run_id"] for item in bindings] != sorted(
            item["run_id"] for item in bindings
        ):
            raise FirewallError("PowerDecision.prediction_bindings must be sorted by run_id")
        if len({item["run_id"] for item in bindings}) != len(bindings):
            raise FirewallError("PowerDecision.prediction_bindings repeat a run_id")
        power_result = _validated_empirical_power_result(raw["power_result"])
        evidence_dir = _nonempty_string(
            raw["development_power_evidence_dir"],
            "PowerDecision.development_power_evidence_dir",
        )
        if not Path(evidence_dir).is_absolute():
            raise FirewallError(
                "PowerDecision.development_power_evidence_dir must be absolute"
            )
        power_method_id = _nonempty_string(
            raw["power_method_id"], "PowerDecision.power_method_id"
        )
        if power_method_id != EMPIRICAL_POWER_METHOD_ID:
            raise FirewallError(
                "PowerDecision.power_method_id is not the production method"
            )
        n_predictions = raw["n_predictions"]
        if (
            isinstance(n_predictions, bool)
            or not isinstance(n_predictions, int)
            or n_predictions < 1
        ):
            raise FirewallError("PowerDecision.n_predictions must be a positive integer")
        return cls(
            selection_lock_id=_sha256(
                raw["selection_lock_id"], "PowerDecision.selection_lock_id"
            ),
            selection_lock_manifest_sha256=_sha256(
                raw["selection_lock_manifest_sha256"],
                "PowerDecision.selection_lock_manifest_sha256",
            ),
            task_id=_nonempty_string(raw["task_id"], "PowerDecision.task_id"),
            dataset_ids=_string_tuple(raw["dataset_ids"], "PowerDecision.dataset_ids"),
            dataset_bindings_sha256=_sha256(
                raw["dataset_bindings_sha256"],
                "PowerDecision.dataset_bindings_sha256",
            ),
            split_id=_nonempty_string(raw["split_id"], "PowerDecision.split_id"),
            row_id_field=_nonempty_string(
                raw["row_id_field"], "PowerDecision.row_id_field"
            ),
            unit_id_field=_nonempty_string(
                raw["unit_id_field"], "PowerDecision.unit_id_field"
            ),
            unit_id_namespace=_nonempty_string(
                raw["unit_id_namespace"], "PowerDecision.unit_id_namespace"
            ),
            biological_unit=_nonempty_string(
                raw["biological_unit"], "PowerDecision.biological_unit"
            ),
            table_schema_sha256=_sha256(
                raw["table_schema_sha256"], "PowerDecision.table_schema_sha256"
            ),
            source_join_key_sha256=_sha256(
                raw["source_join_key_sha256"],
                "PowerDecision.source_join_key_sha256",
            ),
            n_predictions=n_predictions,
            prediction_bindings=bindings,
            prediction_set_sha256=_sha256(
                raw["prediction_set_sha256"], "PowerDecision.prediction_set_sha256"
            ),
            development_power_evidence_dir=evidence_dir,
            development_power_evidence_sha256=_sha256(
                raw["development_power_evidence_sha256"],
                "PowerDecision.development_power_evidence_sha256",
            ),
            development_power_evidence_manifest_sha256=_sha256(
                raw["development_power_evidence_manifest_sha256"],
                "PowerDecision.development_power_evidence_manifest_sha256",
            ),
            power_method_id=power_method_id,
            power_result=power_result,
            power_result_sha256=_sha256(
                raw["power_result_sha256"], "PowerDecision.power_result_sha256"
            ),
            power_decision_sha256=_sha256(
                raw["power_decision_sha256"], "PowerDecision.power_decision_sha256"
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": POWER_DECISION_SCHEMA_VERSION,
            "selection_lock_id": self.selection_lock_id,
            "selection_lock_manifest_sha256": self.selection_lock_manifest_sha256,
            "task_id": self.task_id,
            "dataset_ids": list(self.dataset_ids),
            "dataset_bindings_sha256": self.dataset_bindings_sha256,
            "split_id": self.split_id,
            "row_id_field": self.row_id_field,
            "unit_id_field": self.unit_id_field,
            "unit_id_namespace": self.unit_id_namespace,
            "biological_unit": self.biological_unit,
            "table_schema_sha256": self.table_schema_sha256,
            "source_join_key_sha256": self.source_join_key_sha256,
            "n_predictions": self.n_predictions,
            "prediction_bindings": [dict(item) for item in self.prediction_bindings],
            "prediction_set_sha256": self.prediction_set_sha256,
            "development_power_evidence_dir": self.development_power_evidence_dir,
            "development_power_evidence_sha256": self.development_power_evidence_sha256,
            "development_power_evidence_manifest_sha256": (
                self.development_power_evidence_manifest_sha256
            ),
            "power_method_id": self.power_method_id,
            "power_result": dict(self.power_result),
            "power_result_sha256": self.power_result_sha256,
            "power_decision_sha256": self.power_decision_sha256,
        }


def _power_decision_identity(decision: PowerDecision) -> dict[str, Any]:
    identity = decision.to_dict()
    identity.pop("power_decision_sha256")
    return identity


def _commit_dir(state: Path, commit_sha256: str) -> Path:
    root = _configured_path(state, "evaluator-state root").resolve()
    return root / "prediction_commits" / _sha256(
        commit_sha256, "prediction_commit_sha256"
    )


def _power_dir(state: Path, power_sha256: str) -> Path:
    root = _configured_path(state, "evaluator-state root").resolve()
    return root / "power_decisions" / _sha256(
        power_sha256, "power_decision_sha256"
    )


def _complete_prediction_set(
    commits: tuple[PredictionCommit, ...],
) -> tuple[SelectionLock, tuple[Mapping[str, str], ...]]:
    first = commits[0]
    if any(commit.task_id != first.task_id for commit in commits):
        raise FirewallError("a PowerDecision may bind only one task")
    if any(commit.dataset_ids != first.dataset_ids for commit in commits):
        raise FirewallError("prediction commits disagree on sealed datasets")
    if any(
        commit.dataset_bindings_sha256 != first.dataset_bindings_sha256
        for commit in commits
    ):
        raise FirewallError("prediction commits disagree on dataset bindings")
    if any(commit.split_id != first.split_id for commit in commits):
        raise FirewallError("selected and baseline predictions use different sealed splits")
    if any(commit.row_id_field != first.row_id_field for commit in commits):
        raise FirewallError("selected and baseline predictions use different row identifiers")
    if any(commit.unit_id_field != first.unit_id_field for commit in commits):
        raise FirewallError("selected and baseline predictions use different unit identifiers")
    if any(commit.unit_id_namespace != first.unit_id_namespace for commit in commits):
        raise FirewallError("selected and baseline predictions use different unit namespaces")
    if any(commit.biological_unit != first.biological_unit for commit in commits):
        raise FirewallError("selected and baseline predictions use different biological units")
    if any(commit.table_schema_sha256 != first.table_schema_sha256 for commit in commits):
        raise FirewallError("selected and baseline predictions use different table schemas")
    if any(
        commit.source_join_key_sha256 != first.source_join_key_sha256
        for commit in commits
    ):
        raise FirewallError("selected and baseline predictions use different source join keys")
    if any(commit.n_predictions != first.n_predictions for commit in commits):
        raise FirewallError("selected and baseline prediction counts differ")
    if any(commit.selection_lock_id != first.selection_lock_id for commit in commits):
        raise FirewallError("prediction commits use different SelectionLocks")
    if any(
        commit.selection_lock_manifest_sha256 != first.selection_lock_manifest_sha256
        for commit in commits
    ):
        raise FirewallError("prediction commits use different SelectionLock manifests")
    lock = _verify_selection_dir(Path(first.selection_lock_dir))
    decision = _task_decision(lock, first.task_id)
    selected_run_ids = {
        str(record["run_id"]) for record in decision.get("selected_runs", ())
    }
    baseline_run_ids = {
        str(record["run_id"]) for record in decision.get("baseline_runs", ())
    }
    secondary_run_ids: set[str] = set()
    secondary = decision.get("variant_secondary_evaluation")
    if first.task_id == _VARIANT_TASK_ID:
        if not isinstance(secondary, Mapping):
            raise FirewallError(
                "variant power decision lacks its locked secondary evaluation"
            )
        raw_comparators = secondary.get("comparators")
        if not isinstance(raw_comparators, (list, tuple)) or not raw_comparators:
            raise FirewallError("variant secondary comparator roster is invalid")
        secondary_run_ids = {
            str(record.get("run_id"))
            for record in raw_comparators
            if isinstance(record, Mapping)
        }
        if len(secondary_run_ids) != len(raw_comparators):
            raise FirewallError("variant secondary comparator roster repeats a run")
    elif secondary is not None:
        raise FirewallError("nonvariant task carries a variant secondary evaluation")
    if secondary_run_ids.intersection(selected_run_ids | baseline_run_ids):
        raise FirewallError("secondary comparator runs overlap signed primary runs")
    expected_roles = {
        run_id: (
            "selected_and_baseline"
            if run_id in selected_run_ids and run_id in baseline_run_ids
            else "selected" if run_id in selected_run_ids else "baseline"
        )
        for run_id in selected_run_ids | baseline_run_ids
    }
    expected_roles.update(
        {run_id: "secondary_comparator" for run_id in secondary_run_ids}
    )
    observed_roles = {commit.run_id: commit.model_role for commit in commits}
    if len(observed_roles) != len(commits) or observed_roles != expected_roles:
        raise FirewallError(
            f"power decision for {first.task_id} requires every locked selected and "
            "baseline run and every secondary comparator run"
        )
    if len({commit.prediction_commit_sha256 for commit in commits}) != len(commits):
        raise FirewallError("prediction commit set contains duplicates")
    by_run = sorted(commits, key=lambda item: item.run_id)
    bindings = tuple(
        {
            "run_id": commit.run_id,
            "model_role": commit.model_role,
            "prediction_commit_sha256": commit.prediction_commit_sha256,
            "prediction_bundle_sha256": commit.prediction_bundle_sha256,
        }
        for commit in by_run
    )
    return lock, bindings


def freeze_power_decision(
    *,
    evaluator_state_dir: str | Path,
    prediction_commit_sha256s: Iterable[str],
    development_power_evidence_dir: str | Path,
) -> PowerDecision:
    """Bind a complete sealed prediction set to verified development power."""

    state = _configured_path(evaluator_state_dir, "evaluator-state root").resolve()
    commit_ids = tuple(prediction_commit_sha256s)
    if not commit_ids or len(set(commit_ids)) != len(commit_ids):
        raise FirewallError(
            "prediction_commit_sha256s must be a non-empty set without duplicates"
        )
    commits = tuple(
        verify_prediction_commit(_commit_dir(state, commit_id), reverify_sources=True)
        for commit_id in commit_ids
    )
    lock, bindings = _complete_prediction_set(commits)
    first = commits[0]
    evidence_dir = _configured_path(
        development_power_evidence_dir, "development-power evidence"
    ).resolve()
    evidence = verify_development_power_evidence(
        evidence_dir,
        reverify_sources=True,
    )
    if (
        evidence.selection_lock_id != lock.lock_id
        or evidence.selection_lock_manifest_sha256
        != first.selection_lock_manifest_sha256
        or evidence.task_id != first.task_id
    ):
        raise FirewallError(
            "DevelopmentPowerEvidence does not bind this SelectionLock task"
        )
    selected_run_ids = tuple(
        sorted(
            binding["run_id"]
            for binding in bindings
            if binding["model_role"] in {"selected", "selected_and_baseline"}
        )
    )
    baseline_run_ids = tuple(
        sorted(
            binding["run_id"]
            for binding in bindings
            if binding["model_role"] in {"baseline", "selected_and_baseline"}
        )
    )
    if (
        selected_run_ids != evidence.selected_run_ids
        or baseline_run_ids != evidence.baseline_run_ids
    ):
        raise FirewallError(
            "DevelopmentPowerEvidence does not bind the complete finalist run set"
        )
    identity = {
        "schema_version": POWER_DECISION_SCHEMA_VERSION,
        "selection_lock_id": lock.lock_id,
        "selection_lock_manifest_sha256": first.selection_lock_manifest_sha256,
        "task_id": first.task_id,
        "dataset_ids": list(first.dataset_ids),
        "dataset_bindings_sha256": first.dataset_bindings_sha256,
        "split_id": first.split_id,
        "row_id_field": first.row_id_field,
        "unit_id_field": first.unit_id_field,
        "unit_id_namespace": first.unit_id_namespace,
        "biological_unit": first.biological_unit,
        "table_schema_sha256": first.table_schema_sha256,
        "source_join_key_sha256": first.source_join_key_sha256,
        "n_predictions": first.n_predictions,
        "prediction_bindings": [dict(item) for item in bindings],
        "prediction_set_sha256": canonical_hash(bindings),
        "development_power_evidence_dir": evidence_dir.as_posix(),
        "development_power_evidence_sha256": evidence.evidence_sha256,
        "development_power_evidence_manifest_sha256": sha256_file(
            evidence_dir / "ARTIFACTS.json"
        ),
        "power_method_id": evidence.power_method_id,
        "power_result": dict(evidence.power_result),
        "power_result_sha256": evidence.power_result_sha256,
    }
    payload = {**identity, "power_decision_sha256": canonical_hash(identity)}
    decision = PowerDecision.from_dict(payload)
    decision_dir = _freeze_state_document(
        state_root=state,
        collection="power_decisions",
        identifier=decision.power_decision_sha256,
        filename="power_decision.json",
        payload=decision.to_dict(),
        artifact_class="power_decision",
    )
    return verify_power_decision(decision_dir, reverify_sources=True)


def verify_power_decision(
    path: str | Path, *, reverify_sources: bool = True
) -> PowerDecision:
    root = _configured_path(path, "PowerDecision").resolve()
    try:
        verify_frozen_tree(root)
    except (ArtifactError, OSError, ValueError, KeyError, TypeError) as error:
        raise FirewallError(f"invalid frozen power_decision: {error}") from error
    payload = _read_json_object(root / "power_decision.json", "PowerDecision")
    decision = PowerDecision.from_dict(payload)
    _verify_state_manifest(
        root,
        artifact_class="power_decision",
        artifact_id=decision.power_decision_sha256,
    )
    if canonical_hash(_power_decision_identity(decision)) != decision.power_decision_sha256:
        raise FirewallError("PowerDecision identity hash mismatch")
    if canonical_hash(decision.power_result) != decision.power_result_sha256:
        raise FirewallError("PowerDecision power-result hash mismatch")
    if canonical_hash(decision.prediction_bindings) != decision.prediction_set_sha256:
        raise FirewallError("PowerDecision prediction-set hash mismatch")
    if not reverify_sources:
        return decision
    if root.parent.name != "power_decisions":
        raise FirewallError("PowerDecision is outside an evaluator-state power_decisions tree")
    state = root.parent.parent
    commits = tuple(
        verify_prediction_commit(
            _commit_dir(state, binding["prediction_commit_sha256"]),
            reverify_sources=True,
        )
        for binding in decision.prediction_bindings
    )
    lock, bindings = _complete_prediction_set(commits)
    first = commits[0]
    evidence_dir = _configured_path(
        decision.development_power_evidence_dir,
        "PowerDecision development-power evidence",
    )
    evidence = verify_development_power_evidence(
        evidence_dir,
        reverify_sources=True,
    )
    try:
        evidence_manifest_sha256 = sha256_file(evidence_dir / "ARTIFACTS.json")
    except OSError as error:
        raise FirewallError(f"DevelopmentPowerEvidence manifest is unreadable: {error}") from error
    if tuple(bindings) != decision.prediction_bindings:
        raise FirewallError("PowerDecision prediction bindings changed")
    if (
        lock.lock_id != decision.selection_lock_id
        or first.selection_lock_manifest_sha256
        != decision.selection_lock_manifest_sha256
        or first.task_id != decision.task_id
        or first.dataset_ids != decision.dataset_ids
        or first.dataset_bindings_sha256 != decision.dataset_bindings_sha256
        or first.split_id != decision.split_id
        or first.row_id_field != decision.row_id_field
        or first.unit_id_field != decision.unit_id_field
        or first.unit_id_namespace != decision.unit_id_namespace
        or first.biological_unit != decision.biological_unit
        or first.table_schema_sha256 != decision.table_schema_sha256
        or first.source_join_key_sha256 != decision.source_join_key_sha256
        or first.n_predictions != decision.n_predictions
        or evidence.evidence_sha256
        != decision.development_power_evidence_sha256
        or evidence_manifest_sha256
        != decision.development_power_evidence_manifest_sha256
        or evidence.selection_lock_id != decision.selection_lock_id
        or evidence.selection_lock_manifest_sha256
        != decision.selection_lock_manifest_sha256
        or evidence.task_id != decision.task_id
        or evidence.power_method_id != decision.power_method_id
        or evidence.power_result != decision.power_result
        or evidence.power_result_sha256 != decision.power_result_sha256
    ):
        raise FirewallError("PowerDecision source bindings do not reverify")
    selected_run_ids = tuple(
        sorted(
            binding["run_id"]
            for binding in bindings
            if binding["model_role"] in {"selected", "selected_and_baseline"}
        )
    )
    baseline_run_ids = tuple(
        sorted(
            binding["run_id"]
            for binding in bindings
            if binding["model_role"] in {"baseline", "selected_and_baseline"}
        )
    )
    if (
        selected_run_ids != evidence.selected_run_ids
        or baseline_run_ids != evidence.baseline_run_ids
    ):
        raise FirewallError("PowerDecision finalist run evidence changed")
    return decision


_OUTCOME_TASK_BINDING_FIELDS = frozenset({"task_id", "dataset_ids"})


def _parse_outcome_task_binding(value: object, label: str) -> dict[str, Any]:
    raw = _strict_mapping(value, fields=_OUTCOME_TASK_BINDING_FIELDS, label=label)
    return {
        "task_id": _nonempty_string(raw["task_id"], f"{label}.task_id"),
        "dataset_ids": list(_string_tuple(raw["dataset_ids"], f"{label}.dataset_ids")),
    }


@dataclass(frozen=True, slots=True)
class SealedOutcomeBundle:
    bundle_id: str
    bundle_sha256: str
    joint_bundle: bool
    task_bindings: tuple[Mapping[str, Any], ...]
    unit_id_fields: Mapping[str, str]
    artifacts: tuple[ArtifactRef, ...]

    _FIELDS = frozenset(
        {
            "schema_version",
            "bundle_id",
            "bundle_sha256",
            "joint_bundle",
            "task_bindings",
            "unit_id_fields",
            "artifacts",
        }
    )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SealedOutcomeBundle":
        raw = _strict_mapping(value, fields=cls._FIELDS, label="SealedOutcomeBundle")
        if raw["schema_version"] != SEALED_OUTCOME_BUNDLE_SCHEMA_VERSION:
            raise FirewallError(
                "unsupported SealedOutcomeBundle schema_version: "
                f"{raw['schema_version']}"
            )
        if raw["joint_bundle"] is not True:
            raise FirewallError("sealed outcomes must be consumed with joint_bundle=true")
        bindings_raw = raw["task_bindings"]
        if not isinstance(bindings_raw, (list, tuple)) or not bindings_raw:
            raise FirewallError("SealedOutcomeBundle.task_bindings must not be empty")
        bindings = tuple(
            _parse_outcome_task_binding(item, f"task_bindings[{index}]")
            for index, item in enumerate(bindings_raw)
        )
        task_ids = [str(item["task_id"]) for item in bindings]
        if task_ids != sorted(task_ids) or len(set(task_ids)) != len(task_ids):
            raise FirewallError(
                "SealedOutcomeBundle.task_bindings must have unique sorted task IDs"
            )
        unit_fields = _canonical_mapping(
            raw["unit_id_fields"], "SealedOutcomeBundle.unit_id_fields"
        )
        if set(unit_fields) != set(task_ids):
            raise FirewallError("unit_id_fields must cover outcome tasks exactly")
        for task_id, field_name in unit_fields.items():
            _nonempty_string(field_name, f"unit_id_fields[{task_id!r}]")
        artifacts_raw = raw["artifacts"]
        if not isinstance(artifacts_raw, (list, tuple)) or not artifacts_raw:
            raise FirewallError("SealedOutcomeBundle.artifacts must not be empty")
        try:
            artifacts = tuple(ArtifactRef.from_dict(item) for item in artifacts_raw)
        except ContractError as error:
            raise FirewallError(f"invalid sealed outcome artifact: {error}") from error
        if len({artifact.path for artifact in artifacts}) != len(artifacts):
            raise FirewallError("sealed outcome artifact paths must be unique")
        bundle = cls(
            bundle_id=_nonempty_string(
                raw["bundle_id"], "SealedOutcomeBundle.bundle_id"
            ),
            bundle_sha256=_sha256(
                raw["bundle_sha256"], "SealedOutcomeBundle.bundle_sha256"
            ),
            joint_bundle=True,
            task_bindings=bindings,
            unit_id_fields=unit_fields,
            artifacts=artifacts,
        )
        if canonical_hash(bundle.identity_payload()) != bundle.bundle_sha256:
            raise FirewallError("SealedOutcomeBundle identity hash mismatch")
        return bundle

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": SEALED_OUTCOME_BUNDLE_SCHEMA_VERSION,
            "bundle_id": self.bundle_id,
            "joint_bundle": self.joint_bundle,
            "task_bindings": [dict(item) for item in self.task_bindings],
            "unit_id_fields": dict(self.unit_id_fields),
            "artifacts": [artifact.to_dict() for artifact in self.artifacts],
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            **self.identity_payload(),
            "bundle_sha256": self.bundle_sha256,
        }

    def validate_artifacts(self, root: str | Path) -> tuple[Path, ...]:
        return tuple(
            artifact.validate(root, require_relative=True)
            for artifact in self.artifacts
        )


_ARTIFACT_RECORD_FIELDS = frozenset({"path", "sha256", "size_bytes"})
_FROZEN_MANIFEST_FIELDS = frozenset(
    {"schema_version", "metadata", "artifacts"}
)
_COMPLETE_FIELDS = frozenset(
    {"schema_version", "manifest_sha256", "artifact_count"}
)


def _preflight_sealed_outcome_bundle(
    path: str | Path,
) -> tuple[Path, SealedOutcomeBundle, str]:
    """Verify outcome metadata without hashing or opening outcome artifacts."""

    root = _configured_path(path, "sealed outcome bundle").resolve()
    if not root.is_dir():
        raise FirewallError(f"sealed outcome bundle must be a directory: {root}")
    manifest_path = root / "ARTIFACTS.json"
    complete_path = root / "COMPLETE"
    manifest = _strict_mapping(
        _read_json_object(manifest_path, "sealed outcome ARTIFACTS.json"),
        fields=_FROZEN_MANIFEST_FIELDS,
        label="sealed outcome ARTIFACTS.json",
    )
    complete = _strict_mapping(
        _read_json_object(complete_path, "sealed outcome COMPLETE"),
        fields=_COMPLETE_FIELDS,
        label="sealed outcome COMPLETE",
    )
    if manifest["schema_version"] != "masld-bench-artifacts-v1":
        raise FirewallError("unsupported sealed outcome artifact-manifest schema")
    if complete["schema_version"] != "masld-bench-complete-v1":
        raise FirewallError("unsupported sealed outcome completion schema")
    manifest_sha256 = sha256_file(manifest_path)
    if complete["manifest_sha256"] != manifest_sha256:
        raise FirewallError("sealed outcome COMPLETE does not match ARTIFACTS.json")
    records_raw = manifest["artifacts"]
    if not isinstance(records_raw, list):
        raise FirewallError("sealed outcome artifact manifest must contain an array")
    records: dict[str, dict[str, Any]] = {}
    for index, item in enumerate(records_raw):
        record = _strict_mapping(
            item,
            fields=_ARTIFACT_RECORD_FIELDS,
            label=f"sealed outcome artifacts[{index}]",
        )
        record_path = _nonempty_string(
            record["path"], f"sealed outcome artifacts[{index}].path"
        )
        _sha256(record["sha256"], f"sealed outcome artifacts[{index}].sha256")
        if (
            isinstance(record["size_bytes"], bool)
            or not isinstance(record["size_bytes"], int)
            or record["size_bytes"] < 0
        ):
            raise FirewallError(
                f"sealed outcome artifacts[{index}].size_bytes must be nonnegative"
            )
        if record_path in records:
            raise FirewallError(f"sealed outcome manifest repeats {record_path}")
        records[record_path] = record
    if complete["artifact_count"] != len(records):
        raise FirewallError("sealed outcome completion artifact count mismatch")

    document_path = root / "sealed_outcome_bundle.json"
    document_record = records.get("sealed_outcome_bundle.json")
    if document_record is None:
        raise FirewallError("sealed outcome manifest omits sealed_outcome_bundle.json")
    try:
        document_changed = (
            document_path.stat().st_size != document_record["size_bytes"]
            or sha256_file(document_path) != document_record["sha256"]
        )
    except OSError as error:
        raise FirewallError(f"sealed outcome bundle document is unreadable: {error}") from error
    if document_changed:
        raise FirewallError("sealed outcome bundle document changed")
    bundle = SealedOutcomeBundle.from_dict(
        _read_json_object(document_path, "SealedOutcomeBundle")
    )
    metadata = manifest["metadata"]
    if not isinstance(metadata, Mapping):
        raise FirewallError("sealed outcome manifest metadata must be an object")
    if metadata.get("artifact_class") != "sealed_outcome_bundle":
        raise FirewallError("artifact tree is not classified as a sealed outcome bundle")
    if metadata.get("bundle_sha256") != bundle.bundle_sha256:
        raise FirewallError("sealed outcome manifest and bundle identity disagree")

    expected_paths = {"sealed_outcome_bundle.json"}
    for artifact in bundle.artifacts:
        configured = Path(artifact.path)
        if configured.is_absolute() or artifact.path == "sealed_outcome_bundle.json":
            raise FirewallError("sealed outcome artifacts must be relative data paths")
        try:
            resolved = _configured_path(
                root / configured, "sealed outcome artifact"
            ).resolve(strict=True)
            resolved.relative_to(root)
        except (OSError, ValueError) as error:
            raise FirewallError(
                f"sealed outcome artifact escapes or is missing: {artifact.path}"
            ) from error
        if not resolved.is_file():
            raise FirewallError(f"sealed outcome artifact is not a file: {artifact.path}")
        if resolved.stat().st_size != artifact.size_bytes:
            raise FirewallError(
                f"sealed outcome artifact size changed: {artifact.path}"
            )
        record = records.get(artifact.path)
        if record is None or (
            record["sha256"] != artifact.sha256
            or record["size_bytes"] != artifact.size_bytes
        ):
            raise FirewallError(
                f"sealed outcome manifest does not bind artifact {artifact.path}"
            )
        expected_paths.add(artifact.path)
    if set(records) != expected_paths:
        raise FirewallError("sealed outcome manifest contains undeclared artifacts")
    observed_paths = {
        member.relative_to(root).as_posix()
        for member in root.rglob("*")
        if member.is_file() and member.name not in {"ARTIFACTS.json", "COMPLETE"}
    }
    if observed_paths != expected_paths:
        raise FirewallError("sealed outcome directory inventory changed")
    return root, bundle, manifest_sha256


def verify_sealed_outcome_bundle(path: str | Path) -> SealedOutcomeBundle:
    root, bundle, _ = _preflight_sealed_outcome_bundle(path)
    try:
        verify_frozen_tree(root)
        bundle.validate_artifacts(root)
    except (ArtifactError, ContractError, OSError, ValueError, KeyError, TypeError) as error:
        raise FirewallError(f"sealed outcome artifact verification failed: {error}") from error
    return bundle


@dataclass(frozen=True, slots=True)
class OutcomeConsumption:
    consumption_sha256: str
    selection_lock_id: str
    selection_lock_manifest_sha256: str
    outcome_bundle_id: str
    outcome_bundle_sha256: str
    outcome_manifest_sha256: str
    task_bindings: tuple[Mapping[str, Any], ...]
    power_bindings: tuple[Mapping[str, Any], ...]
    terminal_policy: str
    marker_path: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": OUTCOME_CONSUMPTION_SCHEMA_VERSION,
            "consumption_sha256": self.consumption_sha256,
            "selection_lock_id": self.selection_lock_id,
            "selection_lock_manifest_sha256": self.selection_lock_manifest_sha256,
            "outcome_bundle_id": self.outcome_bundle_id,
            "outcome_bundle_sha256": self.outcome_bundle_sha256,
            "outcome_manifest_sha256": self.outcome_manifest_sha256,
            "task_bindings": [dict(item) for item in self.task_bindings],
            "power_bindings": [dict(item) for item in self.power_bindings],
            "terminal_policy": self.terminal_policy,
            "marker_path": self.marker_path,
        }


_OUTCOME_POWER_BINDING_FIELDS = frozenset(
    {
        "task_id",
        "power_decision_sha256",
        "power_result_sha256",
        "prediction_set_sha256",
        "prediction_bindings",
        "n_join_units",
        "n_join_rows",
        "unit_set_sha256",
        "row_set_sha256",
    }
)
_OUTCOME_CONSUMPTION_FIELDS = frozenset(
    {
        "schema_version",
        "consumption_sha256",
        "selection_lock_id",
        "selection_lock_manifest_sha256",
        "outcome_bundle_id",
        "outcome_bundle_sha256",
        "outcome_manifest_sha256",
        "task_bindings",
        "power_bindings",
        "terminal_policy",
    }
)


def _parse_outcome_power_binding(value: object, label: str) -> dict[str, Any]:
    raw = _strict_mapping(value, fields=_OUTCOME_POWER_BINDING_FIELDS, label=label)
    predictions_raw = raw["prediction_bindings"]
    if not isinstance(predictions_raw, (list, tuple)) or not predictions_raw:
        raise FirewallError(f"{label}.prediction_bindings must not be empty")
    predictions = tuple(
        _parse_prediction_binding(item, f"{label}.prediction_bindings[{index}]")
        for index, item in enumerate(predictions_raw)
    )
    run_ids = [item["run_id"] for item in predictions]
    if run_ids != sorted(run_ids) or len(run_ids) != len(set(run_ids)):
        raise FirewallError(
            f"{label}.prediction_bindings must have unique sorted run IDs"
        )
    counts: dict[str, int] = {}
    for field_name in ("n_join_units", "n_join_rows"):
        item = raw[field_name]
        if isinstance(item, bool) or not isinstance(item, int) or item < 1:
            raise FirewallError(f"{label}.{field_name} must be a positive integer")
        counts[field_name] = item
    if counts["n_join_units"] > counts["n_join_rows"]:
        raise FirewallError(f"{label}.n_join_units cannot exceed n_join_rows")
    return {
        "task_id": _nonempty_string(raw["task_id"], f"{label}.task_id"),
        "power_decision_sha256": _sha256(
            raw["power_decision_sha256"], f"{label}.power_decision_sha256"
        ),
        "power_result_sha256": _sha256(
            raw["power_result_sha256"], f"{label}.power_result_sha256"
        ),
        "prediction_set_sha256": _sha256(
            raw["prediction_set_sha256"], f"{label}.prediction_set_sha256"
        ),
        "prediction_bindings": [dict(item) for item in predictions],
        **counts,
        "unit_set_sha256": _sha256(
            raw["unit_set_sha256"], f"{label}.unit_set_sha256"
        ),
        "row_set_sha256": _sha256(
            raw["row_set_sha256"], f"{label}.row_set_sha256"
        ),
    }


def verify_outcome_consumption(path: str | Path) -> OutcomeConsumption:
    """Verify the immutable one-time marker without reopening sealed outcomes."""

    configured = _configured_path(path, "outcome-consumption marker")
    try:
        marker = configured.resolve(strict=True)
    except OSError as error:
        raise FirewallError(f"outcome-consumption marker is unavailable: {error}") from error
    if not marker.is_file():
        raise FirewallError("outcome-consumption marker must be a regular file")
    raw = _strict_mapping(
        _read_json_object(marker, "OutcomeConsumption"),
        fields=_OUTCOME_CONSUMPTION_FIELDS,
        label="OutcomeConsumption",
    )
    if raw["schema_version"] != OUTCOME_CONSUMPTION_SCHEMA_VERSION:
        raise FirewallError(
            f"unsupported OutcomeConsumption schema_version: {raw['schema_version']}"
        )
    task_bindings_raw = raw["task_bindings"]
    if not isinstance(task_bindings_raw, (list, tuple)) or not task_bindings_raw:
        raise FirewallError("OutcomeConsumption.task_bindings must not be empty")
    task_bindings = tuple(
        _parse_outcome_task_binding(item, f"task_bindings[{index}]")
        for index, item in enumerate(task_bindings_raw)
    )
    power_bindings_raw = raw["power_bindings"]
    if not isinstance(power_bindings_raw, (list, tuple)) or not power_bindings_raw:
        raise FirewallError("OutcomeConsumption.power_bindings must not be empty")
    power_bindings = tuple(
        _parse_outcome_power_binding(item, f"power_bindings[{index}]")
        for index, item in enumerate(power_bindings_raw)
    )
    task_ids = [str(item["task_id"]) for item in task_bindings]
    power_task_ids = [str(item["task_id"]) for item in power_bindings]
    if (
        task_ids != sorted(task_ids)
        or len(task_ids) != len(set(task_ids))
        or power_task_ids != task_ids
    ):
        raise FirewallError(
            "OutcomeConsumption task and power bindings must cover the same "
            "unique sorted task IDs"
        )
    terminal_policy = _nonempty_string(
        raw["terminal_policy"], "OutcomeConsumption.terminal_policy"
    )
    if terminal_policy != TERMINAL_OUTCOME_POLICY:
        raise FirewallError("OutcomeConsumption terminal policy changed")
    outcome_bundle_id = _nonempty_string(
        raw["outcome_bundle_id"], "OutcomeConsumption.outcome_bundle_id"
    )
    expected_name = (
        canonical_hash({"sealed_outcome_bundle_id": outcome_bundle_id})
        + ".consumed.json"
    )
    if marker.name != expected_name:
        raise FirewallError("OutcomeConsumption marker filename is not canonical")
    identity = {
        "schema_version": OUTCOME_CONSUMPTION_SCHEMA_VERSION,
        "selection_lock_id": _sha256(
            raw["selection_lock_id"], "OutcomeConsumption.selection_lock_id"
        ),
        "selection_lock_manifest_sha256": _sha256(
            raw["selection_lock_manifest_sha256"],
            "OutcomeConsumption.selection_lock_manifest_sha256",
        ),
        "outcome_bundle_id": outcome_bundle_id,
        "outcome_bundle_sha256": _sha256(
            raw["outcome_bundle_sha256"], "OutcomeConsumption.outcome_bundle_sha256"
        ),
        "outcome_manifest_sha256": _sha256(
            raw["outcome_manifest_sha256"],
            "OutcomeConsumption.outcome_manifest_sha256",
        ),
        "task_bindings": [dict(item) for item in task_bindings],
        "power_bindings": [dict(item) for item in power_bindings],
        "terminal_policy": terminal_policy,
    }
    consumption_sha256 = _sha256(
        raw["consumption_sha256"], "OutcomeConsumption.consumption_sha256"
    )
    if canonical_hash(identity) != consumption_sha256:
        raise FirewallError("OutcomeConsumption identity hash mismatch")
    return OutcomeConsumption(
        consumption_sha256=consumption_sha256,
        selection_lock_id=identity["selection_lock_id"],
        selection_lock_manifest_sha256=identity[
            "selection_lock_manifest_sha256"
        ],
        outcome_bundle_id=outcome_bundle_id,
        outcome_bundle_sha256=identity["outcome_bundle_sha256"],
        outcome_manifest_sha256=identity["outcome_manifest_sha256"],
        task_bindings=task_bindings,
        power_bindings=power_bindings,
        terminal_policy=terminal_policy,
        marker_path=marker.as_posix(),
    )


def _expected_joint_bindings(
    lock: SelectionLock, bundle: SealedOutcomeBundle
) -> dict[str, tuple[str, ...]]:
    bundle_datasets = {
        dataset_id
        for binding in bundle.task_bindings
        for dataset_id in binding["dataset_ids"]
    }
    expected: dict[str, tuple[str, ...]] = {}
    for decision in lock.task_decisions:
        sealed_ids = tuple(decision.get("sealed_dataset_ids", ()))
        if bundle_datasets.intersection(sealed_ids):
            expected[str(decision["task_id"])] = sealed_ids
    return expected


def _trees_overlap(first: Path, second: Path) -> bool:
    try:
        first.relative_to(second)
    except ValueError:
        pass
    else:
        return True
    try:
        second.relative_to(first)
    except ValueError:
        return False
    return True


def _artifact_rows(
    path: Path,
    artifact: ArtifactRef,
    *,
    expected_tsv_fields: tuple[str, ...] | None = None,
) -> list[Mapping[str, Any]]:
    try:
        if artifact.media_type == "application/json":
            payload = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(payload, Mapping) or set(payload) != {"rows"}:
                raise FirewallError(
                    f"join artifact {path} must contain exactly one rows array"
                )
            rows = payload["rows"]
            if not isinstance(rows, list):
                raise FirewallError(f"join artifact {path} rows must be an array")
            if any(not isinstance(row, Mapping) for row in rows):
                raise FirewallError(f"join artifact {path} rows must be objects")
            return [dict(row) for row in rows]
        if artifact.media_type == "text/tab-separated-values":
            with path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if (
                    reader.fieldnames is None
                    or len(reader.fieldnames) != len(set(reader.fieldnames))
                    or any(not field for field in reader.fieldnames)
                ):
                    raise FirewallError(
                        f"join artifact {path} has no unique non-empty header"
                    )
                if (
                    expected_tsv_fields is not None
                    and tuple(reader.fieldnames) != expected_tsv_fields
                ):
                    raise FirewallError(
                        f"join artifact {path} header must be exactly "
                        + "\t".join(expected_tsv_fields)
                    )
                return [dict(row) for row in reader]
    except FirewallError:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, csv.Error) as error:
        raise FirewallError(f"cannot read join artifact {path}: {error}") from error
    raise FirewallError(
        "exact outcome joining supports only application/json or "
        "text/tab-separated-values artifacts"
    )


def _join_identity(
    *,
    path: Path,
    artifact: ArtifactRef,
    row_id_field: str,
    unit_id_field: str,
    expected_rows: int,
    expected_tsv_fields: tuple[str, ...] | None = None,
) -> dict[str, Any]:
    rows = _artifact_rows(
        path,
        artifact,
        expected_tsv_fields=expected_tsv_fields,
    )
    if len(rows) != expected_rows:
        raise FirewallError(
            f"join artifact {path} has {len(rows)} rows; expected {expected_rows}"
        )
    row_ids: list[str] = []
    unit_ids: list[str] = []
    for index, row in enumerate(rows):
        if row_id_field not in row or unit_id_field not in row:
            raise FirewallError(
                f"join artifact {path} row {index} lacks its row or unit identifier"
            )
        row_id = row[row_id_field]
        unit_id = row[unit_id_field]
        if (
            not isinstance(row_id, str)
            or not row_id.strip()
            or "\x00" in row_id
            or not isinstance(unit_id, str)
            or not unit_id.strip()
            or "\x00" in unit_id
        ):
            raise FirewallError(
                f"join artifact {path} row {index} has an invalid row or unit identifier"
            )
        if row_id_field.endswith("_hash") and not is_sha256(row_id):
            raise FirewallError(
                f"join artifact {path} row {index} has a non-SHA-256 row identifier"
            )
        if unit_id_field.endswith("_hash") and not is_sha256(unit_id):
            raise FirewallError(
                f"join artifact {path} row {index} has a non-SHA-256 unit identifier"
            )
        row_ids.append(row_id)
        unit_ids.append(unit_id)
    if len(set(row_ids)) != len(row_ids):
        raise FirewallError(
            f"join artifact {path} must contain unique row identifiers"
        )
    sorted_units = sorted(set(unit_ids))
    sorted_rows = sorted(zip(row_ids, unit_ids, strict=True))
    return {
        "n_join_units": len(sorted_units),
        "n_join_rows": len(rows),
        "unit_set_sha256": canonical_hash(
            {"unit_id_field": unit_id_field, "unit_ids": sorted_units}
        ),
        "row_set_sha256": canonical_hash(
            {
                "row_id_field": row_id_field,
                "unit_id_field": unit_id_field,
                "row_unit_pairs": [
                    {"row_id": row_id, "unit_id": unit_id}
                    for row_id, unit_id in sorted_rows
                ],
            }
        ),
    }


def _prediction_join_identity(state: Path, power: PowerDecision) -> dict[str, Any]:
    expected_fields = _PREDICTION_TABLE_FIELDS.get(power.task_id)
    if expected_fields is None:
        raise FirewallError(
            f"task {power.task_id} has no registered sealed prediction-table schema"
        )
    if (power.row_id_field, power.unit_id_field) != expected_fields[:2]:
        raise FirewallError(
            f"PowerDecision {power.task_id} row/unit fields are not task-native"
        )
    if power.table_schema_sha256 != canonical_hash(
        {"format": "tsv", "fields": list(expected_fields)}
    ):
        raise FirewallError(
            f"PowerDecision {power.task_id} does not bind its task-native table schema"
        )
    observed: dict[str, Any] | None = None
    for binding in power.prediction_bindings:
        commit = verify_prediction_commit(
            _commit_dir(state, binding["prediction_commit_sha256"]),
            reverify_sources=True,
        )
        if (
            commit.standardized_table.media_type != "text/tab-separated-values"
            or commit.row_ids.media_type != "text/tab-separated-values"
        ):
            raise FirewallError(
                "exact prediction joining requires canonical TSV standardized tables "
                "and row-ID inventories"
            )
        try:
            table_path = commit.standardized_table.validate(
                Path(commit.prediction_bundle_path).parent,
                require_relative=True,
            )
            row_ids_path = commit.row_ids.validate(
                Path(commit.prediction_bundle_path).parent,
                require_relative=True,
            )
        except ContractError as error:
            raise FirewallError(
                f"prediction join artifact failed validation: {error}"
            ) from error
        table_identity = _join_identity(
            path=table_path,
            artifact=commit.standardized_table,
            row_id_field=power.row_id_field,
            unit_id_field=power.unit_id_field,
            expected_rows=power.n_predictions,
            expected_tsv_fields=expected_fields,
        )
        inventory_rows = _artifact_rows(
            row_ids_path,
            commit.row_ids,
            expected_tsv_fields=(power.row_id_field, power.unit_id_field),
        )
        inventory_order = [
            (str(row[power.row_id_field]), str(row[power.unit_id_field]))
            for row in inventory_rows
        ]
        if inventory_order != sorted(inventory_order):
            raise FirewallError(
                f"PredictionCommit {commit.run_id} row inventory is not canonical"
            )
        inventory_identity = _join_identity(
            path=row_ids_path,
            artifact=commit.row_ids,
            row_id_field=power.row_id_field,
            unit_id_field=power.unit_id_field,
            expected_rows=power.n_predictions,
        )
        if inventory_identity != table_identity:
            raise FirewallError(
                f"PredictionCommit {commit.run_id} row inventory differs from its table"
            )
        identity = table_identity
        if observed is None:
            observed = identity
        elif observed != identity:
            raise FirewallError(
                f"selected and baseline prediction row sets differ for {power.task_id}"
            )
    if observed is None:
        raise FirewallError(f"PowerDecision {power.task_id} has no predictions")
    return observed


def _outcome_artifact_for_task(
    bundle: SealedOutcomeBundle, task_id: str
) -> ArtifactRef:
    matches = [
        artifact for artifact in bundle.artifacts if artifact.role == f"outcome:{task_id}"
    ]
    if len(matches) != 1:
        raise FirewallError(
            f"sealed outcome task {task_id} needs exactly one role=outcome:{task_id} artifact"
        )
    return matches[0]


def authorize_one_time_outcome_join(
    *,
    evaluator_state_dir: str | Path,
    power_decision_sha256s: Iterable[str],
    sealed_outcome_bundle_dir: str | Path,
    consumption_ledger_dir: str | Path,
) -> OutcomeConsumption:
    """Atomically consume a whole sealed bundle after every task is ready."""

    state = _configured_path(evaluator_state_dir, "evaluator-state root").resolve()
    power_ids = tuple(power_decision_sha256s)
    if not power_ids or len(set(power_ids)) != len(power_ids):
        raise FirewallError(
            "power_decision_sha256s must be a non-empty set without duplicates"
        )
    powers = tuple(
        verify_power_decision(_power_dir(state, power_id), reverify_sources=True)
        for power_id in power_ids
    )
    if len({power.task_id for power in powers}) != len(powers):
        raise FirewallError("outcome authorization has duplicate power decisions per task")
    if len({power.selection_lock_id for power in powers}) != 1:
        raise FirewallError("joint outcome tasks must share one SelectionLock")
    if len({power.selection_lock_manifest_sha256 for power in powers}) != 1:
        raise FirewallError("joint outcome tasks use different SelectionLock manifests")

    first_binding = powers[0].prediction_bindings[0]
    first_commit = verify_prediction_commit(
        _commit_dir(state, first_binding["prediction_commit_sha256"]),
        reverify_sources=True,
    )
    lock = _verify_selection_dir(Path(first_commit.selection_lock_dir))
    outcome_root, outcome_bundle, outcome_manifest_sha256 = (
        _preflight_sealed_outcome_bundle(sealed_outcome_bundle_dir)
    )
    if _trees_overlap(outcome_root, state):
        raise FirewallError("sealed outcome bundle must be outside evaluator state")
    outcome_bindings = {
        str(binding["task_id"]): tuple(binding["dataset_ids"])
        for binding in outcome_bundle.task_bindings
    }
    expected_bindings = _expected_joint_bindings(lock, outcome_bundle)
    if outcome_bindings != expected_bindings:
        raise FirewallError(
            "sealed outcome bundle must jointly cover every locked task sharing its datasets"
        )
    powers_by_task = {power.task_id: power for power in powers}
    if set(powers_by_task) != set(outcome_bindings):
        raise FirewallError(
            "every task in the joint outcome bundle needs one frozen power decision"
        )
    join_identities: dict[str, dict[str, Any]] = {}
    for task_id, dataset_ids in outcome_bindings.items():
        power = powers_by_task[task_id]
        if power.dataset_ids != dataset_ids:
            raise FirewallError(
                f"PowerDecision datasets do not match outcome binding for {task_id}"
            )
        if outcome_bundle.unit_id_fields[task_id] != power.unit_id_field:
            raise FirewallError(
                f"prediction and outcome unit identifiers do not match for {task_id}"
            )
        if power.power_result.get("passed") is not True:
            raise FirewallError(
                f"underpowered task {task_id} must be downgraded before outcomes are opened"
            )
        _outcome_artifact_for_task(outcome_bundle, task_id)
        join_identities[task_id] = _prediction_join_identity(state, power)

    power_bindings = tuple(
        {
            "task_id": power.task_id,
            "power_decision_sha256": power.power_decision_sha256,
            "power_result_sha256": power.power_result_sha256,
            "prediction_set_sha256": power.prediction_set_sha256,
            "prediction_bindings": [dict(item) for item in power.prediction_bindings],
            **join_identities[power.task_id],
        }
        for power in sorted(powers, key=lambda item: item.task_id)
    )
    identity = {
        "schema_version": OUTCOME_CONSUMPTION_SCHEMA_VERSION,
        "selection_lock_id": lock.lock_id,
        "selection_lock_manifest_sha256": powers[0].selection_lock_manifest_sha256,
        "outcome_bundle_id": outcome_bundle.bundle_id,
        "outcome_bundle_sha256": outcome_bundle.bundle_sha256,
        "outcome_manifest_sha256": outcome_manifest_sha256,
        "task_bindings": [dict(item) for item in outcome_bundle.task_bindings],
        "power_bindings": [dict(item) for item in power_bindings],
        "terminal_policy": TERMINAL_OUTCOME_POLICY,
    }
    consumption_sha256 = canonical_hash(identity)
    ledger = _configured_path(
        consumption_ledger_dir, "consumption-ledger root"
    ).resolve()
    if not ledger.is_dir():
        raise FirewallError(
            "consumption_ledger_dir must be a pre-existing independent directory"
        )
    if _trees_overlap(ledger, state):
        raise FirewallError("consumption ledger must be outside evaluator state")
    if _trees_overlap(ledger, outcome_root):
        raise FirewallError("consumption ledger must be outside the sealed outcome bundle")
    bundle_key = canonical_hash({"sealed_outcome_bundle_id": outcome_bundle.bundle_id})
    marker = ledger / f"{bundle_key}.consumed.json"
    marker_payload = {
        **identity,
        "consumption_sha256": consumption_sha256,
    }
    try:
        write_json_exclusive(marker, marker_payload, mode=0o440)
    except ArtifactError as error:
        raise FirewallError(
            f"sealed outcome bundle {outcome_bundle.bundle_id} was already consumed"
        ) from error

    # The atomic marker is deliberately written before any outcome artifact is
    # opened or hashed. A corrupt bundle therefore remains consumed and cannot
    # be repaired and replayed after partial exposure.
    for power in powers:
        verify_power_decision(
            _power_dir(state, power.power_decision_sha256),
            reverify_sources=True,
        )
    verified_bundle = verify_sealed_outcome_bundle(outcome_root)
    if verified_bundle.bundle_sha256 != outcome_bundle.bundle_sha256:
        raise FirewallError("sealed outcome bundle changed during atomic consumption")
    for task_id, expected_identity in join_identities.items():
        artifact = _outcome_artifact_for_task(verified_bundle, task_id)
        try:
            outcome_path = artifact.validate(outcome_root, require_relative=True)
        except ContractError as error:
            raise FirewallError(
                f"sealed outcome join artifact failed validation: {error}"
            ) from error
        observed_identity = _join_identity(
            path=outcome_path,
            artifact=artifact,
            row_id_field=powers_by_task[task_id].row_id_field,
            unit_id_field=verified_bundle.unit_id_fields[task_id],
            expected_rows=powers_by_task[task_id].n_predictions,
            expected_tsv_fields=_OUTCOME_TABLE_FIELDS[task_id],
        )
        if observed_identity != expected_identity:
            raise FirewallError(
                f"sealed outcome unit or row set differs from predictions for {task_id}"
            )
    return verify_outcome_consumption(marker)


__all__ = [
    "CampaignPhase",
    "DEFAULT_RESOURCE_AUTHORITIES",
    "DEVELOPMENT_POWER_EVIDENCE_SCHEMA_VERSION",
    "DevelopmentPowerEvidence",
    "EMPIRICAL_POWER_METHOD_ID",
    "FirewallError",
    "OutcomeConsumption",
    "PowerDecision",
    "PredictionCommit",
    "SEALED_OUTCOME_BUNDLE_SCHEMA_VERSION",
    "SealedArtifact",
    "SealedOutcomeBundle",
    "UNIT_TEST_RESAMPLING_POLICY",
    "assert_dataset_access",
    "assert_resource_unchanged",
    "authorize_one_time_outcome_join",
    "commit_predictions",
    "freeze_development_power_evidence",
    "freeze_power_decision",
    "resource_snapshot",
    "verify_development_power_evidence",
    "verify_outcome_consumption",
    "verify_power_decision",
    "verify_prediction_commit",
    "verify_sealed_outcome_bundle",
]
