"""Freeze reproducible, reviewable benchmark campaign plans."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
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
    write_json_exclusive,
    write_text_exclusive,
)
from .contracts import (
    ArtifactRef,
    ContractError,
    DatasetManifest,
    ModelExecutionContract,
    ModelManifest,
    RunSpec,
    TaskSpec,
)
from .firewall import resource_snapshot
from .hashing import canonicalize
from .registry import Registry, RegistryError
from .slurm import ResourceProfile, load_resource_profiles, render_sbatch, resource_totals


class PlanningError(RuntimeError):
    """Raised when a campaign cannot be safely frozen."""


def _reject_configured_symlinks(path: str | Path, label: str) -> Path:
    try:
        return reject_symlink_components(Path(path), label=label)
    except ArtifactError as error:
        raise PlanningError(f"{label} path rejected: {error}") from error


_PAIRING_LEVELS = frozenset(
    {
        "same_molecule",
        "same_cell",
        "same_nucleus",
        "same_section",
        "adjacent_section",
        "same_sample_different_aliquot",
        "same_donor_different_tissue",
        "same_study_unpaired",
    }
)
_MISSING_STATES = frozenset(
    {
        "observed",
        "structurally_missing",
        "not_applicable",
        "below_qc",
        "unavailable_permission",
        "join_unresolved",
        "withheld_sealed",
        "derivable_not_processed",
    }
)
_IDENTIFIER = re.compile(r"^[a-z0-9][a-z0-9._-]{0,127}$")
_SOURCE_OUTPUT_ROOTS = frozenset(
    {"candidates", "executions", "releases", "results", "selections"}
)
_SOURCE_IGNORED_PARTS = frozenset({".git", ".pytest_cache", "__pycache__"})
_UNRESOLVED_SENTINELS = frozenset({"", "UNKNOWN", "UNRESOLVED", "TBD", "LATEST"})
_UNRESOLVED_PROSE_KEYS = frozenset({"blockers", "description", "notes"})
_SCIENTIFIC_WAVE_ACTIONS: Mapping[str, tuple[str, ...]] = {
    "smoke": ("prepare", "fit", "predict"),
    "frozen_screen": ("prepare", "fit", "predict"),
    "full_specialist_screen": ("prepare", "fit", "predict"),
    "adaptation": ("prepare", "fit", "predict"),
    "error_audit": ("predict",),
    "conditional_model": ("prepare", "fit", "predict"),
    "prediction_first_stress": ("prepare", "predict"),
    "sealed_inference": ("predict",),
    "terminal_reporting": ("export",),
    "fnih_activation_audit": ("probe",),
}


@dataclass(frozen=True, slots=True)
class RegistryDocument:
    path: Path
    payload: Mapping[str, Any]
    sha256: str


def _load_toml(path: Path) -> Mapping[str, Any]:
    try:
        with path.open("rb") as handle:
            payload = tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise PlanningError(f"cannot load TOML document {path}: {error}") from error
    if not isinstance(payload, Mapping):
        raise PlanningError(f"TOML root must be a table: {path}")
    return payload


def _documents(directory: Path) -> tuple[RegistryDocument, ...]:
    if not directory.is_dir():
        raise PlanningError(f"registry directory does not exist: {directory}")
    return tuple(
        RegistryDocument(path, _load_toml(path), sha256_file(path))
        for path in sorted(directory.glob("*.toml"))
    )


def _stable_file_record(path: Path, *, relative_to: Path) -> dict[str, Any]:
    """Hash one regular file and reject a concurrent replacement."""

    if path.is_symlink():
        raise PlanningError(f"provenance trees may not contain symlinks: {path}")
    if not path.is_file():
        raise PlanningError(f"provenance tree member is not a regular file: {path}")
    try:
        before = path.stat(follow_symlinks=False)
        digest = sha256_file(path)
        after = path.stat(follow_symlinks=False)
    except OSError as error:
        raise PlanningError(f"cannot hash provenance tree member {path}: {error}") from error
    identity_before = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    )
    identity_after = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    )
    if identity_before != identity_after or path.is_symlink():
        raise PlanningError(f"provenance tree member changed while hashing: {path}")
    return {
        "path": path.relative_to(relative_to).as_posix(),
        "sha256": digest,
        "size_bytes": after.st_size,
    }


def _tree_file_records(
    root: Path,
    *,
    suffix: str | None = None,
    excluded_parts: frozenset[str] = frozenset(),
    excluded_roots: frozenset[str] = frozenset(),
) -> list[dict[str, Any]]:
    """Return a deterministic regular-file manifest without following links."""

    records: list[dict[str, Any]] = []
    pending = [root]
    while pending:
        directory = pending.pop()
        try:
            members = sorted(directory.iterdir(), key=lambda path: path.name)
        except OSError as error:
            raise PlanningError(f"cannot enumerate provenance directory {directory}: {error}") from error
        for path in members:
            relative = path.relative_to(root)
            if (
                (relative.parts and relative.parts[0] in excluded_roots)
                or any(part in excluded_parts for part in relative.parts)
            ):
                continue
            if path.is_symlink():
                raise PlanningError(f"provenance trees may not contain symlinks: {path}")
            if path.is_dir():
                pending.append(path)
                continue
            if not path.is_file():
                raise PlanningError(f"provenance tree contains a non-regular member: {path}")
            if suffix is None or path.suffix == suffix:
                records.append(_stable_file_record(path, relative_to=root))
    return sorted(records, key=lambda record: str(record["path"]))


def _config_registry_snapshot(config_root: str | Path) -> dict[str, Any]:
    """Bind every regular configuration file by relative path, digest, and size."""

    source = _reject_configured_symlinks(config_root, "configuration root")
    try:
        root = source.resolve(strict=True)
    except OSError as error:
        raise PlanningError(f"configuration root cannot be resolved: {source}") from error
    if not root.is_dir():
        raise PlanningError(f"configuration root is not a directory: {root}")
    files = _tree_file_records(root)
    if not files:
        raise PlanningError(f"configuration tree has no registry documents: {root}")
    return {
        "schema_version": "masld-bench-registry-file-snapshot-v1",
        "files": files,
    }


def _require_identifier(value: Any, label: str) -> str:
    identifier = str(value or "")
    if not _IDENTIFIER.fullmatch(identifier):
        raise PlanningError(
            f"{label} must match {_IDENTIFIER.pattern!r}; path syntax is forbidden"
        )
    return identifier


def _runtimes(config_root: Path) -> dict[str, dict[str, Any]]:
    runtimes: dict[str, dict[str, Any]] = {}
    for document in _documents(config_root / "runtimes"):
        runtime = dict(document.payload)
        runtime_id = _require_identifier(runtime.get("runtime_id"), "runtime_id")
        if runtime_id in runtimes:
            raise PlanningError(f"duplicate runtime identifier: {runtime_id}")
        prefix = runtime.get("command_prefix")
        if not isinstance(prefix, list) or not prefix or any(
            not isinstance(item, str) or not item.strip() for item in prefix
        ):
            runtime["admission_blocking"] = True
            blockers = list(runtime.get("blockers", []))
            blockers.append("runtime command_prefix is unresolved")
            runtime["blockers"] = blockers
        runtime["registry_path"] = document.path.relative_to(config_root).as_posix()
        runtime["registry_sha256"] = document.sha256
        runtimes[runtime_id] = runtime
    return runtimes


def _models(config_root: Path) -> tuple[dict[str, Any], ...]:
    models: list[dict[str, Any]] = []
    for document in _documents(config_root / "models"):
        family_id = _require_identifier(document.payload.get("family_id"), "family_id")
        raw_models = document.payload.get("models")
        if not isinstance(raw_models, list) or not raw_models:
            raise PlanningError(f"model family has no [[models]] entries: {document.path}")
        for raw_model in raw_models:
            if not isinstance(raw_model, Mapping):
                raise PlanningError(f"invalid model entry in {document.path}")
            model = dict(raw_model)
            model["model_id"] = _require_identifier(model.get("model_id"), "model_id")
            model["family_id"] = family_id
            model["family_modality"] = document.payload.get("modality")
            model["registry_path"] = document.path.relative_to(config_root).as_posix()
            model["registry_sha256"] = document.sha256
            models.append(model)
    identifiers = [model["model_id"] for model in models]
    duplicates = sorted({item for item in identifiers if identifiers.count(item) > 1})
    if duplicates:
        raise PlanningError(f"duplicate model identifiers: {', '.join(duplicates)}")
    return tuple(sorted(models, key=lambda model: model["model_id"]))


def _tasks(config_root: Path) -> dict[str, dict[str, Any]]:
    tasks: dict[str, dict[str, Any]] = {}
    for document in _documents(config_root / "tasks"):
        task = dict(document.payload)
        task_id = _require_identifier(task.get("task_id"), "task_id")
        if task_id in tasks:
            raise PlanningError(f"duplicate task identifier: {task_id}")
        task["registry_path"] = document.path.relative_to(config_root).as_posix()
        task["registry_sha256"] = document.sha256
        tasks[task_id] = task
    return tasks


def _datasets(config_root: Path) -> dict[str, dict[str, Any]]:
    datasets: dict[str, dict[str, Any]] = {}
    for document in _documents(config_root / "datasets"):
        dataset = dict(document.payload)
        dataset_id = _require_identifier(dataset.get("dataset_id"), "dataset_id")
        if dataset_id in datasets:
            raise PlanningError(f"duplicate dataset identifier: {dataset_id}")
        dataset["registry_path"] = document.path.relative_to(config_root).as_posix()
        dataset["registry_sha256"] = document.sha256
        datasets[dataset_id] = dataset
    return datasets


def _splits(config_root: Path) -> dict[str, dict[str, Any]]:
    splits: dict[str, dict[str, Any]] = {}
    for document in _documents(config_root / "splits"):
        split = dict(document.payload)
        split_id = _require_identifier(split.get("split_id"), "split_id")
        if split_id in splits:
            raise PlanningError(f"duplicate split identifier: {split_id}")
        split["registry_path"] = document.path.relative_to(config_root).as_posix()
        split["registry_sha256"] = document.sha256
        splits[split_id] = split
    return splits


def _resource_document(config_root: Path) -> Mapping[str, Any]:
    path = config_root / "resources.toml"
    payload = dict(_load_toml(path))
    if "profiles" not in payload:
        nested = payload.get("cluster")
        if isinstance(nested, Mapping) and isinstance(nested.get("profiles"), Mapping):
            payload["profiles"] = nested["profiles"]
    return payload


def validate_registry_tree(config_root: str | Path) -> dict[str, Any]:
    """Validate cross-document admission invariants without opening any dataset."""

    registry_snapshot = _config_registry_snapshot(config_root)
    registry_snapshot_sha256 = canonical_hash(registry_snapshot)
    config = _reject_configured_symlinks(
        config_root, "configuration root"
    ).resolve(strict=True)
    from .registry import Registry

    strict_registry = Registry.load(config, validate_references=False)
    models = _models(config)
    tasks = _tasks(config)
    datasets = _datasets(config)
    profiles = load_resource_profiles(_resource_document(config))
    splits = _splits(config)
    split_ids = set(splits)
    runtimes = _runtimes(config)
    runtime_ids = set(runtimes)
    for runtime_id, runtime in runtimes.items():
        profile_id = runtime.get("resource_profile")
        if profile_id is not None and str(profile_id) not in profiles:
            raise PlanningError(f"runtime {runtime_id} names unknown profile {profile_id}")
    for dataset_id, dataset in datasets.items():
        pairing = {str(item) for item in dataset.get("pairing_levels", [])}
        invalid_pairing = sorted(pairing.difference(_PAIRING_LEVELS))
        if invalid_pairing:
            raise PlanningError(
                f"dataset {dataset_id} has invalid pairing levels: {', '.join(invalid_pairing)}"
            )
        missing_default = str(dataset.get("modality_status_default", ""))
        if missing_default not in _MISSING_STATES:
            raise PlanningError(
                f"dataset {dataset_id} has invalid modality_status_default: {missing_default}"
            )
        if dataset.get("role") in {
            "project_sealed_final",
            "conditional_bulk_seal",
            "conditional_future_seal",
            "intervention_stress_test",
            "withheld_sealed",
            "conditional_sealed",
            "conditional_future_holdout",
        } and bool(dataset.get("automatic_download", False)):
            raise PlanningError(f"sealed dataset {dataset_id} may not enable automatic download")
    datasets_casefold = {key.casefold(): value for key, value in datasets.items()}
    if "gse244832" in datasets_casefold and any(
        level in datasets_casefold["gse244832"].get("pairing_levels", [])
        for level in ("same_cell", "same_nucleus", "same_molecule")
    ):
        raise PlanningError("GSE244832 must not encode cell-to-cell or nucleus pairing")
    if "gse296875" in datasets_casefold and "same_nucleus" not in datasets_casefold["gse296875"].get(
        "pairing_levels", []
    ):
        raise PlanningError("GSE296875 must retain its same-nucleus RNA-ATAC topology")
    for task_id, task in tasks.items():
        split_id = str(task.get("split_id", ""))
        if split_id not in split_ids:
            raise PlanningError(f"task {task_id} names unknown split {split_id}")
        for dataset_id in _task_dataset_ids(task):
            if dataset_id not in datasets:
                raise PlanningError(f"task {task_id} names unknown dataset {dataset_id}")
        missingness_policy = str(task.get("missingness_policy", ""))
        if missingness_policy not in {"explicit_mask", "explicit_state_and_mask"}:
            raise PlanningError(
                f"task {task_id} must use explicit missing-state masks, not {missingness_policy!r}"
            )
    dispositions = []
    for model in models:
        model_contract = strict_registry.models[model["model_id"]]
        open_champion_decision = _open_champion_license_decision(
            model, model_contract.execution
        )
        disposition, blockers = _model_disposition(
            model, open_champion_decision=open_champion_decision
        )
        dispositions.append(
            {
                "model_id": model["model_id"],
                "disposition": disposition,
                "blockers": list(blockers),
                "open_champion_license_decision": open_champion_decision,
            }
        )
    return {
        "schema_version": "masld-bench-registry-validation-v1",
        "datasets": len(datasets),
        "model_families": len({model["family_id"] for model in models}),
        "models": len(models),
        "tasks": len(tasks),
        "splits": len(split_ids),
        "runtimes": len(runtime_ids),
        "resource_profiles": len(profiles),
        "registry_snapshot": registry_snapshot,
        "registry_snapshot_sha256": registry_snapshot_sha256,
        "core_registry_contract_sha256": strict_registry.snapshot_sha256,
        "model_dispositions": dispositions,
    }


def _git_output(args: list[str], cwd: Path) -> bytes | None:
    process = subprocess.run(
        args,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return process.stdout if process.returncode == 0 else None


def source_lock(package_root: str | Path) -> dict[str, Any]:
    """Record commit, dirty patch, and package tree as one source identity."""

    source = _reject_configured_symlinks(package_root, "package root")
    try:
        root = source.resolve(strict=True)
    except OSError as error:
        raise PlanningError(f"package root cannot be resolved: {source}") from error
    if not root.is_dir():
        raise PlanningError(f"package root is not a directory: {root}")
    repository_text = _git_output(["git", "rev-parse", "--show-toplevel"], root)
    if not repository_text:
        raise PlanningError("source locking requires a Git repository")
    repository = _reject_configured_symlinks(
        Path(repository_text.decode().strip()), "Git repository root"
    ).resolve(strict=True)
    try:
        relative = root.relative_to(repository).as_posix()
    except ValueError as error:
        raise PlanningError("package root is outside its reported Git repository") from error
    head_payload = _git_output(["git", "rev-parse", "HEAD"], repository)
    if not head_payload or not head_payload.decode().strip():
        raise PlanningError("source locking requires a resolved Git commit")
    head = head_payload.decode().strip()

    from hashlib import sha256

    dirty_components: list[dict[str, Any]] = []
    dirty_payload = b""
    for label, command in (
        ("unstaged", ["git", "diff", "--binary", "--", relative]),
        ("staged", ["git", "diff", "--cached", "--binary", "--", relative]),
    ):
        payload = _git_output(command, repository)
        if payload is None:
            raise PlanningError(f"cannot capture {label} Git patch for {root}")
        dirty_components.append(
            {
                "component": label,
                "sha256": sha256(payload).hexdigest(),
                "size_bytes": len(payload),
            }
        )
        dirty_payload += label.encode("utf-8") + b"\0" + payload + b"\0"

    records = _tree_file_records(
        root,
        excluded_parts=_SOURCE_IGNORED_PARTS,
        excluded_roots=_SOURCE_OUTPUT_ROOTS,
    )
    payload = {
        "schema_version": "masld-bench-source-lock-v2",
        "git_commit": head,
        "dirty_patch_sha256": sha256(dirty_payload).hexdigest(),
        "dirty_patch_components": dirty_components,
        "tree_sha256": canonical_hash(records),
        "tree_files": len(records),
        "tree_manifest": records,
    }
    result = dict(payload)
    result["source_lock_sha256"] = canonical_hash(payload)
    return result


_OPEN_CODE_LICENSES = frozenset(
    {
        "MIT",
        "BSD-2-Clause",
        "BSD-3-Clause",
        "Apache-2.0",
        "GPL-2.0-only",
        "GPL-2.0-or-later",
        "GPL-3.0-only",
        "GPL-3.0-or-later",
        "LGPL-3.0-only",
        "LGPL-3.0-or-later",
        "project_owned",
    }
)
_REDISTRIBUTABLE_WEIGHT_LICENSES = frozenset(
    {
        "MIT",
        "BSD-2-Clause",
        "BSD-3-Clause",
        "Apache-2.0",
        "CC-BY-4.0",
        "CC-BY-SA-4.0",
        "CC0-1.0",
        "GPL-2.0-only",
        "GPL-2.0-or-later",
        "GPL-3.0-only",
        "GPL-3.0-or-later",
        "LGPL-3.0-only",
        "LGPL-3.0-or-later",
        "project_owned",
        "not_applicable_classical",
    }
)
_OPEN_MANIFEST_LICENSES = frozenset(
    {
        "open",
        "open_source",
        "permissive",
        "redistributable",
        "project_owned",
        "MIT",
        "BSD-2-Clause",
        "BSD-3-Clause",
        "Apache-2.0",
        "CC-BY-4.0",
        "CC-BY-SA-4.0",
        "GPL-2.0",
        "GPL-2.0-only",
        "GPL-2.0-or-later",
        "GPL-3.0",
        "GPL-3.0-only",
        "GPL-3.0-or-later",
        "LGPL-3.0",
        "LGPL-3.0-only",
        "LGPL-3.0-or-later",
    }
)


def _open_champion_license_decision(
    model: Mapping[str, Any], execution: ModelExecutionContract
) -> dict[str, Any]:
    top_level = str(model.get("license_status", "UNRESOLVED"))
    exposure = str(model.get("exposure_status", "unknown"))
    if execution.ready and top_level not in (
        _OPEN_MANIFEST_LICENSES
        | {
            "research_only",
            "noncommercial",
            "restricted",
            "redistribution_prohibited",
            "source_terms_apply",
        }
    ):
        raise PlanningError(
            f"ready model {model.get('model_id')} has an unrecognized manifest license_status"
        )
    checks = {
        "manifest_license_open": top_level in _OPEN_MANIFEST_LICENSES,
        "code_use_open": execution.code_license in _OPEN_CODE_LICENSES,
        "weights_redistributable": (
            execution.weights_license in _REDISTRIBUTABLE_WEIGHT_LICENSES
        ),
        "derivative_weights_redistributable": (
            execution.derivative_weights_license
            in _REDISTRIBUTABLE_WEIGHT_LICENSES
        ),
        "exposure_eligible": exposure
        in {"clean_declared", "target_label_unexposed"},
    }
    reasons = sorted(key for key, passed in checks.items() if not passed)
    eligible = execution.ready and not reasons
    return {
        "schema_version": "masld-bench-open-champion-license-decision-v1",
        "audited": execution.ready,
        "model_license_status": top_level,
        "code_license": execution.code_license,
        "weights_license": execution.weights_license,
        "derivative_weights_license": execution.derivative_weights_license,
        "exposure_status": exposure,
        **checks,
        "eligible": eligible,
        "ineligibility_reasons": reasons if execution.ready else [
            "model_execution_contract_not_ready"
        ],
    }


def _model_disposition(
    model: Mapping[str, Any], *, open_champion_decision: Mapping[str, Any]
) -> tuple[str, tuple[str, ...]]:
    blockers = tuple(str(item) for item in model.get("blockers", []))
    status = str(model.get("status", "unknown"))
    blocking = bool(model.get("admission_blocking", False))
    required_exact = {
        "checkpoint_revision": model.get("checkpoint_revision"),
        "checkpoint_sha256": model.get("checkpoint_sha256"),
        "license_status": model.get("license_status"),
        "exposure_status": model.get("exposure_status"),
    }
    unresolved_sentinels = {"", "UNKNOWN", "UNRESOLVED", "TBD", "LATEST", "NONE"}
    missing = tuple(
        key
        for key, value in required_exact.items()
        if value is None or str(value).strip().upper() in unresolved_sentinels
    )
    if missing:
        blockers += tuple(f"missing_{key}" for key in missing)
        blocking = True
    implementation_type = str(model.get("implementation_type", ""))
    checkpoint_sha256 = str(model.get("checkpoint_sha256", ""))
    checkpoint_not_applicable = checkpoint_sha256.lower() in {
        "not_applicable",
        "not_applicable_pretraining",
        "trained_within_fold",
    }
    if (
        not checkpoint_not_applicable
        and not blocking
        and (
            len(checkpoint_sha256) != 64
            or any(character not in "0123456789abcdefABCDEF" for character in checkpoint_sha256)
        )
    ):
        blockers += ("checkpoint_sha256_not_exact",)
        blocking = True
    if checkpoint_not_applicable and implementation_type not in {
        "classical_baseline",
        "graph_baseline",
        "train_locally",
        "sequence_baseline",
        "fixture",
    }:
        blockers += ("pretrained_checkpoint_cannot_be_not_applicable",)
        blocking = True
    if status in {"blocked", "blocked_terms", "deferred"}:
        blocking = True
    if blocking:
        return "blocked", tuple(sorted(set(blockers or (f"status_{status}",))))
    champion_eligible = bool(open_champion_decision.get("eligible")) and status not in {
        "restricted_comparator"
    }
    return ("eligible_open" if champion_eligible else "comparator_only"), tuple(sorted(set(blockers)))


def _task_dataset_ids(
    task: Mapping[str, Any],
    partitions: Iterable[str] = (
        "datasets_train",
        "datasets_development",
        "datasets_sealed",
    ),
) -> tuple[str, ...]:
    result: list[str] = []
    for key in partitions:
        value = task.get(key, [])
        if isinstance(value, str):
            result.append(value)
        elif isinstance(value, Iterable):
            result.extend(str(item) for item in value)
    return tuple(dict.fromkeys(result))


def _run_dataset_ids(task: Mapping[str, Any], wave: str) -> tuple[str, ...]:
    if wave == "admission":
        return _task_dataset_ids(task)
    if wave == "fnih_activation_audit":
        return tuple(
            dataset_id
            for dataset_id in _task_dataset_ids(task, ("datasets_sealed",))
            if dataset_id == "fnih_86_liver"
        )
    if wave == "terminal_reporting":
        return ()
    if wave == "prediction_first_stress":
        return _task_dataset_ids(task, ("datasets_development",))
    if wave == "sealed_inference":
        return _task_dataset_ids(task, ("datasets_sealed",))
    return _task_dataset_ids(task, ("datasets_train", "datasets_development"))


def _dataset_access_metadata(
    *,
    task: TaskSpec,
    run_dataset_ids: tuple[str, ...],
    dataset_locks: Mapping[str, Mapping[str, Any]],
    wave: str,
) -> dict[str, Any]:
    run_set = set(run_dataset_ids)
    fit_dataset_ids = tuple(
        dataset_id for dataset_id in task.datasets_train if dataset_id in run_set
    )
    development_prediction_dataset_ids = tuple(
        dataset_id
        for dataset_id in task.datasets_development
        if dataset_id in run_set
        and dataset_locks[dataset_id].get("prediction_first_policy") is None
    )
    prediction_first_stress_dataset_ids = tuple(
        dataset_id
        for dataset_id in task.datasets_development
        if dataset_id in run_set
        and dataset_locks[dataset_id].get("prediction_first_policy") is not None
    )
    sealed_prediction_dataset_ids = tuple(
        dataset_id for dataset_id in task.datasets_sealed if dataset_id in run_set
    )
    routing: dict[str, dict[str, Any]] = {}
    for dataset_id in run_dataset_ids:
        lock = dataset_locks[dataset_id]
        if dataset_id in fit_dataset_ids:
            partition = "fit"
        elif dataset_id in development_prediction_dataset_ids:
            partition = "development_prediction"
        elif dataset_id in prediction_first_stress_dataset_ids:
            partition = "prediction_first_stress"
        elif dataset_id in sealed_prediction_dataset_ids:
            partition = "sealed_prediction"
        else:
            raise PlanningError(
                f"dataset {dataset_id} has no explicit task partition in {task.task_id}"
            )
        routing[dataset_id] = {
            "task_partition": partition,
            "registry_role": lock["role"],
            "label_visibility": lock["label_visibility"],
            "prediction_first_policy": lock["prediction_first_policy"],
            "prediction_first_policy_sha256": lock[
                "prediction_first_policy_sha256"
            ],
            "fit_allowed": partition == "fit",
            "model_selection_allowed": partition == "development_prediction",
            "prediction_first_outcome_artifact_included": False,
            "dataset_view_id": (
                lock.get("dataset_view", {}).get("view_id")
                if isinstance(lock.get("dataset_view"), Mapping)
                else None
            ),
            "dataset_view_contract_sha256": (
                lock.get("dataset_view", {}).get("contract_sha256")
                if isinstance(lock.get("dataset_view"), Mapping)
                else None
            ),
        }
    return {
        "fit_dataset_ids": list(fit_dataset_ids),
        "development_prediction_dataset_ids": list(
            development_prediction_dataset_ids
        ),
        "prediction_first_stress_dataset_ids": list(
            prediction_first_stress_dataset_ids
        ),
        "sealed_prediction_dataset_ids": list(sealed_prediction_dataset_ids),
        "dataset_routing": routing,
        "dataset_routing_wave": wave,
    }


def _adaptation_regimes(campaign: Mapping[str, Any], model_id: str) -> tuple[str, ...]:
    wave = str(campaign.get("wave", ""))
    execution = campaign.get("execution", {})
    if wave == "admission":
        return ("admission",)
    raw = execution.get("adaptation_regimes")
    regimes: list[str]
    if isinstance(raw, list) and raw:
        regimes = [str(item) for item in raw]
    else:
        single = str(execution.get("adaptation_regime", ""))
        regimes = (
            ["common_lane", "native_lane"]
            if single == "common_and_native"
            else [single or wave]
        )
    if model_id == "corgi_regular" and isinstance(execution.get("corgi_ladder"), list):
        regimes.extend(f"corgi_{item}" for item in execution["corgi_ladder"])
    return tuple(dict.fromkeys(regimes))


def _model_has_runnable_task_contract(
    model: ModelManifest,
    *,
    task_id: str,
    actions: Iterable[str],
    adaptation_regimes: Iterable[str],
) -> bool:
    """Return whether a registered model can emit a run for this exact task."""

    execution = model.execution
    if (
        model.admission_blocking
        or model.status.value in {"blocked", "blocked_terms", "deferred"}
        or not execution.ready
    ):
        return False
    supported_tasks = execution.executable_tasks or model.supported_tasks
    if task_id not in supported_tasks:
        return False
    if not set(actions).issubset(execution.supported_actions):
        return False
    return bool(
        set(adaptation_regimes).intersection(
            execution.supported_adaptation_regimes
        )
    )


def _resolved_model_hyperparameters(
    model_id: str,
    *,
    default: Mapping[str, Any],
    by_model: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Resolve one exact model schema; model tables replace rather than merge."""

    return dict(by_model.get(model_id, default))


def _common_head_ids(campaign: Mapping[str, Any]) -> tuple[str, ...]:
    execution = campaign.get("execution", {})
    raw = execution.get("common_head_ids")
    if raw is None:
        return ()
    if not isinstance(raw, list) or not raw:
        raise PlanningError(
            "execution.common_head_ids must be a non-empty unique array"
        )
    head_ids = tuple(
        _require_identifier(item, "execution.common_head_ids item") for item in raw
    )
    if any(head_id != head_id.lower() for head_id in head_ids):
        raise PlanningError("execution.common_head_ids must use lowercase identifiers")
    if len(set(head_ids)) != len(head_ids):
        raise PlanningError("execution.common_head_ids must not contain duplicates")
    return head_ids


def _contains_unresolved(value: Any) -> bool:
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.upper() in _UNRESOLVED_SENTINELS:
            return True
        return bool(
            re.search(
                r"(?:^|[^A-Z0-9])(?:UNKNOWN|UNRESOLVED|TBD|LATEST)(?:$|[^A-Z0-9])",
                stripped.upper(),
            )
        )
    if isinstance(value, Mapping):
        return any(
            _contains_unresolved(item)
            for key, item in value.items()
            if key not in _UNRESOLVED_PROSE_KEYS
        )
    if isinstance(value, (list, tuple)):
        return any(_contains_unresolved(item) for item in value)
    return False


def _verified_artifact(
    artifact: ArtifactRef,
    *,
    config_root: Path,
    expected_role: str,
    label: str,
) -> ArtifactRef:
    if artifact.role != expected_role:
        raise PlanningError(
            f"{label} must use ArtifactRef.role={expected_role!r}, got {artifact.role!r}"
        )
    configured = Path(artifact.path)
    candidate = configured if configured.is_absolute() else config_root / configured
    if candidate.is_symlink():
        raise PlanningError(f"{label} may not be a symlink: {candidate}")
    try:
        resolved = artifact.validate(None if configured.is_absolute() else config_root)
    except ContractError as error:
        raise PlanningError(f"{label} failed artifact verification: {error}") from error
    return ArtifactRef(
        path=resolved.as_posix(),
        sha256=artifact.sha256,
        size_bytes=artifact.size_bytes,
        media_type=artifact.media_type,
        role=artifact.role,
    )


def _campaign_prerequisite_artifacts(
    campaign: Mapping[str, Any], *, config_root: Path, scientific: bool
) -> tuple[ArtifactRef, ...]:
    raw_prerequisites = campaign.get("prerequisites", [])
    if not isinstance(raw_prerequisites, list) or any(
        not isinstance(item, str) or not _IDENTIFIER.fullmatch(item)
        for item in raw_prerequisites
    ):
        if scientific:
            raise PlanningError(
                "scientific campaign prerequisites must be unique identifier strings"
            )
        return ()
    prerequisites = tuple(raw_prerequisites)
    if len(set(prerequisites)) != len(prerequisites):
        raise PlanningError("campaign prerequisites must not contain duplicates")
    raw_receipts = campaign.get("prerequisite_receipts", {})
    if not isinstance(raw_receipts, Mapping):
        raise PlanningError("campaign prerequisite_receipts must be a table")
    if not scientific and not prerequisites:
        if raw_receipts:
            raise PlanningError("admission campaign has receipts but no prerequisites")
        return ()
    if set(raw_receipts) != set(prerequisites):
        raise PlanningError(
            "scientific campaign prerequisite_receipts must cover every and only "
            "declared prerequisite"
        )
    artifacts: list[ArtifactRef] = []
    for prerequisite in sorted(prerequisites):
        try:
            artifact = ArtifactRef.from_dict(raw_receipts[prerequisite])
        except ContractError as error:
            raise PlanningError(
                f"invalid prerequisite receipt for {prerequisite}: {error}"
            ) from error
        verified = _verified_artifact(
            artifact,
            config_root=config_root,
            expected_role=f"prerequisite:{prerequisite}",
            label=f"prerequisite receipt {prerequisite}",
        )
        if scientific:
            try:
                payload = json.loads(Path(verified.path).read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as error:
                raise PlanningError(
                    f"prerequisite {prerequisite} is not a readable receipt: {error}"
                ) from error
            if not isinstance(payload, Mapping):
                raise PlanningError(
                    f"prerequisite {prerequisite} receipt must be a JSON object"
                )
            schema_version = payload.get("schema_version")
            if schema_version == "masld-bench-run-execution-receipt-v1":
                try:
                    from .campaign import CampaignError, verify_run_execution_attempt

                    receipt = verify_run_execution_attempt(Path(verified.path).parent)
                except CampaignError as error:
                    raise PlanningError(
                        f"prerequisite {prerequisite} execution receipt failed "
                        f"recursive verification: {error}"
                    ) from error
                if receipt.get("status") != "succeeded":
                    raise PlanningError(
                        f"prerequisite {prerequisite} execution did not succeed"
                    )
            elif schema_version == "masld-bench-selection-lock-v1":
                try:
                    from .selection import SelectionError, verify_selection_lock

                    verify_selection_lock(Path(verified.path).parent).require_locked()
                except (SelectionError, ContractError) as error:
                    raise PlanningError(
                        f"prerequisite {prerequisite} SelectionLock failed "
                        f"recursive verification: {error}"
                    ) from error
            else:
                raise PlanningError(
                    f"prerequisite {prerequisite} has unsupported receipt schema "
                    f"{schema_version!r}; prose or placeholder artifacts are forbidden"
                )
        artifacts.append(verified)
    return tuple(artifacts)


def _runtime_lock(
    runtime_id: str,
    *,
    runtimes: Mapping[str, Mapping[str, Any]],
    profiles: Mapping[str, ResourceProfile],
) -> dict[str, Any]:
    if runtime_id not in runtimes:
        raise PlanningError(f"unknown runtime_id: {runtime_id}")
    runtime = runtimes[runtime_id]
    profile_id = _require_identifier(
        runtime.get("resource_profile"), f"runtime {runtime_id} resource_profile"
    )
    if profile_id not in profiles:
        raise PlanningError(f"runtime {runtime_id} names unknown profile {profile_id}")
    profiles[profile_id].validate(for_scheduling=True)
    prefix = runtime.get("command_prefix")
    if not isinstance(prefix, list) or not prefix or any(
        not isinstance(item, str) or not item.strip() for item in prefix
    ):
        raise PlanningError(f"runtime {runtime_id} has no exact command_prefix")
    unresolved = any(
        _contains_unresolved(value)
        for value in (
            runtime.get("environment_lock"),
            runtime.get("container_digest"),
            runtime.get("modules", []),
        )
    )
    if bool(runtime.get("admission_blocking", False)) or unresolved:
        reason = "; ".join(str(item) for item in runtime.get("blockers", []))
        raise PlanningError(
            f"runtime {runtime_id} is admission-blocking: "
            + (reason or "unresolved runtime identity")
        )
    return {
        "runtime_id": runtime_id,
        "registry_path": runtime["registry_path"],
        "registry_sha256": runtime["registry_sha256"],
        "environment_lock": runtime.get("environment_lock"),
        "container_digest": runtime.get("container_digest"),
        "modules": list(runtime.get("modules", [])),
        "command_prefix": list(prefix),
        "resource_profile": profile_id,
    }


def _adapter_actions(campaign: Mapping[str, Any]) -> tuple[str, ...]:
    wave = str(campaign.get("wave", ""))
    if wave == "admission":
        return ("probe",)
    expected = _SCIENTIFIC_WAVE_ACTIONS.get(wave)
    configured = campaign.get("execution", {}).get("actions")
    if expected is not None:
        if configured is not None and tuple(configured) != expected:
            raise PlanningError(
                f"campaign actions may not weaken the locked {wave} action sequence"
            )
        return expected
    if not isinstance(configured, list) or not configured:
        raise PlanningError(
            f"unsupported scientific wave {wave!r} requires explicit execution.actions"
        )
    actions = tuple(str(item) for item in configured)
    unknown = sorted(set(actions).difference(RunSpec._ADAPTER_ACTIONS))
    if unknown or len(set(actions)) != len(actions):
        raise PlanningError("execution.actions must be unique supported adapter actions")
    return actions


def _folds(campaign: Mapping[str, Any]) -> tuple[int, ...]:
    raw = campaign.get("execution", {}).get("folds", [0])
    if not isinstance(raw, list) or not raw:
        raise PlanningError("execution.folds must be a non-empty integer array")
    if any(isinstance(item, bool) or not isinstance(item, int) or item < 0 for item in raw):
        raise PlanningError("execution.folds must contain non-negative integers")
    if len(set(raw)) != len(raw):
        raise PlanningError("execution.folds must not contain duplicates")
    return tuple(raw)


def _campaign_dataset_views(
    campaign: Mapping[str, Any],
    *,
    registry: Registry,
    selected_task_ids: tuple[str, ...],
    wave: str,
) -> dict[str, Any]:
    """Resolve purpose-limited smoke views by parent dataset identifier."""

    raw = campaign.get("dataset_view_ids", [])
    if not isinstance(raw, list) or any(not isinstance(item, str) for item in raw):
        raise PlanningError("campaign dataset_view_ids must be a string array")
    if len(set(raw)) != len(raw):
        raise PlanningError("campaign dataset_view_ids must not contain duplicates")
    if raw and wave != "smoke":
        raise PlanningError("dataset views are restricted to the smoke wave")
    selected: dict[str, Any] = {}
    selected_train_ids = {
        dataset_id
        for task_id in selected_task_ids
        for dataset_id in registry.tasks[task_id].datasets_train
    }
    nontrain_ids = {
        dataset_id
        for task_id in selected_task_ids
        for dataset_id in (
            *registry.tasks[task_id].datasets_development,
            *registry.tasks[task_id].datasets_sealed,
        )
    }
    for view_id in raw:
        if view_id not in registry.dataset_views:
            raise PlanningError(f"campaign names unknown dataset view: {view_id}")
        view = registry.dataset_views[view_id]
        if wave not in view.allowed_waves:
            raise PlanningError(
                f"dataset view {view_id} does not permit campaign wave {wave}"
            )
        if view.parent_dataset_id not in selected_train_ids:
            raise PlanningError(
                f"dataset view {view_id} parent is not a selected task training dataset"
            )
        if view.parent_dataset_id in nontrain_ids:
            raise PlanningError(
                f"dataset view {view_id} parent also occurs in a prediction partition"
            )
        if view.parent_dataset_id in selected:
            raise PlanningError(
                f"campaign selects multiple views for parent {view.parent_dataset_id}"
            )
        selected[view.parent_dataset_id] = view
    if selected:
        cell_budget = campaign.get("execution", {}).get("cell_budget")
        if (
            isinstance(cell_budget, bool)
            or not isinstance(cell_budget, int)
            or cell_budget < 1
        ):
            raise PlanningError(
                "dataset-view smoke campaigns must freeze a positive execution.cell_budget"
            )
        if any(view.row_count != cell_budget for view in selected.values()):
            raise PlanningError(
                "every selected dataset view row_count must equal execution.cell_budget"
            )
    return selected


def build_plan(
    campaign_path: str | Path,
    *,
    config_root: str | Path,
    package_root: str | Path,
) -> dict[str, Any]:
    campaign_source = _reject_configured_symlinks(campaign_path, "campaign file")
    try:
        campaign_file = campaign_source.resolve(strict=True)
    except OSError as error:
        raise PlanningError(f"campaign file cannot be resolved: {campaign_source}") from error
    try:
        validation = validate_registry_tree(config_root)
    except RegistryError as error:
        raise PlanningError(f"registry validation failed: {error}") from error
    config = _reject_configured_symlinks(
        config_root, "configuration root"
    ).resolve(strict=True)
    registry_snapshot = validation["registry_snapshot"]
    registry_snapshot_sha256 = str(validation["registry_snapshot_sha256"])
    core_registry_contract_sha256 = str(
        validation["core_registry_contract_sha256"]
    )
    try:
        campaign_relative = campaign_file.relative_to(config).as_posix()
    except ValueError as error:
        raise PlanningError("campaign file must be inside the configuration root") from error
    campaign_file_sha256 = sha256_file(campaign_file)
    campaign_records = {
        str(record["path"]): record for record in registry_snapshot["files"]
    }
    campaign_record = campaign_records.get(campaign_relative)
    if campaign_record is None:
        raise PlanningError("campaign file is absent from the full configuration snapshot")
    if campaign_record["sha256"] != campaign_file_sha256:
        raise PlanningError("campaign file changed after the configuration snapshot")
    campaign = dict(_load_toml(campaign_file))
    campaign_id = _require_identifier(campaign.get("campaign_id"), "campaign_id")
    wave = _require_identifier(campaign.get("wave"), "campaign.wave")
    scientific = wave != "admission"
    review_only = campaign.get("review_only", False)
    if not isinstance(review_only, bool):
        raise PlanningError("campaign.review_only must be boolean")
    if review_only and wave != "admission":
        raise PlanningError("review-only campaigns must use wave='admission'")
    campaign_status = str(campaign.get("status", "ready"))
    if campaign_status.startswith("blocked_"):
        raise PlanningError(
            f"campaign {campaign_id} is a template blocked by prerequisites: {campaign_status}"
        )
    if scientific and campaign_status != "ready":
        raise PlanningError(
            f"scientific campaign {campaign_id} must set status='ready'"
        )
    if scientific and _contains_unresolved(campaign):
        raise PlanningError(
            f"scientific campaign {campaign_id} contains unresolved control values"
        )
    conditional_model_spec_binding: dict[str, str] | None = None
    if wave == "conditional_model":
        configured_spec = campaign.get("conditional_model_spec_path")
        configured_sha256 = campaign.get("conditional_model_spec_sha256")
        if not isinstance(configured_spec, str) or not configured_spec:
            raise PlanningError(
                "conditional-model campaign requires conditional_model_spec_path"
            )
        try:
            # Local import avoids a module cycle: context_spec recursively
            # verifies a conditional decision through tournament, while
            # tournament imports load_frozen_plan from this module.
            from .context_spec import ContextSpecError, verify_context_model_spec

            spec = verify_context_model_spec(configured_spec)
            spec_root = _reject_configured_symlinks(
                configured_spec, "conditional model spec"
            ).resolve(strict=True)
        except (ContextSpecError, OSError, ValueError) as error:
            raise PlanningError(f"invalid conditional model spec: {error}") from error
        if configured_sha256 != spec["context_model_spec_id"]:
            raise PlanningError(
                "conditional_model_spec_sha256 differs from the verified spec identity"
            )
        conditional_model_spec_binding = {
            "path": spec_root.as_posix(),
            "manifest_sha256": sha256_file(spec_root / "ARTIFACTS.json"),
            "document_sha256": sha256_file(spec_root / "context_model_spec.json"),
            "context_model_spec_id": str(spec["context_model_spec_id"]),
        }
    if bool(campaign.get("allow_arrays", False)):
        raise PlanningError("array campaigns require a separate explicit implementation and review")
    if bool(campaign.get("allow_sealed_labels", False)):
        raise PlanningError("campaign planning may never read project-sealed labels")
    try:
        strict_registry = Registry.load(config, validate_references=False)
    except RegistryError as error:
        raise PlanningError(f"strict registry validation failed: {error}") from error
    models = _models(config)
    tasks = _tasks(config)
    datasets = _datasets(config)
    splits = _splits(config)
    strict_models = strict_registry.models
    strict_tasks = strict_registry.tasks
    strict_datasets = strict_registry.datasets
    strict_splits = strict_registry.splits
    variant_capability_records: dict[str, dict[str, Any]] = {}
    variant_capability_registry_sha256: str | None = None
    if strict_registry.variant_capabilities is not None:
        variant_capability_document = canonicalize(
            strict_registry.variant_capabilities
        )
        variant_capability_records = {
            str(record["model_id"]): canonicalize(record)
            for record in variant_capability_document["models"]
        }
        variant_capability_registry_sha256 = canonical_hash(
            variant_capability_document
        )
    selected_task_ids = tuple(str(item) for item in campaign.get("task_ids", []))
    if not selected_task_ids or len(set(selected_task_ids)) != len(selected_task_ids):
        raise PlanningError("campaign task_ids must be a non-empty unique array")
    unknown_tasks = sorted(set(selected_task_ids).difference(tasks))
    if unknown_tasks:
        raise PlanningError(f"campaign names unknown tasks: {', '.join(unknown_tasks)}")
    selected_dataset_views = _campaign_dataset_views(
        campaign,
        registry=strict_registry,
        selected_task_ids=selected_task_ids,
        wave=wave,
    )
    selected_statuses = frozenset(
        str(item) for item in campaign.get("selection", {}).get("model_statuses", [])
    )
    default_profile = _require_identifier(
        campaign.get("execution", {}).get("default_resource_profile"),
        "execution.default_resource_profile",
    )
    default_runtime_id = _require_identifier(
        campaign.get("execution", {}).get("runtime_id"), "execution.runtime_id"
    )
    raw_seeds = campaign.get("execution", {}).get("seeds", [])
    if not isinstance(raw_seeds, list) or not raw_seeds or any(
        isinstance(seed, bool) or not isinstance(seed, int) or seed < 0
        for seed in raw_seeds
    ):
        raise PlanningError("campaign must freeze non-negative integer seeds")
    seeds = tuple(raw_seeds)
    if len(set(seeds)) != len(seeds):
        raise PlanningError("campaign seeds must not contain duplicates")
    folds = _folds(campaign)
    actions = _adapter_actions(campaign)
    resources_path = config / "resources.toml"
    profiles = load_resource_profiles(_resource_document(config))
    if default_profile not in profiles:
        raise PlanningError(f"unknown default resource profile: {default_profile}")
    profiles[default_profile].validate(for_scheduling=True)
    runtimes = _runtimes(config)
    default_runtime_lock = _runtime_lock(
        default_runtime_id, runtimes=runtimes, profiles=profiles
    )
    if default_runtime_lock["resource_profile"] != default_profile:
        raise PlanningError(
            f"runtime {default_runtime_id} requires profile "
            f"{default_runtime_lock['resource_profile']}, "
            f"not {default_profile}"
        )
    prerequisite_artifacts = _campaign_prerequisite_artifacts(
        campaign, config_root=config, scientific=scientific
    )

    package_source = _reject_configured_symlinks(package_root, "package root")
    locked_source = source_lock(package_source)
    package = package_source.resolve(strict=True)
    source_lock_sha256 = str(locked_source["source_lock_sha256"])

    dataset_locks: dict[str, dict[str, Any]] = {}
    dataset_artifacts: dict[str, ArtifactRef] = {}
    dataset_authority_artifacts: dict[str, tuple[ArtifactRef, ...]] = {}
    dataset_view_supporting_artifacts: dict[str, tuple[ArtifactRef, ...]] = {}
    for task_id in selected_task_ids:
        for dataset_id in strict_tasks[task_id].dataset_ids:
            if dataset_id not in datasets:
                raise PlanningError(f"task {task_id} names unknown dataset {dataset_id}")
            dataset = datasets[dataset_id]
            contract = strict_datasets[dataset_id]
            activation_sha256 = canonical_hash(contract.activation.to_dict())
            artifact: ArtifactRef | None = None
            if contract.activation.ready:
                assert contract.activation.artifact_manifest is not None
                artifact = _verified_artifact(
                    contract.activation.artifact_manifest,
                    config_root=config,
                    expected_role=f"dataset:{dataset_id}",
                    label=f"dataset activation {dataset_id}",
                )
                dataset_artifacts[dataset_id] = artifact
                authority_artifacts = tuple(
                    _verified_artifact(
                        authority_artifact,
                        config_root=config,
                        expected_role=str(authority_artifact.role),
                        label=(
                            f"dataset activation authority {dataset_id} "
                            f"{authority_artifact.role}"
                        ),
                    )
                    for authority_artifact in contract.activation.evidence_artifacts
                )
                dataset_authority_artifacts[dataset_id] = authority_artifacts
            dataset_locks[dataset_id] = {
                "registry_path": dataset["registry_path"],
                "registry_sha256": dataset["registry_sha256"],
                "role": dataset.get("role"),
                "status": dataset.get("status"),
                "label_visibility": dataset.get("split", {}).get("label_visibility"),
                "admission_blocking": bool(dataset.get("admission_blocking", False)),
                "blockers": list(dataset.get("blockers", [])),
                "activation_ready": contract.activation.ready,
                "activation_sha256": activation_sha256,
                "prediction_first_policy": (
                    contract.prediction_first_outcome.to_dict()
                    if contract.prediction_first_outcome is not None
                    else None
                ),
                "prediction_first_policy_sha256": (
                    canonical_hash(contract.prediction_first_outcome.to_dict())
                    if contract.prediction_first_outcome is not None
                    else None
                ),
                "artifact_manifest": artifact.to_dict() if artifact is not None else None,
                "authority_artifacts": [
                    item.to_dict()
                    for item in dataset_authority_artifacts.get(dataset_id, ())
                ],
            }

    for dataset_id, view in selected_dataset_views.items():
        view_manifest = _verified_artifact(
            view.artifact_manifest,
            config_root=config,
            expected_role=f"dataset_view_manifest:{view.view_id}",
            label=f"dataset view manifest {view.view_id}",
        )
        view_authorities = tuple(
            _verified_artifact(
                artifact,
                config_root=config,
                expected_role=str(artifact.role),
                label=f"dataset view authority {view.view_id} {artifact.role}",
            )
            for artifact in view.authority_artifacts
        )
        view_supporting = tuple(
            _verified_artifact(
                artifact,
                config_root=config,
                expected_role=str(artifact.role),
                label=f"dataset view artifact {view.view_id} {artifact.role}",
            )
            for artifact in (
                view.complete_receipt,
                view.data_artifact,
                view.selection_artifact,
                view.ontology_artifact,
                view.derivation_artifact,
            )
        )
        parent_lock = dataset_locks[dataset_id]
        parent_lock["parent_activation_ready"] = parent_lock["activation_ready"]
        parent_lock["parent_activation_sha256"] = parent_lock["activation_sha256"]
        parent_lock["activation_ready"] = True
        parent_lock["activation_sha256"] = canonical_hash(view.to_dict())
        parent_lock["activation_source"] = "dataset_view"
        parent_lock["artifact_manifest"] = view_manifest.to_dict()
        parent_lock["authority_artifacts"] = [
            artifact.to_dict() for artifact in view_authorities
        ]
        parent_lock["view_supporting_artifacts"] = [
            artifact.to_dict() for artifact in view_supporting
        ]
        parent_lock["dataset_view"] = {
            "view_id": view.view_id,
            "parent_dataset_id": view.parent_dataset_id,
            "purpose": view.purpose,
            "allowed_waves": list(view.allowed_waves),
            "row_count": view.row_count,
            "biological_unit_count": view.biological_unit_count,
            "contract_sha256": canonical_hash(view.to_dict()),
            "registry_path": (
                config / "dataset_views" / f"{view.view_id}.toml"
            ).relative_to(config).as_posix(),
            "registry_sha256": sha256_file(
                config / "dataset_views" / f"{view.view_id}.toml"
            ),
        }
        dataset_artifacts[dataset_id] = view_manifest
        dataset_authority_artifacts[dataset_id] = view_authorities
        dataset_view_supporting_artifacts[dataset_id] = view_supporting

    task_dispositions: list[dict[str, Any]] = []
    task_eligible: dict[str, bool] = {}
    task_run_dataset_ids: dict[str, tuple[str, ...]] = {}
    task_evaluator_artifacts: dict[str, ArtifactRef] = {}
    for task_id in selected_task_ids:
        task = tasks[task_id]
        task_contract = strict_tasks[task_id]
        split_contract = strict_splits[task_contract.split_id]
        split_registry = splits[task_contract.split_id]
        task_parameters = task_contract.to_dict()["evaluator_parameters"]
        raw_roster_authority = (
            task_parameters.get("roster_authority")
            if isinstance(task_parameters, Mapping)
            else None
        )
        evaluator_authority: ArtifactRef | None = None
        if raw_roster_authority is not None:
            try:
                evaluator_authority = _verified_artifact(
                    ArtifactRef.from_dict(raw_roster_authority),
                    config_root=config,
                    expected_role=f"task_evaluator_roster:{task_id}",
                    label=f"task evaluator roster {task_id}",
                )
            except ContractError as error:
                raise PlanningError(
                    f"task {task_id} evaluator roster artifact is invalid: {error}"
                ) from error
            task_evaluator_artifacts[task_id] = evaluator_authority
        blockers: list[str] = []
        omitted_datasets: list[dict[str, str]] = []
        if scientific and task.get("status") in {"blocked", "deferred"}:
            blockers.append(f"task_status_{task.get('status')}")
        if scientific and (
            task_contract.primary_evaluator_id == "UNRESOLVED"
            or _contains_unresolved(task_contract.evaluator_parameters)
        ):
            blockers.append("primary_evaluator_roster_authority_unresolved")
        sequence_reference_required = "dna_sequence" in task_contract.input_modalities
        if task_id == "rna_conditioned_atac" and sequence_reference_required:
            capabilities = strict_registry.rna_atac_capabilities
            selection_for_reference = campaign.get("selection", {})
            selected_for_reference = (
                selection_for_reference.get("model_ids")
                if isinstance(selection_for_reference, Mapping)
                else None
            )
            if capabilities is not None and isinstance(selected_for_reference, list):
                selected_reference_ids = {
                    str(model_id) for model_id in selected_for_reference
                }
                if bool(
                    selection_for_reference.get("include_mandatory_baselines", False)
                ):
                    selected_reference_ids.update(task_contract.baseline_model_ids)
                selected_reference_ids = {
                    model_id
                    for model_id in selected_reference_ids
                    if model_id in strict_models
                    and _model_has_runnable_task_contract(
                        strict_models[model_id],
                        task_id=task_id,
                        actions=actions,
                        adaptation_regimes=(
                            ("native_lane",)
                            if strict_models[model_id].family_id
                            == "mandatory_baselines"
                            else _adaptation_regimes(campaign, model_id)
                        ),
                    )
                }
                sequence_reference_required = bool(
                    selected_reference_ids.intersection(
                        capabilities["sequence_required_models"]
                    )
                )
        if (
            scientific
            and sequence_reference_required
            and not strict_registry.reference_bundle.sequence_extraction_ready
        ):
            blockers.append("sequence_reference_bundle_not_ready")
        full_candidate_dataset_ids = _run_dataset_ids(task_contract.to_dict(), wave)
        candidate_dataset_ids = full_candidate_dataset_ids
        admitted_dataset_ids: list[str] = []
        if selected_dataset_views:
            candidate_dataset_ids = tuple(
                dataset_id
                for dataset_id in full_candidate_dataset_ids
                if dataset_id in selected_dataset_views
            )
            omitted_datasets.extend(
                {
                    "dataset_id": dataset_id,
                    "reason": f"smoke_campaign_view_scope:{dataset_id}",
                }
                for dataset_id in full_candidate_dataset_ids
                if dataset_id not in selected_dataset_views
            )
            if not candidate_dataset_ids:
                blockers.append("no_selected_dataset_view_for_task")
        for dataset_id in candidate_dataset_ids:
            dataset_contract = strict_datasets[dataset_id]
            dataset_view = selected_dataset_views.get(dataset_id)
            if not scientific:
                admitted_dataset_ids.append(dataset_id)
                continue
            is_fit_dataset = dataset_id in task_contract.datasets_train
            if (
                is_fit_dataset
                and dataset_contract.prediction_first_outcome is not None
            ):
                blockers.append(
                    f"prediction_first_dataset_declared_for_fit:{dataset_id}"
                )
                continue
            prediction_first_policy = dataset_contract.prediction_first_outcome
            if wave == "prediction_first_stress" and prediction_first_policy is None:
                omitted_datasets.append(
                    {
                        "dataset_id": dataset_id,
                        "reason": f"not_prediction_first_stress:{dataset_id}",
                    }
                )
                continue
            if (
                prediction_first_policy is not None
                and wave != "prediction_first_stress"
            ):
                omitted_datasets.append(
                    {
                        "dataset_id": dataset_id,
                        "reason": (
                            "prediction_first_requires_dedicated_post_selection_wave:"
                            f"{dataset_id}"
                        ),
                    }
                )
                continue
            if (
                prediction_first_policy is not None
                and prediction_first_policy.prediction_universe_sha256
                == "UNRESOLVED"
            ):
                omitted_datasets.append(
                    {
                        "dataset_id": dataset_id,
                        "reason": (
                            f"prediction_universe_not_frozen:{dataset_id}"
                        ),
                    }
                )
                continue
            if not dataset_contract.activation.ready and dataset_view is None:
                reason = f"dataset_activation_not_ready:{dataset_id}"
                if is_fit_dataset or wave == "fnih_activation_audit":
                    blockers.append(reason)
                else:
                    omitted_datasets.append(
                        {"dataset_id": dataset_id, "reason": reason}
                    )
                continue
            if wave == "sealed_inference":
                if not bool(campaign.get("allow_sealed_features", False)):
                    blockers.append(f"sealed_features_not_authorized:{dataset_id}")
                    continue
            elif (
                dataset_view is None and dataset_contract.admission_blocking
            ) or dataset_contract.status.value != "available":
                reason = f"dataset_not_admitted:{dataset_id}"
                if is_fit_dataset or wave == "fnih_activation_audit":
                    blockers.append(reason)
                else:
                    omitted_datasets.append(
                        {"dataset_id": dataset_id, "reason": reason}
                    )
                    continue
            admitted_dataset_ids.append(dataset_id)
        if wave in {
            "sealed_inference",
            "fnih_activation_audit",
            "prediction_first_stress",
        } and not (
            admitted_dataset_ids
        ):
            blockers.append("no_eligible_prediction_dataset")
        task_run_dataset_ids[task_id] = tuple(admitted_dataset_ids)
        task_eligible[task_id] = not blockers
        disposition = "blocked"
        if not blockers:
            disposition = "eligible_partial" if omitted_datasets else "eligible"
        task_dispositions.append(
            {
                "task_id": task_id,
                "status": task.get("status"),
                "disposition": disposition,
                "blockers": sorted(set(blockers)),
                "run_dataset_ids": list(admitted_dataset_ids),
                "omitted_datasets": sorted(
                    omitted_datasets,
                    key=lambda item: (item["dataset_id"], item["reason"]),
                ),
                "registry_path": task["registry_path"],
                "registry_sha256": task["registry_sha256"],
                "contract_sha256": canonical_hash(task_contract.to_dict()),
                "split_id": split_contract.split_id,
                "split_registry_path": split_registry["registry_path"],
                "split_registry_sha256": split_registry["registry_sha256"],
                "split_contract": split_contract.to_dict(),
                "split_contract_sha256": canonical_hash(split_contract.to_dict()),
                "datasets_train": list(task_contract.datasets_train),
                "datasets_development": list(task_contract.datasets_development),
                "datasets_sealed": list(task_contract.datasets_sealed),
                "dataset_views": [
                    dataset_locks[dataset_id]["dataset_view"]
                    for dataset_id in admitted_dataset_ids
                    if "dataset_view" in dataset_locks[dataset_id]
                ],
                "baseline_model_ids": list(task_contract.baseline_model_ids),
                "uncertainty_method": task_contract.uncertainty_method,
                "resampling_units": list(task_contract.resampling_units),
                "bootstrap_replicates": task_contract.bootstrap_replicates,
                "multiplicity_family": task_contract.multiplicity_family,
                "multiplicity_method": task_contract.multiplicity_method,
                "promotion_gate_id": task_contract.promotion_gate_id,
                "promotion_gate_config_sha256": (
                    task_contract.promotion_gate_config_sha256
                ),
                "primary_evaluator_id": task_contract.primary_evaluator_id,
                "evaluator_parameters": task_contract.to_dict()[
                    "evaluator_parameters"
                ],
                "evaluator_contract_sha256": canonical_hash(
                    {
                        "primary_evaluator_id": task_contract.primary_evaluator_id,
                        "evaluator_parameters": task_contract.to_dict()[
                            "evaluator_parameters"
                        ],
                    }
                ),
                "evaluator_authority": (
                    evaluator_authority.to_dict()
                    if evaluator_authority is not None
                    else None
                ),
                "variant_capability_registry_sha256": (
                    variant_capability_registry_sha256
                    if task_id == "variant_to_regulation"
                    else None
                ),
            }
        )
    if scientific and not any(task_eligible.values()):
        details = "; ".join(
            f"{item['task_id']}={','.join(item['blockers'])}"
            for item in task_dispositions
        )
        raise PlanningError(
            "scientific campaign has zero eligible tasks after independent task "
            f"admission: {details}"
        )

    dispositions: list[dict[str, Any]] = []
    runs: list[dict[str, Any]] = []
    runtime_locks: dict[str, dict[str, Any]] = {
        default_runtime_id: default_runtime_lock
    }
    selection = campaign.get("selection", {})
    if not isinstance(selection, Mapping):
        raise PlanningError("campaign selection must be a table")
    configured_census_sha256 = selection.get("census_model_ids_sha256")
    observed_census_sha256 = canonical_hash(
        sorted(str(model["model_id"]) for model in models)
    )
    if configured_census_sha256 is not None:
        if (
            not isinstance(configured_census_sha256, str)
            or not re.fullmatch(r"[0-9a-f]{64}", configured_census_sha256)
        ):
            raise PlanningError(
                "selection.census_model_ids_sha256 must be a lowercase SHA-256"
            )
        if configured_census_sha256 != observed_census_sha256:
            raise PlanningError(
                "selection.census_model_ids_sha256 differs from the frozen model census"
            )
    requested_model_ids_raw = selection.get("model_ids")
    requested_model_ids: set[str] | None = None
    if requested_model_ids_raw is not None:
        if not isinstance(requested_model_ids_raw, list) or any(
            not isinstance(item, str) or item not in strict_models
            for item in requested_model_ids_raw
        ):
            raise PlanningError("selection.model_ids must name known models")
        if len(set(requested_model_ids_raw)) != len(requested_model_ids_raw):
            raise PlanningError("selection.model_ids must not contain duplicates")
        requested_model_ids = set(requested_model_ids_raw)
    if bool(selection.get("include_mandatory_baselines", False)):
        baseline_ids = {
            model_id
            for task_id in selected_task_ids
            for model_id in strict_tasks[task_id].baseline_model_ids
        }
        if requested_model_ids is not None:
            requested_model_ids |= baseline_ids
    if review_only:
        all_model_ids = {str(model["model_id"]) for model in models}
        all_model_statuses = {str(model.get("status")) for model in models}
        if set(selected_task_ids) != set(tasks):
            raise PlanningError(
                "review-only admission must include every registered task"
            )
        if requested_model_ids is not None and requested_model_ids != all_model_ids:
            raise PlanningError(
                "review-only admission model_ids must include every registered model"
            )
        if selected_statuses != all_model_statuses:
            raise PlanningError(
                "review-only admission model_statuses must cover the complete census"
            )

    execution_config = campaign.get("execution", {})
    raw_hyperparameters = execution_config.get("hyperparameters", {})
    raw_hyperparameters_by_model = execution_config.get(
        "hyperparameters_by_model", {}
    )
    raw_adapter_request = execution_config.get("adapter_request", {})
    if not isinstance(raw_hyperparameters, Mapping):
        raise PlanningError("execution.hyperparameters must be a table")
    if not isinstance(raw_hyperparameters_by_model, Mapping):
        raise PlanningError("execution.hyperparameters_by_model must be a table")
    unknown_parameter_models = sorted(
        set(str(model_id) for model_id in raw_hyperparameters_by_model).difference(
            strict_models
        )
    )
    if unknown_parameter_models:
        raise PlanningError(
            "execution.hyperparameters_by_model names unknown models: "
            + ", ".join(unknown_parameter_models)
        )
    if any(
        not isinstance(model_id, str) or not isinstance(parameters, Mapping)
        for model_id, parameters in raw_hyperparameters_by_model.items()
    ):
        raise PlanningError(
            "execution.hyperparameters_by_model must map model IDs to tables"
        )
    if not isinstance(raw_adapter_request, Mapping):
        raise PlanningError("execution.adapter_request must be a table")
    common_head_ids = _common_head_ids(campaign)
    if common_head_ids and (
        "common_head_id" in raw_hyperparameters
        or any(
            "common_head_id" in parameters
            for parameters in raw_hyperparameters_by_model.values()
        )
    ):
        raise PlanningError(
            "configured hyperparameters.common_head_id conflicts with "
            "execution.common_head_ids"
        )
    forbidden_request_fields = {
        "schema_version",
        "run_id",
        "action",
        "run_spec",
        "prior_action_outputs",
        "fit_dataset_ids",
        "development_prediction_dataset_ids",
        "prediction_first_stress_dataset_ids",
        "sealed_prediction_dataset_ids",
        "dataset_routing",
        "action_dataset_ids",
        "full_input_inventory_sha256",
        "action_input_inventory_sha256",
    }
    if forbidden_request_fields.intersection(raw_adapter_request):
        raise PlanningError(
            "execution.adapter_request may not override protocol identity fields"
        )
    requested_adapter_timeout = execution_config.get("adapter_timeout_seconds")
    if requested_adapter_timeout is not None and (
        isinstance(requested_adapter_timeout, bool)
        or not isinstance(requested_adapter_timeout, int)
        or requested_adapter_timeout < 1
    ):
        raise PlanningError("execution.adapter_timeout_seconds must be a positive integer")

    for model in models:
        if requested_model_ids is not None and model["model_id"] not in requested_model_ids:
            continue
        if selected_statuses and str(model.get("status")) not in selected_statuses:
            continue
        model_contract = strict_models[model["model_id"]]
        execution = model_contract.execution
        open_champion_decision = _open_champion_license_decision(model, execution)
        disposition, blockers = _model_disposition(
            model, open_champion_decision=open_champion_decision
        )
        execution_sha256 = canonical_hash(execution.to_dict())
        execution_blockers = list(blockers)
        model_authority_artifacts: tuple[ArtifactRef, ...] = ()
        if execution.ready:
            model_authority_artifacts = tuple(
                _verified_artifact(
                    authority_artifact,
                    config_root=config,
                    expected_role=str(authority_artifact.role),
                    label=(
                        f"model execution authority {model_contract.model_id} "
                        f"{authority_artifact.role}"
                    ),
                )
                for authority_artifact in execution.evidence_artifacts
            )
        if not execution.ready:
            execution_blockers.append("model_execution_contract_not_ready")
        unsupported_actions = sorted(set(actions).difference(execution.supported_actions))
        if unsupported_actions:
            execution_blockers.append(
                "unsupported_adapter_actions:" + ",".join(unsupported_actions)
            )
        requested_adaptation_regimes = _adaptation_regimes(
            campaign, model["model_id"]
        )
        if model_contract.family_id == "mandatory_baselines" and scientific:
            requested_adaptation_regimes = ("native_lane",)
        omitted_adaptation_regimes: tuple[str, ...] = ()
        executable_adaptation_regimes = requested_adaptation_regimes
        if scientific:
            supported_adaptation_regimes = frozenset(
                execution.supported_adaptation_regimes
            )
            omitted_adaptation_regimes = tuple(
                regime
                for regime in requested_adaptation_regimes
                if regime not in supported_adaptation_regimes
            )
            executable_adaptation_regimes = tuple(
                regime
                for regime in requested_adaptation_regimes
                if regime in supported_adaptation_regimes
            )
            if not executable_adaptation_regimes:
                execution_blockers.append("no_supported_adaptation_regime")
        variant_capability = variant_capability_records.get(
            str(model["model_id"])
        )
        dispositions.append(
            {
                "model_id": model["model_id"],
                "family_id": model["family_id"],
                "status": model.get("status"),
                "disposition": disposition,
                "champion_eligible": disposition == "eligible_open",
                "open_champion_eligible": bool(
                    open_champion_decision["eligible"]
                )
                and disposition == "eligible_open",
                "open_champion_license_decision": open_champion_decision,
                "open_champion_license_decision_sha256": canonical_hash(
                    open_champion_decision
                ),
                "blockers": sorted(set(execution_blockers)),
                "registry_path": model["registry_path"],
                "registry_sha256": model["registry_sha256"],
                "execution_ready": execution.ready,
                "execution_sha256": execution_sha256,
                "runtime_id": execution.runtime_id,
                "supported_actions": list(execution.supported_actions),
                "executable_tasks": list(execution.executable_tasks),
                "supported_adaptation_regimes": list(
                    execution.supported_adaptation_regimes
                ),
                "requested_adaptation_regimes": list(
                    requested_adaptation_regimes
                ),
                "omitted_adaptation_regimes": list(
                    omitted_adaptation_regimes
                ),
                "executable_adaptation_regimes": list(
                    executable_adaptation_regimes
                ),
                "common_head_ids": (
                    list(common_head_ids)
                    if "common_lane" in executable_adaptation_regimes
                    else []
                ),
                "authority_artifacts": [
                    artifact.to_dict() for artifact in model_authority_artifacts
                ],
                "variant_capability": variant_capability,
                "variant_capability_sha256": (
                    canonical_hash(variant_capability)
                    if variant_capability is not None
                    else None
                ),
            }
        )
        if review_only:
            continue
        if (
            not execution.ready
            or unsupported_actions
            or (scientific and not executable_adaptation_regimes)
        ):
            continue
        if scientific and disposition == "blocked":
            continue
        model_runtime_lock = _runtime_lock(
            execution.runtime_id, runtimes=runtimes, profiles=profiles
        )
        runtime_locks[execution.runtime_id] = model_runtime_lock
        adapter_command = (
            *model_runtime_lock["command_prefix"],
            *execution.adapter_command,
        )
        if execution.environment_artifact is None:
            raise PlanningError(
                f"ready model {model_contract.model_id} has no environment artifact"
            )
        environment_artifact = _verified_artifact(
            execution.environment_artifact,
            config_root=config,
            expected_role=f"environment:{execution.runtime_id}",
            label=f"model environment {model_contract.model_id}",
        )
        profile = profiles[str(model_runtime_lock["resource_profile"])]
        shutdown_buffer_seconds = max(60, min(900, profile.wall_seconds // 10))
        action_budget_seconds = profile.wall_seconds - shutdown_buffer_seconds
        adapter_timeout_seconds = (
            requested_adapter_timeout
            if requested_adapter_timeout is not None
            else action_budget_seconds // len(actions)
        )
        if (
            adapter_timeout_seconds < 1
            or adapter_timeout_seconds * len(actions) > action_budget_seconds
        ):
            raise PlanningError(
                f"adapter timeout for {model_contract.model_id} exceeds the "
                f"{profile.wall_time} Slurm wall time after a "
                f"{shutdown_buffer_seconds}-second receipt-freeze buffer"
            )
        dispositions[-1]["environment_artifact"] = environment_artifact.to_dict()
        dispositions[-1]["capture_r_session"] = execution.capture_r_session
        dispositions[-1]["adapter_timeout_seconds"] = adapter_timeout_seconds
        checkpoint: ArtifactRef | None = None
        if execution.checkpoint_artifact is not None:
            checkpoint = _verified_artifact(
                execution.checkpoint_artifact,
                config_root=config,
                expected_role=f"checkpoint:{model_contract.model_id}",
                label=f"model checkpoint {model_contract.model_id}",
            )
        supported = frozenset(
            execution.executable_tasks or model_contract.supported_tasks
        )
        for task_id in selected_task_ids:
            if task_id not in supported:
                continue
            if not task_eligible[task_id]:
                continue
            task_variant_capability: dict[str, Any] | None = None
            if task_id == "variant_to_regulation":
                task_variant_capability = variant_capability_records.get(
                    model_contract.model_id
                )
                if task_variant_capability is None:
                    raise PlanningError(
                        f"variant model {model_contract.model_id} has no capability record"
                    )
            task_contract = strict_tasks[task_id]
            run_dataset_ids = task_run_dataset_ids[task_id]
            dataset_access = _dataset_access_metadata(
                task=task_contract,
                run_dataset_ids=run_dataset_ids,
                dataset_locks=dataset_locks,
                wave=wave,
            )
            immutable_inputs = {
                dataset_id: dataset_locks[dataset_id]["activation_sha256"]
                for dataset_id in run_dataset_ids
            }
            inputs = [
                *prerequisite_artifacts,
                environment_artifact,
                *model_authority_artifacts,
            ]
            evaluator_authority = task_evaluator_artifacts.get(task_id)
            if evaluator_authority is not None:
                inputs.append(evaluator_authority)
            for dataset_id in run_dataset_ids:
                inputs.extend(dataset_authority_artifacts.get(dataset_id, ()))
                inputs.extend(dataset_view_supporting_artifacts.get(dataset_id, ()))
            if scientific:
                inputs.extend(dataset_artifacts[dataset_id] for dataset_id in run_dataset_ids)
            inputs.sort(key=lambda item: (item.role or "", item.path))
            dataset_input_owner: dict[str, str] = {}
            for dataset_id in run_dataset_ids:
                owned_artifacts = [
                    *dataset_authority_artifacts.get(dataset_id, ()),
                    *dataset_view_supporting_artifacts.get(dataset_id, ()),
                ]
                if scientific:
                    owned_artifacts.append(dataset_artifacts[dataset_id])
                for artifact in owned_artifacts:
                    assert artifact.role is not None
                    dataset_input_owner[artifact.role] = dataset_id
            input_owner_by_role = {
                str(artifact.role): dataset_input_owner.get(str(artifact.role))
                for artifact in inputs
            }
            full_input_inventory_sha256 = canonical_hash(
                [artifact.to_dict() for artifact in inputs]
            )
            regimes = executable_adaptation_regimes
            for seed in seeds:
                for fold in folds:
                    for adaptation_regime in regimes:
                        head_ids: tuple[str | None, ...] = (
                            tuple(common_head_ids)
                            if adaptation_regime == "common_lane" and common_head_ids
                            else (None,)
                        )
                        for common_head_id in head_ids:
                            run_hyperparameters = _resolved_model_hyperparameters(
                                str(model["model_id"]),
                                default=raw_hyperparameters,
                                by_model=raw_hyperparameters_by_model,
                            )
                            if common_head_id is not None:
                                run_hyperparameters["common_head_id"] = common_head_id
                            run_spec = RunSpec(
                                schema_version="masld-bench-run-v1",
                                campaign_id=campaign_id,
                                task_id=task_id,
                                model_id=model["model_id"],
                                dataset_ids=run_dataset_ids,
                                split_id=task_contract.split_id,
                                seed=seed,
                                stage=wave,
                                adaptation_regime=adaptation_regime,
                                fold=fold,
                                action=actions,
                                adapter_command=adapter_command,
                                immutable_inputs=immutable_inputs,
                                inputs=tuple(inputs),
                                checkpoint=checkpoint,
                                hyperparameters=run_hyperparameters,
                                runtime_id=execution.runtime_id,
                                runtime_registry_sha256=model_runtime_lock[
                                    "registry_sha256"
                                ],
                                resource_profile=model_runtime_lock[
                                    "resource_profile"
                                ],
                                code_lock_sha256=source_lock_sha256,
                                metadata={
                                    "wave": wave,
                                    "model_registry_sha256": model[
                                        "registry_sha256"
                                    ],
                                    "model_execution_sha256": execution_sha256,
                                    "open_champion_eligible": dispositions[-1][
                                        "open_champion_eligible"
                                    ],
                                    "open_champion_license_decision_sha256": (
                                        dispositions[-1][
                                            "open_champion_license_decision_sha256"
                                        ]
                                    ),
                                    "task_registry_sha256": tasks[task_id][
                                        "registry_sha256"
                                    ],
                                    "task_contract_sha256": canonical_hash(
                                        task_contract.to_dict()
                                    ),
                                    **dataset_access,
                                    "split_registry_sha256": splits[
                                        task_contract.split_id
                                    ]["registry_sha256"],
                                    "split_contract": strict_splits[
                                        task_contract.split_id
                                    ].to_dict(),
                                    "split_contract_sha256": canonical_hash(
                                        strict_splits[
                                            task_contract.split_id
                                        ].to_dict()
                                    ),
                                    "primary_evaluator_id": (
                                        task_contract.primary_evaluator_id
                                    ),
                                    "evaluator_parameters": task_contract.to_dict()[
                                        "evaluator_parameters"
                                    ],
                                    "evaluator_contract_sha256": canonical_hash(
                                        {
                                            "primary_evaluator_id": (
                                                task_contract.primary_evaluator_id
                                            ),
                                            "evaluator_parameters": (
                                                task_contract.to_dict()[
                                                    "evaluator_parameters"
                                                ]
                                            ),
                                        }
                                    ),
                                    "variant_capability": task_variant_capability,
                                    "variant_capability_sha256": (
                                        canonical_hash(task_variant_capability)
                                        if task_variant_capability is not None
                                        else None
                                    ),
                                    "variant_capability_registry_sha256": (
                                        variant_capability_registry_sha256
                                        if task_variant_capability is not None
                                        else None
                                    ),
                                    "primary_endpoint_scoring_allowed": (
                                        bool(
                                            task_variant_capability[
                                                "primary_eligible"
                                            ]
                                        )
                                        if task_variant_capability is not None
                                        else None
                                    ),
                                    "registry_snapshot_sha256": (
                                        registry_snapshot_sha256
                                    ),
                                    "campaign_file_sha256": campaign_file_sha256,
                                    "prerequisite_sha256s": [
                                        artifact.sha256
                                        for artifact in prerequisite_artifacts
                                    ],
                                    "adapter_request": dict(raw_adapter_request),
                                    "adapter_timeout_seconds": adapter_timeout_seconds,
                                    "common_head_id": common_head_id,
                                    "capture_r_session": execution.capture_r_session,
                                    "input_owner_by_role": input_owner_by_role,
                                    "full_input_inventory_sha256": (
                                        full_input_inventory_sha256
                                    ),
                                },
                            )
                            identity = run_spec.identity_payload
                            run = dict(identity)
                            run["run_id"] = run_spec.run_id
                            try:
                                reloaded = RunSpec.from_dict(identity)
                            except ContractError as error:
                                raise PlanningError(
                                    f"planner emitted an invalid RunSpec: {error}"
                                ) from error
                            if reloaded.run_id != run["run_id"]:
                                raise PlanningError(
                                    "RunSpec round-trip changed run identity"
                                )
                            runs.append(run)

    runs.sort(key=lambda item: item["run_id"])
    if review_only:
        disposition_ids = {str(item["model_id"]) for item in dispositions}
        expected_ids = {str(model["model_id"]) for model in models}
        if disposition_ids != expected_ids:
            raise PlanningError(
                "review-only admission dispositions must cover the complete model census"
            )
        if runs:
            raise PlanningError("review-only admission must not emit executable runs")
    if scientific:
        if not runs:
            raise PlanningError(
                "scientific campaign has no executable runs after contract admission"
            )
        run_pairs = {(run["task_id"], run["model_id"]) for run in runs}
        runnable_baseline_pairs = {
            (task_id, model_id)
            for task_id in selected_task_ids
            if task_eligible[task_id]
            for model_id in strict_tasks[task_id].baseline_model_ids
            if _model_has_runnable_task_contract(
                strict_models[model_id],
                task_id=task_id,
                actions=actions,
                adaptation_regimes=("native_lane",),
            )
        }
        for task_id in selected_task_ids:
            if task_eligible[task_id] and not any(
                pair[0] == task_id for pair in runnable_baseline_pairs
            ):
                raise PlanningError(
                    f"scientific campaign task {task_id} has no runnable mandatory baseline"
                )
        missing_baseline_runs = sorted(
            f"{task_id}:{model_id}"
            for task_id, model_id in runnable_baseline_pairs
            if (task_id, model_id) not in run_pairs
        )
        selected_baseline_model_ids = {
            model_id
            for task_id in selected_task_ids
            if task_eligible[task_id]
            for model_id in strict_tasks[task_id].baseline_model_ids
        }
        incremental_baseline_smoke = (
            wave == "smoke"
            and selection.get("include_mandatory_baselines") is False
            and isinstance(requested_model_ids_raw, list)
            and bool(requested_model_ids_raw)
            and set(requested_model_ids_raw).issubset(selected_baseline_model_ids)
            and all(
                (task_id, model_id) in run_pairs
                and model_id in strict_tasks[task_id].baseline_model_ids
                for task_id, model_id in run_pairs
            )
            and all(
                any(
                    (task_id, model_id) in run_pairs
                    for model_id in strict_tasks[task_id].baseline_model_ids
                )
                for task_id in selected_task_ids
                if task_eligible[task_id]
            )
        )
        if missing_baseline_runs and not incremental_baseline_smoke:
            raise PlanningError(
                "scientific campaign cannot drop mandatory baseline runs: "
                + ", ".join(missing_baseline_runs)
            )
    repository_text = _git_output(["git", "rev-parse", "--show-toplevel"], package)
    if not repository_text:
        raise PlanningError("campaign freezing requires a Git repository for the Resource firewall")
    repository_root = Path(repository_text.decode().strip())
    totals = resource_totals(
        profiles, (str(run["resource_profile"]) for run in runs)
    )
    resource_document = _resource_document(config)
    cluster = resource_document.get("cluster", {})
    aggregate_memory_cap = int(cluster.get("aggregate_memory_cap_gb", 7000))
    if int(totals["memory_gb_if_all_concurrent"]) > aggregate_memory_cap:
        raise PlanningError(
            "campaign could exceed the nslab aggregate memory cap if submitted "
            f"together: {totals['memory_gb_if_all_concurrent']} GB > "
            f"{aggregate_memory_cap} GB"
        )
    if _config_registry_snapshot(config) != registry_snapshot:
        raise PlanningError("configuration tree changed while the plan was being built")
    if sha256_file(campaign_file) != campaign_file_sha256:
        raise PlanningError("campaign file changed while the plan was being built")
    if source_lock(package) != locked_source:
        raise PlanningError("source tree changed while the plan was being built")
    control_plane_python = Path(sys.executable).resolve(strict=True)
    control_plane_lock = {
        "python_executable": control_plane_python.as_posix(),
        "python_executable_sha256": sha256_file(control_plane_python),
        "python_executable_size_bytes": control_plane_python.stat().st_size,
        "python_version": sys.version,
        "package_root": package.as_posix(),
        "pythonpath": (package / "src").resolve(strict=True).as_posix(),
        "module": "masld_bench.cli",
        "source_lock_sha256": source_lock_sha256,
    }
    plan_payload: dict[str, Any] = {
        "schema_version": "masld-bench-plan-v1",
        "campaign": campaign,
        "campaign_file": {
            "path": campaign_relative,
            "sha256": campaign_file_sha256,
        },
        "config_root": config.as_posix(),
        "package_root": package.as_posix(),
        "repository_root": repository_root.resolve().as_posix(),
        "resources_file": {
            "path": resources_path.relative_to(config).as_posix(),
            "sha256": sha256_file(resources_path),
        },
        "runtime_lock": default_runtime_lock,
        "runtime_locks": dict(sorted(runtime_locks.items())),
        "control_plane_lock": control_plane_lock,
        "retry_policy": {
            "max_attempts": 1,
            "retryable_states": [],
            "checkpoint_resume_required": True,
            "software_correction": "new_campaign_revision",
            "reason": (
                "No standardized immutable resume-checkpoint ArtifactRef contract "
                "is implemented; repeated execution fails closed."
            ),
        },
        "conditional_model_spec_binding": conditional_model_spec_binding,
        "prerequisite_artifacts": [
            artifact.to_dict() for artifact in prerequisite_artifacts
        ],
        "registry_snapshot": registry_snapshot,
        "registry_snapshot_sha256": registry_snapshot_sha256,
        "core_registry_contract_sha256": core_registry_contract_sha256,
        "source_lock": locked_source,
        "source_lock_sha256": source_lock_sha256,
        "resource_firewall": resource_snapshot(repository_root),
        "dataset_locks": dict(sorted(dataset_locks.items())),
        "model_dispositions": dispositions,
        "model_census_sha256": observed_census_sha256,
        "task_dispositions": task_dispositions,
        "admission_review": (
            {
                "review_only": True,
                "model_count": len(models),
                "model_ids": sorted(str(model["model_id"]) for model in models),
                "task_count": len(tasks),
                "task_ids": sorted(tasks),
                "expected_run_count": 0,
                "expected_job_count": 0,
            }
            if review_only
            else None
        ),
        "runs": runs,
        "dag": {
            "prerequisites": list(campaign.get("prerequisites", [])),
            "run_nodes": [run["run_id"] for run in runs],
            "edges": [
                {"from": str(prerequisite), "to": run["run_id"]}
                for prerequisite in campaign.get("prerequisites", [])
                for run in runs
            ],
        },
        "resource_profiles": {key: value.as_dict() for key, value in sorted(profiles.items())},
        "resource_totals": totals,
        "safety": {
            "arrays": False,
            "sealed_features": bool(campaign.get("allow_sealed_features", False)),
            "sealed_labels": False,
            "submit_enabled": bool(campaign.get("submit_enabled", False)),
        },
    }
    plan_payload["plan_sha256"] = canonical_hash(plan_payload)
    return plan_payload


def _review_text(plan: Mapping[str, Any]) -> str:
    totals = plan["resource_totals"]
    omitted_task_lines = []
    for task in plan.get("task_dispositions", []):
        if task.get("disposition") == "blocked":
            omitted_task_lines.append(
                f"- {task['task_id']}: " + ", ".join(task.get("blockers", []))
            )
        elif task.get("omitted_datasets"):
            details = ", ".join(
                f"{item['dataset_id']} ({item['reason']})"
                for item in task["omitted_datasets"]
            )
            omitted_task_lines.append(
                f"- {task['task_id']} partial dataset omission: {details}"
            )
    lines = [
        f"Campaign: {plan['campaign']['campaign_id']}",
        f"Plan SHA-256: {plan['plan_sha256']}",
        f"Jobs: {totals['jobs']}",
        f"Requested CPU-hours: {totals['cpu_hours_requested']}",
        f"Requested GPU-hours: {totals['gpu_hours_requested']}",
        f"Requested memory-GB-hours: {totals['memory_gb_hours_requested']}",
        f"CPUs if all jobs are concurrent: {totals['cpus_if_all_concurrent']}",
        f"GPUs if all jobs are concurrent: {totals['gpus_if_all_concurrent']}",
        f"Memory GB if all jobs are concurrent: {totals['memory_gb_if_all_concurrent']}",
        f"Arrays: {plan['safety']['arrays']}",
        f"Sealed features: {plan['safety']['sealed_features']}",
        f"Sealed labels: {plan['safety']['sealed_labels']}",
        f"Submission enabled: {plan['safety']['submit_enabled']}",
        "",
        "Blocked tasks and dataset omissions:",
        *(omitted_task_lines or ["- None"]),
        "",
        "Exact distinct sbatch headers are in sbatch_headers.txt.",
        "The frozen DAG is in dag.json.",
        "Approve the campaign SHA-256 emitted after the complete tree is frozen.",
        "Nothing has been submitted.",
        "",
    ]
    return "\n".join(lines)


def freeze_campaign(
    campaign_path: str | Path,
    *,
    config_root: str | Path,
    package_root: str | Path,
    output_root: str | Path,
) -> Path:
    package = _reject_configured_symlinks(
        package_root, "package root"
    ).resolve(strict=True)
    requested_output_root = _reject_configured_symlinks(
        output_root, "campaign output root"
    ).resolve()
    try:
        output_relative = requested_output_root.relative_to(package)
    except ValueError:
        output_relative = None
    if output_relative is not None and (
        not output_relative.parts
        or output_relative.parts[0] not in _SOURCE_OUTPUT_ROOTS
    ):
        raise PlanningError(
            "an output root inside the package must use a provenance-excluded "
            f"top-level directory: {', '.join(sorted(_SOURCE_OUTPUT_ROOTS))}"
        )
    plan = build_plan(campaign_path, config_root=config_root, package_root=package_root)
    campaign_id = str(plan["campaign"]["campaign_id"])
    target_root = requested_output_root
    target_root.mkdir(parents=True, exist_ok=True)
    target = target_root / f"{campaign_id}--{plan['plan_sha256'][:16]}"
    if target.exists():
        raise PlanningError(f"candidate already exists and is immutable: {target}")
    staging = Path(tempfile.mkdtemp(prefix=f".{campaign_id}.", dir=target_root))
    try:
        write_json_exclusive(staging / "plan.json", plan)
        write_json_exclusive(staging / "dag.json", plan["dag"])
        write_text_exclusive(staging / "plan.sha256", f"{plan['plan_sha256']}  plan.payload\n")
        write_text_exclusive(staging / "REVIEW.txt", _review_text(plan))
        profiles = {
            key: ResourceProfile.from_serialized(key, value)
            for key, value in plan["resource_profiles"].items()
        }
        script_headers: set[str] = set()
        scripts = staging / "jobs"
        execution_root = target_root / "executions" / target.name
        for run in plan["runs"]:
            run_id = str(run["run_id"])
            runtime = plan["runtime_locks"].get(str(run["runtime_id"]))
            if not isinstance(runtime, Mapping):
                raise PlanningError(
                    f"run {run_id} has no frozen model-specific runtime lock"
                )
            if runtime.get("registry_sha256") != run.get("runtime_registry_sha256"):
                raise PlanningError(
                    f"run {run_id} runtime registry binding changed"
                )
            control_plane = plan["control_plane_lock"]
            script = render_sbatch(
                job_name=f"mb-{run_id[:12]}",
                profile=profiles[str(run["resource_profile"])],
                command=(
                    str(control_plane["python_executable"]),
                    "-m",
                    str(control_plane["module"]),
                    "run",
                    "execute",
                    "--candidate",
                    target.as_posix(),
                    "--run-id",
                    run_id,
                    "--output-root",
                    execution_root.as_posix(),
                ),
                stdout_path=(execution_root / "logs" / f"{run_id}.out").as_posix(),
                stderr_path=(execution_root / "logs" / f"{run_id}.err").as_posix(),
                environment={
                    "PYTHONHASHSEED": str(run["seed"]),
                    "PYTHONPATH": str(control_plane["pythonpath"]),
                },
                modules=tuple(str(item) for item in runtime["modules"]),
                working_directory=str(control_plane["package_root"]),
            )
            write_text_exclusive(scripts / f"{run_id}.sbatch", script, mode=0o750)
            script_headers.add("\n".join(line for line in script.splitlines() if line.startswith("#")))
        write_text_exclusive(
            staging / "sbatch_headers.txt",
            "\n\n".join(sorted(script_headers)) + ("\n" if script_headers else ""),
        )
        freeze_tree(
            staging,
            {
                "campaign_id": campaign_id,
                "plan_sha256": plan["plan_sha256"],
                "artifact_class": "candidate_plan",
            },
        )
        try:
            publish_directory_noreplace(staging, target)
        except ArtifactError as error:
            raise PlanningError(f"candidate publication failed: {error}") from error
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return target


def load_frozen_plan(candidate: str | Path) -> dict[str, Any]:
    candidate_path = Path(candidate)
    plan_path = candidate_path / "plan.json"
    try:
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise PlanningError(f"cannot read candidate plan: {error}") from error
    claimed = str(plan.pop("plan_sha256", ""))
    actual = canonical_hash(plan)
    plan["plan_sha256"] = claimed
    if not claimed or claimed != actual:
        raise PlanningError(f"candidate plan hash mismatch: expected {claimed}, observed {actual}")
    return plan
