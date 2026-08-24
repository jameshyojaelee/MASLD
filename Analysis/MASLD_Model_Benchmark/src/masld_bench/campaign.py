"""Execute reviewed campaign jobs without weakening plan immutability."""

from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from typing import Any, Mapping

from .adapters import (
    AdapterAction,
    AdapterError,
    AdapterReceipt,
    SubprocessAdapter,
    expected_action_dataset_ids,
    withheld_input_roles_for_action,
)
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
    RunExecutionReceipt,
    RunSpec,
    RunState,
)
from .firewall import FirewallError, assert_resource_unchanged
from .planner import PlanningError, load_frozen_plan, source_lock
from .runtime import capture_runtime_lock


class CampaignError(RuntimeError):
    """Raised when a frozen campaign or run fails its execution contract."""


def _resolve_existing_directory(path: str | Path, label: str) -> Path:
    configured = Path(path)
    try:
        lexical = reject_symlink_components(configured, label=label)
        resolved = lexical.resolve(strict=True)
    except ArtifactError as error:
        raise CampaignError(f"{label} path rejected: {error}") from error
    except OSError as error:
        raise CampaignError(f"{label} is missing or unreadable: {configured}") from error
    if not resolved.is_dir():
        raise CampaignError(f"{label} is not a directory: {resolved}")
    return resolved


def _resolve_output_directory(path: str | Path, label: str) -> Path:
    configured = Path(path)
    try:
        lexical = reject_symlink_components(configured, label=label)
        resolved = lexical.resolve()
    except ArtifactError as error:
        raise CampaignError(f"{label} path rejected: {error}") from error
    except OSError as error:
        raise CampaignError(f"{label} is unreadable: {configured}") from error
    if resolved.exists() and not resolved.is_dir():
        raise CampaignError(f"{label} is not a directory: {resolved}")
    return resolved


def _resolve_existing_file(path: str | Path, label: str) -> Path:
    configured = Path(path)
    try:
        lexical = reject_symlink_components(configured, label=label)
        resolved = lexical.resolve(strict=True)
    except ArtifactError as error:
        raise CampaignError(f"{label} path rejected: {error}") from error
    except OSError as error:
        raise CampaignError(f"{label} is missing or unreadable: {configured}") from error
    if not resolved.is_file():
        raise CampaignError(f"{label} is not a file: {resolved}")
    return resolved


def _assert_hash(path: Path, expected: str, label: str) -> None:
    if not path.is_file():
        raise CampaignError(f"{label} is missing: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise CampaignError(f"{label} checksum mismatch: {path}")


def _next_attempt(run_root: Path) -> tuple[int, Path]:
    existing = []
    if run_root.is_dir():
        for path in run_root.glob("attempt-*"):
            try:
                existing.append(int(path.name.removeprefix("attempt-")))
            except ValueError:
                continue
    number = max(existing, default=0) + 1
    return number, run_root / f"attempt-{number:03d}"


def _read_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CampaignError(f"cannot read {label}: {error}") from error
    if not isinstance(value, Mapping):
        raise CampaignError(f"{label} must contain a JSON object")
    return dict(value)


def _candidate_run(
    candidate: Path, run_id: str
) -> tuple[dict[str, Any], dict[str, Any], RunSpec, str]:
    candidate = _resolve_existing_directory(candidate, "candidate")
    try:
        manifest = verify_frozen_tree(candidate)
    except ArtifactError as error:
        raise CampaignError(f"candidate is not an intact immutable tree: {error}") from error
    metadata = manifest.get("metadata")
    if not isinstance(metadata, Mapping) or metadata.get("artifact_class") != "candidate_plan":
        raise CampaignError("candidate tree is not a frozen campaign plan")
    plan = load_frozen_plan(candidate)
    if metadata.get("campaign_id") != plan.get("campaign", {}).get("campaign_id"):
        raise CampaignError("candidate manifest campaign identity does not match plan")
    if metadata.get("plan_sha256") != plan.get("plan_sha256"):
        raise CampaignError("candidate manifest plan identity does not match plan")
    matching = [run for run in plan.get("runs", []) if run.get("run_id") == run_id]
    if len(matching) != 1:
        raise CampaignError(f"run_id must identify exactly one frozen run: {run_id}")
    raw_run = dict(matching[0])
    claimed_run_id = raw_run.pop("run_id", None)
    try:
        run_spec = RunSpec.from_dict(raw_run)
    except ContractError as error:
        raise CampaignError(f"frozen run violates RunSpec: {error}") from error
    if claimed_run_id != run_id or run_spec.run_id != run_id:
        raise CampaignError("frozen run identity does not match its RunSpec")
    return manifest, plan, run_spec, sha256_file(candidate / "ARTIFACTS.json")


def _assert_config_snapshot(plan: Mapping[str, Any]) -> None:
    config_root = _resolve_existing_directory(
        Path(str(plan.get("config_root", ""))), "frozen configuration root"
    )
    snapshot = plan.get("registry_snapshot")
    if not isinstance(snapshot, Mapping) or snapshot.get("schema_version") != (
        "masld-bench-registry-file-snapshot-v1"
    ):
        raise CampaignError("frozen registry snapshot is invalid")
    raw_expected = snapshot.get("files")
    if not isinstance(raw_expected, list):
        raise CampaignError("frozen registry snapshot has no file inventory")
    expected: dict[str, Mapping[str, Any]] = {}
    for value in raw_expected:
        if not isinstance(value, Mapping):
            raise CampaignError("frozen registry file inventory contains a non-object")
        relative = str(value.get("path", ""))
        if relative in expected:
            raise CampaignError(f"duplicate frozen registry path: {relative}")
        expected[relative] = value
    observed_paths: set[str] = set()
    for path in sorted(config_root.rglob("*")):
        if path.is_symlink():
            raise CampaignError(f"configuration tree contains a symlink: {path}")
        if path.is_dir():
            continue
        if not path.is_file():
            raise CampaignError(
                f"configuration tree contains a non-regular member: {path}"
            )
        relative = path.relative_to(config_root).as_posix()
        observed_paths.add(relative)
    if set(expected) != observed_paths:
        changed = sorted(set(expected) ^ observed_paths)
        raise CampaignError(
            "configuration TOML membership changed after review: " + ", ".join(changed)
        )
    for relative, record in expected.items():
        path = _resolve_existing_file(
            config_root / relative, f"registry file {relative}"
        )
        expected_size = record.get("size_bytes")
        if (
            isinstance(expected_size, bool)
            or not isinstance(expected_size, int)
            or path.stat().st_size != expected_size
        ):
            raise CampaignError(f"registry file size changed after review: {relative}")
        _assert_hash(path, str(record.get("sha256", "")), f"registry file {relative}")
    if canonical_hash(snapshot) != plan.get("registry_snapshot_sha256"):
        raise CampaignError("frozen registry snapshot hash is invalid")


def _validate_external_artifact(artifact: ArtifactRef, label: str) -> Path:
    if not Path(artifact.path).is_absolute():
        raise CampaignError(f"{label} must use an absolute frozen artifact path")
    try:
        return artifact.validate()
    except ContractError as error:
        raise CampaignError(f"{label} failed artifact verification: {error}") from error


def _environment_artifact(run_spec: RunSpec) -> ArtifactRef:
    expected_role = f"environment:{run_spec.runtime_id}"
    matches = [artifact for artifact in run_spec.inputs if artifact.role == expected_role]
    if len(matches) != 1:
        raise CampaignError(
            "run must bind exactly one environment artifact with role "
            f"{expected_role}"
        )
    return matches[0]


def _assert_adapter_environment_receipt(
    adapter_receipt: AdapterReceipt,
    *,
    environment_artifact: ArtifactRef,
    runtime_id: str,
    fit_dataset_ids: tuple[str, ...],
) -> None:
    if adapter_receipt.metadata.get("environment_artifact_sha256") != (
        environment_artifact.sha256
    ):
        raise CampaignError(
            "adapter receipt does not attest the frozen environment artifact"
        )
    if adapter_receipt.metadata.get("runtime_id") != runtime_id:
        raise CampaignError("adapter receipt runtime_id differs from the RunSpec")
    if adapter_receipt.metadata.get("fit_dataset_ids") != list(fit_dataset_ids):
        raise CampaignError(
            "adapter receipt does not attest the frozen fit dataset identifiers"
        )


def _expected_dataset_access(
    *,
    run_spec: RunSpec,
    task_entry: Mapping[str, Any],
    dataset_locks: Mapping[str, Any],
) -> dict[str, Any]:
    run_ids = tuple(run_spec.dataset_ids)
    run_set = set(run_ids)
    train = tuple(str(item) for item in task_entry.get("datasets_train", []))
    development = tuple(
        str(item) for item in task_entry.get("datasets_development", [])
    )
    sealed = tuple(str(item) for item in task_entry.get("datasets_sealed", []))
    fit_ids = tuple(item for item in train if item in run_set)
    development_ids = tuple(
        item
        for item in development
        if item in run_set
        and dataset_locks[item].get("prediction_first_policy") is None
    )
    stress_ids = tuple(
        item
        for item in development
        if item in run_set
        and dataset_locks[item].get("prediction_first_policy") is not None
    )
    sealed_ids = tuple(item for item in sealed if item in run_set)
    routing: dict[str, dict[str, Any]] = {}
    for dataset_id in run_ids:
        lock = dataset_locks.get(dataset_id)
        if not isinstance(lock, Mapping):
            raise CampaignError(f"run names an unlocked dataset: {dataset_id}")
        if dataset_id in fit_ids:
            partition = "fit"
        elif dataset_id in development_ids:
            partition = "development_prediction"
        elif dataset_id in stress_ids:
            partition = "prediction_first_stress"
        elif dataset_id in sealed_ids:
            partition = "sealed_prediction"
        else:
            raise CampaignError(
                f"run dataset has no frozen task partition: {dataset_id}"
            )
        routing[dataset_id] = {
            "task_partition": partition,
            "registry_role": lock.get("role"),
            "label_visibility": lock.get("label_visibility"),
            "prediction_first_policy": lock.get("prediction_first_policy"),
            "prediction_first_policy_sha256": lock.get(
                "prediction_first_policy_sha256"
            ),
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
        "fit_dataset_ids": list(fit_ids),
        "development_prediction_dataset_ids": list(development_ids),
        "prediction_first_stress_dataset_ids": list(stress_ids),
        "sealed_prediction_dataset_ids": list(sealed_ids),
        "dataset_routing": routing,
        "dataset_routing_wave": run_spec.stage,
    }


def _assert_run_bindings(
    *, plan: Mapping[str, Any], run_spec: RunSpec
) -> tuple[Path, str]:
    _assert_config_snapshot(plan)
    if run_spec.code_lock_sha256 != plan.get("source_lock_sha256"):
        raise CampaignError("run is not bound to the frozen source lock")
    package_root = _resolve_existing_directory(
        Path(str(plan.get("package_root", Path(__file__).resolve().parents[2]))),
        "frozen package root",
    )
    try:
        observed_source = source_lock(package_root)
    except PlanningError as error:
        raise CampaignError(f"cannot revalidate source lock: {error}") from error
    if observed_source != plan.get("source_lock"):
        raise CampaignError("source tree differs from the reviewed campaign plan")
    control_plane = plan.get("control_plane_lock")
    if not isinstance(control_plane, Mapping):
        raise CampaignError("plan has no frozen control-plane runtime")
    if (
        control_plane.get("package_root") != package_root.as_posix()
        or control_plane.get("source_lock_sha256") != run_spec.code_lock_sha256
        or control_plane.get("module") != "masld_bench.cli"
    ):
        raise CampaignError("control-plane runtime binding changed")
    configured_control_python = Path(
        str(control_plane.get("python_executable", ""))
    )
    if not configured_control_python.is_absolute():
        raise CampaignError("frozen control-plane Python is unavailable")
    control_python = _resolve_existing_file(
        configured_control_python, "frozen control-plane Python"
    )
    expected_python_size = control_plane.get("python_executable_size_bytes")
    if (
        isinstance(expected_python_size, bool)
        or not isinstance(expected_python_size, int)
        or control_python.stat().st_size != expected_python_size
    ):
        raise CampaignError("frozen control-plane Python size changed")
    if sha256_file(control_python) != control_plane.get("python_executable_sha256"):
        raise CampaignError("frozen control-plane Python checksum changed")
    if control_plane.get("python_version") != sys.version:
        raise CampaignError("control-plane Python version differs from the frozen plan")
    pythonpath = _resolve_existing_directory(
        Path(str(control_plane.get("pythonpath", ""))),
        "frozen control-plane PYTHONPATH",
    )
    if pythonpath != package_root / "src":
        raise CampaignError("frozen control-plane PYTHONPATH is unavailable")

    runtime_locks = plan.get("runtime_locks")
    if not isinstance(runtime_locks, Mapping):
        raise CampaignError("plan has no model-specific runtime locks")
    runtime = runtime_locks.get(run_spec.runtime_id)
    if not isinstance(runtime, Mapping):
        raise CampaignError("run runtime_id has no frozen runtime lock")
    if runtime.get("registry_sha256") != run_spec.runtime_registry_sha256:
        raise CampaignError("run runtime registry binding changed")
    if runtime.get("resource_profile") != run_spec.resource_profile:
        raise CampaignError("run resource profile differs from its runtime contract")
    if list(run_spec.adapter_command[: len(runtime.get("command_prefix", []))]) != list(
        runtime.get("command_prefix", [])
    ):
        raise CampaignError("run adapter command does not begin with its runtime prefix")

    metadata = run_spec.metadata
    model_entries = [
        item
        for item in plan.get("model_dispositions", [])
        if item.get("model_id") == run_spec.model_id
    ]
    task_entries = [
        item
        for item in plan.get("task_dispositions", [])
        if item.get("task_id") == run_spec.task_id
    ]
    if len(model_entries) != 1 or len(task_entries) != 1:
        raise CampaignError("run model/task disposition is not unique")
    model_entry, task_entry = model_entries[0], task_entries[0]
    if metadata.get("model_registry_sha256") != model_entry.get("registry_sha256"):
        raise CampaignError("run model registry binding changed")
    if metadata.get("model_execution_sha256") != model_entry.get("execution_sha256"):
        raise CampaignError("run model execution binding changed")
    if metadata.get("open_champion_eligible") != model_entry.get(
        "open_champion_eligible"
    ):
        raise CampaignError("run open-champion eligibility binding changed")
    if metadata.get("open_champion_license_decision_sha256") != model_entry.get(
        "open_champion_license_decision_sha256"
    ):
        raise CampaignError("run open-champion license decision binding changed")
    if metadata.get("capture_r_session") != model_entry.get("capture_r_session"):
        raise CampaignError("run R-session capture binding changed")
    if metadata.get("adapter_timeout_seconds") != model_entry.get(
        "adapter_timeout_seconds"
    ):
        raise CampaignError("run adapter-timeout binding changed")
    if metadata.get("task_registry_sha256") != task_entry.get("registry_sha256"):
        raise CampaignError("run task registry binding changed")
    if metadata.get("task_contract_sha256") != task_entry.get("contract_sha256"):
        raise CampaignError("run task contract binding changed")
    if run_spec.split_id != task_entry.get("split_id"):
        raise CampaignError("run split_id differs from the task split contract")
    if metadata.get("split_registry_sha256") != task_entry.get(
        "split_registry_sha256"
    ):
        raise CampaignError("run split registry binding changed")
    if canonical_hash(metadata.get("split_contract")) != canonical_hash(
        task_entry.get("split_contract")
    ):
        raise CampaignError("run parsed split contract changed")
    if metadata.get("split_contract_sha256") != task_entry.get(
        "split_contract_sha256"
    ):
        raise CampaignError("run parsed split contract hash changed")
    if metadata.get("primary_evaluator_id") != task_entry.get(
        "primary_evaluator_id"
    ):
        raise CampaignError("run primary evaluator binding changed")
    if canonical_hash(metadata.get("evaluator_parameters")) != canonical_hash(
        task_entry.get("evaluator_parameters")
    ):
        raise CampaignError("run evaluator parameter binding changed")
    if metadata.get("evaluator_contract_sha256") != task_entry.get(
        "evaluator_contract_sha256"
    ):
        raise CampaignError("run evaluator contract hash changed")
    if metadata.get("registry_snapshot_sha256") != plan.get("registry_snapshot_sha256"):
        raise CampaignError("run registry snapshot binding changed")
    if metadata.get("campaign_file_sha256") != plan.get("campaign_file", {}).get("sha256"):
        raise CampaignError("run campaign-file binding changed")

    dataset_locks = plan.get("dataset_locks")
    if not isinstance(dataset_locks, Mapping):
        raise CampaignError("plan has no dataset locks")
    if list(run_spec.dataset_ids) != task_entry.get("run_dataset_ids"):
        raise CampaignError("run dataset identifiers differ from task admission")
    expected_dataset_access = _expected_dataset_access(
        run_spec=run_spec,
        task_entry=task_entry,
        dataset_locks=dataset_locks,
    )
    for key, expected_value in expected_dataset_access.items():
        if canonical_hash(metadata.get(key)) != canonical_hash(expected_value):
            raise CampaignError(f"run {key} binding changed")
    scientific = run_spec.stage != "admission"
    for dataset_id in run_spec.dataset_ids:
        lock = dataset_locks.get(dataset_id)
        if not isinstance(lock, Mapping):
            raise CampaignError(f"run names an unlocked dataset: {dataset_id}")
        if run_spec.immutable_inputs.get(dataset_id) != lock.get("activation_sha256"):
            raise CampaignError(f"run dataset activation binding changed: {dataset_id}")
        if scientific and lock.get("activation_ready") is not True:
            raise CampaignError(f"scientific run uses an unactivated dataset: {dataset_id}")
        view = lock.get("dataset_view")
        if view is not None:
            if (
                not isinstance(view, Mapping)
                or view.get("parent_dataset_id") != dataset_id
                or view.get("purpose") != "compatibility_smoke"
                or view.get("allowed_waves") != ["smoke"]
                or run_spec.stage != "smoke"
                or lock.get("activation_source") != "dataset_view"
            ):
                raise CampaignError(
                    f"run dataset view binding is invalid: {dataset_id}"
                )

    planned_prerequisites = plan.get("prerequisite_artifacts", [])
    if not isinstance(planned_prerequisites, list):
        raise CampaignError("plan prerequisite artifact inventory is invalid")
    expected_inputs = [ArtifactRef.from_dict(item) for item in planned_prerequisites]
    dataset_input_owner: dict[str, str] = {}
    raw_environment = model_entry.get("environment_artifact")
    if not isinstance(raw_environment, Mapping):
        raise CampaignError("run model has no frozen environment artifact")
    expected_inputs.append(ArtifactRef.from_dict(raw_environment))
    raw_model_authorities = model_entry.get("authority_artifacts")
    if not isinstance(raw_model_authorities, list):
        raise CampaignError("run model authority artifact inventory is invalid")
    expected_inputs.extend(
        ArtifactRef.from_dict(item) for item in raw_model_authorities
    )
    raw_evaluator_authority = task_entry.get("evaluator_authority")
    if raw_evaluator_authority is not None:
        if not isinstance(raw_evaluator_authority, Mapping):
            raise CampaignError("run task evaluator authority is invalid")
        expected_inputs.append(ArtifactRef.from_dict(raw_evaluator_authority))
    for dataset_id in run_spec.dataset_ids:
        raw_dataset_authorities = dataset_locks[dataset_id].get(
            "authority_artifacts"
        )
        if not isinstance(raw_dataset_authorities, list):
            raise CampaignError(
                f"run dataset authority artifact inventory is invalid: {dataset_id}"
            )
        expected_inputs.extend(
            ArtifactRef.from_dict(item) for item in raw_dataset_authorities
        )
        for item in raw_dataset_authorities:
            artifact = ArtifactRef.from_dict(item)
            assert artifact.role is not None
            dataset_input_owner[artifact.role] = dataset_id
        raw_view_supporting = dataset_locks[dataset_id].get(
            "view_supporting_artifacts", []
        )
        if not isinstance(raw_view_supporting, list):
            raise CampaignError(
                f"run dataset view artifact inventory is invalid: {dataset_id}"
            )
        expected_inputs.extend(
            ArtifactRef.from_dict(item) for item in raw_view_supporting
        )
        for item in raw_view_supporting:
            artifact = ArtifactRef.from_dict(item)
            assert artifact.role is not None
            dataset_input_owner[artifact.role] = dataset_id
    if scientific:
        for dataset_id in run_spec.dataset_ids:
            raw_artifact = dataset_locks[dataset_id].get("artifact_manifest")
            if not isinstance(raw_artifact, Mapping):
                raise CampaignError(
                    f"scientific run dataset has no artifact manifest: {dataset_id}"
                )
            expected_inputs.append(ArtifactRef.from_dict(raw_artifact))
            dataset_artifact = ArtifactRef.from_dict(raw_artifact)
            assert dataset_artifact.role is not None
            dataset_input_owner[dataset_artifact.role] = dataset_id
    if sorted(
        [item.to_dict() for item in expected_inputs], key=lambda item: str(item)
    ) != sorted(
        [item.to_dict() for item in run_spec.inputs], key=lambda item: str(item)
    ):
        raise CampaignError("run input artifacts differ from the reviewed plan")
    expected_input_owner_by_role = {
        str(artifact.role): dataset_input_owner.get(str(artifact.role))
        for artifact in expected_inputs
    }
    if metadata.get("input_owner_by_role") != expected_input_owner_by_role:
        raise CampaignError("run input ownership binding changed")
    if metadata.get("full_input_inventory_sha256") != canonical_hash(
        [artifact.to_dict() for artifact in run_spec.inputs]
    ):
        raise CampaignError("run full input inventory binding changed")
    for index, artifact in enumerate(run_spec.inputs):
        _validate_external_artifact(artifact, f"run input {index}")
    environment_artifact = _environment_artifact(run_spec)
    if environment_artifact.to_dict() != ArtifactRef.from_dict(raw_environment).to_dict():
        raise CampaignError("run environment artifact differs from the model contract")
    if run_spec.checkpoint is not None:
        _validate_external_artifact(run_spec.checkpoint, "run checkpoint")

    repository_value = plan.get("repository_root")
    repository_root = _resolve_existing_directory(
        Path(str(repository_value))
        if repository_value is not None
        else _repository_root(package_root),
        "frozen repository root",
    )
    observed_firewall = assert_resource_unchanged(
        plan["resource_firewall"], repository_root=repository_root
    )
    return repository_root, str(observed_firewall["snapshot_sha256"])


def _failure_state(error: BaseException) -> tuple[RunState, str]:
    if isinstance(error, subprocess.TimeoutExpired):
        return RunState.FAILED_RESOURCE, "resource"
    if isinstance(error, AdapterError):
        text = str(error).lower()
        invalid_markers = ("receipt", "artifact", "inventory", "schema", "run_id")
        if any(marker in text for marker in invalid_markers):
            return RunState.INVALID_OUTPUT, "software"
        return RunState.FAILED_SOFTWARE, "software"
    if isinstance(error, (ArtifactError, ContractError, FirewallError, CampaignError)):
        return RunState.INVALID_OUTPUT, "software"
    return RunState.FAILED_SOFTWARE, "software"


def _attempt_manifest_metadata(
    *,
    run_id: str,
    attempt: int,
    status: RunState,
    plan_sha256: str,
    candidate_manifest_sha256: str,
) -> dict[str, Any]:
    return {
        "artifact_class": "run_execution_attempt",
        "run_id": run_id,
        "attempt": attempt,
        "status": status.value,
        "plan_sha256": plan_sha256,
        "candidate_manifest_sha256": candidate_manifest_sha256,
    }


def _publish_minimal_failure_attempt(
    *,
    working_attempt: Path,
    attempt: Path,
    run_spec: RunSpec,
    receipt: RunExecutionReceipt,
    publication_error: BaseException,
) -> tuple[Path, Path]:
    """Quarantine untrusted output and publish a separately built terminal receipt.

    The untrusted tree is moved atomically without inventorying or copying it.  A
    new minimal tree is then constructed from in-memory contract objects only.
    """

    quarantine = attempt.parent / f"quarantine-{attempt.name}-untrusted"
    quarantine_path = working_attempt
    quarantine_error: str | None = None
    try:
        publish_directory_noreplace(working_attempt, quarantine)
        quarantine_path = quarantine
    except (ArtifactError, OSError) as error:
        quarantine_error = f"{type(error).__name__}: {error}"

    publication_failure = (
        f"attempt_publication_preflight={type(publication_error).__name__}: "
        f"{publication_error}"
    )
    original_failure = receipt.failure or (
        f"original_status={receipt.status.value}; output could not be frozen"
    )
    failure_parts = [original_failure, publication_failure]
    if quarantine_error is not None:
        failure_parts.append(f"quarantine_publication={quarantine_error}")
    terminal_receipt = replace(
        receipt,
        action_receipts=(),
        status=RunState.INVALID_OUTPUT,
        failure="; ".join(failure_parts),
        retry_class="software",
    )

    minimal_staging = Path(
        tempfile.mkdtemp(prefix=f".{attempt.name}.terminal.", dir=attempt.parent)
    )
    write_json_exclusive(minimal_staging / "run_spec.json", run_spec.to_dict())
    write_json_exclusive(
        minimal_staging / "failure_context.json",
        {
            "schema_version": "masld-bench-attempt-publication-failure-v1",
            "run_id": run_spec.run_id,
            "attempt": run_spec.attempt,
            "original_status": receipt.status.value,
            "original_failure": receipt.failure,
            "publication_failure": publication_failure,
            "untrusted_quarantine_path": quarantine_path.as_posix(),
            "untrusted_quarantine_published": quarantine_path == quarantine,
            "quarantine_publication_failure": quarantine_error,
            "discarded_action_receipts": [
                dict(action_receipt)
                for action_receipt in receipt.action_receipts
            ],
        },
    )
    write_json_exclusive(
        minimal_staging / "run_execution_receipt.json",
        terminal_receipt.to_dict(),
    )
    freeze_tree(
        minimal_staging,
        _attempt_manifest_metadata(
            run_id=run_spec.run_id,
            attempt=run_spec.attempt,
            status=terminal_receipt.status,
            plan_sha256=terminal_receipt.plan_sha256,
            candidate_manifest_sha256=(
                terminal_receipt.candidate_manifest_sha256
            ),
        ),
    )
    publish_directory_noreplace(minimal_staging, attempt)
    return attempt, quarantine_path


def execute_run(
    *, candidate: str | Path, run_id: str, output_root: str | Path
) -> Path:
    candidate_path = _resolve_existing_directory(candidate, "candidate")
    _, plan, planned_run, candidate_manifest_sha256 = _candidate_run(
        candidate_path, run_id
    )
    retry_policy = plan.get("retry_policy")
    expected_retry_policy = {
        "max_attempts": 1,
        "retryable_states": [],
        "checkpoint_resume_required": True,
        "software_correction": "new_campaign_revision",
        "reason": (
            "No standardized immutable resume-checkpoint ArtifactRef contract "
            "is implemented; repeated execution fails closed."
        ),
    }
    if retry_policy != expected_retry_policy:
        raise CampaignError("frozen retry policy is missing or unsupported")
    execution = _resolve_output_directory(output_root, "execution output root")
    run_root = execution / "runs" / run_id
    run_root.parent.mkdir(parents=True, exist_ok=True)
    try:
        run_root.mkdir(exist_ok=False)
    except FileExistsError as error:
        raise CampaignError(
            "a prior attempt or interrupted staging directory exists; retry requires a future reviewed "
            "resume-checkpoint contract, while software correction requires a new campaign"
        ) from error
    number, attempt = _next_attempt(run_root)
    if number != 1:
        raise CampaignError("initial run attempt numbering is not canonical")
    working_attempt = Path(
        tempfile.mkdtemp(prefix=f".{attempt.name}.", dir=run_root)
    )
    run_spec = replace(
        planned_run,
        output_dir=attempt.as_posix(),
        attempt=number,
    )
    write_json_exclusive(working_attempt / "run_spec.json", run_spec.to_dict())
    action_receipts: list[dict[str, Any]] = []
    status = RunState.SUCCEEDED
    failure: str | None = None
    retry_class: str | None = None
    firewall_before: str | None = None
    firewall_after: str | None = None
    repository_root: Path | None = None
    try:
        include_r = run_spec.metadata.get("capture_r_session")
        if not isinstance(include_r, bool):
            raise CampaignError("run capture_r_session binding is not boolean")
        capture_runtime_lock(
            working_attempt / "runtime", include_r=include_r, include_cuda=True
        )
        repository_root, firewall_before = _assert_run_bindings(
            plan=plan, run_spec=run_spec
        )
        adapter_request = run_spec.metadata.get("adapter_request", {})
        if not isinstance(adapter_request, Mapping):
            raise CampaignError("run adapter_request metadata is not an object")
        timeout_seconds = run_spec.metadata.get("adapter_timeout_seconds")
        if (
            isinstance(timeout_seconds, bool)
            or not isinstance(timeout_seconds, int)
            or timeout_seconds < 1
        ):
            raise CampaignError("run adapter timeout is not a positive integer")
        adapter = SubprocessAdapter(
            run_spec.adapter_command,
            timeout_seconds=timeout_seconds,
        )
        environment_artifact = _environment_artifact(run_spec)
        raw_fit_dataset_ids = run_spec.metadata.get("fit_dataset_ids")
        if not isinstance(raw_fit_dataset_ids, (list, tuple)) or any(
            not isinstance(item, str) for item in raw_fit_dataset_ids
        ):
            raise CampaignError("run fit_dataset_ids metadata is invalid")
        fit_dataset_ids = tuple(raw_fit_dataset_ids)
        prior_outputs: list[dict[str, Any]] = []
        for index, action_text in enumerate(run_spec.action, start=1):
            action = AdapterAction(action_text)
            output_relative = Path("adapter_actions") / f"{index:03d}-{action.value}"
            output = working_attempt / output_relative
            action_dataset_ids = expected_action_dataset_ids(
                action=action,
                stage=run_spec.stage,
                run_dataset_ids=run_spec.dataset_ids,
                fit_dataset_ids=fit_dataset_ids,
                development_prediction_dataset_ids=run_spec.metadata[
                    "development_prediction_dataset_ids"
                ],
                prediction_first_stress_dataset_ids=run_spec.metadata[
                    "prediction_first_stress_dataset_ids"
                ],
                sealed_prediction_dataset_ids=run_spec.metadata[
                    "sealed_prediction_dataset_ids"
                ],
            )
            owners = run_spec.metadata["input_owner_by_role"]
            if not isinstance(owners, Mapping):
                raise CampaignError("run input ownership binding is invalid")
            allowed_dataset_ids = frozenset(action_dataset_ids)
            available_roles = frozenset(str(role) for role in owners)
            withheld_roles = withheld_input_roles_for_action(
                request=adapter_request,
                action=action,
                available_roles=available_roles,
            )
            scoped_artifacts = [
                artifact
                for artifact in run_spec.inputs
                if (
                    owners.get(str(artifact.role)) is None
                    or owners.get(str(artifact.role)) in allowed_dataset_ids
                )
                and str(artifact.role) not in withheld_roles
            ]
            scoped_inputs = [artifact.to_dict() for artifact in scoped_artifacts]
            scoped_run_spec = run_spec.to_dict()
            scoped_run_spec["inputs"] = scoped_inputs
            request = dict(adapter_request)
            request.update(
                {
                    "schema_version": "masld-bench-adapter-request-v1",
                    "run_id": run_id,
                    "action": action.value,
                    "run_spec": scoped_run_spec,
                    "prior_action_outputs": prior_outputs,
                    "fit_dataset_ids": list(fit_dataset_ids),
                    "development_prediction_dataset_ids": list(
                        run_spec.metadata["development_prediction_dataset_ids"]
                    ),
                    "prediction_first_stress_dataset_ids": list(
                        run_spec.metadata["prediction_first_stress_dataset_ids"]
                    ),
                    "sealed_prediction_dataset_ids": list(
                        run_spec.metadata["sealed_prediction_dataset_ids"]
                    ),
                    "dataset_routing": run_spec.metadata["dataset_routing"],
                    "action_dataset_ids": list(action_dataset_ids),
                    "full_input_inventory_sha256": run_spec.metadata[
                        "full_input_inventory_sha256"
                    ],
                    "action_input_inventory_sha256": canonical_hash(scoped_inputs),
                }
            )
            request_path = (
                working_attempt / "requests" / f"{index:03d}-{action.value}.json"
            )
            write_json_exclusive(request_path, request)
            for input_index, artifact in enumerate(scoped_artifacts):
                _validate_external_artifact(
                    artifact,
                    f"{action.value} scoped input {input_index} before execution",
                )
            if run_spec.checkpoint is not None:
                _validate_external_artifact(
                    run_spec.checkpoint,
                    f"{action.value} checkpoint before execution",
                )
            try:
                adapter_receipt = adapter.invoke(action, request_path, output)
            finally:
                for input_index, artifact in enumerate(scoped_artifacts):
                    _validate_external_artifact(
                        artifact,
                        f"{action.value} scoped input {input_index} after execution",
                    )
                if run_spec.checkpoint is not None:
                    _validate_external_artifact(
                        run_spec.checkpoint,
                        f"{action.value} checkpoint after execution",
                    )
            if adapter_receipt.status != "complete":
                raise AdapterError(
                    f"adapter action {action.value} returned {adapter_receipt.status}"
                )
            _assert_adapter_environment_receipt(
                adapter_receipt,
                environment_artifact=environment_artifact,
                runtime_id=run_spec.runtime_id,
                fit_dataset_ids=fit_dataset_ids,
            )
            output_manifest_sha256 = sha256_file(output / "ARTIFACTS.json")
            adapter_receipt_sha256 = sha256_file(output / "adapter_receipt.json")
            action_record = {
                "action": action.value,
                "status": adapter_receipt.status,
                "output_path": output_relative.as_posix(),
                "adapter_receipt_sha256": adapter_receipt_sha256,
                "output_manifest_sha256": output_manifest_sha256,
            }
            action_receipts.append(action_record)
            prior_outputs.append(dict(action_record))

        _assert_config_snapshot(plan)
        if source_lock(Path(str(plan["package_root"]))) != plan["source_lock"]:
            raise CampaignError("source tree changed during adapter execution")
        for index, artifact in enumerate(run_spec.inputs):
            _validate_external_artifact(artifact, f"run input {index} after execution")
        if run_spec.checkpoint is not None:
            _validate_external_artifact(run_spec.checkpoint, "run checkpoint after execution")
        assert repository_root is not None
        observed_after = assert_resource_unchanged(
            plan["resource_firewall"], repository_root=repository_root
        )
        firewall_after = str(observed_after["snapshot_sha256"])
    except Exception as error:  # freeze a terminal receipt for every begun attempt
        status, retry_class = _failure_state(error)
        failure = f"{type(error).__name__}: {error}"
        if repository_root is not None and firewall_before is not None:
            try:
                observed_after = assert_resource_unchanged(
                    plan["resource_firewall"], repository_root=repository_root
                )
                firewall_after = str(observed_after["snapshot_sha256"])
            except Exception as post_error:
                failure += f"; postflight={type(post_error).__name__}: {post_error}"
                status = RunState.INVALID_OUTPUT
    receipt = RunExecutionReceipt(
        schema_version="masld-bench-run-execution-receipt-v1",
        candidate_path=candidate_path.as_posix(),
        campaign_id=run_spec.campaign_id,
        plan_sha256=str(plan["plan_sha256"]),
        candidate_manifest_sha256=candidate_manifest_sha256,
        run_id=run_id,
        attempt=number,
        wave=run_spec.stage,
        task_id=run_spec.task_id,
        model_id=run_spec.model_id,
        seed=run_spec.seed,
        fold=run_spec.fold,
        adaptation_regime=run_spec.adaptation_regime,
        runtime_id=run_spec.runtime_id,
        runtime_registry_sha256=run_spec.runtime_registry_sha256,
        adapter_actions=run_spec.action,
        action_receipts=tuple(action_receipts),
        resource_firewall_before_sha256=firewall_before,
        resource_firewall_after_sha256=firewall_after,
        status=status,
        failure=failure,
        retry_class=retry_class,
    )
    try:
        write_json_exclusive(
            working_attempt / "run_execution_receipt.json", receipt.to_dict()
        )
        freeze_tree(
            working_attempt,
            _attempt_manifest_metadata(
                run_id=run_id,
                attempt=number,
                status=status,
                plan_sha256=str(plan["plan_sha256"]),
                candidate_manifest_sha256=candidate_manifest_sha256,
            ),
        )
        publish_directory_noreplace(working_attempt, attempt)
    except (ArtifactError, OSError) as error:
        try:
            terminal_attempt, quarantine = _publish_minimal_failure_attempt(
                working_attempt=working_attempt,
                attempt=attempt,
                run_spec=run_spec,
                receipt=receipt,
                publication_error=error,
            )
        except (ArtifactError, OSError, ContractError) as terminal_error:
            raise CampaignError(
                "run-attempt publication failed and the minimal terminal receipt "
                f"could not be published; untrusted staging retained at "
                f"{working_attempt}: {terminal_error}"
            ) from terminal_error
        raise CampaignError(
            "run output failed immutable publication; minimal terminal receipt: "
            f"{terminal_attempt}; untrusted output retained at {quarantine}"
        ) from error
    if status is not RunState.SUCCEEDED:
        raise CampaignError(f"run failed; immutable receipt: {attempt}")
    return attempt


def execute_admission_run(
    *, candidate: str | Path, run_id: str, output_root: str | Path
) -> Path:
    """Compatibility wrapper that refuses to execute a non-admission run."""

    candidate_path = _resolve_existing_directory(candidate, "candidate")
    _, _, run_spec, _ = _candidate_run(candidate_path, run_id)
    if run_spec.stage != "admission":
        raise CampaignError("execute_admission_run may run only the admission wave")
    return execute_run(candidate=candidate_path, run_id=run_id, output_root=output_root)


def verify_run_execution_attempt(
    path: str | Path, *, require_succeeded: bool = True
) -> dict[str, Any]:
    attempt = _resolve_existing_directory(path, "run execution attempt")
    try:
        manifest = verify_frozen_tree(attempt)
    except ArtifactError as error:
        raise CampaignError(f"run attempt is not an intact immutable tree: {error}") from error
    metadata = manifest.get("metadata")
    expected_metadata_fields = {
        "artifact_class",
        "run_id",
        "attempt",
        "status",
        "plan_sha256",
        "candidate_manifest_sha256",
    }
    if not isinstance(metadata, Mapping) or set(metadata) != expected_metadata_fields:
        raise CampaignError("run-attempt manifest metadata is incomplete or ambiguous")
    if metadata.get("artifact_class") != "run_execution_attempt":
        raise CampaignError("artifact tree is not a run execution attempt")
    try:
        receipt = RunExecutionReceipt.from_dict(
            _read_json_object(
                attempt / "run_execution_receipt.json", "run execution receipt"
            )
        )
    except ContractError as error:
        raise CampaignError(f"invalid run execution receipt: {error}") from error
    candidate = Path(receipt.candidate_path)
    _, plan, planned_run, candidate_manifest_sha256 = _candidate_run(
        candidate, receipt.run_id
    )
    if candidate_manifest_sha256 != receipt.candidate_manifest_sha256:
        raise CampaignError("execution receipt candidate manifest binding changed")
    if receipt.plan_sha256 != plan.get("plan_sha256"):
        raise CampaignError("execution receipt plan binding changed")
    firewall_sha256 = plan.get("resource_firewall", {}).get("snapshot_sha256")
    if (
        receipt.resource_firewall_before_sha256 is not None
        and receipt.resource_firewall_before_sha256 != firewall_sha256
    ) or (
        receipt.resource_firewall_after_sha256 is not None
        and receipt.resource_firewall_after_sha256 != firewall_sha256
    ):
        raise CampaignError("execution receipt Resource firewall binding changed")
    if metadata != {
        "artifact_class": "run_execution_attempt",
        "run_id": receipt.run_id,
        "attempt": receipt.attempt,
        "status": receipt.status.value,
        "plan_sha256": receipt.plan_sha256,
        "candidate_manifest_sha256": receipt.candidate_manifest_sha256,
    }:
        raise CampaignError("run-attempt manifest metadata differs from its receipt")
    if attempt.name != f"attempt-{receipt.attempt:03d}":
        raise CampaignError("run attempt directory name differs from its receipt")
    try:
        executed_run = RunSpec.from_dict(
            _read_json_object(attempt / "run_spec.json", "executed RunSpec")
        )
    except ContractError as error:
        raise CampaignError(f"invalid executed RunSpec: {error}") from error
    if executed_run.run_id != receipt.run_id:
        raise CampaignError("executed RunSpec run identity differs from its receipt")
    if executed_run.identity_payload != planned_run.identity_payload:
        raise CampaignError("executed RunSpec differs from the frozen candidate run")
    if executed_run.output_dir != attempt.as_posix() or executed_run.attempt != receipt.attempt:
        raise CampaignError("executed RunSpec attempt location/number is inconsistent")
    for index, artifact in enumerate(executed_run.inputs):
        _validate_external_artifact(artifact, f"executed run input {index}")
    environment_artifact = _environment_artifact(executed_run)
    if executed_run.checkpoint is not None:
        _validate_external_artifact(executed_run.checkpoint, "executed run checkpoint")
    receipt_bindings = {
        "campaign_id": executed_run.campaign_id,
        "wave": executed_run.stage,
        "task_id": executed_run.task_id,
        "model_id": executed_run.model_id,
        "seed": executed_run.seed,
        "fold": executed_run.fold,
        "adaptation_regime": executed_run.adaptation_regime,
        "runtime_id": executed_run.runtime_id,
        "runtime_registry_sha256": executed_run.runtime_registry_sha256,
        "adapter_actions": list(executed_run.action),
    }
    receipt_payload = receipt.to_dict()
    for key, expected in receipt_bindings.items():
        if receipt_payload.get(key) != expected:
            raise CampaignError(f"execution receipt {key} binding changed")
    for action_record in receipt.action_receipts:
        output = attempt / str(action_record["output_path"])
        try:
            output_manifest = verify_frozen_tree(output)
        except ArtifactError as error:
            raise CampaignError(
                f"adapter action output is not an intact frozen tree: {error}"
            ) from error
        if sha256_file(output / "ARTIFACTS.json") != action_record[
            "output_manifest_sha256"
        ]:
            raise CampaignError("adapter action output manifest binding changed")
        if sha256_file(output / "adapter_receipt.json") != action_record[
            "adapter_receipt_sha256"
        ]:
            raise CampaignError("adapter action receipt binding changed")
        try:
            adapter_receipt = AdapterReceipt.from_mapping(
                _read_json_object(output / "adapter_receipt.json", "adapter receipt")
            )
        except AdapterError as error:
            raise CampaignError(f"invalid adapter action receipt: {error}") from error
        if (
            adapter_receipt.run_id != receipt.run_id
            or adapter_receipt.action.value != action_record["action"]
            or adapter_receipt.status != "complete"
        ):
            raise CampaignError("adapter action receipt identity changed")
        _assert_adapter_environment_receipt(
            adapter_receipt,
            environment_artifact=environment_artifact,
            runtime_id=executed_run.runtime_id,
            fit_dataset_ids=tuple(
                str(item)
                for item in executed_run.metadata.get("fit_dataset_ids", ())
            ),
        )
        adapter_receipt.validate_artifacts(output)
        output_metadata = output_manifest.get("metadata")
        if output_metadata != {
            "artifact_class": "adapter_output",
            "action": action_record["action"],
            "run_id": receipt.run_id,
        }:
            raise CampaignError("adapter output manifest metadata changed")
    if require_succeeded and receipt.status is not RunState.SUCCEEDED:
        raise CampaignError(
            f"run execution attempt is terminal but not successful: {receipt.status.value}"
        )
    return receipt_payload


def _repository_root(path: Path) -> Path:
    process = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        cwd=path,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if process.returncode != 0:
        raise CampaignError("cannot resolve repository root for Resource firewall")
    return Path(process.stdout.strip())


def submission_commands(
    *,
    candidate: str | Path,
    approved_campaign_sha256: str,
    execute: bool,
    submission_ledger_root: str | Path | None = None,
) -> list[str]:
    """Return or execute the exact scripts covered by the reviewed tree hash."""

    candidate_path = _resolve_existing_directory(candidate, "candidate")
    try:
        manifest = verify_frozen_tree(candidate_path)
    except ArtifactError as error:
        raise CampaignError(f"candidate is not an intact immutable tree: {error}") from error
    plan = load_frozen_plan(candidate_path)
    campaign_sha256 = sha256_file(candidate_path / "ARTIFACTS.json")
    if approved_campaign_sha256 != campaign_sha256:
        raise CampaignError(
            "approved campaign SHA-256 does not match the frozen candidate tree"
        )
    if plan["safety"].get("arrays"):
        raise CampaignError("array submission is not implemented or authorized")
    if execute and not plan["safety"].get("submit_enabled"):
        raise CampaignError("this frozen campaign has submission disabled")
    expected_relative = {
        f"jobs/{run['run_id']}.sbatch" for run in plan.get("runs", [])
    }
    manifested_relative = {
        str(item["path"])
        for item in manifest.get("artifacts", [])
        if str(item.get("path", "")).startswith("jobs/")
        and str(item.get("path", "")).endswith(".sbatch")
    }
    if manifested_relative != expected_relative:
        raise CampaignError(
            "frozen job script set does not exactly match planned run identifiers"
        )
    scripts = [candidate_path / relative for relative in sorted(expected_relative)]
    commands = [f"sbatch --parsable {path.as_posix()}" for path in scripts]
    if execute:
        if submission_ledger_root is None:
            raise CampaignError(
                "production submission requires a stable submission_ledger_root"
            )
        ledger_root = _resolve_output_directory(
            submission_ledger_root, "submission ledger root"
        )
        ledger_root.mkdir(parents=True, exist_ok=True)
        campaign_id = str(plan["campaign"]["campaign_id"])
        ledger = ledger_root / f"{campaign_id}--{campaign_sha256[:16]}"
        claim = ledger_root / f".{ledger.name}.submission-claim.json"
        try:
            write_json_exclusive(
                claim,
                {
                    "schema_version": "masld-bench-submission-claim-v1",
                    "campaign_id": campaign_id,
                    "candidate_manifest_sha256": campaign_sha256,
                    "target": ledger.as_posix(),
                },
            )
        except ArtifactError as error:
            raise CampaignError(
                "this frozen campaign already has an initial-submission ledger"
            ) from error
        staging_ledger = Path(
            tempfile.mkdtemp(prefix=f".{ledger.name}.", dir=ledger_root)
        )
        write_json_exclusive(
            staging_ledger / "submission_request.json",
            {
                "schema_version": "masld-bench-submission-request-v1",
                "campaign_id": campaign_id,
                "candidate_path": candidate_path.as_posix(),
                "candidate_manifest_sha256": campaign_sha256,
                "approved_campaign_sha256": approved_campaign_sha256,
                "commands": commands,
            },
        )
        execution_root = candidate_path.parent / "executions" / candidate_path.name
        (execution_root / "logs").mkdir(parents=True, exist_ok=True)
        records: list[dict[str, Any]] = []
        failure: str | None = None
        for index, path in enumerate(scripts, start=1):
            command = ["sbatch", "--parsable", path.as_posix()]
            try:
                process = subprocess.run(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=False,
                )
                returncode: int | None = process.returncode
                stdout = process.stdout
                stderr = process.stderr
            except OSError as error:
                returncode = None
                stdout = ""
                stderr = f"{type(error).__name__}: {error}"
            job_id: str | None = None
            if returncode == 0:
                candidate_job_id = stdout.strip().split(";", 1)[0]
                if re.fullmatch(r"[0-9]+", candidate_job_id):
                    job_id = candidate_job_id
                else:
                    failure = (
                        f"sbatch returned an invalid parsable job ID for {path.name}: "
                        f"{stdout.strip()!r}"
                    )
            else:
                failure = (
                    f"sbatch failed for {path.name} with exit {returncode}: "
                    f"{stderr.strip()}"
                )
            record = {
                "index": index,
                "script": path.as_posix(),
                "command": command,
                "returncode": returncode,
                "stdout": stdout,
                "stderr": stderr,
                "job_id": job_id,
                "status": "submitted" if failure is None else "failed",
            }
            records.append(record)
            write_json_exclusive(
                staging_ledger / "records" / f"{index:05d}.json", record
            )
            if failure is not None:
                break
        status = "succeeded" if failure is None else "partial_failure"
        write_json_exclusive(
            staging_ledger / "submission_receipt.json",
            {
                "schema_version": "masld-bench-submission-receipt-v1",
                "campaign_id": campaign_id,
                "candidate_manifest_sha256": campaign_sha256,
                "status": status,
                "submitted_job_ids": [
                    record["job_id"]
                    for record in records
                    if record["job_id"] is not None
                ],
                "records_written": len(records),
                "planned_jobs": len(scripts),
                "failure": failure,
            },
        )
        freeze_tree(
            staging_ledger,
            {
                "artifact_class": "campaign_submission_ledger",
                "campaign_id": campaign_id,
                "candidate_manifest_sha256": campaign_sha256,
                "status": status,
            },
        )
        try:
            publish_directory_noreplace(staging_ledger, ledger)
        except (ArtifactError, OSError) as error:
            raise CampaignError(
                "submission-ledger publication failed; claim and staging are "
                f"retained at {claim} and {staging_ledger}: {error}"
            ) from error
        if failure is not None:
            raise CampaignError(f"submission stopped; immutable ledger: {ledger}")
    return commands
