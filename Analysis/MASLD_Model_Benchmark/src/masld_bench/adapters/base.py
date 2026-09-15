"""Environment-neutral adapter interface and subprocess implementation."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import StrEnum
import json
from pathlib import Path
import re
import subprocess
from typing import Any, Mapping, Sequence

from ..artifacts import (
    canonical_hash,
    freeze_tree,
    reject_symlink_components,
    sha256_file,
    verify_frozen_tree,
    write_text_exclusive,
)


class AdapterError(RuntimeError):
    """Raised when a model adapter violates the exchange protocol."""


class AdapterAction(StrEnum):
    PROBE = "probe"
    PREPARE = "prepare"
    FIT = "fit"
    PREDICT = "predict"
    EXPORT = "export"


_FORBIDDEN_EVALUATION_KEYS = frozenset(
    {
        "metrics",
        "benchmark_metrics",
        "macro_f1",
        "auprc",
        "brier",
        "p_value",
        "q_value",
        "promotion_decision",
    }
)
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_RECEIPT_FIELDS = frozenset(
    {"schema_version", "action", "run_id", "status", "artifacts", "metadata"}
)
_RECEIPT_REQUIRED = _RECEIPT_FIELDS - {"metadata"}
_ARTIFACT_FIELDS = frozenset({"path", "sha256", "size_bytes"})
_PRIOR_OUTPUT_FIELDS = frozenset(
    {
        "action",
        "status",
        "output_path",
        "adapter_receipt_sha256",
        "output_manifest_sha256",
    }
)


def _reject_evaluation_fields(value: Any, path: str = "receipt") -> None:
    if isinstance(value, Mapping):
        forbidden = _FORBIDDEN_EVALUATION_KEYS.intersection(str(key) for key in value)
        if forbidden:
            raise AdapterError(
                f"adapters may not return benchmark evaluation fields at {path}: "
                + ", ".join(sorted(forbidden))
            )
        for key, child in value.items():
            _reject_evaluation_fields(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _reject_evaluation_fields(child, f"{path}[{index}]")


def _validate_dataset_access_request(
    request: Mapping[str, Any], action: AdapterAction
) -> None:
    run_spec = request.get("run_spec")
    if not isinstance(run_spec, Mapping):
        raise AdapterError("adapter request must bind a RunSpec object")
    raw_run_ids = run_spec.get("dataset_ids")
    if not isinstance(raw_run_ids, list) or any(
        not isinstance(item, str) for item in raw_run_ids
    ):
        raise AdapterError("adapter RunSpec dataset_ids must be a string array")
    run_ids = tuple(raw_run_ids)
    list_fields = (
        "fit_dataset_ids",
        "development_prediction_dataset_ids",
        "prediction_first_stress_dataset_ids",
        "sealed_prediction_dataset_ids",
        "action_dataset_ids",
    )
    partitions: dict[str, tuple[str, ...]] = {}
    for field_name in list_fields:
        raw = request.get(field_name)
        if not isinstance(raw, list) or any(not isinstance(item, str) for item in raw):
            raise AdapterError(f"adapter request {field_name} must be a string array")
        values = tuple(raw)
        if len(set(values)) != len(values):
            raise AdapterError(f"adapter request {field_name} contains duplicates")
        partitions[field_name] = values
    partition_ids = (
        *partitions["fit_dataset_ids"],
        *partitions["development_prediction_dataset_ids"],
        *partitions["prediction_first_stress_dataset_ids"],
        *partitions["sealed_prediction_dataset_ids"],
    )
    if len(set(partition_ids)) != len(partition_ids) or set(partition_ids) != set(
        run_ids
    ):
        raise AdapterError(
            "adapter request dataset partitions must uniquely cover RunSpec.dataset_ids"
        )
    routing = request.get("dataset_routing")
    if not isinstance(routing, Mapping) or set(routing) != set(run_ids):
        raise AdapterError("adapter request dataset_routing must cover every run dataset")
    expected_action_ids = expected_action_dataset_ids(
        action=action,
        stage=str(run_spec.get("stage", "")),
        run_dataset_ids=run_ids,
        fit_dataset_ids=partitions["fit_dataset_ids"],
        development_prediction_dataset_ids=partitions[
            "development_prediction_dataset_ids"
        ],
        prediction_first_stress_dataset_ids=partitions[
            "prediction_first_stress_dataset_ids"
        ],
        sealed_prediction_dataset_ids=partitions[
            "sealed_prediction_dataset_ids"
        ],
    )
    if partitions["action_dataset_ids"] != expected_action_ids:
        raise AdapterError(
            f"{action.value} action_dataset_ids differ from the frozen action scope"
        )
    if action is AdapterAction.FIT:
        for dataset_id in partitions["action_dataset_ids"]:
            route = routing.get(dataset_id)
            if (
                not isinstance(route, Mapping)
                or route.get("task_partition") != "fit"
                or route.get("fit_allowed") is not True
                or route.get("prediction_first_policy") is not None
            ):
                raise AdapterError(
                    f"fit action received a non-fit or prediction-first dataset: {dataset_id}"
                )
    _validate_scoped_input_inventory(
        request=request,
        run_spec=run_spec,
        action_dataset_ids=expected_action_ids,
        action=action,
    )


def expected_action_dataset_ids(
    *,
    action: AdapterAction,
    stage: str,
    run_dataset_ids: Sequence[str],
    fit_dataset_ids: Sequence[str],
    development_prediction_dataset_ids: Sequence[str],
    prediction_first_stress_dataset_ids: Sequence[str],
    sealed_prediction_dataset_ids: Sequence[str],
) -> tuple[str, ...]:
    """Return the only dataset roster an adapter action may receive."""

    run_ids = tuple(run_dataset_ids)
    fit = frozenset(fit_dataset_ids)
    prediction = frozenset(
        (
            *development_prediction_dataset_ids,
            *prediction_first_stress_dataset_ids,
            *sealed_prediction_dataset_ids,
        )
    )
    if action is AdapterAction.PROBE:
        allowed = frozenset(run_ids)
    elif action is AdapterAction.PREPARE:
        allowed = fit if fit else prediction
    elif action is AdapterAction.FIT:
        allowed = fit
    elif action is AdapterAction.PREDICT:
        if prediction:
            allowed = prediction
        elif stage == "smoke":
            allowed = fit
        else:
            allowed = frozenset()
    elif action is AdapterAction.EXPORT:
        allowed = frozenset()
    else:  # pragma: no cover - AdapterAction is exhaustive
        raise AdapterError(f"unsupported adapter action: {action}")
    return tuple(dataset_id for dataset_id in run_ids if dataset_id in allowed)


def _validate_scoped_input_inventory(
    *,
    request: Mapping[str, Any],
    run_spec: Mapping[str, Any],
    action_dataset_ids: Sequence[str],
    action: AdapterAction,
) -> None:
    inputs = run_spec.get("inputs")
    metadata = run_spec.get("metadata")
    if not isinstance(inputs, list) or any(not isinstance(item, Mapping) for item in inputs):
        raise AdapterError("adapter RunSpec inputs must be an object array")
    if not isinstance(metadata, Mapping):
        raise AdapterError("adapter RunSpec metadata must be an object")
    owners = metadata.get("input_owner_by_role")
    if not isinstance(owners, Mapping) or any(
        not isinstance(role, str)
        or (owner is not None and not isinstance(owner, str))
        for role, owner in owners.items()
    ):
        raise AdapterError("adapter RunSpec input_owner_by_role is invalid")
    run_ids = run_spec.get("dataset_ids")
    if not isinstance(run_ids, list):
        raise AdapterError("adapter RunSpec dataset_ids must be a string array")
    if any(owner is not None and owner not in run_ids for owner in owners.values()):
        raise AdapterError("input_owner_by_role names a dataset outside the RunSpec")
    roles = [item.get("role") for item in inputs]
    if any(not isinstance(role, str) for role in roles) or len(set(roles)) != len(roles):
        raise AdapterError("scoped RunSpec inputs require unique string roles")
    allowed = frozenset(action_dataset_ids)
    expected_roles = {
        role for role, owner in owners.items() if owner is None or owner in allowed
    }
    expected_roles.difference_update(
        withheld_input_roles_for_action(
            request=request,
            action=action,
            available_roles=frozenset(owners),
        )
    )
    if set(roles) != expected_roles:
        raise AdapterError("scoped RunSpec inputs expose an incorrect action inventory")
    full_hash = metadata.get("full_input_inventory_sha256")
    if not isinstance(full_hash, str) or not _SHA256.fullmatch(full_hash):
        raise AdapterError("RunSpec full input inventory hash is invalid")
    if request.get("full_input_inventory_sha256") != full_hash:
        raise AdapterError("adapter request full input inventory binding changed")
    action_hash = request.get("action_input_inventory_sha256")
    if action_hash != canonical_hash(inputs):
        raise AdapterError("adapter request action input inventory binding changed")


def withheld_input_roles_for_action(
    *,
    request: Mapping[str, Any],
    action: AdapterAction,
    available_roles: frozenset[str],
) -> frozenset[str]:
    """Return prospectively frozen roles hidden from one adapter action.

    This supports outcome-separated workflows in which ``prepare`` may read a
    paired development object and materialize training targets plus query
    inputs, while later ``fit`` and ``predict`` requests must not receive the
    original outcome-bearing object.  The omission is part of the reviewed
    adapter request and therefore enters the campaign hash.
    """

    raw = request.get("withheld_input_roles_by_action", {})
    if not isinstance(raw, Mapping):
        raise AdapterError("withheld_input_roles_by_action must be an object")
    unknown_actions = sorted(
        str(key) for key in raw if str(key) not in {item.value for item in AdapterAction}
    )
    if unknown_actions:
        raise AdapterError(
            "withheld_input_roles_by_action names unsupported actions: "
            + ", ".join(unknown_actions)
        )
    raw_roles = raw.get(action.value, [])
    if not isinstance(raw_roles, (list, tuple)) or any(
        not isinstance(role, str) or not role for role in raw_roles
    ):
        raise AdapterError(
            f"withheld_input_roles_by_action.{action.value} must be a string array"
        )
    roles = tuple(raw_roles)
    if len(set(roles)) != len(roles) or roles != tuple(sorted(roles)):
        raise AdapterError(
            f"withheld_input_roles_by_action.{action.value} must be sorted and unique"
        )
    unknown_roles = sorted(set(roles).difference(available_roles))
    if unknown_roles:
        raise AdapterError(
            "withheld input roles are absent from the full RunSpec inventory: "
            + ", ".join(unknown_roles)
        )
    environment_roles = sorted(
        role for role in roles if role.startswith("environment:")
    )
    if environment_roles:
        raise AdapterError("adapter actions may not withhold their environment lock")
    return frozenset(roles)


def _stable_file_snapshot(path: Path, label: str) -> tuple[int, int, int, int, int, str]:
    try:
        lexical = reject_symlink_components(path, label=label)
        if not lexical.is_file():
            raise AdapterError(f"{label} is not a regular file: {lexical}")
        before = lexical.stat()
        digest = sha256_file(lexical)
        after = lexical.stat()
        reject_symlink_components(lexical, label=label)
    except AdapterError:
        raise
    except (OSError, RuntimeError) as error:
        raise AdapterError(f"cannot bind {label}: {error}") from error
    identity = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
        digest,
    )
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    ) != identity[:-1]:
        raise AdapterError(f"{label} changed while it was hashed")
    return identity


def _validate_output_inventory(
    receipt: "AdapterReceipt", output: Path, *, frozen: bool = False
) -> None:
    declared_paths = {
        Path(str(artifact["path"])).as_posix() for artifact in receipt.artifacts
    }
    observed_paths = {
        path.relative_to(output).as_posix()
        for path in output.rglob("*")
        if path.is_file() or path.is_symlink()
    }
    expected_paths = declared_paths | {
        "adapter_receipt.json",
        "adapter_stdout.log",
        "adapter_stderr.log",
    }
    if frozen:
        expected_paths.update({"ARTIFACTS.json", "COMPLETE"})
    if observed_paths != expected_paths:
        unexpected = sorted(observed_paths - expected_paths)
        missing = sorted(expected_paths - observed_paths)
        details = []
        if unexpected:
            details.append("unexpected=" + ",".join(unexpected))
        if missing:
            details.append("missing=" + ",".join(missing))
        raise AdapterError(
            "adapter output inventory does not match its receipt: "
            + "; ".join(details)
        )


def _validate_prior_action_outputs(
    request: Path, request_payload: Mapping[str, Any]
) -> None:
    raw_prior = request_payload.get("prior_action_outputs", [])
    if not isinstance(raw_prior, list) or any(
        not isinstance(item, Mapping) for item in raw_prior
    ):
        raise AdapterError("prior_action_outputs must be an object array")
    if not raw_prior:
        return
    run_spec = request_payload.get("run_spec")
    if not isinstance(run_spec, Mapping):
        raise AdapterError("prior outputs require a bound RunSpec")
    actions = run_spec.get("action")
    if not isinstance(actions, list) or any(not isinstance(item, str) for item in actions):
        raise AdapterError("prior outputs require a RunSpec action array")
    current_action = str(request_payload.get("action", ""))
    current_index = len(raw_prior) + 1
    if (
        len(actions) < current_index
        or actions[current_index - 1] != current_action
        or [str(item.get("action", "")) for item in raw_prior]
        != actions[: current_index - 1]
    ):
        raise AdapterError("prior action sequence differs from the RunSpec")
    if (
        request.parent.name != "requests"
        or request.name != f"{current_index:03d}-{current_action}.json"
    ):
        raise AdapterError("request with prior outputs is outside the canonical layout")
    try:
        attempt_root = reject_symlink_components(
            request.parent.parent, label="adapter attempt root"
        ).resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise AdapterError("adapter attempt root is invalid") from error
    if not attempt_root.is_dir():
        raise AdapterError("adapter attempt root is not a directory")
    expected_run_id = str(request_payload.get("run_id", ""))
    for index, item in enumerate(raw_prior, start=1):
        if set(item) != _PRIOR_OUTPUT_FIELDS:
            raise AdapterError(
                f"prior action output {index} has an invalid field inventory"
            )
        action_text = str(item["action"])
        if item["status"] != "complete":
            raise AdapterError(f"prior action output {index} is not complete")
        expected_relative = Path("adapter_actions") / f"{index:03d}-{action_text}"
        relative = Path(str(item["output_path"]))
        if relative != expected_relative or relative.is_absolute() or ".." in relative.parts:
            raise AdapterError(f"prior action output {index} path is not canonical")
        try:
            output = reject_symlink_components(
                attempt_root / relative, label=f"prior action output {index}"
            ).resolve(strict=True)
            output.relative_to(attempt_root)
            manifest = verify_frozen_tree(output)
        except (OSError, RuntimeError, ValueError) as error:
            raise AdapterError(
                f"prior action output {index} is not a verified frozen tree: {error}"
            ) from error
        manifest_path = output / "ARTIFACTS.json"
        receipt_path = output / "adapter_receipt.json"
        if (
            item["output_manifest_sha256"] != sha256_file(manifest_path)
            or item["adapter_receipt_sha256"] != sha256_file(receipt_path)
        ):
            raise AdapterError(f"prior action output {index} binding changed")
        if manifest.get("metadata") != {
            "artifact_class": "adapter_output",
            "action": action_text,
            "run_id": expected_run_id,
        }:
            raise AdapterError(f"prior action output {index} metadata changed")
        try:
            receipt_raw = json.loads(receipt_path.read_text(encoding="utf-8"))
            receipt = AdapterReceipt.from_mapping(receipt_raw)
        except (OSError, json.JSONDecodeError) as error:
            raise AdapterError(
                f"prior action output {index} receipt is invalid: {error}"
            ) from error
        if (
            receipt.action.value != action_text
            or receipt.run_id != expected_run_id
            or receipt.status != "complete"
        ):
            raise AdapterError(f"prior action output {index} receipt binding changed")
        receipt.validate_artifacts(output)
        _validate_output_inventory(receipt, output, frozen=True)


@dataclass(frozen=True, slots=True)
class AdapterReceipt:
    schema_version: str
    action: AdapterAction
    run_id: str
    status: str
    artifacts: tuple[Mapping[str, Any], ...]
    metadata: Mapping[str, Any]

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "AdapterReceipt":
        if not isinstance(value, Mapping):
            raise AdapterError("adapter receipt must be an object")
        unknown = sorted(set(value).difference(_RECEIPT_FIELDS))
        missing = sorted(_RECEIPT_REQUIRED.difference(value))
        if unknown:
            raise AdapterError(
                "adapter receipt has unknown fields: " + ", ".join(unknown)
            )
        if missing:
            raise AdapterError(
                "adapter receipt is missing fields: " + ", ".join(missing)
            )
        _reject_evaluation_fields(value)
        if value.get("schema_version") != "masld-bench-adapter-receipt-v1":
            raise AdapterError("unsupported adapter receipt schema")
        status = str(value.get("status", ""))
        if status not in {"complete", "failed", "ineligible"}:
            raise AdapterError(f"invalid adapter receipt status: {status}")
        raw_artifacts = value.get("artifacts", [])
        if not isinstance(raw_artifacts, list):
            raise AdapterError("adapter receipt artifacts must be a list")
        if any(not isinstance(item, Mapping) for item in raw_artifacts):
            raise AdapterError("adapter receipt artifacts must contain objects")
        for index, item in enumerate(raw_artifacts):
            if set(item) != _ARTIFACT_FIELDS:
                raise AdapterError(
                    f"adapter receipt artifact {index} must contain exactly: "
                    + ", ".join(sorted(_ARTIFACT_FIELDS))
                )
        artifact_paths = [str(item["path"]) for item in raw_artifacts]
        if len(set(artifact_paths)) != len(artifact_paths):
            raise AdapterError("adapter receipt artifact paths must be unique")
        metadata = value.get("metadata", {})
        if not isinstance(metadata, Mapping):
            raise AdapterError("adapter receipt metadata must be an object")
        run_id = str(value.get("run_id", ""))
        if not _SHA256.fullmatch(run_id):
            raise AdapterError("adapter receipt run_id must be a full lowercase SHA-256")
        try:
            action = AdapterAction(str(value["action"]))
        except ValueError as error:
            raise AdapterError(f"invalid adapter receipt action: {value['action']!r}") from error
        return cls(
            schema_version=str(value["schema_version"]),
            action=action,
            run_id=run_id,
            status=status,
            artifacts=tuple(dict(item) for item in raw_artifacts),
            metadata=dict(metadata),
        )

    def validate_artifacts(self, root: str | Path) -> None:
        configured_root = Path(root)
        try:
            base = reject_symlink_components(
                configured_root, label="adapter output root"
            ).resolve(strict=True)
        except (OSError, RuntimeError) as error:
            raise AdapterError(f"adapter output root is invalid: {configured_root}") from error
        if not base.is_dir():
            raise AdapterError(f"adapter output root is not a directory: {base}")
        for artifact in self.artifacts:
            relative = Path(str(artifact.get("path", "")))
            if relative.is_absolute() or ".." in relative.parts:
                raise AdapterError(f"artifact path escapes output root: {relative}")
            configured_path = base / relative
            if configured_path.is_symlink():
                raise AdapterError(f"adapter artifact may not be a symlink: {relative}")
            path = configured_path.resolve()
            try:
                path.relative_to(base)
            except ValueError as error:
                raise AdapterError(
                    f"adapter artifact resolves outside output root: {relative}"
                ) from error
            if not path.is_file():
                raise AdapterError(f"adapter artifact is missing: {relative}")
            expected = str(artifact.get("sha256", ""))
            if not _SHA256.fullmatch(expected):
                raise AdapterError(f"adapter artifact has invalid SHA-256: {relative}")
            expected_size = artifact.get("size_bytes")
            if (
                isinstance(expected_size, bool)
                or not isinstance(expected_size, int)
                or expected_size < 0
            ):
                raise AdapterError(f"adapter artifact has invalid size: {relative}")
            if path.stat().st_size != expected_size:
                raise AdapterError(f"adapter artifact size mismatch: {relative}")
            observed = sha256_file(path)
            if expected != observed:
                raise AdapterError(f"adapter artifact checksum mismatch: {relative}")


class Adapter(ABC):
    """Five-action interface implemented inside each isolated environment."""

    @abstractmethod
    def invoke(
        self,
        action: AdapterAction,
        request_path: str | Path,
        output_dir: str | Path,
    ) -> AdapterReceipt:
        raise NotImplementedError

    def probe(self, request_path: str | Path, output_dir: str | Path) -> AdapterReceipt:
        return self.invoke(AdapterAction.PROBE, request_path, output_dir)

    def prepare(self, request_path: str | Path, output_dir: str | Path) -> AdapterReceipt:
        return self.invoke(AdapterAction.PREPARE, request_path, output_dir)

    def fit(self, request_path: str | Path, output_dir: str | Path) -> AdapterReceipt:
        return self.invoke(AdapterAction.FIT, request_path, output_dir)

    def predict(self, request_path: str | Path, output_dir: str | Path) -> AdapterReceipt:
        return self.invoke(AdapterAction.PREDICT, request_path, output_dir)

    def export(self, request_path: str | Path, output_dir: str | Path) -> AdapterReceipt:
        return self.invoke(AdapterAction.EXPORT, request_path, output_dir)


class SubprocessAdapter(Adapter):
    """Call a pinned environment without importing its dependencies."""

    def __init__(self, command_prefix: Sequence[str], *, timeout_seconds: int = 3600):
        if not command_prefix:
            raise AdapterError("adapter command prefix is empty")
        if any(
            not isinstance(item, str) or not item or "\x00" in item
            for item in command_prefix
        ):
            raise AdapterError("adapter command prefix must contain non-empty NUL-free strings")
        if (
            isinstance(timeout_seconds, bool)
            or not isinstance(timeout_seconds, int)
            or timeout_seconds < 1
        ):
            raise AdapterError("adapter timeout must be a positive integer")
        self._command_prefix = tuple(command_prefix)
        self._timeout_seconds = timeout_seconds

    def invoke(
        self,
        action: AdapterAction,
        request_path: str | Path,
        output_dir: str | Path,
    ) -> AdapterReceipt:
        configured_request = Path(request_path)
        configured_output = Path(output_dir)
        try:
            request = reject_symlink_components(
                configured_request, label="adapter request"
            ).resolve(strict=True)
            output = reject_symlink_components(
                configured_output, label="adapter output"
            )
        except (OSError, RuntimeError) as error:
            raise AdapterError("adapter request or output path is invalid") from error
        if not request.is_file():
            raise AdapterError(f"adapter request does not exist: {request}")
        request_snapshot = _stable_file_snapshot(request, "adapter request")
        try:
            request_payload = json.loads(request.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise AdapterError(f"invalid adapter request: {error}") from error
        if not isinstance(request_payload, Mapping):
            raise AdapterError("adapter request must be a JSON object")
        if request_payload.get("schema_version") != "masld-bench-adapter-request-v1":
            raise AdapterError("unsupported adapter request schema")
        expected_run_id = str(request_payload.get("run_id", ""))
        if not _SHA256.fullmatch(expected_run_id):
            raise AdapterError("adapter request run_id must be a full lowercase SHA-256")
        requested_action = request_payload.get("action")
        if str(requested_action) != action.value:
            raise AdapterError("adapter request action does not match invoked action")
        _validate_dataset_access_request(request_payload, action)
        _validate_prior_action_outputs(request, request_payload)
        if output.exists():
            raise AdapterError(f"adapter output must be a new directory: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        command = [
            *self._command_prefix,
            "--action",
            action.value,
            "--request",
            request.as_posix(),
            "--output",
            output.as_posix(),
        ]
        try:
            process = subprocess.run(
                command,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=self._timeout_seconds,
            )
        except subprocess.TimeoutExpired as error:
            if not output.exists():
                output.mkdir(parents=True, exist_ok=False)
            reject_symlink_components(output, label="adapter output")
            stdout = error.stdout or ""
            stderr = error.stderr or ""
            if isinstance(stdout, bytes):
                stdout = stdout.decode("utf-8", errors="replace")
            if isinstance(stderr, bytes):
                stderr = stderr.decode("utf-8", errors="replace")
            write_text_exclusive(output / "adapter_stdout.log", stdout)
            write_text_exclusive(output / "adapter_stderr.log", stderr)
            if _stable_file_snapshot(request, "adapter request") != request_snapshot:
                raise AdapterError("adapter request changed during subprocess execution")
            _validate_prior_action_outputs(request, request_payload)
            raise
        if _stable_file_snapshot(request, "adapter request") != request_snapshot:
            raise AdapterError("adapter request changed during subprocess execution")
        _validate_prior_action_outputs(request, request_payload)
        if not output.exists():
            output.mkdir(parents=True, exist_ok=False)
        reject_symlink_components(output, label="adapter output")
        if output.is_symlink() or not output.is_dir():
            raise AdapterError("adapter output is not a new non-symlink directory")
        write_text_exclusive(output / "adapter_stdout.log", process.stdout)
        write_text_exclusive(output / "adapter_stderr.log", process.stderr)
        if process.returncode != 0:
            raise AdapterError(
                f"adapter failed with exit {process.returncode}: {process.stderr[-4000:]}"
            )
        receipt_path = output / "adapter_receipt.json"
        try:
            receipt = AdapterReceipt.from_mapping(
                json.loads(receipt_path.read_text(encoding="utf-8"))
            )
        except (OSError, json.JSONDecodeError) as error:
            raise AdapterError(f"invalid adapter receipt: {error}") from error
        if receipt.action is not action:
            raise AdapterError("adapter receipt action does not match the request")
        if receipt.run_id != expected_run_id:
            raise AdapterError("adapter receipt run_id does not match the request")
        receipt.validate_artifacts(output)
        _validate_output_inventory(receipt, output)
        freeze_tree(
            output,
            {
                "artifact_class": "adapter_output",
                "action": action.value,
                "run_id": expected_run_id,
            },
        )
        verify_frozen_tree(output)
        return receipt
