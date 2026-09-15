"""Fail-closed assembly and verification of open benchmark releases."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import ipaddress
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tempfile
from typing import Any, Callable, Mapping, Sequence
from urllib.parse import urlsplit

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
from .hashing import is_sha256
from .planner import PlanningError, load_frozen_plan
from .selection import SelectionError, verify_selection_lock
from .tournament import (
    TournamentError,
    verify_selection_candidate_ledger,
    verify_terminal_evaluation_authorization,
)


class ReleaseError(RuntimeError):
    """Raised when publication output files violate the open-release rules."""


SPEC_SCHEMA_VERSION = "masld-bench-open-release-spec-v1"
RELEASE_SCHEMA_VERSION = "masld-bench-open-release-v1"
AUTHORIZATION_SCHEMA_VERSION = "masld-bench-terminal-authorization-v1"
TIMESTAMP_SCHEMA_VERSION = "masld-bench-external-timestamp-v1"
RELEASE_SCOPE = "research_checkpoint_no_api_portal_or_clinical_use"

REQUIRED_RELEASE_CLASSES = frozenset(
    {
        "source",
        "environment",
        "model_card",
        "frozen_predictions",
        "benchmark_table",
        "fetch_instructions",
    }
)
_IDENTIFIER = re.compile(r"^[a-z0-9][a-z0-9._-]{0,127}$")
_DESTINATION_PART = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
_UTC_TIMESTAMP = re.compile(
    r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d{1,6})?Z$"
)
_DANGEROUS_SUFFIXES = frozenset(
    {
        ".bcf",
        ".ckpt",
        ".cram",
        ".fastq",
        ".fq",
        ".h5",
        ".h5ad",
        ".loom",
        ".mtx",
        ".npy",
        ".npz",
        ".parquet",
        ".pickle",
        ".pkl",
        ".pt",
        ".pth",
        ".rds",
        ".safetensors",
        ".vcf",
    }
)
_ALLOWED_SUFFIXES = {
    "source": frozenset(
        {".c", ".cpp", ".json", ".md", ".py", ".r", ".sh", ".toml", ".txt", ".yaml", ".yml"}
    ),
    "environment": frozenset({".json", ".lock", ".toml", ".txt", ".yaml", ".yml"}),
    "model_card": frozenset({".json", ".md", ".txt"}),
    "frozen_predictions": frozenset({".csv", ".jsonl", ".tsv"}),
    "benchmark_table": frozenset({".csv", ".tsv"}),
    "fetch_instructions": frozenset({".json", ".md", ".txt"}),
}
_SOURCE_SPECIAL_NAMES = frozenset({"Dockerfile", "Makefile", "Snakefile"})
_PROHIBITED_DATA_TOKENS = frozenset(
    {
        "checkpoint",
        "checkpoints",
        "donor",
        "donors",
        "genotype",
        "genotypes",
        "label",
        "labels",
        "outcome",
        "outcomes",
        "participant",
        "participants",
        "processed",
        "raw",
        "state_dict",
        "weights",
    }
)
_PREDICTION_FIXED_FIELDS = frozenset(
    {"unit_hash", "run_id", "task_id", "model_id", "seed", "fold", "split_id", "class_id"}
)
_PREDICTION_PREFIXES = ("prediction", "probability", "score", "uncertainty", "logit")
_BENCHMARK_FORBIDDEN_FIELDS = frozenset(
    {
        "donor_id",
        "genotype",
        "individual_id",
        "label",
        "outcome",
        "participant_id",
        "patient_id",
        "sample_id",
        "subject_id",
        "unit_hash",
    }
)

TimestampProviderVerifier = Callable[..., datetime]
_TIMESTAMP_PROVIDER_VERIFIERS: dict[str, TimestampProviderVerifier] = {}


def _strict_mapping(value: object, fields: frozenset[str], label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ReleaseError(f"{label} must be an object")
    if any(not isinstance(key, str) for key in value):
        raise ReleaseError(f"{label} keys must be strings")
    missing = sorted(fields.difference(value))
    unknown = sorted(set(value).difference(fields))
    if missing:
        raise ReleaseError(f"{label} missing fields: {', '.join(missing)}")
    if unknown:
        raise ReleaseError(f"{label} has unknown fields: {', '.join(unknown)}")
    return dict(value)


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ReleaseError(f"cannot read {label}: {error}") from error
    if not isinstance(value, Mapping):
        raise ReleaseError(f"{label} must contain a JSON object")
    return dict(value)


def _string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value.strip() or "\x00" in value:
        raise ReleaseError(f"{label} must be a non-empty string")
    return value


def _sha256(value: object, label: str) -> str:
    if not is_sha256(value):
        raise ReleaseError(f"{label} must be a lowercase SHA-256")
    return str(value)


def _size(value: object, label: str, *, allow_zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ReleaseError(f"{label} must be an integer")
    if value < (0 if allow_zero else 1):
        raise ReleaseError(f"{label} must be {'non-negative' if allow_zero else 'positive'}")
    return value


def _absolute_dir(value: object, label: str) -> Path:
    text = _string(value, label)
    configured = Path(text)
    if not configured.is_absolute() or ".." in configured.parts:
        raise ReleaseError(f"{label} must be an absolute normalized directory")
    try:
        configured = reject_symlink_components(configured, label=label)
    except ArtifactError as error:
        raise ReleaseError(str(error)) from error
    if not configured.is_dir():
        raise ReleaseError(f"{label} must be a non-symlink directory")
    return configured.resolve()


def _absolute_file(value: object, label: str) -> Path:
    text = _string(value, label)
    configured = Path(text)
    if not configured.is_absolute() or ".." in configured.parts:
        raise ReleaseError(f"{label} must be an absolute normalized file")
    try:
        configured = reject_symlink_components(configured, label=label)
    except ArtifactError as error:
        raise ReleaseError(str(error)) from error
    if not configured.is_file():
        raise ReleaseError(f"{label} must be a non-symlink regular file")
    return configured.resolve()


def _absolute_protected_path(value: object, label: str) -> Path:
    """Validate a protected source root that may be a file or directory."""

    text = _string(value, label)
    configured = Path(text)
    if not configured.is_absolute() or ".." in configured.parts:
        raise ReleaseError(f"{label} must be an absolute normalized path")
    try:
        configured = reject_symlink_components(configured, label=label)
    except ArtifactError as error:
        raise ReleaseError(str(error)) from error
    if not configured.is_dir() and not configured.is_file():
        raise ReleaseError(
            f"{label} must be a non-symlink directory or regular file"
        )
    return configured.resolve()


def _frozen_binding(value: object, label: str) -> tuple[dict[str, Any], Path, dict[str, Any]]:
    binding = _strict_mapping(
        value,
        frozenset({"path", "manifest_sha256"}),
        label,
    )
    root = _absolute_dir(binding["path"], f"{label}.path")
    expected = _sha256(binding["manifest_sha256"], f"{label}.manifest_sha256")
    try:
        manifest = verify_frozen_tree(root)
        observed = sha256_file(root / "ARTIFACTS.json")
    except (ArtifactError, OSError, ValueError, KeyError, TypeError) as error:
        raise ReleaseError(f"invalid frozen {label}: {error}") from error
    if observed != expected:
        raise ReleaseError(f"{label} manifest SHA-256 mismatch")
    return binding, root, manifest


def _manifest_records(manifest: Mapping[str, Any], label: str) -> dict[str, dict[str, Any]]:
    records = manifest.get("artifacts")
    if not isinstance(records, list):
        raise ReleaseError(f"{label} frozen manifest has no artifact inventory")
    result: dict[str, dict[str, Any]] = {}
    for raw in records:
        if not isinstance(raw, Mapping):
            raise ReleaseError(f"{label} frozen manifest entry is not an object")
        path = str(raw.get("path", ""))
        if path in result:
            raise ReleaseError(f"{label} frozen manifest repeats {path}")
        result[path] = dict(raw)
    return result


def _document_in_tree(
    *,
    root: Path,
    manifest: Mapping[str, Any],
    filename: str,
    expected_sha256: object,
    label: str,
) -> dict[str, Any]:
    records = _manifest_records(manifest, label)
    record = records.get(filename)
    if record is None:
        raise ReleaseError(f"{label} manifest omits {filename}")
    expected = _sha256(expected_sha256, f"{label}.document_sha256")
    document = root / filename
    if sha256_file(document) != expected or record.get("sha256") != expected:
        raise ReleaseError(f"{label} document SHA-256 mismatch")
    if document.stat().st_size != record.get("size_bytes"):
        raise ReleaseError(f"{label} document size mismatch")
    return _load_json(document, label)


def _canonical_string_array(value: object, label: str, *, allow_empty: bool = False) -> list[str]:
    if not isinstance(value, list):
        raise ReleaseError(f"{label} must be an array")
    result = [_string(item, f"{label}[{index}]") for index, item in enumerate(value)]
    if not allow_empty and not result:
        raise ReleaseError(f"{label} must not be empty")
    if result != sorted(result) or len(result) != len(set(result)):
        raise ReleaseError(f"{label} must be unique and canonically sorted")
    return result


def _prediction_bindings(
    value: object, label: str, *, allow_empty: bool = False
) -> list[dict[str, str]]:
    if not isinstance(value, list) or (not value and not allow_empty):
        raise ReleaseError(
            f"{label} must be {'an array' if allow_empty else 'a non-empty array'}"
        )
    result: list[dict[str, str]] = []
    for index, item in enumerate(value):
        raw = _strict_mapping(
            item,
            frozenset({"run_id", "prediction_commit_sha256"}),
            f"{label}[{index}]",
        )
        result.append(
            {
                "run_id": _sha256(raw["run_id"], f"{label}[{index}].run_id"),
                "prediction_commit_sha256": _sha256(
                    raw["prediction_commit_sha256"],
                    f"{label}[{index}].prediction_commit_sha256",
                ),
            }
        )
    run_ids = [item["run_id"] for item in result]
    if run_ids != sorted(run_ids) or len(run_ids) != len(set(run_ids)):
        raise ReleaseError(f"{label} must be unique and sorted by run_id")
    return result


def _optional_sha256(value: object, label: str) -> str | None:
    return None if value is None else _sha256(value, label)


_DECISION_FIELDS = frozenset(
    {
        "decision_sha256",
        "task_id",
        "selected_model_id",
        "model_disposition",
        "task_disposition",
        "sealed",
        "open_champion",
        "gate_passed",
        "promoted",
        "terminal",
        "prediction_bindings",
        "power_decision_sha256",
        "outcome_bundle_sha256",
        "consumption_sha256",
    }
)


def _terminal_decisions(value: object) -> list[dict[str, Any]]:
    if not isinstance(value, list) or not value:
        raise ReleaseError("terminal_authorization.decisions must be a non-empty array")
    decisions: list[dict[str, Any]] = []
    for index, item in enumerate(value):
        raw = _strict_mapping(item, _DECISION_FIELDS, f"terminal decision {index}")
        decision = {
            "task_id": _string(raw["task_id"], f"terminal decision {index}.task_id"),
            "selected_model_id": _string(
                raw["selected_model_id"], f"terminal decision {index}.selected_model_id"
            ),
            "model_disposition": _string(
                raw["model_disposition"], f"terminal decision {index}.model_disposition"
            ),
            "task_disposition": _string(
                raw["task_disposition"], f"terminal decision {index}.task_disposition"
            ),
            "sealed": raw["sealed"],
            "open_champion": raw["open_champion"],
            "gate_passed": raw["gate_passed"],
            "promoted": raw["promoted"],
            "terminal": raw["terminal"],
            "prediction_bindings": _prediction_bindings(
                raw["prediction_bindings"],
                f"terminal decision {index}.prediction_bindings",
                allow_empty=True,
            ),
            "power_decision_sha256": _optional_sha256(
                raw["power_decision_sha256"],
                f"terminal decision {index}.power_decision_sha256",
            ),
            "outcome_bundle_sha256": _optional_sha256(
                raw["outcome_bundle_sha256"],
                f"terminal decision {index}.outcome_bundle_sha256",
            ),
            "consumption_sha256": _optional_sha256(
                raw["consumption_sha256"],
                f"terminal decision {index}.consumption_sha256",
            ),
        }
        for field in (
            "sealed",
            "open_champion",
            "gate_passed",
            "promoted",
            "terminal",
        ):
            if not isinstance(decision[field], bool):
                raise ReleaseError(f"terminal decision {index}.{field} must be boolean")
        claimed = _sha256(raw["decision_sha256"], f"terminal decision {index}.decision_sha256")
        if canonical_hash(decision) != claimed:
            raise ReleaseError(f"terminal decision {index} identity hash mismatch")
        decisions.append({"decision_sha256": claimed, **decision})
    task_ids = [item["task_id"] for item in decisions]
    if task_ids != sorted(task_ids) or len(task_ids) != len(set(task_ids)):
        raise ReleaseError("terminal decisions must be unique and sorted by task_id")
    return decisions


_AUTH_FIELDS = frozenset(
    {
        "schema_version",
        "authorization_id",
        "candidate_manifest_sha256",
        "plan_sha256",
        "selection_lock_id",
        "selection_lock_manifest_sha256",
        "selection_lock_sha256",
        "terminal_policy",
        "decisions",
    }
)


def _verify_authorization(
    value: object,
    *,
    candidate_manifest_sha256: str,
    plan_sha256: str,
    lock: Any,
    selection_lock_dir: Path,
    selection_manifest_sha256: str,
    selection_lock_sha256: str,
) -> tuple[dict[str, Any], dict[str, str], tuple[Path, ...]]:
    """Recursively verify the terminal checks and every authority they bind."""

    binding = _strict_mapping(
        value,
        frozenset({"path", "manifest_sha256", "document_sha256"}),
        "terminal_authorization",
    )
    _, root, manifest = _frozen_binding(
        {"path": binding["path"], "manifest_sha256": binding["manifest_sha256"]},
        "terminal_authorization",
    )
    document = _document_in_tree(
        root=root,
        manifest=manifest,
        filename="terminal_evaluation_authorization.json",
        expected_sha256=binding["document_sha256"],
        label="terminal_authorization",
    )
    try:
        verified = verify_terminal_evaluation_authorization(
            root, reverify_sources=True
        )
    except (TournamentError, OSError, ValueError) as error:
        raise ReleaseError(
            f"terminal authorization failed recursive verification: {error}"
        ) from error
    protected_raw = verified.pop("_protected_source_roots", None)
    if verified != document:
        raise ReleaseError("terminal authorization verifier returned a different document")
    selection_binding = verified.get("selection_lock_binding")
    if not isinstance(selection_binding, Mapping):
        raise ReleaseError("terminal authorization lacks a SelectionLock binding")
    expected_bindings = {
        "candidate_manifest_sha256": candidate_manifest_sha256,
        "plan_sha256": plan_sha256,
        "selection_lock_id": lock.lock_id,
        "selection_lock_manifest_sha256": selection_manifest_sha256,
        "terminal_policy": lock.terminal_policy,
    }
    if any(verified.get(key) != expected for key, expected in expected_bindings.items()):
        raise ReleaseError("terminal authorization does not bind the verified candidate and SelectionLock")
    expected_selection_binding = {
        "path": selection_lock_dir.as_posix(),
        "manifest_sha256": selection_manifest_sha256,
        "document_sha256": selection_lock_sha256,
        "selection_lock_id": lock.lock_id,
    }
    if dict(selection_binding) != expected_selection_binding:
        raise ReleaseError("terminal authorization SelectionLock binding changed")
    decisions = verified.get("decisions")
    if not isinstance(decisions, list) or not decisions:
        raise ReleaseError("terminal authorization decisions must be a non-empty array")
    if any(
        not isinstance(decision, Mapping) or decision.get("terminal") is not True
        for decision in decisions
    ):
        raise ReleaseError("terminal authorization contains a nonterminal task decision")
    claimed = _sha256(
        verified.get("authorization_id"), "authorization.authorization_id"
    )
    metadata = manifest.get("metadata")
    if metadata != {
        "artifact_class": "terminal_evaluation_authorization",
        "authorization_id": claimed,
    }:
        raise ReleaseError("terminal authorization frozen-tree metadata mismatch")
    if set(_manifest_records(manifest, "terminal_authorization")) != {
        "terminal_evaluation_authorization.json"
    }:
        raise ReleaseError("terminal authorization tree contains undeclared artifacts")
    if not isinstance(protected_raw, list) or any(
        not isinstance(item, str) for item in protected_raw
    ):
        raise ReleaseError("terminal authorization verifier omitted protected source roots")
    protected_roots = tuple(
        _absolute_protected_path(
            item, f"terminal_authorization.protected_source_roots[{index}]"
        )
        for index, item in enumerate(protected_raw)
    )
    return verified, {
        "path": root.as_posix(),
        "manifest_sha256": str(binding["manifest_sha256"]),
        "document_sha256": str(binding["document_sha256"]),
    }, protected_roots


def _plan_index(plan: Mapping[str, Any], field: str, key: str) -> dict[str, Mapping[str, Any]]:
    raw = plan.get(field)
    if not isinstance(raw, list):
        raise ReleaseError(f"candidate plan {field} must be an array")
    result: dict[str, Mapping[str, Any]] = {}
    for index, item in enumerate(raw):
        if not isinstance(item, Mapping):
            raise ReleaseError(f"candidate plan {field}[{index}] must be an object")
        identifier = _string(item.get(key), f"candidate plan {field}[{index}].{key}")
        if identifier in result:
            raise ReleaseError(f"candidate plan repeats {field} identity {identifier}")
        result[identifier] = item
    return result


def _validate_decisions(
    decisions: list[dict[str, Any]], *, lock: Any, plan: Mapping[str, Any]
) -> list[dict[str, Any]]:
    lock_tasks = {str(item["task_id"]): item for item in lock.task_decisions}
    if set(item["task_id"] for item in decisions) != set(lock_tasks):
        raise ReleaseError("terminal authorization must cover every and only locked task")
    model_dispositions = _plan_index(plan, "model_dispositions", "model_id")
    task_dispositions = _plan_index(plan, "task_dispositions", "task_id")
    runs = _plan_index(plan, "runs", "run_id")
    champions: list[dict[str, Any]] = []
    for decision in decisions:
        task_id = decision["task_id"]
        locked = lock_tasks[task_id]
        model_id = str(locked["selected_model_id"])
        selected_run_ids = sorted(str(item["run_id"]) for item in locked["selected_runs"])
        baseline_model_id = str(locked["baseline_model_id"])
        baseline_run_ids = sorted(str(item["run_id"]) for item in locked["baseline_runs"])
        observed_run_ids = sorted(
            str(item["run_id"]) for item in decision["prediction_bindings"]
        )
        # The variant power decision is required by firewall._complete_prediction_set
        # to bind every secondary comparator run in addition to the selected and
        # baseline runs, so the expected set must include them or a correctly
        # formed release is rejected.  Note that `locked` came out of a
        # SelectionLock, so its nested objects are mappingproxy and its arrays
        # are tuples: read them through Mapping/Sequence, never dict/list.
        secondary_run_ids: list[str] = []
        secondary = locked.get("variant_secondary_evaluation")
        if isinstance(secondary, Mapping):
            comparators = secondary.get("comparators")
            if not isinstance(comparators, (str, bytes)) and isinstance(
                comparators, Sequence
            ):
                secondary_run_ids = [
                    str(item["run_id"])
                    for item in comparators
                    if isinstance(item, Mapping)
                ]
        expected_run_ids = sorted(
            (*selected_run_ids, *baseline_run_ids, *secondary_run_ids)
        )
        if decision["selected_model_id"] != model_id:
            raise ReleaseError(f"terminal decision does not bind selected model for {task_id}")
        if decision["sealed"] and observed_run_ids != expected_run_ids:
            raise ReleaseError(f"terminal decision does not bind the required runs for {task_id}")
        if decision.get("terminal") is not True:
            raise ReleaseError(f"terminal decision is not terminal for {task_id}")
        if decision["sealed"] and (
            decision["power_decision_sha256"] is None
            or decision["outcome_bundle_sha256"] is None
            or decision.get("consumption_sha256") is None
        ):
            raise ReleaseError(f"opened terminal decision lacks power/outcome authority for {task_id}")
        if not decision["sealed"] and (
            decision["outcome_bundle_sha256"] is not None
            or decision.get("consumption_sha256") is not None
        ):
            raise ReleaseError(f"unopened terminal decision binds sealed outcomes for {task_id}")
        # Secondary comparator runs are bound by the variant power decision and
        # carry their own model_id, so map each to the model the selection record names for
        # it rather than to the selected or baseline model.
        secondary_models: dict[str, str] = {}
        if isinstance(secondary, Mapping):
            comparators = secondary.get("comparators")
            if not isinstance(comparators, (str, bytes)) and isinstance(
                comparators, Sequence
            ):
                secondary_models = {
                    str(item["run_id"]): str(item["model_id"])
                    for item in comparators
                    if isinstance(item, Mapping)
                }
        expected_models = {
            **{run_id: model_id for run_id in selected_run_ids},
            **{run_id: baseline_model_id for run_id in baseline_run_ids},
            **secondary_models,
        }
        for run_id in observed_run_ids:
            run = runs.get(run_id)
            if (
                run is None
                or run.get("task_id") != task_id
                or run.get("model_id") != expected_models.get(run_id)
            ):
                raise ReleaseError(f"candidate run disposition mismatch for {task_id}/{run_id}")
        model = model_dispositions.get(model_id)
        task = task_dispositions.get(task_id)
        if model is None or task is None:
            raise ReleaseError(f"candidate dispositions are incomplete for {task_id}/{model_id}")
        thresholds = locked.get("thresholds")
        promotion = thresholds.get("promotion") if isinstance(thresholds, Mapping) else None
        if (
            not isinstance(promotion, Mapping)
            or decision["model_disposition"] != promotion.get("promotion_mode")
        ):
            raise ReleaseError(f"terminal promotion mode mismatch for {task_id}")
        _string(decision["task_disposition"], f"terminal disposition for {task_id}")
        if decision["open_champion"] is not bool(locked.get("open_champion", False)):
            raise ReleaseError(f"terminal open-champion flag mismatch for {task_id}")
        promotable = (
            decision["sealed"] is True
            and decision["gate_passed"] is True
            and decision["open_champion"] is True
            and model.get("disposition") == "eligible_open"
            and model.get("champion_eligible") is True
            and task.get("disposition") == "eligible"
            and bool(locked.get("sealed_dataset_ids"))
        )
        if decision["promoted"] and not promotable:
            raise ReleaseError(f"terminal decision illegally promotes {task_id}/{model_id}")
        if decision["promoted"] is not (decision["task_disposition"] == "promoted"):
            raise ReleaseError(f"terminal promotion flags disagree for {task_id}")
        if decision["promoted"]:
            champions.append(
                {
                    "task_id": task_id,
                    "model_id": model_id,
                    "decision_sha256": decision["decision_sha256"],
                    "run_ids": selected_run_ids,
                }
            )
    return champions


def _identity_sets(authorization: Mapping[str, Any]) -> dict[str, str]:
    decisions = authorization["decisions"]
    prediction_identity = [
        {
            "task_id": item["task_id"],
            "prediction_bindings": item["prediction_bindings"],
        }
        for item in decisions
    ]
    power_identity = [
        {"task_id": item["task_id"], "power_decision_sha256": item["power_decision_sha256"]}
        for item in decisions
    ]
    outcome_identity = [
        {"task_id": item["task_id"], "outcome_bundle_sha256": item["outcome_bundle_sha256"]}
        for item in decisions
    ]
    return {
        "prediction_set_sha256": canonical_hash(prediction_identity),
        "power_set_sha256": canonical_hash(power_identity),
        "outcome_set_sha256": canonical_hash(outcome_identity),
        "terminal_decisions_sha256": canonical_hash(decisions),
    }


def _external_locator(value: object) -> str:
    locator = _string(value, "timestamp_receipt.external_locator")
    parsed = urlsplit(locator)
    host = (parsed.hostname or "").lower()
    local_host = False
    try:
        local_host = not ipaddress.ip_address(host).is_global
    except ValueError:
        pass
    if (
        parsed.scheme != "https"
        or not host
        or "." not in host
        or parsed.username is not None
        or parsed.password is not None
        or host in {"localhost", "127.0.0.1", "::1"}
        or host.endswith((".local", ".internal", ".invalid", ".test"))
        or local_host
    ):
        raise ReleaseError("timestamp receipt requires a non-local HTTPS external_locator")
    return locator


_RECEIPT_FIELDS = frozenset(
    {
        "schema_version",
        "receipt_id",
        "phase",
        "provider",
        "provider_receipt_id",
        "external_locator",
        "timestamp_utc",
        "bound_identity",
        "bound_identity_sha256",
        "proof",
    }
)


def _verify_provider_timestamp_proof(
    *,
    provider: str,
    provider_receipt_id: str,
    external_locator: str,
    timestamp: datetime,
    bound_identity_sha256: str,
    proof_path: Path,
    proof_sha256: str,
) -> None:
    """Dispatch to an explicit cryptographic/provider API verifier.

    The production registry is intentionally empty until a provider protocol,
    trust root, and verifier are pinned.  Merely naming a provider, supplying
    local bytes, or using a public-looking HTTPS URL is never evidence that an
    identity was externally timestamped.
    """

    verifier = _TIMESTAMP_PROVIDER_VERIFIERS.get(provider.casefold())
    if verifier is None:
        raise ReleaseError(
            "no provider-specific external timestamp verifier is configured "
            f"for {provider!r}; release remains fail-closed"
        )
    try:
        attested_time = verifier(
            provider_receipt_id=provider_receipt_id,
            external_locator=external_locator,
            bound_identity_sha256=bound_identity_sha256,
            proof_path=proof_path,
            proof_sha256=proof_sha256,
        )
    except ReleaseError:
        raise
    except Exception as error:
        raise ReleaseError(
            f"external timestamp provider verification failed for {provider!r}: {error}"
        ) from error
    if not isinstance(attested_time, datetime) or attested_time.tzinfo is None:
        raise ReleaseError(
            "external timestamp provider verifier did not return an aware datetime"
        )
    verified_utc = attested_time.astimezone(timezone.utc)
    if verified_utc != timestamp:
        raise ReleaseError(
            "timestamp_utc does not equal the provider-verified attestation time"
        )


def _verify_receipt(
    value: object,
    *,
    phase: str,
    expected_identity: Mapping[str, str],
) -> tuple[dict[str, Any], dict[str, str], datetime]:
    binding = _strict_mapping(
        value,
        frozenset({"path", "manifest_sha256", "document_sha256"}),
        f"receipts.{phase}",
    )
    _, root, manifest = _frozen_binding(
        {"path": binding["path"], "manifest_sha256": binding["manifest_sha256"]},
        f"receipts.{phase}",
    )
    raw = _strict_mapping(
        _document_in_tree(
            root=root,
            manifest=manifest,
            filename="timestamp_receipt.json",
            expected_sha256=binding["document_sha256"],
            label=f"receipts.{phase}",
        ),
        _RECEIPT_FIELDS,
        f"receipts.{phase}",
    )
    if raw["schema_version"] != TIMESTAMP_SCHEMA_VERSION or raw["phase"] != phase:
        raise ReleaseError(f"timestamp receipt phase/schema mismatch for {phase}")
    provider = _string(raw["provider"], f"receipts.{phase}.provider")
    if provider.strip().lower() in {"local", "self", "project", "internal"}:
        raise ReleaseError("timestamp provider must be independent and non-local")
    timestamp_text = _string(raw["timestamp_utc"], f"receipts.{phase}.timestamp_utc")
    if not _UTC_TIMESTAMP.fullmatch(timestamp_text):
        raise ReleaseError("timestamp_utc must be an RFC-3339 UTC timestamp ending in Z")
    timestamp = datetime.fromisoformat(timestamp_text[:-1] + "+00:00").astimezone(timezone.utc)
    identity = _strict_mapping(
        raw["bound_identity"], frozenset(expected_identity), f"receipts.{phase}.bound_identity"
    )
    if identity != dict(expected_identity):
        raise ReleaseError(f"timestamp receipt {phase} binds the wrong identities")
    identity_sha = _sha256(
        raw["bound_identity_sha256"], f"receipts.{phase}.bound_identity_sha256"
    )
    if canonical_hash(identity) != identity_sha:
        raise ReleaseError(f"timestamp receipt {phase} identity digest mismatch")
    proof = _strict_mapping(
        raw["proof"], frozenset({"path", "sha256", "size_bytes"}), f"receipts.{phase}.proof"
    )
    proof_path = _safe_relative(proof["path"], f"receipts.{phase}.proof.path")
    if proof_path.as_posix() == "timestamp_receipt.json":
        raise ReleaseError("timestamp provider proof must be separate from its wrapper")
    proof_sha = _sha256(proof["sha256"], f"receipts.{phase}.proof.sha256")
    proof_size = _size(proof["size_bytes"], f"receipts.{phase}.proof.size_bytes")
    proof_file = root.joinpath(*proof_path.parts)
    records = _manifest_records(manifest, f"receipts.{phase}")
    record = records.get(proof_path.as_posix())
    if (
        record is None
        or not proof_file.is_file()
        or proof_file.is_symlink()
        or proof_file.stat().st_size != proof_size
        or sha256_file(proof_file) != proof_sha
        or record.get("sha256") != proof_sha
        or record.get("size_bytes") != proof_size
    ):
        raise ReleaseError(f"timestamp receipt {phase} provider proof mismatch")
    if set(records) != {"timestamp_receipt.json", proof_path.as_posix()}:
        raise ReleaseError(f"timestamp receipt {phase} tree contains undeclared artifacts")
    provider_receipt_id = _string(
        raw["provider_receipt_id"], f"receipts.{phase}.provider_receipt_id"
    )
    external_locator = _external_locator(raw["external_locator"])
    _verify_provider_timestamp_proof(
        provider=provider,
        provider_receipt_id=provider_receipt_id,
        external_locator=external_locator,
        timestamp=timestamp,
        bound_identity_sha256=identity_sha,
        proof_path=proof_file,
        proof_sha256=proof_sha,
    )
    normalized = {
        "schema_version": TIMESTAMP_SCHEMA_VERSION,
        "phase": phase,
        "provider": provider,
        "provider_receipt_id": provider_receipt_id,
        "external_locator": external_locator,
        "timestamp_utc": timestamp_text,
        "bound_identity": identity,
        "bound_identity_sha256": identity_sha,
        "proof": {"path": proof_path.as_posix(), "sha256": proof_sha, "size_bytes": proof_size},
    }
    receipt_id = _sha256(raw["receipt_id"], f"receipts.{phase}.receipt_id")
    if canonical_hash(normalized) != receipt_id:
        raise ReleaseError(f"timestamp receipt {phase} canonical receipt ID mismatch")
    metadata = manifest.get("metadata")
    if metadata != {"artifact_class": "external_timestamp_receipt", "receipt_id": receipt_id}:
        raise ReleaseError(f"timestamp receipt {phase} frozen-tree metadata mismatch")
    return (
        {"receipt_id": receipt_id, **normalized},
        {
            "manifest_sha256": str(binding["manifest_sha256"]),
            "document_sha256": str(binding["document_sha256"]),
        },
        timestamp,
    )


def _safe_relative(value: object, label: str) -> PurePosixPath:
    text = _string(value, label)
    if "\\" in text or text.startswith("/") or text.endswith("/"):
        raise ReleaseError(f"{label} must be a normalized relative POSIX path")
    path = PurePosixPath(text)
    if (
        not path.parts
        or any(part in {"", ".", ".."} for part in path.parts)
        or len(path.parts) > 8
        or len(text) > 240
        or any(not _DESTINATION_PART.fullmatch(part) for part in path.parts)
    ):
        raise ReleaseError(f"{label} is not a safe release destination")
    if path.as_posix() != text:
        raise ReleaseError(f"{label} must be canonically normalized")
    return path


def _validate_text(path: Path, label: str) -> None:
    try:
        with path.open("r", encoding="utf-8", errors="strict", newline="") as handle:
            while chunk := handle.read(1024 * 1024):
                if "\x00" in chunk:
                    raise ReleaseError(f"{label} contains binary NUL bytes")
    except UnicodeDecodeError as error:
        raise ReleaseError(f"{label} must be UTF-8 text, not a binary artifact") from error


def _normalized_field(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def _validate_predictions(path: Path) -> None:
    suffix = path.suffix.lower()
    if suffix == ".jsonl":
        row_count = 0
        with path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                try:
                    row = json.loads(line)
                except json.JSONDecodeError as error:
                    raise ReleaseError(f"invalid frozen prediction JSONL line {line_number}") from error
                if not isinstance(row, Mapping):
                    raise ReleaseError("frozen prediction JSONL rows must be objects")
                _validate_prediction_row(dict(row), f"line {line_number}")
                row_count += 1
        if not row_count:
            raise ReleaseError("frozen predictions must contain at least one row")
        return
    delimiter = "\t" if suffix == ".tsv" else ","
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ReleaseError("frozen predictions require a unique header")
        row_count = 0
        for row_count, row in enumerate(reader, start=1):
            if None in row:
                raise ReleaseError("frozen prediction row has more values than its header")
            _validate_prediction_row(dict(row), f"row {row_count}")
        if not row_count:
            raise ReleaseError("frozen predictions must contain at least one row")


def _validate_prediction_row(row: Mapping[str, Any], label: str) -> None:
    fields = {_normalized_field(str(key)) for key in row}
    unknown = sorted(
        field
        for field in fields
        if field not in _PREDICTION_FIXED_FIELDS
        and not field.startswith(_PREDICTION_PREFIXES)
    )
    if unknown:
        raise ReleaseError("frozen predictions contain forbidden/non-output fields: " + ", ".join(unknown))
    if "unit_hash" not in fields or not any(field.startswith(_PREDICTION_PREFIXES) for field in fields):
        raise ReleaseError("frozen predictions require unit_hash and prediction output fields")
    unit_key = next(key for key in row if _normalized_field(str(key)) == "unit_hash")
    if not is_sha256(row[unit_key]):
        raise ReleaseError(f"frozen prediction {label} unit_hash must be a lowercase SHA-256")


def _validate_benchmark(path: Path) -> None:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        try:
            header = next(reader)
        except StopIteration as error:
            raise ReleaseError("benchmark table must not be empty") from error
        fields = [_normalized_field(item) for item in header]
        if not fields or any(not field for field in fields) or len(fields) != len(set(fields)):
            raise ReleaseError("benchmark table requires a unique non-empty header")
        forbidden = sorted(set(fields).intersection(_BENCHMARK_FORBIDDEN_FIELDS))
        if forbidden:
            raise ReleaseError("benchmark table contains unit-level or sealed fields: " + ", ".join(forbidden))
        try:
            next(reader)
        except StopIteration as error:
            raise ReleaseError("benchmark table must contain aggregate result rows") from error


def _validate_open_artifact(path: Path, artifact_class: str, destination: PurePosixPath) -> None:
    suffix = destination.suffix.lower()
    if suffix in _DANGEROUS_SUFFIXES:
        raise ReleaseError(f"prohibited data/checkpoint/weight artifact: {destination}")
    allowed = _ALLOWED_SUFFIXES[artifact_class]
    if suffix not in allowed and not (
        artifact_class == "source" and destination.name in _SOURCE_SPECIAL_NAMES
    ):
        raise ReleaseError(f"unsupported {artifact_class} artifact type: {destination}")
    if artifact_class not in {"source", "fetch_instructions"}:
        tokens = {
            token
            for part in destination.parts
            for token in re.split(r"[^a-z0-9]+", part.lower())
            if token
        }
        prohibited = sorted(tokens.intersection(_PROHIBITED_DATA_TOKENS))
        if prohibited:
            raise ReleaseError("artifact destination signals prohibited content: " + ", ".join(prohibited))
    _validate_text(path, f"{artifact_class} artifact {destination}")
    if artifact_class == "frozen_predictions":
        _validate_predictions(path)
    elif artifact_class == "benchmark_table":
        _validate_benchmark(path)


_ARTIFACT_FIELDS = frozenset({"class", "path", "sha256", "size_bytes", "destination"})


def _stage_artifacts(
    artifacts: object,
    *,
    staging: Path,
    champions: list[dict[str, Any]],
    lock_tasks: Mapping[str, Mapping[str, Any]],
    protected_roots: tuple[Path, ...],
) -> list[dict[str, Any]]:
    if not isinstance(artifacts, list) or not artifacts:
        raise ReleaseError("release artifacts must be a non-empty array")
    inventory: list[dict[str, Any]] = []
    source_paths: set[Path] = set()
    target_paths: set[str] = set()
    document_text: dict[str, list[str]] = {"model_card": [], "fetch_instructions": []}
    for index, item in enumerate(artifacts):
        raw = _strict_mapping(item, _ARTIFACT_FIELDS, f"artifact {index}")
        artifact_class = _string(raw["class"], f"artifact {index}.class")
        if artifact_class not in REQUIRED_RELEASE_CLASSES:
            raise ReleaseError(f"artifact {index} has prohibited or unknown class {artifact_class!r}")
        source = _absolute_file(raw["path"], f"artifact {index}.path")
        for protected in protected_roots:
            try:
                source.relative_to(protected)
            except ValueError:
                continue
            raise ReleaseError(
                f"open artifact may not be sourced from protected frozen tree: {source}"
            )
        if source in source_paths:
            raise ReleaseError(f"release artifact source is repeated: {source}")
        source_paths.add(source)
        destination = _safe_relative(raw["destination"], f"artifact {index}.destination")
        release_path = PurePosixPath("artifacts") / artifact_class / destination
        if release_path.as_posix() in target_paths:
            raise ReleaseError(f"release destination is repeated: {release_path}")
        target_paths.add(release_path.as_posix())
        expected_sha = _sha256(raw["sha256"], f"artifact {index}.sha256")
        expected_size = _size(raw["size_bytes"], f"artifact {index}.size_bytes")
        if source.stat().st_size != expected_size or sha256_file(source) != expected_sha:
            raise ReleaseError(f"release artifact checksum/size mismatch: {source}")
        _validate_open_artifact(source, artifact_class, destination)
        target = staging.joinpath(*release_path.parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        if target.stat().st_size != expected_size or sha256_file(target) != expected_sha:
            raise ReleaseError(f"staged artifact changed while copying: {source}")
        os.chmod(target, 0o444)
        if artifact_class in document_text:
            document_text[artifact_class].append(target.read_text(encoding="utf-8"))
        inventory.append(
            {
                "class": artifact_class,
                "path": release_path.as_posix(),
                "sha256": expected_sha,
                "size_bytes": expected_size,
            }
        )
    inventory.sort(key=lambda item: item["path"])
    classes = {item["class"] for item in inventory}
    missing = sorted(REQUIRED_RELEASE_CLASSES.difference(classes))
    if missing:
        raise ReleaseError("release is missing artifact classes: " + ", ".join(missing))
    model_card = "\n".join(document_text["model_card"]).lower()
    fetch = "\n".join(document_text["fetch_instructions"]).lower()
    if not any(phrase in fetch for phrase in ("not redistributed", "not included in this release", "no external weights")):
        raise ReleaseError("fetch instructions must state that weights are not redistributed")
    for champion in champions:
        model_id = champion["model_id"].lower()
        checkpoint_sha = str(lock_tasks[champion["task_id"]]["checkpoint_sha256"]).lower()
        if model_id not in model_card:
            raise ReleaseError(f"model card does not identify champion {champion['model_id']}")
        if model_id not in fetch or checkpoint_sha not in fetch:
            raise ReleaseError(f"fetch instructions do not bind champion {champion['model_id']} checkpoint")
    return inventory


def _readonly_tree(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        os.chmod(path, 0o555 if path.is_dir() else 0o444)
    os.chmod(root, 0o555)


def stage_release(*, spec_path: str | Path, output_root: str | Path) -> Path:
    """Validate every upstream authority, then atomically stage a read-only release."""

    try:
        configured_spec = reject_symlink_components(
            Path(spec_path), label="release specification path"
        )
        resolved_spec = configured_spec.resolve(strict=True)
    except (ArtifactError, OSError) as error:
        raise ReleaseError(f"invalid release specification path: {error}") from error
    if not resolved_spec.is_file():
        raise ReleaseError("release specification path must be a regular file")
    spec = _strict_mapping(
        _load_json(resolved_spec, "release specification"),
        frozenset(
            {
                "schema_version",
                "release_id",
                "candidate",
                "selection_lock",
                "terminal_authorization",
                "receipts",
                "artifacts",
            }
        ),
        "release specification",
    )
    if spec["schema_version"] != SPEC_SCHEMA_VERSION:
        raise ReleaseError("unsupported release specification schema_version")
    release_id = _string(spec["release_id"], "release_id")
    if not _IDENTIFIER.fullmatch(release_id):
        raise ReleaseError(f"release_id must match {_IDENTIFIER.pattern!r}")

    _, candidate_dir, candidate_manifest = _frozen_binding(spec["candidate"], "candidate")
    metadata = candidate_manifest.get("metadata")
    artifact_class = metadata.get("artifact_class") if isinstance(metadata, Mapping) else None
    if artifact_class == "candidate_plan":
        try:
            plan = load_frozen_plan(candidate_dir)
        except (PlanningError, OSError, ValueError, KeyError, TypeError) as error:
            raise ReleaseError(f"invalid frozen candidate plan: {error}") from error
    elif artifact_class == "selection_candidate_ledger":
        try:
            ledger = verify_selection_candidate_ledger(candidate_dir)
        except (TournamentError, OSError, ValueError) as error:
            raise ReleaseError(f"invalid selection candidate ledger: {error}") from error
        plan = {
            "plan_sha256": ledger["ledger_id"],
            "model_dispositions": ledger["model_dispositions"],
            "task_dispositions": ledger["task_dispositions"],
            "runs": ledger["runs"],
        }
    else:
        raise ReleaseError(
            "candidate tree must be a candidate_plan or selection_candidate_ledger"
        )
    candidate_manifest_sha = sha256_file(candidate_dir / "ARTIFACTS.json")

    selection_spec = _strict_mapping(
        spec["selection_lock"],
        frozenset({"path", "manifest_sha256", "lock_sha256"}),
        "selection_lock",
    )
    _, selection_dir, _ = _frozen_binding(
        {"path": selection_spec["path"], "manifest_sha256": selection_spec["manifest_sha256"]},
        "selection_lock",
    )
    expected_lock_sha = _sha256(selection_spec["lock_sha256"], "selection_lock.lock_sha256")
    if sha256_file(selection_dir / "selection_lock.json") != expected_lock_sha:
        raise ReleaseError("SelectionLock document SHA-256 mismatch")
    try:
        lock = verify_selection_lock(selection_dir)
    except SelectionError as error:
        raise ReleaseError(f"invalid frozen SelectionLock: {error}") from error
    selection_manifest_sha = sha256_file(selection_dir / "ARTIFACTS.json")
    if (
        lock.candidate_manifest_sha256 != candidate_manifest_sha
        or lock.plan_sha256 != plan.get("plan_sha256")
    ):
        raise ReleaseError("SelectionLock does not bind the supplied frozen candidate")

    authorization, authorization_binding, authorization_source_roots = _verify_authorization(
        spec["terminal_authorization"],
        candidate_manifest_sha256=candidate_manifest_sha,
        plan_sha256=str(plan["plan_sha256"]),
        lock=lock,
        selection_lock_dir=selection_dir,
        selection_manifest_sha256=selection_manifest_sha,
        selection_lock_sha256=expected_lock_sha,
    )
    champions = _validate_decisions(authorization["decisions"], lock=lock, plan=plan)
    identity_sets = _identity_sets(authorization)

    receipts = _strict_mapping(
        spec["receipts"], frozenset({"pre_unblind", "terminal"}), "receipts"
    )
    pre_identity = {
        "selection_lock_id": lock.lock_id,
        "selection_lock_manifest_sha256": selection_manifest_sha,
        "prediction_set_sha256": identity_sets["prediction_set_sha256"],
        "power_set_sha256": identity_sets["power_set_sha256"],
    }
    pre_receipt, pre_binding, pre_time = _verify_receipt(
        receipts["pre_unblind"], phase="pre_unblind", expected_identity=pre_identity
    )
    terminal_identity = {
        "authorization_id": authorization["authorization_id"],
        "authorization_manifest_sha256": authorization_binding["manifest_sha256"],
        "outcome_set_sha256": identity_sets["outcome_set_sha256"],
        "terminal_decisions_sha256": identity_sets["terminal_decisions_sha256"],
    }
    terminal_receipt, terminal_binding, terminal_time = _verify_receipt(
        receipts["terminal"], phase="terminal", expected_identity=terminal_identity
    )
    if not pre_time < terminal_time:
        raise ReleaseError("pre-unblind receipt must strictly precede terminal receipt")

    try:
        root = reject_symlink_components(Path(output_root), label="release output_root")
    except ArtifactError as error:
        raise ReleaseError(str(error)) from error
    if root.exists() and not root.is_dir():
        raise ReleaseError("output_root must be a non-symlink directory")
    root.mkdir(parents=True, exist_ok=True)
    try:
        root = reject_symlink_components(root, label="release output_root")
    except ArtifactError as error:
        raise ReleaseError(str(error)) from error
    root = root.resolve()
    protected_roots = tuple(
        dict.fromkeys(
            (
                candidate_dir,
                selection_dir,
                _absolute_dir(
                    authorization_binding["path"], "terminal_authorization.path"
                ),
                _absolute_dir(
                    receipts["pre_unblind"]["path"], "receipts.pre_unblind.path"
                ),
                _absolute_dir(
                    receipts["terminal"]["path"], "receipts.terminal.path"
                ),
                *authorization_source_roots,
            )
        )
    )
    for protected in protected_roots:
        try:
            root.relative_to(protected)
        except ValueError:
            continue
        raise ReleaseError("output_root may not be inside a protected frozen tree")
    target = root / release_id
    if target.exists():
        raise ReleaseError(f"release already exists and cannot be overwritten: {target}")
    staging = Path(tempfile.mkdtemp(prefix=f".{release_id}.", dir=root))
    try:
        lock_tasks = {str(item["task_id"]): item for item in lock.task_decisions}
        inventory = _stage_artifacts(
            spec["artifacts"],
            staging=staging,
            champions=champions,
            lock_tasks=lock_tasks,
            protected_roots=protected_roots,
        )
        frozen_sources = (
            (candidate_dir, candidate_manifest_sha),
            (selection_dir, selection_manifest_sha),
            (
                _absolute_dir(
                    authorization_binding["path"], "terminal_authorization.path"
                ),
                authorization_binding["manifest_sha256"],
            ),
            (
                _absolute_dir(
                    receipts["pre_unblind"]["path"], "receipts.pre_unblind.path"
                ),
                pre_binding["manifest_sha256"],
            ),
            (
                _absolute_dir(
                    receipts["terminal"]["path"], "receipts.terminal.path"
                ),
                terminal_binding["manifest_sha256"],
            ),
        )
        for frozen_root, expected_manifest_sha in frozen_sources:
            try:
                verify_frozen_tree(frozen_root)
            except (ArtifactError, OSError, ValueError, KeyError, TypeError) as error:
                raise ReleaseError(
                    f"upstream frozen tree changed during release staging: {frozen_root}: {error}"
                ) from error
            if sha256_file(frozen_root / "ARTIFACTS.json") != expected_manifest_sha:
                raise ReleaseError(
                    f"upstream frozen-tree identity changed during release staging: {frozen_root}"
                )
        reverified_authorization, reverified_binding, reverified_roots = _verify_authorization(
            spec["terminal_authorization"],
            candidate_manifest_sha256=candidate_manifest_sha,
            plan_sha256=str(plan["plan_sha256"]),
            lock=lock,
            selection_lock_dir=selection_dir,
            selection_manifest_sha256=selection_manifest_sha,
            selection_lock_sha256=expected_lock_sha,
        )
        if (
            reverified_authorization != authorization
            or reverified_binding != authorization_binding
            or reverified_roots != authorization_source_roots
        ):
            raise ReleaseError("terminal authorization sources changed during release staging")
        release_payload: dict[str, Any] = {
            "schema_version": RELEASE_SCHEMA_VERSION,
            "release_id": release_id,
            "scope": RELEASE_SCOPE,
            "publication_disposition": (
                "open_champions" if champions else "negative_result_no_champion"
            ),
            "candidate": {
                "manifest_sha256": candidate_manifest_sha,
                "plan_sha256": plan["plan_sha256"],
            },
            "selection_lock": {
                "lock_id": lock.lock_id,
                "manifest_sha256": selection_manifest_sha,
                "lock_sha256": expected_lock_sha,
            },
            "terminal_authorization": {
                "authorization_id": authorization["authorization_id"],
                "manifest_sha256": authorization_binding["manifest_sha256"],
                "document_sha256": authorization_binding["document_sha256"],
                **identity_sets,
            },
            "receipts": {
                "pre_unblind": {
                    "receipt_id": pre_receipt["receipt_id"],
                    **pre_binding,
                    "bound_identity_sha256": pre_receipt["bound_identity_sha256"],
                },
                "terminal": {
                    "receipt_id": terminal_receipt["receipt_id"],
                    **terminal_binding,
                    "bound_identity_sha256": terminal_receipt["bound_identity_sha256"],
                },
            },
            "champions": champions,
            "artifacts": inventory,
            "restricted_weight_policy": "fetch_instructions_only_no_checkpoints_or_weights",
        }
        release_identity_sha = canonical_hash(release_payload)
        manifest = {**release_payload, "release_identity_sha256": release_identity_sha}
        record = write_json_exclusive(staging / "release_manifest.json", manifest, mode=0o444)
        freeze_tree(
            staging,
            {
                "artifact_class": "open_research_release",
                "release_id": release_id,
                "release_identity_sha256": release_identity_sha,
                "release_manifest_sha256": record.sha256,
            },
        )
        _readonly_tree(staging)
        try:
            publish_directory_noreplace(staging, target)
        except ArtifactError as error:
            raise ReleaseError(f"could not publish immutable release: {error}") from error
    except Exception:
        if staging.exists():
            try:
                for path in staging.rglob("*"):
                    if path.exists():
                        os.chmod(path, 0o755 if path.is_dir() else 0o644)
                os.chmod(staging, 0o755)
            except OSError:
                pass
            shutil.rmtree(staging, ignore_errors=True)
        raise
    verify_release(target)
    return target


_RELEASE_FIELDS = frozenset(
    {
        "schema_version",
        "release_id",
        "scope",
        "publication_disposition",
        "candidate",
        "selection_lock",
        "terminal_authorization",
        "receipts",
        "champions",
        "artifacts",
        "restricted_weight_policy",
        "release_identity_sha256",
    }
)


def _verify_release_bindings(manifest: Mapping[str, Any]) -> None:
    candidate = _strict_mapping(
        manifest["candidate"],
        frozenset({"manifest_sha256", "plan_sha256"}),
        "release candidate binding",
    )
    _sha256(candidate["manifest_sha256"], "release candidate manifest_sha256")
    _sha256(candidate["plan_sha256"], "release candidate plan_sha256")
    selection = _strict_mapping(
        manifest["selection_lock"],
        frozenset({"lock_id", "manifest_sha256", "lock_sha256"}),
        "release SelectionLock binding",
    )
    for field in ("lock_id", "manifest_sha256", "lock_sha256"):
        _sha256(selection[field], f"release SelectionLock {field}")
    authorization = _strict_mapping(
        manifest["terminal_authorization"],
        frozenset(
            {
                "authorization_id",
                "manifest_sha256",
                "document_sha256",
                "prediction_set_sha256",
                "power_set_sha256",
                "outcome_set_sha256",
                "terminal_decisions_sha256",
            }
        ),
        "release terminal authorization binding",
    )
    for field, value in authorization.items():
        _sha256(value, f"release terminal authorization {field}")
    receipts = _strict_mapping(
        manifest["receipts"],
        frozenset({"pre_unblind", "terminal"}),
        "release receipt bindings",
    )
    for phase in ("pre_unblind", "terminal"):
        receipt = _strict_mapping(
            receipts[phase],
            frozenset(
                {
                    "receipt_id",
                    "manifest_sha256",
                    "document_sha256",
                    "bound_identity_sha256",
                }
            ),
            f"release {phase} receipt binding",
        )
        for field, value in receipt.items():
            _sha256(value, f"release {phase} receipt {field}")
    champions = manifest["champions"]
    if not isinstance(champions, list):
        raise ReleaseError("release champions must be an array")
    normalized_keys: list[tuple[str, str]] = []
    for index, item in enumerate(champions):
        champion = _strict_mapping(
            item,
            frozenset({"task_id", "model_id", "decision_sha256", "run_ids"}),
            f"release champion {index}",
        )
        task_id = _string(champion["task_id"], f"release champion {index}.task_id")
        model_id = _string(champion["model_id"], f"release champion {index}.model_id")
        _sha256(champion["decision_sha256"], f"release champion {index}.decision_sha256")
        run_ids = _canonical_string_array(
            champion["run_ids"], f"release champion {index}.run_ids"
        )
        for run_index, run_id in enumerate(run_ids):
            _sha256(run_id, f"release champion {index}.run_ids[{run_index}]")
        normalized_keys.append((task_id, model_id))
    if normalized_keys != sorted(normalized_keys) or len(normalized_keys) != len(
        set(normalized_keys)
    ):
        raise ReleaseError("release champions must be unique and canonically sorted")


def verify_release(path: str | Path) -> Mapping[str, Any]:
    """Verify the exact inventory, canonical identity, and read-only release tree."""

    try:
        root = reject_symlink_components(Path(path), label="release path")
    except ArtifactError as error:
        raise ReleaseError(str(error)) from error
    if not root.is_dir():
        raise ReleaseError("release must be a non-symlink directory")
    root = root.resolve()
    try:
        frozen_manifest = verify_frozen_tree(root)
    except (ArtifactError, OSError, ValueError, KeyError, TypeError) as error:
        raise ReleaseError(f"invalid frozen release: {error}") from error
    manifest = _strict_mapping(
        _load_json(root / "release_manifest.json", "release manifest"),
        _RELEASE_FIELDS,
        "release manifest",
    )
    if (
        manifest["schema_version"] != RELEASE_SCHEMA_VERSION
        or manifest["scope"] != RELEASE_SCOPE
        or manifest["restricted_weight_policy"]
        != "fetch_instructions_only_no_checkpoints_or_weights"
    ):
        raise ReleaseError("release manifest policy/schema mismatch")
    release_id = _string(manifest["release_id"], "release manifest release_id")
    if not _IDENTIFIER.fullmatch(release_id) or root.name != release_id:
        raise ReleaseError("release directory and strict release_id disagree")
    claimed_identity = _sha256(
        manifest["release_identity_sha256"], "release_identity_sha256"
    )
    identity = dict(manifest)
    identity.pop("release_identity_sha256")
    if canonical_hash(identity) != claimed_identity:
        raise ReleaseError("release manifest canonical identity mismatch")
    _verify_release_bindings(manifest)
    metadata = frozen_manifest.get("metadata")
    if metadata != {
        "artifact_class": "open_research_release",
        "release_id": release_id,
        "release_identity_sha256": claimed_identity,
        "release_manifest_sha256": sha256_file(root / "release_manifest.json"),
    }:
        raise ReleaseError("release frozen-tree metadata mismatch")
    artifacts = manifest["artifacts"]
    if not isinstance(artifacts, list) or not artifacts:
        raise ReleaseError("release artifact inventory must be non-empty")
    expected_paths = {"release_manifest.json"}
    observed_classes: set[str] = set()
    previous_path = ""
    for index, item in enumerate(artifacts):
        raw = _strict_mapping(
            item, frozenset({"class", "path", "sha256", "size_bytes"}), f"release artifact {index}"
        )
        artifact_class = _string(raw["class"], f"release artifact {index}.class")
        if artifact_class not in REQUIRED_RELEASE_CLASSES:
            raise ReleaseError(f"release contains prohibited artifact class {artifact_class!r}")
        relative = _safe_relative(raw["path"], f"release artifact {index}.path")
        if len(relative.parts) < 3 or relative.parts[:2] != ("artifacts", artifact_class):
            raise ReleaseError("release artifact path does not match its class")
        path_text = relative.as_posix()
        if path_text <= previous_path or path_text in expected_paths:
            raise ReleaseError("release artifacts must be unique and sorted by path")
        previous_path = path_text
        expected_paths.add(path_text)
        observed_classes.add(artifact_class)
        artifact_path = root.joinpath(*relative.parts)
        expected_sha = _sha256(raw["sha256"], f"release artifact {index}.sha256")
        expected_size = _size(raw["size_bytes"], f"release artifact {index}.size_bytes")
        if (
            artifact_path.is_symlink()
            or not artifact_path.is_file()
            or artifact_path.stat().st_size != expected_size
            or sha256_file(artifact_path) != expected_sha
        ):
            raise ReleaseError(f"release artifact checksum/size mismatch: {path_text}")
        destination = PurePosixPath(*relative.parts[2:])
        _validate_open_artifact(artifact_path, artifact_class, destination)
    missing = sorted(REQUIRED_RELEASE_CLASSES.difference(observed_classes))
    if missing:
        raise ReleaseError("release is missing artifact classes: " + ", ".join(missing))
    frozen_paths = set(_manifest_records(frozen_manifest, "release"))
    if frozen_paths != expected_paths:
        raise ReleaseError("release manifest and frozen-tree inventories disagree")
    disposition = manifest["publication_disposition"]
    champions = manifest["champions"]
    expected_disposition = "open_champions" if champions else "negative_result_no_champion"
    if disposition != expected_disposition:
        raise ReleaseError("publication disposition disagrees with derived champions")
    for member in (root, *root.rglob("*")):
        if member.stat().st_mode & 0o222:
            raise ReleaseError(f"release tree is not read-only: {member}")
    return manifest


__all__ = [
    "AUTHORIZATION_SCHEMA_VERSION",
    "RELEASE_SCHEMA_VERSION",
    "REQUIRED_RELEASE_CLASSES",
    "ReleaseError",
    "SPEC_SCHEMA_VERSION",
    "TIMESTAMP_SCHEMA_VERSION",
    "stage_release",
    "verify_release",
]
