#!/usr/bin/env python3
"""Compare two isolated REL01--04 candidate rebuilds before REL05.

This is deliberately a pre-REL05 gate: it compares the frozen REL02--04
product inventories and never reads or creates a REL05 seal.
"""

from __future__ import annotations

import argparse
import importlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
if str(PROGRAM_ROOT) not in sys.path:
    sys.path.insert(0, str(PROGRAM_ROOT))

from release_common import (  # noqa: E402
    CANDIDATE_ID,
    ReleaseContractError,
    atomic_write_tsv,
    canonical_json_bytes,
    parse_nonnegative_int,
    read_tsv_exact,
    require_sha256,
    resolve_project_path,
    sha256_bytes,
    sha256_file,
)
from release_products import TRANSITION_PRODUCT_FIELDS  # noqa: E402
from rebuild_source_spec import load_rebuild_source_spec  # noqa: E402


CLEAN_REBUILD_COMPARISON_FIELDS = (
    "comparison_id",
    "candidate_id",
    "closure_sha256",
    "rebuild_source_spec_sha256",
    "authority_snapshot_spec_sha256",
    "build_a_snapshot_spec_sha256",
    "build_b_snapshot_spec_sha256",
    "authority_project_root",
    "authority_candidate_root",
    "authority_product_manifest",
    "authority_product_manifest_sha256",
    "authority_execution_context",
    "authority_execution_context_sha256",
    "authority_execution_identity",
    "authority_product_transition",
    "authority_product_transition_sha256",
    "retained_namespace_root",
    "comparison_root",
    "build_a_id",
    "build_a_process_id",
    "build_a_project_root",
    "build_a_candidate_root",
    "build_a_product_manifest",
    "build_a_product_manifest_sha256",
    "build_a_execution_context",
    "build_a_execution_context_sha256",
    "build_a_execution_identity",
    "build_a_product_transition",
    "build_a_product_transition_sha256",
    "build_b_id",
    "build_b_process_id",
    "build_b_project_root",
    "build_b_candidate_root",
    "build_b_product_manifest",
    "build_b_product_manifest_sha256",
    "build_b_execution_context",
    "build_b_execution_context_sha256",
    "build_b_execution_identity",
    "build_b_product_transition",
    "build_b_product_transition_sha256",
    "artifact_key",
    "artifact_phase",
    "artifact_role",
    "producer_id",
    "artifact_relpath_authority",
    "artifact_relpath_a",
    "artifact_relpath_b",
    "authority_sha256",
    "authority_bytes",
    "build_a_sha256",
    "build_a_bytes",
    "build_b_sha256",
    "build_b_bytes",
    "equal",
)
BUILD_ID_RE = re.compile(r"^[A-Za-z0-9_.:-]+$")
RETAINED_COMPARISON_RELATIVE = (
    Path("RNA-seq/results/manuscript_release/candidates")
    / CANDIDATE_ID
    / "retained_clean_rebuilds"
)
CANDIDATE_RELATIVE = (
    Path("RNA-seq/results/manuscript_release/candidates") / CANDIDATE_ID
)
ALLOWED_PRODUCERS_BY_PHASE = {
    "REL02": {"rel02", "rel02_build_candidate_tables"},
    "REL03": {"rel03", "rel03_render_candidate_panels"},
    "REL04": {"rel04_build_candidate_manuscript"},
}
EXECUTION_BOUND_PRODUCT_ROLES = {"rel02_state", "rel03_state"}


def require_plain_directory(path: Path, label: str) -> Path:
    """Require one retained, real directory rather than a symlink alias."""
    absolute = path.absolute()
    if absolute.is_symlink() or not absolute.is_dir():
        raise ReleaseContractError(
            f"{label} must be a retained non-symlink directory: {absolute}"
        )
    return absolute.resolve()


def require_no_symlink_components(base: Path, target: Path, label: str) -> None:
    """Reject symlink aliases at every component below a retained namespace."""
    base = base.absolute()
    target = target.absolute()
    try:
        relative = target.relative_to(base)
    except ValueError as error:
        raise ReleaseContractError(
            f"{label} is outside retained namespace {base}: {target}"
        ) from error
    cursor = base
    for component in relative.parts:
        cursor = cursor / component
        if cursor.is_symlink():
            raise ReleaseContractError(
                f"{label} contains a symlinked component: {cursor}"
            )


def retained_paths(
    authority_project: Path,
    comparison_id: str,
) -> tuple[Path, Path, Path, Path, Path, Path, Path, Path]:
    """Resolve the only allowed persistent two-build directory layout."""
    authority = require_plain_directory(authority_project, "authority project root")
    authority_candidate = authority / CANDIDATE_RELATIVE
    namespace = authority / RETAINED_COMPARISON_RELATIVE
    comparison = namespace / comparison_id
    project_a = comparison / "build_a/project"
    project_b = comparison / "build_b/project"
    candidate_a = project_a / CANDIDATE_RELATIVE
    candidate_b = project_b / CANDIDATE_RELATIVE
    for path, label in (
        (authority_candidate, "authority candidate"),
        (namespace, "retained comparison namespace"),
        (comparison, "retained comparison root"),
        (project_a, "retained build-A project"),
        (project_b, "retained build-B project"),
        (candidate_a, "retained build-A candidate"),
        (candidate_b, "retained build-B candidate"),
    ):
        require_no_symlink_components(authority, path, label)
        require_plain_directory(path, label)
    return (
        authority,
        authority_candidate,
        namespace,
        comparison,
        project_a,
        candidate_a,
        project_b,
        candidate_b,
    )


def load_spec(candidate: Path, *, fixture_mode: bool) -> tuple[str, str, str]:
    path = candidate / "manifests/snapshot_spec.json"
    if not path.is_file() or path.is_symlink():
        raise ReleaseContractError(f"clean rebuild lacks a safe snapshot spec: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ReleaseContractError(f"invalid snapshot spec: {path}") from error
    spec = payload.get("spec")
    spec_hash = payload.get("spec_sha256")
    if not isinstance(spec, dict) or not isinstance(spec_hash, str):
        raise ReleaseContractError(f"snapshot spec lacks spec/spec_sha256: {path}")
    require_sha256(spec_hash, "clean-rebuild spec_sha256")
    if sha256_bytes(canonical_json_bytes(spec)) != spec_hash:
        raise ReleaseContractError(f"snapshot spec hash does not reproduce: {path}")
    if spec.get("candidate_id") != CANDIDATE_ID:
        raise ReleaseContractError(f"candidate ID drift in snapshot spec: {path}")
    closure_hash = str(spec.get("closure_sha256") or "")
    require_sha256(closure_hash, "clean-rebuild closure_sha256")
    _, rebuild_hash = load_rebuild_source_spec(candidate, fixture_mode=fixture_mode)
    if spec.get("rebuild_source_spec_sha256") != rebuild_hash:
        raise ReleaseContractError(
            f"snapshot spec does not bind rebuild-source identity: {path}"
        )
    return spec_hash, rebuild_hash, closure_hash


def load_json_file(path: Path, label: str) -> dict[str, object]:
    if not path.is_file() or path.is_symlink():
        raise ReleaseContractError(f"{label} is missing or symlinked: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ReleaseContractError(f"{label} is invalid JSON: {path}") from error
    if not isinstance(payload, dict):
        raise ReleaseContractError(f"{label} must be a JSON object: {path}")
    return payload


def execution_identity(
    candidate: Path,
    label: str,
    fixture_mode: bool = False,
) -> dict[str, str]:
    """Derive an execution identity from the immutable REL01 context."""

    path = candidate / "manifests/execution_context.json"
    payload = load_json_file(path, f"{label} execution context")
    required = {
        "recorded_at_utc",
        "command_line",
        "host",
        "process_id",
        "slurm_job_id",
        "slurm_array_task_id",
        "slurm_step_id",
        "fixture_mode",
    }
    if not required.issubset(payload):
        raise ReleaseContractError(
            f"{label} execution context lacks {sorted(required - set(payload))}"
        )
    if payload["fixture_mode"] is not fixture_mode:
        raise ReleaseContractError(f"{label} clean-rebuild fixture-mode contract drift")
    if "rel01_snapshot_candidate.py" not in str(payload["command_line"]):
        raise ReleaseContractError(
            f"{label} execution context is not a REL01 invocation"
        )
    host = str(payload["host"])
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", host):
        raise ReleaseContractError(f"{label} execution host is unsafe: {host!r}")
    process_id = payload["process_id"]
    if type(process_id) is not int or process_id <= 0:
        raise ReleaseContractError(
            f"{label} execution process_id must be a positive integer"
        )
    timestamp = str(payload["recorded_at_utc"])
    try:
        parsed = datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
    except ValueError as error:
        raise ReleaseContractError(
            f"{label} execution timestamp is invalid: {timestamp!r}"
        ) from error
    if parsed.tzinfo is None:
        raise ReleaseContractError(f"{label} execution timestamp lacks a timezone")
    timestamp_id = parsed.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    slurm_job_id = str(payload["slurm_job_id"])
    slurm_array_task_id = str(payload["slurm_array_task_id"])
    slurm_step_id = str(payload["slurm_step_id"])
    for value, field in (
        (slurm_array_task_id, "SLURM_ARRAY_TASK_ID"),
        (slurm_step_id, "SLURM_STEP_ID"),
    ):
        if value and not BUILD_ID_RE.fullmatch(value):
            raise ReleaseContractError(f"{label} {field} is unsafe: {value!r}")
    if slurm_job_id:
        if not BUILD_ID_RE.fullmatch(slurm_job_id):
            raise ReleaseContractError(
                f"{label} SLURM_JOB_ID is unsafe: {slurm_job_id!r}"
            )
        identity = (
            f"slurm:{slurm_job_id}:"
            f"array:{slurm_array_task_id or 'none'}:"
            f"step:{slurm_step_id or 'none'}"
        )
    else:
        if slurm_array_task_id or slurm_step_id:
            raise ReleaseContractError(
                f"{label} has SLURM task/step identity without SLURM_JOB_ID"
            )
        identity = f"local:{host}:{process_id}:{timestamp_id}"
    digest = sha256_file(path)
    state_path = candidate / "manifests/candidate_state.json"
    state = load_json_file(state_path, f"{label} candidate state")
    if (
        state.get("execution_context_sha256") != digest
        or state.get("canonical_promotion_status") != "not_promoted"
    ):
        raise ReleaseContractError(
            f"{label} candidate state does not bind its REL01 execution context"
        )
    return {
        "path": str(path),
        "sha256": digest,
        "identity": identity,
        "candidate_state_sha256": sha256_file(state_path),
    }


def validate_complete_pre_rel05_tree(
    project: Path,
    candidate: Path,
    label: str,
    fixture_mode: bool = False,
) -> dict[str, object]:
    """Reuse the full read-only REL01--04 validator, excluding REL05 outputs."""

    module = importlib.import_module("rel05_validate_candidate")
    validator = getattr(module, "validate_complete_pre_rel05_candidate", None)
    if not callable(validator):
        raise ReleaseContractError("full pre-REL05 validator is unavailable")
    summary = validator(project, fixture_mode=fixture_mode)
    observed_candidate = summary.get("candidate")
    if not isinstance(observed_candidate, Path) or (
        observed_candidate.resolve() != candidate.resolve()
    ):
        raise ReleaseContractError(
            f"{label} full pre-REL05 validator resolved a different candidate"
        )
    if int(summary.get("n_products", 0)) <= 0 or int(summary.get("n_panels", 0)) <= 0:
        raise ReleaseContractError(f"{label} is not a complete pre-REL05 candidate")
    return summary


def validate_product_transition(
    candidate: Path,
    spec_hash: str,
    execution: dict[str, str],
    label: str,
    fixture_mode: bool = False,
) -> dict[str, str]:
    """Bind deterministic products to this candidate's distinct REL01 run."""

    manifest = candidate / "manifests/rel02_04_product_manifest.tsv"
    manifest_rows = read_tsv_exact(manifest, TRANSITION_PRODUCT_FIELDS)
    path = candidate / "manifests/rel02_04_transition.json"
    transition = load_json_file(path, f"{label} REL02--04 transition")
    expected = {
        "candidate_id": CANDIDATE_ID,
        "transition": "REL02_04_PRODUCTS_FROZEN",
        "from_candidate_state_sha256": execution["candidate_state_sha256"],
        "input_snapshot_spec_sha256": spec_hash,
        "product_manifest_sha256": sha256_file(manifest),
        "n_products": len(manifest_rows),
        "fixture_mode": fixture_mode,
        "canonical_promotion_status": "not_promoted",
        "full_release_pass": False,
    }
    drift = {
        key: {"expected": value, "observed": transition.get(key)}
        for key, value in expected.items()
        if transition.get(key) != value
    }
    if drift:
        raise ReleaseContractError(
            f"{label} REL02--04 transition is not bound to this execution: {drift}"
        )
    return {"path": str(path), "sha256": sha256_file(path)}


def product_inventory(
    project: Path,
    candidate: Path,
    expected_spec_hash: str,
) -> dict[str, dict[str, object]]:
    manifest = candidate / "manifests/rel02_04_product_manifest.tsv"
    rows = read_tsv_exact(manifest, TRANSITION_PRODUCT_FIELDS)
    if not rows:
        raise ReleaseContractError(
            f"clean rebuild product manifest is empty: {manifest}"
        )
    inventory: dict[str, dict[str, object]] = {}
    for row in rows:
        key = row["product_id"]
        if not key or key in inventory:
            raise ReleaseContractError(
                f"blank/duplicate clean-rebuild product ID: {key!r}"
            )
        if row["input_snapshot_spec_sha256"] != expected_spec_hash:
            raise ReleaseContractError(
                f"product {key} points to a different snapshot spec"
            )
        phase = row["phase"]
        producer = row["producer_id"]
        if (
            phase not in ALLOWED_PRODUCERS_BY_PHASE
            or not key.startswith(f"{phase}:")
            or producer not in ALLOWED_PRODUCERS_BY_PHASE[phase]
        ):
            raise ReleaseContractError(
                f"clean-rebuild product phase/producer contract drift: {key}:{phase}:{producer}"
            )
        artifact = resolve_project_path(
            project, row["project_relative_path"], f"clean-rebuild product {key}"
        )
        if not artifact.is_file() or artifact.is_symlink():
            raise ReleaseContractError(
                f"clean-rebuild product is missing or symlinked: {artifact}"
            )
        size = parse_nonnegative_int(row["bytes"], f"clean-rebuild product {key} bytes")
        digest = require_sha256(row["sha256"], f"clean-rebuild product {key} sha256")
        if artifact.stat().st_size != size or sha256_file(artifact) != digest:
            raise ReleaseContractError(f"clean-rebuild product hash/byte drift: {key}")
        if row["artifact_role"] in EXECUTION_BOUND_PRODUCT_ROLES:
            continue
        inventory[key] = {
            "phase": phase,
            "role": row["artifact_role"],
            "producer": producer,
            "path": row["project_relative_path"],
            "sha256": digest,
            "bytes": size,
        }
    return inventory


def derive_clean_rebuild_rows(
    authority_project: Path,
    build_a_id: str,
    build_a_process_id: str,
    build_b_id: str,
    build_b_process_id: str,
    comparison_id: str,
    *,
    fixture_mode: bool = False,
) -> list[dict[str, object]]:
    identifiers = (
        build_a_id,
        build_a_process_id,
        build_b_id,
        build_b_process_id,
        comparison_id,
    )
    if any(not BUILD_ID_RE.fullmatch(value) for value in identifiers):
        raise ReleaseContractError(
            "build/process/comparison IDs must be nonblank stable identifiers"
        )
    if build_a_id == build_b_id or build_a_process_id == build_b_process_id:
        raise ReleaseContractError(
            "clean rebuilds require distinct build and process IDs"
        )
    (
        authority_project,
        authority_candidate,
        retained_namespace,
        comparison_root,
        project_a,
        candidate_a,
        project_b,
        candidate_b,
    ) = retained_paths(authority_project, comparison_id)
    validate_complete_pre_rel05_tree(
        authority_project,
        authority_candidate,
        "authority",
        fixture_mode=fixture_mode,
    )
    validate_complete_pre_rel05_tree(
        project_a,
        candidate_a,
        "retained build A",
        fixture_mode=fixture_mode,
    )
    validate_complete_pre_rel05_tree(
        project_b,
        candidate_b,
        "retained build B",
        fixture_mode=fixture_mode,
    )
    spec_authority, rebuild_authority, closure_authority = load_spec(
        authority_candidate, fixture_mode=fixture_mode
    )
    spec_a, rebuild_a, closure_a = load_spec(candidate_a, fixture_mode=fixture_mode)
    spec_b, rebuild_b, closure_b = load_spec(candidate_b, fixture_mode=fixture_mode)
    if not (
        rebuild_authority == rebuild_a == rebuild_b
        and closure_authority == closure_a == closure_b
    ):
        raise ReleaseContractError(
            "authority and retained builds do not share one rebuild-source identity"
        )
    execution_authority = execution_identity(
        authority_candidate, "authority", fixture_mode=fixture_mode
    )
    execution_a = execution_identity(
        candidate_a, "retained build A", fixture_mode=fixture_mode
    )
    execution_b = execution_identity(
        candidate_b, "retained build B", fixture_mode=fixture_mode
    )
    if (
        len(
            {
                execution_authority["sha256"],
                execution_a["sha256"],
                execution_b["sha256"],
            }
        )
        != 3
        or len(
            {
                execution_authority["identity"],
                execution_a["identity"],
                execution_b["identity"],
            }
        )
        != 3
    ):
        raise ReleaseContractError(
            "authority and retained clean rebuilds require three distinct REL01 executions"
        )
    if (
        build_a_process_id != execution_a["identity"]
        or build_b_process_id != execution_b["identity"]
    ):
        raise ReleaseContractError(
            "CLI process IDs do not equal identities derived from retained execution contexts"
        )
    inventory_authority = product_inventory(
        authority_project, authority_candidate, spec_authority
    )
    inventory_a = product_inventory(project_a, candidate_a, spec_a)
    inventory_b = product_inventory(project_b, candidate_b, spec_b)
    if not (set(inventory_authority) == set(inventory_a) == set(inventory_b)):
        raise ReleaseContractError(
            "authority and retained clean-rebuild artifact-key sets differ"
        )
    manifest_authority = authority_candidate / "manifests/rel02_04_product_manifest.tsv"
    manifest_a = candidate_a / "manifests/rel02_04_product_manifest.tsv"
    manifest_b = candidate_b / "manifests/rel02_04_product_manifest.tsv"
    manifest_authority_sha256 = sha256_file(manifest_authority)
    manifest_a_sha256 = sha256_file(manifest_a)
    manifest_b_sha256 = sha256_file(manifest_b)
    transition_authority = validate_product_transition(
        authority_candidate,
        spec_authority,
        execution_authority,
        "authority",
        fixture_mode=fixture_mode,
    )
    transition_a = validate_product_transition(
        candidate_a,
        spec_a,
        execution_a,
        "retained build A",
        fixture_mode=fixture_mode,
    )
    transition_b = validate_product_transition(
        candidate_b,
        spec_b,
        execution_b,
        "retained build B",
        fixture_mode=fixture_mode,
    )

    rows: list[dict[str, object]] = []
    for key in sorted(inventory_authority):
        authority_item = inventory_authority[key]
        left = inventory_a[key]
        right = inventory_b[key]
        equal = (
            authority_item["phase"] == left["phase"] == right["phase"]
            and authority_item["role"] == left["role"] == right["role"]
            and authority_item["producer"] == left["producer"] == right["producer"]
            and authority_item["path"] == left["path"] == right["path"]
            and authority_item["sha256"] == left["sha256"] == right["sha256"]
            and authority_item["bytes"] == left["bytes"] == right["bytes"]
        )
        if not equal:
            raise ReleaseContractError(
                f"authority/retained clean-rebuild product differs: {key}"
            )
        rows.append(
            {
                "comparison_id": comparison_id,
                "candidate_id": CANDIDATE_ID,
                "closure_sha256": closure_authority,
                "rebuild_source_spec_sha256": rebuild_authority,
                "authority_snapshot_spec_sha256": spec_authority,
                "build_a_snapshot_spec_sha256": spec_a,
                "build_b_snapshot_spec_sha256": spec_b,
                "authority_project_root": str(authority_project),
                "authority_candidate_root": str(authority_candidate),
                "authority_product_manifest": str(manifest_authority),
                "authority_product_manifest_sha256": manifest_authority_sha256,
                "authority_execution_context": execution_authority["path"],
                "authority_execution_context_sha256": execution_authority["sha256"],
                "authority_execution_identity": execution_authority["identity"],
                "authority_product_transition": transition_authority["path"],
                "authority_product_transition_sha256": transition_authority["sha256"],
                "retained_namespace_root": str(retained_namespace),
                "comparison_root": str(comparison_root),
                "build_a_id": build_a_id,
                "build_a_process_id": build_a_process_id,
                "build_a_project_root": str(project_a),
                "build_a_candidate_root": str(candidate_a),
                "build_a_product_manifest": str(manifest_a),
                "build_a_product_manifest_sha256": manifest_a_sha256,
                "build_a_execution_context": execution_a["path"],
                "build_a_execution_context_sha256": execution_a["sha256"],
                "build_a_execution_identity": execution_a["identity"],
                "build_a_product_transition": transition_a["path"],
                "build_a_product_transition_sha256": transition_a["sha256"],
                "build_b_id": build_b_id,
                "build_b_process_id": build_b_process_id,
                "build_b_project_root": str(project_b),
                "build_b_candidate_root": str(candidate_b),
                "build_b_product_manifest": str(manifest_b),
                "build_b_product_manifest_sha256": manifest_b_sha256,
                "build_b_execution_context": execution_b["path"],
                "build_b_execution_context_sha256": execution_b["sha256"],
                "build_b_execution_identity": execution_b["identity"],
                "build_b_product_transition": transition_b["path"],
                "build_b_product_transition_sha256": transition_b["sha256"],
                "artifact_key": key,
                "artifact_phase": str(authority_item["phase"]),
                "artifact_role": str(authority_item["role"]),
                "producer_id": str(authority_item["producer"]),
                "artifact_relpath_authority": str(authority_item["path"]),
                "artifact_relpath_a": str(left["path"]),
                "artifact_relpath_b": str(right["path"]),
                "authority_sha256": str(authority_item["sha256"]),
                "authority_bytes": int(authority_item["bytes"]),
                "build_a_sha256": str(left["sha256"]),
                "build_a_bytes": int(left["bytes"]),
                "build_b_sha256": str(right["sha256"]),
                "build_b_bytes": int(right["bytes"]),
                "equal": "true",
            }
        )
    return rows


def compare_clean_rebuilds(
    authority_project: Path,
    build_a_id: str,
    build_a_process_id: str,
    build_b_id: str,
    build_b_process_id: str,
    comparison_id: str,
    *,
    fixture_mode: bool = False,
) -> list[dict[str, object]]:
    _, authority_candidate, *_ = retained_paths(authority_project, comparison_id)
    output = authority_candidate / "manifests/two_clean_rebuild_comparison.tsv"
    if output.exists() or output.is_symlink():
        raise ReleaseContractError(
            f"refusing to overwrite clean-rebuild report: {output}"
        )
    rows = derive_clean_rebuild_rows(
        authority_project,
        build_a_id,
        build_a_process_id,
        build_b_id,
        build_b_process_id,
        comparison_id,
        fixture_mode=fixture_mode,
    )
    atomic_write_tsv(output, rows, CLEAN_REBUILD_COMPARISON_FIELDS)
    return rows


def revalidate_clean_rebuild_report(
    authority_project: Path,
    report: Path,
    *,
    fixture_mode: bool = False,
) -> list[dict[str, object]]:
    """Recheck all three retained product trees against an existing report."""
    authority = require_plain_directory(authority_project, "authority project root")
    expected_report = (
        authority / CANDIDATE_RELATIVE / "manifests/two_clean_rebuild_comparison.tsv"
    )
    if report.absolute().resolve() != expected_report.resolve():
        raise ReleaseContractError(
            f"clean-rebuild report is not at the authority candidate path: {report}"
        )
    observed = read_tsv_exact(report, CLEAN_REBUILD_COMPARISON_FIELDS)
    if not observed:
        raise ReleaseContractError("clean-rebuild report is empty")
    invariant_fields = (
        "comparison_id",
        "build_a_id",
        "build_a_process_id",
        "build_b_id",
        "build_b_process_id",
    )
    metadata = {field: observed[0][field] for field in invariant_fields}
    if any(
        any(row[field] != metadata[field] for field in invariant_fields)
        for row in observed
    ):
        raise ReleaseContractError(
            "clean-rebuild report metadata varies across products"
        )
    expected = derive_clean_rebuild_rows(
        authority,
        metadata["build_a_id"],
        metadata["build_a_process_id"],
        metadata["build_b_id"],
        metadata["build_b_process_id"],
        metadata["comparison_id"],
        fixture_mode=fixture_mode,
    )
    normalized = [
        {field: str(row[field]) for field in CLEAN_REBUILD_COMPARISON_FIELDS}
        for row in expected
    ]
    if observed != normalized:
        raise ReleaseContractError(
            "clean-rebuild report no longer rederives from authority/build-A/build-B products"
        )
    return expected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority-project-root", type=Path, required=True)
    parser.add_argument("--build-a-id", required=True)
    parser.add_argument(
        "--build-a-process-id",
        required=True,
        help="must equal the identity derived from build A execution_context.json",
    )
    parser.add_argument("--build-b-id", required=True)
    parser.add_argument(
        "--build-b-process-id",
        required=True,
        help="must equal the identity derived from build B execution_context.json",
    )
    parser.add_argument("--comparison-id", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        rows = compare_clean_rebuilds(
            args.authority_project_root,
            args.build_a_id,
            args.build_a_process_id,
            args.build_b_id,
            args.build_b_process_id,
            args.comparison_id,
        )
    except ReleaseContractError as error:
        print(f"CLEAN_REBUILD_COMPARISON_BLOCKED: {error}", file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                "status": "TWO_CLEAN_REBUILDS_IDENTICAL",
                "n_products": len(rows),
                "report": rows[0]["authority_candidate_root"]
                + "/manifests/two_clean_rebuild_comparison.tsv",
                "report_sha256": sha256_file(
                    Path(rows[0]["authority_candidate_root"])
                    / "manifests/two_clean_rebuild_comparison.tsv"
                ),
                "retained_namespace_root": rows[0]["retained_namespace_root"],
                "comparison_root": rows[0]["comparison_root"],
                "rel05_executed": False,
                "canonical_promotion_authorized": False,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
