#!/usr/bin/env python3
"""Materialize one isolated, no-overwrite project for a real clean rebuild.

This does not copy an authority candidate or reproduce unrelated authority
workspace dirtiness.  It starts at the authority's frozen commit, restores
only the exact release producers and source inputs in the deterministic
``rebuild_source_spec``, and verifies all hashes before a READY marker is
written outside the isolated project.  REL01--04 remain separate executions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
if str(PROGRAM_ROOT) not in sys.path:
    sys.path.insert(0, str(PROGRAM_ROOT))

from release_common import (  # noqa: E402
    CANDIDATE_ID,
    CANDIDATE_REL,
    PROTECTED_BASELINE_FIELDS,
    PROTECTED_SCOPE_FIELDS,
    SNAPSHOT_MANIFEST_FIELDS,
    ReleaseContractError,
    atomic_write_json,
    canonical_json_bytes,
    clean_relative_path,
    parse_nonnegative_int,
    read_tsv_exact,
    require_sha256,
    sha256_bytes,
    sha256_file,
    validate_closure,
    verify_protected_baseline,
)
from release_products import BASE_SELECTION_FIELDS, BASE_SNAPSHOT_FIELDS  # noqa: E402
from rebuild_source_spec import load_rebuild_source_spec  # noqa: E402
from rel05_validate_candidate import (  # noqa: E402
    validate_complete_pre_rel05_candidate,
)


MATERIALIZATION_CONTRACT = "plan60_clean_rebuild_materialization_v1"
SAFE_ID_RE = re.compile(r"^[A-Za-z0-9_.:-]+$")
DIRTY_GITLINK_RE = re.compile(
    rb"^\+Subproject commit [0-9a-f]{40}-dirty$", re.MULTILINE
)
RETAINED_ROOT = (
    Path("RNA-seq/results/manuscript_release/candidates")
    / CANDIDATE_ID
    / "retained_clean_rebuilds"
)
PREPARATION_MANIFEST_FIELDS = (
    "artifact_id",
    "relative_path",
    "sha256",
    "bytes",
    "artifact_role",
    "canonical_promotion_authorized",
)


class MaterializationError(ReleaseContractError):
    """Raised before READY when an isolated rebuild cannot be reconstructed."""


def load_json(path: Path, label: str) -> dict[str, object]:
    if not path.is_file() or path.is_symlink():
        raise MaterializationError(f"{label} is missing or symlinked: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise MaterializationError(f"{label} is invalid JSON: {path}") from error
    if not isinstance(payload, dict):
        raise MaterializationError(f"{label} must be a JSON object: {path}")
    return payload


def git_output(
    project: Path, *arguments: str, input_bytes: bytes | None = None
) -> bytes:
    result = subprocess.run(
        ["git", "-C", str(project), *arguments],
        input=input_bytes,
        check=False,
        capture_output=True,
    )
    if result.returncode != 0:
        raise MaterializationError(
            f"git {' '.join(arguments)} failed in {project}: "
            f"{result.stderr.decode('utf-8', errors='replace').strip()}"
        )
    return result.stdout


def authority_spec(
    authority_project: Path,
) -> tuple[Path, dict[str, object], str, dict[str, object], str]:
    candidate = authority_project / CANDIDATE_REL
    payload = load_json(candidate / "manifests/snapshot_spec.json", "authority spec")
    spec = payload.get("spec")
    spec_sha = payload.get("spec_sha256")
    if not isinstance(spec, dict) or not isinstance(spec_sha, str):
        raise MaterializationError("authority snapshot spec lacks spec/spec_sha256")
    require_sha256(spec_sha, "authority spec SHA256")
    if sha256_bytes(canonical_json_bytes(spec)) != spec_sha:
        raise MaterializationError("authority snapshot spec hash does not reproduce")
    if (
        spec.get("candidate_id") != CANDIDATE_ID
        or spec.get("fixture_mode") is not False
    ):
        raise MaterializationError("authority is not the real locked candidate")
    rebuild_spec, rebuild_spec_sha = load_rebuild_source_spec(
        candidate, fixture_mode=False
    )
    if spec.get("rebuild_source_spec_sha256") != rebuild_spec_sha:
        raise MaterializationError(
            "authority snapshot does not bind the rebuild source specification"
        )
    return candidate, spec, spec_sha, rebuild_spec, rebuild_spec_sha


def validate_preparation(
    authority_project: Path,
    preparation: Path,
    candidate: Path,
) -> tuple[Path, Path]:
    expected_parent = (
        authority_project / "scripts/manuscript/program_context_v2/coordinator/prepared"
    ).resolve()
    prepared = preparation.resolve()
    try:
        relative = prepared.relative_to(expected_parent)
    except ValueError as error:
        raise MaterializationError(
            f"preparation is outside the coordinator-owned root: {prepared}"
        ) from error
    if len(relative.parts) != 1 or prepared.is_symlink() or not prepared.is_dir():
        raise MaterializationError(
            "preparation must be one real decision directory below prepared/"
        )
    ready = load_json(prepared / "COORDINATOR_PREPARATION_READY", "preparation READY")
    required = {
        "status": "COORDINATOR_PREPARATION_READY_NOT_EXECUTED_NOT_PROMOTED",
        "rel00_05_executed": False,
        "candidate_root_created": False,
        "canonical_promotion_authorized": False,
    }
    drift = {
        key: {"expected": value, "observed": ready.get(key)}
        for key, value in required.items()
        if ready.get(key) != value
    }
    if drift:
        raise MaterializationError(f"preparation READY contract drift: {drift}")
    manifest = prepared / "coordinator_preparation_manifest.tsv"
    if ready.get("preparation_manifest_sha256") != sha256_file(manifest):
        raise MaterializationError("preparation manifest hash drift")
    manifest_rows = read_tsv_exact(manifest, PREPARATION_MANIFEST_FIELDS)
    manifested_paths: set[str] = set()
    for row in manifest_rows:
        relative_path = clean_relative_path(
            row["relative_path"], f"preparation artifact {row['artifact_id']}"
        )
        if (
            relative_path in manifested_paths
            or row["canonical_promotion_authorized"] != "false"
        ):
            raise MaterializationError(
                f"duplicate/promoting preparation artifact: {relative_path}"
            )
        manifested_paths.add(relative_path)
        artifact = prepared / relative_path
        if (
            not artifact.is_file()
            or artifact.is_symlink()
            or sha256_file(artifact) != row["sha256"]
            or artifact.stat().st_size
            != parse_nonnegative_int(
                row["bytes"], f"preparation artifact {relative_path} bytes"
            )
        ):
            raise MaterializationError(
                f"preparation artifact hash/byte drift: {relative_path}"
            )
    discovered = {
        path.relative_to(prepared).as_posix()
        for path in prepared.rglob("*")
        if path.is_file()
        and path.name
        not in {"coordinator_preparation_manifest.tsv", "COORDINATOR_PREPARATION_READY"}
    }
    if manifested_paths != discovered:
        raise MaterializationError(
            "preparation file set differs from its manifest: "
            f"missing={sorted(manifested_paths - discovered)}; "
            f"extra={sorted(discovered - manifested_paths)}"
        )
    closure = prepared / "workstream_closure.tsv"
    selection = prepared / "base_input_selection.tsv"
    if (
        ready.get("closure_sha256") != sha256_file(closure)
        or ready.get("base_selection_sha256") != sha256_file(selection)
        or sha256_file(closure)
        != sha256_file(candidate / "manifests/workstream_closure.tsv")
        or sha256_file(selection)
        != sha256_file(candidate / "manifests/base_input_selection.tsv")
    ):
        raise MaterializationError(
            "preparation closure/selection does not equal the authority snapshot"
        )
    return closure, selection


def retained_destination(
    authority_project: Path,
    comparison_id: str,
    build_slot: str,
) -> tuple[Path, Path]:
    if not SAFE_ID_RE.fullmatch(comparison_id):
        raise MaterializationError("comparison ID is blank or unsafe")
    if build_slot not in {"build_a", "build_b"}:
        raise MaterializationError("build slot must be build_a or build_b")
    slot_root = authority_project / RETAINED_ROOT / comparison_id / build_slot
    project = slot_root / "project"
    ready = slot_root / "MATERIALIZATION_READY.json"
    if project.exists() or project.is_symlink() or ready.exists() or ready.is_symlink():
        raise MaterializationError(
            f"refusing to overwrite retained build slot: {slot_root}"
        )
    return project, ready


def authority_repository_provenance(
    spec: dict[str, object],
    rebuild_spec: dict[str, object],
) -> dict[str, str]:
    """Validate frozen global Git provenance without replaying authority dirt."""

    repository = spec.get("repository")
    if not isinstance(repository, dict):
        raise MaterializationError("authority spec lacks repository state")
    commit = str(repository.get("repository_commit") or "")
    tracked_diff_sha = str(repository.get("tracked_diff_sha256") or "")
    if not re.fullmatch(r"[0-9a-f]{40,64}", commit):
        raise MaterializationError(f"unsafe authority commit: {commit!r}")
    require_sha256(tracked_diff_sha, "authority tracked-diff SHA256")
    if rebuild_spec.get("repository_commit") != commit:
        raise MaterializationError(
            "rebuild source commit differs from frozen authority provenance"
        )
    return {
        "authority_repository_commit": commit,
        # Provenance only: deliberately excluded from rebuild-source identity
        # and never applied to a retained scientific build.
        "authority_tracked_diff_sha256": tracked_diff_sha,
    }


def copy_frozen_file(
    source: Path,
    destination: Path,
    sha256: str,
    size: int,
    *,
    allow_isolated_overlay: bool = False,
) -> None:
    require_sha256(sha256, f"rehydrated file {destination}")
    if (
        not source.is_file()
        or source.is_symlink()
        or source.stat().st_size != size
        or sha256_file(source) != sha256
    ):
        raise MaterializationError(f"frozen source hash/byte drift: {source}")
    if destination.exists() or destination.is_symlink():
        matches = (
            not destination.is_symlink()
            and destination.is_file()
            and destination.stat().st_size == size
            and sha256_file(destination) == sha256
        )
        if matches:
            return
        if not allow_isolated_overlay:
            raise MaterializationError(
                f"worktree file conflicts with frozen source: {destination}"
            )
        if destination.is_symlink() or not destination.is_file():
            raise MaterializationError(
                f"isolated overlay target is not a regular file: {destination}"
            )
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
    try:
        shutil.copyfile(source, temporary)
        if temporary.stat().st_size != size or sha256_file(temporary) != sha256:
            raise MaterializationError(
                f"rehydrated temporary copy failed verification: {destination}"
            )
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()
    if (
        destination.is_symlink()
        or not destination.is_file()
        or destination.stat().st_size != size
        or sha256_file(destination) != sha256
    ):
        raise MaterializationError(
            f"rehydrated copy failed verification: {destination}"
        )


def build_rehydration_plan(
    authority_project: Path,
    authority_candidate: Path,
    preparation_closure: Path,
    preparation_selection: Path,
) -> tuple[dict[str, tuple[Path, str, int]], str, str]:
    planned: dict[str, tuple[Path, str, int]] = {}

    def add(relative_text: str, source: Path, digest: str, size_value: object) -> None:
        relative = clean_relative_path(relative_text, "rehydrated source path")
        size = parse_nonnegative_int(size_value, f"rehydrated {relative} bytes")
        previous = planned.get(relative)
        item = (source, require_sha256(digest, f"rehydrated {relative}"), size)
        if previous is not None and previous[1:] != item[1:]:
            raise MaterializationError(f"conflicting rehydration plan: {relative}")
        planned[relative] = item

    snapshot_rows = read_tsv_exact(
        authority_candidate / "manifests/input_snapshot_manifest.tsv",
        SNAPSHOT_MANIFEST_FIELDS,
    )
    for row in snapshot_rows:
        add(
            row["source_path_provenance"],
            authority_candidate / row["snapshot_relpath"],
            row["source_sha256"],
            row["source_bytes"],
        )
    base_rows = read_tsv_exact(
        authority_candidate / "manifests/base_input_snapshot_manifest.tsv",
        BASE_SNAPSHOT_FIELDS,
    )
    for row in base_rows:
        add(
            row["source_path_provenance"],
            authority_candidate / row["snapshot_path"],
            row["source_sha256"],
            row["source_bytes"],
        )

    # REL01 requires the BASE selection beside the closure.  It is copied into
    # the authority candidate by REL02 but is not itself a BASE artifact row.
    closure_relative = preparation_closure.relative_to(authority_project).as_posix()
    selection_relative = preparation_selection.relative_to(authority_project).as_posix()
    add(
        closure_relative,
        authority_candidate / "manifests/workstream_closure.tsv",
        sha256_file(preparation_closure),
        preparation_closure.stat().st_size,
    )
    add(
        selection_relative,
        authority_candidate / "manifests/base_input_selection.tsv",
        sha256_file(preparation_selection),
        preparation_selection.stat().st_size,
    )

    scopes_path = authority_candidate / "manifests/protected_scopes.tsv"
    baseline_path = authority_candidate / "manifests/protected_release_baseline.tsv"
    scopes = read_tsv_exact(scopes_path, PROTECTED_SCOPE_FIELDS)
    baseline = read_tsv_exact(baseline_path, PROTECTED_BASELINE_FIELDS)
    roots = {
        row["scope_id"]: clean_relative_path(row["root_path"], "scope root")
        for row in scopes
    }
    for row in baseline:
        root = roots.get(row["scope_id"])
        if root is None:
            raise MaterializationError(
                f"protected baseline has unknown scope: {row['scope_id']}"
            )
        relative = (Path(root) / row["relative_path"]).as_posix()
        add(
            relative,
            authority_project / relative,
            row["sha256"],
            row["bytes"],
        )

    return planned, closure_relative, selection_relative


def rehydrate_snapshot_sources(
    planned: dict[str, tuple[Path, str, int]],
    destination_project: Path,
    closure_relative: str,
    selection_relative: str,
) -> tuple[Path, Path, int]:
    for relative, (source, digest, size) in sorted(planned.items()):
        copy_frozen_file(
            source,
            destination_project / relative,
            digest,
            size,
            allow_isolated_overlay=True,
        )
    return (
        destination_project / closure_relative,
        destination_project / selection_relative,
        len(planned),
    )


def restore_frozen_producers(
    authority_candidate: Path,
    destination_project: Path,
) -> int:
    fields = ("producer_id", "repository_path", "sha256", "bytes", "snapshot_path")
    rows = read_tsv_exact(
        authority_candidate / "manifests/producer_script_manifest.tsv", fields
    )
    for row in rows:
        relative = clean_relative_path(row["repository_path"], "producer path")
        copy_frozen_file(
            authority_candidate / row["snapshot_path"],
            destination_project / relative,
            row["sha256"],
            parse_nonnegative_int(row["bytes"], f"producer {relative} bytes"),
            allow_isolated_overlay=True,
        )
    return len(rows)


def assert_symlink_free(root: Path, label: str) -> None:
    if root.is_symlink():
        raise MaterializationError(f"{label} root is a symlink: {root}")
    for path in root.rglob("*"):
        if path.is_symlink():
            raise MaterializationError(f"{label} contains a symlink: {path}")


def sparse_pattern(relative: str) -> str:
    """Return one exact no-cone sparse-checkout pattern."""

    relative = clean_relative_path(relative, "sparse checkout path")
    if "\n" in relative or "\r" in relative:
        raise MaterializationError(f"unsafe newline in sparse path: {relative!r}")
    escaped = re.sub(r"([\\*?\[\]])", r"\\\1", relative)
    return f"/{escaped}"


def create_sparse_worktree(
    authority: Path,
    destination: Path,
    commit: str,
    required_files: set[str],
) -> None:
    """Create a detached exact-file worktree without tracked symlink spillover."""

    result = subprocess.run(
        [
            "git",
            "-C",
            str(authority),
            "worktree",
            "add",
            "--detach",
            "--no-checkout",
            str(destination),
            commit,
        ],
        check=False,
        capture_output=True,
    )
    if result.returncode != 0:
        raise MaterializationError(
            "git worktree add failed: "
            + result.stderr.decode("utf-8", errors="replace").strip()
        )
    patterns = (
        "\n".join(sparse_pattern(relative) for relative in sorted(required_files))
        + "\n"
    )
    git_output(destination, "sparse-checkout", "init", "--no-cone")
    git_output(
        destination,
        "sparse-checkout",
        "set",
        "--stdin",
        input_bytes=patterns.encode("utf-8"),
    )
    # ``worktree add --no-checkout`` leaves the index/worktree empty; populate
    # it only after the exact no-cone file patterns are installed.
    git_output(destination, "checkout", "--detach", commit)
    assert_symlink_free(destination, "sparse retained source project")


def validate_isolated_repository_state(
    destination: Path,
    expected_commit: str,
    allowed_overlay_paths: set[str],
) -> dict[str, object]:
    """Require a clean commit plus only exact frozen release-source overlays."""

    observed_commit = git_output(destination, "rev-parse", "HEAD").decode().strip()
    if observed_commit != expected_commit:
        raise MaterializationError(
            "isolated worktree commit differs from rebuild source specification"
        )
    patch = git_output(destination, "diff", "--binary", "HEAD", "--")
    if DIRTY_GITLINK_RE.search(patch):
        raise MaterializationError(
            "isolated rebuild contains a dirty gitlink; only frozen file overlays "
            "are permitted"
        )
    changed_raw = git_output(destination, "diff", "--name-only", "-z", "HEAD", "--")
    untracked_raw = git_output(destination, "ls-files", "--others", "-z", "--")
    changed = {
        clean_relative_path(item.decode("utf-8"), "isolated tracked overlay")
        for item in changed_raw.split(b"\0")
        if item
    }
    untracked = {
        clean_relative_path(item.decode("utf-8"), "isolated untracked source")
        for item in untracked_raw.split(b"\0")
        if item
    }
    unexpected = (changed | untracked) - allowed_overlay_paths
    if unexpected:
        raise MaterializationError(
            "isolated rebuild contains paths outside the frozen source overlay: "
            f"{sorted(unexpected)}"
        )
    return {
        "repository_commit": observed_commit,
        "isolated_tracked_diff_sha256": hashlib.sha256(patch).hexdigest(),
        "isolated_tracked_overlay_paths": sorted(changed),
        "isolated_untracked_source_paths": sorted(untracked),
    }


def materialize(
    authority_project: Path,
    preparation_root: Path,
    comparison_id: str,
    build_slot: str,
) -> dict[str, object]:
    authority = authority_project.resolve()
    if (
        authority.is_symlink()
        or not authority.is_dir()
        or not (authority / ".git").exists()
    ):
        raise MaterializationError("authority must be a real Git working tree")
    candidate, spec, spec_sha, rebuild_spec, rebuild_spec_sha = authority_spec(
        authority
    )
    complete = validate_complete_pre_rel05_candidate(authority, fixture_mode=False)
    if complete.get("candidate") != candidate:
        raise MaterializationError(
            "authority full pre-REL05 validation resolved a different candidate"
        )
    closure, selection = validate_preparation(authority, preparation_root, candidate)
    repository_provenance = authority_repository_provenance(spec, rebuild_spec)
    destination, ready_path = retained_destination(authority, comparison_id, build_slot)
    rehydration_plan, closure_relative, selection_relative = build_rehydration_plan(
        authority,
        candidate,
        closure,
        selection,
    )
    producer_rows = read_tsv_exact(
        candidate / "manifests/producer_script_manifest.tsv",
        ("producer_id", "repository_path", "sha256", "bytes", "snapshot_path"),
    )
    producer_paths = {
        clean_relative_path(row["repository_path"], "producer path")
        for row in producer_rows
    }
    sparse_files = {
        *rehydration_plan,
        *producer_paths,
    }

    # No filesystem mutation occurs before every authority/source/spec gate
    # above passes.  A later failure intentionally leaves an inspectable,
    # unsealed worktree; this function never deletes or overwrites it.
    destination.parent.mkdir(parents=True, exist_ok=True)
    commit = str(rebuild_spec["repository_commit"])
    create_sparse_worktree(authority, destination, commit, sparse_files)
    n_producers = restore_frozen_producers(candidate, destination)
    destination_closure, destination_selection, n_sources = rehydrate_snapshot_sources(
        rehydration_plan,
        destination,
        closure_relative,
        selection_relative,
    )
    assert_symlink_free(destination, "rehydrated retained source project")
    isolated_repository = validate_isolated_repository_state(
        destination,
        commit,
        sparse_files,
    )
    validate_closure(destination, destination_closure)
    selection_rows = read_tsv_exact(destination_selection, BASE_SELECTION_FIELDS)
    for row in selection_rows:
        source = destination / clean_relative_path(
            row["source_path"], f"BASE source {row['artifact_id']}"
        )
        if (
            not source.is_file()
            or source.is_symlink()
            or sha256_file(source) != row["source_sha256"]
            or source.stat().st_size
            != parse_nonnegative_int(
                row["source_bytes"], f"BASE source {row['artifact_id']} bytes"
            )
        ):
            raise MaterializationError(
                f"rehydrated BASE source drift: {row['artifact_id']}"
            )
    verify_protected_baseline(
        destination,
        destination / closure.parent.relative_to(authority) / "protected_scopes.tsv",
        destination
        / closure.parent.relative_to(authority)
        / "protected_release_baseline.tsv",
    )
    ready = {
        "contract_version": MATERIALIZATION_CONTRACT,
        "status": "SOURCE_PROJECT_READY_FOR_FRESH_REL01_NOT_EXECUTED",
        "candidate_id": CANDIDATE_ID,
        "comparison_id": comparison_id,
        "build_slot": build_slot,
        "project_root": str(destination),
        "closure_path": str(destination_closure),
        "base_selection_path": str(destination_selection),
        "base_selection_sha256": sha256_file(destination_selection),
        "authority_spec_sha256": spec_sha,
        "rebuild_source_spec_sha256": rebuild_spec_sha,
        **repository_provenance,
        **isolated_repository,
        "n_rehydrated_sources": n_sources,
        "n_frozen_producers": n_producers,
        "rel01_04_executed": False,
        "candidate_tree_copied": False,
        "canonical_promotion_authorized": False,
    }
    atomic_write_json(ready_path, ready)
    return ready


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority-project-root", type=Path, required=True)
    parser.add_argument("--preparation-root", type=Path, required=True)
    parser.add_argument("--comparison-id", required=True)
    parser.add_argument("--build-slot", choices=("build_a", "build_b"), required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        result = materialize(
            args.authority_project_root,
            args.preparation_root,
            args.comparison_id,
            args.build_slot,
        )
    except (MaterializationError, ReleaseContractError) as error:
        print(f"CLEAN_REBUILD_MATERIALIZATION_BLOCKED: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
