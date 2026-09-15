#!/usr/bin/env python3
"""Shared dependency-free contracts for the program-context v2 release scaffold."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass
from datetime import date
from pathlib import Path, PurePosixPath
from typing import Iterable, Mapping, Sequence


CANDIDATE_ID = "program-context-v2-candidate-2026-08-07"
CANDIDATE_REL = Path("RNA-seq/results/manuscript_release/candidates") / CANDIDATE_ID
REQUIRED_WORKSTREAMS = ("PLAN13", "PLAN20", "PLAN30", "PLAN40", "PLAN50")
WORKSTREAM_LANES = {
    "PLAN13": "SP-INT",
    "PLAN20": "HS-V2",
    "PLAN30": "GEN",
    "PLAN40": "HLF",
    "PLAN50": "PASS",
}
ALLOWED_TERMINAL_STATES = {
    "accepted_main",
    "accepted_supplement",
    "skipped_source_gate",
    "rejected_qc",
}
CLOSURE_FIELDS = (
    "workstream_id",
    "terminal_state",
    "gate_id",
    "gate_verdict",
    "artifact_manifest",
    "manifest_sha256",
    "include_main",
    "include_supplement",
    "figure_role",
    "allowed_wording",
    "prohibited_wording",
    "exclusion_reason",
    "adjudicator",
    "date",
)
ARTIFACT_FIELDS = (
    "artifact_id",
    "source_path",
    "source_sha256",
    "source_bytes",
    "snapshot_relpath",
    "artifact_role",
    "downstream_read_path",
)
PROTECTED_SCOPE_FIELDS = (
    "scope_id",
    "protection_class",
    "root_path",
)
PROTECTED_BASELINE_FIELDS = (
    "scope_id",
    "relative_path",
    "sha256",
    "bytes",
)
SNAPSHOT_MANIFEST_FIELDS = (
    "workstream_id",
    "artifact_id",
    "artifact_role",
    "source_path_provenance",
    "source_sha256",
    "source_bytes",
    "snapshot_relpath",
    "snapshot_sha256",
    "snapshot_bytes",
)
DOWNSTREAM_FIELDS = (
    "consumer_id",
    "workstream_id",
    "artifact_id",
    "snapshot_path",
    "sha256",
    "bytes",
)
SCIENTIFIC_REGISTRY_FIELDS = (
    "check_id",
    "status",
    "report_path",
    "report_sha256",
    "warnings_count",
    "adjudication_path",
    "adjudication_sha256",
)
REQUIRED_SCIENTIFIC_CHECKS = (
    "hotspot_registry_consistency",
    "spatial_assay_native_consistency",
    "numbers_claim_ledger",
    "figure_panel_manifest",
    "passport_semantics",
    "rebuild_determinism",
    "null_compatibility",
    "manuscript_retired_claims",
    "biological_unit_audit",
    "multiplicity_universe_audit",
)
# One foundational producer universe for both the scientific-report builder
# and REL05 verification.  Ordering is deterministic for serialized bindings.
SCIENTIFIC_REPORT_PRODUCER_PATHS = (
    "scripts/manuscript/program_context_v2/coordinator/build_scientific_validation_registry.py",
    "scripts/manuscript/program_context_v2/coordinator/coordinator_contract.py",
    "scripts/manuscript/program_context_v2/fibrosis_candidate_contract.py",
    "scripts/manuscript/program_context_v2/rel02_build_candidate_tables.py",
    "scripts/manuscript/program_context_v2/rel03_render_candidate_panels.py",
    "scripts/manuscript/program_context_v2/rel04_build_candidate_manuscript.py",
    "scripts/manuscript/program_context_v2/rel05_validate_candidate.py",
    "scripts/manuscript/program_context_v2/rebuild_source_spec.py",
    "scripts/manuscript/program_context_v2/release_common.py",
    "scripts/manuscript/program_context_v2/release_products.py",
)
# Stable logical producer IDs recorded by the sealed Plan 50 payload manifest.
# These are contract identifiers, independent of a frozen copy's physical
# execution path.  The payload is allowed to name exactly this producer set.
PLAN50_MANIFEST_PRODUCER_PATHS = (
    "scripts/portal/adjudicate_gen_ensembl_identity.py",
    "scripts/portal/build_evidence_passport_bundle.py",
    "scripts/portal/build_evidence_passport_ui.py",
    "scripts/portal/prepare_passport_attestations.py",
    "scripts/portal/prepare_passport_selection.py",
    "scripts/portal/render_evidence_passport_ui.py",
    "scripts/portal/validate_evidence_passports.py",
)
# REL01 additionally freezes the scripts and wrappers needed to rebuild and
# independently test Plan 50 in the source-resolving project environment.
# The list is explicit so filename discovery cannot expand the trust boundary.
PLAN50_FROZEN_PRODUCER_PATHS = (
    *PLAN50_MANIFEST_PRODUCER_PATHS,
    "scripts/portal/generate_evidence_passports.py",
    "scripts/portal/run_gene_catalog.sbatch",
    "scripts/portal/run_gen_identity_adjudication.sbatch",
    "scripts/portal/test_evidence_passport_bundle.py",
)
REQUIRED_CANDIDATE_DIRS = (
    "inputs/BASE",
    "inputs/HS-V2",
    "inputs/GEN",
    "inputs/SP-INT",
    "inputs/HLF",
    "inputs/PASS",
    "tables",
    "figure_sources",
    "manifests/upstream",
    "claims",
    "environments",
    "logs",
)
FORBIDDEN_WRITE_RELS = (
    Path("figures/main"),
    Path("docs/PAPER.md"),
    Path("docs/STATUS.md"),
    Path("docs/RESULTS.md"),
    Path("docs/ROADMAP.md"),
    Path("docs/manuscript/release"),
    Path("docs/manuscript/draft"),
    Path("docs/manuscript/METHODS.md"),
    Path("docs/manuscript/working"),
    Path("RNA-seq/results/manuscript_release/2026-07-10-r1"),
    Path("RNA-seq/results/manuscript_release/2026-07-15-r2"),
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class ReleaseContractError(RuntimeError):
    """Raised whenever a release boundary cannot be proven."""


@dataclass(frozen=True)
class Artifact:
    workstream_id: str
    artifact_id: str
    source_path_text: str
    source_path: Path
    source_sha256: str
    source_bytes: int
    snapshot_relpath: str
    artifact_role: str
    downstream_read_path: str


@dataclass(frozen=True)
class WorkstreamClosure:
    row: dict[str, str]
    manifest_path_text: str
    manifest_path: Path
    manifest_sha256: str
    artifacts: tuple[Artifact, ...]


@dataclass(frozen=True)
class ClosureBundle:
    closure_path: Path
    closure_sha256: str
    workstreams: tuple[WorkstreamClosure, ...]

    @property
    def artifacts(self) -> tuple[Artifact, ...]:
        return tuple(
            artifact
            for workstream in self.workstreams
            for artifact in workstream.artifacts
        )


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def sha256_file(path: Path, chunk_bytes: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_json_bytes(payload: object) -> bytes:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def parse_nonnegative_int(value: str, context: str) -> int:
    try:
        result = int(value)
    except ValueError as error:
        raise ReleaseContractError(f"{context} is not an integer: {value!r}") from error
    if result < 0:
        raise ReleaseContractError(f"{context} is negative: {result}")
    return result


def parse_bool(value: str, context: str) -> bool:
    if value == "true":
        return True
    if value == "false":
        return False
    raise ReleaseContractError(f"{context} must be exactly true or false: {value!r}")


def require_sha256(value: str, context: str) -> str:
    if not SHA256_RE.fullmatch(value):
        raise ReleaseContractError(f"{context} is not lowercase SHA256: {value!r}")
    return value


def read_tsv_exact(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ReleaseContractError(
            f"required regular TSV is missing or symlinked: {path}"
        )
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise ReleaseContractError(
                f"TSV schema mismatch for {path}: expected {tuple(fields)}, "
                f"observed {tuple(reader.fieldnames or ())}"
            )
        return [dict(row) for row in reader]


def atomic_write_tsv(
    path: Path,
    rows: Iterable[Mapping[str, object]],
    fields: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def atomic_write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def resolve_project_path(project_root: Path, path_text: str, context: str) -> Path:
    project = project_root.resolve()
    raw = Path(path_text)
    candidate = raw if raw.is_absolute() else project / raw
    # Normalize ``..`` without resolving symlinks first.  Calling ``resolve()``
    # immediately would erase whether the manifest named a symlink and would let
    # a link masquerade as an read-only regular-file input.
    lexical = Path(os.path.abspath(candidate))
    if not is_relative_to(lexical, project):
        raise ReleaseContractError(f"{context} escapes project root: {path_text}")
    cursor = project
    for part in lexical.relative_to(project).parts:
        cursor /= part
        if cursor.is_symlink():
            raise ReleaseContractError(
                f"{context} must use a non-symlink path component: {cursor}"
            )
    resolved = lexical.resolve()
    if not is_relative_to(resolved, project):
        raise ReleaseContractError(f"{context} escapes project root: {path_text}")
    return resolved


def project_relative(project_root: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(project_root.resolve()).as_posix()
    except ValueError as error:
        raise ReleaseContractError(f"path escapes project root: {path}") from error


def clean_relative_path(value: str, context: str) -> str:
    if "\\" in value:
        raise ReleaseContractError(f"{context} contains a backslash: {value!r}")
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or not path.parts
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise ReleaseContractError(f"{context} is not a clean relative path: {value!r}")
    return path.as_posix()


def candidate_root(project_root: Path) -> Path:
    # Keep the lexical path so a pre-existing symlink at the candidate location
    # cannot disappear through ``resolve()`` before the guard inspects it.
    return project_root.resolve() / CANDIDATE_REL


def assert_candidate_root(project_root: Path, proposed: Path) -> Path:
    project = project_root.resolve()
    expected = Path(os.path.abspath(candidate_root(project)))
    observed = Path(os.path.abspath(proposed))
    if observed != expected:
        raise ReleaseContractError(
            f"candidate root must be exactly {expected}; observed {observed}"
        )
    cursor = project
    for part in observed.relative_to(project).parts:
        cursor /= part
        if cursor.is_symlink():
            raise ReleaseContractError(
                f"candidate root must use non-symlink path components: {cursor}"
            )
    for relative in FORBIDDEN_WRITE_RELS:
        forbidden = Path(os.path.abspath(project / relative))
        if is_relative_to(observed, forbidden) or is_relative_to(forbidden, observed):
            raise ReleaseContractError(
                f"candidate root overlaps forbidden output: {forbidden}"
            )
    return observed


def allowed_workstream_roots(
    project_root: Path, workstream_id: str
) -> tuple[Path, ...]:
    project = project_root.resolve()
    multimodal = (
        project / "Analysis/Multimodal_Program_Projection/candidates" / CANDIDATE_ID
    )
    roots = {
        "PLAN13": (
            project / "Analysis/Spatial/candidates" / CANDIDATE_ID,
            multimodal / "spatial_context",
            multimodal / "spatial_context_semantic_v2_2026-08-08",
        ),
        "PLAN20": (multimodal / "hotspot",),
        "PLAN30": (multimodal / "genetics_context",),
        "PLAN40": (multimodal / "myojin_hlf",),
        "PLAN50": (
            project / "RNA-seq/results/evidence_passports/candidates" / CANDIDATE_ID,
        ),
    }
    return tuple(path.resolve() for path in roots[workstream_id])


def require_within_any(path: Path, roots: Sequence[Path], context: str) -> None:
    if not any(is_relative_to(path.resolve(), root.resolve()) for root in roots):
        raise ReleaseContractError(
            f"{context} is outside its workstream-owned roots: {path}"
        )


def validate_terminal_inclusion(row: Mapping[str, str]) -> None:
    workstream = row["workstream_id"]
    state = row["terminal_state"]
    include_main = parse_bool(row["include_main"], f"{workstream} include_main")
    include_supplement = parse_bool(
        row["include_supplement"], f"{workstream} include_supplement"
    )
    if state == "accepted_main" and not include_main:
        raise ReleaseContractError(f"{workstream} accepted_main must include main")
    if state == "accepted_supplement" and (include_main or not include_supplement):
        raise ReleaseContractError(
            f"{workstream} accepted_supplement must be supplement-only"
        )
    if state == "skipped_source_gate" and include_main:
        raise ReleaseContractError(f"{workstream} source-gate skip cannot enter main")
    if state == "rejected_qc" and (include_main or include_supplement):
        raise ReleaseContractError(f"{workstream} rejected_qc cannot enter figures")
    if (
        state in {"skipped_source_gate", "rejected_qc"}
        and not row["exclusion_reason"].strip()
    ):
        raise ReleaseContractError(f"{workstream} terminal exclusion lacks a reason")


def validate_artifact_manifest(
    project_root: Path,
    workstream_id: str,
    manifest_path: Path,
) -> tuple[Artifact, ...]:
    rows = read_tsv_exact(manifest_path, ARTIFACT_FIELDS)
    if not rows:
        raise ReleaseContractError(f"empty terminal artifact manifest: {manifest_path}")
    artifacts = []
    artifact_ids = set()
    snapshot_paths = set()
    lane = WORKSTREAM_LANES[workstream_id]
    roots = allowed_workstream_roots(project_root, workstream_id)
    for row in rows:
        artifact_id = row["artifact_id"].strip()
        if not artifact_id or artifact_id in artifact_ids:
            raise ReleaseContractError(
                f"blank or duplicate artifact_id in {manifest_path}: {artifact_id!r}"
            )
        if artifact_id == "terminal_artifact_manifest":
            raise ReleaseContractError(
                f"reserved artifact_id in {manifest_path}: {artifact_id}"
            )
        artifact_ids.add(artifact_id)
        source = resolve_project_path(
            project_root,
            row["source_path"],
            f"{workstream_id} {artifact_id} source",
        )
        require_within_any(source, roots, f"{workstream_id} {artifact_id} source")
        if not source.is_file() or source.is_symlink():
            raise ReleaseContractError(
                f"artifact source must be a regular non-symlink file: {source}"
            )
        source_hash = require_sha256(
            row["source_sha256"], f"{workstream_id} {artifact_id} source_sha256"
        )
        source_bytes = parse_nonnegative_int(
            row["source_bytes"], f"{workstream_id} {artifact_id} source_bytes"
        )
        if source.stat().st_size != source_bytes or sha256_file(source) != source_hash:
            raise ReleaseContractError(
                f"artifact source hash/byte mismatch: {workstream_id} {artifact_id}"
            )
        snapshot = clean_relative_path(
            row["snapshot_relpath"],
            f"{workstream_id} {artifact_id} snapshot_relpath",
        )
        if not snapshot.startswith(f"inputs/{lane}/"):
            raise ReleaseContractError(
                f"{workstream_id} artifact must snapshot beneath inputs/{lane}/: {snapshot}"
            )
        if snapshot in snapshot_paths:
            raise ReleaseContractError(
                f"duplicate snapshot path within {workstream_id}: {snapshot}"
            )
        snapshot_paths.add(snapshot)
        downstream = clean_relative_path(
            row["downstream_read_path"],
            f"{workstream_id} {artifact_id} downstream_read_path",
        )
        if downstream != snapshot:
            raise ReleaseContractError(
                f"downstream reads must name the copied snapshot, not a live path: "
                f"{workstream_id} {artifact_id}"
            )
        role = row["artifact_role"].strip()
        if not role:
            raise ReleaseContractError(
                f"blank artifact role: {workstream_id} {artifact_id}"
            )
        artifacts.append(
            Artifact(
                workstream_id=workstream_id,
                artifact_id=artifact_id,
                source_path_text=row["source_path"],
                source_path=source,
                source_sha256=source_hash,
                source_bytes=source_bytes,
                snapshot_relpath=snapshot,
                artifact_role=role,
                downstream_read_path=downstream,
            )
        )
    terminal_count = sum(
        artifact.artifact_role == "terminal_verdict" for artifact in artifacts
    )
    if terminal_count != 1:
        raise ReleaseContractError(
            f"{workstream_id} requires exactly one terminal_verdict artifact; "
            f"observed {terminal_count}"
        )
    return tuple(artifacts)


def validate_closure(project_root: Path, closure_path: Path) -> ClosureBundle:
    project = project_root.resolve()
    closure = resolve_project_path(project, str(closure_path), "closure matrix")
    rows = read_tsv_exact(closure, CLOSURE_FIELDS)
    counts = {workstream: 0 for workstream in REQUIRED_WORKSTREAMS}
    for row in rows:
        workstream = row["workstream_id"]
        if workstream not in counts:
            raise ReleaseContractError(
                f"unexpected workstream in closure: {workstream}"
            )
        counts[workstream] += 1
    if any(count != 1 for count in counts.values()) or len(rows) != len(counts):
        raise ReleaseContractError(
            "closure must contain exactly one row for Plans 13, 20, 30, 40, and 50; "
            f"observed {counts}"
        )
    closures = []
    global_snapshot_paths = set()
    for row in sorted(
        rows, key=lambda item: REQUIRED_WORKSTREAMS.index(item["workstream_id"])
    ):
        workstream = row["workstream_id"]
        state = row["terminal_state"]
        if state not in ALLOWED_TERMINAL_STATES:
            raise ReleaseContractError(
                f"unsupported terminal state for {workstream}: {state!r}"
            )
        for field in (
            "gate_id",
            "gate_verdict",
            "figure_role",
            "allowed_wording",
            "prohibited_wording",
            "adjudicator",
            "date",
        ):
            if not row[field].strip():
                raise ReleaseContractError(
                    f"{workstream} closure field is blank: {field}"
                )
        try:
            date.fromisoformat(row["date"])
        except ValueError as error:
            raise ReleaseContractError(
                f"{workstream} date is not ISO YYYY-MM-DD: {row['date']!r}"
            ) from error
        validate_terminal_inclusion(row)
        manifest = resolve_project_path(
            project,
            row["artifact_manifest"],
            f"{workstream} artifact_manifest",
        )
        require_within_any(
            manifest,
            allowed_workstream_roots(project, workstream),
            f"{workstream} artifact_manifest",
        )
        expected_manifest_hash = require_sha256(
            row["manifest_sha256"], f"{workstream} manifest_sha256"
        )
        if sha256_file(manifest) != expected_manifest_hash:
            raise ReleaseContractError(f"{workstream} terminal manifest hash mismatch")
        artifacts = validate_artifact_manifest(project, workstream, manifest)
        for artifact in artifacts:
            if artifact.snapshot_relpath in global_snapshot_paths:
                raise ReleaseContractError(
                    f"snapshot destination collision: {artifact.snapshot_relpath}"
                )
            global_snapshot_paths.add(artifact.snapshot_relpath)
        included = parse_bool(
            row["include_main"], f"{workstream} include_main"
        ) or parse_bool(row["include_supplement"], f"{workstream} include_supplement")
        if workstream == "PLAN20" and included:
            if row["gate_verdict"].upper() != "READY":
                raise ReleaseContractError(
                    "Plan 20 inclusion requires gate_verdict READY"
                )
            if (
                sum(artifact.artifact_role == "ready_seal" for artifact in artifacts)
                != 1
            ):
                raise ReleaseContractError(
                    "Plan 20 inclusion requires exactly one ready_seal artifact"
                )
        if workstream == "PLAN40" and not row["figure_role"].strip():
            raise ReleaseContractError(
                "Plan 40 requires the frozen figure-role verdict"
            )
        closures.append(
            WorkstreamClosure(
                row=dict(row),
                manifest_path_text=row["artifact_manifest"],
                manifest_path=manifest,
                manifest_sha256=expected_manifest_hash,
                artifacts=artifacts,
            )
        )
    return ClosureBundle(
        closure_path=closure,
        closure_sha256=sha256_file(closure),
        workstreams=tuple(closures),
    )


def read_protected_scopes(
    project_root: Path, scopes_path: Path
) -> list[dict[str, str]]:
    rows = read_tsv_exact(scopes_path, PROTECTED_SCOPE_FIELDS)
    if not rows:
        raise ReleaseContractError("protected scope registry is empty")
    ids = set()
    classes = set()
    roots = set()
    candidate = candidate_root(project_root)
    for row in rows:
        scope_id = row["scope_id"].strip()
        protection_class = row["protection_class"]
        if not scope_id or scope_id in ids:
            raise ReleaseContractError(
                f"blank or duplicate protected scope: {scope_id!r}"
            )
        ids.add(scope_id)
        if protection_class not in {"named_manuscript_release", "program_v1"}:
            raise ReleaseContractError(
                f"invalid protection class for {scope_id}: {protection_class}"
            )
        classes.add(protection_class)
        root = resolve_project_path(
            project_root, row["root_path"], f"protected scope {scope_id}"
        )
        if not root.is_dir() or root.is_symlink():
            raise ReleaseContractError(
                f"protected scope must be a non-symlink directory: {root}"
            )
        if root in roots:
            raise ReleaseContractError(f"duplicate protected root: {root}")
        roots.add(root)
        if is_relative_to(root, candidate) or is_relative_to(candidate, root):
            raise ReleaseContractError(
                f"protected scope overlaps candidate root: {root}"
            )
    if classes != {"named_manuscript_release", "program_v1"}:
        raise ReleaseContractError(
            "protected scopes require both named_manuscript_release and program_v1"
        )
    release_parent = project_root.resolve() / "RNA-seq/results/manuscript_release"
    discovered_named = {
        path.resolve()
        for path in release_parent.iterdir()
        if path.is_dir() and path.name != "candidates" and not path.is_symlink()
    }
    registered_named = {
        resolve_project_path(project_root, row["root_path"], row["scope_id"])
        for row in rows
        if row["protection_class"] == "named_manuscript_release"
    }
    if discovered_named != registered_named:
        raise ReleaseContractError(
            "protected scope registry does not cover every pre-existing named "
            f"manuscript release: discovered={sorted(map(str, discovered_named))}, "
            f"registered={sorted(map(str, registered_named))}"
        )
    return rows


def discover_scope_files(root: Path) -> list[Path]:
    files = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ReleaseContractError(f"protected scope contains a symlink: {path}")
        if path.is_file():
            files.append(path)
    if not files:
        raise ReleaseContractError(f"protected scope contains no files: {root}")
    return files


def capture_protected_baseline(
    project_root: Path,
    scopes_path: Path,
    output_path: Path,
) -> list[dict[str, object]]:
    scopes = read_protected_scopes(project_root, scopes_path)
    rows = []
    for scope in scopes:
        root = resolve_project_path(
            project_root, scope["root_path"], f"protected scope {scope['scope_id']}"
        )
        for path in discover_scope_files(root):
            rows.append(
                {
                    "scope_id": scope["scope_id"],
                    "relative_path": path.relative_to(root).as_posix(),
                    "sha256": sha256_file(path),
                    "bytes": path.stat().st_size,
                }
            )
    atomic_write_tsv(output_path, rows, PROTECTED_BASELINE_FIELDS)
    return rows


def verify_protected_baseline(
    project_root: Path,
    scopes_path: Path,
    baseline_path: Path,
) -> dict[str, object]:
    scopes = read_protected_scopes(project_root, scopes_path)
    baseline = read_tsv_exact(baseline_path, PROTECTED_BASELINE_FIELDS)
    by_scope: dict[str, dict[str, dict[str, str]]] = {
        row["scope_id"]: {} for row in scopes
    }
    for row in baseline:
        scope_id = row["scope_id"]
        if scope_id not in by_scope:
            raise ReleaseContractError(
                f"baseline references unknown protected scope: {scope_id}"
            )
        relative = clean_relative_path(
            row["relative_path"], f"protected {scope_id} relative_path"
        )
        if relative in by_scope[scope_id]:
            raise ReleaseContractError(
                f"duplicate protected path in {scope_id}: {relative}"
            )
        require_sha256(row["sha256"], f"protected {scope_id} {relative} sha256")
        parse_nonnegative_int(row["bytes"], f"protected {scope_id} {relative} bytes")
        by_scope[scope_id][relative] = row
    checked_files = 0
    for scope in scopes:
        scope_id = scope["scope_id"]
        root = resolve_project_path(
            project_root, scope["root_path"], f"protected scope {scope_id}"
        )
        discovered = {
            path.relative_to(root).as_posix(): path
            for path in discover_scope_files(root)
        }
        if set(discovered) != set(by_scope[scope_id]):
            raise ReleaseContractError(
                f"protected scope membership changed for {scope_id}: "
                f"baseline_only={sorted(set(by_scope[scope_id]) - set(discovered))}, "
                f"current_only={sorted(set(discovered) - set(by_scope[scope_id]))}"
            )
        for relative, path in discovered.items():
            row = by_scope[scope_id][relative]
            if (
                path.stat().st_size != int(row["bytes"])
                or sha256_file(path) != row["sha256"]
            ):
                raise ReleaseContractError(
                    f"protected release drift: {scope_id}/{relative}"
                )
            checked_files += 1
    return {
        "n_scopes": len(scopes),
        "n_files": checked_files,
        "scope_registry_sha256": sha256_file(scopes_path),
        "baseline_sha256": sha256_file(baseline_path),
    }


def repository_state(project_root: Path, fixture_mode: bool) -> dict[str, str]:
    if fixture_mode:
        return {
            "repository_commit": "synthetic-fixture",
            "tracked_diff_sha256": sha256_bytes(b""),
        }
    git_dir = project_root.resolve() / ".git"
    if not git_dir.exists():
        raise ReleaseContractError("real candidate snapshot requires a Git repository")
    commit = subprocess.run(
        ["git", "-C", str(project_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
    ).stdout.strip()
    diff = subprocess.run(
        ["git", "-C", str(project_root), "diff", "--binary", "HEAD", "--"],
        check=True,
        capture_output=True,
    ).stdout
    return {
        "repository_commit": commit.decode("ascii"),
        "tracked_diff_sha256": sha256_bytes(diff),
    }


def assert_no_symlinks(root: Path) -> None:
    if root.is_symlink():
        raise ReleaseContractError(f"candidate root is a symlink: {root}")
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ReleaseContractError(f"candidate contains a symlink: {path}")
