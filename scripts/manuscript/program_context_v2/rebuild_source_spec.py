#!/usr/bin/env python3
"""Deterministic scientific source identity for retained clean rebuilds.

The ordinary REL01 snapshot specification intentionally records the authority
workspace's global tracked diff and the runtime environment.  Those fields are
valuable provenance, but unrelated workspace dirtiness and execution context
must not define whether two scientific rebuilds used the same frozen sources.

This module defines the narrower identity used by retained A/B rebuilds:
repository commit, signed closure/artifacts, BASE selection, protected-source
contracts, and the exact current release-producer bytes frozen by REL01.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping, Sequence

from release_common import (
    CANDIDATE_ID,
    ReleaseContractError,
    canonical_json_bytes,
    parse_nonnegative_int,
    read_tsv_exact,
    require_sha256,
    sha256_bytes,
)
from release_products import BASE_SELECTION_FIELDS


REBUILD_SOURCE_CONTRACT = "plan60_rebuild_source_v1"
REBUILD_SOURCE_FILENAME = "manifests/rebuild_source_spec.json"
PRODUCER_FIELDS = (
    "producer_id",
    "repository_path",
    "sha256",
    "bytes",
    "snapshot_path",
)


def _normalized_rows(
    rows: Sequence[Mapping[str, object]],
    fields: Sequence[str],
    label: str,
) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    for index, row in enumerate(rows):
        if not all(field in row for field in fields):
            raise ReleaseContractError(f"{label} row {index} lacks required fields")
        normalized.append({field: row[field] for field in fields})
    return normalized


def build_rebuild_source_spec(
    *,
    repository_commit: str,
    closure_sha256: str,
    workstreams: Sequence[Mapping[str, object]],
    artifacts: Sequence[Mapping[str, object]],
    base_selection_sha256: str,
    base_selection_rows: Sequence[Mapping[str, object]],
    producers: Sequence[Mapping[str, object]],
    protected_scope_registry_sha256: str,
    protected_baseline_sha256: str,
    fixture_mode: bool,
) -> dict[str, object]:
    """Build the path-portable source identity shared by authority/A/B."""

    if not repository_commit.strip():
        raise ReleaseContractError("rebuild source requires a repository commit")
    require_sha256(closure_sha256, "rebuild closure SHA256")
    require_sha256(base_selection_sha256, "rebuild BASE selection SHA256")
    require_sha256(
        protected_scope_registry_sha256,
        "rebuild protected-scope registry SHA256",
    )
    require_sha256(protected_baseline_sha256, "rebuild protected baseline SHA256")
    normalized_workstreams = sorted(
        _normalized_rows(
            workstreams,
            (
                "workstream_id",
                "terminal_state",
                "gate_id",
                "gate_verdict",
                "manifest_sha256",
                "include_main",
                "include_supplement",
                "figure_role",
                "manifest_path_provenance",
            ),
            "rebuild workstream",
        ),
        key=lambda row: str(row["workstream_id"]),
    )
    normalized_artifacts = sorted(
        _normalized_rows(
            artifacts,
            (
                "workstream_id",
                "artifact_id",
                "artifact_role",
                "source_sha256",
                "source_bytes",
                "snapshot_relpath",
            ),
            "rebuild artifact",
        ),
        key=lambda row: (
            str(row["workstream_id"]),
            str(row["artifact_id"]),
        ),
    )
    normalized_selection = sorted(
        _normalized_rows(
            base_selection_rows,
            BASE_SELECTION_FIELDS,
            "rebuild BASE selection",
        ),
        key=lambda row: str(row["artifact_id"]),
    )
    normalized_producers = sorted(
        _normalized_rows(producers, PRODUCER_FIELDS, "rebuild producer"),
        key=lambda row: str(row["repository_path"]),
    )
    if (
        not normalized_workstreams
        or not normalized_artifacts
        or not normalized_selection
    ):
        raise ReleaseContractError(
            "rebuild source contract cannot contain empty inputs"
        )
    producer_paths: set[str] = set()
    producer_ids: set[str] = set()
    for row in normalized_producers:
        producer_id = str(row["producer_id"])
        repository_path = str(row["repository_path"])
        if (
            not producer_id
            or producer_id in producer_ids
            or not repository_path
            or repository_path in producer_paths
        ):
            raise ReleaseContractError("rebuild producer IDs/paths must be unique")
        producer_ids.add(producer_id)
        producer_paths.add(repository_path)
        require_sha256(str(row["sha256"]), f"rebuild producer {producer_id}")
        parse_nonnegative_int(row["bytes"], f"rebuild producer {producer_id} bytes")
    for row in normalized_artifacts:
        require_sha256(
            str(row["source_sha256"]),
            f"rebuild artifact {row['workstream_id']}:{row['artifact_id']}",
        )
        parse_nonnegative_int(
            row["source_bytes"],
            f"rebuild artifact {row['workstream_id']}:{row['artifact_id']} bytes",
        )
    selection_ids: set[str] = set()
    for row in normalized_selection:
        artifact_id = str(row["artifact_id"])
        if not artifact_id or artifact_id in selection_ids:
            raise ReleaseContractError(
                "rebuild BASE selection artifact IDs must be unique"
            )
        selection_ids.add(artifact_id)
        require_sha256(
            str(row["source_sha256"]),
            f"rebuild BASE source {artifact_id}",
        )
        parse_nonnegative_int(
            row["source_bytes"], f"rebuild BASE source {artifact_id} bytes"
        )
    spec = {
        "contract_version": REBUILD_SOURCE_CONTRACT,
        "candidate_id": CANDIDATE_ID,
        "repository_commit": repository_commit,
        "closure_sha256": closure_sha256,
        "workstreams": normalized_workstreams,
        "artifacts": normalized_artifacts,
        "base_selection_sha256": base_selection_sha256,
        "base_selection": normalized_selection,
        "producer_scripts": normalized_producers,
        "protected_scope_registry_sha256": protected_scope_registry_sha256,
        "protected_baseline_sha256": protected_baseline_sha256,
        "fixture_mode": fixture_mode,
    }
    return {
        "rebuild_source_spec": spec,
        "rebuild_source_spec_sha256": sha256_bytes(canonical_json_bytes(spec)),
    }


def validate_rebuild_source_payload(
    payload: Mapping[str, object],
    *,
    fixture_mode: bool,
) -> tuple[dict[str, object], str]:
    spec = payload.get("rebuild_source_spec")
    digest = payload.get("rebuild_source_spec_sha256")
    if not isinstance(spec, dict) or not isinstance(digest, str):
        raise ReleaseContractError(
            "rebuild source payload lacks rebuild_source_spec/hash"
        )
    require_sha256(digest, "rebuild source spec SHA256")
    if sha256_bytes(canonical_json_bytes(spec)) != digest:
        raise ReleaseContractError("rebuild source spec hash does not reproduce")
    if (
        spec.get("contract_version") != REBUILD_SOURCE_CONTRACT
        or spec.get("candidate_id") != CANDIDATE_ID
        or spec.get("fixture_mode") is not fixture_mode
    ):
        raise ReleaseContractError("rebuild source identity contract drift")
    rebuilt = build_rebuild_source_spec(
        repository_commit=str(spec.get("repository_commit") or ""),
        closure_sha256=str(spec.get("closure_sha256") or ""),
        workstreams=spec.get("workstreams")
        if isinstance(spec.get("workstreams"), list)
        else [],
        artifacts=spec.get("artifacts")
        if isinstance(spec.get("artifacts"), list)
        else [],
        base_selection_sha256=str(spec.get("base_selection_sha256") or ""),
        base_selection_rows=(
            spec.get("base_selection")
            if isinstance(spec.get("base_selection"), list)
            else []
        ),
        producers=(
            spec.get("producer_scripts")
            if isinstance(spec.get("producer_scripts"), list)
            else []
        ),
        protected_scope_registry_sha256=str(
            spec.get("protected_scope_registry_sha256") or ""
        ),
        protected_baseline_sha256=str(spec.get("protected_baseline_sha256") or ""),
        fixture_mode=fixture_mode,
    )
    if rebuilt != dict(payload):
        raise ReleaseContractError("rebuild source payload is not canonical")
    return spec, digest


def load_rebuild_source_spec(
    candidate: Path,
    *,
    fixture_mode: bool,
) -> tuple[dict[str, object], str]:
    path = candidate / REBUILD_SOURCE_FILENAME
    if not path.is_file() or path.is_symlink():
        raise ReleaseContractError(f"rebuild source spec is missing or unsafe: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ReleaseContractError(f"invalid rebuild source JSON: {path}") from error
    if not isinstance(payload, dict):
        raise ReleaseContractError("rebuild source payload must be an object")
    return validate_rebuild_source_payload(payload, fixture_mode=fixture_mode)


def read_base_selection(path: Path) -> tuple[list[dict[str, str]], str]:
    if not path.is_file() or path.is_symlink():
        raise ReleaseContractError(f"BASE selection is missing or unsafe: {path}")
    from release_common import sha256_file

    return read_tsv_exact(path, BASE_SELECTION_FIELDS), sha256_file(path)
