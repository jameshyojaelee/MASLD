#!/usr/bin/env python3
"""Stage coordinator-selected base hero inputs as an immutable snapshot extension."""

from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

from fibrosis_candidate_contract import (
    FIBROSIS_VALIDATION_ROLE,
    validate_fibrosis_selection_rows,
    validate_frozen_fibrosis_bundle,
)
from rel01_snapshot_candidate import (
    copy_verified,
    materialize_immutable_json,
    materialize_immutable_tsv,
)
from rel05_validate_candidate import validate_snapshot_base
from release_common import (
    CANDIDATE_ID,
    ReleaseContractError,
    assert_no_symlinks,
    clean_relative_path,
    is_relative_to,
    parse_nonnegative_int,
    project_relative,
    read_tsv_exact,
    require_sha256,
    resolve_project_path,
    sha256_file,
    verify_protected_baseline,
)
from release_products import (
    BASE_DOWNSTREAM_FIELDS,
    BASE_SELECTION_FIELDS,
    BASE_SNAPSHOT_FIELDS,
)


REQUIRED_BASE_ROLES = {
    "fibrosis_transition_raw",
    FIBROSIS_VALIDATION_ROLE,
    "cohort_overview",
}


def validate_selection(
    project_root: Path,
    candidate: Path,
    selection_path: Path,
    expected_sha256: str,
    fixture_mode: bool,
) -> tuple[list[dict[str, str]], str]:
    require_sha256(expected_sha256, "base input selection SHA256")
    selection = resolve_project_path(
        project_root, str(selection_path), "base input selection"
    )
    if is_relative_to(selection, candidate):
        raise ReleaseContractError(
            "base input selection must be coordinator-owned before candidate copy"
        )
    observed_hash = sha256_file(selection)
    if observed_hash != expected_sha256:
        raise ReleaseContractError("base input selection SHA256 mismatch")
    rows = read_tsv_exact(selection, BASE_SELECTION_FIELDS)
    if not rows:
        raise ReleaseContractError("base input selection is empty")
    selection_ids = {row["selection_id"].strip() for row in rows}
    coordinators = {row["coordinator"].strip() for row in rows}
    dates = {row["date"].strip() for row in rows}
    if len(selection_ids) != 1 or "" in selection_ids:
        raise ReleaseContractError("base selection requires one nonblank selection_id")
    if len(coordinators) != 1 or "" in coordinators:
        raise ReleaseContractError("base selection requires one named coordinator")
    if len(dates) != 1:
        raise ReleaseContractError("base selection requires one signing date")
    try:
        date.fromisoformat(next(iter(dates)))
    except ValueError as error:
        raise ReleaseContractError(
            "base selection date must be ISO YYYY-MM-DD"
        ) from error
    artifact_ids = set()
    destinations = set()
    roles = set()
    for row in rows:
        artifact_id = row["artifact_id"].strip()
        if not artifact_id or artifact_id in artifact_ids:
            raise ReleaseContractError(
                f"blank or duplicate base artifact_id: {artifact_id!r}"
            )
        artifact_ids.add(artifact_id)
        roles.add(row["artifact_role"])
        destination = clean_relative_path(
            row["snapshot_relpath"], f"base {artifact_id} snapshot_relpath"
        )
        if not destination.startswith("inputs/BASE/"):
            raise ReleaseContractError(
                f"base input must snapshot beneath inputs/BASE/: {artifact_id}"
            )
        if destination in destinations:
            raise ReleaseContractError(
                f"duplicate base snapshot destination: {destination}"
            )
        destinations.add(destination)
        require_sha256(row["source_sha256"], f"base {artifact_id} source SHA256")
        parse_nonnegative_int(row["source_bytes"], f"base {artifact_id} source bytes")
        if not row["allowed_wording"].strip() or not row["prohibited_wording"].strip():
            raise ReleaseContractError(
                f"base artifact lacks wording boundary: {artifact_id}"
            )
    missing_roles = REQUIRED_BASE_ROLES - roles
    if missing_roles:
        raise ReleaseContractError(
            f"base input selection lacks required roles: {sorted(missing_roles)}"
        )
    validate_fibrosis_selection_rows(
        project_root,
        rows,
        fixture_mode=fixture_mode,
        verify_external_sources=False,
    )
    return rows, observed_hash


def stage_base_inputs(
    project_root: Path,
    selection_path: Path,
    selection_sha256: str,
    fixture_mode: bool = False,
) -> dict[str, object]:
    base = validate_snapshot_base(project_root, fixture_mode=fixture_mode)
    project = base["project"]
    candidate = base["candidate"]
    if not isinstance(project, Path) or not isinstance(candidate, Path):
        raise ReleaseContractError("internal base-transition path type drift")
    state_path = candidate / "manifests/candidate_state.json"
    state_hash_before = sha256_file(state_path)
    protected_before = base["protected"]
    rows, selection_hash = validate_selection(
        project, candidate, selection_path, selection_sha256, fixture_mode
    )
    selection = resolve_project_path(
        project, str(selection_path), "base input selection"
    )
    prepared_inputs = []
    for row in sorted(rows, key=lambda item: item["snapshot_relpath"]):
        artifact_id = row["artifact_id"]
        source = resolve_project_path(
            project, row["source_path"], f"base {artifact_id} source"
        )
        if is_relative_to(source, candidate):
            raise ReleaseContractError(
                f"base source cannot point into candidate: {artifact_id}"
            )
        expected_hash = row["source_sha256"]
        expected_bytes = int(row["source_bytes"])
        if (
            not source.is_file()
            or source.is_symlink()
            or source.stat().st_size != expected_bytes
            or sha256_file(source) != expected_hash
        ):
            raise ReleaseContractError(
                f"base source failed pre-write hash/byte verification: {artifact_id}"
            )
        destination = candidate / row["snapshot_relpath"]
        if destination.exists() or destination.is_symlink():
            if (
                destination.is_symlink()
                or not destination.is_file()
                or destination.stat().st_size != expected_bytes
                or sha256_file(destination) != expected_hash
            ):
                raise ReleaseContractError(
                    f"existing immutable base destination differs: {artifact_id}"
                )
        prepared_inputs.append(
            (row, source, destination, expected_hash, expected_bytes)
        )

    copied_selection = candidate / "manifests/base_input_selection.tsv"
    copy_verified(
        selection,
        copied_selection,
        selection_hash,
        selection.stat().st_size,
        candidate,
    )

    snapshot_rows = []
    downstream_rows = []
    for row, source, destination, expected_hash, expected_bytes in prepared_inputs:
        artifact_id = row["artifact_id"]
        copy_verified(
            source,
            destination,
            expected_hash,
            expected_bytes,
            candidate,
        )
        snapshot_rows.append(
            {
                "artifact_id": artifact_id,
                "artifact_role": row["artifact_role"],
                "source_path_provenance": project_relative(project, source),
                "source_sha256": expected_hash,
                "source_bytes": expected_bytes,
                "snapshot_path": row["snapshot_relpath"],
                "snapshot_sha256": sha256_file(destination),
                "snapshot_bytes": destination.stat().st_size,
                "allowed_wording": row["allowed_wording"],
                "prohibited_wording": row["prohibited_wording"],
            }
        )
        downstream_rows.append(
            {
                "consumer_id": f"BASE:{artifact_id}",
                "artifact_id": artifact_id,
                "snapshot_path": row["snapshot_relpath"],
                "sha256": expected_hash,
                "bytes": expected_bytes,
            }
        )

    snapshot_manifest = candidate / "manifests/base_input_snapshot_manifest.tsv"
    downstream_manifest = candidate / "manifests/base_downstream_inputs.tsv"
    materialize_immutable_tsv(
        snapshot_manifest,
        snapshot_rows,
        BASE_SNAPSHOT_FIELDS,
        candidate,
        "base input snapshot manifest",
    )
    materialize_immutable_tsv(
        downstream_manifest,
        downstream_rows,
        BASE_DOWNSTREAM_FIELDS,
        candidate,
        "base downstream input manifest",
    )
    transition = {
        "contract_version": 1,
        "candidate_id": CANDIDATE_ID,
        "transition": "BASE_INPUTS_FROZEN",
        "from_candidate_state_sha256": state_hash_before,
        "input_snapshot_spec_sha256": base["spec_hash"],
        "selection_sha256": selection_hash,
        "selection_snapshot_sha256": sha256_file(copied_selection),
        "base_snapshot_manifest_sha256": sha256_file(snapshot_manifest),
        "base_downstream_manifest_sha256": sha256_file(downstream_manifest),
        "n_base_inputs": len(snapshot_rows),
        "fixture_mode": fixture_mode,
        "scientific_assembly_performed": False,
        "canonical_promotion_status": "not_promoted",
    }
    transition_path = candidate / "manifests/base_input_transition.json"
    materialize_immutable_json(
        transition_path,
        transition,
        candidate,
        "base input transition",
    )
    if sha256_file(state_path) != state_hash_before:
        raise ReleaseContractError(
            "base transition altered the immutable candidate state"
        )
    protected_after = verify_protected_baseline(
        project,
        candidate / "manifests/protected_scopes.tsv",
        candidate / "manifests/protected_release_baseline.tsv",
    )
    if protected_before != protected_after:
        raise ReleaseContractError(
            "protected release verification changed during base-input transition"
        )
    observed_base = {
        path.relative_to(candidate).as_posix()
        for path in (candidate / "inputs/BASE").rglob("*")
        if path.is_file()
    }
    expected_base = {row["snapshot_path"] for row in snapshot_rows}
    if observed_base != expected_base:
        raise ReleaseContractError(
            "base input directory does not exactly match signed selection"
        )
    validate_frozen_fibrosis_bundle(candidate, fixture_mode=fixture_mode)
    assert_no_symlinks(candidate)
    return {
        "status": "REL02_BASE_INPUTS_FROZEN",
        "candidate_id": CANDIDATE_ID,
        "selection_sha256": selection_hash,
        "transition_sha256": sha256_file(transition_path),
        "n_base_inputs": len(snapshot_rows),
        "real_candidate_assembled": not fixture_mode,
        "canonical_promotion_status": "not_promoted",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--base-input-selection", type=Path, required=True)
    parser.add_argument("--selection-sha256", required=True)
    parser.add_argument("--fixture-mode", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        result = stage_base_inputs(
            args.project_root,
            args.base_input_selection,
            args.selection_sha256,
            fixture_mode=args.fixture_mode,
        )
    except ReleaseContractError as error:
        raise SystemExit(f"REL02_BASE_BLOCKED: {error}") from error
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
