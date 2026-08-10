#!/usr/bin/env python3
"""Validate immutable snapshot infrastructure without promoting or assembling it."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from datetime import date
from pathlib import Path

from release_common import (
    ALLOWED_TERMINAL_STATES,
    ARTIFACT_FIELDS,
    CANDIDATE_ID,
    CANDIDATE_REL,
    CLOSURE_FIELDS,
    DOWNSTREAM_FIELDS,
    PLAN50_MANIFEST_PRODUCER_PATHS,
    REQUIRED_WORKSTREAMS,
    SCIENTIFIC_REPORT_PRODUCER_PATHS,
    REQUIRED_CANDIDATE_DIRS,
    REQUIRED_SCIENTIFIC_CHECKS,
    SCIENTIFIC_REGISTRY_FIELDS,
    SNAPSHOT_MANIFEST_FIELDS,
    WORKSTREAM_LANES,
    Artifact,
    ClosureBundle,
    WorkstreamClosure,
    ReleaseContractError,
    assert_candidate_root,
    assert_no_symlinks,
    atomic_write_json,
    candidate_root,
    canonical_json_bytes,
    clean_relative_path,
    is_relative_to,
    parse_nonnegative_int,
    parse_bool,
    read_tsv_exact,
    require_sha256,
    resolve_project_path,
    sha256_bytes,
    sha256_file,
    validate_terminal_inclusion,
    verify_protected_baseline,
)
from rel01_snapshot_candidate import (
    ENVIRONMENT_INDEX_FIELDS,
    MANIFEST_REGISTRY_FIELDS,
    PLAN50_RELOCATED_VERIFICATION_REL,
    PLAN50_TERMINAL_PROVENANCE_FIELDS,
    PYTHON_PACKAGE_FIELDS,
    PRODUCER_SCRIPT_FIELDS,
    verify_relocated_plan50_bundle,
)
from publish_plan60_workstream_handoff import PLAN50_PAYLOAD_MANIFEST_FIELDS
from adapter_contract import (
    ADAPTER_CONTRACT,
    ADAPTER_CONTRACT_SHA256,
    ADAPTER_PROVENANCE_FIELDS,
    ADAPTER_PROVENANCE_VERSION,
    EXPECTED_ADAPTER_IDS,
    RECURSIVE_PRODUCER_FIELDS,
    require_exact_adapter_binding_paths,
)
from coordinator.compare_clean_rebuilds import revalidate_clean_rebuild_report
from fibrosis_candidate_contract import (
    TRUE_KLEINER_TRANSITIONS,
    validate_frozen_fibrosis_bundle,
)
from rebuild_source_spec import load_rebuild_source_spec
from release_products import (
    ANALYSIS_MANIFEST_FIELDS,
    BASE_SNAPSHOT_FIELDS,
    CLAIM_FIELDS,
    EXPECTED_FIGURES,
    EXPECTED_SUPPLEMENTARY_FIGURES,
    EXPECTED_SUPPLEMENTARY_PANELS,
    EXPECTED_SUPPLEMENTARY_SECTIONS,
    FIGURE_SOURCE_MANIFEST_FIELDS,
    JOURNAL_BRANCH,
    MANUSCRIPT_MANIFEST_FIELDS,
    MYOJIN_ROLE,
    NUMBERS_FIELDS,
    SOURCE_DEPENDENCY_FIELDS,
    SOURCE_ROW_FIELDS,
    TRANSITION_PRODUCT_FIELDS,
    candidate_figure_root,
    candidate_manuscript_root,
    contains_nmf_discrete_class_claim,
    validate_nmf_source_language,
    load_frozen_input_index,
    read_json_object,
    validate_source_dependency_contract,
)


ALLOWED_CANDIDATE_PATH_PREFIXES = {
    "inputs",
    "tables",
    "figure_sources",
    "manifests",
    "claims",
    "environments",
    "logs",
}
MUTABLE_PROJECT_PREFIXES = (
    "Analysis/",
    "RNA-seq/",
    "figures/",
    "docs/",
    "scripts/",
    "data/",
    "GWAS/",
    "archive/",
    "config/",
    "streamlit_deg_explorer/",
    "Cas13_Library_Design/",
)
SCIENTIFIC_REPORT_CONTRACT = "plan60_real_scientific_validation_v2"


def verify_frozen_adapter_provenance(
    candidate: Path,
    state: dict[str, object],
    workstream_snapshot_rows: list[dict[str, str]],
) -> dict[str, object]:
    provenance_path = candidate_file(
        candidate,
        "inputs/BASE/adapter_provenance.tsv",
        "adapter provenance",
    )
    recursive_path = candidate_file(
        candidate,
        "inputs/BASE/recursive_release_producers.tsv",
        "recursive producer manifest",
    )
    rows = read_tsv_exact(provenance_path, ADAPTER_PROVENANCE_FIELDS)
    recursive_rows = read_tsv_exact(recursive_path, RECURSIVE_PRODUCER_FIELDS)
    frozen_rows = read_tsv_exact(
        candidate / "manifests/producer_script_manifest.tsv",
        PRODUCER_SCRIPT_FIELDS,
    )
    if (
        len(rows) != len(EXPECTED_ADAPTER_IDS)
        or {row["artifact_id"] for row in rows} != EXPECTED_ADAPTER_IDS
    ):
        raise ReleaseContractError(
            "frozen adapter provenance does not contain exactly seven adapters"
        )
    recursive_by_path = {row["repository_path"]: row for row in recursive_rows}
    frozen_by_path = {row["repository_path"]: row for row in frozen_rows}
    if len(recursive_by_path) != len(recursive_rows) or len(frozen_by_path) != len(
        frozen_rows
    ):
        raise ReleaseContractError("duplicate producer path in frozen manifests")
    recursive_sha = sha256_file(recursive_path)
    base_snapshot_rows = read_tsv_exact(
        candidate / "manifests/base_input_snapshot_manifest.tsv",
        BASE_SNAPSHOT_FIELDS,
    )
    frozen_sources: dict[tuple[str, str, int], list[Path]] = defaultdict(list)
    for snapshot_row in [*workstream_snapshot_rows, *base_snapshot_rows]:
        source_path = snapshot_row["source_path_provenance"]
        source_sha = require_sha256(
            snapshot_row["source_sha256"], f"frozen source {source_path}"
        )
        source_bytes = parse_nonnegative_int(
            snapshot_row["source_bytes"], f"frozen source {source_path} bytes"
        )
        snapshot_field = (
            "snapshot_relpath"
            if "snapshot_relpath" in snapshot_row
            else "snapshot_path"
        )
        frozen_path = candidate_file(
            candidate,
            snapshot_row[snapshot_field],
            f"frozen source {source_path}",
        )
        if (
            sha256_file(frozen_path) != source_sha
            or frozen_path.stat().st_size != source_bytes
        ):
            raise ReleaseContractError(
                f"frozen adapter-input source drift: {source_path}"
            )
        frozen_sources[(source_path, source_sha, source_bytes)].append(frozen_path)
    for row in rows:
        artifact_id = row["artifact_id"]
        if (
            row["contract_version"] != ADAPTER_PROVENANCE_VERSION
            or row["adapter_contract_sha256"] != ADAPTER_CONTRACT_SHA256
            or row["canonical_promotion_authorized"].lower() != "false"
            or row["producer_manifest_relative_path"]
            != "inputs/BASE/recursive_release_producers.tsv"
            or row["producer_manifest_sha256"] != recursive_sha
        ):
            raise ReleaseContractError(
                f"frozen adapter provenance contract drift: {artifact_id}"
            )
        output_relative = Path(row["output_relative_path"])
        if (
            output_relative.is_absolute()
            or len(output_relative.parts) != 1
            or output_relative.as_posix()
            != ADAPTER_CONTRACT[artifact_id]["output_relative_path"]
        ):
            raise ReleaseContractError(
                f"unsafe frozen adapter output path: {artifact_id}"
            )
        output = candidate_file(
            candidate,
            f"inputs/BASE/{output_relative.as_posix()}",
            f"adapter output {artifact_id}",
        )
        if sha256_file(output) != require_sha256(
            row["output_sha256"], f"adapter output {artifact_id}"
        ) or output.stat().st_size != parse_nonnegative_int(
            row["output_bytes"], f"adapter output {artifact_id} bytes"
        ):
            raise ReleaseContractError(
                f"frozen adapter output hash/byte drift: {artifact_id}"
            )
        try:
            producer_bindings = json.loads(row["producer_bindings_json"])
            input_bindings = json.loads(row["input_bindings_json"])
        except json.JSONDecodeError as error:
            raise ReleaseContractError(
                f"invalid frozen adapter provenance JSON: {artifact_id}"
            ) from error
        if not isinstance(producer_bindings, list) or not producer_bindings:
            raise ReleaseContractError(
                f"frozen adapter lacks producer bindings: {artifact_id}"
            )
        if not isinstance(input_bindings, list) or not input_bindings:
            raise ReleaseContractError(
                f"frozen adapter lacks input bindings: {artifact_id}"
            )
        require_exact_adapter_binding_paths(
            artifact_id,
            (
                str(binding.get("source_path", ""))
                for binding in input_bindings
                if isinstance(binding, dict)
            ),
            (
                str(binding.get("repository_path", ""))
                for binding in producer_bindings
                if isinstance(binding, dict)
            ),
        )
        for binding in input_bindings:
            if not isinstance(binding, dict) or set(binding) != {
                "source_path",
                "sha256",
                "bytes",
            }:
                raise ReleaseContractError(
                    f"frozen adapter input schema drift: {artifact_id}"
                )
            source_path = str(binding["source_path"])
            source_sha = require_sha256(
                str(binding["sha256"]),
                f"frozen adapter input {artifact_id}:{source_path}",
            )
            source_bytes = parse_nonnegative_int(
                binding["bytes"],
                f"frozen adapter input {artifact_id}:{source_path} bytes",
            )
            matches = frozen_sources.get((source_path, source_sha, source_bytes), [])
            if not matches:
                raise ReleaseContractError(
                    "adapter input is not represented by an immutable frozen snapshot: "
                    f"{artifact_id}:{source_path}"
                )
        for binding in producer_bindings:
            if not isinstance(binding, dict) or set(binding) != {
                "producer_id",
                "repository_path",
                "sha256",
                "bytes",
            }:
                raise ReleaseContractError(
                    f"frozen adapter producer schema drift: {artifact_id}"
                )
            repository_path = str(binding["repository_path"])
            recursive = recursive_by_path.get(repository_path)
            frozen = frozen_by_path.get(repository_path)
            if recursive is None or frozen is None:
                raise ReleaseContractError(
                    f"frozen adapter producer is not in both manifests: {repository_path}"
                )
            binding_sha = require_sha256(
                str(binding["sha256"]), f"adapter producer {repository_path}"
            )
            binding_bytes = parse_nonnegative_int(
                binding["bytes"], f"adapter producer {repository_path} bytes"
            )
            if (
                str(binding["producer_id"]) != recursive["producer_id"]
                or binding_sha != recursive["sha256"]
                or binding_bytes
                != parse_nonnegative_int(
                    recursive["bytes"], f"recursive producer {repository_path} bytes"
                )
                or binding_sha != frozen["sha256"]
                or binding_bytes
                != parse_nonnegative_int(
                    frozen["bytes"], f"frozen producer {repository_path} bytes"
                )
            ):
                raise ReleaseContractError(
                    f"frozen adapter producer binding drift: {repository_path}"
                )
    observed_summary = {
        "contract_version": ADAPTER_PROVENANCE_VERSION,
        "adapter_contract_sha256": ADAPTER_CONTRACT_SHA256,
        "n_adapters": len(rows),
        "base_input_selection_sha256": sha256_file(
            candidate_file(
                candidate,
                "manifests/base_input_selection.tsv",
                "frozen BASE selection",
            )
        ),
        "adapter_provenance_sha256": sha256_file(provenance_path),
        "recursive_producer_manifest_sha256": recursive_sha,
    }
    if state.get("adapter_provenance") != observed_summary:
        raise ReleaseContractError("candidate-state adapter provenance summary drift")
    return observed_summary


def load_json(path: Path, context: str) -> dict[str, object]:
    if not path.is_file() or path.is_symlink():
        raise ReleaseContractError(f"{context} is missing or symlinked: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ReleaseContractError(f"{context} is invalid JSON: {path}") from error
    if not isinstance(payload, dict):
        raise ReleaseContractError(f"{context} must be a JSON object: {path}")
    return payload


def candidate_file(candidate: Path, relative_text: str, context: str) -> Path:
    relative = clean_relative_path(relative_text, context)
    if relative.split("/", 1)[0] not in ALLOWED_CANDIDATE_PATH_PREFIXES:
        raise ReleaseContractError(
            f"{context} does not use a candidate-owned prefix: {relative}"
        )
    path = candidate / relative
    try:
        path.resolve().relative_to(candidate.resolve())
    except ValueError as error:
        raise ReleaseContractError(
            f"{context} escapes candidate: {relative}"
        ) from error
    return path


def verify_plan50_relocated_verification(
    candidate: Path,
    state: dict[str, object],
    frozen_producers: list[dict[str, str]],
    fixture_mode: bool,
) -> dict[str, object] | None:
    """Re-derive the Plan 50 producer and relocated-verification bindings."""

    verification_path = candidate / PLAN50_RELOCATED_VERIFICATION_REL
    if fixture_mode:
        if verification_path.exists() or verification_path.is_symlink():
            raise ReleaseContractError(
                "fixture snapshot unexpectedly contains Plan 50 relocated verification"
            )
        if (
            state.get("plan50_producer_binding") is not None
            or state.get("plan50_relocated_verification_sha256") is not None
        ):
            raise ReleaseContractError(
                "fixture candidate state claims production Plan 50 verification"
            )
        return None

    record = load_json(verification_path, "relocated Plan 50 verification")
    if state.get("plan50_relocated_verification_sha256") != sha256_file(
        verification_path
    ):
        raise ReleaseContractError(
            "candidate-state relocated Plan 50 verification hash drift"
        )
    expected_keys = {
        "contract_version",
        "status",
        "bundle_path",
        "validator_path",
        "validator_sha256",
        "seal_sha256",
        "payload_manifest_sha256",
        "terminal_provenance_sha256",
        "verification_result_sha256",
        "artifact_count",
    }
    if set(record) != expected_keys:
        raise ReleaseContractError("relocated Plan 50 verification schema drift")
    if (
        record["contract_version"] != 1
        or record["status"] != "pass"
        or record["bundle_path"] != "inputs/PASS"
        or record["validator_path"]
        != "environments/producers/portal/validate_evidence_passports.py"
    ):
        raise ReleaseContractError("relocated Plan 50 verification semantics drift")

    bundle = candidate / "inputs/PASS"
    validator = candidate_file(
        candidate, str(record["validator_path"]), "frozen Plan 50 validator"
    )
    if (
        not bundle.is_dir()
        or bundle.is_symlink()
        or not validator.is_file()
        or validator.is_symlink()
    ):
        raise ReleaseContractError("relocated Plan 50 bundle/validator is unsafe")
    payload_manifest = bundle / "passport_release_manifest.tsv"
    terminal_provenance = bundle / "passport_terminal_provenance.tsv"
    seal = bundle / "PASS06_VALIDATED"
    observed_hashes = {
        "validator_sha256": sha256_file(validator),
        "seal_sha256": sha256_file(seal),
        "payload_manifest_sha256": sha256_file(payload_manifest),
        "terminal_provenance_sha256": sha256_file(terminal_provenance),
    }
    for field, observed in observed_hashes.items():
        if require_sha256(str(record[field]), field) != observed:
            raise ReleaseContractError(f"relocated Plan 50 verification {field} drift")
    require_sha256(
        str(record["verification_result_sha256"]),
        "relocated Plan 50 verification result SHA256",
    )
    artifact_count = sum(path.is_file() for path in bundle.rglob("*"))
    if (
        parse_nonnegative_int(
            record["artifact_count"], "relocated Plan 50 artifact count"
        )
        != artifact_count
    ):
        raise ReleaseContractError(
            "relocated Plan 50 verification artifact count drift"
        )

    producer_index = {row["repository_path"]: row for row in frozen_producers}
    if len(producer_index) != len(frozen_producers):
        raise ReleaseContractError("duplicate frozen producer repository path")
    payload_rows = read_tsv_exact(payload_manifest, PLAN50_PAYLOAD_MANIFEST_FIELDS)
    observed_producers = {row["producer"] for row in payload_rows}
    if observed_producers != set(PLAN50_MANIFEST_PRODUCER_PATHS):
        raise ReleaseContractError("relocated Plan 50 payload producer universe drift")

    def verify_producer(row: dict[str, str], context: str) -> None:
        producer = row["producer"]
        frozen = producer_index.get(producer)
        if frozen is None or require_sha256(
            row["producer_sha256"], f"{context} producer SHA256"
        ) != require_sha256(frozen["sha256"], f"frozen producer {producer} SHA256"):
            raise ReleaseContractError(
                f"{context} producer differs from the frozen copy: {producer}"
            )

    for row in payload_rows:
        verify_producer(row, f"relocated Plan 50 payload {row['relative_path']}")
    terminal_rows = read_tsv_exact(
        terminal_provenance, PLAN50_TERMINAL_PROVENANCE_FIELDS
    )
    if len(terminal_rows) != 2 or {row["artifact_role"] for row in terminal_rows} != {
        "payload_manifest",
        "plan60_handoff",
    }:
        raise ReleaseContractError("relocated Plan 50 terminal provenance role drift")
    for row in terminal_rows:
        verify_producer(
            row,
            f"relocated Plan 50 terminal provenance {row['artifact_role']}",
        )

    expected_binding = {
        "payload_manifest_sha256": sha256_file(payload_manifest),
        "terminal_provenance_sha256": sha256_file(terminal_provenance),
        "manifest_producer_count": len(observed_producers),
    }
    if state.get("plan50_producer_binding") != expected_binding:
        raise ReleaseContractError(
            "candidate-state Plan 50 frozen-producer binding drift"
        )
    replayed = verify_relocated_plan50_bundle(candidate)
    if replayed != record:
        raise ReleaseContractError(
            "REL05 frozen Plan 50 replay differs from the REL01 verification record"
        )
    return record


def validate_frozen_closure(candidate: Path, spec: dict[str, object]) -> ClosureBundle:
    """Validate closure semantics using only files copied into the snapshot."""

    closure_path = candidate / "manifests/workstream_closure.tsv"
    rows = read_tsv_exact(closure_path, CLOSURE_FIELDS)
    if sha256_file(closure_path) != spec.get("closure_sha256"):
        raise ReleaseContractError("frozen closure hash differs from snapshot spec")
    counts = {workstream: 0 for workstream in REQUIRED_WORKSTREAMS}
    for row in rows:
        workstream = row["workstream_id"]
        if workstream not in counts:
            raise ReleaseContractError(
                f"unexpected workstream in frozen closure: {workstream}"
            )
        counts[workstream] += 1
    if any(count != 1 for count in counts.values()) or len(rows) != len(counts):
        raise ReleaseContractError(
            f"frozen closure must contain each required workstream once: {counts}"
        )

    workstreams: list[WorkstreamClosure] = []
    global_snapshots: set[str] = set()
    for row in sorted(
        rows, key=lambda item: REQUIRED_WORKSTREAMS.index(item["workstream_id"])
    ):
        workstream = row["workstream_id"]
        if row["terminal_state"] not in ALLOWED_TERMINAL_STATES:
            raise ReleaseContractError(
                f"unsupported frozen terminal state for {workstream}: "
                f"{row['terminal_state']!r}"
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
                    f"frozen {workstream} closure field is blank: {field}"
                )
        try:
            date.fromisoformat(row["date"])
        except ValueError as error:
            raise ReleaseContractError(
                f"frozen {workstream} date is not ISO YYYY-MM-DD"
            ) from error
        validate_terminal_inclusion(row)
        manifest_path_text = clean_relative_path(
            row["artifact_manifest"],
            f"frozen {workstream} artifact-manifest provenance",
        )
        if (
            workstream == "PLAN13"
            and not bool(spec.get("fixture_mode"))
            and "spatial_context_semantic_v2_2026-08-08/" not in manifest_path_text
        ):
            raise ReleaseContractError(
                "Plan13 production closure must originate from the sealed "
                "semantic-v2 workstream; legacy spatial_context is non-consumable"
            )
        manifest_hash = require_sha256(
            row["manifest_sha256"], f"frozen {workstream} manifest SHA256"
        )
        copied_manifest = candidate / f"manifests/upstream/{workstream}.tsv"
        if (
            not copied_manifest.is_file()
            or copied_manifest.is_symlink()
            or sha256_file(copied_manifest) != manifest_hash
        ):
            raise ReleaseContractError(
                f"frozen upstream terminal manifest drift: {workstream}"
            )
        artifact_rows = read_tsv_exact(copied_manifest, ARTIFACT_FIELDS)
        if not artifact_rows:
            raise ReleaseContractError(
                f"frozen terminal artifact manifest is empty: {workstream}"
            )
        artifact_ids: set[str] = set()
        artifacts: list[Artifact] = []
        lane = WORKSTREAM_LANES[workstream]
        for artifact_row in artifact_rows:
            artifact_id = artifact_row["artifact_id"].strip()
            if (
                not artifact_id
                or artifact_id == "terminal_artifact_manifest"
                or artifact_id in artifact_ids
            ):
                raise ReleaseContractError(
                    f"blank, reserved, or duplicate frozen artifact ID: "
                    f"{workstream}:{artifact_id!r}"
                )
            artifact_ids.add(artifact_id)
            source_path_text = clean_relative_path(
                artifact_row["source_path"],
                f"frozen {workstream}:{artifact_id} source provenance",
            )
            if (
                workstream == "PLAN13"
                and not bool(spec.get("fixture_mode"))
                and "spatial_context_semantic_v2_2026-08-08/" not in source_path_text
            ):
                raise ReleaseContractError(
                    "Plan13 artifact bypasses semantic-v2 evidence-state remediation: "
                    f"{artifact_id}"
                )
            source_hash = require_sha256(
                artifact_row["source_sha256"],
                f"frozen {workstream}:{artifact_id} source SHA256",
            )
            source_bytes = parse_nonnegative_int(
                artifact_row["source_bytes"],
                f"frozen {workstream}:{artifact_id} source bytes",
            )
            snapshot = clean_relative_path(
                artifact_row["snapshot_relpath"],
                f"frozen {workstream}:{artifact_id} snapshot path",
            )
            if not snapshot.startswith(f"inputs/{lane}/"):
                raise ReleaseContractError(
                    f"frozen artifact escaped inputs/{lane}/: {workstream}:{artifact_id}"
                )
            if snapshot in global_snapshots:
                raise ReleaseContractError(
                    f"duplicate frozen snapshot destination: {snapshot}"
                )
            global_snapshots.add(snapshot)
            downstream = clean_relative_path(
                artifact_row["downstream_read_path"],
                f"frozen {workstream}:{artifact_id} downstream path",
            )
            if downstream != snapshot:
                raise ReleaseContractError(
                    f"frozen downstream path differs from snapshot: "
                    f"{workstream}:{artifact_id}"
                )
            role = artifact_row["artifact_role"].strip()
            if not role:
                raise ReleaseContractError(
                    f"blank frozen artifact role: {workstream}:{artifact_id}"
                )
            artifacts.append(
                Artifact(
                    workstream_id=workstream,
                    artifact_id=artifact_id,
                    source_path_text=source_path_text,
                    source_path=candidate / snapshot,
                    source_sha256=source_hash,
                    source_bytes=source_bytes,
                    snapshot_relpath=snapshot,
                    artifact_role=role,
                    downstream_read_path=downstream,
                )
            )
        if sum(item.artifact_role == "terminal_verdict" for item in artifacts) != 1:
            raise ReleaseContractError(
                f"frozen {workstream} requires exactly one terminal verdict"
            )
        included = parse_bool(
            row["include_main"], f"frozen {workstream} include_main"
        ) or parse_bool(
            row["include_supplement"], f"frozen {workstream} include_supplement"
        )
        if workstream == "PLAN20" and included:
            if (
                row["gate_verdict"].upper() != "READY"
                or sum(item.artifact_role == "ready_seal" for item in artifacts) != 1
            ):
                raise ReleaseContractError(
                    "frozen Plan20 inclusion requires READY and one ready seal"
                )
        workstreams.append(
            WorkstreamClosure(
                row=dict(row),
                manifest_path_text=manifest_path_text,
                manifest_path=copied_manifest,
                manifest_sha256=manifest_hash,
                artifacts=tuple(artifacts),
            )
        )
    return ClosureBundle(
        closure_path=closure_path,
        closure_sha256=sha256_file(closure_path),
        workstreams=tuple(workstreams),
    )


def verify_snapshot_manifest(
    candidate: Path,
) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    path = candidate / "manifests/input_snapshot_manifest.tsv"
    rows = read_tsv_exact(path, SNAPSHOT_MANIFEST_FIELDS)
    if not rows:
        raise ReleaseContractError("input snapshot manifest is empty")
    by_snapshot = {}
    for row in rows:
        relative = clean_relative_path(
            row["snapshot_relpath"], "snapshot manifest snapshot_relpath"
        )
        if relative in by_snapshot:
            raise ReleaseContractError(f"duplicate frozen snapshot path: {relative}")
        destination = candidate_file(candidate, relative, "frozen snapshot path")
        if not destination.is_file() or destination.is_symlink():
            raise ReleaseContractError(
                f"frozen snapshot is missing or symlinked: {destination}"
            )
        snapshot_hash = require_sha256(
            row["snapshot_sha256"], f"snapshot {relative} sha256"
        )
        snapshot_bytes = parse_nonnegative_int(
            row["snapshot_bytes"], f"snapshot {relative} bytes"
        )
        if (
            destination.stat().st_size != snapshot_bytes
            or sha256_file(destination) != snapshot_hash
        ):
            raise ReleaseContractError(f"frozen snapshot drift: {relative}")
        source_hash = require_sha256(
            row["source_sha256"], f"source provenance {relative} sha256"
        )
        source_bytes = parse_nonnegative_int(
            row["source_bytes"], f"source provenance {relative} bytes"
        )
        if source_hash != snapshot_hash or source_bytes != snapshot_bytes:
            raise ReleaseContractError(
                f"source/snapshot integrity fields differ: {relative}"
            )
        provenance = clean_relative_path(
            row["source_path_provenance"], f"source provenance {relative}"
        )
        candidate_marker = CANDIDATE_REL.as_posix()
        if provenance == candidate_marker or provenance.startswith(
            f"{candidate_marker}/"
        ):
            raise ReleaseContractError(
                f"source provenance cannot point back into the candidate: {relative}"
            )
        by_snapshot[relative] = row
    return rows, by_snapshot


def verify_spec_and_snapshot_contract(
    project_root: Path,
    candidate: Path,
    closure,
    spec: dict[str, object],
    snapshot_rows: list[dict[str, str]],
) -> None:
    if spec.get("contract_version") != 1:
        raise ReleaseContractError("unsupported or missing snapshot contract version")
    expected_workstreams = [
        {
            "workstream_id": workstream.row["workstream_id"],
            "terminal_state": workstream.row["terminal_state"],
            "gate_id": workstream.row["gate_id"],
            "gate_verdict": workstream.row["gate_verdict"],
            "manifest_sha256": workstream.manifest_sha256,
            "include_main": workstream.row["include_main"],
            "include_supplement": workstream.row["include_supplement"],
            "figure_role": workstream.row["figure_role"],
            "manifest_path_provenance": workstream.manifest_path_text,
        }
        for workstream in closure.workstreams
    ]
    if spec.get("workstreams") != expected_workstreams:
        raise ReleaseContractError(
            "snapshot specification workstreams differ from the signed closure"
        )
    expected_artifacts = [
        {
            "workstream_id": artifact.workstream_id,
            "artifact_id": artifact.artifact_id,
            "artifact_role": artifact.artifact_role,
            "source_sha256": artifact.source_sha256,
            "source_bytes": artifact.source_bytes,
            "snapshot_relpath": artifact.snapshot_relpath,
        }
        for artifact in sorted(
            closure.artifacts, key=lambda item: item.snapshot_relpath
        )
    ]
    if spec.get("artifacts") != expected_artifacts:
        raise ReleaseContractError(
            "snapshot specification artifacts differ from terminal manifests"
        )

    expected: dict[tuple[str, str], dict[str, object]] = {
        ("REL00", "workstream_closure"): {
            "artifact_role": "closure_matrix",
            "snapshot_relpath": "manifests/workstream_closure.tsv",
            "sha256": spec.get("closure_sha256"),
        },
        ("REL01", "protected_scopes"): {
            "artifact_role": "protected_release_registry",
            "snapshot_relpath": "manifests/protected_scopes.tsv",
            "sha256": spec.get("protected_scope_registry_sha256"),
        },
        ("REL01", "protected_release_baseline"): {
            "artifact_role": "protected_release_baseline",
            "snapshot_relpath": "manifests/protected_release_baseline.tsv",
            "sha256": spec.get("protected_baseline_sha256"),
        },
    }
    for workstream in closure.workstreams:
        key = (workstream.row["workstream_id"], "terminal_artifact_manifest")
        expected[key] = {
            "artifact_role": "upstream_terminal_manifest",
            "snapshot_relpath": (
                f"manifests/upstream/{workstream.row['workstream_id']}.tsv"
            ),
            "sha256": workstream.manifest_sha256,
            "source_path_provenance": workstream.manifest_path_text,
        }
        for artifact in workstream.artifacts:
            expected[(artifact.workstream_id, artifact.artifact_id)] = {
                "artifact_role": artifact.artifact_role,
                "snapshot_relpath": artifact.snapshot_relpath,
                "sha256": artifact.source_sha256,
                "bytes": artifact.source_bytes,
                "source_path_provenance": artifact.source_path_text,
            }

    observed: dict[tuple[str, str], dict[str, str]] = {}
    for row in snapshot_rows:
        key = (row["workstream_id"], row["artifact_id"])
        if key in observed:
            raise ReleaseContractError(
                f"duplicate snapshot contract identity: {key[0]}:{key[1]}"
            )
        observed[key] = row
    if set(observed) != set(expected):
        raise ReleaseContractError(
            "input snapshot manifest does not exactly match the closure contract: "
            f"missing={sorted(set(expected) - set(observed))}, "
            f"extra={sorted(set(observed) - set(expected))}"
        )
    for key, contract in expected.items():
        row = observed[key]
        if row["artifact_role"] != contract["artifact_role"]:
            raise ReleaseContractError(f"snapshot role drift: {key[0]}:{key[1]}")
        if row["snapshot_relpath"] != contract["snapshot_relpath"]:
            raise ReleaseContractError(f"snapshot path drift: {key[0]}:{key[1]}")
        if row["source_sha256"] != contract["sha256"]:
            raise ReleaseContractError(f"snapshot source hash drift: {key[0]}:{key[1]}")
        if "bytes" in contract and int(row["source_bytes"]) != contract["bytes"]:
            raise ReleaseContractError(f"snapshot source byte drift: {key[0]}:{key[1]}")
        if (
            "source_path_provenance" in contract
            and row["source_path_provenance"] != contract["source_path_provenance"]
        ):
            raise ReleaseContractError(
                f"snapshot source provenance drift: {key[0]}:{key[1]}"
            )


def verify_downstream_firewall(
    candidate: Path,
    snapshot_rows: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    path = candidate / "manifests/downstream_inputs.tsv"
    rows = read_tsv_exact(path, DOWNSTREAM_FIELDS)
    if not rows:
        raise ReleaseContractError("downstream input manifest is empty")
    consumers = set()
    for row in rows:
        consumer = row["consumer_id"]
        if not consumer or consumer in consumers:
            raise ReleaseContractError(
                f"blank or duplicate downstream consumer_id: {consumer!r}"
            )
        consumers.add(consumer)
        relative = clean_relative_path(
            row["snapshot_path"], f"{consumer} downstream snapshot_path"
        )
        if not relative.startswith("inputs/"):
            raise ReleaseContractError(
                f"downstream consumer is not pinned to candidate inputs/: {consumer}"
            )
        if relative not in snapshot_rows:
            raise ReleaseContractError(
                f"downstream consumer references an unmanifested snapshot: {relative}"
            )
        source_row = snapshot_rows[relative]
        if (
            row["sha256"] != source_row["snapshot_sha256"]
            or row["bytes"] != source_row["snapshot_bytes"]
            or row["workstream_id"] != source_row["workstream_id"]
            or row["artifact_id"] != source_row["artifact_id"]
        ):
            raise ReleaseContractError(
                f"downstream manifest differs from frozen snapshot: {consumer}"
            )
        destination = candidate_file(
            candidate, relative, f"{consumer} downstream snapshot"
        )
        if not destination.is_file() or destination.is_symlink():
            raise ReleaseContractError(
                f"downstream snapshot is missing or symlinked: {destination}"
            )
    manifested_inputs = {row["snapshot_path"] for row in rows}
    observed_inputs = {
        path.relative_to(candidate).as_posix()
        for path in (candidate / "inputs").rglob("*")
        if path.is_file()
        and not path.relative_to(candidate).as_posix().startswith("inputs/BASE/")
    }
    if manifested_inputs != observed_inputs:
        raise ReleaseContractError(
            "candidate inputs contain unregistered or missing files: "
            f"manifest_only={sorted(manifested_inputs - observed_inputs)}, "
            f"candidate_only={sorted(observed_inputs - manifested_inputs)}"
        )
    return rows


def verify_registered_downstream_manifests(candidate: Path) -> None:
    registry = read_tsv_exact(
        candidate / "manifests/downstream_manifest_registry.tsv",
        MANIFEST_REGISTRY_FIELDS,
    )
    if not registry:
        raise ReleaseContractError("downstream manifest registry is empty")
    ids = set()
    for row in registry:
        manifest_id = row["manifest_id"]
        if not manifest_id or manifest_id in ids:
            raise ReleaseContractError(
                f"blank or duplicate downstream manifest_id: {manifest_id!r}"
            )
        ids.add(manifest_id)
        if row["consumer_scope"] != "candidate_only":
            raise ReleaseContractError(
                f"downstream manifest is not candidate-only: {manifest_id}"
            )
        manifest = candidate_file(
            candidate, row["manifest_path"], f"{manifest_id} manifest_path"
        )
        expected_hash = require_sha256(row["sha256"], f"{manifest_id} registry sha256")
        expected_bytes = parse_nonnegative_int(
            row["bytes"], f"{manifest_id} registry bytes"
        )
        if (
            not manifest.is_file()
            or manifest.is_symlink()
            or manifest.stat().st_size != expected_bytes
            or sha256_file(manifest) != expected_hash
        ):
            raise ReleaseContractError(
                f"registered downstream manifest drift: {manifest_id}"
            )
        with manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ReleaseContractError(
                    f"registered downstream manifest lacks a header: {manifest_id}"
                )
            path_columns = [
                field
                for field in reader.fieldnames
                if field.endswith(("_path", "_file", "_input", "_source"))
                or field in {"path", "file", "input", "source"}
            ]
            for manifest_row in reader:
                for field, raw_value in manifest_row.items():
                    value = raw_value.strip()
                    if not value:
                        continue
                    if value.startswith(("/", "file://", "../")) or value.startswith(
                        MUTABLE_PROJECT_PREFIXES
                    ):
                        raise ReleaseContractError(
                            "registered downstream manifest contains a mutable/live "
                            f"reference: {manifest_id} {field}={value!r}"
                        )
                for field in path_columns:
                    value = manifest_row[field].strip()
                    if not value:
                        continue
                    candidate_file(
                        candidate,
                        value,
                        f"{manifest_id} {field}",
                    )


def verify_scientific_registry(
    candidate: Path,
    spec_hash: str,
    fixture_mode: bool,
) -> dict[str, int]:
    path = candidate / "manifests/scientific_validation_registry.tsv"
    rows = read_tsv_exact(path, SCIENTIFIC_REGISTRY_FIELDS)
    counts = {check_id: 0 for check_id in REQUIRED_SCIENTIFIC_CHECKS}
    for row in rows:
        check_id = row["check_id"]
        if check_id not in counts:
            raise ReleaseContractError(
                f"unexpected scientific validation check: {check_id}"
            )
        counts[check_id] += 1
    if any(count != 1 for count in counts.values()) or len(rows) != len(counts):
        raise ReleaseContractError(
            "scientific validation registry must contain every required check "
            f"exactly once: {counts}"
        )
    warning_total = 0
    for row in rows:
        check_id = row["check_id"]
        if row["status"] != "pass":
            raise ReleaseContractError(
                f"scientific validation check is not pass: {check_id}={row['status']}"
            )
        report = candidate_file(
            candidate, row["report_path"], f"{check_id} report_path"
        )
        report_hash = require_sha256(row["report_sha256"], f"{check_id} report_sha256")
        if (
            not report.is_file()
            or report.is_symlink()
            or sha256_file(report) != report_hash
        ):
            raise ReleaseContractError(f"scientific report drift: {check_id}")
        if not fixture_mode:
            payload = load_json(report, f"{check_id} scientific report")
            expected_keys = {
                "contract",
                "candidate_id",
                "check_id",
                "status",
                "input_snapshot_spec_sha256",
                "validator_producer_binding",
                "metrics",
                "warnings",
                "canonical_promotion_authorized",
            }
            if set(payload) != expected_keys:
                raise ReleaseContractError(
                    f"scientific report schema drift: {check_id}"
                )
            if (
                payload["contract"] != SCIENTIFIC_REPORT_CONTRACT
                or payload["candidate_id"] != CANDIDATE_ID
                or payload["check_id"] != check_id
                or payload["status"] != "pass"
                or payload["input_snapshot_spec_sha256"] != spec_hash
                or payload["canonical_promotion_authorized"] is not False
                or not isinstance(payload["metrics"], dict)
                or not isinstance(payload["warnings"], list)
            ):
                raise ReleaseContractError(
                    f"scientific report contract/specification drift: {check_id}"
                )
            binding = payload["validator_producer_binding"]
            if not isinstance(binding, dict) or set(binding) != {
                "producer_manifest_snapshot_path",
                "producer_manifest_sha256",
                "producer_bindings",
            }:
                raise ReleaseContractError(
                    f"scientific report producer-binding schema drift: {check_id}"
                )
            producer_manifest = candidate_file(
                candidate,
                str(binding["producer_manifest_snapshot_path"]),
                f"{check_id} producer manifest",
            )
            producer_manifest_sha = require_sha256(
                str(binding["producer_manifest_sha256"]),
                f"{check_id} producer manifest SHA256",
            )
            if sha256_file(producer_manifest) != producer_manifest_sha:
                raise ReleaseContractError(
                    f"scientific producer manifest drift: {check_id}"
                )
            producer_rows = read_tsv_exact(producer_manifest, RECURSIVE_PRODUCER_FIELDS)
            producer_by_path = {row["repository_path"]: row for row in producer_rows}
            bindings = binding["producer_bindings"]
            if not isinstance(bindings, list) or {
                str(item.get("repository_path", ""))
                for item in bindings
                if isinstance(item, dict)
            } != set(SCIENTIFIC_REPORT_PRODUCER_PATHS):
                raise ReleaseContractError(
                    f"scientific report producer set drift: {check_id}"
                )
            for item in bindings:
                if not isinstance(item, dict) or set(item) != {
                    "producer_id",
                    "repository_path",
                    "sha256",
                    "bytes",
                }:
                    raise ReleaseContractError(
                        f"scientific report producer schema drift: {check_id}"
                    )
                repository_path = str(item["repository_path"])
                frozen = producer_by_path.get(repository_path)
                if frozen is None or (
                    str(item["producer_id"]) != frozen["producer_id"]
                    or require_sha256(
                        str(item["sha256"]),
                        f"{check_id} producer {repository_path}",
                    )
                    != frozen["sha256"]
                    or parse_nonnegative_int(
                        item["bytes"], f"{check_id} producer {repository_path} bytes"
                    )
                    != parse_nonnegative_int(
                        frozen["bytes"], f"frozen producer {repository_path} bytes"
                    )
                ):
                    raise ReleaseContractError(
                        f"scientific report producer binding drift: {check_id}:{repository_path}"
                    )
        warnings = parse_nonnegative_int(
            row["warnings_count"], f"{check_id} warnings_count"
        )
        warning_total += warnings
        adjudication_path = row["adjudication_path"].strip()
        adjudication_hash = row["adjudication_sha256"].strip()
        if warnings:
            if not adjudication_path:
                raise ReleaseContractError(
                    f"warnings lack written adjudication: {check_id}"
                )
            adjudication = candidate_file(
                candidate,
                adjudication_path,
                f"{check_id} adjudication_path",
            )
            require_sha256(adjudication_hash, f"{check_id} adjudication_sha256")
            if (
                not adjudication.is_file()
                or adjudication.is_symlink()
                or sha256_file(adjudication) != adjudication_hash
            ):
                raise ReleaseContractError(f"warning adjudication drift: {check_id}")
        elif adjudication_path or adjudication_hash:
            raise ReleaseContractError(
                f"zero-warning check carries an unexplained adjudication: {check_id}"
            )
    return {"n_checks": len(rows), "n_warnings": warning_total}


def validate_rel02_04_products(
    project: Path,
    candidate: Path,
    spec_hash: str,
    fixture_mode: bool,
) -> dict[str, int]:
    load_frozen_input_index(project, candidate, spec_hash, fixture_mode)
    transition_path = candidate / "manifests/rel02_04_transition.json"
    transition = read_json_object(transition_path, "REL-02--04 transition")
    expected_transition = {
        "contract_version": 1,
        "candidate_id": CANDIDATE_ID,
        "transition": "REL02_04_PRODUCTS_FROZEN",
        "from_candidate_state_sha256": sha256_file(
            candidate / "manifests/candidate_state.json"
        ),
        "input_snapshot_spec_sha256": spec_hash,
        "journal_branch": JOURNAL_BRANCH,
        "figure_count": 5,
        "myojin_role": MYOJIN_ROLE,
        "fixture_mode": fixture_mode,
        "canonical_promotion_status": "not_promoted",
        "full_release_pass": False,
    }
    for key, value in expected_transition.items():
        if transition.get(key) != value:
            raise ReleaseContractError(
                f"REL-02--04 transition drift for {key}: {transition.get(key)!r}"
            )
    transition_hash_paths = {
        "base_input_transition_sha256": candidate
        / "manifests/base_input_transition.json",
        "rel02_state_sha256": candidate / "manifests/rel02_state.json",
        "rel03_state_sha256": candidate / "manifests/rel03_state.json",
        "manuscript_manifest_sha256": candidate / "manifests/manuscript_manifest.tsv",
    }
    for field, path in transition_hash_paths.items():
        if transition.get(field) != sha256_file(path):
            raise ReleaseContractError(f"REL-02--04 transition hash drift: {field}")
    product_manifest_path = candidate / "manifests/rel02_04_product_manifest.tsv"
    if transition.get("product_manifest_sha256") != sha256_file(product_manifest_path):
        raise ReleaseContractError("REL-02--04 product manifest hash drift")
    products = read_tsv_exact(product_manifest_path, TRANSITION_PRODUCT_FIELDS)
    if transition.get("n_products") != len(products):
        raise ReleaseContractError("REL-02--04 product count drift")
    product_ids = set()
    product_paths = set()
    allowed_roots = (
        candidate,
        candidate_figure_root(project),
        candidate_manuscript_root(project),
    )
    for row in products:
        product_id = row["product_id"]
        if not product_id or product_id in product_ids:
            raise ReleaseContractError(
                f"blank or duplicate transition product: {product_id}"
            )
        product_ids.add(product_id)
        path = resolve_project_path(
            project,
            row["project_relative_path"],
            f"transition product {product_id}",
        )
        if not any(is_relative_to(path, root) for root in allowed_roots):
            raise ReleaseContractError(
                f"transition product escaped candidate-only roots: {product_id}"
            )
        expected_hash = require_sha256(row["sha256"], f"{product_id} SHA256")
        expected_bytes = parse_nonnegative_int(row["bytes"], f"{product_id} bytes")
        if (
            not path.is_file()
            or path.is_symlink()
            or path.stat().st_size != expected_bytes
            or sha256_file(path) != expected_hash
        ):
            raise ReleaseContractError(f"transition product drift: {product_id}")
        if row["input_snapshot_spec_sha256"] != spec_hash:
            raise ReleaseContractError(
                f"transition product uses a different input snapshot: {product_id}"
            )
        product_paths.add(path)

    numbers = read_tsv_exact(candidate / "tables/numbers_ledger.tsv", NUMBERS_FIELDS)
    claims = read_tsv_exact(candidate / "claims/claim_ledger.tsv", CLAIM_FIELDS)
    dependencies = read_tsv_exact(
        candidate / "claims/source_dependency_ledger.tsv",
        SOURCE_DEPENDENCY_FIELDS,
    )
    number_by_id = {row["number_id"]: row for row in numbers}
    claim_by_id = {row["claim_id"]: row for row in claims}
    dependency_by_claim = {row["claim_id"]: row for row in dependencies}
    if (
        len(number_by_id) != len(numbers)
        or len(claim_by_id) != len(claims)
        or len(dependency_by_claim) != len(dependencies)
    ):
        raise ReleaseContractError(
            "number, claim, or source-dependency ledger contains duplicate IDs"
        )
    if set(dependency_by_claim) != set(claim_by_id):
        raise ReleaseContractError(
            "source-dependency ledger does not contain exactly one row per claim"
        )
    for claim_id, dependency in dependency_by_claim.items():
        if (
            dependency["source_dependence_class"]
            != claim_by_id[claim_id]["source_dependence"]
        ):
            raise ReleaseContractError(
                f"claim/source-dependency class mismatch: {claim_id}"
            )
        for field in (
            "discovery_source",
            "evaluation_source",
            "reuse_detail",
            "independence_boundary",
        ):
            if not dependency[field].strip():
                raise ReleaseContractError(
                    f"source-dependency ledger has blank {field}: {claim_id}"
                )
        validate_source_dependency_contract(
            dependency["source_dependence_class"],
            dependency["discovery_source"],
            dependency["evaluation_source"],
            dependency["reuse_detail"],
            dependency["independence_boundary"],
            f"source-dependency ledger claim {claim_id}",
        )
    for number in numbers:
        if number["claim_id"] not in claim_by_id:
            raise ReleaseContractError(
                f"number lacks a claim-ledger row: {number['number_id']}"
            )
        source = candidate / number["source_artifact"]
        if (
            not number["source_artifact"].startswith("inputs/")
            or not source.is_file()
            or source.is_symlink()
        ):
            raise ReleaseContractError(
                f"number reads a mutable/non-snapshot source: {number['number_id']}"
            )

    figure_manifest = read_tsv_exact(
        candidate / "manifests/figure_source_manifest.tsv",
        FIGURE_SOURCE_MANIFEST_FIELDS,
    )
    main_figures = {
        row["figure_id"]
        for row in figure_manifest
        if row["figure_id"] in EXPECTED_FIGURES
    }
    if main_figures != set(EXPECTED_FIGURES):
        raise ReleaseContractError(
            f"figure manifest does not contain locked five figures: {main_figures}"
        )
    supplementary_panels = {
        (row["figure_id"], row["panel_id"])
        for row in figure_manifest
        if row["figure_id"].startswith("FigureS")
    }
    expected_supplementary_panels = {
        (figure_id, panel[0])
        for figure_id, panels in EXPECTED_SUPPLEMENTARY_PANELS.items()
        for panel in panels
    }
    if supplementary_panels != expected_supplementary_panels:
        raise ReleaseContractError(
            "figure manifest does not contain exactly the locked supplementary panels: "
            f"observed={sorted(supplementary_panels)}"
        )
    if {figure_id for figure_id, _ in supplementary_panels} != set(
        EXPECTED_SUPPLEMENTARY_FIGURES
    ):
        raise ReleaseContractError("supplementary figure identity drift")
    if transition.get("n_panels") != len(figure_manifest):
        raise ReleaseContractError("REL-02--04 transition panel count drift")
    if any(row["figure_id"] == "Figure6" for row in figure_manifest):
        raise ReleaseContractError("five-figure branch contains Figure6")
    panel_ids = set()
    plotted_number_ids = set()
    all_source_rows: list[dict[str, str]] = []
    for panel in figure_manifest:
        identity = (panel["figure_id"], panel["panel_id"])
        if identity in panel_ids:
            raise ReleaseContractError(f"duplicate figure panel: {identity}")
        panel_ids.add(identity)
        if panel["normal_font_pt"] != "6" or panel["format"] != "PDF":
            raise ReleaseContractError(f"figure typography/format drift: {identity}")
        if panel["figure_id"] == "Figure4" and panel["figure_order"] != "4":
            raise ReleaseContractError("Figure 4 order drift")
        source = candidate / panel["source_table_path"]
        if (
            not panel["source_table_path"].startswith("figure_sources/")
            or sha256_file(source) != panel["source_table_sha256"]
        ):
            raise ReleaseContractError(f"figure source table drift: {identity}")
        source_rows = read_tsv_exact(source, SOURCE_ROW_FIELDS)
        if not source_rows:
            raise ReleaseContractError(f"empty figure source table: {identity}")
        for row in source_rows:
            if (
                row["figure_id"] != panel["figure_id"]
                or row["panel_id"] != panel["panel_id"]
            ):
                raise ReleaseContractError(f"figure/source identity drift: {identity}")
            number_id = f"N::{row['record_id']}"
            number = number_by_id.get(number_id)
            if number is None:
                raise ReleaseContractError(
                    f"plotted number lacks ledger row: {number_id}"
                )
            expected_number = {
                "claim_id": row["claim_id"],
                "number_role": row["number_role"],
                "raw_value": row["value"],
                "display_value": row["display_value"],
                "figure_id": row["figure_id"],
                "panel_id": row["panel_id"],
                "source_artifact": row["source_snapshot_path"],
                "source_row_filter": f"record_id={row['record_id']}",
            }
            if any(number[field] != value for field, value in expected_number.items()):
                raise ReleaseContractError(
                    f"figure/number ledger mismatch: {number_id}"
                )
            plotted_number_ids.add(number_id)
            all_source_rows.append(row)
        pdf = resolve_project_path(project, panel["pdf_path"], f"{identity} PDF")
        if not is_relative_to(pdf, candidate_figure_root(project)):
            raise ReleaseContractError(f"PDF outside candidate misc root: {identity}")
        if (
            pdf.suffix.lower() != ".pdf"
            or sha256_file(pdf) != panel["pdf_sha256"]
            or pdf.stat().st_size != int(panel["pdf_bytes"])
        ):
            raise ReleaseContractError(f"PDF manifest mismatch: {identity}")
        payload = pdf.read_bytes()
        if b"/BaseFont /Helvetica" not in payload or b"/F1 6 Tf" not in payload:
            raise ReleaseContractError(
                f"PDF typography is not normal 6-point: {identity}"
            )
    if plotted_number_ids != {
        number["number_id"] for number in numbers if number["figure_id"]
    }:
        raise ReleaseContractError("figure and numbers-ledger plotted universes differ")
    rows_by_claim: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in all_source_rows:
        rows_by_claim[row["claim_id"]].append(row)
        source = candidate / row["source_snapshot_path"]
        source_hash = require_sha256(
            row["source_sha256"], f"{row['record_id']} frozen-source SHA256"
        )
        if (
            not row["source_snapshot_path"].startswith("inputs/")
            or not source.is_file()
            or source.is_symlink()
            or sha256_file(source) != source_hash
        ):
            raise ReleaseContractError(
                f"figure row does not resolve to its frozen source: {row['record_id']}"
            )
    if set(rows_by_claim) != set(claim_by_id):
        raise ReleaseContractError(
            "claim ledger and figure-source claim universes differ"
        )
    for claim_id, claim_rows in rows_by_claim.items():
        claim = claim_by_id[claim_id]
        expected_number_ids = sorted(f"N::{row['record_id']}" for row in claim_rows)
        expected_artifact_ids = sorted({row["artifact_id"] for row in claim_rows})
        statuses = sorted({row["evidence_status"] for row in claim_rows})
        expected_status = statuses[0] if len(statuses) == 1 else "mixed"
        if (
            sorted(item for item in claim["number_ids"].split(";") if item)
            != expected_number_ids
            or sorted(item for item in claim["source_artifact_ids"].split(";") if item)
            != expected_artifact_ids
            or claim["claim_status"] != expected_status
        ):
            raise ReleaseContractError(
                f"claim ledger does not reproduce its source rows: {claim_id}"
            )
        dependency = dependency_by_claim[claim_id]
        source_contract = {
            "source_dependence_class": {row["source_dependence"] for row in claim_rows},
            "discovery_source": {row["discovery_sources"] for row in claim_rows},
            "evaluation_source": {row["evaluation_sources"] for row in claim_rows},
            "reuse_detail": {row["reuse_detail"] for row in claim_rows},
            "independence_boundary": {
                row["independence_boundary"] for row in claim_rows
            },
        }
        drift = {
            field: values
            for field, values in source_contract.items()
            if len(values) != 1
        }
        if drift or any(
            dependency[field] != next(iter(values))
            for field, values in source_contract.items()
        ):
            raise ReleaseContractError(
                f"source-dependency ledger does not reproduce claim rows: "
                f"{claim_id}; drift={drift}"
            )
    figure_tree = {
        path for path in candidate_figure_root(project).rglob("*") if path.is_file()
    }
    manifested_pdfs = {
        resolve_project_path(project, row["pdf_path"], "manifested panel PDF")
        for row in figure_manifest
    }
    if figure_tree != manifested_pdfs or any(
        path.suffix.lower() in {".png", ".svg"} for path in figure_tree
    ):
        raise ReleaseContractError(
            "candidate figure tree must contain exactly manifested individual PDFs"
        )

    yak_rows = [row for row in all_source_rows if row["group_id"] == "yak_module8"]
    if yak_rows:
        roles = {row["number_role"] for row in yak_rows}
        if roles != {
            "descriptive_median_slope_direction",
            "inferential_signed_stouffer_direction",
        }:
            raise ReleaseContractError(
                "Yak module 8 direction summaries were collapsed"
            )
    for row in all_source_rows:
        if row["evidence_status"] == "tested_negative":
            if not row["negative_adequacy_criterion"].strip() or not any(
                token in row["allowed_wording"].lower()
                for token in ("adequate-negative", "equivalence criterion")
            ):
                raise ReleaseContractError(
                    f"tested_negative lacks adequacy criterion: {row['record_id']}"
                )
        elif row["negative_adequacy_criterion"].strip():
            raise ReleaseContractError(
                f"negative adequacy criterion/status mismatch: {row['record_id']}"
            )
        combined = " ".join(
            row[field]
            for field in (
                "label",
                "category",
                "model_contrast",
                "effect_unit",
                "allowed_wording",
                "claim_text",
            )
        ).lower()
        if "broad validation" in combined:
            raise ReleaseContractError(
                f"broad-validation overclaim: {row['record_id']}"
            )
        if any(token in combined for token in ("vu array", "cosmx array")) and (
            "donor unknown" not in row["biological_unit"].lower()
            or "technical" not in row["unit"].lower()
        ):
            raise ReleaseContractError(
                f"array technical unit misreported as biological n: {row['record_id']}"
            )
        effect_text = row["effect_unit"].lower()
        if "moran" in combined and any(
            token in effect_text and f"not {token}" not in effect_text
            for token in ("sampling se", "sampling ci")
        ):
            raise ReleaseContractError(
                f"matched-null dispersion mislabeled as sampling uncertainty: {row['record_id']}"
            )
        if row["section_id"] == "supplementary_myojin":
            if (
                row["figure_id"] != "FigureS1"
                or row["evidence_status"] not in {"indeterminate", "untestable"}
                or not any(
                    token in combined
                    for token in ("assay-specific non-support", "nonconfirmatory")
                )
            ):
                raise ReleaseContractError("Myojin supplement-only semantic drift")
            if row["evidence_status"] == "untestable" and (
                "untestable" not in row["display_value"].lower()
                or row["p_value"].strip()
                or row["q_value"].strip()
                or not any(
                    token in row["unit"].lower()
                    for token in ("not an effect", "not a zero effect")
                )
                or "not applicable" not in row["effect_unit"].lower()
                or row["value"] != "1"
                or "testability_indicator" not in row["number_role"]
            ):
                raise ReleaseContractError(
                    "Myojin untestable row is encoded as an observed effect"
                )
        if row["section_id"] == "supplementary_nmf":
            if row["figure_id"] != "FigureS2":
                raise ReleaseContractError("continuous NMF supplement panel drift")
            validate_nmf_source_language(row)
            if contains_nmf_discrete_class_claim(combined):
                raise ReleaseContractError("continuous NMF supplement semantic drift")

    observed_supplementary_sections = {
        row["section_id"]
        for row in all_source_rows
        if row["section_id"].startswith("supplementary_")
    }
    if observed_supplementary_sections != set(EXPECTED_SUPPLEMENTARY_SECTIONS):
        raise ReleaseContractError(
            "supplementary source section identity drift: "
            f"{sorted(observed_supplementary_sections)}"
        )

    fibrosis_contract = validate_frozen_fibrosis_bundle(
        candidate, fixture_mode=fixture_mode
    )
    fibrosis_summary = fibrosis_contract["summary"]
    if not isinstance(fibrosis_summary, dict):
        raise ReleaseContractError("sealed fibrosis summary is not a mapping")
    stage_expected = {
        contrast: {
            "deg_count": str(fibrosis_summary[contrast]["n_deg_05"]),
            "sample_count": str(spec["n_samples"]),
            "cohort_count": str(len(spec["cohorts"])),
        }
        for contrast, spec in TRUE_KLEINER_TRANSITIONS.items()
    }
    stage_observed: dict[str, dict[str, str]] = defaultdict(dict)
    for row in all_source_rows:
        if row["model_contrast"] in stage_expected:
            stage_observed[row["model_contrast"]][row["number_role"]] = row["value"]
            if "cross-sectional" not in row["allowed_wording"].lower():
                raise ReleaseContractError("stage count lacks cross-sectional wording")
    if stage_observed != stage_expected:
        raise ReleaseContractError(
            f"sealed true-Kleiner stage-count source drift: {dict(stage_observed)}"
        )

    manuscript_manifest_rows = read_tsv_exact(
        candidate / "manifests/manuscript_manifest.tsv", MANUSCRIPT_MANIFEST_FIELDS
    )
    if len(manuscript_manifest_rows) != 1:
        raise ReleaseContractError("manuscript manifest requires exactly one document")
    manuscript_manifest = manuscript_manifest_rows[0]
    if (
        manuscript_manifest["journal_branch"] != JOURNAL_BRANCH
        or manuscript_manifest["figure_count"] != "5"
        or manuscript_manifest["myojin_role"] != MYOJIN_ROLE
    ):
        raise ReleaseContractError("manuscript branch metadata drift")
    manuscript = resolve_project_path(
        project, manuscript_manifest["candidate_path"], "candidate manuscript"
    )
    if not is_relative_to(manuscript, candidate_manuscript_root(project)):
        raise ReleaseContractError("manuscript escaped candidate-only root")
    if sha256_file(manuscript) != manuscript_manifest[
        "sha256"
    ] or manuscript.stat().st_size != int(manuscript_manifest["bytes"]):
        raise ReleaseContractError("manuscript manifest mismatch")
    manuscript_text = manuscript.read_text(encoding="utf-8")
    marked_numbers = set(re.findall(r"NUMBER:(N::[A-Za-z0-9_.:-]+)", manuscript_text))
    marked_claims = set(re.findall(r"CLAIM:([A-Za-z0-9_.:-]+)", manuscript_text))
    expected_manuscript_numbers = {
        row["number_id"] for row in numbers if row["manuscript_included"] == "true"
    }
    if marked_numbers != expected_manuscript_numbers:
        raise ReleaseContractError(
            "manuscript number markers differ from numbers-ledger inclusion set"
        )
    if marked_numbers != set(manuscript_manifest["number_ids"].split(";")):
        raise ReleaseContractError("manuscript manifest number IDs drift")
    if marked_claims != set(manuscript_manifest["claim_ids"].split(";")):
        raise ReleaseContractError("manuscript manifest claim IDs drift")
    for line in manuscript_text.splitlines():
        ids = re.findall(r"NUMBER:(N::[A-Za-z0-9_.:-]+)", line)
        for number_id in ids:
            if number_by_id[number_id]["display_value"] not in line:
                raise ReleaseContractError(
                    f"manuscript displayed number differs from ledger: {number_id}"
                )
        visible = re.sub(r"<!--.*?-->", "", line)
        # k4/k6 are prespecified NMF model labels, not reported numerical
        # results.  Remove only those exact tokens before enforcing the
        # unledgered-number gate; any other digit on the same line remains
        # visible to the validator.
        visible_numeric = re.sub(r"\bk(?:4|6)\b", "k", visible, flags=re.IGNORECASE)
        if (
            re.search(r"\d", visible_numeric)
            and not ids
            and not re.match(r"^\d+\. ", visible_numeric)
        ):
            raise ReleaseContractError(
                f"unledgered manuscript number: {visible.strip()!r}"
            )
    abstract = manuscript_text.split("## Abstract", 1)[1].split("## Introduction", 1)[0]
    if any(
        name.lower() in abstract.lower()
        for name in ("yakubovsky", "gse287826", "myojin")
    ):
        raise ReleaseContractError("abstract names a processed external dataset")
    if "Figure 6" in manuscript_text or "Figure6" in manuscript_text:
        raise ReleaseContractError("five-figure manuscript contains Figure 6")
    manuscript_tree = {
        path for path in candidate_manuscript_root(project).rglob("*") if path.is_file()
    }
    if manuscript_tree != {manuscript}:
        raise ReleaseContractError(
            "candidate manuscript tree contains unmanifested products"
        )
    analysis_rows = read_tsv_exact(
        candidate / "manifests/analysis_release_manifest.tsv",
        ANALYSIS_MANIFEST_FIELDS,
    )
    expected_product_paths = {
        candidate / row["candidate_path"]
        for row in analysis_rows
        if row["producer_id"] == "rel02_build_candidate_tables"
    }
    expected_product_paths.update(
        {
            candidate / "manifests/rel02_state.json",
            candidate / "manifests/figure_source_manifest.tsv",
            candidate / "manifests/rel03_state.json",
            candidate / "manifests/manuscript_manifest.tsv",
            manuscript,
            *manifested_pdfs,
        }
    )
    if product_paths != expected_product_paths:
        raise ReleaseContractError(
            "REL-02--04 transition manifest does not exactly cover candidate products"
        )
    return {
        "n_products": len(products),
        "n_panels": len(figure_manifest),
        "n_numbers": len(numbers),
        "n_claims": len(claims),
        "n_manuscript_numbers": len(marked_numbers),
    }


def validate_snapshot_base(
    project_root: Path,
    fixture_mode: bool = False,
) -> dict[str, object]:
    project = project_root.resolve()
    candidate = assert_candidate_root(project, candidate_root(project))
    if not candidate.is_dir() or candidate.is_symlink():
        raise ReleaseContractError(f"candidate snapshot does not exist: {candidate}")
    assert_no_symlinks(candidate)
    for relative in REQUIRED_CANDIDATE_DIRS:
        path = candidate / relative
        if not path.is_dir() or path.is_symlink():
            raise ReleaseContractError(
                f"required candidate directory is missing or symlinked: {relative}"
            )
    specification = load_json(
        candidate / "manifests/snapshot_spec.json", "snapshot spec"
    )
    spec = specification.get("spec")
    spec_hash = specification.get("spec_sha256")
    if not isinstance(spec, dict) or not isinstance(spec_hash, str):
        raise ReleaseContractError("snapshot spec lacks spec/spec_sha256")
    require_sha256(spec_hash, "snapshot spec_sha256")
    if sha256_bytes(canonical_json_bytes(spec)) != spec_hash:
        raise ReleaseContractError("snapshot spec hash does not reproduce")
    if spec.get("candidate_id") != CANDIDATE_ID:
        raise ReleaseContractError("snapshot candidate ID drift")
    if bool(spec.get("fixture_mode")) != fixture_mode:
        raise ReleaseContractError(
            "fixture-mode mismatch between validator and frozen snapshot"
        )
    rebuild_spec, rebuild_spec_hash = load_rebuild_source_spec(
        candidate, fixture_mode=fixture_mode
    )
    if (
        specification.get("rebuild_source_spec") != rebuild_spec
        or specification.get("rebuild_source_spec_sha256") != rebuild_spec_hash
        or spec.get("rebuild_source_spec_sha256") != rebuild_spec_hash
    ):
        raise ReleaseContractError(
            "snapshot specification does not bind the rebuild source identity"
        )
    expected_rebuild_projection = {
        "repository_commit": str(
            spec.get("repository", {}).get("repository_commit", "")
        )
        if isinstance(spec.get("repository"), dict)
        else "",
        "closure_sha256": spec.get("closure_sha256"),
        "workstreams": sorted(
            spec.get("workstreams", []), key=lambda row: str(row["workstream_id"])
        ),
        "artifacts": sorted(
            spec.get("artifacts", []),
            key=lambda row: (
                str(row["workstream_id"]),
                str(row["artifact_id"]),
            ),
        ),
        "producer_scripts": sorted(
            spec.get("producer_scripts", []),
            key=lambda row: str(row["repository_path"]),
        ),
        "protected_scope_registry_sha256": spec.get("protected_scope_registry_sha256"),
        "protected_baseline_sha256": spec.get("protected_baseline_sha256"),
        "fixture_mode": fixture_mode,
    }
    for key, expected in expected_rebuild_projection.items():
        if rebuild_spec.get(key) != expected:
            raise ReleaseContractError(
                f"rebuild source differs from snapshot source for {key}"
            )
    state = load_json(candidate / "manifests/candidate_state.json", "candidate state")
    expected_state = {
        "candidate_id": CANDIDATE_ID,
        "state": "SNAPSHOT_FROZEN",
        "spec_sha256": spec_hash,
        "rebuild_source_spec_sha256": rebuild_spec_hash,
        "scientific_assembly_performed": False,
        "canonical_promotion_status": "not_promoted",
    }
    for key, expected in expected_state.items():
        if state.get(key) != expected:
            raise ReleaseContractError(
                f"candidate state drift for {key}: {state.get(key)!r}"
            )
    closure = validate_frozen_closure(candidate, spec)
    protected_scopes_path = candidate / "manifests/protected_scopes.tsv"
    protected_baseline_path = candidate / "manifests/protected_release_baseline.tsv"
    if sha256_file(protected_scopes_path) != spec.get(
        "protected_scope_registry_sha256"
    ):
        raise ReleaseContractError(
            "frozen protected-scope registry differs from snapshot spec"
        )
    if sha256_file(protected_baseline_path) != spec.get("protected_baseline_sha256"):
        raise ReleaseContractError(
            "frozen protected baseline differs from snapshot spec"
        )
    repository = load_json(
        candidate / "manifests/repository_state.json", "repository state"
    )
    if repository != spec.get("repository"):
        raise ReleaseContractError("repository state differs from snapshot spec")
    producers = read_tsv_exact(
        candidate / "manifests/producer_script_manifest.tsv",
        PRODUCER_SCRIPT_FIELDS,
    )
    frozen_producers = spec.get("producer_scripts")
    if not isinstance(frozen_producers, list):
        raise ReleaseContractError("snapshot spec lacks producer-script inventory")
    normalized_frozen_producers = [
        {field: str(row[field]) for field in PRODUCER_SCRIPT_FIELDS}
        for row in frozen_producers
        if isinstance(row, dict)
        and all(field in row for field in PRODUCER_SCRIPT_FIELDS)
    ]
    if (
        len(normalized_frozen_producers) != len(frozen_producers)
        or producers != normalized_frozen_producers
    ):
        raise ReleaseContractError(
            "producer script hashes differ from the frozen snapshot specification"
        )
    producer_paths: set[str] = set()
    for producer in producers:
        snapshot_path = clean_relative_path(
            producer["snapshot_path"],
            f"producer {producer['producer_id']} snapshot path",
        )
        if not snapshot_path.startswith("environments/producers/"):
            raise ReleaseContractError(
                f"producer copy escaped environments/producers/: {snapshot_path}"
            )
        if snapshot_path in producer_paths:
            raise ReleaseContractError(
                f"duplicate frozen producer path: {snapshot_path}"
            )
        producer_paths.add(snapshot_path)
        producer_copy = candidate / snapshot_path
        expected_hash = require_sha256(
            producer["sha256"], f"producer {producer['producer_id']} SHA256"
        )
        expected_bytes = parse_nonnegative_int(
            producer["bytes"], f"producer {producer['producer_id']} bytes"
        )
        if (
            not producer_copy.is_file()
            or producer_copy.is_symlink()
            or producer_copy.stat().st_size != expected_bytes
            or sha256_file(producer_copy) != expected_hash
        ):
            raise ReleaseContractError(
                f"frozen producer copy drift: {producer['producer_id']}"
            )
    observed_producer_paths = {
        path.relative_to(candidate).as_posix()
        for path in (candidate / "environments/producers").rglob("*")
        if path.is_file()
    }
    if observed_producer_paths != producer_paths:
        raise ReleaseContractError(
            "frozen producer tree does not exactly match its manifest"
        )
    runtime_environment = load_json(
        candidate / "environments/runtime_environment.json", "runtime environment"
    )
    if runtime_environment != spec.get("runtime_environment"):
        raise ReleaseContractError("runtime-environment snapshot differs from spec")
    python_packages = read_tsv_exact(
        candidate / "environments/python_packages.tsv", PYTHON_PACKAGE_FIELDS
    )
    frozen_packages = spec.get("python_packages")
    if not isinstance(frozen_packages, list) or python_packages != [
        {field: str(row[field]) for field in PYTHON_PACKAGE_FIELDS}
        for row in frozen_packages
        if isinstance(row, dict)
        and all(field in row for field in PYTHON_PACKAGE_FIELDS)
    ]:
        raise ReleaseContractError("Python package export differs from snapshot spec")
    environment_index_path = candidate / "environments/environment_index.tsv"
    environment_rows = read_tsv_exact(environment_index_path, ENVIRONMENT_INDEX_FIELDS)
    expected_environment_paths = {
        "environments/runtime_environment.json",
        "environments/python_packages.tsv",
    }
    if {row["snapshot_path"] for row in environment_rows} != expected_environment_paths:
        raise ReleaseContractError(
            "environment index has missing or extra runtime products"
        )
    for row in environment_rows:
        environment_file = candidate_file(
            candidate, row["snapshot_path"], f"environment {row['environment_id']}"
        )
        if sha256_file(environment_file) != require_sha256(
            row["sha256"], f"environment {row['environment_id']} SHA256"
        ) or environment_file.stat().st_size != parse_nonnegative_int(
            row["bytes"], f"environment {row['environment_id']} bytes"
        ):
            raise ReleaseContractError(
                f"frozen environment index drift: {row['environment_id']}"
            )
    execution_path = candidate / "manifests/execution_context.json"
    execution = load_json(execution_path, "execution context")
    if bool(execution.get("fixture_mode")) != fixture_mode:
        raise ReleaseContractError("execution-context fixture mode drift")
    for workstream in closure.workstreams:
        copied = candidate / f"manifests/upstream/{workstream.row['workstream_id']}.tsv"
        if sha256_file(copied) != workstream.manifest_sha256:
            raise ReleaseContractError(
                f"copied upstream terminal manifest drift: {workstream.row['workstream_id']}"
            )
    protected = verify_protected_baseline(
        project,
        protected_scopes_path,
        protected_baseline_path,
    )
    snapshot_rows, snapshot_by_path = verify_snapshot_manifest(candidate)
    verify_spec_and_snapshot_contract(
        project,
        candidate,
        closure,
        spec,
        snapshot_rows,
    )
    downstream_rows = verify_downstream_firewall(candidate, snapshot_by_path)
    verify_registered_downstream_manifests(candidate)
    snapshot_path = candidate / "manifests/input_snapshot_manifest.tsv"
    downstream_path = candidate / "manifests/downstream_inputs.tsv"
    registry_path = candidate / "manifests/downstream_manifest_registry.tsv"
    if state.get("snapshot_manifest_sha256") != sha256_file(snapshot_path):
        raise ReleaseContractError("candidate-state snapshot manifest hash drift")
    if state.get("downstream_manifest_sha256") != sha256_file(downstream_path):
        raise ReleaseContractError("candidate-state downstream manifest hash drift")
    if state.get("downstream_registry_sha256") != sha256_file(registry_path):
        raise ReleaseContractError("candidate-state downstream registry hash drift")
    if state.get("repository_state_sha256") != sha256_file(
        candidate / "manifests/repository_state.json"
    ):
        raise ReleaseContractError("candidate-state repository hash drift")
    rebuild_source_path = candidate / "manifests/rebuild_source_spec.json"
    if (
        state.get("rebuild_source_spec_file_sha256") != sha256_file(rebuild_source_path)
        or state.get("rebuild_source_spec_sha256") != rebuild_spec_hash
    ):
        raise ReleaseContractError("candidate-state rebuild source identity drift")
    if state.get("producer_script_manifest_sha256") != sha256_file(
        candidate / "manifests/producer_script_manifest.tsv"
    ):
        raise ReleaseContractError("candidate-state producer-script hash drift")
    if state.get("environment_index_sha256") != sha256_file(
        candidate / "environments/environment_index.tsv"
    ):
        raise ReleaseContractError("candidate-state environment-index hash drift")
    if state.get("execution_context_sha256") != sha256_file(execution_path):
        raise ReleaseContractError("candidate-state execution-context hash drift")
    plan50_relocated_verification = verify_plan50_relocated_verification(
        candidate,
        state,
        normalized_frozen_producers,
        fixture_mode,
    )
    adapter_provenance = None
    if not fixture_mode:
        adapter_provenance = verify_frozen_adapter_provenance(
            candidate,
            state,
            snapshot_rows,
        )
    return {
        "project": project,
        "candidate": candidate,
        "spec": spec,
        "spec_hash": spec_hash,
        "rebuild_source_spec": rebuild_spec,
        "rebuild_source_spec_sha256": rebuild_spec_hash,
        "state": state,
        "closure": closure,
        "protected": protected,
        "snapshot_rows": snapshot_rows,
        "snapshot_by_path": snapshot_by_path,
        "downstream_rows": downstream_rows,
        "adapter_provenance": adapter_provenance,
        "plan50_relocated_verification": plan50_relocated_verification,
    }


def validate_candidate(
    project_root: Path,
    fixture_mode: bool = False,
    write_report: bool = True,
) -> dict[str, object]:
    base = validate_snapshot_base(project_root, fixture_mode=fixture_mode)
    candidate = base["candidate"]
    if not isinstance(candidate, Path):
        raise ReleaseContractError("internal candidate validation type drift")
    project = base["project"]
    if not isinstance(project, Path):
        raise ReleaseContractError("internal project validation type drift")
    products = validate_rel02_04_products(
        project,
        candidate,
        str(base["spec_hash"]),
        fixture_mode,
    )
    if fixture_mode:
        clean_rebuild_rows = []
    else:
        clean_rebuild_rows = revalidate_clean_rebuild_report(
            project,
            candidate / "manifests/two_clean_rebuild_comparison.tsv",
        )
    scientific = verify_scientific_registry(
        candidate,
        str(base["spec_hash"]),
        fixture_mode,
    )
    assert_no_symlinks(candidate)
    status = "SYNTHETIC_PASS" if fixture_mode else "REL05_INFRASTRUCTURE_PASS"
    report = {
        "status": status,
        "candidate_id": CANDIDATE_ID,
        "spec_sha256": base["spec_hash"],
        "n_closed_workstreams": len(base["closure"].workstreams),
        "n_snapshot_records": len(base["snapshot_rows"]),
        "n_downstream_inputs": len(base["downstream_rows"]),
        "n_protected_scopes": base["protected"]["n_scopes"],
        "n_protected_files": base["protected"]["n_files"],
        "n_scientific_checks": scientific["n_checks"],
        "n_adjudicated_warnings": scientific["n_warnings"],
        "n_transition_products": products["n_products"],
        "n_candidate_panels": products["n_panels"],
        "n_ledger_numbers": products["n_numbers"],
        "n_manuscript_numbers": products["n_manuscript_numbers"],
        "n_clean_rebuild_products": len(clean_rebuild_rows),
        "two_clean_rebuilds_revalidated": not fixture_mode,
        "adapter_producer_bindings_revalidated": not fixture_mode,
        "candidate_products_frozen": True,
        "scientific_assembly_performed": False,
        "canonical_promotion_authorized": False,
        "full_release_pass": False,
        "note": (
            "Infrastructure consistency passed. Full scientific release PASS and "
            "promotion remain separate, explicitly authorized operations."
        ),
    }
    if write_report:
        atomic_write_json(candidate / "manifests/rel05_consistency_report.json", report)
        markdown = "\n".join(
            [
                f"# Candidate consistency scaffold: {CANDIDATE_ID}",
                "",
                f"**Status:** `{status}`",
                "",
                f"- Frozen snapshot records: {len(base['snapshot_rows'])}",
                f"- Candidate-only downstream inputs: {len(base['downstream_rows'])}",
                f"- Protected release files reverified: {base['protected']['n_files']}",
                f"- Registered scientific validation reports: {scientific['n_checks']}",
                f"- Written warning adjudications: {scientific['n_warnings']}",
                f"- Candidate-only REL-02--04 products: {products['n_products']}",
                (
                    "- Retained clean-rebuild products independently revalidated: "
                    f"{len(clean_rebuild_rows)}"
                ),
                f"- Individual PDF panels: {products['n_panels']}",
                f"- Ledger-linked manuscript numbers: {products['n_manuscript_numbers']}",
                "- Scientific assembly performed: no",
                "- Canonical promotion authorized: no",
                "- Full release PASS: no",
                "",
            ]
        )
        report_path = candidate / "candidate_consistency_report.md"
        report_path.write_text(markdown, encoding="utf-8")
    return report


def validate_complete_pre_rel05_candidate(
    project_root: Path,
    fixture_mode: bool = False,
) -> dict[str, object]:
    """Read-only proof that REL01--04 are complete before REL05 artifacts exist.

    This intentionally excludes the scientific-validation registry, the
    two-clean-rebuild report, and REL05 reports/seals.  All immutable snapshot,
    adapter-provenance, staged-input, table, panel, manuscript, and exact
    REL02--04 product-tree checks are reused from the terminal validator.
    """

    base = validate_snapshot_base(project_root, fixture_mode=fixture_mode)
    candidate = base["candidate"]
    project = base["project"]
    if not isinstance(candidate, Path) or not isinstance(project, Path):
        raise ReleaseContractError("internal pre-REL05 candidate type drift")
    products = validate_rel02_04_products(
        project,
        candidate,
        str(base["spec_hash"]),
        fixture_mode,
    )
    assert_no_symlinks(candidate)
    return {
        "project": project,
        "candidate": candidate,
        "spec_hash": str(base["spec_hash"]),
        "closure_sha256": str(base["spec"]["closure_sha256"]),
        "n_snapshot_records": len(base["snapshot_rows"]),
        "n_products": products["n_products"],
        "n_panels": products["n_panels"],
        "n_numbers": products["n_numbers"],
        "n_claims": products["n_claims"],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument(
        "--fixture-mode",
        action="store_true",
        help="Validate an isolated synthetic project root; never a real candidate.",
    )
    parser.add_argument(
        "--no-write-report",
        action="store_true",
        help="Run read-only validation without writing candidate reports.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        report = validate_candidate(
            args.project_root,
            fixture_mode=args.fixture_mode,
            write_report=not args.no_write_report,
        )
    except ReleaseContractError as error:
        raise SystemExit(f"REL05_BLOCKED: {error}") from error
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
