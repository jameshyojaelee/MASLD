#!/usr/bin/env python3
"""Copy a closed Plan 60 input set into the exact read-only candidate root."""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import io
import json
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping

from release_common import (
    CANDIDATE_ID,
    DOWNSTREAM_FIELDS,
    PLAN50_FROZEN_PRODUCER_PATHS,
    PLAN50_MANIFEST_PRODUCER_PATHS,
    REQUIRED_CANDIDATE_DIRS,
    SNAPSHOT_MANIFEST_FIELDS,
    ClosureBundle,
    ReleaseContractError,
    assert_candidate_root,
    assert_no_symlinks,
    candidate_root,
    canonical_json_bytes,
    parse_nonnegative_int,
    project_relative,
    read_tsv_exact,
    repository_state,
    require_sha256,
    resolve_project_path,
    sha256_bytes,
    sha256_file,
    validate_closure,
    verify_protected_baseline,
)
from release_products import BASE_SELECTION_FIELDS
from adapter_contract import (
    ADAPTER_CONTRACT,
    ADAPTER_CONTRACT_SHA256,
    ADAPTER_PROVENANCE_FIELDS,
    ADAPTER_PROVENANCE_VERSION,
    EXPECTED_ADAPTER_IDS,
    NMF_FREEZER_PRODUCER_PATH,
    RECURSIVE_PRODUCER_FIELDS,
    require_exact_adapter_binding_paths,
)
from rebuild_source_spec import (
    REBUILD_SOURCE_FILENAME,
    build_rebuild_source_spec,
    read_base_selection,
)
from publish_plan60_workstream_handoff import PLAN50_PAYLOAD_MANIFEST_FIELDS


MANIFEST_REGISTRY_FIELDS = (
    "manifest_id",
    "manifest_path",
    "sha256",
    "bytes",
    "consumer_scope",
)
PRODUCER_SCRIPT_FIELDS = (
    "producer_id",
    "repository_path",
    "sha256",
    "bytes",
    "snapshot_path",
)
PYTHON_PACKAGE_FIELDS = (
    "distribution",
    "version",
    "location",
    "metadata_path",
)
ENVIRONMENT_INDEX_FIELDS = (
    "environment_id",
    "snapshot_path",
    "sha256",
    "bytes",
    "scope",
)
SAFE_ENVIRONMENT_KEYS = (
    "CONDA_DEFAULT_ENV",
    "CONDA_PREFIX",
    "MAMBA_ROOT_PREFIX",
    "SLURM_CLUSTER_NAME",
    "SLURM_ARRAY_TASK_ID",
    "SLURM_JOB_ID",
    "SLURM_JOB_NAME",
    "SLURM_STEP_ID",
)
PLAN50_TERMINAL_PROVENANCE_FIELDS = (
    "analysis_release_id",
    "bundle_uri",
    "artifact_role",
    "relative_path",
    "sha256",
    "bytes",
    "producer",
    "producer_sha256",
)
PLAN50_RELOCATED_VERIFICATION_REL = "manifests/plan50_relocated_verification.json"


def validate_plan50_frozen_producer_binding(
    closure: ClosureBundle,
    frozen_producers: list[Mapping[str, object]],
) -> dict[str, object]:
    """Bind every sealed Plan 50 producer hash to the REL01 frozen copy."""

    indexed_producers: dict[str, Mapping[str, object]] = {}
    for row in frozen_producers:
        repository_path = str(row["repository_path"])
        if repository_path in indexed_producers:
            raise ReleaseContractError(
                f"duplicate frozen producer repository path: {repository_path}"
            )
        indexed_producers[repository_path] = row
    missing = set(PLAN50_FROZEN_PRODUCER_PATHS) - set(indexed_producers)
    if missing:
        raise ReleaseContractError(
            f"REL01 frozen producer universe lacks Plan 50 producers: {sorted(missing)}"
        )

    plan50 = [
        workstream
        for workstream in closure.workstreams
        if workstream.row["workstream_id"] == "PLAN50"
    ]
    if len(plan50) != 1:
        raise ReleaseContractError("closure must contain exactly one PLAN50 workstream")

    def one_artifact(filename: str) -> Path:
        matches = [
            artifact.source_path
            for artifact in plan50[0].artifacts
            if artifact.source_path.name == filename
        ]
        if len(matches) != 1:
            raise ReleaseContractError(
                f"PLAN50 handoff must contain exactly one {filename}; found {len(matches)}"
            )
        return matches[0]

    payload_manifest = one_artifact("passport_release_manifest.tsv")
    payload_rows = read_tsv_exact(payload_manifest, PLAN50_PAYLOAD_MANIFEST_FIELDS)
    observed = {row["producer"] for row in payload_rows}
    if observed != set(PLAN50_MANIFEST_PRODUCER_PATHS):
        raise ReleaseContractError(
            "sealed PLAN50 payload producer universe drift: "
            f"observed={sorted(observed)}"
        )

    def validate_row(row: Mapping[str, str], context: str) -> None:
        producer = row["producer"]
        if producer not in PLAN50_MANIFEST_PRODUCER_PATHS:
            raise ReleaseContractError(
                f"{context} producer is outside the Plan 50 allowlist: {producer}"
            )
        frozen = indexed_producers[producer]
        sealed_hash = require_sha256(
            row["producer_sha256"], f"{context} producer SHA256"
        )
        if sealed_hash != str(frozen["sha256"]):
            raise ReleaseContractError(
                f"{context} producer hash differs from REL01 frozen bytes: {producer}"
            )

    for row in payload_rows:
        validate_row(row, f"PLAN50 payload {row['relative_path']}")

    terminal_provenance = one_artifact("passport_terminal_provenance.tsv")
    terminal_rows = read_tsv_exact(
        terminal_provenance, PLAN50_TERMINAL_PROVENANCE_FIELDS
    )
    if len(terminal_rows) != 2 or {row["artifact_role"] for row in terminal_rows} != {
        "payload_manifest",
        "plan60_handoff",
    }:
        raise ReleaseContractError(
            "PLAN50 terminal provenance must contain exactly the two terminal roles"
        )
    for row in terminal_rows:
        validate_row(row, f"PLAN50 terminal provenance {row['artifact_role']}")
    return {
        "payload_manifest_sha256": sha256_file(payload_manifest),
        "terminal_provenance_sha256": sha256_file(terminal_provenance),
        "manifest_producer_count": len(observed),
    }


def verify_relocated_plan50_bundle(candidate: Path) -> dict[str, object]:
    """Replay the frozen read-only Plan 50 validator on the relocated snapshot."""

    bundle = candidate / "inputs/PASS"
    validator_rel = "environments/producers/portal/validate_evidence_passports.py"
    validator = candidate / validator_rel
    if not bundle.is_dir() or bundle.is_symlink():
        raise ReleaseContractError(
            f"relocated Plan 50 bundle is missing or unsafe: {bundle}"
        )
    if not validator.is_file() or validator.is_symlink():
        raise ReleaseContractError(
            f"frozen Plan 50 validator is missing or unsafe: {validator}"
        )

    def bundle_state() -> list[tuple[str, str, int]]:
        rows: list[tuple[str, str, int]] = []
        for path in sorted(bundle.rglob("*")):
            if path.is_symlink():
                raise ReleaseContractError(
                    f"relocated Plan 50 bundle contains a symlink: {path}"
                )
            if path.is_file():
                rows.append(
                    (
                        path.relative_to(bundle).as_posix(),
                        sha256_file(path),
                        path.stat().st_size,
                    )
                )
        return rows

    before = bundle_state()
    environment = dict(os.environ)
    environment.update(
        {
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONNOUSERSITE": "1",
            "PYTHONHASHSEED": "0",
        }
    )
    command = [
        sys.executable,
        "-B",
        str(validator),
        "--bundle",
        str(bundle),
        "--verify-sealed",
    ]
    completed = subprocess.run(
        command,
        cwd=validator.parent,
        env=environment,
        text=True,
        capture_output=True,
        timeout=1200,
        check=False,
    )
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        raise ReleaseContractError(
            "frozen relocated Plan 50 validation failed: " + detail[-4000:]
        )
    try:
        result = json.loads(completed.stdout.strip())
    except json.JSONDecodeError as error:
        raise ReleaseContractError(
            "frozen relocated Plan 50 validator emitted invalid JSON"
        ) from error
    if (
        not isinstance(result, list)
        or len(result) != 1
        or not isinstance(result[0], dict)
        or set(result[0]) != {"check_id", "status", "detail"}
        or result[0]["check_id"] != "PASS06_POST_SEAL_READ_ONLY"
        or result[0]["status"] != "pass"
        or not isinstance(result[0]["detail"], str)
        or not result[0]["detail"].strip()
    ):
        raise ReleaseContractError(
            "frozen relocated Plan 50 validator result schema or verdict drift"
        )
    after = bundle_state()
    if before != after:
        raise ReleaseContractError(
            "frozen relocated Plan 50 verification mutated the sealed bundle"
        )
    return {
        "contract_version": 1,
        "status": "pass",
        "bundle_path": "inputs/PASS",
        "validator_path": validator_rel,
        "validator_sha256": sha256_file(validator),
        "seal_sha256": sha256_file(bundle / "PASS06_VALIDATED"),
        "payload_manifest_sha256": sha256_file(
            bundle / "passport_release_manifest.tsv"
        ),
        "terminal_provenance_sha256": sha256_file(
            bundle / "passport_terminal_provenance.tsv"
        ),
        "verification_result_sha256": sha256_bytes(canonical_json_bytes(result)),
        "artifact_count": len(before),
    }


def producer_script_inventory() -> list[dict[str, object]]:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parents[2]
    rows = []
    producer_ids: set[str] = set()
    for path in sorted(script_dir.rglob("*")):
        if (
            not path.is_file()
            or path.is_symlink()
            or "__pycache__" in path.parts
            or path.suffix not in {".py", ".R", ".sh", ".sbatch"}
        ):
            continue
        relative = path.relative_to(script_dir).as_posix()
        repository_path = f"scripts/manuscript/program_context_v2/{relative}"
        producer_id = repository_path.replace("/", "__")
        if producer_id in producer_ids:
            raise ReleaseContractError(
                f"recursive producer inventory has an ID collision: {relative}"
            )
        producer_ids.add(producer_id)
        if not path.is_file() or path.is_symlink():
            raise ReleaseContractError(
                f"release producer script is missing or symlinked: {path}"
            )
        rows.append(
            {
                "producer_id": producer_id,
                "repository_path": repository_path,
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "snapshot_path": f"environments/producers/{relative}",
            }
        )
    for repository_path in PLAN50_FROZEN_PRODUCER_PATHS:
        path = project_root / repository_path
        if not path.is_file() or path.is_symlink():
            raise ReleaseContractError(
                f"Plan 50 producer script is missing or symlinked: {path}"
            )
        producer_id = repository_path.replace("/", "__")
        if producer_id in producer_ids:
            raise ReleaseContractError(
                f"external producer inventory has an ID collision: {repository_path}"
            )
        producer_ids.add(producer_id)
        rows.append(
            {
                "producer_id": producer_id,
                "repository_path": repository_path,
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "snapshot_path": f"environments/producers/portal/{path.name}",
            }
        )
    external_path = project_root / NMF_FREEZER_PRODUCER_PATH
    external_id = NMF_FREEZER_PRODUCER_PATH.replace("/", "__")
    if not external_path.is_file() or external_path.is_symlink():
        raise ReleaseContractError(
            f"external NMF producer is missing or symlinked: {external_path}"
        )
    if external_id in producer_ids:
        raise ReleaseContractError(
            f"external producer inventory has an ID collision: {external_id}"
        )
    producer_ids.add(external_id)
    rows.append(
        {
            "producer_id": external_id,
            "repository_path": NMF_FREEZER_PRODUCER_PATH,
            "sha256": sha256_file(external_path),
            "bytes": external_path.stat().st_size,
            "snapshot_path": (
                "environments/producers/external/"
                "519_freeze_nmf_continuous_supplement.py"
            ),
        }
    )
    return rows


def producer_source_path(row: Mapping[str, object]) -> Path:
    repository_path = str(row["repository_path"])
    permitted = (
        repository_path.startswith("scripts/manuscript/program_context_v2/")
        or repository_path in PLAN50_FROZEN_PRODUCER_PATHS
        or repository_path == NMF_FREEZER_PRODUCER_PATH
    )
    if not permitted:
        raise ReleaseContractError(
            f"producer path is outside the frozen producer allowlist: {repository_path}"
        )
    project_root = Path(__file__).resolve().parents[3]
    source = project_root / repository_path
    if not source.is_file() or source.is_symlink():
        raise ReleaseContractError(f"producer source is missing or symlinked: {source}")
    return source


def environment_inventory(
    project_root: Path,
) -> tuple[dict[str, object], list[dict[str, str]]]:
    project = project_root.resolve()

    def portable_path(value: str) -> str:
        if not value:
            return value
        path = Path(value)
        if not path.is_absolute():
            return value
        resolved = path.resolve()
        try:
            relative = resolved.relative_to(project)
        except ValueError:
            return str(resolved)
        return "${PROJECT_ROOT}/" + relative.as_posix()

    runtime = {
        "python_executable": portable_path(sys.executable),
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "sys_path": [portable_path(path) for path in sys.path],
    }
    packages: list[dict[str, str]] = []
    for distribution in importlib.metadata.distributions():
        name = str(distribution.metadata.get("Name") or "").strip()
        version = str(distribution.version).strip()
        if not name or not version:
            continue
        canonical = re.sub(r"[-_.]+", "-", name).lower()
        location = portable_path(str(Path(distribution.locate_file("")).resolve()))
        metadata_path = portable_path(str(getattr(distribution, "_path", "")))
        packages.append(
            {
                "distribution": canonical,
                "version": version,
                "location": location,
                "metadata_path": metadata_path,
            }
        )
    rows = sorted(
        packages,
        key=lambda row: (
            row["distribution"],
            row["version"],
            row["location"],
            row["metadata_path"],
        ),
    )
    return runtime, rows


def build_snapshot_spec(
    project_root: Path,
    closure: ClosureBundle,
    protected_scopes: Path,
    protected_baseline: Path,
    base_selection_path: Path,
    fixture_mode: bool,
) -> dict[str, object]:
    repository = repository_state(project_root, fixture_mode)
    producers = producer_script_inventory()
    runtime_environment, python_packages = environment_inventory(project_root)
    base_selection_rows, base_selection_sha256 = read_base_selection(
        base_selection_path
    )
    spec = {
        "contract_version": 1,
        "candidate_id": CANDIDATE_ID,
        "closure_sha256": closure.closure_sha256,
        "protected_scope_registry_sha256": sha256_file(protected_scopes),
        "protected_baseline_sha256": sha256_file(protected_baseline),
        "repository": repository,
        "producer_scripts": producers,
        "runtime_environment": runtime_environment,
        "python_packages": python_packages,
        "workstreams": [
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
        ],
        "artifacts": [
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
        ],
        "fixture_mode": fixture_mode,
    }
    rebuild_source = build_rebuild_source_spec(
        repository_commit=repository["repository_commit"],
        closure_sha256=closure.closure_sha256,
        workstreams=spec["workstreams"],
        artifacts=spec["artifacts"],
        base_selection_sha256=base_selection_sha256,
        base_selection_rows=base_selection_rows,
        producers=producers,
        protected_scope_registry_sha256=spec["protected_scope_registry_sha256"],
        protected_baseline_sha256=spec["protected_baseline_sha256"],
        fixture_mode=fixture_mode,
    )
    spec["rebuild_source_spec_sha256"] = rebuild_source["rebuild_source_spec_sha256"]
    return {
        "spec": spec,
        "spec_sha256": sha256_bytes(canonical_json_bytes(spec)),
        **rebuild_source,
    }


def validate_prepared_adapter_provenance(
    project_root: Path,
    base_selection_path: Path,
    frozen_producers: list[Mapping[str, object]],
) -> dict[str, object]:
    """Bind prepared adapter outputs to inputs and the exact frozen producers.

    This runs before candidate-root creation.  A producer edit after adapter
    preparation therefore fails closed without leaving a partial read-only
    candidate behind.
    """

    selection_rows = read_tsv_exact(base_selection_path, BASE_SELECTION_FIELDS)
    artifacts: dict[str, Mapping[str, str]] = {}
    for artifact in selection_rows:
        artifact_id = artifact["artifact_id"]
        if artifact_id in artifacts:
            raise ReleaseContractError(
                f"duplicate BASE artifact ID in adapter provenance gate: {artifact_id}"
            )
        artifacts[artifact_id] = artifact
    required = EXPECTED_ADAPTER_IDS | {
        "adapter_provenance",
        "recursive_release_producers",
    }
    missing = sorted(required - set(artifacts))
    if missing:
        raise ReleaseContractError(
            f"prepared release lacks adapter-provenance artifacts: {missing}"
        )

    provenance_artifact = artifacts["adapter_provenance"]
    recursive_artifact = artifacts["recursive_release_producers"]
    provenance_source = resolve_project_path(
        project_root,
        provenance_artifact["source_path"],
        "adapter provenance source",
    )
    recursive_source = resolve_project_path(
        project_root,
        recursive_artifact["source_path"],
        "recursive producer source",
    )
    provenance_rows = read_tsv_exact(provenance_source, ADAPTER_PROVENANCE_FIELDS)
    recursive_rows = read_tsv_exact(recursive_source, RECURSIVE_PRODUCER_FIELDS)
    if {row["artifact_id"] for row in provenance_rows} != EXPECTED_ADAPTER_IDS:
        raise ReleaseContractError(
            "adapter provenance must contain exactly the seven presentation adapters"
        )
    if len(provenance_rows) != len(EXPECTED_ADAPTER_IDS):
        raise ReleaseContractError("adapter provenance contains duplicate artifact IDs")

    recursive_by_path: dict[str, Mapping[str, str]] = {}
    for row in recursive_rows:
        repository_path = row["repository_path"]
        if not repository_path or repository_path in recursive_by_path:
            raise ReleaseContractError(
                f"blank/duplicate recursive producer path: {repository_path!r}"
            )
        require_sha256(row["sha256"], f"recursive producer {repository_path}")
        parse_nonnegative_int(row["bytes"], f"recursive producer {repository_path}")
        recursive_by_path[repository_path] = row

    frozen_by_path: dict[str, Mapping[str, object]] = {}
    for row in frozen_producers:
        repository_path = str(row["repository_path"])
        if not repository_path or repository_path in frozen_by_path:
            raise ReleaseContractError(
                f"blank/duplicate REL01 producer path: {repository_path!r}"
            )
        frozen_by_path[repository_path] = row

    # The prepared recursive manifest is the complete temporal dependency
    # universe.  REL01 may freeze additional historical/portal producers, but
    # it must freeze every recursive producer with the identical stable ID and
    # exact bytes; otherwise a non-adapter helper could drift after preparation.
    missing_recursive = sorted(set(recursive_by_path) - set(frozen_by_path))
    if missing_recursive:
        raise ReleaseContractError(
            f"recursive producer absent from REL01 freeze: {missing_recursive}"
        )
    for repository_path, recursive in recursive_by_path.items():
        frozen = frozen_by_path[repository_path]
        if (
            str(frozen["producer_id"]) != recursive["producer_id"]
            or str(frozen["sha256"]) != recursive["sha256"]
            or parse_nonnegative_int(
                frozen["bytes"], f"REL01 producer {repository_path} bytes"
            )
            != parse_nonnegative_int(
                recursive["bytes"], f"recursive producer {repository_path} bytes"
            )
        ):
            raise ReleaseContractError(
                f"recursive producer differs from REL01 freeze: {repository_path}"
            )

    recursive_sha256 = sha256_file(recursive_source)
    for row in provenance_rows:
        artifact_id = row["artifact_id"]
        output = artifacts[artifact_id]
        output_source = resolve_project_path(
            project_root,
            output["source_path"],
            f"adapter output {artifact_id}",
        )
        if (
            row["contract_version"] != ADAPTER_PROVENANCE_VERSION
            or row["adapter_contract_sha256"] != ADAPTER_CONTRACT_SHA256
        ):
            raise ReleaseContractError(
                f"adapter provenance contract drift: {artifact_id}"
            )
        if row["canonical_promotion_authorized"].lower() != "false":
            raise ReleaseContractError(
                f"adapter provenance authorizes promotion: {artifact_id}"
            )
        if (
            row["output_relative_path"] != output_source.name
            or row["output_relative_path"]
            != ADAPTER_CONTRACT[artifact_id]["output_relative_path"]
        ):
            raise ReleaseContractError(
                f"adapter output path binding drift: {artifact_id}"
            )
        expected_output_sha = require_sha256(
            row["output_sha256"], f"adapter output {artifact_id}"
        )
        expected_output_bytes = parse_nonnegative_int(
            row["output_bytes"], f"adapter output {artifact_id} bytes"
        )
        if (
            output["source_sha256"] != expected_output_sha
            or parse_nonnegative_int(
                output["source_bytes"], f"BASE adapter {artifact_id} bytes"
            )
            != expected_output_bytes
            or sha256_file(output_source) != expected_output_sha
            or output_source.stat().st_size != expected_output_bytes
        ):
            raise ReleaseContractError(
                f"adapter output no longer matches provenance: {artifact_id}"
            )
        if (
            row["producer_manifest_relative_path"]
            != "inputs/BASE/recursive_release_producers.tsv"
            or row["producer_manifest_sha256"] != recursive_sha256
        ):
            raise ReleaseContractError(
                f"adapter producer-manifest binding drift: {artifact_id}"
            )
        try:
            producer_bindings = json.loads(row["producer_bindings_json"])
            input_bindings = json.loads(row["input_bindings_json"])
        except json.JSONDecodeError as error:
            raise ReleaseContractError(
                f"invalid adapter provenance JSON: {artifact_id}"
            ) from error
        if not isinstance(producer_bindings, list) or not producer_bindings:
            raise ReleaseContractError(
                f"adapter lacks producer bindings: {artifact_id}"
            )
        seen_producers: set[str] = set()
        for binding in producer_bindings:
            if not isinstance(binding, dict) or set(binding) != {
                "producer_id",
                "repository_path",
                "sha256",
                "bytes",
            }:
                raise ReleaseContractError(
                    f"invalid producer binding schema: {artifact_id}"
                )
            repository_path = str(binding["repository_path"])
            if repository_path in seen_producers:
                raise ReleaseContractError(
                    f"duplicate adapter producer binding: {artifact_id}:{repository_path}"
                )
            seen_producers.add(repository_path)
            recursive = recursive_by_path.get(repository_path)
            frozen = frozen_by_path.get(repository_path)
            if recursive is None or frozen is None:
                raise ReleaseContractError(
                    f"adapter producer absent from recursive/REL01 freeze: {repository_path}"
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
                or binding_sha != str(frozen["sha256"])
                or binding_bytes != int(frozen["bytes"])
            ):
                raise ReleaseContractError(
                    f"adapter producer hash differs from REL01 freeze: {repository_path}"
                )
        if not isinstance(input_bindings, list) or not input_bindings:
            raise ReleaseContractError(f"adapter lacks input bindings: {artifact_id}")
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
        seen_inputs: set[str] = set()
        for binding in input_bindings:
            if not isinstance(binding, dict) or set(binding) != {
                "source_path",
                "sha256",
                "bytes",
            }:
                raise ReleaseContractError(
                    f"invalid adapter input binding schema: {artifact_id}"
                )
            source_text = str(binding["source_path"])
            if source_text in seen_inputs:
                raise ReleaseContractError(
                    f"duplicate adapter input binding: {artifact_id}:{source_text}"
                )
            seen_inputs.add(source_text)
            source = resolve_project_path(
                project_root, source_text, f"adapter input {artifact_id}"
            )
            expected_sha = require_sha256(
                str(binding["sha256"]), f"adapter input {artifact_id}:{source_text}"
            )
            expected_bytes = parse_nonnegative_int(
                binding["bytes"], f"adapter input {artifact_id}:{source_text} bytes"
            )
            if (
                not source.is_file()
                or source.is_symlink()
                or sha256_file(source) != expected_sha
                or source.stat().st_size != expected_bytes
            ):
                raise ReleaseContractError(
                    f"adapter input changed after preparation: {artifact_id}:{source_text}"
                )
    return {
        "contract_version": ADAPTER_PROVENANCE_VERSION,
        "adapter_contract_sha256": ADAPTER_CONTRACT_SHA256,
        "n_adapters": len(provenance_rows),
        "base_input_selection_sha256": sha256_file(base_selection_path),
        "adapter_provenance_sha256": sha256_file(provenance_source),
        "recursive_producer_manifest_sha256": recursive_sha256,
    }


def ensure_parent_is_candidate(destination: Path, candidate: Path) -> None:
    candidate_lexical = Path(os.path.abspath(candidate))
    destination_lexical = Path(os.path.abspath(destination))
    try:
        relative_parent = destination_lexical.parent.relative_to(candidate_lexical)
    except ValueError as error:
        raise ReleaseContractError(
            f"snapshot destination escapes candidate root: {destination}"
        ) from error
    if candidate_lexical.is_symlink() or not candidate_lexical.is_dir():
        raise ReleaseContractError(
            f"snapshot candidate is missing or symlinked: {candidate_lexical}"
        )
    cursor = candidate_lexical
    for part in relative_parent.parts:
        cursor /= part
        if cursor.is_symlink():
            raise ReleaseContractError(
                f"snapshot destination has a symlinked parent: {cursor}"
            )
        if cursor.exists() and not cursor.is_dir():
            raise ReleaseContractError(
                f"snapshot destination parent is not a directory: {cursor}"
            )
        cursor.mkdir(exist_ok=True)


def canonical_tsv_bytes(
    rows: list[Mapping[str, object]], fields: tuple[str, ...]
) -> bytes:
    handle = io.StringIO(newline="")
    writer = csv.DictWriter(
        handle,
        fieldnames=list(fields),
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue().encode("utf-8")


def pretty_json_bytes(payload: object) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")


def materialize_immutable(
    path: Path,
    payload: bytes,
    candidate: Path,
    context: str,
) -> None:
    ensure_parent_is_candidate(path, candidate)
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_file() or path.read_bytes() != payload:
            raise ReleaseContractError(
                f"existing candidate {context} differs; refusing replacement: {path}"
            )
        return
    with tempfile.NamedTemporaryFile(
        mode="wb",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    try:
        try:
            os.link(temporary, path)
        except FileExistsError as error:
            raise ReleaseContractError(
                f"candidate {context} appeared concurrently; refusing overwrite: {path}"
            ) from error
    finally:
        temporary.unlink(missing_ok=True)
    if path.is_symlink() or not path.is_file() or path.read_bytes() != payload:
        raise ReleaseContractError(
            f"candidate {context} failed post-write integrity check: {path}"
        )


def materialize_immutable_tsv(
    path: Path,
    rows: list[Mapping[str, object]],
    fields: tuple[str, ...],
    candidate: Path,
    context: str,
) -> None:
    materialize_immutable(
        path,
        canonical_tsv_bytes(rows, fields),
        candidate,
        context,
    )


def materialize_immutable_json(
    path: Path,
    payload: object,
    candidate: Path,
    context: str,
) -> None:
    materialize_immutable(
        path,
        pretty_json_bytes(payload),
        candidate,
        context,
    )


def copy_verified(
    source: Path,
    destination: Path,
    expected_sha256: str,
    expected_bytes: int,
    candidate: Path,
) -> None:
    if not source.is_file() or source.is_symlink():
        raise ReleaseContractError(
            f"copy source must be a regular non-symlink file: {source}"
        )
    if (
        source.stat().st_size != expected_bytes
        or sha256_file(source) != expected_sha256
    ):
        raise ReleaseContractError(f"copy source changed before snapshot: {source}")
    ensure_parent_is_candidate(destination, candidate)
    if destination.exists() or destination.is_symlink():
        if destination.is_symlink() or not destination.is_file():
            raise ReleaseContractError(
                f"existing snapshot destination is not a regular file: {destination}"
            )
        if (
            destination.stat().st_size != expected_bytes
            or sha256_file(destination) != expected_sha256
        ):
            raise ReleaseContractError(
                f"existing snapshot differs from requested immutable input: {destination}"
            )
        if (
            source.stat().st_size != expected_bytes
            or sha256_file(source) != expected_sha256
        ):
            raise ReleaseContractError(f"copy source changed during snapshot: {source}")
        return
    with tempfile.NamedTemporaryFile(
        mode="wb",
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
    try:
        shutil.copyfile(source, temporary)
        if (
            temporary.stat().st_size != expected_bytes
            or sha256_file(temporary) != expected_sha256
        ):
            raise ReleaseContractError(
                f"snapshot copy failed post-copy integrity check: {source}"
            )
        if (
            source.stat().st_size != expected_bytes
            or sha256_file(source) != expected_sha256
        ):
            raise ReleaseContractError(f"copy source changed during snapshot: {source}")
        try:
            os.link(temporary, destination)
        except FileExistsError as error:
            raise ReleaseContractError(
                f"snapshot destination appeared concurrently; refusing overwrite: {destination}"
            ) from error
    finally:
        if temporary.exists():
            temporary.unlink()
    if destination.is_symlink() or sha256_file(destination) != expected_sha256:
        raise ReleaseContractError(
            f"snapshot destination failed final integrity check: {destination}"
        )


def snapshot_record(
    project_root: Path,
    source: Path,
    destination_relpath: str,
    workstream_id: str,
    artifact_id: str,
    artifact_role: str,
) -> dict[str, object]:
    digest = sha256_file(source)
    size = source.stat().st_size
    return {
        "workstream_id": workstream_id,
        "artifact_id": artifact_id,
        "artifact_role": artifact_role,
        "source_path_provenance": project_relative(project_root, source),
        "source_sha256": digest,
        "source_bytes": size,
        "snapshot_relpath": destination_relpath,
        "snapshot_sha256": digest,
        "snapshot_bytes": size,
    }


def load_existing_spec(path: Path) -> Mapping[str, object]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReleaseContractError(
            f"cannot read existing candidate spec: {path}"
        ) from error


def snapshot_candidate(
    project_root: Path,
    closure_path: Path,
    protected_scopes_path: Path,
    protected_baseline_path: Path,
    fixture_mode: bool = False,
) -> dict[str, object]:
    project = project_root.resolve()
    closure = validate_closure(project, closure_path)
    protected_scopes = resolve_project_path(
        project,
        str(protected_scopes_path),
        "protected scope registry",
    )
    protected_baseline = resolve_project_path(
        project,
        str(protected_baseline_path),
        "protected release baseline",
    )
    protected_before = verify_protected_baseline(
        project, protected_scopes, protected_baseline
    )
    base_selection_path = (
        closure.closure_path.resolve().parent / "base_input_selection.tsv"
    )
    if not base_selection_path.is_file() or base_selection_path.is_symlink():
        raise ReleaseContractError(
            f"prepared BASE selection is missing or unsafe: {base_selection_path}"
        )
    specification = build_snapshot_spec(
        project,
        closure,
        protected_scopes,
        protected_baseline,
        base_selection_path,
        fixture_mode,
    )
    plan50_producer_binding = None
    adapter_provenance = None
    if not fixture_mode:
        plan50_producer_binding = validate_plan50_frozen_producer_binding(
            closure,
            list(specification["spec"]["producer_scripts"]),
        )
        adapter_provenance = validate_prepared_adapter_provenance(
            project,
            base_selection_path,
            list(specification["spec"]["producer_scripts"]),
        )
    candidate = assert_candidate_root(project, candidate_root(project))
    specification_path = candidate / "manifests/snapshot_spec.json"
    if candidate.exists():
        if candidate.is_symlink() or not candidate.is_dir():
            raise ReleaseContractError(
                f"existing candidate is not a regular directory: {candidate}"
            )
        assert_no_symlinks(candidate)
        if not specification_path.is_file() or specification_path.is_symlink():
            raise ReleaseContractError(
                "existing candidate lacks an immutable snapshot spec; refusing reuse"
            )
        existing = load_existing_spec(specification_path)
        if (
            existing != specification
            or specification_path.read_bytes() != pretty_json_bytes(specification)
        ):
            raise ReleaseContractError(
                "existing candidate snapshot spec differs; refusing in-place replacement"
            )
    else:
        candidate.mkdir(parents=True)
        for relative in REQUIRED_CANDIDATE_DIRS:
            (candidate / relative).mkdir(parents=True, exist_ok=True)
        materialize_immutable_json(
            specification_path,
            specification,
            candidate,
            "snapshot specification",
        )
    for relative in REQUIRED_CANDIDATE_DIRS:
        directory = candidate / relative
        if not directory.is_dir() or directory.is_symlink():
            raise ReleaseContractError(
                f"candidate directory is missing or symlinked: {directory}"
            )
    assert_no_symlinks(candidate)

    records = []
    administrative = (
        (
            closure.closure_path,
            "manifests/workstream_closure.tsv",
            "REL00",
            "workstream_closure",
            "closure_matrix",
        ),
        (
            protected_scopes,
            "manifests/protected_scopes.tsv",
            "REL01",
            "protected_scopes",
            "protected_release_registry",
        ),
        (
            protected_baseline,
            "manifests/protected_release_baseline.tsv",
            "REL01",
            "protected_release_baseline",
            "protected_release_baseline",
        ),
    )
    for source, relative, workstream, artifact_id, role in administrative:
        digest = sha256_file(source)
        size = source.stat().st_size
        copy_verified(source, candidate / relative, digest, size, candidate)
        records.append(
            snapshot_record(project, source, relative, workstream, artifact_id, role)
        )

    for workstream in closure.workstreams:
        relative = f"manifests/upstream/{workstream.row['workstream_id']}.tsv"
        copy_verified(
            workstream.manifest_path,
            candidate / relative,
            workstream.manifest_sha256,
            workstream.manifest_path.stat().st_size,
            candidate,
        )
        records.append(
            snapshot_record(
                project,
                workstream.manifest_path,
                relative,
                workstream.row["workstream_id"],
                "terminal_artifact_manifest",
                "upstream_terminal_manifest",
            )
        )

    downstream_rows = []
    for artifact in sorted(closure.artifacts, key=lambda item: item.snapshot_relpath):
        destination = candidate / artifact.snapshot_relpath
        copy_verified(
            artifact.source_path,
            destination,
            artifact.source_sha256,
            artifact.source_bytes,
            candidate,
        )
        records.append(
            {
                "workstream_id": artifact.workstream_id,
                "artifact_id": artifact.artifact_id,
                "artifact_role": artifact.artifact_role,
                "source_path_provenance": project_relative(
                    project, artifact.source_path
                ),
                "source_sha256": artifact.source_sha256,
                "source_bytes": artifact.source_bytes,
                "snapshot_relpath": artifact.snapshot_relpath,
                "snapshot_sha256": sha256_file(destination),
                "snapshot_bytes": destination.stat().st_size,
            }
        )
        downstream_rows.append(
            {
                "consumer_id": f"{artifact.workstream_id}:{artifact.artifact_id}",
                "workstream_id": artifact.workstream_id,
                "artifact_id": artifact.artifact_id,
                "snapshot_path": artifact.downstream_read_path,
                "sha256": artifact.source_sha256,
                "bytes": artifact.source_bytes,
            }
        )

    snapshot_manifest = candidate / "manifests/input_snapshot_manifest.tsv"
    materialize_immutable_tsv(
        snapshot_manifest,
        sorted(records, key=lambda row: str(row["snapshot_relpath"])),
        SNAPSHOT_MANIFEST_FIELDS,
        candidate,
        "input snapshot manifest",
    )
    downstream_manifest = candidate / "manifests/downstream_inputs.tsv"
    materialize_immutable_tsv(
        downstream_manifest,
        sorted(downstream_rows, key=lambda row: str(row["consumer_id"])),
        DOWNSTREAM_FIELDS,
        candidate,
        "downstream input manifest",
    )
    repository_path = candidate / "manifests/repository_state.json"
    materialize_immutable_json(
        repository_path,
        specification["spec"]["repository"],
        candidate,
        "repository state",
    )
    rebuild_source_path = candidate / REBUILD_SOURCE_FILENAME
    materialize_immutable_json(
        rebuild_source_path,
        {
            "rebuild_source_spec": specification["rebuild_source_spec"],
            "rebuild_source_spec_sha256": specification["rebuild_source_spec_sha256"],
        },
        candidate,
        "rebuild source specification",
    )
    producer_path = candidate / "manifests/producer_script_manifest.tsv"
    for producer in specification["spec"]["producer_scripts"]:
        source = producer_source_path(producer)
        destination = candidate / str(producer["snapshot_path"])
        copy_verified(
            source,
            destination,
            str(producer["sha256"]),
            int(producer["bytes"]),
            candidate,
        )
    materialize_immutable_tsv(
        producer_path,
        specification["spec"]["producer_scripts"],
        PRODUCER_SCRIPT_FIELDS,
        candidate,
        "producer script manifest",
    )
    plan50_relocated_verification = None
    if not fixture_mode:
        plan50_relocated_verification = verify_relocated_plan50_bundle(candidate)
        materialize_immutable_json(
            candidate / PLAN50_RELOCATED_VERIFICATION_REL,
            plan50_relocated_verification,
            candidate,
            "relocated Plan 50 verification",
        )
    runtime_environment_path = candidate / "environments/runtime_environment.json"
    materialize_immutable_json(
        runtime_environment_path,
        specification["spec"]["runtime_environment"],
        candidate,
        "runtime environment",
    )
    python_packages_path = candidate / "environments/python_packages.tsv"
    materialize_immutable_tsv(
        python_packages_path,
        specification["spec"]["python_packages"],
        PYTHON_PACKAGE_FIELDS,
        candidate,
        "Python package export",
    )
    environment_index_path = candidate / "environments/environment_index.tsv"
    environment_rows = [
        {
            "environment_id": "release_runtime",
            "snapshot_path": "environments/runtime_environment.json",
            "sha256": sha256_file(runtime_environment_path),
            "bytes": runtime_environment_path.stat().st_size,
            "scope": "REL01-REL05 runtime",
        },
        {
            "environment_id": "release_python_packages",
            "snapshot_path": "environments/python_packages.tsv",
            "sha256": sha256_file(python_packages_path),
            "bytes": python_packages_path.stat().st_size,
            "scope": "installed Python distributions at snapshot time",
        },
    ]
    materialize_immutable_tsv(
        environment_index_path,
        environment_rows,
        ENVIRONMENT_INDEX_FIELDS,
        candidate,
        "environment index",
    )
    execution_path = candidate / "manifests/execution_context.json"
    if execution_path.exists() or execution_path.is_symlink():
        if execution_path.is_symlink() or not execution_path.is_file():
            raise ReleaseContractError(
                f"execution context is not a regular file: {execution_path}"
            )
        execution_context = load_existing_spec(execution_path)
        if execution_path.read_bytes() != pretty_json_bytes(execution_context):
            raise ReleaseContractError(
                "existing execution context is not in its immutable canonical form"
            )
    else:
        execution_context = {
            "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
            "command_line": shlex.join(sys.argv),
            "host": platform.node(),
            "process_id": os.getpid(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
            "slurm_step_id": os.environ.get("SLURM_STEP_ID", ""),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "safe_environment_exports": {
                key: os.environ[key]
                for key in SAFE_ENVIRONMENT_KEYS
                if key in os.environ
            },
            "fixture_mode": fixture_mode,
        }
        materialize_immutable_json(
            execution_path,
            execution_context,
            candidate,
            "execution context",
        )
    registry_path = candidate / "manifests/downstream_manifest_registry.tsv"
    registry_rows = [
        {
            "manifest_id": "downstream_inputs",
            "manifest_path": "manifests/downstream_inputs.tsv",
            "sha256": sha256_file(downstream_manifest),
            "bytes": downstream_manifest.stat().st_size,
            "consumer_scope": "candidate_only",
        }
    ]
    materialize_immutable_tsv(
        registry_path,
        registry_rows,
        MANIFEST_REGISTRY_FIELDS,
        candidate,
        "downstream manifest registry",
    )
    state = {
        "candidate_id": CANDIDATE_ID,
        "state": "SNAPSHOT_FROZEN",
        "spec_sha256": specification["spec_sha256"],
        "snapshot_manifest_sha256": sha256_file(snapshot_manifest),
        "downstream_manifest_sha256": sha256_file(downstream_manifest),
        "downstream_registry_sha256": sha256_file(registry_path),
        "repository_state_sha256": sha256_file(repository_path),
        "rebuild_source_spec_file_sha256": sha256_file(rebuild_source_path),
        "rebuild_source_spec_sha256": specification["rebuild_source_spec_sha256"],
        "producer_script_manifest_sha256": sha256_file(producer_path),
        "environment_index_sha256": sha256_file(environment_index_path),
        "execution_context_sha256": sha256_file(execution_path),
        "adapter_provenance": adapter_provenance,
        "plan50_producer_binding": plan50_producer_binding,
        "plan50_relocated_verification_sha256": (
            sha256_file(candidate / PLAN50_RELOCATED_VERIFICATION_REL)
            if plan50_relocated_verification is not None
            else None
        ),
        "scientific_assembly_performed": False,
        "canonical_promotion_status": "not_promoted",
    }
    materialize_immutable_json(
        candidate / "manifests/candidate_state.json",
        state,
        candidate,
        "candidate state",
    )
    assert_no_symlinks(candidate)
    protected_after = verify_protected_baseline(
        project,
        candidate / "manifests/protected_scopes.tsv",
        candidate / "manifests/protected_release_baseline.tsv",
    )
    if protected_before != protected_after:
        raise ReleaseContractError(
            "protected release verification changed during candidate snapshot"
        )
    return {
        "status": "REL01_SNAPSHOT_FROZEN",
        "candidate_id": CANDIDATE_ID,
        "candidate_root": str(candidate),
        "spec_sha256": specification["spec_sha256"],
        "rebuild_source_spec_sha256": specification["rebuild_source_spec_sha256"],
        "snapshot_manifest_sha256": state["snapshot_manifest_sha256"],
        "n_workstreams": len(closure.workstreams),
        "n_copied_artifacts": len(closure.artifacts),
        "n_protected_files": protected_after["n_files"],
        "scientific_assembly_performed": False,
        "canonical_promotion_status": "not_promoted",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--closure", type=Path, required=True)
    parser.add_argument("--protected-scopes", type=Path, required=True)
    parser.add_argument("--protected-baseline", type=Path, required=True)
    parser.add_argument(
        "--fixture-mode",
        action="store_true",
        help="Use only for isolated synthetic test roots; records synthetic Git state.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        result = snapshot_candidate(
            project_root=args.project_root,
            closure_path=args.closure,
            protected_scopes_path=args.protected_scopes,
            protected_baseline_path=args.protected_baseline,
            fixture_mode=args.fixture_mode,
        )
    except ReleaseContractError as error:
        raise SystemExit(f"REL01_BLOCKED: {error}") from error
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
