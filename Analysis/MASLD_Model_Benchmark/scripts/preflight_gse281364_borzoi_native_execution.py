#!/usr/bin/env python3
"""Audit the native Borzoi execution check without loading weights or outcomes."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Callable, Mapping

from masld_bench.artifacts import (
    ArtifactError,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)


CONFIG_SCHEMA = "masld-bench-gse281364-borzoi-native-execution-gate-v1"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
EXPECTED_GATE_ORDER = (
    "weight_terms",
    "official_checkpoints",
    "native_runtime",
    "reference_parity",
    "native_parity",
)


class BorzoiExecutionError(RuntimeError):
    """Raised when a fixed Borzoi execution authority differs."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise BorzoiExecutionError(f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise BorzoiExecutionError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise BorzoiExecutionError(f"{label} must be a JSON object")
    return value


def project_path(root: Path, relative: str, *, must_exist: bool) -> Path:
    configured = Path(relative)
    if configured.is_absolute() or ".." in configured.parts:
        raise BorzoiExecutionError(f"unsafe project-relative path: {relative}")
    path = reject_symlink_components(root / configured, label="Borzoi authority")
    if must_exist:
        try:
            path.resolve(strict=True).relative_to(root)
        except (OSError, ValueError) as error:
            raise BorzoiExecutionError(
                f"authority is missing or escapes project root: {relative}"
            ) from error
    return path


def validate_fixed_authorities(
    root: Path, config: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    raw = config.get("frozen_authorities")
    if not isinstance(raw, dict) or set(raw) != {
        "screen_contract",
        "checkpoint_contract",
        "exposure_audit",
        "development_crosswalk",
        "native_fixture",
        "admission_sources",
        "project_reference",
    }:
        raise BorzoiExecutionError("fixed authority census differs")
    records: dict[str, dict[str, Any]] = {}
    for name, binding in raw.items():
        if not isinstance(binding, dict):
            raise BorzoiExecutionError(f"invalid fixed authority: {name}")
        expected = str(binding.get("sha256", ""))
        if not SHA256.fullmatch(expected):
            raise BorzoiExecutionError(f"invalid fixed authority SHA-256: {name}")
        path = project_path(root, str(binding.get("path", "")), must_exist=True)
        kind = binding.get("kind")
        if kind == "file":
            if path.is_symlink() or not path.is_file():
                raise BorzoiExecutionError(f"fixed authority is not a file: {name}")
            observed = sha256_file(path)
        elif kind == "tree":
            try:
                verify_frozen_tree(path)
            except ArtifactError as error:
                raise BorzoiExecutionError(
                    f"invalid frozen authority tree {name}: {error}"
                ) from error
            observed = sha256_file(path / "ARTIFACTS.json")
        else:
            raise BorzoiExecutionError(f"unknown fixed authority kind: {name}")
        if observed != expected:
            raise BorzoiExecutionError(f"fixed authority changed: {name}")
        records[name] = {
            "kind": kind,
            "path": str(path.relative_to(root)),
            "sha256": observed,
        }
    return records


def validate_fixture_and_contracts(root: Path, config: Mapping[str, Any]) -> None:
    authorities = config["frozen_authorities"]
    screen = load_json(
        root / authorities["screen_contract"]["path"], label="screen contract"
    )
    checkpoint = load_json(
        root / authorities["checkpoint_contract"]["path"],
        label="checkpoint contract",
    )
    exposure = load_json(
        root / authorities["exposure_audit"]["path"], label="exposure audit"
    )
    crosswalk = load_json(
        root / authorities["development_crosswalk"]["path"],
        label="development crosswalk",
    )
    fixture = load_json(
        root / authorities["native_fixture"]["path"] / "fixture/receipt.json",
        label="fixture receipt",
    )
    release = config.get("official_release_contract", {})
    members = config.get("official_checkpoint_members")
    if (
        screen.get("schema_version")
        != "masld-bench-gse281364-borzoi-native-screen-v1"
        or screen.get("dataset_id") != "gse281364"
        or screen.get("model_id") != "borzoi_ensemble"
        or screen.get("execution_gate", {}).get("weight_terms")
        != "UNDECLARED_BLOCKED"
        or screen.get("execution_gate", {}).get("model_forward_allowed") is not False
        or screen.get("claim_boundary", {}).get("sealed_external_claim_allowed")
        is not False
    ):
        raise BorzoiExecutionError("screen execution firewall differs")
    if (
        fixture.get("schema_version")
        != "masld-bench-gse281364-borzoi-native-fixture-v1"
        or fixture.get("status") != "pass_fixture_checkpoint_execution_blocked"
        or fixture.get("fixture_elements") != 1_033
        or fixture.get("borzoi_long_range_groups") != 239
        or fixture.get("input_length_bp") != 524_288
        or fixture.get("native_human_tracks") != 7_611
        or fixture.get("official_ensemble_members") != 4
        or fixture.get("orientation_predictions_per_allele") != 8
        or fixture.get("checkpoint_bytes_loaded") is not False
        or fixture.get("model_forward_executed") is not False
        or fixture.get("outcomes_read") is not False
        or fixture.get("sealed_outcomes_read") is not False
    ):
        raise BorzoiExecutionError("outcome-blind fixture receipt differs")
    checkpoint_release = checkpoint.get("release", {})
    checkpoint_members = checkpoint.get("canonical_ensemble", {}).get("members")
    if (
        checkpoint_release.get("borzoi_revision")
        != release.get("borzoi_revision")
        or checkpoint_release.get("borzoi_paper_revision")
        != release.get("borzoi_paper_revision")
        or checkpoint.get("weight_license") != "UNDECLARED"
        or checkpoint.get("runtime_contract", {}).get(
            "checkpoint_bound_baskerville_revision"
        )
        != "UNRESOLVED"
        or not isinstance(checkpoint_members, list)
        or not isinstance(members, list)
        or len(checkpoint_members) != 4
        or len(members) != 4
    ):
        raise BorzoiExecutionError("official release or checkpoint contract differs")
    if (
        exposure.get("checkpoint_findings", {})
        .get("borzoi_ensemble", {})
        .get("exposure_state")
        != "target_label_unexposed"
        or exposure.get("license_disposition", {}).get("weights")
        != "UNDECLARED_on_the_official_GCS_objects"
        or crosswalk.get("findings", {}).get("gse281364", {}).get("exposure_state")
        != "clean_declared"
        or crosswalk.get("findings", {}).get("gse289173", {}).get("exposure_state")
        != "target_label_unexposed"
    ):
        raise BorzoiExecutionError("checkpoint contamination or license audit differs")
    for expected, observed in zip(members, checkpoint_members, strict=True):
        for configured_key, checkpoint_key in (
            ("replicate", "replicate"),
            ("name", "name"),
            ("generation", "gcs_generation"),
            ("size_bytes", "size_bytes"),
            ("md5", "md5"),
            ("crc32c_base64", "crc32c_base64"),
        ):
            if expected.get(configured_key) != observed.get(checkpoint_key):
                raise BorzoiExecutionError(
                    f"official checkpoint member differs: {expected.get('replicate')}"
                )


def observe_local_port(config: Mapping[str, Any]) -> dict[str, Any]:
    observation = config.get("local_port_observation")
    if not isinstance(observation, dict):
        raise BorzoiExecutionError("local port observation contract differs")
    environment = reject_symlink_components(
        Path(str(observation.get("environment", ""))), label="local port environment"
    )
    record: dict[str, Any] = {
        "environment": str(environment),
        "observed": environment.is_dir(),
        "modules_imported": False,
        "checkpoint_bytes_opened": False,
        "eligible_as_native_runtime": False,
        "disposition": observation.get("disposition"),
    }
    if not environment.is_dir():
        record["reason"] = "environment_absent"
        return record
    identities: dict[str, str] = {}
    for key, path_key, hash_key in (
        ("conda_history", "conda_history", "conda_history_sha256"),
        ("port_metadata", "port_metadata", "port_metadata_sha256"),
    ):
        relative = Path(str(observation.get(path_key, "")))
        if relative.is_absolute() or ".." in relative.parts:
            raise BorzoiExecutionError(f"unsafe local port path: {path_key}")
        path = reject_symlink_components(environment / relative, label="local port file")
        if not path.is_file() or path.is_symlink():
            raise BorzoiExecutionError(f"registered local port file missing: {path_key}")
        observed = sha256_file(path)
        if observed != observation.get(hash_key):
            raise BorzoiExecutionError(f"registered local port file changed: {path_key}")
        identities[key] = observed
    site_packages = environment / "lib/python3.11/site-packages"
    identities["tensorflow_distribution_count"] = str(
        len(list(site_packages.glob("tensorflow-*.dist-info")))
    )
    identities["baskerville_distribution_count"] = str(
        len(list(site_packages.glob("baskerville-*.dist-info")))
    )
    identities["westminster_distribution_count"] = str(
        len(list(site_packages.glob("westminster-*.dist-info")))
    )
    record["identities"] = identities
    record["reason"] = (
        "registered_borzoi_pytorch_0_5_1_python_3_11_port_lacks_"
        "native_TensorFlow_2_15_and_checkpoint_bound_Baskerville"
    )
    return record


def require_bool(receipt: Mapping[str, Any], key: str, value: bool) -> None:
    if receipt.get(key) is not value:
        raise BorzoiExecutionError(f"gate receipt field differs: {key}")


def validate_weight_terms(
    receipt: Mapping[str, Any], _: Mapping[str, Any]
) -> dict[str, Any]:
    if receipt.get("status") != "pass_authorized_for_research_execution":
        raise BorzoiExecutionError("weight execution authority is absent")
    for key in ("checkpoint_download_allowed", "research_inference_allowed"):
        require_bool(receipt, key, True)
    authority = receipt.get("determination_authority")
    if not isinstance(authority, str) or not authority.strip():
        raise BorzoiExecutionError("weight determination authority is absent")
    return {"determination_authority": authority, "status": receipt["status"]}


def validate_checkpoints(
    receipt: Mapping[str, Any], config: Mapping[str, Any]
) -> dict[str, Any]:
    if receipt.get("status") != "pass_strict_native_checkpoint_bundle":
        raise BorzoiExecutionError("native checkpoint bundle did not pass")
    require_bool(receipt, "external_or_soft_hdf5_links_present", False)
    require_bool(receipt, "strict_tensor_inventory_passed", True)
    observed = receipt.get("members")
    expected = config.get("official_checkpoint_members")
    if not isinstance(observed, list) or not isinstance(expected, list) or len(observed) != 4:
        raise BorzoiExecutionError("native checkpoint member census differs")
    member_files: list[dict[str, Any]] = []
    for wanted, actual in zip(expected, observed, strict=True):
        if not isinstance(actual, dict):
            raise BorzoiExecutionError("native checkpoint member record differs")
        for key in ("replicate", "name", "generation", "size_bytes", "md5", "crc32c_base64"):
            if actual.get(key) != wanted.get(key):
                raise BorzoiExecutionError(f"native checkpoint identity differs: {key}")
        checksum = str(actual.get("sha256", ""))
        if not SHA256.fullmatch(checksum):
            raise BorzoiExecutionError("native checkpoint SHA-256 is unresolved")
        relative = Path(str(actual.get("path", "")))
        if relative.is_absolute() or ".." in relative.parts or not relative.parts:
            raise BorzoiExecutionError("unsafe native checkpoint member path")
        member_files.append(
            {
                "replicate": wanted["replicate"],
                "path": relative.as_posix(),
                "sha256": checksum,
                "size_bytes": wanted["size_bytes"],
            }
        )
    return {"members": member_files, "status": receipt["status"]}


def validate_runtime(
    receipt: Mapping[str, Any], config: Mapping[str, Any]
) -> dict[str, Any]:
    release = config["official_release_contract"]
    if receipt.get("status") != "pass_offline_native_inference_runtime":
        raise BorzoiExecutionError("native runtime did not pass")
    if {
        "documented_release_python_series": release.get(
            "documented_release_python_series"
        ),
        "pyproject_python_requirement": release.get("pyproject_python_requirement"),
        "tensorflow_requirement": release.get("tensorflow_requirement"),
        "h5py_requirement": release.get("h5py_requirement"),
        "numpy_requirement": release.get("numpy_requirement"),
    } != {
        "documented_release_python_series": "3.10",
        "pyproject_python_requirement": ">=3.9",
        "tensorflow_requirement": "~=2.15.0",
        "h5py_requirement": "~=3.10.0",
        "numpy_requirement": "~=1.24.3",
    }:
        raise BorzoiExecutionError("official runtime requirement contract differs")
    required = {
        "borzoi_revision": release["borzoi_revision"],
        "westminster_required": False,
        "network_access": False,
    }
    for key, expected in required.items():
        if receipt.get(key) != expected:
            raise BorzoiExecutionError(f"native runtime field differs: {key}")
    version_requirements = {
        "python_version": (3, 10, 0),
        "tensorflow_version": (2, 15, 0),
        "h5py_version": (3, 10, 0),
        "numpy_version": (1, 24, 3),
    }
    versions: dict[str, str] = {}
    for key, expected in version_requirements.items():
        value = str(receipt.get(key, ""))
        match = re.fullmatch(r"(\d+)\.(\d+)\.(\d+)", value)
        if match is None:
            raise BorzoiExecutionError(f"native runtime version is not exact: {key}")
        observed = tuple(int(part) for part in match.groups())
        if observed[:2] != expected[:2] or observed[2] < expected[2]:
            raise BorzoiExecutionError(f"native runtime version differs: {key}")
        versions[key] = value
    baskerville = str(receipt.get("baskerville_revision", ""))
    if not GIT_SHA.fullmatch(baskerville):
        raise BorzoiExecutionError("checkpoint-bound Baskerville revision is unresolved")
    entrypoint = receipt.get("prediction_entrypoint")
    python = receipt.get("python_executable")
    bound_files: dict[str, dict[str, str]] = {}
    for name, value in (
        ("python_executable", python),
        ("prediction_entrypoint", entrypoint),
    ):
        if not isinstance(value, dict):
            raise BorzoiExecutionError(f"native {name} is absent")
        relative = Path(str(value.get("path", "")))
        checksum = str(value.get("sha256", ""))
        if relative.is_absolute() or ".." in relative.parts or not relative.parts:
            raise BorzoiExecutionError(f"unsafe native {name}")
        if not SHA256.fullmatch(checksum):
            raise BorzoiExecutionError(f"native {name} is not hash-bound")
        bound_files[name] = {"path": relative.as_posix(), "sha256": checksum}
    return {
        "status": receipt["status"],
        "baskerville_revision": baskerville,
        "versions": versions,
        **bound_files,
    }


def validate_reference(
    receipt: Mapping[str, Any], config: Mapping[str, Any]
) -> dict[str, Any]:
    if receipt.get("status") != "pass_primary_contig_inference_parity":
        raise BorzoiExecutionError("inference reference parity did not pass")
    require_bool(receipt, "primary_contig_sequences_identical", True)
    if receipt.get("primary_contig_differences") != []:
        raise BorzoiExecutionError("inference reference has primary-contig differences")
    if (
        receipt.get("benchmark_fasta_sha256")
        != config["reference_scope"]["benchmark_fasta_sha256"]
        or not SHA256.fullmatch(str(receipt.get("official_inference_fasta_sha256", "")))
    ):
        raise BorzoiExecutionError("inference reference identity differs")
    return {
        "status": receipt["status"],
        "official_inference_fasta_sha256": receipt["official_inference_fasta_sha256"],
        "training_reference_required_for_zero_shot": False,
    }


def validate_native_parity(
    receipt: Mapping[str, Any], config: Mapping[str, Any]
) -> dict[str, Any]:
    release = config["official_release_contract"]
    if receipt.get("status") != "pass_all_member_native_numeric_parity":
        raise BorzoiExecutionError("native numeric parity did not pass")
    for key in (
        "strict_restore_passed",
        "forward_passed",
        "reverse_complement_passed",
        "strand_pair_restore_passed",
        "inverse_transform_passed",
        "zero_shift_passed",
        "determinism_passed",
    ):
        require_bool(receipt, key, True)
    if (
        receipt.get("members_passed") != [0, 1, 2, 3]
        or receipt.get("human_tracks") != release["human_tracks"]
        or receipt.get("raw_output_bins") != release["raw_output_bins"]
        or receipt.get("released_prediction_bins")
        != release["released_prediction_bins"]
        or receipt.get("central_training_loss_bins")
        != release["central_training_loss_bins"]
    ):
        raise BorzoiExecutionError("native parity tensor contract differs")
    tolerance = receipt.get("numeric_tolerance")
    if not isinstance(tolerance, dict) or tolerance.get("passed") is not True:
        raise BorzoiExecutionError("native parity tolerance is absent")
    return {"status": receipt["status"], "members_passed": [0, 1, 2, 3]}


GATE_VALIDATORS: dict[
    str, Callable[[Mapping[str, Any], Mapping[str, Any]], dict[str, Any]]
] = {
    "weight_terms": validate_weight_terms,
    "official_checkpoints": validate_checkpoints,
    "native_runtime": validate_runtime,
    "reference_parity": validate_reference,
    "native_parity": validate_native_parity,
}


def validate_execution_gates(
    root: Path, config: Mapping[str, Any]
) -> tuple[dict[str, dict[str, Any]], list[str]]:
    gates = config.get("required_execution_artifacts")
    if not isinstance(gates, dict) or tuple(gates) != EXPECTED_GATE_ORDER:
        raise BorzoiExecutionError("execution gate order or census differs")
    results: dict[str, dict[str, Any]] = {}
    blockers: list[str] = []
    for name in EXPECTED_GATE_ORDER:
        gate = gates[name]
        if not isinstance(gate, dict):
            raise BorzoiExecutionError(f"invalid execution gate: {name}")
        path = project_path(root, str(gate.get("path", "")), must_exist=False)
        expected_hash = str(gate.get("artifacts_sha256", ""))
        record: dict[str, Any] = {
            "path": str(path.relative_to(root)),
            "expected_artifacts_sha256": expected_hash,
            "exists": path.is_dir(),
            "passed": False,
        }
        if expected_hash == "UNRESOLVED":
            blockers.append(f"artifact_hash_unresolved:{name}")
            if not path.is_dir():
                blockers.append(f"artifact_missing:{name}")
            record["status"] = "blocked_unbound"
            results[name] = record
            continue
        if not SHA256.fullmatch(expected_hash):
            raise BorzoiExecutionError(f"invalid execution artifact SHA-256: {name}")
        if not path.is_dir():
            blockers.append(f"artifact_missing:{name}")
            record["status"] = "blocked_missing"
            results[name] = record
            continue
        try:
            verify_frozen_tree(path)
            observed_hash = sha256_file(path / "ARTIFACTS.json")
            if observed_hash != expected_hash:
                raise BorzoiExecutionError("frozen ARTIFACTS.json hash differs")
            receipt = load_json(path / str(gate.get("receipt", "")), label=f"{name} receipt")
            if receipt.get("schema_version") != gate.get("receipt_schema"):
                raise BorzoiExecutionError("receipt schema differs")
            details = GATE_VALIDATORS[name](receipt, config)
            if name == "native_runtime":
                for file_name in ("python_executable", "prediction_entrypoint"):
                    bound = details[file_name]
                    bound_path = path / bound["path"]
                    if bound_path.is_symlink() or not bound_path.is_file():
                        raise BorzoiExecutionError(
                            f"bound native runtime file is absent: {file_name}"
                        )
                    if sha256_file(bound_path) != bound["sha256"]:
                        raise BorzoiExecutionError(
                            f"bound native runtime file changed: {file_name}"
                        )
            if name == "official_checkpoints":
                for member in details["members"]:
                    member_path = path / member["path"]
                    if member_path.is_symlink() or not member_path.is_file():
                        raise BorzoiExecutionError(
                            f"bound checkpoint member is absent: {member['replicate']}"
                        )
                    if (
                        member_path.stat().st_size != member["size_bytes"]
                        or sha256_file(member_path) != member["sha256"]
                    ):
                        raise BorzoiExecutionError(
                            f"bound checkpoint member changed: {member['replicate']}"
                        )
        except (ArtifactError, BorzoiExecutionError) as error:
            blockers.append(f"artifact_invalid:{name}:{error}")
            record["status"] = "blocked_invalid"
            results[name] = record
            continue
        record.update(
            {
                "status": "pass",
                "passed": True,
                "observed_artifacts_sha256": observed_hash,
                "details": details,
            }
        )
        results[name] = record
    return results, blockers


def preflight(*, root: Path, config_path: Path) -> dict[str, Any]:
    root = reject_symlink_components(root, label="project root").resolve(strict=True)
    config_path = reject_symlink_components(
        config_path, label="Borzoi execution config"
    ).resolve(strict=True)
    try:
        config_path.relative_to(root)
    except ValueError as error:
        raise BorzoiExecutionError("execution config escapes project root") from error
    config = load_json(config_path, label="Borzoi execution config")
    if (
        config.get("schema_version") != CONFIG_SCHEMA
        or config.get("dataset_id") != "gse281364"
        or config.get("model_id") != "borzoi_ensemble"
    ):
        raise BorzoiExecutionError("Borzoi execution config identity differs")
    fixed = validate_fixed_authorities(root, config)
    validate_fixture_and_contracts(root, config)
    gates, blockers = validate_execution_gates(root, config)
    policy = config.get("policy_transition")
    if (
        not isinstance(policy, dict)
        or policy.get("checkpoint_acquisition_authorized_by_this_contract") is not False
        or policy.get("dependency_install_authorized_by_this_contract") is not False
    ):
        raise BorzoiExecutionError("execution policy transition contract differs")
    if policy.get("current_model_forward_authorized") is not True:
        blockers.append("policy_transition_closed:model_forward")
    local_port = observe_local_port(config)
    executable = not blockers and all(record["passed"] for record in gates.values())
    runtime_details = gates.get("native_runtime", {}).get("details", {})
    return {
        "schema_version": "masld-bench-gse281364-borzoi-native-execution-preflight-v1",
        "status": "pass_executable" if executable else config["status"],
        "dataset_id": "gse281364",
        "model_id": "borzoi_ensemble",
        "config_sha256": sha256_file(config_path),
        "executable": executable,
        "fixed_authorities": fixed,
        "execution_gates": gates,
        "blockers": blockers,
        "local_port_observation": local_port,
        "prediction_entrypoint": runtime_details.get("prediction_entrypoint"),
        "checkpoint_bytes_opened": False,
        "model_modules_imported": False,
        "model_forward_executed": False,
        "reporter_counts_read": False,
        "reporter_outcomes_read": False,
        "sealed_assets_read": False,
        "metrics_calculated": False,
        "claim_role": config["claim_boundary"]["role"],
        "open_champion_eligible": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--mode", choices=("audit", "require-executable"), default="audit"
    )
    arguments = parser.parse_args()
    if arguments.output.exists() or arguments.output.is_symlink():
        raise BorzoiExecutionError(f"refusing to overwrite {arguments.output}")
    receipt = preflight(root=arguments.root, config_path=arguments.config)
    write_json_exclusive(arguments.output, receipt, mode=0o640)
    print(json.dumps(receipt, sort_keys=True))
    if arguments.mode == "require-executable" and not receipt["executable"]:
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
