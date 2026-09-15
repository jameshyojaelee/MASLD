#!/usr/bin/env python3
"""Admit AlphaGenome v0.8 SDK semantics without contacting or running the model."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import sys
from typing import Any, Mapping
from unittest import mock
import zipfile

from masld_bench.artifacts import (
    ArtifactError,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)


SCHEMA = "masld-bench-alphagenome-sdk-native-task-gate-v1"
RECEIPT_SCHEMA = "masld-bench-alphagenome-sdk-native-task-receipt-v1"
OUTPUT_ENUM = (
    "ATAC",
    "CAGE",
    "DNASE",
    "RNA_SEQ",
    "CHIP_HISTONE",
    "CHIP_TF",
    "SPLICE_SITES",
    "SPLICE_SITE_USAGE",
    "SPLICE_JUNCTIONS",
    "CONTACT_MAPS",
    "PROCAP",
)


class AlphaGenomeSdkGateError(RuntimeError):
    """Raised when an official source or the closed execution policy differs."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise AlphaGenomeSdkGateError(f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AlphaGenomeSdkGateError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise AlphaGenomeSdkGateError(f"{label} must be a JSON object")
    return value


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AlphaGenomeSdkGateError(message)


def message_body(proto: str, name: str) -> str:
    match = re.search(rf"^message {re.escape(name)} \{{(.*?)^\}}", proto, re.M | re.S)
    if match is None:
        raise AlphaGenomeSdkGateError(f"protobuf message is absent: {name}")
    return match.group(1)


def enum_members(proto: str, name: str, *, prefix: str) -> tuple[str, ...]:
    match = re.search(rf"^enum {re.escape(name)} \{{(.*?)^\}}", proto, re.M | re.S)
    if match is None:
        raise AlphaGenomeSdkGateError(f"protobuf enum is absent: {name}")
    return tuple(
        member.removeprefix(prefix)
        for member in re.findall(r"^\s*([A-Z][A-Z0-9_]*)\s*=\s*\d+\s*;", match.group(1), re.M)
        if member != f"{prefix}UNSPECIFIED"
    )


def validate_prior(root: Path, config: Mapping[str, Any]) -> dict[str, str]:
    prior = config["prior_preflight"]
    config_path = root / prior["config_path"]
    artifact_path = root / prior["artifact_path"]
    require(sha256_file(config_path) == prior["config_sha256"], "prior preflight config changed")
    try:
        verify_frozen_tree(artifact_path)
    except ArtifactError as error:
        raise AlphaGenomeSdkGateError(f"prior preflight artifact differs: {error}") from error
    require(
        sha256_file(artifact_path / "ARTIFACTS.json") == prior["artifacts_sha256"],
        "prior preflight artifact identity changed",
    )
    receipt = load_json(artifact_path / "preflight/receipt.json", label="prior receipt")
    require(not receipt["executable"] and not receipt["open_champion_eligible"], "prior firewall opened")
    require(not receipt["api_route"]["api_connection_attempted"], "prior API connection state changed")
    return {
        "config_path": prior["config_path"],
        "config_sha256": prior["config_sha256"],
        "artifact_path": prior["artifact_path"],
        "artifacts_sha256": prior["artifacts_sha256"],
    }


def validate_source_files(source_dir: Path, config: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for relative, expected in config["official_source_files"].items():
        path = source_dir / "github" / relative
        require(path.is_file() and not path.is_symlink(), f"official source absent: {relative}")
        require(path.stat().st_size == expected["size_bytes"], f"official source size differs: {relative}")
        digest = sha256_file(path)
        require(digest == expected["sha256"], f"official source hash differs: {relative}")
        records[relative] = {"sha256": digest, "size_bytes": path.stat().st_size}
    return records


def validate_remote_metadata(source_dir: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    release = load_json(source_dir / "github_release_latest.json", label="GitHub release metadata")
    tag = load_json(source_dir / "github_tag_v0.8.0.json", label="GitHub tag metadata")
    pypi = load_json(source_dir / "pypi_0.8.0.json", label="PyPI release metadata")
    hf = load_json(source_dir / "hf_alphagenome_all_folds.json", label="Hugging Face metadata")
    official = config["official_release"]
    package = config["official_package"]
    checkpoint = config["official_checkpoint"]
    require(
        release.get("tag_name") == official["tag"]
        and release.get("target_commitish") == official["revision"]
        and release.get("published_at") == official["published_at"],
        "latest official GitHub release differs",
    )
    tag_object = tag.get("object") or {}
    require(tag_object.get("type") == "commit" and tag_object.get("sha") == official["revision"], "official tag differs")
    require(
        (pypi.get("info") or {}).get("version") == package["version"]
        and (pypi.get("info") or {}).get("requires_python") == package["requires_python"],
        "PyPI version or Python requirement differs",
    )
    files = {row.get("filename"): row for row in pypi.get("urls", []) if isinstance(row, dict)}
    for kind in ("wheel", "sdist"):
        filename = package[f"{kind}_filename"]
        row = files.get(filename)
        require(row is not None, f"official {kind} is absent")
        require(
            row.get("size") == package[f"{kind}_size_bytes"]
            and (row.get("digests") or {}).get("sha256") == package[f"{kind}_sha256"],
            f"official {kind} identity differs",
        )
        require(
            str(row.get("url", "")).startswith("https://files.pythonhosted.org/"),
            f"official {kind} host differs",
        )
    card = hf.get("cardData") or {}
    require(
        hf.get("sha") == checkpoint["revision"]
        and hf.get("lastModified") == checkpoint["last_modified"]
        and hf.get("gated") == checkpoint["gated"]
        and hf.get("private") is checkpoint["private"],
        "official checkpoint metadata differs",
    )
    require(
        card.get("license_name") == checkpoint["license_name"]
        and card.get("license_link") == checkpoint["model_terms_url"]
        and "non-commercial use only" in str(card.get("extra_gated_prompt", "")),
        "official checkpoint terms marker differs",
    )
    model_card = source_dir / "hf_model_card.md"
    require(
        model_card.is_file()
        and model_card.stat().st_size == checkpoint["model_card_size_bytes"]
        and sha256_file(model_card) == checkpoint["model_card_sha256"],
        "official checkpoint model card differs",
    )
    return {
        "github_release_tag": release["tag_name"],
        "github_release_revision": release["target_commitish"],
        "github_release_published_at": release["published_at"],
        "pypi_version": pypi["info"]["version"],
        "pypi_wheel_sha256": package["wheel_sha256"],
        "pypi_sdist_sha256": package["sdist_sha256"],
        "huggingface_checkpoint_revision": hf["sha"],
        "huggingface_gate": hf["gated"],
        "checkpoint_model_terms_url": card["license_link"],
    }


def safe_extract_wheel(wheel: Path, destination: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    package = config["official_package"]
    require(wheel.is_file() and not wheel.is_symlink(), "official wheel is absent")
    require(wheel.stat().st_size == package["wheel_size_bytes"], "wheel size differs")
    require(sha256_file(wheel) == package["wheel_sha256"], "wheel hash differs")
    destination.mkdir(mode=0o750)
    members = 0
    with zipfile.ZipFile(wheel) as archive:
        for info in archive.infolist():
            relative = PurePosixPath(info.filename)
            mode = info.external_attr >> 16
            require(not relative.is_absolute() and ".." not in relative.parts, "unsafe wheel member path")
            require((mode & 0o170000) != 0o120000, "wheel symlink is forbidden")
            target = destination.joinpath(*relative.parts)
            if info.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(info) as source, target.open("xb") as sink:
                shutil.copyfileobj(source, sink)
            members += 1
    require(members == 71, "wheel member census differs")
    return {"wheel_sha256": package["wheel_sha256"], "regular_file_members": members}


def validate_protocol(source_dir: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    model_proto = (source_dir / "github/src/alphagenome/protos/dna_model.proto").read_text(encoding="utf-8")
    service_proto = (source_dir / "github/src/alphagenome/protos/dna_model_service.proto").read_text(encoding="utf-8")
    observed_outputs = enum_members(model_proto, "OutputType", prefix="OUTPUT_TYPE_")
    require(observed_outputs == tuple(config["native_interface"]["output_types"]), "native output enum differs")
    request_names = config["native_interface"]["request_methods_with_model_version"]
    for name in request_names:
        require(re.search(r"\bstring model_version\s*=", message_body(service_proto, name)) is not None, f"model version is absent from {name}")
    metadata_body = message_body(service_proto, "MetadataRequest")
    require("model_version" not in metadata_body, "metadata request unexpectedly became model-version-specific")
    response_names = (
        "PredictSequenceResponse",
        "PredictIntervalResponse",
        "PredictVariantResponse",
        "ScoreIntervalResponse",
        "ScoreVariantResponse",
        "ScoreIsmVariantResponse",
        "MetadataResponse",
    )
    identity_pattern = re.compile(r"\b(model_version|checkpoint_digest|server_build|model_digest)\b")
    require(all(identity_pattern.search(message_body(service_proto, name)) is None for name in response_names), "response now exposes immutable model identity")
    return {
        "output_types": list(observed_outputs),
        "requests_with_logical_model_version": list(request_names),
        "metadata_request_is_model_version_specific": False,
        "responses_echo_model_version_or_checkpoint_digest": False,
    }


def offline_runtime_probe(sdk_root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    sys.path.insert(0, str(sdk_root))
    for name in tuple(sys.modules):
        if name == "alphagenome" or name.startswith("alphagenome."):
            del sys.modules[name]
    alphagenome = importlib.import_module("alphagenome")
    dna_client = importlib.import_module("alphagenome.models.dna_client")
    dna_model = importlib.import_module("alphagenome.models.dna_model")
    dna_output = importlib.import_module("alphagenome.models.dna_output")
    genome = importlib.import_module("alphagenome.data.genome")
    require(alphagenome.__version__ == config["official_package"]["version"], "runtime SDK version differs")
    versions = tuple(item.name for item in dna_model.ModelVersion)
    outputs = tuple(item.name for item in dna_output.OutputType)
    lengths = tuple(sorted(dna_client.SUPPORTED_SEQUENCE_LENGTHS.values()))
    require(versions == tuple(config["native_interface"]["model_versions"]), "runtime model versions differ")
    require(outputs == tuple(config["native_interface"]["output_types"]), "runtime output types differ")
    require(lengths == tuple(config["native_interface"]["supported_sequence_lengths_bp"]), "runtime sequence lengths differ")

    captured: dict[str, Any] = {}

    class FakeStub:
        def __init__(self, channel: object):
            require(channel is sentinel_channel, "offline fake channel differs")

        def PredictSequence(self, requests: object, *, metadata: object) -> object:
            captured["prediction"] = next(iter(requests))
            captured["prediction_metadata"] = metadata
            return iter(())

        def GetMetadata(self, request: object, *, metadata: object) -> object:
            captured["metadata_request"] = request
            captured["metadata_metadata"] = metadata
            return iter(())

    sentinel_channel = object()
    sentinel_output = object()
    client = dna_client.DnaClient(channel=sentinel_channel, model_version=dna_model.ModelVersion.ALL_FOLDS, metadata=())
    with (
        mock.patch.object(dna_client.dna_model_service_pb2_grpc, "DnaModelServiceStub", FakeStub),
        mock.patch.object(dna_client, "_make_output", return_value=sentinel_output),
        mock.patch.object(dna_client, "construct_output_metadata", return_value=sentinel_output),
    ):
        result = client.predict_sequence(
            "N" * min(lengths),
            requested_outputs=[dna_output.OutputType.RNA_SEQ],
            ontology_terms=["UBERON:0002107"],
        )
        metadata_result = client.output_metadata()
    require(result is sentinel_output and metadata_result is sentinel_output, "offline stub result differs")
    request = captured["prediction"]
    require(request.model_version == "ALL_FOLDS", "ALL_FOLDS was not serialized")
    require(len(request.requested_outputs) == 1, "requested output serialization differs")
    metadata_fields = {field.name for field in captured["metadata_request"].DESCRIPTOR.fields}
    require("model_version" not in metadata_fields, "metadata runtime request gained model version")
    require(captured["prediction_metadata"] == () and captured["metadata_metadata"] == (), "credential metadata unexpectedly present")
    return {
        "sdk_version": alphagenome.__version__,
        "model_versions": list(versions),
        "supported_sequence_lengths_bp": list(lengths),
        "output_types": list(outputs),
        "logical_model_version_serialized": request.model_version,
        "credential_metadata_present": False,
        "grpc_channel_was_fake": True,
        "api_connection_attempted": False,
    }


def audit(root: Path, config_path: Path, source_dir: Path, output: Path) -> dict[str, Any]:
    root = reject_symlink_components(root, label="benchmark root").resolve(strict=True)
    config_path = reject_symlink_components(config_path, label="gate config").resolve(strict=True)
    source_dir = reject_symlink_components(source_dir, label="source directory").resolve(strict=True)
    require(root in config_path.parents and root in source_dir.parents, "gate input escapes benchmark root")
    require(not output.exists() and not output.is_symlink(), "refusing to overwrite gate output")
    config = load_json(config_path, label="gate config")
    require(config.get("schema_version") == SCHEMA, "gate schema differs")
    disposition = config["disposition"]
    require(
        disposition == {
            "sdk_source_and_offline_runtime_probe_allowed": True,
            "api_key_inspection_allowed": False,
            "api_connection_allowed": False,
            "checkpoint_download_allowed": False,
            "model_forward_allowed": False,
            "outcome_access_allowed": False,
            "sealed_access_allowed": False,
            "open_champion_eligible": False,
            "conditional_open_model_training_allowed": False,
            "terminal_after_gate": "sdk_native_interface_admitted_but_model_execution_remains_blocked",
        },
        "closed execution disposition differs",
    )
    prior = validate_prior(root, config)
    source_records = validate_source_files(source_dir, config)
    remote = validate_remote_metadata(source_dir, config)
    protocol = validate_protocol(source_dir, config)
    output.mkdir(parents=True, mode=0o750)
    wheel = source_dir / config["official_package"]["wheel_filename"]
    wheel_receipt = safe_extract_wheel(wheel, output / "sdk", config)
    runtime = offline_runtime_probe(output / "sdk", config)
    receipt = {
        "schema_version": RECEIPT_SCHEMA,
        "status": "pass_sdk_native_interface_admitted_model_execution_blocked",
        "config_sha256": sha256_file(config_path),
        "prior_preflight": prior,
        "official_source_records": source_records,
        "official_remote_metadata": remote,
        "wheel_inventory": wheel_receipt,
        "native_protocol": protocol,
        "offline_runtime_probe": runtime,
        "task_scope": config["task_scope"],
        "sdk_source_admitted": True,
        "native_interface_admitted": True,
        "api_key_inspected": False,
        "api_connection_attempted": False,
        "checkpoint_bytes_downloaded": False,
        "checkpoint_deserialized": False,
        "model_forward_executed": False,
        "observed_outcomes_loaded": False,
        "sealed_assets_loaded": False,
        "metrics_calculated": False,
        "open_champion_eligible": False,
        "conditional_open_model_training_allowed": False,
        "remaining_blockers": [
            "gated_noncommercial_model_terms",
            "no_authorized_checkpoint_or_api_execution_artifact",
            "api_responses_do_not_echo_an_immutable_server_build_identity",
            "metadata_request_is_not_model_version_specific",
            "no_frozen_numeric_repeatability_or_local_API_parity_fixture",
        ],
        "terminal_disposition": disposition["terminal_after_gate"],
    }
    write_json_exclusive(output / "receipt.json", receipt, mode=0o640)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = audit(arguments.root, arguments.config, arguments.source_dir, arguments.output)
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
