#!/usr/bin/env python3
"""Validate small official Borzoi/Enformer sources without loading checkpoint bytes."""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
from hashlib import md5, sha256
import json
from pathlib import Path
from typing import IO, Mapping
import tomllib


class AdmissionSourceError(ValueError):
    """Raised when an admission source differs from its immutable contract."""


def _hash(path: Path, algorithm: str = "sha256") -> str:
    value = sha256() if algorithm == "sha256" else md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _validate_file(path: Path, expected: Mapping[str, object]) -> dict[str, object]:
    if path.is_symlink() or not path.is_file():
        raise AdmissionSourceError(f"source is not a regular file: {path.name}")
    size = path.stat().st_size
    digest = _hash(path)
    if size != expected["size_bytes"] or digest != expected["sha256"]:
        raise AdmissionSourceError(f"source bytes differ: {path.name}")
    return {"path": path.name, "size_bytes": size, "sha256": digest}


def _target_rows(
    path: Path, compressed: bool
) -> tuple[list[str], list[dict[str, str]]]:
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        raw_fields = list(reader.fieldnames or [])
        fields = ["index" if field == "" else field for field in raw_fields]
        rows = []
        for row in reader:
            normalized = {
                ("index" if key == "" else str(key)): str(value)
                for key, value in row.items()
            }
            rows.append(normalized)
    return fields, rows


def _validate_targets(
    path: Path, expected: Mapping[str, object], *, borzoi: bool
) -> dict[str, object]:
    _validate_file(path, expected)
    fields, rows = _target_rows(path, compressed=borzoi)
    if fields != expected["fields"] or len(rows) != expected["rows"]:
        raise AdmissionSourceError(f"target manifest schema differs: {path.name}")
    try:
        indexes = [int(row["index"]) for row in rows]
    except ValueError as exc:
        raise AdmissionSourceError("target indexes differ") from exc
    if indexes != list(range(len(rows))):
        raise AdmissionSourceError("target indexes are not contiguous and ordered")
    prefixes: dict[str, int] = {}
    for row in rows:
        prefix = row["description"].split(":", 1)[0].upper()
        prefixes[prefix] = prefixes.get(prefix, 0) + 1
    payload: dict[str, object] = {
        "path": path.name,
        "rows": len(rows),
        "fields": fields,
        "description_prefix_counts": dict(sorted(prefixes.items())),
        "sha256": _hash(path),
    }
    if borzoi:
        try:
            strand_pair = [int(row["strand_pair"]) for row in rows]
        except ValueError as exc:
            raise AdmissionSourceError("Borzoi strand-pair mapping differs") from exc
        if any(value < 0 or value >= len(rows) for value in strand_pair) or any(
            strand_pair[strand_pair[index]] != index for index in range(len(rows))
        ):
            raise AdmissionSourceError("Borzoi strand-pair mapping is not an involution")
        expected_counts = {
            "ATAC": 232,
            "CAGE": 1276,
            "CHIP": 3886,
            "DNASE": 674,
            "RNA": 1543,
        }
        if prefixes != expected_counts:
            raise AdmissionSourceError("Borzoi assay target counts differ")
        payload["strand_pair_is_involution"] = True
        payload["reverse_complement_target_permutation_sha256"] = sha256(
            ("\n".join(str(value) for value in strand_pair) + "\n").encode("ascii")
        ).hexdigest()
    else:
        if "strand_pair" in fields:
            raise AdmissionSourceError("Enformer target manifest unexpectedly has strand pairs")
        payload["strand_pair_field_present"] = False
        payload["reverse_complement_target_axis"] = "identity_if_RC_adaptation_is_run"
    return payload


def _validate_metadata(path: Path, expected: Mapping[str, object]) -> dict[str, object]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 100_000:
        raise AdmissionSourceError("checkpoint metadata file differs")
    raw = json.loads(path.read_text(encoding="utf-8"))
    decoded_md5 = base64.b64decode(raw.get("md5Hash", "")).hex()
    if (
        raw.get("name") != expected["name"]
        or raw.get("generation") != expected["generation"]
        or int(raw.get("size", -1)) != expected["size_bytes"]
        or raw.get("crc32c") != expected["crc32c_base64"]
        or decoded_md5 != expected["md5_hex"]
    ):
        raise AdmissionSourceError(f"checkpoint object metadata differs: {path.name}")
    custom = raw.get("metadata") or {}
    license_fields = sorted(
        key for key in set(raw) | set(custom) if "licen" in key.lower()
    )
    if license_fields:
        raise AdmissionSourceError(
            "checkpoint metadata now declares terms; reassess before admission"
        )
    return {
        "name": raw["name"],
        "generation": raw["generation"],
        "metageneration": raw.get("metageneration"),
        "size_bytes": int(raw["size"]),
        "md5_hex": decoded_md5,
        "crc32c_base64": raw["crc32c"],
        "license_metadata_present": False,
        "checkpoint_bytes_downloaded": False,
    }


def _validate_audits(
    project_root: Path, config: Mapping[str, object]
) -> dict[str, object]:
    receipts = {}
    for key, binding in config["audit_bindings"].items():
        path = (project_root / binding["path"]).resolve(strict=True)
        if project_root.resolve() not in path.parents or path.is_symlink():
            raise AdmissionSourceError("audit binding escapes the benchmark root")
        digest = _hash(path)
        if digest != binding["sha256"]:
            raise AdmissionSourceError(f"audit binding differs: {key}")
        receipts[key] = {"path": binding["path"], "sha256": digest}
    borzoi = json.loads(
        (project_root / config["audit_bindings"]["borzoi_exposure_audit"]["path"])
        .read_text(encoding="utf-8")
    )
    enformer = json.loads(
        (project_root / config["audit_bindings"]["enformer_exposure_audit"]["path"])
        .read_text(encoding="utf-8")
    )
    registry = tomllib.loads(
        (project_root / config["audit_bindings"]["model_registry"]["path"])
        .read_text(encoding="utf-8")
    )
    registry_by_id = {row["model_id"]: row for row in registry["models"]}
    if (
        borzoi["checkpoint_findings"]["borzoi_ensemble"]["exposure_state"]
        != "target_label_unexposed"
        or "UNDECLARED" not in borzoi["license_disposition"]["weights"]
        or enformer["checkpoint_findings"]["official_native_sonnet_human"][
            "exposure_state"
        ]
        != "target_label_unexposed"
        or "UNDECLARED" not in enformer["license_disposition"]["native_weights"]
        or registry_by_id["borzoi_ensemble"]["status"] != "blocked_terms"
        or registry_by_id["enformer"]["status"] != "restricted_comparator"
        or registry_by_id["borzoi_ensemble"]["admission_blocking"] is not True
        or registry_by_id["enformer"]["admission_blocking"] is not True
    ):
        raise AdmissionSourceError("source/license/exposure disposition differs")
    return {
        "bindings": receipts,
        "gse289173_exposure": {
            "borzoi_ensemble": "target_label_unexposed",
            "enformer": "target_label_unexposed",
        },
        "weight_terms": {
            "borzoi_ensemble": "UNDECLARED_BLOCKED",
            "enformer_native": "UNDECLARED_BLOCKED",
        },
        "registry_status": {
            "borzoi_ensemble": "blocked_terms",
            "enformer": "restricted_comparator",
        },
        "review_passed_before_model_execution": True,
    }


def validate(
    config_path: Path, source_dir: Path, project_root: Path, output: Path
) -> dict[str, object]:
    if output.exists() or config_path.is_symlink() or source_dir.is_symlink():
        raise AdmissionSourceError("source validation request differs")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("schema_version") != "masld-bench-borzoi-enformer-admission-v1":
        raise AdmissionSourceError("admission config schema differs")
    if config["admission_disposition"] != {
        "source_and_fixture_probe_allowed": True,
        "checkpoint_download_allowed": False,
        "checkpoint_deserialization_allowed": False,
        "gpu_forward_allowed": False,
        "reason": "official_weight_terms_and_exact_no_network_native_runtimes_remain_unresolved",
        "open_champion_eligible": False,
    }:
        raise AdmissionSourceError("admission disposition differs")
    audit_receipt = _validate_audits(project_root, config)

    small_sources = {}
    for key, expected in config["source_files"].items():
        small_sources[key] = _validate_file(source_dir / key, expected)
    borzoi_targets = _validate_targets(
        source_dir / "borzoi_targets_human.txt.gz",
        config["models"]["borzoi_ensemble"]["target_manifest"],
        borzoi=True,
    )
    enformer_targets = _validate_targets(
        source_dir / "enformer_targets_human.txt",
        config["models"]["enformer"]["target_manifest"],
        borzoi=False,
    )
    parameter_expected = config["models"]["borzoi_ensemble"][
        "checkpoint_parameter_object"
    ]
    parameter_receipt = _validate_file(
        source_dir / "borzoi_checkpoint_params.json", parameter_expected
    )
    parameters = json.loads(
        (source_dir / "borzoi_checkpoint_params.json").read_text(encoding="utf-8")
    )
    trunk = parameters.get("model", {}).get("trunk", [])
    crops = [layer.get("cropping") for layer in trunk if layer.get("name") == "Cropping1D"]
    if (
        parameters.get("model", {}).get("seq_length") != 524288
        or parameters.get("model", {}).get("head_human", {}).get("units") != 7611
        or crops != [5120]
    ):
        raise AdmissionSourceError("Borzoi checkpoint parameter semantics differ")
    parameter_receipt["central_training_crop_bins_each_side"] = 5120

    checkpoint_metadata = {"borzoi_ensemble": [], "enformer": []}
    for index, expected in enumerate(
        config["models"]["borzoi_ensemble"]["checkpoint_members"]
    ):
        checkpoint_metadata["borzoi_ensemble"].append(
            _validate_metadata(source_dir / f"borzoi_checkpoint_{index}.metadata.json", expected)
        )
    for index, expected in enumerate(config["models"]["enformer"]["checkpoint_objects"]):
        checkpoint_metadata["enformer"].append(
            _validate_metadata(source_dir / f"enformer_checkpoint_{index}.metadata.json", expected)
        )

    if any(path.stat().st_size > 1_100_000 for path in source_dir.iterdir() if path.is_file()):
        raise AdmissionSourceError("unexpected checkpoint-sized object in source artifact")
    output.mkdir(mode=0o750)
    receipt = {
        "schema_version": "masld-bench-borzoi-enformer-source-admission-v1",
        "status": "pass",
        "small_sources": small_sources,
        "target_contracts": {
            "borzoi_ensemble": borzoi_targets,
            "enformer": enformer_targets,
        },
        "borzoi_checkpoint_parameters": parameter_receipt,
        "checkpoint_metadata": checkpoint_metadata,
        "checkpoint_bytes_downloaded": False,
        "checkpoint_bytes_deserialized": False,
        "model_forward_executed": False,
        "observed_outcomes_loaded": False,
        "weight_terms": {
            "borzoi_ensemble": "UNDECLARED_BLOCKED",
            "enformer": "UNDECLARED_BLOCKED",
        },
        "runtime_disposition": {
            "borzoi_ensemble": "blocked_exact_Baskerville_Westminster_and_no_network_runtime_unresolved",
            "enformer": "blocked_native_TF2.4.1_Sonnet2_runtime_not_L40S_admitted",
        },
        "independent_source_license_exposure_review": audit_receipt,
        "terminal_disposition": "source_and_fixture_admitted_checkpoint_execution_blocked",
        "open_champion_eligible": False,
    }
    (output / "source_admission_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    validate(
        arguments.config,
        arguments.source_dir,
        arguments.project_root,
        arguments.output,
    )


if __name__ == "__main__":
    main()
