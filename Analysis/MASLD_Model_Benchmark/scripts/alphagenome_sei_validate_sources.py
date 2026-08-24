#!/usr/bin/env python3
"""Freeze AlphaGenome and Sei source, license, exposure, and checkpoint gates."""

from __future__ import annotations

import argparse
from hashlib import sha1, sha256
import json
from pathlib import Path
from typing import Mapping
import tomllib


class AdmissionSourceError(ValueError):
    """Raised when an official source differs from the admission contract."""


def _hash(path: Path) -> str:
    value = sha256()
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


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def _load_json(path: Path) -> object:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 20_000_000:
        raise AdmissionSourceError(f"remote metadata differs: {path.name}")
    return json.loads(path.read_text(encoding="utf-8"))


def _normalize_oid(value: object) -> str:
    text = str(value or "")
    return text.split(":", 1)[-1]


def _validate_alphagenome_remote(
    source_dir: Path, model: Mapping[str, object]
) -> dict[str, object]:
    api = _load_json(source_dir / "alphagenome_hf_model_api.json")
    tree = _load_json(source_dir / "alphagenome_hf_tree_api.json")
    if not isinstance(api, dict) or not isinstance(tree, list):
        raise AdmissionSourceError("AlphaGenome Hugging Face API schema differs")
    if (
        api.get("sha") != model["checkpoint_revision"]
        or api.get("gated") not in ("auto", True)
        or api.get("private") is True
    ):
        raise AdmissionSourceError("AlphaGenome revision or gate differs")
    by_path = {
        str(row.get("path") or row.get("rfilename")): row
        for row in tree
        if isinstance(row, dict)
    }
    receipts = []
    for expected in model["checkpoint_objects"]:
        row = by_path.get(str(expected["name"]))
        if row is None or int(row.get("size", -1)) != expected["size_bytes"]:
            raise AdmissionSourceError(
                f"AlphaGenome checkpoint inventory differs: {expected['name']}"
            )
        lfs = row.get("lfs") or {}
        public_lfs_identity = _normalize_oid(lfs.get("oid") or lfs.get("sha256"))
        if "lfs_sha256" in expected:
            pointer_identity = _normalize_oid(row.get("oid") or row.get("blobId"))
            if (
                pointer_identity != expected["git_blob_sha1"]
                or int(lfs.get("size", -1)) != expected["size_bytes"]
            ):
                raise AdmissionSourceError("AlphaGenome gated LFS identity differs")
            observed = pointer_identity
            identity_source = (
                "pinned_git_pointer_and_size;_LFS_SHA256_from_bound_prior_audit_not_publicly_reverified"
                if set(public_lfs_identity) == {"*"}
                else "public_LFS_metadata_and_pinned_git_pointer"
            )
            if set(public_lfs_identity) != {"*"} and public_lfs_identity != expected["lfs_sha256"]:
                raise AdmissionSourceError("AlphaGenome public LFS identity differs")
        else:
            observed = _normalize_oid(row.get("oid") or row.get("blobId"))
            if observed != expected["git_blob_sha1"]:
                raise AdmissionSourceError("AlphaGenome git-blob identity differs")
            pointer_identity = observed
            identity_source = "git_blob"
        receipts.append(
            {
                "name": expected["name"],
                "size_bytes": expected["size_bytes"],
                "object_identity": observed,
                "git_pointer_identity": pointer_identity,
                "identity_source": identity_source,
                "bound_prior_audit_lfs_sha256": expected.get("lfs_sha256"),
                "lfs_sha256_publicly_reverified": (
                    "lfs_sha256" not in expected or set(public_lfs_identity) != {"*"}
                ),
                "large_checkpoint_bytes_downloaded": False,
            }
        )
    card = source_dir / "alphagenome_model_card.md"
    card_expected = next(
        row for row in model["checkpoint_objects"] if row["name"] == "README.md"
    )
    if (
        card.is_symlink()
        or not card.is_file()
        or card.stat().st_size != card_expected["size_bytes"]
        or _git_blob_sha1(card) != card_expected["git_blob_sha1"]
    ):
        raise AdmissionSourceError("AlphaGenome model card differs")
    return {
        "checkpoint_revision": api["sha"],
        "gate": api["gated"],
        "inventory_objects": receipts,
        "inventory_size_bytes": sum(int(row["size_bytes"]) for row in receipts),
        "model_card_git_blob_sha1": _git_blob_sha1(card),
        "checkpoint_bytes_downloaded": False,
    }


def _validate_sei_zenodo(
    source_dir: Path, model: Mapping[str, object]
) -> dict[str, object]:
    raw = _load_json(source_dir / "sei_zenodo_record_api.json")
    if not isinstance(raw, dict) or int(raw.get("id", -1)) != 4906997:
        raise AdmissionSourceError("Sei Zenodo record identity differs")
    files = {
        str(row.get("key")): row for row in raw.get("files", []) if isinstance(row, dict)
    }
    record = files.get("sei_model.tar.gz")
    if record is None:
        raise AdmissionSourceError("Sei archive is absent from its Zenodo record")
    checksum = str(record.get("checksum", ""))
    if (
        int(record.get("size", -1)) != model["checkpoint_size_bytes"]
        or checksum.replace("md5:", "") != model["checkpoint_md5"]
    ):
        raise AdmissionSourceError("Sei archive metadata differs")
    license_id = str((raw.get("metadata") or {}).get("license", {}).get("id", ""))
    if license_id.upper().replace("_", "-") != "CC-BY-4.0":
        raise AdmissionSourceError("Sei Zenodo license differs")
    return {
        "record_id": 4906997,
        "filename": "sei_model.tar.gz",
        "size_bytes": int(record["size"]),
        "md5": checksum.replace("md5:", ""),
        "record_license": license_id,
        "checkpoint_bytes_downloaded": False,
        "checkpoint_deserialized": False,
    }


def _registry_row(path: Path, model_id: str) -> Mapping[str, object]:
    raw = tomllib.loads(path.read_text(encoding="utf-8"))
    matches = [row for row in raw["models"] if row["model_id"] == model_id]
    if len(matches) != 1:
        raise AdmissionSourceError(f"registry identity differs: {model_id}")
    return matches[0]


def _validate_audits(
    project_root: Path, config: Mapping[str, object]
) -> dict[str, object]:
    bindings = {}
    for name, expected in config["audit_bindings"].items():
        path = (project_root / expected["path"]).resolve(strict=True)
        if project_root.resolve() not in path.parents or path.is_symlink():
            raise AdmissionSourceError("audit binding escapes benchmark root")
        digest = _hash(path)
        if digest != expected["sha256"]:
            raise AdmissionSourceError(f"audit binding differs: {name}")
        bindings[name] = {"path": expected["path"], "sha256": digest}
    alpha_exposure = json.loads(
        (project_root / config["audit_bindings"]["alphagenome_exposure_audit"]["path"])
        .read_text(encoding="utf-8")
    )
    sei_exposure = json.loads(
        (project_root / config["audit_bindings"]["sei_exposure_audit"]["path"])
        .read_text(encoding="utf-8")
    )
    alpha_registry = _registry_row(
        project_root / config["audit_bindings"]["alphagenome_registry"]["path"],
        "alphagenome",
    )
    sei_registry = _registry_row(
        project_root / config["audit_bindings"]["sei_registry"]["path"], "sei"
    )
    alpha_finding = alpha_exposure["checkpoint_findings"]["alphagenome_all_folds"]
    sei_finding = sei_exposure["checkpoint_findings"]["sei_zenodo_4906997"]
    if (
        alpha_finding["exposure_state"] != "target_label_unexposed"
        or sei_finding["exposure_state"] != "target_label_unexposed"
        or alpha_registry["status"] != "restricted_comparator"
        or sei_registry["status"] != "restricted_comparator"
        or alpha_registry["admission_blocking"] is not True
        or sei_registry["admission_blocking"] is not True
        or "noncommercial" not in alpha_registry["license_status"]
        or "research_only" not in sei_registry["license_status"]
    ):
        raise AdmissionSourceError("source/license/exposure disposition differs")
    return {
        "bindings": bindings,
        "gse289173_exposure": {
            "alphagenome_all_folds": "target_label_unexposed",
            "sei": "target_label_unexposed",
        },
        "registry_status": {
            "alphagenome_all_folds": "restricted_comparator_blocked",
            "sei": "restricted_comparator_pending_archive_and_runtime",
        },
        "review_passed_before_model_execution": True,
    }


def validate(
    config_path: Path, source_dir: Path, project_root: Path, output: Path
) -> dict[str, object]:
    if output.exists() or config_path.is_symlink() or source_dir.is_symlink():
        raise AdmissionSourceError("source validation request differs")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("schema_version") != "masld-bench-alphagenome-sei-admission-v1":
        raise AdmissionSourceError("admission config schema differs")
    expected_disposition = {
        "source_and_fixture_probe_allowed": True,
        "alphagenome_checkpoint_download_allowed": False,
        "alphagenome_gpu_forward_allowed": False,
        "sei_checkpoint_download_allowed_after_source_review": True,
        "sei_checkpoint_deserialization_allowed": False,
        "sei_gpu_forward_allowed": False,
        "observed_outcomes_allowed": False,
        "open_champion_eligible": False,
    }
    if config.get("admission_disposition") != expected_disposition:
        raise AdmissionSourceError("admission disposition differs")
    small_sources = {
        name: _validate_file(source_dir / name, expected)
        for name, expected in config["source_files"].items()
    }
    audits = _validate_audits(project_root, config)
    alpha = _validate_alphagenome_remote(
        source_dir, config["models"]["alphagenome_all_folds"]
    )
    sei = _validate_sei_zenodo(source_dir, config["models"]["sei"])
    if any(
        path.stat().st_size > 4_000_000
        for path in source_dir.iterdir()
        if path.is_file()
    ):
        raise AdmissionSourceError("checkpoint-sized object found in source review")
    output.mkdir(mode=0o750)
    receipt = {
        "schema_version": "masld-bench-alphagenome-sei-source-admission-v1",
        "status": "pass",
        "small_sources": small_sources,
        "remote_checkpoint_metadata": {
            "alphagenome_all_folds": alpha,
            "sei": sei,
        },
        "independent_source_license_exposure_review": audits,
        "checkpoint_bytes_downloaded": False,
        "checkpoint_bytes_deserialized": False,
        "model_forward_executed": False,
        "observed_outcomes_loaded": False,
        "model_dispositions": {
            "alphagenome_all_folds": {
                "terminal_disposition": config["models"]["alphagenome_all_folds"][
                    "terminal_disposition"
                ],
                "next_model_job_allowed": False,
                "open_champion_eligible": False,
            },
            "sei": {
                "terminal_disposition": "source_review_admitted_safe_archive_acquisition_only",
                "safe_archive_acquisition_allowed": True,
                "checkpoint_deserialization_allowed": False,
                "gpu_forward_allowed": False,
                "open_champion_eligible": False,
            },
        },
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
    args = parser.parse_args()
    validate(args.config, args.source_dir, args.project_root, args.output)


if __name__ == "__main__":
    main()
