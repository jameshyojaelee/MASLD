#!/usr/bin/env python3
"""Validate five DNA language-model source, license, and exposure gates."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping
import tomllib


MODEL_IDS = (
    "dnabert2",
    "nucleotide_transformer",
    "hyenadna",
    "caduceus",
    "evo",
)


class DnaLmAdmissionError(ValueError):
    """Raised when a DNA language-model admission source differs."""


def _hash(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _validate_file(path: Path, expected: Mapping[str, object]) -> dict[str, object]:
    if path.is_symlink() or not path.is_file():
        raise DnaLmAdmissionError(f"source is not a regular file: {path.name}")
    size = path.stat().st_size
    digest = _hash(path)
    if size != expected["size_bytes"] or digest != expected["sha256"]:
        raise DnaLmAdmissionError(f"source bytes differ: {path.name}")
    return {"path": path.name, "size_bytes": size, "sha256": digest}


def _load_json(path: Path) -> object:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 20_000_000:
        raise DnaLmAdmissionError(f"remote metadata differs: {path.name}")
    return json.loads(path.read_text(encoding="utf-8"))


def _normalize_oid(value: object) -> str:
    return str(value or "").split(":", 1)[-1]


def _validate_hf_model(
    source_dir: Path, model_id: str, model: Mapping[str, object]
) -> dict[str, object]:
    api = _load_json(source_dir / f"{model_id}_hf_model_api.json")
    tree = _load_json(source_dir / f"{model_id}_hf_tree_api.json")
    if not isinstance(api, dict) or not isinstance(tree, list):
        raise DnaLmAdmissionError(f"Hugging Face API schema differs: {model_id}")
    if (
        api.get("sha") != model["checkpoint_revision"]
        or api.get("private") is True
        or api.get("gated") is True
        or api.get("gated") == "auto"
    ):
        raise DnaLmAdmissionError(f"checkpoint revision or access differs: {model_id}")
    rows = {
        str(row.get("path") or row.get("rfilename")): row
        for row in tree
        if isinstance(row, dict)
    }
    checkpoint = rows.get(str(model["checkpoint_filename"]))
    if checkpoint is None or int(checkpoint.get("size", -1)) != model["checkpoint_size_bytes"]:
        raise DnaLmAdmissionError(f"checkpoint file identity differs: {model_id}")
    lfs = checkpoint.get("lfs") or {}
    lfs_identity = _normalize_oid(lfs.get("oid") or lfs.get("sha256"))
    if (
        lfs_identity != model["checkpoint_sha256"]
        or int(lfs.get("size", -1)) != model["checkpoint_size_bytes"]
    ):
        raise DnaLmAdmissionError(f"checkpoint LFS identity differs: {model_id}")
    card_data = api.get("cardData") or api.get("card_data") or {}
    return {
        "checkpoint_repository": model["checkpoint_repository"],
        "checkpoint_revision": api["sha"],
        "checkpoint_filename": model["checkpoint_filename"],
        "checkpoint_size_bytes": int(checkpoint["size"]),
        "checkpoint_sha256": lfs_identity,
        "checkpoint_git_pointer_sha1": _normalize_oid(checkpoint.get("oid")),
        "repository_gate": api.get("gated", False),
        "card_license": card_data.get("license"),
        "checkpoint_bytes_downloaded": False,
    }


def _registry_rows(path: Path) -> dict[str, Mapping[str, object]]:
    raw = tomllib.loads(path.read_text(encoding="utf-8"))
    rows = {row["model_id"]: row for row in raw["models"]}
    if not set(MODEL_IDS).issubset(rows):
        raise DnaLmAdmissionError("DNA language-model registry census differs")
    return {model_id: rows[model_id] for model_id in MODEL_IDS}


def _validate_audits(
    project_root: Path, config: Mapping[str, object]
) -> dict[str, object]:
    bindings = {}
    for name, expected in config["audit_bindings"].items():
        path = (project_root / expected["path"]).resolve(strict=True)
        if project_root.resolve() not in path.parents or path.is_symlink():
            raise DnaLmAdmissionError("audit binding escapes benchmark root")
        digest = _hash(path)
        if digest != expected["sha256"]:
            raise DnaLmAdmissionError(f"audit binding differs: {name}")
        bindings[name] = {"path": expected["path"], "sha256": digest}
    registry_path = project_root / config["audit_bindings"]["model_registry"]["path"]
    registry = _registry_rows(registry_path)
    expected_registry_status = {
        "dnabert2": "candidate",
        "nucleotide_transformer": "restricted_comparator",
        "hyenadna": "candidate",
        "caduceus": "candidate",
        "evo": "candidate",
    }
    exposure = {}
    for model_id in MODEL_IDS:
        audit_path = project_root / config["audit_bindings"][
            f"{model_id}_exposure_audit"
        ]["path"]
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        finding = audit.get("checkpoint_finding")
        if not isinstance(finding, dict) or finding.get("exposure_state") != "target_label_unexposed":
            raise DnaLmAdmissionError(f"checkpoint exposure differs: {model_id}")
        row = registry[model_id]
        if (
            row["status"] != expected_registry_status[model_id]
            or row["admission_blocking"] is not True
            or row["exposure_status"] != "target_label_unexposed"
        ):
            raise DnaLmAdmissionError(f"registry disposition differs: {model_id}")
        exposure[model_id] = {
            "gse289173_exposure": "target_label_unexposed",
            "reference_sequence_exposure": finding["reference_sequence_exposure"],
            "registry_status": row["status"],
            "admission_blocking_before_runtime": True,
        }
    return {
        "bindings": bindings,
        "models": exposure,
        "review_passed_before_checkpoint_download": True,
        "review_passed_before_model_execution": True,
    }


def validate(
    config_path: Path, source_dir: Path, project_root: Path, output: Path
) -> dict[str, object]:
    if output.exists() or config_path.is_symlink() or source_dir.is_symlink():
        raise DnaLmAdmissionError("source validation request differs")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("schema_version") != "masld-bench-dna-lm-admission-v1":
        raise DnaLmAdmissionError("admission config schema differs")
    if tuple(config.get("models", {})) != MODEL_IDS:
        raise DnaLmAdmissionError("DNA language-model census order differs")
    disposition = config.get("admission_disposition")
    if disposition != {
        "source_and_fixture_probe_allowed": True,
        "safe_checkpoint_download_allowed_after_source_review": True,
        "checkpoint_deserialization_allowed": False,
        "gpu_forward_allowed": False,
        "observed_outcomes_allowed": False,
        "supervised_eqtl_or_ieqtl_head_allowed": False,
    }:
        raise DnaLmAdmissionError("admission disposition differs")
    sources = {
        name: _validate_file(source_dir / name, expected)
        for name, expected in config["source_files"].items()
    }
    checkpoint_metadata = {
        model_id: _validate_hf_model(source_dir, model_id, config["models"][model_id])
        for model_id in MODEL_IDS
    }
    audits = _validate_audits(project_root, config)
    if any(
        path.stat().st_size > 2_000_000
        for path in source_dir.iterdir()
        if path.is_file()
    ):
        raise DnaLmAdmissionError("checkpoint-sized object found in source review")
    model_dispositions = {}
    for model_id in MODEL_IDS:
        model = config["models"][model_id]
        model_dispositions[model_id] = {
            "terminal_source_disposition": model["terminal_source_disposition"],
            "safe_checkpoint_download_allowed": model["safe_download_allowed"],
            "checkpoint_deserialization_allowed": False,
            "gpu_forward_allowed": False,
            "open_champion_eligible_after_all_runtime_and_task_gates": model[
                "open_champion_eligible_after_all_runtime_and_task_gates"
            ],
            "embeddings_biologically_interpretable_without_head": False,
            "supervised_eqtl_or_ieqtl_head_allowed": False,
        }
    output.mkdir(mode=0o750)
    receipt = {
        "schema_version": "masld-bench-dna-lm-source-admission-v1",
        "status": "pass",
        "small_sources": sources,
        "checkpoint_metadata": checkpoint_metadata,
        "independent_source_license_exposure_review": audits,
        "model_dispositions": model_dispositions,
        "checkpoint_bytes_downloaded": False,
        "checkpoint_bytes_deserialized": False,
        "model_forward_executed": False,
        "observed_outcomes_loaded": False,
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
