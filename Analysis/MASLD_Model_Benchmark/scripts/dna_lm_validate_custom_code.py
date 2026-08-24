#!/usr/bin/env python3
"""Validate complete pinned HyenaDNA and Caduceus custom-code bundles."""

from __future__ import annotations

import argparse
import ast
from hashlib import sha1, sha256
import json
from pathlib import Path
from typing import Mapping


EXPECTED_FILES = {
    "hyenadna": (
        "config.json",
        "configuration_hyena.py",
        "modeling_hyena.py",
        "tokenization_hyena.py",
    ),
    "caduceus": (
        "config.json",
        "configuration_caduceus.py",
        "modeling_caduceus.py",
        "modeling_rcps.py",
        "tokenization_caduceus.py",
    ),
}


class CustomCodeAdmissionError(ValueError):
    """Raised when a custom-code bundle differs from its pinned repository tree."""


def _sha256(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def _imports(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=path.name)
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            names.add(("." * node.level) + (node.module or ""))
    return sorted(names)


def validate(source_artifact: Path, code_root: Path, output: Path) -> dict[str, object]:
    if output.exists() or source_artifact.is_symlink() or code_root.is_symlink():
        raise CustomCodeAdmissionError("custom-code validation request differs")
    model_receipt = json.loads(
        (source_artifact / "review" / "source_admission_receipt.json").read_text()
    )
    if (
        model_receipt["status"] != "pass"
        or model_receipt["checkpoint_bytes_downloaded"]
        or model_receipt["model_forward_executed"]
    ):
        raise CustomCodeAdmissionError("upstream source admission differs")
    models: dict[str, object] = {}
    for model_id, filenames in EXPECTED_FILES.items():
        tree_path = source_artifact / "sources" / f"{model_id}_hf_tree_api.json"
        tree = json.loads(tree_path.read_text(encoding="utf-8"))
        tree_rows: Mapping[str, Mapping[str, object]] = {
            str(row.get("path")): row for row in tree if isinstance(row, dict)
        }
        directory = code_root / model_id
        observed_names = sorted(path.name for path in directory.iterdir() if path.is_file())
        if observed_names != sorted(filenames):
            raise CustomCodeAdmissionError(f"custom-code file census differs: {model_id}")
        records = []
        for filename in filenames:
            path = directory / filename
            row = tree_rows.get(filename)
            if (
                row is None
                or path.is_symlink()
                or not path.is_file()
                or path.stat().st_size != int(row.get("size", -1))
                or _git_blob_sha1(path) != row.get("oid")
            ):
                raise CustomCodeAdmissionError(
                    f"custom-code Git identity differs: {model_id}/{filename}"
                )
            record = {
                "filename": filename,
                "size_bytes": path.stat().st_size,
                "git_blob_sha1": _git_blob_sha1(path),
                "sha256": _sha256(path),
            }
            if filename.endswith(".py"):
                record["imports"] = _imports(path)
            records.append(record)
        models[model_id] = {
            "checkpoint_revision": model_receipt["checkpoint_metadata"][model_id][
                "checkpoint_revision"
            ],
            "files": records,
            "custom_code_executed": False,
        }
    output.mkdir(mode=0o750)
    receipt = {
        "schema_version": "masld-bench-dna-lm-custom-code-admission-v1",
        "status": "pass",
        "models": models,
        "custom_code_executed": False,
        "checkpoint_tensor_data_loaded": False,
        "model_forward_executed": False,
        "runtime_network_allowed": False,
        "observed_outcomes_loaded": False,
        "next_allowed_action": "isolated_no_network_import_and_common_6000bp_runtime_probe",
    }
    (output / "custom_code_admission_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-artifact", type=Path, required=True)
    parser.add_argument("--code-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    validate(args.source_artifact, args.code_root, args.output)


if __name__ == "__main__":
    main()
