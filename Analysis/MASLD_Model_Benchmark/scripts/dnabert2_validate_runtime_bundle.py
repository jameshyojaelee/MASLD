#!/usr/bin/env python3
"""Validate the complete revision-pinned DNABERT-2 runtime bundle."""

from __future__ import annotations

import argparse
import ast
from hashlib import sha1, sha256
import json
from pathlib import Path


EXPECTED_FILES = (
    "LICENSE",
    "bert_layers.py",
    "bert_padding.py",
    "config.json",
    "configuration_bert.py",
    "flash_attn_triton.py",
    "tokenizer.json",
    "tokenizer_config.json",
)


class DnabertBundleError(ValueError):
    """Raised when the DNABERT-2 bundle differs from the pinned tree."""


def _sha256(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _blob_sha1(path: Path) -> str:
    payload = path.read_bytes()
    return sha1(f"blob {len(payload)}\0".encode("ascii") + payload).hexdigest()


def _imports(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=path.name)
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            names.add(("." * node.level) + (node.module or ""))
    return sorted(names)


def validate(source: Path, bundle: Path, output: Path) -> dict[str, object]:
    if output.exists() or source.is_symlink() or bundle.is_symlink():
        raise DnabertBundleError("DNABERT-2 bundle request differs")
    upstream = json.loads((source / "review" / "source_admission_receipt.json").read_text())
    if upstream["status"] != "pass" or not upstream[
        "independent_source_license_exposure_review"
    ]["review_passed_before_model_execution"]:
        raise DnabertBundleError("DNABERT-2 upstream admission differs")
    tree = json.loads((source / "sources" / "dnabert2_hf_tree_api.json").read_text())
    identities = {row["path"]: row for row in tree}
    observed = sorted(path.name for path in bundle.iterdir() if path.is_file())
    if observed != sorted(EXPECTED_FILES):
        raise DnabertBundleError("DNABERT-2 runtime file census differs")
    files = []
    for name in EXPECTED_FILES:
        path = bundle / name
        row = identities.get(name)
        if (
            row is None
            or path.is_symlink()
            or path.stat().st_size != row["size"]
            or _blob_sha1(path) != row["oid"]
        ):
            raise DnabertBundleError(f"DNABERT-2 Git identity differs: {name}")
        record = {
            "filename": name,
            "size_bytes": path.stat().st_size,
            "git_blob_sha1": _blob_sha1(path),
            "sha256": _sha256(path),
        }
        if name.endswith(".py"):
            record["imports"] = _imports(path)
        files.append(record)
    output.mkdir(mode=0o750)
    receipt = {
        "schema_version": "masld-bench-dnabert2-runtime-bundle-v1",
        "status": "pass",
        "checkpoint_revision": upstream["checkpoint_metadata"]["dnabert2"][
            "checkpoint_revision"
        ],
        "license": "apache-2.0_from_revision_pinned_LICENSE",
        "files": files,
        "custom_code_executed": False,
        "checkpoint_loaded": False,
        "model_forward_executed": False,
        "observed_outcomes_loaded": False,
        "next_allowed_action": "no_network_import_then_converted_checkpoint_key_load",
    }
    (output / "runtime_bundle_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    validate(args.source, args.bundle, args.output)


if __name__ == "__main__":
    main()
