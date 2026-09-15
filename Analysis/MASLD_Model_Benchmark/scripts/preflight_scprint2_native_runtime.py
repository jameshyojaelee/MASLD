#!/usr/bin/env python3
"""Restricted exact-source restoration preflight for scPRINT-2 small-v2."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
from typing import Any


CHECKPOINT_SHA256 = "2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822"
CHECKPOINT_SIZE = 889_706_878
SOURCE_COMMIT = "7ee2de5aaa3492f326f540f703d44ddeaf1f2fd4"
SOURCE_RECORDS = {
    "LICENSE": "ae49d34bb24814b78a8c3390342986b3df683fe719864f9d6c4569ad27f0ec9a",
    "pyproject.toml": "9a24613df23178ff23694fede78172047a116449b54d84103c4baf411633b83a",
    "scprint2/model/model.py": "7866483f99ab50886c12a18bfedfff57bfef0993fe80b28d706136005f73103a",
    "scprint2/tasks/cell_emb.py": "b08acbdafa16f39e006ecf5441dd0604a0d8938b7bfb35c1618a7b0974e76e66",
    "uv.lock": "088c2bec476ab072b57893d6839dec8014fdddbf8eb1e496559b1ab1338dbd22",
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def source_record(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve(strict=True)),
        "sha256": digest(path),
        "size_bytes": path.stat().st_size,
    }


def run(checkpoint: Path, source: Path) -> dict[str, Any]:
    checkpoint = checkpoint.resolve(strict=True)
    source = source.resolve(strict=True)
    if checkpoint.stat().st_size != CHECKPOINT_SIZE or digest(checkpoint) != CHECKPOINT_SHA256:
        raise RuntimeError("scPRINT-2 checkpoint identity differs")
    records = []
    for relative, expected in SOURCE_RECORDS.items():
        path = source / relative
        record = source_record(path)
        if record["sha256"] != expected:
            raise RuntimeError(f"scPRINT-2 source record differs: {relative}")
        records.append(record)

    import torch
    import scprint2
    from scprint2 import scPRINT2

    runtime = {
        "python": sys.version,
        "torch": torch.__version__,
        "scprint2": importlib.metadata.version("scprint2"),
        "cuda_available": torch.cuda.is_available(),
        "module_path": str(Path(scprint2.__file__).resolve(strict=True)),
    }
    if runtime["torch"].split("+", 1)[0] != "2.8.0" or runtime["scprint2"] != "1.0.3":
        raise RuntimeError("scPRINT-2 native runtime version differs")
    if source not in Path(runtime["module_path"]).parents:
        raise RuntimeError("scPRINT-2 imported outside the pinned source tree")

    payload = torch.load(checkpoint, map_location="cpu", weights_only=True, mmap=True)
    if not isinstance(payload, dict):
        raise RuntimeError("restricted checkpoint payload is not a dictionary")
    state = payload.get("state_dict")
    hyperparameters = payload.get("hyper_parameters")
    if not isinstance(state, dict) or not isinstance(hyperparameters, dict):
        raise RuntimeError("checkpoint lacks state_dict or hyper_parameters")
    genes = hyperparameters.get("genes")
    organisms = hyperparameters.get("organisms")
    if not isinstance(genes, dict) or not isinstance(organisms, list):
        raise RuntimeError("checkpoint multispecies vocabulary contract differs")
    if "NCBITaxon:9606" not in genes or "NCBITaxon:9606" not in organisms:
        raise RuntimeError("checkpoint lacks the registered human vocabulary")

    constructor = dict(hyperparameters)
    instantiator_metadata = constructor.pop("_instantiator", None)
    constructor["precpt_gene_emb"] = None
    constructor["gene_pos_file"] = None
    restoration_error = None
    missing: list[str] = []
    unexpected: list[str] = []
    parameter_count = None
    try:
        model = scPRINT2(**constructor)
        model.on_load_checkpoint(payload)
        incompatibility = model.load_state_dict(state, strict=False)
        missing = sorted(incompatibility.missing_keys)
        unexpected = sorted(incompatibility.unexpected_keys)
        parameter_count = sum(parameter.numel() for parameter in model.parameters())
    except Exception as error:
        restoration_error = f"{type(error).__name__}: {error}"
    status = (
        "pass"
        if restoration_error is None and not missing and not unexpected
        else "failed_strict_restore"
    )

    return {
        "schema_version": "masld-bench-scprint2-native-runtime-preflight-v1",
        "status": status,
        "checkpoint_deserialization": "torch_load_weights_only_true",
        "outcomes_or_labels_read": False,
        "source_commit": SOURCE_COMMIT,
        "checkpoint": {
            "sha256": CHECKPOINT_SHA256,
            "size_bytes": CHECKPOINT_SIZE,
            "state_tensor_count": len(state),
            "organism_count": len(genes),
            "human_gene_vocabulary_size": len(genes["NCBITaxon:9606"]),
        },
        "runtime": runtime,
        "restoration": {
            "missing_keys": missing,
            "unexpected_keys": unexpected,
            "restoration_error": restoration_error,
            "parameter_count": parameter_count,
            "lightning_instantiator_metadata_removed": instantiator_metadata is not None,
            "pretrained_gene_embedding_path_disabled": True,
            "gene_position_path_disabled": True,
        },
        "source_records": records,
        "interpretation": (
            "Exact-source strict restoration only. This does not admit inference, "
            "read benchmark cells, fit a head, or support a performance claim."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    result = run(args.checkpoint, args.source)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
