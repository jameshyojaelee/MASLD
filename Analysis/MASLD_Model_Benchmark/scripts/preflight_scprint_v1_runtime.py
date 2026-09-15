#!/usr/bin/env python3
"""Outcome-blind restoration preflight for the frozen scPRINT v1.5 checkpoint."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
from typing import Any


CHECKPOINT_SHA256 = "a4cf0753270d4ff451a5dbddadd4c59a53e10e4c3e96a8af0db407ad893c36c5"
CHECKPOINT_SIZE = 221_775_592
EXPECTED_SCPRINT_VERSION = "1.6.4"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def source_record(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve(strict=True)),
        "sha256": digest(path),
        "size_bytes": path.stat().st_size,
    }


def run(checkpoint: Path) -> dict[str, Any]:
    checkpoint = checkpoint.resolve(strict=True)
    if checkpoint.stat().st_size != CHECKPOINT_SIZE:
        raise RuntimeError("scPRINT checkpoint size differs from the frozen identity")
    if digest(checkpoint) != CHECKPOINT_SHA256:
        raise RuntimeError("scPRINT checkpoint SHA-256 differs from the frozen identity")

    import torch
    import scprint
    from scprint import scPrint

    package_version = importlib.metadata.version("scprint")
    if package_version != EXPECTED_SCPRINT_VERSION:
        raise RuntimeError(
            f"scPRINT version differs: expected {EXPECTED_SCPRINT_VERSION}, observed {package_version}"
        )
    package_file = Path(scprint.__file__).resolve(strict=True)
    model_source = package_file.parent / "model" / "model.py"
    encoder_source = package_file.parent / "model" / "encoders.py"

    payload = torch.load(
        checkpoint,
        map_location="cpu",
        weights_only=True,
        mmap=True,
    )
    if not isinstance(payload, dict):
        raise RuntimeError("restricted checkpoint payload is not a dictionary")
    state = payload.get("state_dict")
    hyperparameters = payload.get("hyper_parameters")
    if not isinstance(state, dict) or not isinstance(hyperparameters, dict):
        raise RuntimeError("checkpoint lacks state_dict or hyper_parameters")
    genes = hyperparameters.get("genes")
    if not isinstance(genes, list) or len(genes) != 44_756:
        raise RuntimeError("checkpoint gene vocabulary differs")

    constructor = dict(hyperparameters)
    constructor["precpt_gene_emb"] = None
    raw_classes = constructor.get("classes")
    class_mapping_reconstruction: dict[str, int] | None = None
    if isinstance(raw_classes, list):
        label_decoders = constructor.get("label_decoders")
        if not isinstance(label_decoders, dict):
            raise RuntimeError(
                "checkpoint classes is a list but label_decoders is unavailable"
            )
        class_mapping_reconstruction = {}
        for name in raw_classes:
            decoder = label_decoders.get(name)
            if not isinstance(name, str) or not isinstance(decoder, dict) or not decoder:
                raise RuntimeError(
                    "checkpoint classes cannot be reconstructed from label_decoders"
                )
            output_bias = state.get(f"cls_decoders.{name}.out_layer.bias")
            if output_bias is None or getattr(output_bias, "ndim", None) != 1:
                raise RuntimeError(
                    f"checkpoint has no one-dimensional classifier bias for {name}"
                )
            class_mapping_reconstruction[name] = int(output_bias.shape[0])
        constructor["classes"] = class_mapping_reconstruction
    original_hierarchy = constructor.get("labels_hierarchy")
    # The released constructor materializes hierarchy tensors with obsolete
    # global category offsets and crashes before state restoration. These
    # tensors are used only by supervised training losses, never frozen
    # encoder inference, and are not checkpoint parameters. Disable them for
    # this restoration probe rather than altering any learned tensor.
    constructor["labels_hierarchy"] = {}

    restoration_error = None
    missing: list[str] = []
    unexpected: list[str] = []
    try:
        model = scPrint(**constructor)
        incompatibility = model.load_state_dict(state, strict=False)
        missing = sorted(incompatibility.missing_keys)
        unexpected = sorted(incompatibility.unexpected_keys)
    except Exception as error:  # frozen terminal diagnostic, never silently repaired
        restoration_error = f"{type(error).__name__}: {error}"
    status = (
        "pass"
        if restoration_error is None and not missing and not unexpected
        else "failed_strict_restore"
    )

    return {
        "schema_version": "masld-bench-scprint-v1-runtime-preflight-v1",
        "status": status,
        "outcomes_or_labels_read": False,
        "checkpoint": {
            "sha256": CHECKPOINT_SHA256,
            "size_bytes": CHECKPOINT_SIZE,
            "tensor_count": len(state),
            "gene_vocabulary_size": len(genes),
        },
        "runtime": {
            "python": sys.version,
            "torch": torch.__version__,
            "scprint": package_version,
            "cuda_available": torch.cuda.is_available(),
        },
        "architecture": {
            "d_model": hyperparameters.get("d_model"),
            "nhead": hyperparameters.get("nhead"),
            "nlayers": hyperparameters.get("nlayers"),
            "transformer": hyperparameters.get("transformer"),
            "cell_emb_style": hyperparameters.get("cell_emb_style"),
            "checkpoint_classes_container": type(raw_classes).__name__,
            "class_mapping_reconstruction": class_mapping_reconstruction,
            "hierarchy_disabled_for_inference_only": bool(original_hierarchy),
            "missing_keys": missing,
            "restoration_error": restoration_error,
            "unexpected_keys": unexpected,
        },
        "source_records": [
            source_record(package_file),
            source_record(model_source),
            source_record(encoder_source),
        ],
        "interpretation": (
            "Strict architecture restoration only; this does not admit the runtime, "
            "run inference, fit a common head, or support a performance claim."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    result = run(args.checkpoint)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
