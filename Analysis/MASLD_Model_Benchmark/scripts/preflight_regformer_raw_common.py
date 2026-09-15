#!/usr/bin/env python3
"""Fail-closed CPU strict-restore preflight for raw RegFormer 512d embeddings."""

from __future__ import annotations

import argparse
import ast
from collections.abc import Mapping
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys
from typing import Any, Optional


SOURCE_COMMIT = "9618f1f4bb7b61c9b236d3ea8a6339c236f2da1b"
CHECKPOINT_SHA256 = "8e8e751a0b05caa4e6e4d6b5c0f8ee847ae607757f0431b4bee90cd9ef26a434"
CHECKPOINT_SIZE = 197_677_201
VOCAB_SHA256 = "acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f"
GRAPH_SHA256 = "f0a530d760c82ef9e264eec3f3f5a230f295fd3c3e3a6b347e2f84a784807397"
CONFIG_SHA256 = "f599274fed2ff2d80a96d2e318e40941ad4c67adb08946b72a990e37b25c57ea"
INVENTORY_SHA256 = "7686e84eabb467fc4b906600b521c902e496d2fa6279d7cef13465fa663b7cef"
INSPECTION_SHA256 = "28217cb286650ca4dd1e767870ef8a5e6496cc63e44410cfde87170ad2b459df"
RUNTIME_ARTIFACTS_SHA256 = "3d59792d8a13ca10a2df3d6386bc266eb0e108def6279188793d1cfee4fc1167"
OUTPUT_DIMENSION = 512
EXPECTED_STATE_TENSORS = 139
EXPECTED_STATE_NUMEL_WITH_TIED_DUPLICATE = 80_483_842
EXPECTED_VOCAB_BEFORE_MASK = 60_697
EXPECTED_VOCAB_AFTER_MASK = 60_698
EXPECTED_SOURCE_HASHES = {
    "LICENSE": "847f9048bf110be3e0df25ae50769b707be769b90eb18c8f201e67511da48111",
    "README.md": "e2aa5a3e2a9bd7b6d51cd500aef24ca185f759cbe3cd5e3761b453c394dfe868",
    "requirements.txt": "3b4061611d51b15b5f4f0f54e3b1066d0355104290ef474d45ec57de9dd39aa3",
    "downstream_task/regformer_emb.py": (
        "1fe56a1349cfd63593ca368408701ebb877deb75076037b2408234ed5390e0f8"
    ),
    "regformer/model/BiMamba.py": (
        "ebbc667e01751069f4b840128404540ab630d21c4d0c33c9a3d78ba7d6996016"
    ),
    "regformer/model/dsbn.py": "666e443b530b8f4c7a26759c18b2c3125d96bb8ae55b7d91f26599ae7f84fb9b",
    "regformer/model/layers.py": "e861d8d6241de4c0ab439fd9b6932f0fef10e2c6ab80d922e7284753c8f0a12c",
    "regformer/model/mambaLM.py": (
        "2392323c20f10d79b8ce1e2967e34e718cc0703984d1fccd65ef4258e3ccd06e"
    ),
    "regformer/utils/utils.py": "844e449ca8d2429b555f5d70f9772360e53fd9d067a770d4ffa17b34988f0406",
}
EXPECTED_RUNTIME = {
    "causal-conv1d": "1.2.0.post2",
    "mamba-ssm": "1.2.0.post1",
    "torch": "2.2.0",
}
REQUIRED_CONFIG_LITERALS = {
    "append_cls": "false",
    "batch_size": "4",
    "bimamba_type": "'none'",
    "cell_emb_style": "'avg-pool'",
    "data_is_raw": "false",
    "do_pretrain": "true",
    "do_train": "true",
    "dropout": "0.2",
    "graph_sort": "true",
    "include_zero_gene": "false",
    "input_emb_style": "'continuous'",
    "input_style": "'binned'",
    "layer_emb": "true",
    "layer_size": "512",
    "max_seq_len": "1200",
    "MVC": "true",
    "n_bins": "51",
    "nlayers": "10",
    "pre_norm": "false",
    "token_emb_freeze": "false",
}


class RegFormerPreflightError(RuntimeError):
    """Raised when an activation prerequisite is not exact."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def require_file(
    path: Path, expected_sha256: str, *, size: Optional[int] = None
) -> None:
    if not path.is_file() or path.is_symlink():
        raise RegFormerPreflightError(f"required regular file is absent: {path}")
    if size is not None and path.stat().st_size != size:
        raise RegFormerPreflightError(f"file size differs: {path}")
    if sha256_file(path) != expected_sha256:
        raise RegFormerPreflightError(f"file SHA-256 differs: {path}")


def git_output(source: Path, *arguments: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(source), *arguments],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def validate_source(source: Path) -> dict[str, Any]:
    if source.is_symlink() or not source.is_dir():
        raise RegFormerPreflightError("source checkout is missing or a symlink")
    try:
        commit = git_output(source, "rev-parse", "HEAD")
        dirty = git_output(source, "status", "--porcelain=v1", "--untracked-files=all")
    except (OSError, subprocess.CalledProcessError) as error:
        raise RegFormerPreflightError("source is not an auditable Git checkout") from error
    if commit != SOURCE_COMMIT or dirty:
        raise RegFormerPreflightError("source revision differs or checkout is dirty")
    for relative, expected in EXPECTED_SOURCE_HASHES.items():
        require_file(source / relative, expected)
    embedding_audit = audit_released_embedding_source(
        (source / "downstream_task/regformer_emb.py").read_text(encoding="utf-8")
    )
    return {
        "repository": "BGIResearch/RegFormer",
        "commit": commit,
        "core_file_sha256": dict(sorted(EXPECTED_SOURCE_HASHES.items())),
        "released_embedding_leakage_audit": embedding_audit,
    }


def audit_released_embedding_source(text: str) -> dict[str, Any]:
    tree = ast.parse(text)
    raw_pool_lines: list[int] = []
    refine_lines: list[int] = []
    label_reads: list[int] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = ""
            if isinstance(node.func, ast.Name):
                name = node.func.id
            elif isinstance(node.func, ast.Attribute):
                name = node.func.attr
            if name == "_get_cell_emb_from_layer":
                raw_pool_lines.append(node.lineno)
            if name == "refine_embedding":
                refine_lines.append(node.lineno)
        if (
            isinstance(node, ast.Subscript)
            and isinstance(node.value, ast.Name)
            and node.value.id == "batch"
        ):
            key = node.slice
            if hasattr(ast, "Index") and isinstance(key, ast.Index):
                key = key.value
            if isinstance(key, ast.Constant) and key.value == "celltype_labels":
                label_reads.append(node.lineno)
    if len(raw_pool_lines) != 1 or len(refine_lines) != 1 or not label_reads:
        raise RegFormerPreflightError("released embedding leakage signature differs")
    if refine_lines[0] <= raw_pool_lines[0]:
        raise RegFormerPreflightError("released refinement no longer follows raw pooling")
    return {
        "raw_pool_line": raw_pool_lines[0],
        "target_label_read_lines": sorted(label_reads),
        "forbidden_refine_embedding_line": refine_lines[0],
        "benchmark_action": "do_not_import_or_execute_released_embedding_driver",
    }


def parse_literal_assignments(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        result[key.strip()] = value.strip()
    return result


def validate_config(path: Path) -> dict[str, str]:
    require_file(path, CONFIG_SHA256)
    assignments = parse_literal_assignments(path.read_text(encoding="utf-8"))
    observed = {key: assignments.get(key) for key in REQUIRED_CONFIG_LITERALS}
    if observed != REQUIRED_CONFIG_LITERALS:
        raise RegFormerPreflightError("canonical model configuration differs")
    return dict(sorted(REQUIRED_CONFIG_LITERALS.items()))


def validate_vocab(path: Path) -> dict[str, Any]:
    require_file(path, VOCAB_SHA256)
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or len(value) != EXPECTED_VOCAB_BEFORE_MASK:
        raise RegFormerPreflightError("released vocabulary cardinality differs")
    indices = list(value.values())
    if (
        any(isinstance(item, bool) or not isinstance(item, int) for item in indices)
        or set(indices) != set(range(EXPECTED_VOCAB_BEFORE_MASK))
        or value.get("<pad>") != 60_694
        or value.get("<cls>") != 60_695
        or "<mask>" in value
    ):
        raise RegFormerPreflightError("released vocabulary indices differ")
    value["<mask>"] = len(value)
    if len(value) != EXPECTED_VOCAB_AFTER_MASK:
        raise RegFormerPreflightError("runtime vocabulary cardinality differs")
    return {
        "released_size": EXPECTED_VOCAB_BEFORE_MASK,
        "runtime_size_after_released_mask_append": EXPECTED_VOCAB_AFTER_MASK,
        "pad_id": value["<pad>"],
        "cls_id": value["<cls>"],
        "mask_id": value["<mask>"],
        "mapping": value,
    }


def inspection_schema(path: Path) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    require_file(path, INSPECTION_SHA256)
    value = json.loads(path.read_text(encoding="utf-8"))
    containers = value.get("containers")
    if (
        value.get("status") != "pass"
        or value.get("sha256") != CHECKPOINT_SHA256
        or value.get("size_bytes") != CHECKPOINT_SIZE
        or value.get("tensor_count") != EXPECTED_STATE_TENSORS
        or value.get("tensor_numel") != EXPECTED_STATE_NUMEL_WITH_TIED_DUPLICATE
        or not isinstance(containers, list)
        or len(containers) != 1
        or containers[0].get("length") != EXPECTED_STATE_TENSORS
        or containers[0].get("path") != "checkpoint"
        or containers[0].get("type") != "collections.OrderedDict"
        or not isinstance(containers[0].get("keys"), list)
        or len(containers[0]["keys"]) != EXPECTED_STATE_TENSORS
    ):
        raise RegFormerPreflightError("trusted weights-only inspection differs")
    schema: dict[str, dict[str, Any]] = {}
    for item in value.get("tensors", []):
        raw_path = item.get("path", "")
        if not raw_path.startswith("checkpoint."):
            raise RegFormerPreflightError("inspection tensor path differs")
        key = raw_path[len("checkpoint.") :]
        if key in schema:
            raise RegFormerPreflightError("inspection repeats a state key")
        schema[key] = {"shape": item.get("shape"), "dtype": item.get("dtype")}
    if len(schema) != EXPECTED_STATE_TENSORS:
        raise RegFormerPreflightError("inspection state schema cardinality differs")
    return schema, {
        "inspection_sha256": INSPECTION_SHA256,
        "state_schema_sha256": canonical_sha256(schema),
        "state_tensor_count": len(schema),
        "state_numel_with_tied_duplicate": value["tensor_numel"],
    }


def validate_static_inputs(arguments: argparse.Namespace) -> dict[str, Any]:
    require_file(arguments.runtime_artifacts, RUNTIME_ARTIFACTS_SHA256)
    require_file(arguments.checkpoint, CHECKPOINT_SHA256, size=CHECKPOINT_SIZE)
    require_file(arguments.graph, GRAPH_SHA256)
    require_file(arguments.inventory, INVENTORY_SHA256)
    inventory = json.loads(arguments.inventory.read_text(encoding="utf-8"))
    if (
        inventory.get("checkpoint_deserialized") is not False
        or inventory.get("archive_extracted") is not False
        or inventory.get("container_type") != "PyTorch_ZIP_with_pickle_metadata"
    ):
        raise RegFormerPreflightError("safe checkpoint inventory differs")
    schema, schema_audit = inspection_schema(arguments.inspection)
    config_audit = validate_config(arguments.config)
    vocab_audit = validate_vocab(arguments.vocab)
    source_audit = validate_source(arguments.source)
    return {
        "source": source_audit,
        "checkpoint": {
            "sha256": CHECKPOINT_SHA256,
            "size_bytes": CHECKPOINT_SIZE,
            "safe_inventory_sha256": INVENTORY_SHA256,
            **schema_audit,
        },
        "assets": {
            "vocab_sha256": VOCAB_SHA256,
            "graph_sha256": GRAPH_SHA256,
            "config_sha256": CONFIG_SHA256,
            "runtime_artifacts_sha256": RUNTIME_ARTIFACTS_SHA256,
        },
        "config": config_audit,
        "vocab": {key: value for key, value in vocab_audit.items() if key != "mapping"},
        "state_schema": schema,
        "runtime_vocab": vocab_audit["mapping"],
    }


def strict_restore(arguments: argparse.Namespace, audit: dict[str, Any]) -> dict[str, Any]:
    import numpy as np
    import torch

    observed_runtime = {name: version(name) for name in EXPECTED_RUNTIME}
    if observed_runtime != EXPECTED_RUNTIME or sys.version_info[:2] != (3, 8):
        raise RegFormerPreflightError(
            f"isolated compatibility runtime differs: {observed_runtime!r}"
        )
    sys.path.insert(0, str(arguments.source))
    from regformer.model.mambaLM import MambaModel

    if "regformer.utils.utils" in sys.modules:
        raise RegFormerPreflightError("label-refining utility entered execution path")
    runtime_vocab = audit.pop("runtime_vocab")
    model = MambaModel(
        ntoken=len(runtime_vocab),
        d_model=512,
        nlayers=10,
        device="cpu",
        vocab=runtime_vocab,
        dropout=0.2,
        pad_token="<pad>",
        pad_value=-2,
        do_mvc=True,
        do_dab=False,
        domain_spec_batchnorm=False,
        num_batch_labels=None,
        n_input_bins=51,
        input_emb_style="continuous",
        cell_emb_style="avg-pool",
        pre_norm=False,
        do_pretrain=True,
        topo_graph=True,
        if_bimamba=False,
        bimamba_type="none",
        token_emb_freeze=False,
        only_value_emb=False,
        bin_cls=False,
        bin_nums=51,
        use_transformer=False,
    )
    state = torch.load(
        arguments.checkpoint,
        map_location="cpu",
        mmap=True,
        weights_only=True,
    )
    if not isinstance(state, Mapping) or not all(
        isinstance(key, str) and isinstance(value, torch.Tensor)
        for key, value in state.items()
    ):
        raise RegFormerPreflightError("checkpoint is not a tensor-only state mapping")
    expected_schema = audit.pop("state_schema")
    observed_schema = {
        key: {"shape": list(value.shape), "dtype": str(value.dtype)}
        for key, value in state.items()
    }
    if observed_schema != expected_schema:
        raise RegFormerPreflightError("loaded checkpoint state schema differs")
    incompatibility = model.load_state_dict(state, strict=True)
    if incompatibility.missing_keys or incompatibility.unexpected_keys:
        raise RegFormerPreflightError("strict checkpoint restore was incomplete")
    if model.lm_head.weight is not model.encoder.embedding.weight:
        raise RegFormerPreflightError("released tied embedding identity differs")

    model.eval()
    hidden = torch.linspace(-2.0, 2.0, steps=6 * OUTPUT_DIMENSION).reshape(
        2, 3, OUTPUT_DIMENSION
    )
    padding = torch.tensor([[False, False, True], [False, False, False]])
    with torch.inference_mode():
        first = model._get_cell_emb_from_layer(
            hidden, weights=None, src_key_padding_mask=padding
        )
        permutation = torch.tensor([1, 0])
        permuted = model._get_cell_emb_from_layer(
            hidden[permutation],
            weights=None,
            src_key_padding_mask=padding[permutation],
        )
    if (
        tuple(first.shape) != (2, OUTPUT_DIMENSION)
        or not torch.isfinite(first).all()
        or not torch.equal(first[permutation], permuted)
    ):
        raise RegFormerPreflightError("raw average-pool 512d probe differs")
    if "regformer.utils.utils" in sys.modules:
        raise RegFormerPreflightError("label-refining utility entered execution path")
    parameter_count = sum(parameter.numel() for parameter in model.parameters())
    buffer_count = sum(buffer.numel() for buffer in model.buffers())
    return {
        "python": sys.version,
        "versions": observed_runtime,
        "device": "cpu",
        "strict_restore": True,
        "missing_keys": [],
        "unexpected_keys": [],
        "parameter_count_unique": parameter_count,
        "buffer_count": buffer_count,
        "tied_lm_head_encoder_embedding_identity": True,
        "raw_common_representation": {
            "method": "MambaModel._get_cell_emb_from_layer",
            "style": "avg-pool_then_checkpoint_cell_norm",
            "dimension": OUTPUT_DIMENSION,
            "finite": True,
            "batch_row_order_equivariant": True,
            "label_inputs": [],
            "refine_embedding_imported": False,
            "refine_embedding_called": False,
        },
        "compiled_mamba_imported": True,
        "full_encoder_forward_executed": False,
        "graph_deserialized": False,
        "checkpoint_moved_to_gpu": False,
        "numpy_version": np.__version__,
    }


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--source", type=Path, required=True)
    value.add_argument("--checkpoint", type=Path, required=True)
    value.add_argument("--vocab", type=Path, required=True)
    value.add_argument("--graph", type=Path, required=True)
    value.add_argument("--config", type=Path, required=True)
    value.add_argument("--inventory", type=Path, required=True)
    value.add_argument("--inspection", type=Path, required=True)
    value.add_argument("--runtime-artifacts", type=Path, required=True)
    value.add_argument("--output", type=Path, required=True)
    return value


def main() -> int:
    arguments = parser().parse_args()
    if arguments.output.exists():
        raise RegFormerPreflightError(f"refusing to replace output: {arguments.output}")
    audit = validate_static_inputs(arguments)
    runtime = strict_restore(arguments, audit)
    receipt = {
        "schema_version": "masld-bench-regformer-raw-common-preflight-v1",
        "status": "pass_strict_restore_cpu_compatibility",
        "model_id": "regformer",
        "common_representation": "raw_unrefined_avg_pooled_512d_encoder_embedding",
        "label_refined_embeddings_used": False,
        "observed_labels_loaded": False,
        "sealed_outcomes_loaded": False,
        "histology_loaded": False,
        "source_and_assets": audit,
        "runtime": runtime,
        "activation_disposition": "blocked_pending_canonical_runtime_and_compiled_forward",
        "remaining_blockers": [
            (
                "The probe runtime is Python 3.8/Torch 2.2/Mamba 1.2 "
                "compatibility evidence, not the upstream-declared Python 3.9/"
                "Torch 2.0/Mamba 1.1.1 runtime."
            ),
            (
                "The upstream Torch CUDA 11.6 and DGL CUDA 11.8 declarations "
                "remain internally mismatched."
            ),
            (
                "A compiled full encoder forward, deterministic row-ID gene "
                "subsampling, graph-order parity, and raw-count preprocessing "
                "parity have not passed."
            ),
            (
                "Checkpoint exposure for GSE289173 remains unknown, so "
                "sealed-champion use is ineligible."
            ),
            (
                "The frozen census still records the checkpoint SHA as "
                "unresolved and remains admission-blocking until this evidence "
                "is reviewed and promoted."
            ),
        ],
    }
    with arguments.output.open("x", encoding="utf-8") as handle:
        json.dump(receipt, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
    print(arguments.output.resolve().as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
