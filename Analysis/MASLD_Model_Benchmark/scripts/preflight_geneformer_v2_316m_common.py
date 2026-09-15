#!/usr/bin/env python3
"""Hash-bound CPU/GPU preflight for raw Geneformer V2-316M embeddings."""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import sys
from typing import Any


MODEL_ID = "geneformer_v2_316m"
EMBEDDING_POLICY = "last_hidden_state_mean_over_non_special_nonpadding_tokens"
HIDDEN_SIZE = 1152
GPU_BATCH_SIZE = 32
CHECKPOINT_SHA256 = "965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3"
CHECKPOINT_SIZE = 1_265_455_076
CONFIG_SHA256 = "2cc4af3442644e84af71814a607b61c958d7b4e5ebde6791ab77b5f534ac6f6e"
GENERATION_CONFIG_SHA256 = (
    "bf1ba320993168b0d6ee5ff170df2108454d239b28eaac5948d30ee6993cccc1"
)
BUNDLE_ARTIFACTS_SHA256 = (
    "b269cef53fa92cc7e7ab9c21a2b85de69904cdfaaf5881e9676c57dc2ae02ae0"
)
INVENTORY_SHA256 = "b28a589b9e2886776eaf888ea7e82e13c3ef3b25d365a1e4fc572f84a5b4c276"
TOKENS_SHA256 = "dece0446a19ce0746e860db615e82a53219c03ebd5cd13ff16580954a6faa38c"
TOKEN_SUMMARY_SHA256 = (
    "d7ba1227048487e485bdb80f6b2ee1018f9b9d80b20e9e838beac3f1049930ca"
)
PRIOR_EMBEDDINGS_SHA256 = (
    "361ac7ec520bd6dcc2a9628af0c45f8b8dd7c12b4cbb895998e4b66c851150aa"
)
PRIOR_ROW_ORDER_SHA256 = (
    "dacd7add5a3e30c02f935168905bed26e2a7c326d99efb48d6c89903d5d995c0"
)
ADAPTER_SHA256 = "08e285ba7cedaf3065d71e8c2537888721c9bb823982ad26ee034512739f0ff8"
EXTRACTION_SCRIPT_SHA256 = (
    "a5f2364570129f68e64ae3ab67480ed28372ab52298b29d84b491b91c465ab70"
)
SAFE_TENSOR_COUNT = 298
SAFE_TENSOR_ELEMENTS = 316_354_867
EXPECTED_RUNTIME = {
    "numpy": "2.4.6",
    "safetensors": "0.8.0",
    "tokenizers": "0.22.2",
    "torch": "2.6.0",
    "transformers": "4.56.2",
}
EXPECTED_MLM_ONLY_KEYS = {
    "cls.predictions.bias",
    "cls.predictions.transform.LayerNorm.bias",
    "cls.predictions.transform.LayerNorm.weight",
    "cls.predictions.transform.dense.bias",
    "cls.predictions.transform.dense.weight",
}


class GeneformerPreflightError(RuntimeError):
    """Raised when the compatibility runtime or raw embedding requirement differs."""


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


def require_file(path: Path, expected: str, size: int | None = None) -> None:
    if not path.is_file() or path.is_symlink():
        raise GeneformerPreflightError(f"required file is absent or unsafe: {path}")
    if size is not None and path.stat().st_size != size:
        raise GeneformerPreflightError(f"file size differs: {path}")
    if sha256_file(path) != expected:
        raise GeneformerPreflightError(f"file SHA-256 differs: {path}")


def validate_runtime() -> dict[str, Any]:
    observed = {name: version(name) for name in EXPECTED_RUNTIME}
    normalized = {**observed, "torch": observed["torch"].split("+", 1)[0]}
    if normalized != EXPECTED_RUNTIME or sys.version_info[:3] != (3, 11, 15):
        raise GeneformerPreflightError(f"captured compatibility runtime differs: {observed}")
    import torch

    if torch.version.cuda != "12.4":
        raise GeneformerPreflightError("captured Torch CUDA build differs")
    return {
        "python": sys.version,
        "packages": observed,
        "torch_cuda": torch.version.cuda,
        "runtime_origin": "captured_mutable_decima_environment",
        "runtime_claim": "compatibility_snapshot_not_upstream_canonical_environment",
    }


def validate_assets(arguments: argparse.Namespace) -> dict[str, Any]:
    bundle = arguments.bundle.resolve(strict=True)
    model = bundle / "upstream/Geneformer-V2-316M"
    require_file(bundle / "ARTIFACTS.json", BUNDLE_ARTIFACTS_SHA256)
    require_file(model / "model.safetensors", CHECKPOINT_SHA256, CHECKPOINT_SIZE)
    require_file(model / "config.json", CONFIG_SHA256)
    require_file(model / "generation_config.json", GENERATION_CONFIG_SHA256)
    require_file(arguments.inventory, INVENTORY_SHA256)
    require_file(arguments.tokens, TOKENS_SHA256)
    require_file(arguments.token_summary, TOKEN_SUMMARY_SHA256)
    require_file(arguments.prior_embeddings, PRIOR_EMBEDDINGS_SHA256)
    require_file(arguments.prior_row_order, PRIOR_ROW_ORDER_SHA256)
    require_file(arguments.adapter, ADAPTER_SHA256)
    require_file(arguments.extraction_script, EXTRACTION_SCRIPT_SHA256)

    config = json.loads((model / "config.json").read_text(encoding="utf-8"))
    expected_config = {
        "architectures": ["BertForMaskedLM"],
        "hidden_size": HIDDEN_SIZE,
        "intermediate_size": 4608,
        "max_position_embeddings": 4096,
        "model_type": "bert",
        "num_attention_heads": 18,
        "num_hidden_layers": 18,
        "pad_token_id": 0,
        "torch_dtype": "float32",
        "vocab_size": 20275,
    }
    if any(config.get(key) != value for key, value in expected_config.items()):
        raise GeneformerPreflightError("Geneformer V2-316M architecture differs")
    inventory = json.loads(arguments.inventory.read_text(encoding="utf-8"))
    if (
        inventory.get("status") != "pass"
        or inventory.get("sha256") != CHECKPOINT_SHA256
        or inventory.get("tensor_count") != SAFE_TENSOR_COUNT
        or inventory.get("parameter_and_buffer_elements") != SAFE_TENSOR_ELEMENTS
        or inventory.get("tensor_data_loaded") is not False
    ):
        raise GeneformerPreflightError("safe checkpoint inventory differs")
    schema = {
        item["name"]: {"shape": item["shape"], "dtype": item["dtype"]}
        for item in inventory.get("tensors", [])
    }
    if len(schema) != SAFE_TENSOR_COUNT:
        raise GeneformerPreflightError("safe checkpoint schema cardinality differs")

    summary = json.loads(arguments.token_summary.read_text(encoding="utf-8"))
    if summary != {
        "cells": 1000,
        "cells_truncated": 29,
        "embedding_policy": EMBEDDING_POLICY,
        "genes_dropped": 18531,
        "genes_total": 37533,
        "genes_usable": 19002,
        "max_tokens": 4096,
        "schema_version": "masld-bench-geneformer-tokenization-v1",
    }:
        raise GeneformerPreflightError("prior smoke tokenization summary differs")
    prior = validate_prior_smoke(
        arguments.tokens, arguments.prior_row_order, arguments.prior_embeddings
    )
    return {
        "bundle_artifacts_sha256": BUNDLE_ARTIFACTS_SHA256,
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "checkpoint_size_bytes": CHECKPOINT_SIZE,
        "config_sha256": CONFIG_SHA256,
        "generation_config_sha256": GENERATION_CONFIG_SHA256,
        "inventory_sha256": INVENTORY_SHA256,
        "checkpoint_schema_sha256": canonical_sha256(schema),
        "checkpoint_tensor_count": len(schema),
        "checkpoint_schema": schema,
        "model_directory": model,
        "prior_smoke": prior,
    }


def validate_prior_smoke(
    token_path: Path, row_order_path: Path, embedding_path: Path
) -> dict[str, Any]:
    import numpy as np

    row_order = row_order_path.read_text(encoding="utf-8").splitlines()
    if (
        len(row_order) != 1000
        or row_order != sorted(row_order)
        or len(set(row_order)) != len(row_order)
        or any(len(value) != 64 for value in row_order)
    ):
        raise GeneformerPreflightError("prior embedding row order differs")
    tokens: dict[str, list[int]] = {}
    with token_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            row_id = str(record["row_hash"])
            values = record["tokens"]
            if (
                row_id in tokens
                or not isinstance(values, list)
                or not values
                or len(values) > 4096
                or any(isinstance(value, bool) or not isinstance(value, int) for value in values)
            ):
                raise GeneformerPreflightError("prior rank-token row differs")
            tokens[row_id] = values
    if set(tokens) != set(row_order):
        raise GeneformerPreflightError("prior tokens and embedding rows differ")
    embeddings = np.load(embedding_path, mmap_mode="r", allow_pickle=False)
    if (
        embeddings.shape != (1000, HIDDEN_SIZE)
        or embeddings.dtype != np.dtype("float32")
        or not np.isfinite(embeddings).all()
    ):
        raise GeneformerPreflightError("prior raw embedding matrix differs")
    return {
        "tokens_sha256": TOKENS_SHA256,
        "token_summary_sha256": TOKEN_SUMMARY_SHA256,
        "embeddings_sha256": PRIOR_EMBEDDINGS_SHA256,
        "row_order_sha256": PRIOR_ROW_ORDER_SHA256,
        "rows": len(row_order),
        "embedding_shape": list(embeddings.shape),
        "embedding_dtype": str(embeddings.dtype),
        "labels_read_by_extraction": False,
        "row_order": row_order,
        "tokens": tokens,
    }


def load_prior_compatible_model(model_directory: Path):
    from transformers import AutoModel

    model, loading = AutoModel.from_pretrained(
        model_directory.as_posix(),
        local_files_only=True,
        output_loading_info=True,
    )
    missing = set(loading.get("missing_keys", []))
    unexpected = set(loading.get("unexpected_keys", []))
    mismatched = loading.get("mismatched_keys", [])
    expected_pooler = {"pooler.dense.bias", "pooler.dense.weight"}
    if missing not in (set(), expected_pooler) or unexpected not in (
        set(),
        EXPECTED_MLM_ONLY_KEYS,
    ) or mismatched:
        raise GeneformerPreflightError(
            f"prior-compatible AutoModel restore differs: {loading}"
        )
    return model, {
        "missing_keys": sorted(missing),
        "unexpected_keys": sorted(unexpected),
        "mismatched_keys": mismatched,
        "unused_pooler_excluded_from_embedding": True,
        "unused_mlm_head_excluded_from_embedding": True,
    }


def cpu_strict_restore(model_directory: Path, schema: dict[str, Any]) -> dict[str, Any]:
    import torch
    from transformers import BertForMaskedLM

    full, full_loading = BertForMaskedLM.from_pretrained(
        model_directory.as_posix(),
        local_files_only=True,
        output_loading_info=True,
    )
    if (
        full_loading.get("missing_keys")
        or full_loading.get("unexpected_keys")
        or full_loading.get("mismatched_keys")
    ):
        raise GeneformerPreflightError(f"full checkpoint restore differs: {full_loading}")
    auto, auto_loading = load_prior_compatible_model(model_directory)
    full_encoder = full.bert.state_dict()
    auto_encoder = {
        key: value for key, value in auto.state_dict().items() if not key.startswith("pooler.")
    }
    if set(full_encoder) != set(auto_encoder) or any(
        not torch.equal(full_encoder[key], auto_encoder[key]) for key in full_encoder
    ):
        raise GeneformerPreflightError("strict and prior-compatible encoder states differ")
    expected_encoder_schema = {
        name[len("bert.") :]: value
        for name, value in schema.items()
        if name.startswith("bert.")
    }
    observed_encoder_schema = {
        key: {"shape": list(value.shape), "dtype": "F32"}
        for key, value in full_encoder.items()
    }
    if observed_encoder_schema != expected_encoder_schema:
        raise GeneformerPreflightError("restored encoder schema differs")
    hidden = torch.linspace(-1.0, 1.0, steps=6 * HIDDEN_SIZE).reshape(
        2, 3, HIDDEN_SIZE
    )
    attention = torch.tensor([[1, 1, 0], [1, 1, 1]], dtype=torch.long)
    pooled = raw_mean_pool(hidden, attention)
    if pooled.shape != (2, HIDDEN_SIZE) or not torch.isfinite(pooled).all():
        raise GeneformerPreflightError("raw mean-pool compatibility probe differs")
    return {
        "full_mlm_strict_restore": True,
        "prior_auto_model_restore": auto_loading,
        "strict_to_prior_encoder_parameter_identity": True,
        "encoder_state_tensor_count": len(full_encoder),
        "raw_embedding_policy": EMBEDDING_POLICY,
        "raw_embedding_dimension": HIDDEN_SIZE,
        "labels_loaded": False,
        "label_refinement_used": False,
        "full_encoder_forward_executed": False,
    }


def raw_mean_pool(hidden: Any, attention: Any) -> Any:
    mask = attention.unsqueeze(-1).to(hidden.dtype)
    return (hidden * mask).sum(dim=1) / mask.sum(dim=1)


def gpu_one_batch(
    model_directory: Path,
    prior: dict[str, Any],
    prior_embeddings_path: Path,
    output_embeddings: Path,
) -> dict[str, Any]:
    import numpy as np
    import torch

    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise GeneformerPreflightError("exactly one CUDA device is required")
    if "L40S" not in torch.cuda.get_device_name(0):
        raise GeneformerPreflightError("one-batch compatibility probe requires L40S")
    torch.manual_seed(20260824)
    torch.cuda.manual_seed_all(20260824)
    torch.use_deterministic_algorithms(True)
    model, loading = load_prior_compatible_model(model_directory)
    model.eval().to("cuda")
    rows = prior["row_order"][:GPU_BATCH_SIZE]
    sequences = [prior["tokens"][row] for row in rows]
    width = max(len(values) for values in sequences)
    ids = torch.zeros((len(rows), width), dtype=torch.long, device="cuda")
    attention = torch.zeros_like(ids)
    for index, values in enumerate(sequences):
        ids[index, : len(values)] = torch.tensor(values, dtype=torch.long, device="cuda")
        attention[index, : len(values)] = 1
    with torch.inference_mode():
        first_hidden = model(input_ids=ids, attention_mask=attention).last_hidden_state
        first = raw_mean_pool(first_hidden, attention).cpu().numpy().astype("float32")
        second_hidden = model(input_ids=ids, attention_mask=attention).last_hidden_state
        second = raw_mean_pool(second_hidden, attention).cpu().numpy().astype("float32")
    prior_embeddings = np.load(prior_embeddings_path, mmap_mode="r", allow_pickle=False)[
        :GPU_BATCH_SIZE
    ]
    repeat_max_abs = float(np.max(np.abs(first - second)))
    prior_max_abs = float(np.max(np.abs(first - prior_embeddings)))
    prior_mean_abs = float(np.mean(np.abs(first - prior_embeddings)))
    if (
        first.shape != (GPU_BATCH_SIZE, HIDDEN_SIZE)
        or not np.isfinite(first).all()
        or repeat_max_abs > 1e-6
        or prior_max_abs > 1e-5
        or prior_mean_abs > 1e-6
    ):
        raise GeneformerPreflightError(
            "frozen-runtime batch does not reproduce prior raw embeddings"
        )
    np.save(output_embeddings, first, allow_pickle=False)
    return {
        "device": torch.cuda.get_device_name(0),
        "cuda": torch.version.cuda,
        "batch_rows": GPU_BATCH_SIZE,
        "batch_width": width,
        "embedding_shape": list(first.shape),
        "embedding_sha256": sha256(first.tobytes()).hexdigest(),
        "prior_batch_sha256": sha256(prior_embeddings.tobytes()).hexdigest(),
        "repeat_max_abs": repeat_max_abs,
        "prior_max_abs": prior_max_abs,
        "prior_mean_abs": prior_mean_abs,
        "prior_smoke_raw_embeddings_reproduced": True,
        "auto_model_restore": loading,
        "labels_loaded": False,
        "label_refinement_used": False,
    }


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--mode", choices=("cpu", "gpu"), required=True)
    value.add_argument("--bundle", type=Path, required=True)
    value.add_argument("--inventory", type=Path, required=True)
    value.add_argument("--tokens", type=Path, required=True)
    value.add_argument("--token-summary", type=Path, required=True)
    value.add_argument("--prior-embeddings", type=Path, required=True)
    value.add_argument("--prior-row-order", type=Path, required=True)
    value.add_argument("--adapter", type=Path, required=True)
    value.add_argument("--extraction-script", type=Path, required=True)
    value.add_argument("--output", type=Path, required=True)
    value.add_argument("--output-embeddings", type=Path)
    return value


def main() -> int:
    arguments = parser().parse_args()
    if arguments.output.exists():
        raise GeneformerPreflightError(f"refusing to replace output: {arguments.output}")
    runtime = validate_runtime()
    assets = validate_assets(arguments)
    schema = assets.pop("checkpoint_schema")
    model_directory = assets.pop("model_directory")
    prior = assets["prior_smoke"]
    if arguments.mode == "cpu":
        mode_result = cpu_strict_restore(model_directory, schema)
        status = "pass_cpu_runtime_and_strict_restore"
    else:
        if arguments.output_embeddings is None or arguments.output_embeddings.exists():
            raise GeneformerPreflightError("GPU output embedding path is absent or exists")
        mode_result = gpu_one_batch(
            model_directory,
            prior,
            arguments.prior_embeddings,
            arguments.output_embeddings,
        )
        status = "pass_l40s_prior_batch_reproduction"
    assets["prior_smoke"] = {
        key: value
        for key, value in prior.items()
        if key not in {"row_order", "tokens"}
    }
    receipt = {
        "schema_version": "masld-bench-geneformer-v2-316m-runtime-preflight-v1",
        "status": status,
        "mode": arguments.mode,
        "model_id": MODEL_ID,
        "embedding_policy": EMBEDDING_POLICY,
        "raw_unrefined_embeddings_only": True,
        "labels_loaded": False,
        "sealed_outcomes_loaded": False,
        "runtime": runtime,
        "assets": assets,
        "result": mode_result,
    }
    with arguments.output.open("x", encoding="utf-8") as handle:
        json.dump(receipt, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
    print(arguments.output.resolve().as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
