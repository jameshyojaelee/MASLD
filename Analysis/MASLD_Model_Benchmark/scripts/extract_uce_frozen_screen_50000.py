#!/usr/bin/env python3
"""Extract UCE 4L or 33L embeddings only; common heads are a separate lane."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
from typing import Any, Mapping

from scripts.uce_frozen_screen_50000_common import (
    ARCH,
    MODEL_IDS,
    UCEActivationError,
    load_activation_contract,
    load_checkpoint_model,
    sha256_file,
    verify_record,
)


def load_fixture(
    root: Path,
    *,
    expected_artifacts_sha256: str | None = None,
    expected_rows: int = 50_000,
    expected_donors: int = 102,
    expected_studies: int = 7,
) -> dict[str, Any]:
    import numpy as np

    if expected_artifacts_sha256 is not None:
        from masld_bench.artifacts import verify_frozen_tree

        if sha256_file(root / "ARTIFACTS.json") != expected_artifacts_sha256:
            raise UCEActivationError("UCE fixture ARTIFACTS SHA-256 differs")
        metadata = verify_frozen_tree(root).get("metadata", {})
        if (
            metadata.get("artifact_class") != "uce_frozen_screen_50000_fixture"
            or metadata.get("dataset_view_id")
            != "resource_atlas_frozen_screen_50000_v1"
            or metadata.get("rows") != expected_rows
            or metadata.get("evaluation_labels_read") is not False
            or metadata.get("sealed_outcomes_read") is not False
        ):
            raise UCEActivationError("UCE frozen fixture metadata differs")
    fixture = root / "fixture"
    receipt = json.loads((fixture / "receipt.json").read_text(encoding="utf-8"))
    if (
        receipt.get("status") != "pass_outcome_blind_uce_fixture"
        or receipt.get("rows") != expected_rows
        or receipt.get("donors") != expected_donors
        or receipt.get("studies") != expected_studies
        or receipt.get("outer_folds") != 5
        or receipt.get("row_roster_preserved") is not True
        or receipt.get("cells_removed") != 0
        or receipt.get("cells_truncated") != 0
        or receipt.get("truncation_policy")
        != "fail_if_sentence_exceeds_1536;never_truncate"
        or receipt.get("evaluation_label_columns_read") != []
        or receipt.get("symbol_source_observation_fields_read") != []
        or receipt.get("sealed_outcomes_read") is not False
        or receipt.get("downstream_head_fit") is not False
        or receipt.get("head_fitting_delegated_to_hardened_common_lane") is not True
    ):
        raise UCEActivationError("UCE label-blind fixture receipt differs")
    tokens = np.load(fixture / "tokens.npy", mmap_mode="r", allow_pickle=False)
    lengths = np.load(fixture / "content_lengths.npy", allow_pickle=False)
    outer_folds = np.load(fixture / "outer_folds.npy", allow_pickle=False)
    row_ids = (fixture / "row_ids.txt").read_text(encoding="utf-8").splitlines()
    row_order_sha = sha256(("\n".join(row_ids) + "\n").encode("utf-8")).hexdigest()
    if (
        tokens.shape != (expected_rows, ARCH["pad_length"])
        or tokens.dtype != np.dtype("int32")
        or lengths.shape != (expected_rows,)
        or outer_folds.shape != (expected_rows,)
        or len(row_ids) != expected_rows
        or len(set(row_ids)) != expected_rows
        or set(map(int, outer_folds)) != set(range(5))
        or np.any(lengths < 1)
        or np.any(lengths > ARCH["pad_length"])
        or int(np.max(tokens)) >= 145469
        or int(np.min(tokens)) < 0
        or row_order_sha != receipt.get("row_order_sha256")
        or sha256_file(fixture / "tokens.npy") != receipt.get("tokens_sha256")
    ):
        raise UCEActivationError("UCE token fixture arrays differ")
    return {
        "receipt": receipt,
        "tokens": tokens,
        "lengths": lengths.astype(np.int16, copy=False),
        "outer_folds": outer_folds.astype(np.int8, copy=False),
        "row_ids": row_ids,
    }


def embed_batch(model: Any, tokens: Any, lengths: Any, *, device: str):
    import numpy as np
    import torch

    token_tensor = torch.as_tensor(np.asarray(tokens, dtype=np.int64), device=device)
    length_tensor = torch.as_tensor(np.asarray(lengths, dtype=np.int64), device=device)
    positions = torch.arange(token_tensor.shape[1], device=device).unsqueeze(0)
    mask = (positions < length_tensor.unsqueeze(1)).to(dtype=torch.float32)
    with torch.inference_mode():
        values = model(token_tensor.transpose(0, 1), mask)
    result = values.detach().cpu().numpy().astype(np.float32, copy=False)
    if (
        result.ndim != 2
        or result.shape[1] != ARCH["output_dim"]
        or not np.isfinite(result).all()
    ):
        raise UCEActivationError("UCE embedding batch shape or values differ")
    norms = np.linalg.norm(result, axis=1)
    if not np.allclose(norms, 1.0, rtol=1e-5, atol=1e-6):
        raise UCEActivationError("UCE embeddings are not unit L2 normalized")
    return result


def extract(
    contract: Mapping[str, Any],
    fixture_root: Path,
    output: Path,
    *,
    model_id: str,
    expected_fixture_artifacts_sha256: str,
    batch_size: int,
    device: str = "cuda",
    expected_rows: int = 50_000,
) -> dict[str, Any]:
    import numpy as np
    import torch

    if output.exists() or model_id not in MODEL_IDS or batch_size < 1:
        raise UCEActivationError("UCE extraction output, model, or batch contract differs")
    if device != "cuda" or not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise UCEActivationError("UCE full extraction requires exactly one CUDA device")
    device_name = torch.cuda.get_device_name(0)
    if "L40S" not in device_name.upper():
        raise UCEActivationError("UCE activation is admitted only on the frozen L40S runtime")
    fixture = load_fixture(
        fixture_root,
        expected_artifacts_sha256=expected_fixture_artifacts_sha256,
        expected_rows=expected_rows,
    )
    checkpoint_record = contract["checkpoints"][model_id]
    checkpoint = verify_record(contract, checkpoint_record, f"{model_id} checkpoint")
    all_tokens = verify_record(
        contract, contract["model_files"]["all_tokens"], "UCE all-token table"
    )
    model, loading = load_checkpoint_model(
        checkpoint, all_tokens, model_id=model_id, device=device
    )
    torch.cuda.reset_peak_memory_stats()

    stage = output.with_name(f".{output.name}.{os.getpid()}.staging")
    if stage.exists():
        raise UCEActivationError("UCE extraction staging path exists")
    embedding_dir = stage / "embeddings"
    embedding_dir.mkdir(parents=True)
    raw_path = embedding_dir / ".common_embeddings.npy.incomplete"
    embeddings = np.lib.format.open_memmap(
        raw_path,
        mode="w+",
        dtype="float32",
        shape=(expected_rows, ARCH["output_dim"]),
    )
    for start in range(0, expected_rows, batch_size):
        stop = min(start + batch_size, expected_rows)
        embeddings[start:stop] = embed_batch(
            model,
            fixture["tokens"][start:stop],
            fixture["lengths"][start:stop],
            device=device,
        )
    embeddings.flush()
    repeated = embed_batch(
        model,
        fixture["tokens"][: min(4, expected_rows)],
        fixture["lengths"][: min(4, expected_rows)],
        device=device,
    )
    if not np.allclose(
        repeated,
        np.asarray(embeddings[: repeated.shape[0]]),
        rtol=0.0,
        atol=1e-6,
    ):
        raise UCEActivationError("UCE deterministic repeat probe differs")
    norms = np.linalg.norm(embeddings, axis=1)
    if not np.allclose(norms, 1.0, rtol=1e-5, atol=1e-6):
        raise UCEActivationError("UCE full embeddings are not unit normalized")
    final_npz = embedding_dir / "common_embeddings.npz"
    temporary_npz = embedding_dir / ".common_embeddings.npz.incomplete"
    with temporary_npz.open("wb") as handle:
        np.savez(
            handle,
            embeddings=np.asarray(embeddings),
            outer_folds=fixture["outer_folds"],
            row_ids=np.asarray(fixture["row_ids"]),
        )
    os.rename(temporary_npz, final_npz)
    del embeddings
    raw_path.unlink()
    with np.load(final_npz, allow_pickle=False) as bundle:
        if set(bundle.files) != {"embeddings", "outer_folds", "row_ids"}:
            raise UCEActivationError("UCE common-head NPZ fields differ")
        if (
            bundle["embeddings"].shape != (expected_rows, ARCH["output_dim"])
            or bundle["outer_folds"].tolist() != fixture["outer_folds"].tolist()
            or bundle["row_ids"].astype(str).tolist() != fixture["row_ids"]
        ):
            raise UCEActivationError("UCE common-head NPZ content differs")
    receipt = {
        "schema_version": "masld-bench-uce-frozen-screen-extraction-v1",
        "status": "pass_raw_outcome_blind_extraction",
        "model_id": model_id,
        "checkpoint_sha256": checkpoint_record["sha256"],
        "checkpoint_layers": checkpoint_record["layers"],
        "exposure_status": checkpoint_record["exposure_status"],
        "sealed_champion_eligible": checkpoint_record["sealed_champion_eligible"],
        "known_encoder_seen_development_studies": contract["exposure"][
            "known_encoder_seen_development_studies"
        ],
        "source_stratified_diagnostics_required": True,
        "dataset_view_id": contract["dataset_view_id"],
        "split_id": contract["split_id"],
        "rows": expected_rows,
        "embedding_width": ARCH["output_dim"],
        "embedding_dtype": "float32",
        "policies": {
            "common": "l2_normalized_cls_decoder_output_from_order_invariant_official_format_sentences"
        },
        "embedding_sha256": sha256_file(final_npz),
        "fixture_artifacts_sha256": expected_fixture_artifacts_sha256,
        "fixture_row_order_sha256": fixture["receipt"]["row_order_sha256"],
        "strict_weights_only_restore": loading["strict_weights_only_restore"],
        "shared_token_table_equal": loading["shared_token_table_equal"],
        "loading": loading,
        "batch_size": batch_size,
        "device": device_name,
        "cuda": torch.version.cuda,
        "cuda_peak_allocated_bytes": int(torch.cuda.max_memory_allocated()),
        "cuda_peak_reserved_bytes": int(torch.cuda.max_memory_reserved()),
        "deterministic_repeat_probe_rows": int(repeated.shape[0]),
        "deterministic_repeat_atol": 1e-6,
        "raw_unrefined_embeddings_only": True,
        "downstream_head_fit": False,
        "head_fitting_delegated_to_hardened_common_lane": True,
        "evaluation_label_columns_read": [],
        "histology_columns_read": [],
        "sealed_outcomes_read": False,
    }
    (embedding_dir / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    os.rename(stage, output)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--fixture-root", type=Path, required=True)
    parser.add_argument("--expected-fixture-artifacts-sha256", required=True)
    parser.add_argument("--model-id", choices=MODEL_IDS, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, required=True)
    arguments = parser.parse_args()
    contract = load_activation_contract(arguments.contract, arguments.project_root)
    result = extract(
        contract,
        arguments.fixture_root,
        arguments.output,
        model_id=arguments.model_id,
        expected_fixture_artifacts_sha256=arguments.expected_fixture_artifacts_sha256,
        batch_size=arguments.batch_size,
    )
    print(
        json.dumps(
            {key: result[key] for key in ("status", "model_id", "rows", "embedding_width")}
        )
    )


if __name__ == "__main__":
    main()
