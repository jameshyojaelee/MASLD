#!/usr/bin/env python3
"""Extract frozen SCimilarity v1.1 embeddings without labels or reference search."""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
import os
from pathlib import Path
from typing import Any


ROWS = 50_000
INPUT_DIMENSION = 28_231
LATENT_DIMENSION = 128
ENCODER_SHA256 = "a08f023788bdbded191e2d9344e5d4240ffc7d7505eca4e0b7a37526b858d19f"
NN_MODELS_SHA256 = "5b3edf1a57ae653202d9655e049b516ec31cfa14e81daaab7c96b647fec6c110"
EMBEDDING_POLICY = "released_encoder_unit_l2_latent_without_reference_index"


class SCimilarityExtractionError(RuntimeError):
    """Raised when the frozen fixture or embedding requirement differs."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_fixture(
    root: Path,
    *,
    expected_rows: int = ROWS,
    expected_dimension: int = INPUT_DIMENSION,
) -> dict[str, Any]:
    import numpy as np
    from scipy import sparse

    receipt = json.loads((root / "fixture_receipt.json").read_text(encoding="utf-8"))
    if (
        receipt.get("status") != "pass_outcome_blind_fixture"
        or receipt.get("model_id") != "scimilarity_v1_1"
        or receipt.get("rows") != expected_rows
        or receipt.get("input_dimension") != expected_dimension
        or receipt.get("evaluation_label_columns_read")
        or receipt.get("histology_columns_read")
        or receipt.get("sealed_outcomes_read") is not False
        or receipt.get("released_reference_index_used") is not False
        or receipt.get("feature_symbol_source_observations_read") is not False
    ):
        raise SCimilarityExtractionError("outcome-blind SCimilarity fixture receipt differs")
    matrix = sparse.load_npz(root / "normalized_counts.npz").tocsr()
    observed_mask = np.load(root / "observed_gene_mask.npy", allow_pickle=False)
    study_observability = np.load(root / "study_gene_observability.npy", allow_pickle=False)
    outer_folds = np.load(root / "outer_folds.npy", allow_pickle=False)
    rows = (root / "embedding_row_order.txt").read_text(encoding="utf-8").splitlines()
    if (
        matrix.shape != (expected_rows, expected_dimension)
        or matrix.dtype != np.dtype("float32")
        or matrix.nnz < 1
        or not np.isfinite(matrix.data).all()
        or np.any(matrix.data < 0.0)
        or observed_mask.shape != (expected_dimension,)
        or observed_mask.dtype != np.dtype("bool")
        or int(observed_mask.sum()) != receipt.get("gene_overlap")
        or study_observability.shape != (receipt.get("studies"), expected_dimension)
        or study_observability.dtype != np.dtype("bool")
        or receipt.get("native_model_accepts_observability_mask") is not False
        or receipt.get("aggregate_exposure_status") != "encoder_seen"
        or outer_folds.shape != (expected_rows,)
        or outer_folds.dtype != np.dtype("int8")
        or set(map(int, outer_folds)) != set(range(5))
        or len(rows) != expected_rows
        or len(set(rows)) != expected_rows
    ):
        raise SCimilarityExtractionError("SCimilarity matrix/mask/row fixture differs")
    row_hash = sha256(("\n".join(rows) + "\n").encode("utf-8")).hexdigest()
    if row_hash != receipt.get("row_order_sha256"):
        raise SCimilarityExtractionError("SCimilarity fixture row order hash differs")
    return {
        "receipt": receipt,
        "matrix": matrix,
        "observed_mask": observed_mask,
        "study_observability": study_observability,
        "outer_folds": outer_folds,
        "rows": rows,
    }


def _load_encoder(checkpoint: Path, *, device: str):
    import torch

    if sha256_file(checkpoint) != ENCODER_SHA256:
        raise SCimilarityExtractionError("SCimilarity encoder checkpoint hash differs")
    spec = importlib.util.find_spec("scimilarity.nn_models")
    if spec is None or spec.origin is None:
        raise SCimilarityExtractionError("SCimilarity architecture module is unavailable")
    nn_models_path = Path(spec.origin).resolve(strict=True)
    if sha256_file(nn_models_path) != NN_MODELS_SHA256:
        raise SCimilarityExtractionError("SCimilarity architecture source hash differs")
    from scimilarity.nn_models import Encoder

    model = Encoder(
        n_genes=INPUT_DIMENSION,
        latent_dim=LATENT_DIMENSION,
        hidden_dim=[1024, 1024, 1024],
        dropout=0.5,
        input_dropout=0.4,
    )
    value = torch.load(checkpoint, map_location="cpu", weights_only=True, mmap=True)
    if not isinstance(value, dict) or set(value) != {"state_dict"}:
        raise SCimilarityExtractionError("SCimilarity checkpoint container differs")
    incompatibility = model.load_state_dict(value["state_dict"], strict=True)
    if incompatibility.missing_keys or incompatibility.unexpected_keys:
        raise SCimilarityExtractionError("SCimilarity checkpoint strict restore differs")
    model.eval().to(device)
    parameter_count = sum(parameter.numel() for parameter in model.parameters())
    if parameter_count != 31_146_115:
        raise SCimilarityExtractionError("SCimilarity encoder parameter count differs")
    return model, {"nn_models_path": str(nn_models_path), "parameter_count": parameter_count}


def embed_batch(model: Any, matrix: Any, *, device: str):
    import numpy as np
    import torch

    values = torch.as_tensor(matrix.toarray(), dtype=torch.float32, device=device)
    with torch.inference_mode():
        result = model(values).detach().cpu().numpy().astype(np.float32, copy=False)
    if result.ndim != 2 or result.shape[1] != LATENT_DIMENSION or not np.isfinite(result).all():
        raise SCimilarityExtractionError("SCimilarity embedding batch differs")
    norms = np.linalg.norm(result, axis=1)
    if not np.allclose(norms, 1.0, rtol=1e-5, atol=1e-6):
        raise SCimilarityExtractionError("SCimilarity embeddings are not unit L2 normalized")
    return result


def extract(
    fixture_root: Path,
    checkpoint: Path,
    output: Path,
    *,
    expected_rows: int = ROWS,
    batch_size: int = 512,
    device: str = "cuda",
) -> dict[str, Any]:
    import numpy as np
    import torch

    if output.exists() or batch_size < 1 or device not in {"cpu", "cuda"}:
        raise SCimilarityExtractionError("output/batch/device contract differs")
    if device == "cuda" and (not torch.cuda.is_available() or torch.cuda.device_count() != 1):
        raise SCimilarityExtractionError("exactly one CUDA device is required")
    fixture = load_fixture(fixture_root, expected_rows=expected_rows)
    model, loading = _load_encoder(checkpoint, device=device)
    output.mkdir(parents=True)
    temporary = output / ".embeddings.npy.incomplete"
    embeddings = np.lib.format.open_memmap(
        temporary,
        mode="w+",
        dtype="float32",
        shape=(expected_rows, LATENT_DIMENSION),
    )
    if device == "cuda":
        torch.cuda.reset_peak_memory_stats()
    for start in range(0, expected_rows, batch_size):
        stop = min(start + batch_size, expected_rows)
        embeddings[start:stop] = embed_batch(
            model, fixture["matrix"][start:stop], device=device
        )
    embeddings.flush()
    del embeddings
    os.rename(temporary, output / ".embeddings.npy.complete")
    completed = np.load(output / ".embeddings.npy.complete", mmap_mode="r", allow_pickle=False)
    representation = output / "common_embeddings.npz"
    np.savez_compressed(
        representation,
        embeddings=np.asarray(completed),
        outer_folds=fixture["outer_folds"],
        row_ids=np.asarray(fixture["rows"]),
    )
    np.save(output / "observed_gene_mask.npy", fixture["observed_mask"], allow_pickle=False)
    np.save(
        output / "study_gene_observability.npy",
        fixture["study_observability"],
        allow_pickle=False,
    )
    del completed
    (output / ".embeddings.npy.complete").unlink()
    with np.load(representation, allow_pickle=False) as representation_data:
        if set(representation_data.files) != {"embeddings", "outer_folds", "row_ids"}:
            raise SCimilarityExtractionError("common-head NPZ field contract differs")
        observed = representation_data["embeddings"]
        observed_folds = representation_data["outer_folds"]
        observed_rows = representation_data["row_ids"]
    if (
        observed.shape != (expected_rows, LATENT_DIMENSION)
        or not np.isfinite(observed).all()
        or observed_folds.shape != (expected_rows,)
        or observed_rows.shape != (expected_rows,)
    ):
        raise SCimilarityExtractionError("SCimilarity embedding output differs")
    receipt = {
        "schema_version": "masld-bench-scimilarity-v1-1-frozen-screen-extraction-v1",
        "status": "pass_outcome_blind_embeddings",
        "model_id": "scimilarity_v1_1",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "rows": expected_rows,
        "embedding_width": LATENT_DIMENSION,
        "embedding_dtype": str(observed.dtype),
        "policies": {
            "common": {
                "embedding_policy": EMBEDDING_POLICY,
                "input_normalization": fixture["receipt"]["normalization"],
                "released_reference_index_used": False,
            }
        },
        "representation_sha256": sha256_file(representation),
        "fixture_row_order_sha256": fixture["receipt"]["row_order_sha256"],
        "checkpoint_sha256": ENCODER_SHA256,
        "exposure_status": "encoder_seen",
        "development_exposure_by_study": fixture["receipt"]["exposure_by_study"],
        "study_gene_observability_order": fixture["receipt"]["study_gene_observability_order"],
        "study_gene_observability_is_global": fixture["receipt"]["study_gene_observability_is_global"],
        "observed_gene_mask_sha256": sha256_file(output / "observed_gene_mask.npy"),
        "study_gene_observability_sha256": sha256_file(output / "study_gene_observability.npy"),
        "native_model_accepts_observability_mask": False,
        "native_fixed_vocabulary_exception": fixture["receipt"]["native_fixed_vocabulary_exception"],
        "architecture_source_sha256": NN_MODELS_SHA256,
        "strict_weights_only_restore": True,
        "batch_size": batch_size,
        "device": torch.cuda.get_device_name(0) if device == "cuda" else "cpu",
        "cuda": torch.version.cuda if device == "cuda" else None,
        "cuda_peak_allocated_bytes": int(torch.cuda.max_memory_allocated()) if device == "cuda" else 0,
        "cuda_peak_reserved_bytes": int(torch.cuda.max_memory_reserved()) if device == "cuda" else 0,
        "loading": loading,
        "raw_unrefined_embeddings_only": True,
        "released_reference_index_used": False,
        "downstream_head_fit": False,
        "evaluation_label_columns_read": [],
        "histology_read": False,
        "sealed_outcomes_read": False,
        "sealed_champion_eligible": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    artifacts = {
        name: {"sha256": sha256_file(output / name), "size_bytes": (output / name).stat().st_size}
        for name in (
            "common_embeddings.npz",
            "observed_gene_mask.npy",
            "receipt.json",
            "study_gene_observability.npy",
        )
    }
    (output / "ARTIFACTS.json").write_text(
        json.dumps(
            {
                "artifacts": artifacts,
                "schema_version": "masld-bench-artifact-manifest-v1",
                "status": "pass",
            },
            sort_keys=True,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-rows", type=int, default=ROWS)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    arguments = parser.parse_args()
    result = extract(
        arguments.fixture_root,
        arguments.checkpoint,
        arguments.output,
        expected_rows=arguments.expected_rows,
        batch_size=arguments.batch_size,
        device=arguments.device,
    )
    print(json.dumps({key: result[key] for key in ("status", "rows", "embedding_width")}))


if __name__ == "__main__":
    main()
