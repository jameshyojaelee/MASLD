#!/usr/bin/env python
"""scGPT common-lane adapter for the frozen cell-state task.

Mirrors ``geneformer_common_lane.py``: a frozen encoder, a trainable common
head fit on outer-training donors only, and byte-identical output column
naming.  The encoder is scGPT (``bowang-lab/scGPT`` at the pinned revision
``cebd6fae655b9c585a4807daa3ac31bb764f06b4``), reimplemented here in plain
PyTorch rather than imported from the upstream ``scgpt`` package, which is
not installed in this environment (it depends on flash-attn) and is not
needed: the checkpoint's ``fast_transformer_backend="flash"`` attention is
mathematically ordinary scaled-dot-product self-attention over a fused QKV
projection (``self_attn.Wqkv`` / ``self_attn.out_proj``); FlashAttention is a
kernel-level optimization of the identical computation, not a different
architecture.  The tokenization/embedding contract below (rank-like 51-bin
digitization of nonzero, in-vocabulary genes; CLS-prepend; sample-without-
replacement truncation to 1200 tokens; L2-normalized CLS hidden state) is
taken from ``config/artifacts/models/scgpt/checkpoints.json``'s
``input_contract`` and cross-checked against the pinned repository's
``scgpt/preprocess.py`` (``binning``/``_digitize``), ``scgpt/data_collator.py``
(``DataCollator``, ``n_bins=51``, ``keep_first_n_tokens=1``), and
``scgpt/tasks/cell_emb.py`` (``embed_data`` / ``get_batch_cell_embeddings``),
read directly from GitHub at the pinned revision.

Restricted deserialization.  The upstream ``best_model.pt`` is untrusted and
is never handed to a bare, unrestricted ``pickle``/``torch.load``.  It is
first verified byte-for-byte against the acquisition's own ``SHA256SUMS``
manifest (the analogue of ``checkpoint_preflight.py``'s hash-then-load
discipline; there is no scGPT bundle-staging module, so this adapter performs
that verification itself), and only then loaded with
``torch.load(..., weights_only=True)`` -- PyTorch's own restricted unpickler,
which admits only tensor/storage reconstruction and builtin containers, no
arbitrary global or class construction.  The Geneformer common-lane adapter
already relies on this same mechanism for its own ``head_state.pt``.  Loading
this checkpoint under ``weights_only=True`` was verified to succeed with the
allowlist untouched: the file is a plain ``collections.OrderedDict`` of
``torch.float32`` tensors, so nothing was widened.

Two upstream steps use an *unseeded* global RNG (the per-cell bin-tie-break
randomization in ``_digitize``, and ``torch.randperm`` truncation sampling for
cells with more than 1199 expressed in-vocabulary genes).  Per the frozen
``sampling_warning`` in the registry contract, this adapter replaces both
draws with a per-row generator seeded from the row's own ``row_hash``, so
tokenization is deterministic and invariant to cell processing order; native
(non-common-lane) scGPT runs are expected to differ and are reported
separately.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from hvg_pca_logistic import (  # type: ignore[import-not-found]
    ScientificAdapterError,
    _artifact_record,
    _input_by_role,
    _require_sha256,
    _sha256_file,
    _validate_file_artifact,
    _write_json,
    _write_tsv,
    fold_index,
    join_hash,
)

TASK_ID = "cell_state_mapping"
DATASET_ID = "resource_atlas_current"
RUNTIME_ID = "gpu_scgpt"
RECEIPT_SCHEMA = "masld-bench-adapter-receipt-v1"
ADAPTER_ID = "scgpt_common_lane_v1"
MODEL_IDS = ("scgpt_continual",)
HEAD_IDS = ("linear", "two_layer_mlp")
EMBEDDING_POLICY = "cls_token_hidden_state_l2_normalized"

# Frozen architecture contract -- from the acquisition-time exchange with the
# team lead and confirmed against config/artifacts/models/scgpt/checkpoints.json.
EMBEDDING_SIZE = 512
D_HID = 512
N_LAYERS = 12
N_HEADS = 8
N_EXPRESSION_BINS = 51
MAX_SEQUENCE_LENGTH = 1200
PAD_TOKEN = "<pad>"
PAD_VALUE = -2.0
CLS_TOKEN = "<cls>"
VOCAB_SIZE = 60697
PAD_ID = 60694
CLS_ID = 60695

GENE_SYMBOL_COLUMN = "source_feature_id"
_REQUIRED_CHECKPOINT_FILES = ("args.json", "vocab.json", "best_model.pt")
# Submodules present in the raw upstream checkpoint that this adapter never
# loads: the masked-value/MLM decoder, the classification and MVC heads, and
# an unused ``flag_encoder``.  None of them appear in scGPT's own zero-shot
# embedding path (``scgpt.tasks.cell_emb.get_batch_cell_embeddings``, which
# calls ``model._encode`` directly and never touches ``model.decoder``,
# ``model.cls_decoder``, ``model.mvc_decoder``, or a ``flag_encoder`` --
# confirmed absent from every class in ``scgpt/model/model.py`` and
# ``scgpt/model/generation_model.py`` at the pinned revision).
_DROPPED_CHECKPOINT_PREFIXES = (
    "flag_encoder.",
    "decoder.",
    "cls_decoder.",
    "mvc_decoder.",
)


def _select_device(requested: str) -> str:
    import torch

    if requested == "cpu":
        return "cpu"
    if torch.cuda.is_available():
        return "cuda"
    if requested == "cuda":
        raise ScientificAdapterError("CUDA was required but is unavailable")
    return "cpu"


def _parse_sha256sums(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            raise ScientificAdapterError("scGPT checksum manifest line is malformed")
        digest, relative = parts[0].strip().lower(), parts[1].strip()
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ScientificAdapterError("scGPT checksum manifest digest is not a SHA-256")
        mapping[relative] = digest
    if not mapping:
        raise ScientificAdapterError("scGPT checksum manifest is empty")
    return mapping


def _verify_checkpoint_file(path: Path, expected_sha256: str, label: str) -> Path:
    if not path.is_absolute() or path.is_symlink() or not path.is_file():
        raise ScientificAdapterError(f"{label} must be an absolute regular file")
    if _sha256_file(path) != expected_sha256:
        raise ScientificAdapterError(f"{label} SHA-256 differs from the acquisition manifest")
    return path


def _validate_architecture(args: Mapping[str, Any], vocab: Mapping[str, Any]) -> None:
    expected = {
        "embsize": EMBEDDING_SIZE,
        "d_hid": D_HID,
        "nlayers": N_LAYERS,
        "nheads": N_HEADS,
        "n_bins": N_EXPRESSION_BINS,
        "max_seq_len": MAX_SEQUENCE_LENGTH,
        "pad_token": PAD_TOKEN,
        "pad_value": int(PAD_VALUE),
    }
    if any(args.get(field) != value for field, value in expected.items()):
        raise ScientificAdapterError(
            "scGPT checkpoint args.json differs from the frozen architecture contract"
        )
    if len(vocab) != VOCAB_SIZE:
        raise ScientificAdapterError(
            "scGPT vocabulary size differs from the frozen architecture contract"
        )


def _verified_checkpoint(request: Mapping[str, Any], model_id: str) -> dict[str, Any]:
    """Verify the acquired-and-hashed scGPT bundle and load its plain-data files.

    The request carries one input with role ``model_checkpoint_bundle:<model_id>``
    pointing at the acquisition's ``SHA256SUMS`` file.  Every required checkpoint
    file is re-verified against it directly (size is implied by the digest;
    symlinks are rejected) before anything is read.
    """

    role = f"model_checkpoint_bundle:{model_id}"
    manifest_artifact = _input_by_role(request, role)
    manifest_path = _validate_file_artifact(manifest_artifact, "scGPT checksum manifest")
    if manifest_path.name != "SHA256SUMS":
        raise ScientificAdapterError(
            "scGPT checkpoint role must reference the acquisition SHA256SUMS manifest"
        )
    checksums = _parse_sha256sums(manifest_path)
    root = manifest_path.parent
    paths: dict[str, Path] = {}
    for name in _REQUIRED_CHECKPOINT_FILES:
        relative = f"{model_id}/{name}"
        expected = checksums.get(relative)
        if expected is None:
            raise ScientificAdapterError(f"scGPT checksum manifest is missing {relative}")
        paths[name] = _verify_checkpoint_file(root / relative, expected, f"scGPT {relative}")

    try:
        args = json.loads(paths["args.json"].read_text(encoding="utf-8"))
        vocab = json.loads(paths["vocab.json"].read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificAdapterError(f"scGPT checkpoint metadata is invalid: {error}") from error
    if not isinstance(args, Mapping) or not isinstance(vocab, Mapping) or not vocab:
        raise ScientificAdapterError("scGPT args.json/vocab.json are invalid")
    _validate_architecture(args, vocab)
    vocab = {str(key): int(value) for key, value in vocab.items()}
    if vocab.get(PAD_TOKEN) != PAD_ID or vocab.get(CLS_TOKEN) != CLS_ID:
        raise ScientificAdapterError(
            "scGPT vocabulary special-token ids differ from the frozen contract"
        )
    return {
        "weights_path": paths["best_model.pt"],
        "args": dict(args),
        "vocab": vocab,
        "pad_id": PAD_ID,
        "cls_id": CLS_ID,
    }


def _row_rng(row_hash: str):
    import numpy as np

    seed = int.from_bytes(bytes.fromhex(row_hash)[:8], "big")
    return np.random.default_rng(seed)


def _digitize(values: Any, bins: Any, rng: Any) -> Any:
    import numpy as np

    left = np.digitize(values, bins)
    right = np.digitize(values, bins, right=True)
    rands = rng.random(len(values))
    digits = rands * (right - left) + left
    return np.ceil(digits).astype(np.int64)


def _rank_bin(values: Any, n_bins: int, rng: Any) -> Any:
    """Reproduce upstream ``scgpt.preprocess.binning`` for one cell's nonzero
    values, with the tie-break draw seeded per row instead of upstream's
    unseeded global RNG (see module docstring)."""

    import numpy as np

    edges = np.quantile(values, np.linspace(0.0, 1.0, n_bins - 1))
    return _digitize(values, edges, rng)


def _truncate_indices(n_usable: int, max_len: int, rng: Any) -> Any:
    import numpy as np

    limit = max_len - 1  # position 0 is always reserved for <cls>
    if n_usable <= limit:
        return np.arange(n_usable)
    return np.sort(rng.choice(n_usable, size=limit, replace=False))


def _encode_rows(
    matrix: Any,
    gene_symbols: Sequence[str],
    vocab: Mapping[str, int],
    *,
    cls_id: int,
    n_bins: int,
    max_len: int,
    row_hashes: Sequence[str],
) -> tuple[list[list[int]], list[list[float]], dict[str, Any]]:
    import numpy as np

    token_id_by_column = np.array(
        [vocab.get(symbol, -1) for symbol in gene_symbols], dtype=np.int64
    )
    genes_total = len(gene_symbols)
    genes_in_vocab = int((token_id_by_column >= 0).sum())

    tokens: list[list[int]] = []
    values: list[list[float]] = []
    truncated = 0
    nonzero_total = 0
    nonzero_in_vocab = 0
    for row_index in range(matrix.shape[0]):
        row = matrix.getrow(row_index)
        columns = row.indices
        data = row.data
        keep = data > 0.0
        columns = columns[keep]
        data = data[keep]
        if data.size and (not np.all(np.isfinite(data)) or np.any(data < 0.0)):
            raise ScientificAdapterError("Atlas matrix has non-finite or negative counts")
        nonzero_total += int(columns.size)
        row_token_ids = token_id_by_column[columns]
        in_vocab = row_token_ids >= 0
        nonzero_in_vocab += int(in_vocab.sum())
        row_token_ids = row_token_ids[in_vocab]
        row_values = data[in_vocab].astype(np.float64)
        if row_token_ids.size == 0:
            raise ScientificAdapterError(
                "a cell has no genes usable by the scGPT vocabulary"
            )
        rng = _row_rng(row_hashes[row_index])
        bins = _rank_bin(row_values, n_bins, rng)
        if row_token_ids.size > max_len - 1:
            truncated += 1
            keep_idx = _truncate_indices(row_token_ids.size, max_len, rng)
            row_token_ids = row_token_ids[keep_idx]
            bins = bins[keep_idx]
        tokens.append([int(cls_id)] + [int(v) for v in row_token_ids])
        values.append([0.0] + [float(v) for v in bins])

    summary = {
        "schema_version": "masld-bench-scgpt-tokenization-v1",
        "embedding_policy": EMBEDDING_POLICY,
        "gene_symbol_column": GENE_SYMBOL_COLUMN,
        "within_cell_expression_transform": "51-bin rank-like binning of nonzero values",
        "sampling_contract": "preserve_cls_then_sample_without_replacement_to_1200_tokens",
        "sampling_determinism": "per_row_hash_seeded_rng_replaces_unseeded_upstream_rng",
        "genes_total": genes_total,
        "genes_in_vocab": genes_in_vocab,
        "genes_in_vocab_fraction": genes_in_vocab / genes_total if genes_total else 0.0,
        "nonzero_genes_total": nonzero_total,
        "nonzero_genes_in_vocab": nonzero_in_vocab,
        "nonzero_genes_in_vocab_fraction": (
            nonzero_in_vocab / nonzero_total if nonzero_total else 0.0
        ),
        "max_tokens": max_len,
        "cells_truncated": truncated,
        "cells": len(tokens),
    }
    return tokens, values, summary


def _build_model(checkpoint: Mapping[str, Any], device: str):
    import torch
    from torch import nn

    class GeneEncoder(nn.Module):
        def __init__(self, vocab_size: int, d_model: int, pad_id: int) -> None:
            super().__init__()
            self.embedding = nn.Embedding(vocab_size, d_model, padding_idx=pad_id)
            self.enc_norm = nn.LayerNorm(d_model)

        def forward(self, x: Any) -> Any:
            return self.enc_norm(self.embedding(x))

    class ContinuousValueEncoder(nn.Module):
        def __init__(self, d_model: int, max_value: float = 512.0) -> None:
            super().__init__()
            self.linear1 = nn.Linear(1, d_model)
            self.linear2 = nn.Linear(d_model, d_model)
            self.norm = nn.LayerNorm(d_model)
            self.max_value = max_value

        def forward(self, x: Any) -> Any:
            x = x.unsqueeze(-1).clamp(max=self.max_value)
            x = torch.relu(self.linear1(x))
            x = self.linear2(x)
            return self.norm(x)

    class FlashStyleSelfAttention(nn.Module):
        """Ordinary scaled-dot-product multi-head attention over a fused QKV
        projection -- the exact math FlashAttention accelerates, run here as
        a plain (non-fused-kernel) computation since flash-attn is not
        installed in this environment."""

        def __init__(self, d_model: int, n_heads: int) -> None:
            super().__init__()
            if d_model % n_heads != 0:
                raise ScientificAdapterError(
                    "scGPT embedding size must divide evenly by head count"
                )
            self.n_heads = n_heads
            self.head_dim = d_model // n_heads
            self.Wqkv = nn.Linear(d_model, 3 * d_model)
            self.out_proj = nn.Linear(d_model, d_model)

        def forward(self, x: Any, key_padding_mask: Any) -> Any:
            batch, length, dim = x.shape
            qkv = (
                self.Wqkv(x)
                .view(batch, length, 3, self.n_heads, self.head_dim)
                .permute(2, 0, 3, 1, 4)
            )
            query, key, value = qkv[0], qkv[1], qkv[2]
            scores = torch.matmul(query, key.transpose(-2, -1)) / (self.head_dim**0.5)
            if key_padding_mask is not None:
                scores = scores.masked_fill(
                    key_padding_mask[:, None, None, :], float("-inf")
                )
            weights = torch.softmax(scores, dim=-1)
            out = torch.matmul(weights, value)
            out = out.transpose(1, 2).reshape(batch, length, dim)
            return self.out_proj(out)

    class ScGptEncoderLayer(nn.Module):
        """Post-norm encoder layer, matching upstream
        ``FlashTransformerEncoderLayer``'s default ``norm_scheme="post"``."""

        def __init__(self, d_model: int, n_heads: int, d_hid: int) -> None:
            super().__init__()
            self.self_attn = FlashStyleSelfAttention(d_model, n_heads)
            self.linear1 = nn.Linear(d_model, d_hid)
            self.linear2 = nn.Linear(d_hid, d_model)
            self.norm1 = nn.LayerNorm(d_model)
            self.norm2 = nn.LayerNorm(d_model)

        def forward(self, x: Any, key_padding_mask: Any) -> Any:
            x = x + self.self_attn(x, key_padding_mask)
            x = self.norm1(x)
            x = x + self.linear2(torch.relu(self.linear1(x)))
            x = self.norm2(x)
            return x

    class ScGptTransformerEncoder(nn.Module):
        def __init__(self, d_model: int, n_heads: int, d_hid: int, n_layers: int) -> None:
            super().__init__()
            self.layers = nn.ModuleList(
                ScGptEncoderLayer(d_model, n_heads, d_hid) for _ in range(n_layers)
            )

        def forward(self, x: Any, key_padding_mask: Any) -> Any:
            for layer in self.layers:
                x = layer(x, key_padding_mask)
            return x

    class ScGptEncoderStack(nn.Module):
        def __init__(
            self,
            vocab_size: int,
            d_model: int,
            n_heads: int,
            d_hid: int,
            n_layers: int,
            pad_id: int,
        ) -> None:
            super().__init__()
            self.encoder = GeneEncoder(vocab_size, d_model, pad_id)
            self.value_encoder = ContinuousValueEncoder(d_model)
            self.transformer_encoder = ScGptTransformerEncoder(d_model, n_heads, d_hid, n_layers)

        def forward(self, gene_ids: Any, values: Any, key_padding_mask: Any) -> Any:
            total = self.encoder(gene_ids) + self.value_encoder(values)
            return self.transformer_encoder(total, key_padding_mask)

    model = ScGptEncoderStack(
        vocab_size=len(checkpoint["vocab"]),
        d_model=EMBEDDING_SIZE,
        n_heads=N_HEADS,
        d_hid=D_HID,
        n_layers=N_LAYERS,
        pad_id=checkpoint["pad_id"],
    )
    state_dict = torch.load(
        checkpoint["weights_path"], map_location="cpu", weights_only=True
    )
    if not isinstance(state_dict, Mapping):
        raise ScientificAdapterError("scGPT checkpoint is not a plain state dict")
    result = model.load_state_dict(dict(state_dict), strict=False)
    if result.missing_keys:
        raise ScientificAdapterError(
            "scGPT checkpoint is missing required weights: "
            + ", ".join(sorted(result.missing_keys)[:5])
        )
    unexpected = [
        key
        for key in result.unexpected_keys
        if not key.startswith(_DROPPED_CHECKPOINT_PREFIXES)
    ]
    if unexpected:
        raise ScientificAdapterError(
            "scGPT checkpoint has unrecognized weights: " + ", ".join(sorted(unexpected)[:5])
        )
    model.eval()
    model.to(device)
    return model


def _pad_batch(
    token_chunk: Sequence[Sequence[int]],
    value_chunk: Sequence[Sequence[float]],
    pad_id: int,
    pad_value: float,
):
    import torch

    width = max(len(item) for item in token_chunk)
    ids = torch.full((len(token_chunk), width), pad_id, dtype=torch.long)
    vals = torch.full((len(token_chunk), width), pad_value, dtype=torch.float32)
    mask = torch.ones((len(token_chunk), width), dtype=torch.bool)  # True = padding
    for index, (item_ids, item_vals) in enumerate(zip(token_chunk, value_chunk)):
        length = len(item_ids)
        ids[index, :length] = torch.tensor(item_ids, dtype=torch.long)
        vals[index, :length] = torch.tensor(item_vals, dtype=torch.float32)
        mask[index, :length] = False
    return ids, vals, mask


def _embed_rows(
    checkpoint: Mapping[str, Any],
    token_lists: Sequence[Sequence[int]],
    value_lists: Sequence[Sequence[float]],
    device: str,
    batch_size: int,
):
    import torch

    model = _build_model(checkpoint, device)
    outputs = []
    with torch.no_grad():
        for start in range(0, len(token_lists), batch_size):
            token_chunk = token_lists[start : start + batch_size]
            value_chunk = value_lists[start : start + batch_size]
            ids, vals, mask = _pad_batch(
                token_chunk, value_chunk, checkpoint["pad_id"], PAD_VALUE
            )
            hidden = model(ids.to(device), vals.to(device), mask.to(device))
            cls_hidden = hidden[:, 0, :]
            normalized = cls_hidden / cls_hidden.norm(p=2, dim=1, keepdim=True).clamp_min(1e-12)
            outputs.append(normalized.detach().cpu())
    return torch.cat(outputs, dim=0)


def _build_head(head_id: str, dimensions: int, classes: int, seed: int):
    import torch
    from torch import nn

    torch.manual_seed(seed)
    if head_id == "linear":
        return nn.Linear(dimensions, classes)
    if head_id == "two_layer_mlp":
        return nn.Sequential(
            nn.Linear(dimensions, 256), nn.ReLU(), nn.Linear(256, classes)
        )
    raise ScientificAdapterError(f"unsupported common head: {head_id}")


def _validate_request(request: Mapping[str, Any], action: str):
    if request.get("schema_version") != "masld-bench-adapter-request-v1":
        raise ScientificAdapterError("unsupported adapter request schema")
    if request.get("action") != action:
        raise ScientificAdapterError("request action differs")
    _require_sha256(request.get("run_id"), "run_id")
    run_spec = request.get("run_spec")
    if not isinstance(run_spec, Mapping):
        raise ScientificAdapterError("request has no RunSpec")
    model_id = str(run_spec.get("model_id", ""))
    if model_id not in MODEL_IDS:
        raise ScientificAdapterError("RunSpec names an unsupported scGPT model")
    for field, expected in (
        ("task_id", TASK_ID),
        ("split_id", "donor_outer"),
        ("adaptation_regime", "common_lane"),
        ("runtime_id", RUNTIME_ID),
    ):
        if run_spec.get(field) != expected:
            raise ScientificAdapterError(f"RunSpec {field} differs from {expected}")
    parameters = run_spec.get("hyperparameters")
    if not isinstance(parameters, Mapping):
        raise ScientificAdapterError("RunSpec carries no frozen hyperparameters")
    if parameters.get("embedding_policy") != EMBEDDING_POLICY:
        raise ScientificAdapterError("embedding policy differs from the frozen policy")
    head_id = str(parameters.get("common_head_id", ""))
    if head_id not in HEAD_IDS:
        raise ScientificAdapterError(f"unsupported common head: {head_id!r}")
    roster = parameters.get("cell_state_roster")
    if not isinstance(roster, (list, tuple)) or not roster:
        raise ScientificAdapterError("frozen cell-state roster is missing")
    roster = tuple(str(item) for item in roster)
    if list(roster) != sorted(set(roster)):
        raise ScientificAdapterError("cell-state roster must be sorted and unique")
    return parameters, roster, model_id, head_id


def _receipt(*, action, request, output, artifacts, extra_metadata) -> None:
    receipt = {
        "schema_version": RECEIPT_SCHEMA,
        "action": action,
        "run_id": request["run_id"],
        "status": "complete",
        "artifacts": [_artifact_record(p, relative_to=output) for p in artifacts],
        "metadata": {
            "adapter": ADAPTER_ID,
            "runtime_id": RUNTIME_ID,
            "fit_dataset_ids": list(request["fit_dataset_ids"]),
            **dict(extra_metadata),
        },
    }
    _write_json(output / "adapter_receipt.json", receipt)


def _atlas_path(request: Mapping[str, Any]) -> Path:
    role = "dataset_view_data:resource_atlas_geneformer_smoke_1000_v1"
    return _validate_file_artifact(_input_by_role(request, role), "Atlas smoke H5AD")


def _load_atlas(path: Path) -> Any:
    import anndata

    adata = anndata.read_h5ad(path)
    if not {"donor_id", "broad_label"}.issubset(adata.obs.columns):
        raise ScientificAdapterError("Atlas H5AD lacks required observation fields")
    if GENE_SYMBOL_COLUMN not in adata.var.columns:
        raise ScientificAdapterError("Atlas H5AD lacks the gene-symbol feature axis")
    if len(set(map(str, adata.obs_names))) != adata.n_obs:
        raise ScientificAdapterError("Atlas row identifiers are not unique")
    return adata


def _counts_csr(adata: Any):
    from scipy import sparse

    matrix = adata.X
    return matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)


def _split_rows(adata: Any, parameters: Mapping[str, Any], run_spec: Mapping[str, Any]):
    outer_folds = int(parameters["outer_folds"])
    fold = int(run_spec["fold"])
    if not 0 <= fold < outer_folds:
        raise ScientificAdapterError("outer-fold contract differs")
    namespace = str(parameters["join_namespace"])
    rows = []
    for row_id, donor in zip(
        map(str, adata.obs_names), map(str, adata.obs["donor_id"])
    ):
        assigned = fold_index(donor, seed=20260821, outer_folds=outer_folds)
        rows.append(
            {
                "row_hash": join_hash(namespace, "row", row_id),
                "unit_hash": join_hash(namespace, "unit", donor),
                "fold": assigned,
                "held_out": assigned == fold,
            }
        )
    rows.sort(key=lambda item: item["row_hash"])
    if len({r["row_hash"] for r in rows}) != len(rows):
        raise ScientificAdapterError("row identifiers collide")
    return rows


def _prior_output(request: Mapping[str, Any], action: str) -> Path:
    outputs = request.get("prior_action_outputs")
    if not isinstance(outputs, Mapping) or action not in outputs:
        raise ScientificAdapterError(f"request lacks the frozen {action} output")
    return Path(str(outputs[action])).resolve(strict=True)


def _read_split(path: Path) -> list[dict[str, Any]]:
    import csv

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [
        {
            "row_hash": r["row_hash"],
            "unit_hash": r["unit_hash"],
            "held_out": r["held_out"] == "true",
        }
        for r in rows
    ]


def _read_tokens(path: Path) -> tuple[dict[str, list[int]], dict[str, list[float]]]:
    tokens: dict[str, list[int]] = {}
    values: dict[str, list[float]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            row_hash = str(record["row_hash"])
            tokens[row_hash] = [int(v) for v in record["tokens"]]
            values[row_hash] = [float(v) for v in record["values"]]
    return tokens, values


def _labels_by_hash(request, parameters, roster):
    adata = _load_atlas(_atlas_path(request))
    namespace = str(parameters["join_namespace"])
    index = {label: position for position, label in enumerate(roster)}
    return {
        join_hash(namespace, "row", str(name)): index[str(label)]
        for name, label in zip(adata.obs_names, adata.obs["broad_label"])
    }


def prepare(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    parameters, roster, model_id, _ = _validate_request(request, "prepare")
    adata = _load_atlas(_atlas_path(request))
    observed = set(map(str, adata.obs["broad_label"]))
    if observed != set(roster):
        raise ScientificAdapterError("Atlas label set differs from the frozen roster")
    rows = _split_rows(adata, parameters, request["run_spec"])
    split_path = output / "split_rows.tsv"
    _write_tsv(
        split_path,
        ("row_hash", "unit_hash", "fold", "held_out"),
        (
            {**row, "held_out": "true" if row["held_out"] else "false"}
            for row in rows
        ),
    )

    checkpoint = _verified_checkpoint(request, model_id)
    gene_symbols = [str(g) for g in adata.var[GENE_SYMBOL_COLUMN]]
    matrix = _counts_csr(adata)
    order = {str(name): index for index, name in enumerate(adata.obs_names)}
    namespace = str(parameters["join_namespace"])
    row_hash_by_index = [""] * adata.n_obs
    for name, index in order.items():
        row_hash_by_index[index] = join_hash(namespace, "row", name)

    tokens, values, summary = _encode_rows(
        matrix,
        gene_symbols,
        checkpoint["vocab"],
        cls_id=checkpoint["cls_id"],
        n_bins=N_EXPRESSION_BINS,
        max_len=MAX_SEQUENCE_LENGTH,
        row_hashes=row_hash_by_index,
    )
    by_hash = {
        row_hash_by_index[index]: (tokens[index], values[index])
        for index in range(adata.n_obs)
    }
    tokens_path = output / "rank_value_tokens.jsonl"
    with tokens_path.open("w", encoding="utf-8") as handle:
        for row in rows:
            row_tokens, row_values = by_hash[row["row_hash"]]
            handle.write(
                json.dumps(
                    {
                        "row_hash": row["row_hash"],
                        "tokens": row_tokens,
                        "values": row_values,
                    },
                    sort_keys=True,
                )
                + "\n"
            )
    summary_path = output / "tokenization_summary.json"
    _write_json(summary_path, summary)
    _receipt(
        action="prepare",
        request=request,
        output=output,
        artifacts=[split_path, tokens_path, summary_path],
        extra_metadata={
            "model_id": model_id,
            "genes_in_vocab_fraction": summary["genes_in_vocab_fraction"],
            "nonzero_genes_in_vocab_fraction": summary["nonzero_genes_in_vocab_fraction"],
            "cells_truncated": summary["cells_truncated"],
        },
    )


def fit(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import torch

    parameters, roster, model_id, head_id = _validate_request(request, "fit")
    prepare_dir = _prior_output(request, "prepare")
    split = _read_split(prepare_dir / "split_rows.tsv")
    tokens_by_hash, values_by_hash = _read_tokens(prepare_dir / "rank_value_tokens.jsonl")
    labels = _labels_by_hash(request, parameters, roster)

    # Only outer-training donors may influence any fitted state.
    training = [row for row in split if not row["held_out"]]
    if not training:
        raise ScientificAdapterError("fold has no training donors")
    held_units = {row["unit_hash"] for row in split if row["held_out"]}
    if held_units & {row["unit_hash"] for row in training}:
        raise ScientificAdapterError("a donor appears in both training and held-out")

    checkpoint = _verified_checkpoint(request, model_id)
    device = _select_device(str(parameters.get("device", "auto")))
    seed = int(request["run_spec"]["seed"])
    torch.manual_seed(seed)

    embeddings = _embed_rows(
        checkpoint,
        [tokens_by_hash[row["row_hash"]] for row in training],
        [values_by_hash[row["row_hash"]] for row in training],
        device,
        int(parameters.get("embedding_batch_size", 8)),
    )
    targets = torch.tensor(
        [labels[row["row_hash"]] for row in training], dtype=torch.long
    )

    # Donor- and class-balanced weighting, matching the Geneformer common lane.
    from collections import Counter

    donor_counts = Counter(row["unit_hash"] for row in training)
    class_counts = Counter(int(t) for t in targets)
    weights = torch.tensor(
        [
            1.0
            / (donor_counts[row["unit_hash"]] * class_counts[int(target)])
            for row, target in zip(training, targets)
        ],
        dtype=torch.float32,
    )
    weights = weights / weights.sum() * len(weights)

    head = _build_head(head_id, embeddings.shape[1], len(roster), seed).to(device)
    optimizer = torch.optim.Adam(
        head.parameters(), lr=float(parameters.get("head_learning_rate", 1e-3))
    )
    epochs = int(parameters.get("head_epochs", 60))
    features = embeddings.to(device)
    target_device = targets.to(device)
    weight_device = weights.to(device)
    losses = []
    for _ in range(epochs):
        optimizer.zero_grad()
        logits = head(features)
        per_row = torch.nn.functional.cross_entropy(
            logits, target_device, reduction="none"
        )
        loss = (per_row * weight_device).mean()
        loss.backward()
        optimizer.step()
        losses.append(float(loss.detach().cpu()))

    state_path = output / "head_state.pt"
    torch.save(
        {
            "head_id": head_id,
            "state_dict": {k: v.cpu() for k, v in head.state_dict().items()},
            "embedding_dimensions": int(embeddings.shape[1]),
            "roster": list(roster),
            "seed": seed,
        },
        state_path,
    )
    summary_path = output / "fit_summary.json"
    _write_json(
        summary_path,
        {
            "schema_version": "masld-bench-scgpt-fit-v1",
            "model_id": model_id,
            "common_head_id": head_id,
            "embedding_policy": EMBEDDING_POLICY,
            "encoder_frozen": True,
            "device": device,
            "training_rows": len(training),
            "training_donors": len(donor_counts),
            "held_out_donors": len(held_units),
            "epochs": epochs,
            "first_loss": losses[0],
            "final_loss": losses[-1],
        },
    )
    _receipt(
        action="fit",
        request=request,
        output=output,
        artifacts=[state_path, summary_path],
        extra_metadata={
            "model_id": model_id,
            "common_head_id": head_id,
            "encoder_frozen": True,
            "device": device,
            "training_donors": len(donor_counts),
        },
    )


def predict(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import torch

    parameters, roster, model_id, head_id = _validate_request(request, "predict")
    prepare_dir = _prior_output(request, "prepare")
    fit_dir = _prior_output(request, "fit")
    split = _read_split(prepare_dir / "split_rows.tsv")
    tokens_by_hash, values_by_hash = _read_tokens(prepare_dir / "rank_value_tokens.jsonl")

    held = [row for row in split if row["held_out"]]
    if not held:
        raise ScientificAdapterError("fold has no held-out rows to predict")

    state = torch.load(fit_dir / "head_state.pt", map_location="cpu", weights_only=True)
    if state["head_id"] != head_id or list(state["roster"]) != list(roster):
        raise ScientificAdapterError("fitted head does not match this run's contract")

    checkpoint = _verified_checkpoint(request, model_id)
    device = _select_device(str(parameters.get("device", "auto")))
    embeddings = _embed_rows(
        checkpoint,
        [tokens_by_hash[row["row_hash"]] for row in held],
        [values_by_hash[row["row_hash"]] for row in held],
        device,
        int(parameters.get("embedding_batch_size", 8)),
    )
    if embeddings.shape[1] != int(state["embedding_dimensions"]):
        raise ScientificAdapterError("embedding width differs from the fitted head")

    head = _build_head(head_id, embeddings.shape[1], len(roster), int(state["seed"]))
    head.load_state_dict(state["state_dict"])
    head.eval()
    with torch.no_grad():
        probabilities = torch.softmax(head(embeddings), dim=1).cpu()

    prediction_rows = []
    probability_rows = []
    for position, row in enumerate(held):
        values = [float(v) for v in probabilities[position]]
        # Deterministic tie-break identical to the sealed evaluator.
        winner = min(range(len(roster)), key=lambda i: (-values[i], roster[i]))
        prediction_rows.append(
            {
                "row_hash": row["row_hash"],
                "unit_hash": row["unit_hash"],
                "predicted_class": roster[winner],
            }
        )
        probability_rows.append(
            {
                "row_hash": row["row_hash"],
                "unit_hash": row["unit_hash"],
                **{
                    f"probability::{label}": repr(values[index])
                    for index, label in enumerate(roster)
                },
            }
        )
    prediction_rows.sort(key=lambda item: item["row_hash"])
    probability_rows.sort(key=lambda item: item["row_hash"])

    predictions_path = output / "predictions.tsv"
    _write_tsv(
        predictions_path, ("row_hash", "unit_hash", "predicted_class"), prediction_rows
    )
    row_ids_path = output / "row_ids.tsv"
    _write_tsv(
        row_ids_path,
        ("row_hash", "unit_hash"),
        ({"row_hash": r["row_hash"], "unit_hash": r["unit_hash"]} for r in prediction_rows),
    )
    probabilities_path = output / "class_probabilities.tsv"
    _write_tsv(
        probabilities_path,
        ("row_hash", "unit_hash", *(f"probability::{label}" for label in roster)),
        probability_rows,
    )
    _receipt(
        action="predict",
        request=request,
        output=output,
        artifacts=[predictions_path, row_ids_path, probabilities_path],
        extra_metadata={
            "model_id": model_id,
            "common_head_id": head_id,
            "device": device,
            "predicted_rows": len(prediction_rows),
            "class_roster": list(roster),
        },
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--action", required=True, choices=("prepare", "fit", "predict"))
    parser.add_argument("--request", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args(argv)
    try:
        request = json.loads(arguments.request.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificAdapterError("adapter request is invalid JSON") from error
    if not isinstance(request, Mapping):
        raise ScientificAdapterError("adapter request must be an object")
    if arguments.output.exists():
        raise ScientificAdapterError(f"adapter output already exists: {arguments.output}")
    arguments.output.mkdir(parents=True, exist_ok=False)
    {"prepare": prepare, "fit": fit, "predict": predict}[arguments.action](
        arguments.request.resolve(strict=True), request, arguments.output
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
