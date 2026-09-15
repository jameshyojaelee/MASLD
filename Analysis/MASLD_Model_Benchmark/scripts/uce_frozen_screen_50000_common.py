#!/usr/bin/env python3
"""Extraction-only UCE requirements shared by the 50k fixture and encoder runners."""

from __future__ import annotations

from hashlib import sha256
import io
import json
import math
from pathlib import Path
import pickle
import re
from typing import Any, Mapping


MODEL_IDS = ("uce_4l", "uce_33l")
N_LAYERS = {"uce_4l": 4, "uce_33l": 33}
ARCH = {
    "chromosome_token_offset": 143574,
    "cls_token_index": 3,
    "d_hid": 5120,
    "d_model": 1280,
    "dropout": 0.05,
    "n_heads": 20,
    "output_dim": 1280,
    "pad_length": 1536,
    "sample_size": 1024,
    "token_dim": 5120,
}
PAD_TOKEN_INDEX = 0
CHROMOSOME_RIGHT_TOKEN_INDEX = 2
PE_EMBEDDING_ROWS = 145469
TOKENIZATION_POLICY = "uce_common_order_invariant_v1"
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


class UCEActivationError(RuntimeError):
    """Raised when a read-only UCE activation requirement differs."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def require_sha256(value: Any, label: str) -> str:
    text = str(value)
    if _SHA256.fullmatch(text) is None:
        raise UCEActivationError(f"{label} is not a lowercase SHA-256")
    return text


def load_activation_contract(path: Path, project_root: Path) -> dict[str, Any]:
    contract = json.loads(path.read_text(encoding="utf-8"))
    if (
        contract.get("schema_version")
        != "masld-bench-uce-frozen-screen-activation-v1"
        or contract.get("dataset_view_id")
        != "resource_atlas_frozen_screen_50000_v1"
        or contract.get("split_id") != "resource_atlas_study_outer_5fold_v1"
        or contract.get("rows") != 50_000
        or contract.get("donors") != 102
        or contract.get("studies") != 7
        or contract.get("repository_revision")
        != "9c416007be15ad6753dc84af4468c1dc10421ab9"
        or contract.get("licenses")
        != {
            "code": "MIT",
            "code_evidence": "https://raw.githubusercontent.com/snap-stanford/UCE/9c416007be15ad6753dc84af4468c1dc10421ab9/LICENSE",
            "weights": "CC-BY-4.0",
            "weights_evidence": "https://api.figshare.com/v2/articles/24320806",
        }
    ):
        raise UCEActivationError("UCE activation identity or license contract differs")
    if set(contract.get("checkpoints", {})) != set(MODEL_IDS):
        raise UCEActivationError("UCE checkpoint roster differs")
    if (
        contract["checkpoints"]["uce_4l"].get("layers") != 4
        or contract["checkpoints"]["uce_33l"].get("layers") != 33
        or contract["checkpoints"]["uce_4l"].get("sha256")
        == contract["checkpoints"]["uce_33l"].get("sha256")
    ):
        raise UCEActivationError("UCE checkpoints are not distinct 4L and 33L identities")
    for model_id, record in contract["checkpoints"].items():
        if (
            record.get("exposure_status") != "target_label_unexposed"
            or record.get("sealed_champion_eligible") is not True
            or record.get("layers") != N_LAYERS[model_id]
        ):
            raise UCEActivationError(f"{model_id} exposure or layer contract differs")
        require_sha256(record.get("sha256"), f"{model_id} checkpoint")
    tokenization = contract.get("tokenization", {})
    required_tokenization = {
        "gene_min_cells": 10,
        "cell_min_genes_before_protein_mapping": 25,
        "sample_size": ARCH["sample_size"],
        "pad_length": ARCH["pad_length"],
        "cls_token_index": ARCH["cls_token_index"],
        "pad_token_index": PAD_TOKEN_INDEX,
        "chromosome_right_token_index": CHROMOSOME_RIGHT_TOKEN_INDEX,
        "chromosome_token_offset": ARCH["chromosome_token_offset"],
        "common_lane_seed": 20260824,
        "truncation_policy": "fail_if_sentence_exceeds_1536;never_truncate",
        "output_policy": "l2_normalized_cls_decoder_output",
    }
    if any(tokenization.get(key) != value for key, value in required_tokenization.items()):
        raise UCEActivationError("UCE tokenization contract differs")
    resolved = dict(contract)
    resolved["_project_root"] = project_root.resolve(strict=True)
    return resolved


def contract_path(contract: Mapping[str, Any], record: Mapping[str, Any]) -> Path:
    root = Path(contract["_project_root"])
    relative = Path(str(record.get("path", "")))
    if relative.is_absolute() or relative == Path("."):
        raise UCEActivationError("contract artifact path is not a relative path")
    try:
        path = (root / relative).resolve(strict=True)
        path.relative_to(root)
    except (OSError, ValueError) as error:
        raise UCEActivationError("contract artifact escapes the project root") from error
    if not path.is_file() and not path.is_dir():
        raise UCEActivationError("contract artifact is not a regular file or directory")
    return path


def verify_record(contract: Mapping[str, Any], record: Mapping[str, Any], label: str) -> Path:
    path = contract_path(contract, record)
    expected = require_sha256(record.get("sha256"), label)
    if not path.is_file() or sha256_file(path) != expected:
        raise UCEActivationError(f"{label} SHA-256 differs")
    size = record.get("size_bytes")
    if size is not None and path.stat().st_size != int(size):
        raise UCEActivationError(f"{label} size differs")
    return path


class _NoGlobalsUnpickler(pickle.Unpickler):
    def find_class(self, module: str, name: str) -> Any:
        raise UCEActivationError(
            f"pickle global/class construction is forbidden: {module}.{name}"
        )


def load_restricted_plain_mapping(path: Path) -> dict[str, Any]:
    stream = io.BytesIO(path.read_bytes())
    try:
        value = _NoGlobalsUnpickler(stream).load()
    except UCEActivationError:
        raise
    except Exception as error:
        raise UCEActivationError("restricted mapping decode failed") from error
    if stream.read(1) or not isinstance(value, dict) or not value:
        raise UCEActivationError("restricted mapping is not one plain nonempty dictionary")
    return value


def derive_cell_seed(row_id: str, run_seed: int) -> int:
    material = f"{TOKENIZATION_POLICY}\0{run_seed}\0{row_id}".encode("utf-8")
    return int.from_bytes(sha256(material).digest()[:8], "big")


def tokenize_cell(
    counts: Any,
    token_rows: Any,
    chromosome_codes: Any,
    genomic_starts: Any,
    *,
    row_id: str,
    run_seed: int,
) -> tuple[Any, int, int]:
    """Sample one official-format sentence without ever truncating it."""

    import numpy as np

    values = np.asarray(counts, dtype=np.float64)
    if values.ndim != 1 or values.size < 1 or np.any(values <= 0.0):
        raise UCEActivationError("UCE tokenization requires positive expressed-gene counts")
    weights = np.log1p(values)
    total = float(weights.sum())
    if not np.isfinite(total) or total <= 0.0:
        raise UCEActivationError("UCE sampling weights are invalid")
    rng = np.random.default_rng(derive_cell_seed(row_id, run_seed))
    choice = rng.choice(values.size, size=ARCH["sample_size"], replace=True, p=weights / total)
    chosen_chromosomes = chromosome_codes[choice]
    order = np.argsort(chosen_chromosomes, kind="stable")
    choice = choice[order]
    chosen_chromosomes = chosen_chromosomes[order]
    chosen_starts = genomic_starts[choice]
    sentence = np.full(ARCH["pad_length"], PAD_TOKEN_INDEX, dtype=np.int32)
    sentence[0] = ARCH["cls_token_index"]
    position = 1
    unique_chromosomes = np.unique(chosen_chromosomes)
    rng.shuffle(unique_chromosomes)
    for chromosome in unique_chromosomes:
        members = np.flatnonzero(chosen_chromosomes == chromosome)
        start_order = np.argsort(chosen_starts[members], kind="stable")
        block = token_rows[choice[members[start_order]]]
        required = 1 + block.size + 1
        if position + required > ARCH["pad_length"]:
            raise UCEActivationError("UCE sentence exceeded 1536 tokens; truncation is forbidden")
        sentence[position] = int(chromosome) + ARCH["chromosome_token_offset"]
        position += 1
        sentence[position : position + block.size] = block
        position += int(block.size)
        sentence[position] = CHROMOSOME_RIGHT_TOKEN_INDEX
        position += 1
    return sentence, position, int(np.unique(choice).size)


def _full_block(in_features: int, out_features: int, dropout: float):
    from torch import nn

    return nn.Sequential(
        nn.Linear(in_features, out_features, bias=True),
        nn.LayerNorm(out_features),
        nn.GELU(),
        nn.Dropout(p=dropout),
    )


def build_model(nlayers: int):
    """Independent port of the pinned MIT UCE encoder, with no task head."""

    import torch
    from torch import nn

    if nlayers not in set(N_LAYERS.values()):
        raise UCEActivationError("UCE layer count differs")
    d_model = ARCH["d_model"]
    dropout = ARCH["dropout"]

    class PositionalEncoding(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.dropout = nn.Dropout(p=dropout)
            positions = torch.arange(ARCH["pad_length"]).unsqueeze(1)
            divisor = torch.exp(
                torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model)
            )
            values = torch.zeros(ARCH["pad_length"], 1, d_model)
            values[:, 0, 0::2] = torch.sin(positions * divisor)
            values[:, 0, 1::2] = torch.cos(positions * divisor)
            self.register_buffer("pe", values)

        def forward(self, values):
            return self.dropout(values + self.pe[: values.size(0)])

    class UCETransformer(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.d_model = d_model
            self.pos_encoder = PositionalEncoding()
            self.encoder = nn.Sequential(
                nn.Linear(ARCH["token_dim"], d_model), nn.GELU(), nn.LayerNorm(d_model)
            )
            layer = nn.TransformerEncoderLayer(
                d_model, ARCH["n_heads"], ARCH["d_hid"], dropout
            )
            self.transformer_encoder = nn.TransformerEncoder(layer, nlayers)
            self.decoder = nn.Sequential(
                _full_block(d_model, 1024, dropout),
                _full_block(1024, ARCH["output_dim"], dropout),
                _full_block(ARCH["output_dim"], ARCH["output_dim"], dropout),
                nn.Linear(ARCH["output_dim"], ARCH["output_dim"]),
            )
            self.binary_decoder = nn.Sequential(
                _full_block(ARCH["output_dim"] + 1280, 2048, dropout),
                _full_block(2048, 512, dropout),
                _full_block(512, 128, dropout),
                nn.Linear(128, 1),
            )
            self.gene_embedding_layer = nn.Sequential(
                nn.Linear(ARCH["token_dim"], d_model), nn.GELU(), nn.LayerNorm(d_model)
            )
            self.pe_embedding = nn.Embedding(PE_EMBEDDING_ROWS, ARCH["token_dim"])

        def forward(self, tokens, mask):
            source = self.pe_embedding(tokens)
            source = nn.functional.normalize(source, dim=2)
            source = self.encoder(source) * math.sqrt(self.d_model)
            source = self.pos_encoder(source)
            # The pinned evaluator passes this 0/1-derived float mask directly.
            encoded = self.transformer_encoder(source, src_key_padding_mask=(1 - mask))
            decoded = self.decoder(encoded)
            return nn.functional.normalize(decoded[0, :, :], dim=1)

    return UCETransformer()


def inspect_state_dict(state_dict: Mapping[str, Any], model_id: str) -> dict[str, Any]:
    import torch

    if model_id not in MODEL_IDS or not isinstance(state_dict, Mapping):
        raise UCEActivationError("UCE state-dict identity differs")
    if not state_dict or any(not isinstance(value, torch.Tensor) for value in state_dict.values()):
        raise UCEActivationError("UCE checkpoint is not a tensor-only state dictionary")
    layers = {
        int(match.group(1))
        for key in state_dict
        if (match := re.match(r"transformer_encoder\.layers\.(\d+)\.", key))
    }
    if layers != set(range(N_LAYERS[model_id])):
        raise UCEActivationError("UCE checkpoint transformer layer roster differs")
    token_table = state_dict.get("pe_embedding.weight")
    if token_table is None or tuple(token_table.shape) != (
        PE_EMBEDDING_ROWS,
        ARCH["token_dim"],
    ):
        raise UCEActivationError("UCE checkpoint token table shape differs")
    return {
        "tensor_count": len(state_dict),
        "total_numel": sum(int(value.numel()) for value in state_dict.values()),
        "layers": len(layers),
        "token_table_shape": list(token_table.shape),
    }


def load_checkpoint_model(
    checkpoint: Path,
    all_tokens_path: Path,
    *,
    model_id: str,
    device: str,
) -> tuple[Any, dict[str, Any]]:
    import torch

    state_dict = torch.load(
        checkpoint, map_location="cpu", weights_only=True, mmap=True
    )
    inspection = inspect_state_dict(state_dict, model_id)
    token_table = torch.load(
        all_tokens_path, map_location="cpu", weights_only=True, mmap=True
    )
    if not isinstance(token_table, torch.Tensor) or tuple(token_table.shape) != (
        PE_EMBEDDING_ROWS,
        ARCH["token_dim"],
    ):
        raise UCEActivationError("UCE all-token table shape differs")
    checkpoint_tokens = state_dict["pe_embedding.weight"]
    for start in range(0, PE_EMBEDDING_ROWS, 4096):
        stop = min(start + 4096, PE_EMBEDDING_ROWS)
        if not torch.equal(token_table[start:stop], checkpoint_tokens[start:stop]):
            raise UCEActivationError("UCE shared token table differs from checkpoint weights")
    del token_table
    model = build_model(N_LAYERS[model_id])
    incompatibility = model.load_state_dict(state_dict, strict=True, assign=True)
    if incompatibility.missing_keys or incompatibility.unexpected_keys:
        raise UCEActivationError("UCE strict checkpoint restore differs")
    model.eval().to(device)
    inspection["strict_weights_only_restore"] = True
    inspection["shared_token_table_equal"] = True
    return model, inspection
