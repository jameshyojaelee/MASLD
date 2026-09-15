#!/usr/bin/env python
"""UCE (Universal Cell Embedding) common-lane adapter for the frozen cell-state task.

The encoder is frozen.  Only the common head is trained, and it is trained on
outer-training donors alone, so the held-out donors of the fold never
influence any fitted state.  Embeddings follow the registered policy exactly
(``config/artifacts/models/uce/checkpoints.json``): the L2-normalized CLS
decoder output of the pretrained UCE transformer.

UCE keys genes by symbol, not Ensembl ID, and needs a per-gene ESM2 protein
embedding plus a genomic chromosome/start position to build each cell's
"sentence" of sampled gene tokens.  The architecture (``_UCETransformerModel``
below), the special-token indices, and the sampling/tokenization algorithm are
ported from the pinned upstream repository (``snap-stanford/UCE`` revision
``9c416007be15ad6753dc84af4468c1dc10421ab9``, MIT licensed; see
``model.py``/``eval_data.py``/``evaluate.py``/``eval_single_anndata.py`` whose
SHA-256 hashes are recorded in ``config/artifacts/models/uce/checkpoints.json
source_records``).  Two upstream behaviors are deliberately reproduced as-is
rather than "fixed": (1) the transformer's padding mask is passed to
``nn.TransformerEncoder`` as a plain 0/1 float tensor (an additive bias) as
upstream's own ``evaluate.py`` does, not converted to a boolean mask; (2) the
sampled-gene sort order uses NumPy's default unstable-quicksort tie-breaking
semantics, made ``stable`` here only for adapter-to-adapter reproducibility,
which does not change upstream's own tie-breaking requirements since ties are
broken on (chromosome id, genomic start) which are the same sort keys either
way.

Gene sampling and chromosome-order shuffling are the one place the official
UCE evaluator mutates NumPy's *global* random state (see the registry's
``sampling_warning``).  This adapter never touches global NumPy state: each
cell's sampling draws come from a local ``numpy.random.Generator`` seeded
deterministically from ``(join_namespace, row_hash, run seed)``, so the same
seed reproduces byte-identical embeddings regardless of batch composition or
processing order.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import io
import json
import math
import pickle
from pathlib import Path
from typing import Any, Mapping, Sequence

from hvg_pca_logistic import (  # type: ignore[import-not-found]
    ScientificAdapterError,
    _artifact_record,
    _input_by_role,
    _require_sha256,
    _validate_file_artifact,
    _write_json,
    _write_tsv,
    fold_index,
    join_hash,
)


class _NoGlobalsUnpickler(pickle.Unpickler):
    """Deny every pickle global/class reference outright.

    ``species_offsets.pkl`` is expected to unpickle to nothing but a small
    plain ``str -> int`` dict, so no pickle opcode should ever need to
    construct a class or call a function.  This mirrors the restricted-
    unpickling discipline in ``checkpoint_preflight.py``
    (``load_restricted_plain_mapping``, which additionally allow-lists two
    numpy scalar reconstructors for Geneformer's dictionaries); it is
    duplicated in this narrower, stricter form here rather than imported,
    because this adapter is invoked as a standalone script from the
    ``adapters/`` directory and ``checkpoint_preflight.py`` lives one
    package level up, outside that script's default import path.
    """

    def find_class(self, module: str, name: str) -> Any:
        raise ScientificAdapterError(
            f"pickle global/class construction is forbidden: {module}.{name}"
        )


def load_restricted_plain_mapping(payload: bytes) -> dict[str, Any]:
    stream = io.BytesIO(payload)
    try:
        value = _NoGlobalsUnpickler(stream).load()
    except ScientificAdapterError:
        raise
    except Exception as error:
        raise ScientificAdapterError(
            f"auxiliary pickle cannot be decoded as restricted plain data: {error}"
        ) from error
    if stream.read(1):
        raise ScientificAdapterError(
            "auxiliary pickle contains trailing content after the first object"
        )
    if not isinstance(value, dict) or not value:
        raise ScientificAdapterError(
            "auxiliary pickle must contain one nonempty plain dictionary"
        )
    return value

TASK_ID = "cell_state_mapping"
DATASET_ID = "resource_atlas_current"
RUNTIME_ID = "gpu_uce"
RECEIPT_SCHEMA = "masld-bench-adapter-receipt-v1"
ADAPTER_ID = "uce_common_lane_v1"
MODEL_IDS = ("uce_4l", "uce_33l")
HEAD_IDS = ("linear", "two_layer_mlp")
EMBEDDING_POLICY = "l2_normalized_cls_decoder_output"

# Frozen architecture requirements, matching
# config/artifacts/models/uce/checkpoints.json:architecture_contract exactly.
N_LAYERS = {"uce_4l": 4, "uce_33l": 33}
CHECKPOINT_FILENAMES = {
    "uce_4l": "4layer_model.torch",
    "uce_33l": "33l_8ep_1024t_1280.torch",
}
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
# Special-token indices below CHROM_TOKEN_OFFSET, from the pinned upstream
# eval_single_anndata.py argparse defaults (--pad_token_idx, --chrom_token_right_idx).
PAD_TOKEN_INDEX = 0
CHROM_TOKEN_RIGHT_INDEX = 2
# Row count of the pe_embedding table shipped inside every UCE checkpoint
# (143,574 real gene/ESM2 rows + 1,895 chromosome-identity rows).
PE_EMBEDDING_ROWS = 145469
HUMAN_SPECIES_KEY = "human"
# The four hyperparameters the RunSpec must restate verbatim; every other
# architecture constant above is intrinsic to model_id, not caller-supplied.
_REQUIRED_ARCH_HYPERPARAMETERS = (
    "pad_length",
    "sample_size",
    "cls_token_index",
    "chromosome_token_offset",
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


def _full_block(in_features: int, out_features: int, dropout: float):
    from torch import nn

    return nn.Sequential(
        nn.Linear(in_features, out_features, bias=True),
        nn.LayerNorm(out_features),
        nn.GELU(),
        nn.Dropout(p=dropout),
    )


def _make_positional_encoding(d_model: int, dropout: float, max_len: int):
    """Port of upstream model.py PositionalEncoding, as a module instance."""

    import torch
    from torch import nn

    class _PositionalEncoding(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.dropout = nn.Dropout(p=dropout)
            position = torch.arange(max_len).unsqueeze(1)
            div_term = torch.exp(
                torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model)
            )
            pe = torch.zeros(max_len, 1, d_model)
            pe[:, 0, 0::2] = torch.sin(position * div_term)
            pe[:, 0, 1::2] = torch.cos(position * div_term)
            self.register_buffer("pe", pe)

        def forward(self, x):
            return self.dropout(x + self.pe[: x.size(0)])

    return _PositionalEncoding()


def _build_model(nlayers: int):
    """Port of upstream model.py TransformerModel, extended with the fixed
    pe_embedding lookup + normalization that upstream's evaluate.py performs
    just before calling the model, so that ``forward`` takes raw token ids.
    Every submodule the checkpoint's state_dict defines is reconstructed here
    (including ``binary_decoder``/``gene_embedding_layer``, unused by the CLS
    embedding path) so the checkpoint loads under ``strict=True`` with no
    silently-ignored or silently-missing keys.
    """

    import torch
    from torch import nn

    d_model = ARCH["d_model"]
    dropout = ARCH["dropout"]

    class _UCETransformerModel(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.d_model = d_model
            self.pos_encoder = _make_positional_encoding(
                d_model, dropout, ARCH["pad_length"]
            )
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
            """tokens: LongTensor[seq_len, batch]; mask: FloatTensor[batch, seq_len]
            with 1 for valid (non-padding) positions, 0 for padding — exactly
            upstream evaluate.py's requirements, including passing the float
            ``(1 - mask)`` directly as ``src_key_padding_mask`` (an additive
            attention bias upstream never converted to a boolean mask).
            """

            src = self.pe_embedding(tokens)
            src = nn.functional.normalize(src, dim=2)
            src = self.encoder(src) * math.sqrt(self.d_model)
            src = self.pos_encoder(src)
            output = self.transformer_encoder(src, src_key_padding_mask=(1 - mask))
            gene_output = self.decoder(output)
            embedding = gene_output[0, :, :]
            return nn.functional.normalize(embedding, dim=1)

    return _UCETransformerModel()


def _load_checkpoint_model(checkpoint_path: Path, all_tokens_path: Path, model_id: str, device: str):
    """Restricted-load a UCE checkpoint and its shared token table.

    Both files are plain ``torch.save``d tensors/state-dicts (verified above
    to contain no custom classes), so ``weights_only=True`` — PyTorch's own
    restricted unpickler that permits only tensor/collection reconstruction —
    is sufficient; no allowlist widening is needed for either file.
    """

    import torch

    model = _build_model(N_LAYERS[model_id])
    state_dict = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    if not isinstance(state_dict, Mapping):
        raise ScientificAdapterError("UCE checkpoint is not a plain state dict")
    model.load_state_dict(state_dict, strict=True)

    all_tokens = torch.load(all_tokens_path, map_location="cpu", weights_only=True)
    if not hasattr(all_tokens, "shape") or tuple(all_tokens.shape) != (
        PE_EMBEDDING_ROWS,
        ARCH["token_dim"],
    ):
        raise ScientificAdapterError("UCE all_tokens.torch shape differs from the frozen contract")
    if not torch.equal(all_tokens, state_dict["pe_embedding.weight"]):
        raise ScientificAdapterError(
            "UCE all_tokens.torch does not match the checkpoint's own token embedding table"
        )

    model.eval()
    model.to(device)
    return model


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
        raise ScientificAdapterError("RunSpec names an unsupported UCE model")
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
    for field in _REQUIRED_ARCH_HYPERPARAMETERS:
        if parameters.get(field) != ARCH[field]:
            raise ScientificAdapterError(
                f"UCE architecture hyperparameter {field!r} differs from the frozen contract"
            )
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


def _checkpoint_path(request: Mapping[str, Any], model_id: str) -> Path:
    role = f"model_checkpoint_weights:{model_id}"
    path = _validate_file_artifact(_input_by_role(request, role), "UCE checkpoint weights")
    if path.name != CHECKPOINT_FILENAMES[model_id]:
        raise ScientificAdapterError(f"UCE checkpoint file name differs for {model_id}")
    return path


def _all_tokens_path(request: Mapping[str, Any]) -> Path:
    role = "model_checkpoint_all_tokens:uce_shared"
    return _validate_file_artifact(_input_by_role(request, role), "UCE token embedding table")


def _species_offsets_path(request: Mapping[str, Any]) -> Path:
    role = "model_checkpoint_species_offsets:uce_shared"
    return _validate_file_artifact(_input_by_role(request, role), "UCE species offsets")


def _species_chrom_path(request: Mapping[str, Any]) -> Path:
    role = "model_checkpoint_species_chrom:uce_shared"
    return _validate_file_artifact(_input_by_role(request, role), "UCE species chromosome table")


def _human_protein_embeddings_path(request: Mapping[str, Any]) -> Path:
    role = "model_checkpoint_protein_embeddings:uce_shared_human"
    return _validate_file_artifact(
        _input_by_role(request, role), "UCE human protein embeddings"
    )


def _atlas_path(request: Mapping[str, Any]) -> Path:
    role = "dataset_view_data:resource_atlas_geneformer_smoke_1000_v1"
    return _validate_file_artifact(_input_by_role(request, role), "Atlas smoke H5AD")


def _load_atlas(path: Path) -> Any:
    import anndata

    adata = anndata.read_h5ad(path)
    if not {"donor_id", "broad_label"}.issubset(adata.obs.columns):
        raise ScientificAdapterError("Atlas H5AD lacks required observation fields")
    if "source_feature_id" not in adata.var.columns:
        raise ScientificAdapterError("Atlas H5AD lacks the gene-symbol feature axis")
    if len(set(map(str, adata.obs_names))) != adata.n_obs:
        raise ScientificAdapterError("Atlas row identifiers are not unique")
    return adata


def _counts(adata: Any):
    import numpy as np

    matrix = adata.X
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


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


def _gene_universe(adata: Any, request: Mapping[str, Any]):
    """Build the per-gene (token row, chromosome code, genomic start) arrays
    UCE's own tokenizer needs, restricted to genes with both a human ESM2
    protein embedding and a chromosome/start record.  The chromosome
    categorical coding is fit over the *entire* multi-species table before
    subsetting to human, exactly matching upstream data_utils.get_spec_chrom_csv
    + adata_path_to_prot_chrom_starts, because the chromosome-identity tokens
    the pretrained model was built with are indexed by that shared coding.
    """

    import numpy as np
    import pandas as pd
    import torch

    offsets = load_restricted_plain_mapping(_species_offsets_path(request).read_bytes())
    human_offset = offsets.get(HUMAN_SPECIES_KEY)
    if not isinstance(human_offset, int):
        raise ScientificAdapterError("species offsets table lacks an integer human offset")

    table = pd.read_csv(_species_chrom_path(request))
    if not {"gene_symbol", "chromosome", "start", "species"}.issubset(table.columns):
        raise ScientificAdapterError("species chromosome table lacks required columns")
    categories = pd.Categorical(table["species"] + "_" + table["chromosome"])
    human_mask = table["species"] == HUMAN_SPECIES_KEY
    human_symbols = table.loc[human_mask, "gene_symbol"].str.upper()
    if human_symbols.duplicated().any():
        raise ScientificAdapterError("human chromosome table has duplicate gene symbols")
    human_codes = pd.Series(categories.codes[human_mask.to_numpy()], index=human_symbols)
    human_starts = pd.Series(table.loc[human_mask, "start"].to_numpy(), index=human_symbols)
    chrom_by_symbol = dict(zip(human_symbols, human_codes))
    start_by_symbol = dict(zip(human_symbols, human_starts))

    protein_embeddings = torch.load(
        _human_protein_embeddings_path(request), map_location="cpu", weights_only=True
    )
    if not isinstance(protein_embeddings, Mapping) or not protein_embeddings:
        raise ScientificAdapterError("human protein embedding table is empty or malformed")
    ordered_symbols = [str(symbol).upper() for symbol in protein_embeddings.keys()]
    if len(set(ordered_symbols)) != len(ordered_symbols):
        raise ScientificAdapterError("human protein embedding keys are not unique under uppercasing")
    local_index_by_symbol = {symbol: index for index, symbol in enumerate(ordered_symbols)}

    gene_symbols = [str(g).upper() for g in adata.var["source_feature_id"]]
    n_genes = len(gene_symbols)
    token_row = np.full(n_genes, -1, dtype=np.int64)
    chrom_code = np.full(n_genes, -1, dtype=np.int64)
    start = np.zeros(n_genes, dtype=np.int64)
    usable = np.zeros(n_genes, dtype=bool)
    for index, symbol in enumerate(gene_symbols):
        if symbol not in local_index_by_symbol or symbol not in chrom_by_symbol:
            continue
        token_row[index] = local_index_by_symbol[symbol] + human_offset
        chrom_code[index] = int(chrom_by_symbol[symbol])
        start[index] = int(start_by_symbol[symbol])
        usable[index] = True

    usable_count = int(usable.sum())
    summary = {
        "schema_version": "masld-bench-uce-gene-universe-v1",
        "gene_identifier_source": "obs.var.source_feature_id",
        "genes_total": n_genes,
        "genes_usable": usable_count,
        "genes_dropped": n_genes - usable_count,
        "overlap_fraction": (usable_count / n_genes) if n_genes else 0.0,
        "human_species_offset": int(human_offset),
        "chromosome_categories_total": int(len(categories.categories)),
    }
    if usable_count == 0:
        raise ScientificAdapterError(
            "no dataset gene maps to a UCE human protein embedding and chromosome "
            f"position; overlap_fraction=0/{n_genes}"
        )
    return token_row, chrom_code, start, usable, summary


def _derive_cell_seed(namespace: str, row_hash: str, run_seed: int) -> int:
    """A per-cell seed for a LOCAL numpy.random.Generator, never global state.

    This is the adapter's answer to the registry's sampling_warning: the
    official evaluator advances NumPy's shared global RNG per cell, making its
    embeddings depend on batch order.  Seeding an independent Generator per
    cell from (namespace, row_hash, run seed) makes this adapter's embeddings
    depend only on the cell's own identity and the run's seed.
    """

    material = f"{namespace}\0{ADAPTER_ID}\0{run_seed}\0{row_hash}".encode("utf-8")
    return int.from_bytes(sha256(material).digest()[:8], "big")


def _tokenize_cell(counts, token_row, chrom_code, start, rng) -> tuple[list[int], int]:
    """Port of upstream eval_data.py sample_cell_sentences, for one cell."""

    import numpy as np

    weights = np.log1p(np.asarray(counts, dtype=np.float64))
    weight_total = float(weights.sum())
    if weight_total <= 0.0:
        raise ScientificAdapterError("a cell has no usable counts over the UCE gene universe")
    weights = weights / weight_total

    choice = rng.choice(weights.shape[0], size=ARCH["sample_size"], replace=True, p=weights)
    chosen_chrom = chrom_code[choice]
    order = np.argsort(chosen_chrom, kind="stable")
    choice = choice[order]
    chosen_chrom = chosen_chrom[order]
    chosen_start = start[choice]

    pad_length = ARCH["pad_length"]
    ordered = np.full(pad_length, ARCH["cls_token_index"], dtype=np.int64)
    position = 1
    unique_chroms = np.unique(chosen_chrom)
    rng.shuffle(unique_chroms)
    for chrom in unique_chroms:
        ordered[position] = int(chrom) + ARCH["chromosome_token_offset"]
        position += 1
        members = np.flatnonzero(chosen_chrom == chrom)
        start_order = np.argsort(chosen_start[members], kind="stable")
        block = token_row[choice[members[start_order]]]
        ordered[position : position + block.shape[0]] = block
        position += int(block.shape[0])
        ordered[position] = CHROM_TOKEN_RIGHT_INDEX
        position += 1
    if position > pad_length:
        raise ScientificAdapterError("UCE cell sentence exceeded the frozen pad length")
    ordered[position:] = PAD_TOKEN_INDEX
    return [int(v) for v in ordered], position


def _embed_rows(model, device: str, tokens_by_row, row_order: Sequence[str], batch_size: int):
    import torch

    pad_length = ARCH["pad_length"]
    outputs = []
    with torch.no_grad():
        for start_index in range(0, len(row_order), batch_size):
            chunk = row_order[start_index : start_index + batch_size]
            token_batch = torch.zeros(len(chunk), pad_length, dtype=torch.long)
            mask_batch = torch.zeros(len(chunk), pad_length, dtype=torch.float32)
            for index, row_hash in enumerate(chunk):
                tokens, content_length = tokens_by_row[row_hash]
                token_batch[index] = torch.tensor(tokens, dtype=torch.long)
                mask_batch[index, :content_length] = 1.0
            token_batch = token_batch.permute(1, 0).to(device)  # [seq_len, batch]
            mask_batch = mask_batch.to(device)
            embedding = model(token_batch, mask_batch)
            outputs.append(embedding.detach().cpu())
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


def _read_tokens(path: Path) -> dict[str, tuple[list[int], int]]:
    result: dict[str, tuple[list[int], int]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            result[str(record["row_hash"])] = (
                [int(v) for v in record["tokens"]],
                int(record["content_length"]),
            )
    return result


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


def _prior_output(request: Mapping[str, Any], action: str) -> Path:
    outputs = request.get("prior_action_outputs")
    if not isinstance(outputs, Mapping) or action not in outputs:
        raise ScientificAdapterError(f"request lacks the frozen {action} output")
    return Path(str(outputs[action])).resolve(strict=True)


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

    token_row, chrom_code, start, usable, gene_summary = _gene_universe(adata, request)
    usable_columns = usable.nonzero()[0]
    counts = _counts(adata)[:, usable_columns]
    usable_token_row = token_row[usable_columns]
    usable_chrom_code = chrom_code[usable_columns]
    usable_start = start[usable_columns]

    namespace = str(parameters["join_namespace"])
    run_seed = int(request["run_spec"]["seed"])
    cell_index_by_row_hash = {
        join_hash(namespace, "row", str(name)): index
        for index, name in enumerate(adata.obs_names)
    }

    import numpy as np

    tokens_path = output / "cell_sentence_tokens.jsonl"
    truncated = 0
    with tokens_path.open("w", encoding="utf-8") as handle:
        for row in rows:
            row_hash = row["row_hash"]
            cell_index = cell_index_by_row_hash[row_hash]
            rng = np.random.Generator(
                np.random.PCG64(_derive_cell_seed(namespace, row_hash, run_seed))
            )
            tokens, content_length = _tokenize_cell(
                counts[cell_index], usable_token_row, usable_chrom_code, usable_start, rng
            )
            if content_length >= ARCH["pad_length"]:
                truncated += 1
            handle.write(
                json.dumps(
                    {
                        "row_hash": row_hash,
                        "tokens": tokens,
                        "content_length": content_length,
                    },
                    sort_keys=True,
                )
                + "\n"
            )

    summary = {**gene_summary, "cells": len(rows), "cells_at_pad_length": truncated}
    summary_path = output / "tokenization_summary.json"
    _write_json(summary_path, summary)
    _receipt(
        action="prepare",
        request=request,
        output=output,
        artifacts=[split_path, tokens_path, summary_path],
        extra_metadata={
            "model_id": model_id,
            "genes_usable": gene_summary["genes_usable"],
            "genes_dropped": gene_summary["genes_dropped"],
            "gene_overlap_fraction": gene_summary["overlap_fraction"],
        },
    )


def fit(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import torch

    parameters, roster, model_id, head_id = _validate_request(request, "fit")
    prepare_dir = _prior_output(request, "prepare")
    split = _read_split(prepare_dir / "split_rows.tsv")
    tokens_by_row = _read_tokens(prepare_dir / "cell_sentence_tokens.jsonl")
    labels = _labels_by_hash(request, parameters, roster)

    training = [row for row in split if not row["held_out"]]
    if not training:
        raise ScientificAdapterError("fold has no training donors")
    held_units = {row["unit_hash"] for row in split if row["held_out"]}
    if held_units & {row["unit_hash"] for row in training}:
        raise ScientificAdapterError("a donor appears in both training and held-out")

    device = _select_device(str(parameters.get("device", "auto")))
    seed = int(request["run_spec"]["seed"])
    torch.manual_seed(seed)

    model = _load_checkpoint_model(
        _checkpoint_path(request, model_id), _all_tokens_path(request), model_id, device
    )
    training_row_hashes = [row["row_hash"] for row in training]
    embeddings = _embed_rows(
        model, device, tokens_by_row, training_row_hashes,
        int(parameters.get("embedding_batch_size", 4)),
    )
    del model
    if device == "cuda":
        torch.cuda.empty_cache()
    targets = torch.tensor(
        [labels[row_hash] for row_hash in training_row_hashes], dtype=torch.long
    )

    from collections import Counter

    donor_by_hash = {row["row_hash"]: row["unit_hash"] for row in training}
    donor_counts = Counter(donor_by_hash[h] for h in training_row_hashes)
    class_counts = Counter(int(t) for t in targets)
    weights = torch.tensor(
        [
            1.0 / (donor_counts[donor_by_hash[row_hash]] * class_counts[int(target)])
            for row_hash, target in zip(training_row_hashes, targets)
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
            "schema_version": "masld-bench-uce-fit-v1",
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
    tokens_by_row = _read_tokens(prepare_dir / "cell_sentence_tokens.jsonl")

    held = [row for row in split if row["held_out"]]
    if not held:
        raise ScientificAdapterError("fold has no held-out rows to predict")

    state = torch.load(fit_dir / "head_state.pt", map_location="cpu", weights_only=True)
    if state["head_id"] != head_id or list(state["roster"]) != list(roster):
        raise ScientificAdapterError("fitted head does not match this run's contract")

    device = _select_device(str(parameters.get("device", "auto")))
    model = _load_checkpoint_model(
        _checkpoint_path(request, model_id), _all_tokens_path(request), model_id, device
    )
    held_row_hashes = [row["row_hash"] for row in held]
    embeddings = _embed_rows(
        model, device, tokens_by_row, held_row_hashes,
        int(parameters.get("embedding_batch_size", 4)),
    )
    del model
    if device == "cuda":
        torch.cuda.empty_cache()
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
