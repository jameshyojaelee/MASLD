#!/usr/bin/env python3
"""Extract outcome-blind pooled embeddings from one frozen DNA-language model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
import sys
from typing import Callable

import numpy as np
import torch
from safetensors.torch import load_file
from transformers import EsmConfig, EsmForMaskedLM, EsmTokenizer, PreTrainedTokenizerFast

from scripts.gse281364_dna_lm_native_contract import (
    ALLELES,
    MODELS,
    NativeContractError,
    dnabert_encoding,
    load_config,
    nt_phase_spans,
    nt_pool_indices,
    overlapping_indices,
    read_fixture,
)


TOKEN_IDS = {"A": 7, "C": 8, "G": 9, "T": 10, "N": 11}


def require_cuda() -> torch.device:
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda" or torch.cuda.device_count() != 1:
        raise NativeContractError("exactly one CUDA device is required")
    return device


def load_dnabert(
    root: Path, record: dict[str, object], device: torch.device
) -> tuple[Callable[[str, int, int], np.ndarray], int, dict[str, object]]:
    bundle = root / record["runtime_bundle_path"] / "bundle"
    sys.path.insert(0, str(bundle.parent))
    from bundle import bert_layers
    from bundle.configuration_bert import BertConfig

    bert_layers.flash_attn_qkvpacked_func = None
    config = BertConfig.from_json_file(str(bundle / "config.json"))
    tokenizer = PreTrainedTokenizerFast(
        tokenizer_file=str(bundle / "tokenizer.json"),
        unk_token="[UNK]",
        cls_token="[CLS]",
        sep_token="[SEP]",
        pad_token="[PAD]",
        mask_token="[MASK]",
    )
    if config.hidden_size != 768 or tokenizer.vocab_size != 4_096:
        raise NativeContractError("DNABERT-2 model/tokenizer identity differs")
    model = bert_layers.BertForMaskedLM(config)
    checkpoint = root / record["weights_path"] / record["checkpoint_relative_path"]
    state = load_file(str(checkpoint), device="cpu")
    incompatible = model.load_state_dict(state, strict=False)
    if (
        incompatible.missing_keys != ["cls.predictions.decoder.weight"]
        or incompatible.unexpected_keys
        or model.cls.predictions.decoder.weight
        is not model.bert.embeddings.word_embeddings.weight
    ):
        raise NativeContractError("DNABERT-2 checkpoint restore differs")
    del state
    model.to(device).eval()

    def embed(sequence: str, pool_start0: int, pool_end0: int) -> np.ndarray:
        _, offsets = dnabert_encoding(tokenizer, sequence)
        encoded = tokenizer(
            sequence,
            add_special_tokens=True,
            return_attention_mask=True,
            truncation=False,
            return_tensors="pt",
        )
        indices = overlapping_indices(offsets, pool_start0, pool_end0)
        if not indices:
            raise NativeContractError("DNABERT-2 pool has no tokens")
        with torch.inference_mode():
            hidden, _ = model.bert(
                input_ids=encoded["input_ids"].to(device),
                attention_mask=encoded["attention_mask"].to(device),
                output_all_encoded_layers=False,
            )
            pooled = hidden[0, indices].float().mean(dim=0)
        return pooled.cpu().numpy()

    return embed, 768, {
        "attention_backend": "upstream_explicit_pytorch_fallback",
        "archived_triton_kernel_executed": False,
        "checkpoint_tied_decoder_alias_validated": True,
    }


def load_nucleotide_transformer(
    root: Path, record: dict[str, object], device: torch.device
) -> tuple[Callable[[str, int, int], np.ndarray], int, dict[str, object]]:
    source = root / record["source_path"] / "sources"
    config = EsmConfig.from_json_file(str(source / "nt_hf_config"))
    tokenizer_settings = json.loads(
        (source / "nt_hf_tokenizer_config").read_text(encoding="utf-8")
    )
    tokenizer = EsmTokenizer(
        str(source / "nt_hf_vocab"), eos_token=tokenizer_settings["eos_token"]
    )
    if (
        config.hidden_size != 1_280
        or config.max_position_embeddings != 1_002
        or tokenizer_settings.get("model_max_length") != 1_000
        or tokenizer.vocab_size != 4_107
        or tokenizer.cls_token_id != 3
        or tokenizer.eos_token_id is not None
    ):
        raise NativeContractError("Nucleotide Transformer model/tokenizer identity differs")
    model = EsmForMaskedLM(config)
    checkpoint = root / record["weights_path"] / record["checkpoint_relative_path"]
    state = load_file(str(checkpoint), device="cpu")
    legacy = state.pop("esm.embeddings.position_ids", None)
    if (
        legacy is None
        or tuple(legacy.shape) != (1, 1_002)
        or not torch.equal(legacy, torch.arange(1_002, dtype=torch.int64).unsqueeze(0))
    ):
        raise NativeContractError("Nucleotide Transformer legacy position IDs differ")
    incompatible = model.load_state_dict(state, strict=True)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        raise NativeContractError("Nucleotide Transformer checkpoint restore differs")
    del state
    model.to(device).eval()

    def embed(sequence: str, pool_start0: int, pool_end0: int) -> np.ndarray:
        phased = [sequence[start:end] for start, end in nt_phase_spans(len(sequence))]
        encoded = tokenizer(
            phased,
            add_special_tokens=True,
            padding=False,
            truncation=False,
            return_tensors="pt",
        )
        if encoded["input_ids"].shape != (6, 682) or int(encoded["input_ids"].max()) >= 4_105:
            raise NativeContractError("Nucleotide Transformer phase tensor differs")
        with torch.inference_mode():
            hidden = model.esm(
                input_ids=encoded["input_ids"].to(device),
                attention_mask=encoded["attention_mask"].to(device),
                return_dict=True,
            ).last_hidden_state
            phase_values = [
                hidden[phase, nt_pool_indices(phase, pool_start0, pool_end0)]
                .float()
                .mean(dim=0)
                for phase in range(6)
            ]
            pooled = torch.stack(phase_values).mean(dim=0)
        return pooled.cpu().numpy()

    return embed, 1_280, {
        "phase_count": 6,
        "per_phase_context_bp": 4_086,
        "tokens_per_phase_including_cls": 682,
        "legacy_nonpersistent_position_ids_validated": True,
    }


def load_hyenadna(
    root: Path, record: dict[str, object], device: torch.device
) -> tuple[Callable[[str, int, int], np.ndarray], int, dict[str, object]]:
    code = root / record["code_path"] / "code/hyenadna"
    sys.path.insert(0, str(code.parent))
    from hyenadna.configuration_hyena import HyenaConfig
    from hyenadna.modeling_hyena import HyenaDNAForCausalLM

    config = HyenaConfig.from_json_file(str(code / "config.json"))
    if config.d_model != 256 or config.max_seq_len != 1_000_002:
        raise NativeContractError("HyenaDNA model identity differs")
    model = HyenaDNAForCausalLM(config)
    checkpoint = root / record["weights_path"] / record["checkpoint_relative_path"]
    state = load_file(str(checkpoint), device="cpu")
    incompatible = model.load_state_dict(state, strict=False)
    aliases = {
        f"hyena.backbone.layers.{layer}.mixer.filter_fn.implicit_filter.{position}.freq"
        for layer in range(config.n_layer)
        for position in (3, 5)
    }
    shared = all(
        model.hyena.backbone.layers[layer].mixer.filter_fn.implicit_filter[1]
        is model.hyena.backbone.layers[layer].mixer.filter_fn.implicit_filter[3]
        is model.hyena.backbone.layers[layer].mixer.filter_fn.implicit_filter[5]
        for layer in range(config.n_layer)
    )
    if set(incompatible.missing_keys) != aliases or incompatible.unexpected_keys or not shared:
        raise NativeContractError("HyenaDNA checkpoint restore differs")
    del state
    model.to(device).eval()

    def embed(sequence: str, pool_start0: int, pool_end0: int) -> np.ndarray:
        ids = torch.tensor(
            [[TOKEN_IDS[base] for base in sequence] + [1]], dtype=torch.long, device=device
        )
        with torch.inference_mode():
            hidden = model.hyena(input_ids=ids, return_dict=True).last_hidden_state
            pooled = hidden[0, pool_start0:pool_end0].float().mean(dim=0)
        return pooled.cpu().numpy()

    return embed, 256, {
        "tokens_including_terminal_sep": 4_097,
        "shared_activation_modules_identity_validated": True,
        "safetensors_omitted_shared_parameter_aliases": sorted(aliases),
    }


def extract(
    *,
    root: Path,
    config_path: Path,
    model_id: str,
    output: Path,
    seed: int,
) -> dict[str, object]:
    if model_id not in MODELS or output.exists() or config_path.is_symlink():
        raise NativeContractError("embedding extraction request differs")
    config = load_config(config_path)
    fixture = root / config["fixture"]["path"]
    rows, fasta = read_fixture(fixture)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.use_deterministic_algorithms(True)
    device = require_cuda()
    record = config["models"][model_id]
    loaders = {
        "dnabert2": load_dnabert,
        "nucleotide_transformer": load_nucleotide_transformer,
        "hyenadna": load_hyenadna,
    }
    embed, width, native = loaders[model_id](root, record, device)
    expected_width = int(record["hidden_width"])
    if width != expected_width:
        raise NativeContractError("model hidden width differs")
    output.mkdir(parents=True)
    embeddings = np.empty((len(rows), len(ALLELES), width), dtype=np.float32)
    first_input: tuple[str, int, int] | None = None
    first_output: np.ndarray | None = None
    for row_index, row in enumerate(rows):
        for allele_index, allele in enumerate(ALLELES):
            sequence = fasta[f"{row['fasta_record_prefix']}|{allele}"]
            value = embed(sequence, int(row["pool_start0"]), int(row["pool_end0"]))
            if value.shape != (width,) or not np.isfinite(value).all():
                raise NativeContractError("pooled embedding differs")
            embeddings[row_index, allele_index] = value
            if first_input is None:
                first_input = (sequence, int(row["pool_start0"]), int(row["pool_end0"]))
                first_output = value
    assert first_input is not None and first_output is not None
    repeated = embed(*first_input)
    reproducibility_max_abs_diff = float(np.max(np.abs(repeated - first_output)))
    if reproducibility_max_abs_diff > 1e-6:
        raise NativeContractError("repeated pooled embedding differs")
    np.savez_compressed(
        output / "allele_embeddings.npz",
        fixture_ids=np.asarray([row["fixture_id"] for row in rows]),
        outer_locus_sequence_group_ids=np.asarray(
            [row["outer_locus_sequence_group_id"] for row in rows]
        ),
        outer_folds=np.asarray([int(row["outer_fold"]) for row in rows], dtype=np.int8),
        allele_order=np.asarray(ALLELES),
        embeddings=embeddings,
    )
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-gse281364-dna-lm-allele-embeddings-v1",
        "status": "pass_outcome_blind_embedding_extraction",
        "dataset_id": "gse281364",
        "model_id": model_id,
        "fixture_artifacts_sha256": config["fixture"]["artifacts_sha256"],
        "checkpoint_sha256": record["checkpoint_sha256"],
        "successful_probe_artifacts_sha256": record[
            "successful_probe_artifacts_sha256"
        ],
        "license": record["license"],
        "restricted_comparator": record["restricted_comparator"],
        "open_champion_eligible_after_task_gates": record[
            "open_champion_eligible_after_task_gates"
        ],
        "elements": len(rows),
        "outer_locus_sequence_groups": len(rows),
        "outer_folds": 5,
        "alleles": list(ALLELES),
        "sequence_or_phase_forwards": len(rows)
        * len(ALLELES)
        * (6 if model_id == "nucleotide_transformer" else 1)
        + (6 if model_id == "nucleotide_transformer" else 1),
        "hidden_width": width,
        "pool_start0": 1_280,
        "pool_end0": 2_816,
        "pool_length_bp": 1_536,
        "reverse_complement_aggregation_deferred_to_fold_projection": True,
        "native_contract": native,
        "reproducibility_max_abs_diff": reproducibility_max_abs_diff,
        "checkpoint_loaded_with_safetensors": True,
        "runtime_network_allowed": False,
        "reporter_counts_read": False,
        "reporter_outcomes_read": False,
        "sealed_labels_read": False,
        "downstream_head_fit": False,
        "device": str(device),
        "torch": torch.__version__,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--model-id", choices=MODELS, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260824)
    arguments = parser.parse_args()
    extract(
        root=arguments.root,
        config_path=arguments.config,
        model_id=arguments.model_id,
        output=arguments.output,
        seed=arguments.seed,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
