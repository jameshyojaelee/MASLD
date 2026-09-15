#!/usr/bin/env python3
"""Mean-pool DNA-language-model embeddings for matched H3K27ac liver region pairs.

Every model is constructed from the frozen vendored custom code and its `config.json`,
then restored from a local safetensors checkpoint.  Nothing calls `AutoConfig`,
`AutoModel`, or `AutoTokenizer`, so no Hugging Face network access or hub cache is
needed, and the DNA tokenizers are plain character maps rather than the compiled
`tokenizers` extension.

`--mode random` builds the identical architecture and skips the checkpoint restore, so
the random-init twin differs from the pretrained arm only in the learned weights.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
import sys
import time

import numpy as np
import torch
from safetensors.torch import load_file


MODELS = ("hyenadna", "caduceus", "dnabert2")
MODES = ("pretrained", "random")
TOKEN_IDS = {"A": 7, "C": 8, "G": 9, "T": 10, "N": 11}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
CODE_ROOT = "executions/dna-lm-custom-code-21064626/code"
CHECKPOINTS = {
    "hyenadna": "executions/dna-lm-open-safetensors-21064556/checkpoints/hyenadna.model.safetensors",
    "caduceus": "executions/dna-lm-open-safetensors-21064556/checkpoints/caduceus.model.safetensors",
    "dnabert2": "executions/dna-lm-weights-only-dnabert2-21064812/converted/dnabert2.model.safetensors",
}
WIDTHS = {"hyenadna": 256, "caduceus": 256, "dnabert2": 768}


class RegionEmbeddingError(RuntimeError):
    """Raised when a model identity, checkpoint restore, or pooling requirement differs."""


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def build_hyenadna(root: Path, mode: str, seed: int, device: torch.device):
    code = root / CODE_ROOT / "hyenadna"
    sys.path.insert(0, str(code.parent))
    from hyenadna.configuration_hyena import HyenaConfig
    from hyenadna.modeling_hyena import HyenaDNAForCausalLM

    config = HyenaConfig.from_json_file(str(code / "config.json"))
    if config.d_model != 256 or config.max_seq_len != 1_000_002:
        raise RegionEmbeddingError("HyenaDNA model identity differs")
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    model = HyenaDNAForCausalLM(config)
    restore: dict[str, object] = {"checkpoint_loaded": False}
    if mode == "pretrained":
        state = load_file(str(root / CHECKPOINTS["hyenadna"]), device="cpu")
        incompatible = model.load_state_dict(state, strict=False)
        aliases = {
            f"hyena.backbone.layers.{layer}.mixer.filter_fn.implicit_filter.{position}.freq"
            for layer in range(config.n_layer)
            for position in (3, 5)
        }
        if set(incompatible.missing_keys) != aliases or incompatible.unexpected_keys:
            raise RegionEmbeddingError("HyenaDNA checkpoint restore differs")
        del state
        restore = {"checkpoint_loaded": True, "shared_parameter_aliases": len(aliases)}
    model.to(device).eval()

    def embed(batch: list[str]) -> np.ndarray:
        length = len(batch[0])
        ids = torch.tensor(
            [[TOKEN_IDS[base] for base in sequence] + [1] for sequence in batch],
            dtype=torch.long,
            device=device,
        )
        with torch.inference_mode():
            hidden = model.hyena(input_ids=ids, return_dict=True).last_hidden_state
            if hidden.shape[1] != length + 1:
                raise RegionEmbeddingError("HyenaDNA token count differs")
            # Pool the base positions only; the terminal SEP is not sequence.
            return hidden[:, :length].float().mean(dim=1).cpu().numpy()

    return embed, restore


def build_caduceus(root: Path, mode: str, seed: int, device: torch.device):
    code = root / CODE_ROOT / "caduceus"
    sys.path.insert(0, str(code.parent))
    from caduceus.configuration_caduceus import CaduceusConfig
    from caduceus.modeling_caduceus import CaduceusForMaskedLM

    config = CaduceusConfig.from_json_file(str(code / "config.json"))
    if config.d_model != 256 or config.n_layer != 16:
        raise RegionEmbeddingError("Caduceus model identity differs")
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    model = CaduceusForMaskedLM(config)
    restore: dict[str, object] = {"checkpoint_loaded": False}
    if mode == "pretrained":
        state = load_file(str(root / CHECKPOINTS["caduceus"]), device="cpu")
        incompatible = model.load_state_dict(state, strict=False)
        aliases = {
            f"caduceus.backbone.layers.{layer}.mixer.submodule.mamba_rev.{projection}.weight"
            for layer in range(16)
            for projection in ("in_proj", "out_proj")
        } | {"lm_head.lm_head.weight"}
        if set(incompatible.missing_keys) != aliases or incompatible.unexpected_keys:
            raise RegionEmbeddingError("Caduceus checkpoint restore differs")
        del state
        restore = {"checkpoint_loaded": True, "shared_parameter_aliases": len(aliases)}
    model.to(device).eval()

    def embed(batch: list[str]) -> np.ndarray:
        length = len(batch[0])
        ids = torch.tensor(
            [[TOKEN_IDS[base] for base in sequence] for sequence in batch],
            dtype=torch.long,
            device=device,
        )
        with torch.inference_mode():
            hidden = model.caduceus(input_ids=ids, return_dict=True).last_hidden_state
            if hidden.shape[1:] != (length, 512):
                raise RegionEmbeddingError("Caduceus hidden shape differs")
            # Reverse-complement parameter sharing: fold the two 256-wide halves back
            # onto one strand-symmetric representation, as the frozen probe does.
            realigned = (
                hidden[..., :256] + torch.flip(hidden[..., 256:], dims=(-2, -1))
            ) / 2
            return realigned.float().mean(dim=1).cpu().numpy()

    return embed, restore


def build_dnabert2(root: Path, mode: str, seed: int, device: torch.device):
    bundle = root / "executions/dnabert2-runtime-bundle-21064834/bundle"
    sys.path.insert(0, str(bundle.parent))
    from bundle import bert_layers
    from bundle.configuration_bert import BertConfig
    from transformers import PreTrainedTokenizerFast

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
        raise RegionEmbeddingError("DNABERT-2 model/tokenizer identity differs")
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    model = bert_layers.BertForMaskedLM(config)
    restore: dict[str, object] = {"checkpoint_loaded": False}
    if mode == "pretrained":
        state = load_file(str(root / CHECKPOINTS["dnabert2"]), device="cpu")
        incompatible = model.load_state_dict(state, strict=False)
        if (
            incompatible.missing_keys != ["cls.predictions.decoder.weight"]
            or incompatible.unexpected_keys
            or model.cls.predictions.decoder.weight
            is not model.bert.embeddings.word_embeddings.weight
        ):
            raise RegionEmbeddingError("DNABERT-2 checkpoint restore differs")
        del state
        restore = {"checkpoint_loaded": True, "tied_decoder_alias_validated": True}
    model.to(device).eval()

    def embed(batch: list[str]) -> np.ndarray:
        encoded = tokenizer(
            batch,
            add_special_tokens=True,
            return_attention_mask=True,
            padding=True,
            truncation=False,
            return_tensors="pt",
        )
        ids = encoded["input_ids"].to(device)
        mask = encoded["attention_mask"].to(device)
        with torch.inference_mode():
            hidden, _ = model.bert(
                input_ids=ids, attention_mask=mask, output_all_encoded_layers=False
            )
            # Mean over real sequence tokens only: drop [CLS], [SEP], and padding.
            weight = mask.clone()
            weight[:, 0] = 0
            weight[ids == tokenizer.sep_token_id] = 0
            weight = weight.unsqueeze(-1).float()
            total = weight.sum(dim=1).clamp(min=1.0)
            return ((hidden * weight).sum(dim=1) / total).float().cpu().numpy()

    return embed, restore


BUILDERS = {"hyenadna": build_hyenadna, "caduceus": build_caduceus, "dnabert2": build_dnabert2}


def extract(*, root: Path, pairs: Path, model_id: str, mode: str, output: Path, seed: int, batch: int) -> dict:
    if model_id not in MODELS or mode not in MODES:
        raise RegionEmbeddingError("model or mode differs")
    if output.exists():
        raise RegionEmbeddingError(f"refusing to overwrite: {output}")
    with np.load(pairs, allow_pickle=False) as data:
        sequences = data["seq"].astype(str)
        carried = {key: data[key] for key in ("label", "chrom", "region") if key in data.files}
        extra = {key: data[key] for key in ("row_index", "width") if key in data.files}
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1 or set("".join(sequences[:64])) - set("ACGT"):
        raise RegionEmbeddingError("region sequence geometry differs")
    if set(carried) != {"label", "chrom", "region"}:
        raise RegionEmbeddingError("carried metadata missing")

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RegionEmbeddingError("CUDA is required")
    embed, restore = BUILDERS[model_id](root, mode, seed, device)

    width = WIDTHS[model_id]
    forward = np.empty((len(sequences), width), dtype=np.float32)
    # Caduceus is strand-symmetric by construction (reverse-complement parameter
    # sharing), so a second pass would be redundant there.
    wants_rc = model_id != "caduceus"
    reverse = np.empty((len(sequences), width), dtype=np.float32) if wants_rc else None
    start = time.time()
    for index in range(0, len(sequences), batch):
        chunk = [str(value) for value in sequences[index : index + batch]]
        forward[index : index + len(chunk)] = embed(chunk)
        if wants_rc:
            reverse[index : index + len(chunk)] = embed([reverse_complement(s) for s in chunk])
        if index % (batch * 100) == 0:
            print(f"{model_id}/{mode} {index}/{len(sequences)} {time.time() - start:.0f}s", flush=True)
    if not np.isfinite(forward).all() or (wants_rc and not np.isfinite(reverse).all()):
        raise RegionEmbeddingError("non-finite embedding")

    payload = {"emb": forward, **carried, **extra}
    if wants_rc:
        payload["emb_rc"] = reverse
        payload["emb_rcmean"] = ((forward + reverse) / 2.0).astype(np.float32)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **payload)

    receipt = {
        "schema_version": "masld-bench-chromatin-region-dna-lm-embeddings-v1",
        "status": "pass_outcome_blind_region_embedding_extraction",
        "model_id": model_id,
        "mode": mode,
        "checkpoint_restore": restore,
        "sequences": int(len(sequences)),
        "sequence_length_bp": int(next(iter(lengths))),
        "embedding_width": width,
        "pooling": "mean_over_sequence_positions",
        "arrays": sorted(payload),
        "reverse_complement_pass": bool(wants_rc),
        "strand_symmetric_by_construction": model_id == "caduceus",
        "seed": seed,
        "batch_size": batch,
        "device": str(device),
        "torch": torch.__version__,
        "elapsed_seconds": round(time.time() - start, 1),
        "hub_network_used": False,
        "labels_read_for_fitting": False,
        "model_fit": False,
    }
    output.with_suffix(".receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--model", required=True, choices=MODELS)
    parser.add_argument("--mode", required=True, choices=MODES)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260905)
    parser.add_argument("--batch", type=int, default=32)
    arguments = parser.parse_args()
    extract(
        root=arguments.root.resolve(strict=True),
        pairs=arguments.pairs.resolve(strict=True),
        model_id=arguments.model,
        mode=arguments.mode,
        output=arguments.output,
        seed=arguments.seed,
        batch=arguments.batch,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
