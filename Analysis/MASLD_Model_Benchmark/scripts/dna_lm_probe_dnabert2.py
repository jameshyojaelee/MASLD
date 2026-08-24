#!/usr/bin/env python3
"""Run exact DNABERT-2 code and weights on the common outcome-free fixture."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import random
import sys

import numpy as np
import torch
from safetensors.torch import load_file
from transformers import PreTrainedTokenizerFast


class Dnabert2ProbeError(ValueError):
    """Raised when the exact DNABERT-2 runtime contract differs."""


def _sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _sha256_text(value: str) -> str:
    return sha256(value.encode("ascii")).hexdigest()


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
                if not name or name in records:
                    raise Dnabert2ProbeError("DNABERT-2 FASTA header differs")
                records[name] = []
            elif name is None:
                raise Dnabert2ProbeError("DNABERT-2 FASTA sequence precedes header")
            else:
                records[name].append(line)
    values = {key: "".join(parts) for key, parts in records.items()}
    if any(len(value) != 6000 or set(value) - set("ACGT") for value in values.values()):
        raise Dnabert2ProbeError("DNABERT-2 FASTA alphabet or length differs")
    return values


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["model_id"] == "dnabert2" and row["context_id"] == "common_6000"
        ]
    if len(rows) != 3 or len({row["genomic_fold"] for row in rows}) != 3:
        raise Dnabert2ProbeError("DNABERT-2 common fixture census differs")
    return rows


def _tokenize(
    tokenizer: PreTrainedTokenizerFast, sequence: str
) -> tuple[torch.Tensor, torch.Tensor, list[tuple[int, int]]]:
    encoded = tokenizer(
        sequence,
        add_special_tokens=True,
        return_attention_mask=True,
        return_offsets_mapping=True,
        truncation=False,
        return_tensors="pt",
    )
    ids = encoded["input_ids"]
    mask = encoded["attention_mask"]
    offsets = [(int(start), int(end)) for start, end in encoded["offset_mapping"][0]]
    if (
        ids.shape[0] != 1
        or ids.shape != mask.shape
        or len(offsets) != ids.shape[1]
        or offsets[0] != (0, 0)
        or offsets[-1] != (0, 0)
        or any(token_id >= 4096 for token_id in ids.flatten().tolist())
    ):
        raise Dnabert2ProbeError("DNABERT-2 BPE tokenization differs")
    cursor = 0
    for start, end in offsets[1:-1]:
        if start != cursor or end <= start:
            raise Dnabert2ProbeError("DNABERT-2 BPE offsets do not tile sequence")
        cursor = end
    if cursor != len(sequence):
        raise Dnabert2ProbeError("DNABERT-2 BPE offsets do not cover sequence")
    return ids, mask, offsets


def _affected_interval(
    ref_offsets: list[tuple[int, int]],
    alt_offsets: list[tuple[int, int]],
    variant_index0: int,
) -> tuple[int, int]:
    start, end = variant_index0, variant_index0 + 1
    changed = True
    while changed:
        changed = False
        for offsets in (ref_offsets, alt_offsets):
            for token_start, token_end in offsets:
                if token_end > start and token_start < end:
                    new_start = min(start, token_start)
                    new_end = max(end, token_end)
                    if (new_start, new_end) != (start, end):
                        start, end = new_start, new_end
                        changed = True
    return start, end


def _indices_overlapping(
    offsets: list[tuple[int, int]], start: int, end: int
) -> list[int]:
    return [
        index
        for index, (token_start, token_end) in enumerate(offsets)
        if token_end > start and token_start < end
    ]


def _allele_masked_log_probability(
    model: torch.nn.Module,
    ids: torch.Tensor,
    attention_mask: torch.Tensor,
    token_indices: list[int],
    mask_token_id: int,
    device: torch.device,
) -> float:
    batch_ids = ids.repeat(len(token_indices), 1)
    batch_mask = attention_mask.repeat(len(token_indices), 1)
    labels = torch.full_like(batch_ids, -100)
    targets = []
    for batch_index, token_index in enumerate(token_indices):
        target = int(batch_ids[batch_index, token_index])
        if target <= 4:
            raise Dnabert2ProbeError("DNABERT-2 affected target is special")
        targets.append(target)
        labels[batch_index, token_index] = target
        batch_ids[batch_index, token_index] = mask_token_id
    with torch.inference_mode():
        output = model(
            input_ids=batch_ids.to(device),
            attention_mask=batch_mask.to(device),
            labels=labels.to(device),
            return_dict=True,
        )
        logits = output.logits.float()
        selected = torch.log_softmax(
            logits[
                torch.arange(len(token_indices), device=device),
                torch.tensor(token_indices, device=device),
            ],
            dim=-1,
        )
        values = selected[
            torch.arange(len(token_indices), device=device),
            torch.tensor(targets, device=device),
        ]
    return float(values.sum().cpu())


def _allele_features(
    model: torch.nn.Module,
    ids: torch.Tensor,
    attention_mask: torch.Tensor,
    offsets: list[tuple[int, int]],
    pool_start0: int,
    pool_end0: int,
    device: torch.device,
) -> np.ndarray:
    pool_indices = _indices_overlapping(offsets, pool_start0, pool_end0)
    if not pool_indices:
        raise Dnabert2ProbeError("DNABERT-2 pool has no tokens")
    with torch.inference_mode():
        hidden, _ = model.bert(
            input_ids=ids.to(device),
            attention_mask=attention_mask.to(device),
            output_all_encoded_layers=False,
        )
        pooled = hidden[0, pool_indices].float().mean(dim=0)
    return pooled.cpu().numpy()


def _score_pair(
    model: torch.nn.Module,
    tokenizer: PreTrainedTokenizerFast,
    reference: str,
    alternative: str,
    variant_index0: int,
    device: torch.device,
) -> tuple[float, dict[str, object]]:
    ref_ids, ref_mask, ref_offsets = _tokenize(tokenizer, reference)
    alt_ids, alt_mask, alt_offsets = _tokenize(tokenizer, alternative)
    start, end = _affected_interval(ref_offsets, alt_offsets, variant_index0)
    ref_indices = _indices_overlapping(ref_offsets, start, end)
    alt_indices = _indices_overlapping(alt_offsets, start, end)
    if not ref_indices or not alt_indices or not (start <= variant_index0 < end):
        raise Dnabert2ProbeError("DNABERT-2 affected BPE interval differs")
    ref_value = _allele_masked_log_probability(
        model, ref_ids, ref_mask, ref_indices, tokenizer.mask_token_id, device
    )
    alt_value = _allele_masked_log_probability(
        model, alt_ids, alt_mask, alt_indices, tokenizer.mask_token_id, device
    )
    width = end - start
    return (alt_value - ref_value) / width, {
        "affected_start0": start,
        "affected_end0": end,
        "affected_bp": width,
        "reference_affected_tokens": len(ref_indices),
        "alternative_affected_tokens": len(alt_indices),
        "reference_total_tokens": int(ref_ids.shape[1]),
        "alternative_total_tokens": int(alt_ids.shape[1]),
    }


def run(
    code_root: Path,
    checkpoint: Path,
    fixture: Path,
    output: Path,
    seed: int,
) -> dict[str, object]:
    if output.exists() or code_root.is_symlink() or checkpoint.is_symlink() or fixture.is_symlink():
        raise Dnabert2ProbeError("DNABERT-2 runtime request differs")
    output.mkdir(mode=0o750)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.use_deterministic_algorithms(True)

    sys.path.insert(0, str(code_root.parent))
    from bundle import bert_layers
    from bundle.configuration_bert import BertConfig

    # The archived upstream Triton kernel targets the 2023 Triton API and
    # fails closed under the admitted Triton 2.3 runtime. Use the unmodified
    # model's explicit PyTorch fallback instead of patching that kernel.
    bert_layers.flash_attn_qkvpacked_func = None
    BertForMaskedLM = bert_layers.BertForMaskedLM
    config = BertConfig.from_json_file(str(code_root / "config.json"))
    tokenizer = PreTrainedTokenizerFast(
        tokenizer_file=str(code_root / "tokenizer.json"),
        unk_token="[UNK]",
        cls_token="[CLS]",
        sep_token="[SEP]",
        pad_token="[PAD]",
        mask_token="[MASK]",
    )
    if config.vocab_size != 4096 or tokenizer.vocab_size != 4096 or tokenizer.mask_token_id != 4:
        raise Dnabert2ProbeError("DNABERT-2 model/tokenizer config differs")
    model = BertForMaskedLM(config)
    state = load_file(str(checkpoint), device="cpu")
    incompatible = model.load_state_dict(state, strict=False)
    tied_alias = "cls.predictions.decoder.weight"
    tied_identity = (
        model.cls.predictions.decoder.weight
        is model.bert.embeddings.word_embeddings.weight
    )
    if (
        incompatible.missing_keys != [tied_alias]
        or incompatible.unexpected_keys
        or not tied_identity
    ):
        raise Dnabert2ProbeError(f"DNABERT-2 checkpoint keys differ: {incompatible!r}")
    del state
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise Dnabert2ProbeError("DNABERT-2 production probe requires CUDA")
    model.to(device).eval()

    fasta = _read_fasta(fixture / "fixture" / "common_6000.alleles.fa.gz")
    rows = _read_manifest(fixture / "fixture" / "sequence_manifest.tsv")
    score_rows = []
    feature_ids = []
    features = []
    allele_delta_ids = []
    allele_deltas = []
    first_pair: tuple[str, str, int] | None = None
    first_value: float | None = None
    for row in rows:
        fixture_id = row["fixture_id"]
        expected_hashes = {
            "REF": row["reference_sequence_sha256"],
            "ALT": row["alternative_sequence_sha256"],
            "REF_RC": row["reverse_complement_reference_sha256"],
            "ALT_RC": row["reverse_complement_alternative_sha256"],
        }
        sequences = {}
        allele_features = {}
        for allele in ("REF", "ALT", "REF_RC", "ALT_RC"):
            record_id = f"{fixture_id}|common_6000|{allele}"
            sequence = fasta.get(record_id)
            if sequence is None or _sha256_text(sequence) != expected_hashes[allele]:
                raise Dnabert2ProbeError("DNABERT-2 FASTA identity differs")
            sequences[allele] = sequence
            ids, mask, offsets = _tokenize(tokenizer, sequence)
            value = _allele_features(
                model,
                ids,
                mask,
                offsets,
                int(row["pool_start0"]),
                int(row["pool_end0"]),
                device,
            )
            allele_features[allele] = value
            feature_ids.append(record_id)
            features.append(value)
        forward, forward_geometry = _score_pair(
            model,
            tokenizer,
            sequences["REF"],
            sequences["ALT"],
            int(row["forward_variant_index0"]),
            device,
        )
        reverse, reverse_geometry = _score_pair(
            model,
            tokenizer,
            sequences["REF_RC"],
            sequences["ALT_RC"],
            int(row["reverse_complement_variant_index0"]),
            device,
        )
        if first_pair is None:
            first_pair = (
                sequences["REF"],
                sequences["ALT"],
                int(row["forward_variant_index0"]),
            )
            first_value = forward
        score_rows.append(
            {
                "fixture_id": fixture_id,
                "genomic_fold": int(row["genomic_fold"]),
                "forward_alt_minus_ref_per_affected_bp": forward,
                "reverse_alt_minus_ref_per_affected_bp": reverse,
                "rc_averaged_alt_minus_ref_per_affected_bp": (forward + reverse) / 2.0,
                "forward_geometry_json": json.dumps(forward_geometry, sort_keys=True),
                "reverse_geometry_json": json.dumps(reverse_geometry, sort_keys=True),
            }
        )
        allele_delta_ids.append(fixture_id)
        allele_deltas.append(
            (
                allele_features["ALT"]
                - allele_features["REF"]
                + allele_features["ALT_RC"]
                - allele_features["REF_RC"]
            )
            / 2.0
        )

    assert first_pair is not None and first_value is not None
    repeated, _ = _score_pair(model, tokenizer, *first_pair, device)
    reproducibility_abs_diff = abs(repeated - first_value)
    if reproducibility_abs_diff > 1e-6:
        raise Dnabert2ProbeError("DNABERT-2 repeated forward differs")
    with (output / "native_scores.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(score_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(score_rows)
    np.savez_compressed(
        output / "pooled_features.npz",
        feature_ids=np.asarray(feature_ids),
        features=np.stack(features),
        allele_delta_ids=np.asarray(allele_delta_ids),
        allele_deltas=np.stack(allele_deltas),
    )
    receipt = {
        "schema_version": "masld-bench-dnabert2-runtime-probe-v1",
        "status": "pass",
        "seed": seed,
        "checkpoint_sha256": _sha256_file(checkpoint),
        "model_parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "fixture_count": len(rows),
        "feature_width": int(features[0].shape[0]),
        "pool_length_bp": 1536,
        "common_input_length_bp": 6000,
        "native_score": "allele_specific_BPE_coordinate_mapped_masked_pseudo_log_likelihood_ALT_minus_REF_per_affected_bp_with_RC_average",
        "feature_contract": "1536bp_coordinate_mapped_mean_pool_REF_ALT_and_RC_plus_RC_averaged_allele_delta",
        "open_champion_eligible_after_task_gates": True,
        "attention_backend": "upstream_explicit_pytorch_fallback",
        "archived_triton_kernel_executed": False,
        "model_source_code_modified": False,
        "tied_embedding_decoder_identity_validated": tied_identity,
        "embeddings_biologically_meaningful_without_trained_head": False,
        "head_fit": False,
        "observed_outcomes_loaded": False,
        "sealed_outcomes_loaded": False,
        "runtime_network_allowed": False,
        "reproducibility_abs_diff": reproducibility_abs_diff,
        "device": str(device),
        "torch": torch.__version__,
    }
    (output / "runtime_probe_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--code-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260824)
    args = parser.parse_args()
    run(args.code_root, args.checkpoint, args.fixture, args.output, args.seed)


if __name__ == "__main__":
    main()
