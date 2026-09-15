#!/usr/bin/env python3
"""Run Nucleotide Transformer on the six-phase outcome-free fixture."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import random

import numpy as np
import torch
from safetensors.torch import load_file
from transformers import EsmConfig, EsmForMaskedLM, EsmTokenizer


class NucleotideTransformerProbeError(ValueError):
    """Raised when the NT fixture or exact runtime requirement differs."""


PHASE_CONTEXT_BP = 5994


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
                    raise NucleotideTransformerProbeError("FASTA header differs")
                records[name] = []
            elif name is None:
                raise NucleotideTransformerProbeError("FASTA sequence precedes header")
            else:
                records[name].append(line)
    values = {key: "".join(parts) for key, parts in records.items()}
    if any(len(value) != 6000 or set(value) - set("ACGT") for value in values.values()):
        raise NucleotideTransformerProbeError("NT common FASTA alphabet or length differs")
    return values


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["model_id"] == "nucleotide_transformer"
            and row["context_id"] == "common_6000"
        ]
    if len(rows) != 3 or len({row["genomic_fold"] for row in rows}) != 3:
        raise NucleotideTransformerProbeError("NT common fixture census differs")
    return rows


def _tokenize_phases(
    tokenizer: EsmTokenizer, sequence: str
) -> tuple[torch.Tensor, torch.Tensor, list[list[str]]]:
    # 999 non-overlapping 6-mers plus <cls> exactly matches the published
    # 1,000-token training length. Each fixed phase trims only fixture edges.
    phased = [sequence[phase : phase + PHASE_CONTEXT_BP] for phase in range(6)]
    encoded = tokenizer(
        phased,
        add_special_tokens=True,
        padding=True,
        truncation=False,
        return_tensors="pt",
    )
    if int(encoded["input_ids"].max()) >= 4105:
        raise NucleotideTransformerProbeError("NT input uses tokenizer-only EOS/BOS IDs")
    tokens = [
        tokenizer.convert_ids_to_tokens(
            encoded["input_ids"][phase][encoded["attention_mask"][phase].bool()].tolist()
        )
        for phase in range(6)
    ]
    for phase, observed in enumerate(tokens):
        expected = ["<cls>"] + [
            sequence[start : start + 6]
            for start in range(phase, phase + PHASE_CONTEXT_BP, 6)
        ]
        if observed != expected or len(observed) != 1000:
            raise NucleotideTransformerProbeError("NT six-phase tokenization differs")
    return encoded["input_ids"], encoded["attention_mask"], tokens


def _pool_hidden(
    hidden: torch.Tensor,
    tokens: list[list[str]],
    pool_start0: int,
    pool_end0: int,
) -> np.ndarray:
    outputs = []
    for phase in range(6):
        positions = []
        base = phase
        for token_index, token in enumerate(tokens[phase][1:], start=1):
            width = len(token)
            if base < pool_end0 and base + width > pool_start0:
                positions.append(token_index)
            base += width
        if not positions:
            raise NucleotideTransformerProbeError("NT pool has no tokens")
        outputs.append(hidden[phase, positions].float().mean(dim=0).cpu().numpy())
    return np.stack(outputs)


def _unmasked_features(
    model: EsmForMaskedLM,
    tokenizer: EsmTokenizer,
    sequence: str,
    pool_start0: int,
    pool_end0: int,
    device: torch.device,
) -> np.ndarray:
    ids, mask, tokens = _tokenize_phases(tokenizer, sequence)
    with torch.inference_mode():
        hidden = model.esm(
            input_ids=ids.to(device),
            attention_mask=mask.to(device),
            return_dict=True,
        ).last_hidden_state
    return _pool_hidden(hidden, tokens, pool_start0, pool_end0)


def _masked_phase_scores(
    model: EsmForMaskedLM,
    tokenizer: EsmTokenizer,
    reference: str,
    alternative: str,
    device: torch.device,
) -> np.ndarray:
    ref_ids, ref_mask, ref_tokens = _tokenize_phases(tokenizer, reference)
    alt_ids, alt_mask, alt_tokens = _tokenize_phases(tokenizer, alternative)
    if not torch.equal(ref_mask, alt_mask):
        raise NucleotideTransformerProbeError("NT allele attention masks differ")
    positions = []
    ref_targets = []
    alt_targets = []
    masked = ref_ids.clone()
    for phase in range(6):
        active = ref_mask[phase].bool()
        differing = torch.nonzero(
            ref_ids[phase, active] != alt_ids[phase, active], as_tuple=False
        ).flatten()
        if len(differing) != 1:
            raise NucleotideTransformerProbeError("NT allele token difference is not unique")
        position = int(differing.item())
        if len(ref_tokens[phase][position]) != 6 or len(alt_tokens[phase][position]) != 6:
            raise NucleotideTransformerProbeError("NT variant is not inside one 6-mer")
        positions.append(position)
        ref_targets.append(int(ref_ids[phase, position]))
        alt_targets.append(int(alt_ids[phase, position]))
        masked[phase, position] = tokenizer.mask_token_id
    for phase, position in enumerate(positions):
        ref_copy = ref_ids[phase].clone()
        alt_copy = alt_ids[phase].clone()
        ref_copy[position] = tokenizer.mask_token_id
        alt_copy[position] = tokenizer.mask_token_id
        if not torch.equal(ref_copy, alt_copy):
            raise NucleotideTransformerProbeError("NT masked allele contexts differ")
    with torch.inference_mode():
        logits = model(
            input_ids=masked.to(device),
            attention_mask=ref_mask.to(device),
            return_dict=True,
        ).logits.float()
        phase_index = torch.arange(6, device=device)
        selected = torch.log_softmax(
            logits[phase_index, torch.tensor(positions, device=device)], dim=-1
        )
        ref_values = selected[
            phase_index, torch.tensor(ref_targets, device=device)
        ]
        alt_values = selected[
            phase_index, torch.tensor(alt_targets, device=device)
        ]
    return (alt_values - ref_values).cpu().numpy()


def run(
    config_path: Path,
    vocab_path: Path,
    tokenizer_config_path: Path,
    checkpoint: Path,
    fixture: Path,
    output: Path,
    seed: int,
) -> dict[str, object]:
    if output.exists() or any(
        path.is_symlink()
        for path in (
            config_path,
            vocab_path,
            tokenizer_config_path,
            checkpoint,
            fixture,
        )
    ):
        raise NucleotideTransformerProbeError("NT runtime request differs")
    output.mkdir(mode=0o750)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.use_deterministic_algorithms(True)
    config = EsmConfig.from_json_file(str(config_path))
    tokenizer_settings = json.loads(tokenizer_config_path.read_text(encoding="utf-8"))
    if tokenizer_settings != {
        "clean_up_tokenization_spaces": True,
        "eos_token": None,
        "model_max_length": 1000,
        "tokenizer_class": "EsmTokenizer",
    }:
        raise NucleotideTransformerProbeError("NT tokenizer config differs")
    tokenizer = EsmTokenizer(
        str(vocab_path), eos_token=tokenizer_settings["eos_token"]
    )
    if (
        config.model_type != "esm"
        or config.vocab_size != 4105
        or tokenizer.vocab_size != 4107
        or tokenizer.mask_token_id != 2
        or tokenizer.cls_token_id != 3
    ):
        raise NucleotideTransformerProbeError("NT model/tokenizer config differs")
    model = EsmForMaskedLM(config)
    state = load_file(str(checkpoint), device="cpu")
    legacy_position_ids = state.pop("esm.embeddings.position_ids", None)
    if (
        legacy_position_ids is None
        or tuple(legacy_position_ids.shape) != (1, config.max_position_embeddings)
        or legacy_position_ids.dtype != torch.int64
        or not torch.equal(
            legacy_position_ids,
            torch.arange(config.max_position_embeddings, dtype=torch.int64).unsqueeze(0),
        )
    ):
        raise NucleotideTransformerProbeError("NT legacy position-id buffer differs")
    incompatible = model.load_state_dict(state, strict=True)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        raise NucleotideTransformerProbeError(
            f"NT checkpoint keys differ: {incompatible!r}"
        )
    del state
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise NucleotideTransformerProbeError("NT production probe requires CUDA")
    model.to(device).eval()

    fasta = _read_fasta(fixture / "fixture" / "common_6000.alleles.fa.gz")
    rows = _read_manifest(fixture / "fixture" / "sequence_manifest.tsv")
    score_rows = []
    feature_ids = []
    features = []
    allele_delta_ids = []
    allele_deltas = []
    first_masked: tuple[str, str] | None = None
    first_scores: np.ndarray | None = None
    for row in rows:
        fixture_id = row["fixture_id"]
        sequences = {}
        expected_hashes = {
            "REF": row["reference_sequence_sha256"],
            "ALT": row["alternative_sequence_sha256"],
            "REF_RC": row["reverse_complement_reference_sha256"],
            "ALT_RC": row["reverse_complement_alternative_sha256"],
        }
        allele_features = {}
        for allele in ("REF", "ALT", "REF_RC", "ALT_RC"):
            record_id = f"{fixture_id}|common_6000|{allele}"
            sequence = fasta.get(record_id)
            if sequence is None or _sha256_text(sequence) != expected_hashes[allele]:
                raise NucleotideTransformerProbeError("NT FASTA identity differs")
            sequences[allele] = sequence
            value = _unmasked_features(
                model,
                tokenizer,
                sequence,
                int(row["pool_start0"]),
                int(row["pool_end0"]),
                device,
            )
            allele_features[allele] = value
            for phase in range(6):
                feature_ids.append(f"{record_id}|phase{phase}")
                features.append(value[phase])
        forward = _masked_phase_scores(
            model, tokenizer, sequences["REF"], sequences["ALT"], device
        )
        reverse = _masked_phase_scores(
            model, tokenizer, sequences["REF_RC"], sequences["ALT_RC"], device
        )
        if first_masked is None:
            first_masked = (sequences["REF"], sequences["ALT"])
            first_scores = forward
        record: dict[str, object] = {
            "fixture_id": fixture_id,
            "genomic_fold": int(row["genomic_fold"]),
        }
        for phase in range(6):
            record[f"forward_phase{phase}_alt_minus_ref_log_probability"] = float(
                forward[phase]
            )
            record[f"reverse_phase{phase}_alt_minus_ref_log_probability"] = float(
                reverse[phase]
            )
        record["forward_six_phase_mean"] = float(forward.mean())
        record["reverse_six_phase_mean"] = float(reverse.mean())
        record["rc_averaged_six_phase_mean"] = float(
            (forward.mean() + reverse.mean()) / 2.0
        )
        score_rows.append(record)
        delta = (
            allele_features["ALT"].mean(axis=0)
            - allele_features["REF"].mean(axis=0)
            + allele_features["ALT_RC"].mean(axis=0)
            - allele_features["REF_RC"].mean(axis=0)
        ) / 2.0
        allele_delta_ids.append(fixture_id)
        allele_deltas.append(delta)

    assert first_masked is not None and first_scores is not None
    repeated = _masked_phase_scores(model, tokenizer, *first_masked, device)
    reproducibility_max_abs_diff = float(np.max(np.abs(repeated - first_scores)))
    if reproducibility_max_abs_diff > 1e-6:
        raise NucleotideTransformerProbeError("NT repeated forward differs")

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
        "schema_version": "masld-bench-nucleotide-transformer-runtime-probe-v1",
        "status": "pass",
        "seed": seed,
        "checkpoint_sha256": _sha256_file(checkpoint),
        "model_parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "fixture_count": len(rows),
        "sequence_phase_forward_count": 36 * len(rows) + 6,
        "feature_width": int(features[0].shape[0]),
        "phase_count": 6,
        "model_vocab_size": 4105,
        "tokenizer_vocab_size": 4107,
        "tokenizer_extra_ids_excluded_from_inputs": [4105, 4106],
        "source_tokenizer_eos_token": tokenizer_settings["eos_token"],
        "pool_length_bp": 1536,
        "common_input_length_bp": 6000,
        "effective_per_phase_context_bp": PHASE_CONTEXT_BP,
        "native_score": "six_phase_masked_reconstruction_ALT_minus_REF_with_reverse_complement_average",
        "feature_contract": "phase_specific_1536bp_mean_pool_plus_phase_and_RC_averaged_allele_delta",
        "restricted_comparator": True,
        "legacy_nonpersistent_position_ids_validated": True,
        "open_champion_eligible": False,
        "embeddings_biologically_meaningful_without_trained_head": False,
        "head_fit": False,
        "observed_outcomes_loaded": False,
        "sealed_outcomes_loaded": False,
        "runtime_network_allowed": False,
        "reproducibility_max_abs_diff": reproducibility_max_abs_diff,
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
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--vocab", type=Path, required=True)
    parser.add_argument("--tokenizer-config", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260824)
    args = parser.parse_args()
    run(
        args.config,
        args.vocab,
        args.tokenizer_config,
        args.checkpoint,
        args.fixture,
        args.output,
        args.seed,
    )


if __name__ == "__main__":
    main()
