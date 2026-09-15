#!/usr/bin/env python3
"""Run the exact HyenaDNA checkpoint on the outcome-free common fixture."""

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


TOKEN_IDS = {"A": 7, "C": 8, "G": 9, "T": 10, "N": 11}


class HyenaProbeError(ValueError):
    """Raised when the runtime probe does not meet its frozen fixture requirements."""


def _sha256_text(value: str) -> str:
    return sha256(value.encode("ascii")).hexdigest()


def _sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
                if not name or name in records:
                    raise HyenaProbeError("FASTA header differs")
                records[name] = []
            elif name is None:
                raise HyenaProbeError("FASTA sequence precedes header")
            else:
                records[name].append(line)
    values = {key: "".join(parts) for key, parts in records.items()}
    if any(set(value) - set(TOKEN_IDS) for value in values.values()):
        raise HyenaProbeError("FASTA contains unsupported bases")
    return values


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["model_id"] == "hyenadna" and row["context_id"] == "common_6000"
        ]
    if len(rows) != 3 or len({row["genomic_fold"] for row in rows}) != 3:
        raise HyenaProbeError("HyenaDNA common fixture census differs")
    return rows


def _tokenize(sequence: str) -> torch.Tensor:
    # This exactly reproduces the included character tokenizer: one token per
    # input base followed by [SEP]=1. No padding or truncation is permitted.
    return torch.tensor([TOKEN_IDS[base] for base in sequence] + [1], dtype=torch.long)


def _predict(
    model: torch.nn.Module,
    sequence: str,
    mutation_index0: int,
    pool_start0: int,
    pool_end0: int,
    device: torch.device,
) -> tuple[float, np.ndarray]:
    ids = _tokenize(sequence).unsqueeze(0).to(device)
    with torch.inference_mode():
        hidden = model.hyena(input_ids=ids, return_dict=True).last_hidden_state
        logits = model.lm_head(hidden).float()
        target_positions = torch.arange(mutation_index0, len(sequence), device=device)
        predictor_positions = target_positions - 1
        selected_logits = logits[0, predictor_positions]
        selected_targets = ids[0, target_positions]
        log_prob = torch.log_softmax(selected_logits, dim=-1)
        mean_suffix_log_likelihood = log_prob.gather(
            1, selected_targets.unsqueeze(1)
        ).mean()
        pooled = hidden[0, pool_start0:pool_end0].float().mean(dim=0)
    return float(mean_suffix_log_likelihood.cpu()), pooled.cpu().numpy()


def run(
    code_root: Path,
    checkpoint: Path,
    fixture: Path,
    output: Path,
    seed: int,
) -> dict[str, object]:
    if output.exists() or code_root.is_symlink() or checkpoint.is_symlink() or fixture.is_symlink():
        raise HyenaProbeError("HyenaDNA runtime request differs")
    output.mkdir(mode=0o750)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.use_deterministic_algorithms(True)

    sys.path.insert(0, str(code_root.parent))
    from hyenadna.configuration_hyena import HyenaConfig
    from hyenadna.modeling_hyena import HyenaDNAForCausalLM

    config_path = code_root / "config.json"
    config = HyenaConfig.from_json_file(str(config_path))
    if config.model_type != "hyenadna" or config.max_seq_len != 1_000_002:
        raise HyenaProbeError("HyenaDNA config differs")
    model = HyenaDNAForCausalLM(config)
    state = load_file(str(checkpoint), device="cpu")
    incompatible = model.load_state_dict(state, strict=False)
    tied_alias_keys = {
        f"hyena.backbone.layers.{layer}.mixer.filter_fn.implicit_filter.{position}.freq"
        for layer in range(config.n_layer)
        for position in (3, 5)
    }
    tied_modules_are_identical = all(
        model.hyena.backbone.layers[layer].mixer.filter_fn.implicit_filter[1]
        is model.hyena.backbone.layers[layer].mixer.filter_fn.implicit_filter[3]
        is model.hyena.backbone.layers[layer].mixer.filter_fn.implicit_filter[5]
        for layer in range(config.n_layer)
    )
    if (
        set(incompatible.missing_keys) != tied_alias_keys
        or incompatible.unexpected_keys
        or not tied_modules_are_identical
    ):
        raise HyenaProbeError("HyenaDNA checkpoint keys differ")
    del state
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise HyenaProbeError("HyenaDNA production probe requires CUDA")
    model.to(device).eval()

    fasta = _read_fasta(fixture / "fixture" / "common_6000.alleles.fa.gz")
    rows = _read_manifest(fixture / "fixture" / "sequence_manifest.tsv")
    score_rows: list[dict[str, object]] = []
    feature_ids: list[str] = []
    features: list[np.ndarray] = []
    allele_delta_ids: list[str] = []
    allele_deltas: list[np.ndarray] = []
    first_reproducibility_input: tuple[str, int, int, int] | None = None
    first_reproducibility_output: tuple[float, np.ndarray] | None = None

    for row in rows:
        fixture_id = row["fixture_id"]
        mutation_indices = {
            "REF": int(row["forward_variant_index0"]),
            "ALT": int(row["forward_variant_index0"]),
            "REF_RC": int(row["reverse_complement_variant_index0"]),
            "ALT_RC": int(row["reverse_complement_variant_index0"]),
        }
        expected_hashes = {
            "REF": row["reference_sequence_sha256"],
            "ALT": row["alternative_sequence_sha256"],
            "REF_RC": row["reverse_complement_reference_sha256"],
            "ALT_RC": row["reverse_complement_alternative_sha256"],
        }
        observed: dict[str, tuple[float, np.ndarray]] = {}
        for allele in ("REF", "ALT", "REF_RC", "ALT_RC"):
            record_id = f"{fixture_id}|common_6000|{allele}"
            sequence = fasta.get(record_id)
            if sequence is None or len(sequence) != 6000:
                raise HyenaProbeError("HyenaDNA FASTA record differs")
            if _sha256_text(sequence) != expected_hashes[allele]:
                raise HyenaProbeError("HyenaDNA FASTA hash differs")
            prediction = _predict(
                model,
                sequence,
                mutation_indices[allele],
                int(row["pool_start0"]),
                int(row["pool_end0"]),
                device,
            )
            observed[allele] = prediction
            feature_ids.append(record_id)
            features.append(prediction[1])
            if first_reproducibility_input is None:
                first_reproducibility_input = (
                    sequence,
                    mutation_indices[allele],
                    int(row["pool_start0"]),
                    int(row["pool_end0"]),
                )
                first_reproducibility_output = prediction
        forward_delta = observed["ALT"][0] - observed["REF"][0]
        reverse_delta = observed["ALT_RC"][0] - observed["REF_RC"][0]
        averaged_delta = (forward_delta + reverse_delta) / 2.0
        feature_delta = (
            observed["ALT"][1]
            - observed["REF"][1]
            + observed["ALT_RC"][1]
            - observed["REF_RC"][1]
        ) / 2.0
        score_rows.append(
            {
                "fixture_id": fixture_id,
                "genomic_fold": int(row["genomic_fold"]),
                "forward_suffix_bases": 6000 - int(row["forward_variant_index0"]),
                "reverse_suffix_bases": 6000 - int(row["reverse_complement_variant_index0"]),
                "forward_alt_minus_ref_mean_log_likelihood": forward_delta,
                "reverse_complement_alt_minus_ref_mean_log_likelihood": reverse_delta,
                "rc_averaged_alt_minus_ref_mean_log_likelihood": averaged_delta,
            }
        )
        allele_delta_ids.append(fixture_id)
        allele_deltas.append(feature_delta)

    assert first_reproducibility_input is not None
    assert first_reproducibility_output is not None
    repeated = _predict(model, *first_reproducibility_input, device)
    reproducibility_score_abs_diff = abs(repeated[0] - first_reproducibility_output[0])
    reproducibility_feature_max_abs_diff = float(
        np.max(np.abs(repeated[1] - first_reproducibility_output[1]))
    )
    if reproducibility_score_abs_diff > 1e-7 or reproducibility_feature_max_abs_diff > 1e-6:
        raise HyenaProbeError("HyenaDNA repeated forward differs")

    score_path = output / "native_scores.tsv"
    with score_path.open("w", newline="", encoding="utf-8") as handle:
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
        "schema_version": "masld-bench-hyenadna-runtime-probe-v1",
        "status": "pass",
        "seed": seed,
        "checkpoint_sha256": _sha256_file(checkpoint),
        "model_parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "fixture_count": len(rows),
        "sequence_forward_count": len(feature_ids) + 1,
        "feature_width": int(features[0].shape[0]),
        "pool_length_bp": 1536,
        "common_input_length_bp": 6000,
        "native_score": "length_normalized_affected_suffix_autoregressive_ALT_minus_REF_with_reverse_complement_average",
        "feature_contract": "mean_pool_1536bp_REF_ALT_and_RC_plus_RC_averaged_allele_delta",
        "embeddings_biologically_meaningful_without_trained_head": False,
        "head_fit": False,
        "observed_outcomes_loaded": False,
        "sealed_outcomes_loaded": False,
        "runtime_network_allowed": False,
        "checkpoint_loaded_with_safetensors": True,
        "safetensors_omitted_shared_parameter_aliases": sorted(tied_alias_keys),
        "shared_activation_modules_identity_validated": tied_modules_are_identical,
        "custom_code_execution": "exact_AST_audited_revision_pinned_bundle",
        "reproducibility_score_abs_diff": reproducibility_score_abs_diff,
        "reproducibility_feature_max_abs_diff": reproducibility_feature_max_abs_diff,
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
