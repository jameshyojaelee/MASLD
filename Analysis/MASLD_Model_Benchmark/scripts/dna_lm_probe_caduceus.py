#!/usr/bin/env python3
"""Run the exact Caduceus checkpoint on the outcome-free common fixture."""

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
MASK_TOKEN_ID = 3


class CaduceusProbeError(ValueError):
    """Raised when the exact runtime or fixture contract differs."""


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
                    raise CaduceusProbeError("Caduceus FASTA header differs")
                records[name] = []
            elif name is None:
                raise CaduceusProbeError("Caduceus FASTA sequence precedes header")
            else:
                records[name].append(line)
    values = {key: "".join(parts) for key, parts in records.items()}
    if any(len(value) != 6000 or set(value) - set(TOKEN_IDS) for value in values.values()):
        raise CaduceusProbeError("Caduceus FASTA alphabet or length differs")
    return values


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["model_id"] == "caduceus" and row["context_id"] == "common_6000"
        ]
    if len(rows) != 3 or len({row["genomic_fold"] for row in rows}) != 3:
        raise CaduceusProbeError("Caduceus common fixture census differs")
    return rows


def _encode(sequence: str, device: torch.device) -> torch.Tensor:
    return torch.tensor(
        [TOKEN_IDS[base] for base in sequence], dtype=torch.long, device=device
    )


def _realigned_features(
    model: torch.nn.Module,
    sequences: list[str],
    pool_start0: int,
    pool_end0: int,
    device: torch.device,
) -> np.ndarray:
    ids = torch.stack([_encode(sequence, device) for sequence in sequences])
    with torch.inference_mode():
        hidden = model.caduceus(input_ids=ids, return_dict=True).last_hidden_state
        d_model = model.config.d_model
        if tuple(hidden.shape) != (len(sequences), 6000, 2 * d_model):
            raise CaduceusProbeError("Caduceus RCPS hidden-state shape differs")
        forward = hidden[..., :d_model]
        reverse = torch.flip(hidden[..., d_model:], dims=(-2, -1))
        realigned = (forward + reverse) / 2.0
        pooled = realigned[:, pool_start0:pool_end0].float().mean(dim=1)
    return pooled.cpu().numpy()


def _masked_score(
    model: torch.nn.Module,
    reference: str,
    alternative: str,
    variant_index0: int,
    device: torch.device,
) -> float:
    reference_ids = _encode(reference, device)
    alternative_ids = _encode(alternative, device)
    differing = torch.nonzero(reference_ids != alternative_ids, as_tuple=False).flatten()
    if differing.tolist() != [variant_index0]:
        raise CaduceusProbeError("Caduceus allele difference is not the expected SNV")
    masked_reference = reference_ids.clone()
    masked_alternative = alternative_ids.clone()
    masked_reference[variant_index0] = MASK_TOKEN_ID
    masked_alternative[variant_index0] = MASK_TOKEN_ID
    if not torch.equal(masked_reference, masked_alternative):
        raise CaduceusProbeError("Caduceus masked allele contexts differ")
    with torch.inference_mode():
        logits = model(input_ids=masked_reference.unsqueeze(0), return_dict=True).logits
        log_probabilities = torch.log_softmax(logits[0, variant_index0], dim=-1)
        score = (
            log_probabilities[alternative_ids[variant_index0]]
            - log_probabilities[reference_ids[variant_index0]]
        )
    return float(score.cpu())


def _input_gradient_probe(
    model: torch.nn.Module,
    reference: str,
    alternative: str,
    variant_index0: int,
    device: torch.device,
) -> tuple[bool, float]:
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    reference_ids = _encode(reference, device).unsqueeze(0)
    alternative_ids = _encode(alternative, device).unsqueeze(0)
    masked_ids = reference_ids.clone()
    masked_ids[0, variant_index0] = MASK_TOKEN_ID
    inputs_embeds = model.get_input_embeddings()(masked_ids).detach().requires_grad_(True)
    logits = model(inputs_embeds=inputs_embeds, return_dict=True).logits
    value = (
        logits[0, variant_index0, alternative_ids[0, variant_index0]]
        - logits[0, variant_index0, reference_ids[0, variant_index0]]
    )
    value.backward()
    gradient = inputs_embeds.grad
    finite = gradient is not None and bool(torch.isfinite(gradient).all().item())
    maximum = float(gradient.detach().abs().max().cpu()) if gradient is not None else 0.0
    return finite, maximum


def run(
    code_root: Path,
    checkpoint: Path,
    fixture: Path,
    output: Path,
    seed: int,
) -> dict[str, object]:
    if output.exists() or any(
        path.is_symlink() for path in (code_root, checkpoint, fixture)
    ):
        raise CaduceusProbeError("Caduceus runtime request differs")
    output.mkdir(mode=0o750)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False

    sys.path.insert(0, str(code_root.parent))
    from caduceus.configuration_caduceus import CaduceusConfig
    from caduceus.modeling_caduceus import CaduceusForMaskedLM

    config = CaduceusConfig.from_json_file(str(code_root / "config.json"))
    if (
        config.model_type != "caduceus"
        or config.d_model != 256
        or config.n_layer != 16
        or config.vocab_size != 16
        or not config.rcps
        or not config.bidirectional
        or not config.bidirectional_weight_tie
        or config.bidirectional_strategy != "add"
    ):
        raise CaduceusProbeError("Caduceus checkpoint config differs")
    model = CaduceusForMaskedLM(config)
    state = load_file(str(checkpoint), device="cpu")
    expected_aliases = {
        f"caduceus.backbone.layers.{layer}.mixer.submodule.mamba_rev.{projection}.weight"
        for layer in range(config.n_layer)
        for projection in ("in_proj", "out_proj")
    } | {"lm_head.lm_head.weight"}
    incompatible = model.load_state_dict(state, strict=False)
    if set(incompatible.missing_keys) != expected_aliases or incompatible.unexpected_keys:
        raise CaduceusProbeError(
            f"Caduceus checkpoint keys differ: {incompatible!r}"
        )
    del state
    tied_aliases_identical = all(
        layer.mixer.submodule.mamba_rev.in_proj.weight
        is layer.mixer.submodule.mamba_fwd.in_proj.weight
        and layer.mixer.submodule.mamba_rev.out_proj.weight
        is layer.mixer.submodule.mamba_fwd.out_proj.weight
        for layer in model.caduceus.backbone.layers
    ) and model.lm_head.weight is model.get_input_embeddings().weight
    if not tied_aliases_identical:
        raise CaduceusProbeError("Caduceus tied parameter identities differ")

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise CaduceusProbeError("Caduceus production probe requires CUDA")
    model.to(device).eval()
    fasta = _read_fasta(fixture / "fixture" / "common_6000.alleles.fa.gz")
    rows = _read_manifest(fixture / "fixture" / "sequence_manifest.tsv")
    score_rows: list[dict[str, object]] = []
    feature_ids: list[str] = []
    features: list[np.ndarray] = []
    allele_delta_ids: list[str] = []
    allele_deltas: list[np.ndarray] = []
    first_masked: tuple[str, str, int] | None = None
    first_masked_score: float | None = None
    first_unmasked: tuple[str, int, int] | None = None
    first_unmasked_feature: np.ndarray | None = None
    external_rc_feature_max_abs_differences: list[float] = []

    for row in rows:
        fixture_id = row["fixture_id"]
        expected_hashes = {
            "REF": row["reference_sequence_sha256"],
            "ALT": row["alternative_sequence_sha256"],
            "REF_RC": row["reverse_complement_reference_sha256"],
            "ALT_RC": row["reverse_complement_alternative_sha256"],
        }
        sequences: dict[str, str] = {}
        record_ids: list[str] = []
        for allele in ("REF", "ALT", "REF_RC", "ALT_RC"):
            record_id = f"{fixture_id}|common_6000|{allele}"
            sequence = fasta.get(record_id)
            if sequence is None or _sha256_text(sequence) != expected_hashes[allele]:
                raise CaduceusProbeError("Caduceus FASTA identity differs")
            sequences[allele] = sequence
            record_ids.append(record_id)
        values = _realigned_features(
            model,
            [sequences[allele] for allele in ("REF", "ALT", "REF_RC", "ALT_RC")],
            int(row["pool_start0"]),
            int(row["pool_end0"]),
            device,
        )
        if first_unmasked is None:
            first_unmasked = (
                sequences["REF"],
                int(row["pool_start0"]),
                int(row["pool_end0"]),
            )
            first_unmasked_feature = values[0].copy()
        for record_id, value in zip(record_ids, values):
            feature_ids.append(record_id)
            features.append(value)
        forward_score = _masked_score(
            model,
            sequences["REF"],
            sequences["ALT"],
            int(row["forward_variant_index0"]),
            device,
        )
        reverse_score = _masked_score(
            model,
            sequences["REF_RC"],
            sequences["ALT_RC"],
            int(row["reverse_complement_variant_index0"]),
            device,
        )
        if first_masked is None:
            first_masked = (
                sequences["REF"],
                sequences["ALT"],
                int(row["forward_variant_index0"]),
            )
            first_masked_score = forward_score
        score_rows.append(
            {
                "fixture_id": fixture_id,
                "genomic_fold": int(row["genomic_fold"]),
                "forward_alt_minus_ref_masked_log_probability": forward_score,
                "reverse_complement_alt_minus_ref_masked_log_probability": reverse_score,
                "rc_averaged_alt_minus_ref_masked_log_probability": (
                    forward_score + reverse_score
                )
                / 2.0,
                "rc_score_absolute_difference": abs(forward_score - reverse_score),
            }
        )
        allele_delta_ids.append(fixture_id)
        allele_deltas.append((values[1] - values[0] + values[3] - values[2]) / 2.0)
        external_rc_feature_max_abs_differences.extend(
            [
                float(np.max(np.abs(values[0] - values[2]))),
                float(np.max(np.abs(values[1] - values[3]))),
            ]
        )

    assert first_masked is not None and first_masked_score is not None
    assert first_unmasked is not None and first_unmasked_feature is not None
    repeated_masked = _masked_score(model, *first_masked, device)
    repeated_unmasked = _realigned_features(
        model,
        [first_unmasked[0]],
        first_unmasked[1],
        first_unmasked[2],
        device,
    )[0]
    masked_repeat_abs_difference = abs(repeated_masked - first_masked_score)
    feature_repeat_max_abs_difference = float(
        np.max(np.abs(repeated_unmasked - first_unmasked_feature))
    )
    (output / "numeric_repeat_diagnostics.json").write_text(
        json.dumps(
            {
                "masked_repeat_abs_difference": masked_repeat_abs_difference,
                "feature_repeat_max_abs_difference": feature_repeat_max_abs_difference,
                "masked_tolerance": 1e-5,
                "feature_tolerance": 1e-5,
            },
            sort_keys=True,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    if masked_repeat_abs_difference > 1e-5 or feature_repeat_max_abs_difference > 1e-5:
        raise CaduceusProbeError("Caduceus deterministic repeat differs")
    input_gradient_finite, input_gradient_max_abs = _input_gradient_probe(
        model, *first_masked, device
    )
    if not input_gradient_finite or input_gradient_max_abs <= 0.0:
        raise CaduceusProbeError("Caduceus input-gradient backward differs")

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
        "schema_version": "masld-bench-caduceus-runtime-probe-v1",
        "status": "pass",
        "seed": seed,
        "checkpoint_sha256": _sha256_file(checkpoint),
        "learned_parameter_count": sum(parameter.numel() for parameter in model.parameters()),
        "checkpoint_parameter_and_buffer_elements": 7_725_344,
        "checkpoint_tensor_count": 276,
        "fixture_count": len(rows),
        "common_input_length_bp": 6000,
        "pool_length_bp": 1536,
        "feature_width": int(features[0].shape[0]),
        "native_score": "RCPS_masked_language_model_ALT_minus_REF_with_external_reverse_complement_average",
        "feature_contract": "RCPS_internal_channels_realigned_and_averaged_then_1536bp_mean_pool",
        "tied_safetensors_aliases": sorted(expected_aliases),
        "tied_alias_parameter_identity_validated": tied_aliases_identical,
        "max_external_rc_score_absolute_difference": max(
            float(row["rc_score_absolute_difference"]) for row in score_rows
        ),
        "max_external_rc_feature_absolute_difference": max(
            external_rc_feature_max_abs_differences
        ),
        "deterministic_masked_repeat_abs_difference": masked_repeat_abs_difference,
        "deterministic_feature_repeat_max_abs_difference": feature_repeat_max_abs_difference,
        "input_gradient_backward_finite": input_gradient_finite,
        "input_gradient_max_abs": input_gradient_max_abs,
        "compiled_mamba_forward_executed": True,
        "checkpoint_loaded_with_safetensors": True,
        "embeddings_biologically_meaningful_without_trained_head": False,
        "head_fit": False,
        "observed_outcomes_loaded": False,
        "sealed_outcomes_loaded": False,
        "runtime_network_allowed": False,
        "model_and_code_terms": "apache-2.0_open_comparator_and_task_candidate",
        "champion_eligible_after_task_and_external_evaluation_gates": True,
        "device": str(device),
        "torch": torch.__version__,
        "torch_cuda_build": torch.version.cuda,
        "gpu": torch.cuda.get_device_name(0),
        "maximum_gpu_memory_bytes": torch.cuda.max_memory_allocated(),
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
    arguments = parser.parse_args()
    run(
        arguments.code_root,
        arguments.checkpoint,
        arguments.fixture,
        arguments.output,
        arguments.seed,
    )


if __name__ == "__main__":
    main()
