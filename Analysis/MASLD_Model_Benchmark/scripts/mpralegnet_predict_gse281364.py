#!/usr/bin/env python3
"""Run MPRALegNet HepG2 on outcome-blind GSE281364 230-bp alleles."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
import platform

import numpy as np
from safetensors.torch import load_file
import torch
from torch import nn

from scripts.mpralegnet_probe_native import (
    COMPLEMENT,
    _load_module,
    _one_hot,
    _read_fasta as read_probe_fasta,
    _sha256_file,
)


EXPECTED_CHECKPOINT_SHA256 = "470dc7bfd3f0912c91f307a5f4019598072b8786294c46fc501b826c030cbe39"
EXPECTED_ELEMENTS = 4_359
EXPECTED_GROUPS = 1_033
EXPECTED_SEQUENCES = EXPECTED_ELEMENTS * 4
ALLELES = ("REF", "ALT", "REF_RC", "ALT_RC")
CPU_GPU_PARITY_ATOL = 1.0e-4


class MPRALegNetPredictionError(ValueError):
    """Raised when MPRALegNet execution or input provenance differs."""


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    pieces: list[str] = []
    with gzip.open(path, "rt", encoding="ascii", errors="strict") as handle:
        for line in handle:
            value = line.strip()
            if value.startswith(">"):
                if name is not None:
                    if name in records:
                        raise MPRALegNetPredictionError("duplicate FASTA identifier")
                    records[name] = "".join(pieces)
                name, pieces = value[1:], []
            elif name is None:
                raise MPRALegNetPredictionError("FASTA sequence precedes identifier")
            else:
                pieces.append(value.upper())
    if name is not None:
        if name in records:
            raise MPRALegNetPredictionError("duplicate FASTA identifier")
        records[name] = "".join(pieces)
    return records


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    if (
        len(rows) != EXPECTED_ELEMENTS
        or len({row["fixture_id"] for row in rows}) != EXPECTED_ELEMENTS
        or len({row["element_id"] for row in rows}) != EXPECTED_ELEMENTS
        or len({row["outer_locus_sequence_group_id"] for row in rows}) != EXPECTED_GROUPS
        or {int(row["outer_fold"]) for row in rows} != set(range(5))
        or any(row["cell_context"] != "HepG2" for row in rows)
        or any(row["input_lane"] != "GRCh38p14_230bp_genomic_window_transfer" for row in rows)
    ):
        raise MPRALegNetPredictionError("fixture manifest census or context differs")
    return rows


def validate_fixture_sequences(
    rows: list[dict[str, str]], records: dict[str, str]
) -> list[str]:
    ordered_names: list[str] = []
    expected_names: set[str] = set()
    for row in rows:
        fixture_id = row["fixture_id"]
        names = [f"{fixture_id}|{allele}" for allele in ALLELES]
        expected_names.update(names)
        forward_ref, forward_alt, reverse_ref, reverse_alt = (
            records.get(name, "") for name in names
        )
        if (
            any(len(sequence) != 230 for sequence in (forward_ref, forward_alt, reverse_ref, reverse_alt))
            or set("".join((forward_ref, forward_alt, reverse_ref, reverse_alt))) - set("ACGT")
            or forward_ref[115] != row["ref"]
            or forward_alt[115] != row["alt"]
            or reverse_ref != forward_ref.translate(COMPLEMENT)[::-1]
            or reverse_alt != forward_alt.translate(COMPLEMENT)[::-1]
            or reverse_ref[114] != row["ref"].translate(COMPLEMENT)
            or reverse_alt[114] != row["alt"].translate(COMPLEMENT)
            or _sha256_file_bytes(forward_ref) != row["reference_sequence_sha256"]
            or _sha256_file_bytes(forward_alt) != row["alternative_sequence_sha256"]
            or _sha256_file_bytes(reverse_ref) != row["reference_reverse_complement_sha256"]
            or _sha256_file_bytes(reverse_alt) != row["alternative_reverse_complement_sha256"]
        ):
            raise MPRALegNetPredictionError("fixture allele sequence or orientation differs")
        ordered_names.extend(names)
    if len(records) != EXPECTED_SEQUENCES or set(records) != expected_names:
        raise MPRALegNetPredictionError("fixture FASTA census differs")
    return ordered_names


def _sha256_file_bytes(value: str) -> str:
    from hashlib import sha256

    return sha256(value.encode()).hexdigest()


def load_model(model_source: Path, config_path: Path, checkpoint: Path) -> nn.Module:
    if _sha256_file(checkpoint) != EXPECTED_CHECKPOINT_SHA256:
        raise MPRALegNetPredictionError("checkpoint checksum differs")
    config = json.loads(config_path.read_text())
    if config.get("activation") != "SiLU" or config.get("in_ch") != 4:
        raise MPRALegNetPredictionError("model config differs")
    module = _load_module(model_source)
    constructor = {
        key: value for key, value in config.items() if key not in {"model_type", "activation"}
    }
    model = module.LegNet(**constructor, activation=nn.SiLU)
    state = load_file(str(checkpoint), device="cpu")
    if len(state) != 134 or sum(tensor.numel() for tensor in state.values()) != 1_330_548:
        raise MPRALegNetPredictionError("checkpoint tensor census differs")
    if any(not key.startswith("model.") for key in state):
        raise MPRALegNetPredictionError("checkpoint state prefix differs")
    normalized = {key.removeprefix("model."): value for key, value in state.items()}
    incompatible = model.load_state_dict(normalized, strict=True)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        raise MPRALegNetPredictionError("strict checkpoint restore differs")
    model.eval()
    if model.training or any(module.training for module in model.modules()):
        raise MPRALegNetPredictionError("BatchNorm eval-mode contract differs")
    return model


def score_sequences(
    model: nn.Module,
    records: dict[str, str],
    ordered_names: list[str],
    batch_size: int,
) -> tuple[np.ndarray, float]:
    try:
        device = next(model.parameters()).device
    except StopIteration as error:
        raise MPRALegNetPredictionError("model has no parameters") from error
    output = np.empty(len(ordered_names), dtype=np.float32)
    repeat_max_abs = 0.0
    with torch.inference_mode():
        for start in range(0, len(ordered_names), batch_size):
            end = min(start + batch_size, len(ordered_names))
            inputs = torch.from_numpy(
                np.stack([_one_hot(records[name]) for name in ordered_names[start:end]])
            ).to(device=device)
            scores = model(inputs)
            if scores.shape != (end - start,) or not torch.isfinite(scores).all():
                raise MPRALegNetPredictionError("prediction shape or finiteness differs")
            if start == 0:
                repeat = model(inputs)
                repeat_max_abs = float((scores - repeat).abs().max())
                if repeat_max_abs != 0.0:
                    raise MPRALegNetPredictionError("deterministic repeat differs")
            output[start:end] = scores.detach().cpu().numpy()
    return output, repeat_max_abs


def aggregate_allele_scores(values: np.ndarray) -> tuple[float, float, float]:
    if values.shape != (4,) or not np.isfinite(values).all():
        raise MPRALegNetPredictionError("allele score vector differs")
    reference = float((values[0] + values[2]) / 2.0)
    alternative = float((values[1] + values[3]) / 2.0)
    return reference, alternative, alternative - reference


def cpu_gpu_parity(
    model: nn.Module,
    probe_manifest: Path,
    probe_fasta: Path,
    probe_scores: Path,
) -> float:
    records_6000 = read_probe_fasta(probe_fasta)
    with probe_manifest.open(encoding="utf-8", newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["model_id"] == "hyenadna" and row["context_id"] == "common_6000"
        ]
    sequences: dict[str, str] = {}
    names: list[str] = []
    for row in sorted(rows, key=lambda value: int(value["genomic_fold"])):
        fixture_id = row["fixture_id"]
        reference = records_6000[f"{fixture_id}|common_6000|REF"][2885:3115]
        alternative = records_6000[f"{fixture_id}|common_6000|ALT"][2885:3115]
        for allele, sequence in (
            ("REF", reference),
            ("ALT", alternative),
            ("REF_RC", reference.translate(COMPLEMENT)[::-1]),
            ("ALT_RC", alternative.translate(COMPLEMENT)[::-1]),
        ):
            name = f"{fixture_id}|{allele}"
            names.append(name)
            sequences[name] = sequence
    observed, _ = score_sequences(model, sequences, names, batch_size=len(names))
    with np.load(probe_scores, allow_pickle=False) as source:
        expected_names = source["names"].astype(str).tolist()
        expected = source["scores"].astype(np.float32)
    if names != expected_names or expected.shape != observed.shape:
        raise MPRALegNetPredictionError("runtime parity fixture differs")
    maximum = float(np.max(np.abs(observed - expected)))
    return maximum


def predict(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists():
        raise MPRALegNetPredictionError("prediction output exists")
    if not 1 <= arguments.batch_size <= 1024:
        raise MPRALegNetPredictionError("batch size differs")
    if not 1 <= arguments.threads <= 16:
        raise MPRALegNetPredictionError("thread count differs")
    torch.manual_seed(20260824)
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(arguments.threads)
    if arguments.device == "cuda":
        if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
            raise MPRALegNetPredictionError("exactly one visible CUDA device is required")
        gpu_name = torch.cuda.get_device_name(0)
        if "L40S" not in gpu_name:
            raise MPRALegNetPredictionError("visible CUDA device is not an L40S")
        torch.cuda.manual_seed_all(20260824)
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
        torch.cuda.reset_peak_memory_stats(0)
        device = torch.device("cuda:0")
    else:
        gpu_name = None
        device = torch.device("cpu")
    model = load_model(arguments.model_source, arguments.config, arguments.checkpoint)
    model.to(device)
    running_before = {
        name: value.detach().clone()
        for name, value in model.named_buffers()
        if "running_" in name or "num_batches_tracked" in name
    }
    parity_max_abs = cpu_gpu_parity(
        model, arguments.probe_manifest, arguments.probe_fasta, arguments.probe_scores
    )
    if parity_max_abs > CPU_GPU_PARITY_ATOL:
        raise MPRALegNetPredictionError(
            "runtime-to-frozen-L40S parity tolerance failed: "
            f"{parity_max_abs:.9g} > {CPU_GPU_PARITY_ATOL:.9g}"
        )
    rows = read_manifest(arguments.manifest)
    records = read_fasta(arguments.fasta)
    ordered_names = validate_fixture_sequences(rows, records)
    scores, repeat_max_abs = score_sequences(
        model, records, ordered_names, arguments.batch_size
    )
    running_after = {
        name: value.detach().clone()
        for name, value in model.named_buffers()
        if "running_" in name or "num_batches_tracked" in name
    }
    if running_before.keys() != running_after.keys() or any(
        not torch.equal(running_before[name], running_after[name]) for name in running_before
    ):
        raise MPRALegNetPredictionError("BatchNorm running state changed")
    if arguments.device == "cuda":
        torch.cuda.synchronize(0)
        peak_gpu_memory_bytes = int(torch.cuda.max_memory_allocated(0))
    else:
        peak_gpu_memory_bytes = 0
    arguments.output.mkdir(parents=True)
    np.savez_compressed(
        arguments.output / "native_scores.npz",
        names=np.asarray(ordered_names),
        scores=scores,
    )
    fields = (
        "element_id",
        "outer_locus_sequence_group_id",
        "outer_fold",
        "contig",
        "variant_pos0",
        "ref",
        "alt",
        "cell_context",
        "input_lane",
        "reference_forward_score",
        "alternative_forward_score",
        "reference_reverse_complement_score",
        "alternative_reverse_complement_score",
        "reference_orientation_mean_score",
        "alternative_orientation_mean_score",
        "alt_minus_ref_score",
    )
    with (arguments.output / "allele_scores.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, row in enumerate(rows):
            values = scores[index * 4 : index * 4 + 4]
            reference, alternative, delta = aggregate_allele_scores(values)
            writer.writerow(
                {
                    **{field: row[field] for field in fields[:9]},
                    "reference_forward_score": float(values[0]),
                    "alternative_forward_score": float(values[1]),
                    "reference_reverse_complement_score": float(values[2]),
                    "alternative_reverse_complement_score": float(values[3]),
                    "reference_orientation_mean_score": reference,
                    "alternative_orientation_mean_score": alternative,
                    "alt_minus_ref_score": delta,
                }
            )
    fold_counts = {
        str(fold): sum(int(row["outer_fold"]) == fold for row in rows) for fold in range(5)
    }
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-mpralegnet-gse281364-prediction-v1",
        "status": "pass_outcome_blind_prediction",
        "dataset_id": "gse281364",
        "model_id": "mpralegnet_hepg2_test1_val2",
        "checkpoint_sha256": EXPECTED_CHECKPOINT_SHA256,
        "checkpoint_provenance": "third_party_Hugging_Face_repackaging_author_equivalence_unresolved",
        "elements": len(rows),
        "outer_locus_sequence_groups": len(
            {row["outer_locus_sequence_group_id"] for row in rows}
        ),
        "fold_elements": fold_counts,
        "input_length_bp": 230,
        "forward_variant_index0": 115,
        "reverse_complement_variant_index0": 114,
        "cell_context": "HepG2",
        "input_lane": "GRCh38p14_230bp_genomic_window_transfer",
        "native_output": "uncalibrated_lentiMPRA_reporter_expression_score",
        "allele_effect_sign": "ALT_minus_REF",
        "strand_policy": "mean_forward_and_reverse_complement",
        "device": arguments.device,
        "cpu": platform.processor() or "unreported",
        "gpu": gpu_name,
        "peak_gpu_memory_bytes": peak_gpu_memory_bytes,
        "threads": arguments.threads,
        "batch_size": arguments.batch_size,
        "torch": torch.__version__,
        "deterministic_repeat_max_abs_diff": repeat_max_abs,
        "cpu_to_l40s_parity_max_abs_diff": parity_max_abs,
        "runtime_to_frozen_l40s_parity_max_abs_diff": parity_max_abs,
        "cpu_to_l40s_parity_tolerance": CPU_GPU_PARITY_ATOL,
        "batchnorm_eval_mode": True,
        "batchnorm_running_state_unchanged": True,
        "fixture_artifacts_sha256": _sha256_file(arguments.fixture_artifacts),
        "runtime_probe_artifacts_sha256": _sha256_file(arguments.runtime_probe_artifacts),
        "native_scores_sha256": _sha256_file(arguments.output / "native_scores.npz"),
        "allele_scores_sha256": _sha256_file(arguments.output / "allele_scores.tsv"),
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_outcomes_read": False,
        "model_fitted_or_adapted": False,
        "standalone_champion_eligible": False,
        "terminal_disposition": "secondary_HepG2_reporter_activity_comparator_only",
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-source", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--fixture-artifacts", type=Path, required=True)
    parser.add_argument("--probe-manifest", type=Path, required=True)
    parser.add_argument("--probe-fasta", type=Path, required=True)
    parser.add_argument("--probe-scores", type=Path, required=True)
    parser.add_argument("--runtime-probe-artifacts", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--output", type=Path, required=True)
    predict(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
