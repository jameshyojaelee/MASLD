#!/usr/bin/env python3
"""Validate task-native sequence controls on the frozen outcome-blind fixture."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path

from safetensors.torch import save
import torch

from scripts.gse281364_task_native_sequence_controls import (
    MODEL_IDS,
    allele_pair_activity,
    build_model,
    joint_orientation_augmentation,
    reverse_complement_one_hot,
)


SEEDS = (1103, 2909, 4721, 6673, 8111)


class SequenceControlValidationError(RuntimeError):
    """Raised when a model or fixture requirement differs."""


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = ""
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            value = line.strip()
            if value.startswith(">"):
                name = value[1:]
            elif not name or name in records:
                raise SequenceControlValidationError("smoke FASTA differs")
            else:
                records[name] = value
                name = ""
    return records


def one_hot(sequences: list[str]) -> torch.Tensor:
    output = torch.zeros((len(sequences), 4, 230), dtype=torch.float32)
    for row, sequence in enumerate(sequences):
        if len(sequence) != 230:
            raise SequenceControlValidationError("sequence length differs")
        for index, base in enumerate(sequence):
            if base not in "AGCTN":
                raise SequenceControlValidationError("sequence alphabet differs")
            if base == "N":
                output[row, :, index] = 0.25
            else:
                output[row, "AGCT".index(base), index] = 1.0
    return output


def state_hash(model: torch.nn.Module) -> str:
    state = {key: value.detach().cpu().contiguous() for key, value in sorted(model.state_dict().items())}
    return sha256(save(state)).hexdigest()


def validate(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists():
        raise SequenceControlValidationError("validation output exists")
    fixture_receipt = json.loads(
        (arguments.fixture / "fixture/receipt.json").read_text(encoding="utf-8")
    )
    if (
        fixture_receipt.get("status") != "pass_outcome_blind_task_fixture_freeze"
        or fixture_receipt.get("outcomes_read")
        or fixture_receipt.get("metrics_calculated")
        or fixture_receipt.get("model_training_executed")
    ):
        raise SequenceControlValidationError("fixture receipt differs")
    with (arguments.fixture / "fixture/smoke_manifest.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        smoke = list(csv.DictReader(handle, delimiter="\t"))
    records = read_fasta(arguments.fixture / "fixture/sequence_controls.alleles.fa.gz")
    if len(smoke) != 10 or len(records) != 4132:
        raise SequenceControlValidationError("fixture census differs")
    reference = one_hot([records[f"{row['element_id']}|REF"] for row in smoke])
    alternative = one_hot([records[f"{row['element_id']}|ALT"] for row in smoke])
    if not torch.equal(
        reverse_complement_one_hot(reverse_complement_one_hot(reference)), reference
    ):
        raise SequenceControlValidationError("reverse-complement involution failed")

    rows = []
    parameter_counts: dict[str, int] = {}
    for model_id in MODEL_IDS:
        hashes = []
        for seed in SEEDS:
            model = build_model(model_id, seed=seed)
            model.train()
            reverse_mask = torch.tensor(
                [(index + seed) % 2 == 0 for index in range(len(smoke))],
                dtype=torch.bool,
            )
            augmented_ref, augmented_alt = joint_orientation_augmentation(
                reference, alternative, reverse_mask
            )
            predictions = allele_pair_activity(
                model,
                augmented_ref,
                augmented_alt,
                reverse_complement_average=False,
            )
            target = torch.zeros_like(predictions["alt_minus_ref"])
            loss = torch.nn.functional.smooth_l1_loss(
                predictions["alt_minus_ref"], target, beta=1.0
            )
            loss.backward()
            if not torch.isfinite(loss) or not all(
                parameter.grad is None or torch.isfinite(parameter.grad).all()
                for parameter in model.parameters()
            ):
                raise SequenceControlValidationError("synthetic backward differs")
            model.eval()
            with torch.inference_mode():
                forward = allele_pair_activity(
                    model,
                    reference,
                    alternative,
                    reverse_complement_average=True,
                )
                swapped = allele_pair_activity(
                    model,
                    alternative,
                    reference,
                    reverse_complement_average=True,
                )
            if (
                forward["alt_minus_ref"].shape != (10, 2)
                or not torch.isfinite(forward["alt_minus_ref"]).all()
                or not torch.equal(forward["alt_minus_ref"], -swapped["alt_minus_ref"])
            ):
                raise SequenceControlValidationError("allele sign or output differs")
            digest = state_hash(model)
            hashes.append(digest)
            rows.append(
                {
                    "model_id": model_id,
                    "seed": seed,
                    "state_sha256_after_synthetic_backward": digest,
                    "synthetic_loss_finite": "true",
                    "ref_alt_swap_negates_delta": "true",
                    "metrics_calculated": "false",
                }
            )
            if seed == SEEDS[0]:
                parameter_counts[model_id] = sum(
                    parameter.numel() for parameter in model.parameters()
                )
        if len(set(hashes)) != 5:
            raise SequenceControlValidationError("five seeds are not distinct states")
        first = build_model(model_id, seed=SEEDS[0])
        second = build_model(model_id, seed=SEEDS[0])
        if state_hash(first) != state_hash(second):
            raise SequenceControlValidationError("same-seed initialization differs")
    ratio = max(parameter_counts.values()) / min(parameter_counts.values())
    if ratio > 1.35 or min(parameter_counts.values()) < 200_000 or max(parameter_counts.values()) > 400_000:
        raise SequenceControlValidationError("parameter-count parity or budget differs")

    source = arguments.architecture_source.read_text(encoding="utf-8").lower()
    if any(token in source for token in ("chrombpnet", "scbasset", "tensorflow", ".h5")):
        raise SequenceControlValidationError("local-ATAC component reuse detected")
    arguments.output.mkdir(mode=0o750)
    with (arguments.output / "seed_state_fixture.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    receipt = {
        "schema_version": "masld-bench-gse281364-sequence-control-validation-v1",
        "status": "pass_outcome_blind_architecture_fixture",
        "models": list(MODEL_IDS),
        "parameter_counts": parameter_counts,
        "parameter_count_ratio": ratio,
        "input_shape": [10, 4, 230],
        "output_shape_per_activity_call": [10, 2],
        "five_genuine_seed_states_per_model": True,
        "same_seed_initialization_deterministic": True,
        "shared_ref_alt_encoder": True,
        "reverse_complement_involution": True,
        "ref_alt_swap_negates_delta": True,
        "synthetic_forward_backward": True,
        "local_atac_code_reused": False,
        "local_atac_weights_read": False,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "metrics_calculated": False,
        "model_training_executed": False,
        "predictions_generated": False,
    }
    (arguments.output / "architecture_contract.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--architecture-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    validate(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
