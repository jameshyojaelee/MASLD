#!/usr/bin/env python3
"""Freeze tokenizer and coordinate requirements for the outcome-blind GSE281364 fixture."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from transformers import EsmTokenizer, PreTrainedTokenizerFast

from scripts.gse281364_dna_lm_native_contract import (
    ALLELES,
    INPUT_LENGTH,
    MODELS,
    NativeContractError,
    affected_interval,
    dnabert_encoding,
    load_config,
    nt_pool_indices,
    validate_nt_pair,
    overlapping_indices,
    read_fixture,
)


FIELDS = (
    "fixture_id",
    "outer_locus_sequence_group_id",
    "outer_fold",
    "dnabert_ref_tokens",
    "dnabert_alt_tokens",
    "dnabert_ref_rc_tokens",
    "dnabert_alt_rc_tokens",
    "dnabert_forward_affected_start0",
    "dnabert_forward_affected_end0",
    "dnabert_reverse_affected_start0",
    "dnabert_reverse_affected_end0",
    "dnabert_min_pool_tokens",
    "dnabert_max_pool_tokens",
    "nt_phases",
    "nt_context_bp_per_phase",
    "nt_tokens_per_phase_including_cls",
    "nt_total_edge_trim_bp_per_phase",
    "nt_min_pool_tokens",
    "nt_max_pool_tokens",
    "hyenadna_tokens_including_sep",
)


def validate_runtime_receipts(root: Path, config: dict[str, object]) -> dict[str, object]:
    expected = {
        "dnabert2": (768, False, True),
        "nucleotide_transformer": (1280, True, False),
        "hyenadna": (256, False, True),
    }
    bindings: dict[str, object] = {}
    for model, (width, restricted, champion) in expected.items():
        record = config["models"][model]
        probe = root / record["successful_probe_path"]
        receipt = json.loads(
            (probe / "predictions/runtime_probe_receipt.json").read_text(encoding="utf-8")
        )
        eligible = receipt.get(
            "open_champion_eligible_after_task_gates",
            receipt.get("open_champion_eligible", True),
        )
        if (
            receipt.get("status") != "pass"
            or receipt.get("feature_width") != width
            or receipt.get("observed_outcomes_loaded") is not False
            or receipt.get("sealed_outcomes_loaded") is not False
            or receipt.get("head_fit") is not False
            or bool(record["restricted_comparator"]) != restricted
            or bool(record["open_champion_eligible_after_task_gates"]) != champion
            or bool(eligible) != champion
        ):
            raise NativeContractError(f"{model} successful runtime receipt differs")
        bindings[model] = {
            "hidden_width": width,
            "checkpoint_sha256": record["checkpoint_sha256"],
            "successful_probe_artifacts_sha256": record[
                "successful_probe_artifacts_sha256"
            ],
            "license": record["license"],
            "restricted_comparator": restricted,
            "open_champion_eligible_after_task_gates": champion,
        }
    return bindings


def preflight(
    *,
    root: Path,
    config_path: Path,
    output: Path,
) -> dict[str, object]:
    if output.exists() or config_path.is_symlink():
        raise NativeContractError("native preflight request differs")
    config = load_config(config_path)
    fixture = root / config["fixture"]["path"]
    rows, fasta = read_fixture(fixture)
    dnabert_bundle = root / config["models"]["dnabert2"]["runtime_bundle_path"]
    nt_source = root / config["models"]["nucleotide_transformer"]["source_path"]
    dnabert = PreTrainedTokenizerFast(
        tokenizer_file=str(dnabert_bundle / "bundle/tokenizer.json"),
        unk_token="[UNK]",
        cls_token="[CLS]",
        sep_token="[SEP]",
        pad_token="[PAD]",
        mask_token="[MASK]",
    )
    nt = EsmTokenizer(str(nt_source / "sources/nt_hf_vocab"), eos_token=None)
    if (
        dnabert.vocab_size != 4_096
        or dnabert.mask_token_id != 4
        or nt.vocab_size != 4_107
        or nt.mask_token_id != 2
        or nt.cls_token_id != 3
    ):
        raise NativeContractError("revision-pinned tokenizer identity differs")
    output.mkdir(parents=True)
    manifest_rows: list[dict[str, object]] = []
    all_dnabert_counts: list[int] = []
    all_dnabert_pool_counts: list[int] = []
    for row in rows:
        sequences = {
            allele: fasta[f"{row['fasta_record_prefix']}|{allele}"] for allele in ALLELES
        }
        encodings = {
            allele: dnabert_encoding(dnabert, sequence)
            for allele, sequence in sequences.items()
        }
        forward_interval = affected_interval(
            encodings["REF"][1], encodings["ALT"][1], int(row["forward_variant_index0"])
        )
        reverse_interval = affected_interval(
            encodings["REF_RC"][1],
            encodings["ALT_RC"][1],
            int(row["reverse_complement_variant_index0"]),
        )
        if (
            not forward_interval[0] <= int(row["forward_variant_index0"]) < forward_interval[1]
            or not reverse_interval[0]
            <= int(row["reverse_complement_variant_index0"])
            < reverse_interval[1]
        ):
            raise NativeContractError("DNABERT-2 affected interval excludes variant")
        pool_counts = [
            len(
                overlapping_indices(
                    offsets, int(row["pool_start0"]), int(row["pool_end0"])
                )
            )
            for _, offsets in encodings.values()
        ]
        if not all(pool_counts):
            raise NativeContractError("DNABERT-2 pool token census differs")
        validate_nt_pair(
            nt,
            sequences["REF"],
            sequences["ALT"],
            int(row["forward_variant_index0"]),
        )
        validate_nt_pair(
            nt,
            sequences["REF_RC"],
            sequences["ALT_RC"],
            int(row["reverse_complement_variant_index0"]),
        )
        nt_pool_counts = [
            len(nt_pool_indices(phase, int(row["pool_start0"]), int(row["pool_end0"])))
            for phase in range(6)
        ]
        token_counts = {allele: len(value[0]) for allele, value in encodings.items()}
        all_dnabert_counts.extend(token_counts.values())
        all_dnabert_pool_counts.extend(pool_counts)
        manifest_rows.append(
            {
                "fixture_id": row["fixture_id"],
                "outer_locus_sequence_group_id": row["outer_locus_sequence_group_id"],
                "outer_fold": row["outer_fold"],
                "dnabert_ref_tokens": token_counts["REF"],
                "dnabert_alt_tokens": token_counts["ALT"],
                "dnabert_ref_rc_tokens": token_counts["REF_RC"],
                "dnabert_alt_rc_tokens": token_counts["ALT_RC"],
                "dnabert_forward_affected_start0": forward_interval[0],
                "dnabert_forward_affected_end0": forward_interval[1],
                "dnabert_reverse_affected_start0": reverse_interval[0],
                "dnabert_reverse_affected_end0": reverse_interval[1],
                "dnabert_min_pool_tokens": min(pool_counts),
                "dnabert_max_pool_tokens": max(pool_counts),
                "nt_phases": 6,
                "nt_context_bp_per_phase": 4_086,
                "nt_tokens_per_phase_including_cls": 682,
                "nt_total_edge_trim_bp_per_phase": 10,
                "nt_min_pool_tokens": min(nt_pool_counts),
                "nt_max_pool_tokens": max(nt_pool_counts),
                "hyenadna_tokens_including_sep": INPUT_LENGTH + 1,
            }
        )
    with (output / "tokenization_manifest.tsv").open(
        "x", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    bindings = validate_runtime_receipts(root, config)
    contract: dict[str, object] = {
        "schema_version": "masld-bench-gse281364-dna-lm-token-coordinate-contract-v1",
        "status": "pass_outcome_blind_native_preflight",
        "dataset_id": "gse281364",
        "fixture_artifacts_sha256": config["fixture"]["artifacts_sha256"],
        "elements": len(rows),
        "outer_locus_sequence_groups": len(rows),
        "outer_folds": 5,
        "allele_sequences": len(fasta),
        "models": list(MODELS),
        "dnabert2": {
            "token_count_min": min(all_dnabert_counts),
            "token_count_max": max(all_dnabert_counts),
            "pool_token_count_min": min(all_dnabert_pool_counts),
            "pool_token_count_max": max(all_dnabert_pool_counts),
            "base_offsets_tile_complete_input": True,
            "truncation": False,
        },
        "nucleotide_transformer": {
            "phases": 6,
            "context_bp_per_phase": 4_086,
            "tokens_per_phase_including_cls": 682,
            "total_edge_trim_bp_per_phase": 10,
            "exactly_one_allele_differing_6mer_per_phase": True,
            "truncation": False,
        },
        "hyenadna": {
            "tokens_including_terminal_sep": 4_097,
            "base_token_coordinate_identity": True,
            "truncation": False,
        },
        "runtime_bindings": bindings,
        "normalization_projection": config["normalization_projection"],
        "tokenizers_loaded": True,
        "checkpoint_bytes_loaded": False,
        "model_forward_executed": False,
        "reporter_counts_read": False,
        "reporter_outcomes_read": False,
        "sealed_labels_read": False,
        "downstream_head_fit": False,
    }
    (output / "contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(contract, sort_keys=True))
    return contract


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    preflight(root=arguments.root, config_path=arguments.config, output=arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
