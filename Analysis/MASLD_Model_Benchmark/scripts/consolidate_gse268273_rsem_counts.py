#!/usr/bin/env python3
"""Consolidate outcome-free GSE268273 participant RSEM count vectors."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


class GSE268273ConsolidationError(RuntimeError):
    """Raised when the participant-level quantification matrix is incomplete."""


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE268273ConsolidationError(f"TSV lacks a header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def consolidate(
    *,
    plan_root: Path,
    quant_state_root: Path,
    model_input_root: Path,
    reference_root: Path,
    output: Path,
    bundle_count: int,
) -> dict[str, Any]:
    if output.exists():
        raise GSE268273ConsolidationError(
            f"refusing to overwrite consolidated model input: {output}"
        )
    verify_frozen_tree(plan_root)
    verify_frozen_tree(model_input_root)
    verify_frozen_tree(reference_root)
    participant_fields, participants = read_tsv(
        model_input_root / "participant_axis.tsv"
    )
    candidate_fields, candidates = read_tsv(
        model_input_root / "candidate_gene_axis.tsv"
    )
    forbidden = {"fibrosis", "outcome", "label", "disease", "sex", "nas", "bmi"}
    if (
        forbidden & set(participant_fields)
        or len(participants) != 109
        or len(candidates) != 14_078
        or len({row["row_id"] for row in participants}) != 109
        or sum(int(row["source_feature_count"]) > 1 for row in candidates) != 11
    ):
        raise GSE268273ConsolidationError("frozen outcome-free axes differ")
    vectors: dict[str, np.ndarray] = {}
    for bundle_id in range(bundle_count):
        bundle = quant_state_root / f"bundle_{bundle_id:02d}"
        verify_frozen_tree(bundle)
        bundle_receipt = json.loads(
            (bundle / "bundle_receipt.json").read_text(encoding="utf-8")
        )
        if (
            bundle_receipt.get("status") != "complete"
            or bundle_receipt.get("target_genes") != 14_078
            or bundle_receipt.get("target_duplicate_sum_groups") != 11
            or bundle_receipt.get("labels_accessed") is not False
            or bundle_receipt.get("fit_or_score_performed") is not False
        ):
            raise GSE268273ConsolidationError("quantification bundle receipt differs")
        for participant in sorted((bundle / "participants").iterdir()):
            if not participant.is_dir():
                continue
            verify_frozen_tree(participant)
            receipt = json.loads(
                (participant / "receipt.json").read_text(encoding="utf-8")
            )
            row_id = str(receipt["row_id"])
            if row_id in vectors:
                raise GSE268273ConsolidationError(
                    "participant appears in multiple quantification bundles"
                )
            values = np.load(
                participant / "v49_expected_counts.npy", allow_pickle=False
            )
            if (
                values.shape != (14_078,)
                or values.dtype != np.float64
                or not np.all(np.isfinite(values))
                or np.any(values < 0)
            ):
                raise GSE268273ConsolidationError(
                    "participant v49 expected-count vector differs"
                )
            vectors[row_id] = values
    row_order = [row["row_id"] for row in participants]
    if set(vectors) != set(row_order):
        raise GSE268273ConsolidationError(
            "quantified participant axis is incomplete or unexpected"
        )
    matrix = np.stack([vectors[row_id] for row_id in row_order])
    if matrix.shape != (109, 14_078):
        raise GSE268273ConsolidationError("consolidated count matrix shape differs")
    output.mkdir(parents=True)
    np.save(output / "rna_expected_counts.npy", matrix, allow_pickle=False)
    np.save(
        output / "rna_observed_mask.npy",
        np.ones(109, dtype=np.bool_),
        allow_pickle=False,
    )
    admitted_participants = [
        {
            **row,
            "rna_observation_state": "observed",
            "quantification_measurement": "RSEM_expected_count_raw_count_scale",
        }
        for row in participants
    ]
    write_tsv(
        output / "participant_axis.tsv",
        tuple(admitted_participants[0]),
        admitted_participants,
    )
    shutil.copyfile(
        model_input_root / "candidate_gene_axis.tsv",
        output / "candidate_gene_axis.tsv",
    )
    receipt = {
        "schema_version": "masld-bench-gse268273-raw-rsem-model-input-v1",
        "status": "passed_outcome_free_quantification_training_transform_blocked",
        "participants": 109,
        "target_stable_genes": 14_078,
        "target_one_to_one_genes": 14_067,
        "target_duplicate_sum_groups": 11,
        "shape": [109, 14_078],
        "dtype": "float64",
        "measurement": "RSEM_expected_count_raw_count_scale",
        "normalization_applied": False,
        "measured_zeros_are_values": True,
        "missing_encoded_as_zero": False,
        "rna_observed_mask_all_true": True,
        "labels_included": False,
        "clinical_covariates_included": False,
        "outcomes_accessed": False,
        "fit_or_score_performed": False,
        "source_global_voom_used": False,
        "processed_differential_expression_used": False,
        "next_gate": "apply_a_preselected_training_frozen_transform_without_target_outcomes",
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(
        output,
        {
            "artifact_class": "gse268273_outcome_free_raw_rsem_model_input",
            "participants": 109,
            "target_genes": 14_078,
            "status": "passed_training_transform_blocked",
        },
    )
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--quant-state-root", type=Path, required=True)
    parser.add_argument("--model-input-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bundles", type=int, default=8)
    args = parser.parse_args()
    receipt = consolidate(
        plan_root=args.plan_root,
        quant_state_root=args.quant_state_root,
        model_input_root=args.model_input_root,
        reference_root=args.reference_root,
        output=args.output,
        bundle_count=args.bundles,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
