#!/usr/bin/env python3
"""Finalize the label-free GSE105127 production activation without fitting."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


class GSE105127FinalizationError(RuntimeError):
    """Raised when final no-fit activation evidence is incomplete."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        forbidden = {"phenotype", "label", "outcome", "disease", "fibrosis", "nas", "sex", "age", "bmi"}
        if forbidden & set(reader.fieldnames or ()):
            raise GSE105127FinalizationError("finalization source violates label firewall")
        return list(reader)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise GSE105127FinalizationError("refusing to write empty final authority")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def finalize(
    *, plan_root: Path, reference_root: Path, collapsed_root: Path,
    crosswalk_root: Path, rna_root: Path, output: Path,
) -> dict[str, object]:
    if output.exists():
        raise GSE105127FinalizationError(f"refusing to overwrite final activation: {output}")
    for root in (plan_root, reference_root, crosswalk_root, rna_root):
        verify_frozen_tree(root)
    plan = json.loads((plan_root / "plan.json").read_text(encoding="utf-8"))
    reference = json.loads((reference_root / "receipt.json").read_text(encoding="utf-8"))
    crosswalk = json.loads((crosswalk_root / "receipt.json").read_text(encoding="utf-8"))
    rna = json.loads((rna_root / "receipt.json").read_text(encoding="utf-8"))
    if (
        plan.get("status") != "planned_no_fit"
        or reference.get("status") != "passed"
        or crosswalk.get("status") != "passed_failure_aware_roundtrip"
        or crosswalk.get("all_source_intervals_retained") is not True
        or int(crosswalk.get("mapping_states", {}).get("mapped_unique_cpg", 0)) <= 0
        or rna.get("status") != "passed_label_free_raw_scale"
        or any(value.get("labels_accessed") is not False for value in (plan, reference, crosswalk, rna))
        or any(value.get("fit_or_score_performed") is not False for value in (plan, reference, crosswalk, rna))
    ):
        raise GSE105127FinalizationError("production authority status differs")
    rows = read_rows(plan_root / "rna_rows.tsv")
    rrbs_rows = {row["row_id"]: row for row in read_rows(plan_root / "rrbs_rows.tsv")}
    if len(rows) != 57 or set(rrbs_rows) != {row["row_id"] for row in rows}:
        raise GSE105127FinalizationError("production row axes differ")
    rna_axis = read_rows(rna_root / "row_axis.tsv")
    if (
        len(rna_axis) != 57
        or [int(row["row_index"]) for row in rna_axis] != list(range(57))
        or [row["row_id"] for row in rna_axis]
        != [row["row_id"] for row in sorted(rows, key=lambda value: (value["participant_group_id"], value["zone"]))]
    ):
        raise GSE105127FinalizationError("consolidated RNA row axis differs")
    participants = sorted({row["participant_group_id"] for row in rows})
    if len(participants) != 19:
        raise GSE105127FinalizationError("participant census differs")
    fold_by_participant = {participant: index % 5 for index, participant in enumerate(participants)}
    fold_counts = {fold: sum(value == fold for value in fold_by_participant.values()) for fold in range(5)}
    if sorted(fold_counts.values()) != [3, 4, 4, 4, 4]:
        raise GSE105127FinalizationError("opaque participant fold counts differ")
    ordered_rows = sorted(rows, key=lambda value: (value["participant_group_id"], value["zone"]))
    row_index = {row["row_id"]: index for index, row in enumerate(ordered_rows)}
    model_rows = []
    fold_rows = []
    for row in ordered_rows:
        collapsed = collapsed_root / "participants" / row["row_id"]
        verify_frozen_tree(collapsed)
        collapsed_receipt = json.loads((collapsed / "receipt.json").read_text(encoding="utf-8"))
        if (
            collapsed_receipt.get("status") != "passed"
            or collapsed_receipt.get("row_id") != row["row_id"]
            or collapsed_receipt.get("participant_group_id") != row["participant_group_id"]
            or collapsed_receipt.get("zone") != row["zone"]
            or collapsed_receipt.get("labels_accessed") is not False
            or collapsed_receipt.get("fit_or_score_performed") is not False
        ):
            raise GSE105127FinalizationError("collapsed RRBS receipt violates firewall")
        model_rows.append({
            "row_id": row["row_id"],
            "participant_group_id": row["participant_group_id"],
            "zone": row["zone"],
            "pairing_topology": "adjacent_section",
            "rna_row_index": row_index[row["row_id"]],
            "rrbs_collapsed_path": str(collapsed / "cpg_counts.hg19.tsv.gz"),
            "rrbs_collapsed_artifacts_sha256": sha256_file(collapsed / "ARTIFACTS.json"),
            "cpg_crosswalk_path": str(crosswalk_root / "cpg_crosswalk.tsv.gz"),
            "rna_matrix_root": str(rna_root),
        })
        fold_rows.append({
            "row_id": row["row_id"],
            "participant_group_id": row["participant_group_id"],
            "zone": row["zone"],
            "participant_outer_fold": fold_by_participant[row["participant_group_id"]],
            "fold_assignment_uses_outcomes": False,
            "available_to_model_as_feature": False,
        })
    output.mkdir(parents=True)
    write_tsv(output / "model_input_manifest.tsv", model_rows)
    write_tsv(output / "participant_folds.tsv", fold_rows)
    evaluator = {
        "schema_version": "masld-bench-gse105127-evaluator-contract-v1",
        "status": "frozen_contract_implementation_pending",
        "unit_of_inference": "participant",
        "repeated_measurements": ["zone", "adjacent_section", "library", "CpG"],
        "target": "RRBS_methylated_reads_conditional_on_observed_coverage",
        "primary_metric": "participant_block_macro_beta_binomial_deviance_skill",
        "dispersion_fit": "outer_training_participants_only",
        "missing_CpG": "masked_never_zero",
        "genomic_split": "held_target_chromosome_or_prespecified_block_with_receptive_field_buffer",
        "participant_folds": 5,
        "participant_fold_counts": fold_counts,
        "labels_or_outcomes_used_for_fold_assignment": False,
        "observed_RRBS_query_lane": "separate_from_RNA_only_lane",
        "fit_or_score_performed": False,
    }
    (output / "evaluator_contract.json").write_text(json.dumps(evaluator, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    receipt = {
        "schema_version": "masld-bench-gse105127-production-activation-v1",
        "status": "assays_ready_no_fit_evaluator_implementation_pending",
        "participants": 19,
        "participant_zone_rows": 57,
        "pairing_topology": "adjacent_section",
        "same_section": False,
        "source_reference_frozen": True,
        "chains_frozen": True,
        "strand_pair_scan_complete": True,
        "failure_aware_roundtrip_crosswalk_complete": True,
        "rna_fastq_md5_gzip_verified": True,
        "rna_rsem_star_quantification_complete": True,
        "rna_consolidation_complete": True,
        "participant_folds_frozen_without_outcomes": True,
        "labels_accessed": False,
        "rsem_native_tpm_available": True,
        "post_quantification_normalization_run": False,
        "model_fit_run": False,
        "model_score_run": False,
        "remaining_blockers": [
            "implement_and_fixture_test_the_frozen_coverage_aware_evaluator",
            "freeze_model_receptive_field_specific_genomic_blocks_and_buffers_before_each_campaign",
        ],
        "source_artifacts": {
            "plan": sha256_file(plan_root / "ARTIFACTS.json"),
            "reference": sha256_file(reference_root / "ARTIFACTS.json"),
            "crosswalk": sha256_file(crosswalk_root / "ARTIFACTS.json"),
            "rna": sha256_file(rna_root / "ARTIFACTS.json"),
        },
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    freeze_tree(output, {"artifact_class": "gse105127_production_activation", "status": receipt["status"]})
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--collapsed-root", type=Path, required=True)
    parser.add_argument("--crosswalk-root", type=Path, required=True)
    parser.add_argument("--rna-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(finalize(
        plan_root=args.plan_root, reference_root=args.reference_root,
        collapsed_root=args.collapsed_root, crosswalk_root=args.crosswalk_root,
        rna_root=args.rna_root, output=args.output,
    ), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
