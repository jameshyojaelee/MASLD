#!/usr/bin/env python3
"""Admit a consolidated GSE268273 RSEM matrix for source integrity only.

This stage does not select features, fit preprocessing, open evaluation data,
fit a model, or score predictions. Because it was designed after a documented
evaluator-receipt exposure, its output is not clean/sealed best model eligible
until independently re-derived or audited by an unexposed agent.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


FORBIDDEN_COLUMNS = frozenset(
    {
        "fibrosis",
        "fibrosis_score",
        "nas",
        "nas_score",
        "sex",
        "bmi",
        "outcome",
        "label",
        "disease",
        "clinical_group",
    }
)


class GSE268273SourceIntegrityError(ValueError):
    """Raised when the outcome-free source-integrity requirement differs."""


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE268273SourceIntegrityError(f"TSV lacks a header: {path}")
        return list(reader.fieldnames), list(reader)


def require_frozen_hash(root: Path, expected: str, label: str) -> dict[str, Any]:
    if sha256_file(root / "ARTIFACTS.json") != expected:
        raise GSE268273SourceIntegrityError(f"{label} ARTIFACTS hash differs")
    return verify_frozen_tree(root)


def run(
    *,
    consolidation_root: Path,
    prior_model_input_root: Path,
    plan_root: Path,
    reference_root: Path,
    incident_root: Path,
    output: Path,
    expected_prior_model_input_artifacts_sha256: str,
    expected_plan_artifacts_sha256: str,
    expected_reference_artifacts_sha256: str,
    expected_incident_artifacts_sha256: str,
    expected_participants: int,
    expected_genes: int,
    expected_technical_runs: int,
    expected_fastq_bytes: int,
) -> None:
    if output.exists():
        raise GSE268273SourceIntegrityError("refusing to overwrite source admission")
    prior_manifest = require_frozen_hash(
        prior_model_input_root,
        expected_prior_model_input_artifacts_sha256,
        "prior model-input",
    )
    plan_manifest = require_frozen_hash(
        plan_root, expected_plan_artifacts_sha256, "raw plan"
    )
    reference_manifest = require_frozen_hash(
        reference_root, expected_reference_artifacts_sha256, "RSEM reference"
    )
    incident_manifest = require_frozen_hash(
        incident_root, expected_incident_artifacts_sha256, "incident"
    )
    consolidation_manifest = verify_frozen_tree(consolidation_root)
    source_root = consolidation_root / "model_input"
    source_manifest = verify_frozen_tree(source_root)
    if (
        consolidation_manifest.get("metadata", {}).get("artifact_class")
        != "gse268273_raw_rsem_consolidation_execution"
        or consolidation_manifest.get("metadata", {}).get("status") != "passed"
        or source_manifest.get("metadata", {}).get("artifact_class")
        != "gse268273_outcome_free_raw_rsem_model_input"
        or source_manifest.get("metadata", {}).get("participants")
        != expected_participants
        or source_manifest.get("metadata", {}).get("target_genes") != expected_genes
        or source_manifest.get("metadata", {}).get("status")
        != "passed_training_transform_blocked"
        or prior_manifest.get("metadata", {}).get("labels_included") is not False
        or prior_manifest.get("metadata", {}).get("rna_values_present") is not False
        or plan_manifest.get("metadata", {}).get("participants")
        != expected_participants
        or plan_manifest.get("metadata", {}).get("status") != "planned_not_submitted"
        or reference_manifest.get("metadata", {}).get("target_stable_genes")
        != expected_genes
        or reference_manifest.get("metadata", {}).get("status") != "passed"
        or incident_manifest.get("metadata", {}).get("clean_or_sealed_champion_eligible")
        is not False
        or incident_manifest.get("metadata", {}).get("status")
        != "contained_fail_closed"
    ):
        raise GSE268273SourceIntegrityError("frozen source-integrity metadata differs")

    receipt = json.loads((source_root / "receipt.json").read_text(encoding="utf-8"))
    if (
        receipt.get("schema_version")
        != "masld-bench-gse268273-raw-rsem-model-input-v1"
        or receipt.get("status")
        != "passed_outcome_free_quantification_training_transform_blocked"
        or receipt.get("participants") != expected_participants
        or receipt.get("target_stable_genes") != expected_genes
        or receipt.get("shape") != [expected_participants, expected_genes]
        or receipt.get("dtype") != "float64"
        or receipt.get("measurement") != "RSEM_expected_count_raw_count_scale"
        or receipt.get("normalization_applied") is not False
        or receipt.get("measured_zeros_are_values") is not True
        or receipt.get("missing_encoded_as_zero") is not False
        or receipt.get("rna_observed_mask_all_true") is not True
        or receipt.get("labels_included") is not False
        or receipt.get("clinical_covariates_included") is not False
        or receipt.get("outcomes_accessed") is not False
        or receipt.get("fit_or_score_performed") is not False
        or receipt.get("source_global_voom_used") is not False
        or receipt.get("processed_differential_expression_used") is not False
    ):
        raise GSE268273SourceIntegrityError("consolidated source receipt differs")

    prior_fields, prior_rows = read_tsv(prior_model_input_root / "participant_axis.tsv")
    fields, rows = read_tsv(source_root / "participant_axis.tsv")
    expected_fields = [
        *prior_fields,
        "quantification_measurement",
    ]
    if (
        fields != expected_fields
        or FORBIDDEN_COLUMNS & set(fields)
        or len(rows) != expected_participants
        or len({row["row_id"] for row in rows}) != expected_participants
        or [row["row_id"] for row in rows] != [row["row_id"] for row in prior_rows]
    ):
        raise GSE268273SourceIntegrityError("participant axis schema or order differs")
    for prior, observed in zip(prior_rows, rows, strict=True):
        for field in prior_fields:
            expected = "observed" if field == "rna_observation_state" else prior[field]
            if observed[field] != expected:
                raise GSE268273SourceIntegrityError(
                    "participant source identity differs after quantification"
                )
        if (
            observed["quantification_measurement"]
            != "RSEM_expected_count_raw_count_scale"
        ):
            raise GSE268273SourceIntegrityError("participant measurement differs")
    if (
        sum(int(row["raw_technical_run_count"]) for row in rows)
        != expected_technical_runs
        or sum(int(row["raw_fastq_bytes"]) for row in rows) != expected_fastq_bytes
    ):
        raise GSE268273SourceIntegrityError("participant source totals differ")

    prior_gene_axis = prior_model_input_root / "candidate_gene_axis.tsv"
    gene_axis = source_root / "candidate_gene_axis.tsv"
    if sha256_file(gene_axis) != sha256_file(prior_gene_axis):
        raise GSE268273SourceIntegrityError("candidate gene axis changed")
    gene_fields, genes = read_tsv(gene_axis)
    if (
        len(genes) != expected_genes
        or len({row["stable_gene_id"] for row in genes}) != expected_genes
        or FORBIDDEN_COLUMNS & set(gene_fields)
    ):
        raise GSE268273SourceIntegrityError("candidate gene source axis differs")

    counts = np.load(source_root / "rna_expected_counts.npy", mmap_mode="r", allow_pickle=False)
    mask = np.load(source_root / "rna_observed_mask.npy", allow_pickle=False)
    if (
        counts.shape != (expected_participants, expected_genes)
        or counts.dtype != np.float64
        or mask.shape != (expected_participants,)
        or mask.dtype != np.bool_
        or not np.all(mask)
    ):
        raise GSE268273SourceIntegrityError("source matrix shape, dtype, or mask differs")
    for start in range(0, expected_participants, 16):
        block = counts[start : start + 16]
        if (
            not np.all(np.isfinite(block))
            or np.any(block < 0.0)
            or np.any(block.sum(axis=1) <= 0.0)
        ):
            raise GSE268273SourceIntegrityError("source count values are invalid")

    output.mkdir(mode=0o750)
    admission = {
        "schema_version": "masld-bench-gse268273-raw-source-integrity-admission-v1",
        "status": "passed_source_integrity_only_unexposed_audit_required",
        "source_artifacts_sha256": sha256_file(source_root / "ARTIFACTS.json"),
        "consolidation_artifacts_sha256": sha256_file(
            consolidation_root / "ARTIFACTS.json"
        ),
        "prior_model_input_artifacts_sha256": (
            expected_prior_model_input_artifacts_sha256
        ),
        "plan_artifacts_sha256": expected_plan_artifacts_sha256,
        "reference_artifacts_sha256": expected_reference_artifacts_sha256,
        "incident_artifacts_sha256": expected_incident_artifacts_sha256,
        "participants": expected_participants,
        "genes": expected_genes,
        "technical_runs": expected_technical_runs,
        "raw_fastq_bytes": expected_fastq_bytes,
        "shape": [expected_participants, expected_genes],
        "dtype": "float64",
        "measurement": "RSEM_expected_count_raw_count_scale",
        "participant_axis_sha256": sha256_file(source_root / "participant_axis.tsv"),
        "candidate_gene_axis_sha256": sha256_file(gene_axis),
        "rna_expected_counts_sha256": sha256_file(
            source_root / "rna_expected_counts.npy"
        ),
        "rna_observed_mask_sha256": sha256_file(
            source_root / "rna_observed_mask.npy"
        ),
        "participant_axis_matches_pre_download_source": True,
        "candidate_gene_axis_matches_pre_download_source": True,
        "normalization_applied": False,
        "transform_selected_or_applied": False,
        "feature_selection_performed": False,
        "labels_read": False,
        "clinical_covariates_read": False,
        "evaluation_data_opened": False,
        "fit_or_score_performed": False,
        "source_integrity_only": True,
        "clean_or_sealed_champion_eligible": False,
        "independent_unexposed_rederivation_or_audit_required": True,
    }
    write_json_exclusive(output / "source_integrity_admission.json", admission)
    freeze_tree(
        output,
        {
            "artifact_class": "gse268273_raw_rsem_source_integrity_admission",
            "participants": expected_participants,
            "genes": expected_genes,
            "source_integrity_only": True,
            "clean_or_sealed_champion_eligible": False,
            "independent_unexposed_rederivation_or_audit_required": True,
            "status": "passed_unexposed_audit_required",
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consolidation-root", required=True, type=Path)
    parser.add_argument("--prior-model-input-root", required=True, type=Path)
    parser.add_argument("--plan-root", required=True, type=Path)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--incident-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-prior-model-input-artifacts-sha256", required=True)
    parser.add_argument("--expected-plan-artifacts-sha256", required=True)
    parser.add_argument("--expected-reference-artifacts-sha256", required=True)
    parser.add_argument("--expected-incident-artifacts-sha256", required=True)
    parser.add_argument("--expected-participants", required=True, type=int)
    parser.add_argument("--expected-genes", required=True, type=int)
    parser.add_argument("--expected-technical-runs", required=True, type=int)
    parser.add_argument("--expected-fastq-bytes", required=True, type=int)
    arguments = parser.parse_args()
    run(
        consolidation_root=arguments.consolidation_root,
        prior_model_input_root=arguments.prior_model_input_root,
        plan_root=arguments.plan_root,
        reference_root=arguments.reference_root,
        incident_root=arguments.incident_root,
        output=arguments.output,
        expected_prior_model_input_artifacts_sha256=(
            arguments.expected_prior_model_input_artifacts_sha256
        ),
        expected_plan_artifacts_sha256=arguments.expected_plan_artifacts_sha256,
        expected_reference_artifacts_sha256=(
            arguments.expected_reference_artifacts_sha256
        ),
        expected_incident_artifacts_sha256=arguments.expected_incident_artifacts_sha256,
        expected_participants=arguments.expected_participants,
        expected_genes=arguments.expected_genes,
        expected_technical_runs=arguments.expected_technical_runs,
        expected_fastq_bytes=arguments.expected_fastq_bytes,
    )
    print(json.dumps({"output": arguments.output.as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
