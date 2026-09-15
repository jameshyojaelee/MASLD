#!/usr/bin/env python3
"""Independently audit GSE268273 raw-RSEM source mechanics without outcomes.

This auditor does not import the provisional inclusion implementation.  It
re-derives the count matrix from frozen participant quantifications, checks it
against the pre-download axes and raw plan, and audits the included code surface
for model, endpoint, prediction, or transform execution.  Its receipt grants no
permission to fit, select, score, unblind, or make a best-model claim.
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


PARTICIPANT_FIELDS = (
    "row_index",
    "row_id",
    "cohort_family_id",
    "rna_observation_state",
    "raw_technical_run_count",
    "raw_fastq_bytes",
)
SOURCE_PARTICIPANT_FIELDS = (*PARTICIPANT_FIELDS, "quantification_measurement")
PLAN_PARTICIPANT_FIELDS = (
    "bundle_id",
    "row_id",
    "experiment_accession",
    "biosample_accession",
    "technical_runs",
    "fastq_files",
    "fastq_bytes",
    "effective_library_layout",
    "aggregation_rule",
)
GENE_FIELDS = (
    "feature_index",
    "stable_gene_id",
    "source_feature_count",
    "source_feature_ids",
    "gencode_v49_gene_name",
    "gencode_v49_gene_type",
    "raw_count_aggregation",
    "admission_state",
)
ALLOWED_IMPORT_ROOTS = frozenset(
    {"__future__", "argparse", "ast", "csv", "json", "pathlib", "typing", "numpy", "masld_bench"}
)
FORBIDDEN_CALLS = frozenset(
    {
        "fit",
        "fit_predict",
        "fit_transform",
        "predict",
        "predict_proba",
        "score",
        "train",
        "transform",
    }
)
FORBIDDEN_PATH_SURFACES = (
    "evaluator_only",
    "outcomes.tsv",
    "predictions.tsv",
    "prediction_bundle",
    "config/evaluation",
    "gse268273_fibrosis_ood_transfer_task",
)
SOURCE_RECEIPT_KEYS = frozenset(
    {
        "schema_version",
        "status",
        "participants",
        "target_stable_genes",
        "target_one_to_one_genes",
        "target_duplicate_sum_groups",
        "shape",
        "dtype",
        "measurement",
        "normalization_applied",
        "measured_zeros_are_values",
        "missing_encoded_as_zero",
        "rna_observed_mask_all_true",
        "labels_included",
        "clinical_covariates_included",
        "outcomes_accessed",
        "fit_or_score_performed",
        "source_global_voom_used",
        "processed_differential_expression_used",
        "next_gate",
    }
)


class GSE268273IndependentSourceAuditError(ValueError):
    """Raised when the independently derived source-only requirement differs."""


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE268273IndependentSourceAuditError(f"TSV lacks header: {path}")
        return tuple(reader.fieldnames), list(reader)


def require_pinned_tree(root: Path, expected_sha256: str, label: str) -> dict[str, Any]:
    if sha256_file(root / "ARTIFACTS.json") != expected_sha256:
        raise GSE268273IndependentSourceAuditError(f"{label} ARTIFACTS hash differs")
    return verify_frozen_tree(root)


def require_exact_keys(observed: dict[str, Any], expected: Iterable[str], label: str) -> None:
    if set(observed) != set(expected):
        raise GSE268273IndependentSourceAuditError(f"{label} keys differ")


def audit_python_surface(path: Path, expected_sha256: str) -> dict[str, Any]:
    if sha256_file(path) != expected_sha256:
        raise GSE268273IndependentSourceAuditError("provisional admission source hash differs")
    source = path.read_text(encoding="utf-8")
    lowered = source.lower()
    for forbidden in FORBIDDEN_PATH_SURFACES:
        if forbidden in lowered:
            raise GSE268273IndependentSourceAuditError(
                f"provisional admission opens forbidden surface: {forbidden}"
            )
    tree = ast.parse(source, filename=path.as_posix())
    imports: set[str] = set()
    calls: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name.split(".", 1)[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module.split(".", 1)[0])
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Attribute):
                calls.add(node.func.attr)
            elif isinstance(node.func, ast.Name):
                calls.add(node.func.id)
    if imports - ALLOWED_IMPORT_ROOTS:
        raise GSE268273IndependentSourceAuditError(
            f"provisional admission imports modeling surface: {sorted(imports - ALLOWED_IMPORT_ROOTS)}"
        )
    if calls & FORBIDDEN_CALLS:
        raise GSE268273IndependentSourceAuditError(
            f"provisional admission calls modeling/transform surface: {sorted(calls & FORBIDDEN_CALLS)}"
        )
    return {"imports": sorted(imports), "call_names": sorted(calls)}


def audit_shell_surface(
    path: Path,
    expected_sha256: str,
    label: str,
    *,
    allowed_negative_assertion_mentions: dict[str, int] | None = None,
) -> None:
    if sha256_file(path) != expected_sha256:
        raise GSE268273IndependentSourceAuditError(f"{label} source hash differs")
    lowered = path.read_text(encoding="utf-8").lower()
    allowed_mentions = allowed_negative_assertion_mentions or {}
    for forbidden in FORBIDDEN_PATH_SURFACES:
        observed_mentions = lowered.count(forbidden)
        if observed_mentions != allowed_mentions.get(forbidden, 0):
            raise GSE268273IndependentSourceAuditError(
                f"{label} forbidden-surface mention count differs: {forbidden}"
            )
    active_lines = tuple(line.strip().lower() for line in lowered.splitlines())
    for forbidden in ("#sbatch --array", "#sbatch --partition=gpu"):
        if any(line.startswith(forbidden) for line in active_lines):
            raise GSE268273IndependentSourceAuditError(
                f"{label} adds non-source execution surface: {forbidden}"
            )


def audit_submission_record(path: Path, expected_job_id: int) -> str:
    fields, rows = read_tsv(path)
    expected_fields = ("JobID", "JobName", "Partition", "QOS", "State", "ExitCode", "SubmitLine")
    if fields != expected_fields or len(rows) != 1:
        raise GSE268273IndependentSourceAuditError("admission SLURM record schema differs")
    row = rows[0]
    submit_line = row["SubmitLine"]
    required = (
        str(expected_job_id),
        "model-data-081",
        "cpu",
        "nslab",
        "COMPLETED",
        "0:0",
        "--dependency=afterok:21083342",
        "WRAPPER_VALIDATION_ROOT=",
        "model-validate-209-21088706",
        "WRAPPER_VALIDATION_SHA256=865ea6a736d6c004c2db1a7ba8ff81d40e3a73a69f0cad51375f319b74661924",
        "slurm/admit_gse268273_raw_source_integrity_cpu.sbatch",
    )
    values = (row["JobID"], row["JobName"], row["Partition"], row["QOS"], row["State"], row["ExitCode"], submit_line)
    for item in required:
        if not any(item in value for value in values):
            raise GSE268273IndependentSourceAuditError(
                f"admission SLURM record lacks frozen contract: {item}"
            )
    lowered = submit_line.lower()
    for forbidden in FORBIDDEN_PATH_SURFACES:
        if forbidden in lowered:
            raise GSE268273IndependentSourceAuditError(
                f"admission submission opens forbidden surface: {forbidden}"
            )
    return sha256_file(path)


def run(
    *,
    consolidation_root: Path,
    quant_state_root: Path,
    prior_model_input_root: Path,
    plan_root: Path,
    reference_root: Path,
    admission_execution_root: Path,
    admission_source: Path,
    admission_wrapper: Path,
    core_validation_wrapper: Path,
    production_validation_wrapper: Path,
    submission_record: Path,
    output: Path,
    expected_prior_model_input_artifacts_sha256: str,
    expected_plan_artifacts_sha256: str,
    expected_reference_artifacts_sha256: str,
    expected_admission_source_sha256: str,
    expected_admission_wrapper_sha256: str,
    expected_core_validation_wrapper_sha256: str,
    expected_production_validation_wrapper_sha256: str,
    expected_core_validation_artifacts_sha256: str,
    expected_wrapper_validation_artifacts_sha256: str,
    expected_admission_job_id: int,
    expected_participants: int,
    expected_genes: int,
    expected_technical_runs: int,
    expected_fastq_files: int,
    expected_fastq_bytes: int,
    expected_bundles: int,
) -> None:
    if output.exists():
        raise GSE268273IndependentSourceAuditError("refusing to overwrite audit receipt")

    prior_manifest = require_pinned_tree(
        prior_model_input_root,
        expected_prior_model_input_artifacts_sha256,
        "pre-download model input",
    )
    plan_manifest = require_pinned_tree(plan_root, expected_plan_artifacts_sha256, "raw plan")
    reference_manifest = require_pinned_tree(
        reference_root, expected_reference_artifacts_sha256, "RSEM reference"
    )
    consolidation_manifest = verify_frozen_tree(consolidation_root)
    source_root = consolidation_root / "model_input"
    source_manifest = verify_frozen_tree(source_root)
    admission_manifest = verify_frozen_tree(admission_execution_root)
    admitted_root = admission_execution_root / "admission"
    admitted_manifest = verify_frozen_tree(admitted_root)

    if (
        prior_manifest.get("metadata", {}).get("labels_included") is not False
        or prior_manifest.get("metadata", {}).get("rna_values_present") is not False
        or plan_manifest.get("metadata", {}).get("participants") != expected_participants
        or plan_manifest.get("metadata", {}).get("fastq_files") != expected_fastq_files
        or plan_manifest.get("metadata", {}).get("fastq_bytes") != expected_fastq_bytes
        or reference_manifest.get("metadata", {}).get("target_stable_genes") != expected_genes
        or reference_manifest.get("metadata", {}).get("status") != "passed"
        or consolidation_manifest.get("metadata", {}).get("artifact_class")
        != "gse268273_raw_rsem_consolidation_execution"
        or consolidation_manifest.get("metadata", {}).get("status") != "passed"
        or source_manifest.get("metadata", {}).get("artifact_class")
        != "gse268273_outcome_free_raw_rsem_model_input"
        or source_manifest.get("metadata", {}).get("participants") != expected_participants
        or source_manifest.get("metadata", {}).get("target_genes") != expected_genes
        or source_manifest.get("metadata", {}).get("status")
        != "passed_training_transform_blocked"
        or admission_manifest.get("metadata", {}).get("artifact_class")
        != "gse268273_raw_source_integrity_admission_execution"
        or admission_manifest.get("metadata", {}).get("source_integrity_only") is not True
        or admission_manifest.get("metadata", {}).get("clean_or_sealed_champion_eligible") is not False
        or admission_manifest.get("metadata", {}).get("independent_unexposed_rederivation_or_audit_required") is not True
        or admission_manifest.get("metadata", {}).get("core_validation_artifacts_sha256")
        != expected_core_validation_artifacts_sha256
        or admission_manifest.get("metadata", {}).get("wrapper_validation_artifacts_sha256")
        != expected_wrapper_validation_artifacts_sha256
        or admitted_manifest.get("metadata", {}).get("source_integrity_only") is not True
        or admitted_manifest.get("metadata", {}).get("clean_or_sealed_champion_eligible") is not False
    ):
        raise GSE268273IndependentSourceAuditError("frozen source-only metadata differs")

    source_receipt = json.loads((source_root / "receipt.json").read_text(encoding="utf-8"))
    require_exact_keys(source_receipt, SOURCE_RECEIPT_KEYS, "consolidated source receipt")
    if (
        source_receipt.get("status")
        != "passed_outcome_free_quantification_training_transform_blocked"
        or source_receipt.get("participants") != expected_participants
        or source_receipt.get("target_stable_genes") != expected_genes
        or source_receipt.get("shape") != [expected_participants, expected_genes]
        or source_receipt.get("dtype") != "float64"
        or source_receipt.get("measurement") != "RSEM_expected_count_raw_count_scale"
        or source_receipt.get("normalization_applied") is not False
        or source_receipt.get("measured_zeros_are_values") is not True
        or source_receipt.get("missing_encoded_as_zero") is not False
        or source_receipt.get("rna_observed_mask_all_true") is not True
        or source_receipt.get("labels_included") is not False
        or source_receipt.get("clinical_covariates_included") is not False
        or source_receipt.get("outcomes_accessed") is not False
        or source_receipt.get("fit_or_score_performed") is not False
        or source_receipt.get("source_global_voom_used") is not False
        or source_receipt.get("processed_differential_expression_used") is not False
    ):
        raise GSE268273IndependentSourceAuditError("consolidated source receipt differs")

    prior_fields, prior_rows = read_tsv(prior_model_input_root / "participant_axis.tsv")
    plan_fields, plan_rows = read_tsv(plan_root / "participants.tsv")
    source_fields, source_rows = read_tsv(source_root / "participant_axis.tsv")
    if (
        prior_fields != PARTICIPANT_FIELDS
        or plan_fields != PLAN_PARTICIPANT_FIELDS
        or source_fields != SOURCE_PARTICIPANT_FIELDS
        or len(prior_rows) != expected_participants
        or len(plan_rows) != expected_participants
        or len(source_rows) != expected_participants
    ):
        raise GSE268273IndependentSourceAuditError("participant source schema differs")
    prior_ids = [row["row_id"] for row in prior_rows]
    source_ids = [row["row_id"] for row in source_rows]
    plan_by_id = {row["row_id"]: row for row in plan_rows}
    if (
        len(set(prior_ids)) != expected_participants
        or prior_ids != source_ids
        or set(plan_by_id) != set(prior_ids)
    ):
        raise GSE268273IndependentSourceAuditError("participant identity/order differs")
    for prior, source in zip(prior_rows, source_rows, strict=True):
        expected_source = dict(prior)
        expected_source["rna_observation_state"] = "observed"
        expected_source["quantification_measurement"] = "RSEM_expected_count_raw_count_scale"
        if source != expected_source:
            raise GSE268273IndependentSourceAuditError("participant axis mutated beyond source state")
        planned = plan_by_id[prior["row_id"]]
        if (
            int(prior["raw_technical_run_count"]) != int(planned["technical_runs"])
            or int(prior["raw_fastq_bytes"]) != int(planned["fastq_bytes"])
        ):
            raise GSE268273IndependentSourceAuditError("participant raw plan totals differ")
    if (
        sum(int(row["technical_runs"]) for row in plan_rows) != expected_technical_runs
        or sum(int(row["fastq_files"]) for row in plan_rows) != expected_fastq_files
        or sum(int(row["fastq_bytes"]) for row in plan_rows) != expected_fastq_bytes
    ):
        raise GSE268273IndependentSourceAuditError("raw plan totals differ")

    prior_gene_fields, prior_genes = read_tsv(prior_model_input_root / "candidate_gene_axis.tsv")
    source_gene_fields, source_genes = read_tsv(source_root / "candidate_gene_axis.tsv")
    if (
        prior_gene_fields != GENE_FIELDS
        or source_gene_fields != GENE_FIELDS
        or prior_genes != source_genes
        or len(source_genes) != expected_genes
        or [int(row["feature_index"]) for row in source_genes] != list(range(expected_genes))
        or len({row["stable_gene_id"] for row in source_genes}) != expected_genes
        or sha256_file(source_root / "candidate_gene_axis.tsv")
        != sha256_file(prior_model_input_root / "candidate_gene_axis.tsv")
    ):
        raise GSE268273IndependentSourceAuditError("candidate gene source axis differs")

    counts = np.load(source_root / "rna_expected_counts.npy", mmap_mode="r", allow_pickle=False)
    mask = np.load(source_root / "rna_observed_mask.npy", allow_pickle=False)
    if (
        counts.shape != (expected_participants, expected_genes)
        or counts.dtype != np.float64
        or mask.shape != (expected_participants,)
        or mask.dtype != np.bool_
        or not np.all(mask)
    ):
        raise GSE268273IndependentSourceAuditError("source matrix shape, dtype, or mask differs")

    vectors: dict[str, np.ndarray] = {}
    quant_bundle_artifact_hashes: list[str] = []
    for bundle_id in range(expected_bundles):
        bundle = quant_state_root / f"bundle_{bundle_id:02d}"
        bundle_manifest = verify_frozen_tree(bundle)
        quant_bundle_artifact_hashes.append(sha256_file(bundle / "ARTIFACTS.json"))
        bundle_receipt = json.loads((bundle / "bundle_receipt.json").read_text(encoding="utf-8"))
        if (
            bundle_manifest.get("metadata", {}).get("status") != "passed"
            or bundle_receipt.get("status") != "complete"
            or bundle_receipt.get("target_genes") != expected_genes
            or bundle_receipt.get("labels_accessed") is not False
            or bundle_receipt.get("fit_or_score_performed") is not False
        ):
            raise GSE268273IndependentSourceAuditError("quantification bundle contract differs")
        for participant in sorted((bundle / "participants").iterdir()):
            if not participant.is_dir():
                continue
            verify_frozen_tree(participant)
            receipt = json.loads((participant / "receipt.json").read_text(encoding="utf-8"))
            row_id = str(receipt.get("row_id"))
            if row_id in vectors or row_id not in plan_by_id:
                raise GSE268273IndependentSourceAuditError("quantified participant identity differs")
            if receipt.get("labels_accessed") is not False or receipt.get("fit_or_score_performed") is not False:
                raise GSE268273IndependentSourceAuditError("quantified participant source-only firewall differs")
            vector = np.load(participant / "v49_expected_counts.npy", allow_pickle=False)
            if (
                vector.shape != (expected_genes,)
                or vector.dtype != np.float64
                or not np.all(np.isfinite(vector))
                or np.any(vector < 0.0)
            ):
                raise GSE268273IndependentSourceAuditError("quantified participant vector differs")
            vectors[row_id] = vector
    if set(vectors) != set(prior_ids):
        raise GSE268273IndependentSourceAuditError("quantification participant census differs")
    for index, row_id in enumerate(prior_ids):
        if not np.array_equal(np.asarray(counts[index]), vectors[row_id], equal_nan=False):
            raise GSE268273IndependentSourceAuditError("consolidated matrix differs from quantification")
        if not np.isfinite(counts[index]).all() or np.any(counts[index] < 0.0) or counts[index].sum() <= 0.0:
            raise GSE268273IndependentSourceAuditError("consolidated source values are invalid")

    provisional_surface = audit_python_surface(admission_source, expected_admission_source_sha256)
    audit_shell_surface(admission_wrapper, expected_admission_wrapper_sha256, "admission wrapper")
    audit_shell_surface(
        core_validation_wrapper,
        expected_core_validation_wrapper_sha256,
        "core validation wrapper",
    )
    audit_shell_surface(
        production_validation_wrapper,
        expected_production_validation_wrapper_sha256,
        "production validation wrapper",
        allowed_negative_assertion_mentions={
            "evaluator_only": 1,
            "outcomes.tsv": 1,
            "gse268273_fibrosis_ood_transfer_task": 1,
        },
    )
    submission_record_sha256 = audit_submission_record(submission_record, expected_admission_job_id)

    admitted = json.loads((admitted_root / "source_integrity_admission.json").read_text(encoding="utf-8"))
    if (
        admitted.get("source_artifacts_sha256") != sha256_file(source_root / "ARTIFACTS.json")
        or admitted.get("consolidation_artifacts_sha256")
        != sha256_file(consolidation_root / "ARTIFACTS.json")
        or admitted.get("participants") != expected_participants
        or admitted.get("genes") != expected_genes
        or admitted.get("technical_runs") != expected_technical_runs
        or admitted.get("raw_fastq_bytes") != expected_fastq_bytes
        or admitted.get("source_integrity_only") is not True
        or admitted.get("normalization_applied") is not False
        or admitted.get("transform_selected_or_applied") is not False
        or admitted.get("feature_selection_performed") is not False
        or admitted.get("labels_read") is not False
        or admitted.get("clinical_covariates_read") is not False
        or admitted.get("evaluation_data_opened") is not False
        or admitted.get("fit_or_score_performed") is not False
        or admitted.get("clean_or_sealed_champion_eligible") is not False
    ):
        raise GSE268273IndependentSourceAuditError("provisional admission receipt differs")

    verify_frozen_tree(consolidation_root)
    verify_frozen_tree(source_root)
    verify_frozen_tree(admission_execution_root)
    output.mkdir(mode=0o750)
    receipt = {
        "schema_version": "masld-bench-gse268273-independent-source-integrity-audit-v1",
        "status": "passed_source_integrity_mechanics_only",
        "source_integrity_mechanics_independently_rederived": True,
        "source_integrity_mechanics_eligible": True,
        "eligible_scope": "source_integrity_mechanics_only",
        "participants": expected_participants,
        "genes": expected_genes,
        "technical_runs": expected_technical_runs,
        "fastq_files": expected_fastq_files,
        "raw_fastq_bytes": expected_fastq_bytes,
        "quantification_bundles": expected_bundles,
        "source_artifacts_sha256": sha256_file(source_root / "ARTIFACTS.json"),
        "consolidation_artifacts_sha256": sha256_file(consolidation_root / "ARTIFACTS.json"),
        "admission_execution_artifacts_sha256": sha256_file(admission_execution_root / "ARTIFACTS.json"),
        "quant_bundle_artifacts_sha256": quant_bundle_artifact_hashes,
        "submission_record_sha256": submission_record_sha256,
        "provisional_admission_imports": provisional_surface["imports"],
        "provisional_admission_model_or_transform_calls_found": False,
        "endpoint_decision_encoded": False,
        "model_decision_encoded": False,
        "transform_decision_encoded": False,
        "normalization_applied": False,
        "feature_selection_performed": False,
        "labels_read": False,
        "clinical_covariates_read": False,
        "outcomes_read": False,
        "predictions_read": False,
        "fit_or_score_performed": False,
        "modeling_or_scoring_authorized_by_this_receipt": False,
        "sealed_evaluation_authorized_by_this_receipt": False,
        "clean_or_sealed_champion_eligible_by_this_receipt": False,
        "separate_training_frozen_transform_gate_required": True,
        "separate_modeling_and_evaluation_gates_required": True,
    }
    write_json_exclusive(output / "independent_source_integrity_audit.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse268273_independent_raw_source_integrity_audit",
            "source_integrity_mechanics_independently_rederived": True,
            "source_integrity_mechanics_eligible": True,
            "eligible_scope": "source_integrity_mechanics_only",
            "modeling_or_scoring_authorized": False,
            "clean_or_sealed_champion_eligible": False,
            "status": "passed_source_integrity_mechanics_only",
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consolidation-root", required=True, type=Path)
    parser.add_argument("--quant-state-root", required=True, type=Path)
    parser.add_argument("--prior-model-input-root", required=True, type=Path)
    parser.add_argument("--plan-root", required=True, type=Path)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--admission-execution-root", required=True, type=Path)
    parser.add_argument("--admission-source", required=True, type=Path)
    parser.add_argument("--admission-wrapper", required=True, type=Path)
    parser.add_argument("--core-validation-wrapper", required=True, type=Path)
    parser.add_argument("--production-validation-wrapper", required=True, type=Path)
    parser.add_argument("--submission-record", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-prior-model-input-artifacts-sha256", required=True)
    parser.add_argument("--expected-plan-artifacts-sha256", required=True)
    parser.add_argument("--expected-reference-artifacts-sha256", required=True)
    parser.add_argument("--expected-admission-source-sha256", required=True)
    parser.add_argument("--expected-admission-wrapper-sha256", required=True)
    parser.add_argument("--expected-core-validation-wrapper-sha256", required=True)
    parser.add_argument("--expected-production-validation-wrapper-sha256", required=True)
    parser.add_argument("--expected-core-validation-artifacts-sha256", required=True)
    parser.add_argument("--expected-wrapper-validation-artifacts-sha256", required=True)
    parser.add_argument("--expected-admission-job-id", required=True, type=int)
    parser.add_argument("--expected-participants", required=True, type=int)
    parser.add_argument("--expected-genes", required=True, type=int)
    parser.add_argument("--expected-technical-runs", required=True, type=int)
    parser.add_argument("--expected-fastq-files", required=True, type=int)
    parser.add_argument("--expected-fastq-bytes", required=True, type=int)
    parser.add_argument("--expected-bundles", required=True, type=int)
    arguments = parser.parse_args()
    run(**vars(arguments))
    print(json.dumps({"output": arguments.output.as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
