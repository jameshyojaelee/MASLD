#!/usr/bin/env python3
"""Tests for the independent GSE268273 source-integrity audit."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.hashing import sha256_file
from scripts import audit_gse268273_raw_source_integrity_unexposed as audit


ROOT = Path(__file__).parents[2]


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_fixture(root: Path, *, plan_byte_error: bool = False) -> dict[str, object]:
    participants = 3
    genes = 2
    bundles = 2
    prior = root / "prior"
    plan = root / "plan"
    reference = root / "reference"
    quant = root / "quant"
    consolidation = root / "consolidation"
    source = consolidation / "model_input"
    execution = root / "execution"
    admitted_root = execution / "admission"
    for path in (prior, plan, reference, source, admitted_root):
        path.mkdir(parents=True)

    prior_rows = [
        {
            "row_index": str(index),
            "row_id": f"g268_fixture_{index}",
            "cohort_family_id": "GSE268273",
            "rna_observation_state": "derivable_not_processed",
            "raw_technical_run_count": str(index + 1),
            "raw_fastq_bytes": str((index + 1) * 100),
        }
        for index in range(participants)
    ]
    write_tsv(prior / "participant_axis.tsv", audit.PARTICIPANT_FIELDS, prior_rows)
    gene_rows = [
        {
            "feature_index": str(index),
            "stable_gene_id": f"ENSG{index + 1:011d}",
            "source_feature_count": "1",
            "source_feature_ids": f"ENSG{index + 1:011d}.1",
            "gencode_v49_gene_name": f"GENE{index + 1}",
            "gencode_v49_gene_type": "protein_coding",
            "raw_count_aggregation": "identity",
            "admission_state": "admitted",
        }
        for index in range(genes)
    ]
    write_tsv(prior / "candidate_gene_axis.tsv", audit.GENE_FIELDS, gene_rows)
    write_json(
        prior / "receipt.json",
        {"labels_included": False, "rna_values_present": False},
    )
    np.save(prior / "rna_observed_mask.npy", np.zeros(participants, dtype=np.bool_))
    freeze_tree(
        prior,
        {
            "artifact_class": "gse268273_fibrosis_ood_model_input",
            "labels_included": False,
            "rna_values_present": False,
            "status": "blocked_raw_reprocessing",
        },
    )

    plan_rows = []
    for index, prior_row in enumerate(prior_rows):
        byte_count = int(prior_row["raw_fastq_bytes"]) + (1 if plan_byte_error and index == 1 else 0)
        plan_rows.append(
            {
                "bundle_id": str(index % bundles),
                "row_id": prior_row["row_id"],
                "experiment_accession": f"SRX{index + 1}",
                "biosample_accession": f"SAMN{index + 1}",
                "technical_runs": prior_row["raw_technical_run_count"],
                "fastq_files": prior_row["raw_technical_run_count"],
                "fastq_bytes": str(byte_count),
                "effective_library_layout": "single_end",
                "aggregation_rule": "one_participant_RSEM_library",
            }
        )
    write_tsv(plan / "participants.tsv", audit.PLAN_PARTICIPANT_FIELDS, plan_rows)
    write_json(plan / "receipt.json", {"labels_included": False, "outcomes_accessed": False})
    freeze_tree(
        plan,
        {
            "artifact_class": "gse268273_raw_campaign_plan",
            "participants": participants,
            "fastq_files": sum(int(row["fastq_files"]) for row in plan_rows),
            "fastq_bytes": sum(int(row["fastq_bytes"]) for row in plan_rows),
            "status": "planned_not_submitted",
        },
    )
    write_json(reference / "receipt.json", {"labels_accessed": False, "fit_or_score_performed": False})
    freeze_tree(
        reference,
        {
            "artifact_class": "gse268273_rsem_star_reference",
            "target_stable_genes": genes,
            "status": "passed",
        },
    )

    matrix = np.arange(1, participants * genes + 1, dtype=np.float64).reshape(participants, genes)
    for bundle_id in range(bundles):
        bundle = quant / f"bundle_{bundle_id:02d}"
        (bundle / "participants").mkdir(parents=True)
        member_rows = [row for index, row in enumerate(prior_rows) if index % bundles == bundle_id]
        for row in member_rows:
            index = int(row["row_index"])
            participant = bundle / "participants" / row["row_id"]
            participant.mkdir()
            np.save(participant / "v49_expected_counts.npy", matrix[index], allow_pickle=False)
            write_json(
                participant / "receipt.json",
                {
                    "row_id": row["row_id"],
                    "labels_accessed": False,
                    "fit_or_score_performed": False,
                },
            )
            freeze_tree(
                participant,
                {
                    "artifact_class": "gse268273_participant_rsem_expected_counts",
                    "row_id": row["row_id"],
                    "target_genes": genes,
                    "status": "passed",
                },
            )
        write_json(
            bundle / "bundle_receipt.json",
            {
                "status": "complete",
                "target_genes": genes,
                "labels_accessed": False,
                "fit_or_score_performed": False,
            },
        )
        freeze_tree(
            bundle,
            {
                "artifact_class": "gse268273_rsem_expected_count_bundle",
                "bundle_id": bundle_id,
                "participants": len(member_rows),
                "status": "passed",
            },
        )

    source_rows = [
        {
            **row,
            "rna_observation_state": "observed",
            "quantification_measurement": "RSEM_expected_count_raw_count_scale",
        }
        for row in prior_rows
    ]
    write_tsv(source / "participant_axis.tsv", audit.SOURCE_PARTICIPANT_FIELDS, source_rows)
    write_tsv(source / "candidate_gene_axis.tsv", audit.GENE_FIELDS, gene_rows)
    np.save(source / "rna_expected_counts.npy", matrix, allow_pickle=False)
    np.save(source / "rna_observed_mask.npy", np.ones(participants, dtype=np.bool_), allow_pickle=False)
    write_json(
        source / "receipt.json",
        {
            "schema_version": "masld-bench-gse268273-raw-rsem-model-input-v1",
            "status": "passed_outcome_free_quantification_training_transform_blocked",
            "participants": participants,
            "target_stable_genes": genes,
            "target_one_to_one_genes": genes,
            "target_duplicate_sum_groups": 0,
            "shape": [participants, genes],
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
            "next_gate": "separate_training_frozen_transform_gate",
        },
    )
    freeze_tree(
        source,
        {
            "artifact_class": "gse268273_outcome_free_raw_rsem_model_input",
            "participants": participants,
            "target_genes": genes,
            "status": "passed_training_transform_blocked",
        },
    )
    freeze_tree(
        consolidation,
        {"artifact_class": "gse268273_raw_rsem_consolidation_execution", "status": "passed"},
    )

    admitted = {
        "source_artifacts_sha256": sha256_file(source / "ARTIFACTS.json"),
        "consolidation_artifacts_sha256": sha256_file(consolidation / "ARTIFACTS.json"),
        "participants": participants,
        "genes": genes,
        "technical_runs": sum(int(row["technical_runs"]) for row in plan_rows),
        "raw_fastq_bytes": sum(int(row["fastq_bytes"]) for row in plan_rows),
        "source_integrity_only": True,
        "normalization_applied": False,
        "transform_selected_or_applied": False,
        "feature_selection_performed": False,
        "labels_read": False,
        "clinical_covariates_read": False,
        "evaluation_data_opened": False,
        "fit_or_score_performed": False,
        "clean_or_sealed_champion_eligible": False,
    }
    write_json(admitted_root / "source_integrity_admission.json", admitted)
    freeze_tree(
        admitted_root,
        {
            "artifact_class": "gse268273_raw_rsem_source_integrity_admission",
            "source_integrity_only": True,
            "clean_or_sealed_champion_eligible": False,
            "status": "passed_unexposed_audit_required",
        },
    )
    freeze_tree(
        execution,
        {
            "artifact_class": "gse268273_raw_source_integrity_admission_execution",
            "core_validation_artifacts_sha256": "core-validation",
            "wrapper_validation_artifacts_sha256": "wrapper-validation",
            "source_integrity_only": True,
            "clean_or_sealed_champion_eligible": False,
            "independent_unexposed_rederivation_or_audit_required": True,
            "status": "passed_unexposed_audit_required",
        },
    )

    source_copy = root / "admission.py"
    wrapper_copy = root / "admission.sbatch"
    core_copy = root / "core.sbatch"
    production_copy = root / "production.sbatch"
    shutil.copyfile(ROOT / "scripts/admit_gse268273_raw_rsem_source_integrity.py", source_copy)
    shutil.copyfile(ROOT / "slurm/admit_gse268273_raw_source_integrity_cpu.sbatch", wrapper_copy)
    shutil.copyfile(ROOT / "slurm/validate_gse268273_raw_source_integrity_admission_cpu.sbatch", core_copy)
    shutil.copyfile(ROOT / "slurm/validate_gse268273_raw_source_integrity_production_cpu.sbatch", production_copy)
    submission = root / "submission.tsv"
    write_tsv(
        submission,
        ("JobID", "JobName", "Partition", "QOS", "State", "ExitCode", "SubmitLine"),
        [
            {
                "JobID": "21088737",
                "JobName": "model-data-081",
                "Partition": "cpu",
                "QOS": "nslab",
                "State": "COMPLETED",
                "ExitCode": "0:0",
                "SubmitLine": (
                    "sbatch --dependency=afterok:21083342 --export=ALL,"
                    "WRAPPER_VALIDATION_ROOT=/tmp/model-validate-209-21088706,"
                    "WRAPPER_VALIDATION_SHA256=865ea6a736d6c004c2db1a7ba8ff81d40e3a73a69f0cad51375f319b74661924 "
                    "slurm/admit_gse268273_raw_source_integrity_cpu.sbatch"
                ),
            }
        ],
    )
    return {
        "consolidation_root": consolidation,
        "quant_state_root": quant,
        "prior_model_input_root": prior,
        "plan_root": plan,
        "reference_root": reference,
        "admission_execution_root": execution,
        "admission_source": source_copy,
        "admission_wrapper": wrapper_copy,
        "core_validation_wrapper": core_copy,
        "production_validation_wrapper": production_copy,
        "submission_record": submission,
        "expected_prior_model_input_artifacts_sha256": sha256_file(prior / "ARTIFACTS.json"),
        "expected_plan_artifacts_sha256": sha256_file(plan / "ARTIFACTS.json"),
        "expected_reference_artifacts_sha256": sha256_file(reference / "ARTIFACTS.json"),
        "expected_admission_source_sha256": sha256_file(source_copy),
        "expected_admission_wrapper_sha256": sha256_file(wrapper_copy),
        "expected_core_validation_wrapper_sha256": sha256_file(core_copy),
        "expected_production_validation_wrapper_sha256": sha256_file(production_copy),
        "expected_core_validation_artifacts_sha256": "core-validation",
        "expected_wrapper_validation_artifacts_sha256": "wrapper-validation",
        "expected_admission_job_id": 21088737,
        "expected_participants": participants,
        "expected_genes": genes,
        "expected_technical_runs": sum(int(row["technical_runs"]) for row in plan_rows),
        "expected_fastq_files": sum(int(row["fastq_files"]) for row in plan_rows),
        "expected_fastq_bytes": sum(int(row["fastq_bytes"]) for row in plan_rows),
        "expected_bundles": bundles,
    }


class GSE268273IndependentSourceAuditTests(unittest.TestCase):
    def test_rederives_source_only_mechanics_and_keeps_modeling_blocked(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            arguments = build_fixture(root)
            output = root / "audit"
            audit.run(**arguments, output=output)
            manifest = verify_frozen_tree(output)
            receipt = json.loads((output / "independent_source_integrity_audit.json").read_text())
        self.assertTrue(receipt["source_integrity_mechanics_eligible"])
        self.assertFalse(receipt["modeling_or_scoring_authorized_by_this_receipt"])
        self.assertFalse(receipt["clean_or_sealed_champion_eligible_by_this_receipt"])
        self.assertEqual(manifest["metadata"]["eligible_scope"], "source_integrity_mechanics_only")

    def test_rejects_raw_plan_participant_total_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            arguments = build_fixture(root, plan_byte_error=True)
            output = root / "audit"
            with self.assertRaisesRegex(
                audit.GSE268273IndependentSourceAuditError, "participant raw plan totals"
            ):
                audit.run(**arguments, output=output)

    def test_rejects_model_call_added_to_provisional_surface(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            arguments = build_fixture(root)
            source = arguments["admission_source"]
            assert isinstance(source, Path)
            source.write_text(source.read_text() + "\nobject().fit()\n", encoding="utf-8")
            arguments["expected_admission_source_sha256"] = sha256_file(source)
            with self.assertRaisesRegex(
                audit.GSE268273IndependentSourceAuditError, "modeling/transform surface"
            ):
                audit.run(**arguments, output=root / "audit")

    def test_rejects_evaluator_path_added_to_wrapper(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            arguments = build_fixture(root)
            wrapper = arguments["admission_wrapper"]
            assert isinstance(wrapper, Path)
            wrapper.write_text(wrapper.read_text() + "\nX=evaluator_only\n", encoding="utf-8")
            arguments["expected_admission_wrapper_sha256"] = sha256_file(wrapper)
            with self.assertRaisesRegex(
                audit.GSE268273IndependentSourceAuditError, "forbidden-surface mention"
            ):
                audit.run(**arguments, output=root / "audit")

    def test_rejects_submission_without_frozen_dependency(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            arguments = build_fixture(root)
            submission = arguments["submission_record"]
            assert isinstance(submission, Path)
            text = submission.read_text().replace("afterok:21083342", "afterok:1")
            submission.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(
                audit.GSE268273IndependentSourceAuditError, "lacks frozen contract"
            ):
                audit.run(**arguments, output=root / "audit")


if __name__ == "__main__":
    unittest.main()
