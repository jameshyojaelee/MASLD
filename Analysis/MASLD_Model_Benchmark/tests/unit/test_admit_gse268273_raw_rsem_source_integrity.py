#!/usr/bin/env python3
"""Tests for outcome-free GSE268273 raw RSEM source admission."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.hashing import sha256_file
from scripts import admit_gse268273_raw_rsem_source_integrity as admission


ROOT = Path(__file__).parents[2]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def build_fixture(
    root: Path,
    *,
    forbidden_column: bool = False,
    reverse_participants: bool = False,
    count_problem: str | None = None,
    incident_eligible: bool = False,
) -> dict[str, object]:
    participants = 5
    genes = 3
    prior = root / "prior"
    plan = root / "plan"
    reference = root / "reference"
    incident = root / "incident"
    consolidation = root / "consolidation"
    source = consolidation / "model_input"
    for path in (prior, plan, reference, incident, source):
        path.mkdir(parents=True)

    prior_fields = [
        "row_id",
        "raw_technical_run_count",
        "raw_fastq_bytes",
        "rna_observation_state",
    ]
    prior_rows = [
        {
            "row_id": f"participant-{index + 1}",
            "raw_technical_run_count": str(index + 1),
            "raw_fastq_bytes": str((index + 1) * 100),
            "rna_observation_state": "derivable_not_processed",
        }
        for index in range(participants)
    ]
    write_tsv(prior / "participant_axis.tsv", prior_fields, prior_rows)
    gene_fields = ["stable_gene_id", "source_feature_count"]
    gene_rows = [
        {"stable_gene_id": f"ENSG{index + 1}", "source_feature_count": "1"}
        for index in range(genes)
    ]
    write_tsv(prior / "candidate_gene_axis.tsv", gene_fields, gene_rows)
    freeze_tree(
        prior,
        {
            "artifact_class": "gse268273_outcome_free_model_input",
            "labels_included": False,
            "rna_values_present": False,
            "status": "passed",
        },
    )
    freeze_tree(
        plan,
        {
            "artifact_class": "gse268273_raw_plan",
            "participants": participants,
            "status": "planned_not_submitted",
        },
    )
    freeze_tree(
        reference,
        {
            "artifact_class": "gse268273_rsem_reference",
            "target_stable_genes": genes,
            "status": "passed",
        },
    )
    freeze_tree(
        incident,
        {
            "artifact_class": "gse268273_evaluator_receipt_incident",
            "clean_or_sealed_champion_eligible": incident_eligible,
            "status": "contained_fail_closed",
        },
    )

    source_fields = [*prior_fields, "quantification_measurement"]
    source_rows = [
        {
            **row,
            "rna_observation_state": "observed",
            "quantification_measurement": "RSEM_expected_count_raw_count_scale",
        }
        for row in prior_rows
    ]
    if forbidden_column:
        source_fields.append("fibrosis_score")
        for row in source_rows:
            row["fibrosis_score"] = "withheld"
    if reverse_participants:
        source_rows.reverse()
    write_tsv(source / "participant_axis.tsv", source_fields, source_rows)
    write_tsv(source / "candidate_gene_axis.tsv", gene_fields, gene_rows)
    counts = np.arange(1, participants * genes + 1, dtype=np.float64).reshape(
        participants, genes
    )
    if count_problem == "nan":
        counts[1, 1] = np.nan
    elif count_problem == "negative":
        counts[1, 1] = -1.0
    elif count_problem == "zero_library":
        counts[1, :] = 0.0
    np.save(source / "rna_expected_counts.npy", counts, allow_pickle=False)
    np.save(
        source / "rna_observed_mask.npy",
        np.ones(participants, dtype=np.bool_),
        allow_pickle=False,
    )
    receipt = {
        "schema_version": "masld-bench-gse268273-raw-rsem-model-input-v1",
        "status": "passed_outcome_free_quantification_training_transform_blocked",
        "participants": participants,
        "target_stable_genes": genes,
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
    }
    (source / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
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
        {
            "artifact_class": "gse268273_raw_rsem_consolidation_execution",
            "status": "passed",
        },
    )
    return {
        "consolidation_root": consolidation,
        "prior_model_input_root": prior,
        "plan_root": plan,
        "reference_root": reference,
        "incident_root": incident,
        "expected_prior_model_input_artifacts_sha256": sha256_file(
            prior / "ARTIFACTS.json"
        ),
        "expected_plan_artifacts_sha256": sha256_file(plan / "ARTIFACTS.json"),
        "expected_reference_artifacts_sha256": sha256_file(
            reference / "ARTIFACTS.json"
        ),
        "expected_incident_artifacts_sha256": sha256_file(
            incident / "ARTIFACTS.json"
        ),
        "expected_participants": participants,
        "expected_genes": genes,
        "expected_technical_runs": sum(range(1, participants + 1)),
        "expected_fastq_bytes": 100 * sum(range(1, participants + 1)),
    }


class GSE268273RawSourceIntegrityAdmissionTests(unittest.TestCase):
    def call(self, root: Path, **fixture_options: object) -> Path:
        arguments = build_fixture(root, **fixture_options)
        output = root / "admission"
        admission.run(**arguments, output=output)
        return output

    def test_admits_only_source_integrity_and_freezes_ineligible_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = self.call(Path(temporary))
            manifest = verify_frozen_tree(output)
            receipt = json.loads(
                (output / "source_integrity_admission.json").read_text()
            )
        self.assertTrue(receipt["source_integrity_only"])
        self.assertFalse(receipt["clean_or_sealed_champion_eligible"])
        self.assertTrue(receipt["independent_unexposed_rederivation_or_audit_required"])
        self.assertFalse(receipt["transform_selected_or_applied"])
        self.assertFalse(receipt["fit_or_score_performed"])
        self.assertEqual(manifest["metadata"]["status"], "passed_unexposed_audit_required")

    def test_rejects_forbidden_participant_column(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(
                admission.GSE268273SourceIntegrityError, "participant axis"
            ):
                self.call(Path(temporary), forbidden_column=True)

    def test_rejects_participant_reordering(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(
                admission.GSE268273SourceIntegrityError, "participant axis"
            ):
                self.call(Path(temporary), reverse_participants=True)

    def test_rejects_invalid_count_values(self) -> None:
        for problem in ("nan", "negative", "zero_library"):
            with self.subTest(problem=problem), tempfile.TemporaryDirectory() as temporary:
                with self.assertRaisesRegex(
                    admission.GSE268273SourceIntegrityError, "count values"
                ):
                    self.call(Path(temporary), count_problem=problem)

    def test_rejects_incident_marked_champion_eligible(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(
                admission.GSE268273SourceIntegrityError, "metadata"
            ):
                self.call(Path(temporary), incident_eligible=True)

    def test_source_has_no_evaluator_or_modeling_surface(self) -> None:
        source = (
            ROOT / "scripts/admit_gse268273_raw_rsem_source_integrity.py"
        ).read_text(encoding="utf-8")
        for forbidden in (
            "evaluator_only",
            '"outcomes.tsv"',
            "gse268273_fibrosis_ood_transfer_task",
            "sklearn",
            ".fit(",
            ".predict(",
            ".score(",
        ):
            self.assertNotIn(forbidden, source)


if __name__ == "__main__":
    unittest.main()
