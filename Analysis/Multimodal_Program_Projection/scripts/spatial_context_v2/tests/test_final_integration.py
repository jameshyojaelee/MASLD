#!/usr/bin/env python3

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from contract_lib import ContractError  # noqa: E402
from final_integration_lib import (  # noqa: E402
    HOTFIX_EXPECTED_SHA256,
    INTEGRATED_EFFECT_COLUMNS,
    independent_bh,
    make_verdicts,
    native_spatial_rows,
    require_complete_groups,
    require_no_cross_assay_construct,
    sha256_file,
    verify_hotfix,
)


def write_tsv(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class FinalIntegrationTests(unittest.TestCase):
    def test_complete_family_rejects_missing_and_duplicate(self):
        rows = [
            {"dataset": "d", "assay": "a", "analysis_set_id": "x", "program_uid": "p1"},
            {"dataset": "d", "assay": "a", "analysis_set_id": "x", "program_uid": "p2"},
        ]
        require_complete_groups(rows, {"p1", "p2"})
        with self.assertRaises(ContractError):
            require_complete_groups(rows[:1], {"p1", "p2"})
        with self.assertRaises(ContractError):
            require_complete_groups(rows + [rows[0]], {"p1", "p2"})

    def test_cross_assay_score_and_comparability_fail_closed(self):
        row = {column: "" for column in INTEGRATED_EFFECT_COLUMNS}
        row.update({"cross_assay_comparable": "FALSE", "effect_unit": "native_unit"})
        require_no_cross_assay_construct([row], INTEGRATED_EFFECT_COLUMNS)
        require_no_cross_assay_construct([dict(row, cross_assay_comparable=False)], INTEGRATED_EFFECT_COLUMNS)
        bad = dict(row, cross_assay_comparable="TRUE")
        with self.assertRaises(ContractError):
            require_no_cross_assay_construct([bad], INTEGRATED_EFFECT_COLUMNS)
        with self.assertRaises(ContractError):
            require_no_cross_assay_construct([dict(row, cross_assay_comparable=True)], INTEGRATED_EFFECT_COLUMNS)
        with self.assertRaises(ContractError):
            require_no_cross_assay_construct([row], INTEGRATED_EFFECT_COLUMNS + ("overall_rank",))

    def test_bh_complete_family(self):
        observed = independent_bh({"p1": 0.01, "p2": 0.04})
        self.assertAlmostEqual(observed["p1"], 0.02)
        self.assertAlmostEqual(observed["p2"], 0.04)

    def test_verdict_distinguishes_indeterminate_valid_null_and_terminal_boundaries(self):
        template = {
            "assay": "a",
            "analysis_set_id": "x",
            "source_dependence": "independent",
            "biological_unit": "donor",
            "n_biological": "4",
            "technical_unit": "section",
            "n_technical": "5",
            "effect_unit": "native",
            "evidence_role": "spatial_organization",
        }
        null = [dict(template, dataset="D", evidence_state="tested_negative") for _ in range(2)]
        indeterminate = [dict(template, dataset="I", evidence_state="indeterminate") for _ in range(2)]
        skip = [dict(template, dataset="S", evidence_state="skipped") for _ in range(2)]
        verdicts = {row["dataset"]: row for row in make_verdicts(null + indeterminate + skip)}
        self.assertEqual(verdicts["D"]["dataset_gate"], "valid_null")
        self.assertEqual(verdicts["I"]["dataset_gate"], "indeterminate")
        self.assertEqual(verdicts["D"]["figure4_verdict"], "retain_complete_program_panel")
        self.assertEqual(verdicts["S"]["dataset_gate"], "skipped")
        self.assertEqual(verdicts["S"]["figure4_verdict"], "retain_transparent_skip_rows")

    def test_native_spatial_adapter_uses_centered_moran_and_null_sd(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Path(tmp)
            hotspot = project / "candidate/hotspot"
            native = project / "candidate/native"
            producer = project / "Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/05b_run_visium_v2_semantic_ready.py"
            producer.parent.mkdir(parents=True)
            producer.write_text("# fixture producer\n", encoding="utf-8")
            registry_fields = (
                "cell_type", "module", "program_uid", "membership_sha256", "module_name",
                "primary_direction", "external_test_eligible",
            )
            registry = [
                {"cell_type": "hepatocytes", "module": str(i), "program_uid": f"p{i}", "membership_sha256": f"h{i}", "module_name": f"P{i}", "primary_direction": "positive", "external_test_eligible": "TRUE"}
                for i in (1, 2)
            ]
            write_tsv(hotspot / "program_registry_v2.tsv", registry_fields, registry)
            write_tsv(native / "READY", ("status",), [{"status": "pass_v2_candidate"}])
            write_tsv(
                native / "source_manifest.tsv",
                ("relative_path", "source_role", "bytes", "sha256"),
                [
                    {
                        "relative_path": producer.relative_to(project).as_posix(),
                        "source_role": "candidate_producer",
                        "bytes": producer.stat().st_size,
                        "sha256": sha256_file(producer),
                    }
                ],
            )
            write_tsv(
                native / "execution_manifest.tsv",
                ("parameter", "value"),
                [{"parameter": "candidate_producer_sha256", "value": sha256_file(producer)}],
            )
            universe_fields = ("program_id", "legacy_program_id", "cell_type", "module", "program_name", "membership_sha256", "primary_direction")
            universe = [
                {"program_id": f"p{i}", "legacy_program_id": f"hepatocytes::{i}", "cell_type": "hepatocytes", "module": str(i), "program_name": f"P{i}", "membership_sha256": f"h{i}", "primary_direction": "positive"}
                for i in (1, 2)
            ]
            write_tsv(native / "tested_universe.tsv", universe_fields, universe)
            design_fields = ("dataset", "biological_unit", "n_biological", "technical_unit", "n_technical", "design_status")
            write_tsv(
                native / "native_design_audit.tsv",
                design_fields,
                [
                    {"dataset": "GSE192741", "biological_unit": "donor", "n_biological": "4", "technical_unit": "section", "n_technical": "5", "design_status": "pass"},
                    {"dataset": "Vu_et_al_2025", "biological_unit": "physical_array", "n_biological": "10", "technical_unit": "physical_array", "n_technical": "10", "design_status": "pass"},
                ],
            )
            result_fields = (
                "program_id", "dataset", "program_name", "n_measured", "retained_l1_weight", "testable",
                "residual_moran_i", "residual_null_mean", "residual_null_sd", "residual_moran_z",
                "residual_pvalue", "residual_padj", "sensitivity_sign_agree", "robust",
            )
            results = []
            for dataset in ("GSE192741", "Vu_et_al_2025"):
                for i in (1, 2):
                    results.append(
                        {"program_id": f"p{i}", "dataset": dataset, "program_name": f"P{i}", "n_measured": "10", "retained_l1_weight": "0.5", "testable": "TRUE", "residual_moran_i": str(0.3 + i / 10), "residual_null_mean": "0.1", "residual_null_sd": "0.05", "residual_moran_z": str(4 + i), "residual_pvalue": str(0.01 * i), "residual_padj": "0.02", "sensitivity_sign_agree": "TRUE", "robust": "TRUE"}
                    )
            write_tsv(native / "spatial_program_results.tsv", result_fields, results)
            rows = native_spatial_rows(project, native, hotspot)
            self.assertEqual(len(rows), 4)
            self.assertAlmostEqual(float(rows[0]["estimate"]), 0.3)
            self.assertEqual(rows[0]["std_error"], "")
            self.assertAlmostEqual(float(rows[0]["matched_null_sd"]), 0.05)
            self.assertEqual(rows[0]["uncertainty_semantics"], "matched_gene_null_standard_deviation_not_sampling_standard_error")
            by_dataset = {row["dataset"]: row["source_dependence"] for row in rows}
            self.assertEqual(by_dataset["GSE192741"], "independent")
            self.assertEqual(by_dataset["Vu_et_al_2025"], "source_dependent")
            vu = next(row for row in rows if row["dataset"] == "Vu_et_al_2025")
            self.assertEqual(vu["n_biological"], "")
            self.assertEqual(vu["n_technical"], "10")

    def test_real_hotfix_is_pinned_and_mutation_fails(self):
        project = SCRIPT_ROOT.parents[3]
        hotfix = project / "Analysis/Spatial/candidates/program-context-v2-candidate-2026-08-07/yakubovsky2026/validator_hotfix_manifest.tsv"
        self.assertEqual(sha256_file(hotfix), HOTFIX_EXPECTED_SHA256)
        verify_hotfix(hotfix)
        with tempfile.TemporaryDirectory() as tmp:
            mutated = Path(tmp) / "validator_hotfix_manifest.tsv"
            mutated.write_bytes(hotfix.read_bytes() + b"\n")
            with self.assertRaises(ContractError):
                verify_hotfix(mutated)


if __name__ == "__main__":
    unittest.main()
