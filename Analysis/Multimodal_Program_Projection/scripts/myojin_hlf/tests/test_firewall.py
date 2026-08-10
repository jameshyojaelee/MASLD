#!/usr/bin/env python3
"""Dependency-free regression tests for the sealed-phase Myojin firewall."""

from __future__ import annotations

import csv
import importlib.util
import math
import sys
import unittest
from collections import defaultdict
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SCRIPT_DIR.parents[3]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
CANDIDATE = (
    PROJECT_ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / "program-context-v2-candidate-2026-08-07/myojin_hlf"
)
sys.path.insert(0, str(SCRIPT_DIR))

from myojin_firewall_common import read_tsv, sha256_file, split_gene_header  # noqa: E402


def load_script_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(
        module_name, SCRIPT_DIR / filename
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_validator_module():
    return load_script_module("03_validate_blind_spec.py", "myojin_validator")


class FirewallTests(unittest.TestCase):
    def test_gene_header_normalization_is_conservative(self) -> None:
        self.assertEqual(split_gene_header("ACSL3 (2181)"), ("ACSL3", "2181"))
        self.assertEqual(split_gene_header("ACSL3"), ("ACSL3", ""))
        self.assertEqual(split_gene_header("GENE-WITH-DASH (42)"), ("GENE-WITH-DASH", "42"))

    def test_masked_schema_exact_allowlist(self) -> None:
        validator = load_validator_module()
        with (CANDIDATE / "masked_screen_schema.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            observed = set(next(csv.reader(handle, delimiter="\t")))
        self.assertEqual(observed, validator.MASKED_ALLOWLIST)
        self.assertFalse(observed & validator.FORBIDDEN_VALUE_COLUMNS)

    def test_depmap_hlf_identity_fixture_excludes_hlf_a(self) -> None:
        depmap = load_script_module("01_extract_depmap_hlf.py", "myojin_depmap")
        model = depmap.read_hlf_model(FIXTURES / "Model_minimal.csv")
        self.assertEqual(model["ModelID"], "ACH-000393")
        self.assertEqual(model["CellLineName"], "HLF")
        self.assertEqual(model["OncotreeLineage"], "Liver")

    def test_depmap_wide_fixture_maps_one_default_hlf_row(self) -> None:
        depmap = load_script_module("01_extract_depmap_hlf.py", "myojin_depmap_wide")
        columns, row, profile, default = depmap.read_wide_hlf_row(
            FIXTURES / "wide_hlf_minimal.csv"
        )
        values = depmap.gene_value_map(columns, row)
        self.assertEqual(profile, "PR-HLF")
        self.assertEqual(default, "TRUE")
        self.assertEqual(values["ACSL3"], ("ACSL3 (2181)", "2181", "0.25"))
        self.assertEqual(values["GPAT3"], ("GPAT3 (84803)", "84803", "-0.75"))

    def test_depmap_duplicate_normalized_symbol_fails_closed(self) -> None:
        depmap = load_script_module("01_extract_depmap_hlf.py", "myojin_depmap_duplicate")
        with self.assertRaisesRegex(ValueError, "Duplicate normalized DepMap symbol"):
            depmap.gene_value_map(
                ["ModelID", "ACSL3 (2181)", "ACSL3 (999999)"],
                ["ACH-000393", "0.25", "0.5"],
            )

    def test_depmap_entrez_conflicts_are_retained_but_covariates_blanked(self) -> None:
        depmap = load_script_module("01_extract_depmap_hlf.py", "myojin_depmap_conflict")
        effect_header, effect_row, _, _ = depmap.read_wide_hlf_row(
            FIXTURES / "wide_effect_entrez_conflict.csv"
        )
        expression_header, expression_row, _, _ = depmap.read_wide_hlf_row(
            FIXTURES / "wide_expression_entrez_conflict.csv"
        )
        rows, audit = depmap.merge_covariate_maps(
            depmap.gene_value_map(effect_header, effect_row),
            depmap.gene_value_map(expression_header, expression_row),
        )
        by_symbol = {row["gene_symbol"]: row for row in rows}
        audit_by_symbol = {row["gene_symbol"]: row for row in audit}
        self.assertEqual(by_symbol["ACSL3"]["depmap_entrez_conflict"], "FALSE")
        self.assertNotEqual(by_symbol["ACSL3"]["HLF_TPM"], "")
        self.assertNotEqual(by_symbol["ACSL3"]["HLF_Chronos"], "")
        for symbol in ("MEF2B", "RLN2"):
            self.assertEqual(by_symbol[symbol]["depmap_entrez_conflict"], "TRUE")
            self.assertEqual(by_symbol[symbol]["depmap_mapping_state"], "entrez_conflict_values_blanked")
            self.assertEqual(by_symbol[symbol]["HLF_TPM"], "")
            self.assertEqual(by_symbol[symbol]["HLF_Chronos"], "")
            self.assertEqual(audit_by_symbol[symbol]["covariate_values_blanked"], "TRUE")
            self.assertEqual(audit_by_symbol[symbol]["eligible_for_primary_mapping"], "FALSE")

    def test_low_observed_tpm_is_complete_but_below_expression_floor(self) -> None:
        freezer = load_script_module("02_freeze_blind_spec.py", "myojin_freezer_fixture")
        with (FIXTURES / "covariate_availability.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            state = freezer.depmap_covariate_state(row)
            self.assertEqual(
                state["complete"], row["complete_numeric_covariates"] == "TRUE"
            )
            self.assertEqual(
                state["expression_floor_pass"], row["expression_floor_pass"] == "TRUE"
            )
        low = next(row for row in rows if row["case"] == "low_but_observed")
        low_state = freezer.depmap_covariate_state(low)
        self.assertTrue(low_state["complete"])
        self.assertFalse(low_state["expression_floor_pass"])

    def test_no_outcome_export_markers(self) -> None:
        masked = read_tsv(CANDIDATE / "masked_screen_schema.tsv")
        universe = read_tsv(CANDIDATE / "screen_gene_universe.tsv")
        predictions = read_tsv(CANDIDATE / "prediction_manifest.tsv")
        self.assertEqual(len(masked), 18343)
        self.assertTrue(all(row["outcome_values_exported"] == "FALSE" for row in masked))
        self.assertTrue(all(row["outcome_values_read"] == "FALSE" for row in universe))
        self.assertTrue(all(row["external_outcomes_read"] == "FALSE" for row in predictions))

    def test_workbook_is_authenticated_but_not_joined(self) -> None:
        manifest = read_tsv(CANDIDATE / "source_manifest.tsv")
        self.assertEqual(len(manifest), 1)
        self.assertEqual(
            manifest[0]["sha256"],
            "2663a811676a1abe2778375e6d28f6397eaded1f895d274f3c241733669c1c66",
        )
        self.assertEqual(manifest[0]["gene_records"], "18343")
        self.assertEqual(manifest[0]["outcome_values_exported_to_blind_phase"], "FALSE")

    def test_program_weights_and_expected_direction(self) -> None:
        rows = [
            row
            for row in read_tsv(CANDIDATE / "prediction_manifest.tsv")
            if row["object_type"] == "hepatocyte_program_gene"
        ]
        sums = defaultdict(float)
        for row in rows:
            sums[row["object_id"]] += float(row["frozen_weight"])
            self.assertEqual(
                row["expected_direction"], "positive_KO_LFC_counter_state_compatible"
            )
        self.assertEqual(len(sums), 2)
        self.assertTrue(
            all(math.isclose(value, 1.0, rel_tol=0, abs_tol=1e-12) for value in sums.values())
        )

    def test_genetics_only_direction_is_not_invented(self) -> None:
        rows = [
            row
            for row in read_tsv(CANDIDATE / "prediction_manifest.tsv")
            if row["object_type"] == "evidence_class_gene"
            and row["primary_evidence_class"] == "genetic_only"
        ]
        self.assertGreater(len(rows), 0)
        self.assertTrue(all(row["expected_direction"] == "not_specified" for row in rows))

    def test_known_hit_set_is_frozen(self) -> None:
        rows = read_tsv(CANDIDATE / "known_hit_exclusion.tsv")
        self.assertEqual(
            {row["canonical_symbol"] for row in rows},
            {"ACSL3", "INSIG1", "NF2", "CASP8", "GPAT3", "RNF213"},
        )
        self.assertTrue(all(row["included_in_primary"] == "TRUE" for row in rows))
        self.assertTrue(
            all(row["remove_all_known_hits_sensitivity"] == "TRUE" for row in rows)
        )

    def test_blocked_gate_is_fail_closed(self) -> None:
        gates = {row["gate"]: row for row in read_tsv(CANDIDATE / "gate_status.tsv")}
        if gates["version_matched_depmap_triplet"]["passed"] == "FALSE":
            self.assertTrue((CANDIDATE / "BLOCKED").is_file())
            self.assertFalse((CANDIDATE / "SEALED").exists())
            self.assertFalse((CANDIDATE / "FIREWALL_READY").exists())
            self.assertEqual(len(read_tsv(CANDIDATE / "detectability_covariates.tsv")), 0)

    def test_spec_contains_fixed_magnitude_and_draw_rules(self) -> None:
        text = (CANDIDATE / "blind_analysis_spec.yaml").read_text(encoding="utf-8")
        for required in (
            "accepted_draws: 100000",
            "null_draws: 10000",
            "class_minimum_standardized_pairwise_effect: 0.20",
            "program_minimum_signed_standardized_effect: 0.20",
            "covariate_coverage_definition: finite_numeric_TPM_and_Chronos_before_expression_floor",
            "expression_floor_role: primary_eligibility_and_testability_only_not_missingness",
            "coefficient_confidence_interval: HC3_robust_Wald_95pct",
            "standardized_pairwise_effect: adjusted_class_coefficient_divided_by_full_model_residual_SD",
            "draw_rule: one_control_per_program_member_without_replacement_within_draw",
            "matched_control_weight: inherit_the_corresponding_program_member_frozen_weight",
            "use_goodsgrna_rank_p_or_FDR_as_detectability_covariates",
            "reconstruct_vehicle_fitness_by_subtracting_MAGeCK_contrasts",
        ):
            self.assertIn(required, text)

    def test_raw_workbook_permissions_are_read_only(self) -> None:
        raw = CANDIDATE / "source/raw/Myojin_Supplemental_Table_S1_C218.xlsx"
        self.assertEqual(sha256_file(raw), "2663a811676a1abe2778375e6d28f6397eaded1f895d274f3c241733669c1c66")
        self.assertEqual(raw.stat().st_mode & 0o222, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
