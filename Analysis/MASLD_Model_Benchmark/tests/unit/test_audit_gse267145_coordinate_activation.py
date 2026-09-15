from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.audit_gse267145_coordinate_activation import (
    CoordinateActivationError,
    activation_contract,
    chromhmm_grid_and_offset_test,
    coordinate_claims,
    coordinate_structural_evidence,
    reduce_postcondition_gap_test,
    validate_task_contract,
)


ROOT = Path(__file__).resolve().parents[2]


class GSE267145CoordinateActivationTests(unittest.TestCase):
    def test_bed_mention_and_in_bounds_do_not_resolve_semantics(self) -> None:
        observed = coordinate_claims(
            "Peaks were merged with BEDtools and all intervals were in bounds."
        )
        self.assertEqual(observed["zero_based_half_open"], [])
        self.assertEqual(observed["one_based_inclusive"], [])

    def test_only_direct_coordinate_claim_is_recognized(self) -> None:
        observed = coordinate_claims(
            "The deposited intervals use 0-based, half-open coordinates."
        )
        self.assertTrue(observed["zero_based_half_open"])
        self.assertEqual(observed["one_based_inclusive"], [])

    def test_declared_task_is_development_only_and_fail_closed(self) -> None:
        observed = validate_task_contract(
            ROOT / "config/evaluation/paired_bulk_rna_h3k27ac_task.toml",
            ROOT
            / "config/evaluation/paired_bulk_rna_h3k27ac_promotion_gate.json",
        )
        self.assertEqual(observed["task_id"], "paired_bulk_rna_h3k27ac")
        self.assertFalse(observed["champion_eligible"])
        self.assertFalse(observed["external_claim_eligible"])

    def test_activation_separates_observed_pair_and_rna_conditioned_lanes(self) -> None:
        coordinate = {
            "coordinate_semantics_resolved": False,
            "sequence_extraction_allowed": False,
        }
        participants = {
            "biological_unit": "participant",
            "participants": 99,
            "pairing": "same_sample_different_aliquot",
        }
        task = {
            "task_id": "paired_bulk_rna_h3k27ac",
            "champion_eligible": False,
            "external_claim_eligible": False,
        }
        observed = activation_contract(coordinate, participants, task)
        self.assertTrue(observed["task_specific_model_fitting_allowed"])
        self.assertFalse(observed["dataset_wide_or_sequence_activation"])
        self.assertFalse(observed["sequence_extraction_allowed"])
        self.assertEqual(observed["h3_feature_identity"], "opaque_source_region_string")
        self.assertTrue(
            observed["allowed_lanes"]["observed_pair_multiview"]
            ["cannot_support_rna_conditioned_claim"]
        )
        self.assertFalse(observed["single_cell_model_families_allowed"])
        self.assertIn("sequence_extraction", observed["blocked_actions"])
        self.assertFalse(
            observed["source_outcome_firewall"]["model_selection_or_error_selector_allowed"]
        )

    def test_task_validation_rejects_promotable_gate(self) -> None:
        task = ROOT / "config/evaluation/paired_bulk_rna_h3k27ac_task.toml"
        with tempfile.TemporaryDirectory() as temporary:
            gate = Path(temporary) / "gate.json"
            gate.write_text(
                json.dumps({"champion_eligible": True, "claim_mode": "confirmatory"}),
                encoding="utf-8",
            )
            with self.assertRaises(CoordinateActivationError):
                validate_task_contract(task, gate)


def _synthetic_rows_and_breakpoints(convention: str, n: int = 1_200):
    """Rows whose ChromHMM feature columns were written under `convention`,
    against a breakpoint set that is 0-based by construction."""
    breakpoints = {"chr1": set()}
    rows = []
    for i in range(n):
        bed_start = 10_000 + 1_000 * i          # true 0-based segment start
        bed_end = bed_start + 200 * (1 + i % 3)  # 200-bp bins, half-open
        breakpoints["chr1"].add(bed_start)
        breakpoints["chr1"].add(bed_end)
        if convention == "one_based_inclusive":
            start, end = bed_start + 1, bed_end
        else:
            start, end = bed_start, bed_end
        rows.append({"seqnames": "chr1", "start_position": str(start),
                     "end_position": str(end)})
    return rows, breakpoints


def _synthetic_regions(convention: str, n: int = 2_000):
    """Adjacent deposited regions with the gap floor each merge tool leaves."""
    regions, cursor = [], 100_000
    for i in range(n):
        width = 300 + 7 * (i % 40)
        start, end = cursor, cursor + width
        regions.append(f"chr1:{start}-{end}")
        if convention == "one_based_inclusive":
            gap = 2 + (i % 5)          # reduce() floor is 2
        else:
            gap = 1 + (i % 5)          # bedtools merge leaves gap-1 pairs
        cursor = end + gap
    return regions


class GSE267145CoordinateStructuralEvidenceTests(unittest.TestCase):
    def test_structural_tests_return_zero_based_on_a_zero_based_fixture(self) -> None:
        rows, bp = _synthetic_rows_and_breakpoints("zero_based_half_open")
        a = chromhmm_grid_and_offset_test(rows, bp)
        b = reduce_postcondition_gap_test(_synthetic_regions("zero_based_half_open"))
        self.assertEqual(a["verdict"], "zero_based_half_open")
        self.assertEqual(a["offset_minus_one_breakpoint_matches"], 0)
        self.assertEqual(b["verdict"], "zero_based_half_open")
        self.assertGreater(b["n_gap_eq_one"], 0)
        frozen = coordinate_structural_evidence(a, b)
        self.assertEqual(frozen["coordinate_convention_structural"], "zero_based_half_open")

    def test_structural_tests_return_one_based_on_a_one_based_fixture(self) -> None:
        rows, bp = _synthetic_rows_and_breakpoints("one_based_inclusive")
        a = chromhmm_grid_and_offset_test(rows, bp)
        b = reduce_postcondition_gap_test(_synthetic_regions("one_based_inclusive"))
        self.assertEqual(a["verdict"], "one_based_inclusive")
        self.assertEqual(a["offset_zero_breakpoint_matches"], 0)
        self.assertEqual(b["verdict"], "one_based_inclusive")
        self.assertEqual(b["n_gap_eq_one"], 0)
        self.assertEqual(b["min_gap"], 2)

    def test_structural_disagreement_fails_closed(self) -> None:
        rows, bp = _synthetic_rows_and_breakpoints("one_based_inclusive")
        a = chromhmm_grid_and_offset_test(rows, bp)
        b = reduce_postcondition_gap_test(_synthetic_regions("zero_based_half_open"))
        frozen = coordinate_structural_evidence(a, b)
        self.assertFalse(frozen["tests_agree"])
        self.assertEqual(frozen["coordinate_convention_structural"], "unresolved")
        self.assertEqual(frozen["status"], "structural_evidence_conflicts_fail_closed")
        self.assertEqual(frozen["resolution_class"], "none")

    def test_structural_freeze_is_never_recorded_as_a_primary_source_statement(self) -> None:
        rows, bp = _synthetic_rows_and_breakpoints("one_based_inclusive")
        frozen = coordinate_structural_evidence(
            chromhmm_grid_and_offset_test(rows, bp),
            reduce_postcondition_gap_test(_synthetic_regions("one_based_inclusive")),
            direct_claims={"workbook": coordinate_claims("merged with BEDtools, in bounds")},
        )
        self.assertEqual(frozen["status"], "frozen_by_structural_inference")
        self.assertFalse(frozen["primary_source_statement_found"])
        self.assertFalse(frozen["coordinate_semantics_resolved_by_primary_source"])
        self.assertFalse(frozen["bounds_compatibility_used_as_semantic_evidence"])
        self.assertFalse(frozen["tool_default_or_file_extension_used_as_semantic_evidence"])

    def test_primary_source_statement_overrides_structural_evidence(self) -> None:
        rows, bp = _synthetic_rows_and_breakpoints("one_based_inclusive")
        frozen = coordinate_structural_evidence(
            chromhmm_grid_and_offset_test(rows, bp),
            reduce_postcondition_gap_test(_synthetic_regions("one_based_inclusive")),
            direct_claims={"article": coordinate_claims(
                "Intervals are given in 0-based, half-open coordinates.")},
        )
        self.assertTrue(frozen["voided_by_primary_source"])
        self.assertEqual(frozen["coordinate_convention_structural"], "unresolved")
        self.assertEqual(frozen["status"], "voided_by_contradicting_primary_source_statement")


if __name__ == "__main__":
    unittest.main()
