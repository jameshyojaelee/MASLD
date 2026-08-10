#!/usr/bin/env python3
"""Outcome-free synthetic tests for the frozen phase-C statistical engine."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from phase_c_stats import (  # noqa: E402
    add_frozen_strata,
    bh_adjust,
    fit_class_binary,
    fit_class_continuous,
    matched_program_test,
    mechanical_figure_verdict,
    rank_first_bins,
    weighted_score,
)


def synthetic_class_frame() -> pd.DataFrame:
    rng = np.random.default_rng(991)
    rows = []
    classes = ("neither", "genetic_only", "disease_state_only", "convergent")
    effects = dict(
        neither=0.0, genetic_only=0.1, disease_state_only=0.2, convergent=0.3
    )
    for class_index, class_name in enumerate(classes):
        for index in range(40):
            guide = 4 + (index % 5)
            symbol = f"G{class_index}_{index:03d}"
            expression = 1.0 + rng.uniform(0, 5)
            chronos = rng.normal(-0.2, 0.4)
            outcome = effects[class_name] + 0.04 * chronos + rng.normal(0, 0.2)
            rows.append(
                {
                    "gene_symbol": symbol,
                    "primary_evidence_class": class_name,
                    "guide_count_min": guide,
                    "log1p_HLF_TPM": expression,
                    "HLF_Chronos": chronos,
                    "HLF_Chronos_missing": 0,
                    "outcome": outcome,
                    "event": (index + class_index) % 9 == 0,
                }
            )
    return add_frozen_strata(pd.DataFrame(rows))


class PhaseCStatsTest(unittest.TestCase):
    def test_bh_and_rank_bins(self) -> None:
        self.assertTrue(np.allclose(bh_adjust([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04]))
        bins = rank_first_bins([2, 1, 1, 3, 4], ["z", "b", "a", "c", "d"])
        self.assertEqual(bins.tolist(), [3, 2, 1, 4, 5])
        with self.assertRaises(ValueError):
            bh_adjust([])

    def test_frozen_strata_survive_subsetting(self) -> None:
        frame = synthetic_class_frame()
        before = frame.set_index("gene_symbol")["expression_quintile"].to_dict()
        subset = add_frozen_strata(frame[~frame["gene_symbol"].str.endswith("0")])
        after = subset.set_index("gene_symbol")["expression_quintile"].to_dict()
        self.assertTrue(all(after[symbol] == before[symbol] for symbol in after))

    def test_class_continuous_is_deterministic_and_complete(self) -> None:
        frame = synthetic_class_frame()
        rows_a, audit_a = fit_class_continuous(
            frame, "outcome", "fixture", 1234, accepted_draws=199, batch_size=31
        )
        rows_b, audit_b = fit_class_continuous(
            frame, "outcome", "fixture", 1234, accepted_draws=199, batch_size=31
        )
        self.assertEqual(rows_a, rows_b)
        self.assertEqual(audit_a, audit_b)
        self.assertEqual(len(rows_a), 4)
        self.assertEqual(audit_a["accepted_draws"], 199)
        self.assertTrue(all(0 < float(row["permutation_p"]) <= 1 for row in rows_a))

    def test_binary_stratified_null_is_deterministic(self) -> None:
        frame = synthetic_class_frame()
        rows_a, audits_a = fit_class_binary(frame, "event", "fixture", 4321, 257)
        rows_b, audits_b = fit_class_binary(frame, "event", "fixture", 4321, 257)
        self.assertEqual(rows_a, rows_b)
        self.assertEqual(audits_a, audits_b)
        self.assertEqual(len(rows_a), 3)
        self.assertTrue(all(int(row["accepted_draws"]) == 257 for row in rows_a))

    def test_weighted_and_matched_program_score(self) -> None:
        self.assertAlmostEqual(weighted_score([1.0, 3.0], [1.0, 3.0]), 2.5)
        rows = []
        for index in range(48):
            rows.append(
                {
                    "gene_symbol": f"P{index}" if index < 8 else f"C{index}",
                    "primary_evidence_class": "neither",
                    "guide_count_min": 4,
                    "log1p_HLF_TPM": 2.0,
                    "HLF_Chronos": -0.2,
                    "HLF_Chronos_missing": 0,
                    "guide_bin": 4,
                    "expression_quintile": 1,
                    "chronos_missing": 0,
                    "chronos_quintile": 1,
                    "permutation_stratum": "g4_e1_m0",
                    "outcome_z": 1.0 if index < 8 else (index - 28) / 20,
                }
            )
        frame = pd.DataFrame(rows)
        members = [
            {"gene_symbol": f"P{index}", "original_l1_weight": 0.03}
            for index in range(8)
        ]
        result, audit = matched_program_test(
            frame,
            members,
            {f"P{index}" for index in range(8)},
            1,
            "program_fixture",
            "primary",
            "original_weight",
            np.random.SeedSequence(123),
            draws=101,
        )
        self.assertEqual(result["testable"], "TRUE")
        self.assertEqual(result["draws"], 101)
        self.assertEqual(audit["minimum_control_pool_per_cell"], 40)
        self.assertGreater(float(result["signed_effect"]), 0)

    def test_zero_discovery_is_valid_figure_verdict(self) -> None:
        class_rows = [
            {
                "analysis_variant": "primary",
                "contrast": "omnibus_evidence_class",
                "permutation_p": 0.7,
            },
            {
                "analysis_variant": "remove_all_known_hits",
                "contrast": "omnibus_evidence_class",
                "permutation_p": 0.8,
            },
        ]
        for variant in ("primary", "remove_all_known_hits"):
            for contrast in (
                "genetic_only_vs_neither",
                "disease_state_only_vs_neither",
                "convergent_vs_neither",
            ):
                class_rows.append(
                    {
                        "analysis_variant": variant,
                        "contrast": contrast,
                        "BH_q": 0.8,
                        "estimate": 0.0,
                        "standardized_effect": 0.0,
                    }
                )
        program_rows = []
        for variant, mode in (
            ("primary", "original_weight"),
            ("primary", "equal_weight"),
            ("primary", "leave_highest_weight_out"),
            ("remove_all_known_hits", "original_weight"),
        ):
            program_rows.append(
                {
                    "analysis_variant": variant,
                    "weight_mode": mode,
                    "program_uid": "p1",
                    "testable": "TRUE",
                    "BH_q": 0.9,
                    "signed_effect": 0.0,
                }
            )
        verdict = mechanical_figure_verdict(True, class_rows, program_rows)
        self.assertEqual(verdict["main_figure_eligible"], "FALSE")
        self.assertIn("class_branch_failed", verdict["verdict_reason_codes"])
        self.assertIn("program_branch_failed", verdict["verdict_reason_codes"])


if __name__ == "__main__":
    unittest.main()
