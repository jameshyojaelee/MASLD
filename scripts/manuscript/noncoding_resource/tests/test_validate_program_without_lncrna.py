#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys
import unittest

import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import validate_program_without_lncrna as sensitivity  # noqa: E402


class ProgramSensitivityTests(unittest.TestCase):
    def fixture(self):
        content = pd.DataFrame(
            {
                "program_uid": ["p1", "p2"],
                "lncrna_l1_fraction": [0.10, 0.01],
                "requires_leave_all_lncrna_out_sensitivity": ["true", "false"],
            }
        )
        scores = pd.DataFrame(
            {
                "program_uid": ["p1"] * 4,
                "sample_id": ["s1", "s2", "s3", "s4"],
                "original_score": [1.0, 2.0, 3.0, 4.0],
                "without_lncrna_score": [1.1, 2.1, 2.9, 4.1],
            }
        )
        effects = pd.DataFrame(
            {
                "program_uid": ["p1"],
                "original_stage_beta": [0.5],
                "without_lncrna_stage_beta": [0.4],
            }
        )
        return content, scores, effects

    def test_passing_triggered_program(self):
        result = sensitivity.validate_program_sensitivity(*self.fixture())
        self.assertEqual(result["program_uid"].tolist(), ["p1"])
        self.assertTrue(bool(result.loc[0, "sensitivity_passed"]))
        self.assertTrue(bool(result.loc[0, "stage_direction_preserved"]))

    def test_direction_reversal_fails(self):
        content, scores, effects = self.fixture()
        effects.loc[0, "without_lncrna_stage_beta"] = -0.4
        result = sensitivity.validate_program_sensitivity(content, scores, effects)
        self.assertFalse(bool(result.loc[0, "sensitivity_passed"]))

    def test_trigger_drift_is_rejected(self):
        content, scores, effects = self.fixture()
        content.loc[0, "requires_leave_all_lncrna_out_sensitivity"] = "false"
        with self.assertRaisesRegex(sensitivity.ProgramSensitivityError, "5% rule"):
            sensitivity.validate_program_sensitivity(content, scores, effects)

    def test_incomplete_triggered_family_is_rejected(self):
        content, scores, effects = self.fixture()
        content.loc[1, "lncrna_l1_fraction"] = 0.08
        content.loc[1, "requires_leave_all_lncrna_out_sensitivity"] = "true"
        with self.assertRaisesRegex(
            sensitivity.ProgramSensitivityError, "every triggered"
        ):
            sensitivity.validate_program_sensitivity(content, scores, effects)


if __name__ == "__main__":
    unittest.main(verbosity=2)
