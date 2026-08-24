#!/usr/bin/env python3
"""Unit tests for explicit MultiVI state-to-row-sum admission."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = ROOT / "scripts/multivi_peakvi_no_outcome_probe.py"
SPEC = importlib.util.spec_from_file_location("multivi_peakvi_probe", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
module = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = module
SPEC.loader.exec_module(module)


class ExplicitStateBridgeTests(unittest.TestCase):
    def test_wrapper_does_not_precreate_script_owned_output(self) -> None:
        wrapper = (ROOT / "slurm/multivi_peakvi_run_no_outcome_probe.sbatch").read_text()
        self.assertNotIn('install -d -m 0750 "${STAGE}/probe"', wrapper)

    def test_observed_and_structural_missing_match_upstream(self) -> None:
        self.assertEqual(
            module.validate_state_row_sums(
                [3.0, 0.0, 2.0],
                ["observed", "structurally_missing", "observed"],
            ),
            [True, False, True],
        )

    def test_observed_all_zero_rejected(self) -> None:
        with self.assertRaises(module.MultiomeProbeError):
            module.validate_state_row_sums([0.0], ["observed"])

    def test_structurally_missing_nonzero_rejected(self) -> None:
        with self.assertRaises(module.MultiomeProbeError):
            module.validate_state_row_sums([1.0], ["structurally_missing"])

    def test_unadmitted_state_rejected_before_tensor_mapping(self) -> None:
        with self.assertRaises(module.MultiomeProbeError):
            module.validate_state_row_sums([0.0], ["below_qc"])

    def test_axis_and_negative_values_rejected(self) -> None:
        with self.assertRaises(module.MultiomeProbeError):
            module.validate_state_row_sums([1.0], [])
        with self.assertRaises(module.MultiomeProbeError):
            module.validate_state_row_sums([-1.0], ["observed"])


if __name__ == "__main__":
    unittest.main()
