from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_ROOT))


def load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_ROOT / filename)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ADJ = load("plan45_stage_a_terminal", "41_adjudicate_stage_a_and_freeze_stage_b.py")


class StageATerminalTests(unittest.TestCase):
    def test_bh_tie_order_is_deterministic(self) -> None:
        observed = ADJ.bh_adjust([0.01, 0.01, 0.9])
        self.assertAlmostEqual(observed[0], 0.015)
        self.assertAlmostEqual(observed[1], 0.015)
        self.assertAlmostEqual(observed[2], 0.9)

    def test_gate_requires_direction_backgrounds_and_q(self) -> None:
        row = {
            "result_status": "valid_estimate",
            "q_value": "0.01",
            "expected_direction_agreement": "true",
            "n_backgrounds": "2",
        }
        self.assertTrue(ADJ.row_pass(row, 0.05, 2))
        for field, value in [
            ("q_value", "0.05"),
            ("expected_direction_agreement", "false"),
            ("n_backgrounds", "1"),
            ("result_status", "untestable_qc"),
        ]:
            altered = dict(row)
            altered[field] = value
            self.assertFalse(ADJ.row_pass(altered, 0.05, 2))

    def test_naive_timestamp_fails_closed(self) -> None:
        with self.assertRaises(RuntimeError):
            ADJ.utc("2026-08-10T10:00:00")


if __name__ == "__main__":
    unittest.main()
