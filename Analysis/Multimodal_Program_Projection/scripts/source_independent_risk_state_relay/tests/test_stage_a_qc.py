from __future__ import annotations

import importlib.util
import datetime as dt
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


QC = load("plan45_blinded_qc", "35_seal_blinded_stage_a_qc.py")
UNBLIND = load("plan45_unblinding", "37_authorize_stage_a_unblinding.py")


class StageAQCTests(unittest.TestCase):
    def test_bh_is_monotone_in_sorted_p_order(self) -> None:
        observed = QC.bh_adjust([0.01, 0.04, 0.03])
        self.assertAlmostEqual(observed[0], 0.03)
        self.assertAlmostEqual(observed[1], 0.04)
        self.assertAlmostEqual(observed[2], 0.04)

    def test_qc_operator_is_explicit(self) -> None:
        self.assertTrue(QC.metric_pass(10, "ge", 10))
        self.assertTrue(QC.metric_pass(0.05, "le", 0.1))
        self.assertFalse(QC.metric_pass(9, "ge", 10))
        with self.assertRaises(RuntimeError):
            QC.metric_pass(1, "chosen_after_qc", 1)

    def test_transfer_pair_is_frozen_from_the_actual_arm_universe(self) -> None:
        transwell_arms = {
            "SRC_CTRL", "SRC_RISK", "MOS_CTRL", "MOS_RISK", "RECIP_CTRL",
            "RECIP_RISK", "TRW_CTRL", "TRW_RISK",
        }
        self.assertEqual(
            UNBLIND.pair_map(transwell_arms)["transfer"],
            ("TRW_CTRL", "TRW_RISK"),
        )
        conditioned_arms = (
            transwell_arms - {"TRW_CTRL", "TRW_RISK"}
        ) | {"CM_CTRL", "CM_RISK"}
        self.assertEqual(
            UNBLIND.pair_map(conditioned_arms)["transfer"],
            ("CM_CTRL", "CM_RISK"),
        )

    def test_missing_risk_control_pair_fails_closed(self) -> None:
        incomplete = {
            "SRC_CTRL", "SRC_RISK", "MOS_CTRL", "MOS_RISK", "RECIP_CTRL",
            "RECIP_RISK", "TRW_CTRL",
        }
        with self.assertRaises(RuntimeError):
            UNBLIND.pair_map(incomplete)

    def test_unblinding_timestamp_requires_timezone(self) -> None:
        observed = UNBLIND.utc("2026-08-10T12:00:00+00:00")
        self.assertEqual(observed.tzinfo, dt.timezone.utc)
        with self.assertRaises(RuntimeError):
            UNBLIND.utc("2026-08-10T12:00:00")


if __name__ == "__main__":
    unittest.main()
