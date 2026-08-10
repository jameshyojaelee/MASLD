from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_ROOT))
SPEC = importlib.util.spec_from_file_location(
    "plan45_stage_a_execution", SCRIPT_ROOT / "31_freeze_stage_a_execution.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class StageAExecutionTests(unittest.TestCase):
    def test_time_normalization(self) -> None:
        self.assertEqual(MODULE.hours("2", "days"), 48.0)
        self.assertEqual(MODULE.hours("24", "hours"), 24.0)
        with self.assertRaises(RuntimeError):
            MODULE.hours("1", "weeks")

    def test_chronic_schedule_cannot_be_replaced_by_acute_loading(self) -> None:
        MODULE.validate_chronic_schedule(
            {"early_cis": 24.0, "intermediate_mediator": 120.0, "late_relay": 336.0},
            336.0,
        )
        with self.assertRaises(RuntimeError):
            MODULE.validate_chronic_schedule(
                {"early_cis": 6.0, "intermediate_mediator": 24.0, "late_relay": 48.0},
                48.0,
            )

    def test_transfer_route_is_explicit(self) -> None:
        self.assertEqual(MODULE.transfer_arms("transwell"), ("TRW_CTRL", "TRW_RISK"))
        self.assertEqual(
            MODULE.transfer_arms("conditioned_medium"), ("CM_CTRL", "CM_RISK")
        )
        with self.assertRaises(RuntimeError):
            MODULE.transfer_arms("both_after_looking")

    def test_risk_control_must_share_batch_within_biological_unit(self) -> None:
        base = {
            "target_uid": "target1", "guide_id": "guide1",
            "biological_unit_id": "unit1", "assay_id": "rna",
            "time_role": "late_relay", "challenge_role": "primary_chronic",
        }
        good = [
            {**base, "arm_id": "MOS_CTRL", "batch_id": "batch1"},
            {**base, "arm_id": "MOS_RISK", "batch_id": "batch1"},
        ]
        MODULE.validate_blocked_pairing(good, "transwell")
        bad = [
            {**base, "arm_id": "MOS_CTRL", "batch_id": "batch1"},
            {**base, "arm_id": "MOS_RISK", "batch_id": "batch2"},
        ]
        with self.assertRaises(RuntimeError):
            MODULE.validate_blocked_pairing(bad, "transwell")

    def test_contract_placeholders_resolve_without_adding_arms(self) -> None:
        active = {"MOS_CTRL", "MOS_RISK", "TRW_CTRL", "TRW_RISK"}
        observed = MODULE.resolve_contract_values(
            "MOS_CTRL;{TRANSFER_RISK}", active, {"early_cis", "late_relay"},
            "transwell",
        )
        self.assertEqual(observed, {"MOS_CTRL", "TRW_RISK"})

    def test_expected_design_uses_biological_units_not_technical_replicates(self) -> None:
        targets_guides = {("target1", "guide1"), ("target1", "guide2")}
        units = {
            ("unit1", "background1", "diff1"),
            ("unit2", "background2", "diff1"),
        }
        assays = [{
            "assay_id": "rna",
            "required_arm_ids": "MOS_CTRL;MOS_RISK",
            "required_time_roles": "early_cis;late_relay",
            "required_challenge_roles": "basal;primary_chronic",
        }]
        arm_times = {
            "MOS_CTRL": {"early_cis", "late_relay"},
            "MOS_RISK": {"early_cis", "late_relay"},
        }
        keys = MODULE.expected_design_keys(targets_guides, units, assays, arm_times)
        self.assertEqual(len(keys), 2 * 2 * 2 * 2 * 2)

    def test_arm_time_intersection_prevents_invalid_transfer_early_cell(self) -> None:
        assays = [{
            "assay_id": "rna",
            "required_arm_ids": "TRW_CTRL;TRW_RISK",
            "required_time_roles": "early_cis;late_relay",
            "required_challenge_roles": "basal",
        }]
        keys = MODULE.expected_design_keys(
            {("target1", "guide1")},
            {("unit1", "background1", "diff1")},
            assays,
            {"TRW_CTRL": {"intermediate_mediator", "late_relay"},
             "TRW_RISK": {"intermediate_mediator", "late_relay"}},
        )
        self.assertEqual({key[6] for key in keys}, {"late_relay"})


if __name__ == "__main__":
    unittest.main()
