from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
AUDIT_SOURCE = ROOT / "scripts/audit_sequence_task_native_five_seed_live_dispatch.py"
DISPATCH_SOURCE = ROOT / "scripts/gpu_bundle_dispatcher.py"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


AUDIT = load(AUDIT_SOURCE, "five_seed_live_dispatch_audit")
DISPATCH = load(DISPATCH_SOURCE, "five_seed_dispatcher_policy")


class FiveSeedLiveDispatchTest(unittest.TestCase):
    def test_contract_has_distinct_v1_and_v2_dispatcher_identities(self) -> None:
        contract = AUDIT.load_contract(
            ROOT / "config/campaigns/sequence_task_native_five_seed_live_dispatch_reconciliation_20260825.json"
        )
        v1 = contract["queue"]["v1_dispatcher_ids"]
        v2 = contract["queue"]["v2_dispatcher_ids"]
        self.assertEqual(len(v1), 30)
        self.assertEqual(len(v2), 30)
        self.assertFalse(set(v1) & set(v2))
        self.assertFalse(contract["gpu_policy"]["manual_gpu_submission_allowed"])

    def test_corrected_controller_import_and_terminal_paths_execute(self) -> None:
        contract = AUDIT.load_contract(
            ROOT / "config/campaigns/sequence_task_native_five_seed_live_dispatch_reconciliation_20260825.json"
        )
        result = AUDIT.validate_controller(ROOT, contract)
        self.assertTrue(result["counter_import_complete"])
        self.assertTrue(result["pass_terminal_summary_executed"])
        self.assertTrue(result["failure_terminal_summary_executed"])

    def test_dispatcher_cap_allows_no_second_pending_job(self) -> None:
        self.assertTrue(DISPATCH.submission_allowed(
            ["RUNNING", "RUNNING", "RUNNING"],
            cap=5,
            max_running=4,
            max_pending=1,
        ))
        self.assertFalse(DISPATCH.submission_allowed(
            ["PENDING"],
            cap=5,
            max_running=4,
            max_pending=1,
        ))
        self.assertFalse(DISPATCH.submission_allowed(
            ["RUNNING", "RUNNING", "RUNNING", "RUNNING"],
            cap=5,
            max_running=4,
            max_pending=1,
        ))

    def test_v1_queue_is_disabled_and_has_no_claim(self) -> None:
        contract = AUDIT.load_contract(
            ROOT / "config/campaigns/sequence_task_native_five_seed_live_dispatch_reconciliation_20260825.json"
        )
        queue = ROOT / "config/campaigns/gpu_bundle_queue"
        claims = ROOT / "executions/gpu-bundle-dispatch-state/claims"
        rows = []
        for path in queue.glob("*.json"):
            value = json.loads(path.read_text(encoding="utf-8"))
            if value.get("bundle_id") in contract["queue"]["v1_dispatcher_ids"]:
                rows.append(value)
        self.assertEqual(len(rows), 30)
        self.assertTrue(all(row["enabled"] is False for row in rows))
        self.assertTrue(all(not (claims / f"{row['bundle_id']}.json").exists() for row in rows))


if __name__ == "__main__":
    unittest.main()
