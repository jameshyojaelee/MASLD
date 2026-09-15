from __future__ import annotations

import json
from pathlib import Path
import unittest

from scripts.prime_all_ready_gpu_queue import load_selected, plan_chains


ROOT = Path(__file__).resolve().parents[2]
QUEUE = ROOT / "config/campaigns/gpu_bundle_queue"
STATE = ROOT / "executions/gpu-bundle-dispatch-state"
CONFIG = ROOT / "config/campaigns/gpu_queue_prime_all_ready_20260825.json"


class PrimeAllReadyGPUQueueTests(unittest.TestCase):
    def test_frozen_live_snapshot_has_exact_resource_census(self) -> None:
        config = json.loads(CONFIG.read_text(encoding="utf-8"))
        items = load_selected(ROOT, QUEUE, STATE, config)
        self.assertEqual(len(items), 40)
        self.assertEqual(sum(item["logical_tasks"] for item in items), 604)
        self.assertEqual(round(sum(item["_hours"] for item in items)), 1658)

    def test_four_dependency_chains_have_only_three_new_roots(self) -> None:
        items = [{"bundle_id": f"bundle-{index}", "_hours": 48.0} for index in range(20)]
        plan = plan_chains(items, "21083052")
        self.assertEqual(len(plan), len(items))
        self.assertEqual(sum(row["dependency_job_id"] is None for row in plan), 3)
        self.assertEqual({row["lane"] for row in plan}, {0, 1, 2, 3})
        self.assertEqual(len({row["bundle_id"] for row in plan}), len(items))


if __name__ == "__main__":
    unittest.main()
