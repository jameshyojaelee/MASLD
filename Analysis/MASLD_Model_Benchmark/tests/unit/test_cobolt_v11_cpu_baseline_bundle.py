"""Focused checks on the v11 CPU control executor.

The v11 executor is a derived copy of the reviewed v9 one, so these tests target
the parts the derivation actually changed: the single reviewed binding, the
150-run roster selected out of the frozen v11 plan, the capacity guard, and the
single-allocation shape of the run wrapper.
"""

from __future__ import annotations

import json
from pathlib import Path
import re
import unittest

import scripts.run_cobolt_v11_cpu_baseline_bundle as bundle

ROOT = Path(__file__).resolve().parents[2]
CANDIDATE = ROOT / "candidates/v1-rna-atac-cobolt-5seed-smoke-v11--e17be3dab9ba48cb"
RUN_WRAPPER = ROOT / "slurm/run_cobolt_v11_cpu_controls.sbatch"


class CoboltV11CpuBundleTests(unittest.TestCase):
    def test_single_reviewed_binding_is_self_consistent(self) -> None:
        self.assertEqual(bundle.CAMPAIGN_ID, "v1-rna-atac-cobolt-5seed-smoke-v11")
        self.assertEqual(list(bundle.CANDIDATE_BINDINGS), [bundle.CANDIDATE_SHA256])
        binding = bundle.CANDIDATE_BINDINGS[bundle.CANDIDATE_SHA256]
        self.assertEqual(binding["plan_sha256"], bundle.PLAN_SHA256)
        self.assertEqual(binding["stability_sha256"], bundle.STABILITY_SHA256)
        self.assertEqual(binding["resource_snapshot_sha256"], bundle.RESOURCE_SNAPSHOT_SHA256)
        self.assertEqual(binding["stability_mode"], "direct_candidate")
        for digest in (
            bundle.CANDIDATE_SHA256,
            bundle.PLAN_SHA256,
            bundle.STABILITY_SHA256,
            bundle.RESOURCE_SNAPSHOT_SHA256,
            bundle.SOURCE_LOCK_SHA256,
        ):
            self.assertRegex(digest, r"^[0-9a-f]{64}$")

    def test_binding_matches_the_frozen_v11_plan(self) -> None:
        plan = json.loads((CANDIDATE / "plan.json").read_text(encoding="utf-8"))
        self.assertEqual(plan["plan_sha256"], bundle.PLAN_SHA256)
        self.assertEqual(plan["source_lock_sha256"], bundle.SOURCE_LOCK_SHA256)
        self.assertEqual(
            plan["resource_firewall"]["snapshot_sha256"], bundle.RESOURCE_SNAPSHOT_SHA256
        )
        self.assertEqual(plan["campaign"]["campaign_id"], bundle.CAMPAIGN_ID)

    def test_selects_exactly_the_150_run_cpu_control_roster(self) -> None:
        plan = json.loads((CANDIDATE / "plan.json").read_text(encoding="utf-8"))
        runs = bundle.select_cpu_runs(plan, CANDIDATE)
        self.assertEqual(len(runs), bundle.EXPECTED_RUNS)
        self.assertEqual(
            sorted({str(run["model_id"]) for run in runs}), sorted(bundle.EXPECTED_MODELS)
        )
        self.assertEqual(
            sorted({int(run["seed"]) for run in runs}), sorted(bundle.EXPECTED_SEEDS)
        )
        self.assertEqual(
            sorted({int(run["fold"]) for run in runs}), sorted(bundle.EXPECTED_FOLDS)
        )
        for run in runs:
            self.assertEqual(run["resource_profile"], bundle.EXPECTED_PROFILE)
            self.assertEqual(run["task_id"], bundle.EXPECTED_TASK)
        # The 25 GPU fits belong to the B6K bundle and must never appear here.
        self.assertNotIn("cobolt", {str(run["model_id"]) for run in runs})

    def test_capacity_guard_is_bounded_by_cpu_and_memory(self) -> None:
        bundle._check_capacity(workers=12, allocated_cpus=16, allocated_memory_gb=200)
        with self.assertRaises(bundle.CpuBundleError):
            bundle._check_capacity(workers=0, allocated_cpus=16, allocated_memory_gb=200)
        with self.assertRaises(bundle.CpuBundleError):
            bundle._check_capacity(
                workers=bundle.MAX_WORKERS + 1, allocated_cpus=64, allocated_memory_gb=800
            )
        with self.assertRaises(bundle.CpuBundleError):
            bundle._check_capacity(workers=12, allocated_cpus=4, allocated_memory_gb=200)
        with self.assertRaises(bundle.CpuBundleError):
            bundle._check_capacity(workers=12, allocated_cpus=16, allocated_memory_gb=16)

    def test_run_wrapper_is_one_allocation_without_nested_sbatch(self) -> None:
        text = RUN_WRAPPER.read_text(encoding="utf-8")
        body = "\n".join(line for line in text.splitlines() if not line.startswith("#"))
        self.assertIsNone(re.search(r"(^|[;&|]\s*)sbatch(?:\s|$)", body))
        self.assertNotIn("--array", text)
        self.assertIn("--partition=cpu", text)
        self.assertIn(bundle.CANDIDATE_SHA256, text)
        self.assertIn(bundle.STABILITY_SHA256, text)


if __name__ == "__main__":
    unittest.main()
