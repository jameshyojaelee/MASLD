from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from scripts.gpu_dependency_dispatcher import (
    GPUDependencyDispatcherError,
    _submit,
    assert_no_incomplete_dependency_claims,
    audit_live_campaign,
    legacy_dispatcher_active,
    load_policy,
    plan_items,
    validate_item_route,
)


ROOT = Path(__file__).resolve().parents[2]
POLICY = ROOT / "config/campaigns/gpu_dependency_dispatcher_policy_20260825.json"


def scheduler_row(job_id: int, *, reason: str) -> dict[str, str]:
    return {
        "job_id": str(job_id),
        "user": "test-user",
        "state": "PENDING",
        "account": "nslab",
        "qos": "nslab",
        "partition": "gpu",
        "job_name": "model-training-700",
        "gres": "gres/gpu:l40s:1",
        "reason": reason,
        "dependency": "(null)" if reason != "Dependency" else "afterany:1",
    }


class GPUDependencyDispatcherTests(unittest.TestCase):
    def test_frozen_policy_rederives_four_bootstrap_tails(self) -> None:
        policy = load_policy(POLICY, ROOT)
        self.assertEqual(policy["lane_count"], 4)
        self.assertIsNone(policy["pending_job_limit"])
        self.assertIsNone(policy["submission_batch_limit"])
        self.assertEqual(
            [row["tail_job_id"] for row in policy["bootstrap"]["lane_tails"]],
            ["21100330", "21100328", "21100331", "21100329"],
        )

    def test_accepts_only_generic_nslab_l40s_nonarray_wrapper(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            wrapper = Path(temporary) / "job.sbatch"
            wrapper.write_text(
                "#!/bin/bash\n"
                "#SBATCH --job-name=model-training-901\n"
                "#SBATCH --partition=gpu\n"
                "#SBATCH --account=nslab\n"
                "#SBATCH --qos=nslab\n"
                "#SBATCH --gres=gpu:l40s:1\n"
                "#SBATCH --time=2-00:00:00\n",
                encoding="utf-8",
            )
            item = validate_item_route(
                {"bundle_id": "bundle-001", "_wrapper": wrapper}
            )
            self.assertEqual(item["_job_name"], "model-training-901")
            self.assertEqual(item["_hours"], 48.0)
            wrapper.write_text(
                wrapper.read_text(encoding="utf-8").replace(
                    "model-training-901", "scbasset-training"
                ),
                encoding="utf-8",
            )
            with self.assertRaises(GPUDependencyDispatcherError):
                validate_item_route({"bundle_id": "bundle-002", "_wrapper": wrapper})

    def test_unlimited_pending_plan_has_at_most_four_new_roots(self) -> None:
        lane_state = {
            "lanes": [
                {
                    "lane": lane,
                    "tail_job_id": str(100 + lane),
                    "tail_bundle_id": f"tail-{lane}",
                    "requested_gpu_hours": 0.0,
                }
                for lane in range(4)
            ]
        }
        items = [
            {"bundle_id": f"model-training-{index:03d}", "_hours": 8.0}
            for index in range(300)
        ]
        fresh = plan_items(items, lane_state, {lane: False for lane in range(4)})
        self.assertEqual(len(fresh), 300)
        self.assertEqual(sum(row["dependency_job_id"] is None for row in fresh), 4)
        chained = plan_items(items, lane_state, {lane: True for lane in range(4)})
        self.assertEqual(len(chained), 300)
        self.assertEqual(sum(row["dependency_job_id"] is None for row in chained), 0)
        self.assertEqual({row["lane"] for row in fresh}, {0, 1, 2, 3})

    def test_live_campaign_audit_allows_four_heads_and_unlimited_dependencies(self) -> None:
        rows = [scheduler_row(1000 + index, reason="(Priority)" if index < 4 else "(Dependency)") for index in range(304)]
        census = audit_live_campaign(rows, {row["job_id"] for row in rows})
        self.assertEqual(census["visible_campaign_jobs"], 304)
        self.assertEqual(census["runnable_or_running_campaign_jobs"], 4)
        self.assertEqual(census["dependency_pending_campaign_jobs"], 300)
        rows[4]["reason"] = "Priority"
        with self.assertRaises(GPUDependencyDispatcherError):
            audit_live_campaign(rows, {row["job_id"] for row in rows})

    def test_incomplete_v2_claim_fails_closed_without_rejecting_legacy_claims(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            state = Path(temporary)
            claims = state / "claims"
            claims.mkdir()
            (claims / "legacy.json").write_text(
                '{"schema_version":"masld-bench-gpu-bundle-dispatch-claim-v1",'
                '"status":"submission_failed"}\n',
                encoding="utf-8",
            )
            assert_no_incomplete_dependency_claims(state)
            (claims / "v2.json").write_text(
                '{"schema_version":"masld-bench-gpu-dependency-dispatch-claim-v1",'
                '"status":"claimed_before_submission"}\n',
                encoding="utf-8",
            )
            with self.assertRaises(GPUDependencyDispatcherError):
                assert_no_incomplete_dependency_claims(state)

    def test_live_campaign_audit_rejects_innovation_or_non_l40s(self) -> None:
        row = scheduler_row(1000, reason="Priority")
        row["qos"] = "innovation"
        with self.assertRaises(GPUDependencyDispatcherError):
            audit_live_campaign([row], {"1000"})

    def test_live_legacy_submission_outside_four_lanes_fails_closed(self) -> None:
        row = scheduler_row(1000, reason="(Priority)")
        with self.assertRaises(GPUDependencyDispatcherError):
            audit_live_campaign([row], {"1000"}, set())
        row["qos"] = "nslab"
        row["gres"] = "gres/gpu:a100:1"
        with self.assertRaises(GPUDependencyDispatcherError):
            audit_live_campaign([row], {"1000"})

    def test_terminal_legacy_dispatcher_is_not_a_scheduler_failure(self) -> None:
        terminal = subprocess.CompletedProcess(
            [],
            1,
            stdout="",
            stderr="slurm_load_jobs error: Invalid job id specified\n",
        )
        with patch(
            "scripts.gpu_dependency_dispatcher.subprocess.run",
            return_value=terminal,
        ):
            self.assertFalse(legacy_dispatcher_active("21075961"))
        scheduler_error = subprocess.CompletedProcess(
            [], 2, stdout="", stderr="scheduler unavailable\n"
        )
        with patch(
            "scripts.gpu_dependency_dispatcher.subprocess.run",
            return_value=scheduler_error,
        ):
            with self.assertRaises(GPUDependencyDispatcherError):
                legacy_dispatcher_active("21075961")

    def test_submission_is_explicit_afterany_nslab_l40s_and_nonarray(self) -> None:
        item = {
            "bundle_id": "model-training-901",
            "_job_name": "model-training-901",
            "_queue_sha256": "a" * 64,
            "_wrapper": Path("wrapper.sbatch"),
            "exports": {"BUNDLE_ID": "bundle-901"},
        }
        completed = subprocess.CompletedProcess([], 0, stdout="123456;nslab\n", stderr="")
        with patch("scripts.gpu_dependency_dispatcher.subprocess.run", return_value=completed) as run:
            self.assertEqual(_submit(item, "123455"), "123456")
        command = run.call_args.args[0]
        self.assertIn("--partition=gpu", command)
        self.assertIn("--account=nslab", command)
        self.assertIn("--qos=nslab", command)
        self.assertIn("--gres=gpu:l40s:1", command)
        self.assertIn("--dependency=afterany:123455", command)
        self.assertFalse(any(value.startswith("--array") for value in command))
        self.assertFalse(any("innovation" in value.lower() for value in command))


if __name__ == "__main__":
    unittest.main()
