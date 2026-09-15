from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
import tempfile
import unittest

import scripts.run_cobolt_v9_cpu_baseline_bundle as bundle


ROOT = Path(__file__).resolve().parents[2]
CANDIDATE = (
    ROOT
    / "candidates"
    / "v1-rna-atac-cobolt-5seed-smoke-v9--7031e0dfbf76bc4d"
)
CURRENT_CANDIDATE = (
    ROOT
    / "candidates"
    / "v1-rna-atac-cobolt-5seed-smoke-v9--edf243d20e2b8f58"
)
RUN_WRAPPER = ROOT / "slurm/run_cobolt_v9_cpu_baseline_bundle.sbatch"
CURRENT_RUN_WRAPPER = ROOT / "slurm/run_cobolt_v9_current_cpu_baseline_bundle.sbatch"
VALIDATE_WRAPPER = ROOT / "slurm/validate_cobolt_v9_cpu_baseline_bundle.sbatch"


class CoboltV9CpuBaselineBundleTests(unittest.TestCase):
    def test_frozen_candidate_has_exact_unsealed_cpu_roster(self) -> None:
        self.assertEqual(
            bundle.STABILITY_SHA256,
            "b0f3ba334165703d9381f4dd04b36daac645a9f68719bfdf61f729cb4bdc9799",
        )
        self.assertNotEqual(
            bundle.CANDIDATE_SHA256,
            bundle.STABILITY_PREFLIGHT_CANDIDATE_SHA256,
        )
        self.assertEqual(
            bundle.SOURCE_LOCK_SHA256,
            "810da04de9fec741ca626406d9c4a9eac6cb3b3c890d990dc8367a9523fabfe4",
        )
        plan = json.loads((CANDIDATE / "plan.json").read_text(encoding="utf-8"))
        runs = bundle.select_cpu_runs(plan, CANDIDATE)
        self.assertEqual(len(runs), 150)
        self.assertEqual(
            Counter(str(run["model_id"]) for run in runs),
            Counter({model_id: 25 for model_id in bundle.EXPECTED_MODELS}),
        )
        self.assertTrue(
            all(run["metadata"]["sealed_prediction_dataset_ids"] == [] for run in runs)
        )
        self.assertTrue(all(run["dataset_ids"] == ["gse296875"] for run in runs))

    def test_capacity_is_bounded_by_cpu_and_memory(self) -> None:
        bundle._check_capacity(workers=12, allocated_cpus=16, allocated_memory_gb=200)
        with self.assertRaisesRegex(bundle.CpuBundleError, "less memory"):
            bundle._check_capacity(
                workers=12, allocated_cpus=16, allocated_memory_gb=191
            )
        with self.assertRaisesRegex(bundle.CpuBundleError, "between 1 and 12"):
            bundle._check_capacity(
                workers=13, allocated_cpus=16, allocated_memory_gb=208
            )

    def test_current_stability_passed_candidate_has_same_cpu_roster(self) -> None:
        plan = json.loads((CURRENT_CANDIDATE / "plan.json").read_text(encoding="utf-8"))
        runs = bundle.select_cpu_runs(plan, CURRENT_CANDIDATE)
        self.assertEqual(len(runs), 150)
        self.assertEqual(
            bundle.CANDIDATE_BINDINGS[bundle.CURRENT_CANDIDATE_SHA256],
            {
                "plan_sha256": bundle.CURRENT_PLAN_SHA256,
                "resource_snapshot_sha256": bundle.CURRENT_RESOURCE_SNAPSHOT_SHA256,
                "stability_sha256": bundle.CURRENT_STABILITY_SHA256,
                "stability_mode": "direct_candidate",
            },
        )

    def test_resume_reuses_only_one_verified_success(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            execution = Path(temporary)
            succeeded = "a" * 64
            missing = "b" * 64
            attempt = execution / "runs" / succeeded / "attempt-001"
            attempt.mkdir(parents=True)
            for name in ("run_execution_receipt.json", "ARTIFACTS.json", "COMPLETE"):
                (attempt / name).write_text(name, encoding="utf-8")

            def tree_verifier(path: Path) -> dict[str, object]:
                self.assertEqual(path, attempt)
                return {}

            def attempt_verifier(
                path: Path, *, require_succeeded: bool
            ) -> dict[str, object]:
                self.assertEqual(path, attempt)
                self.assertTrue(require_succeeded)
                return {"run_id": succeeded, "status": "succeeded"}

            observed = bundle.classify_existing_successes(
                [{"run_id": succeeded}, {"run_id": missing}],
                execution,
                attempt_verifier=attempt_verifier,
                tree_verifier=tree_verifier,
            )
            self.assertEqual(set(observed), {succeeded})
            self.assertEqual(observed[succeeded]["run_status"], "succeeded")

    def test_resume_fails_closed_on_interrupted_attempt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            execution = Path(temporary)
            run_id = "c" * 64
            run_root = execution / "runs" / run_id
            run_root.mkdir(parents=True)
            (run_root / ".attempt-001.interrupted").mkdir()
            with self.assertRaisesRegex(
                bundle.CpuBundleError, "zero new runs started"
            ):
                bundle.classify_existing_successes(
                    [{"run_id": run_id}], execution
                )

    def test_slurm_surfaces_are_single_allocation_without_nested_sbatch(self) -> None:
        run_text = RUN_WRAPPER.read_text(encoding="utf-8")
        current_run_text = CURRENT_RUN_WRAPPER.read_text(encoding="utf-8")
        validate_text = VALIDATE_WRAPPER.read_text(encoding="utf-8")
        for line in (
            "#SBATCH --job-name=model-training-213",
            "#SBATCH --partition=cpu",
            "#SBATCH --account=nslab",
            "#SBATCH --qos=nslab",
            "#SBATCH --cpus-per-task=16",
            "#SBATCH --mem=200G",
            "#SBATCH --time=18:00:00",
            "#SBATCH --signal=B:USR1@3600",
            "--workers 12",
            "--allocated-cpus 16",
            "--allocated-memory-gb 200",
        ):
            self.assertIn(line, run_text)
        for line in (
            "#SBATCH --job-name=model-training-216",
            "#SBATCH --partition=cpu",
            "#SBATCH --account=nslab",
            "#SBATCH --qos=nslab",
            "#SBATCH --cpus-per-task=16",
            "#SBATCH --mem=200G",
            "#SBATCH --time=18:00:00",
            "#SBATCH --signal=B:USR1@3600",
            "--workers 12",
            "--allocated-cpus 16",
            "--allocated-memory-gb 200",
        ):
            self.assertIn(line, current_run_text)
        for text in (run_text, current_run_text, validate_text):
            self.assertNotIn("#SBATCH --array", text)
            body = "\n".join(
                line for line in text.splitlines() if not line.startswith("#")
            )
            self.assertIsNone(bundle.NESTED_SBATCH.search(body))
        source = (
            ROOT / "scripts/run_cobolt_v9_cpu_baseline_bundle.py"
        ).read_text(encoding="utf-8")
        for value in (
            "ThreadPoolExecutor(max_workers=arguments.workers)",
            "publish_directory_noreplace(staging, target)",
            "verify_run_execution_attempt",
            "fcntl.LOCK_EX | fcntl.LOCK_NB",
            '"outcome_records_parsed_or_joined": False',
        ):
            self.assertIn(value, source)


if __name__ == "__main__":
    unittest.main()
