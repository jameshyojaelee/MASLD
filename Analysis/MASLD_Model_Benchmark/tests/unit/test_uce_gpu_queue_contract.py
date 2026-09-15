from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts.gpu_bundle_dispatcher import dispatch_once, load_item


ROOT = Path(__file__).resolve().parents[2]
WRAPPER = ROOT / "slurm/run_uce_frozen_screen_50000.sbatch"
QUEUE = ROOT / "config/campaigns/gpu_bundle_queue/053-model-work-208.json"
VALIDATION = ROOT / "executions/model-work-208-queue-validation-v1/ARTIFACTS.json"


class UCEGPUQueueContractTest(unittest.TestCase):
    def test_exact_priority_53_sequential_bundle(self) -> None:
        value = json.loads(QUEUE.read_text(encoding="utf-8"))
        self.assertEqual(value["bundle_id"], "model-work-208")
        self.assertEqual(value["priority"], 53)
        self.assertTrue(value["enabled"])
        self.assertEqual(value["logical_tasks"], 2)
        self.assertEqual(value["family"], "cell_foundation")
        self.assertEqual(value["exports"]["MODEL_IDS"], "uce_4l uce_33l")
        self.assertEqual(
            value["exports"]["MASLD_GPU_DISPATCHER_AUTHORITY"],
            "sequence_gpu_dispatcher_v2",
        )
        self.assertEqual(
            value["wrapper_sha256"], sha256(WRAPPER.read_bytes()).hexdigest()
        )
        item = load_item(QUEUE, ROOT.resolve(strict=True))
        self.assertEqual(item["_ready"], VALIDATION.is_file())

    def test_wrapper_is_exact_nslab_allocation_and_generic_name(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        for expected in (
            "#SBATCH --job-name=model-work-208",
            "#SBATCH --partition=gpu",
            "#SBATCH --account=nslab",
            "#SBATCH --qos=nslab",
            "#SBATCH --gres=gpu:l40s:1",
            "#SBATCH --cpus-per-task=12",
            "#SBATCH --mem=128G",
            "#SBATCH --time=120:00:00",
        ):
            self.assertIn(expected, text)
        self.assertNotIn("innovation", text.lower())
        self.assertNotIn("#SBATCH --array", text)
        self.assertFalse(any(line.strip().startswith("sbatch ") for line in text.splitlines()))

    def test_wrapper_has_model_level_resume_and_separate_immutable_roots(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        self.assertIn('MODEL_IDS:-"uce_4l uce_33l"', text)
        self.assertIn('MODEL_IDS}" != "uce_4l uce_33l"', text)
        self.assertIn("model-work-208-${MODEL_ID}", text)
        self.assertIn("Skipping exact completed immutable model chunk", text)
        self.assertIn("verify_frozen_tree(target)", text)
        self.assertIn("publish_directory_noreplace", text)
        self.assertIn("staging-${SLURM_JOB_ID}", text)
        self.assertIn("failed-${SLURM_JOB_ID}", text)
        self.assertIn('"common_embeddings.npz"', (
            ROOT / "scripts/extract_uce_frozen_screen_50000.py"
        ).read_text(encoding="utf-8"))

    def test_bundle_is_extraction_only_without_metric_or_head_runner(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        self.assertIn("scripts/extract_uce_frozen_screen_50000.py", text)
        for forbidden in (
            "fit_predict_common_cell_heads",
            "score_cell_baselines",
            "evaluate_",
            "metric",
            "optimizer.step",
        ):
            self.assertNotIn(forbidden, text)
        self.assertIn('receipt.get("downstream_head_fit") is not False', text)

    def test_shared_gate_blocks_queue_submission_at_four_running(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            queue_dir = root / "queue"
            state = root / "state"
            queue_dir.mkdir()
            wrapper = root / "bundle.sbatch"
            wrapper.write_text(
                "#!/bin/bash\n#SBATCH --job-name=model-work-208\n"
                "#SBATCH --partition=gpu\n#SBATCH --account=nslab\n"
                "#SBATCH --qos=nslab\n#SBATCH --gres=gpu:l40s:1\ntrue\n",
                encoding="utf-8",
            )
            required = root / "ready.txt"
            required.write_text("ready\n", encoding="utf-8")
            payload = {
                "schema_version": "masld-bench-gpu-bundle-queue-item-v1",
                "bundle_id": "model-work-208",
                "priority": 53,
                "enabled": True,
                "wrapper_path": "bundle.sbatch",
                "wrapper_sha256": sha256(wrapper.read_bytes()).hexdigest(),
                "exports": {},
                "required_paths": ["ready.txt"],
                "logical_tasks": 2,
                "family": "cell_foundation",
            }
            (queue_dir / "053.json").write_text(
                json.dumps(payload) + "\n", encoding="utf-8"
            )
            with (
                patch(
                    "scripts.gpu_bundle_dispatcher.gpu_job_states",
                    return_value=["RUNNING"] * 4,
                ),
                patch("scripts.gpu_bundle_dispatcher.submit") as submit,
            ):
                jobs = dispatch_once(
                    root=root,
                    queue=queue_dir,
                    state=state,
                    user="test-user",
                    cap=5,
                    max_running=4,
                    max_pending=1,
                )
            self.assertEqual(jobs, [])
            submit.assert_not_called()
            self.assertEqual(list((state / "claims").glob("*.json")), [])


if __name__ == "__main__":
    unittest.main()
