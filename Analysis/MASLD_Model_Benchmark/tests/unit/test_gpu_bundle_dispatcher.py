from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

from scripts.gpu_bundle_dispatcher import (
    GPUDispatcherError,
    load_item,
    submission_allowed,
)


ROOT = Path(__file__).resolve().parents[2]


class GPUBundleDispatcherTests(unittest.TestCase):
    @staticmethod
    def _queue_payload(wrapper: Path, **updates: object) -> dict[str, object]:
        payload: dict[str, object] = {
            "schema_version": "masld-bench-gpu-bundle-queue-item-v1",
            "bundle_id": "bundle-001",
            "priority": 1,
            "enabled": True,
            "wrapper_path": "bundle.sbatch",
            "wrapper_sha256": sha256(wrapper.read_bytes()).hexdigest(),
            "exports": {"BUNDLE_ID": "bundle-001"},
            "required_paths": ["input.COMPLETE"],
            "logical_tasks": 8,
            "family": "sequence",
        }
        payload.update(updates)
        return payload

    def test_validates_frozen_nslab_nonarray_wrapper(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            wrapper = root / "bundle.sbatch"
            wrapper.write_text(
                "#!/bin/bash\n#SBATCH --partition=gpu\n#SBATCH --qos=nslab\n"
                "#SBATCH --gres=gpu:l40s:1\ntrue\n"
            )
            required = root / "input.COMPLETE"
            required.write_text("done\n")
            queue = root / "item.json"
            queue.write_text(json.dumps(self._queue_payload(wrapper)) + "\n")
            item = load_item(queue, root.resolve())
            self.assertTrue(item["_ready"])

    def test_rejects_required_path_outside_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            wrapper = root / "bundle.sbatch"
            wrapper.write_text(
                "#!/bin/bash\n#SBATCH --partition=gpu\n#SBATCH --qos=nslab\n"
                "#SBATCH --gres=gpu:l40s:1\ntrue\n"
            )
            queue = root / "item.json"
            queue.write_text(
                json.dumps(self._queue_payload(wrapper, required_paths=["../escape.COMPLETE"]))
                + "\n"
            )
            with self.assertRaises(GPUDispatcherError):
                load_item(queue, root.resolve())

    def test_rejects_arrays_and_innovation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            wrapper = root / "bundle.sbatch"
            wrapper.write_text(
                "#!/bin/bash\n#SBATCH --partition=gpu\n#SBATCH --qos=nslab\n"
                "#SBATCH --gres=gpu:l40s:1\n#SBATCH --array=0-4\n--qos=innovation\n"
            )
            queue = root / "item.json"
            queue.write_text(
                json.dumps(
                    self._queue_payload(
                        wrapper,
                        bundle_id="bundle-002",
                        exports={},
                        required_paths=[],
                        logical_tasks=2,
                        family="integration",
                    )
                )
                + "\n"
            )
            with self.assertRaises(GPUDispatcherError):
                load_item(queue, root.resolve())

    def test_allows_only_four_running_plus_one_pending_policy(self) -> None:
        policy = {"cap": 5, "max_running": 4, "max_pending": 1}
        self.assertFalse(submission_allowed(["RUNNING"] * 4, **policy))
        self.assertTrue(submission_allowed(["RUNNING"] * 3, **policy))
        self.assertFalse(submission_allowed(["RUNNING"] * 4 + ["PENDING"], **policy))
        self.assertFalse(submission_allowed(["RUNNING"] * 3 + ["PENDING"], **policy))
        self.assertFalse(submission_allowed(["RUNNING"] * 5, **policy))
        self.assertFalse(submission_allowed(["RUNNING"] * 3 + ["CONFIGURING"], **policy))
        self.assertTrue(submission_allowed(["RUNNING"] * 2 + ["CONFIGURING"], **policy))

    def test_rejects_policy_with_more_than_one_pending_slot(self) -> None:
        with self.assertRaises(GPUDispatcherError):
            submission_allowed([], cap=5, max_running=4, max_pending=2)

    def test_next_sequence_queue_is_generic_nslab_and_diverse(self) -> None:
        wrapper = ROOT / "slurm/run_sequence_gpu_bundle.sbatch"
        wrapper_text = wrapper.read_text()
        self.assertIn("#SBATCH --job-name=model-training-400", wrapper_text)
        self.assertIn("#SBATCH --account=nslab", wrapper_text)
        self.assertIn("#SBATCH --qos=nslab", wrapper_text)
        self.assertNotIn("#SBATCH --job-name=masld", wrapper_text.lower())
        expected_exports = [
            "sequence_gpu_bundle_016_sequence_cnn_control-ready_36823fe74c",
            "sequence_gpu_bundle_023_sequence_transformer_control-ready_ea978e9249",
            "sequence_gpu_bundle_004_bpnet-ready_71d22f6453",
            "sequence_gpu_bundle_010_chrombpnet-ready_a7f761b0ef",
        ]
        queue_paths = [
            ROOT / f"config/campaigns/gpu_bundle_queue/{number:03d}-model-training-{355 + number}.json"
            for number in range(46, 50)
        ]
        values = [json.loads(path.read_text()) for path in queue_paths]
        self.assertEqual([value["priority"] for value in values], list(range(46, 50)))
        self.assertEqual(
            [value["exports"]["BUNDLE_ID"] for value in values], expected_exports
        )
        self.assertEqual(
            {value["bundle_id"] for value in values},
            {f"model-training-{number}" for number in range(401, 405)},
        )
        self.assertEqual({value["logical_tasks"] for value in values}, {8})
        self.assertEqual({value["family"] for value in values}, {"model_training"})
        self.assertEqual(
            {value["wrapper_sha256"] for value in values},
            {sha256(wrapper.read_bytes()).hexdigest()},
        )
        self.assertTrue(
            all("executions/sequence-gpu-bundles-next-validation-21076887/ARTIFACTS.json"
                in value["required_paths"] for value in values)
        )

    def test_frozen_screen_queue_is_generic_nslab_and_fail_closed(self) -> None:
        wrapper = ROOT / "slurm/run_model_training_500.sbatch"
        wrapper_text = wrapper.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --job-name=model-training-500", wrapper_text)
        self.assertIn("#SBATCH --partition=gpu", wrapper_text)
        self.assertIn("#SBATCH --account=nslab", wrapper_text)
        self.assertIn("#SBATCH --qos=nslab", wrapper_text)
        self.assertNotIn("--qos=innovation", wrapper_text)
        self.assertNotIn("#SBATCH --array", wrapper_text)
        value = json.loads(
            (
                ROOT
                / "config/campaigns/gpu_bundle_queue/050-model-training-500.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(value["bundle_id"], "model-training-500")
        self.assertEqual(value["priority"], 50)
        self.assertEqual(value["logical_tasks"], 2)
        self.assertEqual(value["family"], "model_training")
        self.assertEqual(
            value["wrapper_sha256"], sha256(wrapper.read_bytes()).hexdigest()
        )
        self.assertEqual(
            value["required_paths"],
            [
                "executions/model-data-051-21077382/ARTIFACTS.json",
                "executions/model-check-032-21077481/ARTIFACTS.json",
            ],
        )

    def test_geneformer_one_batch_queue_is_next_and_cpu_hash_bound(self) -> None:
        wrapper = ROOT / "slurm/preflight_geneformer_v2_316m_one_batch_l40s.sbatch"
        queue = (
            ROOT
            / "config/campaigns/gpu_bundle_queue/051-model-probe-501.json"
        )
        runtime = (
            ROOT
            / "executions/geneformer-v2-316m-compat-runtime-21077474/ARTIFACTS.json"
        )
        validation = (
            ROOT
            / "executions/geneformer-v2-316m-preflight-validation-21077937/ARTIFACTS.json"
        )
        wrapper_text = wrapper.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --job-name=model-probe-043", wrapper_text)
        self.assertIn("#SBATCH --partition=gpu", wrapper_text)
        self.assertIn("#SBATCH --account=nslab", wrapper_text)
        self.assertIn("#SBATCH --qos=nslab", wrapper_text)
        self.assertIn("#SBATCH --gres=gpu:l40s:1", wrapper_text)
        self.assertNotIn("innovation", wrapper_text.lower())
        self.assertNotIn("#SBATCH --array", wrapper_text)
        self.assertEqual(
            sha256(wrapper.read_bytes()).hexdigest(),
            "75c3defae651c0331ab10591a5a0a49e3bff2c9a9935aac738a31cde342637f2",
        )
        self.assertEqual(
            sha256(runtime.read_bytes()).hexdigest(),
            "c034ee71fb2c67f04f11e91b41ea8c9bf6a7b9df3a4bb72378e5b2696058c307",
        )
        self.assertEqual(
            sha256(validation.read_bytes()).hexdigest(),
            "7e936a9ac9c39e90bcb538f5b0562198f62ab20a4b4c3987c53d6c8d233cc063",
        )
        value = json.loads(queue.read_text(encoding="utf-8"))
        self.assertEqual(value["bundle_id"], "model-probe-501")
        self.assertEqual(value["priority"], 51)
        self.assertTrue(value["enabled"])
        self.assertEqual(value["logical_tasks"], 1)
        self.assertEqual(value["family"], "cell_foundation")
        self.assertEqual(value["exports"], {})
        self.assertEqual(
            value["wrapper_sha256"], sha256(wrapper.read_bytes()).hexdigest()
        )
        self.assertEqual(
            value["required_paths"],
            [
                "executions/geneformer-v2-316m-compat-runtime-21077474/ARTIFACTS.json",
                "executions/geneformer-v2-316m-preflight-validation-21077937/ARTIFACTS.json",
            ],
        )
        item = load_item(queue, ROOT.resolve(strict=True))
        self.assertTrue(item["_ready"])


if __name__ == "__main__":
    unittest.main()
