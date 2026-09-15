from __future__ import annotations

from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class SequenceTaskNativeFiveSeedExecutionCorrectionTest(unittest.TestCase):
    def test_controller_diff_is_exact_import_and_summary_extraction(self) -> None:
        old = (ROOT / "scripts/run_sequence_task_native_five_seed_bundle.py").read_text().splitlines()
        new = (ROOT / "scripts/run_sequence_task_native_five_seed_bundle_v2.py").read_text().splitlines()
        expected = list(old)
        expected.insert(expected.index("import csv"), "from collections import Counter")
        main_index = expected.index("def main() -> None:")
        expected[main_index:main_index] = [
            "def terminal_summary(receipts: Mapping[str, Mapping[str, Any]]) -> tuple[Counter[str], str]:",
            '    """Summarize one synthetic or real completed bundle without reading predictions."""',
            '    counts = Counter(row["status"] for row in receipts.values())',
            '    terminal_status = "passed" if set(counts) <= PASS_STATES else "completed_with_fit_failures"',
            "    return counts, terminal_status",
            "",
            "",
        ]
        old_summary = [
            '    counts = Counter(row["status"] for row in receipts.values())',
            '    terminal_status = "passed" if set(counts) <= PASS_STATES else "completed_with_fit_failures"',
        ]
        summary_index = next(
            index for index in range(len(expected) - 2, -1, -1)
            if expected[index : index + 2] == old_summary
        )
        expected[summary_index : summary_index + 2] = [
            "    counts, terminal_status = terminal_summary(receipts)"
        ]
        self.assertEqual(new, expected)

    def test_v2_controller_imports_counter(self) -> None:
        source = ROOT / "scripts/run_sequence_task_native_five_seed_bundle_v2.py"
        spec = importlib.util.spec_from_file_location("sequence_bundle_v2", source)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertEqual(module.Counter(["passed", "passed"]), {"passed": 2})
        counts, status = module.terminal_summary({
            "fit-a": {"status": "passed"},
            "fit-b": {"status": "recovered_passed"},
        })
        self.assertEqual(counts, {"passed": 1, "recovered_passed": 1})
        self.assertEqual(status, "passed")
        failed_counts, failed_status = module.terminal_summary({
            "fit-a": {"status": "passed"},
            "fit-b": {"status": "failed"},
        })
        self.assertEqual(failed_counts, {"passed": 1, "failed": 1})
        self.assertEqual(failed_status, "completed_with_fit_failures")

    def test_correction_binds_exact_sources(self) -> None:
        correction = json.loads(
            (ROOT / "config/campaigns/sequence_task_native_five_seed_execution_correction_20260825.json").read_text()
        )
        for label in ("v1_controller", "v2_controller", "v1_bundle_wrapper", "v2_bundle_wrapper"):
            binding = correction[label]
            self.assertEqual(digest(ROOT / binding["path"]), binding["sha256"])
        self.assertFalse(correction["diff_contract"]["fit_inputs_changed"])
        self.assertFalse(correction["diff_contract"]["model_outputs_changed"])
        self.assertFalse(correction["queue_disposition"]["manual_gpu_submission_allowed"])

    def test_v1_queue_is_disabled_and_unclaimed(self) -> None:
        queue = ROOT / "config/campaigns/gpu_bundle_queue"
        state = ROOT / "executions/gpu-bundle-dispatch-state/claims"
        paths = sorted(queue.glob("0[7-9][0-9]-model-training-7*.json"))
        self.assertEqual(len(paths), 30)
        for path in paths:
            item = json.loads(path.read_text())
            self.assertFalse(item["enabled"])
            self.assertFalse((state / f"{item['bundle_id']}.json").exists())


if __name__ == "__main__":
    unittest.main()
