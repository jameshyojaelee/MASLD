from __future__ import annotations

from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "scripts/build_sequence_task_native_five_seed_rectangle.py"
SPEC = importlib.util.spec_from_file_location("sequence_rectangle", SOURCE)
assert SPEC is not None and SPEC.loader is not None
RECTANGLE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RECTANGLE)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class SequenceTaskNativeFiveSeedRectangleTest(unittest.TestCase):
    def test_contract_is_exact_full_rectangle(self) -> None:
        contract = json.loads(
            (ROOT / "config/campaigns/sequence_task_native_five_seed_rectangle_20260825.json").read_text()
        )
        self.assertEqual(contract["fixed_seeds"], [20260824, 20260825, 20260826, 20260827, 20260828])
        self.assertEqual(len(contract["lineages"]), 5)
        self.assertEqual(contract["diagonal_outer_folds"], list(range(5)))
        self.assertEqual(contract["fit_count_contract"]["total_fits"], 500)
        self.assertEqual(contract["fit_count_contract"]["total_prediction_views"], 625)
        self.assertTrue(contract["fit_count_contract"]["chrombpnet_full_and_nobias_share_one_fit"])
        self.assertFalse(contract["claims"]["cross_model_ranking_allowed"])

    def test_bundle_sizes_keep_jobs_long_but_checkpointable(self) -> None:
        self.assertEqual(RECTANGLE.balanced_sizes(105, 16, 8), [15] * 7)
        self.assertEqual(RECTANGLE.balanced_sizes(113, 16, 8), [15] + [14] * 7)
        self.assertEqual(sum(RECTANGLE.balanced_sizes(438, 16, 8)), 438)

    def test_old_bound_sources_are_byte_identical(self) -> None:
        expected = {
            "slurm/bpnet_train_predict_hepatocyte_donor0_genomic0.sbatch": "d22d3a840014921b7b8db1f801d37081f6e7092691228246c07cde7243c0386d",
            "slurm/train_chrombpnet_hepatocyte_donor0_genomic0_seed20260824.sbatch": "81966ffb9d63f3f4ca902fd1cf9ee641ac30dbff2701fdafc1d24ae9402f7edd",
            "scripts/sequence_control_production_driver.sh": "d0e256e09adb6416f7269493652fb5e1230e38a7e48c40a1ccba167fb9e103c3",
        }
        for relative, value in expected.items():
            self.assertEqual(digest(ROOT / relative), value, relative)

    def test_revisioned_sources_admit_only_the_five_fixed_seeds(self) -> None:
        for relative in (
            "slurm/bpnet_train_predict_five_seed_v2.sbatch",
            "slurm/chrombpnet_train_predict_five_seed_v2.sbatch",
            "scripts/sequence_control_production_driver_five_seed_v2.sh",
        ):
            text = (ROOT / relative).read_text()
            self.assertIn("20260824|20260825|20260826|20260827|20260828", text)
            self.assertNotIn("--qos=innovation", text)

    def test_gpu_wrapper_is_generic_nslab_and_not_arrayed(self) -> None:
        text = (ROOT / "slurm/run_sequence_task_native_five_seed_bundle.sbatch").read_text()
        self.assertIn("#SBATCH --job-name=model-training-700", text)
        self.assertIn("#SBATCH --partition=gpu", text)
        self.assertIn("#SBATCH --account=nslab", text)
        self.assertIn("#SBATCH --qos=nslab", text)
        self.assertIn("#SBATCH --gres=gpu:l40s:1", text)
        self.assertIn("#SBATCH --time=48:00:00", text)
        self.assertNotIn("#SBATCH --array", text)
        self.assertNotIn("--qos=innovation", text)


if __name__ == "__main__":
    unittest.main()
