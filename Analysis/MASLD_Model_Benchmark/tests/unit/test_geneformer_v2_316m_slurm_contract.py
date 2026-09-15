from __future__ import annotations

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
CPU = ROOT / "slurm/build_geneformer_v2_316m_compat_runtime_cpu.sbatch"
GPU = ROOT / "slurm/preflight_geneformer_v2_316m_one_batch_l40s.sbatch"


class GeneformerV2316MSlurmContractTest(unittest.TestCase):
    def test_cpu_runtime_gate_is_cpu_only_and_install_free(self) -> None:
        text = CPU.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=cpu", text)
        self.assertIn("#SBATCH --qos=nslab", text)
        self.assertNotIn("#SBATCH --gres=gpu", text)
        self.assertNotIn("--array=", text)
        self.assertNotIn("pip install", text)
        self.assertNotIn("conda install", text)
        self.assertNotIn("micromamba create", text)
        self.assertIn("--mode cpu", text)

    def test_gpu_preflight_is_hash_bound_generic_and_not_self_submitting(self) -> None:
        text = GPU.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --job-name=model-probe-043", text)
        self.assertIn("#SBATCH --partition=gpu", text)
        self.assertIn("#SBATCH --account=nslab", text)
        self.assertIn("#SBATCH --qos=nslab", text)
        self.assertIn("#SBATCH --gres=gpu:l40s:1", text)
        self.assertNotIn("innovation", text.lower())
        self.assertNotIn("--array=", text)
        self.assertNotIn("sbatch ", text)
        self.assertNotIn("pip install", text)
        self.assertNotIn("conda install", text)
        self.assertIn(
            "RUNTIME_ARTIFACTS_SHA256="
            "c034ee71fb2c67f04f11e91b41ea8c9bf6a7b9df3a4bb72378e5b2696058c307",
            text,
        )
        self.assertIn(
            "RUNTIME_ARCHIVE_SHA256="
            "a4fdc1f1d0b67c4753354c5a04ccc05c84a3bf7442d8655bc8ed6493a2487775",
            text,
        )
        self.assertIn("--mode gpu", text)
        self.assertIn("--output-embeddings", text)


if __name__ == "__main__":
    unittest.main()
