from __future__ import annotations

from pathlib import Path
import unittest


class GSE296875TranscriptFormerGPUWrapperTest(unittest.TestCase):
    def test_generic_nslab_l40s_nonarray_dispatch_contract(self) -> None:
        root = Path(__file__).parents[2]
        text = (root / "slurm/run_gse296875_transcriptformer_query_embeddings.sbatch").read_text(
            encoding="utf-8"
        )
        self.assertIn("#SBATCH --job-name=model-work-335", text)
        self.assertIn("#SBATCH --partition=gpu", text)
        self.assertIn("#SBATCH --qos=nslab", text)
        self.assertIn("#SBATCH --gres=gpu:l40s:1", text)
        self.assertNotIn("#SBATCH --array", text)
        self.assertNotIn("--qos=innovation", text)
        self.assertIn("MASLD_GPU_DISPATCHER_AUTHORITY", text)
        self.assertNotIn("${CHECKPOINT}/ARTIFACTS.json", text)
        for digest in (
            "eff027e7393f92a32cadcbb2e13a917e24b62b71341571bf40f38ff6ec6d3021",
            "2413d512e8fdd550389bc574dc465357e943e23b4411cab6d23423af08095c2f",
            "676640be87ba4ba7ea7af8d9dc470a8ae8be02768989478027dec5fe61f465e2",
            "0c405e4ead45a4b8350d8e874f834273efdd3f7ad4669b52d3e3727fb4fe70af",
        ):
            self.assertIn(digest, text)


if __name__ == "__main__":
    unittest.main()
