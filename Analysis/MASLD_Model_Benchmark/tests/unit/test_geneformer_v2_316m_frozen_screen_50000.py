from __future__ import annotations

import csv
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

from scripts import build_geneformer_v2_316m_frozen_screen_50000 as builder
from scripts import extract_geneformer_v2_316m_frozen_screen_50000 as extractor


ROOT = Path(__file__).resolve().parents[2]
REAL_FIXTURE = ROOT / "executions/geneformer-v2-316m-frozen-screen-fixture-21078339"
GPU_WRAPPER = ROOT / "slurm/run_geneformer_v2_316m_frozen_screen_50000.sbatch"


class GeneformerV2316MFrozenScreen50000Test(unittest.TestCase):
    def _inputs(self, root: Path, *, forbidden: bool = False) -> tuple[Path, ...]:
        import anndata
        import numpy as np
        import pandas as pd
        from scipy import sparse

        row_ids = [f"cell-{index}" for index in range(5)]
        obs_values: dict[str, object] = {
            "row_id": row_ids,
            "donor_id": [f"donor-{index}" for index in range(5)],
            "dataset_id": [f"study-{index}" for index in range(5)],
            "outer_fold": [0, 1, 2, 3, 4],
            "assay": ["unknown"] * 5,
        }
        if forbidden:
            obs_values["broad_label"] = ["hepatocyte"] * 5
        obs = pd.DataFrame(obs_values, index=pd.Index(row_ids, name="row_id_index"))
        var = pd.DataFrame(index=pd.Index(["ENSG000001", "ENSG000002", "ENSG000003"]))
        matrix = sparse.csr_matrix(
            np.asarray(
                [[10, 5, 0], [0, 2, 8], [4, 0, 1], [1, 1, 1], [0, 3, 2]],
                dtype=np.float32,
            )
        )
        source = root / "unlabeled.h5ad"
        anndata.AnnData(X=matrix, obs=obs, var=var).write_h5ad(source)
        split = root / "split.tsv"
        with split.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["row_id", "donor_id", "dataset", "outer_fold"])
            for index, row_id in enumerate(row_ids):
                writer.writerow([row_id, f"donor-{index}", f"study-{index}", index])
        medians = root / "median.json"
        medians.write_text(
            json.dumps({"ENSG000001": 1.0, "ENSG000002": 10.0, "ENSG000003": 2.0}),
            encoding="utf-8",
        )
        token_values = {
            f"ENSG{index:06d}": index for index in range(1, 20_276)
        }
        tokens = root / "tokens.json"
        tokens.write_text(json.dumps(token_values), encoding="utf-8")
        config = root / "config.json"
        config.write_text(
            json.dumps(
                {
                    "architectures": ["BertForMaskedLM"],
                    "hidden_size": 1152,
                    "max_position_embeddings": 4096,
                    "vocab_size": 20275,
                }
            ),
            encoding="utf-8",
        )
        return source, split, medians, tokens, config

    def test_builds_outcome_blind_tokens_in_source_row_order(self) -> None:
        import numpy as np

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, split, medians, tokens, config = self._inputs(root)
            output = root / "fixture"
            receipt = builder.build(
                source,
                split,
                medians,
                tokens,
                config,
                output,
                expected_source_sha256=builder.sha256_file(source),
                expected_split_sha256=builder.sha256_file(split),
                expected_rows=5,
                expected_donors=5,
                expected_studies=5,
            )
            self.assertEqual(receipt["status"], "pass_outcome_blind_fixture")
            self.assertEqual(receipt["evaluation_label_columns_read"], [])
            self.assertEqual(receipt["histology_columns_read"], [])
            self.assertFalse(receipt["sealed_outcomes_read"])
            self.assertEqual(
                (output / "embedding_row_order.txt").read_text().splitlines(),
                [f"cell-{index}" for index in range(5)],
            )
            token_ids = np.load(output / "token_ids.npy", allow_pickle=False)
            offsets = np.load(output / "token_offsets.npy", allow_pickle=False)
            self.assertEqual(offsets.shape, (6,))
            self.assertEqual(token_ids[offsets[0] : offsets[1]].tolist(), [1, 2])
            loaded = extractor.load_fixture(output, expected_rows=5)
            self.assertEqual(loaded["rows"], [f"cell-{index}" for index in range(5)])

    def test_rejects_a_fixture_that_exposes_labels(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, split, medians, tokens, config = self._inputs(root, forbidden=True)
            with self.assertRaisesRegex(builder.GeneformerFixtureError, "firewall"):
                builder.build(
                    source,
                    split,
                    medians,
                    tokens,
                    config,
                    root / "fixture",
                    expected_source_sha256=builder.sha256_file(source),
                    expected_split_sha256=builder.sha256_file(split),
                    expected_rows=5,
                    expected_donors=5,
                    expected_studies=5,
                )

    def test_rejects_a_study_split_join_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, split, medians, tokens, config = self._inputs(root)
            text = split.read_text(encoding="utf-8").replace("donor-3", "wrong-donor")
            split.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(builder.GeneformerFixtureError, "join"):
                builder.build(
                    source,
                    split,
                    medians,
                    tokens,
                    config,
                    root / "fixture",
                    expected_source_sha256=builder.sha256_file(source),
                    expected_split_sha256=builder.sha256_file(split),
                    expected_rows=5,
                    expected_donors=5,
                    expected_studies=5,
                )

    def test_real_fixture_is_row_stable_and_outcome_blind(self) -> None:
        self.assertEqual(
            builder.sha256_file(REAL_FIXTURE / "ARTIFACTS.json"),
            "bb94045a46c7231d50a17ffd921a146ad6858b39bcfcda516da61230fcfb9bf6",
        )
        value = extractor.load_fixture(REAL_FIXTURE / "fixture")
        self.assertEqual(len(value["rows"]), 50_000)
        self.assertEqual(value["receipt"]["evaluation_label_columns_read"], [])
        self.assertEqual(value["receipt"]["histology_columns_read"], [])
        self.assertFalse(value["receipt"]["sealed_outcomes_read"])

    def test_gpu_wrapper_is_generic_nslab_and_fail_closed(self) -> None:
        text = GPU_WRAPPER.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --job-name=model-training-503", text)
        self.assertIn("#SBATCH --partition=gpu", text)
        self.assertIn("#SBATCH --account=nslab", text)
        self.assertIn("#SBATCH --qos=nslab", text)
        self.assertIn("#SBATCH --gres=gpu:l40s:1", text)
        self.assertNotIn("innovation", text.lower())
        self.assertNotIn("#SBATCH --array", text)
        self.assertFalse(
            any(line.strip().startswith("sbatch ") for line in text.splitlines())
        )
        self.assertIn("ONE_BATCH_ARTIFACTS_SHA256", text)
        self.assertIn("SCREEN_PREFLIGHT_ARTIFACTS_SHA256", text)
        self.assertIn("sequence_gpu_dispatcher_v2", text)
        self.assertIn(
            "FIXTURE_SHA256=bb94045a46c7231d50a17ffd921a146ad6858b39bcfcda516da61230fcfb9bf6",
            text,
        )
        queue = ROOT / "config/campaigns/gpu_bundle_queue"
        references = [
            path
            for path in queue.glob("*.json")
            if "run_geneformer_v2_316m_frozen_screen_50000.sbatch"
            in path.read_text(encoding="utf-8")
        ]
        self.assertEqual(
            references,
            [queue / "052-model-training-503.json"],
        )
        value = json.loads(references[0].read_text(encoding="utf-8"))
        self.assertEqual(value["bundle_id"], "model-training-503")
        self.assertEqual(value["priority"], 52)
        self.assertTrue(value["enabled"])
        self.assertEqual(value["logical_tasks"], 1)
        self.assertEqual(value["family"], "cell_foundation")
        self.assertEqual(
            value["wrapper_sha256"],
            sha256(GPU_WRAPPER.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            value["exports"]["ONE_BATCH_ARTIFACTS_SHA256"],
            "050263e5cf675fb3b421a751cf06b57f294a322d7128b1ff63a2f6f633764f68",
        )
        self.assertEqual(
            value["exports"]["SCREEN_PREFLIGHT_ARTIFACTS_SHA256"],
            "c6d9b1c6f8ad22b7f7fc15d4ba8d18c6c60fbd2681e8668cfd70871949dd4f94",
        )
        self.assertIn(
            "executions/model-check-075-21081348/ARTIFACTS.json",
            value["required_paths"],
        )


if __name__ == "__main__":
    unittest.main()
