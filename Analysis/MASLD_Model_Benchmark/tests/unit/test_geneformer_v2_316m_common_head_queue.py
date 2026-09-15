#!/usr/bin/env python3
"""Admission checks for the Geneformer 316M common-head dependency chain."""

from __future__ import annotations

import json
from pathlib import Path
import unittest

from masld_bench.artifacts import verify_frozen_tree
from masld_bench.hashing import sha256_file
from scripts.fit_predict_common_cell_heads_study_50000 import derive_embedding_identity
from scripts.gpu_bundle_dispatcher import load_item


ROOT = Path(__file__).parents[2]
QUEUE = ROOT / "config/campaigns/gpu_bundle_queue/058-model-training-506.json"
CHAIN = ROOT / "config/campaigns/geneformer_v2_316m_common_head_chain_v1.json"
WRAPPER = ROOT / "slurm/run_common_cell_head_bundle_50000.sbatch"
AUDITOR = ROOT / "slurm/audit_aggregate_common_cell_heads_study_50000.sbatch"
SCORER = ROOT / "slurm/score_common_cell_heads_study_50000.sbatch"
EMBEDDINGS = ROOT / "executions/model-data-053-21088295/embeddings"


class GeneformerCommonHeadQueueTests(unittest.TestCase):
    def test_prepared_embedding_identity_is_exact(self) -> None:
        manifest = verify_frozen_tree(EMBEDDINGS)
        self.assertEqual(
            sha256_file(EMBEDDINGS / "ARTIFACTS.json"),
            "a675daa170020f4e12eac2eb6dcf046cac4304511789805f802e56d08900f86d",
        )
        self.assertFalse(manifest["metadata"]["evaluation_labels_read"])
        self.assertFalse(manifest["metadata"]["sealed_outcomes_read"])
        identity, path = derive_embedding_identity(
            EMBEDDINGS, Path("embeddings/common_embeddings.npz")
        )
        self.assertEqual(
            identity,
            {
                "model_id": "geneformer_v2_316m",
                "representation_id": "common",
                "checkpoint_sha256": "965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3",
                "exposure_status": "target_label_unexposed",
                "embedding_width": 1152,
                "embedding_relative": "embeddings/common_embeddings.npz",
                "sealed_champion_eligible": True,
            },
        )
        self.assertEqual(path.stat().st_size, 238_850_774)

    def test_queue_is_one_dispatcher_only_nslab_bundle(self) -> None:
        value = json.loads(QUEUE.read_text(encoding="utf-8"))
        self.assertEqual(value["bundle_id"], "model-training-506")
        self.assertEqual(value["priority"], 58)
        self.assertTrue(value["enabled"])
        self.assertEqual(value["logical_tasks"], 30)
        self.assertEqual(value["family"], "cell_foundation_common_head")
        self.assertEqual(value["wrapper_sha256"], sha256_file(WRAPPER))
        self.assertEqual(
            value["exports"],
            {
                "EMBEDDINGS_ARTIFACTS_SHA256": "a675daa170020f4e12eac2eb6dcf046cac4304511789805f802e56d08900f86d",
                "EMBEDDINGS_RELATIVE": "embeddings/common_embeddings.npz",
                "EMBEDDINGS_ROOT": EMBEDDINGS.as_posix(),
                "HEAD_IDS": "linear two_layer_mlp",
                "MASLD_GPU_DISPATCHER_AUTHORITY": "sequence_gpu_dispatcher_v2",
                "OUTER_FOLDS": "0 1 2 3 4",
                "SCREEN_SEEDS": "1103 1201 1301",
            },
        )
        self.assertIn(
            "executions/geneformer-v2-316m-common-head-queue-validation-v1/ARTIFACTS.json",
            value["required_paths"],
        )
        load_item(QUEUE, ROOT.resolve(strict=True))
        text = WRAPPER.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=gpu", text)
        self.assertIn("#SBATCH --qos=nslab", text)
        self.assertIn("#SBATCH --gres=gpu:l40s:1", text)
        self.assertNotIn("#SBATCH --array", text)
        self.assertNotIn("innovation", text.lower())
        self.assertNotIn("masld", next(line for line in text.splitlines() if "--job-name" in line).lower())

    def test_downstream_chain_is_hash_bound_and_no_seal(self) -> None:
        value = json.loads(CHAIN.read_text(encoding="utf-8"))
        self.assertEqual(value["model_id"], "geneformer_v2_316m")
        self.assertEqual(value["logical_gpu_shards"], 30)
        self.assertFalse(value["evaluation_labels_read_before_scoring"])
        self.assertFalse(value["histology_read"])
        self.assertFalse(value["sealed_outcomes_read"])
        self.assertFalse(value["project_sealed_external_evaluation_performed"])
        stages = {stage["stage_id"]: stage for stage in value["stages"]}
        self.assertEqual(
            stages["bundled_common_head_fit_predict"]["wrapper_sha256"],
            sha256_file(WRAPPER),
        )
        self.assertEqual(
            stages["label_free_shard_audit"]["wrapper_sha256"],
            sha256_file(AUDITOR),
        )
        self.assertEqual(
            stages["development_label_join_and_scoring"]["wrapper_sha256"],
            sha256_file(SCORER),
        )
        self.assertIn(
            "EXPECTED_MODEL_ID=geneformer_v2_316m",
            stages["label_free_shard_audit"]["required_exports"],
        )
        self.assertIn(
            "EXPECTED_EXPOSURE_STATUS=target_label_unexposed",
            stages["development_label_join_and_scoring"]["required_exports"],
        )


if __name__ == "__main__":
    unittest.main()
