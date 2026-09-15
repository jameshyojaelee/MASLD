from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]


class Scprint2EmbeddingContractTest(unittest.TestCase):
    def test_extractor_is_label_blind_and_neighbor_free(self) -> None:
        text = (ROOT / "scripts/scprint2_extract_no_neighbor_embeddings.py").read_text()
        fixture = (ROOT / "scripts/build_scprint2_no_neighbor_fixture.py").read_text()
        self.assertIn('get_knn_cells=False', text)
        self.assertIn('knn_cells=None', text)
        self.assertIn('knn_cells_info=None', text)
        self.assertIn('source_prefilter_depth_used', text)
        self.assertIn(
            '"released_selection_default_no_neighbor": ("random expr", 2000)',
            text,
        )
        self.assertIn('"random_checkpoint_max": ("random expr", 3200)', text)
        self.assertIn('output_obs["source_total_count"] = source_total_counts', fixture)
        self.assertIn('dtype=np.str_', text)
        self.assertIn('"source_depth_preserved_for_inference": True', fixture)
        self.assertIn('model.load_state_dict(state, strict=False)', text)
        self.assertIn('incompatibility.missing_keys', text)
        self.assertNotIn('broad_label', text)
        self.assertNotIn('GSE289173', text)

    def test_gpu_wrapper_uses_frozen_nslab_policy(self) -> None:
        text = (ROOT / "slurm/model_work_302.sbatch").read_text()
        self.assertIn('#SBATCH --partition=gpu', text)
        self.assertIn('#SBATCH --qos=nslab', text)
        self.assertIn('#SBATCH --gres=gpu:l40s:1', text)
        self.assertNotIn('#SBATCH --array', text)
        self.assertNotIn('--qos=innovation', text)
        self.assertIn('--fixture-root "${fixture}"', text)
        self.assertIn('neighbor_tensors_used', text)

    def test_common_evaluator_accepts_a_frozen_model_identity_and_width(self) -> None:
        text = (ROOT / "scripts/evaluate_transcriptformer_common_lane.py").read_text()
        self.assertIn('parser.add_argument("--model-id"', text)
        self.assertIn('parser.add_argument("--representation-id"', text)
        self.assertIn('parser.add_argument("--expected-width"', text)
        self.assertIn('"--row-id-source"', text)
        self.assertIn('"model_id": args.model_id', text)
        self.assertIn('"representation_id": args.representation_id', text)
        self.assertIn('(1000, args.expected_width)', text)
        self.assertIn('np.array_equal(activation_folds, source_folds)', text)
        self.assertNotIn('allow_pickle=True', text)

    def test_common_evaluation_wrapper_is_frozen_and_seal_blind(self) -> None:
        text = (ROOT / "slurm/evaluate_scprint2_common_lane.sbatch").read_text()
        self.assertIn('#SBATCH --partition=cpu', text)
        self.assertIn('#SBATCH --qos=nslab', text)
        self.assertIn(
            'embedding_artifact_sha256="3c57bcf5f9c821e395966809033b4319971432102714c0051f2c8d449e654183"',
            text,
        )
        self.assertIn('--expected-width 424', text)
        self.assertIn('--row-id-source "${row_id_source}"', text)
        self.assertIn('row_id_source_sha256="951b2615c437c86991bf039b885c453c51d9af3b3d39a111cf33a914d8489635"', text)
        self.assertIn('--bootstrap-resamples 10000', text)
        self.assertIn('common\n  random_checkpoint_max\n  released_selection_default_no_neighbor', text)
        self.assertIn('evaluation["sealed_outcomes_read"]', text)
        self.assertIn('evaluation["sealed_champion_eligible"]', text)
        self.assertNotIn('GSE289173', text)


if __name__ == "__main__":
    unittest.main()
