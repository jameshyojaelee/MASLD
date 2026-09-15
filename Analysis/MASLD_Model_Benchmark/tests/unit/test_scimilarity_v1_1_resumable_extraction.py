from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from scripts import extract_scimilarity_v1_1_frozen_screen_50000_resumable as resumable


class SCimilarityV11ResumableExtractionTest(unittest.TestCase):
    def test_atomic_chunks_resume_without_labels_or_heads(self) -> None:
        import numpy as np
        from scipy import sparse

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture_root = root / "fixture"
            fixture_root.mkdir()
            rows = [f"cell-{index}" for index in range(5)]
            receipt = {
                "row_order_sha256": resumable.sha256(
                    ("\n".join(rows) + "\n").encode("utf-8")
                ).hexdigest(),
                "normalization": "per_cell_total_10000_over_retained_model_genes_then_natural_log1p",
                "exposure_by_study": {f"study-{index}": "clean_declared" for index in range(7)},
                "study_gene_observability_order": [f"study-{index}" for index in range(7)],
                "study_gene_observability_is_global": True,
                "native_fixed_vocabulary_exception": "synthetic fixed-vocabulary fixture",
            }
            (fixture_root / "fixture_receipt.json").write_text(
                json.dumps(receipt), encoding="utf-8"
            )
            fixture = {
                "receipt": receipt,
                "matrix": sparse.csr_matrix(np.eye(5, 4, dtype=np.float32)),
                "observed_mask": np.ones(4, dtype=bool),
                "study_observability": np.ones((7, 4), dtype=bool),
                "outer_folds": np.arange(5, dtype=np.int8),
                "rows": rows,
            }

            def fake_embed(_model, matrix, *, device):
                self.assertEqual(device, "cpu")
                values = np.zeros((matrix.shape[0], 128), dtype=np.float32)
                values[:, 0] = 1.0
                return values

            checkpoint = root / "encoder.ckpt"
            checkpoint.write_bytes(b"synthetic")
            state = root / "state"
            first_output = root / "first-output"
            final_output = root / "final-output"
            with (
                mock.patch.object(resumable, "load_fixture", return_value=fixture),
                mock.patch.object(
                    resumable,
                    "_load_encoder",
                    return_value=(object(), {"parameter_count": 31_146_115}),
                ),
                mock.patch.object(resumable, "embed_batch", side_effect=fake_embed),
            ):
                first = resumable.extract_resumable(
                    fixture_root,
                    checkpoint,
                    state,
                    first_output,
                    expected_rows=5,
                    expected_dimension=4,
                    batch_size=1,
                    chunk_size=2,
                    device="cpu",
                    stop_after_new_chunks=1,
                )
                self.assertEqual(first["status"], "paused_after_atomic_chunk")
                self.assertFalse(first_output.exists())
                second = resumable.extract_resumable(
                    fixture_root,
                    checkpoint,
                    state,
                    final_output,
                    expected_rows=5,
                    expected_dimension=4,
                    batch_size=1,
                    chunk_size=2,
                    device="cpu",
                )
            self.assertEqual(
                second["status"], "pass_raw_outcome_blind_resumable_embeddings"
            )
            self.assertEqual(second["chunks_reused_this_attempt"], 1)
            self.assertFalse(second["released_reference_index_used"])
            self.assertFalse(second["downstream_head_fit"])
            self.assertFalse(second["common_head_fit"])
            self.assertFalse(second["metrics_run"])
            self.assertEqual(second["evaluation_label_columns_read"], [])
            self.assertFalse(second["sealed_outcomes_read"])
            with np.load(final_output / "raw_embeddings.npz", allow_pickle=False) as data:
                self.assertEqual(
                    set(data.files), {"embeddings", "outer_folds", "row_ids"}
                )
                self.assertEqual(data["embeddings"].shape, (5, 128))


if __name__ == "__main__":
    unittest.main()
