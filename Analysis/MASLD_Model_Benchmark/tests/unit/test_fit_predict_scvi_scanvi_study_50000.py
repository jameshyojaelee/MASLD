#!/usr/bin/env python3
"""Contract checks for the frozen 50k scVI/scANVI screen."""

from __future__ import annotations

from pathlib import Path
import unittest

import numpy as np
from scipy import sparse

from scripts import fit_predict_scvi_scanvi_study_50000 as screen
from masld_bench.adapters.scvi_scanvi_inductive import UNLABELED_CATEGORY


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts/fit_predict_scvi_scanvi_study_50000.py"


class FitPredictSCVISCANVIStudy50000Tests(unittest.TestCase):
    def test_frozen_hyperparameters_and_seeds(self) -> None:
        self.assertEqual(screen.SEEDS, (20260824, 20260825, 20260826))
        self.assertEqual(screen.N_TOP_HVG, 4_000)
        self.assertEqual(screen.N_LATENT, 30)
        self.assertEqual(screen.N_HIDDEN, 128)
        self.assertEqual(screen.N_LAYERS, 2)
        self.assertEqual(screen.SCVI_EPOCHS, 200)
        self.assertEqual(screen.SCANVI_EPOCHS, 100)
        self.assertEqual(screen.BATCH_SIZE, 512)
        self.assertEqual(screen.MODEL_IDS, ("scvi_baseline", "scanvi_baseline"))

    def test_model_facing_query_contains_no_truth(self) -> None:
        genes = [f"ENSG{index:011d}" for index in range(20)]
        training_labels = np.asarray(list(screen.ROSTER), dtype=str)
        training, query = screen._training_and_query_objects(
            raw_train=sparse.csr_matrix(np.ones((5, 20), dtype=np.float32)),
            raw_query=sparse.csr_matrix(np.ones((3, 20), dtype=np.float32)),
            gene_ids=genes,
            train_row_ids=np.asarray([f"train-{index}" for index in range(5)]),
            train_donors=np.asarray([f"donor-{index}" for index in range(5)]),
            train_studies=np.asarray(["training-study"] * 5),
            train_labels=training_labels,
            query_row_ids=np.asarray([f"query-{index}" for index in range(3)]),
            query_donors=np.asarray(["held-donor"] * 3),
            query_studies=np.asarray(["held-study"] * 3),
        )
        self.assertEqual(training.shape, (5, 20))
        self.assertEqual(query.shape, (3, 20))
        self.assertEqual(set(map(str, query.obs["broad_label"])), {UNLABELED_CATEGORY})
        self.assertNotIn("held-study", set(map(str, training.obs["study"])))

    def test_training_and_query_studies_must_be_disjoint(self) -> None:
        with self.assertRaises(screen.SCVIFrozenScreenError):
            screen._training_and_query_objects(
                raw_train=sparse.csr_matrix(np.ones((5, 10))),
                raw_query=sparse.csr_matrix(np.ones((2, 10))),
                gene_ids=[f"ENSG{index:011d}" for index in range(10)],
                train_row_ids=[f"train-{index}" for index in range(5)],
                train_donors=[f"donor-{index}" for index in range(5)],
                train_studies=["same-study"] * 5,
                train_labels=list(screen.ROSTER),
                query_row_ids=["query-0", "query-1"],
                query_donors=["held-donor", "held-donor"],
                query_studies=["same-study", "same-study"],
            )

    def test_source_has_no_query_adaptation_or_common_head(self) -> None:
        text = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("load_query_data", text)
        self.assertNotIn("prepare_query_anndata", text)
        self.assertNotIn("LogisticRegression", text)
        self.assertNotIn("MLPClassifier", text)
        self.assertIn("labels[train]", text)
        self.assertNotIn("labels[query]", text)
        self.assertIn("early_stopping=False", text)
        self.assertIn("train_size=1.0", text)
        self.assertIn("set(studies[train]) & set(studies[query])", text)
        self.assertIn("set(inner_by_outer[fold]) != set(donors[train])", text)


if __name__ == "__main__":
    unittest.main()
