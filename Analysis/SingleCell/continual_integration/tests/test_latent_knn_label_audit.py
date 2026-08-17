import unittest
from pathlib import Path

import numpy as np

from masld_cl.config import load_config
from masld_cl.latent_knn_label_audit import (
    _fit_predict_knn,
    _majority_vote,
    load_latent_knn_label_audit_policy,
)


class TestLatentKNNLabelAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(__file__).resolve().parents[1]
        cls.config = load_config(cls.root / "config_v1.json")
        _, cls.policy = load_latent_knn_label_audit_policy(
            cls.config,
            cls.root / "reference" / "latent_knn_label_audit_policy_v25.json",
        )

    def test_majority_vote_has_lexicographic_tie_break(self):
        classes = np.asarray(["A", "B", "C"], dtype=object)
        codes = np.asarray([[0, 1], [1, 1], [2, 0]])
        self.assertEqual(_majority_vote(codes, classes).tolist(), ["A", "B", "A"])

    def test_fit_accepts_no_query_label_argument_and_is_seeded(self):
        rng = np.random.default_rng(17)
        reference = np.vstack([
            rng.normal(-2, 0.1, size=(40, 3)),
            rng.normal(2, 0.1, size=(40, 3)),
        ]).astype(np.float32)
        labels = np.asarray(["A"] * 40 + ["B"] * 40)
        query = np.asarray([[-2, -2, -2], [2, 2, 2]], dtype=np.float32)
        classifier = dict(self.policy["classifier"])
        classifier.update({"neighbors": 30, "n_jobs": 1})
        first, first_identity = _fit_predict_knn(reference, labels, query, classifier)
        second, second_identity = _fit_predict_knn(reference, labels, query, classifier)
        self.assertEqual(first.tolist(), ["A", "B"])
        self.assertTrue(np.array_equal(first, second))
        self.assertEqual(first_identity, second_identity)


if __name__ == "__main__":
    unittest.main()
