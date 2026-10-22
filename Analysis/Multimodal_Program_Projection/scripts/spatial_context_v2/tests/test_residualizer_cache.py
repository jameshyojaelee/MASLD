#!/usr/bin/env python3
"""Check complete spatial-null equivalence when nuisance inverses are reused."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy import sparse

SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from visium_rerun_lib import build_paths, cached_residualizer, load_legacy_engine


class ResidualizerCacheSemantics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine = load_legacy_engine(build_paths())
        cls.original = staticmethod(cls.engine.residualize.__wrapped__)

    def test_repeated_mutated_and_evicted_designs(self):
        rng = np.random.default_rng(42)
        design = np.column_stack([np.ones(24), rng.normal(size=(24, 2))])
        values = rng.normal(size=(24, 7))
        cached = cached_residualizer(self.original)
        with patch.object(np.linalg, "pinv", wraps=np.linalg.pinv) as pinv:
            for _ in range(2):
                np.testing.assert_array_equal(cached(values, design), self.original(values, design))
            # The baseline evaluations also call pinv: three calls means one cached fit.
            self.assertEqual(pinv.call_count, 3)
            design[:, 1] *= 2
            np.testing.assert_array_equal(cached(values, design), self.original(values, design))
            self.assertEqual(pinv.call_count, 5)
        # The cache stays bounded while revisiting evicted designs remains exact.
        for offset in range(6):
            changed = design.copy()
            changed[:, 1] += offset
            np.testing.assert_array_equal(cached(values, changed), self.original(values, changed))
        np.testing.assert_array_equal(cached(values, design), self.original(values, design))

    def test_missing_scores_and_rank_deficient_designs(self):
        rng = np.random.default_rng(73)
        values = rng.normal(size=(24, 3))
        values[3, 1] = np.nan
        for design in (
            np.zeros((24, 3)),
            np.ones((24, 3)),
            np.asfortranarray(np.column_stack([np.ones(24), np.arange(24)])),
        ):
            cached = cached_residualizer(self.original)
            for score in (values, values[:, 0], values.astype(np.float32)):
                np.testing.assert_array_equal(cached(score, design), self.original(score, design))

    def test_complete_nulls_and_random_stream_are_identical(self):
        engine = self.engine
        rng = np.random.default_rng(20261021)
        norm = sparse.csr_matrix(rng.poisson(2, (36, 30)).astype(np.float32))
        obs = pd.DataFrame({
            "sample_id": np.repeat(["H35a", "H35b", "H36"], 12),
            "individual": np.repeat(["H35", "H35", "H36"], 12),
        })
        graphs = {}
        for number, sample in enumerate(["H35a", "H35b", "H36"]):
            graphs[sample] = [(np.arange(number * 12, (number + 1) * 12),
                               np.arange(12), np.roll(np.arange(12), -1))]
        design = np.column_stack([np.ones(36), rng.normal(size=(36, 3))])
        zonation = np.column_stack([design, rng.normal(size=36)])
        weights = np.array([0.5, 0.3, 0.2])
        measured = pd.DataFrame({"gene_symbol": ["G00", "G01", "G02"]})
        pool = pd.DataFrame({
            "gene": [f"G{i:02d}" for i in range(30)],
            "biotype": ["protein_coding"] * 30,
            "mt": [False] * 30, "ribo": [False] * 30,
            "mean_bin": [5] * 30, "detect_bin": [5] * 30,
            "lineage_corr_bin": [5] * 30,
        })
        outputs, states = [], []
        installed = engine.residualize
        try:
            for residualizer in (self.original, cached_residualizer(self.original)):
                engine.residualize = residualizer
                stream = np.random.default_rng(42)
                matched, _, _ = engine.build_matched_sets(
                    measured, pool, {"G00", "G01", "G02"}, 73, stream
                )
                primary = engine.null_moran(
                    norm, matched, weights, design, zonation, graphs, obs,
                    "GSE192741", chunk_size=16,
                )
                equal = engine.null_residual_moran(
                    norm, matched, np.repeat(1 / 3, 3), design, graphs, obs,
                    "GSE192741", chunk_size=16,
                )
                leave = engine.null_residual_moran(
                    norm, matched[:, 1:], weights[1:] / weights[1:].sum(),
                    design, graphs, obs, "GSE192741", chunk_size=16,
                )
                outputs.append((*primary, equal, leave))
                states.append(stream.bit_generator.state)
        finally:
            engine.residualize = installed
        for original, cached in zip(*outputs):
            np.testing.assert_array_equal(original, cached)
        self.assertEqual(states[0], states[1])


if __name__ == "__main__":
    unittest.main(verbosity=2)
