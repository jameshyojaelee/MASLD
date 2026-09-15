from __future__ import annotations

import json
from pathlib import Path
import unittest

import numpy as np

from scripts.gse281364_dna_lm_native_contract import (
    affected_interval,
    apply_projection,
    fit_projection,
    head_features,
    load_config,
    nt_phase_spans,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/gse281364_dna_lm_native_contract.json"


class GSE281364DnaLmNativeContractTests(unittest.TestCase):
    def test_nt_phase_geometry_has_fixed_ten_base_edge_trim(self) -> None:
        spans = nt_phase_spans()
        self.assertEqual(spans, [(phase, phase + 4_086) for phase in range(6)])
        self.assertEqual(
            {start + (4_096 - end) for start, end in spans},
            {10},
        )
        self.assertTrue(all(start <= 2_048 < end for start, end in spans))

    def test_dnabert_affected_interval_closes_across_both_tokenizations(self) -> None:
        reference = [(0, 0), (0, 2), (2, 5), (5, 8), (0, 0)]
        alternative = [(0, 0), (0, 3), (3, 6), (6, 8), (0, 0)]
        self.assertEqual(affected_interval(reference, alternative, 4), (0, 8))

    def test_projection_fit_excludes_held_out_fold(self) -> None:
        random = np.random.default_rng(20260824)
        folds = np.repeat(np.arange(5), 4)
        embeddings = random.normal(size=(20, 4, 8))
        first = fit_projection(embeddings, folds, held_out_fold=4, width=4)
        altered = embeddings.copy()
        altered[folds == 4] += 10_000
        second = fit_projection(altered, folds, held_out_fold=4, width=4)
        for field in ("mean", "standard_deviation", "components", "whitening_scale"):
            np.testing.assert_allclose(first[field], second[field], atol=0, rtol=0)
        transformed = apply_projection(embeddings, first)
        self.assertEqual(transformed.shape, (20, 4, 4))

    def test_head_blocks_are_strand_averaged_in_frozen_order(self) -> None:
        values = np.asarray([[[1.0], [3.0], [5.0], [7.0]]], dtype=np.float32)
        observed = head_features(values)
        np.testing.assert_array_equal(observed, [[3.0, 5.0, 2.0, 2.0]])

    def test_config_binds_exact_models_and_head_width(self) -> None:
        config = load_config(CONFIG)
        self.assertEqual(
            list(config["models"]),
            ["dnabert2", "nucleotide_transformer", "hyenadna"],
        )
        self.assertEqual(config["normalization_projection"]["head_input_width"], 1024)
        self.assertFalse(any(json.loads(CONFIG.read_text())["firewall"].values()))


if __name__ == "__main__":
    unittest.main()
