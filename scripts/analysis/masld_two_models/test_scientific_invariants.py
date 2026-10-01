#!/usr/bin/env python3
"""Focused checks of protected evaluation and signed-effect invariants."""
import unittest
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import allele_combination as allele
import chromatin_acquisition as chromatin
import variant_effect_query as variant_query


class ScientificInvariants(unittest.TestCase):
    def test_held_chromatin_cannot_change_predictions(self):
        rng = np.random.default_rng(67)
        x = rng.normal(size=(24, 45))
        y = rng.normal(size=(24, 30))
        tr = np.arange(18)
        held = np.arange(18, 24)
        selected = np.arange(12)
        cis = np.tile(np.arange(10), (len(selected), 1))
        original = chromatin.fit_profile(x, y, tr, held, cis, selected)
        perturbed = y.copy()
        perturbed[held] = rng.normal(size=(len(held), y.shape[1]))*1000
        changed = chromatin.fit_profile(x, perturbed, tr, held, cis, selected)
        for a, b in zip(original, changed):
            self.assertTrue(np.array_equal(a, b))

    def test_allele_reversal_preserves_zero_and_changes_direction(self):
        rng = np.random.default_rng(101)
        native = rng.normal(size=70)
        adapter = native*.8+rng.normal(scale=.2, size=70)
        effect = native*.6+adapter*.3+rng.normal(scale=.1, size=70)
        weight = allele.stack_fit(native, adapter, effect)
        reversed_weight = allele.stack_fit(-native, -adapter, -effect)
        self.assertTrue(np.allclose(weight, reversed_weight, atol=1e-12, rtol=1e-12))
        self.assertEqual(float(np.array([0.0, 0.0])@weight), 0.0)
        self.assertTrue(np.allclose(np.column_stack([-native, -adapter])@reversed_weight,
                                    -(np.column_stack([native, adapter])@weight)))

    def test_variant_query_requires_exact_source_peak_and_alleles(self):
        identity = "chr6:42911813:C:T"
        labels = pd.read_csv(variant_query.LABELS, sep="\t")
        row = labels.loc[labels.lead_variant_id.eq(identity)]
        self.assertEqual(len(row), 1)
        row = row.iloc[0]
        interval = f"{row.chr}:{int(row.peak_start_hg38)}-{int(row.peak_stop_hg38)}"
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/"predictions.tsv"
            pd.DataFrame([{"key": identity.removeprefix("chr"), "fold": row.heldout_fold,
                           "ref": row.ref, "alt": row.alt, "beta_alt": row.beta_alt,
                           "native2k": 0.1, "fixed_head": 0.1, "adapter": 0.1,
                           "native1m": 0.1, "combined": 0.1}]).to_csv(path, sep="\t", index=False)
            result = variant_query.query("GRCh38:"+identity, interval, "caQTL", path)
            self.assertEqual(result["source_peak_id"], row.peak_id)
            self.assertIsNone(result["target_gene"])
            wrong_interval = f"{row.chr}:{int(row.peak_start_hg38)+1}-{int(row.peak_stop_hg38)}"
            with self.assertRaises(ValueError):
                variant_query.query("GRCh38:"+identity, wrong_interval, "caQTL", path)
            with self.assertRaises(ValueError):
                variant_query.query("GRCh38:chr6:42911813:T:C", interval, "caQTL", path)


if __name__ == "__main__":
    unittest.main()
