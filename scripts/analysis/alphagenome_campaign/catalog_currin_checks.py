#!/usr/bin/env python3
"""Scientific identity checks for the Currin specimen extension, synthetic only."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from catalog_build import Catalog, sha
from catalog_currin import RECIPES, copy_specimen, matched_uncertainty, validate_predictions
from catalog_query import query


class CurrinChecks(unittest.TestCase):
    def setUp(self):
        self.labels = pd.DataFrame(dict(lead_variant_id=[f"chr1:{100+i}:A:G" for i in range(5)],
            peak_id=[f"peak{i}" for i in range(5)], beta_source=np.arange(5)/10,
            beta_alt=np.arange(5)/10, heldout_fold=np.arange(5), block_1mb=["chr1:0"]*5))
        self.source = self.labels.rename(columns={"lead_variant_id":"variant_id", "peak_id":"target_id",
            "beta_source":"beta_nominal"}).assign(varbeta=.04, effect_se=.2,
            uncertainty_state="sqrt_source_nominal_varbeta_matched_by_variant_and_peak")
        self.pred = self.labels.assign(key=self.labels.lead_variant_id.str.removeprefix("chr"))[
            ["key", "block_1mb", "heldout_fold", "beta_alt"]].copy()
        for name in RECIPES:
            self.pred[name] = self.pred.beta_alt + .1

    def test_uncertainty_requires_same_target(self):
        self.source.loc[0,"target_id"] = "wrong_peak"
        result = matched_uncertainty(self.labels, self.source)
        self.assertTrue(np.isnan(result.loc[0,"verified_se"]))
        self.assertEqual(result.loc[0,"se_status"],"missing_or_nonunique_exact_variant_peak")
        np.testing.assert_allclose(result.loc[1:,"verified_se"],.2)

    def test_unknown_and_conflicting_uncertainty_stay_missing(self):
        self.source.loc[0,"beta_nominal"] = 2
        self.source.loc[1,"varbeta"] = -1
        self.source.loc[2,"effect_se"] = .7
        self.source.loc[3,"uncertainty_state"] = "unknown"
        result = matched_uncertainty(self.labels,self.source)
        self.assertTrue(result.loc[:3,"verified_se"].isna().all())
        self.assertEqual(result.loc[4,"verified_se"],.2)

    def test_duplicate_source_is_not_replication(self):
        result = matched_uncertainty(self.labels,pd.concat([self.source,self.source.iloc[:1]]))
        self.assertEqual(len(result),5)
        self.assertTrue(np.isnan(result.loc[0,"verified_se"]))

    def test_fold_and_allele_effect_identity(self):
        validate_predictions(self.labels,self.pred)
        self.pred.loc[0,"heldout_fold"] = 1
        with self.assertRaises(AssertionError):
            validate_predictions(self.labels,self.pred)
        self.pred.loc[0,"heldout_fold"] = 0
        self.pred.loc[1,"beta_alt"] *= -1
        with self.assertRaises(AssertionError):
            validate_predictions(self.labels,self.pred)

    def test_incomplete_prediction_is_not_zero(self):
        self.pred.loc[0,RECIPES[0]] = np.nan
        with self.assertRaisesRegex(ValueError,"Incomplete"):
            validate_predictions(self.labels,self.pred)

    def test_copy_keeps_original_and_separate_prediction_uncertainty(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp)/"base.sqlite"
            original=Catalog(base)
            original.db.execute("INSERT INTO source VALUES ('fixture','synthetic','fixture','synthetic','not_scientific_evidence')")
            variant=original.entity("variant","GRCh38:chr1:100:A:G")
            original.evidence("fixture","measured","caQTL","measured","caQTL",.3,"beta","ALT_dosage",[(variant,"variant")],se=.2)
            original.db.commit()
            original.db.close()
            before=sha(base)
            copied=Path(tmp)/"copy.sqlite"
            catalog=copy_specimen(base,copied)
            catalog.evidence("fixture","predicted","caQTL","predicted","caQTL",.4,"beta","ALT_dosage",[(variant,"variant")])
            catalog.db.commit()
            catalog.db.close()
            self.assertEqual(sha(base),before)
            self.assertEqual(query(base,"variant","GRCh38:chr1:100:A:G")["total_evidence"],1)
            records=query(copied,"variant","GRCh38:chr1:100:A:G")["evidence"]
            self.assertEqual(len(records),2)
            self.assertIsNone(next(row for row in records if row["kind"]=="predicted")["se"])
            self.assertEqual(next(row for row in records if row["kind"]=="measured")["se"],.2)
            with self.assertRaises(FileExistsError):
                copy_specimen(base,copied)


if __name__ == "__main__":
    unittest.main()
