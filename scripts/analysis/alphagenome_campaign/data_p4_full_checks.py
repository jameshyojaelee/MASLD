#!/usr/bin/env python3
"""Compute-node smoke for worker serialization and P4 minimum-support rules."""
import csv
import gzip
import json
from multiprocessing import get_context
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

import data_p4 as p4
import data_p4_full as full
import data_p4_families as families


class FullReductionInvariants(unittest.TestCase):
    def test_original_minimum_each_side(self):
        d=pd.DataFrame({"x":np.arange(399,dtype=float)})
        result=families.measured_contrast(d,"x",np.arange(399)<200)
        self.assertEqual(result["status"],"indeterminate_minimum_group")
        self.assertIsNone(result["p_nominal"])

    def test_distinct_matched_control_support(self):
        strata=np.array(["shared"]*200+["shared"]*5+["unmatched"]*300)
        group=np.arange(len(strata))<200
        self.assertEqual(families.matching_support(strata,group),(200,5))

    def test_U2_actual_donor_axes_and_full_libraries(self):
        path=families.CTX / "counts/GSE244832/hepatocyte.peaks.tsv"
        selected=pd.read_csv(path,sep="\t",nrows=2).peak_coordinate.tolist()
        cov=families.u2_signal_covariates("GSE244832","hepatocyte",set(selected))
        self.assertEqual(set(cov.region_key),set(selected))
        self.assertTrue(np.isfinite(cov.signal_mean).all())
        self.assertTrue(np.isfinite(cov.signal_sd).all())
        self.assertTrue((cov.n_donors>=2).all())
        self.assertTrue((cov.source_count_columns>2).all())

    def test_two_worker_actual_archive(self):
        pilot=p4.ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/p4-21771964/PILOT_region_features.tsv.gz"
        with gzip.open(pilot,"rt") as handle:
            reader=csv.DictReader(handle,delimiter="\t")
            selected=[next(reader),next(reader)]
        tasks=[(str(p4.SAT / "raw/saturation"),r["region_key"],r["universe"],r["chrom"],int(r["start0"]),int(r["end"])) for r in selected]
        p4.OLD.group_columns=p4.group_columns
        with get_context("fork").Pool(2) as pool:
            results=list(pool.imap(full.process,tasks))
        for (row,top,state),expected in zip(results,selected):
            self.assertEqual(state,"ok")
            self.assertEqual(row["region_key"],expected["region_key"])
            self.assertTrue(set(row).issubset(full.feature_columns()))
            for column in ("liver_h3k27ac_rel_mean","liver_h3k27ac_adult_verified_rel_mean","liver_atac_rel_share_best50"):
                np.testing.assert_allclose(row[column],float(expected[column]),rtol=1e-12,atol=1e-12)


if __name__=="__main__":
    unittest.main()
