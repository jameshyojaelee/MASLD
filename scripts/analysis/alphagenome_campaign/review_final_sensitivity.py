#!/usr/bin/env python3
"""Controlled numerical checks for whole-chromosome uncertainty sensitivity."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from model_final import describe_both


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run numerical checks on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    rows = []
    for chrom in range(1, 11):
        for site in range(3):
            y = chrom*.13 + site*.17
            rows.append(dict(key=f"{chrom}:{100+site}:A:G", heldout_fold=(chrom-1)//2,
                block_1mb=f"chr{chrom}:{site//2}", beta_alt=y, perfect=y,
                same_as_anchor=y+1, anchor=y+1))
    frame = pd.DataFrame(rows)
    untouched = frame.copy(deep=True)
    describe_both(frame,["perfect","same_as_anchor","anchor"],["anchor"],"synthetic",args.out,100,30)
    pd.testing.assert_frame_equal(frame,untouched)
    original = pd.read_csv(args.out/"synthetic_performance.tsv",sep="\t").set_index("arm")
    sensitivity = pd.read_csv(args.out/"synthetic_chromosome_sensitivity_performance.tsv",sep="\t").set_index("arm")
    assert original.blocks.eq(20).all() and sensitivity.blocks.eq(10).all()
    assert original.n.eq(30).all() and sensitivity.n.eq(30).all()
    assert sensitivity.uncertainty_unit.eq("paired_whole_chromosomes_within_fold").all()
    for column in ("RMSE_beta_units","MAE_beta_units","pooled_signed_spearman","macro_spearman"):
        np.testing.assert_allclose(original[column],sensitivity[column],atol=1e-12,rtol=0)
    for name in ("synthetic","synthetic_chromosome_sensitivity"):
        comparison = pd.read_csv(args.out/(name+"_paired_contrasts.tsv"),sep="\t").set_index("arm")
        assert comparison.planned_family_contrasts.eq(2).all()
        np.testing.assert_allclose(comparison.loc["perfect",["MSE_improvement","MSE_improvement_low95","MSE_improvement_high95"]].to_numpy(dtype=float),1,atol=1e-12,rtol=0)
        np.testing.assert_allclose(comparison.loc["same_as_anchor",["MSE_improvement","MSE_improvement_low95","MSE_improvement_high95"]].to_numpy(dtype=float),0,atol=0,rtol=0)
    assert (args.out/"synthetic_matched_predictions.tsv.gz").is_file()
    assert not (args.out/"synthetic_chromosome_sensitivity_matched_predictions.tsv.gz").exists()
    result = dict(status="pass", rows=30, folds=5, chromosomes=10, historical_bins=20,
        identical_point_estimates=True, constant_error_difference_interval=True,
        identical_predictor_zero_interval=True, source_unchanged=True, no_duplicate_prediction_export=True,
        interpretation="checks_resampling_implementation; does_not_establish_22_independent_biological_participants",
        protected_outcomes_read=False)
    (args.out/"checks.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
