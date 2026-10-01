#!/usr/bin/env python3
"""Independent unique-pair direction reconstruction from the unchanged table."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/"GWAS/finemapping/results/alphagenome_campaign/week1-20260915"
SOURCE=ROOT/"GWAS/finemapping/results/alphagenome_atlas/p3a-splice-direction-20260915T141415Z/tables/splice_direction_pairs.tsv"


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Compute-node check required")
    args.out.mkdir(parents=True,exist_ok=False)
    original=pd.read_csv(SOURCE,sep="\t")
    rows=original.loc[original.state.eq("scored")].drop_duplicates().reset_index(drop=True)
    assert len(original)==657 and original.state.eq("scored").sum()==637 and len(rows)==621
    assert not rows.duplicated(["variant_uid","phenotype_id"]).any()
    chromosome=rows.variant_uid.str.split(":").str[0].to_numpy()
    position=rows.variant_uid.str.split(":").str[1].astype(int).to_numpy()
    np.testing.assert_array_equal(position,rows.position)
    units=np.full(len(rows),"",dtype=object)
    for chrom in np.unique(chromosome):
        indices=np.flatnonzero(chromosome==chrom)
        indices=indices[np.argsort(position[indices],kind="stable")]
        component=np.cumsum(np.r_[True,np.diff(position[indices])>1_000_000])
        units[indices]=[f"{chrom}:{number}" for number in component]
    measured=rows.slope.to_numpy(float)
    predicted=rows.excision_ratio_log2.to_numpy(float)
    assert np.isfinite(measured).all() and np.isfinite(predicted).all() and np.all(predicted!=0)
    agree=(measured>0)==(predicted>0)
    keys=pd.unique(units)
    assert len(keys)==247
    size=np.array([(units==key).sum() for key in keys])
    success=np.array([agree[units==key].sum() for key in keys])
    draw=np.random.default_rng(123).choice(len(keys),size=(2000,len(keys)),replace=True)
    interval=np.quantile(success[draw].sum(1)/size[draw].sum(1),[.025,.975])
    balance=2*success-size
    p=binomtest(np.sum(balance>0),np.sum(balance!=0),.5).pvalue
    pp,pm=(predicted>0).mean(),(measured>0).mean()
    marginal=pp*pm+(1-pp)*(1-pm)
    rho=spearmanr(np.abs(measured),np.abs(predicted)).statistic
    reviewed=json.loads((BASE/"review/splice-dedup-21773602/comparison.json").read_text())["exact_unique_rows"]
    expected=reviewed["direction"]["excision_ratio_log2"]
    actual=[agree.mean(),marginal,*interval,p,rho]
    np.testing.assert_allclose(actual,[expected["concordance"],expected["marginal_expectation"],expected["block_lo"],expected["block_hi"],
        expected["block_binomial_p"],reviewed["magnitude"]["spearman_abs_slope_vs_abs_ratio"]],rtol=0,atol=1e-12)
    report=dict(status="pass",raw_scored_rows=637,unique_scored_pairs=621,blocks=247,removed_exact_extra_copies=16,
        concordance=float(agree.mean()),marginal_expectation=float(marginal),interval95=interval.tolist(),
        original_nominal_block_majority_p=float(p),magnitude_spearman=float(rho),
        interpretation="weak_direction_below_original_ranking_criteria; loose_event_definition_and_pretraining_exposure_remain",
        old_source_unchanged=True,new_inference=False,protected_outcomes_read=False,
        independently_checked="direction_marginal_expectation_block_bootstrap_p_magnitude; detection_not_reconstructed_here")
    (args.out/"checks.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    main(parser.parse_args())
