#!/usr/bin/env python3
"""Focused measurement, holdout and coordinate-selection checks; compute only."""
import argparse,json,os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import model_native_retention_analysis as A
import model_native_retention_panel as P


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):raise ValueError("Compute allocation required")
    rng=np.random.default_rng(20260915);x=rng.normal(size=(80,11));y=rng.normal(size=80);training=np.arange(40)
    first=A.calibrate(x,y,training);changed=y.copy();changed[40:]+=1e6;second=A.calibrate(x,changed,training)
    for a,b in zip(first,second):np.testing.assert_array_equal(a,b)
    for j in range(x.shape[1]):
        beta=np.linalg.lstsq(np.column_stack([np.ones(40),x[:40,j]]),y[:40],rcond=None)[0]
        np.testing.assert_allclose(first[0][:,j],beta[0]+beta[1]*x[:,j],atol=1e-12)
    inverted=-x[:,0];pred,slope,_,_=A.calibrate(x,inverted,training)
    assert slope[0]<0 and spearmanr(x[:,0],inverted).statistic<-.999 and spearmanr(pred[:,0],inverted).statistic>.999
    class Fasta:
        def fetch(self,chrom,start,end):return ("ACGT"*((end-start+3)//4))[:end-start]
    peaks=pd.DataFrame({"#chr":["chr1"]*600+["chr2"]*600+["chr3"]*20,
        "start":np.arange(1220)*3000+10000,"end":np.arange(1220)*3000+10330,"peakID":[f"peak{i}" for i in range(1220)]})
    one,_,_=P.select(peaks,peaks.peakID.to_numpy(),{"chr1":2,"chr2":1,"chr3":0},Fasta())
    altered=peaks.copy();altered["activity"]=rng.normal(size=len(altered))*1e9
    two,_,_=P.select(altered,peaks.peakID.to_numpy(),{"chr1":2,"chr2":1,"chr3":0},Fasta())
    pd.testing.assert_frame_equal(one,two)
    assert len(one)==1024 and not one.fold.eq(0).any() and one.input_component.nunique()==1024
    # Overlapping windows stay one resampling unit, independent of row order.
    geometry=pd.DataFrame({"chrom":["chr1","chr1","chr2"],"window_start0":[100,1100,100],"window_end0":[2148,3148,2148]})
    block=P.components(geometry);assert block[0]==block[1] and block[0]!=block[2]
    # Inspect the shared donor-resampled target reaching training-only calibration.
    panel=pd.DataFrame({"fold":[2]*20+[1]*20,"input_component":np.arange(40),"chrom":["chr1"]*20+["chr2"]*20})
    observed=rng.normal(size=(40,7));features=rng.normal(size=(40,11));captured=[];original=A.calibrate
    def capture(xx,yy,tr):captured.append(yy.copy());return original(xx,yy,tr)
    A.calibrate=capture
    try:A.bootstrap(features,observed,panel,np.arange(20),np.arange(20,40),"input_component",1)
    finally:A.calibrate=original
    draw=np.random.default_rng(A.SEED).choice(7,7,replace=True)
    np.testing.assert_array_equal(captured[0],observed[:,draw].mean(1))
    result={"calibration_unchanged_by_validation_outcomes":True,"affine_calibration_equals_independent_lstsq":True,
        "negative_slope_cannot_hide_raw_activity_order_reversal":True,"panel_selection_invariant_to_activity_columns":True,
        "overlapping_inputs_share_resampling_component":True,"same_donor_draw_calibration_and_validation":True,
        "all_ten_registered_adaptations_retained":len(A.ADAPTATION_NAMES)==10,"synthetic_only_no_performance_claim":True}
    args.out.write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result,indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--out",type=Path,required=True);main(p.parse_args())
