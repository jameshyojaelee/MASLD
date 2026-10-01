#!/usr/bin/env python3
"""Scientific invariants for nested donor/region comparisons; compute only."""
import argparse,json,os
from pathlib import Path
import numpy as np
import pandas as pd
import model_donor_common as C
import model_donor_nn as N
from model_donor_fit import selection_indices,baseline_pair,baseline_models,contrast_pairs


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):raise ValueError("Compute allocation required")
    rng=np.random.default_rng(20260915);x=rng.normal(size=(99,80));y=rng.normal(size=(99,128));seq=rng.normal(size=(128,12))
    fold=np.arange(99)%5
    reg=pd.DataFrame({"region_role":["train"]*64+["held"]*64,"inner_region_role":["inner_train"]*32+["inner_valid"]*32+["outer_held"]*64})
    for mode in ("trained_regions","held_regions"):
        tr,te,it,iv,rt,rv,irt,irv=selection_indices(fold,reg,0,mode)
        changed=y.copy();changed[te]+=1000
        if mode == "held_regions":changed[:,rv]-=2000
        ri=C.kernel_projection(x,it,np.ones(80,bool));ro=C.kernel_projection(x,tr,np.ones(80,bool))
        one=baseline_pair(y,seq,it,iv,irt,irv,tr,rt,mode);two=baseline_pair(changed,seq,it,iv,irt,irv,tr,rt,mode)
        bi,bo,si,so,_=one
        for a,b in zip(one[:2],two[:2]):np.testing.assert_array_equal(a,b)
        first,selection=baseline_models(y,ri,ro,si,so,bi,bo,it,iv,irt,irv,tr,rt,mode)
        second,selection2=baseline_models(changed,ri,ro,si,so,bi,bo,it,iv,irt,irv,tr,rt,mode)
        assert selection == selection2
        for key in first:np.testing.assert_allclose(first[key],second[key],rtol=0,atol=1e-12)
        if mode == "held_regions":assert not any(k in first for k in ("pcr","full_rna_ridge","pls_svd_ridge","selected_rna_linear"))
        else:assert "selected_rna_linear" in first
        for objective in ("profile","residual","pairwise"):
            original,_=N.train(ri["nn"],si["nn"],y,it,irt,bi,objective,True,(1,3))
            altered,_=N.train(ri["nn"],si["nn"],changed,it,irt,bi,objective,True,(1,3))
            for epoch in original:np.testing.assert_array_equal(original[epoch],altered[epoch])
            choice=lambda predictions:min(predictions,key=lambda e:np.mean((predictions[e][np.ix_(iv,irv)]-y[np.ix_(iv,irv)])**2))
            assert choice(original)==choice(altered)
    donors=np.flatnonzero(fold != 0)
    for epoch in (1,7,22):
        a,b,mask=N.balanced_partners(donors,epoch);a=a[mask>0];b=b[mask>0]
        assert np.all(a != b)
        for donor in donors:
            assert np.sum(a == donor)==2 and np.sum(b == donor)==2
            assert len(set(b[a == donor]))==2
    keep=np.ones(80,bool);keep[:20]=False
    p=C.kernel_projection(x,donors,keep);changed=x.copy();changed[:,:20]+=1e6
    q=C.kernel_projection(changed,donors,keep);np.testing.assert_array_equal(p["nn"],q["nn"])
    held=fold == 0;changed=x.copy();changed[held]+=1000
    q=C.kernel_projection(changed,donors,keep);np.testing.assert_allclose(p["nn"][donors],q["nn"][donors],atol=1e-10)
    gene_chrom=np.array(["unknown"]*10+["chr1"]*20+["chr2"]*50)
    masks={a:C.gene_mask(x,donors,gene_chrom,np.zeros(80),np.ones(80),reg,a,"chr1",np.array([0]))
        for a in ("all_rna","annotated_all_rna","chromosome","random_chromosome")}
    assert masks["all_rna"][0].all() and masks["annotated_all_rna"][0].sum()==70
    for a in ("annotated_all_rna","chromosome","random_chromosome"):assert not masks[a][0][:10].any()
    assert masks["chromosome"][0].sum()==masks["random_chromosome"][0].sum()==50
    keys=[a+"|film_profile" for a in masks]
    contrasts=contrast_pairs(keys,"held_regions")
    assert ("chromosome|film_profile","annotated_all_rna|film_profile","context_exclusion_MSE") in contrasts
    assert ("chromosome|film_profile","all_rna|film_profile","context_exclusion_MSE") not in contrasts
    # Product-kernel formula equals an explicit small Kronecker ridge system.
    xx=rng.normal(size=(5,4));ss=rng.normal(size=(3,2));kx=xx@xx.T;ks=ss@ss.T;yy=rng.normal(size=(5,3));c=.01
    actual=C.product_prediction({"kernel":kx,"training_donors":np.arange(5)},{"kernel":ks,"training_regions":np.arange(3)},yy,c,np.zeros(3))
    full=np.kron(kx,ks+1);penalty=c*np.linalg.eigvalsh(kx).max()*np.linalg.eigvalsh(ks+1).max()
    expected=(full@np.linalg.solve(full+penalty*np.eye(15),(yy-yy.mean(0)).ravel())).reshape(5,3)
    np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-10)
    result={"held_donor_and_region_outcomes_leave_fits_and_inner_selection_unchanged":True,
        "full_RNA_PCR_PLS_train_only":True,"stable_requested_rank_recipes":True,"region_specific_RNA_not_used_for_held_regions":True,
        "two_distinct_balanced_partners_per_training_donor":True,"removed_explicit_genes_cannot_enter_PCA":True,
        "held_donor_inputs_do_not_change_training_projection":True,"product_kernel_equals_direct_system":True,
        "exclusion_anchors_share_known_gene_universe":True,"random_exclusions_dimension_matched":True,
        "replication":"synthetic_invariants_only; no_scientific_performance_claim"}
    args.out.write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result,indent=2))


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--out",type=Path,required=True);main(p.parse_args())
