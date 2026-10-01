#!/usr/bin/env python3
"""Check shared reporter identities against archived 2-kb four-allele features."""
import argparse
import json
from pathlib import Path
import numpy as np


def main(args):
    with np.load(args.archived,allow_pickle=False) as old, np.load(args.current,allow_pickle=False) as new:
        if not np.array_equal(old["allele_order"],new["allele_order"]):
            raise ValueError("Allele representation order differs")
        oi=old["element_ids"].astype(str);ni=new["element_ids"].astype(str)
        lookup={x:i for i,x in enumerate(ni)}
        if not set(oi).issubset(lookup):
            raise ValueError("Some archived reporter identities are absent")
        order=np.array([lookup[x] for x in oi])
        a=old["pooled_centre_bin"];b=new["pooled_centre_bin"][order]
        keep=old["extracted"]&new["extracted"][order]
        a=a[keep].astype(float);b=b[keep].astype(float)
        difference=b-a
        old_delta=0.5*((a[:,1]-a[:,0])+(a[:,3]-a[:,2]))
        new_delta=0.5*((b[:,1]-b[:,0])+(b[:,3]-b[:,2]))
        delta_difference=new_delta-old_delta
        result={"archived":str(args.archived),"current":str(args.current),"common_identities":len(oi),
            "both_extracted":int(keep.sum()),"allele_order":old["allele_order"].tolist(),
            "maximum_absolute_difference":float(np.max(np.abs(difference))),
            "mean_absolute_difference":float(np.mean(np.abs(difference))),
            "relative_L2_difference":float(np.linalg.norm(difference.ravel())/max(np.linalg.norm(a.ravel()),1e-20)),
            "bitwise_equal":bool(np.array_equal(a,b)),
            "numerically_equal_rtol1e-5_atol5e-5":bool(np.allclose(a,b,rtol=1e-5,atol=5e-5)),
            "differences_per_allele":[float(np.max(np.abs(difference[:,i]))) for i in range(4)],
            "signed_effect_feature_maximum_absolute_difference":float(np.max(np.abs(delta_difference))),
            "signed_effect_feature_relative_L2_difference":float(np.linalg.norm(delta_difference.ravel())/max(np.linalg.norm(old_delta.ravel()),1e-20)),
            "signed_effect_feature_global_correlation":float(np.corrcoef(old_delta.ravel(),new_delta.ravel())[0,1]),
            "purpose":"shared_reporter_input_and_variant_bin_equivalence_check; no_fitting"}
    args.out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archived",type=Path,required=True)
    parser.add_argument("--current",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    main(parser.parse_args())
