#!/usr/bin/env python3
"""Materialize exactly the original 128-region paired donor pilot."""
import argparse,json,os
from pathlib import Path
import numpy as np
import pandas as pd
import review_donor as D
from model_donor_common import region_splits


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):raise ValueError("Compute allocation required")
    args.out.mkdir(parents=True,exist_ok=False)
    ids,fold,rna,y,reg,seq,coords=D.load_substrate(2)
    expected=pd.read_csv(args.regions,sep="\t")
    pd.testing.assert_frame_equal(reg[expected.columns],expected,check_dtype=False)
    blocks,inner=region_splits(reg);reg["input_component"]=blocks;reg["inner_region_role"]=inner
    genes=pd.read_csv(D.FIX/"molecular/rna_feature_axis.tsv",sep="\t").stable_gene_id.to_numpy(dtype=str)
    x=D.logcpm(rna)
    np.savez_compressed(args.out/"paired.npz",participant_ids=ids.astype(str),donor_fold=fold,
        rna_log2cpm_full_library=x,target_log2cpm_full_library=y,dinucleotide_features=seq,
        gene_ids=genes,gene_chrom=coords[0].astype(str),gene_start=coords[1],gene_end=coords[2])
    reg.to_csv(args.out/"regions.tsv",sep="\t",index=False)
    (args.out/"receipt.json").write_text(json.dumps({"participants":99,"regions":128,"full_RNA_features":42163,
        "full_H3_library_regions":96460,"RNA_fractional_continuous":True,"normalization":"log2(1+CPM), full source library before every exclusion",
        "exclusion_interpretation":"removes_explicit_gene_features; excluded_genes_still_contribute_to_full_library_denominator",
        "sequence_semantics":"counted_interval_1bp_residual_not_called_peak_boundary",
        "no_H3_derived_predictor_covariates":True,"inner_donors":"next_established_outer_fold; remaining3folds_train",
        "inner_regions":"fixed_hash_half_of_training_tile_and_16kb_overlap_components; no_outerheld_regions_used",
        "all_arrays_MiB":sum(a.nbytes for a in (x,y,seq))/2**20},indent=2)+"\n")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--out",type=Path,required=True);p.add_argument("--regions",type=Path,required=True);main(p.parse_args())
