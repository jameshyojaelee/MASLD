#!/usr/bin/env python3
"""Outcome-free local AlphaGenome counted-interval sequence representations."""
import argparse,json,os,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
import jax,jax.numpy as jnp
import pysam
from model_scalar import pool_weights
from model_donor_common import region_splits
PROJ=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(PROJ/"scripts/analysis/alphagenome_program"))
import i1_common as C
from i1_extract_mpra_embeddings import load_trunk


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);start=time.monotonic();reg=pd.read_csv(args.regions,sep="\t")
    if len(reg) != args.expect:raise ValueError(f"Exactly{args.expect} fixed source regions required")
    if "tile_index" in reg:
        components,inner=region_splits(reg)
    else:
        # The learning-curve panel assigns and checks its own input-overlap components
        # before any outcome is read; recomputing them here would be a second, separate
        # definition of the same split. Re-check the guard on the supplied assignment.
        for column in ("input_component","region_role"):
            if column not in reg:raise ValueError(f"Region panel lacks {column}; cannot verify the outer split")
        components=reg.input_component.to_numpy()
        for block in set(components):
            if len(set(reg.region_role.to_numpy()[components == block])) != 1:
                raise ValueError("Actual input windows span outer region split")
        inner=reg.get("inner_region_role",pd.Series(["unassigned"]*len(reg))).to_numpy()
    if C.gpu_card().get("card_tag") != "l40s":raise ValueError("L40S required")
    length=2048;out=np.full((len(reg),2,3072),np.nan,np.float32)
    from alphagenome_research.model import one_hot_encoder
    encoder=one_hot_encoder.DNAOneHotEncoder();fasta=pysam.FastaFile(C.FASTA_PATH)
    trunk,base,state,device,restore=load_trunk();timings=[]
    complement=str.maketrans("ACGT","TGCA")
    for i,row in enumerate(reg.itertuples()):
        tick=time.monotonic();mid=(int(row.start0)+int(row.end0))//2;lo=mid-length//2
        sequence=fasta.fetch(row.chrom,lo,lo+length).upper()
        if len(sequence) != length or set(sequence)-set("ACGT"):raise ValueError(f"Source region {row.region_key} lacks completeACGT context; no silent population change")
        weights=pool_weights(length,"target",target_start=int(row.start0),target_end=int(row.end0),window_start=lo)
        for strand,seq in enumerate((sequence,sequence.translate(complement)[::-1])):
            x=jnp.asarray(encoder.encode(seq)[None],jnp.float32)
            emb=np.asarray(trunk(base,state,x,jnp.zeros(1,jnp.int32)).get_sequence_embeddings(128),dtype=np.float32)[0]
            w=weights if strand == 0 else weights[::-1]
            out[i,strand]=np.sum(emb*w[:,None],axis=0,dtype=np.float32)
        timings.append(time.monotonic()-tick)
        if not np.isfinite(out[i]).all():raise ValueError("Nonfinite representation")
    np.savez_compressed(args.out/"features.npz",region_key=reg.region_key.to_numpy(dtype=str),strand_features=out,mean_strands=out.mean(1,dtype=np.float32))
    (args.out/"receipt.json").write_text(json.dumps({"regions":len(reg),"length":2048,"strands":2,"pooled_width":3072,
        "input_overlap_components":int(len(set(components))),"component_source":"tile_index" if "tile_index" in reg else "supplied_input_component",
        "file_feature_MiB":out.nbytes/2**20,"seconds":time.monotonic()-start,"restore_seconds":restore,
        "mean_warm_region_seconds":float(np.mean(timings[1:])),"outer_and_inner_input_overlap_components_checked_at_lengths":[2048,16384],
        "whole_counted_interval_weighted128bp_pooling":True,"coordinate_residual_bp":1,"called_peak_boundary_claim":False,
        "sequence_backbone":"local_AlphaGenome_unchanged","outcomes_read":False,"jax":jax.__version__,**C.gpu_card()},indent=2)+"\n")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--regions",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    p.add_argument("--expect",type=int,default=128,help="Region count the panel must hold; 128 is the original fixture.")
    main(p.parse_args())
