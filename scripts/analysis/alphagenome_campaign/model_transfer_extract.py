#!/usr/bin/env python3
"""Local 4,359-reporter representations and fully observed caQTL target pools.

Eligibility is source-native direct/swapped SNVs, reference-matched strict ACGT
sequence and whole measured-peak coverage. No significance filtering or
minimum-p-value target selection occurs. Features use local AlphaGenome;
reporter windows are native genomic flanks, not reconstructed vector sequence.
"""
import argparse
import json
from pathlib import Path
import sys
import time

import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
import pysam

from model_scalar import pool_weights
PROJ=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(PROJ/"scripts/analysis/alphagenome_program"))
import i1_common as C
from i1_extract_mpra_embeddings import load_trunk


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    columns=["element_id","contig","variant_pos1","genomic_ref","genomic_alt","reference_match_start0","reference_match_end0"]
    frame=pd.read_csv(args.reporters,sep="\t",usecols=columns).drop_duplicates()
    if frame.element_id.duplicated().any() or len(frame) != 4359:
        raise ValueError("Source reporter molecular identities are not unique/complete")
    folds=pd.read_csv(args.folds,sep="\t",usecols=["element_id","outer_fold","long_range_block_id"])
    frame=frame.merge(folds,on="element_id",validate="one_to_one")
    if len(frame) != 4359 or frame.groupby("long_range_block_id").outer_fold.nunique().max() != 1:
        raise ValueError("Source fold identity invalid")
    # A deterministic hash order spreads the first 100 timing examples across loci.
    import hashlib
    frame["timing_order"]=[hashlib.sha256(x.encode()).hexdigest() for x in frame.element_id]
    frame=frame.sort_values("timing_order").drop(columns="timing_order").reset_index(drop=True)
    frame.to_csv(args.out/"reporter_manifest.tsv",sep="\t",index=False)
    source=pd.read_csv(args.pairs,sep="\t")
    source=source.loc[source.eligible.astype(str).str.lower().eq("true")].copy()
    if source.duplicated(["element_id","peak"]).any():
        raise ValueError("Repeated variant-target identity requires explicit covariance handling")
    lo=source.variant_pos1-1-16384//2
    target_eligible=(source.peak_start0>=lo)&(source.peak_end0<=lo+16384)&(source.peak_end0>source.peak_start0)
    pairs=source.loc[target_eligible].copy().reset_index(drop=True)
    pairs["pair_index"]=np.arange(len(pairs))
    pairs.to_csv(args.out/"target_labels.tsv.gz",sep="\t",index=False)
    # Whole source blocks must also separate the actual 16-kb input intervals.
    for chrom,group in frame.groupby("contig"):
        ordered=group.sort_values("variant_pos1")
        pos=ordered.variant_pos1.to_numpy();fold=ordered.outer_fold.to_numpy()
        if np.any((np.diff(pos)<16384)&(fold[1:] != fold[:-1])):
            raise ValueError(f"Overlapping 16-kb windows span source folds on {chrom}")
    by_element={e:g for e,g in pairs.groupby("element_id",sort=False)}
    card=C.gpu_card()
    if card.get("card_tag") != "l40s":
        raise ValueError("L40S required")
    from alphagenome_research.model import one_hot_encoder
    encoder=one_hot_encoder.DNAOneHotEncoder()
    fasta=pysam.FastaFile(C.FASTA_PATH)
    trunk,base,state,device,load_seconds=load_trunk()
    org=jnp.zeros((1,),jnp.int32)
    complement=str.maketrans("ACGT","TGCA")
    rejected=[];timings={}
    for length in (2048,16384):
        variant=np.full((len(frame),4,3072),np.nan,dtype=np.float32)
        construct=np.full_like(variant,np.nan)
        target=np.full((len(pairs),4,3072),np.nan,dtype=np.float32) if length == 16384 else None
        keep=np.zeros(len(frame),dtype=bool);times=[]
        for i,row in enumerate(frame.itertuples(index=False)):
            tick=time.monotonic()
            lo,hi=C.window_bounds(row.variant_pos1,length)
            check=C.validate_window(fasta,row.contig,lo,hi,row.variant_pos1,row.genomic_ref)
            if not check["acgt_ok"] or not check["reference_match"]:
                rejected.append({"element_id":row.element_id,"length":length,"reason":"reference_or_ACGT"})
                continue
            ref,alt=C.extract_ref_alt(fasta,row.contig,lo,hi,row.variant_pos1,row.genomic_ref,row.genomic_alt)
            seqs=(ref,alt,ref.translate(complement)[::-1],alt.translate(complement)[::-1])
            cweights=pool_weights(length,"target",target_start=row.reference_match_start0,target_end=row.reference_match_end0,window_start=lo)
            targets=by_element.get(row.element_id)
            for allele,seq in enumerate(seqs):
                x=jnp.asarray(encoder.encode(seq)[None],jnp.float32)
                embedding=np.asarray(trunk(base,state,x,org).get_sequence_embeddings(128),dtype=np.float32)[0]
                if not np.isfinite(embedding).all():
                    raise ValueError("Nonfinite representation")
                reverse=allele>=2
                centre=(length-1-length//2) if reverse else length//2
                variant[i,allele]=embedding[centre//128]
                cw=cweights[::-1] if reverse else cweights
                construct[i,allele]=np.sum(embedding*cw[:,None],axis=0,dtype=np.float32)
                if target is not None and targets is not None:
                    for peak in targets.itertuples(index=False):
                        w=pool_weights(length,"target",target_start=peak.peak_start0,target_end=peak.peak_end0,window_start=lo)
                        w=w[::-1] if reverse else w
                        target[peak.pair_index,allele]=np.sum(embedding*w[:,None],axis=0,dtype=np.float32)
            keep[i]=True;times.append(time.monotonic()-tick)
            if len(times) == 100:
                # At each length, measured sequence + host pooling time plus
                # 2x remaining work must fit the one-hour allocation.
                remaining=2*float(np.mean(times[1:]))*(len(frame)-i-1)
                if length == 2048:
                    remaining += 4*float(np.mean(times[1:]))*len(frame)
                budget=args.max_hours*3600-(time.monotonic()-started)
                receipt={"length":length,"mean_warm_wall_seconds":float(np.mean(times[1:])),
                    "projected_remaining_seconds_with_buffer":remaining,"budget_remaining_seconds":budget}
                (args.out/f"timing_{length}.json").write_text(json.dumps(receipt,indent=2)+"\n")
                if remaining > budget-120:
                    raise RuntimeError("Measured extraction projection exceeds reserved allocation")
            if (i+1)%100 == 0:
                print(json.dumps({"length":length,"attempted":i+1,"total":len(frame),"mean_wall_seconds":float(np.mean(times[1:] or times))}),flush=True)
        np.savez_compressed(args.out/f"reporters_{length}.npz",element_ids=frame.element_id.to_numpy(dtype=str),
            allele_order=np.array(["REF","ALT","REF_RC","ALT_RC"]),pooled_centre_bin=variant,
            pooled_construct=construct,extracted=keep,outer_fold=frame.outer_fold.to_numpy(),
            long_range_block_id=frame.long_range_block_id.to_numpy(dtype=str),input_length_bp=length,gpu_card="l40s")
        if target is not None:
            np.savez_compressed(args.out/"target_pools_16384.npz",pair_index=pairs.pair_index.to_numpy(),
                                pooled_target=target,extracted=np.isfinite(target).all(axis=(1,2)))
        timings[length]={"rows_extracted":int(keep.sum()),"wall_seconds_per_row":float(np.mean(times[1:]))}
    pd.DataFrame(rejected,columns=["element_id","length","reason"]).to_csv(args.out/"excluded.tsv",sep="\t",index=False)
    (args.out/"complete.json").write_text(json.dumps({"reporter_population":len(frame),"source_eligible_pairs":len(source),
        "whole_target_eligible_pairs":len(pairs),"target_variants":pairs.element_id.nunique(),
        "target_blocks":pairs.borzoi_long_range_group_id.nunique(),"timings":timings,
        "elapsed_seconds":time.monotonic()-started,"restore_seconds":load_seconds,
        "overlapping_windows_cross_folds":0,"outcome_filtering":False,
        "strata":"fixed data/transfer_strata_spec.json; no outcome-driven revision",
        "source_restrictions":"internal_representations_not_source_data_redistribution",**card},indent=2)+"\n")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reporters",type=Path,required=True)
    parser.add_argument("--pairs",type=Path,required=True)
    parser.add_argument("--folds",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--max-hours",type=float,default=0.95)
    main(parser.parse_args())
