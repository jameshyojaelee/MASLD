#!/usr/bin/env python3
"""Output-only zero-adapter agreement across independent compiled routes.

No labels or optimizer steps. A common already-trained fixed head is applied
to 32 position-thinned development sequences. Repeated and cache-cleared
compilations establish an observed repeat floor; its size is reported rather
than hidden in a permissive tolerance. Production scalar-only and diagnostic
per-allele routes are measured separately because extra outputs can alter XLA.
"""
import argparse
import json
import os
from pathlib import Path
import sys
import time

import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
import pysam
import model_scalar as M

PROJ=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(PROJ/"scripts/analysis/alphagenome_program"))
import i1_common as C
from i1_extract_mpra_embeddings import load_trunk


def difference(a,b):
    x=np.asarray(a,dtype=float);y=np.asarray(b,dtype=float);d=y-x
    return {"max_absolute":float(np.max(np.abs(d))),"RMS_absolute":float(np.sqrt(np.mean(d*d))),
        "relative_L2":float(np.linalg.norm(d.ravel())/max(np.linalg.norm(x.ravel()),1e-20)),
        "exact":bool(np.array_equal(x,y))}


def main(args):
    args.out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    jax.config.update("jax_enable_compilation_cache",False)
    card=C.gpu_card()
    if card.get("card_tag") != "l40s" or not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Allocated L40S required")
    design={"variants":32,"microbatch":4,"length":2048,"poolings":["variant","symmetric"],
        "ranks":[4,16],"last_blocks":[3,5],"training_steps":0,
        "sequence_selection":"even_position_indices_on_existing_manifest_folds2_to4; no_outcome_columns_read",
        "head":str(args.head),"tolerance":"measured_repeat_max_absolute_plus_8_float32_eps_times_max1_reference_magnitude",
        "strict_tolerance":"8_float32_eps_times_max1_reference_magnitude; reported_separately_from_repeat_floor",
        "repeats":"same_compiled_call_twice_plus_clear_jax_caches_and_recompile_each_route",
        "persistent_compilation_cache_enabled":jax.config.jax_enable_compilation_cache,
        "interpretation":"nonzero_repeat_floor_is_not_proof_small_adaptation_gain_exceeds_numerical_variation"}
    (args.out/"design.json").write_text(json.dumps(design,indent=2)+"\n")
    frame=pd.read_csv(args.manifest,sep="\t")
    if "beta_alt" in frame:
        raise ValueError("Use outcome-free sequence_manifest.tsv")
    frame=frame.loc[frame.heldout_fold.isin([2,3,4])].reset_index(drop=True)
    selected=frame.iloc[np.linspace(0,len(frame)-1,32,dtype=int)].copy()
    fasta=pysam.FastaFile(C.FASTA_PATH)
    from alphagenome_research.model import one_hot_encoder
    encoder=one_hot_encoder.DNAOneHotEncoder();refs=[];alts=[]
    for row in selected.itertuples(index=False):
        lo,hi=C.window_bounds(int(row.pos_hg38),2048)
        check=C.validate_window(fasta,row.chr,lo,hi,int(row.pos_hg38),row.ref)
        if not check["acgt_ok"] or not check["reference_match"]:
            raise ValueError("Position-thinned example ineligible; do not replace silently")
        ref,alt=C.extract_ref_alt(fasta,row.chr,lo,hi,int(row.pos_hg38),row.ref,row.alt)
        refs.append(encoder.encode(ref));alts.append(encoder.encode(alt))
    refs=jnp.asarray(refs,jnp.float32);alts=jnp.asarray(alts,jnp.float32)
    selected.to_csv(args.out/"sequence_examples.tsv",sep="\t",index=False)
    trunk,base,state,device,restore_seconds=load_trunk()
    common_head=M.load_trainable(args.head)["head"]
    base_hash=M.tree_hash(base);state_hash=M.tree_hash(state)
    configurations=[("fixed","frozen",4,3)]+[(f"adapter_r{rank}_last{last}","adapter",rank,last) for rank in (4,16) for last in (3,5)]
    outputs={};repeats=[];rows=[]
    for name,mode,rank,last in configurations:
        trainable,paths=M.initialize(base,mode=mode,rank=rank,last_blocks=last)
        trainable["head"]=common_head
        adapted=M.adapted_params(base,trainable,paths,mode)
        for module,leaf in paths:
            np.testing.assert_array_equal(np.asarray(adapted[module][leaf]),np.asarray(base[module][leaf]))
        del adapted
        def scalar(t,p,s,r,a,w):
            return M.sequence_effect(trunk,p,s,t,paths,mode,r,a,w)
        def diagnostic(t,p,s,r,a,w):
            params=M.adapted_params(p,t,paths,mode);n=len(r)
            embeddings=trunk(params,s,jnp.concatenate([r,a]),jnp.zeros(2*n,dtype=jnp.int32)).get_sequence_embeddings(128).astype(jnp.float32)
            pooled=jnp.sum(embeddings*jnp.concatenate([w,w])[...,None],axis=1,dtype=jnp.float32)
            ref_score=M.shared_head(t["head"],pooled[:n]);alt_score=M.shared_head(t["head"],pooled[n:])
            return {"diagnostic_effect":alt_score-ref_score,"ref_score":ref_score,"alt_score":alt_score,
                    "ref_pooled":pooled[:n],"alt_pooled":pooled[n:],"signed_pooled":pooled[n:]-pooled[:n]}
        scalar_jit=jax.jit(scalar);diagnostic_jit=jax.jit(diagnostic)
        def collect():
            result={}
            for pooling in ("variant","symmetric"):
                w=jnp.asarray(np.repeat(M.pool_weights(2048,pooling)[None],4,axis=0))
                chunks=[]
                for i in range(0,32,4):
                    r,a=refs[i:i+4],alts[i:i+4]
                    prod=np.asarray(scalar_jit(trainable,base,state,r,a,w))
                    diagnostic_values={k:np.asarray(v) for k,v in diagnostic_jit(trainable,base,state,r,a,w).items()}
                    if not np.isfinite(prod).all() or any(not np.isfinite(v).all() for v in diagnostic_values.values()):
                        raise ValueError("Nonfinite output in repeatability/route comparison")
                    chunks.append({"production_effect":prod,**diagnostic_values})
                result[pooling]={key:np.concatenate([c[key] for c in chunks]) for key in chunks[0]}
            return result
        first=collect();same=collect()
        jax.clear_caches()
        recompiled=collect()
        outputs[name]=first
        for pooling in first:
            for output in first[pooling]:
                repeat=difference(first[pooling][output],same[pooling][output])
                recomp=difference(first[pooling][output],recompiled[pooling][output])
                repeats.append({"route":name,"pooling":pooling,"output":output,
                    "same_compilation":repeat,"fresh_compilation":recomp,
                    "repeat_floor_max_absolute":max(repeat["max_absolute"],recomp["max_absolute"])})
            for i,key in enumerate(selected.key):
                rows.append({"route":name,"pooling":pooling,"key":key,
                    "production_effect":float(first[pooling]["production_effect"][i]),
                    "diagnostic_effect":float(first[pooling]["diagnostic_effect"][i]),
                    "ref_score":float(first[pooling]["ref_score"][i]),"alt_score":float(first[pooling]["alt_score"][i])})
        print(json.dumps({"route_completed":name,"elapsed_seconds":time.monotonic()-start}),flush=True)
        if time.monotonic()-start > 13*60:
            raise RuntimeError("Output audit reached reserved wall-clock guard")
    contrasts=[];floor={(r["route"],r["pooling"],r["output"]):r["repeat_floor_max_absolute"] for r in repeats}
    for name in outputs:
        if name == "fixed":continue
        for pooling in outputs[name]:
            for output in outputs[name][pooling]:
                reference=outputs["fixed"][pooling][output]
                observed=difference(reference,outputs[name][pooling][output])
                strict=float(8*np.finfo(np.float32).eps*max(1,float(np.max(np.abs(reference)))))
                repeat_floor=max(floor[("fixed",pooling,output)],floor[(name,pooling,output)])
                contrasts.append({"route":name,"pooling":pooling,"output":output,**observed,
                    "strict_float32_tolerance":strict,"repeat_floor":repeat_floor,"repeat_based_tolerance":repeat_floor+strict,
                    "strict_output_agreement":observed["max_absolute"] <= strict,
                    "within_measured_repeat_floor":observed["max_absolute"] <= repeat_floor+strict})
    np.testing.assert_equal(M.tree_hash(base),base_hash);np.testing.assert_equal(M.tree_hash(state),state_hash)
    (args.out/"repeatability.json").write_text(json.dumps(repeats,indent=2)+"\n")
    pd.DataFrame(contrasts).to_csv(args.out/"route_agreement.tsv",sep="\t",index=False)
    pd.DataFrame(rows).to_csv(args.out/"effect_outputs.tsv",sep="\t",index=False)
    receipt={"elapsed_seconds":time.monotonic()-start,"restore_seconds":restore_seconds,"variants":32,
        "routes":len(configurations),"training_steps":0,"labels_read":False,"frozen_parameters_state_unchanged":True,
        "all_strict_output_agreement":all(r["strict_output_agreement"] for r in contrasts),
        "all_within_measured_repeat_floor":all(r["within_measured_repeat_floor"] for r in contrasts),
        "nonzero_repeat_floor":any(r["repeat_floor_max_absolute"]>0 for r in repeats),
        "jax":jax.__version__,"numpy":np.__version__,**card}
    (args.out/"complete.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--head",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    main(parser.parse_args())
