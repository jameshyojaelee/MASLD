#!/usr/bin/env python3
"""Read-only fixed-backbone native liver-ATAC assessment on a fixed count panel."""
import argparse,hashlib,importlib.metadata,json,os,sys,time
from pathlib import Path
import jax,jax.numpy as jnp
import numpy as np
import pandas as pd
import pysam
import model_scalar as M
import model_native_retention_analysis as A

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts/analysis/alphagenome_program"))
import i1_common as C
from i1_extract_mpra_embeddings import load_trunk

TRACK_NAMES=("UBERON:0001114 ATAC-seq","UBERON:0001115 ATAC-seq","UBERON:0002107 ATAC-seq")
TRACK_INDICES=(129,130,145)
ARMS=A.ADAPTATION_NAMES


def emit(path,value):path.write_text(json.dumps(value,indent=2,default=str)+"\n")


def main(args):
    if not os.environ.get("SLURM_JOB_ID") or C.gpu_card().get("card_tag")!="l40s":raise ValueError("Allocated L40S required")
    args.out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    panel=pd.read_csv(args.fixture/"panel.tsv",sep="\t",dtype={"peak_id":str})
    receipt=json.loads((args.fixture/"targets_receipt.json").read_text())
    if len(panel)!=1024 or receipt["donors"]!=138 or receipt["full_peak_library_regions"]!=349685:raise ValueError("Registered source population differs")
    if C.sha256_file(args.fixture/"panel.tsv")!=receipt["panel_sha256"]:raise ValueError("Registered measured peak map changed after target extraction")
    if receipt["source_VST_or_normalization_factors_used"]:raise ValueError("Source-wide normalization/VST not allowed")
    if panel.peak_id.duplicated().any() or panel[["chrom","start0","end0"]].duplicated().any():raise ValueError("Duplicate measured region identity")
    if set(panel.loc[panel.role=="training","fold"])!={2,3,4} or set(panel.loc[panel.role=="validation","fold"])!={1}:raise ValueError("Role/fold mismatch")
    if panel.fold.eq(0).any() or panel.groupby("input_component").role.nunique().max()!=1:raise ValueError("Held fold0 or overlapping input windows admitted")
    # Reconstruct the assay from raw counts before inference; no prediction-based exclusions.
    measured,donors=A.targets(args.fixture,panel)
    from alphagenome.models import dna_model as api,dna_output
    from alphagenome_research.model import dna_model as local_model,one_hot_encoder
    from alphagenome_research.model.metadata import metadata
    md={o:metadata.load(o) for o in (api.Organism.HOMO_SAPIENS,api.Organism.MUS_MUSCULUS)}
    tracks=md[api.Organism.HOMO_SAPIENS].atac;selected=tracks.iloc[list(TRACK_INDICES)].copy()
    if tuple(selected.name)!=TRACK_NAMES:raise ValueError("Registered liver ATAC channel identities changed")
    for column,want in (("biosample_life_stage","adult"),("biosample_type","tissue"),("data_source","encode"),("strand","."),("Assay title","ATAC-seq")):
        if not selected[column].eq(want).all():raise ValueError(f"Liver ATAC metadata differs at {column}")
    if md[api.Organism.HOMO_SAPIENS].padding[dna_output.OutputType.ATAC][list(TRACK_INDICES)].any():raise ValueError("Padding track admitted")
    selected.insert(0,"native_channel_index",TRACK_INDICES);selected.to_csv(args.out/"native_liver_ATAC_tracks.tsv",sep="\t",index=False)
    encoder=one_hot_encoder.DNAOneHotEncoder();fasta=pysam.FastaFile(C.FASTA_PATH)
    sequences=[];masks=[]
    for row in panel.itertuples(index=False):
        sequence=fasta.fetch(row.chrom,int(row.window_start0),int(row.window_end0)).upper()
        if len(sequence)!=2048 or set(sequence)-set("ACGT") or hashlib.sha256(sequence.encode()).hexdigest()!=row.sequence_sha256:raise ValueError("Registered reference sequence changed")
        lo=int(row.start0-row.window_start0);hi=int(row.end0-row.window_start0)
        if not 0<=lo<hi<=2048:raise ValueError("Measured peak not wholly in input")
        mask=np.zeros(2048,np.float32);mask[lo:hi]=1
        if mask.sum()!=row.end0-row.start0:raise ValueError("1bp count interval integration mismatch")
        sequences.append(encoder.encode(sequence));masks.append(mask)
    sequences=np.asarray(sequences,np.float32);masks=np.asarray(masks,np.float32)
    _,base,state,_,restore_seconds=load_trunk();base_hash=M.tree_hash(base);state_hash=M.tree_hash(state)
    _,full_apply,_,_,_=local_model.create_model(md)
    @jax.jit
    def native(params,stats,x,mask):
        predictions=full_apply(params,stats,x,jnp.zeros(x.shape[0],jnp.int32))["atac"]["predictions_1bp"]
        values=jnp.take(predictions,jnp.asarray(TRACK_INDICES),axis=-1).astype(jnp.float32)
        return jnp.log2(1+jnp.sum(values*mask[:,:,None],axis=1,dtype=jnp.float32))
    configurations=[]
    for name in ARMS:
        root=args.sweep/name;weights=root/"weights.npz";saved=json.loads(Path(str(weights)+".json").read_text());feasibility=json.loads((root/"feasibility.json").read_text())
        config=saved["config"]
        if config["id"]!=name or config["loss"]!="mse" or config["length"]!=2048 or feasibility["completed_steps"]!=5000:raise ValueError("Prespecified complete MSE adaptation differs")
        if Path(saved["checkpoint"])!=C.CHECKPOINT:raise ValueError("Backbone checkpoint identity differs")
        paths=[tuple(p) for p in saved["paths"]]
        if paths!=M.qv_paths(base,config["last_blocks"]):raise ValueError("Adapted query/value paths differ")
        configurations.append({"name":name,"weights":str(weights),"weights_sha256":C.sha256_file(weights),"paths":paths,"mode":config["mode"],"actual_steps":feasibility["completed_steps"]})
    design={"status":"prespecified_measured_native_retention","models":["original"]+list(ARMS),"configuration_records":configurations,
        "sequence_length":2048,"orientation":"reference_forward_only","per_peak_integration":"exact_whole_BED0_interval_sum_of_1bp_native_ATAC_predictions",
        "raw_predictor":"mean_across_three_verified_adult_liver_tracks_of_log2_1plus_integrated_track_signal",
        "native_output_scale":"installed_heads.predict unscaled predictions_1bp using source_metadata_nonzero_mean; pooled_in_float32",
        "population_ascertainment":"source_consensus_counted_peaks; coordinate_and_width_selection_without_new_activity_thresholds",
        "calibration":"model_specific_intercept_slope_using512training_regions_only; no_hyperparameter_selection",
        "target":"donor_mean_log2_1plus_CPM_raw_WASP_filtered_read_counts; full_deposited_peak_library_totals",
        "training_donor_prediction":False,"foundation_source_exposure":"unresolved","experimental_correspondence":"adult_liver_tissue_ATAC; ENCODE_vs_Currin_not_identical_experiments",
        "no_DNase_or_whole_native_function_claim":True,"all_ten_MSE_adaptations_retained":True,
        "bootstrap":"same_donor_resample_in_calibration_and_validation_across_every_model; independent_training_validation_component_samples; affine_refit_only",
        "source_GC_correction_VST_offsets_used":False,"max_hours":args.max_hours,"bootstrap_draws":args.bootstrap,
        "environment":{"python":sys.version,"versions":{p:importlib.metadata.version(p) for p in ("numpy","scipy","pandas","jax","jaxlib","optax")},**C.gpu_card()},
        "panel_sha256":C.sha256_file(args.fixture/"panel.tsv"),"targets_receipt_sha256":C.sha256_file(args.fixture/"targets_receipt.json"),
        "restored_backbone_hash":base_hash,"restored_running_state_hash":state_hash,"checkpoint":str(C.CHECKPOINT)}
    emit(args.out/"design.json",design)
    names=["original"]+list(ARMS);outputs=np.full((len(names),1024,3),np.nan,np.float32);completed=[]
    try:
        for arm,name in enumerate(names):
            tick=time.monotonic()
            if arm==0:params=base
            else:
                config=configurations[arm-1];trainable=M.load_trainable(config["weights"])
                params=M.adapted_params(base,trainable,config["paths"],config["mode"])
                allowed={p[0] for p in config["paths"]}
                if any(params[k] is not base[k] for k in base if k not in allowed):raise ValueError("Native head or other frozen module changed")
            call_times=[]
            for i in range(1024):
                step=time.monotonic()
                values=np.asarray(native(params,state,jnp.asarray(sequences[i:i+1]),jnp.asarray(masks[i:i+1])),np.float32)[0]
                if not np.isfinite(values).all() or np.any(values<0):raise ValueError("Nonfinite or negative native ATAC logsum")
                outputs[arm,i]=values;call_times.append(time.monotonic()-step)
                if arm==0 and i==99:
                    warm=float(np.quantile(call_times[10:],.9));elapsed=time.monotonic()-started
                    # Analysis timing uses synthetic predictor columns; activity never selects a panel/model.
                    probe=np.column_stack([np.linspace(0,1,1024)**(1+j/20) for j in range(11)])
                    tr=np.flatnonzero(panel.role=="training");va=np.flatnonzero(panel.role=="validation")
                    bt=time.monotonic();A.bootstrap(probe,measured,panel,tr,va,"input_component",20);bootstrap20=time.monotonic()-bt
                    projected=2*(elapsed+len(names)*1024*warm+2*args.bootstrap*bootstrap20/20+120)
                    admission={"regions_measured":100,"wall_seconds_including_restore_and_checks":elapsed,"p90_host_seconds_per_region":warm,
                        "bootstrap20_seconds":bootstrap20,"projected_total_seconds_factor2":projected,"projected_GPU_hours":projected/3600,
                        "projected_allocated_CPU_core_hours":projected*4/3600,"admitted":projected<=args.max_hours*3600,
                        "includes":"all11backbones1024regions; raw_target_validation; loading; bootstrap_both_schemes;120s_IO_allowance;factor2"}
                    emit(args.out/"throughput_100regions.json",admission)
                    if not admission["admitted"]:raise RuntimeError("Native-retention projection exceeds reserved allocation")
                if time.monotonic()-started>args.max_hours*3600-300:raise RuntimeError("Native-retention wall guard reserves five minutes for records")
            repeated=np.asarray(native(params,state,jnp.asarray(sequences[:1]),jnp.asarray(masks[:1])),np.float32)[0]
            np.testing.assert_array_equal(repeated,outputs[arm,0])
            completed.append({"model":name,"regions":1024,"seconds":time.monotonic()-tick,"repeat_exact":True})
            emit(args.out/"completed_models.json",completed);print(json.dumps(completed[-1]),flush=True)
        if M.tree_hash(base)!=base_hash or M.tree_hash(state)!=state_hash:raise ValueError("Frozen backbone/state mutated")
        A.run(args.fixture,panel,outputs,names,args.out,args.bootstrap)
    except Exception as exc:
        emit(args.out/"failure.json",{"error":repr(exc),"completed_models":completed,"seconds":time.monotonic()-started});raise
    finally:
        np.savez_compressed(args.out/"native_track_predictions.npz",model_names=np.array(names),peak_id=panel.peak_id.to_numpy(str),native_track_indices=np.array(TRACK_INDICES),log2_track_sums=outputs)
    emit(args.out/"completion.json",{"status":"measured_regional_native_ATAC_comparison_complete","seconds":time.monotonic()-started,"models":11,"regions_per_model":1024,"restore_seconds":restore_seconds,"all_required_models_complete":True})


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--fixture",type=Path,required=True);p.add_argument("--sweep",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    p.add_argument("--max-hours",type=float,default=1.95);p.add_argument("--bootstrap",type=int,default=1000);main(p.parse_args())
