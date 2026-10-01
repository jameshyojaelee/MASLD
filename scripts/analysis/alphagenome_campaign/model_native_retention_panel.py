#!/usr/bin/env python3
"""Choose regional native-ATAC assessment inputs without reading activity values."""
import argparse,hashlib,json,os,sys
from pathlib import Path
import numpy as np
import pandas as pd
import pysam

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts/analysis/alphagenome_program"))
import i1_common as C
FOLDS=ROOT/"GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/fold_manifest.tsv"
SEED=20260915


def components(frame):
    result=np.full(len(frame),-1,int);component=-1
    for _,rows in frame.groupby("chrom",sort=True):
        end=-1
        for i,row in rows.sort_values("window_start0").iterrows():
            if row.window_start0>=end:component+=1
            result[i]=component;end=max(end,int(row.window_end0))
    return result


def select(peaks,source_ids,folds,fasta,n=1024):
    if n!=1024:raise ValueError("Fixed initial panel is1024 regions")
    if peaks.peakID.duplicated().any() or pd.Series(source_ids).duplicated().any():raise ValueError("Duplicate source peak identities")
    frame=peaks.loc[peaks.peakID.isin(source_ids),["#chr","start","end","peakID"]].rename(columns={"#chr":"chrom","start":"start0","end":"end0","peakID":"peak_id"}).copy()
    if len(frame)!=len(source_ids):raise ValueError("Not every count-matrix peak maps exactly to source BED")
    frame["fold"]=frame.chrom.map(folds)
    if frame.fold.isna().any():raise ValueError("Nonautosomal or unassigned count region")
    frame["width"]=frame.end0-frame.start0
    census={"source_count_rows":len(frame),"fold0_excluded":int((frame.fold==0).sum()),"nonpositive_width":int((frame.width<=0).sum())}
    frame=frame.loc[frame.fold.isin([1,2,3,4])&(frame.width>0)&(frame.width<=2048)].copy()
    census["whole_width_eligible_before_sequence"]=len(frame)
    frame["window_start0"]=(frame.start0+frame.end0)//2-1024;frame["window_end0"]=frame.window_start0+2048
    frame["coordinate_hash"]=[hashlib.sha256(f"native-retention|{SEED}|{c}:{s}-{e}|{p}".encode()).hexdigest()
        for c,s,e,p in frame[["chrom","start0","end0","peak_id"]].itertuples(index=False,name=None)]
    selected=[];rejected=[]
    for role,allowed in (("training",[2,3,4]),("validation",[1])):
        candidates=frame.loc[frame.fold.isin(allowed)].sort_values(["coordinate_hash","peak_id"])
        census[role+"_width_candidates"]=len(candidates);admitted=0
        for row in candidates.itertuples(index=False):
            sequence=fasta.fetch(row.chrom,max(0,int(row.window_start0)),max(0,int(row.window_end0))).upper()
            reason=None
            if row.window_start0<0 or len(sequence)!=2048:reason="incomplete_reference_window"
            elif set(sequence)-set("ACGT"):reason="non_ACGT_reference_window"
            if reason:rejected.append({"peak_id":row.peak_id,"role":role,"reason":reason});continue
            record=row._asdict();record["role"]=role;record["sequence_sha256"]=hashlib.sha256(sequence.encode()).hexdigest();selected.append(record)
            admitted+=1
            if admitted==n//2:break
        census[role+"_admitted"]=admitted
        if admitted!=n//2:raise ValueError(f"Only{admitted} eligible {role} rows; no outcome-selected replacement")
    panel=pd.DataFrame(selected).reset_index(drop=True);panel["input_component"]=components(panel)
    for _,rows in panel.groupby("input_component"):
        if rows.role.nunique()!=1:raise ValueError("Actual2kb input overlap across calibration and validation")
    census["sequence_rejections_before_fixed_panel"] = len(rejected)
    return panel,pd.DataFrame(rejected,columns=["peak_id","role","reason"]),census


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):raise ValueError("Compute allocation required")
    args.out.mkdir(parents=True,exist_ok=True)
    peaks=pd.read_csv(args.peaks,sep="\t");ids=pd.read_csv(args.source_axis,sep="\t").peak_id.to_numpy()
    fmap={}
    for row in pd.read_csv(FOLDS,sep="\t").itertuples(index=False):
        for chrom in row.test_chromosomes.split(","):
            if chrom in fmap:raise ValueError("Chromosome assigned twice")
            fmap[chrom]=int(row.fold)
    if set(fmap)!={f"chr{i}" for i in range(1,23)}:raise ValueError("Incomplete chromosome folds")
    panel,rejected,census=select(peaks,ids,fmap,pysam.FastaFile(C.FASTA_PATH))
    panel.to_csv(args.out/"panel.tsv",sep="\t",index=False);rejected.to_csv(args.out/"sequence_rejections.tsv",sep="\t",index=False)
    (args.out/"panel_recipe.json").write_text(json.dumps({"status":"coordinate_panel_fixed_before_activity_extraction","census":census,
        "selection":"ascending_SHA256_of_seed_coordinate_peakID; first512eligible_training_and512validation",
        "training_folds":[2,3,4],"validation_fold":1,"fold0_per_region_outcomes_exported":False,
        "reference":C.FASTA_PATH,"sequence_length":2048,"orientation":"reference_forward_only",
        "actual_input_overlap_components":int(panel.input_component.nunique()),"units":"whole_source_BED0_halfopen_peaks",
        "ascertainment":"coordinate_selected_within_source_ascertained_consensus_counted_peaks; not_genomewide_closed_regions",
        "count_values_accessed_by_panel_selector":False,"source_axis_sha256":C.sha256_file(args.source_axis),
        "peaks_sha256":C.sha256_file(args.peaks),"fold_manifest_sha256":C.sha256_file(FOLDS)},indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--peaks",type=Path,required=True);p.add_argument("--source-axis",type=Path,required=True);p.add_argument("--out",type=Path,required=True);main(p.parse_args())
