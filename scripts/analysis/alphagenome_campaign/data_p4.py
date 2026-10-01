#!/usr/bin/env python3
"""Account for every P4 region and reduce a bounded, life-stage-separated pilot.

The new classification is an amendment informed by prior inspected metadata/results.
Historical groups are retained verbatim. Hosted P4 predictions remain interpretation
outputs and are not adaptation features. A presence census is not an HDF5 integrity
check; the pilot separately records every reduction failure.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
ATLAS=ROOT / "scripts/analysis/alphagenome_atlas"
SAT=ROOT / "GWAS/finemapping/results/alphagenome_atlas/p4-saturation-20260909T191917Z"
META=ROOT / "GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/raw/scorer_metadata"
sys.path.insert(0,str(ATLAS))
spec=importlib.util.spec_from_file_location("original_p4_features", ATLAS / "43_saturation_features.py")
OLD=importlib.util.module_from_spec(spec)
spec.loader.exec_module(OLD)
ORIGINAL_GROUP_COLUMNS=OLD.group_columns
STAGES=("adult_verified","developmental","life_stage_unresolved")
SEED=20260915


def life_stage(row):
    stage=str(row.get("biosample_life_stage","")).strip().lower()
    kind=str(row.get("biosample_type","")).strip().lower()
    if stage in {"embryonic","fetal","newborn","child","pediatric","juvenile","postnatal"}:
        return "developmental"
    if stage == "adult" and kind in {"tissue","primary_cell"}:
        return "adult_verified"
    return "life_stage_unresolved"


def group_columns(scorer, var):
    historical=ORIGINAL_GROUP_COLUMNS(scorer,var)
    n=len(next(iter(var.values()))) if var else 0
    rows=[{key:value[i] for key,value in var.items()} for i in range(n)]
    classes=[OLD.la.track_class(r.get("ontology_curie",""),r.get("biosample_name",""),r.get("biosample_type","")) for r in rows]
    stem={"ATAC":"liver_atac","DNASE":"liver_dnase","CHIP_HISTONE":"liver_h3k27ac","CAGE":"cage_liver","CHIP_TF":"chip_tf_liver"}.get(scorer)
    if stem:
        for stage in STAGES:
            historical[f"{stem}_{stage}"]=np.asarray([i for i,r in enumerate(rows) if classes[i] in {"primary_liver","hepatocyte"} and life_stage(r)==stage and (scorer != "CHIP_HISTONE" or str(r.get("histone_mark","")).upper()=="H3K27AC")],dtype=int)
    return historical


def write_rows(path, rows, columns=None):
    if path.exists():
        raise FileExistsError(path)
    if columns is None:
        columns=list(dict.fromkeys(key for row in rows for key in row))
    opener=gzip.open if str(path).endswith(".gz") else open
    with opener(path,"wt",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter="\t",lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def track_registry(out):
    rows=[]
    for scorer in OLD.SCORERS:
        path=META / f"{scorer}.track_metadata.tsv"
        with path.open() as handle:
            tracks=list(csv.DictReader(handle,delimiter="\t"))
        var={k:[r[k] for r in tracks] for k in tracks[0]} if tracks else {}
        groups=group_columns(scorer,var)
        for i,row in enumerate(tracks):
            selected=[group for group,idx in groups.items() if i in idx]
            if selected:
                rows.append({"scorer":scorer,**row,"new_life_stage_group":life_stage(row),"groups":";".join(selected),"metadata_source":str(path.relative_to(ROOT))})
    write_rows(out / "track_life_stage_registry.tsv",rows)
    return Counter((r["scorer"],group) for r in rows for group in r["groups"].split(";"))


def census(out):
    with gzip.open(SAT / "tables/saturation_universe.tsv.gz","rt") as handle:
        original=list(csv.DictReader(handle,delimiter="\t"))
    merged={}
    for row in original:
        key=(row["universe"],row["chrom"],int(row["start0"]),int(row["end"]))
        if key not in merged:
            merged[key]={**row,"source_universe_rows":1,"lineages":{row["lineage"]},"da_states":{row["u2a_da_state"]} - {""}}
        else:
            merged[key]["source_universe_rows"]+=1
            merged[key]["lineages"].add(row["lineage"])
            merged[key]["da_states"].add(row["u2a_da_state"])
            for flag in ("u2b_sig_either","u2c_genetic_overlap","u2d_program_linked"):
                merged[key][flag]=str(str(merged[key][flag]).lower()=="true" or str(row[flag]).lower()=="true")
    rows=[]
    counts=defaultdict(Counter)
    expected={f"{s}.h5ad" for s in OLD.SCORERS} | {"request.json"}
    started=time.time()
    for i,(key,row) in enumerate(merged.items()):
        universe,chrom,start,end=key
        directory=OLD.region_dir(SAT / "raw/saturation",universe,chrom,start,end)
        missing=sorted(expected)
        state="directory_missing"
        size=0
        if directory.is_dir():
            with os.scandir(directory) as entries:
                files={entry.name:entry for entry in entries if entry.name in expected}
            missing=sorted(expected - files.keys())
            empty=[name for name,entry in files.items() if entry.stat().st_size == 0]
            size=sum(entry.stat().st_size for entry in files.values())
            missing+= [f"empty:{name}" for name in empty]
            state="present_not_integrity_validated" if not missing else "archive_incomplete"
        row={**row,"lineages":";".join(sorted(row["lineages"])),"da_states":";".join(sorted(row["da_states"]-{""})),"archive_state":state,"missing_members":";".join(missing),"archive_bytes":size}
        rows.append(row)
        counts[universe][state]+=1
        if i and i%20000==0:
            print(f"census {i}/{len(merged)} elapsed={time.time()-started:.1f}s",flush=True)
    write_rows(out / "region_completeness.tsv.gz",rows)
    return rows,{k:dict(v) for k,v in counts.items()},len(original)


def reduce_pilot(out,rows,max_regions):
    available=[r for r in rows if r["archive_state"]=="present_not_integrity_validated"]
    selected=np.random.default_rng(SEED).choice(len(available),size=min(max_regions,len(available)),replace=False)
    OLD.group_columns=group_columns
    features,tf,states=[],[],Counter()
    started=time.time()
    for i in selected:
        r=available[int(i)]
        task=(str(SAT / "raw/saturation"),r["region_key"],r["universe"],r["chrom"],int(r["start0"]),int(r["end"]))
        row,top,state=OLD.process_region(task)
        row["archive_state"]=state
        features.append(row)
        tf.extend(top)
        states[state]+=1
    write_rows(out / "PILOT_region_features.tsv.gz",features)
    write_rows(out / "PILOT_tf_top5.tsv.gz",tf)
    elapsed=time.time()-started
    return {"regions":len(features),"states":dict(states),"wall_seconds":elapsed,"seconds_per_region":elapsed/max(len(features),1),"sampling":"seeded random among present archives; not full-U1 or full-U2 results","seed":SEED}


def family_dispositions(rows,counts):
    total_u1=sum(counts.get("U1_h3k27ac",{}).values())
    u1_ok=counts.get("U1_h3k27ac",{}).get("present_not_integrity_validated",0)==total_u1 and total_u1==96460
    u2=[r for r in rows if r["universe"] != "U1_h3k27ac"]
    u2_ok=len(u2)==6930 and all(r["archive_state"]=="present_not_integrity_validated" for r in u2)
    supported=sum("supported" in r["da_states"].split(";") for r in u2)
    source_dependent=sum("source_dependent" in r["da_states"].split(";") for r in u2)
    genetic=sum(str(r["u2c_genetic_overlap"]).lower()=="true" for r in u2)
    families=[]
    for family,ready in (("F1_predictability",u1_ok),("F2_positional",u1_ok),("F3_TF_composition",u1_ok),("F4_DA_grammar",u2_ok),("F5_program_profiles",u2_ok)):
        reason="requires_full_integrity_validated_life_stage_reduction_and_minimum_group_check" if ready else "retrieval_incomplete; full_family_not_tested"
        if family=="F4_DA_grammar":
            reason+=f"; supported={supported}; source_dependent={source_dependent}; supported_arm_under_200_is_indeterminate"
        if family=="F5_program_profiles":
            reason+="; requires_original_program_gene_links_and_size_matched_random_gene_set_null"
        families.append({"family":family,"status":"not_run","reason":reason,"minimum_regions_each_group":200,"genetic_overlap_peaks":genetic})
    return families


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",required=True,type=Path)
    parser.add_argument("--max-regions",type=int,default=100)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    track_counts=track_registry(args.out)
    rows,counts,n_source=census(args.out)
    pilot=reduce_pilot(args.out,rows,args.max_regions)
    dispositions=family_dispositions(rows,counts)
    write_rows(args.out / "F1_F5_dispositions.tsv",dispositions)
    summary={"snapshot_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),"source_universe_rows":n_source,"unique_regions":len(rows),"archive_presence":counts,"pilot":pilot,"track_counts":{f"{k[0]}:{k[1]}":v for k,v in track_counts.items()},"amendment":"life_stage_split_informed_by_previously_inspected_metadata_and_results; historical_groups_retained","hosted_features_training_allowed":False,"full_family_results":False,"retrieval_active_snapshot":True}
    (args.out / "summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":
    main()
