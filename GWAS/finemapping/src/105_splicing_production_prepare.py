#!/usr/bin/env python3
"""Prepare five cohort-specific STAR-SJ to LeafCutter production manifests."""
import csv, os
ROOT=os.environ.get("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INV=os.path.join(ROOT,"GWAS/finemapping/results/seqfunc/disease_splicing/junction_sample_inventory.tsv")
OUT=os.path.join(ROOT,"GWAS/finemapping/results/seqfunc/disease_splicing/production");os.makedirs(OUT,exist_ok=True)
cohorts=["GSE126848","GSE130970","GSE135251","GSE162694","GSE213621"]
rows=[r for r in csv.DictReader(open(INV),delimiter="\t") if r["dataset"] in cohorts and r["group_binary"] in {"Control","Disease"}]
for i,r in enumerate(rows):
    d=os.path.join(OUT,r["dataset"]);os.makedirs(os.path.join(d,"junctions"),exist_ok=True)
    r["task_id"]=str(i);r["junc_path"]=os.path.join(d,"junctions",r["sample_id"]+".junc")
cols=["task_id","sample_id","dataset","group_binary","sj_path","junc_path"]
with open(os.path.join(OUT,"production_manifest.tsv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=cols,delimiter="\t",extrasaction="ignore",lineterminator="\n");w.writeheader();w.writerows(rows)
for ds in cohorts:
    z=[r for r in rows if r["dataset"]==ds]; d=os.path.join(OUT,ds)
    with open(os.path.join(d,"groups.tsv"),"w",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        for r in z:w.writerow([r["sample_id"],r["group_binary"]])
    with open(os.path.join(d,"juncfiles.txt"),"w") as f:
        for r in z:f.write(r["junc_path"]+"\n")
print(f"wrote {len(rows)} samples across {len(cohorts)} primary cohorts")
