#!/usr/bin/env python3
"""Select a balanced 20-sample/four-cohort final-BAM LeafCutter pilot."""
import csv, os
ROOT=os.environ.get("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INV=os.path.join(ROOT,"GWAS/finemapping/results/seqfunc/disease_splicing/junction_sample_inventory.tsv")
OUT=os.path.join(ROOT,"GWAS/finemapping/results/seqfunc/disease_splicing/pilot");os.makedirs(OUT,exist_ok=True)
cohorts=["GSE126848","GSE130970","GSE162694","GSE213621"]
rows=list(csv.DictReader(open(INV),delimiter="\t")); selected=[]
for ds in cohorts:
    z=[r for r in rows if r["dataset"]==ds]
    c=sorted([r for r in z if r["group_binary"]=="Control"],key=lambda r:-int(r["unique_junction_reads"]))[:2]
    d=sorted([r for r in z if r["group_binary"]=="Disease"],key=lambda r:-int(r["unique_junction_reads"]))[:3]
    if len(c)!=2 or len(d)!=3: raise SystemExit(f"pilot balance failed for {ds}")
    selected.extend(c+d)
for i,r in enumerate(selected):
    bam=os.path.join(ROOT,"RNA-seq/Human/Patient_Cohorts/results",r["dataset"],"alignments/star",r["sample_id"],f'{r["sample_id"]}.Aligned.sortedByCoord.out.bam')
    if not os.path.exists(bam): raise SystemExit(f"missing {bam}")
    r["bam_path"]=bam; r["junc_path"]=os.path.join(OUT,"junctions",r["sample_id"]+".junc");r["task_id"]=str(i)
cols=["task_id","sample_id","dataset","group_binary","bam_path","sj_path","junc_path"]
with open(os.path.join(OUT,"pilot_manifest.tsv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=cols,delimiter="\t",extrasaction="ignore");w.writeheader();w.writerows(selected)
with open(os.path.join(OUT,"groups.tsv"),"w",newline="") as f:
    # LeafCutter 2.0.3 reads this file with header=None.
    w=csv.writer(f,delimiter="\t")
    for r in selected:w.writerow([r["sample_id"],r["group_binary"],r["dataset"]])
with open(os.path.join(OUT,"juncfiles.txt"),"w") as f:
    for r in selected:f.write(r["junc_path"]+"\n")
print(f"wrote {len(selected)} pilot samples to {OUT}")
