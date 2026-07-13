#!/usr/bin/env python3
"""Inventory quantitative STAR junction inputs and enforce cohort-level analysis."""
import csv, glob, json, os
from collections import Counter, defaultdict
ROOT=os.environ.get("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT=os.path.join(ROOT,"GWAS/finemapping/results/seqfunc/disease_splicing"); os.makedirs(OUT,exist_ok=True)
meta_path=os.path.join(ROOT,"RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
meta={r["sample_id"]:r for r in csv.DictReader(open(meta_path))}
files=glob.glob(os.path.join(ROOT,"RNA-seq/Human/Patient_Cohorts/results/*/alignments/star/*/*SJ.out.tab"))
rows=[]
for p in files:
    parts=p.split(os.sep); dataset=parts[parts.index("results")+1]; sample=parts[parts.index("star")+1]
    m=meta.get(sample,{})
    # Fast integrity/coverage scan; STAR columns 7/8 are unique/multi-mapping reads.
    n=uniq=annot=0
    with open(p) as f:
        for line in f:
            x=line.rstrip().split("\t")
            if len(x)<9: continue
            n+=1; uniq+=int(x[6]); annot+=int(x[5])==1
    rows.append([sample,dataset,m.get("group_binary",""),m.get("condition",""),m.get("fibrosis_stage",""),p,n,uniq,annot])
with open(os.path.join(OUT,"junction_sample_inventory.tsv"),"w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["sample_id","dataset","group_binary","condition","fibrosis_stage","sj_path","n_junctions","unique_junction_reads","n_annotated_junctions"]); w.writerows(rows)
by=defaultdict(list)
for r in rows: by[r[1]].append(r)
summary=[]
for ds,z in sorted(by.items()):
    groups=Counter(r[2] for r in z); n_ctrl=groups.get("Control",0); n_dis=groups.get("Disease",0)
    summary.append([ds,len(z),n_ctrl,n_dis,sum(r[7] for r in z),n_ctrl>=10 and n_dis>=10,
                    "per-cohort disease_vs_control; meta-analyze junction effects"])
with open(os.path.join(OUT,"cohort_feasibility.tsv"),"w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["dataset","n_samples","n_control","n_disease","unique_junction_reads","primary_gate","design"]); w.writerows(summary)
with open(os.path.join(OUT,"contract.json"),"w") as f: json.dump({"status":"inventory_complete","apply_only_firewall":True,
  "source":"STAR pass1 SJ.out.tab; use unique-read counts and verify against final BAM in pilot",
  "primary":"cohort-specific intron-cluster differential usage followed by random-effects meta-analysis",
  "not_allowed":["raw-count pooling across cohorts","sample-level pseudoreplication","post-hoc expression effect filtering"],
  "blocker":"LeafCutter/regtools are not installed in rnaseq; isolated environment approval required."},f,indent=2)
print(f"wrote {len(rows)} sample inventories to {OUT}")
