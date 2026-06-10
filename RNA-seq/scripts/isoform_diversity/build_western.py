#!/usr/bin/env python3
"""
Harvest the Western_Diet_Datasets for the Cas13 'Western' diet group (and the
GSE246088 HFD arm), mapping GSM->SRR and condition->group/cas13_group with the
M02b rules. Control genotype only (drops KO to avoid the genotype confound).

Outputs:
  western_driver.tsv  (dataset, sample_id=SRR, layout, fastq_r1, fastq_r2)
  western_samples.tsv (sample_id, dataset, group_binary, cas13_group, diet_model, sex)
"""
import csv, glob, json, os, re, sys, collections

PROJ = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WD   = f"{PROJ}/RNA-seq/Mouse/Western_Diet_Datasets"
ID   = f"{PROJ}/RNA-seq/scripts/isoform_diversity"

def gsm2srr(ds):
    m = {}
    p = f"{WD}/{ds}/metadata/gsm_to_srr.tsv"
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            m[r["gsm"]] = r["srr"]
    return m

def abundance_for(ds, srr):
    # reuse existing quant (same vM38 index, 278,396 targets); within-group
    # Control-vs-Disease contrast is strand-cancelling so unstranded is valid here
    p = f"{WD}/{ds}/quant/kallisto/{srr}/abundance.tsv"
    return p if os.path.exists(p) else None

# per-dataset: (gsm_col, condition->(_group_binary, cas13_group, diet_model) or None to drop)
def rule_gse220575(r):
    c = r["condition"]
    return {"Control": ("Control","Western","DIAMOND"),
            "MASH":    ("Disease","Western","DIAMOND")}.get(c)  # HCC -> drop
def rule_gse246088(r):
    if "Control" not in r["genotype"]:   # keep only the control genotype (drop KO)
        return None
    return {"Chow_Control": ("Control","Western","WD"),      # shared control
            "WD_Control":   ("Disease","Western","WD"),
            "HFD_Control":  ("Disease","HFD","HFD")}.get(r["condition"])
def rule_gse305484(r):
    if "KO" in (r.get("genotype","")):   # keep WT only if any KO present
        return None
    return {"Chow": ("Control","Western","WD_Fructose"),
            "WD":   ("Disease","Western","WD_Fructose")}.get(r["condition"])

DATASETS = {
    "GSE220575": ("gsm",          rule_gse220575),
    "GSE246088": ("geo_accession",rule_gse246088),
    "GSE305484": ("gsm_accession",rule_gse305484),
}

drv, smp, drop = [], [], collections.Counter()
for ds,(gsmcol,rule) in DATASETS.items():
    g2s = gsm2srr(ds)
    geno_vals = collections.Counter()
    for r in csv.DictReader(open(f"{WD}/{ds}/metadata/sample_metadata.csv")):
        geno_vals[r.get("genotype","")] += 1
    for r in csv.DictReader(open(f"{WD}/{ds}/metadata/sample_metadata.csv")):
        res = rule(r)
        gsm = r.get(gsmcol,"")
        srr = g2s.get(gsm,"")
        if res is None: drop[f"{ds}:drop"]+=1; continue
        if not srr:     drop[f"{ds}:no_srr"]+=1; continue
        ab = abundance_for(ds, srr)
        if not ab:      drop[f"{ds}:no_abundance"]+=1; continue
        gb, cg, dm = res
        drv.append((ds, srr, ab))
        smp.append((srr, ds, gb, cg, dm, r.get("sex","")))
    sys.stderr.write(f"[western] {ds} genotypes={dict(geno_vals)}\n")

with open(f"{ID}/western_driver.tsv","w") as f:
    f.write("dataset\tsample_id\tabundance_path\n")
    for x in drv: f.write("\t".join(x)+"\n")
with open(f"{ID}/western_samples.tsv","w") as f:
    f.write("sample_id\tdataset\tgroup_binary\tcas13_group\tdiet_model\tsex\n")
    for x in smp: f.write("\t".join(x)+"\n")

c = collections.Counter((s[3], s[2]) for s in smp)   # cas13_group x group_binary
sys.stderr.write(f"[western] {len(smp)} samples; cas13_group x group: {dict(c)}\n")
sys.stderr.write(f"[western] dropped: {dict(drop)}\n")
