#!/usr/bin/env python3
"""Audit held-out chromosome coverage before saturation scoring is submitted."""
import csv, glob, json, os
from collections import Counter
from pathlib import Path

ROOT=Path(os.environ.get("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
SEQ=ROOT/"GWAS/finemapping/results/seqfunc/haplotype_saturation/sequences"
ADULT=ROOT/"GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet"
OUT=ROOT/"GWAS/finemapping/results/seqfunc/haplotype_saturation/fold_preflight"
OUT.mkdir(parents=True,exist_ok=True)
counts=Counter()
for p in sorted(SEQ.glob("genomic_variants.shard???.tsv")):
    with p.open() as f:
        for x in csv.DictReader(f,delimiter="\t"): counts[x["chrom"]]+=1
folds=[]
for fp in sorted(ADULT.glob("hepatocyte_fold*/fold*.json")):
    fold=json.loads(fp.read_text()); base=fp.parent
    model=base/"chrombpnet_full/models/chrombpnet_nobias.h5"
    folds.append({"fold_id":base.name,"fold_json":str(fp),"model":str(model),
                  "model_exists":model.exists() and model.stat().st_size>0,"test":fold.get("test",[])})
rows=[]
for chrom,n in sorted(counts.items(),key=lambda x:(int(x[0].removeprefix("chr")) if x[0].removeprefix("chr").isdigit() else 99)):
    eligible=[f for f in folds if f["model_exists"] and chrom in f["test"]]
    rows.append({"chrom":chrom,"n_unique_genomic_snvs":n,"n_eligible_heldout_models":len(eligible),
      "eligible_fold_ids":";".join(f["fold_id"] for f in eligible),
      "status":"PASS" if eligible else "BLOCKED_no_heldout_model"})
with (OUT/"chromosome_fold_coverage.tsv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
missing=[x["chrom"] for x in rows if x["n_eligible_heldout_models"]==0]
contract={"status":"PASS" if not missing else "BLOCKED","n_models_found":sum(f["model_exists"] for f in folds),
 "n_scored_chromosomes":len(rows),"n_chromosomes_covered":len(rows)-len(missing),
 "missing_heldout_model_chromosomes":missing,"full_array_submission_allowed":not missing,
 "interpretation":"mechanistic_nomination_not_causal_evidence","models":folds}
(OUT/"fold_preflight_contract.json").write_text(json.dumps(contract,indent=2)+"\n")
print(json.dumps(contract,indent=2))
