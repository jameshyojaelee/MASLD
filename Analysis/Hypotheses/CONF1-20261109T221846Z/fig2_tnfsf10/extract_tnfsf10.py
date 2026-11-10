"""Upstream rows for the TNFSF10 locus panel (UKBB_ALT, chr3, hg19). Same sources and filters as
PANEL-DEEP-20261017T022400Z/fig2/recompute/extract_upstream.py, restricted to one window."""
import os, subprocess
import numpy as np, pandas as pd
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; FM = f"{ROOT}/GWAS/finemapping"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "upstream")
C, A, B = "3", 171_900_000, 172_700_000
tier = pd.read_csv(f"{FM}/config/gwas_trait_tier.tsv", sep="\t"); reg = pd.read_csv(f"{FM}/config/gwas_registry.tsv", sep="\t")
t12 = tier[tier.tier.isin([1, 2])].merge(reg[["study_name", "sumstats_path"]], on="study_name"); studies = set(t12.study_name)
for fname in ["combined_finemapping.csv", "combined_finemapping_coloc_targeted.csv"]:
    keep = []
    for ch in pd.read_csv(f"{FM}/results/{fname}", chunksize=2_000_000, dtype={"chromosome": str}, low_memory=False):
        ch = ch[ch.study.isin(studies) & (ch.chromosome.astype(str) == C) & (ch.position >= A) & (ch.position <= B)]
        keep.append(ch)
    pd.concat(keep).to_csv(f"{OUT}/{fname.replace('.csv','')}_GH_windows.tsv.gz", sep="\t", index=False)
    print(fname, "done", flush=True)
path = f"{FM}/" + t12.loc[t12.study_name == "UKBB_ALT", "sumstats_path"].iloc[0]
subprocess.run(f"awk -F'\\t' 'NR==1 || ($1==\"{C}\" && $2>={A} && $2<={B})' {path} > {OUT}/gwas_window_UKBB_ALT_chr{C}.tsv", shell=True, check=True)
subprocess.run(f"awk -F'\\t' 'NR==1 || $14==\"ENSG00000121858\"' {ROOT}/data/broadaway_eqtl/chr{C}_marginal_summary_results.tsv > {OUT}/eqtl_TNFSF10.tsv", shell=True, check=True)
print("done", flush=True)
