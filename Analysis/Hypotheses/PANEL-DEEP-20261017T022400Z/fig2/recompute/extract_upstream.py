"""Extract the upstream rows needed to recompute Figure 2 panels.

Inputs (read-only):
  COLOC r2 master, GWAS registry/tier, combined fine-mapping tables,
  GWAS summary statistics (hg19 reformatted), Broadaway liver eQTL marginals (hg19).
Outputs: reduced tables in ./upstream/ next to this script.
"""
import os, sys, subprocess, json
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
FM = f"{ROOT}/GWAS/finemapping"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "upstream")
os.makedirs(OUT, exist_ok=True)

tier = pd.read_csv(f"{FM}/config/gwas_trait_tier.tsv", sep="\t")
reg = pd.read_csv(f"{FM}/config/gwas_registry.tsv", sep="\t")
t12 = tier[tier.tier.isin([1, 2])].merge(reg[["study_name", "ancestry", "sumstats_path", "ld_panel"]], on="study_name")
t12.to_csv(f"{OUT}/registry_t12.tsv", sep="\t", index=False)
studies = set(t12.study_name)
print("tier1/2 studies", len(studies), flush=True)

# 1. r2 master, all 50 studies kept for reference, but write the 35-study subset
cols = ["gwas_name", "gene", "ensembl", "chr", "PP.H3.abf", "PP.H4.abf", "PP.H3.susie", "PP.H4.susie",
        "n_cs_pairs", "n_cs_pairs_eligible", "method", "top_snp", "top_snp_PP", "ancestry", "ld_panel",
        "eligibility_source"]
m = pd.read_csv(f"{FM}/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv", usecols=cols,
                dtype={"chr": str, "top_snp": str})
print("master rows", len(m), flush=True)
print(m.method.value_counts(dropna=False).to_string(), flush=True)
print(m.eligibility_source.value_counts(dropna=False).to_string(), flush=True)
mt = m[m.gwas_name.isin(studies)].copy()
mt.to_csv(f"{OUT}/r2_master_t12.tsv.gz", sep="\t", index=False)
print("t12 rows", len(mt), flush=True)

# lead variant set: every (study, top_snp) with finite SuSiE PP.H4, or ABF PP.H4 > 0.5
mt["susie_f"] = np.isfinite(pd.to_numeric(mt["PP.H4.susie"], errors="coerce"))
lead = mt[(mt.susie_f) | (pd.to_numeric(mt["PP.H4.abf"], errors="coerce") > 0.5)][["gwas_name", "top_snp"]].dropna().drop_duplicates()
lead[["c", "p"]] = lead.top_snp.str.split(":", expand=True).iloc[:, :2]
lead["p"] = lead.p.astype(int)
lead.to_csv(f"{OUT}/lead_set.tsv", sep="\t", index=False)
print("lead set", len(lead), flush=True)
keyset = set(zip(lead.gwas_name, lead.c.astype(str), lead.p))

# Windows for G/H locus views (hg19)
win = {"15": (60_400_000, 61_300_000), "2": (88_000_000, 88_800_000)}

# 2. fine-mapping records at lead variants (combined_finemapping.csv, 6.5 GB) streamed
for fname in ["combined_finemapping.csv", "combined_finemapping_coloc_targeted.csv"]:
    keep, keepwin = [], []
    for ch in pd.read_csv(f"{FM}/results/{fname}", chunksize=2_000_000, dtype={"chromosome": str}, low_memory=False):
        ch = ch[ch.study.isin(studies)]
        k = [(s, c, p) in keyset for s, c, p in zip(ch.study, ch.chromosome.astype(str), ch.position.astype(int))]
        keep.append(ch[k])
        w = np.zeros(len(ch), bool)
        for c, (a, b) in win.items():
            w |= (ch.chromosome.astype(str) == c).values & (ch.position >= a).values & (ch.position <= b).values
        keepwin.append(ch[w])
    pd.concat(keep).to_csv(f"{OUT}/{fname.replace('.csv','')}_at_leads.tsv.gz", sep="\t", index=False)
    pd.concat(keepwin).to_csv(f"{OUT}/{fname.replace('.csv','')}_GH_windows.tsv.gz", sep="\t", index=False)
    print(fname, "done", flush=True)

# 3. GWAS p at lead variants + G/H windows, per study via awk
rows = []
for s, g in lead.groupby("gwas_name"):
    path = f"{FM}/" + t12.loc[t12.study_name == s, "sumstats_path"].iloc[0]
    posf = f"{OUT}/_pos_{s}.txt"
    g[["c", "p"]].drop_duplicates().to_csv(posf, sep="\t", header=False, index=False)
    outf = f"{OUT}/_gwas_{s}.tsv"
    cmd = (f"awk -F'\\t' 'NR==FNR{{k[$1\":\"$2]=1;next}} FNR==1{{print;next}} ($1\":\"$2) in k' {posf} {path} > {outf}")
    subprocess.run(cmd, shell=True, check=True)
    d = pd.read_csv(outf, sep="\t", dtype={"chromosome": str})
    d["study"] = s
    rows.append(d)
    os.remove(posf); os.remove(outf)
pd.concat(rows).to_csv(f"{OUT}/gwas_p_at_leads.tsv.gz", sep="\t", index=False)
print("gwas p done", flush=True)

for s, c, (a, b) in [("MVP_NAFLD_EUR", "15", win["15"]), ("BBJ_GGT", "15", win["15"]), ("UKBB_ALT", "2", win["2"]),
                     ("UKBB_GGT", "15", win["15"])]:
    path = f"{FM}/" + t12.loc[t12.study_name == s, "sumstats_path"].iloc[0]
    cmd = f"awk -F'\\t' 'NR==1 || ($1==\"{c}\" && $2>={a} && $2<={b})' {path} > {OUT}/gwas_window_{s}_chr{c}.tsv"
    subprocess.run(cmd, shell=True, check=True)
print("gwas windows done", flush=True)

# 4. eQTL marginals for RORA and FABP1 (Broadaway, hg19)
for chrom, ensg, name in [("15", "ENSG00000069667", "RORA"), ("2", "ENSG00000163586", "FABP1")]:
    cmd = (f"awk -F'\\t' 'NR==1 || $14==\"{ensg}\"' {ROOT}/data/broadaway_eqtl/chr{chrom}_marginal_summary_results.tsv"
           f" > {OUT}/eqtl_{name}.tsv")
    subprocess.run(cmd, shell=True, check=True)
print("eqtl done", flush=True)
