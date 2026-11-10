# DIG1 job 69 (written before any MAGMA set result for these sets was read).
# Do the genes that rise in CCL20+ (inflamed) hepatocytes carry liver-trait GWAS signal beyond being disease genes?
# Sets (ENSG via HGNC symbol; ENSG kept as is):
#   STATE_UP: excess > 0.25 in all four single-cell cohorts (jobs 53/56, CCL20 anchor, depth mode), minus FCAMR, TRIM31,
#             UBD, CXCL10, IL32, CCL20, HSD17B13, PNPLA3 (the §41 definition).
#   STATE_DN: excess < -0.25 in all four.
#   STATE_UP_HEP: STATE_UP genes with hepatocyte UMI share >= 0.85 in all three cohorts of job 68a (descriptive).
# MAGMA v1.10 competitive gene-set analysis on the CG0 gene results (run-20261023T022444Z; SNP-wise mean, 0 kb),
# universe = the 12,510 CGB4 genes, conditioned (condition-hide) on P_cont, P_dis (disease |t| strength), AveExpr and
# LOEUF from CGB4 covar.txt; direction=pos.
# Primary: STATE_UP, one-sided, Holm over the 6 liver traits (UKBB ALT, AST; MVP ALT, AST; Pazoki PDFF; Ghouse
# cirrhosis). Called enriched if Holm < 0.05 in a trait and P < 0.05 in its other-cohort enzyme replicate (ALT/AST);
# BMI (Yengo) is the negative-control trait and must have P > 0.05. STATE_DN and STATE_UP_HEP are descriptive.
# Sensitivity: STATE_UP without genes within 500 kb of a lead locus is not run (no clumped loci file); instead the set is
# rerun without genes with gene P < 2.5e-6 in that trait (descriptive).
# Usage: python -I 69_state_magma_sets.py <new_out_dir>
import sys, os, gzip, csv
import pandas as pd
out = sys.argv[1]; assert not os.path.exists(out); os.makedirs(out)
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; RUN = f"{ROOT}/Analysis/Hypotheses/DIG1-20261024T113951Z"
g53 = pd.read_csv(f"{RUN}/out_chemokine_hepatocyte_state/gene_lfc_matched.tsv.gz", sep="\t")
g53 = g53[g53.donors == "masld"].drop(columns="donors")
g56 = pd.read_csv(f"{RUN}/out_atlas_state_replication/gene_lfc_matched.tsv.gz", sep="\t")
SC = pd.concat([g53, g56])
SC = SC[(SC["mode"] == "depth") & (SC.anchor == "CCL20") & SC.mean_excess.notna()]
x = SC.groupby("gene").mean_excess.agg(["count", "min", "max"]); x = x[x["count"] == 4]
DROP = {"FCAMR", "TRIM31", "UBD", "CXCL10", "IL32", "CCL20", "HSD17B13", "PNPLA3"}
up = sorted(set(x.index[x["min"] > 0.25]) - DROP); dn = sorted(set(x.index[x["max"] < -0.25]) - DROP)
hs = pd.read_csv(f"{RUN}/out_hep_share/gene_hep_share.tsv.gz", sep="\t")
hmin = hs.groupby("gene").median_hep_share.agg(["min", "count"]); hepg = set(hmin.index[(hmin["count"] == 3) & (hmin["min"] >= 0.85)])
H = pd.read_csv(f"{ROOT}/data/external/hgnc_complete_set_20260830/hgnc_complete_set.txt", sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
s2e = dict(zip(H.symbol, H.ensembl_gene_id))
def ens(gs): return sorted({g if g.startswith("ENSG") else s2e.get(g) for g in gs} - {None})
sets = {"STATE_UP": ens(up), "STATE_DN": ens(dn), "STATE_UP_HEP": ens([g for g in up if g in hepg])}
cov = pd.read_csv(f"{ROOT}/Analysis/Hypotheses/CGB4-20261022T122538Z/out-20261023T032051Z/covar.txt", sep="\t")
U = set(cov.GENE)
with open(f"{out}/state.sets", "w") as f:
    for k, v in sets.items():
        f.write(k + " " + " ".join(g for g in v if g in U) + "\n")
with open(f"{out}/set_sizes.tsv", "w") as f:
    f.write("set\tsymbols\tensg\tin_universe\n")
    for k, v, s in [("STATE_UP", sets["STATE_UP"], up), ("STATE_DN", sets["STATE_DN"], dn), ("STATE_UP_HEP", sets["STATE_UP_HEP"], [g for g in up if g in hepg])]:
        f.write(f"{k}\t{len(s)}\t{len(v)}\t{sum(g in U for g in v)}\n")
print(open(f"{out}/set_sizes.tsv").read())
