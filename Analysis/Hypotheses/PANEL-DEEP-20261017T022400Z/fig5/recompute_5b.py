"""Independent recomputation of Figure 5B values from upstream analysis outputs.

Inputs (upstream, not the panel's own source table):
  protein: Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv
  RNA:     RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/
           primary_source_controls/deg_results.csv (774-unit pooled fit)
  universe: Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt (July 10,044 genes)
"""
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = BASE / "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5/recompute"
OUT.mkdir(parents=True, exist_ok=True)

prot = pd.read_csv(BASE / "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv", sep="\t")
rna = pd.read_csv(BASE / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
uni = [x.strip() for x in open(BASE / "Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt") if x.strip()]
uni = pd.Series(uni).drop_duplicates()
print("universe", len(uni))
print("protein rows", len(prot), "unique genes", prot.gene.nunique(), "non-NA logFC", prot.logFC.notna().sum())
prot_ok = prot[prot.logFC.notna() & prot.padj.notna()].copy()
print("protein rows with logFC+padj", len(prot_ok), "unique", prot_ok.gene.nunique())
dupp = prot_ok.gene.duplicated(keep=False)
print("duplicated protein gene symbols", dupp.sum())

# RNA: symbol duplicates
rna["symbol"] = rna["symbol"].astype(str)
dup = rna.symbol.duplicated(keep=False)
print("RNA rows", len(rna), "dup symbols rows", dup.sum())

in_uni_prot = prot_ok[prot_ok.gene.isin(set(uni))]
print("universe genes measured in protein", in_uni_prot.gene.nunique())

# Join: RNA by symbol; unique symbols only first
rna_u = rna[~dup].copy()
# Resolve duplicated RNA symbols through the multi-evidence atlas symbol -> Ensembl map
atlas = pd.read_csv(BASE / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv", usecols=["human_symbol", "ensembl_id"])
amap = dict(zip(atlas.human_symbol, atlas.ensembl_id))
rd = rna[dup].copy()
rd["base"] = rd.gene.str.split(".").str[0]
keep = rd[[amap.get(s) == b for s, b in zip(rd.symbol, rd.base)]]
print("duplicated-symbol RNA rows resolved through atlas:", len(keep))
rna_u = pd.concat([rna_u, keep.drop(columns="base")])
m = in_uni_prot.merge(rna_u[["gene", "symbol", "logFC", "padj"]].rename(columns={"gene": "ensg", "logFC": "rna_logFC", "padj": "rna_padj"}),
                      left_on="gene", right_on="symbol", how="left")
print("joined rows", len(m), "with RNA", m.rna_logFC.notna().sum())
missing = m[m.rna_logFC.isna()].gene.tolist()
print("measured-in-protein universe genes without unique RNA symbol:", len(missing), missing[:40])
# Of these, how many have duplicated RNA symbol rows
print("  of which have duplicated RNA symbol:", [g for g in missing if g in set(rna[dup].symbol)])
m = m[m.rna_logFC.notna()].copy()
n = len(m)
rho = stats.spearmanr(m.rna_logFC, m.logFC).correlation
agree = (np.sign(m.rna_logFC) == np.sign(m.logFC)).mean()
slope, intercept, r, p, se = stats.linregress(m.rna_logFC, m.logFC)
sig = m.padj < 0.05
conc = sig & (np.sign(m.rna_logFC) == np.sign(m.logFC))
disc = sig & ~(np.sign(m.rna_logFC) == np.sign(m.logFC))
rna_sig = m.rna_padj < 0.05
res = {
    "n_measured": n, "spearman_rho": rho, "direction_agreement": agree,
    "ols_slope": slope, "ols_intercept": intercept, "ols_r2": r ** 2,
    "n_protein_sig": int(sig.sum()), "n_concordant": int(conc.sum()), "n_discordant": int(disc.sum()),
    "pct_concordant": 100 * conc.sum() / sig.sum(), "pct_discordant": 100 * disc.sum() / sig.sum(),
    "n_not_sig": int((~sig).sum()),
    "sig_rnasig_n": int((sig & rna_sig).sum()), "sig_rnasig_conc": int((conc & rna_sig).sum()),
    "sig_rnans_n": int((sig & ~rna_sig).sum()), "sig_rnans_conc": int((conc & ~rna_sig).sum()),
    "x_min": m.rna_logFC.min(), "x_max": m.rna_logFC.max(), "y_min": m.logFC.min(), "y_max": m.logFC.max(),
}
# binomial for RNA-nonsig agreement
res["binom_p_rnans"] = stats.binomtest(res["sig_rnans_conc"], res["sig_rnans_n"], 0.5).pvalue
# alternative: sign agreement among protein-significant restricted to |RNA log2FC|>0.5 & padj<.05
strong = sig & rna_sig & (m.rna_logFC.abs() > 0.5)
res["sig_rnastrong_n"] = int(strong.sum()); res["sig_rnastrong_conc"] = int((conc & strong).sum())
# Pearson for reference
res["pearson_r"] = stats.pearsonr(m.rna_logFC, m.logFC)[0]
# all measured proteins rho
allm = prot_ok.merge(rna_u[["symbol", "logFC"]].rename(columns={"logFC": "rna_logFC"}), left_on="gene", right_on="symbol")
res["all_measured_n"] = len(allm)
res["all_measured_rho"] = stats.spearmanr(allm.rna_logFC, allm.logFC).correlation
for k, v in res.items():
    print(f"{k}\t{v}")
pd.Series(res).to_csv(OUT / "5b_recomputed_summary.tsv", sep="\t", header=False)
m["class"] = np.where(~sig, "not_significant", np.where(conc, "concordant", "discordant"))
m[["gene", "ensg", "rna_logFC", "rna_padj", "logFC", "padj", "class"]].to_csv(OUT / "5b_recomputed_genes.tsv", sep="\t", index=False)
