#!/usr/bin/env python3
"""A2 Q2: lineage chromatin of the COLOC r2 shared-signal posterior vs expression lineage.

Unit: gene (masses averaged over the gene's primary pairs). Family Q2 (BH over 2):
(i) one-sided Fisher 2x2, (ii) five-lineage agreement permutation.
"""
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from a2common import (LIN, NONHEP, RES, SEED, W, add_bins, bh, boot_p, boot_ratio, ci,
                      matched_expectation, peak_classes)

MIN_POS = 5  # Amendment A5

cv = pd.read_parquet(W / "coloc_vars.parquet")
pos = pd.read_parquet(W / "pos_annot.parquet")
anc = pd.read_parquet(W / "anc_annot.parquet")
edges = pd.read_csv(W / "ld_tertile_edges.tsv", sep="\t")
edges = {r.ancestry: (r.ld_q1, r.ld_q2) for r in edges.itertuples()}
expr = pd.read_csv(RES / "q2_expression_lineage.tsv", sep="\t")

cv["pid"] = cv.gwas_name + "|" + cv.ensembl + "|" + cv.signal_pair_index.astype(str)
cv["allele1"] = cv.allele1.astype(str).str.upper()
cv["allele2"] = cv.allele2.astype(str).str.upper()
tot = cv.groupby("pid").snp_pp_h4.sum()
cv = cv[cv.mapping_status == "mapped"].dropna(subset=["chromosome", "position"]).copy()
cv = cv.astype({"chromosome": np.int64, "position": np.int64})
cv = cv.merge(pos, on=["chromosome", "position"], how="inner")
mapped = cv.groupby("pid").snp_pp_h4.sum() / tot
testable = mapped[mapped >= 0.95].index
print("pairs", len(tot), "testable", len(testable), "min mapped", float(mapped.min()), flush=True)
cv = cv[cv.pid.isin(testable)].copy()
cv["w"] = cv.snp_pp_h4 / cv.pid.map(cv.groupby("pid").snp_pp_h4.sum())
cv = cv.merge(anc[["ancestry", "chromosome", "position", "allele1", "allele2", "maf", "ld"]],
              on=["ancestry", "chromosome", "position", "allele1", "allele2"], how="left")
classes = peak_classes(cv, "in_") + peak_classes(cv, "inb_")
# descriptive: which non-hepatocyte lineages' peaks hold the non-hepatocyte-only posterior
for l in NONHEP:
    cv[f"Nin_{l}"] = (cv.Nonly.astype(bool) & cv[f"in_{l}"].astype(bool)).astype(np.int8)
classes += [f"Nin_{l}" for l in NONHEP]
cv = add_bins(cv, edges).reset_index(drop=True)
cv.to_parquet(W / "q2_coloc_annot.parquet", index=False)

ex = matched_expectation(cv, cv[["pid", "w", "bin_full", "bin_tm", "bin_t"] + classes],
                         "pid", classes)
ob = ex[classes].mul(ex.w, axis=0).groupby(ex.pid).sum().add_prefix("obs_")
ee = ex[[f"exp_{c}" for c in classes]].mul(ex.w, axis=0).groupby(ex.pid).sum()
pair = ob.join(ee).reset_index()
meta = cv.drop_duplicates("pid")[["pid", "gwas_name", "gene", "ensembl", "trait_class",
                                  "ancestry", "pp_h4"]]
pair = meta.merge(pair, on="pid")
pair.to_csv(RES / "q2_pair_level.tsv.gz", sep="\t", index=False)

# ---- gene level ------------------------------------------------------------------------
vals = [c for c in pair.columns if c.startswith(("obs_", "exp_"))]
gene = pair.groupby("gene")[vals].mean()
gene["n_pairs"] = pair.groupby("gene").size()
gene["studies"] = pair.groupby("gene").gwas_name.apply(lambda s: ";".join(sorted(set(s))))
gene["trait_classes"] = pair.groupby("gene").trait_class.apply(lambda s: ";".join(sorted(set(s))))
gene["best_pp_h4"] = pair.groupby("gene").pp_h4.max()
gene = gene.reset_index()


def chrom_call(g, floor, sfx=""):
    H, N, HN = g[f"obs_Honly{sfx}"], g[f"obs_Nonly{sfx}"], g[f"obs_HN{sfx}"]
    hep = (H >= floor) & (H >= 2 * N)
    non = (N >= floor) & (N >= 2 * H)
    mixed = (H >= floor) & (N >= floor) & ~hep & ~non
    shared = (HN >= floor) & ~hep & ~non & ~mixed
    return np.select([hep, non, mixed, shared],
                     ["hepatocyte", "non_hepatocyte", "mixed", "shared_H_and_N"], "not_localized")


for fl in (0.10, 0.25, 0.50):
    gene[f"chrom_call_{fl:.2f}"] = chrom_call(gene, fl)
gene["chrom_call_both_0.25"] = chrom_call(gene, 0.25, "b")
spec_o = gene[[f"obs_spec_{l}" for l in NONHEP]].to_numpy()
gene["nonhep_lineage_top"] = np.array(NONHEP)[spec_o.argmax(axis=1)]
gene["nonhep_lineage_top_mass"] = spec_o.max(axis=1)
gene.loc[gene.nonhep_lineage_top_mass < 0.01, "nonhep_lineage_top"] = "none"
nin = gene[[f"obs_Nin_{l}" for l in NONHEP]].to_numpy()
gene["Nonly_mass_by_lineage"] = [
    ";".join(f"{NONHEP[j]}:{row[j]:.2f}" for j in np.argsort(-row) if row[j] >= 0.01)
    for row in nin]
exc = np.stack([gene[f"obs_spec_{l}"] - gene[f"exp_spec_{l}"] for l in LIN], axis=1)
gene["chrom_argmax5_excess"] = np.array(LIN)[exc.argmax(axis=1)]
gene["chrom_max_excess"] = exc.max(axis=1)
gene = gene.merge(expr, on="gene", how="left")
gene.to_csv(RES / "q2_gene_level.tsv", sep="\t", index=False)

# ---- tests -----------------------------------------------------------------------------
tests, xt = [], []
for fl in ("0.25", "0.10", "0.50", "both_0.25"):
    cc = f"chrom_call_{fl}"
    sub = gene[gene[cc].isin(["hepatocyte", "non_hepatocyte"]) &
               gene.expr_hep_class.isin(["hepatocyte_dominant", "non_hepatocyte_dominant"])]
    t = pd.crosstab(sub[cc], sub.expr_hep_class).reindex(
        index=["hepatocyte", "non_hepatocyte"],
        columns=["hepatocyte_dominant", "non_hepatocyte_dominant"], fill_value=0)
    orr, p = fisher_exact(t.to_numpy(), alternative="greater")
    tests.append({"test": f"fisher_2x2_floor_{fl}", "n_genes": len(sub), "a_hepC_hepE": t.iat[0, 0],
                  "b_hepC_nonE": t.iat[0, 1], "c_nonC_hepE": t.iat[1, 0], "d_nonC_nonE": t.iat[1, 1],
                  "odds_ratio": orr, "p": p})
    full = pd.crosstab(gene[cc], gene.expr_hep_class)
    full.insert(0, "floor", fl)
    xt.append(full)
pd.concat(xt).to_csv(RES / "q2_chromatin_by_expression_crosstab.tsv", sep="\t")

rng = np.random.default_rng(SEED)
sub = gene[(gene.chrom_max_excess > 0) & gene.expr_call5.isin(LIN)]
obs_agree = float((sub.chrom_argmax5_excess == sub.expr_call5).mean())
lab = sub.expr_call5.to_numpy()
perm = np.array([(sub.chrom_argmax5_excess.to_numpy() == rng.permutation(lab)).mean()
                 for _ in range(10_000)])
p_perm = (np.sum(perm >= obs_agree) + 1) / (len(perm) + 1)
tests.append({"test": "five_lineage_agreement_perm", "n_genes": len(sub), "odds_ratio": np.nan,
              "agreement": obs_agree, "null_mean": perm.mean(), "p": p_perm})
pd.crosstab(sub.chrom_argmax5_excess, sub.expr_call5).to_csv(
    RES / "q2_five_lineage_crosstab.tsv", sep="\t")
tests = pd.DataFrame(tests)
fam = tests.test.isin(["fisher_2x2_floor_0.25", "five_lineage_agreement_perm"])
tests["family"] = np.where(fam, "Q2 (BH over 2)", "sensitivity")
tests["q_bh"] = np.nan
tests.loc[fam, "q_bh"] = bh(tests.loc[fam, "p"])
tests.to_csv(RES / "q2_tests.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(tests.round(4).to_string(), flush=True)

# ---- aggregate O/E by trait class (gene bootstrap) --------------------------------------
REPORT = [f"spec_{l}" for l in LIN] + [f"shared_{l}" for l in LIN] + \
    [f"any_{l}" for l in LIN] + ["Honly", "Nonly", "HN", "anypeak"]
rows = []
for tc, sdf in pair.groupby("trait_class"):
    units = sdf.groupby("gene")[[f"obs_{c}" for c in REPORT] + [f"exp_{c}" for c in REPORT]].mean()
    r = boot_ratio(units, REPORT)
    for i, c in enumerate(REPORT):
        lo, hi = ci(r["boot"][:, i])
        npos = int((units[f"obs_{c}"] > 0).sum())
        with np.errstate(divide="ignore"):
            lr = np.log(r["boot"][:, i])
        rows.append({"trait_class": tc, "peak_class": c, "n_genes": len(units), "n_pairs": len(sdf),
                     "share_obs": r["obs"][i], "share_exp": r["exp"][i], "ratio": r["ratio"][i],
                     "ratio_lo": lo, "ratio_hi": hi, "n_genes_obs_pos": npos,
                     "p_boot_descriptive": boot_p(lr) if npos >= MIN_POS else np.nan})
agg = pd.DataFrame(rows)
agg.to_csv(RES / "q2_oe_by_class.tsv", sep="\t", index=False)
print(agg.round(4).to_string(), flush=True)

# ---- discordant gene lists --------------------------------------------------------------
top = (cv.assign(wm=cv.w)
       .sort_values("wm", ascending=False))


def top_variant(gname, cls):
    s = top[(top.gene == gname) & (top[cls] == 1)]
    if s.empty:
        return pd.Series({"top_variant_hg19": "", "top_variant_hg38": "", "top_variant_w": np.nan,
                          "top_variant_study": ""})
    r = s.iloc[0]
    return pd.Series({"top_variant_hg19": f"{r.chromosome}:{r.position}:{r.allele1}:{r.allele2}",
                      "top_variant_hg38": f"chr{r.chromosome}:{r.pos38}",
                      "top_variant_w": r.w, "top_variant_study": r.gwas_name})


keep_cols = ["gene", "studies", "trait_classes", "n_pairs", "best_pp_h4", "obs_Honly", "obs_Nonly",
             "obs_HN", "exp_Honly", "exp_Nonly", "nonhep_lineage_top", "nonhep_lineage_top_mass",
             "Nonly_mass_by_lineage",
             "chrom_call_0.25", "chrom_call_both_0.25", "expr_call6", "stab6", "hep_share6"] + \
    [f"log2cpm_{l}" for l in LIN + ["endothelial"]]
a = gene[(gene["chrom_call_0.25"] == "non_hepatocyte") &
         (gene.expr_hep_class == "hepatocyte_dominant")][keep_cols].copy()
a = pd.concat([a.reset_index(drop=True),
               a.gene.apply(lambda g: top_variant(g, "Nonly")).reset_index(drop=True)], axis=1)
a = a.sort_values("obs_Nonly", ascending=False)
a.to_csv(RES / "q2_nonhep_chromatin_hep_expression.tsv", sep="\t", index=False)
b = gene[(gene["chrom_call_0.25"] == "hepatocyte") &
         (gene.expr_hep_class == "non_hepatocyte_dominant")][keep_cols].copy()
b = pd.concat([b.reset_index(drop=True),
               b.gene.apply(lambda g: top_variant(g, "Honly")).reset_index(drop=True)], axis=1)
b = b.sort_values("obs_Honly", ascending=False)
b.to_csv(RES / "q2_hep_chromatin_nonhep_expression.tsv", sep="\t", index=False)
print("non-hep chromatin / hep expression", len(a), a.gene.tolist()[:30])
print("hep chromatin / non-hep expression", len(b), b.gene.tolist()[:30])
