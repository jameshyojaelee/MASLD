#!/usr/bin/env python3
"""Independent recomputation of Figure 6B, 6C and 6D values from upstream outputs.

Inputs (upstream analysis outputs, not the panel's own source tables):
  COLOC r2 per-study master  GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv
  Study tiers                GWAS/finemapping/config/gwas_trait_tier.tsv
  774-unit bulk fit          RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv
  Catalog identity list      catalog_r9d_gene_table.parquet: ensembl_id, symbol,
                             gen_estimate_from_other_ensembl_id ONLY (identity keys, no class columns)

Rules (catalog_rules_r9.tsv):
  genetics supported iff best eligible Tier-1/2 SuSiE PP.H4 > 0.5
  genetics untestable iff no Tier-1/2 COLOC record, or SuSiE untestable
      (insufficient shared posterior) in >= 1 Tier-1/2 study and not positive
  RNA supported iff BH q < 0.05 and |log2FC| > 0.5 (BH over the 24,115 fitted genes)
  RNA untestable iff gene not in the bulk fit (not detected)
  class: both supported -> concordant; gen sup + RNA indet -> genetically_anchored;
         RNA sup + gen indet -> established_state_associated; both indet -> unresolved;
         either untestable -> untested
The class columns of r9d are read only afterwards, for the comparison.
"""
import csv, os, re, json, collections
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig6/recompute")
MASTER = os.path.join(ROOT, "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv")
TIERS = os.path.join(ROOT, "GWAS/finemapping/config/gwas_trait_tier.tsv")
DEG = os.path.join(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
CAT = os.path.join(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.parquet")
UNTESTABLE = "susie_untestable_insufficient_shared_posterior"
ENSG = re.compile(r"^ENSG\d{11}$")

np.random.seed(20261017)  # no randomness used; set for the record

# ---------------- genetics from per-study COLOC r2 ----------------
tiers = pd.read_csv(TIERS, sep="\t")
tier = dict(zip(tiers.study_name, tiers.tier))
t12 = {s for s, t in tier.items() if t in (1, 2)}
assert len(t12) == 35, len(t12)

m = pd.read_csv(MASTER, usecols=["gwas_name", "gene", "ensembl", "PP.H4.abf", "PP.H4.susie", "PP.H3.susie", "method", "ancestry"],
                low_memory=False)
m = m[m.ensembl.astype(str).str.match(ENSG)]
m12 = m[m.gwas_name.isin(t12)].copy()
m12["PP.H4.susie"] = pd.to_numeric(m12["PP.H4.susie"], errors="coerce")
m12["PP.H4.abf"] = pd.to_numeric(m12["PP.H4.abf"], errors="coerce")


def best_by(df, col):
    x = df.dropna(subset=[col]).sort_values([col, "gwas_name"], ascending=[False, True])
    x = x.drop_duplicates("ensembl")
    return x.set_index("ensembl")[[col, "gwas_name", "ancestry"]]


bs = best_by(m12, "PP.H4.susie").rename(columns={"PP.H4.susie": "susie_best", "gwas_name": "susie_study", "ancestry": "susie_anc"})
ba = best_by(m12, "PP.H4.abf").rename(columns={"PP.H4.abf": "abf_best", "gwas_name": "abf_study", "ancestry": "abf_anc"})
nu = m12[m12.method == UNTESTABLE].groupby("ensembl").size().rename("n_untestable_studies")
n_susie_studies = m12.dropna(subset=["PP.H4.susie"]).groupby("ensembl").size().rename("n_susie_studies")
gen = pd.DataFrame(index=pd.Index(sorted(m12.ensembl.unique()), name="ensembl"))
gen = gen.join(bs).join(ba).join(nu).join(n_susie_studies)
gen["n_untestable_studies"] = gen["n_untestable_studies"].fillna(0).astype(int)
gen["n_susie_studies"] = gen["n_susie_studies"].fillna(0).astype(int)


def gstate(r):
    if pd.notna(r.susie_best) and r.susie_best > 0.5:
        return "positive"
    if r.n_untestable_studies > 0:
        return "untestable"
    return "tested_le_0.5" if pd.notna(r.susie_best) else "no_susie_result"


gen["susie_state"] = gen.apply(gstate, axis=1)
gen.to_csv(os.path.join(OUT, "genetics_tier12_per_gene_recomputed.tsv"), sep="\t")

# ---------------- RNA from 774-unit fit ----------------
deg = pd.read_csv(DEG)
deg["ensembl"] = deg.gene.str.replace(r"\.\d+$", "", regex=True)
assert deg.ensembl.is_unique
# re-derive BH from raw p over all fitted genes
p = deg["P.Value"].values
n = len(p)
order = np.argsort(p)
ranked = p[order] * n / (np.arange(n) + 1)
q = np.minimum.accumulate(ranked[::-1])[::-1]
q = np.minimum(q, 1)
bh = np.empty(n); bh[order] = q
deg["bh_recomputed"] = bh
deg["bh_absdiff"] = (deg.bh_recomputed - deg.padj).abs()
deg["canonical"] = (deg.padj < 0.05) & (deg.logFC.abs() > 0.5)
rna = deg.set_index("ensembl")[["logFC", "SE", "P.Value", "padj", "bh_recomputed", "canonical", "symbol"]]

# ---------------- Catalog identity list (keys only) ----------------
cat = pd.read_parquet(CAT, columns=["ensembl_id", "symbol", "gen_estimate_from_other_ensembl_id"])
cat["gen_key"] = np.where(cat.gen_estimate_from_other_ensembl_id.fillna("") != "",
                          cat.gen_estimate_from_other_ensembl_id, cat.ensembl_id)


def calls(cat, gen_key_col):
    c = cat.copy()
    g = gen.reindex(c[gen_key_col].values)
    c["susie_state"] = g.susie_state.values
    c["susie_best"] = g.susie_best.values
    c["abf_best"] = g.abf_best.values
    c["susie_study"] = g.susie_study.values
    c["abf_study"] = g.abf_study.values
    c["n_untestable_studies"] = g.n_untestable_studies.values
    c["n_susie_studies"] = g.n_susie_studies.values
    c["gen_call"] = np.select(
        [c.susie_state.eq("positive"), c.susie_state.isna() | c.susie_state.eq("untestable")],
        ["supported", "untestable"], "indeterminate")
    r = rna.reindex(c.ensembl_id.values)
    c["rna_log2fc"] = r.logFC.values
    c["rna_q"] = r.padj.values
    c["rna_canonical"] = r.canonical.values
    c["rna_call"] = np.where(r.logFC.isna().values, "untestable",
                             np.where(r.canonical.fillna(False).values.astype(bool), "supported", "indeterminate"))
    cls = []
    for gc, rc in zip(c.gen_call, c.rna_call):
        if gc == "untestable" or rc == "untestable":
            cls.append("untested")
        elif gc == "supported" and rc == "supported":
            cls.append("concordant")
        elif gc == "supported":
            cls.append("genetically_anchored")
        elif rc == "supported":
            cls.append("established_state_associated")
        else:
            cls.append("unresolved")
    c["class_recomputed"] = cls
    c["untested_domain"] = np.select([c.rna_call.eq("untestable"), c.gen_call.eq("untestable")],
                                     ["no_rna_call_not_detected", "no_genetic_call"], "")
    return c


out = {}
for label, key in [("crosswalk", "gen_key"), ("same_id", "ensembl_id")]:
    c = calls(cat, key)
    mat = pd.crosstab(c.gen_call, c.rna_call, margins=True)
    mat = mat.reindex(index=["supported", "indeterminate", "untestable", "All"],
                      columns=["supported", "indeterminate", "untestable", "All"]).fillna(0).astype(int)
    mat.to_csv(os.path.join(OUT, f"fig6b_matrix_recomputed_{label}.tsv"), sep="\t")
    cl = c.class_recomputed.value_counts()
    cl.to_csv(os.path.join(OUT, f"fig6d_class_counts_recomputed_{label}.tsv"), sep="\t")
    ud = c[c.class_recomputed == "untested"].untested_domain.value_counts()
    ud.to_csv(os.path.join(OUT, f"fig6d_untested_split_recomputed_{label}.tsv"), sep="\t")
    out[label] = c
    print("==", label)
    print(mat)
    print(cl)
    print(ud)

c = out["crosswalk"]
c.to_csv(os.path.join(OUT, "catalog_calls_recomputed.tsv.gz"), sep="\t", index=False)

# diagnostics
diag = {}
diag["bulk_fit_genes"] = int(len(deg))
diag["bulk_canonical"] = int(deg.canonical.sum())
diag["bulk_bh_max_absdiff_vs_padj"] = float(deg.bh_absdiff.max())
diag["bulk_q_lt_0.05_absfc_le_0.5_in_fit"] = int(((deg.padj < 0.05) & (deg.logFC.abs() <= 0.5)).sum())
diag["catalog_q_lt_0.05_absfc_le_0.5"] = int(((c.rna_q < 0.05) & (c.rna_log2fc.abs() <= 0.5)).sum())
diag["bulk_ids_without_catalog_record"] = int((~deg.ensembl.isin(cat.ensembl_id)).sum())
diag["bulk_canonical_ids_without_catalog_record"] = deg.loc[deg.canonical & ~deg.ensembl.isin(cat.ensembl_id), ["ensembl", "symbol", "logFC", "padj"]].to_dict("records")
diag["tier12_genes"] = int(len(gen))
diag["tier12_state_counts"] = gen.susie_state.value_counts().to_dict()
diag["catalog_gen_state_counts"] = c.susie_state.fillna("no_tier12_record").value_counts().to_dict()
diag["tier12_ids_without_catalog_record"] = int((~gen.index.isin(c.gen_key)).sum())
diag["universe_union_bulk_tier12"] = int(len(set(deg.ensembl) | set(gen.index)))
diag["catalog_rows"] = int(len(c))
diag["catalog_ids_in_neither_input"] = int((~c.ensembl_id.isin(deg.ensembl) & ~c.gen_key.isin(gen.index)).sum())
u7 = c[c.susie_state.eq("untestable")]
diag["untestable7"] = u7[["symbol", "susie_best", "abf_best", "n_untestable_studies", "n_susie_studies", "rna_call"]].to_dict("records")
gi = c[c.gen_call.eq("indeterminate")]
diag["gen_indet_abf_gt_0.5"] = int((gi.abf_best > 0.5).sum())
diag["gen_indet_abf_ge_0.8"] = int((gi.abf_best >= 0.8).sum())
diag["gen_indet_no_susie"] = int(gi.susie_best.isna().sum())
diag["gen_supported_min_susie"] = float(c[c.gen_call.eq("supported")].susie_best.min())
# per-study SuSiE positives (for the 7 untestable genes, PP.H3)
json.dump(diag, open(os.path.join(OUT, "diagnostics.json"), "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in diag.items() if k not in ("untestable7",)}, indent=1, default=str))
print(pd.DataFrame(diag["untestable7"]))

# ---------------- 6C worked entries ----------------
we = ["GLP1R", "HKDC1", "MTARC1", "THRB", "PNPLA3", "TM6SF2", "GCKR", "APOE"]
w = c[c.symbol.isin(we)][["symbol", "ensembl_id", "gen_call", "susie_best", "susie_study", "abf_best", "abf_study",
                          "n_susie_studies", "n_untestable_studies", "rna_call", "rna_log2fc", "rna_q", "class_recomputed"]]
w.to_csv(os.path.join(OUT, "fig6c_worked_entries_recomputed.tsv"), sep="\t", index=False)
print(w.to_string())
# per-study detail for worked entries
det = m12[m12.ensembl.isin(w.ensembl_id)].merge(w[["ensembl_id", "symbol"]], left_on="ensembl", right_on="ensembl_id")
det = det.sort_values(["symbol", "PP.H4.abf"], ascending=[True, False])
det[["symbol", "gwas_name", "PP.H4.abf", "PP.H4.susie", "PP.H3.susie", "method", "ancestry"]].to_csv(
    os.path.join(OUT, "fig6c_worked_entries_per_study.tsv"), sep="\t", index=False)
