#!/usr/bin/env python3
"""Check of the Gene Catalog r9b candidate (adapted from check_catalog_r9.py) against r8 and the plan source files.

Recomputes every new gene-table column from its source file (full table) and reports
the seeded 20-gene sample per column plus known-answer genes. Writes only to OUT.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20260923
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
R4 = ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z"
CAND = R4 / "catalog_r9c_candidate-20261015T191932Z/build_out"
OUT = Path(sys.argv[1])
HYP = ROOT / "Analysis/Hypotheses"
R8 = ROOT / "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8"
BULK = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
COLOC = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
FROZEN = ROOT / "Analysis/Multimodal_Program_Projection/candidates/gen-corrected-r8-20261009T142027Z/genetics_context/frozen_evidence_classes.tsv"
FIG3E = ROOT / "figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv"
FIG3E_COUNTS = ROOT / "figures/main/fig3_bulk_transcriptomics/source_tables/fig3e_deg_counts.tsv"
LIN = ROOT / "figures/candidates/bulk-context-774-20261014T140904Z_22236311/class_difference/results/gene_classes.tsv.gz"
SIX = ROOT / "figures/candidates/review-r3-fixes-20261015T112148Z/F6/panels_20261015T114537Z/R07_11_six_genes.tsv"
SRC = {
    "H53": HYP / "H53-20261015T112655Z/results/catalog_gene_caqtl_annotation.tsv",
    "H55": HYP / "H55-20261015T112558Z/run-20261015T114811Z/outputs/catalog_column_gtex_postmortem_fat.tsv",
    "H63": HYP / "H63-20261015T125909Z/run-20261015T131038Z/h63_zonation_columns.parquet",
    "H66": HYP / "H66-20261015T124530Z/run-20261015T125923Z/h66_catalog_fields_shipped.tsv.gz",
    "H72": HYP / "H72-20261015T134833Z/run-20261015T142309Z/catalog_sex_effect_state.tsv",
    "H87": HYP / "H87-20261015T130324Z/run-20261015T131921Z/h87_coding_cs_gene_column.tsv",
}
KNOWN = ["PNPLA3", "TM6SF2", "MTARC1", "APOE", "SERPINA1", "GNMT", "HSD17B13", "GCKR", "MBOAT7", "TRIB1",
         "PLIN2", "GLUL", "CYP2E1", "CPS1", "ALB", "APOB", "CYP3A4", "KRT19", "DCN", "PECAM1"]

log = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    log.append(s)


def rtsv(p, sep="\t", **k):
    return pd.read_csv(p, sep=sep, float_precision="round_trip", low_memory=False, **k)


# ------------------------------------------------------------------ load candidate and r8
r9 = pd.read_parquet(CAND / "catalog_r9c_candidate_gene_table.parquet")
r9_tsv = rtsv(CAND / "catalog_r9c_candidate_gene_table.tsv.gz", keep_default_na=False, na_values=[""], dtype=str)
gi = pd.read_parquet(R8 / "passport_gene_index.parquet")
ev = pd.read_parquet(R8 / "passport_evidence_long.parquet")
say("r9 shape", r9.shape, "| tsv shape", r9_tsv.shape, "| r8 gene index", gi.shape, "| r8 evidence", ev.shape)
assert r9.ensembl_id.is_unique and gi.ensembl_id.is_unique
ids = pd.Index(gi.ensembl_id)
r9 = r9.set_index("ensembl_id", drop=False)
say("r9 ensembl set == r8 ensembl set:", set(r9.index) == set(ids), "| parquet vs tsv same column order:",
    list(r9.columns) == list(r9_tsv.columns))
r9 = r9.reindex(ids)
import os
if os.environ.get("NEGCTL") == "1":
    # negative control: plant known defects; each must be detected
    r9 = r9.copy()
    pick = list(ids[:3])
    r9.loc[pick, "rna_log2fc"] = r9.loc[pick, "rna_log2fc"] * 1.001
    r9.loc[pick[:2], "rna_canonical_call"] = ~r9.loc[pick[:2], "rna_canonical_call"].astype(bool)
    r9.loc[pick[:1], "primary_evidence_class"] = "concordant"
    r9.loc[pick[:1], "zon_mrna_call_yakubovsky2026"] = "portal_planted"
    r9.loc[pick[:2], "gen_coloc_best_abf_ancestry"] = "XXX"
    r9.loc[pick[:1], "lineage_measured_in_sc"] = pd.NA
    # r9b planted errors on the changed columns (each must be detected)
    nt_ids = list(r9.index[r9.caqtl_gene_link == "not_tested_in_source"][:2])
    r9.loc[nt_ids, "caqtl_gene_link"] = "not_linked"                       # revert two to the r9 lumping
    ind_id = list(r9.index[r9.caqtl_gene_link == "indeterminate"][:1])
    r9.loc[ind_id, "caqtl_not_tested_reason"] = "no_tier12_eur_pair"        # reason on a non-not-tested gene
    tie_ids = list(r9.index[r9.sc_best_lineage == "tied"][:3])
    r9.loc[tie_ids, "sc_best_lineage"] = r9.loc[tie_ids, "sc_lineages_at_max"].str.split(";").str[0]  # r9 behaviour
    one_ids = list(r9.index[r9.sc_best_lineage.notna() & (r9.sc_best_lineage != "tied")][:1])
    r9.loc[one_ids, "sc_best_lineage"] = "tied"
    r9.loc[tie_ids[:1], "sc_lineages_at_max"] = r9.loc[tie_ids[:1], "sc_lineages_at_max"].str.rsplit(";", n=1).str[0]  # drop one
    hcol = [c for c in r9.columns if c.startswith("h51_gse213621_F3F4")]
    rep_ids = list(r9.index[r9["h51_gse213621_F3F4_vs_F0F1_state"] == "replicated"][:2])
    r9.loc[rep_ids, "h51_gse213621_F3F4_vs_F0F1_state"] = "underpowered"
    r9.loc[rep_ids[:1], "h51_gse213621_F3F4_vs_F0F1_log2fc"] = r9.loc[rep_ids[:1], "h51_gse213621_F3F4_vs_F0F1_log2fc"] * -1
    PLANTED = {"caqtl_gene_link": 2, "caqtl_not_tested_reason": 1, "sc_best_lineage": 4, "sc_lineages_at_max": 1,
               "h51_gse213621_F3F4_vs_F0F1_state": 2, "h51_gse213621_F3F4_vs_F0F1_log2fc": 1, "primary_evidence_class": 1,
               "blocked_on_fig3e_row": 1}

# ------------------------------------------------------------------ (1) r8 unchanged
results = []  # rows: check, column, n_compared, n_diff, note
r8_cols = list(gi.columns)
g8 = gi.set_index("ensembl_id", drop=False)
for c in r8_cols:
    a, b = g8[c], r9[c]
    diff = ~((a.astype(str) == b.astype(str)) | (a.isna() & b.isna()))
    results.append(("1_r8_gene_index", c, len(a), int(diff.sum()), ""))
gen = ev[ev.evidence_domain == "genetics"].set_index("ensembl_id").reindex(ids)
rna = ev[ev.evidence_domain == "transcriptomics"].set_index("ensembl_id").reindex(ids)
for dom, t in [("genetics", gen), ("rna", rna)]:
    for s in ["call_state", "testability_state"]:
        a, b = t[s].astype(str), r9[f"r8_{dom}_{s}"].astype(str)
        results.append(("1_r8_evidence_calls", f"r8_{dom}_{s}", len(a), int((a != b).sum()), ""))
pec_diff = int((g8.primary_evidence_class.astype(str) != r9.primary_evidence_class.astype(str)).sum())
say("(1) primary_evidence_class differences r8 vs r9:", pec_diff)
say("(1) class counts r9:", r9.primary_evidence_class.value_counts().to_dict())
# Fig 6B 3x3
x = pd.crosstab(gen.call_state, rna.call_state)
say("(1) r8 3x3 genetics (rows) x RNA (cols):\n" + x.to_string())

# ------------------------------------------------------------------ recompute new columns
E = pd.DataFrame(index=ids)  # expected values

# bulk 774-unit fit
coef = rtsv(BULK / "coefficient_table.tsv.gz")
coef["eid"] = coef.gene.str.replace(r"\.\d+$", "", regex=True)
assert coef.eid.is_unique
cf = coef.set_index("eid").reindex(ids)
E["rna_log2fc"], E["rna_se"], E["rna_t"] = cf.logFC, cf.SE, cf.t
E["rna_p"], E["rna_bh_q"], E["rna_ave_log2cpm"] = cf["P.Value"], cf.padj, cf.AveExpr
inf = cf.logFC.notna()
E["rna_bh_q_lt_0.05"] = (cf.padj < 0.05).where(inf)
E["rna_canonical_call"] = ((cf.padj < 0.05) & (cf.logFC.abs() > 0.5)).where(inf)
design = rtsv(BULK / "model_design.tsv")
n_u, n_d = len(design), int((design.group_binaryDisease == 1).sum())
say("design units / disease / control:", n_u, n_d, n_u - n_d)
E["rna_n_units"] = pd.Series(n_u, index=ids).where(inf)
E["rna_n_control_units"] = pd.Series(n_u - n_d, index=ids).where(inf)
E["rna_n_disease_units"] = pd.Series(n_d, index=ids).where(inf)
E["rna_unit"] = pd.Series("genetic-identity cluster (not verified donor)", index=ids).where(inf, "")
# independent BH re-derivation over the 24,115-gene family
from statsmodels.stats.multitest import multipletests
bh = multipletests(coef["P.Value"].values, method="fdr_bh")[1]
say("BH over", len(coef), "genes from P.Value vs file padj: max rel diff",
    float(np.nanmax(np.abs(bh - coef.padj.values) / np.maximum(coef.padj.values, 1e-300))))
say("canonical calls (q<0.05,|lfc|>0.5) in full fit:", int(((coef.padj < 0.05) & (coef.logFC.abs() > 0.5)).sum()))

# COLOC r2 Tier-1/2
ct = rtsv(COLOC / "gene_level_coloc_tier12.csv", sep=",")
assert ct.ensembl.is_unique
c2 = ct.set_index("ensembl").reindex(ids)
anc = pd.read_csv(COLOC / "susie_coloc_all_gwas.csv", usecols=["gwas_name", "ancestry"], dtype=str).drop_duplicates()
say("COLOC studies with >1 ancestry label:", int(anc.gwas_name.duplicated().sum()))
amap = anc.set_index("gwas_name").ancestry
E["gen_coloc_best_susie_pp4"] = c2.coloc_best_susie_pp4
E["gen_coloc_best_susie_study"] = c2.susie_driving_gwas
E["gen_coloc_best_susie_ancestry"] = c2.susie_driving_gwas.map(amap)
E["gen_coloc_best_abf_pp4"] = c2.coloc_best_abf_pp4
E["gen_coloc_best_abf_study"] = c2.abf_driving_gwas
E["gen_coloc_best_abf_ancestry"] = c2.abf_driving_gwas.map(amap)
E["gen_susie_state_t12"] = c2.susie_state_t12
# suffix cross-check of ancestry label (independent of the map)
suf = amap.index.to_series().str.extract(r"_(EUR|EAS|AFR|AMR|SAS|CSA|MID)$")[0]
mism = amap[suf.notna() & (suf != amap)]
say("studies whose name suffix disagrees with the ancestry label:", len(mism), mism.to_dict())
say("studies without suffix -> ancestry:", amap[suf.isna()].to_dict())
gest = pd.to_numeric(gen.estimate, errors="coerce")
su = np.isclose(gest, c2.coloc_best_susie_pp4, atol=1e-9, rtol=0)
ab = np.isclose(gest, c2.coloc_best_abf_pp4, atol=1e-9, rtol=0)
E["gen_r8_displayed_pp4_method"] = np.select([gest.isna(), su, ab], ["none", "SuSiE", "ABF"], "value_from_other_ensembl_id")
six = rtsv(SIX)
E["gen_estimate_from_other_ensembl_id"] = pd.Series(ids.map(dict(zip(six.catalog_id, six.ensembl_genetic))), index=ids).fillna("")
# the six r8 estimates should equal the other ID's Tier-1/2 value
ct_i = ct.set_index("ensembl")
for _, r in six.iterrows():
    e8 = float(gest.loc[r.catalog_id])
    o = ct_i.loc[r.ensembl_genetic] if r.ensembl_genetic in ct_i.index else None
    say("  six-gene", r.symbol, r.catalog_id, "r8 est", e8, "| other-ID susie/abf",
        None if o is None else (o.coloc_best_susie_pp4, o.coloc_best_abf_pp4))
E["gen_untestable_reason_r9"] = gen.testability_state.map(
    {"underpowered_source": "no_tier12_coloc_record", "insufficient_shared_posterior": "susie_untestable_in_ge1_tier12_study"}).fillna("")
# genetics call rule re-derived from COLOC
sup_rule = set(ids[(c2.coloc_best_susie_pp4 > 0.5).values])
sup_r8 = set(ids[(gen.call_state == "supported").values])
say("genetics supported (r8) == SuSiE PP.H4>0.5 (COLOC r2):", sup_rule == sup_r8, len(sup_r8), len(sup_rule))
# RNA call rule re-derived
rs = set(ids[(E.rna_canonical_call == True).values])
say("RNA supported (r8) == canonical rule on fit:", rs == set(ids[(rna.call_state == "supported").values]), len(rs))

cls = g8.primary_evidence_class.reindex(ids)
E["class_display_r9"] = cls.map({"concordant": "both supported", "genetically_anchored": "genetically anchored",
                                  "established_state_associated": "established state-associated",
                                  "unresolved": "unresolved", "untested": "untested"})
unt = cls == "untested"
E["untested_missing_domain"] = np.where(unt & (rna.call_state == "untestable"), "no_rna_call_not_detected",
                                        np.where(unt & (gen.call_state == "untestable"), "no_genetic_call", ""))
E["untested_one_domain_supported"] = np.where(unt & (rna.call_state == "supported"), "rna_supported",
                                              np.where(unt & (gen.call_state == "supported"), "genetics_supported", ""))
E["assay_bulk_rna_r9"] = rna.call_state.astype(str)
E["assay_proteomics_r9"] = "not_assessed_at_gene_level"
E["assay_chromatin_r9"] = "not_joined_at_gene_level"
E["assay_spatial_r9"] = "not_joined_at_gene_level"

# H53
h53 = rtsv(SRC["H53"])
assert h53.ensembl.is_unique
say("H53 rows", len(h53), "| H53 genes == r8 genetics supported:", set(h53.ensembl) == sup_r8)
h = h53.set_index("ensembl").reindex(ids)
for c in ["caqtl_region_hg38", "caqtl_pp_h4", "caqtl_lead_ld_r2", "caqtl_three_way_eqtl", "caqtl_abc_predicted_link",
          "caqtl_source_flag"]:
    E[c] = h[c]
E["caqtl_region_state"] = h.caqtl_region_state.fillna("not_evaluated_genetics_not_supported")
E["caqtl_n_gwas_signals"] = h.n_gwas_signals
car = h.caqtl_region_state == "colocalized accessible region"
tw = h.caqtl_three_way_eqtl.astype(str).str.lower()
# r9b link states, written independently of the build: keyed on the H53 state text
rs_ = h.caqtl_region_state.astype(object).fillna("")
link = pd.Series("not_evaluated_genetics_not_supported", index=ids, dtype=object)
link[car & (tw == "true")] = "M1_three_way_supported_linked"
link[car & (tw == "false")] = "M2_r2_pair_only_not_called_linked"
link[rs_ == "different variant (unsupported)"] = "not_linked"
link[rs_.str.startswith("not tested in this source")] = "not_tested_in_source"
link[rs_ == "indeterminate"] = "indeterminate"
link[rs_ == "colocalized, lead LD not confirmed"] = "colocalized_lead_ld_unconfirmed"
say("H53 states not mapped by the checker:", sorted(set(rs_[(rs_ != "") & (link == "not_evaluated_genetics_not_supported")])))
E["caqtl_gene_link"] = link
E["caqtl_not_tested_reason"] = np.where(rs_ == "not tested in this source (no tier-1/2 EUR pair)", "no_tier12_eur_pair",
                                        np.where(rs_ == "not tested in this source", "not_tested_other", ""))
say("(r9b) caqtl_gene_link expected counts:", link.value_counts().to_dict())
say("(r9b) caqtl_not_tested_reason expected counts:", pd.Series(E["caqtl_not_tested_reason"]).value_counts().to_dict())
# independent M1 check from the H53 three-way table
tw_tab = rtsv(HYP / "H53-20261015T112655Z/results/three_way.tsv")
say("H53 three_way.tsv columns", list(tw_tab.columns), "rows", len(tw_tab))
say("H53 region_state counts:", h53.caqtl_region_state.value_counts(dropna=False).to_dict())
say("H53 three_way among colocalized region:", h53[h53.caqtl_region_state == "colocalized accessible region"]
    .caqtl_three_way_eqtl.value_counts(dropna=False).to_dict())

# H55
h55 = rtsv(SRC["H55"])
assert h55.ensembl_id.is_unique
h5 = h55.set_index("ensembl_id").reindex(ids)
for c in ["gtex_postmortem_fat_slope", "gtex_postmortem_fat_ci_low", "gtex_postmortem_fat_ci_high", "gtex_postmortem_fat_q"]:
    E[c] = h5[c]
E["gtex_postmortem_fat_state"] = h5.gtex_postmortem_fat_state.fillna("not_a_canonical_bulk_call")
say("H55 genes == canonical calls of the fit:", set(h55.ensembl_id) == set(coef.eid[(coef.padj < 0.05) & (coef.logFC.abs() > 0.5)]))

# H63
h63 = pd.read_parquet(SRC["H63"]).set_index("ensembl_id").reindex(ids)
say("H63 symbol mismatches vs r8:", int((h63.symbol.astype(str) != g8.symbol.reindex(ids).astype(str)).sum()))
for c in [c for c in h63.columns if c != "symbol"]:
    E[c] = h63[c]

# H66
h66 = rtsv(SRC["H66"]).set_index("ensembl_id").reindex(ids)
say("H66 symbol mismatches vs r8:", int((h66.symbol.astype(str) != g8.symbol.reindex(ids).astype(str)).sum()))
for c in [c for c in h66.columns if c != "symbol"]:
    E[c] = h66[c]
# r9b: ties at the maximum donor detection fraction, recomputed from the per-lineage fractions
h66a = rtsv(HYP / "H66-20261015T124530Z/run-20261015T125923Z/h66_fields_all_including_unshipped.tsv.gz").set_index("ensembl_id").reindex(ids)
fcols = sorted(c for c in h66a.columns if c.startswith("sc_frac_"))
Fv = h66a[fcols].to_numpy(dtype=float)
row = ~np.isnan(Fv).all(axis=1)
mxv = np.where(row, np.nanmax(np.where(np.isnan(Fv), -1, Fv), axis=1), np.nan)
best_l, at_l = [], []
for i in range(len(ids)):
    if not row[i]:
        best_l.append(np.nan); at_l.append(np.nan); continue
    tl = [fcols[j][8:] for j in range(len(fcols)) if Fv[i, j] == mxv[i]]
    if mxv[i] == 0:  # r9c: detected in no lineage
        at_l.append(np.nan); best_l.append("not_detected"); continue
    at_l.append(";".join(tl)); best_l.append("tied" if len(tl) > 1 else tl[0])
E["sc_best_lineage"] = pd.Series(best_l, index=ids, dtype=object)
E["sc_lineages_at_max"] = pd.Series(at_l, index=ids, dtype=object)
near = int(sum(1 for i in range(len(ids)) if row[i] and
               ((np.abs(Fv[i] - mxv[i]) < 1e-5) & (Fv[i] != mxv[i])).any()))
say("(r9b) sc ties:", int((E.sc_best_lineage == "tied").sum()), "of", int(row.sum()),
    "| at fraction 0:", int(((E.sc_best_lineage == "tied") & (mxv == 0)).sum()),
    "| near-ties within 1e-5 but not equal:", near)
say("(r9b) r9 label among tied lineages (alphabetically first):",
    int(sum(1 for i in ids[(E.sc_best_lineage == "tied").values]
            if str(h66.at[i, "sc_best_lineage"]) == E.at[i, "sc_lineages_at_max"].split(";")[0])))
# r9b: H51 F3F4 candidate columns
h51 = rtsv(HYP / "H51-20261015T112617Z/stage2-20261015T164829Z/per_gene_states.tsv.gz")
h51 = h51[h51.contrast == "F3F4"].copy()
h51["eid"] = h51.gene.str.replace(r"\.\d+$", "", regex=True)
assert h51.eid.is_unique
say("(r9b) H51 F3F4 rows", len(h51), "| in Catalog", int(h51.eid.isin(set(ids)).sum()),
    "| catalog labels", h51.catalog.value_counts().to_dict())
hh = h51.set_index("eid").reindex(ids)
pre = "h51_gse213621_F3F4_vs_F0F1"
E[f"{pre}_log2fc"], E[f"{pre}_ci95_lo"], E[f"{pre}_ci95_hi"], E[f"{pre}_bh_q"] = hh.G0_log2FC, hh.G0_lo, hh.G0_hi, hh.G0_BH_q
E[f"{pre}_state"] = hh.state.fillna("not_a_discovery_hit")
E[f"{pre}_state_accession_adjusted"] = hh.acc_state.fillna("not_a_discovery_hit")
h51f2 = rtsv(HYP / "H51-20261015T112617Z/stage2-20261015T164829Z/per_gene_states.tsv.gz")
say("(r9b) H51 F2 Catalog label:", h51f2[h51f2.contrast == "F2"].catalog.unique().tolist())

# H87
h87 = rtsv(SRC["H87"]).set_index("ensembl_id").reindex(ids)
say("H87 primary_evidence_class mismatches vs r8:",
    int((h87.primary_evidence_class.astype(str) != cls.astype(str)).sum()))
for c in [c for c in h87.columns if c.startswith("coding_cs_")]:
    E[c] = h87[c]

# H72
h72 = rtsv(SRC["H72"])
assert h72.ensembl_id.is_unique
say("H72 rows", len(h72), "| == anchored+concordant:",
    set(h72.ensembl_id) == set(ids[cls.isin(["genetically_anchored", "concordant"]).values]))
E["enzyme_sex_effect_state"] = h72.set_index("ensembl_id").state.reindex(ids).fillna("not_in_tested_family")

# Fig 3E
st = rtsv(FIG3E, usecols=["gene_id_base", "logFC", "FDR", "axis", "contrast"])
cnt = rtsv(FIG3E_COUNTS)
for k in ["F1", "F2", "F3", "F4"]:
    s = st[(st.axis == "fibrosis") & (st.contrast == f"{k}_vs_F0")]
    assert s.gene_id_base.is_unique and len(s) == 24115
    hi, lo = int(((s.FDR < 0.05) & (s.logFC > 0)).sum()), int(((s.FDR < 0.05) & (s.logFC < 0)).sum())
    ph = int(cnt[(cnt.axis == "fibrosis") & (cnt.display == k) & (cnt.direction == "Higher")].n.iloc[0])
    pl = int(cnt[(cnt.axis == "fibrosis") & (cnt.display == k) & (cnt.direction == "Lower")].n.iloc[0])
    say(f"Fig 3E {k} vs F0 higher/lower (full 24,115) {hi}/{lo} | installed panel table {ph}/{pl}")
    ss = s.set_index("gene_id_base").reindex(ids)
    E[f"fig3e_{k}_vs_F0_log2fc"] = ss.logFC
    E[f"fig3e_{k}_vs_F0_bh_q"] = ss.FDR
    E[f"fig3e_{k}_vs_F0_call"] = np.select([ss.logFC.isna(), (ss.FDR < 0.05) & (ss.logFC > 0), (ss.FDR < 0.05) & (ss.logFC < 0)],
                                           ["not_in_fit", f"higher_in_{k}", f"lower_in_{k}"], "not_called")
    # does the fig3e fit share the same gene universe as the 774 fit?
    say(f"  {k}: fig3e gene set == 774-fit gene set:", set(s.gene_id_base) == set(coef.eid))

# lineage
lin = rtsv(LIN)
assert lin.g.is_unique
li = lin.set_index("g").reindex(ids)
E["lineage_hep_specificity_log2cpm"] = li.hep_specificity
E["lineage_measured_in_sc"] = li.in_sc

# ------------------------------------------------------------------ (2) compare
new_cols = [c for c in r9.columns if c not in r8_cols]
missing_expect = [c for c in new_cols if c not in E.columns]
say("(2) new columns:", len(new_cols), "| without an independent recomputation:", missing_expect)


def norm(s):
    """Return (numeric series or None, string series) with NA and '' unified."""
    nn = s.dropna()
    if len(nn) and set(map(str, nn.unique())) <= {"True", "False", "true", "false"}:
        return None, s.astype(object).where(s.notna(), "").map(lambda v: str(v).lower())
    num = pd.to_numeric(s, errors="coerce")
    nonblank = s.notna() & (s.astype(str).str.strip() != "") & (s.astype(str) != "<NA>") & (s.astype(str) != "nan")
    if nonblank.any() and num[nonblank].notna().all() and s.dtype != bool and str(s.dtype) != "boolean":
        return num, None
    t = s.astype(object).where(nonblank, "")
    t = t.map(lambda v: str(v).strip().lower() if isinstance(v, (bool, np.bool_)) or str(v) in ("True", "False") else str(v))
    # 1.0/0.0 vs 1/0 style
    t = t.str.replace(r"^(-?\d+)\.0$", r"\1", regex=True)
    return None, t


def compare(c):
    a, b = E[c], r9[c]
    na, sa = norm(a)
    nb, sb = norm(b)
    if na is not None and nb is not None:
        both_na = na.isna() & nb.isna()
        close = np.isclose(na.fillna(0), nb.fillna(0), rtol=1e-9, atol=1e-12) & (na.isna() == nb.isna())
        ok = both_na | close
    else:
        if sa is None:
            sa = a.astype(object).where(a.notna(), "").map(lambda v: str(v)).str.replace(r"^(-?\d+)\.0$", r"\1", regex=True)
        if sb is None:
            sb = b.astype(object).where(b.notna(), "").map(lambda v: str(v)).str.replace(r"^(-?\d+)\.0$", r"\1", regex=True)
        ok = sa == sb
    return ok


rng = np.random.default_rng(SEED)
sample_rows = []
sym = g8.symbol.reindex(ids)
known_ids = [i for i in ids if sym.loc[i] in KNOWN]
say("known-answer genes found:", sorted(sym.loc[known_ids].tolist()))
for c in new_cols:
    if c not in E.columns:
        results.append(("2_column_reproduction", c, 0, 0, "r8 copy: checked in section 1 against evidence_long"))
        continue
    ok = compare(c)
    src_has = E[c].notna() & (E[c].astype(str) != "")
    pool = ids[src_has.values] if src_has.sum() >= 20 else ids
    samp = list(rng.choice(np.asarray(pool), size=min(20, len(pool)), replace=False))
    sel = list(dict.fromkeys(samp + known_ids))
    nd_s = int((~ok.loc[samp]).sum())
    nd_k = int((~ok.loc[known_ids]).sum())
    results.append(("2_column_reproduction", c, len(ids), int((~ok).sum()),
                    f"sample20_diff={nd_s}; known_diff={nd_k}"))
    for i in sel:
        sample_rows.append((c, i, sym.loc[i], "sample" if i in samp else "known_answer",
                            E.at[i, c], r9.at[i, c], bool(ok.loc[i])))
    if (~ok).any():
        bad = ids[(~ok).values][:5]
        for i in bad:
            say(f"  MISMATCH {c} {i} {sym.loc[i]}: expected={E.at[i, c]!r} r9={r9.at[i, c]!r}")
pd.DataFrame(sample_rows, columns=["column", "ensembl_id", "symbol", "draw", "expected_from_source", "r9_value", "match"]
             ).to_csv(OUT / "sample_reproduction.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ (3) H56 absent
h56_hits = []
files = [CAND / "catalog_r9c_candidate_gene_table.parquet"] + sorted(CAND.glob("sidecar_*.tsv*")) + [CAND / "column_dictionary.tsv"]
for f in files:
    cols = pd.read_parquet(f).columns if f.suffix == ".parquet" else pd.read_csv(f, sep="\t", nrows=0).columns
    for col in cols:
        if re.search(r"composition|h56|deconv|cell_mix|celltype_adj", col, re.I):
            h56_hits.append((f.name, col))
dic = pd.read_csv(CAND / "column_dictionary.tsv", sep="\t", dtype=str, keep_default_na=False)
for _, r in dic.iterrows():
    if re.search(r"bulk_composition_sensitivity|\bH56\b", " ".join(r.values)):
        h56_hits.append(("column_dictionary.tsv", f"{r.table}:{r.column}"))
say("(3) H56 / composition hits in column names or dictionary rows:", h56_hits)
results.append(("3_H56_absent", "all files", len(files), len(h56_hits), str(h56_hits)))

# ------------------------------------------------------------------ (4) dictionary vs files
dic_tables = dic.table.unique()
for t in dic_tables:
    if t in ("table",):
        continue
    f = CAND / ("catalog_r9c_candidate_gene_table.parquet" if t == "catalog_r9c_candidate_gene_table" else t)
    if not f.exists():
        results.append(("4_dictionary_columns", t, 0, -1, "file missing"))
        continue
    if f.suffix == ".parquet":
        fc = list(pd.read_parquet(f).columns)
    else:
        fc = list(pd.read_csv(f, sep="\t", nrows=0).columns)
    dc = list(dic[dic.table == t].column)
    extra_d, extra_f = sorted(set(dc) - set(fc)), sorted(set(fc) - set(dc))
    results.append(("4_dictionary_columns", t, len(fc), len(extra_d) + len(extra_f) + int(dc != fc and not extra_d and not extra_f),
                    f"dict_not_in_file={extra_d}; file_not_in_dict={extra_f}; order_same={dc == fc}"))
# allowed-values counts for gene-table columns
gt = dic[dic.table == "catalog_r9c_candidate_gene_table"].set_index("column")
r9t = r9_tsv.set_index("ensembl_id", drop=False)
for c in r9.columns:
    av = gt.at[c, "allowed_values"] if c in gt.index else ""
    s = r9[c]
    m = re.match(r"numeric \[([^,]+), ([^\]]+)\]; n non-missing (\d+)", av)
    if m:
        v = pd.to_numeric(s, errors="coerce")
        lo, hi, n = float(m.group(1)), float(m.group(2)), int(m.group(3))
        good = (np.isclose(v.min(), lo, rtol=1e-5) and np.isclose(v.max(), hi, rtol=1e-5) and int(v.notna().sum()) == n)
        results.append(("4_dictionary_values", c, n, 0 if good else 1,
                        f"dict [{lo},{hi}] n={n}; file [{v.min()},{v.max()}] n={int(v.notna().sum())}"))
        continue
    listed = re.findall(r"\((\d+)\)(?:; |$)", av)
    vc = s.astype(object).where(s.notna(), "<missing>").map(str).str.replace(r"^(-?\d+)\.0$", r"\1", regex=True).value_counts()
    vc_raw = s.astype(object).where(s.notna(), "<missing>").map(str).value_counts()
    found = 0
    for vcx in (vc, vc_raw):
        f_ = sum(1 for val, n in vcx.items() if f"{val} ({n})" in av)
        found = max(found, f_)
    # values listed in dictionary but with a count that matches no value in the file
    results.append(("4_dictionary_values", c, len(listed), max(0, len(listed) - found),
                    f"listed={len(listed)} matched={found} distinct_in_file={len(vc)}"))

# sidecar row/col counts vs build_summary
bs = json.load(open(CAND / "build_summary.json"))
for f, (nr, nc) in bs["sidecars"].items():
    d = pd.read_csv(CAND / f, sep="\t", low_memory=False)
    results.append(("4_sidecar_shape", f, nr, int(d.shape != (nr, nc)), f"summary {nr}x{nc}; file {d.shape}"))

# ------------------------------------------------------------------ (5) universes
fr = rtsv(FROZEN, dtype={"joint_testable": str})
u_sym = int((fr.joint_testable.astype(str).str.lower() == "true").sum())
deg = rtsv(BULK / "deg_results.csv", sep=",")
deg_ids = set(deg.gene.str.replace(r"\.\d+$", "", regex=True))
u_sid = len(deg_ids & set(ct.ensembl))
u_cat = int((cls != "untested").sum())
for name, obs, exp in [("symbol-grouped joint universe", u_sym, 14720), ("stable-ID join bulk x COLOC T1/2", u_sid, 14704),
                       ("Catalog r8 joint", u_cat, 14682), ("Catalog r8 rows", len(gi), 28233), ("r9 rows", len(r9), 28233),
                       ("774 fit genes", len(coef), 24115), ("Catalog genes with fit row", int(inf.sum()), 24001)]:
    results.append(("5_universe", name, exp, int(obs != exp), f"observed {obs}"))
uc = pd.read_csv(CAND / "universe_counts.tsv", sep="\t")
say("(5) candidate universe_counts.tsv:\n" + uc[["universe", "n", "expected"]].to_string())

# blocked_on vs standing rules
bo = pd.read_csv(CAND / "blocked_on.tsv", sep="\t")
say("blocked_on.tsv:\n" + bo.to_string())
dep = dic[dic.table == "catalog_r9c_candidate_gene_table"][["column", "release_state", "blocked_on"]]
say("release_state counts (gene table):", dep.release_state.value_counts().to_dict())
say("blocked_on counts (gene table):", dep.blocked_on.value_counts().to_dict())

# ------------------------------------------------------------------ (6) r9b-specific checks
CLS = {"untested": 13551, "unresolved": 13127, "established_state_associated": 1129, "genetically_anchored": 400, "concordant": 26}
obs_cls = r9.primary_evidence_class.value_counts().to_dict()
results.append(("6_class_counts", "primary_evidence_class", len(r9), int(sum(abs(obs_cls.get(k, 0) - v) for k, v in CLS.items())
                + sum(v for k, v in obs_cls.items() if k not in CLS)), f"observed {obs_cls}"))
results.append(("6_primary_class_unchanged_vs_r8", "primary_evidence_class", len(ids), pec_diff, ""))
# every column the r9b fixes did not target must equal r9 exactly
R9 = R4 / "catalog_r9_candidate-20261015T173124Z/catalog_r9_candidate_gene_table.parquet"
old = pd.read_parquet(R9).set_index("ensembl_id", drop=False).reindex(ids)
CHANGED = {"caqtl_gene_link", "caqtl_not_tested_reason", "sc_best_lineage", "sc_lineages_at_max"} | {c for c in r9.columns if c.startswith("h51_")}
added = [c for c in r9.columns if c not in old.columns]
dropped = [c for c in old.columns if c not in r9.columns]
results.append(("6_columns_added_vs_r9", ";".join(added), len(added), int(set(added) - CHANGED != set()), f"dropped={dropped}"))
results.append(("6_columns_dropped_vs_r9", ";".join(dropped), len(dropped), len(dropped), ""))
for c in [c for c in old.columns if c in r9.columns and c not in CHANGED]:
    a, b = old[c], r9[c]
    diff = ~((a.astype(str) == b.astype(str)) | (a.isna() & b.isna()))
    if diff.any():
        results.append(("6_unchanged_vs_r9", c, len(ids), int(diff.sum()), ""))
results.append(("6_unchanged_vs_r9", "all untargeted columns (rows only written on difference)", len(old.columns) - len(CHANGED & set(old.columns)), 0, ""))
# changed columns: what moved from r9
for c in ["caqtl_gene_link", "sc_best_lineage"]:
    mv = pd.crosstab(old[c].fillna("<missing>"), r9[c].fillna("<missing>"))
    mv = mv.stack()
    mv = mv[mv > 0]
    say(f"(6) r9 -> r9b transitions {c}:", {f"{a} -> {b}": int(n) for (a, b), n in mv.items() if a != b})
# blocked_on and dictionary
if os.environ.get("NEGCTL") == "1":
    bo.loc[bo.columns_or_table.str.startswith("fig3e_"), "blocked_on"] = "H51;H82"
fr_ = bo[bo.columns_or_table.str.startswith("fig3e_")]
results.append(("6_blocked_on_fig3e_row", "blocked_on.tsv", len(fr_),
                int(len(fr_) != 1 or not str(fr_.blocked_on.iloc[0]).startswith("none"))
                + int(not ("RESULT_H51.md" in str(fr_.reason.iloc[0]) and "RESULT_H82.md" in str(fr_.reason.iloc[0])))
                + int(not ("untestable_S5" in str(fr_.reason.iloc[0]) and "F3F4" in str(fr_.reason.iloc[0]))),
                str(fr_.to_dict("records"))))
f3d = dic[(dic.table == "catalog_r9c_candidate_gene_table") & dic.column.str.startswith("fig3e_")]
results.append(("6_dictionary_fig3e_blocked_on", "fig3e_* rows", len(f3d), int(f3d.blocked_on.str.contains(r"H51|H82").sum()), ""))
results.append(("6_dictionary_fig3e_cites_results", "fig3e_* rows", len(f3d),
                int((~(f3d.result.str.contains("RESULT_H51.md") & f3d.result.str.contains("RESULT_H82.md"))).sum()), ""))
results.append(("6_r9c_fig3e_caveat_has_no_H82_numbers", "fig3e_* rows", len(f3d), int(f3d.caveat.str.contains("104-188").sum()), ""))
results.append(("6_r9c_not_detected_has_no_lineage", "sc_best_lineage", int((r9.sc_best_lineage == "not_detected").sum()), int(((r9.sc_best_lineage == "not_detected") & r9.sc_lineages_at_max.notna()).sum()) + int(((r9.sc_best_lineage == "not_detected") & (r9.sc_best_donor_frac.astype(float) != 0)).sum()), ""))
scd = dic[(dic.table == "catalog_r9c_candidate_gene_table") & dic.column.isin(["sc_best_lineage", "sc_lineages_at_max"])]
results.append(("6_dictionary_sc_note", "sc_best_lineage; sc_lineages_at_max", len(scd),
                int((~scd.caveat.str.contains("detection breadth") | ~scd.caveat.str.contains("not the expressing cell type")).sum()) + (2 - len(scd)), ""))
h51d = dic[(dic.table == "catalog_r9c_candidate_gene_table") & dic.column.str.startswith("h51_")]
results.append(("6_dictionary_h51_candidate", "h51_* rows", len(h51d),
                int((~h51d.release_state.str.contains("candidate")).sum()) + int((~h51d.blocked_on.str.contains("user approval")).sum()), ""))
results.append(("6_no_F2_column", "h51_*F2*", 0, int(any("_F2_" in c for c in r9.columns if c.startswith("h51_"))), ""))
if os.environ.get("NEGCTL") == "1":
    say("PLANTED (column: rows):", PLANTED)

res = pd.DataFrame(results, columns=["check", "column_or_item", "n", "n_diff", "note"])
res.to_csv(OUT / "check_results.tsv", sep="\t", index=False)
say("\n=== summary: checks with n_diff != 0 ===")
say(res[res.n_diff != 0].to_string())
(OUT / "check_log.txt").write_text("\n".join(log) + "\n")
