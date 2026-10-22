#!/usr/bin/env python3
"""F6 re-derivation of Gene Catalog findings R07-02/04/05/06/07/11, R15-02/03/05/09.

Reads saved tables only (Catalog r8 reviewed copy, adopted 774-unit bulk fit,
COLOC r2, PXD051911 protein DE). Writes per-finding tables to OUT.
No randomness.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
CAT = ROOT / "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8"
COMP = ROOT / "figures/candidates/catalog-r8-completion-20261014T141645Z_22236527"
BULK = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
COLOC = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
PROT = ROOT / "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv"
GEN = ROOT / "Analysis/Multimodal_Program_Projection/candidates/gen-corrected-r8-20261009T142027Z/genetics_context/frozen_evidence_classes.tsv"
OUT = Path(sys.argv[1])
OUT.mkdir(parents=False, exist_ok=False)

S = {}  # summary numbers
e = pd.read_parquet(CAT / "passport_evidence_long.parquet")
g = pd.read_parquet(CAT / "passport_gene_index.parquet")
gen = e[e.evidence_domain == "genetics"].set_index("ensembl_id")
rna = e[e.evidence_domain == "transcriptomics"].set_index("ensembl_id")
S["n_genes"] = len(g)
S["n_evidence_rows"] = len(e)

# ---------- R07-02: negative rule never evaluated; PP.H3 evidence exists ----------
S["R07_02_negative_rule_na"] = int(e.negative_call_rule_id.isna().sum())
S["R07_02_negative_passed_na"] = int(e.negative_call_passed.isna().sum())
S["R07_02_call_states"] = e.groupby(["evidence_domain", "call_state"]).size().to_dict()
S["R07_02_call_states"] = {f"{k[0]}|{k[1]}": int(v) for k, v in S["R07_02_call_states"].items()}
sp = pd.read_csv(COLOC / "susie_coloc_signal_pairs.tsv.gz", sep="\t")
el = sp[sp.eligible == True].copy()
el["ens"] = el.ensembl.str.split(".").str[0]
gm = el.groupby("ens").agg(max_h4=("PP.H4", "max"), max_h3=("PP.H3", "max"), symbol=("gene", "first"))
S["R07_02_signal_pair_rows"] = len(sp)
S["R07_02_eligible_pairs"] = len(el)
S["R07_02_genes_eligible"] = len(gm)
h3 = gm[(gm.max_h3 > 0.8) & (gm.max_h4 <= 0.5)].copy()
h3["catalog_genetics_call"] = h3.index.map(gen.call_state).fillna("absent_from_catalog")
h3["strict_h3_0.9_h4_0.1"] = (h3.max_h3 > 0.9) & (h3.max_h4 < 0.1)
S["R07_02_h3_gt0.8_h4_le0.5"] = len(h3)
S["R07_02_h3_by_catalog_call"] = h3.catalog_genetics_call.value_counts().to_dict()
S["R07_02_strict"] = int(h3["strict_h3_0.9_h4_0.1"].sum())
h3.sort_values("max_h3", ascending=False).to_csv(OUT / "R07_02_pph3_distinct_signal_genes.tsv", sep="\t")

# ---------- bulk coefficient table (adopted 774-unit fit) ----------
c = pd.read_csv(BULK / "coefficient_table.tsv.gz", sep="\t")
c["ens"] = c.gene.str.split(".").str[0]
c["canonical"] = (c.padj < 0.05) & (c.logFC.abs() > 0.5)
S["bulk_genes"] = len(c)
S["bulk_canonical"] = int(c.canonical.sum())
cb = c.set_index("ens")

# ---------- R07-04: RNA indeterminate with resolved direction ----------
r = rna.join(cb[["logFC", "SE", "t", "P.Value", "padj", "canonical"]], how="left")
ind = r[r.call_state == "indeterminate"]
S["R07_04_rna_indeterminate"] = len(ind)
S["R07_04_rna_indet_q_lt_0.05"] = int((ind.padj < 0.05).sum())
S["R07_04_rna_indet_q_ge_0.05"] = int((ind.padj >= 0.05).sum())
S["R07_04_rna_indet_no_bulk_row"] = int(ind.padj.isna().sum())
gi = g.set_index("ensembl_id")
for cls in ["unresolved", "genetically_anchored"]:
    ids = gi.index[gi.primary_evidence_class == cls]
    sub = r.loc[r.index.intersection(ids)]
    S[f"R07_04_{cls}_n"] = len(ids)
    S[f"R07_04_{cls}_q_lt_0.05"] = int((sub.padj < 0.05).sum())
for sym in ["THRB", "MTARC1", "HKDC1", "GLP1R", "GNMT"]:
    row = r[r.symbol == sym]
    S[f"worked_{sym}_rna"] = row[["call_state", "estimate", "q_value", "logFC", "SE", "padj"]].to_dict("records")

# ---------- R07-05: bulk-tested IDs absent from the Catalog ----------
cat_ids = set(g.ensembl_id)
miss = c[~c.ens.isin(cat_ids)].copy()
S["R07_05_bulk_ids_absent"] = len(miss)
S["R07_05_absent_padj_lt_0.05"] = int((miss.padj < 0.05).sum())
S["R07_05_absent_canonical"] = miss[miss.canonical][["ens", "symbol", "logFC", "padj"]].to_dict("records")
symcount = c.symbol.value_counts()
miss["symbol_shared_with_other_id"] = miss.symbol.map(symcount).fillna(1) > 1
S["R07_05_absent_symbol_shared"] = int(miss.symbol_shared_with_other_id.sum())
miss.sort_values("padj").to_csv(OUT / "R07_05_bulk_ids_absent_from_catalog.tsv", sep="\t", index=False)
S["R07_05_catalog_rna_supported"] = int((rna.call_state == "supported").sum())
pm = pd.read_parquet(CAT / "passport_program_membership.parquet")
idcol = [x for x in pm.columns if "ensembl" in x.lower()]
if idcol:
    pid = pm[idcol[0]].astype(str).str.split(".").str[0]
    S["R07_05_program_member_ids_absent"] = int(pid[~pid.isin(cat_ids)].nunique())
    S["R07_05_program_memberships_absent"] = int((~pid.isin(cat_ids)).sum())
coloc_gene = pd.read_csv(COLOC / "gene_level_coloc_tier12.csv")
cg_ids = set(coloc_gene.ensembl.astype(str).str.split(".").str[0])
S["R07_05_coloc_t12_ids"] = len(cg_ids)
S["R07_05_coloc_t12_ids_absent"] = len(cg_ids - cat_ids)

# ---------- R07-06: concordant has no direction comparison ----------
S["R07_06_genetics_direction"] = gen.direction.value_counts().to_dict()
conc = gi.index[gi.primary_evidence_class == "concordant"]
cc = r.loc[conc, ["symbol", "estimate", "q_value"]].copy()
cc["rna_direction"] = np.where(cc.estimate > 0, "up", "down")
S["R07_06_concordant_n"] = len(conc)
S["R07_06_concordant_rna_up"] = int((cc.estimate > 0).sum())
cc.to_csv(OUT / "R07_06_concordant_genes.tsv", sep="\t")

# ---------- R07-07: untested decomposition, wording ----------
S["R07_07_testability"] = e.groupby(["evidence_domain", "testability_state"]).size()
S["R07_07_testability"] = {f"{k[0]}|{k[1]}": int(v) for k, v in S["R07_07_testability"].items()}
S["R07_07_not_applicable_rows"] = int((e.call_state == "not_applicable").sum())
S["R07_07_source_dependent_true"] = int(e.source_dependent.sum())
S["R07_07_provenance"] = e.provenance_state.value_counts().to_dict()
cross = pd.crosstab(gen.call_state.rename("genetics"), rna.call_state.rename("rna"))
cross.to_csv(OUT / "R07_07_call_matrix.tsv", sep="\t")
untested = gi.index[gi.primary_evidence_class == "untested"]
ut = pd.DataFrame({"genetics_call": gen.loc[untested, "call_state"],
                   "genetics_testability": gen.loc[untested, "testability_state"],
                   "rna_call": rna.loc[untested, "call_state"],
                   "rna_testability": rna.loc[untested, "testability_state"]})
ut_tab = ut.value_counts().rename("n_genes").reset_index()
ut_tab.to_csv(OUT / "R07_07_untested_breakdown.tsv", sep="\t", index=False)
S["R07_07_untested_n"] = len(untested)
wording = e.groupby(["evidence_domain", "call_state", "testability_state"]).allowed_wording.agg(lambda s: " || ".join(sorted(set(s))))
wording.reset_index().to_csv(OUT / "R07_07_allowed_wording_by_state.tsv", sep="\t", index=False)
# genes with no COLOC record at all among 'underpowered_source'
up = gen[gen.testability_state == "underpowered_source"].index
S["R07_07_underpowered_in_any_coloc_t12"] = int(pd.Index(up).isin(list(cg_ids)).sum())
coloc_all = pd.read_csv(COLOC / "susie_coloc_all_gwas.csv", usecols=["gwas_name", "ensembl", "PP.H4.abf", "PP.H4.susie", "ancestry", "method"])
coloc_all["ens"] = coloc_all.ensembl.astype(str).str.split(".").str[0]
all_ids = set(coloc_all.ens)
S["R07_07_underpowered_in_any_coloc_all"] = int(pd.Index(up).isin(list(all_ids)).sum())
S["R07_07_role_hypothesis"] = g.role_hypothesis.value_counts().to_dict()

# ---------- R15-02 + R15-03: per-gene sidecar with method and denominators ----------
cg = coloc_gene.copy()
cg["ens"] = cg.ensembl.astype(str).str.split(".").str[0]
cg = cg.set_index("ens")
gg = gen[["symbol", "estimate", "call_state", "testability_state"]].join(
    cg[["coloc_best_susie_pp4", "coloc_best_abf_pp4", "susie_driving_gwas", "abf_driving_gwas", "abf_driving_trait", "driving_trait"]], how="left")
est = gg.estimate.astype(float)
gg["displayed_equals_susie"] = np.isclose(est, gg.coloc_best_susie_pp4, atol=1e-9)
gg["displayed_equals_abf"] = np.isclose(est, gg.coloc_best_abf_pp4, atol=1e-9)
gg["displayed_method"] = np.select([gg.displayed_equals_susie, gg.displayed_equals_abf, est.isna()],
                                   ["SuSiE", "ABF", "none"], "other")
t = pd.crosstab(gg.call_state, gg.displayed_method)
t.to_csv(OUT / "R15_02_displayed_pp4_method_by_call.tsv", sep="\t")
S["R15_02_method_by_call"] = {f"{a}|{b}": int(t.loc[a, b]) for a in t.index for b in t.columns}
ind_g = gg[gg.call_state == "indeterminate"]
S["R15_02_indet_est_gt_0.5"] = int((ind_g.estimate.astype(float) > 0.5).sum())
S["R15_02_indet_est_ge_0.8"] = int((ind_g.estimate.astype(float) >= 0.8).sum())
S["R15_02_indet_est_gt0.5_with_finite_susie"] = int(((ind_g.estimate.astype(float) > 0.5) & ind_g.coloc_best_susie_pp4.notna()).sum())
S["R15_02_supported_pp4_lt_0.8"] = int((gg[gg.call_state == "supported"].estimate.astype(float) < 0.8).sum())
for sym in ["GNMT", "DGAT2", "CFLAR", "EZR", "THRB", "HKDC1", "MTARC1", "GLP1R"]:
    row = gg[gg.symbol == sym]
    S[f"R15_02_{sym}"] = row[["estimate", "call_state", "coloc_best_susie_pp4", "coloc_best_abf_pp4", "susie_driving_gwas", "abf_driving_gwas"]].to_dict("records")
# best SuSiE study ancestry for the supported gene
anc = coloc_all.dropna(subset=["PP.H4.susie"]).sort_values("PP.H4.susie", ascending=False).drop_duplicates("ens").set_index("ens")
side = pd.DataFrame(index=g.ensembl_id)
side["symbol"] = gi.symbol
side["primary_evidence_class"] = gi.primary_evidence_class
side["genetics_call"] = gen.call_state
side["genetics_testability"] = gen.testability_state
side["coloc_best_susie_pp4"] = cg.coloc_best_susie_pp4
side["coloc_best_susie_study"] = cg.susie_driving_gwas
side["coloc_best_susie_ancestry"] = anc.ancestry
side["coloc_best_abf_pp4"] = cg.coloc_best_abf_pp4
side["coloc_best_abf_study"] = cg.abf_driving_gwas
side["coloc_call_rule"] = "supported iff best eligible SuSiE PP.H4 > 0.5 (ABF not used for the call)"
side["coloc_displayed_estimate_method"] = gg.displayed_method
side["rna_call"] = rna.call_state
side["rna_log2fc"] = cb.logFC
side["rna_se"] = cb.SE
side["rna_t"] = cb.t
side["rna_p"] = cb["P.Value"]
side["rna_bh_q"] = cb.padj
side["rna_direction_resolved_bh05"] = cb.padj < 0.05
side["rna_n_units"] = np.where(side.rna_log2fc.notna(), 774, np.nan)
side["rna_n_control_units"] = np.where(side.rna_log2fc.notna(), 138, np.nan)
side["rna_n_disease_units"] = np.where(side.rna_log2fc.notna(), 636, np.nan)
side["rna_unit"] = "genetic-identity cluster (not verified donor)"
side.reset_index().rename(columns={"index": "ensembl_id"}).to_csv(OUT / "R15_03_gene_effect_sidecar.tsv.gz", sep="\t", index=False)
S["R15_03_null_counts"] = e[["standard_error", "ci_lower", "ci_upper", "p_value", "n_biological_units", "n_technical_units"]].notna().sum().to_dict()
S["R15_03_biological_unit_values"] = e.biological_unit.value_counts().to_dict()
S["R15_03_sidecar_rna_se_filled"] = int(side.rna_se.notna().sum())

# ---------- R15-05: molecular-object assay_applicability ----------
rows = [json.loads(l) for l in open(COMP / "molecular_objects/molecular_object_sidecar.jsonl")]
mo = pd.DataFrame(rows)
aa = mo.assay_applicability.apply(json.loads).apply(pd.Series)
mo = pd.concat([mo[["ensembl_id_base", "symbol", "rna_call", "genetics_call", "primary_evidence_class"]], aa], axis=1)
S["R15_05_records"] = len(mo)
S["R15_05_rna_call_by_bulk_rna"] = {f"{a}|{b}": int(v) for (a, b), v in mo.groupby(["rna_call", "bulk_rna"]).size().items()}
for col in ["proteomics", "chromatin", "spatial"]:
    S[f"R15_05_{col}_values"] = mo[col].value_counts().to_dict()
prot = pd.read_csv(PROT, sep="\t")
S["R15_05_protein_rows"] = len(prot)
S["R15_05_protein_bh_lt_0.05"] = int((prot.padj < 0.05).sum())
pdict = prot.drop_duplicates("gene").set_index("gene")
symdup = g.symbol.value_counts()
mo["protein_measured"] = mo.symbol.isin(pdict.index) & (mo.symbol.map(symdup) == 1)
mo["protein_log2fc"] = mo.symbol.map(pdict.logFC).where(mo.protein_measured)
mo["protein_bh_q"] = mo.symbol.map(pdict.padj).where(mo.protein_measured)
mo["bulk_rna_corrected"] = mo.rna_call
mo["proteomics_corrected"] = np.select([~mo.protein_measured, mo.protein_bh_q < 0.05], ["not_measured", "supported"], "indeterminate")
mo["chromatin_corrected"] = "not_joined_at_gene_level"
mo["spatial_corrected"] = "not_joined_at_gene_level"
S["R15_05_protein_measured_genes"] = int(mo.protein_measured.sum())
S["R15_05_protein_supported_genes"] = int((mo.proteomics_corrected == "supported").sum())
S["R15_05_protein_symbols_unmatched"] = int((~prot.gene.isin(g.symbol)).sum())
S["R15_05_GNMT"] = mo[mo.symbol == "GNMT"][["proteomics", "proteomics_corrected", "protein_log2fc", "protein_bh_q"]].to_dict("records")
mo[["ensembl_id_base", "symbol", "rna_call", "bulk_rna", "bulk_rna_corrected", "proteomics", "proteomics_corrected",
    "protein_log2fc", "protein_bh_q", "chromatin", "chromatin_corrected", "spatial", "spatial_corrected"]].to_csv(
    OUT / "R15_05_assay_applicability_corrected.tsv.gz", sep="\t", index=False)
S["R15_05_protein_by_class"] = {f"{a}|{b}": int(v) for (a, b), v in mo.groupby(["primary_evidence_class", "proteomics_corrected"]).size().items()}

# ---------- R07-11: six genes keyed on bulk ID with a genetics row from another ID ----------
fj = pd.read_csv(COMP / "identity_reconciliation_v3/source_symbol_false_joint_links.tsv", sep="\t", dtype=str)
six = fj[fj.catalog_id.notna() & (fj.catalog_id != "")].copy()
six["bulk_id_has_coloc_t12"] = six.ensembl_bulk.isin(cg_ids)
six["genetic_id_has_coloc_t12"] = six.ensembl_genetic.isin(cg_ids)
six["bulk_id_in_any_coloc"] = six.ensembl_bulk.isin(all_ids)
six["catalog_genetics_call"] = six.catalog_id.map(gen.call_state)
six["catalog_genetics_estimate"] = six.catalog_id.map(gen.estimate)
six["catalog_rna_call"] = six.catalog_id.map(rna.call_state)
six["catalog_class"] = six.catalog_id.map(gi.primary_evidence_class)
six.to_csv(OUT / "R07_11_six_genes.tsv", sep="\t", index=False)
S["R07_11_n"] = len(six)
S["R07_11_rows"] = six[["symbol", "catalog_id", "ensembl_genetic", "bulk_id_has_coloc_t12", "bulk_id_in_any_coloc", "catalog_genetics_call", "catalog_rna_call", "catalog_class"]].to_dict("records")

# ---------- R15-09: program context table ----------
pc = pd.read_parquet(CAT / "passport_program_context.parquet")
pc.to_csv(OUT / "R15_09_program_context_rows.tsv", sep="\t", index=False)
S["R15_09_rows"] = len(pc)
S["R15_09_columns"] = pc.columns.tolist()
pi = pd.read_parquet(CAT / "passport_program_index.parquet")
S["R15_09_program_index_rows"] = len(pi)
S["R15_09_program_index_release"] = pi.filter(like="release").iloc[:, 0].value_counts().to_dict() if pi.filter(like="release").shape[1] else {}

# ---------- Fig 6B/6D corrected source tables ----------
S["class_counts"] = g.primary_evidence_class.value_counts().to_dict()
json.dump(S, open(OUT / "f6_summary.json", "w"), indent=1, default=str)
print(json.dumps(S, indent=1, default=str))
