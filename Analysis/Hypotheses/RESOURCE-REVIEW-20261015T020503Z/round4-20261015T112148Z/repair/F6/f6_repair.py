#!/usr/bin/env python3
"""F6 repair round: re-derive the values the blind checker disputed, from release files.

Items:
  A  R07-02  the seven genetics-untestable genes (insufficient shared posterior)
  B  R07-02  PP.H3 > 0.8 with best PP.H4 <= 0.5, Tier-1/2 studies only (Catalog call basis)
  C  R07-05  why 114 bulk, 61 Tier-1/2 COLOC and 76 program-member IDs have no Catalog record
  D  R15-03  per-gene sidecar rebuilt: full-precision q (Catalog) and p (recomputed from t),
             SuSiE and ABF ancestry of the driving study only
  E  R15-05  molecular-object assay_applicability candidate, proteomics as PXD051911 coverage
             using the H86 protein-row assignment (no protein call; H86 prespec keeps
             protein significance out of the Catalog summary)
  F  R10-03  Myojin class test with the frozen Plan 40 code (stratified Freedman-Lane
             partial F, 100,000 draws, prespecified seed 2026080740); v1 reproduction first

Reads saved tables only. Writes to OUT (new directory). Seed for the permutation is the
prespecified Plan 40 seed; nothing else is random.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
CAT = ROOT / "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8"
COMP = ROOT / "figures/candidates/catalog-r8-completion-20261014T141645Z_22236527"
BULK = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
COLOC = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
TIERCFG = ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv"
H86 = ROOT / "Analysis/Hypotheses/H86-20261015T125018Z/run-20261015T125926Z/protein_row_assignment.tsv"
MY = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf"
WS = Path(__file__).resolve().parent
OUT = Path(sys.argv[1])
OUT.mkdir(parents=False, exist_ok=False)
S = {}


def base(s):
    return s.astype(str).str.replace(r"\..*$", "", regex=True)


def vc(s):
    return {str(k): int(v) for k, v in s.value_counts(dropna=False).items()}


e = pd.read_parquet(CAT / "passport_evidence_long.parquet")
g = pd.read_parquet(CAT / "passport_gene_index.parquet")
gen = e[e.evidence_domain == "genetics"].set_index("ensembl_id")
rna = e[e.evidence_domain == "transcriptomics"].set_index("ensembl_id")
gi = g.set_index("ensembl_id")
cat_ids = set(g.ensembl_id)
tier = pd.read_csv(TIERCFG, sep="\t").set_index("study_name")["tier"]

t12 = pd.read_csv(COLOC / "gene_level_coloc_tier12.csv")
t12["ens"] = base(t12.ensembl)
S_dup = int(t12.ens.duplicated().sum())
T = t12.drop_duplicates("ens").set_index("ens")
glc = pd.read_csv(COLOC / "gene_level_coloc.csv", usecols=["ensembl"])
allg = pd.read_csv(COLOC / "susie_coloc_all_gwas.csv",
                   usecols=["gwas_name", "gene", "ensembl", "PP.H4.susie", "PP.H4.abf", "ancestry", "method"])
allg["ens"] = base(allg.ensembl)
allg["tier"] = allg.gwas_name.map(tier)
study_anc = allg.groupby("gwas_name").ancestry.agg(lambda s: ";".join(sorted(set(s.dropna().astype(str)))))
S["t12_duplicate_base_ids"] = S_dup
S["study_ancestry_unique"] = bool((allg.groupby("gwas_name").ancestry.nunique() <= 1).all())

# ---------------- A: seven genetics-untestable genes ----------------
gu = gen[gen.call_state == "untestable"]
S["A_gen_untestable_by_reason"] = vc(gu.testability_state)
up = gu.index[gu.testability_state == "underpowered_source"]
S["A_underpowered_in_t12"] = int(pd.Index(up).isin(t12.ens).sum())
S["A_underpowered_in_gene_level_coloc"] = int(pd.Index(up).isin(set(base(glc.ensembl))).sum())
S["A_underpowered_in_all_gwas_master"] = int(pd.Index(up).isin(set(allg.ens)).sum())
ins = gu.index[gu.testability_state == "insufficient_shared_posterior"]
rows = []
for ens in ins:
    a = allg[(allg.ens == ens) & allg.tier.isin([1, 2])]
    fin = a.dropna(subset=["PP.H4.susie"]).sort_values("PP.H4.susie", ascending=False)
    unt = a[a.method == "susie_untestable_insufficient_shared_posterior"]
    rows.append(dict(ensembl_id=ens, symbol=gi.loc[ens, "symbol"],
                     t12_susie_state=T.loc[ens, "susie_state_t12"] if ens in T.index else "absent",
                     t12_susie_untestable_studies=T.loc[ens, "susie_untestable_t12_studies"] if ens in T.index else np.nan,
                     n_t12_studies=a.gwas_name.nunique(), n_t12_finite_susie=len(fin),
                     max_finite_susie_pp4=fin["PP.H4.susie"].iloc[0] if len(fin) else np.nan,
                     max_finite_susie_study=fin.gwas_name.iloc[0] if len(fin) else "",
                     untestable_study=";".join(sorted(unt.gwas_name)),
                     best_abf_pp4=T.loc[ens, "coloc_best_abf_pp4"] if ens in T.index else np.nan,
                     catalog_class=gi.loc[ens, "primary_evidence_class"],
                     rna_call=rna.loc[ens, "call_state"]))
seven = pd.DataFrame(rows).sort_values("symbol")
seven.to_csv(OUT / "A_R07_02_seven_susie_untestable_genes.tsv", sep="\t", index=False)
S["A_seven"] = seven.to_dict("records")

# ---------------- B: PP.H3 among eligible SuSiE signal pairs ----------------
sp = pd.read_csv(COLOC / "susie_coloc_signal_pairs.tsv.gz", sep="\t")
sp["tier"] = sp.gwas.map(tier)
sp["ens"] = base(sp.ensembl)
S["B_sp_tier_counts"] = vc(sp.tier)
el = sp[sp.eligible.astype(str).str.upper() == "TRUE"].copy()


def h3(df, label):
    agg = df.groupby("ens").agg(max_h3=("PP.H3", "max"), max_h4=("PP.H4", "max"), symbol=("gene", "first"))
    top = df.sort_values("PP.H3", ascending=False).drop_duplicates("ens").set_index("ens")
    agg["max_h3_study"] = top.gwas
    agg["max_h3_study_tier"] = top.tier
    hit = agg[(agg.max_h3 > 0.8) & (agg.max_h4 <= 0.5)].copy()
    hit["catalog_genetics_call"] = [gen.call_state.get(x, "absent_from_catalog") for x in hit.index]
    hit["catalog_class"] = [gi.primary_evidence_class.get(x, "absent_from_catalog") for x in hit.index]
    hit["strict_h3_gt0.9_h4_lt0.1"] = (hit.max_h3 > 0.9) & (hit.max_h4 < 0.1)
    S[f"B_{label}_genes_with_eligible_pair"] = len(agg)
    S[f"B_{label}_h3gt0.8_h4le0.5"] = len(hit)
    S[f"B_{label}_by_call"] = vc(hit.catalog_genetics_call)
    S[f"B_{label}_strict"] = int(hit["strict_h3_gt0.9_h4_lt0.1"].sum())
    hit.sort_values("max_h3", ascending=False).to_csv(OUT / f"B_R07_02_pph3_{label}.tsv", sep="\t")
    return hit


h_all = h3(el, "all_tiers")
h_12 = h3(el[el.tier.isin([1, 2])], "tier12")
S["B_only_all_tiers"] = sorted(h_all.symbol[~h_all.index.isin(h_12.index)].tolist())
S["B_only_tier12"] = sorted(h_12.symbol[~h_12.index.isin(h_all.index)].tolist())

# ---------------- C: IDs without a Catalog record ----------------
q = pd.read_csv(CAT / "passport_identity_quarantine.tsv", sep="\t", dtype=str)
q_ids = {}
for _, rq in q.iterrows():
    for col in ["normalized_ensembl_candidates", "source_ensembl_text"]:
        for x in str(rq[col]).split(";"):
            x = x.strip().split(".")[0]
            if x and x != "nan":
                q_ids.setdefault(x, set()).add(rq.reason_code)
S["C_quarantine_reason_codes"] = vc(q.reason_code)


def reason(x):
    return ";".join(sorted(q_ids.get(x, []))) or "not_quarantined"


c = pd.read_csv(BULK / "coefficient_table.tsv.gz", sep="\t")
c["ens"] = base(c.gene)
c["canonical"] = (c.padj < 0.05) & (c.logFC.abs() > 0.5)
ab = c[~c.ens.isin(cat_ids)].copy()
ab["quarantine_reason"] = ab.ens.map(reason)
S["C_bulk_absent"] = len(ab)
S["C_bulk_absent_by_reason"] = vc(ab.quarantine_reason)
S["C_bulk_absent_q_lt_0.05"] = int((ab.padj < 0.05).sum())
S["C_bulk_absent_canonical"] = ab[ab.canonical][["ens", "symbol", "logFC", "padj"]].to_dict("records")
S["C_bulk_absent_not_quarantined"] = ab[ab.quarantine_reason == "not_quarantined"][["ens", "symbol", "padj"]].to_dict("records")
sym_n = c.groupby("symbol").ens.nunique()
S["C_bulk_absent_symbol_shared_in_bulk"] = int(ab.symbol.map(sym_n).gt(1).sum())
ab.to_csv(OUT / "C_R07_05_bulk_ids_absent.tsv", sep="\t", index=False)
fj = pd.read_csv(COMP / "identity_reconciliation_v3/source_symbol_false_joint_links.tsv", sep="\t", dtype=str)
r0711_genetic = set(fj.loc[fj.catalog_id.notna() & (fj.catalog_id != ""), "ensembl_genetic"])
ta = t12[~t12.ens.isin(cat_ids)].copy()
ta["quarantine_reason"] = ta.ens.map(reason)
ta["is_R07_11_genetic_id"] = ta.ens.isin(r0711_genetic)
S["C_t12_ids"] = int(t12.ens.nunique())
S["C_t12_absent"] = len(ta)
S["C_t12_absent_by_reason"] = vc(ta.quarantine_reason)
S["C_t12_absent_not_quarantined"] = ta[ta.quarantine_reason == "not_quarantined"][["ens", "is_R07_11_genetic_id", "coloc_best_susie_pp4", "coloc_best_abf_pp4"]].to_dict("records")
ta.to_csv(OUT / "C_R07_05_t12_ids_absent.tsv", sep="\t", index=False)
pm = pd.read_parquet(CAT / "passport_program_membership.parquet")
pm["ens"] = base(pm.ensembl_id)
pa = pm[~pm.ens.isin(cat_ids)].copy()
pa["quarantine_reason"] = pa.ens.map(reason)
S["C_program_member_ids"] = int(pm.ens.nunique())
S["C_program_member_ids_in_catalog"] = int(pm.ens[pm.ens.isin(cat_ids)].nunique())
S["C_program_member_ids_absent"] = int(pa.ens.nunique())
S["C_program_memberships_absent"] = len(pa)
S["C_program_absent_by_reason"] = vc(pa.drop_duplicates("ens").quarantine_reason)
S["C_program_absent_in_bulk"] = int(pa.ens.drop_duplicates().isin(c.ens).sum())
S["C_program_absent_in_any_coloc"] = int(pa.ens.drop_duplicates().isin(set(allg.ens) | set(t12.ens) | set(base(glc.ensembl))).sum())
S["C_catalog_rna_supported"] = int((rna.call_state == "supported").sum())
S["C_bulk_canonical"] = int(c.canonical.sum())

# ---------------- D: sidecar rebuild ----------------
mid = c[(c["P.Value"] > 1e-6) & (c["P.Value"] < 0.5)]


def logerr(df):
    return float(np.max(np.abs(np.log(2 * stats.t.sf(np.abs(mid.t), df)) - np.log(mid["P.Value"]))))


fit = optimize.minimize_scalar(logerr, bounds=(100, 1e5), method="bounded", options={"xatol": 1e-6})
df_total = float(fit.x)
S["D_df_total_fitted"] = df_total
S["D_df_total_max_log_err"] = float(fit.fun)
c["p_full"] = 2 * stats.t.sf(np.abs(c.t.to_numpy()), df_total)
chk = c["P.Value"] >= 1e-10
S["D_p_full_vs_table_max_rel_err_p_ge_1e-10"] = float(((c.p_full - c["P.Value"]).abs() / c["P.Value"])[chk].max())
p = c.p_full.to_numpy()
o = np.argsort(p)
qq = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
q_full = np.empty_like(qq)
q_full[o] = np.minimum(qq, 1)
c["q_full"] = q_full
cq = rna.q_value.astype(float)
jj = pd.DataFrame({"q_cat": cq}).join(c.set_index("ens")[["q_full", "padj"]], how="inner")
S["D_q_full_vs_catalog_max_rel_err"] = float(((jj.q_full - jj.q_cat).abs() / jj.q_cat).max())
S["D_catalog_q_zero"] = int((cq == 0).sum())
S["D_table_padj_zero"] = int((c.padj == 0).sum())
S["D_table_p_zero"] = int((c["P.Value"] == 0).sum())
cb = c.set_index("ens")
side = pd.DataFrame(index=pd.Index(g.ensembl_id, name="ensembl_id"))
side["symbol"] = gi.symbol
side["primary_evidence_class"] = gi.primary_evidence_class
side["genetics_call"] = gen.call_state
side["genetics_testability"] = gen.testability_state
has_s = T.coloc_best_susie_pp4.notna()
has_a = T.coloc_best_abf_pp4.notna()
side["coloc_best_susie_pp4"] = T.coloc_best_susie_pp4
side["coloc_best_susie_study"] = T.susie_driving_gwas.where(has_s)
side["coloc_best_susie_ancestry"] = T.susie_driving_gwas.where(has_s).map(study_anc)
side["coloc_best_abf_pp4"] = T.coloc_best_abf_pp4
side["coloc_best_abf_study"] = T.abf_driving_gwas.where(has_a)
side["coloc_best_abf_ancestry"] = T.abf_driving_gwas.where(has_a).map(study_anc)
side["coloc_call_rule"] = "supported iff best Tier-1/2 SuSiE PP.H4 > 0.5; untestable if any Tier-1/2 study is SuSiE-untestable; ABF not used for the call"
est = gen.estimate.astype(float)
side["coloc_displayed_estimate_method"] = np.select(
    [est.reindex(side.index).isna(),
     np.isclose(est.reindex(side.index), side.coloc_best_susie_pp4, rtol=0, atol=1e-12) & side.coloc_best_susie_pp4.notna(),
     np.isclose(est.reindex(side.index), side.coloc_best_abf_pp4, rtol=0, atol=1e-12) & side.coloc_best_abf_pp4.notna()],
    ["none", "SuSiE", "ABF"], "other_id")
side["rna_call"] = rna.call_state
side["rna_log2fc"] = cb.logFC
side["rna_se"] = cb.SE
side["rna_t"] = cb.t
side["rna_p"] = cb.p_full
side["rna_bh_q"] = cq
side["rna_direction_resolved_bh05"] = side.rna_bh_q < 0.05
tested = side.rna_log2fc.notna()
side["rna_n_units"] = np.where(tested, 774, np.nan)
side["rna_n_control_units"] = np.where(tested, 138, np.nan)
side["rna_n_disease_units"] = np.where(tested, 636, np.nan)
side["rna_unit"] = np.where(tested, "genetic-identity cluster (not verified donor)", "")
side["rna_p_note"] = np.where(tested, f"recomputed as 2*pt(-|t|, {df_total:.4f}); BH q from the Catalog (full precision)", "")
side.reset_index().to_csv(OUT / "R15_03_gene_effect_sidecar_repair.tsv.gz", sep="\t", index=False)
S["D_side_rows"] = len(side)
S["D_side_q_zero"] = int((side.rna_bh_q == 0).sum())
S["D_side_p_zero"] = int((side.rna_p == 0).sum())
S["D_side_rna_tested"] = int(tested.sum())
S["D_side_q_nonnull"] = int(side.rna_bh_q.notna().sum())
S["D_side_susie_anc_without_susie"] = int((side.coloc_best_susie_ancestry.notna() & side.coloc_best_susie_pp4.isna()).sum())
S["D_side_susie_with_anc"] = int(side.coloc_best_susie_ancestry.notna().sum())
S["D_side_susie_values"] = int(side.coloc_best_susie_pp4.notna().sum())
S["D_side_abf_with_anc"] = int(side.coloc_best_abf_ancestry.notna().sum())
S["D_side_abf_values"] = int(side.coloc_best_abf_pp4.notna().sum())
S["D_side_method_by_call"] = {f"{a}|{b}": int(n) for (a, b), n in side.groupby(["genetics_call", "coloc_displayed_estimate_method"]).size().items()}
for sym in ["GNMT", "THRB", "HKDC1", "MTARC1", "GLP1R", "FIRRM", "MAD1L1"]:
    S[f"D_{sym}"] = side[side.symbol == sym].reset_index().to_dict("records")

# ---------------- E: molecular-object assay_applicability with H86 protein coverage ----------------
mo = pd.DataFrame([json.loads(line) for line in open(COMP / "molecular_objects/molecular_object_sidecar.jsonl")])
aa = mo.assay_applicability.apply(lambda x: json.loads(x) if isinstance(x, str) else x).apply(pd.Series)
mo = pd.concat([mo[["ensembl_id_base", "symbol", "rna_call", "genetics_call", "primary_evidence_class"]], aa], axis=1)
S["E_records"] = len(mo)
S["E_rna_call_x_bulk_rna"] = {f"{a}|{b}": int(n) for (a, b), n in mo.groupby(["rna_call", "bulk_rna"]).size().items()}
for col in ["proteomics", "chromatin", "spatial"]:
    S[f"E_{col}_values"] = vc(mo[col])
pr = pd.read_csv(H86, sep="\t")
S["E_h86_routes"] = vc(pr.route)
symdup = g.symbol.value_counts()
pq = pr[pr.route.isin(["direct", "remapped"]) & pr.catalog_symbol.isin(symdup[symdup == 1].index)]
S["E_h86_quantified_catalog_genes"] = int(pq.catalog_symbol.nunique())
S["E_h86_quantified_by_route"] = vc(pq.route)
route = pq.drop_duplicates("catalog_symbol").set_index("catalog_symbol").route
mo["bulk_rna_corrected"] = mo.rna_call
mo["proteomics_corrected"] = np.where(mo.symbol.isin(route.index), "quantified_PXD051911", "not_quantified_PXD051911")
mo["protein_symbol_route"] = mo.symbol.map(route).fillna("")
mo["chromatin_corrected"] = "not_joined_at_gene_level"
mo["spatial_corrected"] = "not_joined_at_gene_level"
S["E_proteomics_corrected"] = vc(mo.proteomics_corrected)
S["E_quantified_by_class"] = {k: int(v) for k, v in mo[mo.proteomics_corrected == "quantified_PXD051911"].primary_evidence_class.value_counts().items()}
mo[["ensembl_id_base", "symbol", "primary_evidence_class", "rna_call", "bulk_rna", "bulk_rna_corrected", "proteomics",
    "proteomics_corrected", "protein_symbol_route", "chromatin", "chromatin_corrected", "spatial", "spatial_corrected"]].to_csv(
    OUT / "R15_05_assay_applicability_corrected_repair.tsv.gz", sep="\t", index=False)

# ---------------- F: Myojin class test, frozen Plan 40 code ----------------
spec = importlib.util.spec_from_file_location("phase_c_stats", WS / "phase_c_stats_copy.py")
pcs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pcs)
sr = pd.read_csv(MY / "screen_results.tsv", sep="\t", dtype=str, keep_default_na=False)
sr = sr[sr.primary_eligible.str.strip().str.lower().isin(["1", "true", "t", "yes", "y"])]
fr = pd.DataFrame({
    "source_row": sr.source_row.astype(int).to_numpy(),
    "gene_symbol": sr.gene_symbol.to_numpy(),
    "primary_evidence_class": sr.primary_evidence_class.to_numpy(),
    "guide_count_min": sr.guide_count_min.astype(int).to_numpy(),
    "log1p_HLF_TPM": [float(x) for x in sr.log1p_HLF_TPM],
    "HLF_Chronos": [np.nan if x == "" else float(x) for x in sr.HLF_Chronos],
    "HLF_Chronos_missing": [int(x == "") for x in sr.HLF_Chronos],
    "PA_vs_vehicle_LFC": [float(x) for x in sr.PA_vs_vehicle_positive_selection_LFC],
    "frozen_guide_bin": sr.analysis_guide_bin.to_numpy(),
    "frozen_expr_q": sr.analysis_expression_quintile.to_numpy(),
})
fr = fr.sort_values("gene_symbol", kind="stable").reset_index(drop=True)
S["F_primary_n"] = len(fr)
prim = pcs.add_frozen_strata(fr.drop(columns=["frozen_guide_bin", "frozen_expr_q"]))
S["F_strata_match_written_bins"] = bool(
    (prim.guide_bin.astype(str).to_numpy() == fr.frozen_guide_bin.to_numpy()).all()
    and (prim.expression_quintile.astype(str).to_numpy() == fr.frozen_expr_q.to_numpy()).all())
SEED, DRAWS = 2026080740, 100_000
res = []
r1, a1 = pcs.fit_class_continuous(prim.copy(), "PA_vs_vehicle_LFC", "v1_reproduction", SEED, accepted_draws=DRAWS)
res += r1
stored = pd.read_csv(MY / "class_effects.tsv", sep="\t")
st = stored[(stored.analysis_variant == "primary") & (stored.contrast == "omnibus_evidence_class")].iloc[0]
S["F_v1_F"] = float(r1[0]["statistic"])
S["F_v1_perm_p"] = float(r1[0]["permutation_p"])
S["F_v1_stored_F"] = float(st.statistic)
S["F_v1_stored_perm_p"] = float(st.permutation_p)
S["F_v1_reproduced_exactly"] = bool(abs(S["F_v1_F"] - S["F_v1_stored_F"]) < 1e-9 and S["F_v1_perm_p"] == S["F_v1_stored_perm_p"])
gu2 = pd.read_csv(MY / "screen_gene_universe.tsv", sep="\t", dtype=str, usecols=["source_row", "gencode_ensembl_base"])
gu2["source_row"] = gu2.source_row.astype(int)
ensmap = gu2.set_index("source_row").gencode_ensembl_base.str.replace(r"\..*$", "", regex=True)
prim["ens"] = prim.source_row.map(ensmap)
prim["r8_class"] = prim.ens.map(gi.primary_evidence_class).fillna("not_in_r8")
S["F_r8_classes"] = vc(prim.r8_class)
S["F_v1_x_r8"] = {f"{a}|{b}": int(n) for (a, b), n in prim.groupby(["primary_evidence_class", "r8_class"]).size().items()}
to_spec = {"unresolved": "neither", "genetically_anchored": "genetic_only",
           "established_state_associated": "disease_state_only", "concordant": "convergent", "untested": "neither"}
# r8_spec: eligible universe = the four spec classes; untested (= not_jointly_testable) excluded,
# strata inherited from the frozen 7,762-gene universe (executor rule for subsets)
r8a = prim[~prim.r8_class.isin(["untested", "not_in_r8"])].copy()
r8a["primary_evidence_class"] = r8a.r8_class.map(to_spec)
r8a = r8a.drop(columns=["ens", "r8_class"]).reset_index(drop=True)
ra, aa_ = pcs.fit_class_continuous(r8a, "PA_vs_vehicle_LFC", "r8_spec_universe_untested_excluded", SEED, accepted_draws=DRAWS)
res += ra
r8b = prim[prim.r8_class != "not_in_r8"].copy()
r8b["primary_evidence_class"] = r8b.r8_class.map(to_spec)
r8b = r8b.drop(columns=["ens", "r8_class"]).reset_index(drop=True)
rb, ab_ = pcs.fit_class_continuous(r8b, "PA_vs_vehicle_LFC", "r8_untested_in_reference", SEED, accepted_draws=DRAWS)
res += rb
rr = pd.DataFrame(res)
rr["contrast"] = rr.contrast.replace({"genetic_only_vs_neither": "genetically_anchored_vs_unresolved",
                                      "disease_state_only_vs_neither": "established_state_associated_vs_unresolved",
                                      "convergent_vs_neither": "concordant_vs_unresolved"}).where(rr.analysis_variant != "v1_reproduction", rr.contrast)
rr.to_csv(OUT / "F_R10_03_myojin_frozen_code_refit.tsv", sep="\t", index=False)
pd.DataFrame([a1, aa_, ab_]).to_csv(OUT / "F_R10_03_permutation_audit.tsv", sep="\t", index=False)
for r in res:
    S[f"F_{r['analysis_variant']}_{r['contrast']}"] = {k: r[k] for k in ["n_total", "n_class", "estimate", "CI95_low", "CI95_high", "standardized_effect", "statistic", "permutation_p", "BH_q", "df_residual"]}
# classical parametric p of the partial F, for reference
for r in res:
    if r["contrast"] == "omnibus_evidence_class":
        S[f"F_{r['analysis_variant']}_partialF_parametric_p"] = float(stats.f.sf(r["statistic"], r["df_numerator"], r["df_residual"]))

json.dump(S, open(OUT / "f6_repair_summary.json", "w"), indent=1, default=str)
print(json.dumps(S, indent=1, default=str))
