#!/usr/bin/env python3
"""Blind re-check of F6 repaired values from release files (r8 Catalog, COLOC r2, 774-unit bulk fit).

Writes recheck.json to the output directory given as argv[1]. Seed 20260923 is set
for completeness; nothing here is random except the frozen Myojin refit, which uses
its own frozen seed 2026080740 from the blind analysis spec.
"""
import gzip
import importlib.util
import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

np.random.seed(20260923)
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
CAT = ROOT / "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8"
COMP = ROOT / "figures/candidates/catalog-r8-completion-20261014T141645Z_22236527"
BULK = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
COLOC = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
TIERCFG = ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv"
H86 = ROOT / "Analysis/Hypotheses/H86-20261015T125018Z/run-20261015T125926Z/protein_row_assignment.tsv"
MY = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf"
FROZEN = ROOT / "Analysis/Multimodal_Program_Projection/scripts/myojin_hlf/phase_c_stats.py"
SPAT = ROOT / "Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11/coverage"
ATAC = ROOT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/da"
GTF = Path("/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz")
PAN = ROOT / "figures/candidates/review-r3-fixes-20261015T112148Z/F6/repair/panels_20261015T173321Z"
OUT = Path(sys.argv[1])
RUN_PERM = len(sys.argv) > 2 and sys.argv[2] == "perm"
R = {}


def b(s):
    return s.astype(str).str.split(".").str[0]


def vc(s):
    return {str(k): int(v) for k, v in s.value_counts(dropna=False).items()}


# ---------- Catalog ----------
e = pd.read_parquet(CAT / "passport_evidence_long.parquet")
g = pd.read_parquet(CAT / "passport_gene_index.parquet")
R["n_evidence_rows"] = len(e)
R["n_genes"] = len(g)
gen = e[e.evidence_domain == "genetics"].set_index("ensembl_id")
rna = e[e.evidence_domain == "transcriptomics"].set_index("ensembl_id")
R["gen_rows"], R["rna_rows"] = len(gen), len(rna)
cls = g.set_index("ensembl_id").primary_evidence_class
R["class_counts"] = vc(cls)
j = pd.DataFrame({"gen": gen.call_state, "rna": rna.call_state, "cls": cls})
R["matrix_gen_x_rna"] = {f"{a}|{c}": int(n) for (a, c), n in j.groupby(["gen", "rna"]).size().items()}
R["matrix_class_map"] = {f"{a}|{c}|{k}": int(n) for (a, c, k), n in j.groupby(["gen", "rna", "cls"]).size().items()}
R["gen_untestable_testability"] = vc(gen.testability_state[gen.call_state == "untestable"])
R["rna_untestable_testability"] = vc(rna.testability_state[rna.call_state == "untestable"])
R["call_state_all_rows"] = vc(e.call_state)
R["not_applicable_rows"] = int((e.call_state == "not_applicable").sum())
R["source_dependent_true_rows"] = int((e.source_dependent == True).sum())  # noqa: E712
R["source_dependent_values"] = vc(e.source_dependent.astype(str))
for col in ["negative_call_rule_id", "negative_decision_boundary", "negative_margin", "negative_call_passed",
            "source_negative_call_rule_id", "source_negative_call_passed"]:
    R[f"nonnull_{col}"] = int(e[col].notna().sum())
for col in ["standard_error", "ci_lower", "ci_upper", "p_value", "n_biological_units", "n_technical_units"]:
    R[f"nonnull_{col}"] = int(e[col].notna().sum())
R["nonnull_q_by_domain"] = {d: int(x.q_value.notna().sum()) for d, x in e.groupby("evidence_domain")}
R["nonnull_estimate_by_domain"] = {d: int(x.estimate.notna().sum()) for d, x in e.groupby("evidence_domain")}
R["biological_unit_values"] = vc(e.biological_unit)
R["tested_negative_rows"] = int(e.call_state.astype(str).str.contains("negative").sum())

# RNA indeterminate with q<0.05 and |logFC|<=0.5 (Catalog values)
ri = rna[rna.call_state == "indeterminate"]
R["rna_indet_q05_absfc_le05_catalog"] = int(((ri.q_value.astype(float) < 0.05) & (ri.estimate.astype(float).abs() <= 0.5)).sum())
R["rna_supported_n"] = int((rna.call_state == "supported").sum())
rs = rna[rna.call_state == "supported"]
R["rna_supported_rule_check"] = int(((rs.q_value.astype(float) < 0.05) & (rs.estimate.astype(float).abs() > 0.5)).sum())

# ---------- bulk 774-unit fit ----------
c = pd.read_csv(BULK / "coefficient_table.tsv.gz", sep="\t")
c["ens"] = b(c.gene)
R["bulk_rows"] = len(c)
R["bulk_dup_base"] = int(c.ens.duplicated().sum())
c["canon"] = (c.padj < 0.05) & (c.logFC.abs() > 0.5)
R["bulk_canonical"] = int(c.canon.sum())
cb = c.set_index("ens")
cg = pd.read_csv(BULK / "cohort_group_counts.tsv", sep="\t")
R["units_control"] = int(cg.new_control.sum())
R["units_disease"] = int(cg.new_disease.sum())
ci = c[~c.canon]
R["bulk_noncanon_q05_absfc_le05"] = int(((ci.padj < 0.05) & (ci.logFC.abs() <= 0.5)).sum())
# restricted to Catalog indeterminate IDs
R["bulk_indet_q05_absfc_le05_in_catalog_indet"] = int(((cb.reindex(ri.index).padj < 0.05) & (cb.reindex(ri.index).logFC.abs() <= 0.5)).sum())
R["rna_untestable_ids_in_bulk"] = int(pd.Index(rna.index[rna.call_state == "untestable"]).isin(c.ens).sum())

# ---------- COLOC r2 ----------
tier = pd.read_csv(TIERCFG, sep="\t").set_index("study_name").tier
t12 = pd.read_csv(COLOC / "gene_level_coloc_tier12.csv")
t12["ens"] = b(t12.ensembl)
R["t12_rows"], R["t12_dup"] = len(t12), int(t12.ens.duplicated().sum())
T = t12.set_index("ens")
glc = pd.read_csv(COLOC / "gene_level_coloc.csv")
glc["ens"] = b(glc.ensembl)
allg = pd.read_csv(COLOC / "susie_coloc_all_gwas.csv",
                   usecols=["gwas_name", "gene", "ensembl", "PP.H4.susie", "PP.H4.abf", "PP.H3.susie", "ancestry", "method"])
allg["ens"] = b(allg.ensembl)
allg["tier"] = allg.gwas_name.map(tier)
R["allg_unmapped_tier_studies"] = sorted(allg.gwas_name[allg.tier.isna()].unique().tolist())
anc = allg.groupby("gwas_name").ancestry.agg(lambda s: sorted(set(s.dropna().astype(str))))
R["study_ancestry_multi"] = {k: v for k, v in anc.items() if len(v) != 1}

gu_up = gen.index[(gen.call_state == "untestable") & (gen.testability_state == "underpowered_source")]
R["underpowered_n"] = len(gu_up)
R["underpowered_in_t12"] = int(pd.Index(gu_up).isin(t12.ens).sum())
R["underpowered_in_gene_level_coloc"] = int(pd.Index(gu_up).isin(glc.ens).sum())
R["underpowered_in_allg"] = int(pd.Index(gu_up).isin(allg.ens).sum())
R["underpowered_in_allg_t12"] = int(pd.Index(gu_up).isin(allg.ens[allg.tier.isin([1, 2])]).sum())
R["underpowered_rna_call"] = vc(rna.call_state.reindex(gu_up))
gu_isp = gen.index[(gen.call_state == "untestable") & (gen.testability_state != "underpowered_source")]
seven = []
for x in gu_isp:
    a = allg[(allg.ens == x) & allg.tier.isin([1, 2])]
    seven.append(dict(ens=x, sym=g.set_index("ensembl_id").symbol.get(x), state=T.susie_state_t12.get(x),
                      n_unt=int(T.susie_untestable_t12_studies.get(x, -1)),
                      unt_studies=sorted(a.gwas_name[a.method.astype(str).str.contains("untestable")].unique().tolist()),
                      rna=rna.call_state.get(x), cls=cls.get(x)))
R["seven"] = seven
R["t12_susie_state_counts"] = vc(T.susie_state_t12)

# PP.H3 from eligible SuSiE signal pairs
sp = pd.read_csv(COLOC / "susie_coloc_signal_pairs.tsv.gz", sep="\t")
sp["ens"] = b(sp.ensembl)
sp["tier"] = sp.gwas.map(tier)
R["sp_unmapped_tier"] = int(sp.tier.isna().sum())
el = sp[sp.eligible.astype(str).str.upper() == "TRUE"]
el12 = el[el.tier.isin([1, 2])]
agg = el12.groupby("ens").agg(h3=("PP.H3", "max"), h4=("PP.H4", "max"), sym=("gene", "first"))
R["pph3_t12_genes_eligible"] = len(agg)
hit = agg[(agg.h3 > 0.8) & (agg.h4 <= 0.5)]
R["pph3_t12_hits"] = len(hit)
R["pph3_t12_hits_by_gen_call"] = vc(pd.Series([gen.call_state.get(x, "absent") for x in hit.index]))
R["pph3_t12_hit_symbols"] = sorted(hit.sym.tolist())
# alternative: best PP.H4 from gene-level SuSiE column
alt = agg.join(T[["coloc_best_susie_pp4"]], how="left")
althit = alt[(alt.h3 > 0.8) & ~(alt.coloc_best_susie_pp4 > 0.5)]
R["pph3_t12_hits_alt_bestsusie"] = len(althit)
R["pph3_t12_genes_eligible_supported"] = int(pd.Series([gen.call_state.get(x) for x in agg.index]).eq("supported").sum())
# check of PP.H3>0.8 & H4<=0.5 within the same pair
pair = el12[(el12["PP.H3"] > 0.8)]
R["pph3_t12_genes_pair_h3gt08"] = int(pair.ens.nunique())

# worked entries
ws = {}
for sym in ["GLP1R", "HKDC1", "MTARC1", "THRB", "GNMT"]:
    ids = g.ensembl_id[g.symbol == sym].tolist()
    for x in ids:
        r = T.loc[x] if x in T.index else None
        a = allg[(allg.ens == x) & allg.tier.isin([1, 2])]
        sus = a.dropna(subset=["PP.H4.susie"]).sort_values("PP.H4.susie", ascending=False)
        abf = a.dropna(subset=["PP.H4.abf"]).sort_values("PP.H4.abf", ascending=False)
        ws[sym] = dict(ens=x, cls=cls.get(x), gen_call=gen.call_state.get(x), gen_est=float(gen.estimate.get(x)) if pd.notna(gen.estimate.get(x)) else None,
                       rna_call=rna.call_state.get(x),
                       t12_susie=None if r is None or pd.isna(r.coloc_best_susie_pp4) else float(r.coloc_best_susie_pp4),
                       t12_susie_study=None if r is None else str(r.susie_driving_gwas),
                       t12_abf=None if r is None or pd.isna(r.coloc_best_abf_pp4) else float(r.coloc_best_abf_pp4),
                       t12_abf_study=None if r is None else str(r.abf_driving_gwas),
                       allg_top_susie=(float(sus["PP.H4.susie"].iloc[0]), sus.gwas_name.iloc[0], sus.ancestry.iloc[0]) if len(sus) else None,
                       allg_top_abf=(float(abf["PP.H4.abf"].iloc[0]), abf.gwas_name.iloc[0], abf.ancestry.iloc[0]) if len(abf) else None,
                       logFC=float(cb.logFC.get(x)) if x in cb.index else None,
                       padj=float(cb.padj.get(x)) if x in cb.index else None,
                       cat_q=float(rna.q_value.get(x)) if pd.notna(rna.q_value.get(x)) else None)
R["worked"] = ws

# ---------- R15-02 displayed PP.H4 provenance ----------
est = gen.estimate.astype(float)
gi = gen.call_state
Tsus = T.coloc_best_susie_pp4.reindex(gen.index)
Tabf = T.coloc_best_abf_pp4.reindex(gen.index)
is_sus = np.isclose(est, Tsus, rtol=0, atol=1e-12) & Tsus.notna()
is_abf = np.isclose(est, Tabf, rtol=0, atol=1e-12) & Tabf.notna() & ~is_sus
meth = pd.Series(np.select([est.isna(), is_sus, is_abf], ["none", "SuSiE", "ABF"], "other"), index=gen.index)
R["displayed_method_by_call"] = {f"{a}|{m}": int(n) for (a, m), n in pd.DataFrame({"a": gi, "m": meth}).groupby(["a", "m"]).size().items()}
ind = gi == "indeterminate"
R["indet_n"] = int(ind.sum())
R["indet_est_gt05"] = int((ind & (est > 0.5)).sum())
R["indet_est_ge08"] = int((ind & (est >= 0.8)).sum())
R["indet_est_gt05_finite_susie_le05"] = int((ind & (est > 0.5) & Tsus.notna() & (Tsus <= 0.5)).sum())
R["indet_est_gt05_by_method"] = vc(meth[ind & (est > 0.5)])
R["supported_n"] = int((gi == "supported").sum())

# sidecar file checks
sc = pd.read_csv(PAN / "R15_03_gene_effect_sidecar_repair.tsv.gz", sep="\t")
R["side_rows"] = len(sc)
R["side_q_zero"] = int((sc.rna_bh_q == 0).sum())
R["side_p_zero"] = int((sc.rna_p == 0).sum())
R["side_q_nonnull"] = int(sc.rna_bh_q.notna().sum())
R["side_susie_n"] = int(sc.coloc_best_susie_pp4.notna().sum())
R["side_susie_anc_n"] = int(sc.coloc_best_susie_ancestry.notna().sum())
R["side_susie_anc_given_value"] = int((sc.coloc_best_susie_pp4.notna() & sc.coloc_best_susie_ancestry.notna()).sum())
R["side_abf_n"] = int(sc.coloc_best_abf_pp4.notna().sum())
R["side_abf_anc_n"] = int(sc.coloc_best_abf_ancestry.notna().sum())
scx = sc.set_index("ensembl_id")
R["side_q_vs_catalog_maxabs"] = float((scx.rna_bh_q - rna.q_value.astype(float).reindex(scx.index)).abs().max())
R["side_susie_vs_t12_maxabs"] = float((scx.coloc_best_susie_pp4 - T.coloc_best_susie_pp4.reindex(scx.index)).abs().max())
R["side_abf_vs_t12_maxabs"] = float((scx.coloc_best_abf_pp4 - T.coloc_best_abf_pp4.reindex(scx.index)).abs().max())
R["side_lfc_vs_bulk_maxabs"] = float((scx.rna_log2fc - cb.logFC.reindex(scx.index)).abs().max())
R["side_se_vs_bulk_maxabs"] = float((scx.rna_se - cb.SE.reindex(scx.index)).abs().max())
# independent p recompute check: single df across genes?
mid = c[(c["P.Value"] > 1e-8) & (c["P.Value"] < 0.9)]
dfs = []
for df in [700, 760, 772, 772.1545, 773, 800, 1e5]:
    dfs.append((df, float(np.max(np.abs(np.log(2 * stats.t.sf(np.abs(mid.t), df)) - np.log(mid["P.Value"]))))))
R["p_df_scan_maxlogerr"] = dfs
pp = 2 * stats.t.sf(np.abs(scx.rna_t.dropna()), 772.1545103460572)
R["side_p_vs_recalc_maxrel"] = float(np.max(np.abs(pp - scx.rna_p.dropna()) / scx.rna_p.dropna()))
R["side_GNMT"] = sc[sc.symbol == "GNMT"].to_dict("records")
R["side_ancestry_values_susie"] = vc(sc.coloc_best_susie_ancestry.dropna())
R["side_ancestry_values_abf"] = vc(sc.coloc_best_abf_ancestry.dropna())
# study->ancestry agreement of the sidecar with all_gwas
st_anc = {k: v[0] for k, v in anc.items() if len(v) == 1}
mm = sc.dropna(subset=["coloc_best_abf_study"])
R["side_abf_anc_mismatch"] = int((mm.coloc_best_abf_study.map(st_anc) != mm.coloc_best_abf_ancestry).sum())
mm = sc.dropna(subset=["coloc_best_susie_study"])
R["side_susie_anc_mismatch"] = int((mm.coloc_best_susie_study.map(st_anc) != mm.coloc_best_susie_ancestry).sum())

# ---------- R07-05 IDs without record ----------
catids = set(g.ensembl_id)
q = pd.read_csv(CAT / "passport_identity_quarantine.tsv", sep="\t", dtype=str)
R["quarantine_reasons"] = vc(q.reason_code)
qmap = {}
for _, rr in q.iterrows():
    for col in ["normalized_ensembl_candidates", "source_ensembl_text"]:
        for x in str(rr[col]).split(";"):
            x = x.strip().split(".")[0]
            if x.startswith("ENSG"):
                qmap.setdefault(x, set()).add(rr.reason_code)
absb = c[~c.ens.isin(catids)].copy()
absb["reason"] = absb.ens.map(lambda x: ";".join(sorted(qmap.get(x, []))) or "none")
R["bulk_absent"] = len(absb)
R["bulk_absent_reason"] = vc(absb.reason)
R["bulk_absent_q05"] = int((absb.padj < 0.05).sum())
R["bulk_absent_canon"] = absb[absb.canon][["ens", "symbol", "logFC", "padj"]].to_dict("records")
symcat = g.groupby("symbol").ensembl_id.apply(list)
nq = absb[absb.reason == "none"]
R["bulk_absent_unquarantined"] = [dict(ens=r.ens, sym=r.symbol, cat_ids_same_symbol=symcat.get(r.symbol, [])) for r in nq.itertuples()]
nsym_bulk = c.groupby("symbol").ens.nunique()
R["bulk_absent_symbol_shared_other_bulk_id"] = int(absb.symbol.map(nsym_bulk).gt(1).sum())
allids_sym = pd.concat([c[["ens", "symbol"]], g.rename(columns={"ensembl_id": "ens"})[["ens", "symbol"]]]).drop_duplicates()
nsym_any = allids_sym.groupby("symbol").ens.nunique()
R["bulk_absent_symbol_shared_any"] = int(absb.symbol.map(nsym_any).gt(1).sum())
absb_t = t12[~t12.ens.isin(catids)].copy()
absb_t["reason"] = absb_t.ens.map(lambda x: ";".join(sorted(qmap.get(x, []))) or "none")
R["t12_ids"] = int(t12.ens.nunique())
R["t12_absent"] = len(absb_t)
R["t12_absent_reason"] = vc(absb_t.reason)
fj = pd.read_csv(COMP / "identity_reconciliation_v3/source_symbol_false_joint_links.tsv", sep="\t", dtype=str)
six = fj[fj.catalog_id.notna() & (fj.catalog_id != "")]
R["six_rows"] = six[["symbol", "ensembl_bulk", "ensembl_genetic", "catalog_id", "catalog_class"]].to_dict("records")
six_gen = set(six.ensembl_genetic)
unq = absb_t[absb_t.reason == "none"]
R["t12_absent_unquarantined"] = [dict(ens=x, in_six=x in six_gen, glc_symbol=str(glc.set_index("ens").gene.get(x)) if "gene" in glc.columns and x in set(glc.ens) else None)
                                 for x in unq.ens]
R["ENSG00000280739_in_catalog"] = dict(present="ENSG00000280739" in catids,
                                       symbol=g.set_index("ensembl_id").symbol.get("ENSG00000280739"))
R["ENSG00000223797_allg_symbols"] = sorted(allg.gene[allg.ens == "ENSG00000223797"].astype(str).unique().tolist())
R["ENSG00000223797_in_bulk"] = bool("ENSG00000223797" in set(c.ens))
pm = pd.read_parquet(CAT / "passport_program_membership.parquet")
pm["ens"] = b(pm.ensembl_id)
R["pm_rows"] = len(pm)
R["pm_programs"] = int(pm.program_uid.nunique())
R["pm_ids"] = int(pm.ens.nunique())
R["pm_ids_in_cat"] = int(pm.ens[pm.ens.isin(catids)].nunique())
R["pm_programs_with_cat_member"] = int(pm.program_uid[pm.ens.isin(catids)].nunique())
pa = pm[~pm.ens.isin(catids)]
R["pm_ids_absent"] = int(pa.ens.nunique())
R["pm_rows_absent"] = len(pa)
anycoloc = set(allg.ens) | set(t12.ens) | set(glc.ens)
R["pm_absent_in_bulk"] = int(pa.ens.drop_duplicates().isin(c.ens).sum())
R["pm_absent_in_coloc"] = int(pa.ens.drop_duplicates().isin(anycoloc).sum())
R["pm_absent_quarantined"] = int(pa.ens.drop_duplicates().map(lambda x: x in qmap).sum())

# GENCODE v49 gene IDs
v49 = set()
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[2] != "gene":
            continue
        m = re.search(r'gene_id "([^"]+)"', f[8])
        v49.add(m.group(1).split(".")[0])
R["v49_genes"] = len(v49)
R["six_bulk_in_v49"] = {r.symbol: r.ensembl_bulk in v49 for r in six.itertuples()}
R["six_genetic_in_v49"] = {r.symbol: r.ensembl_genetic in v49 for r in six.itertuples()}
R["six_bulk_in_any_coloc"] = {r.symbol: r.ensembl_bulk in anycoloc for r in six.itertuples()}
R["six_catalog_gen_est_eq_genetic_abf"] = {
    r.symbol: [None if pd.isna(gen.estimate.get(r.catalog_id)) else float(gen.estimate.get(r.catalog_id)),
               None if r.ensembl_genetic not in T.index else float(T.coloc_best_abf_pp4.get(r.ensembl_genetic))]
    for r in six.itertuples()}
R["six_catalog_class"] = {r.symbol: cls.get(r.catalog_id) for r in six.itertuples()}
R["t12_absent_v49_absent"] = int((~absb_t.ens.isin(v49)).sum())
R["t12_absent_unquarantined_v49"] = {x: (x in v49) for x in unq.ens}

# ---------- R15-05 assay applicability ----------
mo = pd.DataFrame([json.loads(line) for line in open(COMP / "molecular_objects/molecular_object_sidecar.jsonl")])
R["mo_rows"] = len(mo)
aa = mo.assay_applicability.apply(lambda x: json.loads(x) if isinstance(x, str) else x).apply(pd.Series)
mo2 = pd.concat([mo[["ensembl_id_base", "symbol", "rna_call"]], aa], axis=1)
R["mo_bulk_rna_by_rna_call"] = {f"{a}|{b_}": int(n) for (a, b_), n in mo2.groupby(["rna_call", "bulk_rna"]).size().items()}
for col in ["proteomics", "chromatin", "spatial"]:
    R[f"mo_{col}"] = vc(mo2[col])
R["mo_candidate_object"] = vc(mo.candidate_object)
R["mo_candidate_not_adopted"] = vc(mo.candidate_not_adopted.astype(str))
R["mo_keys"] = sorted(mo.columns.tolist())
h = pd.read_csv(H86, sep="\t")
R["h86_route"] = vc(h.route)
gsym = set(g.symbol)
R["h86_direct_in_catalog_symbols"] = int(h[h.route == "direct"].catalog_symbol.isin(gsym).sum())
R["h86_remapped_in_catalog_symbols"] = int(h[h.route == "remapped"].catalog_symbol.isin(gsym).sum())
qset = set(h[h.route.isin(["direct", "remapped"])].catalog_symbol)
R["h86_quantified_catalog_genes"] = int(g.symbol.isin(qset).sum())
R["h86_quantified_unique_symbols"] = len(qset)
R["h86_by_class"] = vc(g[g.symbol.isin(qset)].primary_evidence_class)

# ---------- program context ----------
pc = pd.read_parquet(CAT / "passport_program_context.parquet")
R["pc_rows"] = len(pc)
R["pc_programs"] = sorted(pc.program_label.unique().tolist())
R["pc_datasets"] = sorted(pc.dataset_id.unique().tolist())
R["pc_source_release"] = sorted(pc.source_release_id.unique().tolist())
R["pc_table"] = pc[["program_label", "dataset_id", "assay", "estimate", "call_state", "testability_state", "source_dependent"]].astype(str).to_dict("records")
pi = pd.read_parquet(CAT / "passport_program_index.parquet")
R["program_index_rows"] = len(pi)
ne = pd.read_csv(CAT / "passport_next_experiment.tsv", sep="\t", usecols=["ensembl_id", "experiment_rule_id"])
R["next_experiment_rows"] = len(ne)
R["next_experiment_unique_ids"] = int(ne.ensembl_id.nunique())
R["gene_biotype_mo"] = vc(mo.gene_biotype) if "gene_biotype" in mo.columns else None

# ---------- observability ----------
cov = pd.read_csv(SPAT / "spatial_program_coverage.tsv", sep="\t")
R["spatial_cov"] = {f"{d}|{s}": int(n) for (d, s), n in cov.groupby(["dataset_id", "coverage_status"]).size().items()}
ga = pd.read_csv(SPAT / "gene_axis_audit.tsv", sep="\t")
R["gene_axis"] = ga[["dataset_id", "gene_universe_denominator", "n_detected_genes"]].astype(str).to_dict("records")
da = pd.read_csv(ATAC / "da_peak_results.tsv.gz", sep="\t")
R["atac_rows"] = len(da)
R["atac_states"] = vc(da.evidence_state)
sd = da[da.evidence_state == "source_dependent"]
q1 = sd.qvalue_gse244832 < 0.05
q2 = sd.qvalue_gse281367 < 0.05
R["atac_sd_sig_one_cohort_only"] = int((q1 ^ q2).sum())
R["atac_sd_sig_both"] = int((q1 & q2).sum())
R["atac_sd_sig_none"] = int((~q1 & ~q2).sum())
R["atac_sd_same_sign"] = int((np.sign(sd.logFC_gse244832) == np.sign(sd.logFC_gse281367)).sum())
ls = pd.read_csv(ATAC / "da_lineage_summary.tsv", sep="\t")
R["atac_donors"] = ls[["lineage", "n_normal_gse244832", "n_mash_gse244832", "n_normal_gse281367", "n_mash_gse281367", "n_jointly_testable"]].astype(str).to_dict("records")

(OUT / "recheck_core.json").write_text(json.dumps(R, indent=1, default=str))
print("core done")

# ---------- Myojin refit with frozen code (copied, hash-checked) ----------
if RUN_PERM:
    import hashlib
    cp = OUT / "phase_c_stats_frozen_copy.py"
    shutil.copy(FROZEN, cp)
    P = {"frozen_sha256": hashlib.sha256(cp.read_bytes()).hexdigest()}
    spec = importlib.util.spec_from_file_location("pcs", cp)
    pcs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pcs)
    sr = pd.read_csv(MY / "screen_results.tsv", sep="\t", dtype=str, keep_default_na=False)
    sr = sr[sr.primary_eligible.str.strip().str.upper() == "TRUE"]
    fr = pd.DataFrame({
        "source_row": sr.source_row.astype(int).to_numpy(),
        "gene_symbol": sr.gene_symbol.to_numpy(),
        "primary_evidence_class": sr.primary_evidence_class.to_numpy(),
        "guide_count_min": sr.guide_count_min.astype(int).to_numpy(),
        "log1p_HLF_TPM": sr.log1p_HLF_TPM.astype(float).to_numpy(),
        "HLF_Chronos": pd.to_numeric(sr.HLF_Chronos.replace("", np.nan)).to_numpy(),
        "HLF_Chronos_missing": (sr.HLF_Chronos == "").astype(int).to_numpy(),
        "PA_vs_vehicle_LFC": sr.PA_vs_vehicle_positive_selection_LFC.astype(float).to_numpy(),
    }).sort_values("gene_symbol", kind="stable").reset_index(drop=True)
    P["n_primary"] = len(fr)
    P["chronos_missing"] = int(fr.HLF_Chronos_missing.sum())
    prim = pcs.add_frozen_strata(fr)
    gu = pd.read_csv(MY / "screen_gene_universe.tsv", sep="\t", dtype=str)
    P["gene_universe_cols"] = gu.columns.tolist()
    ensmap = gu.assign(source_row=gu.source_row.astype(int)).set_index("source_row").gencode_ensembl_base.str.split(".").str[0]
    prim["r8"] = prim.source_row.map(ensmap).map(cls).fillna("not_in_r8")
    P["r8_counts"] = vc(prim.r8)
    m = {"unresolved": "neither", "genetically_anchored": "genetic_only", "established_state_associated": "disease_state_only",
         "concordant": "convergent", "untested": "neither"}
    A = prim[~prim.r8.isin(["untested", "not_in_r8"])].copy()
    A["primary_evidence_class"] = A.r8.map(m)
    A = A.drop(columns=["r8"]).reset_index(drop=True)
    Bv = prim[prim.r8 != "not_in_r8"].copy()
    Bv["primary_evidence_class"] = Bv.r8.map(m)
    Bv = Bv.drop(columns=["r8"]).reset_index(drop=True)
    # independent observed F via statsmodels formula (no frozen code)
    for nm, df in [("A", A), ("B", Bv)]:
        X = pd.get_dummies(df.primary_evidence_class)[["genetic_only", "disease_state_only", "convergent"]].astype(float)
        Z = df[["log1p_HLF_TPM", "guide_count_min", "HLF_Chronos"]].astype(float)
        full = sm.add_constant(pd.concat([X, Z], axis=1))
        red = sm.add_constant(Z)
        f1 = sm.OLS(df.PA_vs_vehicle_LFC, full).fit()
        f0 = sm.OLS(df.PA_vs_vehicle_LFC, red).fit()
        F = ((f0.ssr - f1.ssr) / 3) / (f1.ssr / f1.df_resid)
        h3 = sm.OLS(df.PA_vs_vehicle_LFC, full).fit(cov_type="HC3")
        sdr = np.sqrt(f1.ssr / f1.df_resid)
        P[f"indep_{nm}"] = dict(n=len(df), F=F, F_param_p=float(stats.f.sf(F, 3, f1.df_resid)),
                                n_by_class=vc(df.primary_evidence_class), n_strata=int(df.permutation_stratum.nunique()),
                                std=(f1.params[["genetic_only", "disease_state_only", "convergent"]] / sdr).round(4).to_dict(),
                                conv_ci=h3.conf_int().loc["convergent"].round(4).tolist())
    ra, aa_ = pcs.fit_class_continuous(A, "PA_vs_vehicle_LFC", "r8_excl_untested", 2026080740, accepted_draws=100_000)
    rb, ab_ = pcs.fit_class_continuous(Bv, "PA_vs_vehicle_LFC", "r8_untested_in_ref", 2026080740, accepted_draws=100_000)
    pd.DataFrame(ra + rb).to_csv(OUT / "myojin_frozen_refit_recheck.tsv", sep="\t", index=False)
    pd.DataFrame([aa_, ab_]).to_csv(OUT / "myojin_permutation_audit_recheck.tsv", sep="\t", index=False)
    P["frozen_A"] = [{k: r[k] for k in ["contrast", "n_class", "n_reference", "statistic", "permutation_p", "standardized_effect", "CI95_low", "CI95_high"]} for r in ra]
    P["frozen_B"] = [{k: r[k] for k in ["contrast", "n_class", "n_reference", "statistic", "permutation_p", "standardized_effect", "CI95_low", "CI95_high"]} for r in rb]
    P["audit"] = [aa_, ab_]
    (OUT / "recheck_myojin.json").write_text(json.dumps(P, indent=1, default=str))
    print("perm done")
