#!/usr/bin/env python3
"""Blind check of F6 proposed Gene Catalog corrections.

Independent re-derivation from release files (r8 reviewed Catalog, adopted
774-unit bulk fit, COLOC r2, PXD051911 protein DE, molecular-object companion,
GENCODE v49, Myojin HLF screen). Compares with the proposer's candidate tables.
Writes one JSON of numbers plus a few TSVs to OUT.
"""
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
CAT = ROOT / "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8"
CAT_ORIG = ROOT / "RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z"
COMP = ROOT / "figures/candidates/catalog-r8-completion-20261014T141645Z_22236527"
BULK = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
COLOC = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
TIERCFG = ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv"
PROT = ROOT / "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv"
MY = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf"
GTF = Path("/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz")
CAND = ROOT / "figures/candidates/review-r3-fixes-20261015T112148Z/F6/panels_20261015T114537Z"
OUT = Path(sys.argv[1])
OUT.mkdir(parents=False, exist_ok=False)

R = {}


def base(s):
    return s.astype(str).str.replace(r"\..*$", "", regex=True)


def vc(s):
    return {str(k): int(v) for k, v in s.value_counts(dropna=False).items()}


# ---------------- Catalog
e = pd.read_parquet(CAT / "passport_evidence_long.parquet")
g = pd.read_parquet(CAT / "passport_gene_index.parquet")
pm = pd.read_parquet(CAT / "passport_program_membership.parquet")
# reviewed copy equals the original build?
eo = pd.read_parquet(CAT_ORIG / "passport_evidence_long.parquet")
go = pd.read_parquet(CAT_ORIG / "passport_gene_index.parquet")
R["reviewed_equals_original_evidence"] = bool(e.equals(eo))
R["reviewed_equals_original_gene_index"] = bool(g.equals(go))

R["n_evidence_rows"] = len(e)
R["n_genes"] = int(g["ensembl_id"].nunique())
R["domains"] = vc(e["evidence_domain"])
R["assays"] = vc(e["assay"])
for c in ["negative_call_rule_id", "negative_decision_boundary", "negative_margin", "negative_call_passed",
          "source_negative_call_rule_id", "source_negative_call_passed",
          "standard_error", "ci_lower", "ci_upper", "p_value", "n_biological_units", "n_technical_units"]:
    R[f"nonnull_{c}"] = int(e[c].notna().sum())
R["call_state_by_domain"] = {f"{a}|{b}": int(n) for (a, b), n in e.groupby(["evidence_domain", "call_state"]).size().items()}
R["testability_by_domain"] = {f"{a}|{b}": int(n) for (a, b), n in e.groupby(["evidence_domain", "testability_state"]).size().items()}
R["source_dependent"] = vc(e["source_dependent"].astype(str))
R["provenance_state"] = vc(e["provenance_state"])
# any column containing the literal not_applicable
na_cols = {c: int((e[c].astype(str) == "not_applicable").sum()) for c in e.columns if e[c].dtype == object}
R["not_applicable_cells_by_column"] = {k: v for k, v in na_cols.items() if v}
R["direction_by_domain"] = {f"{a}|{b}": int(n) for (a, b), n in e.groupby(["evidence_domain", e["direction"].astype(str)]).size().items()}
R["biological_unit_values"] = vc(e["biological_unit"])
R["primary_class"] = vc(g["primary_evidence_class"])
R["role_hypothesis"] = vc(g["role_hypothesis"])

e["ens"] = base(e["ensembl_id"])
g["ens"] = base(g["ensembl_id"])
gen = e[e["evidence_domain"] == "genetics"].set_index("ens")
rna = e[e["evidence_domain"] == "transcriptomics"].set_index("ens")
assert gen.index.is_unique and rna.index.is_unique
gi = g.set_index("ens")
M = pd.DataFrame({"gen": gen["call_state"], "rna": rna["call_state"], "cls": gi["primary_evidence_class"],
                  "gen_test": gen["testability_state"], "rna_test": rna["testability_state"],
                  "gen_est": gen["estimate"], "rna_est": rna["estimate"], "rna_q": rna["q_value"],
                  "sym": gi["symbol"]})
R["matrix_gen_x_rna"] = {f"{a}|{b}": int(n) for (a, b), n in M.groupby(["gen", "rna"]).size().items()}
R["matrix_cell_class"] = {f"{a}|{b}|{c}": int(n) for (a, b, c), n in M.groupby(["gen", "rna", "cls"]).size().items()}
R["untested_by_test_states"] = {f"{a}|{b}": int(n) for (a, b), n in
                                M[M.cls == "untested"].groupby(["gen_test", "rna_test"]).size().items()}

# ---------------- Bulk (adopted 774-unit fit)
ct = pd.read_csv(BULK / "coefficient_table.tsv.gz", sep="\t")
dr = pd.read_csv(BULK / "deg_results.csv")
ct["ens"] = base(ct["gene"])
dr["ens"] = base(dr["gene"])
R["bulk_rows_coef"] = len(ct)
R["bulk_rows_deg"] = len(dr)
mm = ct.merge(dr[["ens", "logFC", "padj", "SE", "P.Value"]], on="ens", suffixes=("", "_dr"))
R["bulk_coef_vs_deg_max_abs_logFC_diff"] = float((mm.logFC - mm.logFC_dr).abs().max())
R["bulk_coef_vs_deg_max_abs_padj_diff"] = float((mm.padj - mm.padj_dr).abs().max())
R["bulk_coef_vs_deg_max_rel_padj_diff"] = float(((mm.padj - mm.padj_dr).abs() / mm.padj_dr.clip(lower=1e-300)).max())
dr["canon"] = (dr.padj < 0.05) & (dr.logFC.abs() > 0.5)
R["bulk_canonical"] = int(dr.canon.sum())
R["bulk_canonical_up"] = int((dr.canon & (dr.logFC > 0)).sum())
R["bulk_q_lt_0.05_abs_le_0.5"] = int(((dr.padj < 0.05) & (dr.logFC.abs() <= 0.5)).sum())
md = pd.read_csv(BULK / "model_design.tsv", sep="\t")
R["design_n"] = len(md)
R["design_disease"] = int(md["group_binaryDisease"].sum())
R["design_control"] = int((md["group_binaryDisease"] == 0).sum())

d = dr.set_index("ens")
M = M.join(d[["logFC", "SE", "padj", "P.Value", "canon"]], how="left")
# Catalog RNA estimate equals bulk logFC?
both = M[M.logFC.notna() & M.rna_est.notna()]
R["catalog_rna_est_vs_bulk_maxdiff"] = float((both.rna_est - both.logFC).abs().max())
R["catalog_rna_q_vs_bulk_maxreldiff"] = float(((both.rna_q - both.padj).abs() / both.padj.clip(lower=1e-300)).max())
R["catalog_rna_supported_equals_canon"] = {f"{a}|{b}": int(n) for (a, b), n in
                                          M.groupby([M.rna, M.canon.astype(str)]).size().items()}
ind = M[M.rna == "indeterminate"]
R["rna_indet"] = len(ind)
R["rna_indet_q_lt_0.05"] = int((ind.rna_q < 0.05).sum())
R["rna_indet_q_lt_0.05_bulkpadj"] = int((ind.padj < 0.05).sum())
R["rna_indet_q_lt_0.05_max_abs_lfc"] = float(ind.loc[ind.rna_q < 0.05, "rna_est"].abs().max())
R["rna_untestable_has_bulk_row"] = int(M[(M.rna == "untestable")].logFC.notna().sum())
R["unresolved_q_lt_0.05"] = int(((M.cls == "unresolved") & (M.rna_q < 0.05)).sum())
R["anchored_q_lt_0.05"] = int(((M.cls == "genetically_anchored") & (M.rna_q < 0.05)).sum())

# bulk IDs absent from the Catalog
cat_ids = set(g["ens"])
absent = dr[~dr.ens.isin(cat_ids)].copy()
R["bulk_ids_absent_from_catalog"] = len(absent)
R["absent_padj_lt_0.05"] = int((absent.padj < 0.05).sum())
R["absent_canonical"] = absent.loc[absent.canon, ["ens", "symbol", "logFC", "padj"]].to_dict("records")
# symbol sharing, three definitions
sym_bulk_n = dr.groupby("symbol")["ens"].nunique()
R["absent_symbol_shared_in_bulk"] = int(absent.symbol.map(sym_bulk_n).gt(1).sum())
# shared with any ID in bulk + Catalog + COLOC t12
t12 = pd.read_csv(COLOC / "gene_level_coloc_tier12.csv")
allc = pd.read_csv(COLOC / "gene_level_coloc.csv", usecols=["gene", "ensembl"])
sym_map = pd.concat([dr[["symbol", "ens"]], g[["symbol", "ens"]],
                     allc.rename(columns={"gene": "symbol", "ensembl": "ens"})]).dropna().drop_duplicates()
sym_any_n = sym_map.groupby("symbol")["ens"].nunique()
R["absent_symbol_shared_any_source"] = int(absent.symbol.map(sym_any_n).gt(1).sum())
# GEN source quarantine table
q = pd.read_csv(CAT / "passport_identity_quarantine.tsv", sep="\t", dtype=str)
R["quarantine_rows"] = len(q)
R["quarantine_columns"] = list(q.columns)
qsyms = set(q[[c for c in q.columns if "symbol" in c.lower()][0]]) if any("symbol" in c.lower() for c in q.columns) else set()
R["absent_symbol_in_quarantine"] = int(absent.symbol.isin(qsyms).sum())
q_ids = {}
for _, rq in q.iterrows():
    for x in str(rq.normalized_ensembl_candidates).split(";"):
        q_ids.setdefault(x.strip().split(".")[0], set()).add(rq.reason_code)
    for x in str(rq.source_ensembl_text).split(";"):
        q_ids.setdefault(x.strip().split(".")[0], set()).add(rq.reason_code)
R["quarantine_reason_codes"] = vc(q.reason_code)
absent["in_quarantine_ids"] = absent.ens.isin(set(q_ids))
absent["quarantine_reason"] = absent.ens.map(lambda x: ";".join(sorted(q_ids.get(x, []))))
R["absent_in_quarantine_ids"] = int(absent.in_quarantine_ids.sum())
R["absent_quarantine_reason"] = vc(absent.quarantine_reason)
t12_tmp = pd.read_csv(COLOC / "gene_level_coloc_tier12.csv", usecols=["ensembl"])
t12a = t12_tmp[~t12_tmp.ensembl.astype(str).isin(cat_ids)].copy()
t12a["reason"] = t12a.ensembl.astype(str).map(lambda x: ";".join(sorted(q_ids.get(x, []))))
R["t12_absent_quarantine_reason"] = vc(t12a.reason)
absent["sym_n_bulk"] = absent.symbol.map(sym_bulk_n)
absent["sym_n_any"] = absent.symbol.map(sym_any_n)
absent["in_quarantine_symbols"] = absent.symbol.isin(qsyms)
absent.sort_values("padj").to_csv(OUT / "bulk_ids_absent_from_catalog_check.tsv", sep="\t", index=False)
R["absent_not_shared_not_quarantined"] = absent.loc[(absent.sym_n_any <= 1) & ~absent.in_quarantine_symbols,
                                                    ["ens", "symbol", "padj"]].to_dict("records")
R["absent_not_in_quarantine"] = absent.loc[~absent.in_quarantine_symbols, ["ens", "symbol", "padj", "sym_n_any"]].to_dict("records")
R["catalog_rna_supported"] = int((M.rna == "supported").sum())

# COLOC Tier-1/2 IDs absent
t12["ens"] = t12["ensembl"].astype(str)
R["t12_ids"] = int(t12.ens.nunique())
t12_abs = t12[~t12.ens.isin(cat_ids)]
R["t12_ids_absent"] = len(t12_abs)
# program members absent
pm["ens"] = base(pm["ensembl_id"])
pma = pm[~pm.ens.isin(cat_ids)]
R["program_member_ids_absent"] = int(pma.ens.nunique())
R["program_memberships_absent"] = len(pma)
R["program_memberships_total"] = len(pm)

# ---------------- COLOC: what the genetics estimate shows
T = t12.set_index("ens")
G = M.join(T[["coloc_best_susie_pp4", "coloc_best_abf_pp4", "susie_driving_gwas", "abf_driving_gwas",
              "susie_state_t12", "susie_untestable_t12_studies"]], how="left")
s, a, est = G.coloc_best_susie_pp4, G.coloc_best_abf_pp4, G.gen_est
G["meth"] = np.where(est.isna(), "none",
                     np.where(np.isclose(est, s, rtol=0, atol=1e-12) & s.notna(), "SuSiE",
                              np.where(np.isclose(est, a, rtol=0, atol=1e-12) & a.notna(), "ABF", "other")))
R["pp4_method_by_call"] = {f"{x}|{y}": int(n) for (x, y), n in G.groupby(["gen", "meth"]).size().items()}
oth = G[G.meth == "other"]
R["pp4_other_rows"] = oth[["sym", "gen", "gen_est", "coloc_best_susie_pp4", "coloc_best_abf_pp4"]].reset_index().to_dict("records")
ind_g = G[G.gen == "indeterminate"]
R["gen_indet_est_gt_0.5"] = int((ind_g.gen_est > 0.5).sum())
R["gen_indet_est_ge_0.8"] = int((ind_g.gen_est >= 0.8).sum())
R["gen_indet_est_gt_0.5_finite_susie_le_0.5"] = int(((ind_g.gen_est > 0.5) & ind_g.coloc_best_susie_pp4.notna()
                                                      & (ind_g.coloc_best_susie_pp4 <= 0.5)).sum())
R["gen_indet_est_gt_0.5_meth"] = vc(ind_g.loc[ind_g.gen_est > 0.5, "meth"])
R["gen_supported_susie_gt_0.5"] = int(((G.gen == "supported") & (G.coloc_best_susie_pp4 > 0.5)).sum())
R["gen_nonsupported_susie_gt_0.5"] = int(((G.gen != "supported") & (G.coloc_best_susie_pp4 > 0.5)).sum())
R["t12_susie_gt_0.5"] = int((t12.coloc_best_susie_pp4 > 0.5).sum())
R["t12_susie_state"] = vc(t12.susie_state_t12)
# genetics untestable = no t12 record / SuSiE untestable
gu = G[G.gen == "untestable"]
R["gen_untestable_reason"] = vc(gu.gen_test)
R["gen_underpowered_in_t12"] = int(gu[gu.gen_test == "underpowered_source"].index.isin(t12.ens).sum())
allc_ids = set(allc.ensembl.astype(str))
R["gen_underpowered_in_any_gene_level_coloc"] = int(gu[gu.gen_test == "underpowered_source"].index.isin(allc_ids).sum())
R["gen_insufficient_susie_state"] = vc(gu.loc[gu.gen_test == "insufficient_shared_posterior", "susie_state_t12"])
R["gen_insufficient_genes"] = sorted(gu.loc[gu.gen_test == "insufficient_shared_posterior", "sym"].tolist())
for gene in ["GNMT", "DGAT2", "CFLAR", "EZR", "THRB", "HKDC1", "MTARC1", "GLP1R"]:
    r = G[G.sym == gene]
    R[f"worked_{gene}"] = r[["gen", "rna", "cls", "gen_est", "coloc_best_susie_pp4", "coloc_best_abf_pp4",
                             "susie_driving_gwas", "abf_driving_gwas", "susie_state_t12", "rna_est", "rna_q",
                             "logFC", "SE", "padj", "P.Value", "gen_test", "rna_test"]].reset_index().to_dict("records")

# ---------------- PP.H3 among eligible SuSiE signal pairs
tier = pd.read_csv(TIERCFG, sep="\t").set_index("study_name")["tier"]
sp = pd.read_csv(COLOC / "susie_coloc_signal_pairs.tsv.gz", sep="\t")
sp["tier"] = sp.gwas.map(tier)
R["sp_rows"] = len(sp)
R["sp_tier_counts"] = vc(sp.tier)
R["sp_eligible_values"] = vc(sp.eligible.astype(str))
el = sp[sp.eligible.astype(str).str.upper() == "TRUE"].copy()
R["sp_eligible_rows"] = len(el)


def h3_summary(df, label):
    agg = df.groupby("ensembl").agg(max_h3=("PP.H3", "max"), max_h4=("PP.H4", "max"), sym=("gene", "first"))
    hit = agg[(agg.max_h3 > 0.8) & (agg.max_h4 <= 0.5)].copy()
    hit["catalog_gen_call"] = hit.index.map(lambda x: G.gen.get(x, "absent_from_catalog"))
    strict = agg[(agg.max_h3 > 0.9) & (agg.max_h4 < 0.1)]
    R[f"{label}_genes"] = len(agg)
    R[f"{label}_h3gt0.8_h4le0.5"] = len(hit)
    R[f"{label}_h3_by_call"] = vc(hit.catalog_gen_call)
    R[f"{label}_strict_h3gt0.9_h4lt0.1"] = len(strict)
    return hit


hit_all = h3_summary(el, "sp_elig_all")
hit_t12 = h3_summary(el[el.tier.isin([1, 2])], "sp_elig_t12")
hit_t12.to_csv(OUT / "pph3_t12_genes.tsv", sep="\t")
hit_all.to_csv(OUT / "pph3_all_genes.tsv", sep="\t")
# per-pair definition: some eligible pair with H3>0.8 and the gene's max H4<=0.5 is what above uses.

# ---------------- concordant direction
conc = M[M.cls == "concordant"]
R["concordant_n"] = len(conc)
R["concordant_rna_up"] = int((conc.rna_est > 0).sum())
R["concordant_cell"] = vc(conc.gen + "|" + conc.rna)
# class mapping to cells
R["class_cells"] = {c: vc(M.loc[M.cls == c, "gen"] + "|" + M.loc[M.cls == c, "rna"]) for c in M.cls.unique()}

# ---------------- R07-11 six genes
six = ["CTSLP2", "ITFG2-AS1", "MATR3", "PAXX", "POLR2J4", "SUGT1P3"]
fj = pd.read_csv(COMP / "identity_reconciliation_v3/source_symbol_false_joint_links.tsv", sep="\t", dtype=str)
fj = fj[fj.symbol.isin(six)]
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
R["gencode_v49_gene_ids"] = len(v49)
allgw = pd.read_csv(COLOC / "susie_coloc_all_gwas.csv",
                    usecols=["gwas_name", "gene", "ensembl", "PP.H4.abf", "PP.H4.susie", "ancestry", "method"])
allgw_ids = set(allgw.ensembl.astype(str))
rows = []
for _, r in fj.iterrows():
    b, gid = r.ensembl_bulk, r.ensembl_genetic
    cat_gen_est = G.gen_est.get(b, np.nan)
    rec = T.loc[gid] if gid in T.index else None
    rows.append(dict(symbol=r.symbol, bulk=b, genetic=gid, catalog_has_bulk=b in cat_ids,
                     catalog_has_genetic=gid in cat_ids, bulk_in_v49=b in v49, genetic_in_v49=gid in v49,
                     bulk_in_t12=b in set(t12.ens), bulk_in_gene_level=b in allc_ids, bulk_in_all_gwas=b in allgw_ids,
                     genetic_in_t12=gid in set(t12.ens), catalog_gen_call=G.gen.get(b), catalog_cls=G.cls.get(b),
                     catalog_rna_call=G.rna.get(b), catalog_gen_est=cat_gen_est,
                     t12_genetic_susie=None if rec is None else rec.coloc_best_susie_pp4,
                     t12_genetic_abf=None if rec is None else rec.coloc_best_abf_pp4,
                     bulk_tested=b in set(dr.ens)))
six_df = pd.DataFrame(rows)
six_df.to_csv(OUT / "R07_11_six_check.tsv", sep="\t", index=False)
R["R07_11_six"] = rows
R["t12_absent_includes_six_genetic"] = int(t12_abs.ens.isin(fj.ensembl_genetic).sum())

# ---------------- worked-entry ABF studies with ancestry from master
for gene, ens in [("GLP1R", None), ("MTARC1", None), ("HKDC1", None), ("THRB", None), ("GNMT", None)]:
    sub = allgw[(allgw.gene == gene) & allgw.gwas_name.map(tier).isin([1, 2])]
    if len(sub):
        ra = sub.loc[sub["PP.H4.abf"].idxmax()]
        R[f"master_{gene}_best_abf_t12"] = dict(gwas=ra.gwas_name, pp4=float(ra["PP.H4.abf"]), anc=ra.ancestry)
        ss = sub[sub["PP.H4.susie"].notna()]
        if len(ss):
            rs = ss.loc[ss["PP.H4.susie"].idxmax()]
            R[f"master_{gene}_best_susie_t12"] = dict(gwas=rs.gwas_name, pp4=float(rs["PP.H4.susie"]), anc=rs.ancestry)
        R[f"master_{gene}_methods"] = vc(sub.method)
del allgw

# ---------------- proposer's sidecar vs release values
sc = pd.read_csv(CAND / "R15_03_gene_effect_sidecar.tsv.gz", sep="\t")
R["sidecar_columns"] = list(sc.columns)
R["sidecar_rows"] = len(sc)
sc["ens"] = base(sc[sc.columns[0]])
scj = sc.set_index("ens").join(d[["logFC", "SE", "t", "P.Value", "padj"]], rsuffix="_rel")
for col in ["logFC", "SE", "t", "P.Value", "padj"]:
    cands = [c for c in sc.columns if c.lower().replace(".", "_") in
             {col.lower().replace(".", "_"), "bulk_" + col.lower().replace(".", "_")}]
    if cands and f"{col}_rel" in scj:
        x = scj[cands[0]].astype(float)
        y = scj[f"{col}_rel"].astype(float)
        ok = x.notna() & y.notna()
        R[f"sidecar_{col}_col"] = cands[0]
        R[f"sidecar_{col}_maxreldiff"] = float(((x[ok] - y[ok]).abs() / y[ok].abs().clip(lower=1e-300)).max())
        R[f"sidecar_{col}_n_nonnull"] = int(x.notna().sum())
scj2 = sc.set_index("ens").join(T[["coloc_best_susie_pp4", "coloc_best_abf_pp4", "susie_driving_gwas",
                                   "abf_driving_gwas"]], rsuffix="_rel")
for col in ["coloc_best_susie_pp4", "coloc_best_abf_pp4"]:
    cands = [c for c in sc.columns if col.split("_")[2] in c.lower() and "pp4" in c.lower()]
    R[f"sidecar_cands_{col}"] = cands
for c in sc.columns:
    if sc[c].dtype != object and sc[c].nunique() < 5:
        R[f"sidecar_const_{c}"] = vc(sc[c].astype(str))
R["sidecar_head"] = sc.head(3).astype(str).to_dict("records")

# ---------------- molecular-object companion assay_applicability
rows = []
with open(COMP / "molecular_objects/molecular_object_sidecar.jsonl") as fh:
    for line in fh:
        j = json.loads(line)
        aa = j["assay_applicability"]
        aa = json.loads(aa) if isinstance(aa, str) else aa
        rows.append(dict(ens=j["ensembl_id_base"], sym=j["symbol"], rna_call=j["rna_call"], gen_call=j["genetics_call"],
                         cls=j["primary_evidence_class"], **{f"aa_{k}": v for k, v in aa.items()}))
mo = pd.DataFrame(rows)
R["mo_records"] = len(mo)
R["mo_unique_genes"] = int(mo.ens.nunique())
R["mo_rna_x_bulk"] = {f"{a}|{b}": int(n) for (a, b), n in mo.groupby(["rna_call", "aa_bulk_rna"]).size().items()}
for k in ["aa_proteomics", "aa_chromatin", "aa_spatial"]:
    R[f"mo_{k}"] = vc(mo[k])
R["mo_GNMT"] = mo[mo.sym == "GNMT"].to_dict("records")

# ---------------- PXD051911 matching
pr = pd.read_csv(PROT, sep="\t")
R["prot_rows"] = len(pr)
R["prot_q_lt_0.05"] = int((pr.padj < 0.05).sum())
R["prot_unique_gene_strings"] = int(pr.gene.nunique())
cat_sym_n = g.groupby("symbol")["ens"].nunique()
uniq_cat_syms = set(cat_sym_n[cat_sym_n == 1].index)
pr_sym_n = pr.groupby("gene").size()
pr_u = pr[pr.gene.map(pr_sym_n) == 1]
matched = pr_u[pr_u.gene.isin(uniq_cat_syms)]
R["prot_matched_unique_symbol"] = len(matched)
R["prot_matched_q_lt_0.05"] = int((matched.padj < 0.05).sum())
matched_any = pr[pr.gene.isin(set(g.symbol))]
R["prot_matched_any_symbol"] = int(matched_any.gene.nunique())
R["prot_matched_any_q_lt_0.05"] = int(matched_any.loc[matched_any.padj < 0.05, "gene"].nunique())
R["prot_GNMT"] = pr[pr.gene == "GNMT"][["gene", "logFC", "padj"]].to_dict("records")
corr = pd.read_csv(CAND / "R15_05_assay_applicability_corrected.tsv.gz", sep="\t")
R["cand_R15_05_columns"] = list(corr.columns)
for c in corr.columns:
    if corr[c].dtype == object and corr[c].nunique() < 12:
        R[f"cand_R15_05_vc_{c}"] = vc(corr[c])

# ---------------- worked-entry candidate source vs release
R["cand_fig6c_source"] = pd.read_csv(CAND / "fig6c_worked_entries_candidate_source.tsv", sep="\t").astype(str).to_dict("records")
R["cand_fig6b_source"] = pd.read_csv(CAND / "fig6b_call_matrix_3x3_source.tsv", sep="\t").astype(str).to_dict("records")
R["cand_fig6d_source"] = pd.read_csv(CAND / "fig6d_next_experiment_candidate_source.tsv", sep="\t").astype(str).to_dict("records")

# ---------------- Myojin refit with r8 classes (independent code)
sr = pd.read_csv(MY / "screen_results.tsv", sep="\t", low_memory=False)
gu_ = pd.read_csv(MY / "screen_gene_universe.tsv", sep="\t", low_memory=False)
R["myojin_sr_cols_subset"] = [c for c in sr.columns if any(k in c for k in
                              ["LFC", "eligible", "class", "TPM", "guide", "Chronos", "ensembl", "symbol"])]
key = "source_row"
sr = sr.merge(gu_[[key, "gencode_ensembl_base"]], on=key, how="left")
prim = sr[sr["primary_eligible"].astype(str).str.upper() == "TRUE"].copy()
R["myojin_primary_n"] = len(prim)
R["myojin_v1_classes"] = vc(prim["primary_evidence_class"])
prim["ens"] = base(prim["gencode_ensembl_base"])
prim["r8"] = prim.ens.map(gi["primary_evidence_class"]).fillna("not_in_r8")
R["myojin_r8_classes"] = vc(prim["r8"])
Y = "PA_vs_vehicle_positive_selection_LFC"
COV = ["log1p_HLF_TPM", "guide_count_min", "HLF_Chronos"]


def design(df, cls_col, levels):
    X = np.column_stack([np.ones(len(df))] + [(df[cls_col] == l).astype(float).values for l in levels]
                        + [df[c].astype(float).values for c in COV])
    return X


def ols_hc3(y, X):
    XtXi = np.linalg.inv(X.T @ X)
    b = XtXi @ X.T @ y
    r = y - X @ b
    h = np.einsum("ij,jk,ik->i", X, XtXi, X)
    w = (r / (1 - h)) ** 2
    V = XtXi @ (X.T * w) @ X @ XtXi
    rank = np.linalg.matrix_rank(X)
    rss = float(r @ r)
    return b, V, rss, rank


from scipy import stats


def class_fit(df, cls_col, levels, label, nperm=0, seed=2026080740):
    y = df[Y].astype(float).values
    X = design(df, cls_col, levels)
    b, V, rss, rank = ols_hc3(y, X)
    n = len(y)
    rsd = np.sqrt(rss / (n - rank))
    out = []
    idx = list(range(1, 1 + len(levels)))
    for i, l in zip(idx, levels):
        se = np.sqrt(V[i, i])
        z = b[i] / se
        out.append(dict(analysis=label, contrast=l, n=n, n_class=int((df[cls_col] == l).sum()), estimate=b[i], hc3_se=se,
                        ci_low=b[i] - 1.959964 * se, ci_high=b[i] + 1.959964 * se, std_effect=b[i] / rsd,
                        p_hc3_normal=2 * stats.norm.sf(abs(z)), residual_sd=rsd))
    Rm = np.zeros((len(levels), X.shape[1]))
    for k, i in enumerate(idx):
        Rm[k, i] = 1
    rb = Rm @ b
    W = float(rb @ np.linalg.inv(Rm @ V @ Rm.T) @ rb)
    # statsmodels wald_test(scalar=True) with HC3 uses F = W/q with F(q, n-k) reference
    q = len(levels)
    out.append(dict(analysis=label, contrast="omnibus_HC3_wald", n=n, estimate=W,
                    p_chi2=stats.chi2.sf(W, q), p_F=stats.f.sf(W / q, q, n - rank)))
    # prespecified partial F (classical) and Freedman-Lane permutation
    Xr = np.delete(X, idx, axis=1)
    br = np.linalg.lstsq(Xr, y, rcond=None)[0]
    fr = Xr @ br
    er = y - fr
    rss_r = float(er @ er)
    F = ((rss_r - rss) / q) / (rss / (n - rank))
    row = dict(analysis=label, contrast="omnibus_partial_F", n=n, estimate=F, p_F=stats.f.sf(F, q, n - rank))
    if nperm:
        rng = np.random.default_rng(seed)
        Q, _ = np.linalg.qr(X)
        Qr, _ = np.linalg.qr(Xr)
        cnt = 0
        B = 2000
        done = 0
        while done < nperm:
            m = min(B, nperm - done)
            P = np.stack([rng.permutation(n) for _ in range(m)], axis=1)
            Ys = fr[:, None] + er[P]
            rf = Ys - Q @ (Q.T @ Ys)
            rr = Ys - Qr @ (Qr.T @ Ys)
            rssf = (rf ** 2).sum(0)
            rssr = (rr ** 2).sum(0)
            Fs = ((rssr - rssf) / q) / (rssf / (n - rank))
            cnt += int((Fs >= F).sum())
            done += m
        row["perm_p"] = (1 + cnt) / (nperm + 1)
        row["perm_draws"] = nperm
    out.append(row)
    return out


res = []
v1 = prim[prim["primary_evidence_class"].isin(["neither", "genetic_only", "disease_state_only", "convergent"])]
res += class_fit(v1, "primary_evidence_class", ["genetic_only", "disease_state_only", "convergent"], "v1_reproduction", nperm=100000)
lv = ["genetically_anchored", "established_state_associated", "concordant"]
bb = prim[prim.r8.isin(["unresolved", "untested"] + lv)].copy()
R["myojin_r8_ref_plus_untested_n"] = len(bb)
R["myojin_r8_not_in_r8"] = int((prim.r8 == "not_in_r8").sum())
res += class_fit(bb, "r8", lv, "r8_ref_unresolved_plus_untested", nperm=100000)
aa = prim[prim.r8.isin(["unresolved"] + lv)].copy()
res += class_fit(aa, "r8", lv, "r8_ref_unresolved_only")
pd.DataFrame(res).to_csv(OUT / "myojin_r8_refit_check.tsv", sep="\t", index=False)
R["myojin_results"] = pd.DataFrame(res).astype(object).where(pd.notna(pd.DataFrame(res)), None).to_dict("records")
cand = pd.read_csv(CAND / "R10_03_myojin_r8_class_refit.tsv", sep="\t")
R["cand_R10_03"] = cand.astype(str).to_dict("records")


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, float) and np.isnan(o):
        return None
    return o


json.dump(clean(R), open(OUT / "f6_check.json", "w"), indent=1, default=str)
print("done")
