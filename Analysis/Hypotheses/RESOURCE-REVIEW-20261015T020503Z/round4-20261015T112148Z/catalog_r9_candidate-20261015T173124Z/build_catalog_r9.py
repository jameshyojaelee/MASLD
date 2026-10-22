#!/usr/bin/env python3
"""Gene Catalog r9 CANDIDATE build.

Starts from Catalog r8 (gene-catalog-corrected-r8-20261009T142027Z, reviewed copy
catalog-r8-reviewed-20261014T142026Z_22236566). The 11 r8 gene-index columns are copied
unchanged; no r8 gene call or primary_evidence_class is altered. New per-gene columns are
added beside them, keyed on ensembl_id, and per-plan or per-figure tables go to sidecars.

Reads saved tables only. No random draws are made; the seed is set for completeness.
All CSV/TSV numbers are parsed with float_precision='round_trip' (pandas' default fast
parser returns 0 for long fixed-point values such as 0.000000000000000000315506...).

Usage: python build_catalog_r9.py <existing empty output dir>
"""
import hashlib
import json
import random
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20260923
random.seed(SEED)
np.random.seed(SEED)

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
HYP = ROOT / "Analysis/Hypotheses"
R4 = HYP / "RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z"
CAT = ROOT / "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8"
COMP = ROOT / "figures/candidates/catalog-r8-completion-20261014T141645Z_22236527"
BULK = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
BULK_REPRO = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls"
COLOC = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
GEN = ROOT / "Analysis/Multimodal_Program_Projection/candidates/gen-corrected-r8-20261009T142027Z/genetics_context/frozen_evidence_classes.tsv"
FIG3E = ROOT / "figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions"
FIG3E_COUNTS = ROOT / "figures/main/fig3_bulk_transcriptomics/source_tables/fig3e_deg_counts.tsv"
FIG5B = ROOT / "figures/main/fig5_molecular_context/panels/data/fig5b_protein_triage_genes.tsv"
FIG5B_SUM = ROOT / "figures/main/fig5_molecular_context/panels/data/fig5b_protein_triage_summary.tsv"
LIN = ROOT / "figures/candidates/bulk-context-774-20261014T140904Z_22236311/class_difference"
F6P = ROOT / "figures/candidates/review-r3-fixes-20261015T112148Z/F6/panels_20261015T114537Z"
QUAR = CAT / "passport_identity_quarantine.tsv"

P = {  # plan: (run dir, RESULT file, shipped file)
    "H53": (HYP / "H53-20261015T112655Z", HYP / "H53-20261015T112655Z/RESULT_H53.md",
            HYP / "H53-20261015T112655Z/results/catalog_gene_caqtl_annotation.tsv"),
    "H55": (HYP / "H55-20261015T112558Z/run-20261015T114811Z", HYP / "H55-20261015T112558Z/RESULT_H55.md",
            HYP / "H55-20261015T112558Z/run-20261015T114811Z/outputs/catalog_column_gtex_postmortem_fat.tsv"),
    "H63": (HYP / "H63-20261015T125909Z/run-20261015T131038Z", HYP / "H63-20261015T125909Z/run-20261015T131038Z/RESULT_H63.md",
            HYP / "H63-20261015T125909Z/run-20261015T131038Z/h63_zonation_columns.parquet"),
    "H66": (HYP / "H66-20261015T124530Z/run-20261015T125923Z", HYP / "H66-20261015T124530Z/run-20261015T125923Z/RESULT_H66.md",
            HYP / "H66-20261015T124530Z/run-20261015T125923Z/h66_catalog_fields_shipped.tsv.gz"),
    "H78": (HYP / "H78-20261015T130152Z/run-20261015T131130Z", HYP / "H78-20261015T130152Z/run-20261015T131130Z/RESULT_H78.md",
            HYP / "H78-20261015T130152Z/run-20261015T131130Z/catalog_followup_annotation.tsv"),
    "H86": (HYP / "H86-20261015T125018Z/run-20261015T125926Z", HYP / "H86-20261015T125018Z/run-20261015T125926Z/RESULT_H86.md",
            HYP / "H86-20261015T125018Z/run-20261015T125926Z/T1_protein_coverage_by_class.tsv"),
    "H87": (HYP / "H87-20261015T130324Z/run-20261015T131921Z", HYP / "H87-20261015T130324Z/run-20261015T131921Z/RESULT_H87.md",
            HYP / "H87-20261015T130324Z/run-20261015T131921Z/h87_coding_cs_gene_column.tsv"),
    "H72": (HYP / "H72-20261015T134833Z/run-20261015T142309Z", HYP / "H72-20261015T134833Z/run-20261015T142309Z/RESULT_H72.md",
            HYP / "H72-20261015T134833Z/run-20261015T142309Z/catalog_sex_effect_state.tsv"),
    "H74": (HYP / "H74-20261015T124604Z/run-20261015T125815Z", HYP / "H74-20261015T124604Z/run-20261015T125815Z/RESULT_H74.md",
            HYP / "H74-20261015T124604Z/run-20261015T125815Z/locus_table.tsv"),
    "H71": (HYP / "H71-20261015T124239Z/run-20261015T125815Z", HYP / "H71-20261015T124239Z/run-20261015T125815Z/RESULT_H71.md",
            HYP / "H71-20261015T124239Z/run-20261015T125815Z/partA_requalification.tsv"),
    "H59": (HYP / "H59-20261015T112432Z", HYP / "H59-20261015T112432Z/RESULT_H59.md",
            HYP / "H59-20261015T112432Z/PRESPEC_H59.md"),
}
H86_D7 = P["H86"][0] / "D7_deposit_detection_by_class.tsv"
F6_JSON = R4 / "repair_inputs/F6.json"

OUT = Path(sys.argv[1])
assert OUT.is_dir(), OUT
INPUTS = {F6_JSON, R4 / "plan_results.json"} | {v[1] for v in P.values()} | {P["H59"][2]}
INV = []


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rd(p, sep="\t", **k):
    """Numeric-safe reader: round-trip floats; only empty/NA/NaN become missing ('None' stays text)."""
    INPUTS.add(Path(p))
    return pd.read_csv(p, sep=sep, float_precision="round_trip", low_memory=False,
                       keep_default_na=False, na_values=["", "NA", "NaN", "nan", "<NA>"], **k)


def rd_str(p, sep="\t"):
    INPUTS.add(Path(p))
    return pd.read_csv(p, sep=sep, dtype=str, keep_default_na=False)


def rpq(p):
    INPUTS.add(Path(p))
    return pd.read_parquet(p)


def check(name, observed, expected, hard=True, ok=None):
    if ok is None:
        ok = observed == expected
    INV.append({"check": name, "expected": str(expected), "observed": str(observed),
                "passed": bool(ok), "hard": hard})
    print(("PASS " if ok else ("FAIL " if hard else "DIFF ")) + name + f" | observed {observed} | expected {expected}")


def vc(s):
    return {str(k): int(v) for k, v in s.value_counts(dropna=False).items()}


# =========================================================================== r8
g8 = rpq(CAT / "passport_gene_index.parquet")
e8 = rpq(CAT / "passport_evidence_long.parquet")
check("r8 gene index sha256 (reviewed copy = 20261009 build)", sha256(CAT / "passport_gene_index.parquet"),
      "ffd3c4650c77a899e9ebed7327697398c4b6462cfffcd7f0e635cd0d023db06b")
check("r8 gene index rows", len(g8), 28233)
check("r8 ensembl_id unique", int(g8.ensembl_id.duplicated().sum()), 0)
check("r8 evidence rows", len(e8), 56466)
CLS_EXPECTED = {"untested": 13551, "unresolved": 13127, "established_state_associated": 1129,
                "genetically_anchored": 400, "concordant": 26}
check("r8 class counts", vc(g8.primary_evidence_class), CLS_EXPECTED)
gen = e8[e8.evidence_domain == "genetics"].set_index("ensembl_id")
rna = e8[e8.evidence_domain == "transcriptomics"].set_index("ensembl_id")
ids = pd.Index(g8.ensembl_id.astype(str), name="ensembl_id")
check("one genetics and one RNA evidence row per gene", (gen.index.is_unique and rna.index.is_unique
      and set(gen.index) == set(ids) and set(rna.index) == set(ids)), True)
gen = gen.reindex(ids)
rna = rna.reindex(ids)
cls = g8.set_index("ensembl_id").primary_evidence_class.reindex(ids)
sym8 = g8.set_index("ensembl_id").symbol.reindex(ids)
check("r8 genetics call states", vc(gen.call_state), {"indeterminate": 18478, "untestable": 9319, "supported": 436})
check("r8 RNA call states", vc(rna.call_state), {"indeterminate": 22148, "untestable": 4232, "supported": 1853})
mat = pd.crosstab(gen.call_state, rna.call_state)
mat_obs = [int(mat.loc[a, b]) if (a in mat.index and b in mat.columns) else 0
           for a in ["supported", "indeterminate", "untestable"] for b in ["supported", "indeterminate", "untestable"]]
check("Fig 6B 3x3 cells (gen rows x RNA cols)", mat_obs, [26, 400, 10, 1129, 13127, 4222, 698, 8621, 0])

N = pd.DataFrame(index=ids)  # new columns, same row order as r8
N["r8_genetics_call_state"] = gen.call_state.astype(str)
N["r8_genetics_testability_state"] = gen.testability_state.astype(str)
N["r8_rna_call_state"] = rna.call_state.astype(str)
N["r8_rna_testability_state"] = rna.testability_state.astype(str)

# ================================================================ universes
raw = rd_str(GEN)
deg = pd.read_csv(BULK / "deg_results.csv", float_precision="round_trip")
INPUTS.add(BULK / "deg_results.csv")
deg_repro = pd.read_csv(BULK_REPRO / "deg_results.csv", float_precision="round_trip")
INPUTS.add(BULK_REPRO / "deg_results.csv")
t12 = pd.read_csv(COLOC / "gene_level_coloc_tier12.csv", float_precision="round_trip", low_memory=False)
INPUTS.add(COLOC / "gene_level_coloc_tier12.csv")
deg["ens"] = deg.gene.str.split(".").str[0]
t12["ens"] = t12.ensembl.astype(str).str.split(".").str[0]
u_symbol = int((raw.joint_testable == "true").sum())
u_stable = len(set(deg.ens) & set(t12.ens))
u_joint = int((g8.primary_evidence_class != "untested").sum())
UNIV = pd.DataFrame([
    {"universe": "pre-adjudication symbol-grouped joint universe (GEN_RAW frozen classes, joint_testable true)",
     "n": u_symbol, "expected": 14720,
     "source": str(GEN)},
    {"universe": "same-stable-ID join: 774-unit bulk fit IDs intersect COLOC r2 Tier-1/2 IDs",
     "n": u_stable, "expected": 14704,
     "source": f"{BULK / 'deg_results.csv'} ; {COLOC / 'gene_level_coloc_tier12.csv'}"},
    {"universe": "Catalog r8 joint (primary_evidence_class not untested)", "n": u_joint, "expected": 14682,
     "source": str(CAT / "passport_gene_index.parquet")},
    {"universe": "Catalog r8 gene index rows", "n": len(g8), "expected": 28233,
     "source": str(CAT / "passport_gene_index.parquet")},
    {"universe": "Catalog r9 candidate gene table rows", "n": None, "expected": 28233, "source": "this build"},
])
check("universe 14,720 symbol", u_symbol, 14720)
check("universe 14,704 stable-ID", u_stable, 14704)
check("universe 14,682 Catalog joint", u_joint, 14682)
check("adopted vs reproduction bulk fit: same gene IDs", set(deg.gene) == set(deg_repro.gene), True)
dr = deg.set_index("gene").join(deg_repro.set_index("gene")[["logFC", "padj"]], rsuffix="_r")
check("adopted vs reproduction bulk fit: max |logFC diff| < 1e-12",
      float(np.abs(dr.logFC - dr.logFC_r).max()), "< 1e-12", ok=float(np.abs(dr.logFC - dr.logFC_r).max()) < 1e-12)

# ============================================================= RNA (774-unit fit)
coef = rd(BULK / "coefficient_table.tsv.gz")
coef["ens"] = coef.gene.str.split(".").str[0]
check("bulk fit genes", len(coef), 24115)
check("bulk fit stable IDs unique", int(coef.ens.duplicated().sum()), 0)
check("bulk fit P.Value and padj have no zeros at round-trip precision",
      int(((coef["P.Value"] == 0) | (coef.padj == 0)).sum()), 0)
cmp_ = coef.set_index("gene")[["logFC", "padj"]].join(deg.set_index("gene")[["logFC", "padj"]], rsuffix="_d")
check("coefficient_table vs deg_results logFC identical", float(np.abs(cmp_.logFC - cmp_.logFC_d).max()), 0.0,
      ok=float(np.abs(cmp_.logFC - cmp_.logFC_d).max()) < 1e-12)
cb = coef.set_index("ens")
canon = (cb.padj < 0.05) & (cb.logFC.abs() > 0.5)
check("bulk canonical calls (BH q<0.05, |log2FC|>0.5)", int(canon.sum()), 1854)
tested = ids.intersection(cb.index)
check("Catalog genes with a bulk fit row", len(tested), 24001)
est = rna.estimate.astype("float64")
d_lfc = float(np.nanmax(np.abs(est.loc[tested] - cb.logFC.loc[tested])))
check("r8 RNA estimate equals fit logFC (max abs diff)", d_lfc, "< 1e-10", ok=d_lfc < 1e-10)
q8 = rna.q_value.astype("float64")
relq = float(np.nanmax(np.abs(q8.loc[tested] - cb.padj.loc[tested]) / cb.padj.loc[tested]))
check("r8 RNA q equals fit BH q (max relative diff)", relq, "< 1e-9", ok=relq < 1e-9)
check("r8 RNA estimate missing exactly where no fit row",
      int(est.isna().sum()), len(ids) - len(tested))
N["rna_log2fc"] = cb.logFC.reindex(ids)
N["rna_se"] = cb.SE.reindex(ids)
N["rna_t"] = cb.t.reindex(ids)
N["rna_p"] = cb["P.Value"].reindex(ids)
N["rna_bh_q"] = cb.padj.reindex(ids)
N["rna_ave_log2cpm"] = cb.AveExpr.reindex(ids)
N["rna_bh_q_lt_0.05"] = (N.rna_bh_q < 0.05).where(N.rna_bh_q.notna()).astype("boolean")
N["rna_canonical_call"] = canon.reindex(ids).where(N.rna_bh_q.notna()).astype("boolean")
has = N.rna_log2fc.notna()
N["rna_n_units"] = pd.Series(np.where(has, 774, np.nan), index=ids).astype("Int64")
N["rna_n_control_units"] = pd.Series(np.where(has, 138, np.nan), index=ids).astype("Int64")
N["rna_n_disease_units"] = pd.Series(np.where(has, 636, np.nan), index=ids).astype("Int64")
N["rna_unit"] = np.where(has, "genetic-identity cluster (not verified donor)", "")
sup_rna = set(ids[rna.call_state == "supported"])
check("r8 RNA supported == canonical call among Catalog genes", sup_rna == set(tested[canon.reindex(tested).values]), True)
missing_canon = sorted(set(cb.index[canon]) - set(ids))
check("canonical call without a Catalog record", missing_canon, ["ENSG00000278732"])
ind = rna.call_state == "indeterminate"
check("RNA indeterminate with BH q<0.05 (|log2FC|<=0.5)", int((ind & (N.rna_bh_q < 0.05)).sum()), 11301)

design = rd(BULK / "model_design.tsv")
check("774-unit design rows / disease / control", (len(design), int(design.group_binaryDisease.sum()),
      int((design.group_binaryDisease == 0).sum())), (774, 636, 138))

# ======================================================= genetics (COLOC r2)
check("COLOC Tier-1/2 stable IDs unique", int(t12.ens.duplicated().sum()), 0)
ct = t12.set_index("ens")
anc_parts = []
for ch in pd.read_csv(COLOC / "susie_coloc_all_gwas.csv", usecols=["gwas_name", "ancestry"], dtype=str,
                      keep_default_na=False, chunksize=1_000_000):
    anc_parts.append(ch.drop_duplicates())
INPUTS.add(COLOC / "susie_coloc_all_gwas.csv")
amap = pd.concat(anc_parts).drop_duplicates()
check("each COLOC study has one ancestry", int(amap.gwas_name.duplicated().sum()), 0)
amap = amap.set_index("gwas_name").ancestry
N["gen_coloc_best_susie_pp4"] = ct.coloc_best_susie_pp4.reindex(ids)
N["gen_coloc_best_susie_study"] = ct.susie_driving_gwas.reindex(ids)
N["gen_coloc_best_susie_ancestry"] = N.gen_coloc_best_susie_study.map(amap)
N["gen_coloc_best_abf_pp4"] = ct.coloc_best_abf_pp4.reindex(ids)
N["gen_coloc_best_abf_study"] = ct.abf_driving_gwas.reindex(ids)
N["gen_coloc_best_abf_ancestry"] = N.gen_coloc_best_abf_study.map(amap)
N["gen_susie_state_t12"] = ct.susie_state_t12.reindex(ids)
check("SuSiE ancestry filled only where a SuSiE study exists",
      int((N.gen_coloc_best_susie_ancestry.notna() & N.gen_coloc_best_susie_study.isna()).sum()), 0)
check("every SuSiE/ABF study mapped to an ancestry",
      int((N.gen_coloc_best_susie_study.notna() & N.gen_coloc_best_susie_ancestry.isna()).sum()
          + (N.gen_coloc_best_abf_study.notna() & N.gen_coloc_best_abf_ancestry.isna()).sum()), 0)
check("SuSiE PP.H4 present iff SuSiE study present",
      bool((N.gen_coloc_best_susie_pp4.notna() == N.gen_coloc_best_susie_study.notna()).all()), True)
gest = gen.estimate.astype("float64")
m_su = np.isclose(gest, N.gen_coloc_best_susie_pp4, atol=1e-9, rtol=0)
m_ab = np.isclose(gest, N.gen_coloc_best_abf_pp4, atol=1e-9, rtol=0)
meth = pd.Series(np.select([m_su, m_ab, gest.isna()], ["SuSiE", "ABF", "none"], "other"), index=ids)
f6six = rd_str(F6P / "R07_11_six_genes.tsv")
six_ids = set(f6six.catalog_id)
meth[meth == "other"] = "value_from_other_ensembl_id"
N["gen_r8_displayed_pp4_method"] = meth
mx = pd.crosstab(gen.call_state, meth)
check("R15-02 displayed PP.H4 method: supported all SuSiE", int(mx.loc["supported", "SuSiE"]), 436)
check("R15-02 indeterminate ABF / SuSiE / other-ID",
      (int(mx.loc["indeterminate"].get("ABF", 0)), int(mx.loc["indeterminate"].get("SuSiE", 0)),
       int(mx.loc["indeterminate"].get("value_from_other_ensembl_id", 0))), (15231, 3241, 6))
check("other-ID displayed values are the six R07-11 genes", set(ids[meth == "value_from_other_ensembl_id"]) == six_ids, True)
check("genetics supported iff best Tier-1/2 SuSiE PP.H4 > 0.5",
      set(ids[gen.call_state == "supported"]) == set(ids[N.gen_coloc_best_susie_pp4 > 0.5]), True)
N["gen_estimate_from_other_ensembl_id"] = pd.Series(ids.map(dict(zip(f6six.catalog_id, f6six.ensembl_genetic))),
                                                    index=ids).fillna("")
up_ids = ids[gen.testability_state == "underpowered_source"]
isp_ids = ids[gen.testability_state == "insufficient_shared_posterior"]
check("underpowered_source genes with no Tier-1/2 COLOC record", (len(up_ids), int(pd.Index(up_ids).isin(ct.index).sum())), (9312, 0))
check("insufficient_shared_posterior genes: all SuSiE-untestable in >=1 Tier-1/2 study",
      (len(isp_ids), int((ct.susie_state_t12.reindex(isp_ids) == "untestable").sum())), (7, 7))
N["gen_untestable_reason_r9"] = np.select(
    [gen.testability_state == "underpowered_source", gen.testability_state == "insufficient_shared_posterior"],
    ["no_tier12_coloc_record", "susie_untestable_in_ge1_tier12_study"], "")

# ================================================= class display and untested split
DISPLAY = {"concordant": "both supported", "genetically_anchored": "genetically anchored",
           "established_state_associated": "established state-associated", "unresolved": "unresolved",
           "untested": "untested"}
N["class_display_r9"] = cls.map(DISPLAY)
unt = cls == "untested"
both_unt = unt & (gen.call_state == "untestable") & (rna.call_state == "untestable")
check("no gene untestable in both domains", int(both_unt.sum()), 0)
N["untested_missing_domain"] = np.select([unt & (rna.call_state == "untestable"), unt & (gen.call_state == "untestable")],
                                         ["no_rna_call_not_detected", "no_genetic_call"], "")
N["untested_one_domain_supported"] = np.select([unt & (rna.call_state == "supported"), unt & (gen.call_state == "supported")],
                                               ["rna_supported", "genetics_supported"], "")
check("untested split by missing domain", vc(N.untested_missing_domain[unt]),
      {"no_genetic_call": 9319, "no_rna_call_not_detected": 4232})
check("untested with one supported domain (708)", vc(N.untested_one_domain_supported[unt & (N.untested_one_domain_supported != "")]),
      {"rna_supported": 698, "genetics_supported": 10})

# ============================================ corrected assay applicability (R15-05)
mo_vals = {"bulk_rna": {}, "proteomics": {}, "chromatin": {}, "spatial": {}}
mo_path = COMP / "molecular_objects/molecular_object_sidecar.jsonl"
INPUTS.add(mo_path)
with open(mo_path) as fh:
    for line in fh:
        aa = json.loads(json.loads(line)["assay_applicability"])
        for k in mo_vals:
            mo_vals[k][aa[k]] = mo_vals[k].get(aa[k], 0) + 1
N["assay_bulk_rna_r9"] = rna.call_state.astype(str)
N["assay_proteomics_r9"] = "not_assessed_at_gene_level"
N["assay_chromatin_r9"] = "not_joined_at_gene_level"
N["assay_spatial_r9"] = "not_joined_at_gene_level"

# ===================================================================== H53
h53 = rd(P["H53"][2])
check("H53 rows / unique ensembl", (len(h53), int(h53.ensembl.duplicated().sum())), (436, 0))
check("H53 genes == r8 genetics-supported genes", set(h53.ensembl) == set(ids[gen.call_state == "supported"]), True)
h = h53.set_index("ensembl").reindex(ids)
for c in ["caqtl_region_state", "caqtl_region_hg38", "caqtl_pp_h4", "caqtl_lead_ld_r2", "caqtl_three_way_eqtl",
          "caqtl_abc_predicted_link", "caqtl_source_flag"]:
    N[c] = h[c]
N["caqtl_n_gwas_signals"] = h.n_gwas_signals.astype("Int64")
ev = N.caqtl_region_state.notna()
N["caqtl_region_state"] = N.caqtl_region_state.where(ev, "not_evaluated_genetics_not_supported")
tw = N.caqtl_three_way_eqtl.astype(str).str.lower()
car = N.caqtl_region_state == "colocalized accessible region"
N["caqtl_gene_link"] = np.select([car & (tw == "true"), car & (tw == "false"), ev],
                                 ["M1_three_way_supported_linked", "M2_r2_pair_only_not_called_linked", "not_linked"],
                                 "not_evaluated_genetics_not_supported")
check("H53 M1 / M2 gene counts", vc(N.caqtl_gene_link[N.caqtl_gene_link.str.startswith("M")]),
      {"M1_three_way_supported_linked": 108, "M2_r2_pair_only_not_called_linked": 35})

# ===================================================================== H55
h55 = rd(P["H55"][2])
check("H55 rows (canonical calls)", len(h55), 1854)
check("H55 keys not in Catalog", sorted(set(h55.ensembl_id) - set(ids)), ["ENSG00000278732"])
h = h55.set_index("ensembl_id").reindex(ids)
for c in [c for c in h55.columns if c != "ensembl_id"]:
    N[c] = h[c]
N["gtex_postmortem_fat_state"] = N.gtex_postmortem_fat_state.fillna("not_a_canonical_bulk_call")
check("H55 states over the 1,853 Catalog canonical calls",
      vc(N.gtex_postmortem_fat_state[N.gtex_postmortem_fat_state != "not_a_canonical_bulk_call"]),
      {k: v for k, v in vc(h55.set_index("ensembl_id").drop(index="ENSG00000278732").gtex_postmortem_fat_state).items()})

# ===================================================================== H63
h63 = rpq(P["H63"][2])
check("H63 keys == Catalog", (len(h63), set(h63.ensembl_id) == set(ids)), (28233, True))
h = h63.set_index("ensembl_id").reindex(ids)
check("H63 symbol equals r8 symbol", int((h.symbol.astype(object).fillna("").astype(str) != sym8.astype(object).fillna("").astype(str)).sum()), 0)
for c in [c for c in h63.columns if c not in ("ensembl_id", "symbol")]:
    N[c] = h[c]
check("H63 zon_layers_flag counts", vc(N.zon_layers_flag),
      {"neither": 26988, "one_layer_only": 1062, "agree": 162, "disagree": 21})

# ===================================================================== H66
h66 = rd(P["H66"][2])
check("H66 keys == Catalog", (len(h66), set(h66.ensembl_id) == set(ids)), (28233, True))
h = h66.set_index("ensembl_id").reindex(ids)
check("H66 symbol equals r8 symbol", int((h.symbol.astype(object).fillna("").astype(str) != sym8.astype(object).fillna("").astype(str)).sum()), 0, hard=False)
for c in [c for c in h66.columns if c not in ("ensembl_id", "symbol")]:
    N[c] = h[c]
check("H66 liver_proteome_any / plasma_proteome_any true",
      (int((N.liver_proteome_any.astype(str).isin(["1", "1.0", "True", "true"])).sum()),
       int((N.plasma_proteome_any.astype(str).isin(["1", "1.0", "True", "true"])).sum())), (7137, 579), hard=False)

# ===================================================================== H87
h87 = rd(P["H87"][2])
check("H87 keys == Catalog", (len(h87), set(h87.ensembl_id) == set(ids)), (28233, True))
h = h87.set_index("ensembl_id").reindex(ids)
check("H87 primary_evidence_class equals r8", int((h.primary_evidence_class != cls).sum()), 0)
for c in [c for c in h87.columns if c not in ("ensembl_id", "symbol", "primary_evidence_class")]:
    N[c] = h[c]
check("H87 coding_cs_state counts", vc(N.coding_cs_state),
      {"no_coding_cs_member": 28137, "coding_cs_member_pip_lt_0.1": 50, "fine_mapped_coding_pip_ge_0.5": 31,
       "coding_cs_member_pip_0.1_to_0.5": 15})

# ===================================================================== H72
h72 = rd(P["H72"][2])
check("H72 rows == anchored + concordant", (len(h72), set(h72.ensembl_id) == set(ids[cls.isin(["genetically_anchored", "concordant"])])),
      (426, True))
check("H72 catalog_class equals r8", int((h72.set_index("ensembl_id").catalog_class != cls.reindex(h72.ensembl_id).values).sum()), 0)
N["enzyme_sex_effect_state"] = h72.set_index("ensembl_id").state.reindex(ids).fillna("not_in_tested_family")
check("H72 state counts", vc(h72.state), {"underpowered": 271, "not_applicable": 123, "heterogeneous": 13,
                                          "one_scale_only": 13, "scale_dependent": 4, "untestable": 2})

# ================================================================== Fig 3E
st = rd(FIG3E / "stage_extension_all_gene_results.tsv")
check("Fig 3E file sha256", sha256(FIG3E / "stage_extension_all_gene_results.tsv")[:8], "9e50bee2")
check("Fig 3E rows", len(st), 337610)
KEEP = [("fibrosis", "F1_vs_F0", "F1"), ("fibrosis", "F2_vs_F0", "F2"), ("fibrosis", "F3_vs_F0", "F3"),
        ("fibrosis", "F4_vs_F0", "F4"), ("NAS", "NAS1_2_vs_NAS0", "NAS1-2"), ("NAS", "NAS3_4_vs_NAS0", "NAS3-4"),
        ("NAS", "NAS5_8_vs_NAS0", "NAS5-8")]
counts = rd(FIG3E_COUNTS)
long_parts = []
for axis, con, disp in KEEP:
    s = st[(st.axis == axis) & (st.contrast == con)].copy()
    check(f"Fig 3E {con}: one row per gene", (len(s), int(s.gene_id_base.duplicated().sum())), (24115, 0))
    hi = int(((s.FDR < 0.05) & (s.logFC > 0)).sum())
    lo = int(((s.FDR < 0.05) & (s.logFC < 0)).sum())
    exp_hi = int(counts[(counts.axis == ("fibrosis" if axis == "fibrosis" else "NAS")) & (counts.display == disp) & (counts.direction == "Higher")].n.iloc[0])
    exp_lo = int(counts[(counts.axis == ("fibrosis" if axis == "fibrosis" else "NAS")) & (counts.display == disp) & (counts.direction == "Lower")].n.iloc[0])
    check(f"Fig 3E {con} higher/lower at BH q<0.05 (no FC filter) vs installed source table", (hi, lo), (exp_hi, exp_lo))
    s["shown_on_fig3e_panel"] = axis == "fibrosis"
    long_parts.append(s)
    if axis == "fibrosis":
        k = con.split("_")[0]
        ss = s.set_index("gene_id_base").reindex(ids)
        N[f"fig3e_{k}_vs_F0_log2fc"] = ss.logFC
        N[f"fig3e_{k}_vs_F0_bh_q"] = ss.FDR
        N[f"fig3e_{k}_vs_F0_call"] = np.select([ss.logFC.isna(), (ss.FDR < 0.05) & (ss.logFC > 0), (ss.FDR < 0.05) & (ss.logFC < 0)],
                                               ["not_in_fit", f"higher_in_{k}", f"lower_in_{k}"], "not_called")
fig3e_long = pd.concat(long_parts, ignore_index=True)
fig3e_long.insert(0, "ensembl_id", fig3e_long.gene_id_base)
fig3e_long["in_catalog_r8"] = fig3e_long.ensembl_id.isin(set(ids))
man = rd(FIG3E / "stage_extension_sample_manifest.tsv")
fib = man[man.axis == "fibrosis"].drop_duplicates("analysis_unit_id")
check("Fig 3E fibrosis units / cohorts", (len(fib), sorted(fib.dataset.unique())), (387, ["GSE130970", "GSE135251", "GSE162694"]))
check("Fig 3E GSE213621 units", int((man.dataset == "GSE213621").sum()), 0)
f0 = fib[fib.group == "F0"].sample_id
f0d = design.set_index("sample_id").group_binaryDisease.reindex(f0)
check("Fig 3E F0 reference: source controls / F0 disease", (int((f0d == 0).sum()), int((f0d == 1).sum()), int(f0d.isna().sum())), (14, 89, 0))
fib_by_ds = fib.groupby("dataset").size().to_dict()

# ================================================================== Fig 5B
f5 = rd(FIG5B)
check("Fig 5B rows / unique ensembl", (len(f5), int(f5.ensembl_id.duplicated().sum())), (3054, 0))
check("Fig 5B protein-significant (confirmed + discordant)", vc(f5["class"]), {"not_corr": 2843, "confirmed": 188, "discordant": 23})
f5.insert(2, "in_catalog_r8", f5.ensembl_id.isin(set(ids)))
f5["release_note"] = ("Figure 5B source values (PXD051911 liver DIA, fixed 3,054-gene prioritized universe); "
                      "not a Catalog evidence call; proteomics gene_evidence_authorized false")
INPUTS.add(FIG5B_SUM)

# ================================================================ lineage
lin = rd(LIN / "results/gene_classes.tsv.gz")
check("lineage table rows (stable-ID universe) / with hep_specificity", (len(lin), int(lin.hep_specificity.notna().sum())), (14704, 14547))
li = lin.set_index("g").reindex(ids)
N["lineage_hep_specificity_log2cpm"] = li.hep_specificity
N["lineage_measured_in_sc"] = li.in_sc.astype("boolean")
lin_not_in_cat = int((~lin.g.isin(set(ids))).sum())
INPUTS.add(LIN / "results/substrate_provenance.json")
INPUTS.add(LIN / "00_assemble_substrate.py")

# ========================================================= assemble the wide table
assert list(N.index) == list(g8.ensembl_id)
W = pd.concat([g8.reset_index(drop=True), N.reset_index(drop=True)], axis=1)
check("r9 rows", len(W), 28233)
check("r8 columns unchanged in r9 table", W[g8.columns].equals(g8), True)
check("r9 class counts unchanged", vc(W.primary_evidence_class), CLS_EXPECTED)
check("no duplicate column names", int(pd.Index(W.columns).duplicated().sum()), 0)
UNIV.loc[UNIV.universe.str.startswith("Catalog r9"), "n"] = len(W)

# ================================================================ sidecars
SIDE = {}
fig3e_long.to_csv(OUT / "sidecar_fig3e_stage_contrasts_long.tsv.gz", sep="\t", index=False,
                  compression={"method": "gzip", "mtime": 0})
SIDE["sidecar_fig3e_stage_contrasts_long.tsv.gz"] = fig3e_long
f5.to_csv(OUT / "sidecar_fig5b_liver_protein.tsv", sep="\t", index=False)
SIDE["sidecar_fig5b_liver_protein.tsv"] = f5
for plan, name in [("H72", "sidecar_h72_enzyme_sex_effect_detail.tsv"), ("H74", "sidecar_h74_nearest_gene_loci.tsv"),
                   ("H78", "sidecar_h78_followup_annotation.tsv")]:
    shutil.copyfile(P[plan][2], OUT / name)
    INPUTS.add(P[plan][2])
    SIDE[name] = pd.read_csv(OUT / name, sep="\t", float_precision="round_trip", low_memory=False,
                             keep_default_na=False, na_values=["", "NA", "NaN", "nan", "<NA>"])
h78 = SIDE["sidecar_h78_followup_annotation.tsv"]
check("H78 keys == Catalog", (len(h78), set(h78.ensembl_id) == set(ids)), (28233, True))
check("H74 loci", len(SIDE["sidecar_h74_nearest_gene_loci.tsv"]), 259)

h71 = rd_str(P["H71"][2])
q71 = pd.to_numeric(h71.bulk_q.replace("", np.nan))
h71["bulk_q_full_precision_r9"] = h71.ensembl_id.map(cb.padj).map(lambda v: "" if pd.isna(v) else repr(float(v)))
n71_zero = int((q71 == 0).sum())
h71.to_csv(OUT / "sidecar_h71_competitor_list_states.tsv", sep="\t", index=False)
SIDE["sidecar_h71_competitor_list_states.tsv"] = h71
check("H71 rows", len(h71), 1999)

t1 = rd(P["H86"][2])
d7 = rd(H86_D7)
t1 = t1.rename(columns={"row": "class"})
h86 = t1[["class", "n_genes", "n_quantified", "pct_quantified", "n_quantified_direct", "n_quantified_remapped"]].merge(
    d7[["class", "n_detected_in_deposit", "pct_detected", "n_detected_only_in_multigene_group",
        "n_detected_not_quantified", "n_quantified_not_detected"]], on="class", how="left")
h86.to_csv(OUT / "sidecar_h86_protein_coverage_by_class.tsv", sep="\t", index=False)
SIDE["sidecar_h86_protein_coverage_by_class.tsv"] = h86
check("H86 quantified per class", dict(zip(h86["class"][:5], h86.n_quantified[:5].astype(int))),
      {"concordant": 11, "established_state_associated": 173, "genetically_anchored": 203, "unresolved": 5241, "untested": 414})

# R07-05: IDs in a source table with no Catalog record
q = rd_str(QUAR)
qmap = {}
for r in q.itertuples(index=False):
    for x in set((r.normalized_ensembl_candidates + ";" + r.source_ensembl_text).split(";")):
        x = x.strip().split(".")[0]
        if x:
            qmap.setdefault(x, set()).add(r.reason_code)
pm = rpq(CAT / "passport_program_membership.parquet")
pm_ids = pm.ensembl_id.astype(str).str.split(".").str[0]
rows = []
cat_set = set(ids)
for src, idset, extra in [
        ("bulk_774_fit", set(coef.ens), lambda e: {"symbol": cb.symbol.get(e), "bulk_log2fc": cb.logFC.get(e), "bulk_bh_q": cb.padj.get(e),
                                                    "bulk_canonical_call": bool(canon.get(e))}),
        ("coloc_r2_tier12", set(t12.ens), lambda e: {"coloc_best_susie_pp4": ct.coloc_best_susie_pp4.get(e), "coloc_best_abf_pp4": ct.coloc_best_abf_pp4.get(e)}),
        ("program_membership_r8", set(pm_ids), lambda e: {"n_memberships": int((pm_ids == e).sum())})]:
    for e in sorted(idset - cat_set):
        d = {"source": src, "ensembl_id": e, "quarantine_reason": ";".join(sorted(qmap.get(e, []))) or "not_quarantined"}
        d.update(extra(e))
        rows.append(d)
absent = pd.DataFrame(rows)
absent.to_csv(OUT / "sidecar_ids_without_catalog_record.tsv", sep="\t", index=False)
SIDE["sidecar_ids_without_catalog_record.tsv"] = absent
ab_counts = absent.groupby(["source", "quarantine_reason"]).size()
check("R07-05 IDs without record: bulk / COLOC / program", (int((absent.source == "bulk_774_fit").sum()),
      int((absent.source == "coloc_r2_tier12").sum()), int((absent.source == "program_membership_r8").sum())), (114, 61, 76))
check("R07-05 quarantine reasons (checker: bulk 109 multi; COLOC 32 multi, 22 absent v49, 7 none; program 0)",
      {f"{a}|{b}": int(v) for (a, b), v in ab_counts.items()}, "see checker", hard=False, ok=True)

# program context (14 rows): r9 call_state for the Vu Stromal ECM row, selection basis under H59 null
pc = rpq(CAT / "passport_program_context.parquet")
pcs = pc[["program_context_id", "program_uid", "program_label", "dataset_id", "assay", "estimate", "q_value",
          "call_state", "provenance_state", "source_dependent"]].rename(
    columns={"call_state": "call_state_r8", "provenance_state": "provenance_state_r8", "source_dependent": "source_dependent_r8"})
vu = (pc.program_label == "Stromal ECM (IGFBP7)") & (pc.dataset_id == "Vu_et_al_2025")
check("Vu Stromal ECM row: one row, r8 call supported", (int(vu.sum()), pc.loc[vu, "call_state"].iloc[0]), (1, "supported"))
pcs["call_state_r9_candidate"] = np.where(vu, "source_dependent", pcs.call_state_r8)
pcs["selection_basis_r9"] = ("prespecified programs (H59 null: 0 of 117 programs pass the corrective stage rule; "
                             "the 14 rows were the two robust programs of the 2026-08-07 release)")
pcs.to_csv(OUT / "sidecar_program_context_r9.tsv", sep="\t", index=False)
SIDE["sidecar_program_context_r9.tsv"] = pcs

# rules
rules = pd.DataFrame([
    ("genetics_call", "supported iff best eligible Tier-1/2 SuSiE PP.H4 > 0.5; ABF PP.H4 is shown but not used for the call", "r8 rule, unchanged; F6 R15-02"),
    ("rna_call", "supported iff 774-unit limma-voom BH q < 0.05 and |log2FC| > 0.5 (BH over 24,115 genes)", "r8 rule, unchanged; RNA indeterminate includes 11,301 genes with q < 0.05 and |log2FC| <= 0.5 (F6 R07-04)"),
    ("negative_call", "none defined for either domain; tested_negative not assessed", "F6 R07-02; defining a PP.H3-based rule is a new analysis, not done"),
    ("class_concordant_display", "shown as 'both supported'; effect directions are not compared", "F6 R07-06; class value 'concordant' unchanged"),
    ("genetics_untestable_reason", "no_tier12_coloc_record (9,312) replaces r8 'underpowered_source' (no power assessment made); susie_untestable_in_ge1_tier12_study (7)", "F6 R07-07 with checker wording"),
    ("fig3e_stage_call", "BH q < 0.05 within each contrast (BH over 24,115 genes), no fold-change floor; covariates dataset + inferred sex", "Fig 3E installed counts"),
    ("same_id_rule_R07_11", "not applied: six genes keep r8 class (would move unresolved 13,127 -> 13,121, untested 13,551 -> 13,557)", "F6 R07-11; Catalog rule decision"),
], columns=["rule", "definition", "source_note"])
rules.to_csv(OUT / "catalog_rules_r9.tsv", sep="\t", index=False)
SIDE["catalog_rules_r9.tsv"] = rules

# ======================================================== column dictionary
RS = "candidate"
G = {}  # column -> metadata


def meta(cols, **m):
    for c in cols:
        G[c] = m


r8_meta = dict(source="Catalog r8 gene index (copied unchanged)", run_dir=str(CAT), result="docs/RESULTS.md (Gene Catalog r8)",
               source_file=str(CAT / "passport_gene_index.parquet"), unit="text", release_state="r8 (unchanged)",
               blocked_on="", caveat="")
meta(list(g8.columns), **r8_meta)
meta(["r8_genetics_call_state", "r8_genetics_testability_state", "r8_rna_call_state", "r8_rna_testability_state"],
     source="Catalog r8 evidence rows (copied unchanged)", run_dir=str(CAT), result="docs/RESULTS.md (Gene Catalog r8)",
     source_file=str(CAT / "passport_evidence_long.parquet"), unit="state", release_state="r8 (unchanged)", blocked_on="",
     caveat="copied so the 3x3 call matrix can be read from this table; r8 'underpowered_source' means no COLOC record, not a power result")
rna_src = dict(source="F6 R15-03 (corrected): adopted 774-unit bulk refit", run_dir=str(BULK), result=str(F6_JSON) + " (R15-03, R07-04)",
               source_file=str(BULK / "coefficient_table.tsv.gz"), release_state=RS, blocked_on="H79;H81",
               caveat="138 control / 636 disease genetic-identity clusters (not verified donors); covariates dataset + inferred sex; parsed at round-trip precision (F6 sidecar printed q=0 for 736 and p=0 for 956 genes because of pandas fast-parser precision, not the file)")
meta(["rna_log2fc"], unit="log2 fold change, disease vs control", **rna_src)
meta(["rna_se"], unit="log2 fold change SE", **rna_src)
meta(["rna_t"], unit="moderated t", **rna_src)
meta(["rna_p"], unit="p value (two-sided, moderated t)", **rna_src)
meta(["rna_bh_q"], unit="BH q over 24,115 genes", **rna_src)
meta(["rna_ave_log2cpm"], unit="mean log2 CPM over 774 units", **rna_src)
meta(["rna_bh_q_lt_0.05"], unit="boolean", **rna_src)
meta(["rna_canonical_call"], unit="boolean (BH q<0.05 and |log2FC|>0.5)", **rna_src)
meta(["rna_n_units", "rna_n_control_units", "rna_n_disease_units"], unit="count of analysis units", **rna_src)
meta(["rna_unit"], unit="text", **rna_src)
gen_src = dict(source="F6 R15-02/R15-03 (corrected): COLOC r2 Tier-1/2 gene table", run_dir=str(COLOC),
               result=str(F6_JSON) + " (R15-02, R15-03)", source_file=str(COLOC / "gene_level_coloc_tier12.csv"),
               release_state=RS, blocked_on="corrected 50-study COLOC rerun (rebuild)",
               caveat="same-Ensembl-ID join; ancestry mapped from the named study (fixes F6 sidecar, which took ancestry from a different study for 950 genes)")
meta(["gen_coloc_best_susie_pp4", "gen_coloc_best_abf_pp4"], unit="PP.H4 (posterior probability, 0-1)", **gen_src)
meta(["gen_coloc_best_susie_study", "gen_coloc_best_abf_study"], unit="GWAS study name", **gen_src)
meta(["gen_coloc_best_susie_ancestry", "gen_coloc_best_abf_ancestry"], unit="ancestry label of that study (susie_coloc_all_gwas.csv)", **gen_src)
meta(["gen_susie_state_t12"], unit="state", **gen_src)
meta(["gen_r8_displayed_pp4_method"], unit="method of the value in the r8 genetics estimate field", **{**gen_src, "caveat": "r8 shows SuSiE for 436 supported and ABF for 15,231 indeterminate genes in one column"})
meta(["gen_estimate_from_other_ensembl_id"], unit="Ensembl ID whose COLOC record supplies the r8 genetics estimate",
     **{**gen_src, "result": str(F6_JSON) + " (R07-11)", "source_file": str(F6P / "R07_11_six_genes.tsv"),
        "caveat": "six genes; the other ID is absent from GENCODE v49; classes kept (same-ID rule is a Catalog decision)"})
meta(["gen_untestable_reason_r9"], unit="state", **{**gen_src, "result": str(F6_JSON) + " (R07-07 and checker)",
     "source_file": str(CAT / "passport_evidence_long.parquet") + " ; " + str(COLOC / "gene_level_coloc_tier12.csv")})
cls_src = dict(source="F6 R07-06/R07-07 display fields", run_dir=str(F6P), result=str(F6_JSON) + " (R07-06, R07-07)",
               source_file=str(CAT / "passport_gene_index.parquet") + " ; " + str(CAT / "passport_evidence_long.parquet"),
               release_state=RS, blocked_on="corrected 50-study COLOC rerun (rebuild)", caveat="display only; primary_evidence_class unchanged")
meta(["class_display_r9"], unit="text", **cls_src)
meta(["untested_missing_domain", "untested_one_domain_supported"], unit="state (blank if class is not untested)", **cls_src)
aa_src = dict(source="F6 R15-05 (corrected assay applicability)", run_dir=str(F6P), result=str(F6_JSON) + " (R15-05)",
              source_file=str(mo_path), release_state=RS, blocked_on="",
              caveat=("r8 molecular-object values: " + json.dumps(mo_vals) +
                      "; proteomics gene evidence is not authorized (passport_domain_coverage program_context_only), so no per-gene protein call; "
                      "a per-gene protein-quantified column is a user decision (H86 prespec)"))
meta(["assay_bulk_rna_r9", "assay_proteomics_r9", "assay_chromatin_r9", "assay_spatial_r9"], unit="state", **aa_src)
h53_src = dict(source="H53 caQTL region links", run_dir=str(P["H53"][0]), result=str(P["H53"][1]), source_file=str(P["H53"][2]),
               release_state=RS, blocked_on="corrected 50-study COLOC rerun (rebuild)",
               caveat="436 genetics-supported genes evaluated; M1 (three-way) called linked; M2 reached only through the r2 GWAS-eQTL pair, not called linked; Broadaway 2024 not checked; Currin 2025 flag in caqtl_source_flag; GEN_RAW classes unchanged")
meta(["caqtl_region_state", "caqtl_gene_link", "caqtl_source_flag"], unit="state", **h53_src)
meta(["caqtl_region_hg38"], unit="GRCh38 region chr:start-end", **h53_src)
meta(["caqtl_pp_h4"], unit="PP.H4 GWAS-caQTL", **h53_src)
meta(["caqtl_lead_ld_r2"], unit="LD r2 between leads", **h53_src)
meta(["caqtl_three_way_eqtl", "caqtl_abc_predicted_link"], unit="boolean", **h53_src)
meta(["caqtl_n_gwas_signals"], unit="count of liver GWAS signals", **h53_src)
h55_src = dict(source="H55 GTEx postmortem histologic fat", run_dir=str(P["H55"][0]), result=str(P["H55"][1]), source_file=str(P["H55"][2]),
               release_state=RS, blocked_on="",
               caveat="222 GTEx postmortem donors; claim limited to histologic steatosis of unrecorded cause in postmortem liver; NC1 tissue-area rho 0.354 (p 0.0564) close to S7; NC2 ischemic prediction missed; BH over 1,002 measurable genes (family C); not an evidence class")
meta(["gtex_postmortem_fat_slope", "gtex_postmortem_fat_ci_low", "gtex_postmortem_fat_ci_high"],
     unit="log2 CPM change from lowest- to highest-fat donor (exposure = fat-fraction rank scaled to 0-1)", **h55_src)
meta(["gtex_postmortem_fat_q"], unit="BH q over 1,002 genes", **h55_src)
meta(["gtex_postmortem_fat_state"], unit="state", **h55_src)
h63_dict = rd(P["H63"][0] / "h63_zonation_columns_dictionary.tsv").set_index("column").definition
for c in [c for c in h63.columns if c not in ("ensembl_id", "symbol")]:
    G[c] = dict(source="H63 healthy-liver zonation", run_dir=str(P["H63"][0]), result=str(P["H63"][1]), source_file=str(P["H63"][2]),
                unit=str(h63_dict.get(c, "")), release_state=RS, blocked_on="",
                caveat="healthy (non-steatotic) liver only; Yakubovsky 2026 and Weiss 2026 values source-owned; per-class counts need recompute after COLOC rerun")
h66_src = dict(source="H66 follow-up assay observability", run_dir=str(P["H66"][0]), result=str(P["H66"][1]), source_file=str(P["H66"][2]),
               release_state=RS, blocked_on="",
               caveat="not observed in an assay means untestable there, not negative; plasma proteome failed A3 (abundance check) and ships with that failure noted; S7 secondary class description not shipped; joining is a user decision")
for c in [c for c in h66.columns if c not in ("ensembl_id", "symbol")]:
    u = ("fraction of participants/donors" if c.endswith("_frac") else "boolean/indicator" if c.endswith(("_any", "_50", "sc_detected"))
         else "count of genes sharing the group" if c.endswith("_shared_group") else "text")
    G[c] = dict(unit=u, **h66_src)
h87_src = dict(source="H87 fine-mapped coding evidence", run_dir=str(P["H87"][0]), result=str(P["H87"][1]), source_file=str(P["H87"][2]),
               release_state=RS, blocked_on="corrected 50-study COLOC rerun (text)",
               caveat="depends on fine-mapping release noncoding-dna-precoloc-candidate-2026-08-12 (frozen candidate, not promoted); does not change the eQTL-colocalization class; HSD17B13 indel orientation ambiguous")
for c in [c for c in h87.columns if c not in ("ensembl_id", "symbol", "primary_evidence_class")]:
    u = ("normalized SuSiE PIP" if c == "coding_cs_max_pip" else "count" if c.startswith(("coding_cs_n_",)) else
         "GRCh38 position" if c.endswith("pos_grch38") else "text/state")
    G[c] = dict(unit=u, **h87_src)
meta(["enzyme_sex_effect_state"], source="H72 sex-heterogeneous enzyme effects", run_dir=str(P["H72"][0]), result=str(P["H72"][1]),
     source_file=str(P["H72"][2]), unit="state", release_state=RS, blocked_on="corrected 50-study COLOC rerun (rebuild)",
     caveat="UK Biobank participants (ALT 184,990 women, 159,146 men); BH over 281 SNP tests per scale; 4 scale_dependent genes (SAR1A, TYSND1, HKDC1, DYRK2) arise from one sex's estimate near zero; says nothing about sex differences in MASLD; not_in_tested_family = class not anchored or concordant")
f3_src = dict(source="Fig 3E stage contrasts (774-unit refit, stage extension)", run_dir=str(FIG3E),
              result="docs/technical/RESOURCE_REVIEW_AND_HYPOTHESES_2026-10-15.md (section 6.1 gap 1; G2-01, G2-02)",
              source_file=str(FIG3E / "stage_extension_all_gene_results.tsv"), release_state=RS, blocked_on="H51;H82",
              caveat=(f"387 units from 3 cohorts {fib_by_ds}; GSE213621 contributes 0 units; F0 = 103 (14 source controls + 89 F0 disease); "
                      "BH within contrast over 24,115 genes, no FC floor; cross-sectional stage-associated remodeling; G2-01: of 1,338 F1-vs-F0 calls 14 remain without GSE130970 "
                      "(929 without GSE135251, 651 without GSE162694), GNMT and MAT1A lose F0-to-F1 support without GSE130970; G2-02: F1 calls fall to 104-188 under joint library-complexity covariates"))
for k in ["F1", "F2", "F3", "F4"]:
    meta([f"fig3e_{k}_vs_F0_log2fc"], unit=f"log2 fold change, {k} vs F0", **f3_src)
    meta([f"fig3e_{k}_vs_F0_bh_q"], unit="BH q over 24,115 genes", **f3_src)
    meta([f"fig3e_{k}_vs_F0_call"], unit="state at BH q<0.05", **f3_src)
lin_src = dict(source="Donor-level lineage specificity (adopted 2026-10-14 corrected class comparison)", run_dir=str(LIN),
               result="docs/RESULTS.md (class comparison, 281 pairs); " + str(LIN / "VERDICT.md"),
               source_file=str(LIN / "results/gene_classes.tsv.gz"), release_state=RS, blocked_on="",
               caveat=("runs collapsed to donors (Hepatocytes 70, Fibroblasts 56, Cholangiocytes 51, Endothelial 66, Macrophages 35, T cells 39); "
                       "joined by symbol in the source; values for the 14,704 stable-ID universe only; the class test that used it is post hoc and opposite to the prespecified direction; "
                       "no per-gene lineage disease contrast exists in an adopted release"))
meta(["lineage_hep_specificity_log2cpm"], unit="mean log2(CPM+1) in hepatocytes minus mean over 5 non-parenchymal lineages", **lin_src)
meta(["lineage_measured_in_sc"], unit="boolean (detected in any lineage pseudobulk)", **lin_src)

dict_rows = []
for c in W.columns:
    m = G.get(c)
    assert m is not None, f"no dictionary entry for {c}"
    s = W[c]
    if pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s) and s.nunique(dropna=True) > 12:
        allowed = f"numeric [{s.min():.6g}, {s.max():.6g}]; n non-missing {int(s.notna().sum())}"
    else:
        v = s.astype(object).where(s.notna(), "<missing>").astype(str).value_counts()
        allowed = "; ".join(f"{k} ({n})" for k, n in v.head(15).items()) + (f"; ... {len(v) - 15} more" if len(v) > 15 else "")
    dict_rows.append({"table": "catalog_r9_candidate_gene_table", "column": c, **{k: m[k] for k in
                      ["source", "run_dir", "result", "source_file", "unit"]}, "allowed_values": allowed,
                      "release_state": m["release_state"], "blocked_on": m["blocked_on"], "caveat": m["caveat"]})
SIDE_META = {
    "sidecar_fig3e_stage_contrasts_long.tsv.gz": dict(f3_src, unit="per contrast: log2FC, SE, CI, t, p, BH q (FDR) as in source", caveat=f3_src["caveat"] + "; 7 contrasts in the installed Fig 3E source table (4 fibrosis on the panel, 3 NAS)"),
    "sidecar_fig5b_liver_protein.tsv": dict(source="Fig 5B liver protein triage", run_dir=str(FIG5B.parent), result="figures/main/fig5_molecular_context (installed from bulk-context-774-20261014T140904Z_22236311)",
                                            source_file=str(FIG5B), unit="log2FC and BH q (bulk and protein)", release_state=RS, blocked_on="",
                                            caveat="figure source values, selection-conditioned (fixed 3,054-gene prioritized universe); not a Catalog evidence call; proteomics gene_evidence_authorized false"),
    "sidecar_h72_enzyme_sex_effect_detail.tsv": dict(G["enzyme_sex_effect_state"], unit="per-gene test detail (exact copy)"),
    "sidecar_h74_nearest_gene_loci.tsv": dict(source="H74 nearest-gene contradiction", run_dir=str(P["H74"][0]), result=str(P["H74"][1]), source_file=str(P["H74"][2]),
                                              unit="per-locus table (exact copy)", release_state=RS, blocked_on="user adoption of H74; corrected 50-study COLOC rerun",
                                              caveat="a per-locus Catalog flag needs user adoption of the H74 result; not joined to the gene table"),
    "sidecar_h78_followup_annotation.tsv": dict(source="H78 gnomAD v4.1 / DepMap 24Q4 / HPA", run_dir=str(P["H78"][0]), result=str(P["H78"][1]), source_file=str(P["H78"][2]),
                                                unit="per-gene annotation (exact copy)", release_state=RS, blocked_on="",
                                                caveat="destination is a supplementary table; joining to the gene index is a user decision; annotations do not enter any class or route"),
    "sidecar_h71_competitor_list_states.tsv": dict(source="H71 competitor list states (Part A)", run_dir=str(P["H71"][0]), result=str(P["H71"][1]), source_file=str(P["H71"][2]),
                                                   unit="per list gene (copy plus bulk_q_full_precision_r9)", release_state=RS, blocked_on="corrected 50-study COLOC rerun (rebuild)",
                                                   caveat=f"verdict null for both lists; 'orthogonal' wording not supported; H71 bulk_q prints 0 for {n71_zero} rows (parser precision), full-precision q added"),
    "sidecar_h86_protein_coverage_by_class.tsv": dict(source="H86 protein measurability by class", run_dir=str(P["H86"][0]), result=str(P["H86"][1]),
                                                      source_file=str(P["H86"][2]) + " ; " + str(H86_D7), unit="counts and percent per class",
                                                      release_state=RS, blocked_on="corrected 50-study COLOC rerun (recount)",
                                                      caveat="quantified denominators and deposit detection only; significance columns excluded (proteomics program_context_only)"),
    "sidecar_ids_without_catalog_record.tsv": dict(source="F6 R07-05 IDs without a Catalog record", run_dir=str(F6P), result=str(F6_JSON) + " (R07-05 and checker)",
                                                   source_file=f"{BULK / 'coefficient_table.tsv.gz'} ; {COLOC / 'gene_level_coloc_tier12.csv'} ; {CAT / 'passport_program_membership.parquet'} ; {QUAR}",
                                                   unit="one row per source ID", release_state=RS, blocked_on="",
                                                   caveat="adding records is a Catalog rebuild, not done; quarantine reason matched on any listed candidate ID"),
    "sidecar_program_context_r9.tsv": dict(source="F6 R15-09 and H59 destination", run_dir=str(CAT), result=f"{F6_JSON} (R15-09); {P['H59'][2]}:418-420; {P['H59'][1]}",
                                           source_file=str(CAT / "passport_program_context.parquet"), unit="per program-context row", release_state=RS, blocked_on="",
                                           caveat="117-program stage rows not added (H59 columns go to the Fig 4 source table)"),
    "catalog_rules_r9.tsv": dict(source="F6 R07-02/R07-04/R07-06/R07-07/R07-11 rule text", run_dir=str(F6P), result=str(F6_JSON), source_file="", unit="text",
                                 release_state=RS, blocked_on="", caveat=""),
}
for name, df in SIDE.items():
    m = SIDE_META[name]
    for c in df.columns:
        s = df[c]
        if pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s) and s.nunique(dropna=True) > 12:
            allowed = f"numeric [{s.min():.6g}, {s.max():.6g}]; n non-missing {int(s.notna().sum())}"
        else:
            v = s.astype(object).where(s.notna(), "<missing>").astype(str).value_counts()
            allowed = "; ".join(f"{k[:60]} ({n})" for k, n in v.head(10).items()) + (f"; ... {len(v) - 10} more" if len(v) > 10 else "")
        dict_rows.append({"table": name, "column": c, "source": m["source"], "run_dir": m["run_dir"], "result": m["result"],
                          "source_file": m["source_file"], "unit": m["unit"], "allowed_values": allowed,
                          "release_state": m["release_state"], "blocked_on": m["blocked_on"], "caveat": m["caveat"]})
coldict = pd.DataFrame(dict_rows)

blocked = pd.DataFrame([
    ("rna_* effect columns (text using them)", "H79;H81", "H79 ancestry robustness of the 1,854 calls and H81 per-decile detection limit still running (jobs 6190, 6210)"),
    ("fig3e_* columns and sidecar_fig3e_stage_contrasts_long", "H51;H82", "H51 GSE213621 replication running (job 5104); RESULT_H82.md (untestable, S5) appeared 2026-10-15T17:22Z but plan_results.json still lists H82 as running"),
    ("gen_coloc_*, gen_r8_displayed_pp4_method, gen_untestable_reason_r9, class_display_r9, untested_*", "corrected 50-study COLOC rerun", "rebuild after adoption; memberships not frozen"),
    ("caqtl_* (H53)", "corrected 50-study COLOC rerun", "built on COLOC r2 supported set"),
    ("enzyme_sex_effect_state, sidecar_h72", "corrected 50-study COLOC rerun", "H72 result says rebuild on the adopted rerun"),
    ("coding_cs_* text (H87)", "corrected 50-study COLOC rerun", "proposed text waits for the rerun"),
    ("sidecar_h71, sidecar_h74, sidecar_h86 counts, H63 per-class counts", "corrected 50-study COLOC rerun", "class-dependent counts"),
    ("sidecar_h74 as a gene column", "user adoption of H74", "H74 destination"),
], columns=["columns_or_table", "blocked_on", "reason"])

# ================================================================ write outputs
W.to_parquet(OUT / "catalog_r9_candidate_gene_table.parquet", index=False)
W.to_csv(OUT / "catalog_r9_candidate_gene_table.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
rt = pd.read_parquet(OUT / "catalog_r9_candidate_gene_table.parquet")
check("parquet round trip: r8 columns equal", rt[g8.columns].equals(g8), True)
coldict.to_csv(OUT / "column_dictionary.tsv", sep="\t", index=False)
blocked.to_csv(OUT / "blocked_on.tsv", sep="\t", index=False)
UNIV.to_csv(OUT / "universe_counts.tsv", sep="\t", index=False)
summary = {
    "out_dir": str(OUT), "seed": SEED,
    "r9_gene_table_rows": len(W), "r9_gene_table_columns": W.shape[1], "r8_columns": list(g8.columns),
    "new_gene_columns": [c for c in W.columns if c not in g8.columns],
    "class_counts": vc(W.primary_evidence_class),
    "r8_molecular_object_assay_applicability": mo_vals,
    "fig3e_fibrosis_units_by_cohort": fib_by_ds,
    "lineage_ids_not_in_catalog": lin_not_in_cat,
    "h71_bulk_q_zero_rows": n71_zero,
    "ids_without_record_by_reason": {f"{a}|{b}": int(v) for (a, b), v in ab_counts.items()},
    "sidecars": {k: [int(v.shape[0]), int(v.shape[1])] for k, v in SIDE.items()},
}
json.dump(summary, open(OUT / "build_summary.json", "w"), indent=1, default=str)
inv = pd.DataFrame(INV)
inv.to_csv(OUT / "invariants.tsv", sep="\t", index=False)
pd.DataFrame([{"path": str(p), "sha256": sha256(p), "bytes": p.stat().st_size} for p in sorted(INPUTS)]).to_csv(
    OUT / "inputs_sha256.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as fh:
    fh.write(subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True).stdout)

hard_fail = inv[(~inv.passed) & inv.hard]
if len(hard_fail):
    print("HARD FAILURES:\n" + hard_fail.to_string())
    sys.exit(1)


def nrows(p):
    if p.suffix == ".parquet":
        return len(pd.read_parquet(p))
    if p.name.endswith((".tsv", ".tsv.gz")):
        return len(pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False))
    return ""


man_rows = []
for p in sorted(OUT.iterdir()):
    if p.is_file() and p.name != "manifest.tsv":
        man_rows.append({"file": p.name, "rows_excluding_header": nrows(p), "bytes": p.stat().st_size, "sha256": sha256(p)})
pd.DataFrame(man_rows).to_csv(OUT / "manifest.tsv", sep="\t", index=False)
print("DONE", OUT)
