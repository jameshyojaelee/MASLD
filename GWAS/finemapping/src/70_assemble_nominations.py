#!/usr/bin/env python
"""
70_assemble_nominations.py  --  Tier-B seqfunc ASSEMBLER
========================================================
ADDITIVE (new src/70). Combines the 4 direction-INDEPENDENT Tier-B axes into a
single eQTL-absent NOMINATION table with confidence tiers.

CONTEXT (shared): the eQTL-DIRECTION gate FAILED for BOTH Borzoi (on-recipe
logSED auROC 0.559, p=0.082) AND AlphaGenome (auROC 0.561, p=0.066). So this
layer is STRICTLY direction-independent: it produces MAGNITUDE + MECHANISM +
WHICH-GENE nominations for the fine-mapped liver signals whose effector gene
does NOT colocalize with an eQTL. Borzoi |logSED| MAGNITUDE is significant
(Spearman 0.39 vs |beta|, p=3.6e-8) -> magnitude nomination is legitimate.

NOMINATION / mechanism-class ONLY. NEVER wired into the convergence atlas or any
scored channel (circularity guardrail vs COLOC + epigenomic). NOT a direction
claim. Needs orthogonal wet-lab / MPRA validation.

The 4 axes (whichever succeeded; degrades gracefully if an axis file is absent):
  Axis-1 Coding hardening  coding_hardening.tsv        (ESM1b + popEVE + AlphaMissense; CODING subset)
  Axis-2 PoPS which-gene   pops_nominations.tsv        (orthogonal polygenic which-gene)
  Axis-3 ChromBPNet ATAC   chrombpnet_accessibility.tsv(HepG2 accessibility + TF-motif footprint)
  Axis-4 Borzoi magnitude  borzoi_magnitude_nominations.tsv (sequence-magnitude which-gene, PRIMARY spine)

Output
  eqtl_absent_nominations.tsv          -- one row per eQTL-absent fine-mapped locus/gene
  eqtl_absent_nominations.README.md    -- honest summary + top convergent nominations

Row grain: the 81 eQTL-absent regulatory fine-mapped leads (Borzoi magnitude spine)
+ the CODING subset (coding_hardening) appended as var_class=coding, tier=CODING.
"""
from __future__ import annotations
import os, sys, csv, json
import numpy as np
import pandas as pd
try:
    from scipy.stats import binomtest as _binomtest
except Exception:
    _binomtest = None

S = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc"

P_BORZOI = os.path.join(S, "borzoi_magnitude_nominations.tsv")
P_BORZOI_LEADS = os.path.join(S, "borzoi_magnitude_leads.tsv")   # spine fallback (no scores)
P_POPS   = os.path.join(S, "pops_nominations.tsv")
P_CBP    = os.path.join(S, "chrombpnet_accessibility.tsv")
P_CODING = os.path.join(S, "coding_hardening.tsv")
P_CONC   = os.path.join(S, "ag_borzoi_pergene_concordance.tsv")
P_SUBS   = os.path.join(S, "variant_substrate_hg38.tsv")
P_AGMECH = os.path.join(S, "ag_mechanism_profile.tsv")   # Axis-5 (src/72, additive)

# The 7 convergent nominations whose which-gene rested on nearest-TSS (barely
# orthogonal). AlphaGenome 3D-contact linking (src/72) gives each an independent
# which-gene vote that corroborates or overturns the nearest-TSS call.
NEAREST_TSS_DEPENDENT = {"P2RX7", "FCGRT", "TM4SF4", "KLHL8", "CD276", "MTARC1", "ABCB11"}

OUT_TSV  = os.path.join(S, "eqtl_absent_nominations.tsv")
OUT_MD   = os.path.join(S, "eqtl_absent_nominations.README.md")

# ---- accessibility-support thresholds ------------------------------------
CBP_PVAL_SIG = 0.05        # cbp_abs_logfc empirical p (from variant-scorer null)


def log(m): print(f"[70_assemble] {m}", flush=True)


def norm_var(v):
    """Normalise an hg19 variant id 'chr:pos:ref:alt' -> 'pos-less strand-agnostic' key
    'chr:pos' (used for fallback position matching). Returns (full, chrpos)."""
    if v is None or (isinstance(v, float) and np.isnan(v)) or str(v).strip() in ("", "-", "nan"):
        return None, None
    v = str(v).strip()
    parts = v.replace("chr", "").split(":")
    if len(parts) >= 2:
        return v.replace("chr", ""), f"{parts[0]}:{parts[1]}"
    return v, None


def read_tsv(path):
    return pd.read_table(path, dtype=str, keep_default_na=False)


# ==========================================================================
# 1. LOAD AXES (with graceful degradation)
# ==========================================================================
axis_status = {}   # axis -> "green"/"red"

# ---- Axis-4 Borzoi magnitude (SPINE) -------------------------------------
if os.path.exists(P_BORZOI) and os.path.getsize(P_BORZOI) > 0:
    bz = read_tsv(P_BORZOI)
    axis_status["borzoi_magnitude"] = "green"
    log(f"Borzoi magnitude: {len(bz)} leads")
elif os.path.exists(P_BORZOI_LEADS):
    bz = read_tsv(P_BORZOI_LEADS)
    bz = bz.rename(columns={"gene_substrate": "substrate_gene",
                            "ensembl_substrate": "substrate_gene_ensembl"})
    for c in ["best_gene", "best_gene_ensembl", "best_gene_biotype", "borzoi_abs_logsed",
              "magnitude_percentile", "exceeds_p80", "fold_agreement", "hg38_ref", "hg38_alt"]:
        if c not in bz.columns:
            bz[c] = ""
    axis_status["borzoi_magnitude"] = "red"
    log(f"Borzoi magnitude RED (scores absent) -- using {len(bz)} leads as bare spine")
else:
    sys.exit("FATAL: no Borzoi spine (neither nominations nor leads). Cannot assemble.")

# canonical spine keys
bz["credible_variant"] = bz["variant_id_hg19"].astype(str).str.replace("chr", "", regex=False)
bz["chrpos"] = bz["credible_variant"].apply(lambda v: norm_var(v)[1])

# ---- Axis-2 PoPS ---------------------------------------------------------
if os.path.exists(P_POPS) and os.path.getsize(P_POPS) > 0:
    pops = read_tsv(P_POPS)
    axis_status["pops"] = "green"
    log(f"PoPS: {len(pops)} loci")
    pops["credible_variant"] = pops["credible_variant"].astype(str).str.replace("chr", "", regex=False)
    pops["chrpos"] = pops["credible_variant"].apply(lambda v: norm_var(v)[1])
else:
    pops = pd.DataFrame(columns=["credible_variant", "chrpos", "pops_top_gene", "pops_top_ensgid",
                                 "nearest_gene", "nearest_ensgid", "agree", "pops_score"])
    axis_status["pops"] = "red"
    log("PoPS RED (pops_nominations.tsv absent)")

# ---- Axis-3 ChromBPNet ----------------------------------------------------
if os.path.exists(P_CBP) and os.path.getsize(P_CBP) > 0:
    cbp = read_tsv(P_CBP)
    axis_status["chrombpnet"] = "green"
    log(f"ChromBPNet: {len(cbp)} variants")
    cbp["variant_id_hg19"] = cbp["variant_id_hg19"].astype(str).str.replace("chr", "", regex=False)
    for c in ["cbp_abs_logfc", "cbp_jsd", "cbp_abs_logfc_pval", "pos_hg38"]:
        if c in cbp.columns:
            cbp[c] = pd.to_numeric(cbp[c], errors="coerce")
    if "in_hepg2_peak" in cbp.columns:
        cbp["in_hepg2_peak_b"] = cbp["in_hepg2_peak"].astype(str).str.upper().isin(("TRUE", "1", "T"))
    else:
        cbp["in_hepg2_peak_b"] = False
    if "disrupted_tf_motif" not in cbp.columns:
        cbp["disrupted_tf_motif"] = "-"
    # per-variant disruptive flag: in a HepG2 ATAC peak AND disrupts a TF motif
    # AND (empirically significant abs_logfc if pval available)
    has_pval = "cbp_abs_logfc_pval" in cbp.columns and cbp["cbp_abs_logfc_pval"].notna().any()
    sig = (cbp["cbp_abs_logfc_pval"] < CBP_PVAL_SIG) if has_pval else True
    cbp["cbp_disruptive"] = (cbp["in_hepg2_peak_b"] &
                             (cbp["disrupted_tf_motif"].astype(str) != "-") & sig)
else:
    cbp = pd.DataFrame(columns=["variant_id_hg19", "chr", "pos_hg38", "cbp_abs_logfc",
                                "disrupted_tf_motif", "in_hepg2_peak_b", "cbp_disruptive"])
    axis_status["chrombpnet"] = "red"
    log("ChromBPNet RED (chrombpnet_accessibility.tsv absent)")

# ---- Axis-1 Coding hardening ---------------------------------------------
if os.path.exists(P_CODING) and os.path.getsize(P_CODING) > 0:
    cod = read_tsv(P_CODING)
    axis_status["coding"] = "green"
    log(f"Coding hardening: {len(cod)} coding variants / {cod['gene'].nunique()} genes")
else:
    cod = pd.DataFrame(columns=["gene", "variant", "esm1b_llr", "popeve_score",
                                "popeve_severity_rank", "am_pathogenicity", "am_class",
                                "esm_call", "popeve_call", "am_call"])
    axis_status["coding"] = "red"
    log("Coding hardening RED (coding_hardening.tsv absent)")

# ---- AG/Borzoi concordance (context only; on eQTL-PRESENT benchmark) ------
conc = read_tsv(P_CONC) if os.path.exists(P_CONC) else pd.DataFrame()

# ---- AG/Borzoi high-confidence directional tier (both_high_conf) ----------
# Top cross-model tier: Borzoi 4/4 folds agree on sign AND |AlphaGenome gnomAD
# quantile| >= 0.9 AND the two models agree on sign. On the eQTL-PRESENT Broadaway
# benchmark this tier's agreed model direction matched the MEASURED eQTL sign in
# 13/16 = 81% of genes (binom p = 0.011). We fold this in as a per-gene CONFIDENCE
# ANNOTATION (ag_borzoi_highconf) = cross-model directional RELIABILITY context.
# It is NOT a per-locus direction claim for this eQTL-absent set (both blanket
# direction gates FAILED); it flags nominated genes for which the two sequence
# models were independently high-confidence AND sign-agreeing.
def _truthy(s):
    return s.astype(str).str.upper().isin(("TRUE", "1", "T"))
hc_genes, hc_stat, hc_k, hc_n = set(), "", 0, 0
if len(conc) and "both_high_conf" in conc.columns:
    _hc = conc[_truthy(conc["both_high_conf"])].copy()
    hc_genes = {g for g in _hc.get("gene", pd.Series(dtype=str)).astype(str) if g and g != "nan"}
    if "measured_concordant" in _hc.columns and "eqtl_sign_if_measured" in _hc.columns:
        _meas = _hc[_hc["eqtl_sign_if_measured"].astype(str).str.strip().replace({"nan": ""}) != ""]
        hc_n = len(_meas)
        hc_k = int(_truthy(_meas["measured_concordant"]).sum())
        if hc_n:
            _p = (f" (binom p={_binomtest(hc_k, hc_n, 0.5, alternative='greater').pvalue:.3f})"
                  if _binomtest is not None else "")
            hc_stat = f"{hc_k}/{hc_n} = {100*hc_k/hc_n:.0f}%{_p}"
log(f"AG/Borzoi both_high_conf genes: {len(hc_genes)}; measured directional concordance {hc_stat or 'n/a'}")

# ---- substrate coloc lookup (per-row eQTL-absent provenance) --------------
coloc_status = {}   # variant_id_hg19 (chr-stripped) -> "eQTL-absent"/"eQTL-present"
if os.path.exists(P_SUBS):
    for r in csv.DictReader(open(P_SUBS), delimiter="\t"):
        v = str(r.get("variant_id_hg19", "")).replace("chr", "")
        coloc_status[v] = ("eQTL-present" if str(r.get("colocalizes", "")).upper() == "TRUE"
                           else "eQTL-absent")
def coloc_of(v):
    return coloc_status.get(str(v).replace("chr", ""), "unknown")


# ---- Axis-5 AlphaGenome multimodal mechanism (src/72; additive, graceful) ----
# Contact-map which-gene + mechanism-CLASS annotations. NOT a direction claim, NOT
# summed into any tier score -- folded in as per-lead ANNOTATIONS. For the 7
# nearest-TSS-dependent convergent loci the AG-contact call SUPERSEDES the weak
# nearest-TSS which-gene vote (corroborated / overturned, flagged).
agmech_by_var = {}
axis_status["ag_mechanism"] = "red"
if os.path.exists(P_AGMECH) and os.path.getsize(P_AGMECH) > 0:
    _agm = read_tsv(P_AGMECH)
    axis_status["ag_mechanism"] = "green"
    for _, r in _agm.iterrows():
        if str(r.get("var_class", "")) != "regulatory":
            continue
        v = str(r.get("variant_id_hg19", r.get("variant_id", ""))).replace("chr", "")
        agmech_by_var[v] = r
    log(f"Axis-5 AG mechanism: {len(agmech_by_var)} regulatory leads "
        f"(contact which-gene + mechanism-class)")
else:
    log("Axis-5 AG mechanism RED (ag_mechanism_profile.tsv absent) -- columns emitted empty")


# ==========================================================================
# 2. JOIN AXES ONTO THE BORZOI SPINE (per regulatory lead)
# ==========================================================================
# PoPS: exact credible_variant, then chrpos fallback
pops_by_var = {r["credible_variant"]: r for _, r in pops.iterrows()} if len(pops) else {}
pops_by_chrpos = {r["chrpos"]: r for _, r in pops.iterrows()} if len(pops) else {}
# ChromBPNet: index by hg19 variant + a per-chr position list for clump-window scan
cbp_by_var = {r["variant_id_hg19"]: r for _, r in cbp.iterrows()} if len(cbp) else {}
cbp_by_chr = {}
if len(cbp):
    for _, r in cbp.iterrows():
        cbp_by_chr.setdefault(str(r.get("chr", "")).replace("chr", ""), []).append(r)

pops_matched = 0
rows = []
for _, L in bz.iterrows():
    cv = L["credible_variant"]
    chrpos = L["chrpos"]
    best_gene = L.get("best_gene", "") or ""
    best_ens = L.get("best_gene_ensembl", "") or ""

    # ---- PoPS join ----
    pr = pops_by_var.get(cv)
    pops_match = "exact"
    if pr is None and chrpos is not None:
        pr = pops_by_chrpos.get(chrpos)
        pops_match = "chrpos" if pr is not None else "none"
    if pr is None:
        pops_match = "none"
    else:
        pops_matched += 1
    pops_top_gene = (pr["pops_top_gene"] if pr is not None else "")
    pops_top_ens = (pr.get("pops_top_ensgid", "") if pr is not None else "")
    nearest_gene = (pr["nearest_gene"] if pr is not None else "")
    nearest_ens = (pr.get("nearest_ensgid", "") if pr is not None else "")

    # which-gene agreement (ensembl-primary, name fallback)
    def agree(g1, e1, g2, e2):
        if e1 and e2 and str(e1) != "" and str(e2) != "":
            return str(e1).split(".")[0] == str(e2).split(".")[0]
        if g1 and g2:
            return str(g1) == str(g2)
        return False
    agree_pops = bool(best_gene) and agree(best_gene, best_ens, pops_top_gene, pops_top_ens)
    agree_nearest = bool(best_gene) and agree(best_gene, best_ens, nearest_gene, nearest_ens)
    whichgene_agree = agree_pops or agree_nearest

    # ---- ChromBPNet: lead-variant (strict) + clump-window (context) ----
    lead_cbp = cbp_by_var.get(cv)
    cbp_abs = float(lead_cbp["cbp_abs_logfc"]) if lead_cbp is not None and pd.notna(lead_cbp.get("cbp_abs_logfc")) else np.nan
    disr_tf = (lead_cbp["disrupted_tf_motif"] if lead_cbp is not None else "-") or "-"
    in_peak = bool(lead_cbp["in_hepg2_peak_b"]) if lead_cbp is not None else False
    lead_disruptive = bool(lead_cbp["cbp_disruptive"]) if lead_cbp is not None else False
    # clump-window max disruptive (context): any cbp variant within clump span of lead
    locus_cbp_max = np.nan
    locus_disruptive = False
    try:
        span = int(float(L.get("clump_span_bp", 0) or 0))
    except Exception:
        span = 0
    try:
        lead_pos = int(float(L.get("pos_hg38", 0) or 0))
    except Exception:
        lead_pos = 0
    chrom = str(L.get("chr", "")).replace("chr", "")
    if chrom in cbp_by_chr and lead_pos and span >= 0:
        lo, hi = lead_pos - span, lead_pos + span
        for r in cbp_by_chr[chrom]:
            p = r.get("pos_hg38")
            if pd.isna(p):
                continue
            if lo <= float(p) <= hi:
                a = r.get("cbp_abs_logfc")
                if pd.notna(a) and (np.isnan(locus_cbp_max) or float(a) > locus_cbp_max):
                    locus_cbp_max = float(a)
                if bool(r.get("cbp_disruptive")):
                    locus_disruptive = True

    accessibility_support = lead_disruptive   # STRICT: lead credible variant disruptive

    # ---- magnitude axis ----
    try:
        mag_pct = float(L.get("magnitude_percentile", "") or "nan")
    except Exception:
        mag_pct = np.nan
    exceeds_p80 = str(L.get("exceeds_p80", "")).upper() == "TRUE"
    high_magnitude = exceeds_p80 or (pd.notna(mag_pct) and mag_pct >= 80)
    borzoi_has_gene = bool(best_gene) and best_gene not in ("", "-")

    # ---- n axes + tier ----
    axes = {
        "borzoi_magnitude": high_magnitude and borzoi_has_gene,
        "whichgene_agree": whichgene_agree,
        "accessibility": accessibility_support,
    }
    n_axes = sum(bool(v) for v in axes.values())

    if borzoi_has_gene and whichgene_agree and (accessibility_support or high_magnitude):
        tier = "CONVERGENT-NOMINATION"
    elif n_axes >= 2:
        tier = "MULTI-AXIS"                 # 2 axes but not the which-gene-anchored convergent combo
    elif n_axes == 1:
        tier = "SINGLE-AXIS"
    else:
        tier = "NONE"

    # ---- Axis-5 AlphaGenome mechanism annotations (contact which-gene + class) ----
    agm = agmech_by_var.get(cv)
    ag_contact_gene = str(agm.get("ag_contact_gene", "")) if agm is not None else ""
    ag_mechanism_class = str(agm.get("mechanism_class", "")) if agm is not None else ""
    ag_element_class = str(agm.get("ag_element_class", "")) if agm is not None else ""
    ag_contact_agrees_abc = str(agm.get("ag_contact_agrees_abc", "")) if agm is not None else ""
    ag_contact_agrees_borzoi = str(agm.get("ag_contact_agrees_borzoi", "")) if agm is not None else ""
    ag_contact_overturns_nearest = str(agm.get("ag_contact_overturns_nearest", "")) if agm is not None else ""
    ag_tf_chip_disrupted = str(agm.get("ag_tf_chip_disrupted", "")) if agm is not None else ""
    # AlphaGenome 3D-contact relation to the nomination. IMPORTANT (honest finding,
    # src/72 README): the AG contact map is O/E-normalised, so the argmax is a long-range
    # LOOP-PARTNER detector, NOT an effector caller -- it disagrees with measured ABC
    # liver Hi-C on the one testable locus (KLHL8) and with Borzoi ~at chance. Contact is
    # therefore an EXPLORATORY annotation and does NOT supersede the nearest-TSS/PoPS
    # which-gene vote; the existing whichgene_agree/tier logic is UNCHANGED.
    contact_corrob = bool(ag_contact_gene) and bool(best_gene) and ag_contact_gene == best_gene
    if not ag_contact_gene:
        ag_contact_relation = ""
    elif contact_corrob:
        ag_contact_relation = "agrees_nomination"
    else:
        ag_contact_relation = "loop_partner_differs"

    rows.append({
        "locus": L.get("locus_id", ""),
        "credible_variant": cv,
        "nominated_gene": best_gene,
        "nominated_ensembl": best_ens,
        "var_class": "regulatory",
        "coloc_status": coloc_of(cv),
        "lead_pip": L.get("lead_pip", ""),
        "borzoi_abs_logsed": L.get("borzoi_abs_logsed", ""),
        "magnitude_percentile": (f"{mag_pct:.1f}" if pd.notna(mag_pct) else ""),
        "exceeds_p80": str(high_magnitude).upper(),
        "fold_agreement": L.get("fold_agreement", ""),
        "chrombpnet_abs_logfc": (f"{cbp_abs:.5f}" if pd.notna(cbp_abs) else ""),
        "disrupted_tf": disr_tf,
        "in_hepg2_peak": str(in_peak).upper(),
        "locus_cbp_max_abs_logfc": (f"{locus_cbp_max:.5f}" if pd.notna(locus_cbp_max) else ""),
        "accessibility_support": str(accessibility_support).upper(),
        "pops_top_gene": pops_top_gene,
        "agree_with_borzoi_gene": str(agree_pops).upper(),
        "nearest_gene": nearest_gene,
        "nearest_agrees_with_borzoi": str(agree_nearest).upper(),
        "whichgene_agree": str(whichgene_agree).upper(),
        "ag_borzoi_highconf": str(bool(best_gene) and best_gene in hc_genes).upper(),
        "ag_contact_gene": ag_contact_gene,
        "ag_contact_agrees_borzoi": ag_contact_agrees_borzoi,
        "ag_contact_agrees_abc": ag_contact_agrees_abc,
        "ag_contact_overturns_nearest": ag_contact_overturns_nearest,
        "ag_element_class": ag_element_class,
        "ag_tf_chip_disrupted": ag_tf_chip_disrupted,
        "mechanism_class": ag_mechanism_class,
        "ag_contact_relation": ag_contact_relation,
        "coding_severity": "",             # regulatory rows: coding N/A
        "coding_severity_rank": "",
        "substrate_gene": L.get("substrate_gene", ""),
        "pops_join": pops_match,
        "n_axes_supporting": n_axes,
        "confidence_tier": tier,
        "note": L.get("note", ""),
    })

log(f"PoPS matched to spine: {pops_matched}/{len(bz)} leads")


# ==========================================================================
# 3. CODING SUBSET (Axis-1) appended as var_class=coding, tier=CODING
# ==========================================================================
def coding_sev(r):
    """Composite severity label + numeric score (higher = more deleterious)."""
    esm = str(r.get("esm_call", "") or "")
    pop = str(r.get("popeve_call", "") or "")
    am = str(r.get("am_class", "") or "")
    score = 0.0
    esm_s = {"ESM_damaging": 2, "ESM_moderate": 1, "ESM_tolerant": 0}.get(esm, 0)
    am_s = {"likely_pathogenic": 2, "ambiguous": 1, "likely_benign": 0}.get(am, 0)
    pop_s = {"popEVE_severe": 2, "popEVE_moderate": 1, "popEVE_tolerant": 0}.get(pop, 0)
    # continuous corroboration
    try:
        amp = float(r.get("am_pathogenicity", "") or "nan")
    except Exception:
        amp = np.nan
    score = esm_s + am_s + pop_s + (amp if pd.notna(amp) else 0)
    if esm == "ESM_damaging" or am == "likely_pathogenic":
        label = "SEVERE"
    elif esm_s or am_s or pop_s:
        label = "MODERATE"
    else:
        label = "TOLERANT"
    detail = f"ESM={esm or 'NA'};popEVE={pop or 'NA'};AM={am or 'NA'}"
    return label, detail, score

# is this coding gene also a regulatory nominated gene? (cross-axis corroboration)
reg_nom_genes = {r["nominated_gene"] for r in rows if r["nominated_gene"] not in ("", "-")}

coding_rows = []
if axis_status["coding"] == "green":
    for _, r in cod.iterrows():
        label, detail, score = coding_sev(r)
        gene = r.get("gene", "")
        cvar = str(r.get("variant", "")).replace("chr", "")
        n_ax = 1 + (1 if gene in reg_nom_genes else 0)  # coding axis (+regulatory corroboration)
        coding_rows.append({
            "gene": gene, "variant": cvar, "label": label, "detail": detail,
            "score": score, "rank_field": r.get("popeve_severity_rank", ""),
            "corroborates_regulatory": gene in reg_nom_genes, "n_ax": n_ax,
            "var_hg38": r.get("variant_hg38", ""), "prot": r.get("protein_variant", ""),
        })
    # severity rank (1 = most severe)
    coding_rows.sort(key=lambda x: -x["score"])
    for i, cr in enumerate(coding_rows, 1):
        rows.append({
            "locus": cr["var_hg38"] or cr["variant"],
            "credible_variant": cr["variant"],
            "nominated_gene": cr["gene"],
            "nominated_ensembl": "",
            "var_class": "coding",
            "coloc_status": coloc_of(cr["variant"]),
            "lead_pip": "",
            "borzoi_abs_logsed": "", "magnitude_percentile": "", "exceeds_p80": "",
            "fold_agreement": "",
            "chrombpnet_abs_logfc": "", "disrupted_tf": "", "in_hepg2_peak": "",
            "locus_cbp_max_abs_logfc": "", "accessibility_support": "",
            "pops_top_gene": "", "agree_with_borzoi_gene": "",
            "nearest_gene": "", "nearest_agrees_with_borzoi": "", "whichgene_agree": "",
            "ag_borzoi_highconf": str(cr["gene"] in hc_genes).upper(),
            "ag_contact_gene": "", "ag_contact_agrees_borzoi": "",
            "ag_contact_agrees_abc": "", "ag_contact_overturns_nearest": "",
            "ag_element_class": "", "ag_tf_chip_disrupted": "",
            "mechanism_class": "protein-altering", "ag_contact_relation": "",
            "coding_severity": f"{cr['label']} ({cr['detail']}; {cr['prot']})",
            "coding_severity_rank": i,
            "substrate_gene": cr["gene"],
            "pops_join": "",
            "n_axes_supporting": cr["n_ax"],
            "confidence_tier": ("CODING+REG" if cr["corroborates_regulatory"] else "CODING"),
            "note": ("also_regulatory_borzoi_nom" if cr["corroborates_regulatory"] else ""),
        })

# ==========================================================================
# 4. WRITE
# ==========================================================================
COLS = ["locus", "credible_variant", "nominated_gene", "nominated_ensembl", "var_class",
        "coloc_status", "lead_pip", "borzoi_abs_logsed", "magnitude_percentile", "exceeds_p80",
        "fold_agreement", "chrombpnet_abs_logfc", "disrupted_tf", "in_hepg2_peak",
        "locus_cbp_max_abs_logfc", "accessibility_support", "pops_top_gene",
        "agree_with_borzoi_gene", "nearest_gene", "nearest_agrees_with_borzoi",
        "whichgene_agree", "ag_borzoi_highconf",
        "ag_contact_gene", "ag_contact_agrees_borzoi", "ag_contact_agrees_abc",
        "ag_contact_overturns_nearest", "ag_element_class", "ag_tf_chip_disrupted",
        "mechanism_class", "ag_contact_relation",
        "coding_severity", "coding_severity_rank", "substrate_gene",
        "pops_join", "n_axes_supporting", "confidence_tier", "note"]
out = pd.DataFrame(rows)[COLS]
# stable sort: convergent first, then multi/single, then by n_axes and magnitude
tier_order = {"CONVERGENT-NOMINATION": 0, "MULTI-AXIS": 1, "SINGLE-AXIS": 2,
              "CODING+REG": 2, "CODING": 3, "NONE": 4}
out["_t"] = out["confidence_tier"].map(lambda t: tier_order.get(t, 5))
out["_m"] = pd.to_numeric(out["magnitude_percentile"], errors="coerce").fillna(-1)
out = out.sort_values(["_t", "n_axes_supporting", "_m"], ascending=[True, False, False]).drop(columns=["_t", "_m"])
out.to_csv(OUT_TSV, sep="\t", index=False)
log(f"wrote {OUT_TSV}  ({len(out)} rows)")

# ---- tier counts ----
reg = out[out["var_class"] == "regulatory"]
tc = out["confidence_tier"].value_counts().to_dict()
log(f"tier counts: {tc}")

# ==========================================================================
# 5. README
# ==========================================================================
n_reg = len(reg)
conv = reg[reg["confidence_tier"] == "CONVERGENT-NOMINATION"]
multi = reg[reg["confidence_tier"] == "MULTI-AXIS"]
single = reg[reg["confidence_tier"] == "SINGLE-AXIS"]
none = reg[reg["confidence_tier"] == "NONE"]
coding_out = out[out["var_class"] == "coding"]
n_cod_absent = int((coding_out["coloc_status"] == "eQTL-absent").sum())
n_cod_present = int((coding_out["coloc_status"] == "eQTL-present").sum())

# per-axis n
n_bz_scored = int((pd.to_numeric(reg["magnitude_percentile"], errors="coerce").notna()).sum())
n_high_mag = int((reg["exceeds_p80"] == "TRUE").sum())
n_wg = int((reg["whichgene_agree"] == "TRUE").sum())
n_acc = int((reg["accessibility_support"] == "TRUE").sum())
n_pops_join = int((reg["pops_join"] == "exact").sum() + (reg["pops_join"] == "chrpos").sum())

# ---- Axis-5 AlphaGenome mechanism-class + contact loop-partner summary ----
ag5_block = ""
if axis_status.get("ag_mechanism") == "green" and "ag_contact_gene" in reg.columns:
    has_contact = reg[reg["ag_contact_gene"].astype(str).str.strip() != ""]
    n_contact = len(has_contact)
    n_overturn = int((reg["ag_contact_overturns_nearest"] == "TRUE").sum())
    abc_ov = reg[reg["ag_contact_agrees_abc"].isin(["TRUE", "FALSE"])]
    abc_agree = int((abc_ov["ag_contact_agrees_abc"] == "TRUE").sum())
    bz_ov = reg[reg["ag_contact_agrees_borzoi"].isin(["TRUE", "FALSE"])]
    bz_agree = int((bz_ov["ag_contact_agrees_borzoi"] == "TRUE").sum())
    mech_counts = reg[reg["mechanism_class"].astype(str).str.strip() != ""]["mechanism_class"].value_counts().to_dict()
    n_tf_corr = int((reg["ag_tf_chip_disrupted"] == "TRUE").sum())
    # the 7 nearest-TSS-dependent convergent loci: mechanism-class + contact loop-partner
    lines = []
    for g in sorted(NEAREST_TSS_DEPENDENT):
        sub = reg[reg["nominated_gene"] == g]
        if not len(sub):
            continue
        r = sub.iloc[0]
        lines.append(f"| {g} | {str(r['mechanism_class']) or '-'} | {str(r['ag_element_class']) or '-'} | "
                     f"{str(r['ag_tf_chip_disrupted']) or '-'} | {str(r['ag_contact_gene']) or '-'} | "
                     f"{str(r['ag_contact_relation']) or '-'} | {str(r['ag_contact_agrees_abc']) or 'n/a'} |")
    tbl = "\n".join(lines) if lines else "| _(none matched)_ | | | | | | |"
    ag5_block = (
        f"\n## Axis-5 AlphaGenome mechanism-class + contact loop-partner (src/72)\n"
        f"Full-head liver-masked fingerprint (MEAN calibrated quantiles) -> mechanism-CLASS, plus an "
        f"independent 3D-contact annotation. **Class ONLY -- NOT a direction claim, NOT summed into any "
        f"tier; the existing tiers are UNCHANGED.**\n"
        f"- Mechanism-class distribution (regulatory, corroboration-aware): {mech_counts}; "
        f"AG-TF-footprint corroborated by ChromBPNet/motifbreakR in **{n_tf_corr}** leads.\n"
        f"- **Contact which-gene = HONEST NEGATIVE (exploratory only).** AlphaGenome's contact map is "
        f"O/E-normalised, so the argmax is a long-range LOOP-PARTNER detector, NOT an effector caller: "
        f"it overturns nearest-TSS in **{n_overturn}/{n_contact}** leads toward mostly distal/non-coding "
        f"TSS, agrees with Borzoi only **{bz_agree}/{len(bz_ov)}** (~chance) and with MEASURED ABC liver "
        f"Hi-C **{abc_agree}/{len(abc_ov)}** (on KLHL8, the one testable lead, ABC calls the nearest gene "
        f"and AG-contact a 235 kb-distal gene). It PASSES the SORT1 positive control (recovers the true "
        f"123 kb-distal effector over nearest PSRC1) -- so contact works for genuine long-range loops but "
        f"is **not** a reliable blanket which-gene method and does **NOT** supersede the nearest-TSS/PoPS "
        f"call. Retained as the `ag_contact_gene` / `ag_contact_relation` exploratory loop-partner flag.\n\n"
        f"### The 7 nearest-TSS-dependent convergent loci -- AG mechanism annotation\n"
        f"| nominated_gene | mechanism_class | element | TF_corroborated | ag_contact_loop_partner | contact_relation | vs_ABC |\n"
        f"|---|---|---|---|---|---|---|\n{tbl}\n")

# AG/Borzoi magnitude concordance context
ag_ctx = ""
if len(conc):
    try:
        c2 = conc.copy()
        c2["ag"] = pd.to_numeric(c2["ag_logsed"], errors="coerce").abs()
        c2["bz"] = pd.to_numeric(c2["borzoi_ensemble_logsed"], errors="coerce").abs()
        d = c2.dropna(subset=["ag", "bz"])
        if len(d) > 3:
            rho = d["ag"].rank().corr(d["bz"].rank())
            ag_ctx = (f"- Cross-model magnitude corroboration (context, on the eQTL-PRESENT "
                      f"benchmark, n={len(d)}): AlphaGenome |logSED| vs Borzoi |logSED| "
                      f"Spearman rho = {rho:.2f}.\n")
    except Exception:
        pass

# ---- AG/Borzoi high-confidence directional tier fold-in (confidence annotation) ----
reg_hc = reg[reg["ag_borzoi_highconf"] == "TRUE"] if "ag_borzoi_highconf" in reg.columns else reg.iloc[0:0]
cod_hc = coding_out[coding_out["ag_borzoi_highconf"] == "TRUE"] if "ag_borzoi_highconf" in coding_out.columns else coding_out.iloc[0:0]
hc_overlap = sorted({g for g in list(reg_hc["nominated_gene"]) + list(cod_hc["nominated_gene"]) if g and g != "-"})
hc_block = (
    f"\n## AG/Borzoi high-confidence directional tier (cross-model reliability annotation)\n"
    f"- Top cross-model tier `both_high_conf` (Borzoi 4/4 folds sign-agree AND |AlphaGenome "
    f"gnomAD quantile| >= 0.9 AND the two models agree on sign): **{len(hc_genes)} genes**.\n"
    f"- On the eQTL-PRESENT Broadaway benchmark, this tier's agreed model direction matched the "
    f"MEASURED eQTL sign in **{hc_stat or 'n/a'}** -- the directional-reliability floor for any "
    f"sequence-model direction call. (The BLANKET eQTL-direction gate still FAILED: Borzoi auROC "
    f"0.559 / AlphaGenome 0.561; direction is defensible ONLY inside this high-confidence subset.)\n"
    f"- Folded into this table as the per-gene `ag_borzoi_highconf` column. It is cross-model "
    f"directional RELIABILITY context, NOT a per-locus direction claim for the eQTL-absent set.\n"
    f"- eQTL-absent nominations whose gene falls in this high-confidence tier: "
    f"{('**' + ', '.join(hc_overlap) + '**') if hc_overlap else '_none_ (expected: the tier is built on eQTL-PRESENT eGenes; the eQTL-absent set is colocalizes==FALSE by construction).'}\n"
)

def top_table(df, k=20):
    lines = ["| rank | locus | credible_variant | nominated_gene | mag_pctile | pops_top | nearest | disrupted_TF | n_axes |",
             "|---|---|---|---|---|---|---|---|---|"]
    for i, (_, r) in enumerate(df.head(k).iterrows(), 1):
        lines.append(f"| {i} | {r['locus']} | {r['credible_variant']} | {r['nominated_gene']} | "
                     f"{r['magnitude_percentile']} | {r['pops_top_gene']} | {r['nearest_gene']} | "
                     f"{r['disrupted_tf']} | {r['n_axes_supporting']} |")
    return "\n".join(lines)

def top_coding(df, k=15):
    lines = ["| rank | gene | variant | protein/severity | corroborates_regulatory |",
             "|---|---|---|---|---|"]
    for _, r in df.head(k).iterrows():
        lines.append(f"| {r['coding_severity_rank']} | {r['nominated_gene']} | {r['credible_variant']} | "
                     f"{r['coding_severity']} | {'yes' if r['confidence_tier']=='CODING+REG' else 'no'} |")
    return "\n".join(lines)

md = f"""# eQTL-absent seqfunc nominations (Tier-B assembler)

**Direction-INDEPENDENT nominations only.** The eQTL-DIRECTION gate FAILED for BOTH
sequence models on the coloc benchmark: Borzoi on-recipe logSED auROC = 0.559
(p = 0.082, FAIL) and AlphaGenome auROC = 0.561 (p = 0.066, FAIL). We therefore do
NOT claim allelic direction. What survives is **magnitude** (Borzoi |logSED| vs
|beta| Spearman = 0.39, p = 3.6e-8) plus **mechanism** (TF-motif footprint) and
**which-gene** (sequence magnitude + PoPS + nearest-TSS). These are **NOMINATIONS /
mechanism-class labels ONLY** for the ~84% of fine-mapped liver signals with NO
eQTL colocalization. They are explicitly **NOT** wired into the convergence atlas or
any scored channel (circularity guardrail with COLOC + epigenomic), and are **NOT a
direction claim**. Every nomination needs orthogonal wet-lab / MPRA validation.

## Axis inventory (green = usable, red = absent)
| axis | file | status |
|---|---|---|
| Axis-1 Coding hardening (ESM1b + popEVE + AlphaMissense) | coding_hardening.tsv | {axis_status['coding']} |
| Axis-2 PoPS which-gene (orthogonal polygenic) | pops_nominations.tsv | {axis_status['pops']} |
| Axis-3 ChromBPNet HepG2 accessibility + TF footprint | chrombpnet_accessibility.tsv | {axis_status['chrombpnet']} |
| Axis-4 Borzoi sequence-magnitude which-gene (SPINE) | borzoi_magnitude_nominations.tsv | {axis_status['borzoi_magnitude']} |
| Axis-5 AlphaGenome 3D-contact which-gene + mechanism-class | ag_mechanism_profile.tsv | {axis_status.get('ag_mechanism','red')} |

## Row grain & counts
- **{n_reg} eQTL-absent regulatory fine-mapped leads** (Borzoi magnitude spine; top-PIP
  credible variant per 500 kb clump of colocalizes==FALSE regulatory variants).
- **{len(coding_out)} coding effectors** (coding_hardening) appended as var_class=coding,
  tier=CODING, severity-ranked (separate; ESM1b / popEVE / AlphaMissense only). Of these,
  **{n_cod_absent} are eQTL-absent** (in-target) and {n_cod_present} also colocalize with an
  eQTL (`coloc_status=eQTL-present`) — the coding/protein-altering mechanism is valid either
  way, but the eQTL-present ones are strictly outside the eQTL-absent target set (flagged).
- All {n_reg} regulatory spine leads are eQTL-absent regulatory (colocalizes==FALSE), in-target.

## Per-axis firing (regulatory leads, n={n_reg})
- Borzoi magnitude scored: {n_bz_scored}/{n_reg}; high-magnitude (>= p80 of the eQTL-present
  reference |logSED| distribution): **{n_high_mag}**.
- PoPS matched to spine: {n_pops_join}/{n_reg}; which-gene agreement (Borzoi gene == PoPS-top
  OR == nearest-TSS): **{n_wg}**.
- ChromBPNet accessibility-supported (lead credible variant in a HepG2 ATAC peak AND
  disrupts a TF motif AND empirically significant): **{n_acc}**.
{ag_ctx}{hc_block}{ag5_block}
## Confidence tiers (definitions)
- **CONVERGENT-NOMINATION** ({len(conv)}): Borzoi nominates a sequence-magnitude gene AND an
  orthogonal which-gene axis (PoPS or nearest-TSS) agrees AND (accessibility-supported OR
  high magnitude percentile). The strongest, still nomination-only.
- **MULTI-AXIS** ({len(multi)}): 2 axes fire but NOT the which-gene-anchored convergent combo
  (e.g. high magnitude + accessibility, but PoPS/nearest point to a different gene).
- **SINGLE-AXIS** ({len(single)}): exactly one axis fires.
- **NONE** ({len(none)}): Borzoi scored the locus but no axis clears threshold.
- **CODING / CODING+REG** ({len(coding_out)}): coding subset, severity-ranked; CODING+REG =
  the coding gene is ALSO an independent regulatory Borzoi nomination.

## Tier counts
{json.dumps(tc, indent=0)}

## Top convergent nominations (regulatory)
{top_table(conv, 20) if len(conv) else '_none reached CONVERGENT-NOMINATION_'}

## Top multi-axis nominations (regulatory)
{top_table(multi, 15) if len(multi) else '_none_'}

## Coding subset (severity-ranked, top)
{top_coding(coding_out, 15) if len(coding_out) else '_coding axis red_'}

## Honest caveats
- **Nomination-only.** No direction; magnitude/mechanism/which-gene only. The direction
  gate failed, so a nominated gene is NOT asserted to be up/down-regulated by the allele.
- The AG/Borzoi concordance table is measured on the **eQTL-PRESENT** coloc benchmark, not
  on this eQTL-absent target set; it is cross-model magnitude context, not per-lead evidence.
- Accessibility_support uses the **lead credible variant** (conservative, non-cherry-picked);
  a clump-window max (`locus_cbp_max_abs_logfc`) is provided as context but does NOT set the flag.
- PoPS and Borzoi clump the same substrate independently; {n_pops_join}/{n_reg} leads matched
  a PoPS locus (exact variant or chrpos fallback). Unmatched leads carry empty PoPS fields.
- Requires orthogonal wet-lab / MPRA (allelic reporter, base-editing, or CRISPRi at the
  nominated gene) to move from nomination to causal effector.

_Generated by src/70_assemble_nominations.py (additive; not wired into 46d/78/27a)._
"""
with open(OUT_MD, "w") as fh:
    fh.write(md)
log(f"wrote {OUT_MD}")

# ---- machine-readable tier summary for the caller ----
summary = {
    "axis_status": axis_status,
    "n_regulatory_leads": n_reg,
    "n_coding": len(coding_out),
    "tier_counts": tc,
    "per_axis": {"borzoi_high_magnitude": n_high_mag, "whichgene_agree": n_wg,
                 "accessibility_support": n_acc, "pops_matched": n_pops_join,
                 "borzoi_scored": n_bz_scored},
}
print("SUMMARY_JSON " + json.dumps(summary))
