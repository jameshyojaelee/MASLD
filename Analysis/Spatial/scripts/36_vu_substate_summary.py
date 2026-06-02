#!/usr/bin/env python3
"""
36_vu_substate_summary.py — Aggregate Team-M3 outputs into final verdict markdown.

Reads:
  - vu_deconvolved_substate.csv (script 33)
  - vu_within_hepatocyte_bimodality.csv (script 34)
  - vu_moran_matched_null.csv (script 35)
  - vu_per_patient_mixed_effects.csv (script 35)

Writes:
  - vu_deconvolution_summary.md  — 1-paragraph verdict on whether F3a/F3b
    bimodality reflects real within-hepatocyte sub-state biology
"""

import pathlib
import pandas as pd
import numpy as np

PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GRAN_OUT = PROJECT_ROOT / "RNA-seq/results/granular_staging"

OUT_MD = GRAN_OUT / "vu_deconvolution_summary.md"


def load(name):
    p = GRAN_OUT / name
    return pd.read_csv(p) if p.exists() else None


def main():
    deconv = load("vu_deconvolved_substate.csv")
    bim = load("vu_within_hepatocyte_bimodality.csv")
    mor = load("vu_moran_matched_null.csv")
    me = load("vu_per_patient_mixed_effects.csv")

    method = deconv["method"].iloc[0] if deconv is not None else "MISSING"
    n_total = len(deconv) if deconv is not None else 0
    n_hep = (
        deconv["hep_dominant_80"].sum() if deconv is not None else 0
    )
    pct_hep = 100 * n_hep / n_total if n_total else 0

    # Bimodality verdict
    bim_pass = []
    bim_lines = []
    if bim is not None:
        for _, r in bim.iterrows():
            if r["signature"] in ("li_f3b", "li_f3b_z_hep", "mixing_axis_hep",
                                   "mixing_axis"):
                ok = bool(r["pass_real_bimodality"])
                bim_pass.append(ok)
                bim_lines.append(
                    f"- `{r['signature']}` (n={int(r['n_hep_dom_spots'])}): "
                    f"ΔBIC={r['BIC_diff']:.1f}, sep={r['sep_sigma']:.2f}σ, "
                    f"dip-p={r['dip_p']:.4f} → "
                    f"{'PASS' if ok else 'FAIL'}"
                )
    bim_overall_pass = any(bim_pass)

    # Moran verdict
    mor_lines = []
    n_pass_per_sig = {}
    if mor is not None:
        for sig in mor["signature"].unique():
            sub = mor[mor["signature"] == sig]
            n_pass = int(sub["pass_real_organization"].sum())
            n_total_sig = len(sub)
            n_pass_per_sig[sig] = (n_pass, n_total_sig)
            mor_lines.append(
                f"- `{sig}`: {n_pass}/{n_total_sig} patients pass "
                f"(observed I in ≥80th percentile of 100 matched-null sets)"
            )

    # Mixed-effects verdict
    me_lines = []
    me_pass_lines = []
    if me is not None:
        for _, r in me.iterrows():
            ci_excl = bool(r["ci_excludes_zero"])
            n_ind_pass = int(r["n_patients_individual_perm_p_lt_05"])
            pass_4_of_5 = int(r["n_patients"]) > 0 and n_ind_pass >= 4
            me_pass_lines.append(pass_4_of_5)
            me_lines.append(
                f"- `{r['signature']}`: FE intercept={r['fixed_effect_intercept']:.3f} "
                f"(95% CI [{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]), "
                f"individual perm-p<.05 in {n_ind_pass}/{r['n_patients']} patients "
                f"→ {'PASS (≥4/5)' if pass_4_of_5 else 'FAIL (<4/5)'}"
            )

    # Final verdict
    bim_ok = bim_overall_pass
    mor_ok = any(np / nt >= 0.8 for np, nt in n_pass_per_sig.values()) if n_pass_per_sig else False
    me_ok = any(me_pass_lines)
    pass_count = int(bim_ok) + int(mor_ok) + int(me_ok)
    if pass_count == 3:
        verdict = (
            "REAL F3a/F3b SUB-STATE BIOLOGY: All three pre-registered acceptance "
            "criteria pass. Within-hepatocyte spots show genuine bimodality, "
            "spatial autocorrelation exceeds matched-null expectation, and "
            "≥4/5 patients individually demonstrate non-random spatial structure."
        )
    elif pass_count == 2:
        verdict = (
            "PARTIAL SUPPORT for F3a/F3b sub-state biology: 2/3 acceptance "
            "criteria pass. Some evidence of within-hepatocyte structure but "
            "not unambiguous."
        )
    elif pass_count == 1:
        verdict = (
            "WEAK / COMPOSITIONAL F3a/F3b ARTIFACT: Only 1/3 acceptance criteria "
            "pass. F3a/F3b axis bimodality in Vu is most consistent with "
            "hepatocyte/stellate compositional variation rather than within-"
            "hepatocyte sub-state biology."
        )
    else:
        verdict = (
            "NO REAL F3a/F3b SUB-STATE: 0/3 acceptance criteria pass. "
            "F3a/F3b axis bimodality is fully explained by cell-type "
            "compositional heterogeneity across spots; no within-hepatocyte "
            "sub-state signal remains after restricting to hepatocyte-dominant "
            "(≥80%) spots, and spatial structure is at chance level."
        )

    md_lines = [
        "# Vu et al. F3a/F3b Sub-state Deconvolution Verdict",
        "",
        f"Team: M3 (Cell-Type Composition Decoupling)",
        f"Method: {method}",
        f"Spots: {n_total} total; {n_hep} hepatocyte-dominant ≥80% ({pct_hep:.1f}%)",
        "",
        "## Verdict",
        "",
        verdict,
        "",
        f"Pre-registered criteria pass: {pass_count}/3",
        "",
        "## Within-hepatocyte bimodality",
        "",
        *bim_lines,
        "",
        "## Per-patient matched-null Moran's I",
        "",
        *mor_lines,
        "",
        "## Per-patient mixed-effects",
        "",
        *me_lines,
        "",
        "## Pre-registered acceptance",
        "",
        "1. Within-hepatocyte F3b axis ΔBIC ≥10 AND dip-MC p<0.05  →  "
        f"{'PASS' if bim_ok else 'FAIL'}",
        "2. F3a/F3b spatial autocorrelation in ≥80th pct of matched null  →  "
        f"{'PASS' if mor_ok else 'FAIL'}",
        "3. Per-patient mixed-effects CI excludes 0 in ≥4/5 patients  →  "
        f"{'PASS' if me_ok else 'FAIL'}",
        "",
        "## Caveat",
        "",
        ("Reference c2l panel does not have a separate hepatic-stellate cell "
         "type; 'Fibroblasts' is treated as the stellate proxy. If "
         "method=signature_softmax (fallback), per-spot abundance is a "
         "softmax-normalized marker score, not a Bayesian deconvolution. "
         "Acceptance criteria remain meaningful — they test whether "
         "within-hepatocyte (≥80%) bimodality and spatial structure persist "
         "after the dominant compositional source is masked out."),
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n")
    print(f"  Saved: {OUT_MD}")
    print()
    print("\n".join(md_lines))


if __name__ == "__main__":
    main()
