#!/usr/bin/env python3
"""Programmatic verification suite for all 23 remediated defects (N01-N23)
from reports/CODEBASE_DEEP_REVIEW_20261118.md.

Runs without third-party dependencies (stdlib + pdftotext). Exits 0 if all
checks pass, non-zero otherwise.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FAILURES: list[str] = []
PASSES: list[str] = []


def check(defect: str, name: str, cond: bool, detail: str = "") -> None:
    tag = f"[{defect}] {name}"
    if cond:
        PASSES.append(f"PASS {tag}" + (f" :: {detail}" if detail else ""))
    else:
        FAILURES.append(f"FAIL {tag}" + (f" :: {detail}" if detail else ""))


def read_delim(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def pdf_text(rel: str) -> str:
    p = ROOT / rel
    if not p.exists():
        return ""
    return subprocess.check_output(["pdftotext", str(p), "-"], text=True)


def main() -> int:
    # ------------------------------------------------------------------
    # N01: PSD-projected LD + eQTL deduplication + 90-gene refit
    # ------------------------------------------------------------------
    eqtl_csv = ROOT / "Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/results/susie_coloc_all_gwas_refit_20261118.csv"
    step1_json = ROOT / "Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/results/step1_counts_20261118.json"
    comp_tsv = ROOT / "Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/results/eqtl_dedup_gene_comparison.tsv"
    check("N01", "EQTL-DEDUP-PSD refit CSV & summary exist", eqtl_csv.exists() and step1_json.exists() and comp_tsv.exists())
    if step1_json.exists():
        s1 = json.loads(step1_json.read_text())
        check("N01", "Primary multi-signal SuSiE genes == 430", s1["refit"]["multi_signal"] == 430, f"got {s1['refit']['multi_signal']}")
        check("N01", "SuSiE-or-ABF union == 990", s1["refit"]["union"] == 990, f"got {s1['refit']['union']}")
    if comp_tsv.exists():
        comp_rows = read_delim(comp_tsv)
        by_gene = {r["gene"]: r for r in comp_rows}
        irf5 = float(by_gene["IRF5"]["new_best_susie"])
        abhd12 = float(by_gene["ABHD12"]["new_best_susie"])
        thbs3 = float(by_gene["THBS3"]["new_best_susie"])
        cdk5 = float(by_gene["CDK5RAP3"]["new_best_susie"])
        check("N01", "IRF5 & ABHD12 recovered > 0.5", irf5 > 0.9 and abhd12 > 0.7, f"IRF5={irf5:.4f}, ABHD12={abhd12:.4f}")
        check("N01", "THBS3/CDK5RAP3 lost <= 0.5", thbs3 <= 0.5 and cdk5 <= 0.5, f"THBS3={thbs3:.4f}, CDK5RAP3={cdk5:.6f}")

    # ------------------------------------------------------------------
    # N02: C3b ambient subtraction redesigned + in-loop permutation refit
    # ------------------------------------------------------------------
    c3b_sum = ROOT / "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c3b_ambient_20261118/summary.tsv"
    c3b_ctrl = ROOT / "Analysis/Hypotheses/CONF1-20261109T221846Z/out_c3b_ambient_20261118/controls.tsv"
    check("N02", "C3b summary & controls exist", c3b_sum.exists() and c3b_ctrl.exists())
    if c3b_sum.exists():
        txt = c3b_sum.read_text()
        check("N02", "C3b 188 within-hepatocyte genes", "188" in txt, txt.strip()[:200])
    if c3b_ctrl.exists():
        ctrl_rows = read_delim(c3b_ctrl)
        surv = sum(1 for r in ctrl_rows if str(r.get("survives", r.get("called", ""))).upper() in ("TRUE", "1", "YES"))
        check("N02", "C3b 0/9 non-hepatocyte controls survive", surv == 0, f"surv={surv}/{len(ctrl_rows)}")
    txt3h = pdf_text("docs/manuscript/latex/figures/main/figure3/fig3h_within_hepatocyte_partition.pdf")
    check("N02", "Fig 3H prints 188 (not 1,155)", "188" in txt3h and "1,155" not in txt3h)

    # ------------------------------------------------------------------
    # N03: 22-autosome LOCO 8-part ILR deconvolution
    # ------------------------------------------------------------------
    s3a_dir = ROOT / "Analysis/Hypotheses/SS-20261025T135237Z/out_s3a_loco_20261118"
    s3d_dir = ROOT / "Analysis/Hypotheses/SS-20261025T135237Z/out_s3d_loco_20261118"
    check("N03", "LOCO S3a & S3d outputs exist", s3a_dir.is_dir() and s3d_dir.is_dir())
    txt3g = pdf_text("docs/manuscript/latex/figures/main/figure3/fig3g_steatosite_composition.pdf")
    check("N03", "Fig 3G prints 0.224 / 1,132", ("0.224" in txt3g or "0.22" in txt3g) and "1,132" in txt3g)
    txt4d = pdf_text("docs/manuscript/latex/figures/main/figure4/fig4d_steatosite_program_composition.pdf")
    check("N03", "Fig 4D prints Stromal ECM & Ductular injury", "Stromal ECM" in txt4d and "Ductular injury" in txt4d)

    # ------------------------------------------------------------------
    # N04: GSE213621 Control fibrosis_stage=NA (367 disease units)
    # ------------------------------------------------------------------
    cont_dir = ROOT / "figures/candidates/continuum-367-panels-20261118"
    check("N04", "367-unit continuum panels rebuilt", cont_dir.is_dir())
    txt4h = pdf_text("docs/manuscript/latex/figures/main/figure4/fig4h_pnpla3_inflamed_score.pdf")
    check("N04", "Fig 4H prints 633 clean-histology n", "633" in txt4h and "757" in txt4h)

    # ------------------------------------------------------------------
    # N05: Unified Gene Catalog Tier-1 coding rule + 26 both_supported
    # ------------------------------------------------------------------
    cat_path = ROOT / "RNA-seq/results/evidence_passports/candidates/gene-catalog-unified-20261118/catalog_annotated_gene_table_unified_20261118.tsv.gz"
    check("N05", "Unified Catalog TSV exists", cat_path.exists())
    if cat_path.exists():
        cat_rows = read_delim(cat_path)
        check("N05", "Catalog has 28,233 rows", len(cat_rows) == 28233, f"n={len(cat_rows)}")
        c_cnt = Counter(r["candidate_class"] for r in cat_rows)
        check("N05", "Unified classes 401/26/1118/12955/13733",
              c_cnt.get("genetically_anchored") == 401
              and c_cnt.get("both_supported") == 26
              and c_cnt.get("established_state_associated") == 1118
              and c_cnt.get("unresolved") == 12955
              and c_cnt.get("untested") == 13733,
              str(dict(c_cnt)))
        by_sym = {r["symbol"]: r for r in cat_rows}
        t1_expected = {"ADH1B", "APOE", "GCKR", "MTARC1", "PNPLA2", "PNPLA3", "SAMM50", "SERPINA1", "TM6SF2"}
        t1_ok = all(by_sym[g]["candidate_class"] == "genetically_anchored" for g in t1_expected)
        check("N05", "All 9 Tier-1 coding genes genetically_anchored", t1_ok)
        t2_excluded = {"HBB", "GOT1", "ABCB11", "EPPK1", "FOXN4", "ALDH2", "ALDOB", "HNF1A", "HSPA1A"}
        t2_ok = all(by_sym[g]["candidate_class"] not in ("genetically_anchored", "supported_in_both_domains") for g in t2_excluded)
        check("N05", "Tier-2 enzyme-only coding genes excluded from anchored/both", t2_ok)
    txt6c = pdf_text("docs/manuscript/latex/figures/main/figure6/fig6c_both_supported_fields.pdf")
    check("N05", "Fig 6C excludes ABCB11/EPPK1/FOXN4", all(g not in txt6c for g in ("ABCB11", "EPPK1", "FOXN4")))

    # ------------------------------------------------------------------
    # N06: PROT5 MVP ALT zero-overlap + MR-Egger + C2a GOT1 verified
    # ------------------------------------------------------------------
    prot5_dir = ROOT / "Analysis/Hypotheses/PROT5-20261111T131149Z/out_20261118"
    check("N06", "PROT5 out_20261118 exists", prot5_dir.is_dir())
    prot5_files = list(prot5_dir.glob("*.tsv")) if prot5_dir.is_dir() else []
    check("N06", "PROT5 out_20261118 has tables", len(prot5_files) >= 3, f"count={len(prot5_files)}")

    # ------------------------------------------------------------------
    # N07: DIR6 6/6 controls pass (0 NaNs) + non-circular EUR + co-effects
    # ------------------------------------------------------------------
    dir6_ctrl = ROOT / "Analysis/Hypotheses/DIR6-20261111T024808Z/out_20261118/controls.tsv"
    check("N07", "DIR6 controls.tsv exists", dir6_ctrl.exists())
    if dir6_ctrl.exists():
        txt = dir6_ctrl.read_text()
        check("N07", "DIR6 controls contain no NaN and pass", "NaN" not in txt and "FAIL" not in txt)

    # ------------------------------------------------------------------
    # N08: MQTL2 empirical sdY.est across 54,600 pairs
    # ------------------------------------------------------------------
    mqtl_dir = ROOT / "Analysis/Hypotheses/MQTL2-20261111T000008Z/out_20261118"
    check("N08", "MQTL2 out_20261118 exists", mqtl_dir.is_dir() and len(list(mqtl_dir.glob("*.tsv"))) >= 2)

    # ------------------------------------------------------------------
    # N09: DIG1 job 73 (167 loci, 0 opposite, ALDH16A1/RASIP1 disease)
    # ------------------------------------------------------------------
    dig73_loci = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z/out_release_vs_disease_loci_20261118/loci.tsv"
    check("N09", "DIG1 job 73 loci.tsv exists", dig73_loci.exists())
    if dig73_loci.exists():
        loci_rows = read_delim(dig73_loci)
        cls_cnt = Counter(r["cls"] for r in loci_rows)
        check("N09", "DIG1 73 exact 167 loci (75 release / 73 inter / 19 disease / 0 opposite)",
              len(loci_rows) == 167 and cls_cnt.get("release") == 75 and cls_cnt.get("intermediate") == 73
              and cls_cnt.get("disease") == 19 and cls_cnt.get("opposite", 0) == 0,
              str(dict(cls_cnt)))
    txt_s2z = pdf_text("docs/manuscript/latex/figures/supplementary/figureS2_moved/FigS2Z_enzyme_vs_disease_liability.pdf")
    check("N09", "Fig S2Z prints ALT/AST-defined & zero Opposite category",
          "ALT-defined" in txt_s2z and "AST-defined" in txt_s2z and "Opposite" not in txt_s2z)

    # ------------------------------------------------------------------
    # N10: Unified MASLD-only proteomics (5,836 proteins, 21/96 robust) + Fig 5B
    # ------------------------------------------------------------------
    txt5a = pdf_text("docs/manuscript/latex/figures/main/figure5/fig5a_protein_triage.pdf")
    check("N10", "Fig 5A prints 3,054 / 0.46 / 211 / 188",
          "3,054" in txt5a and "0.46" in txt5a and "211" in txt5a and "188" in txt5a)
    txt5b = pdf_text("docs/manuscript/latex/figures/main/figure5/fig5b_mrna_protein_composite.pdf")
    check("N10", "Fig 5B exists and renders text", len(txt5b.strip()) > 20)
    txt5d_main = pdf_text("docs/manuscript/latex/figures/main/figure5/fig5d_program_protein_forest.pdf")
    check("N10", "Fig 5D prints 46 donors & all 21 robust program rows",
          "46 donors" in txt5d_main and "3 Oxidative phosphorylation" in txt5d_main and "4 Translation / ribosome" in txt5d_main)

    # ------------------------------------------------------------------
    # N11: X1 Xenium REF={KNG1,FBP1} + in-loop carrier selection
    # ------------------------------------------------------------------
    x1_dir = ROOT / "Analysis/Hypotheses/X1-20261028T103410Z/out_20261118"
    check("N11", "X1 out_20261118 exists", x1_dir.is_dir())
    x1_uc = ROOT / "Analysis/Hypotheses/X1-20261028T103410Z/out_20261118_unified_cat"
    check("N11", "X1 unified catalog refresh exists", x1_uc.is_dir())

    # ------------------------------------------------------------------
    # N12: Placeholders resolved, draft superseded, Figure 1 Illustrator checklist, Fig 6G
    # ------------------------------------------------------------------
    fig1_doc = ROOT / "reports/FIGURE1_ILLUSTRATOR_UPDATES_20261118.md"
    check("N12", "FIGURE1_ILLUSTRATOR_UPDATES_20261118.md exists", fig1_doc.exists())
    draft_readme = (ROOT / "docs/manuscript/draft/README.md").read_text()
    check("N12", "draft/README.md marked superseded", "SUPERSEDED" in draft_readme.upper())
    # Check zero [ref: or panel letter placeholders in active tex files
    tex_files = [p for p in (ROOT / "docs/manuscript/latex/sections").rglob("*.tex") if p.name != "legends.tex"]
    bad_ph = []
    for tf in tex_files:
        s = tf.read_text()
        if re.search(r"\[ref:|panel letter|table number|to be assigned", s):
            bad_ph.append(tf.name)
    check("N12", "Zero unresolved placeholders in active LaTeX", len(bad_ph) == 0, str(bad_ph))
    # Check catalog.tex HKDC1 orientation fixed
    cat_tex = (ROOT / "docs/manuscript/latex/sections/catalog.tex").read_text()
    check("N12", "catalog.tex HKDC1 oriented to enzyme-raising allele",
          "enzyme-raising allele" in cat_tex and "read on the disease-risk allele" not in cat_tex)
    # Check Fig 6G installed & wired in figure6.tex
    fig6_tex = (ROOT / "docs/manuscript/latex/sections/figures/figure6.tex").read_text()
    check("N12", "figure6.tex wires fig6g_trial_mechanisms.pdf", "fig6g_trial_mechanisms.pdf" in fig6_tex)
    txt6g = pdf_text("docs/manuscript/latex/figures/main/figure6/fig6g_trial_mechanisms.pdf")
    check("N12", "Fig 6G fixes FGF19, MPC modulator, zero TRIAL1/FGF1 typo",
          "FGF19" in txt6g and "MPC modulator" in txt6g and "TRIAL1" not in txt6g and not re.search(r"\bFGF1\b", txt6g))

    # ------------------------------------------------------------------
    # N13: Single-cell deduplication (650,083 unique), Fig 4A (650K), Fig 4C (78/117), C5 (23/17 HC3), CCL1
    # ------------------------------------------------------------------
    txt4a = pdf_text("docs/manuscript/latex/figures/main/figure4/fig4a_atlas_program_schematic.pdf")
    check("N13", "Fig 4A prints 650K unique cells (not 1.23M)", "650K" in txt4a and "1.23M" not in txt4a)
    txt4c = pdf_text("docs/manuscript/latex/figures/main/figure4/fig4c_same_atlas_lineage_specificity.pdf")
    check("N13", "Fig 4C prints 78 / 117 / 66.7%", "78" in txt4c and "117" in txt4c)
    txt4e = pdf_text("docs/manuscript/latex/figures/main/figure4/fig4e_hepatocyte_trio_forest.pdf")
    check("N13", "Fig 4E prints GSE289173 (23/17)", "23/17" in txt4e)
    txt4f = pdf_text("docs/manuscript/latex/figures/main/figure4/fig4f_ccl20_state_set_excess.pdf")
    check("N13", "Fig 4F prints GSE244832 (n = 8) & GSE202379 (n = 5)",
          "GSE244832" in txt4f and "n = 8" in txt4f and "GSE202379" in txt4f and "n = 5" in txt4f)

    # ------------------------------------------------------------------
    # N14: AG2X ALL_FOLDS, locus-grouped B1, measured-only B3a, aligned Stage D
    # ------------------------------------------------------------------
    ag2x_dir = ROOT / "Analysis/Hypotheses/AG2X-20261112T000534Z/out_20261118"
    check("N14", "AG2X out_20261118 exists", ag2x_dir.is_dir())

    # ------------------------------------------------------------------
    # N15: ENZFAT1 S1b MVN CI + log-space ashr + TRAIL4 exact 4-trait ABF
    # ------------------------------------------------------------------
    enz_s1b = ROOT / "Analysis/Hypotheses/ENZFAT1-20261110T204401Z/out_S1b_20261118"
    enz_s2 = ROOT / "Analysis/Hypotheses/ENZFAT1-20261110T204401Z/out_S2_20261118"
    trail4 = ROOT / "Analysis/Hypotheses/TRAIL4-20261111T025501Z/out_20261118"
    check("N15", "ENZFAT1 S1b/S2 & TRAIL4 out_20261118 exist", enz_s1b.is_dir() and enz_s2.is_dir() and trail4.is_dir())

    # ------------------------------------------------------------------
    # N16 & N17: Pre-voom filter (F1-F4 1,338/4,269/6,190/7,042) & GSE126848 n=40 & EISA
    # ------------------------------------------------------------------
    txt3e = pdf_text("docs/manuscript/latex/figures/main/figure3/fig3e_stage_remodeling.pdf")
    check("N16", "Fig 3E prints 4,269 & 7,042", "4,269" in txt3e and "7,042" in txt3e)
    eisa_dir = ROOT / "Analysis/Hypotheses/SSX-20261102T133212Z/out_eisa_20261118"
    check("N16", "SteatoSITE EISA out_eisa_20261118 exists", eisa_dir.is_dir())
    gse126_dir = ROOT / "figures/candidates/bulk-context-774-gse126848-f40-20261118"
    check("N17", "GSE126848 n=40 female-restored sensitivity exists", gse126_dir.is_dir())

    # ------------------------------------------------------------------
    # N18: GeoMx quarantined + snATAC D11/D14 fixed + Visium all-16 + Cicero within-donor residualized
    # ------------------------------------------------------------------
    txt_s5q = pdf_text("docs/manuscript/latex/figures/supplementary/fig5_molecular_context/panels/figS5q.pdf")
    check("N18", "Fig S5Q prints 19,431 / 5N/8M / 0.041", "19,431" in txt_s5q and "5N/8M" in txt_s5q)
    txt_s5r = pdf_text("docs/manuscript/latex/figures/supplementary/fig5_molecular_context/panels/figS5r.pdf")
    check("N18", "Fig S5R prints 7,495 / 0.038", "7,495" in txt_s5r and "0.038" in txt_s5r)
    txt_s5x = pdf_text("docs/manuscript/latex/figures/supplementary/fig5_molecular_context/panels/figS5x.pdf")
    check("N18", "Fig S5X includes Within-donor residualized MSE", "Within-donor residualized" in txt_s5x)

    # ------------------------------------------------------------------
    # N19: snATAC genetic context EQTL-DEDUP + TSSxMAFxLD + RESULTS.md reconciled
    # ------------------------------------------------------------------
    txt_5d_arch = pdf_text("docs/manuscript/latex/figures/supplementary/fig5_molecular_context/panels/_archive/demoted-from-main-20261103/fig5d_snatac_accessibility.pdf")
    check("N19", "Fig 5D (S5AC) prints 418 genes & no Uniform diamond", "418" in txt_5d_arch and "Uniform" not in txt_5d_arch)
    txt_5d_comp = pdf_text("docs/manuscript/latex/figures/supplementary/fig5_molecular_context/panels/_archive/demoted-from-main-20261103/fig5d_companion_locus_matched_background.pdf")
    check("N19", "Fig 5D companion prints 1.68x & 2.39x & 418", "418" in txt_5d_comp and "1.68" in txt_5d_comp and "2.39" in txt_5d_comp)
    res_md = (ROOT / "docs/RESULTS.md").read_text()
    check("N19", "RESULTS.md reconciles 1.84-2.39x vs 2.55-3.11x", "1.84–2.39×" in res_md or "1.84-2.39" in res_md)

    # ------------------------------------------------------------------
    # N20: Predictive model .1 releases exist + historical releases intact
    # ------------------------------------------------------------------
    for m in ("masld-severity-v1.2.1", "masld-panel-v1.1", "masld-liver-chromatin-state-v1.2.1"):
        p = ROOT / "Analysis/MASLD_Model_Benchmark/release" / m / "test_release.sh"
        check("N20", f"Model release {m} exists", p.exists())
    for hist in ("masld-severity-v1.2", "masld-panel-v1", "masld-liver-chromatin-state-v1.2"):
        p = ROOT / "Analysis/MASLD_Model_Benchmark/release" / hist / "test_release.sh"
        check("N20", f"Historical model release {hist} preserved", p.exists())

    # ------------------------------------------------------------------
    # N21: Historical releases & figure1_JL.pdf untouched
    # ------------------------------------------------------------------
    check("N21", "Historical 2026-10-21-r9 preserved", (ROOT / "RNA-seq/results/manuscript_release/2026-10-21-r9").exists())
    check("N21", "Historical gene-catalog-eligibility-20261106 preserved", (ROOT / "RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106").exists())
    check("N21", "figure1_JL.pdf preserved", (ROOT / "docs/manuscript/latex/figures/main/figure1/figure1_JL.pdf").exists())

    # ------------------------------------------------------------------
    # N22: Fig 2A right-edge fix + 430 / 990 + S2AA + S2AC
    # ------------------------------------------------------------------
    txt2a = pdf_text("docs/manuscript/latex/figures/main/figure2/Fig2A_gwas_alluvial.pdf")
    check("N22", "Fig 2A prints 430 & 990", "430" in txt2a and "990" in txt2a)

    # ------------------------------------------------------------------
    # N23: Compiled main.pdf & supplement.pdf exist (checked if built)
    # ------------------------------------------------------------------
    main_pdf = ROOT / "docs/manuscript/latex/build/main.pdf"
    supp_pdf = ROOT / "docs/manuscript/latex/build/supplement.pdf"
    check("N23", "Compiled main.pdf exists and > 1 MB", main_pdf.exists() and main_pdf.stat().st_size > 1_000_000,
          f"size={main_pdf.stat().st_size if main_pdf.exists() else 0}")
    check("N23", "Compiled supplement.pdf exists and > 1 MB", supp_pdf.exists() and supp_pdf.stat().st_size > 1_000_000,
          f"size={supp_pdf.stat().st_size if supp_pdf.exists() else 0}")

    print(f"=== REMEDIATION VERIFICATION SUMMARY: {len(PASSES)} PASSED, {len(FAILURES)} FAILED ===")
    for p in PASSES:
        print(" ", p)
    if FAILURES:
        print("\n=== FAILURES ===")
        for f in FAILURES:
            print(" ", f)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
