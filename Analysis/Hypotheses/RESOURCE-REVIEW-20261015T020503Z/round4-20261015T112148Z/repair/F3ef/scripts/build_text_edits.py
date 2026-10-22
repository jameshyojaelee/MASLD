"""F3ef repair: text edits written against the current docs; validates each old_text (exactly once),
checks for em dashes and evidence tags in new_text, and writes text_edits_repair.json."""
import json, re
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
OUT = ROOT + "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F3ef/text_edits_repair.json"
REP = "figures/candidates/review-r3-fixes-20261015T112148Z/F3ef/repair"
LEG = "docs/manuscript/05_figure_legends.md"; F3 = "docs/manuscript/draft/fig3.md"; RES = "docs/RESULTS.md"
STA = "docs/STATUS.md"; MET = "docs/manuscript/METHODS.md"; PAP = "docs/PAPER.md"

E = {}
E["R03-03"] = [
 (LEG, "All panels use the corrected 774-unit substrate. Per-cohort membership,",
  "All panels draw on the corrected 774-unit substrate; the stage contrasts in E and F use the 387 units with an exact F0–F4 fibrosis stage in GSE130970, GSE135251 and GSE162694 (F0/F1/F2/F3/F4 = 103/102/88/69/25). Per-cohort membership,"),
 (LEG, "**(A–B)** Five-cohort coverage (326, 207, 140, 71 and 30 units) and documented\nfibrosis-by-NAS support.",
  "**(A–B)** Five-cohort coverage (326, 207, 140, 71 and 30 units) and documented\nfibrosis-by-NAS support. GSE213621 records fibrosis only as source bins (F0–F1, F2, F3–F4) and enters no exact-stage contrast."),
 (F3, "On the installed 774-unit stage analysis, Figure 3E contains 1,338 / 4,266 / 6,190 /\n7,043 BH-significant genes at F1 / F2 / F3 / F4 versus F0.",
  "Figure 3E contains 1,338 / 4,266 / 6,190 / 7,043 BH-significant genes at F1 / F2 / F3 / F4 versus F0, from the 387 units with an exact F0–F4 fibrosis stage in GSE130970, GSE135251 and GSE162694."),
 (F3, "F0-centered NMF activity. Each contrast uses eligible cohorts and fixed cohort",
  "F0-centered NMF activity. Each contrast uses the three cohorts with an exact F0–F4 fibrosis stage (GSE130970, GSE135251 and GSE162694; 387 units); GSE213621 records only F0–F1, F2 and F3–F4 bins and is not used. Each model has fixed cohort"),
 (MET, "Adjacent fibrosis stages and adjacent NAS groups were compared\ncross-sectionally in the selected five-cohort fragment-count substrate.",
  "Adjacent fibrosis stages and adjacent NAS groups were compared\ncross-sectionally in the selected five-cohort fragment-count substrate. Fibrosis contrasts used the three cohorts with an exact F0–F4 fibrosis stage (GSE130970, GSE135251 and GSE162694; 387 units); GSE213621 deposits only F0–F1, F2 and F3–F4 bins and was not used in these contrasts."),
]
E["G2-01"] = [
 (LEG, "**(E)** F1–F4-versus-F0 counts are 1,338 / 4,266 / 6,190 / 7,043 at BH q <0.05\nin each complete 24,115-gene family, without a fold-change filter. This rule\ndiffers from the canonical disease/control gate. The F0 reference contains\n103 units (14 source controls and 89 disease biopsies).",
  "**(E)** F1–F4-versus-F0 counts are 1,338 / 4,266 / 6,190 / 7,043 at BH q <0.05\nin each complete 24,115-gene family, without a fold-change filter, from 387 units in GSE130970, GSE135251 and GSE162694 with fixed cohort and inferred-sex terms. This rule\ndiffers from the canonical disease/control gate. The F0 reference contains\n103 units (14 source controls and 89 disease biopsies). The F1 count rests mainly on GSE130970: without that cohort 14 genes remain (155 units), whereas omitting any one cohort leaves 1,934 to 9,546 genes for F2–F4. In a cohort-by-stage model, 95.8% of the 1,338 F1 effects have the pooled sign in all three cohorts. The pooled F4 count averages cohort effects of different size: 4,047 of the 7,043 genes have a cohort-by-stage interaction at BH q <0.05 (2-df moderated F over 24,115 genes), and 75.4% have the pooled sign in all three cohorts."),
 (F3, "depend on sample size, effect and uncertainty and cannot locate remodeling onset.",
  "depend on sample size, effect and uncertainty and cannot locate remodeling onset. The F1 count rests mainly on GSE130970: without it 14 of the 1,338 genes remain (155 units), and fitted alone the cohorts give 252 (GSE130970, 50 units), 3 (GSE135251, 90 units) and 0 (GSE162694, 65 units) genes. In a cohort-by-stage model, 95.8% of the 1,338 effects have the pooled sign in all three cohorts, so the cohorts mostly agree on the direction of the F0-to-F1 differences, but only GSE130970 gives many calls on its own. The F2, F3 and F4 versus F0 counts stay between 1,934 and 9,546 when any one cohort is omitted. The later adjacent steps do depend on single cohorts: F1 to F2 gives 0 genes without GSE162694 (133 units), and F2 to F3 gives 2 without GSE135251 (54 units). Effect sizes also differ between cohorts at F4: 4,047 of the 7,043 F4-versus-F0 genes have a cohort-by-stage interaction at BH q <0.05."),
 (RES, "the same totals are in `figures/candidates/fig3e-s3-774-20261011T165944Z/source_tables/fig3e_deg_counts.tsv`.",
  "the same totals are in `figures/candidates/fig3e-s3-774-20261011T165944Z/source_tables/fig3e_deg_counts.tsv`. All fibrosis contrasts use GSE130970, GSE135251 and GSE162694 only (387 units); GSE213621 stage is binned and not used. Refitted without GSE130970 / GSE135251 / GSE162694 (review repair job 7550, same design; the base fits reproduce the release exactly), the BH q<0.05 counts are F0→F1 14 / 929 / 651 (155 / 115 / 140 units), F1→F2 240 / 205 / 0 (154 / 93 / 133), F2→F3 464 / 2 / 280 (138 / 54 / 122), F3→F4 102 / 6,565 / 508 (82 / 29 / 77), F2 vs F0 3,293 / 4,105 / 1,934, F3 vs F0 4,804 / 3,994 / 4,301 and F4 vs F0 4,800 / 9,546 / 4,830. All 14 F0→F1 calls without GSE130970 are among the 1,338, and 4,596 of the 7,043 F4 vs F0 calls remain without it. F0→F1 fitted in each cohort alone gives 252 / 3 / 0 (50 / 90 / 65 units). Of the 1,338 F0→F1 calls, 95.8% have the pooled sign in all three cohort coefficients of a cohort-by-stage model (95.4% in separate single-cohort fits). F4 vs F0: 4,047 of the 7,043 calls have a cohort-by-stage interaction at BH q<0.05 (2-df moderated F over 24,115 genes), and 75.4% have the pooled sign in all three cohorts (review job 1032). Tables: `" + REP + "/fibrosis_counts_cohort_omission.tsv` and `figures/candidates/review-r3-fixes-20261015T112148Z/F3ef/fig3e_cohort_heterogeneity.tsv`."),
 (STA, "4,266, 6,190 and 7,043 for F1 to F4, and versus NAS0 5,338, 9,337 and 9,841 for\nNAS 1-2, 3-4 and 5-8.",
  "4,266, 6,190 and 7,043 for F1 to F4 (387 units in GSE130970, GSE135251 and GSE162694; the F1 count is 14 without GSE130970), and versus NAS0 5,338, 9,337 and 9,841 for\nNAS 1-2, 3-4 and 5-8 (NAS0 is 45 source controls and 1 disease unit, so these are mainly disease-versus-control counts)."),
 (MET, "Non-estimable branches remain explicit. Early-stage GNMT and MAT1A claims that\ndo not survive these comparisons are not retained.",
  "Non-estimable branches remain explicit. Early-stage GNMT and MAT1A claims that\ndo not survive these comparisons are not retained. As a further sensitivity, every fibrosis contrast was refitted without each cohort in turn with the same design, and F1 versus F0 was fitted in each cohort alone (`~ inferred_sex + group`). Cohort heterogeneity of F1 and F4 versus F0 was tested with a cohort-by-stage interaction (`~ inferred_sex + dataset + dataset:group`), using a 2-df moderated F test with BH correction over 24,115 genes."),
 (RES, "All-gene F0→F1 effects correlate at Spearman rho 0.950 between the all-F0 and\ndisease-only arms (`f0_arm_effect_agreement.tsv`).",
  "All-gene F0→F1 effects correlate at Spearman rho 0.950 between the all-F0 and\ndisease-only arms (`f0_arm_effect_agreement.tsv`). The all-F0 F0→F1 count rests mainly on GSE130970: 14 of the 1,338 calls remain without it (155 units)."),
 (STA, "814 with a control-status term; GNMT F1→F2 has BH q 0.0562 and is no longer below 0.05.",
  "814 with a control-status term; GNMT F1→F2 has BH q 0.0562 and is no longer below 0.05. The all-F0 F0→F1 count rests mainly on GSE130970 (14 of 1,338 calls without it, 155 units), and GNMT and MAT1A F0→F1 have q 0.207 and 0.0892 without GSE130970."),
]
E["G2-03"] = [
 (LEG, "**(F)** Fifty-four genes with BH FDR <0.05 in at least one adjacent fibrosis\ncontrast and membership in the native 436-gene primary SuSiE set of the corrected\nCOLOC r2 input. Four complete 24,115-gene contrast families give 216 displayed\neffects and 59 significant marks (30 / 13 / 11 / 5). Colors show\nhigher-minus-lower stage log2FC, clipped at ±0.688 for display; black dots mark\nBH significance.",
  "**(F)** Fifty-four genes with BH FDR <0.05 in at least one adjacent fibrosis\ncontrast and membership in the native 436-gene primary SuSiE set of the corrected\nCOLOC r2 input. FDR is controlled within each contrast, not over the union of\nfour contrasts; with the budget split across contrasts (BH q <0.0125 in each),\nwhich bounds the FDR over all four 24,115-gene families at 5%, 15 of the 54 genes remain. Four complete 24,115-gene contrast families give 216 displayed\neffects and 59 significant marks (30 / 13 / 11 / 5). Colors show\nhigher-minus-lower stage log2FC, clipped at ±0.688 for display; black dots mark\nBH significance. Twenty-nine genes are selected only through an F0→F1 mark. None of the 30 F0→F1 marks keeps BH q <0.05 without GSE130970, and 18 keep it with disease-only F0; none of the 13 F1→F2 marks keeps it without GSE162694."),
 (F3, "keeping the measurements separate. Fifty-four genes belong to the native 436-gene\nprimary SuSiE set and pass BH FDR <0.05 in at least one adjacent fibrosis\ncontrast. Four complete 24,115-gene families give 216 displayed effects and 59\nsignificant marks (30 / 13 / 11 / 5).",
  "keeping the measurements separate. Fifty-four genes belong to the native 436-gene\nprimary SuSiE set and pass BH FDR <0.05 in at least one adjacent fibrosis\ncontrast. Each contrast is its own BH family, so the union carries no joint 5%\nFDR; 15 genes pass when the budget is split across the four contrasts (BH\nq <0.0125 in each, which bounds the FDR over all four 24,115-gene families at 5%; 11 / 2 / 2 / 1 marks). Four complete 24,115-gene families give 216 displayed effects and 59\nsignificant marks (30 / 13 / 11 / 5). Twenty-nine of the 54 genes are selected only through an F0-to-F1 mark, and those marks share the cohort dependence of the F1 count: none of the 30 keeps BH q <0.05 without GSE130970, and 18 keep it with disease-only F0. Later marks also depend on single cohorts: none of the 13 F1-to-F2 marks keeps q <0.05 without GSE162694, and 1 of the 59 marks (OMG, F2 to F3) keeps it after each single-cohort omission. Each omission also removes 12% to 69% of a contrast's units, which lowers power."),
 (RES, "(65 before): 30, 13, 11 and 5 for F0→F1, F1→F2, F2→F3 and F3→F4 (32, 21, 7 and 5\nbefore).",
  "(65 before): 30, 13, 11 and 5 for F0→F1, F1→F2, F2→F3 and F3→F4 (32, 21, 7 and 5\nbefore). Twenty-nine of the 54 genes are selected only through an F0→F1 mark. Of the 30 F0→F1 marks, 18 keep BH q < 0.05 (same sign) with disease-only F0 and 0 without GSE130970 (smallest q 0.103). Marks kept (q < 0.05, same sign) without GSE130970 / GSE135251 / GSE162694, for F0→F1, F1→F2, F2→F3 and F3→F4: 0/6/6/0, 15/3/1/2 and 12/0/3/4; 1 of the 59 (OMG, F2→F3) is kept under all three omissions (review repair job 7550; `" + REP + "/fig3f_marks_cohort_omission.tsv`). Reselecting with disease-only F0 gives 45 genes and 50 marks (21 / 13 / 11 / 5): 42 of the 54 stay, 12 leave (ACVR2B, HSPA4, MGAT2, MTARC2, MTOR, NUDT4, PABPC4, PFKFB2, PMF1, RC3H2, RPL8, WIPI1) and FCHO2, HPS6 and SEC14L2 enter."),
 (STA, "157 / 94 units; ten genes leave and three enter (list in `docs/RESULTS.md`).",
  "157 / 94 units; ten genes leave and three enter (list in `docs/RESULTS.md`). Twenty-nine of the 54 genes are selected only through an F0→F1 mark; of the 30 F0→F1 marks, 18 keep BH q < 0.05 with disease-only F0 and none without GSE130970."),
]
E["R03-06"] = [
 (LEG, "change. The full fibrosis/NAS, composition and F0-reference sensitivities are\nshown in Figure S3.",
  "change. The full fibrosis/NAS and composition views, the F0-reference arms (disease-only F0, n=89: 783 / 3,946 / 6,459 / 6,929; all F0 with a control-status term, n=103: 814 / 3,679 / 6,166 / 6,640) and the F1-versus-F0 counts without each cohort and in each cohort alone are shown in Figure S3."),
 (LEG, "endpoint DEG agreement, the full fibrosis/NAS remodeling panel, and hepatocyte\nand non-parenchymal composition.",
  "endpoint DEG agreement, the full fibrosis/NAS remodeling panel, hepatocyte\nand non-parenchymal composition, F1–F4-versus-F0 counts under three F0 references (all 103 F0 units; 89 disease-only F0 units; all 103 F0 units with a control-status term), and F1-versus-F0 counts without each cohort (14 / 929 / 651 genes without GSE130970 / GSE135251 / GSE162694; 155 / 115 / 140 units) and in each cohort alone (252 / 3 / 0 genes; 50 / 90 / 65 units)."),
]
E["R12-04"] = [
 (LEG, "NAS counts use BH q <0.05 without a fold-change filter, and continuous\nNAS activity is centered on NAS0 (46 units with fixed scores).",
  "NAS counts use BH q <0.05 without a fold-change filter. The NAS0 reference (46 units) contains 45 source controls and 1 disease biopsy, so the NAS 1–2, 3–4 and 5–8 versus NAS0 counts (5,338, 9,337 and 9,841) mainly contrast disease with controls rather than activity grades. Among disease biopsies, NAS 1–2 to 3–4 gives 3 genes (65 and 117 units; NAS 1–2 includes 2 GSE130970 controls) and NAS 3–4 to 5–8 gives 4,305 genes (117 and 167 units). Continuous NAS activity is centered on NAS0 (46 units with fixed scores)."),
 (RES, "| NAS versus NAS0 | 5–8 vs 0 | 46 / 167 | 3 | 4,571 / 5,270 |",
  "| NAS versus NAS0 | 5–8 vs 0 | 46 / 167 | 3 | 4,571 / 5,270 |\n\nNAS0 is 45 source controls and 1 GSE162694 disease unit (GSE130970 4, GSE135251 10 and GSE162694 31 controls), so the versus-NAS0 rows and the NAS 0→1–2 step mainly contrast disease with controls. NAS 1–2 has 63 disease units and 2 GSE130970 controls; NAS 3–4 (117) and 5–8 (167) are all disease. Among disease biopsies, NAS 1–2→3–4 gives 3 genes and NAS 3–4→5–8 gives 4,305 (1,954 higher, 2,351 lower)."),
]
E["G2-04"] = [
 (F3, "control-status term. GNMT decreases from F0 to F1 across these three fits\n(q=0.00913 / 0.0426 / 0.0430); its F1-to-F2 q is 0.0562 and is not a current\nstage call. MAT1A decreases from F0 to F1 (q=2.45e-4 / 9.5e-4 / 1.6e-3), and\nCYP2C19 decreases from F1 to F2 (q=0.0477).",
  "control-status term. All three genes have lower abundance at F2, F3 and F4 than at F0, with BH q <0.001 under every F0 reference and without GSE130970; omitting GSE135251 or GSE162694 leaves every effect negative, with q ≤0.034 except MAT1A at F4 without GSE135251 (q=0.059; 11 F4 units). The differences at earlier steps depend on single cohorts. GNMT and MAT1A have lower abundance at F1 than at F0 under all three F0 references (GNMT q=0.00913 / 0.0426 / 0.0430; MAT1A q=2.45e-4 / 9.5e-4 / 1.6e-3) but not without GSE130970 (q=0.207 and 0.0892). CYP2C19 has lower abundance at F2 than at F1 (q=0.0477), with q=0.052 with a control-status term and q=0.219 and 0.499 without GSE135251 and GSE162694; GNMT's F1-to-F2 q is 0.0562. These earlier differences are not used to place where the lower abundance begins."),
 (RES, "F0), 0.0426 (disease-only) and 0.0430 (control-status term); MAT1A: 2.45e-4,\n9.5e-4 and 1.6e-3 (`f0_arm_hero_genes.tsv`).",
  "F0), 0.0426 (disease-only) and 0.0430 (control-status term); MAT1A: 2.45e-4,\n9.5e-4 and 1.6e-3 (`f0_arm_hero_genes.tsv`). Neither F0→F1 difference holds without GSE130970 (GNMT q 0.207, MAT1A 0.0892); without GSE135251 they are q 0.0429 and 0.00251, and without GSE162694 0.0338 and 0.00631. CYP2C19 F1→F2 (q 0.0477) is q 0.0239 / 0.219 / 0.499 without GSE130970 / GSE135251 / GSE162694. The F2, F3 and F4 versus F0 differences of all three genes are negative in every fit: q < 0.001 under every F0 reference and without GSE130970, and q ≤ 0.034 without GSE135251 or GSE162694 except MAT1A F4 vs F0 without GSE135251 (q 0.0592; 11 F4 units) (`" + REP + "/hero_genes_cohort_omission.tsv`)."),
 (PAP, "   disease map, GNMT joins MAT1A and CYP2C19 in lower multicohort RNA,\n   cross-sectional fibrosis-stage abundance, and liver protein.",
  "   disease map, GNMT joins MAT1A and CYP2C19 in lower multicohort RNA,\n   lower abundance at F2 to F4 than at F0 in three staged cohorts (the GNMT and\n   MAT1A differences at F1 depend on GSE130970), and liver protein."),
 (RES, "F0→F1 -0.428, q=0.00913, which stays supported on disease-only F0 (q=0.0426) and with a control-status term (q=0.0430).",
  "F0→F1 -0.428, q=0.00913; q=0.0426 on disease-only F0 and 0.0430 with a control-status term, but q=0.207 without GSE130970."),
 (RES, "F0→F1 -0.239, q=2.45e-4; q=9.5e-4 on disease-only F0 and 1.6e-3 with a control-status term.",
  "F0→F1 -0.239, q=2.45e-4; q=9.5e-4 on disease-only F0 and 1.6e-3 with a control-status term, but q=0.0892 without GSE130970."),
 (RES, "F1→F2 -0.697, BH q=0.0477 (q=0.052 with a control-status term)",
  "F1→F2 -0.697, BH q=0.0477 (q=0.052 with a control-status term; q=0.219 without GSE135251 and 0.499 without GSE162694)"),
]

out, bad = {}, []
spans = {}
for fid, edits in E.items():
    out[fid] = []
    for f, old, new in edits:
        s = open(ROOT + f).read()
        n = s.count(old)
        i = s.find(old)
        l0 = s[:i].count("\n") + 1; l1 = l0 + old.count("\n")
        if n != 1: bad.append((fid, f, "count", n))
        if "—" in new: bad.append((fid, f, "em dash"))
        if re.search(r"\[(Inference|Unverified|Speculation|Estimate)", new): bad.append((fid, f, "tag"))
        if "Kleiner" in new or "TMM" in new: bad.append((fid, f, "word"))
        spans.setdefault(f, []).append((i, i + len(old), fid))
        out[fid].append(dict(file=f, line=f"{l0}" if l0 == l1 else f"{l0}-{l1}", old_text=old, new_text=new, claim_level=False))
for f, L in spans.items():
    L.sort()
    for a, b in zip(L, L[1:]):
        if b[0] < a[1]: bad.append(("overlap", f, a[2], b[2]))
print("problems:", bad)
json.dump(out, open(OUT, "w"), indent=1, ensure_ascii=False)
for fid, L in out.items():
    for e in L: print(fid, e["file"], e["line"])
