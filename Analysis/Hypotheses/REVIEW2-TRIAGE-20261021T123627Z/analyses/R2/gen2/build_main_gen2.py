# Round-2 genetics edits (writer gen2). Checks: each old occurs once in the round-1 copy
# (or the repo source when no copy exists); no two edits overlap in the same file.
import json, os
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN = "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z"
AN1 = RUN + "/analyses/AN1/gene_level_20261021T165252Z/results.json"
AN1R = RUN + "/analyses/AN1/RESULT.md"
AN5 = RUN + "/analyses/AN5/RESULT.md"
AN5TW = RUN + "/analyses/AN5/threeway/run_20261021T125538Z/threeway_summary.json"
AN5D = RUN + "/analyses/AN5/pph4_1_diagnostic_summary.tsv"
PR = RUN + "/analyses/R2/gen2/pair_rates_out.txt"

FIG1 = "docs/manuscript/draft/fig1.md"
FIG2 = "docs/manuscript/draft/fig2.md"
FIG3 = "docs/manuscript/draft/fig3.md"
DISC = "docs/manuscript/draft/fig5_discussion.md"
SUPP = "docs/manuscript/draft/supplementary_notes.md"
METH = "docs/manuscript/METHODS.md"
LEG = "docs/manuscript/05_figure_legends.md"
PAPER = "docs/PAPER.md"

E = []
def add(id_, file, old, new, items, srcs, reason, decision=None):
    E.append(dict(id=id_, file=file, old=old, new=new, items=items,
                  number_sources=srcs, decision=decision, reason=reason))

# ---------------- 1. symmetric shared-posterior check ----------------
add("gen2-fig1-1", FIG1,
    "Only first-run positives were rechecked, so some of the 3,769 genes tested at PP.H4 ≤ 0.5 could also fail the check; 3,769 is therefore an upper bound. The 26 shared genes and the colocalized genes are unaffected, because every positive was rechecked.",
    "The rerun applied this check only to first-run positives. In a sensitivity analysis, we applied it to the carried negatives and counted a gene as tested if any tier-1/2 study stayed eligible. Then 41 genes became untestable, the universe was 3,736 genes, and the 26 shared genes stood against 21.1 expected (Mantel-Haenszel OR 1.30, 95% CI 0.84 to 2.02; exact *p* = 0.27). The check cannot add or remove a colocalized gene, and 3,769 remains the prespecified universe (Supplementary Note 3).",
    ["B-G1"],
    [AN1 + ":states_rule_b.transitions (41)", AN1 + ":H57.b_sens (N 3736, E_strat 21.146, MH_OR 1.3003, MH_lo 0.8384, MH_hi 2.0168, p_exact 0.2746)"],
    "D1: replace the round-1 upper-bound wording with the measured rule-(b) sensitivity.", "D1")

add("gen2-fig2-1", FIG2,
    "A carried negative that would fail it could become untestable but not positive, so counts of tested genes, including the 3,769 genes of Figure 1B, are upper bounds.",
    "As a sensitivity analysis, we applied the check to the 12,087 carried negative rows of the primary studies: 343 (2.8%) failed on every signal pair and became untestable, 1,744 failed on some pairs and kept their call, and 10,000 passed (Methods). A carried negative that fails the check becomes untestable, not positive, so the colocalized genes do not change.",
    ["B-G1"],
    [AN1 + ":rows (primary 12087, all_sets_fail 343, some_sets_fail 1744, none_fail 10000)"],
    "Replace the upper-bound sentence with the measured row-level result.", "D1")

add("gen2-fig2-2", FIG2,
    "The colocalized pairs passed the shared-posterior check, but most pairs in the denominators are first-run negatives that were not checked, so both rates are lower bounds for fully checked pairs. Before the check, the rates were 787 of 12,132 (6.5%) and 29 of 771 (3.8%).",
    "Most pairs in these denominators are first-run negatives carried forward. When the shared-posterior check was also applied to them, the denominators fell to 11,797 and 752 pairs (6.3% and 3.7% colocalized). In the first run, before any check, the rates were 787 of 12,132 (6.5%) and 29 of 771 (3.8%).",
    ["B-G1"],
    [PR + ":lines 4-5 (sensitivity tier 2 finite 11797 pos 739; tier 1 finite 752 pos 28), computed by " + RUN + "/analyses/R2/gen2/pair_rates.py from AN1 master_reduced_sensitivity.csv; baseline lines reproduce 12121/739 and 771/28"],
    "The per-pair denominators are now measured after the check, so the lower-bound wording is replaced.", "D1")

add("gen2-fig2-3", FIG2,
    "One gene, MAU2, has PP.H4 of 1.0 for both a direct and an enzyme trait and is counted under the enzyme class.",
    "One gene, MAU2, has PP.H4 of 1.0 in two enzyme GWAS and in the MVP European NAFLD GWAS. Exact ties go to the GWAS whose study name sorts first, so MAU2 is counted under the enzyme class (Methods).",
    ["A2.1", "B-G3"],
    [AN5 + ": section 4 (MAU2 PP.H4 exactly 1 in MVP_ALT_EUR, MVP_AST_EUR, MVP_NAFLD_EUR; driver MVP_ALT_EUR)"],
    "State the declared tie rule that decides MAU2.")

add("gen2-fig3-1", FIG3,
    "The seven SuSiE-untestable genes are drawn separately and are not counted as Neither.",
    "The seven SuSiE-untestable genes are drawn as open circles at their best finite PP.H4 from another tier-1/2 study and are not counted as Neither. The panel shows the prespecified state. After the shared-posterior check of carried negatives (Methods), 40 drawn genes keep no finite PP.H4 (37 Neither, 3 Transcriptomics only) and would leave the panel. If a gene counts as tested when any tier-1/2 study stays eligible, no drawn gene is untestable, and the four groups hold 3,119, 191, 400 and 26 genes.",
    ["B-G1"],
    [AN1 + ":fig3g.dropped_no_finite_pp4 (Neither 37, Transcriptomics 3)", AN1 + ":fig3g.sens.by_state (3119/191/400/26)", AN1 + ":H57.b_sens.n_joint_untestable 40",
     RUN + "/analyses/AN1/gene_level_20261021T165252Z/fig3g_states.tsv (seven state_a_base untestable genes have susie_pp4 = coloc_best_susie_pp4_base)"],
    "State how the untestable genes are drawn and what the check does to the panel.", "D1")

add("gen2-leg-1", LEG,
    "Open circles, seven SuSiE-untestable genes, not counted as Neither.",
    "Open circles, seven SuSiE-untestable genes, drawn at their best finite PP.H4 from another tier-1/2 study and not counted as Neither. The panel shows the prespecified state, before the shared-posterior check of carried negatives (Supplementary Note 3).",
    ["B-G1"],
    [RUN + "/analyses/AN1/gene_level_20261021T165252Z/fig3g_states.tsv (susie_pp4 of the seven untestable genes equals coloc_best_susie_pp4_base)"],
    "Say where the open circles sit on the PP.H4 axis.")

add("gen2-disc-1", DISC,
    "Carried negatives were not rechecked, so 3,918 is\nan upper bound. The 3,335 and 3,319 could move in either direction under a\nrecheck, because each row's PP.H3 comes from its best-PP.H4 signal pair.",
    "When the shared-posterior check was also applied\nto the carried negatives, 3,877 genes kept a finite Tier-1/2 SuSiE result and\n3,306 showed this PP.H3 pattern. The 3,306 is a point value within a range of\n3,101 to 3,416, because rows that lose only some signal pairs have an unknown\nPP.H3 after the check.",
    ["B-G1"],
    [AN1 + ":pp3.carried_PP.H3 (finite 3877, pp3 3306)", AN1 + ":pp3.lower_partial_not_counted.pp3 3101", AN1 + ":pp3.upper_partial_counted.pp3 3416"],
    "Replace the upper-bound wording with the measured PP.H3 counts.", "D1")

add("gen2-meth-1", METH,
    "Among eQTL credible sets in the rechecked primary rows that had no colocalizing pair, 435 of 1,233 (35%) kept less than half of their posterior on variants in the GWAS fit, so an unknown number of carried negatives would fail the check, and counts of tested genes are upper bounds. This partial untestability is distinct from a complete absence of finite results. Full-portfolio nonfinite states comprise 921,663 ABF-fallback rows, 5,016 ABF-only rows and the 25 shared-posterior failures; the primary counts are 647,615, 3,504 and 11. Missing or untestable SuSiE results are not negative tests.",
    "This partial untestability is distinct from a complete absence of finite results. Full-portfolio nonfinite states comprise 921,663 ABF-fallback rows, 5,016 ABF-only rows and the 25 shared-posterior failures; the primary counts are 647,615, 3,504 and 11. Missing or untestable SuSiE results are not negative tests.\n\n"
    "As a sensitivity analysis, we applied the same shared-posterior check to the 12,087 primary negative rows carried forward from the first run. For each row, we rebuilt the GWAS-fit variant set from the GWAS harmonization and the saved eQTL fit and computed the posterior share of each eQTL credible set on those variants. On the 816 primary rows of the rerun, this reproduced the eligibility of all 5,362 signal pairs. Of the 12,087 rows, 343 (2.8%, in 201 genes) failed on every signal pair and became untestable, 1,744 failed on some pairs and kept their call, because the best PP.H4 over the remaining pairs cannot exceed the first-run value of at most 0.5. The other 10,000 passed on every pair. Failures were more frequent in East Asian rows (166 of 1,178) than in European rows (129 of 9,826). Finite primary SuSiE rows fell from 12,892 to 12,549, and genes with a finite primary result from 3,918 to 3,877. Rows that passed on every pair were assumed to keep their first-run PP.H4 and PP.H3 and were not refit; rows that failed on some pairs have an unknown PP.H3 after the check. Under the priority rule above, 188 more genes became untestable, and the Figure 1B universe fell from 3,769 to 3,589 genes. When a gene counted as tested if any tier-1/2 study stayed eligible, 41 genes became untestable, and the universe fell from 3,776 to 3,736 genes; this rule also counts the seven genes untestable in the rerun as tested. No positive call changed. Overlap results under both rules are in Supplementary Note 3.",
    ["B-G1"],
    [AN1 + ":rows", AN1 + ":diagnostics.primary_genes_with_all_sets_fail_row 201", AN1 + ":diagnostics.primary_status_by_ancestry (EAS 736+276+166=1178; EUR 8395+1302+129=9826)",
     AN1 + ":diagnostics.finite_susie_rows_master_tier12 12892 / finite_susie_rows_after_tier12 12549", AN1 + ":pp3.baseline.finite 3918 / carried_PP.H3.finite 3877",
     AN1R + ": Validation table (5,362 signal pairs, 0 different; 816 rows)", AN1 + ":states_rule_a.transitions (188); H57.a_base.N 3769; H57.a_sens.N 3589",
     AN1 + ":states_rule_b.transitions (41); H57.b_base.N 3776; H57.b_sens.N 3736"],
    "Replace the 435/1,233 upper-bound statement with the measured symmetric check, both gene rules.", "D1")

add("gen2-meth-2", METH,
    "Because carried negatives were not rechecked for the shared posterior, 3,769 is an upper bound on the genes tested at PP.H4 ≤ 0.5 or positive.",
    "This is the prespecified universe. The universes after the shared-posterior check of carried negatives, 3,589 genes under the priority rule and 3,736 when any eligible tier-1/2 study counts, are sensitivity analyses (Supplementary Note 3).",
    ["B-G1"],
    [AN1 + ":H57.a_sens.N 3589; H57.b_sens.N 3736"],
    "D1: 3,769 stays prespecified; point to the measured sensitivity.", "D1")

add("gen2-supp-1", SUPP,
    "The annotation-adjudicated Gene Catalog has a separate 14,682-gene interface, and we did not substitute it for the universe of the native RNA comparison.",
    "The annotation-adjudicated Gene Catalog has a separate 14,682-gene interface, and we did not substitute it for the universe of the native RNA comparison.\n\n"
    "The colocalization rerun applied the shared-posterior check only to first-run positives. We applied the same check to the 12,087 primary negative rows carried forward from the first run (Methods). Under the adopted priority rule, in which one untestable tier-1/2 study makes a gene without a positive call untestable, 188 more genes became untestable and the universe fell to 3,589 genes. It held 426 colocalized genes and 209 disease-state calls, with 26 in both against 21.4 expected within expression deciles (Mantel-Haenszel odds ratio 1.28, 95% CI 0.83 to 1.99; exact *p* = 0.28). When a gene counts as tested if any tier-1/2 study stays eligible, 41 genes became untestable (40 of them among the 14,704 genes) and the universe was 3,736 genes, with 217 disease-state calls and 26 in both against 21.1 expected (odds ratio 1.30, 95% CI 0.84 to 2.02; exact *p* = 0.27). These are nominal sensitivity values for the one prespecified test, and 3,769 remains the prespecified universe. Under the priority rule, 147 genes drawn in Figure 3G would become untestable (139 Neither, 8 Transcriptomics only). The Gene Catalog was not rebuilt; under the priority rule, 169 of its unresolved genes and 11 established state-associated genes would move to the untested class.",
    ["B-G1"],
    [AN1 + ":states_rule_a.transitions (188)", AN1 + ":H57.a_sens (N 3589, G 426, D 209, O 26, E_strat 21.396, MH_OR 1.2815, MH_lo 0.8262, MH_hi 1.9878, p_exact 0.2770)",
     AN1 + ":states_rule_b.transitions (41); H57.b_sens.n_joint_untestable 40", AN1 + ":H57.b_sens (N 3736, D 217, E_strat 21.146, MH_OR 1.3003, MH_lo 0.8384, MH_hi 2.0168, p_exact 0.2746)",
     AN1 + ":fig3g.sens.by_state_rule_a_untestable (Neither 139, Transcriptomics 8)", AN1 + ":catalog_moves.a.by_primary_evidence_class (unresolved 169, established_state_associated 11)",
     AN1R + ": 'In r9d, all 9,319 genetics-untestable genes have class untested' (basis for the move to untested; Catalog builder not rerun)"],
    "Give both gene-level rules with exact values; adopted rule (a) as required.", "D1")

add("gen2-paper-1", PAPER,
    "or tested at PP.H4 ≤ 0.5; an upper bound, because first-run negatives were\n  not rechecked for the shared posterior), the 426 colocalized and 220 canonical\n  disease-state genes share 26, against 21.21 expected within 10 AveExpr\n  quantile bins (Mantel-Haenszel OR 1.295, 95% CI 0.835 to 2.007; exact\n  two-sided P 0.2755; H57, one prespecified test).",
    "or tested at PP.H4 ≤ 0.5; the prespecified universe), the 426 colocalized and\n  220 canonical disease-state genes share 26, against 21.21 expected within 10\n  AveExpr quantile bins (Mantel-Haenszel OR 1.295, 95% CI 0.835 to 2.007; exact\n  two-sided P 0.2755; H57, one prespecified test). In the sensitivity that\n  applies the shared-posterior check to the carried negatives and counts a gene\n  as tested if any tier-1/2 study stays eligible, 41 genes become\n  untestable, the universe is 3,736, and 26 shared genes stand against 21.15\n  expected (OR 1.300, 95% CI 0.838 to 2.017; exact P 0.2746). Under the\n  adopted priority rule the universe is 3,589 (OR 1.282, 95% CI 0.826 to 1.988;\n  exact P 0.2770).",
    ["B-G1"],
    [AN1 + ":H57.b_sens (N 3736, E_strat 21.146, MH 1.3003 [0.8384, 2.0168], p 0.2746)", AN1 + ":states_rule_b.transitions (41)", AN1 + ":H57.a_sens (N 3589, MH 1.2815 [0.8262, 1.9878], p 0.2770)"],
    "D1: replace the upper-bound parenthetical with the measured sensitivity.", "D1")

# ---------------- 2. three-way lead-LD rule ----------------
add("gen2-fig2-4", FIG2,
    "In at most 108 of them, the same accessible region also colocalizes with the gene's eQTL; this caQTL-eQTL comparison applied no lead-variant LD rule, so 108 is an upper count. The remaining genes, at least 35, are reached only through the GWAS-eQTL pair.",
    "In 108 of them, the same accessible region also colocalizes with the gene's eQTL. This comparison uses no lead-variant LD rule; when the accessibility-QTL and eQTL lead variants must also have r² ≥ 0.8, 100 of the 108 remain. The other 35 genes (43 under the LD rule) are reached only through the GWAS-eQTL pair.",
    ["B-G2"],
    [AN5TW + ":M1_under_r2_ge_0.8 (100)", AN5 + ": section 3 (M2 35 -> 43, arithmetic 143-100)"],
    "Replace 'at most 108 ... upper count' with the measured lead-LD result.")

add("gen2-disc-2", DISC,
    "In at most 108 Catalog genes, the same\nregion also colocalizes with the gene's eQTL (an upper count, because this\ncomparison applied no lead-variant LD rule).",
    "In 108 Catalog genes, the same region also\ncolocalizes with the gene's eQTL, and in 100 when the accessibility-QTL and\neQTL lead variants must have r² ≥ 0.8.",
    ["B-G2"],
    [AN5TW + ":M1_under_r2_ge_0.8 (100)"],
    "Replace the upper-count wording with the measured lead-LD result.")

add("gen2-meth-3", METH,
    "The lead-variant LD rule was not applied to this comparison, so the 108 linked genes are an upper count. The GWAS signal was matched across GWAS refits by its lead variant, so for 14 of the 108 genes every qualifying GWAS-caQTL entry comes from another gene's refit that shares the same GWAS lead. For two genes, NARS2 and MTMR3, the colocalized eQTL SuSiE effect has no marginal eQTL association in the Broadaway summary statistics: the smallest marginal p in the credible set is 0.18 for NARS2, and the one variant in the MTMR3 set has p = 0.096.",
    "This comparison used no lead-variant LD rule and links 108 genes. As a sensitivity, we required LD r² ≥ 0.8 between the caQTL lead, the matched variant with the smallest caQTL p, and the eQTL lead, the variant with the largest posterior in the colocalized eQTL effect, using the PolyFun UK Biobank European LD panel of the eQTL fits. Of the 108 genes, 100 passed, and 107 passed at r² ≥ 0.5. The GWAS signal was matched across GWAS refits by its lead variant, so for 14 of the 108 genes every qualifying GWAS-caQTL entry comes from another gene's refit that shares the same GWAS lead; all 14 pass the LD rule. For two genes, NARS2 and MTMR3, the colocalized eQTL SuSiE effect has no marginal eQTL association in the Broadaway summary statistics: the smallest marginal p in the credible set is 0.18 for NARS2, and the one variant in the MTMR3 set has p = 0.096. Both pass the LD rule (r² 0.904 and 0.989), so the rule does not remove them.",
    ["B-G2"],
    [AN5TW + ":M1_under_r2_ge_0.8 (100), M1_under_r2_ge_0.5 (107)", AN5 + ": section 3 (14 no-own-refit genes all pass; NARS2 r2 0.904, MTMR3 r2 0.989; PolyFun UKB EUR panel; eQTL lead = top-alpha variant of idx2)"],
    "Replace the upper-count statement; keep the 14-gene and NARS2/MTMR3 statements consistent with the LD result.")

add("gen2-paper-2", PAPER,
    "with at most 108 Catalog genes where",
    "with 108 Catalog genes (100 when the caQTL and eQTL\n  lead variants must have r² ≥ 0.8) where",
    ["B-G2"], [AN5TW + ":M1_under_r2_ge_0.8 (100)"],
    "Replace the round-1 'at most 108' with the measured lead-LD result.")

add("gen2-paper-3", PAPER,
    "In at most 108 Catalog genes the same region",
    "In 108 Catalog genes (100 with caQTL-eQTL lead r² ≥ 0.8) the same region",
    ["B-G2"], [AN5TW + ":M1_under_r2_ge_0.8 (100)"],
    "Replace the round-1 'at most 108' with the measured lead-LD result.")

# ---------------- 3. PP.H4 = 1.0 diagnostic ----------------
add("gen2-meth-4", METH,
    "No Catalog class was changed.",
    "No Catalog class was changed. We also refit the 18 tier-1/2 rows with PP.H4 = 1.0 that carry the only support of 14 genes, with the LD diagonal, the residual variance and the variants fitted varied one at a time (Supplementary Note 1).",
    ["C2-new"],
    [AN5 + ": section 2 (18 rows, 14 genes; only rows with PP.H4 > 0.5 for these genes)"],
    "Point Methods to the PP.H4 = 1.0 sensitivity.")

add("gen2-supp-2", SUPP,
    "The Catalog classes use p12 = 5 × 10⁻⁶ and are unchanged.",
    "The Catalog classes use p12 = 5 × 10⁻⁶ and are unchanged.\n\n"
    "Fourteen colocalized genes rest only on tier-1/2 rows with SuSiE PP.H4 = 1.0 (18 rows). We refit these rows with the inputs and seed of the colocalization run, which reproduced every PP.H4. All 18 fits converged, and PP.H4 stayed at 1.0 with 1 × 10⁻⁶ instead of 1 × 10⁻³ added to the LD diagonal and with the residual variance estimated. In 14 of the 18 rows, the GWAS variant that drives the colocalized signal does not fit the LD panel: its z score differs from the value predicted from nearby variants through LD by more than 5.45 standardized units (the χ² 1-df cut at p = 5 × 10⁻⁸), against at most 0.92 in three comparison rows. After variants above this cut were removed, 11 of the 14 genes no longer had PP.H4 > 0.5 (BRCA1 kept no GWAS credible set), and GALK1 (0.550), PKN3 and TMPRSS2 kept it. All 11 are genetically anchored genes supported by MVP European GWAS. We chose this filter after the first diagnostic, so it is a post hoc sensitivity. It may also remove true signals at loci with several strong associations, which we did not test, and it changes no prespecified class.",
    ["C2-new"],
    [AN5 + ": section 2 bottom line and verdict (18 rows, 14 genes; 14 of 18 |z_std_diff| > 5.45; 11 lose support; GALK1 0.550, PKN3, TMPRSS2 survive; all 11 MVP EUR, genetically anchored in r9d)",
     AN5D + " (pp4_primary_refit = pp4_r2_release in all rows; comparator |z_std_diff| 0.92, 0.42, 0.59)"],
    "Add the PP.H4 = 1.0 sensitivity as a QC line; name it post hoc; no class change.")

# ---------------- 4. tie rule ----------------
add("gen2-meth-5", METH,
    "The prespecified high-resolution support threshold was SuSiE PP.H4 > 0.5.",
    "The prespecified high-resolution support threshold was SuSiE PP.H4 > 0.5. Each colocalized gene took the trait class of the tier-1/2 GWAS with its highest SuSiE PP.H4. When PP.H4 was exactly equal in more than one GWAS, the GWAS whose study name sorts first alphabetically was used. Only MAU2 is affected: its PP.H4 is 1 in the MVP European ALT, AST and NAFLD GWAS, so it is counted under the enzyme class, and a rule that preferred direct traits would give 422 enzyme and 14 direct genes. Counts of support from any trait class do not use this rule.",
    ["A2.1", "B-G3"],
    [AN5 + ": section 4 (rebuild_tier12.py:77-79 better(); MAU2 driver MVP_ALT_EUR; direct-first split 422/14; any-class counts do not use the tie rule)"],
    "State the declared tie rule (AN5 part 4).")

# ---------------- checks ----------------
def cur(f):
    p = os.path.join(ROOT, RUN, "staged/round1", f)
    return open(p if os.path.exists(p) else os.path.join(ROOT, f)).read()
probs = []
spans = {}
for e in E:
    t = cur(e["file"])
    n = t.count(e["old"])
    if n != 1:
        probs.append(f"{e['id']}: old count {n}")
        continue
    s = t.index(e["old"]); spans.setdefault(e["file"], []).append((s, s + len(e["old"]), e["id"]))
    if "—" in e["new"]:
        probs.append(f"{e['id']}: em dash")
    if "TREAT" in e["new"] or "causal gene" in e["new"].lower():
        probs.append(f"{e['id']}: banned term")
for f, sp in spans.items():
    sp.sort()
    for a, b in zip(sp, sp[1:]):
        if b[0] < a[1]:
            probs.append(f"overlap {a[2]} {b[2]}")
print("edits", len(E), "problems", probs)
out = dict(edits=E, panels=[], notes=open(os.path.join(ROOT, RUN, "analyses/R2/gen2/notes.txt")).read().strip())
json.dump(out, open(os.path.join(ROOT, RUN, "staged/main_gen2.json"), "w"), indent=1, ensure_ascii=False)
