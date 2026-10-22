# Blind checker for gen2: rewrite the writer's draft edits into the corrected final list.
import json, os, re
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN = ROOT + "/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z"
A1 = "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/gene_level_20261021T165252Z/"
A5 = "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN5/"
CHK = "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/R2/gen2/check_fig3g_counts_out.txt"
PR = "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/R2/gen2/pair_rates_out.txt"

draft = json.load(open(RUN + "/staged/main_gen2.json"))
edits = {e["id"]: e for e in draft.get("edits", draft.get("final_edits"))}
order = [e["id"] for e in draft.get("edits", draft.get("final_edits"))]

def upd(i, new, srcs=None, reason=None):
    e = edits[i]; e["new"] = new
    if srcs is not None: e["number_sources"] = srcs
    if reason: e["reason"] = reason

upd("gen2-fig1-1",
    "The rerun applied this check only to first-run positives. In a sensitivity analysis, we also applied it to the first-run negatives carried forward and counted a gene as tested if any tier-1/2 study stayed eligible. This rule counts the seven genes above as tested, so it starts from 3,776 genes. After the check, 40 of them became untestable, leaving 3,736 genes, and 26 genes were in both sets against 21.1 expected (Mantel-Haenszel OR 1.30, 95% CI 0.84 to 2.02; exact *p* = 0.27). The check cannot add or remove a colocalized gene, and 3,769 remains the prespecified universe (Supplementary Note 3).",
    [A1 + "results.json:H57.b_base.N 3776",
     CHK + ": rule b crosstab, tested_le_0.5 -> untestable 40 among the 3,776",
     A1 + "results.json:H57.b_sens (N 3736, O 26, E_strat 21.146, MH_OR 1.3003, MH_lo 0.8384, MH_hi 2.0168, p_exact 0.2746)"],
    "D1: replace the upper-bound wording with the measured rule-(b) sensitivity; state the 3,776 starting universe so 3,769 to 3,736 is not read as 41 losses.")

upd("gen2-fig2-1",
    "As a sensitivity analysis, we applied the check to the 12,087 carried-forward negative gene-by-GWAS tests in the tier-1/2 studies. In 343 of them (2.8%), every signal pair failed, so the test became untestable; in 1,744, only some signal pairs failed and the call stayed the same; the other 10,000 passed (Methods). A carried negative that fails the check becomes untestable, not positive, so the colocalized genes do not change.",
    [A1 + "results.json:rows (primary 12087, all_sets_fail 343, some_sets_fail 1744, none_fail 10000)"],
    "Replace the upper-bound sentence with the measured result; say gene-by-GWAS tests, not rows.")

upd("gen2-fig2-2",
    "Most pairs in these denominators are first-run negatives carried forward. When the shared-posterior check was also applied to them, the denominators fell to 11,797 and 752 pairs, and the rates rose to 6.3% and 3.7%. In the first run, before any check, the rates were 787 of 12,132 (6.5%) and 29 of 771 (3.8%).",
    [PR + ": sensitivity tier 2 finite 11797 pos 739 (6.26%); tier 1 finite 752 pos 28 (3.72%); rerun by the checker, identical output; computed from " + A1 + "master_reduced_sensitivity.csv"])

upd("gen2-fig2-3",
    "One gene, MAU2, has PP.H4 of 1.0 in the MVP European ALT, AST and NAFLD GWAS. Exact ties go to the GWAS whose study identifier sorts first alphabetically, so MAU2 is counted under the enzyme class (Methods).",
    [A5 + "RESULT.md: section 4 (MAU2 PP.H4.susie exactly 1 in MVP_ALT_EUR, MVP_AST_EUR, MVP_NAFLD_EUR; driver MVP_ALT_EUR)",
     A5 + "pph4_1_rows_all_t12.tsv (MAU2 rows: PP.H4.susie 1.0 only in these three)"])

upd("gen2-fig2-4",
    "In 108 of them, the same accessible region also colocalizes with the gene's eQTL. This comparison uses no lead-variant LD rule. When the accessibility-QTL and eQTL lead variants must also have r² ≥ 0.8, 100 of the 108 remain (107 at r² ≥ 0.5). The other 35 genes (43 under the r² ≥ 0.8 rule) are reached only through the GWAS-eQTL pair.",
    [A5 + "threeway/run_20261021T125538Z/threeway_summary.json:M1_under_r2_ge_0.8 (100), M1_under_r2_ge_0.5 (107)",
     A5 + "RESULT.md: section 3 (M2 35 -> 43, arithmetic 143 - 100)"],
    "Replace 'at most 108 ... upper count' with the measured lead-LD result; add the r² ≥ 0.5 count because the same paragraph uses r² ≥ 0.5 for the GWAS-caQTL step.")

upd("gen2-fig3-1",
    "The seven SuSiE-untestable genes are drawn as open circles at their best finite PP.H4 from another tier-1/2 study and are not counted as Neither. The panel shows states before the shared-posterior check of carried-forward negatives (Methods). After that check, 40 drawn genes, one of them an open circle, keep no finite tier-1/2 PP.H4 and would leave the panel (36 Neither and 3 Transcriptomics only). If a gene counts as tested when any tier-1/2 study stays eligible, no drawn gene is untestable, and the four groups hold 3,119, 191, 400 and 26 genes.",
    [CHK + ": no finite PP.H4 after check 40 (Neither 37, Transcriptomics 3), of which rerun-untestable (open circles) 1, so 36 Neither",
     A1 + "results.json:fig3g.dropped_no_finite_pp4 (Neither 37, Transcriptomics 3)",
     A1 + "results.json:fig3g.sens.by_state (3119/191/400/26)",
     "figures/candidates/panel-deep-fixes-20261017T035633Z/fig3/source_tables/fig3g.tsv (7 'SuSiE untestable' rows drawn at susie_pp4, equal to " + A1 + "fig3g_states.tsv coloc_best_susie_pp4_base)"],
    "State where the open circles sit and what the check does to the drawn genes; the 37 dropped Neither genes include one open circle (ASAP3).")

upd("gen2-leg-1",
    "Open circles, seven SuSiE-untestable genes, drawn at their best finite PP.H4 from another tier-1/2 study and not counted as Neither. States precede the shared-posterior check of carried-forward negatives (Supplementary Note 3).",
    ["figures/candidates/panel-deep-fixes-20261017T035633Z/fig3/source_tables/fig3g.tsv (susie_pp4 of the 7 'SuSiE untestable' rows equals " + A1 + "fig3g_states.tsv coloc_best_susie_pp4_base)"])

upd("gen2-disc-1",
    "When the shared-posterior check was also applied\nto the carried negatives, 3,877 genes kept a finite Tier-1/2 SuSiE result and\n3,306 showed this PP.H3 pattern. This count keeps the first-run PP.H3 of\ngene-study results in which only some signal pairs failed. Their PP.H3 after\nthe check is unknown, so the count lies between 3,101 and 3,416. Under the\nCatalog's rule for untestable genes, 3,148 of the 3,306 would be called\nindeterminate (range 2,948 to 3,252).",
    [A1 + "results.json:pp3.carried_PP.H3 (finite 3877, pp3 3306, indeterminate_rule_a 3148)",
     A1 + "results.json:pp3.lower_partial_not_counted (pp3 3101, indeterminate_rule_a 2948)",
     A1 + "results.json:pp3.upper_partial_counted (pp3 3416, indeterminate_rule_a 3252)"],
    "Replace the upper-bound wording with the measured PP.H3 counts, including the 3,319 indeterminate count the old text also qualified.")

upd("gen2-meth-1",
    "This partial untestability is distinct from a complete absence of finite results. Full-portfolio nonfinite states comprise 921,663 ABF-fallback rows, 5,016 ABF-only rows and the 25 shared-posterior failures; the primary counts are 647,615, 3,504 and 11. Missing or untestable SuSiE results are not negative tests.\n\nAs a sensitivity analysis, we applied the same shared-posterior check to the 12,087 primary negative rows carried forward from the first run. For each row, we rebuilt the GWAS-fit variant set from the GWAS harmonization and the saved eQTL fit and computed the posterior share of each eQTL credible set on those variants. On the 816 primary rows of the rerun, this reproduced the eligibility of all 5,362 signal pairs. Of the 12,087 rows, 343 (2.8%, in 201 genes) failed on every signal pair and became untestable. Another 1,744 failed on some pairs and kept their call, because the best PP.H4 over the remaining pairs cannot exceed the first-run value of at most 0.5. The other 10,000 passed on every pair. Failures were more frequent in East Asian rows (166 of 1,178) than in European rows (129 of 9,826). Finite primary SuSiE rows fell from 12,892 to 12,549, and genes with a finite primary result from 3,918 to 3,877. Rows that passed on every pair were assumed to keep their first-run PP.H4 and PP.H3 and were not refit; rows that failed on some pairs have an unknown PP.H3 after the check. Under the priority rule above, 188 more genes became untestable (180 of them among the 3,769 genes of Figure 1B), and the Figure 1B universe fell to 3,589 genes. Under a rule in which a gene counts as tested if any tier-1/2 study stays eligible, the seven genes untestable in the rerun count as tested, so the universe before the check is 3,776 genes. After the check, 41 genes became untestable under this rule (40 of them among the 3,776), and the universe fell to 3,736 genes. No positive call changed. Overlap results under both rules are in Supplementary Note 3.",
    [A1 + "results.json:rows (primary 12087, all_sets_fail 343, some_sets_fail 1744, none_fail 10000)",
     A1 + "results.json:diagnostics.primary_genes_with_all_sets_fail_row 201",
     A1 + "results.json:diagnostics.primary_status_by_ancestry (EAS 736+276+166=1178; EUR 8395+1302+129=9826)",
     A1 + "results.json:diagnostics.finite_susie_rows_master_tier12 12892 / finite_susie_rows_after_tier12 12549",
     A1 + "results.json:pp3.baseline.finite 3918 / pp3.carried_PP.H3.finite 3877",
     "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/RESULT.md: Validation table (0 of 5,362 signal pairs differ; 816 rows)",
     A1 + "results.json:states_rule_a.transitions (188); H57.a_sens.N 3589",
     CHK + ": rule a crosstab tested_le_0.5 -> untestable 180 among the 3,769; rule b crosstab 40 among the 3,776",
     A1 + "results.json:states_rule_b.transitions (41); H57.b_base.N 3776; H57.b_sens.N 3736"],
    "Replace the 435/1,233 upper-bound statement with the measured check under both gene rules, with the gene counts that reconcile each universe.")

upd("gen2-meth-2",
    "This is the prespecified universe. Supplementary Note 3 repeats the test after the shared-posterior check of carried negatives, on 3,589 genes under the priority rule and on 3,736 genes when a gene counts as tested if any tier-1/2 study stays eligible.",
    [A1 + "results.json:H57.a_sens.N 3589; H57.b_sens.N 3736"])

upd("gen2-supp-1",
    "The annotation-adjudicated Gene Catalog has a separate 14,682-gene interface, and we did not substitute it for the universe of the native RNA comparison.\n\nThe colocalization rerun applied the shared-posterior check only to first-run positives. We applied the same check to the 12,087 primary negative results carried forward from the first run (Methods). Under the priority rule used for the 3,769 genes, in which one untestable tier-1/2 study makes a gene without a positive call untestable, 188 more genes became untestable, 180 of them among the 3,769, and the universe fell to 3,589 genes. It held 426 colocalized genes and 209 disease-state calls, with 26 in both against 21.4 expected within expression deciles (Mantel-Haenszel odds ratio 1.28, 95% CI 0.83 to 1.99; exact *p* = 0.28). When a gene counts as tested if any tier-1/2 study stays eligible, the seven SuSiE-untestable genes count as tested, and the universe is 3,776 genes before the check. After the check, 40 of these genes became untestable and the universe was 3,736 genes, with 217 disease-state calls and 26 in both against 21.1 expected (odds ratio 1.30, 95% CI 0.84 to 2.02; exact *p* = 0.27). These are nominal sensitivity values for the one prespecified test, and 3,769 remains the prespecified universe. In Figure 3G, 40 drawn genes would keep no finite PP.H4 under either rule. Under the priority rule, 141 more drawn genes would become untestable, giving 147 open circles (139 Neither, 8 Transcriptomics only). The Gene Catalog was not rebuilt. Under the priority rule, 169 of its unresolved genes and 11 of its established-state associated genes would move to the untested class.",
    [A1 + "results.json:states_rule_a.transitions (188)",
     CHK + ": rule a crosstab 180 among the 3,769; rule b crosstab 40 among the 3,776; rule a drawn untestable after 147 (139 Neither, 8 Transcriptomics), 6 already untestable, so 141 new; no finite PP.H4 after check 40",
     A1 + "results.json:H57.a_sens (N 3589, G 426, D 209, O 26, E_strat 21.396, MH_OR 1.2815, MH_lo 0.8262, MH_hi 1.9878, p_exact 0.2770)",
     A1 + "results.json:H57.b_base.N 3776",
     A1 + "results.json:H57.b_sens (N 3736, D 217, O 26, E_strat 21.146, MH_OR 1.3003, MH_lo 0.8384, MH_hi 2.0168, p_exact 0.2746)",
     A1 + "results.json:fig3g.sens.by_state_rule_a_untestable (Neither 139, Transcriptomics 8); fig3g.dropped_no_finite_pp4 (37+3)",
     A1 + "results.json:catalog_moves.a.by_primary_evidence_class (unresolved 169, established_state_associated 11)",
     "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/RESULT.md: 'In r9d, all 9,319 genetics-untestable genes have class untested' (basis for the move; Catalog builder not rerun)"],
    "Give both gene rules with values that reconcile (180 of 188; 3,776 start for rule b); 147 open circles include the 6 existing ones; use the paper's class name 'established-state associated'.")

upd("gen2-paper-1",
    "or tested at PP.H4 ≤ 0.5; the prespecified universe), the 426 colocalized and\n  220 canonical disease-state genes share 26, against 21.21 expected within 10\n  AveExpr quantile bins (Mantel-Haenszel OR 1.295, 95% CI 0.835 to 2.007; exact\n  two-sided P 0.2755; H57, one prespecified test). In the sensitivity that\n  applies the shared-posterior check to the carried negatives and counts a gene\n  as tested if any tier-1/2 study stays eligible (3,776 genes before the check,\n  because the seven SuSiE-untestable genes count as tested), 40 genes become\n  untestable, the universe is 3,736, and 26 genes are shared against 21.15\n  expected (OR 1.300, 95% CI 0.838 to 2.017; exact P 0.2746). Under the\n  priority rule of the 3,769 universe it is 3,589 (OR 1.282, 95% CI 0.826 to\n  1.988; exact P 0.2770).",
    [A1 + "results.json:H57.b_base.N 3776",
     CHK + ": rule b crosstab 40 among the 3,776",
     A1 + "results.json:H57.b_sens (N 3736, E_strat 21.146, MH 1.3003 [0.8384, 2.0168], p 0.2746)",
     A1 + "results.json:H57.a_sens (N 3589, MH 1.2815 [0.8262, 1.9878], p 0.2770)"])

upd("gen2-meth-3",
    "This comparison used no lead-variant LD rule and links 108 genes. As a sensitivity analysis, we required LD r² ≥ 0.8 between the caQTL lead and the eQTL lead in the PolyFun UK Biobank European LD panel used for the eQTL fits. The caQTL lead was the matched variant with the smallest caQTL p, and the eQTL lead was the variant with the largest posterior in the colocalized eQTL effect. Of the 108 genes, 100 passed, and 107 passed at r² ≥ 0.5. The GWAS signal was matched across GWAS refits by its lead variant, so for 14 of the 108 genes every qualifying GWAS-caQTL entry comes from another gene's refit that shares the same GWAS lead; all 14 pass the LD rule. For two genes, NARS2 and MTMR3, the colocalized eQTL SuSiE effect has no marginal eQTL association in the Broadaway summary statistics: the smallest marginal p in the credible set is 0.18 for NARS2, and the one variant in the MTMR3 set has p = 0.096. Both pass the LD rule (r² 0.904 and 0.989), so the rule does not remove them.")

upd("gen2-meth-4",
    "No Catalog class was changed. We also refit the 18 tier-1/2 GWAS fits with SuSiE PP.H4 of 1.0 that carry the only support of 14 colocalized genes, varying the LD diagonal, the residual variance and the variants fitted one at a time (Supplementary Note 1).")

upd("gen2-meth-5",
    "The prespecified high-resolution support threshold was SuSiE PP.H4 > 0.5. Each colocalized gene took the trait class of the tier-1/2 GWAS with its highest SuSiE PP.H4. When PP.H4 was exactly equal in more than one GWAS, the GWAS whose study identifier sorts first alphabetically was used. MAU2 is the only gene whose tied maximum spans two trait classes: its PP.H4 is 1 in the MVP European ALT, AST and NAFLD GWAS, so it is counted under the enzyme class, and a rule that preferred direct traits would give 422 enzyme and 14 direct genes. Counts of support from any trait class do not use this rule.",
    [A5 + "RESULT.md: section 4 (rebuild_tier12.py:77-79 better(); 3 genes tied at their maximum, MAU2, WHRN, THBS3, only MAU2 spans two classes; driver MVP_ALT_EUR; direct-first split 422/14; any-class counts do not use the tie rule)"])

upd("gen2-supp-2",
    "The Catalog classes use p12 = 5 × 10⁻⁶ and are unchanged.\n\nFourteen colocalized genes rest only on tier-1/2 GWAS fits with SuSiE PP.H4 of 1.0 (at least 0.999999; 18 gene-by-GWAS tests). We refit these tests with the inputs and seed of the colocalization run, which reproduced every PP.H4. All 18 fits converged, and PP.H4 stayed at 1.0 with 1 × 10⁻⁶ instead of 1 × 10⁻³ added to the LD diagonal and with the residual variance estimated. In 14 of the 18 tests, the GWAS variant that drives the colocalized signal does not fit the LD panel: its z score differs from the value predicted from nearby variants through LD by more than 5.45 standardized units (the χ² 1-df cut at p = 5 × 10⁻⁸), against at most 0.92 in three comparison tests (RORA with MVP European NAFLD, FABP1 and MAU2 with UK Biobank ALT). After variants above this cut were removed, 11 of the 14 genes no longer had PP.H4 > 0.5 (BRCA1 kept no GWAS credible set), and GALK1 (0.550), PKN3 and TMPRSS2 kept it. All 11 are genetically anchored genes supported by MVP European GWAS. We chose this filter after the first diagnostic, so it is a post hoc sensitivity analysis. It may also remove true signals at loci with several strong associations, which we did not test, and it changes no prespecified class.",
    [A5 + "RESULT.md: section 2 (rows with PP.H4.susie >= 0.999999; 18 rows, 14 genes; 14 of 18 |z_std_diff| > 5.45; 11 lose support; GALK1 0.550, PKN3, TMPRSS2 survive; all 11 MVP EUR, genetically anchored in r9d; comparators RORA MVP_NAFLD_EUR, FABP1 UKBB_ALT, MAU2 UKBB_ALT)",
     A5 + "pph4_1_diagnostic_summary.tsv (pp4_primary_refit = pp4_r2_release in the 18 pph4_1 rows; z_std_diff_at_best_hit1 > 5.45 in 14 of 18; comparator 0.92, 0.42, 0.59)",
     A5 + "job_index.tsv (filter (e) run as second-pass job 3143 after job 3141)"])

final = [edits[i] for i in order]

# checks
spans = {}
for e in final:
    p = RUN + "/staged/round1/" + e["file"]
    src = p if os.path.exists(p) else ROOT + "/" + e["file"]
    t = open(src).read()
    assert t.count(e["old"]) == 1, e["id"]
    assert "—" not in e["new"] and "–" not in e["new"], e["id"]
    assert "TREAT" not in e["new"], e["id"]
    i = t.find(e["old"]); spans.setdefault(e["file"], []).append((i, i + len(e["old"]), e["id"]))
for fn in ["main_fig3c.json", "main_misc.json"]:
    fp = RUN + "/staged/" + fn
    if not os.path.exists(fp): continue
    o = json.load(open(fp))
    for e in o.get("final_edits", o.get("edits", [])):
        p = RUN + "/staged/round1/" + e["file"]
        src = p if os.path.exists(p) else ROOT + "/" + e["file"]
        t = open(src).read(); i = t.find(e["old"])
        spans.setdefault(e["file"], []).append((i, i + len(e["old"]), fn + ":" + e["id"]))
for f, s in spans.items():
    s.sort()
    for a, b in zip(s, s[1:]):
        assert b[0] >= a[1], ("overlap", f, a, b)
json.dump({"final_edits": final, "panels": []}, open(RUN + "/staged/main_gen2.json", "w"), ensure_ascii=False, indent=1)
print("wrote", len(final), "edits; all old strings unique; no overlaps with", [k for k in ["main_fig3c.json", "main_misc.json"]])
