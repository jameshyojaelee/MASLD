"""Build and check the bulk text-writer edit list (REVIEW2 triage).

Checks: each `old` occurs exactly once in its file; no two edits overlap in
the same file. Writes staged/edits_bulk.json. Read-only on the manuscript.
"""
import json
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z"
ID = "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified"
SD = "figures/candidates/resource-corrections-20261005-v1/identity_resource_strict_depth"
FIT = "Analysis/Hypotheses/refit-adoption-review-20261011T151011Z/out"
U3C = "figures/candidates/review-r3-fixes-20261015T112148Z/F3cdg/R03-02_gse126848_depth"
H82 = "Analysis/Hypotheses/H82-20261015T125052Z"
H56 = "Analysis/Hypotheses/H56-20261015T112556Z/RESULT_H56.md"
H79 = "Analysis/Hypotheses/H79-20261015T125125Z/RESULT_H79.md"
NMF = "RNA-seq/results/manuscript_release/candidates/fixed-nmf-corrected-scores-20261009T141312Z"
SYS = "RNA-seq/results/histology_anchored_continuum/molecular_systems_stress/hac-system-stress-corrected774-20261009T171017Z/tables/system_redundancy.tsv"
DESIGN = "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/model_design.tsv"

M = "docs/manuscript/METHODS.md"
F3 = "docs/manuscript/draft/fig3.md"
F4 = "docs/manuscript/draft/fig4.md"
SN = "docs/manuscript/draft/supplementary_notes.md"
LG = "docs/manuscript/05_figure_legends.md"
PP = "docs/PAPER.md"

E = []


def add(id_, file, old, new, items, sources, reason, decision=None):
    E.append({"id": id_, "file": file, "old": old, "new": new, "items": items,
              "number_sources": sources, "decision": decision, "reason": reason})


# ---------------- A3.1 identity ----------------
add("bulk-01", M,
    "We did not perform a complete cross-study check for overlapping participants.",
    "Genotypes were called from the RNA-seq reads at 8,660 autosomal sites, keeping calls with read depth of at least 10. "
    "Libraries with fewer than 1,000 called sites were not used. Two libraries were joined into one genetic-identity cluster "
    "when KING-robust kinship (PLINK 2.0) was at least 0.354 over at least 500 shared sites and also in each of two disjoint "
    "halves of the sites (at least 250 shared sites in each half). A cluster that contained an unresolved pair or conflicting "
    "source labels was excluded. Kinship was computed for every pair of the 1,466 libraries from ten cohorts in one joint panel, "
    "so the screen covered pairs within and across studies. All 72 pairs at or above 0.354 were within one cohort. Each pair of "
    "the 774 selected libraries shared at least 598 sites, and no pair reached 0.354. One pair from different cohorts, both "
    "disease units (GSE213621 and GSE130970), had kinship 0.288 (0.294 and 0.282 in the two halves) and opposite homozygous "
    "genotypes at 54 of 2,845 shared sites, against at most 6 in any pair at or above 0.354. This fits first-degree relatives, "
    "not one person, and both libraries were kept. No other pair of selected libraries reached kinship 0.0884 in PLINK. "
    "Because no library pair from different cohorts reached 0.354, the screen has no cross-cohort positive example, and KING "
    "cannot distinguish monozygotic twins from repeated libraries of one person.",
    ["A3.1"],
    [f"{ID}/analysis_definition.json:thresholds (genotype_depth 10, full_called_sites 1000, full_pair_shared_sites 500, each_half_pair_shared_sites 250, same_genetic_identity_king 0.354)",
     f"{ID}/filtered_fingerprint_commands.log:24,26,79 (1466 samples, 8660 variants, 114 relationships + 1073731 filtered = C(1466,2))",
     f"{ID}/analysis_unit_census.tsv all_screened rows (10 cohorts); all_screened_roster.tsv 1,466 rows, 10 datasets (recomputed)",
     f"{ID}/full_filtered_related_pairs.tsv: 72 rows KINSHIP > 0.354, 0 cross-cohort, IBS0 max 6 (recomputed)",
     f"{ID}/full_filtered_related_pairs.tsv row SRR21622926/SRR9036348: NSNP 2845, IBS0 54, KINSHIP 0.28815; half0 0.293532, half1 0.282456",
     "inputs/C6-identity.json critic.challenges[A3.1].rederivation (min shared sites 598 over 299,151 pairs of the 774)",
     f"{ID}/resource_identity_roster.tsv: SRR9036348 group_binary Disease, SRR21622926 Disease",
     "scripts/analysis/fig34_corrections/11_resource_identity.py:128-152, 312-314 (half-split rule, unresolved-component and site-minimum exclusion)",
     f"{ID}/analysis_definition.json:limitations (KING cannot distinguish monozygotic twins)"],
    "METHODS.md:113 says no cross-study check was done; the all-pairs KING screen exists and the clustering rule was never stated.")

add("bulk-02", M,
    "A separate identity sensitivity (seed 20260929) excluded ambiguous calls and every group with conflicting labels, and kept "
    "the lexicographically first eligible library in each concordant same-person component. It changed no donor roster or gene membership.",
    "An earlier identity sensitivity on the previous 844-record roster (seed 20260929) excluded ambiguous calls and every group "
    "with conflicting labels, and kept the lexicographically first eligible library in each concordant same-person component. "
    "It was evaluated only on continuum associations and is not a sensitivity of the 774-unit differential expression. A "
    "strict-depth identity sensitivity repeated the screen described below with read depth of at least 20 and at least 500 shared "
    "sites in each half of the sites. It keeps 738 clusters and is a four-cohort disease contrast, because GSE126848 keeps no "
    "disease unit under it (Supplementary Note 7).",
    ["A3.1"],
    ["docs/technical/FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md:16-31 (844-record identity_sensitivity, continuum associations)",
     f"{SD}/analysis_definition.json:thresholds (genotype_depth 20, each_half_pair_shared_sites 500)",
     "docs/RESULTS.md:86,129 (738 units; GSE126848 keeps 5 controls and no disease unit)"],
    "The 'changed no donor roster or gene membership' sentence describes the old 844-record sensitivity; the 774-unit strict-depth sensitivity changes membership.")

add("bulk-03", F3,
    "These clusters are not verified independent donors.",
    "These clusters are not verified independent donors, but no two selected libraries, within or across cohorts, had genotypes "
    "consistent with one person (Methods).",
    ["A3.1"],
    [f"{ID}/full_filtered_related_pairs.tsv (0 pairs >= 0.354 among the 774; recomputed)"],
    "States the cross-cohort identity result where n = 774 is introduced.")

add("bulk-04", SN,
    "and up calls fall from 1,333 to 1,085 or 1,100.",
    "and up calls fall from 1,333 to 1,085 or 1,100.\n\n"
    "A strict-depth identity screen (read depth of at least 20 and at least 500 shared sites in each half of the sites) keeps 738 "
    "clusters. The primary model refitted on them tests 24,509 genes and gives 1,906 calls (1,344 higher and 562 lower in disease). "
    "Of the 1,854 primary calls, 1,779 are also called there, and none has the opposite sign. GSE126848 keeps 5 controls and no "
    "disease unit under this screen, so the refit is a four-cohort disease contrast.",
    ["A3.1"],
    [f"{SD}/analysis_definition.json:thresholds",
     f"{FIT}/04_call_counts_all_fits.tsv row 738_primary_stored (24509, 1906, 1344, 562)",
     f"{FIT}/03_04_fit_comparisons.tsv row 774_primary vs 738_primary (both 1779, both_opposite_sign 0)",
     "docs/RESULTS.md:129 (GSE126848 keeps 5 controls and no disease unit)"],
    "Reports the unreported strict-depth identity sensitivity, which changes gene membership.")

# ---------------- A3.2 / B-B4 tissue-level ----------------
add("bulk-05", F3,
    "It gave 1,854 canonical calls at Benjamini-Hochberg (BH) padj < 0.05 across all 24,115 genes and |log2FC| > 0.5.",
    "It gave 1,854 canonical calls at Benjamini-Hochberg (BH) padj < 0.05 across all 24,115 genes and |log2FC| > 0.5. "
    "These are tissue-level differences, so a call can reflect a change in cell-type proportions, a change in expression within "
    "cells, or both.",
    ["A3.2", "B-B4"],
    [],
    "H56 prespecified fallback reading assumes the calls are described as tissue-level; the text never said so (C7 critic missed item).")

# ---------------- B-B2 Fig 3C (D3) ----------------
add("bulk-06", F3,
    "The union of per-cohort calls contains 13,761 genes: 9,065 called in one cohort and 4,696 in two or more. Only 2,831 of the "
    "4,696 have the same sign in every cohort that calls them (1,933 in two cohorts, 636 in three, 218 in four and 44 in all five). "
    "The other 1,865 (39.7%) are called up in one cohort and down in another, and 1,632 of these include a GSE126848 call. "
    "Without GSE126848, 346 of 2,565 multi-cohort genes (13.5%) have mixed signs. These 20,054 cohort memberships are a different "
    "set from the 1,854 pooled calls, and shared membership is not replication.",
    "Because GSE126848 calls track sequencing depth (next paragraph), we summarize overlap in the other four cohorts and show "
    "GSE126848 as a separate row. The union of their per-cohort calls contains 10,630 genes: 8,065 called in one cohort and 2,565 "
    "in two or more. Of the 2,565, 2,219 have the same sign in every cohort that calls them (1,686 in two cohorts, 439 in three and "
    "94 in all four). The other 346 (13.5%) are called up in one cohort and down in another. These 13,907 cohort memberships are a "
    "different set from the 1,854 pooled calls, and shared membership is not replication. Adding GSE126848 gives a 13,761-gene "
    "union in which 1,865 of 4,696 multi-cohort genes (39.7%) have mixed signs, and 1,632 of these include a GSE126848 call. "
    "Removing only the 437 GSE126848 calls on genes with no reads in any of its 30 libraries still leaves 1,784 of 4,604 (38.7%) "
    "mixed in a 13,447-gene union.",
    ["B-B2"],
    [f"{U3C}/fig3c_union_under_gse126848_sensitivities.tsv rows 'four cohorts without GSE126848' (13907, 10630, 8065, 2565, 346, 1686, 439, 94) and 'GSE126848 zero-read calls removed' (13447, 4604, 1784)"],
    "Four-cohort union becomes the primary Fig 3C overlap; GSE126848 reported separately with the zero-read-masked union. Filtered-refit union (12,180; 31.5%) held for AN3.",
    "D3")

add("bulk-07", F3,
    "Of the 6,147 calls, 437 are such genes, with zero reads in all 30 libraries.",
    "Of the 6,147 calls, 437 are such genes, with zero reads in all 30 libraries. More broadly, 2,730 of the 6,147 calls are on "
    "genes that the cohort's own expression filter removes, and 2,639 of these are higher in disease.",
    ["B-B2"],
    [f"{U3C}/gse126848_depth_sensitivity_genes.tsv.gz: base_call True and filterByExpr_keep False = 2730; base_logFC > 0 = 2639 (recomputed)"],
    "Shows the low-count problem extends beyond the 437 zero-read genes and is directional.")

add("bulk-08", LG,
    "Middle, genes per overlap class in the 13,761-gene union; darker bars include a GSE126848 call. Bottom, mixed-sign share of "
    "multi-cohort genes, with and without GSE126848.",
    "Middle, genes per overlap class in the 10,630-gene union of the four cohorts other than GSE126848; GSE126848 is shown as a "
    "separate row because its calls track library depth. Bottom, mixed-sign share of multi-cohort genes: 346 of 2,565 (13.5%) in "
    "four cohorts and 1,865 of 4,696 (39.7%) with GSE126848 added.",
    ["B-B2"],
    [f"{U3C}/fig3c_union_under_gse126848_sensitivities.tsv"],
    "Legend for the D3 layout; must match the redrawn Fig 3C panel.",
    "D3")

add("bulk-09", PP,
    "The per-cohort union contains 13,761 genes, including 9,065 called in one cohort and 4,696 in two or more, of which 2,831 "
    "have one sign in every calling cohort (44 in all five) and 1,865 have opposite signs in different cohorts; it is distinct "
    "from the 1,854 pooled calls.",
    "Per-cohort overlap is summarized in the four cohorts other than GSE126848: their union contains 10,630 genes, including 8,065 "
    "called in one cohort and 2,565 in two or more, of which 2,219 have one sign in every calling cohort and 346 (13.5%) have "
    "opposite signs in different cohorts; it is distinct from the 1,854 pooled calls. GSE126848 is reported separately because "
    "its calls track library depth (five-cohort union 13,761 genes; 1,865 of 4,696 multi-cohort genes with mixed signs).",
    ["B-B2"],
    [f"{U3C}/fig3c_union_under_gse126848_sensitivities.tsv"],
    "PAPER.md Fig 3 row matches the D3 primary union.",
    "D3")

add("bulk-10", PP,
    "shared per-cohort membership is not replication, because 1,865 of 4,696 multi-cohort genes change sign between cohorts.",
    "shared per-cohort membership is not replication, because 346 of 2,565 multi-cohort genes in the four cohorts other than "
    "GSE126848 change sign between cohorts.",
    ["B-B2"],
    [f"{U3C}/fig3c_union_under_gse126848_sensitivities.tsv row 'four cohorts without GSE126848'"],
    "PAPER.md claim column matches the D3 primary union.",
    "D3")

# ---------------- B-B3 ----------------
add("bulk-11", F3,
    "When we removed each cohort and refit the full model, 66.7% to 96.8% of calls kept the same sign, and removing GSE213621 "
    "changed the most calls.",
    "When we removed each cohort and refit the full model, 66.7% to 96.8% of the 1,854 calls were still called with the same "
    "sign, and among calls still tested, 99.7% to 100% kept their sign. Removing GSE213621 lost the most calls (618, of which 108 "
    "were no longer tested). Fitted alone, GSE213621 calls 1,654 of the 1,854 pooled calls (89%) with the same sign, and the next "
    "cohort, GSE162694, calls 1,014.",
    ["B-B3"],
    [f"{FIT}/05_loco_recomputed.tsv: pct_kept 66.7-96.8; sign_agree_all_full_calls_tested 0.997136-1.0; GSE213621 lost 618, lost_not_tested 108",
     "figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv + primary_source_controls/deg_results.csv: GSE213621 1,654 same-sign of 1,854; GSE162694 1,014 (recomputed, FDR < 0.05 and |logFC| > 0.5)"],
    "'Kept the same sign' described call retention; sign agreement is 0.997-1.000. Adds the single-cohort dependence.")

add("bulk-12", LG,
    "Full refits that repeat filtering after cohort removal retain 66.7%–96.8% of calls with the same sign.",
    "Full refits that repeat filtering after cohort removal keep 66.7%–96.8% of calls as same-sign calls; among calls still "
    "tested, 99.7%–100% keep their sign.",
    ["B-B3"],
    [f"{FIT}/05_loco_recomputed.tsv sign_agree_all_full_calls_tested"],
    "Same correction in the Figure S3C legend.")

add("bulk-13", PP,
    "Removing one cohort and refitting keeps 88.6% to 96.8% of the 1,854 calls with the same sign for four cohorts and 66.7% "
    "(1,236) for GSE213621, which supplies",
    "Removing one cohort and refitting keeps 88.6% to 96.8% of the 1,854 calls as same-sign calls for four cohorts and 66.7% "
    "(1,236) for GSE213621, which alone calls 1,654 of the 1,854 with the same sign and supplies",
    ["B-B3"],
    [f"{FIT}/05_loco_recomputed.tsv", "cohort_disease_all_gene_results.tsv + deg_results.csv (1,654; recomputed)"],
    "PAPER.md wording matches fig3.md.")

# ---------------- A3.2 / B-B4 F1 heading and paragraph (A6) ----------------
add("bulk-14", F3,
    "## Expression differs from F0 at every fibrosis stage, and the F1 difference depends mainly on one cohort",
    "## Expression differs from F0 at F2 to F4; the F1 difference depends on one cohort and on the F0 reference",
    ["A3.2", "B-B4"],
    [],
    "F1 rests on 1,338 calls that fall to 14 without GSE130970 and to 783 with a disease-only F0. 'Library properties' omitted: under H82, F2 and NAS 1-2 are also 'severe', so that clause would be F1-unique wording.")

add("bulk-15", F3,
    "and 0 (GSE162694, 65 units) genes (Figure S3 [panel letter to be assigned]).",
    "and 0 (GSE162694, 65 units) genes (Figure S3 [panel letter to be assigned]). With F0 restricted to its 89 disease biopsies "
    "(191 units), F1 gives 783 genes (Figure S3E).",
    ["A3.2", "B-B4"],
    [f"{H82}/run-20261015T131530Z/d4_f0_reference.tsv row base_diseaseF0 (191 units, 783 calls)",
     "docs/manuscript/05_figure_legends.md:125 (89 disease-only F0 units, 783)"],
    "Supports the heading's F0-reference clause in the main text.")

add("bulk-16", F3,
    "The Figure 3E counts are unchanged.",
    "The Figure 3E counts are unchanged. As a post hoc description, Supplementary Note 7 gives the calls in each stage contrast "
    "that stay significant with these properties added.",
    ["B-B4"],
    [],
    "Points to the H82 descriptive counts reported as a logged deviation.",
    "D4")

add("bulk-17", M,
    "It did not, so the prespecified stop rule applied and no result from this check is reported.",
    "It fell below in 18 of 20 draws, not in all 20, so the prespecified stop rule applied and the test is not interpreted. As a "
    "post hoc deviation from the prespecified reading, the retained counts for each contrast are reported descriptively "
    "(Supplementary Note 7).",
    ["B-B4"],
    [f"{H82}/RESULT_H82.md Controls PC-3 (18 of 20 below the F1 null 2.5th percentile 0.5119)"],
    "PRESPEC_H82.md:398-399 vs :450 conflict; reporting D1 counts is a deviation needing user approval and a DEVIATIONS.md entry.",
    "D4")

add("bulk-18", SN,
    "### Transfer of grouped fibrosis-stage effects to GSE213621",
    "### Library properties of the stage contrasts\n\n"
    "This section is a post hoc descriptive deviation from the prespecified reading of the library-property check, whose positive "
    "control failed (Methods). With within-cohort log2 library size, Gini coefficient and normalization factor added together, the "
    "calls that stay BH-significant are 99 of 1,338 for F1, 2,113 of 4,266 for F2, 4,797 of 6,190 for F3 and 5,307 of 7,043 for F4 "
    "versus F0; 2,595 of 5,338, 7,296 of 9,337 and 7,083 of 9,841 for NAS 1–2, 3–4 and 5–8 versus NAS0; 3,006 of 4,305 for NAS "
    "5–8 versus 3–4; and 356 of 398, 497 of 512 and 114 of 281 for F1 to F2, F2 to F3 and F3 to F4. These counts are not a test. "
    "The three properties are computed from the same counts and may reflect disease-related RNA composition, so a lower count does "
    "not show that a stage difference is technical.\n\n"
    "### Transfer of grouped fibrosis-stage effects to GSE213621",
    ["B-B4", "A3.2"],
    [f"{H82}/RESULT_H82.md primary table 'Joint retained / released' for all 11 contrasts"],
    "H82 joint-adjustment retention for all 11 contrasts, no F1-unique wording (PRESPEC_H82.md:438-442 graded reading).",
    "D4")

add("bulk-19", M,
    "The stromal-specific positive control did not pass, so no per-gene labels were produced.",
    "The stromal-specific positive control did not pass (median ratio 0.528 in one combination). The uniform-change negative "
    "control also failed in seven of the eight combinations (median ratios 0.484 to 0.695), the within-hepatocyte benchmark "
    "failed, and the calibration null exceeded its limit in one combination (lambda 1.254). No per-gene labels were produced.",
    ["B-B4"],
    [f"{H56} Controls table (PC2 0.528 C7_earlier; NC1 0.484-0.695 in the 7 failing cells, 0.705 pass in C7_earlier; NC2 fail)",
     f"{H56} Diagnostic 9 (A16_earlier mean_lambda 1.254, S11_stop TRUE)"],
    "Methods named only one of four control failures.")

# ---------------- A8 ancestry ----------------
add("bulk-20", F3,
    "Restricting the fit to units called as European genetic similarity (EUR) changed 570 calls (Holm-adjusted *p* = 0.0199) "
    "(Supplementary Note 7).",
    "Restricting the fit to units called as European genetic similarity (EUR) changed 570 calls (Holm-adjusted *p* = 0.0199) "
    "(Supplementary Note 7). In each test, one reference fit changed at least as many calls as observed (one of 200 permuted-PC "
    "fits and one of 100 random subsets), so both p values are limited by the number of reference fits.",
    ["A8"],
    [f"{H79}:21-22 (B 200 and 100; #ref >= obs 1 and 1)"],
    "Holm p 0.0199 reflects the permutation resolution, not a precise p value.")

add("bulk-21", SN,
    "Both Holm-adjusted *p* values are 0.0199 over the two tests.",
    "Both Holm-adjusted *p* values are 0.0199 over the two tests. One of the 200 permuted-PC fits and one of the 100 random "
    "subsets changed at least as many calls as observed, giving one-sided *p* = 2/201 (0.00995) and 2/101 (0.0198). The smallest "
    "attainable one-sided values are 1/201 and 1/101, so more reference fits would be needed to estimate these p values precisely.",
    ["A8"],
    [f"{H79}:21-22 (p 0.00995 and 0.0198, Holm 0.0199)", "docs/manuscript/METHODS.md:151 (p = (1 + k)/(B + 1); B = 200 and 100)"],
    "States the resolution floor of the two permutation nulls.")

# ---------------- B-B5 / B-S6 / B-BM ----------------
add("bulk-22", F4,
    "(466-unit descriptive roster)",
    "(all 466 GSE162694 and GSE213621 units, including their 99 controls; descriptive)",
    ["B-B5"],
    [f"{SYS} row pooled_centred (n 466); GSE162694 140, GSE213621 326",
     f"{DESIGN}: controls GSE162694 31, GSE213621 68 (recomputed)"],
    "Defines the 466 roster; it includes 99 controls, unlike the 367-unit coefficients.")

add("bulk-23", F4,
    "System scores use the 24,115 genes tested in the expression data, and memberships, weights and coordinates are unchanged.",
    "Memberships, weights and coordinates were built on an earlier 23,370-gene expression universe and are unchanged. System "
    "scores use the 24,115 genes tested in the current data. Because 392 genes of the earlier universe are not in these data, "
    "four nodes have no coefficient.",
    ["B-B5", "B-S6"],
    ["inputs/C8-versions.json B-S6 evidence and critic rederivation (bulk/continuum_gene_meta.tsv.gz: 392 only in 23,370, 1,137 only in 24,115, 22,978 shared)",
     "docs/RESULTS.md:1757 (map built on the 23,370-gene universe)",
     "figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/repair/fig4g/fig4g_node_checks.tsv (unmatched_main_nodes 4)"],
    "States the geometry universe and why four nodes are untestable.")

add("bulk-24", M,
    "The system redundancy summary used a 466-unit descriptive roster.",
    "The system redundancy summary used all 466 GSE162694 and GSE213621 units (140 and 326), including their 31 and 68 controls. "
    "It is descriptive. The map was built on the earlier 23,370-gene universe, and 392 of those genes are absent from the "
    "24,115-gene family. As a result, four nodes keep less than 80% of their original membership after signature exclusion and "
    "have no coefficient.",
    ["B-B5", "B-S6"],
    [f"{SYS} row pooled_centred", f"{DESIGN} (31 and 68 controls; recomputed)",
     "inputs/C8-versions.json B-S6 (392 genes; testability.tsv reason less_than_80_percent_original_membership_observed_after_signature_exclusion)"],
    "Defines 466 and the four untestable nodes in Methods.")

add("bulk-25", M,
    "The 139 observed genes from the 145-gene continuum signature were removed before member genes were mapped to GENCODE v49 identifiers.",
    "All 145 continuum signature genes were removed before member genes were mapped to GENCODE v49 identifiers. The map was built "
    "on the earlier 23,370-gene expression universe, in which 139 of the 145 genes were observed; none of the 145 appears in the "
    "final membership.",
    ["B-B1", "B-B5"],
    ["inputs/C8-versions.json B-B1 critic (outcome_blind_gene_membership.tsv.gz 0 signature-gene matches; geometry signature_exclusion input 145 rows)",
     "hac-molecular-layers-20260818T173348Z/inputs/signature_mapping.tsv (139 True / 6 False)"],
    "Exclusion covers all 145 genes; names the universe the geometry was built on.")

add("bulk-26", M,
    "Of the 145 signature genes, 139 were observed; the six unavailable genes are listed in the mapping table.",
    "Of the 145 signature genes, 143 are observed in the 24,115-gene family (DPEP1 and ZACN are not). In the earlier 23,370-gene "
    "family used for the molecular-systems map and the expression-relationship extension, 139 were observed, because KRT23, "
    "PRAMEF10, IQGAP3 and IER3 were also absent. The mapping table lists the unavailable genes.",
    ["B-B1"],
    ["hac-molecular-layers-corrected774-20261009T171017Z/inputs/signature_mapping.tsv (143 True / 2 False); bulk/run_summary.tsv signature_observed 143",
     "hac-molecular-layers-20260818T173348Z/inputs/signature_mapping.tsv (139 / 6)"],
    "Fixes the 139/6 print in the general continuum paragraph; 143 is current.")

add("bulk-27", M,
    "The 139 observed axis genes keep their separate leave-one-out results and 139-gene BH family.",
    "On that 23,370-gene basis, the 139 observed axis genes keep their separate leave-one-out results and 139-gene BH family; in "
    "the 24,115-gene fit this family has 143 genes.",
    ["B-B1"],
    ["docs/manuscript/METHODS.md:197 (143-gene family)"],
    "Separates the two leave-one-out families.")

add("bulk-28", LG,
    "Four untestable nodes are counted at right.",
    "Four nodes are untestable in the 24,115-gene data, because fewer than 80% of their original members remain observed after "
    "signature exclusion; they are counted at right.",
    ["B-S6"],
    ["inputs/C8-versions.json B-S6 (testability.tsv; fig4g_nodes_disease_only.tsv.gz 4 NA continuum_beta)"],
    "Gives the reason for the four untestable nodes.")

add("bulk-29", F3,
    "the sixth axis, k6-P6, is unavailable",
    "the sixth axis, k6-P6, is not shown because, with its fixed weights restricted to the 4,517 retained genes, its "
    "reconstruction of the released scores from the original counts had Spearman 0.9874, below the prespecified 0.99 rule",
    ["B-BM"],
    [f"{NMF}/axis_availability.tsv row k6 P6 full_stage (retained_rows 4517, historical_global_spearman 0.98742, 'historical global Spearman 0.9874 < 0.99')"],
    "Gives the reason k6-P6 is missing (critic-corrected: historical-count reconstruction, not corrected counts).")

add("bulk-30", F3,
    "These axes describe continuous variation, not patient classes.",
    "These axes describe continuous variation, not patient classes. On the corrected counts, the five shown axes agree with the "
    "released scores at Spearman 0.928 (k6-P5) to 0.998 (k6-P1) over 617 shared libraries.",
    ["B-BM"],
    [f"{NMF}/corrected_vs_historical_diagnostics.tsv k6 full_stage spearman_vs_historical_H (P1 0.99773, P5 0.92781; n_overlap 617)"],
    "Tells the reader how closely the displayed k6 axes match the released ones.")

add("bulk-31", LG,
    "Bottom, F0-centered activity of five k6 NMF axes",
    "Bottom, F0-centered activity of five k6 NMF axes (k6-P6 fails the fixed reconstruction rule and is not shown)",
    ["B-BM"],
    [],
    "Gives the k6-P6 reason in the Fig 3E legend.")

# ---------------- C8 new (D6) ----------------
add("bulk-32", M,
    "The signal-level direction analysis used the same 50-study colocalization results as the gene-level analysis.",
    "The signal-level direction analysis used the 50-study colocalization results of the earlier 2026-08-17 release (462 SuSiE "
    "genes), not the corrected rerun that gives the 436 SuSiE genes reported elsewhere in this paper.",
    ["C8-new"],
    ["docs/RESULTS.md:3871-3872 (816 gene-study pairs over 462 genes)",
     "docs/STATUS.md:843-845 ('This extension used the 2026-08-17 release (462 SuSiE and 822 ABF genes)')"],
    "The sentence misstated the COLOC release behind the 5,362 / 22,356,697 / 816 counts.",
    "D6")


# ---------------- checks ----------------
problems = []
texts = {}
for e in E:
    p = ROOT / e["file"]
    t = texts.setdefault(e["file"], p.read_text())
    n = t.count(e["old"])
    if n != 1:
        problems.append(f"{e['id']}: old occurs {n} times in {e['file']}")
for f, t in texts.items():
    spans = sorted((t.find(e["old"]), t.find(e["old"]) + len(e["old"]), e["id"]) for e in E if e["file"] == f)
    for (a0, a1, ai), (b0, b1, bi) in zip(spans, spans[1:]):
        if b0 < a1:
            problems.append(f"overlap {ai} / {bi} in {f}")
for e in E:
    if "—" in e["new"]:
        problems.append(f"{e['id']}: em dash")
    if "TREAT" in e["new"] or "causal" in e["new"].lower():
        problems.append(f"{e['id']}: banned term")
print("\n".join(problems) if problems else "all checks pass", len(E), "edits")

held = [
    "Fig 3C filtered-refit union (12,180 genes; 1,211 of 3,839 mixed, 31.5%; 946 involving GSE126848): held until AN3 saves a producer run.",
    "GTEx area-adjusted fat rho (AN6/§4 A9): held; fig3.md:15 not edited.",
    "GSE126848 binomial-thinning sensitivity (C7 critic): not run; no text.",
    "773-unit close-relative refit (drops SRR9036348): not run; no text.",
    "Redundancy summary on the 367 disease units (could change 'about 14' and 41%): not run; fig4.md:45 keeps the 466-unit values with the roster defined.",
]
out = {"writer": "bulk", "edits": E, "held": held}
(RUN / "staged/edits_bulk.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
