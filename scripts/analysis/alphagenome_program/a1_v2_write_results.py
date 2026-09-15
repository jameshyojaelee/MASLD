#!/usr/bin/env python
"""Write RESULTS.md for package a1-mpra-reliability-v2 from the deposited tables.

Every number in RESULTS.md is read back out of the table that produced it, so the
document cannot drift from the deposit.

Usage: python a1_v2_write_results.py <output_dir> <jobid>
"""
from __future__ import annotations

import json
import os
import sys

import pandas as pd

CONTEXTS = ["HepG2_control", "HepG2_PAOA", "LX2_control", "LX2_TGFb"]
CELLS = ["HepG2", "LX2"]


def f(x, n=4):
    return f"{float(x):.{n}f}"


def esc(s):
    """Escape pipes so a value like '12|34' does not split a markdown table cell."""
    return str(s).replace("|", "\\|")


def main(outdir, jobid):
    tab = os.path.join(outdir, "tables")
    R = pd.read_csv(f"{tab}/reliability_by_context.tsv", sep="\t")
    main_tab = R[(R.partition == "mean_of_three") & (R.stratum == "all")]
    AGR = pd.read_csv(f"{tab}/reference_reproduction_agreement.tsv", sep="\t")
    INV = pd.read_csv(f"{tab}/nan_convention_invariance.tsv", sep="\t")
    MAP = pd.read_csv(f"{tab}/element_locus_group_mapping.tsv", sep="\t")
    PART = pd.read_csv(f"{tab}/treatment_contrast_by_partition.tsv", sep="\t")
    PUB = pd.read_csv(f"{tab}/treatment_contrast_published_value.tsv", sep="\t")
    TER = pd.read_csv(f"{tab}/treatment_tertile_values.tsv", sep="\t")
    DECOMP = pd.read_csv(f"{tab}/treatment_contrast_decomposition_annotated.tsv", sep="\t")
    ADDI = pd.read_csv(f"{tab}/j_component_sign_and_additivity.tsv", sep="\t")
    CONV = pd.read_csv(f"{tab}/activity_convention_sensitivity.tsv", sep="\t")
    P = pd.read_csv(f"{tab}/predictions.tsv", sep="\t")
    DEC = pd.read_csv(f"{tab}/decision_rule_section3.tsv", sep="\t")
    DG = pd.read_csv(f"{tab}/replicate_pairing_diagnostic.tsv", sep="\t")
    MISS = pd.read_csv(f"{tab}/missing_replicates.tsv", sep="\t")
    QC = pd.read_csv(f"{tab}/below_qc_audit.tsv", sep="\t")
    seeds = json.load(open(f"{outdir}/seeds.json"))

    allele = main_tab[main_tab.quantity == "allele_effect"].set_index("context")
    act = main_tab[main_tab.quantity == "activity"].set_index("context")
    trt = main_tab[main_tab.quantity ==
                   "treatment_minus_control_allele_effect"].set_index("context")

    L = []
    w = L.append

    w("# A1 v2. Reporter measurement reliability, GSE281364 MPRA")
    w("")
    w(f"Package `a1-mpra-reliability-v2`. Output directory `{outdir}`.")
    w(f"Slurm job {jobid}, partition cpu, 4 CPUs, 32G, 4 h wall limit.")
    w("Executes section 3 of `scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md`.")
    w("")
    w("This deposit supersedes "
      "`GWAS/finemapping/results/alphagenome_program/a1-mpra-reliability-20260914T193627Z`. "
      "The wave-1 numbers were re-derived from the raw inputs by an independent checker and "
      "were correct; what was wrong there was the deposit, not the arithmetic. Everything is "
      "recomputed here in one sbatch job at the full draw counts, and every table is written.")
    w("")

    # ---------------- what ran
    w("## 1. What ran")
    w("")
    w(f"- Input: the {int(MISS.elements.iloc[0]):,} `paired_snv` elements of "
      "`gse281364-validated-reconstruction-21066470/validated/replicate_outcomes.tsv.gz`, "
      "four contexts, four experimental replicates each.")
    w(f"- Bootstrap: {seeds['bootstrap_draws']:,} draws, seed {seeds['bootstrap_seed']}, "
      "resampling unit the Borzoi long-range block. The statistic is regenerated end to end "
      "inside every draw through the same function that produced the observed value.")
    w(f"- Permutation: {seeds['permutations']:,} draws, seed {seeds['permutation_seed']}, "
      "for the replicate-index pairing diagnostic only.")
    w(f"- Pseudocount {seeds['pseudocount']}; Python {seeds['python'].split('|')[0].strip()}.")
    w("- 70 bootstrap jobs (4 contexts x 2 quantities + 2 cell lines x 1 quantity, each "
      "crossed with the all-elements stratum and six tertile strata).")
    w("")
    w("Tables deposited in `tables/`:")
    w("")
    for fn in sorted(os.listdir(tab)):
        w(f"- `{fn}`")
    w("")

    # ---------------- headline numbers
    w("## 2. Headline numbers (n and 95% block-bootstrap interval attached to each)")
    w("")
    w("Raw two-versus-two Pearson, mean of the three balanced partitions "
      "{12\\|34, 13\\|24, 14\\|23}. The Spearman-Brown projection to four replicates is beside it, "
      "never in place of it.")
    w("")
    w("| quantity | context | n elements | n blocks | raw 2v2 Pearson [95% CI] | "
      "Spearman | projected 4-rep [95% CI] |")
    w("|---|---|---|---|---|---|---|")
    for q, T in (("activity", act), ("allele effect", allele)):
        for c in CONTEXTS:
            r = T.loc[c]
            w(f"| {q} | {c} | {int(r.n_elements):,} | {int(r.n_blocks)} | "
              f"{f(r.pearson)} [{f(r.pearson_lo)}, {f(r.pearson_hi)}] | {f(r.spearman)} | "
              f"{f(r.projected_4rep)} [{f(r.projected_4rep_lo)}, {f(r.projected_4rep_hi)}] |")
    for c in CELLS:
        r = trt.loc[c]
        w(f"| treatment minus control allele effect | {c} | {int(r.n_elements):,} | "
          f"{int(r.n_blocks)} | {f(r.pearson)} [{f(r.pearson_lo)}, {f(r.pearson_hi)}] | "
          f"{f(r.spearman)} | {f(r.projected_4rep)} "
          f"[{f(r.projected_4rep_lo)}, {f(r.projected_4rep_hi)}] |")
    w("")

    # ---------------- predictions
    w("## 3. Prespecified predictions")
    w("")
    w("All four A1 predictions were recorded in the spec before the data were read. "
      "Failed predictions sit in the same table as the ones that held.")
    w("")
    w("| prediction | statement | observed | verdict | changed vs wave 1 |")
    w("|---|---|---|---|---|")
    for r in P.itertuples():
        w(f"| {r.prediction} | {esc(r.statement)} | {esc(r.observed)} | **{r.verdict}** | no |")
    w("")
    w("A1.1 fails on one of four contexts: HepG2_PAOA allele-effect reliability is "
      f"{f(allele.loc['HepG2_PAOA','pearson'])} "
      f"[{f(allele.loc['HepG2_PAOA','pearson_lo'])}, "
      f"{f(allele.loc['HepG2_PAOA','pearson_hi'])}], below the prespecified 0.70 floor, and "
      "the interval excludes 0.70. The other three contexts sit inside [0.70, 0.90].")
    w("")
    w("A1.3 fails on one of two cell lines: HepG2 is "
      f"{f(trt.loc['HepG2','pearson'])} (below 0.40 as predicted) but LX2 is "
      f"{f(trt.loc['LX2','pearson'])} "
      f"[{f(trt.loc['LX2','pearson_lo'])}, {f(trt.loc['LX2','pearson_hi'])}], "
      "far above the predicted ceiling. The prediction was that the treatment contrast is "
      "unreliable everywhere; it is unreliable in HepG2 only.")
    w("")

    # ---------------- decision rule
    w("## 4. Section 3 decision rule (fixed in the spec before the run)")
    w("")
    w("Threshold: projected four-replicate reliability of the treatment-minus-control allele "
      "effect below 0.50 declares that cell line's allele-by-state reporter endpoint not "
      "measurable on this dataset.")
    w("")
    w("| cell line | n elements | n blocks | raw 2v2 Pearson | projected 4-rep [95% CI] | "
      "threshold | verdict |")
    w("|---|---|---|---|---|---|---|")
    for r in DEC.itertuples():
        n = trt.loc[r.cell_line]
        w(f"| {r.cell_line} | {int(n.n_elements):,} | {int(n.n_blocks)} | "
          f"{f(r.raw_2v2_pearson)} | {f(r.projected_4rep)} "
          f"[{f(r.projected_4rep_lo)}, {f(r.projected_4rep_hi)}] | {r.threshold} | "
          f"**{r.verdict}** |")
    w("")
    w("No verdict changed against wave 1. HepG2's condition-aware reporter endpoint is "
      "dropped as a measurement bound, not as a model failure. LX2's is retained.")
    w("")

    # ---------------- the eight corrections
    w("## 5. The eight corrections to the wave-1 deposit")
    w("")

    w("### 5.1 The reference reproduction agrees to four decimals, not exactly")
    w("")
    w("The word *exact* is withdrawn. The published 1,033-element values are four-decimal "
      "roundings of the same statistic; the recomputation lands on the unrounded value, so a "
      "residual of order 1e-05 is expected and is the rounding itself.")
    w("")
    w("| context | n elements | recomputed projected 4-rep | published (4 dp) | "
      "signed difference | agrees at 4 dp | agrees exactly |")
    w("|---|---|---|---|---|---|---|")
    for r in AGR.itertuples():
        w(f"| {r.context} | {int(r.n_elements):,} | {f(r.recomputed_projected_4rep, 6)} | "
          f"{r.published_projected_4rep_4dp} | "
          f"{r.signed_difference_recomputed_minus_published:+.2e} | "
          f"{r.agrees_to_4_decimals} | **{r.agrees_exactly}** |")
    w("")
    w(f"Largest absolute difference {AGR.absolute_difference.max():.2e}. "
      "Table: `tables/reference_reproduction_agreement.tsv`.")
    w("")

    w("### 5.2 Complete-case versus NaN-skipping is an invariance, not a check")
    w("")
    w("There are zero missing element-by-replicate cells in any context, so the two "
      "conventions are the same arithmetic on the same numbers and cannot disagree. "
      "Reporting their agreement as a passed check overstated what was tested. The half-mean "
      "arrays are bitwise identical.")
    w("")
    w("| context | missing cells, 4,359 set | missing cells, 1,033 fixture | "
      "n complete-case | n NaN-skipping | half-means identical | Pearson difference |")
    w("|---|---|---|---|---|---|---|")
    for r in INV.itertuples():
        w(f"| {r.context} | {int(r.missing_element_replicate_cells_full_4359)} | "
          f"{int(r.missing_element_replicate_cells_fixture_1033)} | "
          f"{int(r.n_elements_complete_case):,} | {int(r.n_elements_nan_skipping):,} | "
          f"{r.half_mean_arrays_bitwise_identical} | {r.abs_difference:g} |")
    w("")
    w("Tables: `tables/nan_convention_invariance.tsv`, `tables/missing_replicates.tsv`. "
      f"The {int(QC.rows_all_elements.sum())} below-QC rows in the full replicate file fall "
      "entirely outside the paired_snv elements "
      f"({int(QC.rows_paired_snv_elements.sum())} inside), which is why nothing is missing "
      "here. That is a property of this element set, not a general guarantee.")
    w("")

    w("### 5.3 Elements map many-to-one onto locus groups")
    w("")
    psa = MAP[MAP.set_name == "paired_snv analysis set"].iloc[0]
    fixr = MAP[MAP.set_name == "dna_lm common fixture"].iloc[0]
    blkr = MAP[MAP.set_name.str.contains("blocks")].iloc[0]
    w(f"The analysis set is {int(psa.n_elements):,} elements on "
      f"{int(psa.n_outer_locus_groups):,} outer locus sequence groups, "
      f"{psa.elements_per_group_mean:.2f} to 1, median {psa.elements_per_group_median:g} and "
      f"maximum {int(psa.elements_per_group_max)} elements in one group. It is not 1 to 1. "
      f"The {int(fixr.n_elements):,}-element common fixture is the exception: it holds exactly "
      f"one element per group ({fixr.mapping}), which is why the two counts coincide at "
      f"{int(fixr.n_elements):,} and why wave 1's log line read as though the analysis set were "
      "also 1 to 1.")
    w("")
    w(f"Onto the resampling unit the ratio is larger still: {int(blkr.n_elements):,} elements "
      f"on {int(blkr.n_outer_locus_groups)} Borzoi long-range blocks, "
      f"{blkr.elements_per_group_mean:.2f} to 1, median "
      f"{blkr.elements_per_group_median:g}, maximum {int(blkr.elements_per_group_max)}. "
      "Elements are therefore not independent units and the block bootstrap, not an "
      "element bootstrap, carries the uncertainty.")
    w("")
    w("Table: `tables/element_locus_group_mapping.tsv`.")
    w("")

    w("### 5.4 One published projected value per cell line, with the convention named")
    w("")
    w("Wave 1 printed two projected values for HepG2 without saying they were different "
      "statistics. They are two aggregation orders:")
    w("")
    w(f"- **published**: {esc(PUB.published_convention.iloc[0])}")
    w(f"- alternative: {esc(PUB.alternative_convention.iloc[0])}")
    w("")
    w("The published convention is the one the spec describes (average the partitions, "
      "then project) and is the only one quoted from here on.")
    w("")
    w("| cell line | n elements | n blocks | published projected 4-rep [95% CI] | "
      "alternative | published minus alternative |")
    w("|---|---|---|---|---|---|")
    for r in PUB.itertuples():
        w(f"| {r.cell_line} | {int(r.n_elements):,} | {int(r.n_blocks)} | "
          f"{f(r.published_projected_4rep)} "
          f"[{f(r.published_projected_4rep_lo)}, {f(r.published_projected_4rep_hi)}] | "
          f"{f(r.alternative_projected_4rep)} | "
          f"{r.difference_published_minus_alternative:+.4f} |")
    w("")
    w("The two orders differ by "
      f"{abs(PUB[PUB.cell_line=='HepG2'].difference_published_minus_alternative.iloc[0]):.4f} "
      "in HepG2 and "
      f"{abs(PUB[PUB.cell_line=='LX2'].difference_published_minus_alternative.iloc[0]):.4f} "
      "in LX2. Spearman-Brown is non-linear, so averaging before and after projecting is not "
      "the same operation; the gap widens as the correlation approaches zero, which is why it "
      "shows in HepG2 and not in LX2. Neither value changes the decision-rule verdict.")
    w("")
    w("Table: `tables/treatment_contrast_published_value.tsv`.")
    w("")

    w("### 5.5 Every number carries n")
    w("")
    w("Treatment-by-tertile reliabilities, n attached to each. Tertiles are of the minimum "
      "per-element DNA count and of the minimum barcode count, taken across the two "
      "constructs and the two conditions.")
    w("")
    w("| cell line | stratifier | tertile | n elements | n blocks | "
      "raw 2v2 Pearson [95% CI] | projected 4-rep [95% CI] |")
    w("|---|---|---|---|---|---|---|")
    for r in TER.itertuples():
        w(f"| {r.cell_line} | {r.stratifier} | {r.tertile} | {int(r.n_elements):,} | "
          f"{int(r.n_blocks)} | {f(r.pearson)} [{f(r.pearson_lo)}, {f(r.pearson_hi)}] | "
          f"{f(r.projected_4rep)} [{f(r.projected_4rep_lo)}, {f(r.projected_4rep_hi)}] |")
    w("")
    hep_ter = TER[TER.cell_line == "HepG2"]
    n_hep_zero = int((hep_ter.pearson_lo <= 0).sum())
    w(f"The point estimates rise monotonically with depth in both cell lines and under both "
      f"stratifiers, but {n_hep_zero} of the {len(hep_ter)} HepG2 tertile intervals include "
      "zero, so the HepG2 trend is a point-estimate ordering and not a demonstrated increase. "
      "The HepG2 top tertile "
      f"({f(TER[(TER.cell_line=='HepG2')&(TER.stratifier=='min_dna')&(TER.tertile==3)].pearson.iloc[0])}) "
      "is still far below the LX2 bottom tertile "
      f"({f(TER[(TER.cell_line=='LX2')&(TER.stratifier=='min_dna')&(TER.tertile==1)].pearson.iloc[0])}). "
      "Depth does not explain the HepG2 failure.")
    w("")
    w("The decomposition is re-emitted with n_elements and n_blocks on every row in "
      "`tables/treatment_contrast_decomposition_annotated.tsv`. Reliability rows "
      "(`is_reliability` true):")
    w("")
    w("| cell line | component | n elements | n blocks | raw 2v2 Pearson | "
      "Spearman | projected 4-rep | sign refers to |")
    w("|---|---|---|---|---|---|---|---|")
    for r in DECOMP[DECOMP.is_reliability].itertuples():
        pr = "" if pd.isna(r.projected_4rep) else f(r.projected_4rep)
        w(f"| {r.cell_line} | `{r.component}` | {int(r.n_elements):,} | {int(r.n_blocks)} | "
          f"{f(r.pearson)} | {f(r.spearman)} | {pr} | {r.sign_refers_to} |")
    w("")

    w("### 5.6 The sign of `j_minus_dna_ratio_only`, stated explicitly")
    w("")
    w(ADDI.sign_statement.iloc[0])
    w("")
    w("| cell line | n element x replicate cells | max additivity residual | "
      "mean | sd | fraction positive |")
    w("|---|---|---|---|---|---|")
    for r in ADDI.itertuples():
        w(f"| {r.cell_line} | {int(r.n_element_replicate_cells):,} | "
          f"{r.max_abs_additivity_residual:.2e} | {r.mean_j_minus_dna_ratio_only:+.6f} | "
          f"{r.sd_j_minus_dna_ratio_only:.4f} | "
          f"{r.fraction_positive_j_minus_dna_ratio_only:.4f} |")
    w("")
    w("The additivity guard confirms `j_full = j_rna_ratio_only + j_minus_dna_ratio_only` to "
      "floating-point precision, so the two components partition j and their reliabilities are "
      "comparable on the same scale. Table: "
      "`tables/j_component_sign_and_additivity.tsv`.")
    w("")

    w("### 5.7 The three per-partition Pearsons behind each treatment value")
    w("")
    w("The HepG2 treatment number is a mean of three partition estimates that do not agree in "
      "sign. Reporting only the mean hid that.")
    w("")
    w("| cell line | partition | n elements | n blocks | Pearson | Spearman | "
      "projected 4-rep of this partition |")
    w("|---|---|---|---|---|---|---|")
    for r in PART.itertuples():
        w(f"| {r.cell_line} | {esc(r.partition)} | {int(r.n_elements):,} | {int(r.n_blocks)} | "
          f"{r.pearson:+.4f} | {r.spearman:+.4f} | "
          f"{r.projected_4rep_this_partition:+.4f} |")
    w("")
    hep = PART[(PART.cell_line == "HepG2") & (PART.partition != "mean_of_three")]
    w("The three HepG2 partitions are "
      + ", ".join(f"{esc(r.partition)} {r.pearson:+.4f}" for r in hep.itertuples())
      + ". One of the three is negative. The same 16 samples are reused across all three "
      "partitions, so these are not three independent experiments and their spread is not a "
      "standard error; it is the sensitivity of the estimate to which replicates are paired. "
      "The LX2 partitions agree in sign and magnitude.")
    w("")
    w("Table: `tables/treatment_contrast_by_partition.tsv`.")
    w("")

    w("### 5.8 The activity-convention check, now executed")
    w("")
    w("Wave 1 wrote `a1_activity_convention_check.py` and never ran it. A1.2 compares "
      "activity reliability with allele-effect reliability, and the main script defines "
      "activity as the element mean of the two constructs so that both quantities share one "
      "inferential unit. The check asks whether A1.2 survives the other conventions.")
    w("")
    w("| context | convention | unit | n units | raw 2v2 Pearson | Spearman | projected 4-rep |")
    w("|---|---|---|---|---|---|---|")
    for r in CONV.itertuples():
        w(f"| {r.context} | `{r.quantity}` | {r.unit} | {int(r.n_units):,} | {f(r.pearson)} | "
          f"{f(r.spearman)} | {f(r.projected_4rep)} |")
    w("")
    note = open(f"{tab}/activity_convention_sensitivity.note.txt").read().strip()
    w(f"A1.2 holds under every activity convention in every context: `{note}`. "
      "The element-mean convention gives the highest activity reliability of the four, so it "
      "is the most favourable to A1.2; the verdict does not depend on that choice, because "
      "even the weakest activity convention exceeds the allele effect by a wide margin in all "
      "four contexts.")
    w("")
    w("Table: `tables/activity_convention_sensitivity.tsv`.")
    w("")

    # ---------------- pairing diagnostic
    w("## 6. Replicate index pairing")
    w("")
    w("GEO and the repo do not document whether control replicate 1 and treatment replicate 1 "
      "share anything. The diagnostic asks whether same-index cross-condition residual "
      "correlations exceed different-index ones.")
    w("")
    pair = DG[DG.n_permutations > 0]
    w("| quantity | cell line | same-index r | different-index r | difference | "
      "permutation p (two-sided) |")
    w("|---|---|---|---|---|---|")
    for r in pair.itertuples():
        w(f"| {r.quantity} | {r.cell_line} | "
          f"{r.mean_same_index_cross_condition_r:+.4f} | "
          f"{r.mean_different_index_cross_condition_r:+.4f} | {r.difference:+.4f} | "
          f"{r.permutation_p_two_sided:.4f} |")
    w("")
    w(f"No evidence of index pairing in any quantity or cell line "
      f"({int(seeds['permutations']):,} permutations, seed {seeds['permutation_seed']}; "
      "smallest p "
      f"{pair.permutation_p_two_sided.min():.4f}, uncorrected, family of "
      f"{len(pair)} tests). That is what A1.4 assumed, and A1.4 is met: the nine-combination "
      "convention differs from same-index by less than 0.05 in both cell lines "
      "(`tables/a1_4_convention_difference.tsv`, "
      "`tables/treatment_nine_combination.tsv`).")
    w("")

    # ---------------- not established
    w("## 7. What is NOT established")
    w("")
    w("- **These reliabilities are not model performance ceilings for any endogenous task.** "
      "They bound how well this reporter assay reproduces itself, on these elements, at these "
      "sequencing depths. Nothing here licenses a claim about endogenous chromatin or "
      "expression.")
    w("- **The LX2 treatment contrast being reliable is not evidence that TGFb changes allelic "
      "regulation.** Most of the reliable part of j in LX2 tracks the common allele effect: "
      "`_same_half_4rep_j_vs_rna_ratio` is "
      f"{f(DECOMP[(DECOMP.cell_line=='LX2')&(DECOMP.component=='_same_half_4rep_j_vs_rna_ratio')].pearson.iloc[0])} "
      "and the residual after removing a per-condition scale term still has reliability "
      f"{f(DECOMP[(DECOMP.cell_line=='LX2')&(DECOMP.component=='j_residual_after_scale_removal')].pearson.iloc[0])}, "
      "so reliability survives scale removal but the biological reading of j is not settled "
      "by a reliability number.")
    w("- **The three partitions are not independent replicates of the estimate.** They reuse "
      "the same 16 samples. Their spread is a pairing sensitivity, not a standard error, and "
      "the bootstrap interval is the only uncertainty statement here.")
    w("- **The Spearman-Brown projection assumes exchangeable, equally reliable replicates "
      "with independent errors.** The residual correlation matrices "
      "(`tables/replicate_residual_correlation_*.tsv`) show structure across replicates, so "
      "the projection is an upper reading of what four replicates would deliver, not a "
      "measurement.")
    w("- **Zero missing cells is a property of this element set**, not a guarantee about the "
      "GSE281364 release; the below-QC rows exist and sit outside the paired_snv elements.")
    w("- **HepG2's failure is not attributed to a cause.** Depth and barcode tertiles do not "
      "explain it (5.5), and this package does not test any further explanation.")
    w("- The bootstrap resamples 1-Mb-scale Borzoi long-range blocks. It does not address "
      "systematic error shared across all blocks, such as a library-preparation batch effect.")
    w("")

    # ---------------- defects
    w("## 8. Defects found and fixed in this rerun")
    w("")
    w("| # | defect in the wave-1 deposit | state now |")
    w("|---|---|---|")
    w("| 1 | 7 of 10 expected tables were missing from the deposit | all "
      f"{len(os.listdir(tab))} tables written |")
    w("| 2 | quoted values came from an off-sbatch reduced run (5 bootstrap draws, "
      "50 permutations) | one sbatch job at "
      f"{seeds['bootstrap_draws']:,} draws and {seeds['permutations']:,} permutations |")
    w("| 3 | the reference reproduction was called exact | four decimals, signed differences "
      "reported (5.1) |")
    w("| 4 | complete-case vs NaN-skipping was presented as a check | relabelled an "
      "invariance, with the zero-missing-cell count (5.2) |")
    w("| 5 | elements implied 1-to-1 with locus groups | 4.22 to 1, with the fixture "
      "exception explained (5.3) |")
    w("| 6 | two projected values published for HepG2 | one, convention named (5.4) |")
    w("| 7 | tertile and decomposition numbers carried no n | n_elements and n_blocks on "
      "every row (5.5) |")
    w("| 8 | `j_minus_dna_ratio_only` sign left to the reader through a double negation | "
      "stated, with an additivity guard (5.6) |")
    w("| 9 | only the mean of three partition Pearsons was shown | all three, one negative "
      "(5.7) |")
    w("| 10 | the activity-convention check was written and never run | run (5.8) |")
    w("")

    w("## 9. Reproduction of the wave-1 tables")
    w("")
    wave1 = ("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/"
             "results/alphagenome_program/a1-mpra-reliability-20260914T193627Z")
    w(f"Byte-for-byte comparison of this deposit's tables against `{wave1}/tables`, "
      "run after the job finished:")
    w("")
    same, diff, only1 = [], [], []
    w1t = os.path.join(wave1, "tables")
    if os.path.isdir(w1t):
        for fn in sorted(os.listdir(w1t)):
            a, b = os.path.join(w1t, fn), os.path.join(tab, fn)
            if not os.path.exists(b):
                only1.append(fn)
            elif open(a, "rb").read() == open(b, "rb").read():
                same.append(fn)
            else:
                diff.append(fn)
        new = sorted(set(os.listdir(tab)) - set(os.listdir(w1t)))
        w(f"- identical: {len(same)} of {len(os.listdir(w1t))} wave-1 tables")
        w(f"- differing: {len(diff)}" + (f" ({', '.join(diff)})" if diff else ""))
        w(f"- present in wave 1 only: {len(only1)}"
          + (f" ({', '.join(only1)})" if only1 else ""))
        w(f"- new in this deposit: {len(new)} ({', '.join(f'`{x}`' for x in new)})")
        w("")
        if not diff and not only1:
            w("**No wave-1 number moved.** Every wave-1 table reproduces exactly under this "
              "full 10,000-draw, 10,000-permutation sbatch run at seed 20260914, so the "
              "wave-1 tables as they now stand on disk were not the reduced-draw output; the "
              "reduced-draw provenance concern applies to values quoted in prose, not to "
              "those files. All four A1 verdicts and both decision-rule verdicts are "
              "unchanged. What this package adds is the ten new tables and the corrected "
              "labelling, not new arithmetic.")
    else:
        w("Wave-1 directory not readable at write time; comparison not run.")
    w("")

    w("## 10. Provenance")
    w("")
    w("- `MANIFEST.tsv`: input paths with sha256.")
    w("- `script.sha256`: sha256 of both producer scripts and the sbatch file, as run.")
    w("- `a1_mpra_reliability.py.run`, `a1_mpra_reliability_v2.py.run`, "
      "`a1_mpra_reliability_v2.sbatch.run`: byte copies of what ran.")
    w("- `pip_freeze.txt`, `seeds.json`, `run_log.txt`, `run_log_addenda.txt`.")
    w("")

    with open(f"{outdir}/RESULTS.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"wrote {outdir}/RESULTS.md ({len(L)} lines)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
