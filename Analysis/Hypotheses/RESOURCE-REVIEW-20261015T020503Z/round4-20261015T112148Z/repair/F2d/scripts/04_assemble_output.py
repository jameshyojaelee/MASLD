# Assemble the F2d repair return (entries, check failures, conflicts) from final_text_edits.json.
import json
import sys
from pathlib import Path

OUT = Path(sys.argv[1])
E = json.loads((OUT / "final_text_edits.json").read_text())
C = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F2d/"
R = C + "repair/"
W = ("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/"
     "round4-20261015T112148Z/repair/F2d/")
H58 = ("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H58-20261015T112423Z/jobs/5800/"
       "H58-20261015T113519Z_22252009/")
H57 = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H57-20261015T112421Z/run-20261015T113107Z/"
JOB = ("sbatch job 7530 (slurm 22256138, COMPLETED 0:0): " + W + "scripts/01_rederive.py -> " + W +
       "out/rederive.json")
REN = ("Rscript " + W + "scripts/02_render_repair.R " + C + " " + R + " (light R on saved tables); pdftotext -layout; "
       "pdftoppm -r 150 PNG viewed")
CHK = "python " + W + "scripts/03_build_edits.py " + W + "out (all 32 old_text found exactly once; no overlaps)"

entries = [
 dict(fid="R02-01", status="fixed", reproduced="yes",
      changed_in_repair=True,
      what_changes=(
        "The pooled 23.34% vs 3.87% coding/splice contrast is reported as a pooled descriptive value, not a trait-class "
        "property, in the abstract, intro, Fig 2 Results, Fig 2D legend, Discussion, METHODS, title note, PAPER, STATUS, "
        "RESULTS and both READMEs. Repair changes: (1) fig2.md:103 now says 15 of 19 direct clusters contain a European "
        "direct set (was 16) and 106 of 207 direct credible sets are European; the abstract replaces 'most direct "
        "credible sets are European' with '106 of 207 direct-trait but 82 of 6,208 enzyme credible sets are European'. "
        "(2) The calibration failures are attributed correctly: EUR failed coverage, power and N1a; AMR failed N1b; "
        "EUR+AMR failed N1c and carries the EUR coverage failure; all-non-EUR pooled had no control and is descriptive. "
        "(3) METHODS states the full rule: coverage within twofold, MDE at most 0.195, 95% coverage in at least 180 of "
        "200 and corrected-level coverage in at least 190 of 200 splits, and a call only when the 98.333% lower bound "
        "is above 0. (4) The 8 of 767 lead windows are named as UKBB ALT/AST/GGT; MVP EUR enzyme GWAS were not "
        "fine-mapped. (5) The causal 'because more enzyme signals were tested' clause and the 17.5%/21.1% rates are "
        "removed from the intro and Discussion; they report 12,121 vs 771 tested gene-by-GWAS pairs instead. "
        "(6) Legend ranges use 'to' and U+2212 minus; legend values match the repaired panel in percentage points. "
        "(7) New edits for PAPER 558-560, docs/README.md 23-25 and figures/main/fig2_genetics/README.md 5-12 "
        "(conflict 18) and title.md 22-24 (conflict 5). The ancestry panel is re-rendered with printed differences "
        "and intervals."),
      old_value=("Pooled 23.34% vs 3.87%, ratio 6.026671 [2.574593, 13.510165] presented as a direct-versus-enzyme "
                 "property; first-round text '16 of 19 direct clusters contain a European set'"),
      new_value=("Pooled 23.34% vs 3.87% (descriptive). EUR 29.88% vs 6.92%, difference 0.2296 [0.0550, 0.4236]; "
                 "AMR 6.73% vs 4.54%, 0.0219 [-0.0489, 0.1065]; EUR+AMR 0.1497 [0.0388, 0.2718]; non-EUR pooled "
                 "0.0331 [-0.0383, 0.1274]. No trait-class call. 15 of 19 direct clusters contain a European direct "
                 "set; 106 of 207 direct and 82 of 6,208 enzyme credible sets are European."),
      panels=[R + "Fig2D_coding_by_ancestry_repair.pdf"],
      source_tables=[R + "Fig2D_coding_by_ancestry_repair_source.tsv", H58 + "family_C.tsv", H58 + "G7_coverage.tsv",
                     H58 + "N1_coverage.tsv", W + "out/rederive.json"],
      commands="; ".join([JOB, REN, CHK]),
      note=("Edits apply to the current files; none depends on another edit's new text. The abstract 24-32 edit "
            "(this entry) and the abstract 32-35 edit (R14-03) touch disjoint parts of line 32. The title.md 22-24 edit "
            "supersedes Claims E197 (same old_text); Claims' title line 10-11 edits are compatible with it. Claims "
            "E195/E196 (abstract 24-32, intro 77-84) should be dropped (conflict 8). The intro edit ends at line 84 "
            "'regulation.'; lines 84-89 belong to conflict 3. The pooled row in the ancestry panel uses the adopted "
            "seed-20260817 interval and is labeled; all rows use the same adopted October 13 credible sets, so no "
            "release state is mixed. The PNPLA3/TM6SF2 absence in the UKBB arm was checked in credible_sets.csv (no "
            "UKBB credible set within 0.5 Mb of either gene). N1 split controls were not rerun here; the counts are "
            "from H58 (final). The repair source table no longer contains embedded newlines in row labels."),
      text_edits=E["R02-01"]),
 dict(fid="R14-02", status="fixed", reproduced="yes",
      changed_in_repair=True,
      what_changes=(
        "The direct-subclass values (PDFF 40.14%, NAFLD+NASH 12.44%) move out of the adopted RESULTS table into a "
        "separate labeled paragraph after the bootstrap-table paragraph (RESULTS 2536-2537), so the 'Adopted corrected "
        "result' column holds only adopted values. The paragraph also carries the H58 family C values, calibration "
        "counts and ancestry composition. fig2.md:103 and the Fig 2D legend carry the subclass values through R02-01 "
        "edits. The first-round subclass panel is unchanged."),
      old_value="Two H58 candidate rows inserted inside the adopted RESULTS table (header 'Adopted corrected result')",
      new_value=("PDFF 40.14% (diff 0.3626 [0.1835, 0.5455]); NAFLD+NASH 12.44% (0.0857 [-0.0164, 0.2096]); enzyme "
                 "3.87%; in a separate RESULTS paragraph"),
      panels=[C + "Fig2D_direct_subclass_split.pdf"],
      source_tables=[C + "Fig2D_direct_subclass_source.tsv", W + "out/rederive.json"],
      commands=JOB + "; " + CHK,
      note="The RESULTS note edit also serves R14-06; it is listed once, here.",
      text_edits=E["R14-02"]),
 dict(fid="R14-06", status="fixed", reproduced="yes",
      changed_in_repair=True,
      what_changes=(
        "Concentration (top 5 of 19 direct clusters carry 76.0%; chr11:823,586 22.5%; ratios 4.93 to 1.96 after "
        "dropping direct rows, labeled a biased diagnostic) and the ungated ascertainment values are in fig2.md:103 "
        "(R02-01 edit) and in the separate RESULTS paragraph (R14-02 edit), not in the adopted table. The "
        "ascertainment sentence now names its ungated set (6,467 credible sets, seed 20260817) and drops the causal "
        "'because'; RESULTS notes that the medians include the 52 sets failing the PIP-sum gate."),
      old_value="Concentration row inside the adopted RESULTS table; ascertainment values given without the ungated set",
      new_value=("Top-5 direct clusters 76.0%; matched ratio 5.78 [1.99, 48.05], 27 matched enzyme loci; median n_eff "
                 "21,893 vs 160,000; beta_std 0.0368 vs 0.0183 over 208 and 6,259 ungated sets"),
      panels=[],
      source_tables=[C + "fig2d_direct_cluster_concentration.tsv",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/"
                     "fig2-coordinate-rerun-20261013T213349Z_22218816/ascertainment/ascertainment_control_summary.tsv",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/"
                     "fig2-coordinate-rerun-20261013T213349Z_22218816/ascertainment/power_imbalance.tsv"],
      commands="read ascertainment_control_summary.tsv and power_imbalance.tsv (not re-run); " + CHK,
      note="Text delivered by the R02-01 fig2.md:103 edit and the R14-02 RESULTS edit; no separate edit here.",
      text_edits=[]),
 dict(fid="R02-02", status="fixed", reproduced="yes",
      changed_in_repair=True,
      what_changes=(
        "fig2.md:107 and the companion legend now name the estimator correctly: context posterior summed over all "
        "credible sets of a class divided by summed noncoding posterior (not 'weighting each credible set equally'). "
        "Values: promoter 13.84/14.67, accessible 9.26/8.44, ABC 0.43/3.70, no context 77.90/77.69, unresolved 0/0. "
        "The companion is re-rendered with the estimator printed on its axis. RESULTS 1837-1840 rows are labeled "
        "pre-repair or withdrawn in the label cell (no fourth cell), and the following paragraph is replaced with the "
        "repaired reading. RESULTS 2539 says 'almost all of it in non-European studies' (EUR 0.27 of 3,652 summed "
        "unresolved PIP). RESULTS 2566 and STATUS 75-76 point to the repair companion."),
      old_value=("Adopted unresolved 52.11%/60.43% (summed) and 36.49%/77.85% (cluster-weighted); first-round text "
                 "called the summed estimator 'weighting each credible set equally'"),
      new_value=("Repaired, summed over sets: promoter 13.84%/14.67%, accessible 9.26%/8.44%, ABC 0.43%/3.70%, no "
                 "context 77.90%/77.69%, unresolved 0%/0%. Equal-weight alternative 12.11/14.71, 8.43/8.47, "
                 "0.37/3.67, 80.36/77.62 (not printed)."),
      panels=[R + "Fig2_companion_noncoding_context_repair.pdf"],
      source_tables=[R + "Fig2_companion_noncoding_context_repair_source.tsv", H58 + "context_descriptive.tsv",
                     H58 + "family_X.tsv", H58 + "cause_split_noncoding_rows.tsv", W + "out/rederive.json"],
      commands="; ".join([JOB, REN, CHK]),
      note=("RESULTS 1837-1840 and the companion use the same estimator (summed over sets), confirmed in "
            "rederive.json context.adopted.summed_over_sets (52.113/60.431)."),
      text_edits=E["R02-02"]),
 dict(fid="R01-10", status="fixed", reproduced="yes",
      changed_in_repair=True,
      what_changes=(
        "Text and companion now use the same named estimator (summed over all credible sets of a class): ABC "
        "0.43%/3.70%, unresolved 0%/0%. The cluster-weighted repaired differences stay in RESULTS 2539 labeled "
        "'cluster-weighted'. Edits are listed under R02-02."),
      old_value="fig2.md:107 described the companion estimator as 'weighting each credible set equally'",
      new_value="'Summed over all credible sets of each class, as in the companion panel': ABC 0.43%/3.70%",
      panels=[R + "Fig2_companion_noncoding_context_repair.pdf"],
      source_tables=[R + "Fig2_companion_noncoding_context_repair_source.tsv", W + "out/rederive.json"],
      commands=JOB,
      note="Covered by the R02-02 fig2.md:107 and legend 92-100 edits.",
      text_edits=[]),
 dict(fid="R02-03", status="fixed", reproduced="yes",
      changed_in_repair=True,
      what_changes=(
        "The Fig 2A support panel is re-rendered with all 20 cells (adding ALT SAS, AST SAS, NAFLD AFR, NAFLD EAS) "
        "and printed Genes / In window / Windows columns plus a total row 436 / 100 / 1,665. Legend 34-36 and "
        "fig2.md:101 state the coverage with both measures named: UKBB ALT/AST/GGT fine-mapped at 6, 9 and 15 loci "
        "that match 8 of 767 lead windows by locus identifier (23 by position within 500 kb); 'GWAS Catalog' and "
        "'catalog lead variants' are removed; 94 genes are supported only by unfine-mapped MVP EUR strata; the UKBB "
        "arm has no credible set at PNPLA3 or TM6SF2; AMR coverage 10/37 vs 84/271."),
      old_value="Panel omitted 4 zero cells (windows summed to 1,646) and printed no totals; text said 'catalog lead variants'",
      new_value="100 of 436 multi-signal genes in a same-study fine-mapped window (+/-250 kb); 94 MVP-EUR-only; 1,665 windows",
      panels=[R + "Fig2A_finemapped_support_by_cell_repair.pdf"],
      source_tables=[R + "Fig2A_finemapped_support_by_cell_repair_source.tsv", H58 + "FIG2A_cells.tsv",
                     H58 + "FIG2A_support.tsv"],
      commands="; ".join([REN, "pandas check of credible_sets.csv: UKBB loci 6/9/15, none near PNPLA3 or TM6SF2", CHK]),
      note="fig2.md:17 is not edited; the legend and fig2.md:101 carry the qualification.",
      text_edits=E["R02-03"]),
 dict(fid="R14-03", status="fixed", reproduced="partly",
      changed_in_repair=True,
      what_changes=(
        "The causal claim 'the enzyme majority reflects more enzyme signals tested' and the 17.5%/21.1% per-eligible-"
        "pair rates are removed from the abstract, intro, fig2.md:105, Discussion and PAPER. The text keeps the "
        "counts (423/13, 22 genes with direct-trait colocalization) and adds the unselected denominators: 12,121 "
        "enzyme vs 771 direct SuSiE-tested gene-by-GWAS pairs, 739 (6.1%) and 28 (3.6%) positive, descriptive only. "
        "fig2.md:105 and PAPER state that the 171/3,844 rerun pairs were selected on first-pass colocalization and "
        "that 124 of the 171 direct pairs come from MVP EUR NAFLD, which has no fine-mapped credible sets."),
      old_value="Per eligible pair 17.5% direct vs 21.1% enzyme; 'because far more enzyme signals were tested'; '22-fold'",
      new_value=("Tested gene-by-GWAS pairs 771 direct (28 positive, 3.6%) vs 12,121 enzyme (739, 6.1%); nominal "
                 "Fisher OR 0.58, p 0.0036, pairs not independent, not reported as a test; 15.7-fold more tested "
                 "pairs (20.8-fold more credible-set pairs); 22 genes with direct-trait PP.H4 > 0.5"),
      panels=[],
      source_tables=["/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/"
                     "susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/"
                     "susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/"
                     "susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv",
                     W + "out/rederive.json"],
      commands=JOB + "; " + CHK,
      note=("Counts reproduced (423/13, 22 genes, 171/30 and 3,844/810); the per-pair-rate inference was not and is "
            "withdrawn. The discussion and intro parts are inside the R02-01 edits. The 'small overlap' wording in "
            "fig2.md:105 and PAPER 572-574 becomes 'the overlap' to agree with H57 (no directional wording)."),
      text_edits=E["R14-03"]),
 dict(fid="R13-07", status="fixed", reproduced="yes",
      changed_in_repair=True,
      what_changes=(
        "'only 13' becomes '13' and '22 genes have any direct-trait colocalization' is added in the abstract, intro, "
        "fig2.md:105 and Discussion (R02-01 and R14-03 edits). 'Coding plus splice mass' is removed by the R02-01 "
        "legend 60-66 edit. The first-round placeholder edit E334 (legend, new_text not text) is dropped."),
      old_value="'for only 13 from a direct disease or liver-fat trait'; placeholder legend edit E334",
      new_value="'for 13 from a direct disease or liver-fat trait; 22 genes have any direct-trait colocalization'",
      panels=[],
      source_tables=[W + "out/rederive.json"],
      commands=JOB,
      note="No own edits; delivered by R02-01 and R14-03 edits.",
      text_edits=[]),
 dict(fid="R13-08", status="fixed", reproduced="yes",
      changed_in_repair=True,
      what_changes=(
        "The overlap statement now uses H57 (final): 26 shared against 21.2 expected for genes of the same expression "
        "level among the 3,769 genes evaluated by multi-signal colocalization (MH OR 1.29, 95% 0.84 to 2.01; exact "
        "p 0.28, one prespecified test). S6 fired, so the text uses no directional wording and says expression was "
        "not fully controlled. The withdrawn '26 vs 33.46 expected' (which counted 10,928 unevaluated genes as "
        "negatives) is removed from fig2.md:86 and :116. Abstract opening and section heading edits are unchanged "
        "from round 1."),
      old_value="'overlap at about the rate expected for sets of these sizes (26 shared against 33.46 expected)'",
      new_value="26 shared against 21.2 expected (10 AveExpr bins), O/E 1.23, MH OR 1.29 [0.84, 2.01], exact p 0.28",
      panels=[],
      source_tables=[H57 + "RESULT_H57.md", H57 + "results.json"],
      commands="read RESULT_H57.md (H57 final; not re-run); " + CHK,
      note=("PAPER 575-580 still prints '26 shared against 33.4 expected ... odds ratio 0.76'; that span belongs to F1 "
            "and is not edited here."),
      text_edits=E["R13-08"]),
]

check_failures = [
 dict(fid="R02-01", value="fig2.md 103 new text: '16 of 19 direct clusters contain a European set'",
      who_was_right="checker",
      final_value=("15 of 19 pooled direct clusters contain a European direct set; 16 contain any European set (direct "
                   "or enzyme) and 16 is also the EUR-only clustering count. 106 of 207 direct credible sets are "
                   "European; 6,126 of 6,208 enzyme sets are non-European."),
      evidence=("sbatch 7530 (slurm 22256138) " + W + "scripts/01_rederive.py: out/rederive.json "
                "direct_clusters_with_EUR_direct_set 15, direct_clusters_with_any_EUR_set 16, "
                "EUR_direct_clusters_within_EUR_only 16, direct_sets_by_ancestry EUR 106 AMR 95 AFR 5 EAS 1.")),
 dict(fid="R14-03", value=("Claim: per-signal colocalization rate is similar (17.5% vs 21.1% of eligible pairs), so the "
                           "enzyme majority reflects 22-fold more enzyme signals tested"),
      who_was_right="checker",
      final_value=("The 171/3,844 eligible pairs come only from the 816 rerun rows (29 direct, 787 enzyme), selected "
                   "on first-pass colocalization; every not-rerun tested row has PP.H4 <= 0.498. Per SuSiE-tested "
                   "gene-by-GWAS row: direct 28/771 = 3.6%, enzyme 739/12,121 = 6.1% (nominal Fisher OR 0.58, "
                   "p 0.0036, rows not independent). Per tested credible-set pair 30/2,035 = 1.5% vs 810/42,382 = 1.9%. "
                   "Fold differences 15.72 (rows) and 20.83 (pairs); 22 = 3,844/171. Causal clause withdrawn."),
      evidence=("sbatch 7530 (slurm 22256138) 01_rederive.py on susie_coloc_all_gwas.csv, susie_coloc_signal_pairs.tsv.gz "
                "and gwas_trait_tier.tsv (35 main tier-1/2 studies): rederive.json tested_rows, rerun_rows_by_tier, "
                "max_PPH4_not_rerun 0.498, fold_tested_rows 15.72, fold_tested_cs_pairs 20.83, "
                "eligible_direct_pairs_by_study MVP_NAFLD_EUR 124.")),
]

conflicts = [
 dict(index=5, resolution=(
   "Dropped the parenthesis with pre-repair values. F2d writes title.md 22-24 against the current file (supersedes "
   "Claims E197, same old_text): keeps E197's structure (definition, 19 direct-trait regions, no causal variant), adds "
   "'pooled across ancestries ... confounded with GWAS ancestry and fine-mapping coverage ... no trait-class call in "
   "any ancestry stratum', and replaces 'Most noncoding probability (74.9% direct, 93.7% enzyme) has no deposited "
   "regulatory context' with 'After position re-annotation, most noncoding probability in both classes lies outside "
   "the deposited promoter, open-chromatin and ABC annotations'. Edit listed under R02-01.")),
 dict(index=6, resolution=(
   "Causal clauses and 17.5%/21.1% removed everywhere F2d writes (abstract 32-35, intro 77-84, fig2.md 105, "
   "Discussion 105-121, PAPER 567-574). Counts kept (423 of 436; 13; 22 genes). Added the unselected denominator: "
   "12,121 enzyme vs 771 direct SuSiE-tested gene-by-GWAS pairs. Per-pair rates 6.1% and 3.6% appear only in "
   "fig2.md 105 and PAPER, labeled descriptive with pairs not independent; fig2.md 105 and PAPER say the 171/3,844 "
   "rerun pairs were selected on first-pass colocalization.")),
 dict(index=8, resolution=(
   "Installed F2d's abstract 24-32 and intro 77-84 text; Claims E195/E196 should be dropped. No ratio printed in the "
   "abstract or intro. The credible-set-pooled clause (15.81% vs 3.26%) is not added, so 'pooled' has one meaning "
   "(across ancestries).")),
 dict(index=16, resolution=(
   "fig2.md 103 prints '106 of 207 are European and 15 of the 19 direct clusters contain a European direct set'. The "
   "abstract replaces 'most direct credible sets are European' with '106 of 207 direct-trait but 82 of 6,208 enzyme "
   "credible sets are European'; legend 60-66, PAPER 545 and STATUS 44 use the same counts.")),
 dict(index=18, resolution=(
   "Added three edits under R02-01: PAPER 558-560 ('None of the 57 comparisons stratifies by ancestry or trait or has "
   "matched fine-mapping coverage, so they are correlated sensitivities of the pooled value, not replications of a "
   "trait-class difference', plus the H58 withdrawal of the unresolved-context fractions); docs/README.md 23-25 and "
   "figures/main/fig2_genetics/README.md 5-12 (pooled qualifier with EUR 29.88% vs 6.92%, AMR 6.73% vs 4.54%, no "
   "trait-class call). The Fig 2 README also says the installed companion's 52.11%/60.43% unresolved fractions are "
   "withdrawn (repaired 0%) and names the estimator as summed over all credible sets.")),
]

res = dict(ws="F2d", entries=entries, check_failures_resolved=check_failures, conflicts_resolved=conflicts)
(OUT / "structured_output.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
print(json.dumps(res, ensure_ascii=False))
