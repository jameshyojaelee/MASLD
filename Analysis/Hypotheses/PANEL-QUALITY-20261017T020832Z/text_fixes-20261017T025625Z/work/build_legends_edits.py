"""Build legends_edits.json for docs/manuscript/05_figure_legends.md (staged only, never applied).
Checks that each `old` occurs exactly once in the live file and attaches the re-read values
from legends_reread.tsv (produced by reread_legends.py)."""
import json
import os

import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
WORK = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-QUALITY-20261017T020832Z/text_fixes-20261017T025625Z/work")
FILE = "docs/manuscript/05_figure_legends.md"
live = open(os.path.join(ROOT, FILE), encoding="utf-8").read()
rr = pd.read_csv(os.path.join(WORK, "legends_reread.tsv"), sep="\t", dtype=str)

E = []


def edit(eid, old, new, reason, source, numbers):
    E.append(dict(id=eid, file=FILE, old=old, new=new, reason=reason, source=source, numbers_in_new=numbers))


edit("L-01",
     "transcription in the native 774-unit interface, restricted to the 3,769 genes",
     "transcription in the native 774-unit interface (774 genetic-identity clusters, 138 control and 636 disease; not verified donors), restricted to the 3,769 genes",
     "Fig 1B: name the unit. The 774 units are genetic-identity clusters, not verified persons.",
     "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/arm_deg_counts.tsv (row primary_source_controls: n_samples, n_control, n_disease); primary_unit_roster.tsv (unit_evidence)",
     "774; 138; 636")
edit("L-02",
     "Brackets bind the\nendpoint bands to their separate multi-signal COLOC (436 genes) and\nsingle-signal COLOC-only (560 genes) denominators.",
     "Brackets group the\nendpoint bands under their separate multi-signal COLOC (436 genes) and\nsingle-signal COLOC-only (560 genes) denominators. Each gene is counted once, in the band of its largest PP.H4 (SuSiE for multi-signal genes, ABF for single-signal-only genes): above 0.9 (244 multi-signal and 67 single-signal-only genes), 0.8 to 0.9 inclusive (89 and 89) and above 0.5 but below 0.8 (103 and 404).",
     "Fig 2A: replace 'bind' with a literal phrase and define the PP.H4 bands with their counts.",
     "figures/main/fig2_genetics/panels/fig2A_gwas_cascade_source.tsv (record_type posterior_band, metric coloc_genes, value); band rule in scripts/figures/fig1_gwas_creative_options.py (su > 0.9; su >= 0.8; else; per-gene maximum)",
     "436; 560; 0.9; 244; 67; 0.8; 0.9; 89; 89; 0.5; 0.8; 103; 404")
edit("L-03",
     "Absence from this panel is not a negative colocalization call.",
     "The dashed line marks PP.H4 = 0.5. Labels are limited to a fixed five-gene roster (HSD17B13, TM6SF2, MBOAT7, GCKR and RORA); TM6SF2, MBOAT7 and GCKR are not among the 327 plotted genes, so only HSD17B13 and RORA are labeled. Absence from this panel is not a negative colocalization call.",
     "Fig 2B: describe the dashed line and the label rule.",
     "figures/candidates/review-r3-fixes-20261015T112148Z/F2abce/repair/Fig2B_PIP_vs_SuSiE-coloc_sameRow_susie_repair_source.csv (327 rows; is_label); render script fig2b_render_repair.R (roster, geom_hline yintercept 0.5 on y = pp4p). Installed PDF is byte-identical to the repair PDF.",
     "0.5; 327")
edit("L-04",
     "are RORA and FABP1. Both use EUR GWAS and do not provide evidence for\npopulation-dependent biology.",
     "are RORA and FABP1. Both use EUR GWAS and do not provide evidence for\npopulation-dependent biology. The FABP1 PIP track is the stored COLOC-targeted SuSiE fine-mapping of UKBB ALT (1,135 variants). Neither stored fine-mapping set has PIPs for the MVP European NAFLD window at RORA, so its PIP track is a single-effect SuSiE fit (L = 1, 95% coverage) run for display on the European GWAS z-scores with the European LD reference (2,282 variants); it is not the COLOC SuSiE fit.",
     "Fig 2G-H: state the PIP track source. Source contradicts the number check: the RORA track is an inline display fit, not the COLOC SuSiE fit.",
     "figures/main/fig2_genetics/panels/Fig2G_RORA_NAFLD_source.csv and Fig2H_FABP1_ALT_source.csv (pip_source, n_pip_variants, gwas); scripts/figures/figS09_locus_zoom_disease.R (stored-PIP fallback order; susie_rss L = 1, coverage = 0.95)",
     "1,135; 1; 95; 2,282")
edit("L-05",
     "give a direction odds ratio of 0.132 (95% CI\n0.090–0.193).",
     "give a conditional direction odds ratio of 0.132 (30 discordant pairs with the lncRNA member up against 227 with the protein-coding member up; Wald 95% CI on the log odds ratio\n0.090–0.193).",
     "Fig 3G: state the lncRNA OR interval type and its discordant-pair counts.",
     "figures/main/fig3_bulk_transcriptomics/source_tables/fig3_lncrna_role/summary.tsv (matched_direction: conditional_odds_ratio, discordant_lncrna_up, discordant_coding_up); scripts/figures/fig3_lncrna_role.R (matched_ci Wald on log scale)",
     "0.132; 30; 227; 95; 0.090; 0.193")
edit("L-06",
     "three-cohort random-effects meta-analysis (REML, Knapp-Hartung); the y axis",
     "three-cohort random-effects meta-analysis (REML, Knapp-Hartung) of within-cohort OLS fits of the program score on sex, fibrosis stage and NAS with HC3 standard errors; the y axis",
     "Fig 4E: state the covariates and the HC3 within-cohort standard errors.",
     "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_three_cohorts/cohort_effects.tsv (se_hc3; n_model_columns 4); producers/program_axis_map_774.R (sex_final, fibrosis_stage, nas_score). Installed 4E is the review-r3 three-cohort render (n_cohorts 3).",
     "none")
edit("L-07",
     "scores across nine prespecified fixed-projection windows in all five Resource\ncohorts.",
     "scores across nine prespecified fixed-projection windows (each spanning 0.2 of the within-cohort fixed-projection percentile, centered at 0.1 to 0.9) in all five Resource\ncohorts (774 genetic-identity clusters: 30 GSE126848, 71 GSE130970, 207 GSE135251, 140 GSE162694 and 326 GSE213621; per window 6, 13 to 15, 41 to 42, 28 and 65 to 66, respectively).",
     "Fig 4F: give unit counts per cohort and per window.",
     "figures/candidates/fig3-fig4-layout-774-20261013T075956Z/panels/fig4/fig4f_continuum_program_trajectories_source.tsv (window bounds; n_participants by dataset); RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-corrected774-20261011T152106Z/projection/participant_scores.tsv (unique sample_id by dataset; fixed_projection axis_percentile). Installed 4F is byte-identical to the layout-774 candidate.",
     "0.2; 0.1; 0.9; 774; 30; 71; 207; 140; 326; 6; 13; 15; 41; 42; 28; 65; 66")
edit("L-08",
     "The 774-unit substrate tests 24,115 genes; memberships, weights and coordinates\nare unchanged.",
     "Pathway and program scores use the 24,115-gene expression matrix of the 774-unit fit; memberships, weights and coordinates\nare unchanged from the geometry built on the 23,370-gene, 844-record expression universe.",
     "Fig 4G: separate the score universe from the geometry universe. Source contradicts the number check: 24,115 is correct for the scores; 23,370 applies only to the geometry.",
     "RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z/inputs/frozen_contract.json (gene_family_size), read by fig4h_disease_only.py; 774 merged_dge.rds 24115 x 774; molecular_systems_atlas/provenance/geometry_input_manifest.tsv expression_universe_only = F_five merged_dge.rds 23370 x 844",
     "24,115; 774; 23,370; 844")
edit("L-09",
     "the adjacent bar includes protein-significant genes only and shows",
     "the adjacent bar includes protein-significant genes only (limma BH-adjusted P <0.05, BH over all 6,071 quantified proteins) and shows",
     "Fig 5B: state the protein significance rule and its BH family.",
     "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv (6,071 rows; padj = BH of pvalue); scripts/figures/fig5b_protein_triage.R (protein_padj < 0.05)",
     "0.05; 6,071")
edit("L-10",
     "fibrosis 0-1 and NAS 0-1, so it combines disease status with severity.",
     "fibrosis 0-1 and NAS 0-1, so it combines disease status with severity. Histology q values are BH-adjusted over 125 tests (25 proteins by 5 histology features). Process NES and q come from preranked fgsea on the adjusted protein moderated t statistic (gene sets of 8 to 600 proteins), with BH over all 1,145 tested sets (1,091 Reactome, 50 Hallmark and 4 curated).",
     "Fig 5C: name the correction families and the enrichment test.",
     "figures/candidates/bulk-context-774-20261014T140904Z_22236311/proteomics/panel4c_histology_partial.tsv (125 rows, correction_family); protein_pathway_enrichment.tsv (1,145 sets, family, padj = BH of pval); Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R (fgsea on t, minSize 8, maxSize 600)",
     "125; 25; 5; 8; 600; 1,145; 1,091; 50; 4")
edit("L-11",
     "(gene-bootstrap 95% intervals above 1)",
     "(percentile 95% intervals from 2,000 gene resamples, seed 20260923, all above 1)",
     "Fig 5D: give bootstrap draws, interval type and seed.",
     "figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/fig5d_locus_background_summary.tsv (bootstrap; liver_enzyme ratio_lo); F5bcd/scripts/fig5d_summaries.py (default_rng(20260923), 2000 resamples, percentile 2.5/97.5)",
     "95; 2,000; 20260923; 1")
edit("L-12",
     "On the endogenous variant-to-gene |z| task, adding either pretrained\nembedding to distance lowers the score (HyenaDNA −0.022, Nucleotide Transformer\n−0.029; intervals exclude zero).",
     "On the endogenous variant-to-gene |z| task (held-out LD fold 4: 111,426 variant-gene rows in 55 LD blocks), adding either pretrained\nembedding to distance lowers the held-out Spearman against the distance ridge (HyenaDNA −0.022, Nucleotide Transformer\n−0.029; paired 95% percentile intervals from 2,000 LD-block bootstrap resamples, seed 20260831, exclude zero).",
     "Fig 7B: state held-out fold, n, metric and paired bootstrap. Seed is 20260831 in the producing file.",
     "Analysis/MASLD_Model_Benchmark/executions/dnalm_heads-20260831T183647Z/logs/03_intervals.json (seed, bootstrap_resamples, bootstrap_unit, primary_held_out_fold, heldout rows/blocks, paired_deltas abs_z m1_dist_emb MINUS b1_dist); fig7_panels_rebuild-20260901T012245Z/tables/panel_7b_gains.tsv (task_metric)",
     "4; 111,426; 55; −0.022; −0.029; 95; 2,000; 20260831")
edit("L-13",
     "not interpreted as clean-only skill.\nGeneformer V2 has zero qualified clean studies among the seven displayed\ncohorts and therefore no clean-only delta. **(D)**",
     "not interpreted as clean-only skill. **(D)**",
     "Fig 7C: Geneformer V2 is not drawn in 7C; its exposure statement moves to 7F (L-17).",
     "Analysis/MASLD_Model_Benchmark/executions/fig7-source-exposure-correction-20261014/tables/panel_7c_exposure.tsv (no geneformer row); pdftotext of installed 7C has no 'Geneformer'",
     "none")
edit("L-14",
     "instrument's Spearman, with its 95% interval, beside each cohort's own",
     "instrument's Spearman, with its 95% percentile bootstrap interval (units resampled within stage, 10,000 draws, seed 20260923), beside each cohort's own",
     "Fig 7D: name the interval method, draws and seed.",
     "figures/candidates/review-r3-fixes-20261015T112148Z/F7b/repair/fig7_transfer_primary_units.tsv (ci_B, ci_seed, ci_method); installed 7D is byte-identical to the F7b repair PDF",
     "95; 10,000; 20260923")
edit("L-15",
     "viral hepatitis). A nonlinear model with no biological",
     "viral hepatitis). The signed-eQTL block uses 111,426 held-out variant-gene rows in 55 LD blocks. The LD-block join uses 644,556 training-fold variant-gene rows in 208 LD blocks for the block-left-in correlation and the 15,192 of them whose gene also has rows in another block for the block-held-out correlation. The routing block uses 521 bulk RNA samples from five cohorts; no donor key exists, so samples are the unit. A nonlinear model with no biological",
     "Fig 7D: give the n of each block.",
     "dnalm_heads-20260831T183647Z/logs/03_intervals.json (heldout rows/blocks); training_matrix.parquet + splits/variant_splits.tsv recomputed as in scripts/03_intervals.py (r 0.86687 / 0.02101 reproduce the recorded values); routed_multihead_eval-20260829T230000Z/results.json (n_samples, unit_of_inference, donor_key_exists, cohorts)",
     "111,426; 55; 644,556; 208; 15,192; 521")
edit("L-16",
     "each held-out study. Per-study deltas are not drawn in this panel; the gains\ndrawn are all-study values.",
     "each held-out study. Per-study deltas (each model minus the liver-specific scVI encoder carrying the same head, within one held-out study) are not drawn in this panel; the gains\ndrawn are all-study gains over HVG + 50-dimensional PCA of raw counts carrying the identical head (102 donors, 50,000 cells, seven studies held out by study; 10,000-resample donor cluster bootstrap).",
     "Fig 7F: name the reference for the per-study deltas and for the drawn gains, and give n.",
     "figures/candidates/review-r3-fixes-20261015T112148Z/F7a/source_tables/panel_7f_gain.tsv (block b block_metric, unit_and_n); encoder-field-score-21314196/evaluation/paired_deltas_by_held_out_study.tsv (reference_id scvi_liver_latent); panel_7f_exposure_grid.tsv (donors, rows per study); installed 7F x axis reads 'gain over HVG + PCA'",
     "50; 102; 50,000; 10,000")
edit("L-17",
     "reported. Geneformer V2, with no clean held-out study, is drawn hollow gray in\nthe gain forest:",
     "reported. Geneformer V2 has no qualified clean study among the seven held-out studies (three encoder-seen, four with unresolved exposure) and therefore no clean-only delta; it is drawn hollow gray in\nthe gain forest:",
     "Fig 7F: receives the Geneformer exposure statement moved from 7C (L-13).",
     "figures/candidates/review-r3-fixes-20261015T112148Z/F7a/source_tables/panel_7f_exposure_grid.tsv (Geneformer study_exposure_state: encoder_seen 3, unknown 4); panel_7f_gain.tsv block c Geneformer gain NaN",
     "3; 4")
edit("L-18",
     "and on those records it matches rather than beats a ridge fitted inside each\ncohort.",
     "and on those records it shows no demonstrated advantage over a ridge fitted inside each\ncohort.",
     "Fig 7 text: intervals that include zero do not establish equivalence, so 'matches' overstates the result.",
     "Analysis/MASLD_Model_Benchmark/executions/external-validation-expansion-20260901T190000Z/results.json (vs_native_paired_delta and ci95 for GSE213621 and GSE276114; GSE268273 no interval); docs/RESULTS.md wording 'No demonstrated advantage over the recorded transductive in-cohort ridge'",
     "none")

bad = []
for e in E:
    e["old_count"] = live.count(e["old"])
    if e["old_count"] != 1:
        bad.append((e["id"], e["old_count"]))
    sub = rr[rr["edit_id"] == e["id"]]
    e["value_reread"] = "; ".join(f"{c}={v}" for c, v in zip(sub["column_or_row"], sub["value_reread"]))
    for banned in ("—", "frozen", "progression", "cascade"):
        assert banned not in e["new"], (e["id"], banned)
json.dump(E, open(os.path.join(WORK, "legends_edits.json"), "w"), indent=1, ensure_ascii=False)
print(f"{len(E)} edits; old_count != 1: {bad}")
