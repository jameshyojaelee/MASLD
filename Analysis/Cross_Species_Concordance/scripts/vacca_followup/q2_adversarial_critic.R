#!/usr/bin/env Rscript
# q2_adversarial_critic.R
# Consolidate the attack-relevant Q2 numbers and quantify the two key confounds:
#   (1) cross-species concordance vs HUMAN-ONLY baseline (does mouse add anything?)
#   (2) circularity: the "continuous concordance" rankings put human LFC on BOTH
#       sides (our_stat = translatability*sign(human_lfc) ; Vacca = human progLFC),
#       so high RRHO/Spearman is partly human-vs-human.
suppressPackageStartupMessages({ library(data.table) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
D    <- file.path(BASE,"Analysis/Cross_Species_Concordance/results/vacca_benchmark")

sm  <- fread(file.path(D,"q2_threshold_free_summary.csv"))
get <- function(m) sm[metric==m, value]

cat("===== CONFOUND 1: cross-species vs human-only =====\n")
cat(sprintf("AUROC  S_geom(cross-sp)=%.3f  vs  human-only=%.3f  -> human-only is %s\n",
  get("AUROC_S_geom"), get("AUROC_h_signlp(human-only)"),
  ifelse(get("AUROC_h_signlp(human-only)")>get("AUROC_S_geom"),"HIGHER (mouse hurts AUROC)","lower")))
cat(sprintf("fgsea  S_geom NES=%.2f  vs  human-only NES=%.2f  -> human-only is %s\n",
  get("fgsea_NES_S_geom"), get("control_humanOnly_fgsea_NES"),
  ifelse(get("control_humanOnly_fgsea_NES")>get("fgsea_NES_S_geom"),"HIGHER","lower")))
cc <- fread(file.path(D,"q2_continuous_concordance_vs_vacca_prog.csv"))
cat(sprintf("Spearman vs Vacca prog: hm_agreement=%.3f  human_ss_only=%.3f  mouse_ss_only=%.3f\n",
  cc[metric=="hm_agreement"&vacca_target=="vacca_prog_lfc_signed",spearman],
  cc[metric=="human_ss_only",spearman], cc[metric=="mouse_ss_only",spearman]))

cat("\n===== CONFOUND 2: pure unsigned cross-species concordance (no human-LFC leakage) =====\n")
fg <- fread(file.path(D,"q2_threshold_free_fgsea.csv"))
cat("These rankings put NO human direction on the axis -> the honest threshold-free test:\n")
print(fg[ranking %in% c("n_concordant","translatability_score","signed_nconc","mean_h_lfc"),
         .(ranking, NES=round(NES,2), padj=signif(padj,2))])
cat(sprintf("  -> n_concordant (purest cross-species count) NES=%.2f ; mean_h_lfc(HUMAN ALONE) NES=%.2f\n",
  fg[ranking=="n_concordant",NES], fg[ranking=="mean_h_lfc",NES]))

cat("\n===== CONFOUND 3: cutoff grid -- is 2.8x a floor or a ceiling? =====\n")
gr <- fread(file.path(D,"q2_cutoff_bottleneck_grid.csv"))
canon <- gr[human_abs_logFC==0 & human_padj==0.05 & mouse_padj==0.05 & n_concordant_req==3]
cat(sprintf("canonical: core=%d overlap=%d pct_in_vacca=%.1f%% OR=%.2f\n",
  canon$core_size, canon$overlap_vacca951, canon$pct_core_in_vacca, canon$fisher_OR))
cat(sprintf("OR range across %d-cell grid: %.2f -> %s ; but pct_core_in_vacca NEVER exceeds %.0f%% until core<%d\n",
  nrow(gr), min(gr$fisher_OR), ifelse(any(is.infinite(gr$fisher_OR)),"Inf",sprintf("%.1f",max(gr$fisher_OR,na.rm=TRUE))),
  max(gr[fisher_OR<Inf & core_size>500, pct_core_in_vacca]),
  min(gr[pct_core_in_vacca>50, core_size])))
cat("Tightening cutoffs RAISES OR (2.7->Inf) but only by SHRINKING the core to a handful;\n")
cat("at any core size >500 genes the % that are Vacca genes stays ~11-13%.\n")

cat("\n===== DIRECTION-CONCORDANCE: the cleanest threshold-free human-side stat =====\n")
th <- fread(file.path(D,"q2_direction_concordance_threshold.csv"))
cat(sprintf("ALL %d shared genes, NO cutoff: human dir-concordance with Vacca = %.1f%% (p~0)\n",
  th[min_abs_both_threshold==0,n], th[min_abs_both_threshold==0,concordance_pct]))
cat("  -> but this is HUMAN-vs-HUMAN (our human LFC vs Vacca human LFC); NOT cross-species.\n")
