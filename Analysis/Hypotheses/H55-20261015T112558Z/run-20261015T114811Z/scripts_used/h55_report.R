# H55 report: writes RESULT_H55.md from outputs/h55_results.rds and draws supplementary panels.
# Usage: Rscript h55_report.R <run_dir> <result_md_path>
suppressMessages({library(data.table); library(ggplot2)})
a <- commandArgs(trailingOnly = TRUE); RUN <- a[1]; RES_MD <- a[2]
OUT <- file.path(RUN, "outputs"); DG <- file.path(RUN, "diagnostics"); FIG <- file.path(RUN, "figures")
r <- readRDS(file.path(OUT, "h55_results.rds"))
f3 <- function(x) formatC(x, digits = 3, format = "f"); g3 <- function(x) formatC(x, digits = 3, format = "g")
L <- c()
add <- function(...) L <<- c(L, paste0(...))
add("# RESULT H55: GTEx postmortem liver fat versus Resource disease-state genes")
add("")
add("Mode: ", r$mode, ". Permutation draws ", r$B, ", bootstrap resamples ", r$BBOOT, ", calibration draws ", r$BCAL, ". Seed 20260923.")
add("Prespecification: `prespec/PRESPEC_H55.md` (sha256 293114cd56cb76a6906480ce4b81d63c87fa41869a6dd013cb3d47e397510a0a), verified with all listed input hashes before the run.")
add("")
add("## Verdict")
add("")
add("**", r$verdict, "**")
add("")
add("## Gates and diagnostics (reported before the primary)")
add("")
if (!is.null(r$S1)) add("- S1 morphometry: ", r$S1$slides, " slides, ", r$S1$ok, " ok, ", r$S1$error, " errors; model donors with a fat fraction ", r$S1$model_with_fat, " (stop if errors > 10 or < 210). Status by center: `diagnostics/d1_morph_status_by_center.tsv`.")
if (!is.null(r$S2)) add("- S2 determinism: ", r$S2$n, " pilot slides, maximum relative difference ", g3(r$S2$max_rel), " (stop if > 1e-9).")
if (!is.null(r$S3)) add("- S3 / PC1 exposure validity: Spearman of fat fraction with note ordinal ", f3(r$S3$rho), " (n ", r$S3$n, "; stop if < 0.5); with stated % ", f3(r$S3$pct[2]), " (n ", r$S3$pct[1], "). Within center: ",
                        paste(sprintf("%s %.3f (n %d)", r$S3$by_center$SMCENTER, r$S3$by_center$rho, r$S3$by_center$n), collapse = ", "), ".")
if (!is.null(r$S3b)) add("- S3b: median fat fraction of ", r$S3b$n, " no-mention model donors ", formatC(r$S3b$median, digits = 5, format = "f"),
                         " (stop if > 0.0335; prediction <= 0.0056: ", ifelse(r$S3b$prediction_met, "met", "not met"), "). Quartiles by center: `diagnostics/d4_no_mention_unspecified_fat.tsv`.")
if (!is.null(r$S4)) add("- S4 gene family: canonical calls ", r$S4$canonical, ", in GTEx v26 ", r$S4$in_gtex, ", family recomputed ", r$S4$family, ", equal to frozen list: ", r$S4$equal, ".")
if (!is.null(r$d6)) add("- Diagnostic 6: 1 - R^2 of fat rank on covariates ", f3(r$d6), ifelse(r$d6 < 0.5, " (below 0.5: reduced power).", "."))
if (!is.null(r$d7)) add("- Diagnostic 7: quality-weight range ", f3(r$d7$range[1]), " to ", f3(r$d7$range[2]), "; donors with weight < 0.3: ", r$d7$n_below_0.3, ". Trend plot `diagnostics/d7_voom_trend.pdf`.")
if (!is.null(r$S5)) add("- S5 identity check: max |Freedman-Lane slope - lmFit coefficient| ", g3(r$S5), " (stop if > 1e-6). Moderated-t identity ", g3(r$modt_identity), ".")
if (!is.null(r$primary)) add("- Diagnostic 9: Freedman-Lane null mean ", f3(r$primary$null["mean"]), ", SD ", f3(r$primary$null["sd"]), ", 2.5-97.5% [", f3(r$primary$null["q025"]), ", ", f3(r$primary$null["q975"]), "]",
                             ifelse(r$d9_flag, " (flag: SD outside 0.15-0.30, planning value 0.215).", " (planning value 0.215)."))
if (!is.null(r$famC)) add("- Diagnostic 10: Family C calibration, fraction of ", min(r$BCAL, r$B), " null draws with any q < 0.05: ", f3(r$famC$cal_frac),
                          ifelse(r$famC$use_empirical, " (> 0.07: Catalog column uses per-gene Freedman-Lane empirical p).", " (<= 0.07: Catalog column uses limma moderated p)."))
if (!is.null(r$stop_rule)) { add(""); add("Stop rule ", r$stop_rule, " fired: ", r$stop_msg, ". The analysis ends here; no Catalog column is written and no text changes.") }
if (!is.null(r$primary)) {
  add(""); add("## Controls"); add("")
  add("- PC2 PLIN2: slope ", f3(r$PC2$slope), " log2 CPM from lowest- to highest-fat donor, limma moderated p ", g3(r$PC2$p), "; prediction (slope > 0, p < 0.05) ", ifelse(r$PC2$pass, "met", "not met"), ".")
  add("- PC3 composite (197 genes): rho ", f3(r$PC3$rho), ", Freedman-Lane p ", g3(r$PC3$p["p_two"]), "; S6 ", ifelse(r$PC3$pass, "does not fire", "FIRES"), "; predicted range 0.3-0.8 ", ifelse(r$PC3$in_predicted_range, "met", "not met"), ".")
  add("- NC1 tissue area: rho ", f3(r$NC1$rho), ", Freedman-Lane p ", g3(r$NC1$p["p_two"]), "; S7 ", ifelse(r$NC1$fail, "FIRES", "does not fire"), "; prediction (|rho| < 0.2, p > 0.05) ", ifelse(r$NC1$prediction_met, "met", "not met"), ".")
  add("- NC2 ischemic time (joint model): rho ", f3(r$NC2$rho_isch), " versus rho_fat ", f3(r$Q2$rho_fat), "; prediction (-0.2 to 0.3 and below rho_fat) ", ifelse(r$NC2$prediction_met, "met", "not met"), ".")
  add(""); add("## Primary (Family A, one test, two-sided, alpha 0.05)"); add("")
  add("Spearman rho across ", r$primary$n_genes, " family genes between bulk disease-vs-control logFC and GTEx fat-rank slope, ", r$primary$n_donors, " donors, residual df ", r$df_res, ": ",
      "**rho ", f3(r$primary$rho), "**, Freedman-Lane p two-sided ", g3(r$primary$p["p_two"]), " (upper ", g3(r$primary$p["p_upper"]), ", lower ", g3(r$primary$p["p_lower"]), "; ", r$B, " within-stratum donor permutations; floor 2/", r$B + 1, ").")
  add(""); add("## Family B qualifiers (Holm over 2, tested only if Family A rejects)"); add("")
  add("- Q1 Hardy 1-4 (n ", r$Q1$n, "): rho ", f3(r$Q1$rho), ", Freedman-Lane p two-sided ", g3(r$Q1$p["p_two"]), ", Holm p ", g3(r$Q1$p_holm), "; ",
      ifelse(r$Q1$tested, ifelse(r$Q1$supports, "supports the wording", "does not support the wording"), "not tested (Family A not rejected); values descriptive"), ".")
  add("- Q2 fat versus ischemic time: rho_fat ", f3(r$Q2$rho_fat), ", rho_isch ", f3(r$Q2$rho_isch), ", D ", f3(r$Q2$D), " (bootstrap 95% ", f3(r$Q2$D_boot_ci[1]), " to ", f3(r$Q2$D_boot_ci[2]),
      "), bootstrap one-sided p ", g3(r$Q2$p_boot), ", Holm p ", g3(r$Q2$p_holm), "; joint-permutation p (descriptive) ", g3(r$Q2$p_joint_perm["p_two"]), "; resamples with a dropped constant column ", r$Q2$n_resample_col_dropped, "; ",
      ifelse(r$Q2$tested, ifelse(r$Q2$supports, "supports the wording", "does not support the wording"), "not tested (Family A not rejected); values descriptive"), ".")
  add(""); add("## Family C (Gene Catalog column; BH over ", r$primary$n_genes, " family genes)"); add("")
  add("Genes at q < 0.05: ", r$famC$n_q05, if (!is.na(r$famC$agree)) paste0("; fraction with slope sign equal to bulk sign ", f3(r$famC$agree)) else "", ".")
  sc <- r$state_counts; add("State counts over 1,854 calls: ", paste(sprintf("%s %d", sc$state, sc$N), collapse = "; "), ".")
  add("Catalog column file written: ", r$catalog_column_written, ifelse(r$catalog_column_written, " (`outputs/catalog_column_gtex_postmortem_fat.tsv`; annotation only, for the next Catalog build).", "."))
  add(""); add("## Descriptive versions (no alpha)"); add("")
  add("| Version | rho | Freedman-Lane p | donors | genes |"); add("|---|---|---|---|---|")
  for (i in seq_len(nrow(r$desc))) add("| ", r$desc$version[i], " | ", f3(r$desc$rho[i]), " | ", g3(r$desc$p_two[i]), " | ", r$desc$n_donors[i], " | ", r$desc$n_genes[i], " |")
  if (isTRUE(r$tnfa_wording)) add("", "Wording rule: the family without TNFA/acute-phase genes has p >= 0.05 while the primary passes, so the alignment depends on acute-phase and immediate-early genes.")
  add(""); add("## Reading (prespecified rule)"); add("")
  rd <- if (grepl("^not interpretable", r$verdict)) "The primary was computed and recorded but is not interpretable because a control failed. No Catalog column and no text change." else
    if (grepl("^null", r$verdict)) "Alignment is not detected for the 1,002 genes. Power was 0.86 at full and 0.54 at half transfer of the biopsy steatosis effect (MDE about rho 0.60 [Estimate]); the result does not exclude alignment of that size or smaller. The calls are not relabeled as cohort-specific." else
    if (grepl("^pass", r$verdict)) paste0("Among 1,002 disease-state calls measurable in GTEx liver, genes change with image-measured histologic fat in ", r$primary$n_donors,
      " postmortem donors not selected for liver disease in the same direction (rho ", f3(r$primary$rho), ") as between disease and control biopsies. ",
      ifelse(r$Q1$supports, "Also in donors who did not die on a ventilator. ", "Q1 not supported: the alignment is concentrated in ventilator deaths and independence from critical illness is not shown. "),
      ifelse(r$Q2$supports, "More than with ischemic time. ", "Q2 not supported: the result is described as aligned with a general postmortem liver-stress axis. "),
      "The claim is histologic steatosis of unrecorded cause in postmortem liver; it does not assign any biopsy-cohort call to steatosis.") else
    "A significant negative rho: reported as opposite direction."
  add(rd)
  add(""); add("Outputs: `outputs/per_gene_1854_calls.tsv`, `outputs/family_slopes.tsv`, `outputs/descriptive_versions.tsv`, `outputs/null_rho_primary_pc3_nc1.tsv.gz`, permutation and bootstrap index matrices (`outputs/*.rds`), figures in `figures/`.")
}
writeLines(L, RES_MD)

## Supplementary panels (individual PDFs, 6 pt, no bold)
if (!is.null(r$primary)) {
  dir.create(FIG, showWarnings = FALSE)
  th <- theme_classic(base_size = 6) + theme(text = element_text(size = 6, colour = "black"), axis.text = element_text(size = 6, colour = "black"),
                                              plot.title = element_text(size = 6, face = "plain"), legend.position = "none", axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25))
  fs <- fread(file.path(OUT, "family_slopes.tsv")); fs[, dir := ifelse(bulk_logFC > 0, "up", "down")]
  p1 <- ggplot(fs, aes(bulk_logFC, slope, colour = dir)) + geom_hline(yintercept = 0, linewidth = 0.2, colour = "#9E9E9E") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "#9E9E9E") +
    geom_point(size = 0.4, alpha = 0.6, stroke = 0) + scale_colour_manual(values = c(up = "#C9265E", down = "#1565C0")) +
    labs(x = "Bulk disease vs control logFC", y = "GTEx fat-rank slope (log2 CPM)") + th
  ggsave(file.path(FIG, "H55_A_scatter_bulk_vs_gtex_slope.pdf"), p1, width = 2.4, height = 2.2, useDingbats = FALSE)
  nl <- fread(file.path(OUT, "null_rho_primary_pc3_nc1.tsv.gz"))
  p2 <- ggplot(nl, aes(rho_null)) + geom_histogram(bins = 60, fill = "#9E9E9E", colour = NA) + geom_vline(xintercept = r$primary$rho, colour = "#C9265E", linewidth = 0.4) +
    labs(x = "Spearman rho (permuted donors)", y = "Permutations") + th
  ggsave(file.path(FIG, "H55_B_null_histogram.pdf"), p2, width = 2.4, height = 2.0, useDingbats = FALSE)
  d <- r$desc; h0 <- d[version == "hardy0_only"]
  hs <- data.table(set = factor(c("All donors", "Hardy 1-4", "Hardy 0"), levels = c("Hardy 0", "Hardy 1-4", "All donors")),
                   rho = c(r$primary$rho, r$Q1$rho, h0$rho), lo = c(r$primary$p["null_q025"], r$Q1$p["null_q025"], h0$null_q025),
                   hi = c(r$primary$p["null_q975"], r$Q1$p["null_q975"], h0$null_q975), n = c(r$primary$n_donors, r$Q1$n, h0$n_donors))
  hs[, lab := paste0(set, " (n ", n, ")")]; hs[, lab := factor(lab, levels = lab[order(as.integer(set))])]
  p3 <- ggplot(hs, aes(y = lab)) + geom_vline(xintercept = 0, linewidth = 0.2, colour = "#9E9E9E") + geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = "#9E9E9E", linewidth = 0.6) +
    geom_point(aes(x = rho), size = 1.2, colour = "black") + labs(x = "Spearman rho (gray: permutation 95% range)", y = NULL) + th
  ggsave(file.path(FIG, "H55_C_hardy_split.pdf"), p3, width = 2.6, height = 1.3, useDingbats = FALSE)
  cs <- data.table(ctl = c("PC3 steatosis genes (197)", "NC1 slide tissue area", "Ischemic time (joint model)", "Fat (joint model)"),
                   rho = c(r$PC3$rho, r$NC1$rho, r$Q2$rho_isch, r$Q2$rho_fat), lo = c(r$PC3$p["null_q025"], r$NC1$p["null_q025"], NA, NA), hi = c(r$PC3$p["null_q975"], r$NC1$p["null_q975"], NA, NA))
  cs[, ctl := factor(ctl, levels = rev(ctl))]
  p4 <- ggplot(cs, aes(y = ctl)) + geom_vline(xintercept = 0, linewidth = 0.2, colour = "#9E9E9E") + geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = "#9E9E9E", linewidth = 0.6, na.rm = TRUE) +
    geom_point(aes(x = rho), size = 1.2, colour = "black") + labs(x = "Spearman rho (gray: permutation 95% range)", y = NULL) + th
  ggsave(file.path(FIG, "H55_D_controls.pdf"), p4, width = 2.8, height = 1.4, useDingbats = FALSE)
  message("Caption: (A) bulk disease-vs-control logFC vs GTEx fat-rank slope for ", r$primary$n_genes, " family genes (magenta up, blue down in biopsies), ", r$primary$n_donors,
          " donors; rho ", f3(r$primary$rho), ", Freedman-Lane p ", g3(r$primary$p["p_two"]), ". (B) permutation null of rho, observed in magenta. (C) rho by Hardy group. (D) controls; joint-model points have no permutation range.")
}
cat("report written:", RES_MD, "\n")
