#!/usr/bin/env Rscript
# H17 step 3: prespecified sensitivities and the prespecified decision rule.
# Usage: Rscript 03_sensitivity_and_call.R <HYP_DIR> <OUT_DIR>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
hyp_dir <- args[1]; out_dir <- args[2]
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "h17_lib.R"))
stopifnot(file.exists(file.path(out_dir, "STEP1_DONE")))
wt <- function(x, f) fwrite(x, file.path(out_dir, f), sep = "\t")
rd <- function(f) fread(file.path(out_dir, f))

de_A <- rd("ipsc_de.tsv"); de_B <- rd("invivo_de.tsv")
st <- readRDS(file.path(out_dir, "state_step1.rds"))

concord <- function(dA, dB, label, n_perm = 10000L) {
  U <- intersect(dA$gene, dB$gene)
  a <- dA[match(U, gene)]; b <- dB[match(U, gene)]
  sa <- sign(a$logFC); sb <- sign(b$logFC)
  ia <- which(a$q < Q_THR); ib <- which(b$q < Q_THR)
  obs <- if (length(ia)) mean(sa[ia] == sb[ia]) else NA_real_
  nul <- if (length(ia)) sign_perm_null(sa[ia], sb, ia, n_perm) else NA_real_
  k <- sum(a$q < Q_THR & b$q < Q_THR)
  data.table(analysis = label, n_U = length(U), n_sig_ipsc = length(ia), n_sig_invivo = length(ib),
             S1_agreement = obs, S1_p_upper = if (length(ia)) p_upper(obs, nul) else NA_real_,
             S1_p_lower = if (length(ia)) p_lower(obs, nul) else NA_real_,
             S2_spearman = cor(a$logFC, b$logFC, method = "spearman"),
             n_sig_both = k,
             S3_hypergeom_p = phyper(k - 1, length(ia), length(U) - length(ia), length(ib), lower.tail = FALSE),
             n_core = sum(a$q < Q_THR & b$q < Q_THR & sa == sb))
}

sens <- list()
# 1. Steatosis-only in vivo arm
B1 <- load_invivo("steatosis_only")
wt(B1$meta[, .N, by = .(dataset, group)], "sens1_steatosis_only_arm_counts.tsv")
if (all(B1$meta[, .N, by = group]$N >= 3) && uniqueN(B1$meta$group) == 2) {
  g1 <- filter_genes(B1$counts, B1$meta$group)
  de1 <- fit_full(B1$counts, g1, invivo_design(B1$meta), "groupsteatotic")
  wt(de1, "sens1_invivo_steatosis_only_de.tsv")
  sens[[1]] <- concord(de_A, de1, "sens1_invivo_steatosis_only")
} else sens[[1]] <- data.table(analysis = "sens1_invivo_steatosis_only: untestable (<3 per arm)")
rm(B1); invisible(gc())

# 2. iPSC without non-responders
A <- load_ipsc(hyp_dir)
nr <- st$nonresp
if (length(nr) > 0) {
  keep <- !(as.character(A$sheet$participant) %in% nr)
  sh <- A$sheet[keep]; sh[, participant := droplevels(participant)]
  cnt <- A$counts[, keep]
  g2 <- filter_genes(cnt, sh$treatment)
  de2 <- fit_full(cnt, g2, ipsc_design(sh$participant, sh$treatment), "treatmentHepFA")
  wt(de2, "sens2_ipsc_without_nonresponders_de.tsv")
  sens[[2]] <- concord(de2, de_B, sprintf("sens2_ipsc_without_%d_nonresponders", length(nr)))
} else sens[[2]] <- data.table(analysis = "sens2_ipsc_without_nonresponders: not run (0 non-responders)")
rm(A); invisible(gc())

# 3. in vivo with inferred sex
B <- load_invivo("steatotic")
if (!anyNA(B$meta$sex) && uniqueN(B$meta$sex) == 2) {
  de3 <- fit_full(B$counts, st$genes_B, invivo_design(B$meta, add_sex = TRUE), "groupsteatotic")
  wt(de3, "sens3_invivo_sex_adjusted_de.tsv")
  sens[[3]] <- concord(de_A, de3, "sens3_invivo_sex_adjusted")
} else sens[[3]] <- data.table(analysis = sprintf(
  "sens3_invivo_sex_adjusted: not run (sex missing for %d donors or one sex only)", sum(is.na(B$meta$sex))))
wt(rbindlist(sens, fill = TRUE), "sensitivity_concordance.tsv")

# ---------------------------------------------------------------- decision rule
conc <- rd("concordance_primary.tsv")
pc <- rd("positive_control_ipsc.tsv")
crit <- data.table(criterion = c("positive_control_pass", "S1_agreement_gt_0.5_and_signperm_p_lt_0.05",
                                 "S1_gt_matchedset_q95", "S1_gt_invivo_labelperm_q95",
                                 "S2_bootstrap_CI_lower_gt_0", "core_nonempty_and_S3_p_lt_0.05"))
labq95 <- bootci <- c(NA_real_, NA_real_)
if (file.exists(file.path(out_dir, "STEP2_DONE"))) {
  lp <- rd("negctrl_invivo_label_permutation.tsv")
  labq95 <- quantile(lp$S1_agreement, 0.95)
  bs <- rd("bootstrap_S2.tsv")
  bootci <- quantile(bs$S2_spearman, c(0.025, 0.975))
}
crit[, value := c(
  as.character(pc$overall_pass[1]),
  sprintf("%.4f (p=%.4g)", conc$S1_agreement, conc$S1_p_upper),
  sprintf("%.4f vs %.4f", conc$S1_agreement, conc$S1_matchedset_q95),
  sprintf("%.4f vs %.4f", conc$S1_agreement, labq95),
  sprintf("rho %.4f [%.4f, %.4f]", conc$S2_spearman, bootci[1], bootci[2]),
  sprintf("core %d, p=%.4g", conc$n_core, conc$S3_hypergeom_p))]
crit[, met := c(
  isTRUE(pc$overall_pass[1]),
  conc$S1_agreement > 0.5 && conc$S1_p_upper < 0.05,
  conc$S1_agreement > conc$S1_matchedset_q95,
  !is.na(labq95) && conc$S1_agreement > labq95,
  !is.na(bootci[1]) && bootci[1] > 0,
  conc$n_core > 0 && conc$S3_hypergeom_p < 0.05)]
call <- if (!crit$met[1]) "untestable" else if (all(crit$met)) "yes" else "no"
opposite <- conc$S1_p_lower < 0.05
crit <- rbind(crit, data.table(criterion = c("CALL", "S1_opposite_direction_lower_tail_p_lt_0.05"),
                               value = c(call, sprintf("p_lower=%.4g", conc$S1_p_lower)),
                               met = c(NA, opposite)))
wt(crit, "decision.tsv")
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
writeLines("ok", file.path(out_dir, "STEP3_DONE"))
message("call: ", call)
