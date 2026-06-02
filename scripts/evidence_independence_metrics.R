#!/usr/bin/env Rscript
# Evidence source independence metrics for the multi-evidence atlas
# Goal: Replace MI=0.001 with concrete, intuitive overlap statistics

suppressPackageStartupMessages({
  library(data.table)
})

cat("======================================================================\n")
cat("EVIDENCE SOURCE INDEPENDENCE ANALYSIS\n")
cat("======================================================================\n\n")

# ── 1. Load atlas ──────────────────────────────────────────────────
atlas <- fread("RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
cat(sprintf("Atlas: %s genes x %s columns\n\n", format(nrow(atlas), big.mark=","), ncol(atlas)))

# ── 2. Define binary evidence layers ──────────────────────────────

# S1: Human bulk transcriptomic (dream mega-analysis)
atlas[, s1_deg := !is.na(dream_padj) & dream_padj < 0.1]

# S3: Genetic causal — multiple definitions for robustness
# (a) SuSiE-COLOC best PP4 > 0.5 (fine-mapped, 24 EUR GWAS)
atlas[, s3_susie := !is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 > 0.5]
# (b) Any COLOC source with PP4 > 0.5 (broadest definition)
atlas[, s3_any_coloc := n_coloc_sources >= 1]
# (c) Causal robustness = "Robust" (multi-method validated)
atlas[, s3_robust := causal_robustness == "Robust"]

# S2: Cross-species (Conserved from mouse concordance)
atlas[, s2_mouse := !is.na(is_conserved) & is_conserved == TRUE]

# Mouse DEG (any mouse significance)
atlas[, s2_mouse_deg := !is.na(mouse_meta_padj) & mouse_meta_padj < 0.1]

# Deconvolution attribution
atlas[, s_hepatocyte := attribution_class == "Hepatocyte_intrinsic"]

cat(sprintf("Layer sizes:\n"))
cat(sprintf("  S1: DEGs (padj<0.1)           : %6d genes\n", sum(atlas$s1_deg)))
cat(sprintf("  S3: SuSiE COLOC (PP4>0.5)     : %6d genes\n", sum(atlas$s3_susie)))
cat(sprintf("  S3: Any COLOC source           : %6d genes\n", sum(atlas$s3_any_coloc)))
cat(sprintf("  S3: Causal robust              : %6d genes\n", sum(atlas$s3_robust)))
cat(sprintf("  S2: Conserved             : %6d genes\n", sum(atlas$s2_mouse)))
cat(sprintf("  S2: Mouse DEG (padj<0.1)       : %6d genes\n", sum(atlas$s2_mouse_deg)))
cat(sprintf("  Hepatocyte-intrinsic           : %6d genes\n", sum(atlas$s_hepatocyte)))

# ── 3. Helper function ────────────────────────────────────────────

compute_overlap <- function(set_a, set_b, name_a, name_b,
                            cont_a=NULL, cont_b=NULL) {
  cat(sprintf("\n%s\n", paste(rep("-", 70), collapse="")))
  cat(sprintf("%s  vs  %s\n", name_a, name_b))
  cat(sprintf("%s\n", paste(rep("-", 70), collapse="")))

  n_a <- sum(set_a, na.rm=TRUE)
  n_b <- sum(set_b, na.rm=TRUE)
  both <- sum(set_a & set_b, na.rm=TRUE)
  either <- sum(set_a | set_b, na.rm=TRUE)
  a_only <- sum(set_a & !set_b, na.rm=TRUE)
  b_only <- sum(!set_a & set_b, na.rm=TRUE)
  neither <- sum(!set_a & !set_b, na.rm=TRUE)

  jaccard <- if(either > 0) both / either else NA
  pct_a_in_b <- if(n_a > 0) 100 * both / n_a else NA
  pct_b_in_a <- if(n_b > 0) 100 * both / n_b else NA

  cat(sprintf("\n  %-25s: %6d genes\n", name_a, n_a))
  cat(sprintf("  %-25s: %6d genes\n", name_b, n_b))
  cat(sprintf("  Both                   : %6d genes\n", both))
  cat(sprintf("  %s only         : %6d genes\n", substr(name_a,1,20), a_only))
  cat(sprintf("  %s only         : %6d genes\n", substr(name_b,1,20), b_only))
  cat(sprintf("  Neither                : %6d genes\n", neither))

  cat(sprintf("\n  Jaccard index     = %.4f  (0 = disjoint, 1 = identical sets)\n", jaccard))
  cat(sprintf("  %% of %s in %s = %.1f%%\n", name_a, name_b, pct_a_in_b))
  cat(sprintf("  %% of %s in %s = %.1f%%\n", name_b, name_a, pct_b_in_a))

  # Fisher's exact test
  mat <- matrix(c(both, a_only, b_only, neither), nrow=2)
  ft <- fisher.test(mat)
  cat(sprintf("\n  Fisher's exact: OR = %.2f [%.2f, %.2f], p = %s\n",
              ft$estimate, ft$conf.int[1], ft$conf.int[2],
              format.pval(ft$p.value, digits=3)))
  if(ft$estimate > 1 & ft$p.value < 0.05) {
    cat(sprintf("  --> ENRICHED (%.1fx), but overlap is still small\n", ft$estimate))
  } else if(ft$estimate < 1 & ft$p.value < 0.05) {
    cat(sprintf("  --> DEPLETED (%.1fx less likely)\n", 1/ft$estimate))
  } else {
    cat(sprintf("  --> No significant enrichment or depletion\n"))
  }

  # Continuous correlation if available
  if(!is.null(cont_a) & !is.null(cont_b)) {
    # Use ALL genes that have both values (including zero PP4)
    both_have <- !is.na(cont_a) & !is.na(cont_b)
    if(sum(both_have) > 10) {
      sp <- cor.test(cont_a[both_have], cont_b[both_have], method="spearman")
      cat(sprintf("\n  Spearman rho = %.4f (n=%s genes with both values), p = %s\n",
                  sp$estimate, format(sum(both_have), big.mark=","),
                  format.pval(sp$p.value, digits=3)))
    }
  }

  invisible(list(n_a=n_a, n_b=n_b, both=both, jaccard=jaccard, or=ft$estimate, pval=ft$p.value,
                 pct_a_in_b=pct_a_in_b, pct_b_in_a=pct_b_in_a))
}

# ── 4. THE KEY COMPARISON: S1 (transcriptomic) vs S3 (genetic) ───
cat("\n\n")
cat("======================================================================\n")
cat("SECTION 1: TRANSCRIPTOMIC vs GENETIC CAUSAL (the headline comparison)\n")
cat("======================================================================\n")

res_main <- compute_overlap(
  atlas$s1_deg, atlas$s3_susie,
  "DEGs", "COLOC PP4>0.5",
  cont_a = atlas$dream_tstat,
  cont_b = atlas$coloc_susie_best_pp4
)

# Also with broader "any COLOC" definition
res_any <- compute_overlap(
  atlas$s1_deg, atlas$s3_any_coloc,
  "DEGs", "Any COLOC source",
  cont_a = atlas$dream_tstat,
  cont_b = atlas$coloc_susie_best_pp4
)

# Also with robust causal
res_robust <- compute_overlap(
  atlas$s1_deg, atlas$s3_robust,
  "DEGs", "Causal-robust"
)

# ── 5. OTHER KEY PAIRS ───────────────────────────────────────────
cat("\n\n")
cat("======================================================================\n")
cat("SECTION 2: OTHER PAIRWISE COMPARISONS\n")
cat("======================================================================\n")

res_s1_cc <- compute_overlap(
  atlas$s1_deg, atlas$s2_mouse,
  "DEGs", "Conserved"
)

res_s1_mdeg <- compute_overlap(
  atlas$s1_deg, atlas$s2_mouse_deg,
  "DEGs (human)", "DEGs (mouse)"
)

res_s3_cc <- compute_overlap(
  atlas$s3_susie, atlas$s2_mouse,
  "COLOC PP4>0.5", "Conserved"
)

res_s3_mdeg <- compute_overlap(
  atlas$s3_susie, atlas$s2_mouse_deg,
  "COLOC PP4>0.5", "Mouse DEGs"
)

res_s1_hep <- compute_overlap(
  atlas$s1_deg, atlas$s_hepatocyte,
  "DEGs", "Hepatocyte-intrinsic"
)

res_s3_hep <- compute_overlap(
  atlas$s3_susie, atlas$s_hepatocyte,
  "COLOC PP4>0.5", "Hepatocyte-intrinsic"
)

# ── 6. CONTINUOUS SPEARMAN MATRIX ─────────────────────────────────
cat("\n\n")
cat("======================================================================\n")
cat("SECTION 3: CONTINUOUS CORRELATIONS (Spearman rho)\n")
cat("======================================================================\n\n")

# For the 6,200 genes tested in SuSiE-COLOC, compute continuous correlations
has_coloc <- !is.na(atlas$coloc_susie_best_pp4)
has_dream <- !is.na(atlas$dream_tstat)
has_mouse <- !is.na(atlas$mouse_meta_logFC)
has_twas  <- !is.na(atlas$twas_z) & atlas$twas_z != 0

cat(sprintf("Genes tested in SuSiE-COLOC: %d\n", sum(has_coloc)))
cat(sprintf("Genes with dream t-stat: %d\n", sum(has_dream)))
cat(sprintf("Genes with mouse meta logFC: %d\n", sum(has_mouse)))
cat(sprintf("Genes with TWAS z: %d\n\n", sum(has_twas)))

# dream_tstat vs coloc_susie_best_pp4
idx <- has_coloc & has_dream
sp1 <- cor.test(atlas$dream_tstat[idx], atlas$coloc_susie_best_pp4[idx], method="spearman")
cat(sprintf("dream_tstat vs coloc_susie_best_pp4  (n=%s): rho = %+.4f, p = %s\n",
            format(sum(idx), big.mark=","), sp1$estimate, format.pval(sp1$p.value, digits=3)))

# dream_tstat vs |coloc_susie_best_pp4| (test if DE magnitude predicts COLOC)
sp1b <- cor.test(abs(atlas$dream_tstat[idx]), atlas$coloc_susie_best_pp4[idx], method="spearman")
cat(sprintf("|dream_tstat| vs coloc_susie_best_pp4 (n=%s): rho = %+.4f, p = %s\n",
            format(sum(idx), big.mark=","), sp1b$estimate, format.pval(sp1b$p.value, digits=3)))

# dream_tstat vs mouse_meta_logFC
idx2 <- has_dream & has_mouse
sp2 <- cor.test(atlas$dream_tstat[idx2], atlas$mouse_meta_logFC[idx2], method="spearman")
cat(sprintf("dream_tstat vs mouse_meta_logFC      (n=%s): rho = %+.4f, p = %s\n",
            format(sum(idx2), big.mark=","), sp2$estimate, format.pval(sp2$p.value, digits=3)))

# coloc_susie_best_pp4 vs mouse_meta_logFC (among genes with both)
idx3 <- has_coloc & has_mouse
sp3 <- cor.test(atlas$coloc_susie_best_pp4[idx3], abs(atlas$mouse_meta_logFC[idx3]), method="spearman")
cat(sprintf("coloc_susie_best_pp4 vs |mouse_LFC|  (n=%s): rho = %+.4f, p = %s\n",
            format(sum(idx3), big.mark=","), sp3$estimate, format.pval(sp3$p.value, digits=3)))

# TWAS vs dream
idx4 <- has_twas & has_dream
sp4 <- cor.test(abs(atlas$twas_z[idx4]), abs(atlas$dream_tstat[idx4]), method="spearman")
cat(sprintf("|twas_z| vs |dream_tstat|             (n=%s): rho = %+.4f, p = %s\n",
            format(sum(idx4), big.mark=","), sp4$estimate, format.pval(sp4$p.value, digits=3)))

# TWAS vs COLOC
idx5 <- has_twas & has_coloc
sp5 <- cor.test(abs(atlas$twas_z[idx5]), atlas$coloc_susie_best_pp4[idx5], method="spearman")
cat(sprintf("|twas_z| vs coloc_susie_best_pp4      (n=%s): rho = %+.4f, p = %s\n",
            format(sum(idx5), big.mark=","), sp5$estimate, format.pval(sp5$p.value, digits=3)))

# ── 7. JACCARD + FISHER MATRICES ─────────────────────────────────
cat("\n\n")
cat("======================================================================\n")
cat("SECTION 4: PAIRWISE JACCARD INDEX MATRIX\n")
cat("(0 = completely disjoint sets, 1 = identical sets)\n")
cat("======================================================================\n\n")

layers <- list(
  "S1: Human DEGs" = atlas$s1_deg,
  "S2: Conserved" = atlas$s2_mouse,
  "S3: COLOC PP4>0.5" = atlas$s3_susie,
  "Hepatocyte-intr." = atlas$s_hepatocyte
)

# Separate list for summary stats (only truly independent layers)
indep_layers <- list(
  "S1: Human DEGs" = atlas$s1_deg,
  "S2: Conserved" = atlas$s2_mouse,
  "S3: COLOC PP4>0.5" = atlas$s3_susie
)

jmat <- matrix(NA, length(layers), length(layers),
               dimnames=list(names(layers), names(layers)))
for(i in seq_along(layers)) {
  for(j in seq_along(layers)) {
    both <- sum(layers[[i]] & layers[[j]], na.rm=TRUE)
    either <- sum(layers[[i]] | layers[[j]], na.rm=TRUE)
    jmat[i,j] <- if(either > 0) both / either else NA
  }
}

cat(sprintf("%-22s", ""))
for(n in colnames(jmat)) cat(sprintf("%18s", n))
cat("\n")
for(i in 1:nrow(jmat)) {
  cat(sprintf("%-22s", rownames(jmat)[i]))
  for(j in 1:ncol(jmat)) {
    if(i == j) { cat(sprintf("%18s", "1.000")); next }
    cat(sprintf("%18.4f", jmat[i,j]))
  }
  cat("\n")
}

cat("\n")
cat("======================================================================\n")
cat("SECTION 5: PAIRWISE FISHER ODDS-RATIO MATRIX\n")
cat("(>1 = enriched co-occurrence, <1 = depleted)\n")
cat("======================================================================\n\n")

or_mat <- matrix(NA, length(layers), length(layers),
                 dimnames=list(names(layers), names(layers)))
pv_mat <- matrix(NA, length(layers), length(layers),
                 dimnames=list(names(layers), names(layers)))

for(i in seq_along(layers)) {
  for(j in seq_along(layers)) {
    if(i == j) { or_mat[i,j] <- Inf; pv_mat[i,j] <- 0; next }
    a <- layers[[i]]; b <- layers[[j]]
    mat <- matrix(c(sum(a&b, na.rm=T), sum(a&!b, na.rm=T),
                     sum(!a&b, na.rm=T), sum(!a&!b, na.rm=T)), nrow=2)
    ft <- fisher.test(mat)
    or_mat[i,j] <- ft$estimate
    pv_mat[i,j] <- ft$p.value
  }
}

cat(sprintf("%-22s", ""))
for(n in colnames(or_mat)) cat(sprintf("%18s", n))
cat("\n")
for(i in 1:nrow(or_mat)) {
  cat(sprintf("%-22s", rownames(or_mat)[i]))
  for(j in 1:ncol(or_mat)) {
    if(i == j) { cat(sprintf("%18s", "---")); next }
    stars <- ifelse(pv_mat[i,j] < 1e-10, "***",
             ifelse(pv_mat[i,j] < 1e-5, "**",
             ifelse(pv_mat[i,j] < 0.05, "*", "ns")))
    cat(sprintf("%14.2f %3s", or_mat[i,j], stars))
  }
  cat("\n")
}
cat("(*** p<1e-10, ** p<1e-5, * p<0.05, ns = not significant)\n")

# ── 8. CONVERGENCE PYRAMID ───────────────────────────────────────
cat("\n\n")
cat("======================================================================\n")
cat("SECTION 6: CONVERGENCE PYRAMID\n")
cat("(Using S1=DEG, S2=Conserved, S3=COLOC PP4>0.5)\n")
cat("======================================================================\n\n")

atlas[, n_top3 := as.integer(s1_deg) + as.integer(s2_mouse) + as.integer(s3_susie)]

cat("Distribution of evidence sources per gene (top 3 independent layers):\n")
src_tab <- table(atlas$n_top3)
for(i in names(src_tab)) {
  cat(sprintf("  %s source(s): %6d genes (%5.1f%%)\n", i, src_tab[i], 100*src_tab[i]/nrow(atlas)))
}

triple <- atlas[n_top3 == 3]
cat(sprintf("\nGenes with ALL THREE sources (DEG + Conserved + COLOC): %d\n", nrow(triple)))
if(nrow(triple) > 0) {
  setorder(triple, -dream_tstat)
  cat("  Genes: ", paste(triple$human_symbol, collapse=", "), "\n")
}

# Also show the 2-source combinations
cat("\n2-source combinations:\n")
cat(sprintf("  DEG + COLOC (no mouse)       : %4d genes\n",
            sum(atlas$s1_deg & atlas$s3_susie & !atlas$s2_mouse)))
cat(sprintf("  DEG + Conserved (no COLOC): %4d genes\n",
            sum(atlas$s1_deg & atlas$s2_mouse & !atlas$s3_susie)))
cat(sprintf("  COLOC + Conserved (no DEG): %4d genes\n",
            sum(atlas$s3_susie & atlas$s2_mouse & !atlas$s1_deg)))

# ── 9. ASYMMETRY CHECK ───────────────────────────────────────────
cat("\n\n")
cat("======================================================================\n")
cat("SECTION 7: ASYMMETRY — WHY EACH LAYER MATTERS\n")
cat("======================================================================\n\n")

# COLOC genes that are NOT DEGs — these are genetically causal but
# transcriptomically silent (or below threshold)
coloc_not_deg <- atlas[s3_susie == TRUE & s1_deg == FALSE]
cat(sprintf("COLOC genes that are NOT DEGs: %d / %d (%.1f%%)\n",
            nrow(coloc_not_deg), sum(atlas$s3_susie),
            100*nrow(coloc_not_deg)/sum(atlas$s3_susie)))
cat("  (These genes have genetic causal evidence but no transcriptomic change)\n")
cat("  Top 10 by PP4: ", paste(head(coloc_not_deg[order(-coloc_susie_best_pp4)]$human_symbol, 10),
                               collapse=", "), "\n\n")

# DEGs that are NOT COLOC — transcriptomically active but not genetically causal
deg_not_coloc <- atlas[s1_deg == TRUE & s3_susie == FALSE]
cat(sprintf("DEGs that are NOT COLOC: %d / %d (%.1f%%)\n",
            nrow(deg_not_coloc), sum(atlas$s1_deg),
            100*nrow(deg_not_coloc)/sum(atlas$s1_deg)))
cat("  (These genes are differentially expressed but lack genetic causal support)\n\n")

# Conserved that are NOT COLOC
cc_not_coloc <- atlas[s2_mouse == TRUE & s3_susie == FALSE]
cat(sprintf("Conserved genes NOT in COLOC: %d / %d (%.1f%%)\n",
            nrow(cc_not_coloc), sum(atlas$s2_mouse),
            100*nrow(cc_not_coloc)/sum(atlas$s2_mouse)))
cat("  (Cross-species conserved but not genetically driven in humans)\n\n")

# COLOC that are NOT Conserved
coloc_not_cc <- atlas[s3_susie == TRUE & s2_mouse == FALSE]
cat(sprintf("COLOC genes NOT in Conserved: %d / %d (%.1f%%)\n",
            nrow(coloc_not_cc), sum(atlas$s3_susie),
            100*nrow(coloc_not_cc)/sum(atlas$s3_susie)))
cat("  (Genetically causal in humans but not conserved in mouse)\n")

# ── 10. PUBLISHABLE SUMMARY ──────────────────────────────────────
cat("\n\n")
cat("======================================================================\n")
cat("PUBLISHABLE SUMMARY (copy-paste ready)\n")
cat("======================================================================\n\n")

cat(sprintf(paste0(
  "Of %s human DEGs (dream padj<0.1), only %d (%.1f%%) also had genetic\n",
  "causal support via COLOC (PP.H4>0.5 across 24 GWAS; Jaccard index = %.3f),\n",
  "and the Spearman correlation between dream t-statistics and COLOC posterior\n",
  "probabilities was rho = %.3f (n=%s genes), confirming that transcriptomic\n",
  "dysregulation and genetic causality nominate largely independent gene sets.\n"),
  format(res_main$n_a, big.mark=","),
  res_main$both,
  res_main$pct_a_in_b,
  res_main$jaccard,
  sp1$estimate,
  format(sum(has_coloc & has_dream), big.mark=",")
))

cat(sprintf(paste0(
  "\nConversely, %.1f%% of COLOC genes were also DEGs, and the two sets showed\n",
  "modest but significant enrichment (Fisher OR = %.2f, p = %s), indicating\n",
  "that convergent genes exist but represent a small, high-confidence subset.\n"),
  res_main$pct_b_in_a,
  res_main$or,
  format.pval(res_main$pval, digits=3)
))

cat(sprintf(paste0(
  "\nCross-species conservation added a third orthogonal axis: only %d of %d\n",
  "Conserved genes (%.1f%%) had COLOC support (Jaccard = %.3f), and\n",
  "only %d of %s COLOC genes (%.1f%%) were cross-species conserved.\n"),
  res_s3_cc$both, res_s1_cc$n_b,
  res_s3_cc$pct_b_in_a, res_s3_cc$jaccard,
  res_s3_cc$both, format(res_s3_cc$n_a, big.mark=","),
  res_s3_cc$pct_a_in_b
))

# Compute max Jaccard among truly independent layers only
jmat_indep <- matrix(NA, length(indep_layers), length(indep_layers),
                     dimnames=list(names(indep_layers), names(indep_layers)))
for(i in seq_along(indep_layers)) {
  for(j in seq_along(indep_layers)) {
    both_ij <- sum(indep_layers[[i]] & indep_layers[[j]], na.rm=TRUE)
    either_ij <- sum(indep_layers[[i]] | indep_layers[[j]], na.rm=TRUE)
    jmat_indep[i,j] <- if(either_ij > 0) both_ij / either_ij else NA
  }
}
max_jaccard_indep <- max(jmat_indep[lower.tri(jmat_indep)], na.rm=TRUE)

cat(sprintf(paste0(
  "\nAll pairwise Jaccard indices between the three independent evidence\n",
  "layers were <=%.3f, and only %d genes (%.2f%%) achieved convergence\n",
  "across all three sources (DEG + COLOC + Conserved), supporting\n",
  "the thesis that evidence convergence — not any single analysis — most\n",
  "reliably identifies therapeutic targets.\n"),
  max_jaccard_indep,
  nrow(triple),
  100*nrow(triple)/nrow(atlas)
))

cat("\nDone.\n")
