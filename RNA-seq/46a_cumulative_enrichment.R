#!/usr/bin/env Rscript
# 46a_cumulative_enrichment.R — Cumulative Enrichment + Permutation Complementarity
#
# Replaces ML target prediction (45b) with unsupervised enrichment analysis.
# Proves multi-source integration adds value via:
#   1. Cumulative enrichment curves (multi-source vs single-source rankings)
#   2. Permutation-based complementarity test (p-values per source)
#   3. Incremental greedy source addition (optimal ordering + bootstrap 95% CIs)
#
# Usage: Rscript 46a_cumulative_enrichment.R
# Compute: login node OK (~2-5 min)

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ME   <- file.path(BASE, "RNA-seq/results/multi_evidence")
DRUG <- file.path(BASE, "RNA-seq/results/drug_repurposing")

# ── Load atlas ──────────────────────────────────────────────────────────────
atlas <- fread(file.path(ME, "multi_evidence_atlas.csv"))
N <- nrow(atlas)
cat("Atlas loaded:", N, "genes\n")

# ── 1. Define external validation sets ──────────────────────────────────────
# V1: DGIdb approved druggable genes
dgidb <- fread(file.path(DRUG, "dgidb_drug_gene_interactions.csv"))
v1 <- intersect(unique(dgidb[approved == TRUE, gene]), atlas$human_symbol)

# V2: OpenTargets known drug targets
ot <- fread(file.path(DRUG, "opentargets_known_drugs.csv"))
v2 <- intersect(unique(ot$gene), atlas$human_symbol)

# V3: Conserved
v3 <- atlas[is_conserved == TRUE, human_symbol]

# V4: Strong clinical drug targets
clin <- fread(file.path(DRUG, "clinical_drug_validation_table.csv"))
v4 <- intersect(unique(clin[atlas_support == "Strong", target_gene]), atlas$human_symbol)

# V5: Strong + Moderate clinical drug targets
v5 <- intersect(unique(clin[atlas_support %in% c("Strong", "Moderate"), target_gene]),
                atlas$human_symbol)

val_sets <- list(V1_DGIdb = v1, V2_OpenTargets = v2, V3_Conserved = v3,
                 V4_Strong = v4, V5_StrongMod = v5)
for (vn in names(val_sets))
  cat(sprintf("  %s: %d genes\n", vn, length(val_sets[[vn]])))

# ── 2. Compute continuous source metrics ────────────────────────────────────
# Higher values = stronger evidence. NA → 0.
safe_neglog10 <- function(x) { x[is.na(x)] <- 1; -log10(pmax(x, 1e-300)) }

# S1: Human bulk RNA-seq
atlas[, s1 := safe_neglog10(dream_padj)]
atlas[is.na(dream_padj), s1 := 0]

# S2: Genetic causal evidence (max across all COLOC PP4 + TWAS flags)
# (mr_pval and sceqtl_twas_best_fdr removed 2026-04-22 — MR ditched from paper)
coloc_cols <- grep("_coloc_pp4$", names(atlas), value = TRUE)
atlas[, s2 := 0]
if (length(coloc_cols) > 0)
  atlas[, s2 := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = coloc_cols]
atlas[is.na(s2), s2 := 0]
atlas[!is.na(twas_pval) & twas_pval < 0.05, s2 := pmax(s2, 0.5)]
if ("ieqtl_disease_interaction" %in% names(atlas))
  atlas[ieqtl_disease_interaction == TRUE, s2 := pmax(s2, 0.5)]

# S3: Essentiality (more negative chronos = more essential)
atlas[, s3 := fifelse(!is.na(essentiality_chronos), pmax(0, -essentiality_chronos), 0)]

# S4: Epigenomic (DA peaks + SCENIC + cross-species promoter)
atlas[, s4 := 0]
if ("mouse_da_padj" %in% names(atlas))
  atlas[!is.na(mouse_da_padj) & mouse_da_padj < 1,
        s4 := safe_neglog10(mouse_da_padj)]
if ("hepatocyte_da_padj" %in% names(atlas))
  atlas[!is.na(hepatocyte_da_padj) & hepatocyte_da_padj < 1,
        s4 := pmax(s4, safe_neglog10(hepatocyte_da_padj))]
if ("scenic_regulon_activity_diff" %in% names(atlas))
  atlas[!is.na(scenic_regulon_activity_diff),
        s4 := s4 + abs(scenic_regulon_activity_diff)]
if ("cross_species_promoter_conserved" %in% names(atlas))
  atlas[cross_species_promoter_conserved == TRUE, s4 := s4 + 1]

# S5: Spatial transcriptomics
atlas[, s5 := 0]
if (all(c("spatial_morans_i", "spatial_is_svg") %in% names(atlas)))
  atlas[spatial_is_svg == TRUE & !is.na(spatial_morans_i), s5 := spatial_morans_i]

# S6: Single-cell
atlas[, s6 := 0]
if ("sc_hepatocyte_logFC" %in% names(atlas))
  atlas[!is.na(sc_hepatocyte_logFC), s6 := abs(sc_hepatocyte_logFC)]
if ("sc_n_celltypes_sig" %in% names(atlas))
  atlas[!is.na(sc_n_celltypes_sig), s6 := s6 + sc_n_celltypes_sig]
if ("liana_n_diff_interactions" %in% names(atlas))
  atlas[!is.na(liana_n_diff_interactions), s6 := s6 + liana_n_diff_interactions]

# Compute percentile ranks (0-1) for each source
src_names <- paste0("s", 1:6)
pr_names  <- paste0("pr", 1:6)
for (i in 1:6)
  atlas[, (pr_names[i]) := rank(get(src_names[i]), ties.method = "average") / .N]

# R_multi: sum of 6 percentile ranks
atlas[, r_multi := rowSums(.SD), .SDcols = pr_names]

cat("Source metrics computed\n")

# ── 3. Fold enrichment function ────────────────────────────────────────────
fold_enrich <- function(rank_idx, is_val, k) {
  # rank_idx: integer indices sorted best-first; is_val: boolean vector
  n_hit <- sum(is_val[rank_idx[seq_len(k)]])
  base_rate <- sum(is_val) / length(is_val)
  if (base_rate == 0) return(NA_real_)
  (n_hit / k) / base_rate
}

# ── 4. Cumulative enrichment curves ────────────────────────────────────────
k_vals <- c(50, 100, 200, 500, 1000, 2000, 5000)

# Precompute rank indices for all strategies
strategies <- list(
  R_multi = order(-atlas$r_multi),
  R_S1    = order(-atlas$s1),
  R_S2    = order(-atlas$s2),
  R_S4    = order(-atlas$s4),
  R_count = order(-atlas$layers_active, -atlas$r_multi)
)

enrichment_out <- rbindlist(lapply(names(strategies), function(sn) {
  ridx <- strategies[[sn]]
  rbindlist(lapply(names(val_sets), function(vn) {
    is_v <- atlas$human_symbol %in% val_sets[[vn]]
    rbindlist(lapply(k_vals, function(k) {
      fe <- fold_enrich(ridx, is_v, k)
      data.table(ranking = sn, validation = vn, k = k,
                 fold_enrichment = fe,
                 n_hits = sum(is_v[ridx[seq_len(k)]]),
                 n_validation = sum(is_v))
    }))
  }))
}))

fwrite(enrichment_out, file.path(ME, "cumulative_enrichment.csv"))
cat("Enrichment results saved:", nrow(enrichment_out), "rows\n")

# ── 5. Permutation-based complementarity test ──────────────────────────────
N_PERM <- 1000
K_TEST <- 200
cat("Running permutation test (", N_PERM, "permutations × 6 sources × ",
    length(val_sets), "validation sets)...\n", sep = "")

pr_mat <- as.matrix(atlas[, ..pr_names])  # N × 6 matrix

comp_results <- rbindlist(lapply(1:6, function(si) {
  cat("  Source S", si, "...\n", sep = "")
  # Base score: sum of all pranks EXCEPT si
  base_score <- rowSums(pr_mat[, -si, drop = FALSE])
  si_vals    <- pr_mat[, si]
  all_score  <- base_score + si_vals

  # Rank indices
  ridx_all     <- order(-all_score)
  ridx_without <- order(-base_score)

  rbindlist(lapply(names(val_sets), function(vn) {
    is_v <- atlas$human_symbol %in% val_sets[[vn]]
    if (sum(is_v) == 0) return(NULL)

    fe_all     <- fold_enrich(ridx_all,     is_v, K_TEST)
    fe_without <- fold_enrich(ridx_without, is_v, K_TEST)
    obs_imp    <- fe_all - fe_without

    # Null distribution: permute source si, recompute enrichment with permuted
    null_imp <- vapply(seq_len(N_PERM), function(p) {
      perm_score <- base_score + sample(si_vals)
      perm_ridx  <- order(-perm_score)
      fold_enrich(perm_ridx, is_v, K_TEST) - fe_without
    }, numeric(1))

    pval <- (sum(null_imp >= obs_imp) + 1) / (N_PERM + 1)

    data.table(
      source = paste0("S", si),
      validation = vn,
      enrichment_all = fe_all,
      enrichment_without = fe_without,
      improvement = obs_imp,
      null_mean = mean(null_imp),
      null_sd   = sd(null_imp),
      perm_pvalue = pval
    )
  }))
}))

fwrite(comp_results, file.path(ME, "complementarity_pvalues.csv"))
cat("Complementarity p-values saved:", nrow(comp_results), "rows\n")

# ── 6. Incremental greedy source addition ──────────────────────────────────
cat("Running incremental source addition...\n")
N_BOOT <- 500

addition_out <- rbindlist(lapply(names(val_sets), function(vn) {
  is_v <- atlas$human_symbol %in% val_sets[[vn]]
  if (sum(is_v) == 0) return(NULL)

  selected <- integer(0)
  available <- 1:6

  rbindlist(lapply(1:6, function(step) {
    remaining <- setdiff(available, selected)

    # Try adding each remaining source; pick best
    best_src <- NA_integer_
    best_fe  <- -Inf
    for (src in remaining) {
      cols <- c(selected, src)
      score <- rowSums(pr_mat[, cols, drop = FALSE])
      ridx  <- order(-score)
      fe    <- fold_enrich(ridx, is_v, K_TEST)
      if (!is.na(fe) && fe > best_fe) { best_fe <- fe; best_src <- src }
    }
    selected <<- c(selected, best_src)

    # Bootstrap 95% CI
    cols <- selected
    boot_fe <- vapply(seq_len(N_BOOT), function(b) {
      idx <- sample.int(N, N, replace = TRUE)
      score_b <- rowSums(pr_mat[idx, cols, drop = FALSE])
      ridx_b  <- order(-score_b)
      fold_enrich(ridx_b, is_v[idx], K_TEST)
    }, numeric(1))

    data.table(
      validation = vn, step = step,
      source_added = paste0("S", best_src),
      sources_selected = paste(paste0("S", selected), collapse = "+"),
      fold_enrichment = best_fe,
      ci_lower = quantile(boot_fe, 0.025, na.rm = TRUE),
      ci_upper = quantile(boot_fe, 0.975, na.rm = TRUE)
    )
  }))
}))

fwrite(addition_out, file.path(ME, "source_addition_order.csv"))
cat("Source addition order saved\n")

# ── Summary ─────────────────────────────────────────────────────────────────
cat("\n=== RESULTS SUMMARY ===\n")
cat("\nFold enrichment at top-200 for V1 (DGIdb approved drugs):\n")
print(enrichment_out[validation == "V1_DGIdb" & k == 200,
                      .(ranking, fold_enrichment = round(fold_enrichment, 2))])

cat("\nComplementarity p-values (V1_DGIdb):\n")
print(comp_results[validation == "V1_DGIdb",
                    .(source, improvement = round(improvement, 3), perm_pvalue)])

cat("\nIncremental source addition (V1_DGIdb):\n")
print(addition_out[validation == "V1_DGIdb",
                    .(step, source_added, FE = round(fold_enrichment, 2),
                      CI = paste0("[", round(ci_lower, 2), ", ", round(ci_upper, 2), "]"))])

sig_sources <- comp_results[perm_pvalue < 0.05, .N, by = source]
cat("\nSources with p < 0.05 in at least one validation set:\n")
print(comp_results[perm_pvalue < 0.05, .(source, validation, perm_pvalue)])

cat("\nDone.\n")
