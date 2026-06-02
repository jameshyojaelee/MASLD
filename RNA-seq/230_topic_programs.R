#!/usr/bin/env Rscript
# =============================================================================
# 230_topic_programs.R
# Gap #16: Latent transcriptomic program discovery via NMF topic modeling
#
# Extends the k=2 subtyping NMF (Script 44) to higher k values (5,7,10) to
# discover co-regulated gene programs (topics). For each program: identifies
# top loading genes, runs Hallmark pathway enrichment, tests association with
# fibrosis stage, and identifies switch-like programs activated at F2.
#
# Inputs:
#   - merged_dge.rds (1,444 QC-passing samples)
#   - unified_metadata.csv + sample_qc_report.csv
#
# Outputs:
#   - figures/supplementary/figS08_subtyping_convergence/figS_topic_programs.pdf (4-panel figure)
#   - RNA-seq/results/topic_programs/ (companion CSVs)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)
  library(edgeR)
  library(limma)
  library(NMF)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(grid)
  library(gridExtra)
  library(RColorBrewer)
})

# Ensure dplyr verbs take priority
select <- dplyr::select
filter <- dplyr::filter

# --- Configuration ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

source(file.path(BASE, "scripts/figures/publication_theme.R"))

dge_path  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
meta_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
qc_path   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv")

out_dir <- file.path(BASE, "RNA-seq/results/topic_programs")
fig_dir <- file.path(BASE, "figures/supplementary/figS08_subtyping_convergence")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Parameters
K_VALUES  <- c(5, 7, 10)
NMF_RUNS  <- 30
N_TOP_GENES <- 5000
TOP_N_PER_PROGRAM <- 50
SEED <- 42

cat("=== 230: Latent Transcriptomic Program Discovery (NMF Topic Modeling) ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. Load and Prepare Data (all samples, not just disease)
# =============================================================================
cat("--- Loading data ---\n")

dge <- readRDS(dge_path)
cat(sprintf("  DGEList: %d genes x %d samples\n", nrow(dge), ncol(dge)))

meta <- fread(meta_path)
qc   <- fread(qc_path)

meta <- meta %>%
  left_join(qc %>% select(sample_id, pass_technical), by = "sample_id")

# Use ALL QC-passing samples (control + disease) for program discovery
pass_samples <- meta %>%
  filter(pass_technical == TRUE) %>%
  pull(sample_id)

available <- intersect(pass_samples, colnames(dge))
cat(sprintf("  QC-passing samples in DGEList: %d\n", length(available)))

dge_sub  <- dge[, available]
meta_sub <- meta %>% filter(sample_id %in% available)

# =============================================================================
# 2. Normalize, batch-correct, select variable genes
# =============================================================================
cat("\n--- Normalizing and selecting variable genes ---\n")

dge_sub <- calcNormFactors(dge_sub, method = "TMM")
v <- voom(dge_sub, design = NULL, plot = FALSE)
logcpm <- v$E
cat(sprintf("  Log-CPM: %d genes x %d samples\n", nrow(logcpm), ncol(logcpm)))

# Remove batch effect (dataset)
logcpm_corrected <- removeBatchEffect(logcpm, batch = meta_sub$dataset)
cat("  Batch correction complete\n")

# Select top variable genes by IQR
gene_iqr <- apply(logcpm_corrected, 1, IQR)
top_genes <- names(sort(gene_iqr, decreasing = TRUE))[1:min(N_TOP_GENES, length(gene_iqr))]
mat <- logcpm_corrected[top_genes, ]
cat(sprintf("  Selected %d most variable genes\n", nrow(mat)))

# Shift to non-negative for NMF
mat_nn <- mat - apply(mat, 1, min)
mat_nn[mat_nn < 0] <- 0

# =============================================================================
# 3. Run NMF at k=5, 7, 10
# =============================================================================
cache_path <- file.path(out_dir, "nmf_topic_cache.rds")

if (file.exists(cache_path)) {
  cat("\n--- Loading cached NMF topic results ---\n")
  cached <- readRDS(cache_path)
  nmf_results <- cached$nmf_results
  metrics     <- cached$metrics
  cat("  Loaded k =", paste(names(nmf_results), collapse = ", "), "\n")
} else {
  cat(sprintf("\n--- Running NMF for k = %s (%d runs each) ---\n",
              paste(K_VALUES, collapse = ", "), NMF_RUNS))

  nmf_results <- list()
  metrics <- data.frame()

  for (k in K_VALUES) {
    cat(sprintf("  NMF k=%d ...\n", k))
    set.seed(SEED)
    res <- nmf(mat_nn, rank = k, nrun = NMF_RUNS,
               method = "brunet", seed = "random",
               .options = list(verbose = FALSE))

    nmf_results[[as.character(k)]] <- res

    coph <- cophcor(res)
    disp <- dispersion(res)
    sil_mean <- mean(silhouette(res, what = "consensus")[, "sil_width"])

    metrics <- bind_rows(metrics, data.frame(
      k = k, cophenetic = coph, dispersion = disp, silhouette = sil_mean
    ))
    cat(sprintf("    cophenetic=%.3f, dispersion=%.3f, silhouette=%.3f\n",
                coph, disp, sil_mean))
  }

  saveRDS(list(nmf_results = nmf_results, metrics = metrics), cache_path)
  cat("  Cached to", cache_path, "\n")
}

cat("\nNMF metrics:\n")
print(metrics, row.names = FALSE)

# Select best k: highest cophenetic among tested values
best_k <- metrics$k[which.max(metrics$cophenetic)]
cat(sprintf("\nBest k = %d (cophenetic = %.3f)\n",
            best_k, metrics$cophenetic[metrics$k == best_k]))

nmf_best <- nmf_results[[as.character(best_k)]]

# =============================================================================
# 4. Extract program loadings and sample scores
# =============================================================================
cat("\n--- Extracting program gene loadings and sample scores ---\n")

W <- basis(nmf_best)  # genes x k (gene loadings per program)
H <- coef(nmf_best)   # k x samples (program activity per sample)

colnames(W) <- paste0("Program_", seq_len(best_k))
rownames(H) <- paste0("Program_", seq_len(best_k))

# Top genes per program
program_genes <- list()
all_top_genes_df <- data.frame()

for (p in seq_len(best_k)) {
  pname <- paste0("Program_", p)
  loadings <- sort(W[, p], decreasing = TRUE)
  top <- head(loadings, TOP_N_PER_PROGRAM)
  program_genes[[pname]] <- names(top)

  all_top_genes_df <- bind_rows(all_top_genes_df, data.frame(
    program = pname,
    gene = names(top),
    loading = as.numeric(top),
    rank = seq_along(top)
  ))
}

cat(sprintf("  Extracted top %d genes per program (k=%d)\n",
            TOP_N_PER_PROGRAM, best_k))

# Sample program scores (transpose H)
sample_scores <- as.data.frame(t(H))
sample_scores$sample_id <- colnames(H)

# Merge with metadata
sample_scores <- sample_scores %>%
  left_join(meta_sub %>% select(sample_id, fibrosis_stage, group_binary, dataset),
            by = "sample_id")

# =============================================================================
# 5. Pathway enrichment per program (Hallmark)
# =============================================================================
cat("\n--- Pathway enrichment per program (fgsea, Hallmark) ---\n")

hallmark <- msigdbr(species = "Homo sapiens",
                    collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# For each program, rank ALL genes by loading and run fgsea
pathway_results <- data.frame()

for (p in seq_len(best_k)) {
  pname <- paste0("Program_", p)
  ranks <- sort(W[, p], decreasing = TRUE)

  # Intersect gene names with pathway gene universe
  fres <- fgsea(pathways = hallmark, stats = ranks,
                minSize = 15, maxSize = 500, nPermSimple = 10000)

  fres$program <- pname
  pathway_results <- bind_rows(pathway_results, as.data.frame(fres) %>%
    select(program, pathway, pval, padj, NES, size))
}

# Top pathways per program
top_pathways <- pathway_results %>%
  group_by(program) %>%
  slice_min(padj, n = 5) %>%
  ungroup()

cat("  Top pathways per program:\n")
for (p in unique(top_pathways$program)) {
  cat(sprintf("  %s:\n", p))
  sub <- top_pathways %>% filter(program == p)
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("    %s (NES=%.2f, padj=%.2e)\n",
                gsub("HALLMARK_", "", sub$pathway[i]), sub$NES[i], sub$padj[i]))
  }
}

# =============================================================================
# 6. Fibrosis stage association + switch detection
# =============================================================================
cat("\n--- Testing fibrosis stage association ---\n")

# Only samples with fibrosis staging
fib_scores <- sample_scores %>%
  filter(!is.na(fibrosis_stage)) %>%
  mutate(fibrosis_stage = paste0("F", as.integer(fibrosis_stage)),
         fibrosis_stage = factor(fibrosis_stage,
                                 levels = c("F0", "F1", "F2", "F3", "F4")))

fibrosis_assoc <- data.frame()
switch_scores  <- data.frame()

for (p in seq_len(best_k)) {
  pname <- paste0("Program_", p)

  # Kruskal-Wallis across fibrosis stages
  kw <- tryCatch(
    kruskal.test(fib_scores[[pname]] ~ fib_scores$fibrosis_stage),
    error = function(e) list(statistic = NA, p.value = NA)
  )

  # Spearman correlation with ordinal fibrosis
  fib_num <- as.integer(fib_scores$fibrosis_stage) - 1  # 0-4
  sp <- tryCatch(
    cor.test(fib_scores[[pname]], fib_num, method = "spearman"),
    error = function(e) list(estimate = NA, p.value = NA)
  )

  # Mean score per stage
  stage_means <- fib_scores %>%
    group_by(fibrosis_stage) %>%
    summarise(mean_score = mean(.data[[pname]]), .groups = "drop")

  # Switch detection: ratio of (F2-F4 mean) / (F0-F1 mean)
  all_stages <- c("F0", "F1", "F2", "F3", "F4")
  means_vec <- setNames(
    sapply(all_stages, function(s) {
      v <- stage_means$mean_score[as.character(stage_means$fibrosis_stage) == s]
      if (length(v) == 0) NA_real_ else v[1]
    }),
    all_stages
  )

  early_mean   <- mean(means_vec[c("F0", "F1")], na.rm = TRUE)
  late_mean    <- mean(means_vec[c("F2", "F3", "F4")], na.rm = TRUE)
  switch_ratio <- late_mean / (early_mean + 1e-10)

  # Largest consecutive jump (only non-NA pairs)
  valid_vec  <- means_vec[!is.na(means_vec)]
  gaps       <- if (length(valid_vec) >= 2) diff(valid_vec) else numeric(0)
  if (length(gaps) > 0) {
    max_gap_idx   <- which.max(abs(gaps))
    max_gap_stage <- names(gaps)[max_gap_idx]
    max_gap_value <- as.numeric(gaps[max_gap_idx])
  } else {
    max_gap_stage <- NA_character_
    max_gap_value <- NA_real_
  }

  # Wilcoxon F0-F1 vs F2-F4
  early_vals <- fib_scores %>%
    filter(fibrosis_stage %in% c("F0", "F1")) %>% pull(!!sym(pname))
  late_vals <- fib_scores %>%
    filter(fibrosis_stage %in% c("F2", "F3", "F4")) %>% pull(!!sym(pname))
  wt <- tryCatch(
    wilcox.test(early_vals, late_vals),
    error = function(e) list(p.value = NA)
  )

  fibrosis_assoc <- bind_rows(fibrosis_assoc, data.frame(
    program = pname,
    kw_pvalue = as.numeric(kw$p.value)[1],
    spearman_rho = as.numeric(sp$estimate)[1],
    spearman_pvalue = as.numeric(sp$p.value)[1],
    switch_ratio = switch_ratio,
    max_gap_stage = max_gap_stage,
    max_gap_value = max_gap_value,
    wilcox_f2_split_pvalue = as.numeric(wt$p.value)[1]
  ))

  switch_scores <- bind_rows(switch_scores, stage_means %>%
    mutate(program = pname))
}

# Identify switch programs: max gap at F1->F2 transition AND significant Wilcoxon
switch_programs <- fibrosis_assoc %>%
  filter(max_gap_stage == "F2" & wilcox_f2_split_pvalue < 0.05)

cat(sprintf("\n  Programs with F2 switch: %d / %d\n",
            nrow(switch_programs), best_k))
if (nrow(switch_programs) > 0) {
  cat("  Switch programs:", paste(switch_programs$program, collapse = ", "), "\n")
}

cat("\nFibrosis association summary:\n")
print(fibrosis_assoc %>%
  select(program, kw_pvalue, spearman_rho, max_gap_stage, wilcox_f2_split_pvalue) %>%
  mutate(across(where(is.numeric), ~ signif(., 3))),
  row.names = FALSE)

# =============================================================================
# 7. Save companion CSVs
# =============================================================================
cat("\n--- Saving results ---\n")

fwrite(all_top_genes_df,
       file.path(out_dir, "topic_program_genes.csv"))
fwrite(pathway_results,
       file.path(out_dir, "topic_program_pathways.csv"))
fwrite(fibrosis_assoc,
       file.path(out_dir, "topic_program_fibrosis_assoc.csv"))
fwrite(switch_scores,
       file.path(out_dir, "topic_program_stage_scores.csv"))
fwrite(sample_scores %>% select(sample_id, starts_with("Program_"),
                                 fibrosis_stage, group_binary, dataset),
       file.path(out_dir, "topic_program_sample_scores.csv"))
fwrite(metrics,
       file.path(out_dir, "topic_program_nmf_metrics.csv"))

cat("  Saved 6 CSVs to", out_dir, "\n")

# =============================================================================
# 8. Four-panel figure
# =============================================================================
cat("\n--- Generating 4-panel figure ---\n")

# ---- Panel (a): Program activity heatmap (programs x fibrosis stages) ----
heatmap_data <- switch_scores %>%
  pivot_wider(names_from = fibrosis_stage, values_from = mean_score) %>%
  as.data.frame() %>%
  { rownames(.) <- .$program; .[, colnames(.) != "program"] }

# Scale rows for visualization
heatmap_mat <- as.matrix(heatmap_data[, c("F0", "F1", "F2", "F3", "F4")])
heatmap_scaled <- t(scale(t(heatmap_mat)))

# Order programs by spearman_rho (strongest fibrosis association first)
prog_order <- fibrosis_assoc %>%
  arrange(desc(abs(spearman_rho))) %>%
  pull(program)
heatmap_scaled <- heatmap_scaled[prog_order, ]

# Panel (a): ggplot2 tile heatmap (avoids pheatmap scaling issues)
heatmap_scaled[!is.finite(heatmap_scaled)] <- 0
hm_df <- as.data.frame(heatmap_scaled) %>%
  tibble::rownames_to_column("program") %>%
  pivot_longer(-program, names_to = "stage", values_to = "score") %>%
  mutate(stage = factor(stage, levels = c("F0","F1","F2","F3","F4")),
         program = factor(program, levels = rev(prog_order)))

pa <- ggplot(hm_df, aes(x = stage, y = program, fill = score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                       midpoint = 0, name = "Scaled\nscore") +
  labs(title = "a  Program activity by fibrosis stage", x = NULL, y = NULL) +
  theme_masld() +
  theme(axis.text = element_text(size = 6), legend.key.height = unit(0.4, "cm"))

# ---- Panel (b): Top genes per program (dot plot) ----
# Show top 10 per program
dot_data <- all_top_genes_df %>%
  filter(rank <= 10) %>%
  mutate(program = factor(program, levels = rev(prog_order)),
         gene = reorder(gene, loading))

pb <- ggplot(dot_data, aes(x = loading, y = program, size = loading)) +
  geom_point(color = masld_colors$up, alpha = 0.7) +
  facet_wrap(~ program, scales = "free_y", ncol = 1) +
  geom_text(aes(label = gene), size = 1.8, hjust = -0.1, vjust = 0.5) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  labs(title = "b  Top loading genes per program", x = "NMF loading", y = NULL) +
  theme_masld() +
  theme(strip.text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# ---- Panel (c): Program-pathway correspondence ----
# Significant pathways only
sig_pw <- pathway_results %>%
  filter(padj < 0.05) %>%
  mutate(pathway_short = gsub("HALLMARK_", "", pathway),
         pathway_short = gsub("_", " ", pathway_short),
         pathway_short = tolower(pathway_short))

# Pivot to heatmap
if (nrow(sig_pw) > 0) {
  # Top 3 pathways per program
  top_pw <- sig_pw %>%
    group_by(program) %>%
    slice_min(padj, n = 3) %>%
    ungroup()

  pw_mat <- top_pw %>%
    select(program, pathway_short, NES) %>%
    pivot_wider(names_from = program, values_from = NES, values_fill = 0) %>%
    as.data.frame() %>%
    { rownames(.) <- .$pathway_short; .[, colnames(.) != "pathway_short"] }

  pw_mat <- as.matrix(pw_mat[, intersect(prog_order, colnames(pw_mat))])

  pc_grob <- pheatmap(pw_mat,
    cluster_rows = nrow(pw_mat) >= 2,
    cluster_cols = ncol(pw_mat) >= 2,
    color = colorRampPalette(c("#1565C0", "white", "#C2185B"))(100),
    border_color = "white",
    main = "c  Program-pathway correspondence (NES)",
    fontsize = 7, fontsize_row = 5, fontsize_col = 6,
    silent = TRUE
  )
} else {
  pc_grob <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No significant Hallmark pathways",
             size = 3, color = "grey50") +
    labs(title = "c  Program-pathway correspondence") +
    theme_void() + theme_masld()
}

# ---- Panel (d): Switch program identification ----
# Line plot: program score vs fibrosis stage, highlight switch programs
line_data <- switch_scores %>%
  left_join(fibrosis_assoc %>% select(program, max_gap_stage, wilcox_f2_split_pvalue),
            by = "program") %>%
  mutate(is_switch = max_gap_stage == "F2" & wilcox_f2_split_pvalue < 0.05,
         fibrosis_stage = factor(fibrosis_stage,
                                 levels = c("F0", "F1", "F2", "F3", "F4")))

# Normalize each program to 0-1 for comparable visualization
line_data <- line_data %>%
  group_by(program) %>%
  mutate(score_norm = (mean_score - min(mean_score)) /
           (max(mean_score) - min(mean_score) + 1e-10)) %>%
  ungroup()

pd <- ggplot(line_data, aes(x = fibrosis_stage, y = score_norm,
                             group = program, color = is_switch)) +
  geom_line(aes(linewidth = is_switch), alpha = 0.8) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  annotate("text", x = 2.7, y = 0.95, label = "F2 switch",
           size = 2.2, color = "grey40", hjust = 0) +
  scale_color_manual(values = c("TRUE" = masld_colors$up,
                                "FALSE" = "grey60"),
                     labels = c("TRUE" = "F2 switch", "FALSE" = "Other"),
                     name = NULL) +
  scale_linewidth_manual(values = c("TRUE" = 0.8, "FALSE" = 0.3), guide = "none") +
  labs(title = "d  Switch program identification",
       x = "Fibrosis stage", y = "Normalized program score") +
  theme_masld() +
  theme(legend.position = c(0.15, 0.85))

# ---- Assemble 4-panel figure ----
fig_path <- file.path(fig_dir, "figS_topic_programs.pdf")

pdf(fig_path, width = 10, height = 10)

# Layout: 2x2
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2,
  widths  = unit(c(0.5, 0.5), "npc"),
  heights = unit(c(0.5, 0.5), "npc")
)))

# Panel a (top-left) — ggplot
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pa, vp = viewport(width = 0.95, height = 0.95))
popViewport()

# Panel b (top-right)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pb, vp = viewport(width = 0.95, height = 0.95))
popViewport()

# Panel c (bottom-left) — ggplot or pheatmap grob
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
if (inherits(pc_grob, "gg")) {
  print(pc_grob, vp = viewport(width = 0.95, height = 0.95))
} else {
  grid.draw(pc_grob$gtable)
}
popViewport()

# Panel d (bottom-right)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
print(pd, vp = viewport(width = 0.95, height = 0.95))
popViewport()

dev.off()

cat(sprintf("\n  Figure saved: %s\n", fig_path))
cat("\n=== 230: Topic Programs COMPLETE ===\n")
cat("End:", format(Sys.time()), "\n")
