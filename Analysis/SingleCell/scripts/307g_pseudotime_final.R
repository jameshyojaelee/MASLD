#!/usr/bin/env Rscript
#' 307g: Pseudotime Final Figure — 6 panels, one insight each.
#'
#' a: UMAPs (trajectory structure)
#' b: Pseudotime shift per cell type (which types respond most)
#' c: Fibrosis concordance heatmap (UP only, 4 rows)
#' d: NAS concordance heatmap (UP only, 3 rows)
#' e: Sequential fibrosis activation in macrophages
#' f: Top trajectory genes per cell type (dot plot)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime")
FIG_DIR <- file.path(BASE, "figures/supplementary/figS02_progression")

theme_src <- file.path(BASE, "scripts/figures/publication_theme.R")
if (file.exists(theme_src)) source(theme_src)

CELL_TYPES <- c("Hepatocytes", "Macrophages", "Fibroblasts",
                "Endothelial_cells", "Cholangiocytes")
CT_LABELS <- c(Hepatocytes = "Hepatocytes", Macrophages = "Macrophages",
               Fibroblasts = "Fibroblasts", Endothelial_cells = "Endothelial",
               Cholangiocytes = "Cholangiocytes")
CT_ORDER <- c("Hepatocytes", "Macrophages", "Fibroblasts", "Endothelial", "Cholangiocytes")

cond_colors <- c(Healthy = "#2196F3", MASLD = "#E91E63")

cat("=== 307g: Final pseudotime figure ===\n")

# =========================================================================
# Load all data
# =========================================================================
cpt <- fread(file.path(PT_DIR, "consensus_pseudotime_all.csv"))
names(cpt)[1] <- "cell"

meta_list <- lapply(CELL_TYPES, function(ct) {
  f <- file.path(PT_DIR, paste0(ct, "_metadata.csv"))
  if (!file.exists(f)) return(NULL)
  m <- fread(f); names(m)[1] <- "cell"; m$ct <- ct; m
})
meta_all <- rbindlist(meta_list, fill = TRUE)
dat <- merge(meta_all, cpt[, .(cell, consensus_pseudotime)], by = "cell", all.x = TRUE)
dat[, ct_label := CT_LABELS[ct]]
dat[, ct_label := factor(ct_label, levels = CT_ORDER)]

conc <- fread(file.path(PT_DIR, "bulk_sc_concordance.csv"))

scores <- fread(file.path(PT_DIR, "bulk_signature_scores.csv"))
names(scores)[1] <- "cell"
scores_pt <- merge(scores, cpt[, .(cell, consensus_pseudotime)], by = "cell")

# =========================================================================
# Panel A: Compact UMAPs
# =========================================================================
cat("Panel a...\n")
set.seed(42)
dp <- dat[!is.na(consensus_pseudotime), .SD[sample(.N, min(.N, 12000))], by = ct]
dp <- dp[sample(.N)]

pa <- ggplot(dp, aes(x = UMAP_1, y = UMAP_2, color = consensus_pseudotime)) +
  geom_point(size = 0.02, alpha = 0.5, shape = 16) +
  scale_color_viridis_c(option = "magma", direction = -1, name = "Pseudotime",
                         breaks = c(0, 1), labels = c("Early", "Late")) +
  facet_wrap(~ct_label, nrow = 1, scales = "free") +
  theme_masld(base_size = 7) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size = 7, face = "bold"),
        legend.position = "right",
        legend.key.width = unit(0.15, "cm"), legend.key.height = unit(0.6, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        panel.spacing = unit(0.05, "cm"),
        plot.margin = margin(1, 1, 1, 1))

# =========================================================================
# Panel B: Pseudotime shift per cell type (ranked)
# =========================================================================
cat("Panel b...\n")
dat_hm <- dat[condition %in% c("Healthy", "MASLD") & !is.na(consensus_pseudotime)]

shift <- dat_hm[, .(
  median_healthy = median(consensus_pseudotime[condition == "Healthy"]),
  median_masld = median(consensus_pseudotime[condition == "MASLD"]),
  n_healthy = sum(condition == "Healthy"),
  n_masld = sum(condition == "MASLD")
), by = ct_label]
shift[, delta := median_masld - median_healthy]
shift[, ct_label := factor(ct_label, levels = ct_label[order(delta)])]

# Wilcoxon p-values
shift[, pval := sapply(as.character(ct_label), function(ct) {
  h <- dat_hm[ct_label == ct & condition == "Healthy", consensus_pseudotime]
  m <- dat_hm[ct_label == ct & condition == "MASLD", consensus_pseudotime]
  wilcox.test(h, m)$p.value
})]
shift[, sig := ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "ns")))]

pb <- ggplot(shift, aes(x = ct_label, y = delta, fill = delta > 0)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sig), vjust = ifelse(shift$delta > 0, -0.3, 1.3), size = 3) +
  geom_text(aes(label = sprintf("%.2f", delta)),
            vjust = ifelse(shift$delta > 0, 1.5, -0.5), size = 2, color = "white") +
  scale_fill_manual(values = c("TRUE" = "#C2185B", "FALSE" = "#1565C0")) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  coord_flip() +
  theme_masld(base_size = 7) +
  theme(plot.margin = margin(2, 5, 2, 2)) +
  labs(x = NULL, y = "Δ median pseudotime\n(MASLD − Healthy)",
       title = "Disease pseudotime shift")

# =========================================================================
# Panel C: Fibrosis UP concordance (4 rows × 5 cols)
# =========================================================================
cat("Panel c...\n")
conc_pal <- conc[pseudotime_method == "palantir"]
if (nrow(conc_pal) == 0) conc_pal <- conc[pseudotime_method == "dpt"]

# Fibrosis UP only
fib_up <- conc_pal[grepl("F[0-4].*_UP$", signature)]
fib_up[, sig_clean := gsub("transition_", "", signature)]
fib_up[, sig_clean := gsub("_to_", " → ", sig_clean)]
fib_up[, sig_clean := gsub("_UP", "", sig_clean)]
fib_up[, sig_clean := factor(sig_clean,
                               levels = rev(c("F0 → F1", "F1 → F2", "F2 → F3", "F3 → F4")))]
fib_up[, ct_label := CT_LABELS[cell_type]]
fib_up[, ct_label := factor(ct_label, levels = CT_ORDER)]

pc <- ggplot(fib_up, aes(x = ct_label, y = sig_clean, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho),
                color = abs(spearman_rho) > 0.25),
            size = 2, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "gray30")) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                        midpoint = 0, limits = c(-0.5, 0.5), name = "ρ") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.5),
        legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.5, "cm")) +
  labs(x = NULL, y = NULL, title = "Fibrosis signature concordance")

# =========================================================================
# Panel D: NAS UP concordance (3 rows × 5 cols)
# =========================================================================
cat("Panel d...\n")
nas_up <- conc_pal[grepl("NAS.*_UP$", signature)]
nas_up[, sig_clean := gsub("transition_", "", signature)]
nas_up[, sig_clean := gsub("_to_", " → ", sig_clean)]
nas_up[, sig_clean := gsub("_UP", "", sig_clean)]
nas_up[, sig_clean := gsub("NAS01", "NAS 0-1", sig_clean)]
nas_up[, sig_clean := gsub("NAS24", "NAS 2-4", sig_clean)]
nas_up[, sig_clean := gsub("NAS5", "NAS 5", sig_clean)]
nas_up[, sig_clean := gsub("NAS68", "NAS 6-8", sig_clean)]
nas_up[, sig_clean := factor(sig_clean,
                               levels = rev(c("NAS 0-1 → NAS 2-4", "NAS 2-4 → NAS 5",
                                               "NAS 5 → NAS 6-8")))]
nas_up[, ct_label := CT_LABELS[cell_type]]
nas_up[, ct_label := factor(ct_label, levels = CT_ORDER)]

pd <- ggplot(nas_up, aes(x = ct_label, y = sig_clean, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho),
                color = abs(spearman_rho) > 0.25),
            size = 2, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "gray30")) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                        midpoint = 0, limits = c(-0.5, 0.5), name = "ρ") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.5),
        legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.5, "cm")) +
  labs(x = NULL, y = NULL, title = "NAS signature concordance")

# =========================================================================
# Panels E/F: Sequential activation — Hepatocytes + Macrophages
# Fibrosis (e) and NAS (f) programs along pseudotime
# =========================================================================
cat("Panels e/f...\n")

fib_up_cols <- grep("^sig_transition_F.*_UP$", names(scores_pt), value = TRUE)
nas_up_cols <- grep("^sig_transition_NAS.*_UP$", names(scores_pt), value = TRUE)

fib_colors <- c("F0 → F1" = "#90CAF9", "F1 → F2" = "#42A5F5",
                "F2 → F3" = "#E91E63", "F3 → F4" = "#880E4F")
nas_colors <- c("NAS 0-1 → 2-4" = "#90CAF9", "NAS 2-4 → 5" = "#AB47BC",
                "NAS 5 → 6-8" = "#880E4F")

# Helper: bin + aggregate for a cell type and set of signature columns
bin_aggregate <- function(ct_name, sig_cols, n_bins = 25) {
  sub <- scores_pt[cell_type == ct_name & !is.na(consensus_pseudotime)]
  brks <- unique(quantile(sub$consensus_pseudotime,
                           probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(brks) < 3) return(NULL)
  sub[, pt_bin := as.integer(cut(consensus_pseudotime, breaks = brks,
                                  include.lowest = TRUE, labels = FALSE))]
  midpoints <- sub[, .(pt_mid = mean(consensus_pseudotime)), by = pt_bin]
  agg <- list()
  for (col in sig_cols) {
    binned <- sub[, .(mean_score = mean(get(col), na.rm = TRUE)), by = pt_bin]
    binned <- merge(binned, midpoints, by = "pt_bin")
    binned$signature <- col
    binned$cell_type <- ct_name
    agg[[col]] <- binned
  }
  rbindlist(agg)
}

# Fibrosis: both cell types
fib_df <- rbindlist(lapply(c("Hepatocytes", "Macrophages"), bin_aggregate, sig_cols = fib_up_cols))
fib_df[, sig_label := gsub("sig_transition_|_UP", "", signature)]
fib_df[, sig_label := gsub("_to_", " → ", sig_label)]
fib_df[, sig_label := factor(sig_label, levels = c("F0 → F1", "F1 → F2", "F2 → F3", "F3 → F4"))]

pe <- ggplot(fib_df, aes(x = pt_mid, y = mean_score, color = sig_label)) +
  geom_smooth(method = "loess", span = 0.4, se = TRUE, alpha = 0.12, linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  scale_color_manual(values = fib_colors, name = NULL) +
  facet_wrap(~cell_type, nrow = 1) +
  theme_masld(base_size = 7) +
  theme(legend.position = "right", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        strip.text = element_text(size = 8, face = "bold")) +
  labs(x = "Pseudotime →", y = "Fibrosis\nsignature score",
       title = "Fibrosis program activation along pseudotime")

# NAS: both cell types
nas_df <- rbindlist(lapply(c("Hepatocytes", "Macrophages"), bin_aggregate, sig_cols = nas_up_cols))
nas_df[, sig_label := gsub("sig_transition_|_UP", "", signature)]
nas_df[, sig_label := gsub("_to_", " → ", sig_label)]
nas_df[, sig_label := gsub("NAS01", "NAS 0-1", sig_label)]
nas_df[, sig_label := gsub("NAS24", "NAS 2-4", sig_label)]
nas_df[, sig_label := gsub("NAS5 → NAS68", "NAS 5 → 6-8", sig_label)]
nas_df[, sig_label := gsub("NAS24", "NAS 2-4", sig_label)]  # catch remaining
nas_df[, sig_label := gsub("NAS5", "NAS 5", sig_label)]
nas_df[, sig_label := gsub("NAS 5 → NAS 6-8", "NAS 5 → 6-8", sig_label)]
nas_df[, sig_label := gsub("NAS 0-1 → NAS 2-4", "NAS 0-1 → 2-4", sig_label)]
nas_df[, sig_label := gsub("NAS 2-4 → NAS 5", "NAS 2-4 → 5", sig_label)]
nas_df[, sig_label := factor(sig_label, levels = c("NAS 0-1 → 2-4", "NAS 2-4 → 5", "NAS 5 → 6-8"))]

pf_nas <- ggplot(nas_df, aes(x = pt_mid, y = mean_score, color = sig_label)) +
  geom_smooth(method = "loess", span = 0.4, se = TRUE, alpha = 0.12, linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  scale_color_manual(values = nas_colors, name = NULL) +
  facet_wrap(~cell_type, nrow = 1) +
  theme_masld(base_size = 7) +
  theme(legend.position = "right", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        strip.text = element_text(size = 8, face = "bold")) +
  labs(x = "Pseudotime →", y = "NAS\nsignature score",
       title = "NAS program activation along pseudotime")

# =========================================================================
# Panel G: Top trajectory genes per cell type (dot plot)
# =========================================================================
cat("Panel g...\n")

# Comprehensive housekeeping / non-informative gene exclusion
hk_pattern <- paste0(
  "^RPL|^RPS|^MT-|^ENSG|^LINC|",        # ribosomal, mitochondrial, ENSG, lincRNA
  "^HLA-[A-C]$|^B2M$|",                   # MHC-I (ubiquitous)
  "^ACTB$|^ACTG1$|^GAPDH$|^TUBA|^TUBB|", # cytoskeletal / metabolic housekeeping
  "^EEF|^EIF|^UBA52$|^FAU$|^UBB$|^UBC$|", # translation / ubiquitin
  "^FTH1$|^FTL$|^RACK1$|",                # ferritin / ribosome-associated
  "^ATP5|^NDUF|^COX[0-9]|^UQCR|",        # electron transport chain
  "^ALDOA$|^ALDOB$|^ENO1$|^PKM$|^LDHA$|^LDHB$|", # glycolysis
  "^TMSB4X$|^TMSB10$|^PTMA$|^PYURF$|",   # thymosin / small peptides
  "^HSP|^DNAJ|",                           # heat shock
  "^S100A[0-9]|^CALM[0-9]|^VIM$"          # calcium binding / ubiquitous structural
)

corr_list <- list()
for (ct in CELL_TYPES) {
  f <- file.path(PT_DIR, paste0("pseudotime_corr_", ct, ".csv"))
  if (!file.exists(f)) next
  c <- fread(f)
  c$ct <- ct
  c <- c[!grepl(hk_pattern, gene) & padj < 0.05]
  # Top 2 positive + top 1 negative per cell type
  top_pos <- head(c[spearman_rho > 0][order(-abs(spearman_rho))], 2)
  top_neg <- head(c[spearman_rho < 0][order(-abs(spearman_rho))], 1)
  corr_list[[ct]] <- rbind(top_pos, top_neg)
}
gene_df <- rbindlist(corr_list)
gene_df[, ct_label := CT_LABELS[ct]]
gene_df[, ct_label := factor(ct_label, levels = CT_ORDER)]
gene_df[, direction := ifelse(spearman_rho > 0, "Up", "Down")]

# Order genes by cell type then by rho
gene_df[, gene := factor(gene, levels = rev(unique(gene[order(ct_label, -spearman_rho)])))]

pg <- ggplot(gene_df, aes(x = ct_label, y = gene, size = abs(spearman_rho),
                            color = direction)) +
  geom_point() +
  scale_size_continuous(range = c(2, 6), name = "|ρ|",
                         breaks = c(0.2, 0.3, 0.4, 0.5)) +
  scale_color_manual(values = c(Up = "#C2185B", Down = "#1565C0"), name = "Direction") +
  theme_masld(base_size = 7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.5),
        axis.text.y = element_text(face = "italic", size = 6),
        legend.key.size = unit(0.3, "cm")) +
  labs(x = NULL, y = NULL, title = "Top trajectory genes")

# =========================================================================
# Assemble: 4 rows
# Row 1: a (UMAPs, full width)
# Row 2: b | c | d (shift + heatmaps)
# Row 3: e | f (fibrosis + NAS virtual staging, full width each half)
# Row 4: g (gene dot plot, full width)
# =========================================================================
cat("Assembling...\n")

row2 <- pb | pc | pd
row3 <- pe | pf_nas

fig <- pa / row2 / row3 / pg +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

fig <- fig + plot_layout(heights = c(0.9, 1.1, 1.0, 1.0))

out <- file.path(FIG_DIR, "figS_pseudotime_final.pdf")
save_fig(fig, out, width = 7.09, height = 10)
cat("Saved:", out, "\n")

cat("=== 307g COMPLETE ===\n")
