#!/usr/bin/env Rscript
# Figures for library vs Integrated DEG overlap comparison
# Only library genes with human orthologs (excl. mouse-only) are used.
# All plots operate on the same set of unique human Ensembl IDs.
# Outputs to results/library/

library(data.table)
library(ggplot2)

outdir <- "results/library"

# --- Load data ---
lib_core <- fread("results/library/final_core_degs.csv")
dream <- fread("RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv")
dream[, ensembl_id := sub("\\.[0-9]+$", "", gene)]

# --- Build unique set of human orthologs from library (excl. mouse-only) ---
lib <- lib_core[has_human_data == TRUE]
lib_human_ids <- unique(unlist(strsplit(lib$human_ortholog_ids, ";")))
lib_human_ids <- lib_human_ids[lib_human_ids != "" & !is.na(lib_human_ids)]
n_lib <- length(lib_human_ids)

cat(sprintf("Library: %d unique human orthologs (from %d mouse genes with human data)\n",
            n_lib, nrow(lib)))

# Shared theme
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "bottom"
  )

cutoffs <- c(0, 0.25, 0.5, 0.75, 1.0)

# =============================================
# Fig 1: Overlap bar chart at LFC cutoffs (UP only)
# =============================================
overlap_data <- rbindlist(lapply(cutoffs, function(lfc_cut) {
  dream_up <- dream[padj < 0.1 & logFC >= lfc_cut, ensembl_id]
  overlap <- intersect(lib_human_ids, dream_up)
  lib_only <- setdiff(lib_human_ids, dream_up)
  dream_only <- setdiff(dream_up, lib_human_ids)
  data.table(
    cutoff = lfc_cut,
    category = factor(c("Overlap", "Library only", "Integrated only"),
                      levels = c("Integrated only", "Overlap", "Library only")),
    count = c(length(overlap), length(lib_only), length(dream_only))
  )
}))

p1 <- ggplot(overlap_data, aes(x = factor(cutoff), y = count, fill = category)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c("Integrated only" = "#4393C3", "Overlap" = "#D6604D",
                                "Library only" = "#FDDBC7"), name = NULL) +
  labs(x = "log2FC cutoff (upregulated, padj < 0.1)",
       y = "Number of genes",
       title = "Library vs Integrated upregulated DEG overlap",
       subtitle = sprintf("Library: %s human orthologs | Integrated: 10-cohort Disease vs Control",
                           format(n_lib, big.mark = ","))) +
  theme_pub +
  theme(legend.position = "top")

ggsave(file.path(outdir, "library_vs_dream_overlap_bars.pdf"), p1,
       width = 7, height = 5)

# =============================================
# Fig 2: Percent of library recovered (UP only)
#         with Integrated DEG count annotated
# =============================================
pct_data <- rbindlist(lapply(cutoffs, function(lfc_cut) {
  dream_up <- dream[padj < 0.1 & logFC >= lfc_cut, ensembl_id]
  n_dream_up <- length(dream_up)
  hits <- length(intersect(lib_human_ids, dream_up))
  data.table(cutoff = lfc_cut,
             pct = 100 * hits / n_lib,
             n = hits,
             n_dream = n_dream_up)
}))

p2 <- ggplot(pct_data, aes(x = cutoff, y = pct)) +
  geom_line(linewidth = 1, color = "#2166AC") +
  geom_point(size = 2.5, color = "#2166AC") +
  geom_text(aes(label = n), vjust = 2, size = 3, color = "#2166AC") +
  geom_text(aes(x = cutoff, y = 100,
                label = format(n_dream, big.mark = ",")),
            vjust = -0.3, size = 2.8, color = "grey30", fontface = "italic") +
  scale_x_continuous(breaks = cutoffs) +
  scale_y_continuous(limits = c(0, 112)) +
  labs(x = "log2FC cutoff (upregulated, padj < 0.1)",
       y = "% of library recovered",
       title = "Library recovery by Integrated upregulated DEGs",
       subtitle = sprintf("Library: %s human orthologs | Integrated DEG counts in italic above",
                           format(n_lib, big.mark = ","))) +
  theme_pub

ggsave(file.path(outdir, "library_recovery_by_lfc.pdf"), p2,
       width = 7, height = 5)

# =============================================
# Fig 3: Library enrichment (UP only)
# =============================================
enrich_data <- rbindlist(lapply(cutoffs, function(lfc_cut) {
  dream_up <- dream[padj < 0.1 & logFC >= lfc_cut, ensembl_id]
  n_dream <- length(dream_up)
  ov <- length(intersect(lib_human_ids, dream_up))
  data.table(
    cutoff = lfc_cut,
    pct_dream = 100 * ov / n_dream
  )
}))

p3 <- ggplot(enrich_data, aes(x = cutoff, y = pct_dream)) +
  geom_line(linewidth = 1, color = "#2166AC") +
  geom_point(size = 2.5, color = "#2166AC") +
  scale_x_continuous(breaks = cutoffs) +
  labs(x = "log2FC cutoff (upregulated)", y = "% of Integrated DEGs in library",
       title = "Library enrichment among Integrated upregulated DEGs",
       subtitle = sprintf("Library: %s human orthologs",
                           format(n_lib, big.mark = ","))) +
  theme_pub

ggsave(file.path(outdir, "library_enrichment_in_dream.pdf"), p3,
       width = 6, height = 4)

# =============================================
# Fig 4: Direction concordance (human ortholog level)
# =============================================
dream_sig <- dream[padj < 0.1]
setkey(dream_sig, ensembl_id)

dir_status <- sapply(lib_human_ids, function(eid) {
  hit <- dream_sig[ensembl_id == eid]
  if (nrow(hit) == 0) return("no_match")
  if (hit$logFC[1] > 0) "up" else "down"
})

dir_counts <- c(
  up = sum(dir_status == "up"),
  down = sum(dir_status == "down"),
  no_match = sum(dir_status == "no_match")
)

dir_data <- data.table(
  direction = factor(c("Concordant\n(Up in both)", "Discordant\n(Down in Integrated)", "No sig. match\nin Integrated"),
                     levels = c("Concordant\n(Up in both)", "Discordant\n(Down in Integrated)", "No sig. match\nin Integrated")),
  n = c(dir_counts["up"], dir_counts["down"], dir_counts["no_match"])
)
dir_data[, pct := round(100 * n / n_lib, 1)]

p4 <- ggplot(dir_data, aes(x = direction, y = n, fill = direction)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Concordant\n(Up in both)" = "#D6604D",
                                "Discordant\n(Down in Integrated)" = "#4393C3",
                                "No sig. match\nin Integrated" = "#BABABA"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Number of genes",
       title = "Direction concordance of library genes in Integrated",
       subtitle = sprintf("Library: %s human orthologs | padj < 0.1 in Integrated",
                           format(n_lib, big.mark = ","))) +
  theme_pub +
  guides(fill = "none")

ggsave(file.path(outdir, "direction_concordance.pdf"), p4,
       width = 6, height = 5)

# =============================================
# Fig 5: Status breakdown (human ortholog level)
# =============================================
dream_up_hi <- dream[padj < 0.1 & logFC >= 0.5, ensembl_id]
dream_up_lo <- dream[padj < 0.1 & logFC >= 0 & logFC < 0.5, ensembl_id]

status <- sapply(lib_human_ids, function(eid) {
  if (eid %in% dream_up_hi) return("Up & LFC >= 0.5")
  if (eid %in% dream_up_lo) return("Up & LFC < 0.5")
  if (eid %in% dream$ensembl_id) return("Down or NS")
  return("Not in Integrated")
})

lvls <- c("Up & LFC >= 0.5", "Up & LFC < 0.5", "Down or NS", "Not in Integrated")
counts <- table(factor(status, levels = lvls))
status_data <- data.table(
  status = factor(names(counts), levels = lvls),
  n = as.integer(counts)
)
status_data[, pct := round(100 * n / sum(n), 1)]

p5 <- ggplot(status_data, aes(x = status, y = n, fill = status)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(
    "Up & LFC >= 0.5" = "#D6604D",
    "Up & LFC < 0.5" = "#F4A582",
    "Down or NS" = "#92C5DE",
    "Not in Integrated" = "#BABABA"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Number of genes",
       title = "Status of library genes in Integrated analysis",
       subtitle = sprintf("Library: %s human orthologs",
                           format(n_lib, big.mark = ","))) +
  theme_pub +
  theme(axis.text.x = element_text(size = 9)) +
  guides(fill = "none")

ggsave(file.path(outdir, "library_gene_status.pdf"), p5,
       width = 7, height = 5)

cat("All figures saved to", outdir, "\n")
