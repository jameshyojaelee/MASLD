#!/usr/bin/env Rscript
# ============================================================================
# Gap #10: Systematic mRNA-protein concordance analysis
# ============================================================================
# Computes concordance between dream mega-analysis mRNA logFC and protein logFC
# across 2 proteomics datasets (PXD051911, PXD052937).
# NOTE: GSE276114 excluded — confirmed bulk RNA-seq, not proteomics.
#
# Outputs:
#   figures/supplementary/figS05_epigenomic_spatial/figS_mrna_protein_concordance.pdf
#   Analysis/Proteomics/results/mrna_protein_concordance_summary.csv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

cat("=== Gap #10: mRNA-protein concordance ===\n")

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
conc  <- fread(file.path(BASE, "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv"))

cat(sprintf("Atlas: %d genes, %d with protein data (n_prot_datasets > 0)\n",
            nrow(atlas), sum(!is.na(atlas$n_prot_datasets) & atlas$n_prot_datasets > 0)))
cat(sprintf("Concordance file: %d rows across %d unique genes\n",
            nrow(conc), length(unique(conc$gene))))

# ---------------------------------------------------------------------------
# 2. Focus on disease_vs_control comparator — PXD052937 (plasma DIA-MS) only
#
#    REWIRED 2026-04-22 (T0.6 modality correction): previously mixed
#    GSE276114 (bulk RNA-seq, wrongly flagged as proteomics) with
#    PXD051911 / PXD052937 via `!duplicated(paste0(gene, dream_comparator))`,
#    which silently let the RNA-seq row dominate the mRNA-vs-protein
#    concordance because it had the smallest p-values. That was circular
#    (RNA vs RNA). We now explicitly select PXD052937 (plasma DIA-MS;
#    Sourianarayanane et al.) as the primary proteomics comparator.
#
#    PXD051911 (liver DIA-MS; Boel et al.) is kept as a tissue-proteomics
#    secondary; GSE276114 fibrosis is excluded here because it is RNA-seq.
# ---------------------------------------------------------------------------
# Primary: PXD052937 plasma DIA-MS, MASLD vs Normal
dvc <- copy(conc[dataset == "PXD052937" &
                 dream_comparator == "disease_vs_control"])
# Secondary tissue comparator: PXD051911 liver DIA-MS, MASLD vs No_MASLD
tissue <- copy(conc[dataset == "PXD051911" &
                    dream_comparator == "disease_vs_control"])
# Fibrosis comparator: PXD051911 MASH-vs-MASL (real liver proteomics)
# — replaces the prior GSE276114 comparison which was RNA-vs-RNA.
fib <- copy(conc[dataset == "PXD051911_mash_vs_masl"])

# Map UniProt -> HGNC symbol for PXD052937 rows
candidates_f <- file.path(BASE,
  "data/PXD052937/20220805_163627_Plasma_liver2/Candidates.tsv")
uniprot_map <- data.table()
if (file.exists(candidates_f)) {
  cand <- fread(candidates_f, sep = "\t")
  map_raw <- unique(cand[, .(ProteinGroups, Genes)])
  map_raw[, protein_id := sub(";.*", "", ProteinGroups)]
  map_raw[, gene_name  := sub(";.*", "", Genes)]
  uniprot_map <- unique(
    map_raw[gene_name != "" & !is.na(gene_name),
            .(protein_id, gene_name)],
    by = "protein_id")
  cat(sprintf("Loaded UniProt->HGNC map for PXD052937: %d entries\n",
              nrow(uniprot_map)))
}
if (nrow(dvc) > 0 && nrow(uniprot_map) > 0) {
  setnames(dvc, "gene", "protein_id")
  dvc <- merge(dvc, uniprot_map, by = "protein_id", all.x = TRUE)
  dvc[, gene := ifelse(!is.na(gene_name) & gene_name != "",
                        gene_name, protein_id)]
  dvc[, gene_name := NULL]
  setorder(dvc, gene, protein_padj)
  dvc <- dvc[!duplicated(gene)]
}

cat(sprintf("\n[Primary]   PXD052937 plasma DIA-MS (MASLD vs Normal): %d genes with both mRNA + protein\n", nrow(dvc)))
cat(sprintf("[Secondary] PXD051911 liver DIA-MS (MASLD vs No_MASLD): %d genes\n", nrow(tissue)))
cat(sprintf("[Fibrosis]  PXD051911 liver (MASH vs MASL): %d genes\n", nrow(fib)))

# Merge atlas annotations
dvc <- merge(dvc, atlas[, .(human_symbol, is_conserved, bulk_padj,
                             n_leading_edge_pathways, top_pathways,
                             human_consensus_tier, primary_category,
                             n_prot_datasets)],
             by.x = "gene", by.y = "human_symbol", all.x = TRUE,
             suffixes = c("", ".atlas"))

# ---------------------------------------------------------------------------
# 3. Compute concordance statistics
# ---------------------------------------------------------------------------
# DEG status (dream padj < 0.05, |logFC| > 0.3)
dvc[, is_deg := !is.na(dream_padj) & dream_padj < 0.05 & abs(dream_logFC) > 0.3]
dvc[, is_prot_sig := !is.na(protein_padj) & protein_padj < 0.05]
dvc[, both_sig := is_deg & is_prot_sig]
dvc[, direction_match := sign(dream_logFC) == sign(protein_logFC)]
dvc[, is_conserved := ifelse(is.na(is_conserved), FALSE, is_conserved)]

# Overall concordance
overall_rho <- cor(dvc$dream_logFC, dvc$protein_logFC, method = "spearman",
                   use = "complete.obs")
overall_dir <- mean(dvc$direction_match, na.rm = TRUE) * 100

cat(sprintf("\n--- Disease vs Control concordance ---\n"))
cat(sprintf("Overall Spearman rho: %.3f\n", overall_rho))
cat(sprintf("Overall direction concordance: %.1f%%\n", overall_dir))

# Stratified
strat <- rbind(
  data.table(stratum = "All genes", n = nrow(dvc),
             rho = overall_rho,
             direction_pct = overall_dir),
  data.table(stratum = "DEGs (padj<0.05, |LFC|>0.3)", n = sum(dvc$is_deg),
             rho = cor(dvc[is_deg == TRUE]$dream_logFC,
                       dvc[is_deg == TRUE]$protein_logFC,
                       method = "spearman", use = "complete.obs"),
             direction_pct = mean(dvc[is_deg == TRUE]$direction_match, na.rm = TRUE) * 100),
  data.table(stratum = "Non-DEGs", n = sum(!dvc$is_deg),
             rho = cor(dvc[is_deg == FALSE]$dream_logFC,
                       dvc[is_deg == FALSE]$protein_logFC,
                       method = "spearman", use = "complete.obs"),
             direction_pct = mean(dvc[is_deg == FALSE]$direction_match, na.rm = TRUE) * 100),
  data.table(stratum = "Conserved", n = sum(dvc$is_conserved),
             rho = if (sum(dvc$is_conserved) > 10)
               cor(dvc[is_conserved == TRUE]$dream_logFC,
                   dvc[is_conserved == TRUE]$protein_logFC,
                   method = "spearman", use = "complete.obs") else NA_real_,
             direction_pct = if (sum(dvc$is_conserved) > 10)
               mean(dvc[is_conserved == TRUE]$direction_match, na.rm = TRUE) * 100 else NA_real_),
  data.table(stratum = "Both significant", n = sum(dvc$both_sig),
             rho = if (sum(dvc$both_sig) > 10)
               cor(dvc[both_sig == TRUE]$dream_logFC,
                   dvc[both_sig == TRUE]$protein_logFC,
                   method = "spearman", use = "complete.obs") else NA_real_,
             direction_pct = if (sum(dvc$both_sig) > 10)
               mean(dvc[both_sig == TRUE]$direction_match, na.rm = TRUE) * 100 else NA_real_),
  data.table(stratum = "Protein sig only", n = sum(dvc$is_prot_sig & !dvc$is_deg),
             rho = if (sum(dvc$is_prot_sig & !dvc$is_deg) > 10)
               cor(dvc[is_prot_sig == TRUE & is_deg == FALSE]$dream_logFC,
                   dvc[is_prot_sig == TRUE & is_deg == FALSE]$protein_logFC,
                   method = "spearman", use = "complete.obs") else NA_real_,
             direction_pct = if (sum(dvc$is_prot_sig & !dvc$is_deg) > 10)
               mean(dvc[is_prot_sig == TRUE & is_deg == FALSE]$direction_match, na.rm = TRUE) * 100 else NA_real_)
)

# Secondary: PXD051911 liver DIA-MS (MASLD vs No_MASLD) — tissue-proteomics
if (nrow(tissue) > 50) {
  t_rho <- cor(tissue$dream_logFC, tissue$protein_logFC, method = "spearman",
               use = "complete.obs")
  t_dir <- mean(sign(tissue$dream_logFC) == sign(tissue$protein_logFC),
                na.rm = TRUE) * 100
  strat <- rbind(strat,
    data.table(stratum = "PXD051911 liver (MASLD vs ctrl)", n = nrow(tissue),
               rho = t_rho, direction_pct = t_dir))
  cat(sprintf("PXD051911 liver Spearman rho: %.3f, direction: %.1f%%\n",
              t_rho, t_dir))
}

# Fibrosis: PXD051911 MASH vs MASL
if (nrow(fib) > 50) {
  fib_rho <- cor(fib$dream_logFC, fib$protein_logFC, method = "spearman",
                 use = "complete.obs")
  fib_dir <- mean(sign(fib$dream_logFC) == sign(fib$protein_logFC), na.rm = TRUE) * 100
  strat <- rbind(strat,
    data.table(stratum = "PXD051911 (MASH vs MASL)", n = nrow(fib),
               rho = fib_rho, direction_pct = fib_dir))
  cat(sprintf("PXD051911 MASH-vs-MASL Spearman rho: %.3f, direction: %.1f%%\n",
              fib_rho, fib_dir))
}

# Pathway-stratified concordance
if ("top_pathways" %in% names(dvc)) {
  # Extract first pathway for each gene
  dvc[, first_pathway := sub(";.*", "", top_pathways)]
  dvc[first_pathway == "" | is.na(first_pathway), first_pathway := "None"]

  pathway_strat <- dvc[, .(
    n = .N,
    rho = if (.N > 20) cor(dream_logFC, protein_logFC, method = "spearman",
                           use = "complete.obs") else NA_real_,
    direction_pct = mean(direction_match, na.rm = TRUE) * 100
  ), by = first_pathway][order(-n)]

  # Keep top pathways with enough genes
  pathway_strat <- pathway_strat[n >= 20 & !is.na(rho)]
  cat(sprintf("\nPathway-stratified concordance (%d pathways with n>=20):\n",
              nrow(pathway_strat)))
  print(head(pathway_strat, 10))
}

cat("\n--- Stratified summary ---\n")
print(strat)

# ---------------------------------------------------------------------------
# 4. Identify discordant genes
# ---------------------------------------------------------------------------
# Discordant = opposite sign AND both have non-trivial effect
dvc[, discordance_score := abs(dream_logFC) + abs(protein_logFC)]
discordant <- dvc[direction_match == FALSE & abs(dream_logFC) > 0.2 & abs(protein_logFC) > 0.2]
discordant <- discordant[order(-discordance_score)]

cat(sprintf("\nDiscordant genes (opp. sign, both |LFC|>0.2): %d\n", nrow(discordant)))
cat("Top 20 discordant:\n")
print(discordant[1:min(20, nrow(discordant)),
                 .(gene, dream_logFC, protein_logFC, dream_padj, protein_padj,
                   is_deg, is_conserved)])

# ---------------------------------------------------------------------------
# 5. Figures
# ---------------------------------------------------------------------------
outdir <- file.path(BASE, "figures/supplementary/figS05_epigenomic_spatial")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "panels"), recursive = TRUE, showWarnings = FALSE)

# Significance category for coloring
dvc[, sig_cat := fifelse(
  both_sig, "Both significant",
  fifelse(is_deg & !is_prot_sig, "mRNA DEG only",
  fifelse(!is_deg & is_prot_sig, "Protein sig only",
  "Neither significant")))]

sig_cat_colors <- c(
  "Both significant"  = masld_colors$up,
  "mRNA DEG only"     = "#42A5F5",
  "Protein sig only"  = "#F48FB1",
  "Neither significant" = masld_colors$ns
)

# --- Panel (a): mRNA vs protein logFC scatter ---
# Label top concordant and discordant genes
label_genes <- c(
  head(dvc[both_sig == TRUE][order(-abs(dream_logFC))]$gene, 8),
  head(discordant$gene, 7)
)
dvc[, label := fifelse(gene %in% label_genes, gene, "")]

pa <- ggplot(dvc, aes(x = dream_logFC, y = protein_logFC, color = sig_cat)) +
  rasterize_layer(geom_point(alpha = 0.4, size = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  geom_smooth(data = dvc[both_sig == TRUE | (is_deg & abs(protein_logFC) > 0.1)],
              method = "lm", se = TRUE, color = "black", linewidth = 0.4, alpha = 0.15,
              inherit.aes = FALSE,
              aes(x = dream_logFC, y = protein_logFC)) +
  geom_text_repel(aes(label = label), size = 1.8, max.overlaps = 20,
                  segment.size = 0.2, show.legend = FALSE) +
  scale_color_manual(values = sig_cat_colors, name = "Significance") +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("rho = %.3f\nn = %s", overall_rho, format(nrow(dvc), big.mark = ",")),
           hjust = 1.1, vjust = -0.5, size = 2.2, fontface = "italic") +
  labs(x = "mRNA logFC (dream mega-analysis)",
       y = "Protein logFC",
       title = "mRNA vs protein fold-change concordance") +
  theme_masld() +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA))

# --- Panel (b): Direction concordance bar chart ---
bar_data <- strat[!is.na(direction_pct) & n >= 10]
bar_data[, stratum := factor(stratum, levels = rev(bar_data$stratum))]

pb <- ggplot(bar_data, aes(x = stratum, y = direction_pct, fill = direction_pct)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%\n(n=%s)", direction_pct, format(n, big.mark = ","))),
            hjust = -0.05, size = 1.8) +
  scale_fill_gradient(low = "#42A5F5", high = masld_colors$up, guide = "none") +
  coord_flip(ylim = c(0, 100)) +
  labs(x = NULL, y = "Direction concordance (%)",
       title = "mRNA-protein direction agreement") +
  theme_masld()

# --- Panel (c): Top 20 discordant genes ---
disc_top <- discordant[!duplicated(gene)][1:min(20, .N)]
disc_top[, gene := factor(gene, levels = rev(unique(gene)))]

# Melt for paired bar chart
disc_melt <- melt(disc_top, id.vars = "gene",
                  measure.vars = c("dream_logFC", "protein_logFC"),
                  variable.name = "modality", value.name = "logFC")
disc_melt[, modality := fifelse(modality == "dream_logFC", "mRNA", "Protein")]

pc <- ggplot(disc_melt, aes(x = gene, y = logFC, fill = modality)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("mRNA" = masld_colors$down, "Protein" = masld_colors$up),
                    name = "Modality") +
  coord_flip() +
  labs(x = NULL, y = "logFC", title = "Top discordant genes (opposite mRNA vs protein)") +
  theme_masld() +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0))

# --- Compose ---
fig <- pa + pb + pc +
  plot_layout(ncol = 3, widths = c(1.2, 1, 1)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(size = 9, face = "bold")))

save_fig(fig, file.path(outdir, "figS_mrna_protein_concordance.pdf"),
         width = fig_full_width, height = 4.5)
cat(sprintf("\nFigure saved: %s\n",
            file.path(outdir, "figS_mrna_protein_concordance.pdf")))

# ---------------------------------------------------------------------------
# 6. Companion CSV
# ---------------------------------------------------------------------------
# Primary output — PXD052937-only (so downstream never re-mixes RNA-seq)
out_csv <- file.path(BASE, "Analysis/Proteomics/results/mrna_protein_concordance_summary.csv")

summary_out <- dvc[, .(gene, dream_logFC, dream_padj, protein_logFC, protein_padj,
                       dataset, dream_comparator,
                       direction_match, is_deg, is_prot_sig, both_sig,
                       is_conserved, sig_cat, discordance_score)]
fwrite(summary_out, out_csv)
cat(sprintf("Summary CSV saved (PXD052937 primary): %s (%d rows)\n",
            out_csv, nrow(summary_out)))

# Dedicated PXD052937 copy (stable path consumed by Fig 4 caption text)
fwrite(summary_out,
       file.path(BASE, "Analysis/Proteomics/results/pxd052937_mrna_protein_concordance.csv"))

# Also save the stratified table
strat_csv <- file.path(BASE, "Analysis/Proteomics/results/mrna_protein_concordance_stratified.csv")
fwrite(strat, strat_csv)
cat(sprintf("Stratified CSV saved: %s\n", strat_csv))

# Save discordant genes
disc_csv <- file.path(BASE, "Analysis/Proteomics/results/mrna_protein_discordant_genes.csv")
fwrite(discordant[, .(gene, dream_logFC, dream_padj, protein_logFC, protein_padj,
                       dataset, is_deg, is_conserved, discordance_score)],
       disc_csv)
cat(sprintf("Discordant genes CSV saved: %s (%d genes)\n", disc_csv, nrow(discordant)))

cat("\n=== Gap #10 complete ===\n")
