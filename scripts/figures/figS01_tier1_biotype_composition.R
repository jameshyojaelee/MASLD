#!/usr/bin/env Rscript
# KEY MESSAGE: canonical Tier 1 DEGs (TREAT FDR<0.05 at lfc=0.25; n~1,918,
# recomputed live) are dominated by protein-coding and lncRNA; the remainder are
# technical-noise-prone biotypes
# (processed pseudogenes, immune V(D)J recombination genes, small RNAs).
#
# Panel for figS01_qc_validation — biotype composition of canonical Tier 1 DEGs.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source(file.path(Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "scripts/figures/load_figure_data.R"))

dr <- load_dream_results()
md <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))

t1 <- dr[is_dream_deg(dr)]
t1[, ensembl_base := tstrsplit(gene, ".", fixed = TRUE, keep = 1L)]
m  <- merge(t1,
            md[, .(ensembl_base, gene_biotype, gene_name, chromosome)],
            by = "ensembl_base", all.x = TRUE)

# Roll-up categories (mutually exclusive, priority order matters)
classify <- function(symbol, biotype, chr) {
  symbol  <- ifelse(is.na(symbol),  "",            symbol)
  biotype <- ifelse(is.na(biotype), "unannotated", biotype)
  is_tcr  <- grepl("^TR[ABGD][VJDC]", symbol)
  is_ig   <- grepl("^IG[HKL][VJDC]",  symbol) | biotype %in% c("IG_V_gene","IG_J_gene","IG_C_gene","IG_D_gene")
  is_hla  <- grepl("^HLA-", symbol)
  is_mt   <- chr %in% "chrM" | biotype %in% c("Mt_tRNA","Mt_rRNA")
  is_pseudo_processed   <- biotype %in% c("processed_pseudogene","unprocessed_pseudogene",
                                          "rRNA_pseudogene","transcribed_unitary_pseudogene")
  is_pseudo_transcribed <- biotype %in% c("transcribed_processed_pseudogene",
                                          "transcribed_unprocessed_pseudogene")
  is_pc   <- biotype == "protein_coding"
  is_lnc  <- biotype == "lncRNA"
  is_mir  <- biotype == "miRNA"
  is_small <- biotype %in% c("snoRNA","snRNA","scaRNA","misc_RNA","sRNA","scRNA",
                             "vault_RNA","Y_RNA","rRNA","ribozyme")
  is_tec  <- biotype == "TEC"
  ifelse(is_tcr,                "TCR V(D)J",
  ifelse(is_ig,                 "Immunoglobulin V(D)J",
  ifelse(is_hla,                "HLA (alt-contig)",
  ifelse(is_mt,                 "Mitochondrial",
  ifelse(is_pseudo_transcribed, "Transcribed pseudogene",
  ifelse(is_pseudo_processed,   "Processed pseudogene",
  ifelse(is_mir,                "miRNA (primary transcript)",
  ifelse(is_small,              "Other small RNA",
  ifelse(is_tec,                "TEC (unclassified)",
  ifelse(is_lnc,                "lncRNA",
  ifelse(is_pc,                 "protein-coding",
                                paste0("Other (", biotype, ")"))))))))))))
}

m[, category   := classify(symbol, gene_biotype, chromosome)]
m[, direction  := ifelse(bulk_logFC > 0, "Up", "Down")]

tally <- m[, .(N = .N), by = .(category, direction)]
order_tally <- m[, .(N = .N), by = category][order(N)]
tally[, category := factor(category, levels = order_tally$category)]
totals <- m[, .(N = .N), by = category][order(-N)]
totals[, pct := round(100 * N / sum(N), 1)]

cat("=== Tier 1 biotype composition (n=", nrow(m), ") ===\n", sep = "")
print(totals)

# --- Plot: horizontal bar, sorted by N, split Up/Down ---
clean_pct <- merge(order_tally, totals[, .(category, pct)], by = "category")
clean_pct[, label := sprintf("%d  (%.1f%%)", N, pct)]

p <- ggplot(tally, aes(x = N, y = category, fill = direction)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  geom_text(data = clean_pct,
            aes(x = N, y = category, label = label),
            inherit.aes = FALSE,
            hjust = -0.05, size = 6 / ggplot2::.pt) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = c(Up = "#D6604D", Down = "#4393C3"),
                    breaks = c("Up", "Down")) +
  labs(x = "Number of Tier 1 DEGs (TREAT FDR<0.05, lfc=0.25)",
       y = NULL,
       fill = "Direction") +
  theme(legend.position = "top")

message(sprintf(
  "[caption] Biotype composition of canonical Tier 1 DEGs. Total n = %s (Up = %s, Down = %s); annotation: GENCODE v49",
  format(nrow(m), big.mark = ","),
  format(sum(m$direction == "Up"),   big.mark = ","),
  format(sum(m$direction == "Down"), big.mark = ",")))

out_pdf <- file.path(FIGS01_DIR, "figS01_tier1_biotype_composition.pdf")
out_csv <- file.path(FIGS01_DIR, "figS01_tier1_biotype_composition.csv")
ggsave(out_pdf, p, width = 7.0, height = 4.5, useDingbats = FALSE)
fwrite(totals, out_csv)
cat("Wrote:\n  ", out_pdf, "\n  ", out_csv, "\n", sep = "")
