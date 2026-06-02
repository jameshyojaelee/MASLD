#!/usr/bin/env Rscript
# 342b2_katsuda_finalize.R - take Katsuda DE results + probe annotation + rat-human orthologs,
# produce a rat-derived polyploid signature in human gene symbols.

suppressPackageStartupMessages({library(data.table)})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
setwd(ROOT)
OUT_DIR <- "data/ploidy_signatures"

# Load DE results (probe-level, both contrasts)
de <- fread(file.path(OUT_DIR, "katsuda2019_de_full.tsv"))
cat("DE rows:", nrow(de), "  contrasts:", paste(unique(de$contrast), collapse=", "), "\n")

# Probe -> rat gene_symbol
p2s <- fread("data/external/katsuda2019/GPL22740_probe2symbol.tsv")
# Rat -> human ortholog map
orth <- fread("data/external/orthologs/rat_human_orthologs.tsv.gz")
one2one <- orth[orthology_type == "ortholog_one2one"][!is.na(human_external_gene_name) & human_external_gene_name != ""]
# A few rat symbols have duplicate orthologs (e.g., E2f1 has E2F1 + NECAB3); keep first
one2one <- unique(one2one, by = "rat_external_gene_name")
r2h <- one2one[, .(rat_symbol = rat_external_gene_name, human_symbol = human_external_gene_name)]
cat("rat->human one2one orthologs:", nrow(r2h), "\n")

# Annotate DE table
de <- merge(de, p2s, by = "PROBE_ID", by.x = "PROBE_ID", by.y = "PROBE_ID", all.x = TRUE)
de <- merge(de, r2h, by.x = "GENE_SYMBOL", by.y = "rat_symbol", all.x = TRUE)
cat("DE rows with human ortholog:", sum(!is.na(de$human_symbol)), "of", nrow(de), "\n")
fwrite(de, file.path(OUT_DIR, "katsuda2019_de_full.tsv"), sep = "\t")  # overwrite with annotations

# Use the binary contrast (polyploid 4c+8c vs 2c)
de_bin <- de[contrast == "polyploid_vs_2c" & !is.na(human_symbol)]
# When multiple probes map to same human gene, take the one with strongest |t|
setorderv(de_bin, c("human_symbol", "t"), order = c(1L, -1L))
# For each gene take row with max |t|
de_bin[, abs_t := abs(t)]
setorder(de_bin, human_symbol, -abs_t)
de_bin_unique <- de_bin[, .SD[1], by = human_symbol]
cat("genes after per-gene reduction:", nrow(de_bin_unique), "\n")

# Top N by Wald-like t-statistic (signed), separately up and down
N_SIG <- 100
sig_up <- de_bin_unique[t > 0][order(-t)][1:N_SIG]
sig_dn <- de_bin_unique[t < 0][order(t)][1:N_SIG]
sig_up <- sig_up[!is.na(human_symbol)]
sig_dn <- sig_dn[!is.na(human_symbol)]

cat("\nKatsuda signature:\n")
cat("  up (4c+8c > 2c):", nrow(sig_up), "  t range:", round(range(sig_up$t),2), "\n")
cat("  down (2c > 4c+8c):", nrow(sig_dn), "  t range:", round(range(sig_dn$t),2), "\n")
cat("Top 10 up:", paste(head(sig_up$human_symbol, 10), collapse=", "), "\n")
cat("Top 10 down:", paste(head(sig_dn$human_symbol, 10), collapse=", "), "\n")

fwrite(sig_up[, .(rat_symbol = GENE_SYMBOL, human_symbol, PROBE_ID, logFC, t, P.Value, adj.P.Val)],
       file.path(OUT_DIR, "katsuda2019_polyploid_up.tsv"), sep = "\t")
fwrite(sig_dn[, .(rat_symbol = GENE_SYMBOL, human_symbol, PROBE_ID, logFC, t, P.Value, adj.P.Val)],
       file.path(OUT_DIR, "katsuda2019_polyploid_down.tsv"), sep = "\t")

# Mechanism marker overlap
mm <- fread(file.path(OUT_DIR, "mechanism_markers.tsv"))
mm_human <- toupper(mm$gene)
ovl_up <- intersect(mm_human, toupper(sig_up$human_symbol))
ovl_dn <- intersect(mm_human, toupper(sig_dn$human_symbol))
cat("\nMechanism marker overlap:\n")
cat("  up:  ", paste(ovl_up, collapse=", "), "\n")
cat("  down:", paste(ovl_dn, collapse=", "), "\n")
cat("  total:", length(ovl_up)+length(ovl_dn), "of", length(mm_human), "markers\n")
