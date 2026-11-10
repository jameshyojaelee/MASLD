# Fig 5B check: protein vs bulk mRNA change on the prioritized universe,
# using my independent 774-unit refit and a unique-symbol GENCODE v49 map.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]; MYDEG <- args[2]
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
PROT <- file.path(R, "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv")
UNI <- file.path(R, "Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt")
GM <- file.path(R, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
PROD <- file.path(R, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
logf <- file(file.path(OUT, "fig5b_log.txt"), "w")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, logf); flush(logf) }

prot <- fread(PROT)
say("protein rows: ", nrow(prot), "; padj<0.05: ", sum(prot$padj < 0.05, na.rm = TRUE))
uni <- trimws(readLines(UNI)); uni <- uni[nzchar(uni)]
gm <- fread(GM)
run <- function(deg, tag) {
  deg[, base := sub("\\..*", "", gene)]
  deg[, sym := gm$gene_name[match(base, gm$ensembl_base)]]
  u <- deg[!is.na(sym) & nzchar(sym)]
  dup <- u[, .N, by = sym][N > 1, sym]
  u <- u[!sym %in% dup]
  m <- merge(prot[, .(sym = gene, p_lfc = logFC, p_padj = padj)], u[, .(sym, r_lfc = logFC, r_padj = padj)], by = "sym")
  m <- m[is.finite(p_lfc) & is.finite(r_lfc)]
  d <- m[sym %in% uni]
  sig <- d[p_padj < 0.05]
  conc <- sum(sign(sig$p_lfc) == sign(sig$r_lfc))
  rs <- sig[r_padj < 0.05]
  say(sprintf("[%s] dup symbols dropped=%d measured all=%d rho_all=%.4f | universe measured=%d rho=%.6f dir=%.1f%% | protein-sig=%d concordant=%d discordant=%d | RNA padj<0.05 among sig: %d of which concordant %d; RNA not sig: %d concordant %d",
              tag, length(dup), nrow(m), cor(m$p_lfc, m$r_lfc, method = "spearman"), nrow(d),
              cor(d$p_lfc, d$r_lfc, method = "spearman"), 100 * mean(sign(d$p_lfc) == sign(d$r_lfc)),
              nrow(sig), conc, nrow(sig) - conc, nrow(rs), sum(sign(rs$p_lfc) == sign(rs$r_lfc)),
              nrow(sig) - nrow(rs), conc - sum(sign(rs$p_lfc) == sign(rs$r_lfc))))
  d
}
mine <- fread(MYDEG)
d1 <- run(mine, "my refit")
prod <- fread(PROD)
d2 <- run(prod[, .(gene, logFC, padj)], "producer deg_results")
say("universe genes in protein table but not matched (my map): ", length(setdiff(intersect(prot$gene, uni), d1$sym)))
close(logf)
