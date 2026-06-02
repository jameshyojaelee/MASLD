#!/usr/bin/env Rscript
#SBATCH --job-name=net_283_dcoloc
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=logs/net_283_dcoloc_%j.out
#SBATCH --error=logs/net_283_dcoloc_%j.err
#
# Script 283: Build D-COLOC edges -- gene pairs where both endpoints share
#             a SuSiE-COLOC locus (PP.H4 > 0.5) either via the same top SNP
#             (strict) or within 500 kb on the same chromosome (relaxed).

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Paths ----------------------------------------------------------------
PROJ    <- Sys.getenv("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
COLOC   <- file.path(PROJ, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv")
ATLAS   <- file.path(PROJ, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
NODES   <- file.path(PROJ, "RNA-seq/results/network/network_nodes.csv")
OUT     <- file.path(PROJ, "RNA-seq/results/network/edges_d_coloc.csv")

PP4_MIN    <- 0.5
WINDOW_BP  <- 500000L

# ---- Load -----------------------------------------------------------------
cat("[283] Loading SuSiE-COLOC, atlas, and node set...\n")
coloc <- fread(COLOC)
nodes <- fread(NODES)
atlas <- fread(ATLAS, select = c("human_symbol", "ensembl_id"))

V <- unique(nodes[[1]])
cat(sprintf("[283] Node set V: %d genes\n", length(V)))

# Expected columns: gene, ensembl, chr, gwas_name, PP.H4.abf, top_snp
stopifnot(all(c("gene", "chr", "gwas_name", "PP.H4.abf", "top_snp") %in% colnames(coloc)))

coloc <- coloc[!is.na(PP.H4.abf) & PP.H4.abf > PP4_MIN]
cat(sprintf("[283] %d gene-GWAS records with PP.H4 > %g (unique genes: %d)\n",
            nrow(coloc), PP4_MIN, length(unique(coloc$gene))))

# Parse top_snp position when of the form "chr:pos"
coloc[, snp_chr := tstrsplit(top_snp, ":", fixed = TRUE, keep = 1L)]
coloc[, snp_pos := suppressWarnings(as.integer(
        tstrsplit(top_snp, ":", fixed = TRUE, keep = 2L)[[1]]))]

# ---- Strict pairs: same top_snp within a GWAS -----------------------------
cat("[283] Building STRICT pairs (same top_snp)...\n")
strict_pairs <- coloc[, {
  if (.N >= 2) {
    g <- gene
    p <- PP.H4.abf
    idx <- combn(.N, 2)
    ga <- g[idx[1, ]]; gb <- g[idx[2, ]]
    pa <- p[idx[1, ]]; pb <- p[idx[2, ]]
    # canonical ordering
    swap <- ga > gb
    tmp <- ga[swap]; ga[swap] <- gb[swap]; gb[swap] <- tmp
    tmp <- pa[swap]; pa[swap] <- pb[swap]; pb[swap] <- tmp
    list(gene_a = ga, gene_b = gb,
         pp4_a = pa, pp4_b = pb,
         pp4_min = pmin(pa, pb),
         top_snp_a = top_snp[1], top_snp_b = top_snp[1],
         locus_id = paste0(gwas_name[1], "|", top_snp[1]),
         distance_bp = 0L,
         mode = "strict",
         gwas_name = gwas_name[1])
  } else NULL
}, by = .(gwas_name, top_snp)]

cat(sprintf("[283] Strict pairs: %d\n", nrow(strict_pairs)))

# ---- Relaxed pairs: same chromosome within WINDOW_BP, same GWAS -----------
cat("[283] Building RELAXED pairs (within 500 kb, same chr, same GWAS)...\n")
coloc_pos <- coloc[!is.na(snp_pos) & !is.na(chr)]
# Work per GWAS x chr
relaxed_list <- coloc_pos[, {
  if (.N < 2) { NULL } else {
    ord <- order(snp_pos)
    g <- gene[ord]; p <- PP.H4.abf[ord]; s <- top_snp[ord]; po <- snp_pos[ord]
    n <- length(g)
    out_a <- character(); out_b <- character()
    out_pa <- numeric(); out_pb <- numeric()
    out_sa <- character(); out_sb <- character()
    out_d <- integer()
    for (i in seq_len(n - 1)) {
      j <- i + 1L
      while (j <= n && (po[j] - po[i]) <= WINDOW_BP) {
        if (g[i] != g[j]) {
          # canonical order
          if (g[i] <= g[j]) {
            out_a <- c(out_a, g[i]); out_b <- c(out_b, g[j])
            out_pa <- c(out_pa, p[i]); out_pb <- c(out_pb, p[j])
            out_sa <- c(out_sa, s[i]); out_sb <- c(out_sb, s[j])
          } else {
            out_a <- c(out_a, g[j]); out_b <- c(out_b, g[i])
            out_pa <- c(out_pa, p[j]); out_pb <- c(out_pb, p[i])
            out_sa <- c(out_sa, s[j]); out_sb <- c(out_sb, s[i])
          }
          out_d <- c(out_d, po[j] - po[i])
        }
        j <- j + 1L
      }
    }
    if (length(out_a) == 0) NULL else
    list(gene_a = out_a, gene_b = out_b,
         pp4_a = out_pa, pp4_b = out_pb,
         pp4_min = pmin(out_pa, out_pb),
         top_snp_a = out_sa, top_snp_b = out_sb,
         locus_id = paste0(gwas_name[1], "|", chr[1], ":",
                           pmin(out_d, WINDOW_BP)),
         distance_bp = out_d,
         mode = "relaxed",
         gwas_name = gwas_name[1])
  }
}, by = .(gwas_name, chr)]

cat(sprintf("[283] Relaxed pairs: %d\n", nrow(relaxed_list)))

# ---- Combine & deduplicate ------------------------------------------------
all_pairs <- rbind(
  strict_pairs[, .(gene_a, gene_b, pp4_a, pp4_b, pp4_min,
                   top_snp_a, top_snp_b, locus_id, distance_bp, mode, gwas_name)],
  relaxed_list[, .(gene_a, gene_b, pp4_a, pp4_b, pp4_min,
                   top_snp_a, top_snp_b, locus_id, distance_bp, mode, gwas_name)]
)

cat(sprintf("[283] Combined pairs: %d\n", nrow(all_pairs)))

# Filter: both endpoints in V
all_pairs <- all_pairs[gene_a %in% V & gene_b %in% V]
cat(sprintf("[283] After V filter: %d\n", nrow(all_pairs)))

# Deduplicate: for each (gene_a, gene_b), keep the record with highest pp4_min
# and aggregate GWAS / locus / mode info.
setorder(all_pairs, gene_a, gene_b, -pp4_min)
dedup <- all_pairs[, {
  keep <- which.max(pp4_min)
  list(
    pp4_min    = pp4_min[keep],
    pp4_a      = pp4_a[keep],
    pp4_b      = pp4_b[keep],
    top_snp_a  = top_snp_a[keep],
    top_snp_b  = top_snp_b[keep],
    gwas_list  = paste(unique(gwas_name), collapse = ";"),
    locus_id   = locus_id[keep],
    mode       = if (any(mode == "strict")) "strict" else "relaxed",
    n_shared_loci = length(unique(locus_id))
  )
}, by = .(gene_a, gene_b)]

cat(sprintf("[283] Final D-COLOC edges: %d\n", nrow(dedup)))

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(dedup, OUT)
cat(sprintf("[283] Wrote %s\n", OUT))
cat(sprintf("  strict-dominant edges: %d\n", sum(dedup$mode == "strict")))
cat(sprintf("  relaxed-only edges:    %d\n", sum(dedup$mode == "relaxed")))
cat(sprintf("  median pp4_min:        %.3f\n", median(dedup$pp4_min)))
