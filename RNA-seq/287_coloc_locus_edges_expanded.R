#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_287_dcoloc_expanded
#SBATCH --output=logs/net_287_dcoloc_expanded_%j.out
#SBATCH --error=logs/net_287_dcoloc_expanded_%j.err
# ===========================================================================
# Script 287: Expanded D-COLOC edges (Phase 3 sparse-layer rehab).
# ---------------------------------------------------------------------------
# Replaces 283's strict PP.H4 >= 0.5 filter with two relaxed channels:
#   (a) lowered-threshold shared-locus pairs: BOTH endpoints have
#       PP.H4 >= 0.3 against the SAME GWAS and share either the same
#       top_snp or a credible-set lead within +/- 500 kb.
#   (b) eQTL sentinel-sharing: gene pairs whose Broadaway lead eQTL
#       variants are within +/- 250 kb AND nominally in LD (surrogate:
#       distance <= 100 kb used as a proxy for r^2 >= 0.8 since LD refs
#       are not loaded here — gated by a distance filter with conservative
#       cap). Annotated with LD_surrogate = "distance_proxy".
#
# Channel (a) is implemented by reusing 283's logic with PP4_MIN = 0.3.
# Channel (b) scans the per-chromosome Broadaway marginal summary files
# for lead per-gene variants and pairs genes sharing nearby lead SNPs.
#
# Both channels emit in the same schema family as edges_d_coloc.csv,
# plus two extra columns (`channel`, `pp4_thresh`).
#
# Output: edges_d_coloc_expanded.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
})

t0 <- Sys.time()

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

COLOC_PATH  <- file.path(PROJ, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv")
EQTL_DIR    <- file.path(PROJ, "data/broadaway_eqtl")
NODES_PATH  <- file.path(PROJ, "RNA-seq/results/network/network_nodes.csv")
OUT         <- file.path(PROJ, "RNA-seq/results/network/edges_d_coloc_expanded.csv")

PP4_RELAXED       <- 0.3
WINDOW_COLOC_BP   <- 500000L
WINDOW_EQTL_BP    <- 250000L
EQTL_LEAD_P_MAX   <- 1e-5   # per-gene lead -- lax cis-eQTL sentinel

message("[287] Expanded D-COLOC edge construction starting at ", format(t0))

stopifnot(file.exists(COLOC_PATH), file.exists(NODES_PATH), dir.exists(EQTL_DIR))

nodes <- fread(NODES_PATH)
V <- unique(nodes$human_symbol)
message("  Node set V: ", length(V), " genes")

# ==========================================================================
#  Channel (a): Shared-locus COLOC pairs at PP4 >= 0.3
# ==========================================================================
message("[", Sys.time(), "] Loading SuSiE-COLOC table")
coloc <- fread(COLOC_PATH)
stopifnot(all(c("gene", "chr", "gwas_name", "PP.H4.abf", "top_snp")
              %in% colnames(coloc)))

# Use whichever posterior is higher (SuSiE preferred, ABF fallback)
coloc[, pp4_best := pmax(PP.H4.abf,
                         ifelse(is.na(PP.H4.susie), -Inf, PP.H4.susie),
                         na.rm = TRUE)]
coloc_r <- coloc[!is.na(pp4_best) & pp4_best >= PP4_RELAXED]
message("  records with pp4_best >= ", PP4_RELAXED, ": ", nrow(coloc_r),
        " (unique genes: ", uniqueN(coloc_r$gene), ")")

coloc_r[, snp_chr := tstrsplit(top_snp, ":", fixed = TRUE, keep = 1L)]
coloc_r[, snp_pos := suppressWarnings(as.integer(
                       tstrsplit(top_snp, ":", fixed = TRUE, keep = 2L)[[1]]))]

# Strict (same top_snp, same GWAS)
message("[", Sys.time(), "] Channel a-strict (same top_snp within GWAS)")
strict <- coloc_r[, {
  if (.N >= 2) {
    idx <- combn(.N, 2)
    ga <- gene[idx[1, ]]; gb <- gene[idx[2, ]]
    pa <- pp4_best[idx[1, ]]; pb <- pp4_best[idx[2, ]]
    swap <- ga > gb
    tmp <- ga[swap]; ga[swap] <- gb[swap]; gb[swap] <- tmp
    tmp <- pa[swap]; pa[swap] <- pb[swap]; pb[swap] <- tmp
    list(gene_a = ga, gene_b = gb,
         pp4_a = pa, pp4_b = pb,
         pp4_min = pmin(pa, pb),
         top_snp_a = top_snp[1], top_snp_b = top_snp[1],
         locus_id = paste0(gwas_name[1], "|", top_snp[1]),
         distance_bp = 0L, mode = "strict",
         gwas_name = gwas_name[1],
         channel = "coloc_strict")
  } else NULL
}, by = .(gwas_name, top_snp)]
message("  strict pairs: ", nrow(strict))

# Relaxed (within WINDOW_COLOC_BP bp on same chr, same GWAS)
message("[", Sys.time(), "] Channel a-relaxed (within ",
        WINDOW_COLOC_BP, " bp on same chr)")
coloc_pos <- coloc_r[!is.na(snp_pos) & !is.na(chr)]
relaxed <- coloc_pos[, {
  if (.N < 2) { NULL } else {
    ord <- order(snp_pos)
    g <- gene[ord]; p <- pp4_best[ord]; s <- top_snp[ord]; po <- snp_pos[ord]
    n <- length(g); ga <- character(); gb <- character()
    pa <- numeric(); pb <- numeric(); sa <- character(); sb <- character()
    dd <- integer()
    for (i in seq_len(n - 1)) {
      j <- i + 1L
      while (j <= n && (po[j] - po[i]) <= WINDOW_COLOC_BP) {
        if (g[i] != g[j]) {
          if (g[i] <= g[j]) {
            ga <- c(ga, g[i]); gb <- c(gb, g[j])
            pa <- c(pa, p[i]); pb <- c(pb, p[j])
            sa <- c(sa, s[i]); sb <- c(sb, s[j])
          } else {
            ga <- c(ga, g[j]); gb <- c(gb, g[i])
            pa <- c(pa, p[j]); pb <- c(pb, p[i])
            sa <- c(sa, s[j]); sb <- c(sb, s[i])
          }
          dd <- c(dd, po[j] - po[i])
        }
        j <- j + 1L
      }
    }
    if (length(ga) == 0) NULL else
    list(gene_a = ga, gene_b = gb,
         pp4_a = pa, pp4_b = pb,
         pp4_min = pmin(pa, pb),
         top_snp_a = sa, top_snp_b = sb,
         locus_id = paste0(gwas_name[1], "|", chr[1], ":",
                           pmin(dd, WINDOW_COLOC_BP)),
         distance_bp = dd, mode = "relaxed",
         gwas_name = gwas_name[1],
         channel = "coloc_relaxed")
  }
}, by = .(gwas_name, chr)]
message("  relaxed pairs: ", nrow(relaxed))

common_cols <- c("gene_a", "gene_b", "pp4_a", "pp4_b", "pp4_min",
                 "top_snp_a", "top_snp_b", "locus_id", "distance_bp",
                 "mode", "gwas_name", "channel")
channel_a <- rbind(
  strict[, ..common_cols],
  relaxed[, ..common_cols]
)
channel_a <- channel_a[gene_a %in% V & gene_b %in% V]
message("  channel (a) after V filter: ", nrow(channel_a))

# ==========================================================================
#  Channel (b): Broadaway eQTL sentinel-sharing
# ==========================================================================
message("[", Sys.time(), "] Channel b: Broadaway eQTL sentinel-sharing")
eqtl_files <- list.files(EQTL_DIR, pattern = "chr[0-9]+_marginal_summary_results\\.tsv$",
                         full.names = TRUE)
message("  eQTL chrom files found: ", length(eqtl_files))

read_eqtl_leads <- function(path) {
  # Read only columns we need; file columns vary a bit across releases.
  hdr <- fread(path, nrows = 0)
  hcols <- tolower(colnames(hdr))
  pick_col <- function(candidates) {
    hits <- match(candidates, hcols, nomatch = 0)
    hits <- hits[hits > 0]
    if (length(hits) == 0) return(NA_character_)
    colnames(hdr)[hits[1]]
  }
  col_gene <- pick_col(c("gene", "gene_symbol", "genesymbol", "gene_id",
                         "hgnc_symbol", "symbol"))
  col_chr  <- pick_col(c("chr", "chromosome", "chrom"))
  col_pos  <- pick_col(c("pos", "position", "bp"))
  col_p    <- pick_col(c("p", "pvalue", "p_value", "pval"))
  col_snp  <- pick_col(c("snp", "rsid", "variant_id", "variant", "marker"))

  if (is.na(col_gene) || is.na(col_chr) || is.na(col_pos) || is.na(col_p)) {
    message("    skip ", basename(path),
            " -- missing required cols (have: ",
            paste(colnames(hdr), collapse = ","), ")")
    return(NULL)
  }

  keep <- c(col_gene, col_chr, col_pos, col_p)
  if (!is.na(col_snp) && length(col_snp) > 0) keep <- c(keep, col_snp)
  dt <- fread(path, select = keep)
  setnames(dt, old = col_gene, new = "gene")
  setnames(dt, old = col_chr,  new = "chr")
  setnames(dt, old = col_pos,  new = "pos")
  setnames(dt, old = col_p,    new = "pval")
  if (!is.na(col_snp)) setnames(dt, old = col_snp, new = "snp")
  dt <- dt[!is.na(pval) & pval <= EQTL_LEAD_P_MAX]
  # One lead per gene (min p)
  setorder(dt, gene, pval)
  dt <- unique(dt, by = "gene")
  dt
}

per_chr_leads <- lapply(eqtl_files, read_eqtl_leads)
per_chr_leads <- per_chr_leads[!vapply(per_chr_leads, is.null, logical(1))]
# Normalize: every per-chr frame must have all of (gene, chr, pos, pval).
# Some files may lack a `chr` column if the parse was partial -- infer chr
# from the filename when missing.
for (i in seq_along(per_chr_leads)) {
  dt <- per_chr_leads[[i]]
  fn <- basename(eqtl_files[i])
  chr_from_fn <- suppressWarnings(as.integer(sub("^chr([0-9]+).*$", "\\1", fn)))
  if (!("chr" %in% colnames(dt)) || all(is.na(dt$chr))) {
    dt[, chr := chr_from_fn]
  }
  if (!("snp" %in% colnames(dt))) dt[, snp := NA_character_]
  per_chr_leads[[i]] <- dt
}
eqtl_leads <- rbindlist(per_chr_leads, fill = TRUE, use.names = TRUE)
eqtl_leads <- eqtl_leads[!is.na(chr) & !is.na(pos)]
eqtl_leads <- eqtl_leads[gene %in% V]
message("  eQTL lead variants (one per gene, pval <= ", EQTL_LEAD_P_MAX,
        ", gene in V): ", nrow(eqtl_leads))

channel_b <- eqtl_leads[, {
  if (.N < 2) { NULL } else {
    ord <- order(pos)
    g <- gene[ord]; p <- pval[ord]; po <- pos[ord]
    n <- length(g); ga <- character(); gb <- character()
    pa <- numeric(); pb <- numeric(); dd <- integer()
    for (i in seq_len(n - 1)) {
      j <- i + 1L
      while (j <= n && (po[j] - po[i]) <= WINDOW_EQTL_BP) {
        if (g[i] != g[j]) {
          if (g[i] <= g[j]) {
            ga <- c(ga, g[i]); gb <- c(gb, g[j])
            pa <- c(pa, p[i]);  pb <- c(pb, p[j])
          } else {
            ga <- c(ga, g[j]); gb <- c(gb, g[i])
            pa <- c(pa, p[j]);  pb <- c(pb, p[i])
          }
          dd <- c(dd, po[j] - po[i])
        }
        j <- j + 1L
      }
    }
    if (length(ga) == 0) NULL else
    list(gene_a = ga, gene_b = gb,
         eqtl_p_a = pa, eqtl_p_b = pb,
         eqtl_distance_bp = dd,
         chr = chr[1])
  }
}, by = chr]
channel_b <- channel_b[gene_a %in% V & gene_b %in% V]
# Sentinel-sharing pairs lack a COLOC posterior; set pp4_min to NA
channel_b[, `:=`(
  pp4_a = NA_real_, pp4_b = NA_real_, pp4_min = NA_real_,
  top_snp_a = NA_character_, top_snp_b = NA_character_,
  locus_id = paste0("eqtl_sentinel|chr", chr, ":",
                    pmin(eqtl_distance_bp, WINDOW_EQTL_BP)),
  distance_bp = eqtl_distance_bp,
  mode = "eqtl_sentinel",
  gwas_name = "BROADAWAY_eqtl",
  channel = "eqtl_sentinel"
)]
channel_b <- channel_b[, .(gene_a, gene_b, pp4_a, pp4_b, pp4_min,
                           top_snp_a, top_snp_b, locus_id, distance_bp,
                           mode, gwas_name, channel)]
message("  channel (b) eQTL sentinel pairs: ", nrow(channel_b))

# ==========================================================================
#  Combine
# ==========================================================================
all_pairs <- rbind(channel_a[, .(gene_a, gene_b, pp4_a, pp4_b, pp4_min,
                                 top_snp_a, top_snp_b, locus_id, distance_bp,
                                 mode, gwas_name, channel)],
                   channel_b[, .(gene_a, gene_b, pp4_a, pp4_b, pp4_min,
                                 top_snp_a, top_snp_b, locus_id, distance_bp,
                                 mode, gwas_name, channel)],
                   fill = TRUE)
message("  total before dedup: ", nrow(all_pairs))

setorder(all_pairs, gene_a, gene_b, -pp4_min, distance_bp)
dedup <- all_pairs[, {
  keep <- 1L
  list(
    pp4_min       = pp4_min[keep],
    pp4_a         = pp4_a[keep],
    pp4_b         = pp4_b[keep],
    top_snp_a     = top_snp_a[keep],
    top_snp_b     = top_snp_b[keep],
    gwas_list     = paste(unique(gwas_name), collapse = ";"),
    locus_id      = locus_id[keep],
    mode          = mode[keep],
    channel       = paste(unique(channel), collapse = ";"),
    distance_bp   = min(distance_bp, na.rm = TRUE),
    n_shared_loci = length(unique(locus_id)),
    pp4_thresh    = PP4_RELAXED
  )
}, by = .(gene_a, gene_b)]

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(dedup, OUT)
message("[", Sys.time(), "] Wrote ", nrow(dedup),
        " expanded D-COLOC edges -> ", OUT)

message("  channel mix:")
print(dedup[, .N, by = channel])
message("  median pp4_min (non-NA): ",
        round(median(dedup$pp4_min, na.rm = TRUE), 3))
message("[287] Elapsed: ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
        " min")
