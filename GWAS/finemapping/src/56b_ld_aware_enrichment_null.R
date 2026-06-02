#!/usr/bin/env Rscript
# 56b_ld_aware_enrichment_null.R
# LD-aware permutation null for cell-type ATAC peak enrichment of fine-mapped GWAS variants.
#
# Rationale
# ---------
# Script 57's current enrichment uses a genome-coverage Fisher null (background = peak
# fraction of GRCh38). This ignores that the 8,329 finemapped variants are clustered in
# LD blocks: nearby variants share peak overlap status, which inflates the effective
# sample size and produces artifactually small p-values / inflated ORs.
# (Macrophages: 1.47x, p=5.6e-9; Hepatocytes: 1.18x, p=3.9e-4 under Fisher null.)
#
# Method
# ------
# Berisa-Pickrell (2016) EUR LD blocks (1,704 blocks) + 1000 Genomes EUR variant pool
# (~8.5M variants, 379 samples). For each observed variant we record (LD block,
# MAF-bin). Each permutation draws, per observed variant, one replacement from the
# same (block, MAF-bin) cell of the 1KG pool, then re-intersects with each cell-type
# peak BED. N=10,000 perms, BiocParallel MulticoreParam workers=16. Empirical p =
# (1 + #(n_perm >= n_obs)) / (N + 1). log2 fold enrichment = log2(n_obs / median(n_perm)).
# Outputs `enrichment_ld_null.csv` + `enrichment_null_perm_matrix.rds`.
#
# Inputs
# ------
#   GWAS/finemapping/results/gwas_atac/variant_overlap_summary.csv  (8,329 unique vars)
#   GWAS/finemapping/data/ld_ref/1kg_eur/approx_LD_blocks.txt       (1,704 EUR blocks)
#   GWAS/finemapping/data/ld_ref/1kg_eur/chr{1..22}_eur.{bed,bim,fam}  (1KG EUR pool)
#   Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2/*.bed
#   data/broadaway_eqtl/hg19ToHg38.over.chain                       (UCSC liftover)
#
# Outputs
# -------
#   GWAS/finemapping/results/gwas_atac/enrichment_ld_null.csv
#     cell_type, n_obs, median_null, p_emp, p_emp_bh, log2FE, ci_lo, ci_hi
#   GWAS/finemapping/results/gwas_atac/enrichment_null_perm_matrix.rds
#     list(perm_matrix = matrix[n_ct, N], n_obs = named vec, cell_types = char vec)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(BiocParallel)
})

set.seed(42)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE_DIR, "GWAS/finemapping")
ATAC_DIR <- file.path(BASE_DIR, "Analysis/ATAC/Human_Multiome")
OUT_DIR  <- file.path(FM_DIR, "results/gwas_atac")
LD_DIR   <- file.path(FM_DIR, "data/ld_ref/1kg_eur")
CHAIN_F  <- file.path(BASE_DIR, "data/broadaway_eqtl/hg19ToHg38.over.chain")
PEAK_DIR <- file.path(ATAC_DIR, "results/label_transfer/cell_type_peak_sets_v2")

# FIX (review B6/56b): default to the session tempdir() instead of a hardcoded personal
# /scratch path that may not exist or be writable on every node.
TMP_DIR <- file.path(Sys.getenv("TMPDIR", unset = tempdir()),
                     paste0("ld_null_56b_", Sys.getpid()))
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Parameters ──────────────────────────────────────────────────────────────
N_PERM   <- as.integer(Sys.getenv("N_PERM", 10000))
N_MAF_BINS <- 10L
N_WORKERS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 16))

cat("==============================================================\n")
cat("56b_ld_aware_enrichment_null.R\n")
cat("==============================================================\n")
cat("N_PERM       :", N_PERM, "\n")
cat("N_MAF_BINS   :", N_MAF_BINS, "\n")
cat("N_WORKERS    :", N_WORKERS, "\n")
cat("TMP_DIR      :", TMP_DIR, "\n\n")

# ── 1. Load observed overlap summary ─────────────────────────────────────────
cat("--- 1. Loading observed variants ---\n")
ov <- fread(file.path(OUT_DIR, "variant_overlap_summary.csv"))
cat("Loaded", nrow(ov), "observed unique finemapped variants\n")
cat("Columns:", paste(colnames(ov), collapse = ", "), "\n")

stopifnot(all(c("chromosome", "position", "chr_hg38", "pos_hg38",
                "cell_types_overlapping") %in% colnames(ov)))
ov[, chromosome := as.integer(chromosome)]
ov <- ov[!is.na(chromosome) & chromosome >= 1 & chromosome <= 22]
cat("After autosomal filter:", nrow(ov), "\n\n")

# ── 2. Load LD blocks ────────────────────────────────────────────────────────
cat("--- 2. Loading LD blocks ---\n")
blocks <- fread(file.path(LD_DIR, "approx_LD_blocks.txt"))
setnames(blocks, c("chr", "start", "stop"))
blocks[, chr := as.integer(chr)]
blocks <- blocks[!is.na(chr) & chr >= 1 & chr <= 22]
blocks[, block_id := paste0(chr, "_", start, "_", stop)]
setkey(blocks, chr, start, stop)
cat("Loaded", nrow(blocks), "LD blocks\n\n")

# Function: assign chr+pos to LD block. Uses data.table foverlaps for vectorised speed.
assign_block <- function(chr_vec, pos_vec) {
  dt <- data.table(chr = as.integer(chr_vec), pos = as.integer(pos_vec))
  dt[, start := pos]; dt[, end := pos]
  setkey(dt, chr, start, end)
  hits <- foverlaps(dt, blocks, by.x = c("chr", "start", "end"),
                    by.y = c("chr", "start", "stop"),
                    nomatch = NA, mult = "first")
  hits$block_id
}

# ── 3. Load 1KG EUR pool from .bim, lift to hg38, compute MAF via plink2 ─────
cat("--- 3. Building 1KG EUR variant pool with MAF + hg38 + block ---\n")

# Step 3a: read all .bim files; combine into one data.table
read_chr_bim <- function(chr) {
  bim_f <- file.path(LD_DIR, paste0("chr", chr, "_eur.bim"))
  if (!file.exists(bim_f)) return(NULL)
  dt <- fread(bim_f, header = FALSE,
              col.names = c("chr", "snp_id", "cm", "pos_hg19", "A1", "A2"))
  dt[, c("cm") := NULL]
  dt
}

bim_list <- lapply(1:22, read_chr_bim)
bim <- rbindlist(bim_list)
bim[, chr := as.integer(chr)]
cat("Loaded", nrow(bim), "1KG EUR variants from bim files\n")

# Step 3b: compute MAF via plink2 --freq per chromosome
# Output: .afreq with ID CHROM POS REF ALT ALT_FREQS OBS_CT
plink_bin <- Sys.getenv("PLINK2_BIN", "")
if (!nzchar(plink_bin) || !file.exists(plink_bin)) {
  # Use module system
  cat("Attempting plink2 via module load...\n")
  plink_bin <- "plink2"
}

compute_maf_one_chr <- function(chr) {
  out_prefix <- file.path(TMP_DIR, paste0("chr", chr, "_freq"))
  in_prefix  <- file.path(LD_DIR, paste0("chr", chr, "_eur"))
  cmd <- paste("plink2 --bfile", shQuote(in_prefix),
               "--freq --silent --threads 2 --out", shQuote(out_prefix))
  rc <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (rc != 0) {
    # Try plink1.9 fallback
    cmd2 <- paste("plink --bfile", shQuote(in_prefix),
                  "--freq --silent --threads 2 --out", shQuote(out_prefix))
    rc <- system(cmd2, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (rc != 0) stop("plink/plink2 frequency computation failed for chr ", chr)
    f <- paste0(out_prefix, ".frq")
    af <- fread(f, header = TRUE)
    return(data.table(chr = as.integer(chr), snp_id = af$SNP, maf = af$MAF))
  }
  af_f <- paste0(out_prefix, ".afreq")
  if (!file.exists(af_f)) stop("Expected ", af_f)
  af <- fread(af_f, header = TRUE)
  # plink2 columns: #CHROM ID REF ALT PROVISIONAL_REF? ALT_FREQS OBS_CT
  # We compute MAF = min(ALT_FREQS, 1-ALT_FREQS)
  alt_col <- if ("ALT_FREQS" %in% colnames(af)) "ALT_FREQS" else "ALT_FREQ"
  id_col  <- if ("ID" %in% colnames(af)) "ID" else "SNP"
  maf <- pmin(af[[alt_col]], 1 - af[[alt_col]])
  data.table(chr = as.integer(chr), snp_id = af[[id_col]], maf = maf)
}

cat("Running plink2 --freq for chr 1..22 ...\n")
t0 <- Sys.time()
freq_list <- lapply(1:22, function(chr) {
  cat(sprintf("  chr%-2d ... ", chr))
  tryCatch({
    dt <- compute_maf_one_chr(chr)
    cat(nrow(dt), "variants\n")
    dt
  }, error = function(e) {
    cat("FAILED: ", conditionMessage(e), "\n")
    NULL
  })
})
freq_dt <- rbindlist(freq_list)
cat("Total freq variants:", nrow(freq_dt), "in", round(difftime(Sys.time(), t0, units = "secs"), 1), "sec\n")

# Step 3c: join bim + freq by (chr, snp_id), drop NA/monomorphic
setkey(bim, chr, snp_id)
setkey(freq_dt, chr, snp_id)
pool <- bim[freq_dt, nomatch = NULL]
pool <- pool[!is.na(maf) & maf > 0]
cat("After freq join + drop monomorphic:", nrow(pool), "variants\n")

# Step 3d: assign LD blocks
cat("Assigning LD blocks to pool ...\n")
pool[, block_id := assign_block(chr, pos_hg19)]
n_drop_no_block <- sum(is.na(pool$block_id))
pool <- pool[!is.na(block_id)]
cat("Dropped", n_drop_no_block, "variants without LD block; kept", nrow(pool), "\n")

# Step 3e: MAF bins (global deciles, computed across pool)
maf_breaks <- quantile(pool$maf, probs = seq(0, 1, length.out = N_MAF_BINS + 1),
                       na.rm = TRUE)
maf_breaks[1]  <- 0
maf_breaks[N_MAF_BINS + 1] <- 0.5 + 1e-9
pool[, maf_bin := as.integer(cut(maf, breaks = maf_breaks, include.lowest = TRUE,
                                  labels = FALSE))]
cat("MAF bin breaks:\n"); print(round(maf_breaks, 4))
cat("Pool MAF bin sizes:\n"); print(table(pool$maf_bin))

# Step 3f: liftover pool hg19 -> hg38
cat("\nLifting 1KG pool hg19 -> hg38 ...\n")
chain <- import.chain(CHAIN_F)
gr_pool <- GRanges(seqnames = paste0("chr", pool$chr),
                   ranges = IRanges(start = pool$pos_hg19, width = 1),
                   strand = "*")
mcols(gr_pool)$pool_idx <- seq_len(nrow(pool))
lifted_pool <- liftOver(gr_pool, chain)
lp_n <- lengths(lifted_pool)
keep_pool <- which(lp_n == 1)
gr_pool_hg38 <- unlist(lifted_pool[keep_pool])
pool_idx_keep <- mcols(gr_pool_hg38)$pool_idx
pool[, c("chr_hg38", "pos_hg38") := list(NA_character_, NA_integer_)]
pool$chr_hg38[pool_idx_keep] <- as.character(seqnames(gr_pool_hg38))
pool$pos_hg38[pool_idx_keep] <- start(gr_pool_hg38)
pool <- pool[!is.na(pos_hg38)]
cat("Pool after liftover:", nrow(pool), "(",
    round(100 * length(keep_pool) / length(gr_pool), 1), "% lifted)\n\n")

# ── 4. Assign observed variants to (block, MAF bin) ─────────────────────────
cat("--- 4. Assigning observed variants to (block, MAF bin) ---\n")

ov[, block_id := assign_block(chromosome, position)]
n_ov_no_block <- sum(is.na(ov$block_id))
cat("Observed variants without LD block:", n_ov_no_block,
    "/", nrow(ov), "(dropped)\n")
ov <- ov[!is.na(block_id)]

# Join MAF: prefer exact pos match in 1KG pool. Fallback = nearest pool variant in same block.
pool_pos_key <- pool[, .(chr, pos_hg19, maf, maf_bin)]
setkey(pool_pos_key, chr, pos_hg19)

# Exact match by (chr, pos_hg19)
ov_match <- pool_pos_key[ov, on = c(chr = "chromosome", pos_hg19 = "position"),
                          mult = "first"]
ov$maf_bin <- ov_match$maf_bin
ov$maf     <- ov_match$maf
n_exact <- sum(!is.na(ov$maf_bin))
cat("Exact pos match in 1KG:", n_exact, "/", nrow(ov), "\n")

# Fallback for unmatched: median MAF of variants in same block
need <- which(is.na(ov$maf_bin))
if (length(need) > 0) {
  block_med <- pool[, .(med_bin = as.integer(median(maf_bin))), by = block_id]
  blk_lookup <- setNames(block_med$med_bin, block_med$block_id)
  ov$maf_bin[need] <- blk_lookup[ov$block_id[need]]
  # Any still NA (block not in pool)? Drop.
  n_still_na <- sum(is.na(ov$maf_bin))
  if (n_still_na > 0) {
    cat("Dropping", n_still_na, "observed variants without recoverable MAF bin\n")
    ov <- ov[!is.na(maf_bin)]
  } else {
    cat("Filled", length(need), "via block-median MAF bin\n")
  }
}

cat("Final observed variant count:", nrow(ov), "\n\n")

# ── 5. Build per-(block, MAF bin) pool lookup ───────────────────────────────
cat("--- 5. Building (block, MAF-bin) -> pool indices lookup ---\n")
pool[, idx := seq_len(.N)]
pool_key <- paste(pool$block_id, pool$maf_bin, sep = "|")
pool_idx_by_cell <- split(pool$idx, pool_key)

# Check cell coverage for observed
ov[, key := paste(block_id, maf_bin, sep = "|")]
ov_cells <- unique(ov$key)
missing_cells <- setdiff(ov_cells, names(pool_idx_by_cell))
cat("Unique (block, bin) cells for observed:", length(ov_cells),
    "; missing in pool:", length(missing_cells), "\n")

# Fallback for missing cells: same block, any MAF bin
if (length(missing_cells) > 0) {
  for (cell in missing_cells) {
    blk <- sub("\\|.*", "", cell)
    cand <- pool$idx[pool$block_id == blk]
    if (length(cand) == 0) cand <- pool$idx  # ultra-fallback: any variant
    pool_idx_by_cell[[cell]] <- cand
  }
  cat("Filled", length(missing_cells), "missing cells with block-only fallback\n")
}

# Pre-extract per observed variant the list of candidate pool indices.
# We then flatten to a single concatenated integer vector + per-obs start offset
# so each permutation can draw all 8329 indices in one vectorised step.
ov_cand <- pool_idx_by_cell[ov$key]
ov_cand_lens <- lengths(ov_cand)
zero_cand <- which(ov_cand_lens == 0)
if (length(zero_cand) > 0) {
  cat("WARNING:", length(zero_cand),
      "observed variants have empty candidate pool; using full-pool fallback.\n")
  for (k in zero_cand) ov_cand[[k]] <- seq_len(nrow(pool))
  ov_cand_lens <- lengths(ov_cand)
}
cat("Mean candidate pool size per obs variant:",
    round(mean(ov_cand_lens), 1), "; median:", median(ov_cand_lens),
    "; min:", min(ov_cand_lens), "; max:", max(ov_cand_lens), "\n")

ov_cand_concat <- unlist(ov_cand, use.names = FALSE)
ov_cand_offset <- c(0L, cumsum(ov_cand_lens)[-length(ov_cand_lens)])
cat("Concatenated candidate vector length:", length(ov_cand_concat), "\n\n")

# ── 6. Load cell-type peaks ─────────────────────────────────────────────────
cat("--- 6. Loading cell-type ATAC peaks (hg38) ---\n")
peak_files <- list.files(PEAK_DIR, pattern = "_peaks\\.bed$", full.names = TRUE)
stopifnot(length(peak_files) > 0)
cell_types <- gsub("_peaks\\.bed$", "", basename(peak_files))

peak_gr_list <- list()
for (i in seq_along(peak_files)) {
  ct <- cell_types[i]
  bed <- fread(peak_files[i], header = FALSE)
  if (ncol(bed) == 3) setnames(bed, c("chr", "start", "end"))
  else { bed <- bed[, 2:4]; setnames(bed, c("chr", "start", "end")) }
  peak_gr_list[[ct]] <- GRanges(seqnames = bed$chr,
                                ranges = IRanges(start = bed$start, end = bed$end))
  cat("  ", ct, ":", length(peak_gr_list[[ct]]), "peaks\n")
}

# ── 7. Compute observed overlap counts (re-confirm) ─────────────────────────
cat("\n--- 7. Observed overlap counts (per cell type) ---\n")
gr_obs <- GRanges(seqnames = ov$chr_hg38,
                  ranges = IRanges(start = ov$pos_hg38, width = 1))
n_obs_vec <- sapply(peak_gr_list, function(gr) {
  length(unique(queryHits(findOverlaps(gr_obs, gr))))
})
print(n_obs_vec)

# Sanity-check vs upstream cell_types_overlapping field
cell_types_count_check <- sapply(cell_types, function(ct)
  sum(grepl(ct, ov$cell_types_overlapping)))
cat("\nUpstream cell_types_overlapping field counts:\n")
print(cell_types_count_check)
cat("(Small differences possible if upstream used different variants subset.)\n\n")

# ── 8. Permutation loop (BiocParallel) ──────────────────────────────────────
cat("--- 8. Permutation loop:", N_PERM, "perms x", length(cell_types),
    "cell types ---\n")

# Pre-build pool hg38 GRanges (used to fetch hg38 positions for permuted indices)
pool_chr_hg38 <- pool$chr_hg38
pool_pos_hg38 <- pool$pos_hg38

# Worker function: one batch of permutations (vectorised draw)
one_batch <- function(perm_ids, ov_cand_concat, ov_cand_offset, ov_cand_lens,
                      pool_chr_hg38, pool_pos_hg38,
                      peak_gr_list, cell_types, seed_base) {
  n_obs <- length(ov_cand_lens)
  out_mat <- matrix(0L, nrow = length(cell_types), ncol = length(perm_ids),
                    dimnames = list(cell_types, NULL))
  for (j in seq_along(perm_ids)) {
    pid <- perm_ids[j]
    set.seed(seed_base + pid)
    # one-shot draw: index_within_cand ~ Unif{1, lk}; flat index = offset + idx_within
    u <- runif(n_obs)
    idx_within <- pmin(as.integer(ceiling(u * ov_cand_lens)), ov_cand_lens)
    idx_within[idx_within < 1L] <- 1L
    flat_idx <- ov_cand_offset + idx_within
    draws <- ov_cand_concat[flat_idx]
    gr_p <- GRanges(seqnames = pool_chr_hg38[draws],
                    ranges = IRanges(start = pool_pos_hg38[draws], width = 1))
    for (ct in cell_types) {
      hits <- findOverlaps(gr_p, peak_gr_list[[ct]])
      out_mat[ct, j] <- length(unique(queryHits(hits)))
    }
  }
  out_mat
}

# Chunk perms into batches; one batch per worker (BiocParallel)
batch_size <- max(1L, ceiling(N_PERM / (N_WORKERS * 4L)))  # 4 batches per worker
batches <- split(seq_len(N_PERM),
                 ceiling(seq_len(N_PERM) / batch_size))
cat("Permutation batches:", length(batches), "of size ~", batch_size, "\n")

bp <- MulticoreParam(workers = N_WORKERS, RNGseed = 42, progressbar = TRUE)
t_perm <- Sys.time()
batch_results <- bplapply(batches, function(b) {
  one_batch(b, ov_cand_concat, ov_cand_offset, ov_cand_lens,
            pool_chr_hg38, pool_pos_hg38,
            peak_gr_list, cell_types, seed_base = 1000L)
}, BPPARAM = bp)
cat("Permutation wall time:", format(difftime(Sys.time(), t_perm)), "\n")

perm_mat <- do.call(cbind, batch_results)
dimnames(perm_mat) <- list(cell_types, NULL)
cat("perm_mat dim:", paste(dim(perm_mat), collapse = " x "), "\n\n")

# ── 9. Empirical p, log2FE, CI per cell type; BH-correct ────────────────────
cat("--- 9. Empirical statistics ---\n")
res <- data.table(cell_type = cell_types)
res[, n_obs := n_obs_vec[cell_type]]
res[, median_null := apply(perm_mat, 1, median)]
res[, mean_null   := apply(perm_mat, 1, mean)]
res[, p_emp := mapply(function(ct, n_o) {
  (1 + sum(perm_mat[ct, ] >= n_o)) / (N_PERM + 1)
}, cell_type, n_obs)]
res[, log2FE := log2(pmax(n_obs, 1) / pmax(median_null, 1))]
res[, ci_lo := log2(pmax(n_obs, 1) /
                    pmax(apply(perm_mat, 1, quantile, probs = 0.975, na.rm = TRUE), 1))]
res[, ci_hi := log2(pmax(n_obs, 1) /
                    pmax(apply(perm_mat, 1, quantile, probs = 0.025, na.rm = TRUE), 1))]
res[, p_emp_bh := p.adjust(p_emp, method = "BH")]
setorder(res, p_emp)

print(res)
fwrite(res, file.path(OUT_DIR, "enrichment_ld_null.csv"))
cat("\nWrote enrichment_ld_null.csv\n")

saveRDS(list(perm_matrix = perm_mat,
             n_obs = n_obs_vec,
             cell_types = cell_types,
             n_perm = N_PERM,
             n_maf_bins = N_MAF_BINS,
             n_obs_used = nrow(ov)),
        file.path(OUT_DIR, "enrichment_null_perm_matrix.rds"))
cat("Wrote enrichment_null_perm_matrix.rds\n")

# ── 10. Compare to existing Fisher result ───────────────────────────────────
cat("\n--- 10. Comparison to existing Fisher null ---\n")
fisher_f <- file.path(OUT_DIR, "enrichment_statistics.csv")
if (file.exists(fisher_f)) {
  fish <- fread(fisher_f)
  cmp <- merge(res[, .(cell_type, n_obs, median_null, log2FE,
                       p_emp, p_emp_bh)],
               fish[, .(cell_type, fisher_or, fisher_padj,
                        fold_enrichment_fisher = fold_enrichment)],
               by = "cell_type")
  cmp[, log2_fisher_or := log2(fisher_or)]
  setorder(cmp, p_emp)
  print(cmp)
  fwrite(cmp, file.path(OUT_DIR, "enrichment_ld_vs_fisher.csv"))
  cat("Wrote enrichment_ld_vs_fisher.csv\n")
}

# Clean tmp
unlink(TMP_DIR, recursive = TRUE)
cat("\n--- Done ---\n")
