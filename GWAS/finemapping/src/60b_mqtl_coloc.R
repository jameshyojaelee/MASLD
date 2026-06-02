#!/usr/bin/env Rscript
# 60b_mqtl_coloc.R
# ---------------------------------------------------------------------------
# Colocalize MASLD-specific GWAS (credible-set loci) vs plasma metabolite-QTL
# (Chen 2023) / lipid-species-QTL (Ottensmann 2023) with coloc.abf().
#
# COMPLEMENTARY (supporting) layer — NOT a primary genetic-causal tier.
#
# Design:
#   - Both metabolite/lipid sumstats are GRCh38 (hg38); MASLD GWAS are hg19.
#     -> lift the metabolite per-locus SNPs hg38->hg19 (chain in data/broadaway_eqtl).
#   - For each metabolite file (this shard) x each MASLD credible-set locus:
#       read metabolite SNPs in the hg38 window, lift to hg19, merge with the
#       MASLD GWAS locus, run coloc.abf (GWAS type per registry; metabolite=quant).
#   - GWAS-vs-metabolite-GWAS coloc: both sides are summary stats. ABF only
#     (first pass; MR is banned — this is colocalization).
#
# Usage (SLURM array): SHARD index over the combined download list.
#   QTL_SOURCE = chen | ottensmann   (which arm)
#   SHARD_ID / N_SHARDS  (0-indexed stride over metabolite files in this arm)
#
# Output: data/external/{chen2023_mqtl,lipidqtl}/coloc_out/<source>_shard<ID>.tsv
# ---------------------------------------------------------------------------
suppressMessages({
  library(data.table); library(coloc)
  library(rtracklayer); library(GenomicRanges)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

QTL_SOURCE <- tolower(Sys.getenv("QTL_SOURCE", unset = "chen"))
SHARD_ID   <- as.integer(Sys.getenv("SHARD_ID",
                unset = if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else "0"))
N_SHARDS   <- as.integer(Sys.getenv("N_SHARDS",
                unset = if (length(commandArgs(TRUE)) >= 2) commandArgs(TRUE)[2] else "1"))

stopifnot(QTL_SOURCE %in% c("chen", "ottensmann"))

# --- paths per arm ---
if (QTL_SOURCE == "chen") {
  QTL_DIR    <- file.path(BASE_DIR, "data/external/chen2023_mqtl/sumstats")
  MAP_FILE   <- file.path(BASE_DIR, "data/external/chen2023_mqtl/chen2023_biomarker_download_list.tsv")
  OUT_DIR    <- file.path(BASE_DIR, "data/external/chen2023_mqtl/coloc_out")
  QTL_N_DEF  <- 8000L  # fallback per-trait N (Chen has per-trait N in download list)
} else {
  QTL_DIR    <- file.path(BASE_DIR, "data/external/lipidqtl/sumstats")
  MAP_FILE   <- file.path(BASE_DIR, "data/external/lipidqtl/ottensmann2023_accession_trait_map.tsv")
  OUT_DIR    <- file.path(BASE_DIR, "data/external/lipidqtl/coloc_out")
  QTL_N_DEF  <- 7174L
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LOCI_FILE  <- file.path(BASE_DIR, "data/external/chen2023_mqtl/masld_credset_loci.tsv")
CHAIN_FILE <- file.path(BASE_DIR, "data/broadaway_eqtl/hg38ToHg19.over.chain")

# --- coloc params (match 06_susie_coloc.R) ---
COLOC_P1 <- 1e-4; COLOC_P2 <- 1e-4; COLOC_P12 <- 5e-6
MIN_SNPS <- 50L          # min overlapping SNPs at a locus
PP4_KEEP <- 0.5          # report rows with PP.H4 > this (plus a thinned summary of all)

cat("============================================================\n")
cat("60b mQTL COLOC | source =", QTL_SOURCE, "| shard", SHARD_ID, "of", N_SHARDS, "\n")
cat("coloc", as.character(packageVersion("coloc")), "\n")
cat("============================================================\n")

# --- chain (hg38 -> hg19) ---
if (!file.exists(CHAIN_FILE)) {
  gz <- paste0(CHAIN_FILE, ".gz")
  if (file.exists(gz)) system(paste("gunzip -k", gz)) else stop("chain file missing: ", CHAIN_FILE)
}
chain <- import.chain(CHAIN_FILE)

# --- MASLD loci ---
loci <- fread(LOCI_FILE)
cat("MASLD credible-set loci:", nrow(loci), "across", uniqueN(loci$study), "GWAS\n")

# Cache GWAS sumstats per study (load once). Path comes from a study->path lookup
# (avoid data.table column/var name collision by precomputing the map).
gwas_path_map <- setNames(loci$sumstats_path, loci$study)
gwas_path_map <- gwas_path_map[!duplicated(names(gwas_path_map))]
gwas_cache <- new.env()
load_gwas <- function(study_name) {
  if (!is.null(gwas_cache[[study_name]])) return(gwas_cache[[study_name]])
  p <- gwas_path_map[[study_name]]
  full <- file.path(BASE_DIR, "GWAS/finemapping", p)
  g <- fread(full)
  setnames(g, tolower(names(g)))
  gwas_cache[[study_name]] <- g
  g
}

# --- metabolite file list for this shard ---
map <- fread(MAP_FILE)
acc_col <- "accession"
all_acc <- map[[acc_col]]
files <- file.path(QTL_DIR, paste0(all_acc, ".tsv.gz"))
present <- file.exists(files)
acc_present <- all_acc[present]; files <- files[present]
trait_lookup <- setNames(map$trait, map[[acc_col]])
n_lookup <- if ("Nnum" %in% names(map)) setNames(map$Nnum, map[[acc_col]]) else
            if ("sample" %in% names(map)) setNames(as.numeric(gsub("[^0-9]","", sub(" .*","", gsub(",","",map$sample)))), map[[acc_col]]) else NULL

# stride assignment
idx <- which((seq_along(acc_present) - 1L) %% N_SHARDS == SHARD_ID)
cat("Files present:", length(acc_present), "| this shard handles:", length(idx), "\n")
if (length(idx) == 0) { cat("Nothing to do for this shard.\n"); quit(save="no", status=0) }

# --- hg38 windows for prefiltering (lift the 90 MASLD hg19 loci -> hg38 once) ---
# Metabolite files are hg38; loci are hg19. Rather than lift ~15M variants per
# file, we lift the small set of locus windows hg19->hg38 and prefilter the
# metabolite file to those hg38 windows, then lift only the survivors back.
hg19to38_chain_file <- file.path(BASE_DIR, "data/broadaway_eqtl/hg19ToHg38.over.chain")
if (!file.exists(hg19to38_chain_file)) {
  gz <- paste0(hg19to38_chain_file, ".gz")
  if (file.exists(gz)) system(paste("gunzip -k", gz)) else stop("hg19ToHg38 chain missing")
}
chain_fwd <- import.chain(hg19to38_chain_file)
win_gr19 <- GRanges(seqnames = paste0("chr", loci$chromosome),
                    ranges = IRanges(start = loci$win_start, end = loci$win_end))
win_lift <- unlist(liftOver(win_gr19, chain_fwd))
hg38_win <- data.table(chr = as.integer(gsub("chr","",as.character(seqnames(win_lift)))),
                       start = start(win_lift), end = end(win_lift))
# Pad windows generously (+/-50kb) so liftOver edge gaps do not drop SNPs.
hg38_win[, start := pmax(1, start - 5e4)][, end := end + 5e4]
hg38_win <- hg38_win[!is.na(chr)]
setkey(hg38_win, chr, start, end)
cat("hg38 prefilter windows:", nrow(hg38_win), "across chr", paste(sort(unique(hg38_win$chr)), collapse=","), "\n")

in_hg38_window <- function(d) {
  # keep metabolite SNPs (hg38) falling in ANY MASLD locus hg38 window.
  # Vectorized interval join via data.table::foverlaps (fast for ~15M rows).
  dd <- d[chr %in% hg38_win$chr]
  if (nrow(dd) == 0) return(d[0])
  dd[, idx := .I]
  dd[, `:=`(start = pos, end = pos)]
  setkey(dd, chr, start, end)
  ov <- foverlaps(dd, hg38_win, by.x = c("chr","start","end"),
                  by.y = c("chr","start","end"), nomatch = NULL, which = TRUE)
  hit_idx <- unique(dd$idx[ov$xid])
  out <- dd[idx %in% hit_idx]
  out[, c("idx","start","end") := NULL]
  out
}

# Helper: read metabolite SNPs within an hg38 window from a tabix-less gz file.
# Files are chr-sorted; we read fully once per file (they're ~per-metabolite,
# manageable) then subset per locus. Column names per arm differ.
read_qtl <- function(f) {
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  setnames(d, tolower(names(d)))
  # harmonize column names
  cn <- names(d)
  ren <- function(from, to) { for (x in from) if (x %in% cn && !(to %in% names(d))) setnames(d, x, to) }
  ren(c("chromosome","chr","#chrom","chrom"), "chr")
  ren(c("base_pair_location","position","pos","bp"), "pos")
  ren(c("effect_allele","ea","alt"), "ea")
  ren(c("other_allele","oa","nea","ref"), "oa")
  ren(c("beta","effect"), "beta")
  ren(c("standard_error","se","standarderror"), "se")
  ren(c("p_value","pval","pvalue","p"), "pval")
  ren(c("neg_log_10_p_value","neg_log10_p","log10p"), "neglog10p")
  ren(c("n","sample_size"), "n_snp")
  if (!"pval" %in% names(d) && "neglog10p" %in% names(d)) d[, pval := 10^(-neglog10p)]
  needed <- c("chr","pos","ea","oa","beta","se")
  if (!all(needed %in% names(d))) return(NULL)
  d[, chr := as.integer(gsub("chr","",chr))]
  d <- d[!is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(se) & se > 0]
  d
}

results <- list()

for (k in idx) {
  acc <- acc_present[k]; f <- files[k]
  trait <- trait_lookup[[acc]]
  qtl_n <- if (!is.null(n_lookup) && !is.na(n_lookup[[acc]])) as.integer(n_lookup[[acc]]) else QTL_N_DEF
  qtl <- read_qtl(f)
  if (is.null(qtl)) { cat("  [skip]", acc, "(unreadable/columns)\n"); next }

  # Prefilter (hg38) to MASLD locus windows BEFORE lifting -> lift only survivors.
  n_full <- nrow(qtl)
  qtl <- in_hg38_window(qtl)
  if (nrow(qtl) == 0) { cat("  [skip]", acc, "(0 SNPs in MASLD windows)\n"); next }

  # lift the (small) prefiltered metabolite SNPs hg38 -> hg19, keep 1:1
  gr <- GRanges(seqnames = paste0("chr", qtl$chr),
                ranges = IRanges(start = qtl$pos, width = 1))
  mcols(gr) <- qtl[, .(ea, oa, beta, se, pval)]
  lifted <- liftOver(gr, chain)
  nmap <- lengths(lifted)
  keep <- which(nmap == 1)
  if (length(keep) == 0) { cat("  [skip]", acc, "(0 lifted)\n"); next }
  glift <- unlist(lifted[keep])
  qh19 <- as.data.table(mcols(glift))
  qh19[, chr := as.integer(gsub("chr","",as.character(seqnames(glift))))]
  qh19[, pos := start(glift)]
  setkey(qh19, chr, pos)

  n_loci_hit <- 0L
  for (li in seq_len(nrow(loci))) {
    L <- loci[li]
    qsub <- qh19[chr == L$chromosome & pos >= L$win_start & pos <= L$win_end]
    if (nrow(qsub) < MIN_SNPS) next

    g <- load_gwas(L$study)
    gsub <- g[chromosome == L$chromosome & position >= L$win_start & position <= L$win_end]
    if (nrow(gsub) < MIN_SNPS) next

    m <- merge(
      qsub[, .(pos, q_ea = ea, q_oa = oa, q_beta = beta, q_se = se)],
      gsub[, .(pos = position, g_a1 = allele1, g_a2 = allele2, g_beta = beta, g_se = se)],
      by = "pos"
    )
    if (nrow(m) < MIN_SNPS) next
    # allele harmonize (orient metabolite to GWAS effect allele)
    m[, match := (q_ea == g_a1 & q_oa == g_a2)]
    m[, flip  := (q_ea == g_a2 & q_oa == g_a1)]
    m <- m[match | flip]
    m[flip == TRUE, q_beta := -q_beta]
    # drop strand-ambiguous
    m <- m[!((g_a1 %in% c("A","T") & g_a2 %in% c("A","T")) |
             (g_a1 %in% c("C","G") & g_a2 %in% c("C","G")))]
    m <- m[!duplicated(pos)]
    if (nrow(m) < MIN_SNPS) next

    snpid <- paste(L$chromosome, m$pos, sep = ":")
    d1 <- list(beta = m$g_beta, varbeta = m$g_se^2, N = L$N_tot,
               type = L$coloc_type, snp = snpid)
    if (L$coloc_type == "cc" && !is.na(L$N_cases) && L$N_cases > 0) d1$s <- L$N_cases / L$N_tot
    if (L$coloc_type == "quant") d1$sdY <- 1
    d2 <- list(beta = m$q_beta, varbeta = m$q_se^2, N = qtl_n,
               type = "quant", sdY = 1, snp = snpid)

    ab <- tryCatch(suppressMessages(suppressWarnings(
      coloc.abf(d1, d2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12))),
      error = function(e) NULL)
    if (is.null(ab)) next
    pp <- ab$summary[paste0("PP.H", 0:4, ".abf")]
    top_snp <- NA_character_; top_pp <- NA_real_
    if (!is.null(ab$results)) {
      ti <- which.max(ab$results$SNP.PP.H4)
      top_snp <- ab$results$snp[ti]; top_pp <- ab$results$SNP.PP.H4[ti]
    }
    n_loci_hit <- n_loci_hit + 1L
    results[[length(results) + 1]] <- data.table(
      source = QTL_SOURCE, accession = acc, metabolite = trait,
      gwas = L$study, locus_id = L$locus_id, chr = L$chromosome,
      win_start = L$win_start, win_end = L$win_end, cs_lead_pos = L$cs_lead_pos,
      n_snps = nrow(m),
      PP.H0 = pp[1], PP.H1 = pp[2], PP.H2 = pp[3], PP.H3 = pp[4], PP.H4 = pp[5],
      top_snp = top_snp, top_snp_PP = top_pp
    )
  }
  cat(sprintf("  [%d/%d] %s | %s | full=%d prefilt=%d lifted=%d | locus tests=%d\n",
              k, length(acc_present), acc, substr(trait, 1, 40),
              n_full, nrow(qtl), length(keep), n_loci_hit))
}

if (length(results) == 0) {
  cat("No coloc results for this shard (no locus with >= MIN_SNPS overlap).\n")
  # still write empty file so the combiner sees the shard ran
  fwrite(data.table(source=character(), accession=character(), metabolite=character(),
                    gwas=character(), locus_id=character(), chr=integer(),
                    win_start=integer(), win_end=integer(), cs_lead_pos=integer(),
                    n_snps=integer(), PP.H0=double(), PP.H1=double(), PP.H2=double(),
                    PP.H3=double(), PP.H4=double(), top_snp=character(), top_snp_PP=double()),
         file.path(OUT_DIR, sprintf("%s_shard%d.tsv", QTL_SOURCE, SHARD_ID)), sep="\t")
  quit(save = "no", status = 0)
}

res <- rbindlist(results, fill = TRUE)
# keep all H4>0.1 rows (for downstream), thin the rest to save space
res_keep <- res[PP.H4 > 0.1]
out <- file.path(OUT_DIR, sprintf("%s_shard%d.tsv", QTL_SOURCE, SHARD_ID))
fwrite(res_keep, out, sep = "\t")
cat("\n============================================================\n")
cat("Shard", SHARD_ID, ":", nrow(res), "locus-metabolite tests |",
    "PP.H4>0.5:", sum(res$PP.H4 > 0.5, na.rm=TRUE),
    "| >0.8:", sum(res$PP.H4 > 0.8, na.rm=TRUE), "\n")
cat("Wrote", nrow(res_keep), "rows (PP.H4>0.1) to", out, "\n")
cat("============================================================\n")
