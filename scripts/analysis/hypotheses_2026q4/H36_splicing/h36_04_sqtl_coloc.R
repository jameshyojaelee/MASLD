#!/usr/bin/env Rscript
# H36 step 4: GTEx liver splicing QTL (eQTL Catalogue r7 QTD000270, leafcutter, .cc file:
# introns with a credible set) colocalized with the liver GWAS at the 274 H27 loci,
# and mapping of GTEx introns to the junctions of the stage test. PRESPEC_H36.md section 7.
#   coloc.abf (coloc 5.2.3), p1 = 1e-4, p2 = 1e-4, p12 = 5e-6, call PP.H4 > 0.5.
#   An intron is tested when min QTL p in the window < 1e-4 and >= 50 variants are shared
#   after hg19 -> hg38 liftOver of the GWAS, matched on position + unordered allele pair.
#   QTL: type quant, sdY = 1, N = max(an) / 2. GWAS: registry type/N/cases, window = lead +/- 500 kb.
# Usage: Rscript h36_04_sqtl_coloc.R <repo> <build_dir> <qtl_dir> <out_dir>
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer); library(coloc); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
REPO <- args[1]; BUILD <- args[2]; QDIR <- args[3]; OUT <- args[4]
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
TABIX <- "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/tabix"
CHAIN <- file.path(REPO, "data/broadaway_eqtl/hg19ToHg38.over.chain")
FM <- file.path(REPO, "GWAS/finemapping")
QFILE <- file.path(QDIR, "QTD000270.cc.tsv.gz")
LOCI <- file.path(REPO, "Analysis/Hypotheses/H27-20261008T230909Z/results/loci_h27.tsv")
WIN <- 500000L; P_QTL_MAX <- 1e-4; MIN_SNPS <- 50L; P1 <- 1e-4; P2 <- 1e-4; P12 <- 5e-6
QCOLS <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf",
           "pvalue", "beta", "se", "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")
set.seed(20260923)
t0 <- Sys.time()
elapsed <- function() round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
akey <- function(a, b) ifelse(a < b, paste(a, b, sep = "_"), paste(b, a, sep = "_"))

# ---- A. GTEx introns and the coordinate offset to our junctions -----------------
hdr <- readLines(pipe(sprintf("zcat %s | head -1", shQuote(QFILE))))
skip <- if (grepl("^molecular_trait_id", hdr)) 1L else 0L
introns <- fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NR>%d && !seen[$1]++ {print $1\"\\t\"$17}'", shQuote(QFILE), skip),
                 header = FALSE, col.names = c("id", "gene_id"), sep = "\t")
sp <- tstrsplit(introns$id, ":", fixed = TRUE)
introns[, `:=`(chrom = paste0("chr", sp[[1]]), s = as.integer(sp[[2]]), e = as.integer(sp[[3]]), clu = sp[[4]])]
introns[, strand := sub(".*_", "", clu)]
cat("GTEx liver leafcutter introns with a credible set:", nrow(introns), "\n")

jall <- fread(file.path(BUILD, "junctions_all_ge5_persons.tsv.gz"))
jk <- fread(file.path(BUILD, "junctions.tsv.gz"))
jall[, key := paste(chr, start, end, sep = ":")]
jk[, key := paste(chr, start, end, sep = ":")]
ann_keys <- jall[annotated == 1L, key]
offs <- CJ(ds = -1:1, de = -1:1)
offs[, c("n_match_annotated", "n_match_any") := {
  r <- lapply(seq_len(.N), function(i) {
    k <- paste(introns$chrom, introns$s + ds[i], introns$e + de[i], sep = ":")
    c(sum(k %in% ann_keys), sum(k %in% jall$key))
  })
  list(vapply(r, `[`, numeric(1), 1L), vapply(r, `[`, numeric(1), 2L))
}]
offs[, share_any := n_match_any / nrow(introns)]
fwrite(offs, file.path(OUT, "sqtl_offset_check.tsv"), sep = "\t")
print(offs)
best <- offs[which.max(n_match_annotated)]
cat(sprintf("best offset: start %+d, end %+d; %.1f%% of introns match a junction (rule: >= 80%%)\n",
            best$ds, best$de, 100 * best$share_any))
map_ok <- best$share_any >= 0.80
introns[, key := paste(chrom, s + best$ds, e + best$de, sep = ":")]
m_all <- jall[, .(key, strand_ours = strand, annotated_ours = annotated, reads_ours = reads, npers_ours = npers)]
m_all <- m_all[!duplicated(key)]
introns <- merge(introns, m_all, by = "key", all.x = TRUE, sort = FALSE)
m_k <- jk[, .(key, cl, jr, strand_kept = strand)][!duplicated(key)]
introns <- merge(introns, m_k, by = "key", all.x = TRUE, sort = FALSE)
introns[, map_class := fifelse(!is.na(cl), "kept_junction_in_tested_cluster",
                       fifelse(!is.na(npers_ours), "junction_present_not_in_tested_cluster", "not_in_our_junctions"))]
cat("mapping classes:\n"); print(introns[, .N, by = map_class])
cat("strand agreement among matched:", round(mean(introns[!is.na(strand_ours), strand == strand_ours]), 4), "\n")
fwrite(introns[, .(molecular_trait_id = id, gene_id, chrom, s_gtex = s, e_gtex = e, clu, strand, key,
                   map_class, cl, jr, annotated_ours, reads_ours, npers_ours)],
       file.path(OUT, "sqtl_intron_mapping.tsv.gz"), sep = "\t")
writeLines(c(if (map_ok) "PASS" else "FAIL",
             sprintf("best_offset_start=%d best_offset_end=%d share_any=%.4f", best$ds, best$de, best$share_any)),
           file.path(OUT, "GATE_mapping.txt"))
cat(sprintf("[%s min] mapping %s\n", elapsed(), if (map_ok) "PASS" else "FAIL"))

# ---- B. loci and GWAS windows ---------------------------------------------------
loci <- fread(LOCI)[status == "test"]
loci[, `:=`(wchr = as.integer(lead_chr), wstart = pmax(1L, as.integer(lead_pos) - WIN), wend = as.integer(lead_pos) + WIN)]
registry <- fread(file.path(FM, "config/gwas_registry.tsv"))
af_src <- fread(file.path(FM, "config/gwas_af_sources.tsv"))
liver_path <- function(study) {
  hit <- af_src[study_name == study & file.exists(af_sumstats_path)]
  if (nrow(hit) == 1L) return(hit$af_sumstats_path)
  file.path(FM, registry[study_name == study]$sumstats_path)
}
chain <- import.chain(CHAIN)
lift38 <- function(chr, pos) {
  lo <- liftOver(GRanges(paste0("chr", chr), IRanges(pos, width = 1)), chain)
  ok <- lengths(lo) == 1L; out <- rep(NA_integer_, length(pos)); out[ok] <- start(unlist(lo[ok])); out
}
gw <- list()
for (st in unique(loci$liver_study)) {
  w <- loci[liver_study == st, .(locus_id, chr = wchr, start = wstart, end = wend)]
  p <- liver_path(st)
  h <- names(fread(p, nrows = 0))
  d <- fread(p, select = intersect(c("chromosome", "position", "allele1", "allele2", "beta", "se", "pval"), h))
  d[, chromosome := suppressWarnings(as.integer(chromosome))]
  d <- d[!is.na(chromosome) & !is.na(beta) & !is.na(se) & se > 0 & !is.na(pval)]
  d[, row_id__ := .I]
  hit <- d[w, on = .(chromosome = chr, position >= start, position <= end), nomatch = NULL,
           .(row_id = x.row_id__, locus_id = i.locus_id)]
  g <- d[hit$row_id]; g[, locus_id := hit$locus_id]; g[, row_id__ := NULL]
  rm(d, hit); gc(FALSE)
  up <- unique(g[, .(chromosome, position)])
  up[, pos38 := lift38(chromosome, position)]
  g <- merge(g, up, by = c("chromosome", "position"))[!is.na(pos38)]
  g[, `:=`(allele1 = toupper(allele1), allele2 = toupper(allele2))]
  g[, k := paste(pos38, akey(allele1, allele2), sep = ":")]
  g <- g[!duplicated(g[, .(locus_id, k)])]
  for (lid in unique(g$locus_id)) gw[[lid]] <- g[locus_id == lid, .(k, pos38, beta.g = beta, se.g = se, p.g = pval)]
  cat(sprintf("[%s min] GWAS %s: %d rows in %d loci\n", elapsed(), st, nrow(g), uniqueN(g$locus_id)))
  rm(g, up); gc(FALSE)
}

# ---- C. colocalization per locus ------------------------------------------------
do_locus <- function(i) {
  lc <- loci[i]; g <- gw[[lc$locus_id]]
  stt <- function(status, ...) data.table(locus_id = lc$locus_id, liver_study = lc$liver_study, status = status, ...)
  if (is.null(g) || nrow(g) < MIN_SNPS) return(list(status = stt("too_few_gwas_variants"), res = NULL))
  r <- registry[study_name == lc$liver_study]
  gw_type <- if (r$trait_type == "binary") "cc" else "quant"
  s38 <- min(g$pos38); e38 <- max(g$pos38)
  q <- tryCatch(fread(cmd = sprintf("%s %s %d:%d-%d", TABIX, shQuote(QFILE), lc$wchr, s38, e38), header = FALSE,
                      col.names = QCOLS, sep = "\t", na.strings = c("NA", "")), error = function(e) data.table())
  if (!nrow(q)) return(list(status = stt("no_qtl_rows_in_window", gwas_min_p = min(g$p.g)), res = NULL))
  q <- q[!is.na(beta) & !is.na(se) & se > 0]
  q[, k := paste(position, akey(ref, alt), sep = ":")]
  tr <- q[, .(minp = min(pvalue), n_qtl_variants = .N, an = max(an)), by = molecular_trait_id]
  cand <- tr[minp < P_QTL_MAX]
  out <- list()
  for (t in cand$molecular_trait_id) {
    qt <- q[molecular_trait_id == t, .(k, beta.q = beta, se.q = se)][!duplicated(k)]
    m <- merge(g, qt, by = "k")
    if (nrow(m) < MIN_SNPS) next
    d1 <- list(beta = m$beta.g, varbeta = m$se.g^2, N = r$N_tot, type = gw_type, snp = m$k)
    if (gw_type == "cc") d1$s <- r$N_cases / r$N_tot else d1$sdY <- 1
    d2 <- list(beta = m$beta.q, varbeta = m$se.q^2, N = cand[molecular_trait_id == t]$an / 2, type = "quant", sdY = 1, snp = m$k)
    invisible(capture.output(cr <- suppressWarnings(suppressMessages(coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12)))))
    top <- cr$results[which.max(cr$results$SNP.PP.H4), ]
    out[[t]] <- data.table(locus_id = lc$locus_id, liver_study = lc$liver_study, molecular_trait_id = t,
                           min_qtl_p = cand[molecular_trait_id == t]$minp, nsnps = nrow(m),
                           PP.H0 = cr$summary[["PP.H0.abf"]], PP.H1 = cr$summary[["PP.H1.abf"]],
                           PP.H2 = cr$summary[["PP.H2.abf"]], PP.H3 = cr$summary[["PP.H3.abf"]],
                           PP.H4 = cr$summary[["PP.H4.abf"]], top_snp_hg38 = top$snp, top_snp_pp_h4 = top$SNP.PP.H4,
                           gwas_min_p_window = min(g$p.g))
  }
  list(status = stt("done", gwas_min_p = min(g$p.g), n_introns_window = nrow(tr), n_introns_p_lt_1e4 = nrow(cand),
                    n_tested = length(out)),
       res = rbindlist(out))
}
cat(sprintf("[%s min] colocalization over %d loci with %d workers\n", elapsed(), nrow(loci), NC))
L <- mclapply(seq_len(nrow(loci)), do_locus, mc.cores = NC, mc.preschedule = FALSE)
bad <- vapply(L, function(x) !is.list(x), logical(1))
if (any(bad)) stop("worker failure for loci: ", paste(loci$locus_id[bad], collapse = ","))
status <- rbindlist(lapply(L, `[[`, "status"), fill = TRUE)
res <- rbindlist(lapply(L, `[[`, "res"), fill = TRUE)
res <- merge(res, introns[, .(molecular_trait_id = id, gene_id, clu, strand, map_class, cl, jr, annotated_ours)],
             by = "molecular_trait_id", all.x = TRUE)
setorder(res, locus_id, -PP.H4)
fwrite(status, file.path(OUT, "sqtl_locus_status.tsv"), sep = "\t")
fwrite(res, file.path(OUT, "sqtl_coloc_locus_intron.tsv.gz"), sep = "\t")

summ <- data.table(
  loci_tested = nrow(loci), loci_with_tested_intron = uniqueN(res$locus_id),
  locus_intron_pairs_tested = nrow(res), distinct_introns_tested = uniqueN(res$molecular_trait_id),
  pairs_PP4_gt_0.5 = sum(res$PP.H4 > 0.5), loci_with_PP4_gt_0.5 = uniqueN(res[PP.H4 > 0.5]$locus_id),
  distinct_introns_PP4_gt_0.5 = uniqueN(res[PP.H4 > 0.5]$molecular_trait_id),
  coloc_pairs_mapped_to_tested_cluster = sum(res$map_class == "kept_junction_in_tested_cluster", na.rm = TRUE),
  coloc_PP4_gt_0.5_mapped_to_tested_cluster = sum(res$PP.H4 > 0.5 & res$map_class == "kept_junction_in_tested_cluster", na.rm = TRUE),
  mapping_share_any = best$share_any, offset_start = best$ds, offset_end = best$de)
fwrite(summ, file.path(OUT, "sqtl_coloc_summary.tsv"), sep = "\t")
print(t(summ))
print(status[, .N, by = status])
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_04.txt"))
cat(sprintf("[%s min] done\n", elapsed()))
