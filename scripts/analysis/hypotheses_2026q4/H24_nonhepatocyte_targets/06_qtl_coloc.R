#!/usr/bin/env Rscript
# H24 step 6: colocalization layers. See PRESPEC_H24.md.
#   L0  GTEx v10 liver ge (QTD000266), bulk; descriptive at H24 loci and the
#       coloc machinery control at non-H24 direct loci with release coloc genes.
#   L3  non-liver cell-type eQTL proxies (eQTL Catalogue r8), mapped to lineages.
#   L4  Currin & Mohlke 2025 bulk-liver caQTL (local, hg38), at H24 variant peaks.
#   L5  GTEx v10 liver leafcutter sQTL (QTD000270, .cc file: traits with a
#       credible set only), bulk.
# coloc 5.x coloc.abf, priors p1 = 1e-4, p2 = 1e-4, p12 = 5e-6 (repo values).
# GWAS per physical locus: among the member clusters' direct studies with an
# _af sumstats file, the one with the smallest p in the window (span +/- 500 kb,
# hg19); untestable if that p > 1e-6. GWAS cc: s = N_cases/N_tot; quant: sdY = 1.
# QTL: type quant, sdY = 1 (normalized phenotypes), N = max(an)/2.
# Traits tested if min QTL p in window < 1e-4 and >= 50 shared variants
# (matched on hg38 position and unordered alleles after liftOver of GWAS).
# Remote layers are marked not_run if the eQTL Catalogue FTP is unreachable.
# Usage: Rscript 06_qtl_coloc.R <OUT_DIR>
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer); library(coloc); library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE); OUT <- normalizePath(args[[1]])
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
ENV <- "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
TABIX <- file.path(ENV, "bin/tabix")
CHAIN <- file.path(ROOT, "data/broadaway_eqtl/hg19ToHg38.over.chain")
GTF <- "/gpfs/commons/home/jameslee/reference_genome/cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
SUMS <- file.path(ROOT, "GWAS/finemapping/data/sumstats_af")
REG <- file.path(ROOT, "GWAS/finemapping/config/gwas_registry.tsv")
CURRIN <- file.path(ROOT, "data/external/currin_2025_caqtl")
FTP <- "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats"
WIN <- 500000L; P_GWAS_MAX <- 1e-6; P_QTL_MAX <- 1e-4; MIN_SNPS <- 50L; CURRIN_N <- 138L; PEAK_PAD <- 250L
P1 <- 1e-4; P2 <- 1e-4; P12 <- 5e-6
QCOLS <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf",
           "pvalue", "beta", "se", "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")
res_dir <- file.path(OUT, "results"); dir.create(res_dir, showWarnings = FALSE)
cache <- file.path(OUT, "qtl_cache"); dir.create(cache, showWarnings = FALSE)
tbi_dir <- file.path(cache, "tbi"); dir.create(tbi_dir, showWarnings = FALSE)

DS <- rbindlist(list(
  list("L0", "QTS000015", "QTD000266", "GTEx_v10 liver ge", "bulk_liver", "all"),
  list("L5", "QTS000015", "QTD000270", "GTEx_v10 liver leafcutter", "bulk_liver", "cc"),
  list("L3", "QTS000023", "QTD000404", "PhLiPS HLC ge", "hepatocyte", "all"),
  list("L3", "QTS000021", "QTD000384", "Nedelec_2016 macrophage naive ge", "macrophage", "all"),
  list("L3", "QTS000001", "QTD000001", "Alasoo_2018 macrophage naive ge", "macrophage", "all"),
  list("L3", "QTS000015", "QTD000216", "GTEx_v10 cultured fibroblast ge", "fibroblast", "all"),
  list("L3", "QTS000012", "QTD000100", "GENCORD fibroblast ge", "fibroblast", "all"),
  list("L3", "QTS000038", "QTD000612", "OneK1K CD4_Naive ge", "t_cell", "all"),
  list("L3", "QTS000038", "QTD000615", "OneK1K CD8_Naive ge", "t_cell", "all"),
  list("L3", "QTS000038", "QTD000620", "OneK1K NK ge", "t_cell", "all"),
  list("L3", "QTS000038", "QTD000608", "OneK1K B_naive ge", "b_cell", "all")))
setnames(DS, c("layer", "qts", "qtd", "label", "proxy_lineage", "kind"))
DS[, url := sprintf("%s/%s/%s/%s.%s.tsv.gz", FTP, qts, qtd, qtd, kind)]

# ---- loci ----
loci <- fread(file.path(OUT, "loci_direct_clusters.tsv"))
phys <- fread(file.path(OUT, "inventory_v2/physical_loci.tsv"))
merge_physical <- function(cl) {
  cl <- cl[order(chrom, span_start_hg19)]
  cl[, `:=`(ws = span_start_hg19 - WIN, we = span_end_hg19 + WIN)]
  grp <- integer(nrow(cl)); g <- 0L; cur_chr <- -1L; cur_end <- -Inf
  for (i in seq_len(nrow(cl))) {
    if (cl$chrom[i] != cur_chr || cl$ws[i] > cur_end) { g <- g + 1L; cur_chr <- cl$chrom[i]; cur_end <- cl$we[i] }
    else cur_end <- max(cur_end, cl$we[i])
    grp[i] <- g
  }
  cl[, physical_locus := paste(sort(cluster_id), collapse = "+"), by = .(g = grp)]
  cl[, .(cluster_id, physical_locus)]
}
ctrl <- loci[h24_primary == FALSE & mhc == FALSE & n_direct_coloc_genes > 0]
ctrl_phys <- merge_physical(copy(ctrl))
L <- rbind(phys[set == "H24", .(cluster_id, physical_locus, locus_set = "H24")],
           ctrl_phys[, .(cluster_id, physical_locus, locus_set = ifelse(grepl("^22:3|22:9", physical_locus), "REF_PNPLA3", "CTRL_L0"))])
L <- merge(L, loci[, .(cluster_id, chrom, span_start_hg19, span_end_hg19, direct_studies, direct_coloc_genes)], by = "cluster_id")
LOC <- L[, .(locus_set = locus_set[1], chrom = chrom[1],
             w19_start = max(1L, min(span_start_hg19) - WIN), w19_end = max(span_end_hg19) + WIN,
             studies = paste(sort(unique(unlist(strsplit(direct_studies, ";")))), collapse = ";"),
             release_direct_coloc_genes = paste(sort(unique(unlist(strsplit(direct_coloc_genes[!is.na(direct_coloc_genes)], ";")))), collapse = ";")),
         by = physical_locus]

reg <- fread(REG)
gtf <- import(GTF, feature.type = "gene"); gtf$gene_id <- sub("\\..*$", "", gtf$gene_id)
sym <- setNames(gtf$gene_name, gtf$gene_id)
chain <- import.chain(CHAIN)

read_gwas <- function(study, chr, s, e) {
  f <- file.path(SUMS, paste0(study, "_af.tsv"))
  d <- fread(cmd = sprintf("awk -F'\\t' 'NR==1 || ($1==\"%s\" && $2>=%d && $2<=%d)' %s", chr, s, e, shQuote(f)))
  d[!is.na(beta) & !is.na(se) & se > 0 & !is.na(pval)]
}
lift38 <- function(chr, pos) {
  lo <- liftOver(GRanges(paste0("chr", chr), IRanges(pos, width = 1)), chain)
  ok <- lengths(lo) == 1L; out <- rep(NA_integer_, length(pos)); out[ok] <- start(unlist(lo[ok])); out
}
akey <- function(a, b) ifelse(a < b, paste(a, b, sep = "_"), paste(b, a, sep = "_"))

# ---- network check ----
net_ok <- identical(suppressWarnings(system2("curl", c("-s", "-o", "/dev/null", "-w", "%{http_code}", "-I", "--max-time", "60",
                                                      shQuote(paste0(DS$url[1], ".tbi"))), stdout = TRUE)), "200")
cat("eQTL Catalogue reachable:", net_ok, "\n")

fetch_qtl <- function(ds, chr, s, e, tag) {
  f <- file.path(cache, sprintf("%s__%s.tsv.gz", ds$qtd, tag))
  if (!file.exists(f)) {
    old <- setwd(tbi_dir); on.exit(setwd(old))
    tmp <- paste0(f, ".part")
    rc <- system2("bash", c("-c", shQuote(sprintf("set -o pipefail; %s %s %s:%d-%d | gzip -c > %s",
                                                  TABIX, ds$url, chr, s, e, shQuote(tmp)))))
    if (rc != 0) return(NULL)
    file.rename(tmp, f)
  }
  if (file.size(f) < 30) return(data.table())
  fread(cmd = paste("zcat", shQuote(f)), header = FALSE, col.names = QCOLS, sep = "\t", na.strings = c("NA", ""))
}

coloc_one <- function(g, q, gw_type, gw_N, gw_s, q_N) {
  m <- merge(g, q, by = "k")
  m <- m[!duplicated(k)]
  if (nrow(m) < MIN_SNPS) return(NULL)
  d1 <- list(beta = m$beta.g, varbeta = m$se.g^2, N = gw_N, type = gw_type, snp = m$k)
  if (gw_type == "cc") d1$s <- gw_s else d1$sdY <- 1
  d2 <- list(beta = m$beta.q, varbeta = m$se.q^2, N = q_N, type = "quant", sdY = 1, snp = m$k)
  invisible(capture.output(r <- suppressWarnings(suppressMessages(coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12)))))
  top <- r$results[which.max(r$results$SNP.PP.H4), ]
  data.table(nsnps = nrow(m), PP.H0 = r$summary[["PP.H0.abf"]], PP.H1 = r$summary[["PP.H1.abf"]],
             PP.H2 = r$summary[["PP.H2.abf"]], PP.H3 = r$summary[["PP.H3.abf"]], PP.H4 = r$summary[["PP.H4.abf"]],
             top_snp_hg38 = top$snp, top_snp_pp_h4 = top$SNP.PP.H4)
}

choice <- list(); status <- list(); results <- list(); l4map <- list()
vpk <- fread(file.path(OUT, "inventory_v2/variant_peaks.tsv"))
for (i in seq_len(nrow(LOC))) {
  lc <- LOC[i]; tag <- gsub("[^A-Za-z0-9]+", "_", lc$physical_locus)
  cand <- intersect(strsplit(lc$studies, ";")[[1]], sub("_af\\.tsv$", "", list.files(SUMS, "_af\\.tsv$")))
  best <- NULL; best_p <- Inf; best_d <- NULL
  for (st in sort(cand)) {
    d <- read_gwas(st, lc$chrom, lc$w19_start, lc$w19_end)
    if (nrow(d) && min(d$pval) < best_p) { best_p <- min(d$pval); best <- st; best_d <- d }
  }
  rg <- if (!is.null(best)) reg[study_name == best] else NULL
  choice[[i]] <- data.table(physical_locus = lc$physical_locus, locus_set = lc$locus_set, candidate_studies = paste(cand, collapse = ";"),
                            gwas_study = if (is.null(best)) NA_character_ else best, gwas_min_p = best_p,
                            gwas_ancestry = if (is.null(rg)) NA_character_ else rg$ancestry,
                            gwas_type = if (is.null(rg)) NA_character_ else rg$trait_type)
  layers <- if (lc$locus_set == "CTRL_L0") DS[layer == "L0"] else DS
  if (is.null(best) || best_p > P_GWAS_MAX) {
    status[[i]] <- data.table(physical_locus = lc$physical_locus, layer = c(layers$layer, if (lc$locus_set != "CTRL_L0") "L4"),
                              dataset = c(layers$qtd, if (lc$locus_set != "CTRL_L0") "Currin2025"),
                              status = "untestable_gwas_p", n_traits = 0L)
    next
  }
  gw_type <- if (rg$trait_type == "binary") "cc" else "quant"
  gw_s <- if (gw_type == "cc") rg$N_cases / rg$N_tot else NA_real_
  best_d[, pos38 := lift38(chromosome, position)]
  g <- best_d[!is.na(pos38), .(k = paste(pos38, akey(allele1, allele2), sep = ":"), beta.g = beta, se.g = se, p.g = pval)]
  s38 <- min(best_d$pos38, na.rm = TRUE); e38 <- max(best_d$pos38, na.rm = TRUE)
  st_rows <- list()
  for (j in seq_len(nrow(layers))) {
    ds <- layers[j]
    if (!net_ok) { st_rows[[j]] <- data.table(physical_locus = lc$physical_locus, layer = ds$layer, dataset = ds$qtd, status = "not_run_no_network", n_traits = 0L); next }
    q <- fetch_qtl(ds, lc$chrom, s38, e38, tag)
    if (is.null(q)) { st_rows[[j]] <- data.table(physical_locus = lc$physical_locus, layer = ds$layer, dataset = ds$qtd, status = "fetch_failed", n_traits = 0L); next }
    if (!nrow(q)) { st_rows[[j]] <- data.table(physical_locus = lc$physical_locus, layer = ds$layer, dataset = ds$qtd, status = "no_rows_in_window", n_traits = 0L); next }
    q <- q[!is.na(beta) & !is.na(se) & se > 0]
    q[, k := paste(position, akey(ref, alt), sep = ":")]
    tr <- q[, .(minp = min(pvalue)), by = molecular_trait_id][minp < P_QTL_MAX]
    out <- list()
    for (t in tr$molecular_trait_id) {
      qt <- q[molecular_trait_id == t, .(k, beta.q = beta, se.q = se)]
      cr <- coloc_one(g, qt, gw_type, rg$N_tot, gw_s, max(q[molecular_trait_id == t]$an, na.rm = TRUE) / 2)
      if (is.null(cr)) next
      gid <- q[molecular_trait_id == t]$gene_id[1]
      out[[t]] <- cbind(data.table(locus_set = lc$locus_set, physical_locus = lc$physical_locus, layer = ds$layer, dataset = ds$qtd,
                                   dataset_label = ds$label, proxy_lineage = ds$proxy_lineage, gwas_study = best,
                                   molecular_trait_id = t, gene_id = gid, symbol = unname(sym[gid]),
                                   min_qtl_p = tr[molecular_trait_id == t]$minp), cr)
    }
    results[[paste(i, j)]] <- rbindlist(out)
    st_rows[[j]] <- data.table(physical_locus = lc$physical_locus, layer = ds$layer, dataset = ds$qtd, status = "tested",
                               n_traits = length(out))
  }
  # L4: Currin caQTL at variant peaks (H24 and REF_PNPLA3 loci)
  if (lc$locus_set != "CTRL_L0") {
    f <- list.files(CURRIN, sprintf("^caQTL_hg38_chr%s_with", lc$chrom), full.names = TRUE)
    cf <- file.path(cache, sprintf("Currin2025__%s.tsv.gz", tag))
    if (!file.exists(cf)) {
      cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==1{print;next}{split($1,a,\":\"); if (a[2]>=%d && a[2]<=%d) print}' | gzip -c > %s",
                     shQuote(f), s38, e38, shQuote(paste0(cf, ".part")))
      rc <- system2("bash", c("-c", shQuote(paste("set -o pipefail;", cmd))))
      if (rc != 0) stop("Currin extraction failed for ", lc$physical_locus)
      file.rename(paste0(cf, ".part"), cf)
    }
    cq <- fread(cmd = paste("zcat", shQuote(cf)))
    cq[, c("chr", "pos", "ref", "alt") := tstrsplit(variant, ":")]; cq[, pos := as.integer(pos)]
    # orientation of distance_from_peakCenter: choose the convention that is constant within a peak
    v1 <- cq[, var(pos - distance_from_peakCenter), by = peak]$V1; v2 <- cq[, var(pos + distance_from_peakCenter), by = peak]$V1
    sgn <- if (sum(v1, na.rm = TRUE) <= sum(v2, na.rm = TRUE)) -1L else 1L
    cent <- cq[, .(centre = as.integer(median(pos + sgn * distance_from_peakCenter)), minp = min(pvalue)), by = peak]
    vp <- vpk[physical_locus == lc$physical_locus]
    if (nrow(vp)) {
      pc <- tstrsplit(vp$peak_id, "-")
      vp[, `:=`(ps = as.integer(pc[[2]]), pe = as.integer(pc[[3]]))]
    }
    mp <- rbindlist(lapply(seq_len(nrow(vp)), function(r) {
      h <- cent[centre >= vp$ps[r] - PEAK_PAD & centre <= vp$pe[r] + PEAK_PAD]
      if (nrow(h)) data.table(peak_id = vp$peak_id[r], currin_peak = h$peak, currin_centre = h$centre, currin_min_p = h$minp) else NULL
    }))
    out <- list()
    if (nrow(mp) > 0) {
      l4map[[i]] <- cbind(physical_locus = lc$physical_locus, mp)
      cq[, k := paste(pos, akey(ref, alt), sep = ":")]
      for (r in seq_len(nrow(mp))) {
        if (mp$currin_min_p[r] >= P_QTL_MAX) next
        qt <- cq[peak == mp$currin_peak[r], .(k, beta.q = beta, se.q = sqrt(varbeta))]
        cr <- coloc_one(g, qt, gw_type, rg$N_tot, gw_s, CURRIN_N)
        if (is.null(cr)) next
        out[[r]] <- cbind(data.table(locus_set = lc$locus_set, physical_locus = lc$physical_locus, layer = "L4", dataset = "Currin2025",
                                     dataset_label = "Currin 2025 liver caQTL", proxy_lineage = "bulk_liver", gwas_study = best,
                                     molecular_trait_id = mp$currin_peak[r], gene_id = NA_character_, symbol = NA_character_,
                                     min_qtl_p = mp$currin_min_p[r], gse296875_peak_id = mp$peak_id[r]), cr)
      }
    }
    results[[paste(i, "L4")]] <- rbindlist(out, fill = TRUE)
    st_rows[[length(st_rows) + 1L]] <- data.table(physical_locus = lc$physical_locus, layer = "L4", dataset = "Currin2025",
                                                  status = if (!nrow(vp)) "no_variant_peak" else if (!nrow(mp)) "no_caqtl_peak_at_variant_peak" else "tested",
                                                  n_traits = length(out))
  }
  status[[i]] <- rbindlist(st_rows)
  cat(sprintf("[%s] %s GWAS %s min p %.2e\n", lc$locus_set, lc$physical_locus, best, best_p))
}
choice <- rbindlist(choice); status <- rbindlist(status, fill = TRUE)
results <- rbindlist(results, fill = TRUE); l4map <- rbindlist(l4map, fill = TRUE)
fwrite(LOC, file.path(res_dir, "coloc_loci.tsv"), sep = "\t")
fwrite(choice, file.path(res_dir, "coloc_gwas_choice.tsv"), sep = "\t")
fwrite(status, file.path(res_dir, "coloc_status.tsv"), sep = "\t")
fwrite(results, file.path(res_dir, "coloc_results.tsv"), sep = "\t")
fwrite(l4map, file.path(res_dir, "l4_caqtl_peak_map.tsv"), sep = "\t")
write_json(list(network_ok = net_ok, priors = c(p1 = P1, p2 = P2, p12 = P12), loci = nrow(LOC),
                coloc_rows = nrow(results)), file.path(res_dir, "coloc_summary.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_coloc.txt"))
cat("coloc done\n")
