#!/usr/bin/env Rscript
# H36 step 1: junction filtering and clustering (coverage-only filters; no stage
# labels are used). Implements PRESPEC_H36.md section 3 and Addendum 1.
# Usage: Rscript h36_01_build_clusters.R <repo> <out_dir>
suppressPackageStartupMessages({ library(data.table); library(parallel) })
args <- commandArgs(trailingOnly = TRUE); REPO <- args[1]; OUT <- args[2]
set.seed(20260923)
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
CAND <- file.path(REPO, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z")
READS <- file.path(REPO, "Analysis/Hypotheses/READS-20261009T141017Z/libraries")
meta <- as.data.table(readRDS(file.path(CAND, "primary_source_controls/meta_matched.rds")))
meta <- meta[, .(sample_id, dataset, fibrosis_stage, inferred_sex, analysis_unit_id)]
meta[, pi := .I]
stopifnot(nrow(meta) == 774L, !anyDuplicated(meta$analysis_unit_id))
fwrite(meta, file.path(OUT, "persons.tsv"), sep = "\t")

CHR <- c(paste0("chr", 1:22), "chrX")
STR <- c("+", "-", ".")
read_one <- function(i) {
  f <- file.path(READS, meta$dataset[i], meta$sample_id[i], "junctions.tsv.gz")
  d <- fread(f, select = c("chrom", "intron_start", "intron_end", "strand", "annotated", "n_unique"))
  d <- d[n_unique > 0L & chrom %in% CHR]
  d[, len := intron_end - intron_start + 1L]
  d <- d[len >= 50L & len <= 100000L]
  data.table(pi = i, chr = match(d$chrom, CHR), st = d$intron_start, en = d$intron_end,
             str = match(d$strand, STR), ann = d$annotated, n = d$n_unique)
}
t0 <- Sys.time()
L <- rbindlist(mclapply(seq_len(nrow(meta)), read_one, mc.cores = NC))
cat("read", nrow(L), "junction rows from", nrow(meta), "libraries in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
L[, jn := .GRP, by = .(chr, st, en, str)]
jt <- L[, .(chr = chr[1], st = st[1], en = en[1], str = str[1], ann = max(ann),
            reads = sum(n), npers = .N), by = jn]
cat("distinct junctions:", nrow(jt), "\n")
jt <- jt[npers >= 5L]
cat("junctions in >=5 persons:", nrow(jt), "\n")
L <- L[jn %in% jt$jn]
setkey(jt, jn)

# ---- clusters: shared donor or acceptor, same chromosome and strand ----------
jt[, don := .GRP, by = .(chr, str, st)]
jt[, acc := .GRP, by = .(chr, str, en)]
lab <- seq_len(nrow(jt))
repeat {
  l1 <- data.table(x = lab, g = jt$don)[, m := min(x), by = g]$m
  l2 <- data.table(x = l1, g = jt$acc)[, m := min(x), by = g]$m
  if (identical(l2, lab)) break
  lab <- l2
}
jt[, cl0 := lab]
L[, cl0 := jt$cl0[match(jn, jt$jn)]]

# ---- person inclusion, junction filters, final inclusion (Addendum 1) ---------
pt <- L[, .(tot = sum(n)), by = .(cl0, pi)]
inc1 <- pt[tot >= 20L, .(cl0, pi)]
cl_all <- L[, .(cl_reads = sum(n)), by = cl0]
jt <- merge(jt, cl_all, by = "cl0")
jt[, frac := reads / cl_reads]
Li <- merge(L, inc1, by = c("cl0", "pi"))
jinc <- Li[, .(n_inc_pers = .N), by = jn]
jt <- merge(jt, jinc, by = "jn", all.x = TRUE)
jt[is.na(n_inc_pers), n_inc_pers := 0L]
keep <- jt[frac >= 0.001 & n_inc_pers >= 10L]
setorder(keep, cl0, -reads)
keep[, jr := seq_len(.N), by = cl0]
keep <- keep[jr <= 10L]
keep[, J := .N, by = cl0]
keep <- keep[J >= 2L]
cat("clusters with >=2 kept junctions:", uniqueN(keep$cl0), "| kept junctions:", nrow(keep), "\n")
Lk <- merge(L, keep[, .(jn, cl0, jr)], by = c("jn", "cl0"))
pt2 <- Lk[, .(tot = sum(n)), by = .(cl0, pi)]
inc2 <- pt2[tot >= 20L, .(cl0, pi)]
Lk <- merge(Lk, inc2, by = c("cl0", "pi"))
# staged-person coverage conditions of the test (section 3.4)
mm <- merge(unique(Lk[, .(cl0, pi)]), meta[, .(pi, dataset, stage = fibrosis_stage)], by = "pi")
cs <- mm[!is.na(stage), .(n_staged = .N, n_lo = sum(stage <= 1), n_hi = sum(stage >= 3)), by = cl0]
cs <- cs[n_staged >= 150L & n_lo >= 30L & n_hi >= 30L]
cat("clusters meeting staged-person coverage:", nrow(cs), "\n")
keep <- keep[cl0 %in% cs$cl0]
keep[, cl := match(cl0, sort(unique(cl0)))]
Lk <- merge(Lk[cl0 %in% cs$cl0], keep[, .(cl0, jn, cl)], by = c("cl0", "jn"))
# re-rank junction order by reads within each cluster (jr: 1 = most reads)
long <- Lk[, .(cl, jr, pi, n)]
setkey(long, cl, pi, jr)
saveRDS(long, file.path(OUT, "long.rds"))
cl_tab <- keep[, .(J = .N, chr = CHR[chr[1]], str = STR[str[1]], start_min = min(st), end_max = max(en),
                   cl_reads = cl_reads[1]), by = .(cl, cl0)]
cl_tab <- merge(cl_tab, cs, by = "cl0")
fwrite(cl_tab, file.path(OUT, "clusters.tsv.gz"), sep = "\t")
fwrite(keep[, .(cl, jr, chr = CHR[chr], start = st, end = en, strand = STR[str], annotated = ann,
                reads, frac, n_inc_pers)], file.path(OUT, "junctions.tsv.gz"), sep = "\t")
# global kept-junction list (all clusters with >=5 persons) for coordinate matching
fwrite(jt[, .(chr = CHR[chr], start = st, end = en, strand = STR[str], annotated = ann, reads, npers)],
       file.path(OUT, "junctions_all_ge5_persons.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_01.txt"))
cat("clusters tested later:", nrow(cl_tab), "| long rows:", nrow(long), "\n")
cat("done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
