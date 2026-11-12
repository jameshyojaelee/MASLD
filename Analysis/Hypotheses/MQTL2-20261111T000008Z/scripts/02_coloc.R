#!/usr/bin/env Rscript
# MQTL2 step 02: coloc.abf of Pan-UKBB EUR ALT (GRCh37) against Chen 2023 metabolite and Ottensmann 2023
# lipid-species GWAS (GRCh38, lifted per SNP to GRCh37) at the 169 ALT loci + PNPLA3 control window.
# Adapted from GWAS/finemapping/src/60b_mqtl_coloc.R (same SNP harmonisation, MIN_SNPS, liftover chains).
# Changes from 60b: ALT GWAS instead of MASLD GWAS; coloc default priors (p12 = 1e-5) as primary and the 60b prior
# (p12 = 5e-6) as sensitivity; metabolite files streamed through a window filter (no whole-file reads); one job
# with forked workers instead of an array; build / effect-allele checks at control SNPs; ALT self-coloc and MVP ALT
# control at GCKR. Prespecification: ../PRESPEC_MQTL2.md.
# Usage: Rscript 02_coloc.R <out_dir> <loci_file> <n_cores>
suppressMessages({
  library(data.table); library(coloc); library(parallel)
  library(rtracklayer); library(GenomicRanges)
})
set.seed(20261110)
setDTthreads(1)  # forked workers; avoid thread oversubscription
args <- commandArgs(TRUE)
OUT <- args[1]; LOCI_FILE <- args[2]; NCORES <- as.integer(args[3])
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/MQTL2-20261111T000008Z"
TESTN <- as.integer(Sys.getenv("MQTL2_TEST_N", "0"))  # >0 only for a pre-submission smoke test
TABIX <- "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/tabix"
EXTRACT <- file.path(RUN, "scripts/extract_windows.sh")
ALT_FILE <- file.path(ROOT, "GWAS/MR_Data/PanUKBB/biomarkers-30620-ALT-both_sexes-irnt.tsv.bgz")
MVP_FILE <- file.path(ROOT, "GWAS/finemapping/data/sumstats/MVP_ALT_EUR_reformatted_hg19.tsv")
ALT_N <- 400822L; MVP_N <- 406803L
MIN_SNPS <- 50L
P1 <- 1e-4; P2 <- 1e-4; P12 <- 1e-5; P12_ALT <- 5e-6
CTRL19 <- data.table(snp = c("rs1260326", "rs738409"), chr = c(2L, 22L), pos19 = c(27730940L, 44324727L),
                     pos38 = c(27508073L, 43928847L))
dir.create(file.path(OUT, "windows"), recursive = TRUE)

cat("MQTL2 02_coloc | coloc", as.character(packageVersion("coloc")), "| cores", NCORES, "|", date(), "\n")
loci <- fread(LOCI_FILE, colClasses = list(character = "locus_id"))
cat("loci:", nrow(loci), "\n")

# ---------------- ALT windows (tabix; GRCh37; effect allele = alt) ----------------
read_alt <- function(chr, s, e) {
  cmd <- sprintf("%s %s %s:%d-%d | cut -f1-4,19,25,31,37,43", TABIX, ALT_FILE, chr, s, e)
  d <- fread(cmd = cmd, header = FALSE, sep = "\t", na.strings = "NA",
             col.names = c("chr", "pos", "ref", "alt", "af", "beta", "se", "nlp", "lowconf"))
  d <- d[lowconf %in% c("false", FALSE) & !is.na(beta) & !is.na(se) & se > 0]
  d[, lowconf := NULL]
  d
}
alt <- lapply(seq_len(nrow(loci)), function(i) read_alt(loci$chromosome[i], loci$win_start[i], loci$win_end[i]))
names(alt) <- loci$locus_id
alt_sum <- rbindlist(lapply(seq_len(nrow(loci)), function(i) {
  d <- alt[[i]]; j <- which.max(d$nlp)
  data.table(locus_id = loci$locus_id[i], n_alt_snps = nrow(d), alt_max_nlp = d$nlp[j], alt_top_pos = d$pos[j])
}))
fwrite(alt_sum, file.path(OUT, "alt_window_summary.tsv"), sep = "\t")
cat("ALT windows read; median SNPs", median(alt_sum$n_alt_snps), "; loci with ALT max -log10p < 7.3:",
    sum(alt_sum$alt_max_nlp < 7.3), "\n")
ctrl_alt <- rbindlist(lapply(seq_len(nrow(CTRL19)), function(i)
  read_alt(CTRL19$chr[i], CTRL19$pos19[i], CTRL19$pos19[i])[, snp := CTRL19$snp[i]]))
print(ctrl_alt)
stopifnot(nrow(ctrl_alt[snp == "rs1260326" & ref == "T" & alt == "C"]) == 1,
          nrow(ctrl_alt[snp == "rs738409" & ref == "C" & alt == "G"]) == 1)

# ---------------- MVP ALT at GCKR (control only; GRCh37; effect = allele1) ----------------
gk <- loci[locus_id == "17"]
mvp_gckr <- fread(cmd = sprintf("awk -F'\\t' 'NR==1 || ($1==%d && $2>=%d && $2<=%d)' %s",
                                gk$chromosome, gk$win_start, gk$win_end, MVP_FILE))
mvp_gckr <- mvp_gckr[!is.na(beta) & !is.na(se) & se > 0]
cat("MVP ALT GCKR window SNPs:", nrow(mvp_gckr), "| min p", min(mvp_gckr$pval), "\n")

# ---------------- hg38 windows for streaming filter ----------------
chain_fwd <- import.chain(file.path(ROOT, "data/broadaway_eqtl/hg19ToHg38.over.chain"))
chain_back <- import.chain(file.path(ROOT, "data/broadaway_eqtl/hg38ToHg19.over.chain"))
w38 <- rbindlist(lapply(seq_len(nrow(loci)), function(i) {
  g <- GRanges(paste0("chr", loci$chromosome[i]), IRanges(loci$win_start[i], loci$win_end[i]))
  l <- unlist(liftOver(g, chain_fwd))
  l <- l[as.character(seqnames(l)) == paste0("chr", loci$chromosome[i])]
  data.table(chr = loci$chromosome[i], start = max(1, min(start(l)) - 5e4), end = max(end(l)) + 5e4)
}))
WIN_FILE <- file.path(OUT, "hg38_windows.tsv")
fwrite(w38, WIN_FILE, sep = "\t", col.names = FALSE)

# ---------------- metabolite file list ----------------
DE <- file.path(ROOT, "data/external")
chen <- fread(file.path(DE, "chen2023_mqtl/chen2023_biomarker_download_list.tsv"))
chen[, `:=`(source = "chen", file = file.path(DE, "chen2023_mqtl/sumstats", paste0(accession, ".tsv.gz")), N = Nnum)]
man <- fread(file.path(RUN, "downloads_manifest.tsv"))
amap <- fread(file.path(DE, "chen2023_mqtl/chen2023_accession_trait_map.tsv"))
amap[, Nnum := as.numeric(sub(" .*", "", gsub(",", "", sample)))]
new <- merge(man[, .(accession)], amap[, .(accession, trait, Nnum)], by = "accession")
if (TESTN == 0) stopifnot(nrow(new) == 5)
new[, `:=`(source = "chen_new", file = file.path(RUN, "downloads", paste0(accession, "_buildGRCh38.tsv.gz")), N = Nnum)]
ott <- fread(file.path(DE, "lipidqtl/ottensmann2023_accession_trait_map.tsv"))
ott[, `:=`(source = "ottensmann", file = file.path(DE, "lipidqtl/sumstats", paste0(accession, ".tsv.gz")),
           N = as.numeric(sub(" .*", "", gsub(",", "", sample))))]
mf <- rbind(chen[, .(source, accession, trait, N, file)], new[, .(source, accession, trait, N, file)],
            ott[, .(source, accession, trait, N, file)])
mf <- mf[file.exists(file)]
if (TESTN > 0) mf <- rbind(head(mf[source == "chen"], TESTN), head(mf[source == "ottensmann"], TESTN))
cat("metabolite files present:", nrow(mf), "\n"); print(mf[, .N, by = source])
fwrite(mf, file.path(OUT, "metabolite_files.tsv"), sep = "\t")

# ---------------- per-file worker ----------------
harmonise_cols <- function(d) {
  setnames(d, tolower(names(d)))
  ren <- function(from, to) for (x in from) if (x %in% names(d) && !(to %in% names(d))) setnames(d, x, to)
  ren(c("chromosome"), "chr"); ren(c("base_pair_location"), "pos")
  ren(c("effect_allele"), "ea"); ren(c("other_allele"), "oa")
  ren(c("standard_error"), "se"); ren(c("effect_allele_frequency"), "eaf")
  if (!"rsid" %in% names(d)) d[, rsid := NA_character_]
  if (!"eaf" %in% names(d)) d[, eaf := NA_real_]
  d[, chr := as.integer(sub("^chr", "", chr))]
  d[, `:=`(ea = toupper(ea), oa = toupper(oa))]
  d[!is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(se) & se > 0, .(chr, pos, ea, oa, eaf, beta, se, rsid)]
}

run_pair <- function(g, q, n_q, n_g, snp_prefix) {
  # g: ALT side with pos, a_eff, a_oth, beta, se; q: metabolite side (GRCh37) with pos, ea, oa, beta, se
  m <- merge(q[, .(pos, q_ea = ea, q_oa = oa, q_beta = beta, q_se = se)],
             g[, .(pos, g_eff = a_eff, g_oth = a_oth, g_beta = beta, g_se = se)], by = "pos")
  m <- m[(q_ea == g_eff & q_oa == g_oth) | (q_ea == g_oth & q_oa == g_eff)]
  m[q_ea == g_oth, q_beta := -q_beta]
  m <- m[!((g_eff %in% c("A", "T") & g_oth %in% c("A", "T")) | (g_eff %in% c("C", "G") & g_oth %in% c("C", "G")))]
  m <- m[!(duplicated(pos) | duplicated(pos, fromLast = TRUE))]
  if (nrow(m) < MIN_SNPS) return(NULL)
  snp <- paste(snp_prefix, m$pos, sep = ":")
  d1 <- list(beta = m$g_beta, varbeta = m$g_se^2, N = n_g, type = "quant", sdY = 1, snp = snp)
  d2 <- list(beta = m$q_beta, varbeta = m$q_se^2, N = n_q, type = "quant", sdY = 1, snp = snp)
  a <- tryCatch(suppressMessages(suppressWarnings(coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12))),
                error = function(e) NULL)
  b <- tryCatch(suppressMessages(suppressWarnings(coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12_ALT))),
                error = function(e) NULL)
  if (is.null(a)) return(NULL)
  pp <- a$summary[paste0("PP.H", 0:4, ".abf")]
  ti <- which.max(a$results$SNP.PP.H4)
  tp <- as.integer(sub(".*:", "", a$results$snp[ti])); mi <- match(tp, m$pos)
  data.table(n_snps = nrow(m), PP.H0 = pp[[1]], PP.H1 = pp[[2]], PP.H2 = pp[[3]], PP.H3 = pp[[4]], PP.H4 = pp[[5]],
             PP.H3_p12_5e6 = if (is.null(b)) NA_real_ else b$summary[["PP.H3.abf"]],
             PP.H4_p12_5e6 = if (is.null(b)) NA_real_ else b$summary[["PP.H4.abf"]],
             top_pos = tp, top_snp_pp = a$results$SNP.PP.H4[ti],
             top_eff = m$g_eff[mi], top_alt_beta = m$g_beta[mi], top_alt_se = m$g_se[mi],
             top_met_beta = m$q_beta[mi], top_met_se = m$q_se[mi])
}

alt_g <- lapply(alt, function(d) d[, .(pos, a_eff = alt, a_oth = ref, beta, se)])

# build check of the ALT locus coordinates: DIG1 lead SNPs must sit in Pan-UKBB (GRCh37) with the DIG1 alleles
lead_ok <- sapply(seq_len(nrow(loci)), function(i) {
  d <- alt[[i]][pos == loci$cs_lead_pos[i]]
  nrow(d) > 0 && any(mapply(function(r, a) setequal(c(r, a), c(loci$lead_a1[i], loci$lead_a2[i])), d$ref, d$alt))
})
cat("DIG1 lead SNPs found in Pan-UKBB GRCh37 with matching alleles:", sum(lead_ok), "of", nrow(loci), "\n")
fwrite(data.table(locus_id = loci$locus_id, lead_found_grch37 = lead_ok), file.path(OUT, "lead_build_check.tsv"), sep = "\t")
stopifnot(mean(lead_ok) >= 0.9)

do_file <- function(k) {
  r <- mf[k]
  shaf <- file.path(OUT, "windows", paste0(r$accession, ".sha256"))
  wf <- file.path(OUT, "windows", paste0(r$accession, ".win.tsv.gz"))
  d <- tryCatch(fread(cmd = sprintf("bash %s %s %s %s", EXTRACT, r$file, WIN_FILE, shaf), sep = "\t",
                      na.strings = c("NA", ""), colClasses = list(character = c("chromosome"))),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(list(chk = data.table(accession = r$accession, status = "unreadable")))
  d <- harmonise_cols(d)
  # build check (raw GRCh38 rows at control SNPs)
  b38 <- rbindlist(lapply(seq_len(nrow(CTRL19)), function(i) {
    x <- d[chr == CTRL19$chr[i] & pos == CTRL19$pos38[i]]
    y <- d[!is.na(rsid) & rsid == CTRL19$snp[i]]
    data.table(snp = CTRL19$snp[i], n_at_pos38 = nrow(x),
               alleles_at_pos38 = paste(unique(paste0(x$ea, "/", x$oa)), collapse = ";"),
               rsid_pos = if (nrow(y)) paste(unique(paste(y$chr, y$pos, sep = ":")), collapse = ";") else NA_character_)
  }))
  # liftover hg38 -> hg19, 1:1
  gr <- GRanges(paste0("chr", d$chr), IRanges(d$pos, width = 1))
  lo <- liftOver(gr, chain_back)
  keep <- which(lengths(lo) == 1)
  gl <- unlist(lo[keep])
  q <- d[keep][, `:=`(chr = as.integer(sub("chr", "", as.character(seqnames(gl)))), pos = start(gl))]
  q <- q[chr == d$chr[keep]]
  fwrite(q, wf, sep = "\t")
  # allele and effect-allele checks at GRCh37 control SNPs
  chk <- rbindlist(lapply(seq_len(nrow(CTRL19)), function(i) {
    x <- q[chr == CTRL19$chr[i] & pos == CTRL19$pos19[i]]
    a <- ctrl_alt[snp == CTRL19$snp[i]]
    x <- x[(ea == a$alt & oa == a$ref) | (ea == a$ref & oa == a$alt)]
    eaf_alt <- if (nrow(x) == 1 && !is.na(x$eaf)) ifelse(x$ea == a$alt, x$eaf, 1 - x$eaf) else NA_real_
    other <- q[chr == CTRL19$chr[i] & pos == CTRL19$pos19[i] & !((ea == a$alt & oa == a$ref) | (ea == a$ref & oa == a$alt))]
    data.table(snp = CTRL19$snp[i], n_match = nrow(x), n_other_alleles_at_pos = nrow(other),
               met_ea = if (nrow(x) == 1) x$ea else NA_character_,
               eaf_of_panukbb_alt = eaf_alt, panukbb_af_alt = a$af,
               met_beta_per_panukbb_alt = if (nrow(x) == 1) ifelse(x$ea == a$alt, x$beta, -x$beta) else NA_real_,
               met_z = if (nrow(x) == 1) x$beta / x$se else NA_real_)
  }))
  chk <- cbind(data.table(accession = r$accession, status = "ok", n_window_rows = nrow(d), n_lifted = nrow(q)),
               merge(chk, b38, by = "snp"))
  # gate before any colocalization: alleles at GCKR rs1260326 and PNPLA3 rs738409 (GRCh37 after liftover)
  # must match Pan-UKBB {ref, alt}; at least one control SNP must be found (confirms GRCh38 input -> GRCh37).
  if (any(chk$n_other_alleles_at_pos > 0 & chk$n_match == 0)) { chk[, status := "allele_fail"]; return(list(chk = chk)) }
  if (all(chk$n_match == 0)) { chk[, status := "build_unconfirmed"]; return(list(chk = chk)) }
  # coloc per locus
  res <- list()
  for (i in seq_len(nrow(loci))) {
    L <- loci[i]
    qs <- q[chr == L$chromosome & pos >= L$win_start & pos <= L$win_end]
    if (nrow(qs) < MIN_SNPS) next
    o <- run_pair(alt_g[[L$locus_id]], qs, r$N, ALT_N, L$chromosome)
    if (!is.null(o)) res[[length(res) + 1]] <- cbind(data.table(gwas = "PanUKBB_ALT_EUR", locus_id = L$locus_id), o)
    if (L$locus_id == "17") {
      g <- mvp_gckr[, .(pos = position, a_eff = toupper(allele1), a_oth = toupper(allele2), beta, se)]
      o <- run_pair(g, qs, r$N, MVP_N, L$chromosome)
      if (!is.null(o)) res[[length(res) + 1]] <- cbind(data.table(gwas = "MVP_ALT_EUR_control", locus_id = "17"), o)
    }
  }
  out <- if (length(res)) cbind(data.table(source = r$source, accession = r$accession, trait = r$trait,
                                           N_met = r$N), rbindlist(res)) else NULL
  list(chk = chk, res = out)
}

t0 <- Sys.time()
R <- mclapply(seq_len(nrow(mf)), function(k) tryCatch(do_file(k), error = function(e)
  list(chk = data.table(accession = mf$accession[k], status = paste("error:", conditionMessage(e))))),
  mc.cores = NCORES, mc.preschedule = FALSE)
cat("worker time:", format(Sys.time() - t0), "\n")
chk <- rbindlist(lapply(R, `[[`, "chk"), fill = TRUE)
res <- rbindlist(lapply(R, `[[`, "res"), fill = TRUE)
fwrite(chk, file.path(OUT, "file_checks.tsv"), sep = "\t")
fwrite(res, file.path(OUT, "coloc_all_tests.tsv.gz"), sep = "\t")
cat("file statuses:\n"); print(chk[, .N, by = status])
if (any(chk$status == "allele_fail")) { print(chk[status == "allele_fail"]); stop("allele alignment failed at a control SNP") }
if (nrow(res) == 0 || mean(chk$status != "ok") > 0.1) { print(chk[status != "ok"]); stop("too many failed files or no tests") }
cat("tests:", nrow(res), "| PanUKBB tests:", nrow(res[gwas == "PanUKBB_ALT_EUR"]),
    "| loci tested:", uniqueN(res[gwas == "PanUKBB_ALT_EUR"]$locus_id), "\n")

# sha256 of metabolite files (computed during the streaming read)
sha <- rbindlist(lapply(seq_len(nrow(mf)), function(k) {
  f <- file.path(OUT, "windows", paste0(mf$accession[k], ".sha256"))
  data.table(file = mf$file[k], sha256 = if (file.exists(f)) readLines(f)[1] else NA_character_)
}))
fwrite(sha, file.path(OUT, "metabolite_input_sha256.tsv"), sep = "\t")

# ---------------- self-colocalization: ALT vs ALT at locus 148 ----------------
g <- alt_g[["148"]]
g <- g[!((a_eff %in% c("A", "T") & a_oth %in% c("A", "T")) | (a_eff %in% c("C", "G") & a_oth %in% c("C", "G")))]
g <- g[!(duplicated(pos) | duplicated(pos, fromLast = TRUE))]
d1 <- list(beta = g$beta, varbeta = g$se^2, N = ALT_N, type = "quant", sdY = 1, snp = as.character(g$pos))
self <- suppressMessages(coloc.abf(d1, d1, p1 = P1, p2 = P2, p12 = P12))$summary
print(self)
fwrite(as.data.table(as.list(self)), file.path(OUT, "control_self_coloc_locus148.tsv"), sep = "\t")
cat("02 done |", date(), "\n")
