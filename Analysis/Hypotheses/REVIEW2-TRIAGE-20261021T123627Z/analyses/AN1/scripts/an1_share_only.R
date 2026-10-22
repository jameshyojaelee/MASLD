#!/usr/bin/env Rscript
# AN1 (REVIEW2 B-G1): share-only symmetric shared-posterior check.
#
# For each (study x gene) row, rebuild the GWAS-fit variant set exactly as
# GWAS/finemapping/src/06_susie_coloc.R builds it (allele harmonization lines
# 395-604 extracted verbatim from that file at run time; eQTL-fit restriction
# line 670; LD-panel matching via get_ld_per_locus), then compute each eQTL
# credible set's posterior share on that set with coloc 5.2.3's convention
# (signal_shared_prop() from src/coloc_signal_eligibility.R, which uses
# coloc:::logbf_to_pp(..., last_is_null = TRUE)).
#
# No SuSiE fit is run on the GWAS side. The GWAS fit's variants are a subset of
# the eQTL-fit variants (06_susie_coloc.R:670), so the GWAS-side share is 1 by
# construction; this is checked per row (gwas_ids_subset_of_eqtl). The GWAS fit
# is assumed to return >= 1 credible set, as it did in the first run (finite
# PP.H4.susie in the master).
#
# LD matching: get_ld_per_locus() reads the .ld matrix only to drop variants
# with NA correlations. For a block set whose .ld files contain no non-numeric
# character (checked with grep), that step cannot fire, so the variant set
# depends only on the .bim files; those rows use ld_variant_positions_bim_only()
# below, a copy of the cached get_ld_per_locus() with the matrix lines removed.
# Block sets with any non-numeric character go through the real
# get_ld_per_locus() (src/perf/get_ld_per_locus_cached.R). A sample of no-NA
# block sets is also run through the real function and must give identical
# positions (identity check).
#
# Usage: Rscript an1_share_only.R <rows.tsv> <out_dir> <chr_list|all> <mode> <n_identity>
suppressPackageStartupMessages({
  library(data.table); library(coloc); library(susieR)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)
ROWS_FILE <- args[1]; OUT <- args[2]; CHR_ARG <- args[3]; MODE <- args[4]
N_IDENTITY <- as.integer(args[5])
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(20261021)

PROJECT_ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FM_DIR <- file.path(PROJECT_ROOT, "GWAS/finemapping")
LDREF <- file.path(FM_DIR, "data/ld_ref")
# Same environment as the r2 eligibility rerun (src/perf/09b_coloc_eligibility_batch.sbatch).
Sys.setenv(MASLD_PROJECT_ROOT = PROJECT_ROOT,
           AFR_LD_DIR = file.path(LDREF, "1kg_afr_gram"),
           AMR_LD_DIR = file.path(LDREF, "1kg_amr_gram"),
           EAS_LD_DIR = file.path(LDREF, "1kg_eas_gram"),
           SAS_LD_DIR = file.path(LDREF, "1kg_sas_gram"),
           EQTL_SUSIE_DIR = file.path(FM_DIR, "results/eqtl_susie_polyfun"))
setwd(FM_DIR)

T0 <- Sys.time()
logf <- function(...) cat(sprintf("[%7.1fs] ", as.numeric(difftime(Sys.time(), T0, units = "secs"))),
                          ..., "\n", sep = "")

sha256 <- function(p) strsplit(system2("sha256sum", p, stdout = TRUE), " ")[[1]][1]
SRC06 <- file.path(FM_DIR, "src/06_susie_coloc.R")
SRC_EL <- file.path(FM_DIR, "src/coloc_signal_eligibility.R")
SRC_FF <- file.path(FM_DIR, "src/finemapping_functions.R")
SRC_LD <- file.path(FM_DIR, "src/perf/get_ld_per_locus_cached.R")
# sha256 recorded by the r2 rerun plan (results/susie_coloc_eligibility_20260923T211525Z/plan/input_manifest.tsv)
stopifnot(sha256(SRC06) == "e22b3ffbef48de9b96e220399275edc0113ec751fc765238ccefd076697114cc",
          sha256(SRC_EL) == "fa7a42c621a55e3dee80bdb8bddad5d5049f697e98f6eb1d4e310aef4fb30675")
manifest <- data.table(role = c("runner_06", "eligibility_helper", "finemapping_functions", "ld_cached", "rows"),
                       path = c(SRC06, SRC_EL, SRC_FF, SRC_LD, ROWS_FILE))
manifest[, sha256 := vapply(path, sha256, "")]

suppressPackageStartupMessages(source(SRC_FF))
source(SRC_LD)       # shadows get_ld_per_locus with the cached production version
source(SRC_EL)

# ---- code extracted verbatim from 06_susie_coloc.R --------------------------
L06 <- readLines(SRC06)
stopifnot(startsWith(L06[54], "PANEL_AF_ROOT <- Sys.getenv"),
          startsWith(L06[88], "strand_verdict_vs_panel <- function"), L06[98] == "}")
eval(parse(text = L06[54:98]))                     # panel AF loader + strand helpers
stopifnot(L06[111] == "COLOC_P1 <- 1e-4", startsWith(L06[121], "SUSIE_L"))
eval(parse(text = L06[111:121]))                   # MIN_SNPS, EQTL_SUSIE_DIR, MIN_TRIPLE_SNPS, ...
HB <- L06[395:604]
stopifnot(trimws(HB[1]) == "eqtl_gene <- eqtl_all[ENSG == gene_id]",
          trimws(HB[length(HB)]) == "if (nrow(merged) < MIN_SNPS) { n_skipped <- n_skipped + 1L; next }")
SKIP <- "{ n_skipped <- n_skipped + 1L; next }"
stopifnot(sum(grepl(SKIP, HB, fixed = TRUE)) == 3)
HB <- gsub(SKIP, "{ return(NULL) }", HB, fixed = TRUE)
stopifnot(!any(grepl("\\bnext\\b", HB)))
build_merged <- function(gene_id, eqtl_all, gwas, pa_gwas, pa_eqtl) NULL
body(build_merged) <- parse(text = c("{", HB, "merged", "}"))[[1]]
logf("06_susie_coloc.R harmonization block extracted (lines 395-604); EQTL_SUSIE_DIR=", EQTL_SUSIE_DIR,
     " MIN_SNPS=", MIN_SNPS, " MIN_TRIPLE_SNPS=", MIN_TRIPLE_SNPS)

# ---- bim-only copy of the cached get_ld_per_locus() --------------------------
# Identical to src/perf/get_ld_per_locus_cached.R except that the .ld matrix is
# never read and the NA-removal step is absent. Valid only for block sets whose
# .ld files have no non-numeric character (checked before use).
ld_variant_positions_bim_only <- function(ss_per_locus, LOCUS, CHR, START, END, ancestry = "EUR") {
  path_to_LD <- find_LD_block(CHR, START, END, ancestry)
  if (is.null(path_to_LD)) return(NULL)
  bim <- data.frame()
  for (file in path_to_LD) {
    bim_file <- paste0(file, ".bim")
    ld_file <- paste0(file, ".ld")
    if (!file.exists(bim_file) || !file.exists(ld_file)) {
      warning(paste("LD files not found:", file))
      next
    }
    BIM <- fread(bim_file)
    bim <- rbind(bim, BIM)
  }
  if (nrow(bim) == 0) return(NULL)
  colnames(bim) <- c("chr", "rsid", "dk", "pos", "alt", "ref")
  bim$SNP <- paste(bim$chr, bim$pos, bim$alt, bim$ref, sep = ":")
  keep_uniq <- !duplicated(bim$SNP)
  if (!all(keep_uniq)) bim <- bim[keep_uniq]
  bim <- data.table::copy(bim)
  ss_per_locus$SNP1 <- paste(ss_per_locus$chromosome, ss_per_locus$position,
                             ss_per_locus$allele1, ss_per_locus$allele2, sep = ":")
  ss_per_locus$SNP2 <- paste(ss_per_locus$chromosome, ss_per_locus$position,
                             ss_per_locus$allele2, ss_per_locus$allele1, sep = ":")
  bim$ss_matched <- bim$SNP %in% c(ss_per_locus$SNP2, ss_per_locus$SNP1)
  n_matched <- sum(bim$ss_matched)
  if (n_matched < 10) return(NULL)
  bim_filtered <- bim[bim$ss_matched, ]
  ss_filtered <- rbind(
    ss_per_locus %>% filter(SNP1 %in% bim_filtered$SNP) %>% rename(SNP = SNP1) %>% select(!c(SNP2)),
    ss_per_locus %>% filter(SNP2 %in% bim_filtered$SNP) %>% mutate(beta = -1 * beta) %>%
      rename(SNP = SNP2, allele1 = allele2, allele2 = allele1) %>% select(!c(SNP1))
  )
  ss_filtered <- ss_filtered[!duplicated(ss_filtered$SNP), ]
  ss_filtered <- ss_filtered[match(bim_filtered$SNP, ss_filtered$SNP), ]
  ss_filtered <- ss_filtered[!is.na(ss_filtered$SNP), ]
  keep <- bim_filtered$SNP %in% ss_filtered$SNP
  if (sum(keep) < 10) return(NULL)
  bim_filtered <- bim_filtered[keep, ]
  ss_filtered <- ss_filtered[match(bim_filtered$SNP, ss_filtered$SNP), ]
  ss_filtered
}

# ---- non-numeric scan of .ld files (cached across jobs) ----------------------
FLAG_DIR <- file.path(dirname(OUT), "ld_nonnumeric_flags")
dir.create(FLAG_DIR, showWarnings = FALSE)
FLAG_FILE <- file.path(FLAG_DIR, paste0("flags_", basename(OUT), ".tsv"))
FLAGS <- new.env()
for (f in list.files(FLAG_DIR, pattern = "^flags_.*\\.tsv$", full.names = TRUE)) {
  x <- tryCatch(fread(f, header = FALSE, col.names = c("path", "code")), error = function(e) NULL)
  if (!is.null(x)) for (k in seq_len(nrow(x))) assign(x$path[k], x$code[k], envir = FLAGS)
}
scan_ld_files <- function(paths) {
  todo <- unique(paths[!vapply(paths, exists, TRUE, envir = FLAGS, inherits = FALSE)])
  if (!length(todo)) return(invisible())
  lst <- tempfile(tmpdir = OUT); writeLines(todo, lst)
  res <- tempfile(tmpdir = OUT)
  t1 <- Sys.time()
  # exit 0: a non-numeric character exists; 1: none; 2: error/missing
  cmd <- sprintf("xargs -a %s -P 8 -I{} sh -c 'LC_ALL=C grep -q -m1 \"[^-0-9.eE+ ]\" \"{}\"; printf \"%%s\\t%%s\\n\" \"{}\" $?' > %s",
                 shQuote(lst), shQuote(res))
  system(cmd)
  x <- fread(res, header = FALSE, col.names = c("path", "code"), sep = "\t")
  stopifnot(setequal(x$path, todo))
  for (k in seq_len(nrow(x))) assign(x$path[k], x$code[k], envir = FLAGS)
  fwrite(x, FLAG_FILE, sep = "\t", append = TRUE, col.names = FALSE)
  logf("scanned ", length(todo), " .ld files in ",
       round(as.numeric(difftime(Sys.time(), t1, units = "secs"))), " s; with non-numeric char: ",
       sum(x$code == 0), "; errors: ", sum(x$code == 2))
  unlink(c(lst, res))
}
block_set_clean <- function(paths) {
  codes <- vapply(paste0(paths, ".ld"), function(p) get(p, envir = FLAGS), 0)
  all(codes == 1)
}

# ---- inputs -------------------------------------------------------------------
registry <- read.delim("config/gwas_registry.tsv", stringsAsFactors = FALSE)
af_src <- fread("config/gwas_af_sources.tsv")
# 2022_36402844_PDFF_EUR_af.tsv was rewritten 2026-10-10 (allele labels swapped);
# the r2 rerun (2026-09-23) read the copy now saved as .before_label_fix_20261010T145917Z.
AF_OVERRIDE <- c("2022_36402844_PDFF_EUR" =
  file.path(FM_DIR, "data/sumstats_af/2022_36402844_PDFF_EUR_af.tsv.before_label_fix_20261010T145917Z"))
gwas_path_of <- function(st) {
  if (st %in% names(AF_OVERRIDE)) return(AF_OVERRIDE[[st]])
  hit <- af_src[study_name == st & file.exists(af_sumstats_path)]
  if (nrow(hit) == 1L) return(hit$af_sumstats_path[[1]])
  registry$sumstats_path[registry$study_name == st]
}

rows <- fread(ROWS_FILE, colClasses = list(character = c("gwas_name", "ensembl", "gene")))
chrs <- if (CHR_ARG == "all") sort(unique(rows$chr)) else as.integer(strsplit(CHR_ARG, ",")[[1]])
rows <- rows[chr %in% chrs]
rows[, row_id := .I]
logf("mode ", MODE, " | rows ", nrow(rows), " | chromosomes ", paste(chrs, collapse = ","),
     " | studies ", uniqueN(rows$gwas_name))
studies <- sort(unique(rows$gwas_name))
st_info <- rbindlist(lapply(studies, function(st) {
  r <- registry[registry$study_name == st, ]
  anc <- r$ancestry; if (length(anc) == 0 || is.na(anc) || anc == "") anc <- "EUR"
  data.table(gwas_name = st, ancestry = anc, sumstats = gwas_path_of(st))
}))
fwrite(st_info, file.path(OUT, "study_inputs.tsv"), sep = "\t")

# GWAS sumstats: read each study once, keep this job's chromosomes.
GW <- list()
for (st in studies) {
  t1 <- Sys.time()
  g <- fread(st_info[gwas_name == st, sumstats], showProgress = FALSE)
  g <- g[chromosome %in% unique(rows[gwas_name == st, chr])]
  GW[[st]] <- split(g, by = "chromosome", keep.by = TRUE)
  rm(g); gc(FALSE)
  logf("GWAS ", st, " loaded (", round(as.numeric(difftime(Sys.time(), t1, units = "secs"))), " s)")
}

EQTL_DIR <- file.path(PROJECT_ROOT, "data/broadaway_eqtl")
P2 <- eval(formals(coloc:::coloc.bf_bf)$p2)
OV <- eval(formals(coloc:::coloc.bf_bf)$overlap.min)
stopifnot(isTRUE(all.equal(OV, COLOC_OVERLAP_MIN)))

row_out <- list(); cs_out <- list(); id_out <- list()
n_identity_done <- 0L

for (chr_num in chrs) {
  rc <- rows[chr == chr_num]
  if (!nrow(rc)) next
  t_chr <- Sys.time()
  eqtl_all <- fread(file.path(EQTL_DIR, paste0("chr", chr_num, "_marginal_summary_results.tsv")),
                    showProgress = FALSE)
  pa_eqtl <- load_panel_af("eur", chr_num)
  PA <- list()
  for (anc in unique(st_info[gwas_name %in% rc$gwas_name, ancestry])) PA[[anc]] <- load_panel_af(anc, chr_num)
  gw_chr <- list()
  for (st in unique(rc$gwas_name)) {
    g <- GW[[st]][[as.character(chr_num)]]
    if (is.null(g)) g <- GW[[st]][[1]][0]
    g <- copy(g)
    g[, merge_key := paste(chromosome, position, sep = ":")]
    gw_chr[[st]] <- g
  }
  logf("chr ", chr_num, ": ", nrow(rc), " rows, eQTL table ", nrow(eqtl_all), " variants")

  # ---- pass 1: harmonize, restrict to eQTL-fit positions, locate LD blocks ----
  P1 <- vector("list", nrow(rc))
  eq_cache <- new.env()
  for (k in seq_len(nrow(rc))) {
    r <- rc[k]; gene_id <- r$ensembl; st <- r$gwas_name
    anc <- st_info[gwas_name == st, ancestry]
    rec <- list(row_id = r$row_id, gwas_name = st, ensembl = gene_id, chr = chr_num, ancestry = anc,
                status = NA_character_, n_snps_rebuilt = NA_integer_, n_gwas_sub = NA_integer_)
    merged <- tryCatch(build_merged(gene_id, eqtl_all, gw_chr[[st]], PA[[anc]], pa_eqtl),
                       error = function(e) structure(conditionMessage(e), class = "err"))
    if (inherits(merged, "err")) { rec$status <- paste0("error_harmonize: ", merged); P1[[k]] <- rec; next }
    if (is.null(merged)) { rec$status <- "prereq_fail_min_snps"; P1[[k]] <- rec; next }
    rec$n_snps_rebuilt <- nrow(merged)
    eqtl_rds <- file.path(EQTL_SUSIE_DIR, paste0("chr", chr_num), paste0(gene_id, "_susie.rds"))
    if (!exists(gene_id, envir = eq_cache, inherits = FALSE)) {
      e <- if (!file.exists(eqtl_rds)) NULL else {
        s <- readRDS(eqtl_rds)
        cs <- s$sets$cs_index
        has_cs <- !(is.null(s$sets$cs) || length(s$sets$cs) == 0)
        list(converged = isTRUE(s$converged), ids = colnames(s$lbf_variable), cs = cs, has_cs = has_cs,
             bf2 = if (has_cs) s$lbf_variable[cs, , drop = FALSE] else NULL)
      }
      assign(gene_id, e, envir = eq_cache)
    }
    e <- get(gene_id, envir = eq_cache)
    if (is.null(e)) { rec$status <- "prereq_fail_no_eqtl_fit"; P1[[k]] <- rec; next }
    if (!e$converged) { rec$status <- "prereq_fail_eqtl_not_converged"; P1[[k]] <- rec; next }
    eqtl_positions <- as.integer(sub("^[0-9]+:", "", e$ids))
    gwas_sub <- merged[eqtl_pos %in% eqtl_positions]                          # 06_susie_coloc.R:670
    rec$n_gwas_sub <- nrow(gwas_sub)
    if (nrow(gwas_sub) < MIN_TRIPLE_SNPS) { rec$status <- "prereq_fail_gwas_sub_lt_50"; P1[[k]] <- rec; next }
    rec$ss_for_ld <- data.frame(chromosome = chr_num, position = gwas_sub$eqtl_pos,
                                allele1 = gwas_sub$gwas_a1, allele2 = gwas_sub$gwas_a2,
                                beta = gwas_sub$gwas_beta, standard_error = gwas_sub$gwas_se)
    rec$START <- min(gwas_sub$eqtl_pos); rec$END <- max(gwas_sub$eqtl_pos)
    paths <- suppressMessages(capture.output(p <- find_LD_block(chr_num, rec$START, rec$END, anc)))
    if (is.null(p)) { rec$status <- "prereq_fail_no_ld_block"; P1[[k]] <- rec; next }
    rec$paths <- p
    rec$key <- paste(c(anc, p), collapse = "|")
    P1[[k]] <- rec
  }
  have_key <- vapply(P1, function(x) !is.null(x$key), TRUE)
  scan_ld_files(unique(paste0(unlist(lapply(P1[have_key], `[[`, "paths")), ".ld")))
  logf("chr ", chr_num, " pass 1 done: ", sum(have_key), " of ", length(P1), " rows reach LD matching")

  # ---- pass 2: LD matching grouped by block set, then eQTL shares -------------
  keys <- unique(vapply(P1[have_key], `[[`, "", "key"))
  for (key in keys) {
    idx <- which(vapply(P1, function(x) identical(x$key, key), TRUE))
    paths <- P1[[idx[1]]]$paths
    clean <- block_set_clean(paths)
    do_identity <- clean && n_identity_done < N_IDENTITY && (MODE == "validate" || runif(1) < 0.05)
    for (k in idx) {
      rec <- P1[[k]]
      ss_ld <- tryCatch({
        if (clean) {
          ld_variant_positions_bim_only(rec$ss_for_ld, rec$ensembl, chr_num, rec$START, rec$END, rec$ancestry)
        } else {
          res <- suppressWarnings(get_ld_per_locus(rec$ss_for_ld, rec$ensembl, chr_num, rec$START, rec$END,
                                                   rec$ancestry))
          if (is.null(res)) NULL else res[[1]]
        }
      }, error = function(e) structure(conditionMessage(e), class = "err"))
      rec$ld_route <- if (clean) "bim_only" else "full_matrix"
      if (inherits(ss_ld, "err")) { rec$status <- paste0("error_ld: ", ss_ld); P1[[k]] <- rec; next }
      if (do_identity && k == idx[1]) {
        ref <- suppressWarnings(get_ld_per_locus(rec$ss_for_ld, rec$ensembl, chr_num, rec$START, rec$END,
                                                 rec$ancestry))
        same <- (is.null(ref) && is.null(ss_ld)) ||
          (!is.null(ref) && !is.null(ss_ld) && identical(as.integer(ref[[1]]$position), as.integer(ss_ld$position)) &&
             identical(ref[[1]]$SNP, ss_ld$SNP))
        id_out[[length(id_out) + 1]] <- data.table(ld_key = key, gwas_name = rec$gwas_name, ensembl = rec$ensembl,
                                                   n_bim_only = if (is.null(ss_ld)) NA_integer_ else nrow(ss_ld),
                                                   n_full = if (is.null(ref)) NA_integer_ else nrow(ref[[1]]),
                                                   identical = same)
        n_identity_done <- n_identity_done + 1L
        logf("identity check ", rec$gwas_name, " ", rec$ensembl, ": ", same)
        rm(ref)
      }
      if (is.null(ss_ld) || nrow(ss_ld) < MIN_TRIPLE_SNPS) {
        rec$status <- "prereq_fail_ld_matched_lt_50"; rec$n_ld <- if (is.null(ss_ld)) 0L else nrow(ss_ld)
        P1[[k]] <- rec; next
      }
      snp_ids <- paste0(chr_num, ":", ss_ld$position)
      e <- get(rec$ensembl, envir = eq_cache)
      common_snps <- intersect(snp_ids, e$ids)
      rec$n_ld <- nrow(ss_ld); rec$n_common <- length(common_snps)
      rec$gwas_ids_subset_of_eqtl <- all(snp_ids %in% e$ids)
      if (length(common_snps) < MIN_TRIPLE_SNPS) { rec$status <- "prereq_fail_common_lt_50"; P1[[k]] <- rec; next }
      if (!e$has_cs) { rec$status <- "prereq_fail_eqtl_no_cs"; P1[[k]] <- rec; next }
      shared <- setdiff(intersect(snp_ids, colnames(e$bf2)), "null")
      sp2 <- signal_shared_prop(e$bf2, shared, P2)
      u2 <- signal_shared_prop_uniform(e$bf2, shared)
      elig <- !(sp2$prop < COLOC_OVERLAP_MIN)
      cs_out[[length(cs_out) + 1]] <- data.table(gwas_name = rec$gwas_name, ensembl = rec$ensembl, chr = chr_num,
                                                 idx2 = e$cs, prop_eqtl = sp2$prop, prop_eqtl_uniform = u2,
                                                 eligible_cs = elig)
      rec$n_eqtl_cs <- length(e$cs); rec$n_eqtl_cs_eligible <- sum(elig %in% TRUE)
      rec$status <- if (rec$n_eqtl_cs_eligible == 0) "untestable_all_pairs_fail" else
        if (rec$n_eqtl_cs_eligible < rec$n_eqtl_cs) "eligible_some_pairs_fail" else "eligible_all_pairs"
      P1[[k]] <- rec
    }
  }
  row_out[[length(row_out) + 1]] <- rbindlist(lapply(P1, function(x) data.table(
    row_id = x$row_id, gwas_name = x$gwas_name, ensembl = x$ensembl, chr = x$chr, ancestry = x$ancestry,
    status = x$status, n_snps_rebuilt = x$n_snps_rebuilt, n_gwas_sub = x$n_gwas_sub,
    n_ld = if (is.null(x$n_ld)) NA_integer_ else x$n_ld,
    n_common = if (is.null(x$n_common)) NA_integer_ else x$n_common,
    gwas_ids_subset_of_eqtl = if (is.null(x$gwas_ids_subset_of_eqtl)) NA else x$gwas_ids_subset_of_eqtl,
    n_eqtl_cs = if (is.null(x$n_eqtl_cs)) NA_integer_ else x$n_eqtl_cs,
    n_eqtl_cs_eligible = if (is.null(x$n_eqtl_cs_eligible)) NA_integer_ else x$n_eqtl_cs_eligible,
    ld_route = if (is.null(x$ld_route)) NA_character_ else x$ld_route,
    ld_key = if (is.null(x$key)) NA_character_ else x$key)), fill = TRUE)
  rm(eqtl_all, gw_chr, PA, pa_eqtl, P1, eq_cache); .LD_CACHE$entry <- NULL; gc(FALSE)
  logf("chr ", chr_num, " done in ", round(as.numeric(difftime(Sys.time(), t_chr, units = "mins")), 1), " min")
  # incremental writes so a time-out keeps finished chromosomes
  fwrite(rbindlist(row_out, fill = TRUE), file.path(OUT, "row_shares.tsv"), sep = "\t")
  fwrite(rbindlist(cs_out), file.path(OUT, "eqtl_cs_shares.tsv"), sep = "\t")
  if (length(id_out)) fwrite(rbindlist(id_out), file.path(OUT, "identity_check.tsv"), sep = "\t")
}

ro <- rbindlist(row_out, fill = TRUE)
ro <- merge(rows, ro[, !c("gwas_name", "ensembl", "chr")], by = "row_id", all.x = TRUE)
fwrite(ro, file.path(OUT, "row_shares.tsv"), sep = "\t")
fwrite(rbindlist(cs_out), file.path(OUT, "eqtl_cs_shares.tsv"), sep = "\t")
if (length(id_out)) fwrite(rbindlist(id_out), file.path(OUT, "identity_check.tsv"), sep = "\t")
fwrite(manifest, file.path(OUT, "input_sha256.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
logf("status counts: ", paste(names(table(ro$status, useNA = "ifany")), table(ro$status, useNA = "ifany"),
                              sep = "=", collapse = "; "))
if (length(id_out) && !all(rbindlist(id_out)$identical)) stop("IDENTITY CHECK FAILED")
logf("DONE")
