#!/usr/bin/env Rscript
# AG2X A3b/A4 input: per-variant PIPs from the production PolyFun-LD Broadaway eQTL SuSiE fits
# (GWAS/finemapping/results/eqtl_susie_polyfun), genes whose top PIP >= 0.5 (0.9 is the primary
# threshold; 0.5 is kept for the sensitivity set), excluding every fit with a duplicated CHR:POS
# identifier (COLOC-IDENTITY audit table; superset of the 90 primary genes named in docs/PAPER.md).
# Joins hg19 alleles, EAF and beta from the marginal eQTL file (gene x position must be a single row),
# then lifts hg19 -> hg38 with data/broadaway_eqtl/hg19ToHg38.over.chain. Read-only on all inputs.
suppressPackageStartupMessages({library(data.table); library(rtracklayer); library(GenomicRanges)})
set.seed(20261111)
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
O <- file.path(R, "Analysis/Hypotheses/AG2X-20261112T000534Z/inputs")
FITS <- file.path(R, "GWAS/finemapping/results/eqtl_susie_polyfun")
AUDIT <- file.path(R, "Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z/audit/all_saved_eqtl_fit_identity_checks.tsv")
EQTL <- file.path(R, "data/broadaway_eqtl")
CHAIN <- file.path(R, "data/broadaway_eqtl/hg19ToHg38.over.chain")

summ <- rbindlist(lapply(1:22, function(ch) fread(file.path(FITS, sprintf("chr%d_summary.tsv", ch)))))
summ <- summ[status == "ok"]
summ <- unique(summ, by = "ensg")
audit <- fread(AUDIT)
counts <- list(fits_ok = nrow(summ), fits_top_pip_ge_0.9 = summ[max_pip >= 0.9, .N],
               fits_top_pip_ge_0.5 = summ[max_pip >= 0.5, .N])
summ <- merge(summ, audit[, .(ensg = ensembl, n_duplicate_ids)], by = "ensg", all.x = TRUE)
counts$top_ge_0.9_without_audit_row <- summ[max_pip >= 0.9 & is.na(n_duplicate_ids), .N]
counts$top_ge_0.9_with_duplicate_ids <- summ[max_pip >= 0.9 & n_duplicate_ids > 0, .N]
counts$top_ge_0.5_with_duplicate_ids <- summ[max_pip >= 0.5 & n_duplicate_ids > 0, .N]
sel <- summ[max_pip >= 0.5 & !is.na(n_duplicate_ids) & n_duplicate_ids == 0]
counts$genes_after_duplicate_exclusion_ge_0.5 <- nrow(sel)
counts$genes_after_duplicate_exclusion_ge_0.9 <- sel[max_pip >= 0.9, .N]

rows <- vector("list", nrow(sel))
for (i in seq_len(nrow(sel))) {
  f <- file.path(FITS, paste0("chr", sel$chr[i]), paste0(sel$ensg[i], "_susie.rds"))
  s <- readRDS(f)
  ids <- colnames(s$lbf_variable)
  stopifnot(length(ids) == length(s$pip), !anyDuplicated(ids))
  rows[[i]] <- data.table(ensg = sel$ensg[i], gene = sel$gene[i], chr = sel$chr[i],
                          pos_hg19 = as.integer(sub("^[0-9]+:", "", ids)), pip = as.numeric(s$pip))
}
pips <- rbindlist(rows)
cat("PIP rows", nrow(pips), "genes", uniqueN(pips$ensg), "\n")

# keep tops (PIP >= 0.5, the gene's maximum) and control candidates (PIP < 0.01)
pips[, is_top := pip == max(pip) & pip >= 0.5, by = ensg]
counts$genes_with_tied_top <- pips[is_top == TRUE, .N, by = ensg][N > 1, .N]
keep <- pips[is_top == TRUE | pip < 0.01]

out <- list()
for (ch in sort(unique(keep$chr))) {
  k <- keep[chr == ch]
  e <- fread(file.path(EQTL, sprintf("chr%d_marginal_summary_results.tsv", ch)),
             select = c("CHR", "POS", "NEA", "EA", "EAF", "Beta", "SE", "PVAL", "N", "ENSG"))
  e <- e[ENSG %in% unique(k$ensg)]
  n_rows <- e[, .N, by = .(ENSG, POS)]
  e <- merge(e, n_rows, by = c("ENSG", "POS"))
  m <- merge(k, e, by.x = c("ensg", "pos_hg19"), by.y = c("ENSG", "POS"), all.x = TRUE)
  out[[as.character(ch)]] <- m
  rm(e); gc(verbose = FALSE)
  cat("chr", ch, nrow(m), "rows\n")
}
d <- rbindlist(out, fill = TRUE)
counts$rows_without_marginal_row <- d[is.na(EA), .N]
counts$rows_multiallelic_position <- d[!is.na(N.y) & N.y > 1, .N]
setnames(d, c("N.x", "N.y"), c("eqtl_n", "n_rows_at_position"), skip_absent = TRUE)

# liftover hg19 -> hg38 (1-based single base)
gr <- GRanges(paste0("chr", d$chr), IRanges(d$pos_hg19, d$pos_hg19), strand = "+")
mcols(gr)$i <- seq_along(gr)
lo <- liftOver(gr, import.chain(CHAIN))
n_hits <- elementNROWS(lo)
d$liftover_hits <- n_hits
ul <- unlist(lo[n_hits == 1])
d$pos_hg38 <- NA_integer_; d$chr_hg38 <- NA_character_; d$strand_hg38 <- NA_character_
idx <- mcols(ul)$i
d$pos_hg38[idx] <- start(ul); d$chr_hg38[idx] <- as.character(seqnames(ul)); d$strand_hg38[idx] <- as.character(strand(ul))
counts$rows_liftover_not_unique <- sum(n_hits != 1)
fwrite(d, file.path(O, "eqtl_pip_candidates_hg38.tsv.gz"), sep = "\t")
write(jsonlite::toJSON(counts, auto_unbox = TRUE, pretty = TRUE), file.path(O, "eqtl_pip_extraction_counts.json"))
capture.output(sessionInfo(), file = file.path(O, "sessionInfo_extract.txt"))
print(counts)
