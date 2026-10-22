#!/usr/bin/env Rscript
# H3 / H31 step 2: one phenotype library per person, frozen Hotspot program scores and the
# person-level count matrix for the carrier signature. Rule in PRESPEC_H3.md.
# No genotype is read here.
# Usage: 02_program_scores.R GENO_DIR OUT_DIR
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})
args <- commandArgs(trailingOnly = TRUE)
geno_dir <- args[[1]]
out <- args[[2]]
stopifnot(!dir.exists(out))
dir.create(out, recursive = TRUE)
set.seed(20260923)
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"

REFIT <- c("GSE130970", "GSE135251", "GSE162694", "GSE213621")
OTHER <- c("GSE174478", "GSE193066", "GSE240729")
PROGRAMS <- c(lipid = "hotspot_hepatocytes_6267e64841a0621c",
              ecm = "hotspot_hepatocytes_f05c535ae5bbc0b9",
              upr = "hotspot_hepatocytes_b75cab32b41701a0")

elig <- fread(file.path(geno_dir, "library_eligibility.tsv"))[eligible == TRUE]
n1 <- fread(file.path(repo, "Analysis/Hypotheses/N1-20261007T210855Z/library_table.tsv"),
            select = c("run", "expr_sex", "input_reads"))
roster <- fread(file.path(repo, "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/resource_identity_roster.tsv"),
                select = c("sample_id", "dataset", "fibrosis_stage", "nas_score", "inferred_sex", "group_binary", "first_biopsy"))
unified <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv",
                 select = c("sample_id", "sex"))
cand <- file.path(repo, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z")
unit774 <- fread(file.path(cand, "primary_unit_roster.tsv"))
dge774 <- readRDS(file.path(cand, "primary_source_controls/merged_dge.rds"))
counts_dir <- file.path(repo, "results/remediation/bg001/bg001-fragment-v211-gencode49-20260807T194845Z/frozen_sets/read_counts")
membership <- fread(file.path(repo, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"))
annotation <- fread(file.path(repo, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"))

base_id <- function(x) sub("\\.[0-9]+.*$", "", x)
annotation[, symbol := toupper(gene_name)]
id2sym <- annotation[symbol != "" & !is.na(symbol), .(symbol = symbol[1]), by = ensembl_base]
memb_all <- membership[program_uid %in% PROGRAMS & !is.na(mapped_symbol) & mapped_symbol != "",
                       .(weight = sum(as.numeric(original_l1_weight))), by = .(program_uid, symbol = toupper(mapped_symbol))]
total_w <- membership[program_uid %in% PROGRAMS, .(total = sum(as.numeric(original_l1_weight))), by = program_uid]

meta <- merge(elig[, .(library, cohort, person_id)], n1[, .(library = run, expr_sex, input_reads)], by = "library", all.x = TRUE)
meta <- merge(meta, roster[, .(library = sample_id, fibrosis_stage, nas_score, inferred_sex, group_binary, first_biopsy)], by = "library", all.x = TRUE)
meta <- merge(meta, unified[, .(library = sample_id, sex_rec = sex)], by = "library", all.x = TRUE)

pick_sex <- function(a, b, c) {
  f <- function(x) { x <- toupper(substr(as.character(x), 1, 1)); ifelse(x %in% c("M", "F"), x, NA_character_) }
  r <- f(a); r[is.na(r)] <- f(b)[is.na(r)]; r[is.na(r)] <- f(c)[is.na(r)]; r
}
meta[, sex := pick_sex(inferred_sex, sex_rec, expr_sex)]
meta[, stage := suppressWarnings(as.numeric(fibrosis_stage))]
meta[is.na(stage) & group_binary == "Control", stage := 0]

score_cohort <- function(m, ref, sel) {
  lc <- cpm(DGEList(m), log = TRUE, prior.count = 1)
  lc <- lc[rowSums(m) > ncol(m), , drop = FALSE]
  sym <- id2sym$symbol[match(base_id(rownames(lc)), id2sym$ensembl_base)]
  ok <- !is.na(sym)
  s <- rowsum(lc[ok, , drop = FALSE], sym[ok], reorder = FALSE)
  s <- s / as.numeric(table(sym[ok])[rownames(s)])
  mu <- rowMeans(s[, ref, drop = FALSE]); sdv <- apply(s[, ref, drop = FALSE], 1, sd)
  keep <- is.finite(sdv) & sdv > 0
  z <- (s[keep, , drop = FALSE] - mu[keep]) / sdv[keep]
  sc <- matrix(NA_real_, length(PROGRAMS), ncol(z), dimnames = list(names(PROGRAMS), colnames(z)))
  tst <- list()
  for (nm in names(PROGRAMS)) {
    p <- PROGRAMS[[nm]]
    mm <- memb_all[program_uid == p & symbol %in% rownames(z)]
    frac <- sum(mm$weight) / total_w[program_uid == p, total]
    test <- nrow(mm) >= 10L && frac >= 0.8
    tst[[nm]] <- data.table(program = nm, program_uid = p, n_observed_genes = nrow(mm), observed_l1_fraction = frac, testable = test)
    if (!test) next
    v <- as.numeric(crossprod(mm$weight / sum(mm$weight), z[mm$symbol, , drop = FALSE]))
    names(v) <- colnames(z)
    sc[nm, ] <- (v - mean(v[ref])) / sd(v[ref])
  }
  list(scores = sc, testability = rbindlist(tst))
}

people <- list(); testab <- list(); cnts <- list()
for (co in c(REFIT, OTHER)) {
  mc <- meta[cohort == co]
  if (co %in% REFIT) {
    ref <- intersect(unit774[sample_id %in% roster[dataset == co, sample_id], sample_id], colnames(dge774$counts))
    m <- dge774$counts[, ref, drop = FALSE]
    src <- "corrected_774"
    ph <- mc[library %in% ref][order(match(library, ref))]
    ph <- ph[!duplicated(person_id)]
  } else {
    ct <- fread(file.path(counts_dir, co, "gene_counts.txt"), skip = 1)
    m0 <- as.matrix(ct[, 7:ncol(ct)])
    colnames(m0) <- sub("^.*?(SRR[0-9]+).*$", "\\1", basename(colnames(m0)))
    rownames(m0) <- ct[[1]]
    src <- "historical_fragment_counts"
    ph <- mc[library %in% colnames(m0)]
    ph[, fb := as.integer(co == "GSE193066" & first_biopsy %in% c("True", "TRUE", TRUE))]
    ph <- ph[order(person_id, -fb, -as.numeric(input_reads), library)]
    ph <- ph[!duplicated(person_id)]
    ref <- ph$library
    m <- m0[, ref, drop = FALSE]
  }
  res <- score_cohort(m, ref, ph$library)
  testab[[co]] <- cbind(cohort = co, res$testability)
  sc <- t(res$scores[, ph$library, drop = FALSE])
  colnames(sc) <- paste0("score_", colnames(sc))
  people[[co]] <- cbind(ph[, .(person_id, cohort, library, stage, sex, group_binary, nas_score = suppressWarnings(as.numeric(nas_score)))],
                        count_source = src, as.data.table(sc))
  mm <- m[, ph$library, drop = FALSE]
  colnames(mm) <- ph$person_id
  rownames(mm) <- base_id(rownames(mm))
  mm <- mm[!duplicated(rownames(mm)), , drop = FALSE]
  cnts[[co]] <- mm
  message(co, ": persons ", nrow(ph), "; reference libraries ", length(ref), "; source ", src)
}
people <- rbindlist(people)
testab <- rbindlist(testab)
fwrite(people, file.path(out, "person_phenotypes.tsv"), sep = "\t")
fwrite(testab, file.path(out, "program_testability.tsv"), sep = "\t")
genes <- Reduce(intersect, lapply(cnts, rownames))
cm <- do.call(cbind, lapply(cnts, function(x) x[genes, , drop = FALSE]))
stopifnot(identical(colnames(cm), people$person_id))
saveRDS(list(counts = cm, cohort = people$cohort, library = people$library), file.path(out, "person_counts.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("persons", nrow(people), "genes", length(genes), "\n")
print(people[, .N, by = .(cohort, count_source)])
print(testab)
