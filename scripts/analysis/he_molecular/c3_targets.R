#!/usr/bin/env Rscript
# C3a: Model C targets and covariates from GTEx v8 liver (PRESPEC_C sections 4-6).
# One liver RNA-seq sample per C1-matched donor. log2(TPM + 1) -> collapse_symbols ->
# score_programs (pinned scorer, coverage 0.80). Per-program ceiling: split-gene-half
# reliability (20 seeded random halvings of member genes, Spearman-Brown). Negative
# controls: unweighted mean z of the two pinned gene lists. Development and sealed
# donors are written to separate files; development fits never read the sealed file.
# Usage: Rscript c3_targets.R <c1_exec_dir> <out_dir>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
c1 <- args[[1]]; out <- args[[2]]
dir.create(file.path(out, "sealed"), recursive = TRUE, showWarnings = FALSE)
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
here <- file.path(repo, "scripts/analysis/he_molecular")
owd <- setwd(here); stopifnot(system2("sha256sum", c("-c", "--quiet", "PRESPEC.sha256")) == 0); setwd(owd)
source(file.path(repo, "scripts/manuscript/program_observability_map/analysis_lib.R"))
gt <- file.path(repo, "data/external/allelic_refs/gtex_v8")
rel <- file.path(repo, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10")
annotation <- fread(file.path(rel, "inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"), select = c("gene_id", "gene_name"))
hot <- file.path(repo, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot")
registry <- fread(file.path(hot, "program_registry_v2.tsv")); membership <- fread(file.path(hot, "program_membership_v2.tsv"))

seal <- fread(file.path(here, "donor_seal.tsv"))
attr <- fread(file.path(gt, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"), sep = "\t", quote = "")
subj <- fread(file.path(gt, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"))
liv <- attr[SMTSD == "Liver" & SMAFRZE == "RNASEQ"][, SUBJID := sub("^(GTEX-[^-]+)-.*$", "\\1", SAMPID)]
liv <- liv[SUBJID %in% seal$SUBJID][order(SUBJID, -SMRIN)][!duplicated(SUBJID)]
cov <- merge(liv[, .(SAMPID, SUBJID, ischemic_min = SMTSISCH, rin = SMRIN, autolysis = SMATSSCR, path_notes = SMPTHNTS)],
             subj[, .(SUBJID, sex = SEX, age_bracket = AGE, hardy = DTHHRDY)], by = "SUBJID")
cov <- merge(cov, seal, by = "SUBJID")

tpm <- fread(file.path(gt, "gene_tpm_2017-06-05_v8_liver.gct.gz"), skip = 2)
keep_cols <- intersect(cov$SAMPID, names(tpm))
m <- as.matrix(tpm[, ..keep_cols]); rownames(m) <- tpm$Name
m <- m[rowMeans(m >= 1) >= 0.10, , drop = FALSE]
lib <- log2(m + 1)
sm <- collapse_symbols(lib, rownames(lib), annotation)
colnames(sm) <- cov$SUBJID[match(colnames(sm), cov$SAMPID)]
meta <- data.table(sample_id = colnames(sm), dataset = "GTEx_liver")
s <- score_programs(sm, meta, membership, registry, 0.80)
scores <- t(s$primary)

# split-gene-half ceiling per program, on the same z-scored symbol matrix
Z <- zscore_rows(sm); Z <- Z[rowSums(is.finite(Z)) == ncol(Z), , drop = FALSE]
mem <- membership[!is.na(mapped_symbol) & mapped_symbol != "", .(weight = sum(as.numeric(original_l1_weight))),
                  by = .(program_uid, gene_symbol = toupper(mapped_symbol))][gene_symbol %in% rownames(Z)]
set.seed(20260923)
ceil <- mem[, {
  if (.N < 6) list(ceiling = NA_real_, n_genes = .N) else {
    r <- replicate(20, {
      h <- sample(rep(1:2, length.out = .N))
      a <- as.numeric(crossprod(weight[h == 1] / sum(weight[h == 1]), Z[gene_symbol[h == 1], , drop = FALSE]))
      b <- as.numeric(crossprod(weight[h == 2] / sum(weight[h == 2]), Z[gene_symbol[h == 2], , drop = FALSE]))
      cor(a, b)
    })
    rr <- mean(r); list(ceiling = 2 * rr / (1 + rr), n_genes = .N)
  }
}, by = program_uid]
neg <- sapply(c(hypoxia = "negative_control_hypoxia_hallmark.txt", immediate_early = "negative_control_immediate_early.txt"),
              function(f) {g <- intersect(toupper(readLines(file.path(here, f))), rownames(Z)); colMeans(Z[g, , drop = FALSE])})

tab <- data.table(SUBJID = rownames(scores), scores)
tab <- merge(cov, tab, by = "SUBJID")
tab <- merge(tab, data.table(SUBJID = rownames(neg), neg_hypoxia = neg[, 1], neg_immediate_early = neg[, 2]), by = "SUBJID")
fwrite(tab[sealed == FALSE], file.path(out, "targets_development.tsv"), sep = "\t")
fwrite(tab[sealed == TRUE], file.path(out, "sealed", "targets_sealed.tsv"), sep = "\t")
Sys.chmod(file.path(out, "sealed"), mode = "0700")
fwrite(ceil[order(-ceiling)], file.path(out, "program_ceilings.tsv"), sep = "\t")
fwrite(s$coverage, file.path(out, "program_coverage.tsv"), sep = "\t")
cat("donors", nrow(tab), "| development", sum(!tab$sealed), "| programs scored", sum(colSums(is.finite(scores)) > 0),
    "| primary family (ceiling >= 0.6):", ceil[ceiling >= 0.6, .N], "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_c3_targets.txt"))
