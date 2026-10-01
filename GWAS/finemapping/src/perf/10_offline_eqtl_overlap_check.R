#!/usr/bin/env Rscript
# 10_offline_eqtl_overlap_check.R — independent check of COLOC signal eligibility
# (review item 2), without refitting SuSiE.
#
# For every tier-1/2 SuSiE-positive gene-study pair replayed by the ATAC v3
# candidate, the replay exported each coloc.susie() signal pair (signal_pairs.tsv:
# hit1, hit2, PP.H4.abf, idx1, idx2) and the SNPs coloc used
# (variant_posteriors.tsv.gz `snp`). Those runs cut both fits to the shared SNPs
# before coloc, so coloc's overlap.min check could not fail. Here the saved eQTL
# fit is reloaded and, for eQTL signal idx2, the share of its FULL posterior on
# the shared SNPs is recomputed with coloc 5.2.3's own convention
# (logbf_to_pp(..., last_is_null = TRUE), via src/coloc_signal_eligibility.R).
#
# GWAS side: 06_susie_coloc.R restricts the GWAS SNPs to eQTL positions before
# fine-mapping (gwas_sub) and names them CHR:POS, so every GWAS SNP is an eQTL
# column and prop_gwas = 1 by construction. This script asserts the part it can
# see offline: the exported shared set is a subset of the eQTL columns.
#
# Env:
#   OUT_TSV        (required) new file; refuses to overwrite
#   REPLAY_ROOT    default: the atac-context-v3-candidate-2026-08-11-r1 replay
#   EQTL_SUSIE_DIR default: results/eqtl_susie_polyfun
#   ONLY_GWAS / ONLY_ENSEMBL  optional, restrict to one study / gene (smoke runs)

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
source(file.path(FM_DIR, "src/coloc_signal_eligibility.R"))

REPLAY_ROOT <- Sys.getenv("REPLAY_ROOT", unset = file.path(BASE_DIR,
  "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1",
  "genetics/replay_execution"))
EQTL_SUSIE_DIR <- Sys.getenv("EQTL_SUSIE_DIR", unset = file.path(FM_DIR, "results/eqtl_susie_polyfun"))
OUT_TSV <- Sys.getenv("OUT_TSV", unset = "")
ONLY_GWAS <- Sys.getenv("ONLY_GWAS", unset = "")
ONLY_ENSEMBL <- Sys.getenv("ONLY_ENSEMBL", unset = "")
if (OUT_TSV == "") stop("OUT_TSV is required")
if (file.exists(OUT_TSV)) stop("Refusing to overwrite ", OUT_TSV)

cat("coloc", as.character(packageVersion("coloc")), "| replay:", REPLAY_ROOT,
    "| eQTL fits:", EQTL_SUSIE_DIR, "\n")
P2 <- eval(formals(coloc:::coloc.bf_bf)$p2)

pair_files <- list.files(REPLAY_ROOT, pattern = "^signal_pairs\\.tsv$",
                         recursive = TRUE, full.names = TRUE)
pair_files <- pair_files[grepl("/exports/", pair_files)]
if (ONLY_GWAS != "")    pair_files <- pair_files[grepl(paste0("/exports/", ONLY_GWAS, "/"), pair_files, fixed = TRUE)]
if (ONLY_ENSEMBL != "") pair_files <- pair_files[grepl(paste0("/", ONLY_ENSEMBL, "/"), pair_files, fixed = TRUE)]
if (!length(pair_files)) stop("No replay signal_pairs.tsv found")
cat("Gene-study pairs:", length(pair_files), "\n")

out <- vector("list", length(pair_files))
for (k in seq_along(pair_files)) {
  gene_dir <- dirname(pair_files[k])
  sig <- fread(pair_files[k], colClasses = list(character = c("gwas_name", "ensembl")))
  snps <- unique(fread(file.path(gene_dir, "variant_posteriors.tsv.gz"), select = "snp")$snp)
  gene_id <- sig$ensembl[1]
  chr_num <- sig$chr[1]
  if (!all(sig$nsnps == length(snps))) stop("nsnps differs from the exported SNP set: ", gene_dir)

  s_eqtl <- readRDS(file.path(EQTL_SUSIE_DIR, paste0("chr", chr_num),
                              paste0(gene_id, "_susie.rds")))
  if (!isTRUE(s_eqtl$converged)) stop("non-converged eQTL fit: ", gene_id)
  eqtl_cols <- colnames(s_eqtl$lbf_variable)
  if (!all(sig$idx2 %in% s_eqtl$sets$cs_index)) stop("idx2 is not an eQTL credible set: ", gene_dir)
  bf <- s_eqtl$lbf_variable[sig$idx2, , drop = FALSE]
  shared <- intersect(eqtl_cols, snps)

  out[[k]] <- data.table(
    batch = sub(".*/(batch_[0-9]+)/exports/.*", "\\1", gene_dir),
    gwas_name = sig$gwas_name, gene = sig$gene, ensembl = gene_id, chr = chr_num,
    idx1 = sig$idx1, idx2 = sig$idx2, hit1 = sig$hit1, hit2 = sig$hit2,
    replay_PP.H4 = sig$PP.H4.abf,
    n_shared_snps = length(snps), n_eqtl_snps = length(eqtl_cols),
    gwas_set_subset_of_eqtl = all(snps %in% eqtl_cols),
    eqtl_last_snp_shared = tail(eqtl_cols, 1) %in% snps,
    prop_eqtl = signal_shared_prop(bf, shared, P2)$prop,
    prop_eqtl_uniform = signal_shared_prop_uniform(bf, shared)
  )
  if (k %% 100 == 0) cat("  ", k, "/", length(pair_files), "\n")
}
res <- rbindlist(out)
res[, eligible := prop_eqtl >= COLOC_OVERLAP_MIN & gwas_set_subset_of_eqtl]
res[, eligible_uniform := prop_eqtl_uniform >= COLOC_OVERLAP_MIN & gwas_set_subset_of_eqtl]

if (!all(res$gwas_set_subset_of_eqtl)) {
  warning(sum(!res$gwas_set_subset_of_eqtl), " signal pairs whose shared SNP set is not within the eQTL columns")
}
dir.create(dirname(OUT_TSV), recursive = TRUE, showWarnings = FALSE)
fwrite(res, OUT_TSV, sep = "\t")

# Gene-study level: the adopted call is the max PP.H4 over all pairs; the
# corrected call uses eligible pairs only.
gene <- res[, .(
  n_pairs = .N,
  n_eligible = sum(eligible),
  adopted_pp4 = max(replay_PP.H4),
  corrected_pp4 = if (any(eligible)) max(replay_PP.H4[eligible]) else NA_real_
), by = .(gwas_name, gene, ensembl, chr)]
plan_file <- file.path(dirname(REPLAY_ROOT), "prepared/replay_plan.tsv")
if (file.exists(plan_file)) {
  plan <- fread(plan_file, select = c("gwas_name", "ensembl", "promoted_pp_h4_susie"))
  gene <- merge(gene, plan, by = c("gwas_name", "ensembl"), all.x = TRUE, sort = FALSE)
  cat("Replay max PP.H4 vs adopted PP.H4.susie, max |diff|:",
      max(abs(gene$adopted_pp4 - gene$promoted_pp_h4_susie), na.rm = TRUE), "\n")
}
gene[, status := fifelse(n_eligible == 0, "untestable_insufficient_shared_posterior",
                  fifelse(corrected_pp4 > 0.5, "retained_pp4_gt_0.5", "eligible_pp4_le_0.5"))]
fwrite(gene, sub("\\.tsv$", "_gene_study.tsv", OUT_TSV), sep = "\t")

cat("\nSignal pairs:", nrow(res), "| eligible:", sum(res$eligible),
    "| ineligible:", sum(!res$eligible),
    "| coloc vs uniform-prior call differ:", sum(res$eligible != res$eligible_uniform), "\n")
cat("Gene-study pairs:", nrow(gene), "| genes:", uniqueN(gene$ensembl), "\n")
print(gene[, .N, by = status])
cat("Genes still SuSiE PP4 > 0.5 in >=1 tier-1/2 study:",
    uniqueN(gene[status == "retained_pp4_gt_0.5", ensembl]), "of", uniqueN(gene$ensembl), "\n")
writeLines(capture.output(sessionInfo()), sub("\\.tsv$", "_sessionInfo.txt", OUT_TSV))
cat("Wrote", OUT_TSV, "\n")
