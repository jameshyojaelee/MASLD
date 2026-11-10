#!/usr/bin/env Rscript
# U9 independent re-derivation of Figure 4E (code review 2026-10-23).
# Written from METHODS.md (Figure 4E paragraph) plus the bulk program-score
# definition (per-cohort gene z-scores, L1-weighted mean of members, >= 80% of
# program L1 weight observed, score re-standardized within cohort).
# Deterministic: no random numbers are drawn. Seed set for the record only.
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(metafor)
})
set.seed(20260923L)
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out <- Sys.getenv("U9_OUT"); stopifnot(nzchar(out)); dir.create(out, recursive = TRUE, showWarnings = FALSE)
P <- file.path(R, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
inp <- c(
  dge = file.path(P, "merged_dge.rds"),
  meta = file.path(P, "meta_matched.rds"),
  ann = file.path(R, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
  registry = file.path(R, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv"),
  membership = file.path(R, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"),
  comp_rebuilt = file.path(R, "figures/candidates/composition-774-20261011T151616Z/composition/accepted_composition_774_plus_unchanged4.tsv.gz"),
  comp_earlier = file.path(R, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/sealed_inputs/nonholdout_composition.tsv.gz"),
  producer_meta = file.path(R, "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_three_cohorts/meta_effects.tsv"),
  producer_roster = file.path(R, "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_three_cohorts/model_roster.tsv")
)
fwrite(data.table(input = names(inp), path = inp,
                  sha256 = sapply(inp, function(f) sub(" .*", "", system2("sha256sum", f, stdout = TRUE)))),
       file.path(out, "input_sha256.tsv"), sep = "\t")

dge <- readRDS(inp[["dge"]])
meta <- as.data.table(readRDS(inp[["meta"]]))
meta <- meta[match(colnames(dge), sample_id)]
stopifnot(identical(meta$sample_id, colnames(dge)), nrow(meta) == 774L)
cat("DGE rows", nrow(dge), "cols", ncol(dge), "norm.factors range", range(dge$samples$norm.factors), "\n")

# ---- symbol matrix: logCPM (prior 1, stored normalization), mean over Ensembl rows per symbol
lcpm <- cpm(dge, log = TRUE, prior.count = 1)
ann <- fread(inp[["ann"]], select = c("gene_id", "gene_name"))
ann[, key := sub("\\..*$", "", gene_id)][, sym := toupper(trimws(gene_name))]
ann <- unique(ann[!is.na(sym) & sym != "", .(key, sym)])
stopifnot(ann[, uniqueN(sym), by = key][V1 > 1, .N] == 0L)
sym <- ann$sym[match(sub("\\..*$", "", rownames(lcpm)), ann$key)]
ok <- !is.na(sym)
cat("Ensembl rows mapped to symbol", sum(ok), "of", length(ok), "\n")
S <- rowsum(lcpm[ok, ], sym[ok])
S <- S / as.numeric(table(sym[ok])[rownames(S)])
rm(lcpm); invisible(gc())
cat("symbols", nrow(S), "\n")

# ---- programs
reg <- fread(inp[["registry"]])
mem <- fread(inp[["membership"]])
stopifnot(nrow(reg) == 117L)
tot <- mem[, .(total = sum(as.numeric(original_l1_weight))), by = program_uid]
# Membership A ("mapped", the rule the producer uses): mapped_symbol only; rows whose
# mapped_symbol is NA (unmapped_or_ambiguous) or an Ensembl ID add weight to the total
# but are never scored.
memA <- mem[!is.na(mapped_symbol) & mapped_symbol != "",
            .(w = sum(as.numeric(original_l1_weight))), by = .(program_uid, s = toupper(trimws(mapped_symbol)))]
# Membership B ("primary", the Figure 5E rule): also score unmapped_or_ambiguous members
# whose symbol names exactly one primary-assembly GENCODE v49 gene, and Ensembl-ID members
# through their GENCODE v49 gene name.
gfull <- fread(inp[["ann"]], select = c("gene_id", "gene_name", "chromosome"))
nprim <- gfull[chromosome %in% paste0("chr", c(1:22, "X", "Y", "M")), .N, by = gene_name]
mem[, sB := fifelse(mapped_symbol_status == "gencode_v49_unique_symbol_confirmed", mapped_symbol,
             fifelse(mapped_symbol_status == "gencode_v49_unambiguous_ensembl_to_symbol",
                     gfull$gene_name[match(sub("\\..*$", "", mapped_symbol), sub("\\..*$", "", gfull$gene_id))],
                     fifelse(canonical_gene %in% nprim[N == 1L, gene_name], canonical_gene, NA_character_)))]
cat("membership B: rows scored", mem[!is.na(sB), .N], "of", nrow(mem), "; rows added over A", mem[!is.na(sB) & (is.na(mapped_symbol) | mapped_symbol_status != "gencode_v49_unique_symbol_confirmed"), .N], "\n")
memB <- mem[!is.na(sB) & sB != "", .(w = sum(as.numeric(original_l1_weight))), by = .(program_uid, s = toupper(trimws(sB)))]

zrows <- function(M) { mu <- rowMeans(M); sdv <- sqrt(rowSums((M - mu)^2) / (ncol(M) - 1)); (M - mu) / sdv }
Zc <- lapply(split(seq_len(nrow(meta)), meta$dataset), function(j) zrows(S[, j, drop = FALSE]))
score_with <- function(mem2, label) {
  cov <- merge(tot, mem2[s %in% rownames(S), .(obs = sum(w)), by = program_uid], all.x = TRUE)
  cov[is.na(obs), obs := 0][, coverage := obs / total][, eligible := coverage >= 0.80]
  cov <- merge(cov, reg[, .(program_uid, cell_type, module_name)], by = "program_uid")
  cat("[", label, "] eligible programs", sum(cov$eligible), "\n"); print(cov[eligible == TRUE, .N, by = cell_type])
  fwrite(cov, file.path(out, paste0("program_coverage_", label, ".tsv")), sep = "\t")
  progs <- cov[eligible == TRUE, program_uid]
  score <- matrix(NA_real_, length(progs), ncol(S), dimnames = list(progs, colnames(S)))
  cc <- list()
  for (ds in names(Zc)) {
    j <- which(meta$dataset == ds); Z <- Zc[[ds]]
    var_ok <- rownames(Z)[apply(is.finite(Z), 1, all)]
    for (p in progs) {
      m <- mem2[program_uid == p & s %in% var_ok]
      ret <- sum(m$w) / cov[program_uid == p, total]
      cc[[length(cc) + 1]] <- data.table(program_uid = p, dataset = ds, retained = ret)
      if (ret < 0.80) next
      sc <- colSums(Z[m$s, , drop = FALSE] * (m$w / sum(m$w)))
      score[p, j] <- (sc - mean(sc)) / sd(sc)
    }
  }
  fwrite(rbindlist(cc), file.path(out, paste0("cohort_coverage_", label, ".tsv")), sep = "\t")
  list(score = score, cov = cov)
}
SC <- list(mapped = score_with(memA, "mapped"), primary = score_with(memB, "primary"))

# ---- composition CLR PCs (within cohort, over all 774-unit cohort records)
clr_pcs <- function(path) {
  comp <- fread(path, check.names = FALSE)
  lin <- setdiff(names(comp)[vapply(comp, is.numeric, logical(1))], c("sample_id", "dataset"))
  stopifnot(length(lin) == 16L)
  comp <- comp[match(meta$sample_id, sample_id)]
  cat(basename(path), ": 774 records matched", sum(!is.na(comp$sample_id)), "\n")
  Fm <- as.matrix(comp[, ..lin]); avail <- rowSums(is.finite(Fm)) == 16L
  Fm[Fm == 0] <- 1e-6; Fm <- Fm / rowSums(Fm); L <- log(Fm); L <- L - rowMeans(L)
  pcs <- matrix(NA_real_, nrow(meta), 3, dimnames = list(meta$sample_id, paste0("pc", 1:3)))
  for (ds in unique(meta$dataset)) {
    j <- which(meta$dataset == ds & avail)
    B <- L[j, , drop = FALSE]; B <- B[, apply(B, 2, sd) > 0, drop = FALSE]
    x <- prcomp(B, center = TRUE, scale. = TRUE)$x[, 1:3]
    pcs[j, ] <- apply(x, 2, function(v) (v - mean(v)) / sd(v))
  }
  pcs
}
PCS <- list(rebuilt = clr_pcs(inp[["comp_rebuilt"]]), earlier = clr_pcs(inp[["comp_earlier"]]))

# ---- roster
cohorts <- c("GSE130970", "GSE135251", "GSE162694")
ros <- meta[dataset %in% cohorts & !is.na(fibrosis_stage) & !is.na(nas_score)]
cat("roster", nrow(ros), "\n"); print(ros[, .N, by = dataset])
prod_ros <- fread(inp[["producer_roster"]])
cat("roster identical to producer roster:", setequal(prod_ros$sample_id, ros$sample_id), "\n")
sexdiff <- ros[sex_final != inferred_sex, .(sample_id, dataset, sex_final, inferred_sex, fibrosis_stage, nas_score)]
fwrite(sexdiff, file.path(out, "sex_disagreements_roster.tsv"), sep = "\t")
cat("sex_final != inferred_sex in roster:", nrow(sexdiff), "; in 774:", meta[sex_final != inferred_sex, .N], "\n")

# ---- per-cohort OLS + HC3 (vectorized over programs)
fit_cohort <- function(ds, sexvar, pcs = NULL, score) {
  d <- ros[dataset == ds]
  X <- cbind(1, male = as.numeric(d[[sexvar]] == "M"),
             fib = d$fibrosis_stage - mean(d$fibrosis_stage),
             nas = d$nas_score - mean(d$nas_score))
  if (!is.null(pcs)) X <- cbind(X, pcs[d$sample_id, , drop = FALSE])
  stopifnot(qr(X)$rank == ncol(X), all(is.finite(X)))
  Y <- t(score[, d$sample_id, drop = FALSE])           # n x programs
  XtXi <- solve(crossprod(X)); B <- XtXi %*% crossprod(X, Y)
  E <- Y - X %*% B; h <- rowSums((X %*% XtXi) * X)
  res <- rbindlist(lapply(c(fib = 3L, nas = 4L), function(k) {
    a <- (X %*% XtXi)[, k]                              # row k of (X'X)^-1 X'
    se <- sqrt(colSums((a^2) * (E^2) / (1 - h)^2))
    data.table(program_uid = colnames(Y), axis = colnames(X)[k], beta = B[k, ], se = se)
  }))
  res[, `:=`(dataset = ds, n = nrow(d), df = nrow(d) - ncol(X))]
  res[!is.finite(beta) | !is.finite(se), `:=`(beta = NA, se = NA)]
  res
}
meta_fit <- function(ce) {
  ce[is.finite(beta) & se > 0, {
    if (.N < 3L) list(k = .N, b = NA_real_, se_kh = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_,
                      dfs = NA_real_, tau2 = NA_real_, p_z = NA_real_) else {
      md <- data.frame(yi = beta, sei = se)
      f <- rma(yi = yi, sei = sei, data = md, method = "REML", test = "knha")
      fz <- rma(yi = yi, sei = sei, data = md, method = "REML", test = "z")
      list(k = .N, b = as.numeric(f$b), se_kh = f$se, lo = f$ci.lb, hi = f$ci.ub, p = f$pval,
           dfs = f$dfs, tau2 = f$tau2, p_z = fz$pval)
    }
  }, by = .(program_uid, axis)]
}
variants <- CJ(membership = c("mapped", "primary"), sexvar = c("sex_final", "inferred_sex"), model = c("unadjusted", "rebuilt", "earlier"), sorted = FALSE)
all_meta <- list(); summ <- list()
for (i in seq_len(nrow(variants))) {
  v <- variants[i]
  pcs <- if (v$model == "unadjusted") NULL else PCS[[v$model]]
  ce <- rbindlist(lapply(cohorts, fit_cohort, sexvar = v$sexvar, pcs = pcs, score = SC[[v$membership]]$score))
  me <- meta_fit(ce)
  me[, `:=`(q = p.adjust(p, "BH"), q_z = p.adjust(p_z, "BH"), family = .N)]
  me[, `:=`(membership = v$membership, sexvar = v$sexvar, model = v$model)]
  all_meta[[i]] <- me
  summ[[i]] <- me[, .(n_programs = .N, family = family[1], k_min = min(k), dfs = paste(unique(dfs), collapse = ","),
                      n_q05 = sum(q < 0.05), min_q = min(q), n_ci_excl0 = sum(lo > 0 | hi < 0),
                      n_q05_z = sum(q_z < 0.05)), by = .(membership, sexvar, model, axis)]
}
all_meta <- rbindlist(all_meta); summ <- rbindlist(summ)
fwrite(all_meta, file.path(out, "meta_effects_all_variants.tsv"), sep = "\t")
fwrite(summ, file.path(out, "summary.tsv"), sep = "\t")
print(summ)

# ---- compare with producer (sex_final arms) and sex sensitivity
pm <- fread(inp[["producer_meta"]])
pm[, model := fifelse(model_kind == "primary", "unadjusted", "rebuilt")]
pm[, axis := fifelse(axis == "fibrosis", "fib", "nas")]
cmp <- merge(all_meta[sexvar == "sex_final" & membership == "mapped"], pm[, .(program_uid = feature_id, axis, model, beta_meta, ci_lower, ci_upper, q_value)],
             by = c("program_uid", "axis", "model"))
cat("producer comparison rows", nrow(cmp), "\n")
print(cmp[, .(max_abs_beta = max(abs(b - beta_meta)), max_abs_q = max(abs(q - q_value)),
              max_abs_lo = max(abs(lo - ci_lower))), by = .(model, axis)])
w <- dcast(all_meta, membership + program_uid + axis + model ~ sexvar, value.var = c("b", "q", "lo", "hi"))
sexcmp <- w[, .(max_abs_dbeta = max(abs(b_inferred_sex - b_sex_final)),
                max_abs_dq = max(abs(q_inferred_sex - q_sex_final)),
                calls_sf = sum(q_sex_final < 0.05), calls_is = sum(q_inferred_sex < 0.05),
                ci_excl_sf = sum(lo_sex_final > 0 | hi_sex_final < 0),
                ci_excl_is = sum(lo_inferred_sex > 0 | hi_inferred_sex < 0),
                ci_flips = sum((lo_sex_final > 0 | hi_sex_final < 0) != (lo_inferred_sex > 0 | hi_inferred_sex < 0)),
                min_q_sf = min(q_sex_final), min_q_is = min(q_inferred_sex)), by = .(membership, model, axis)]
fwrite(sexcmp, file.path(out, "sex_variable_comparison.tsv"), sep = "\t")
print(sexcmp)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("U9_DONE\n")
