#!/usr/bin/env Rscript
# Trio covariance analysis for GNMT-MAT1A-CYP2C19 (20261118 candidate, N10/N18).
#
# Key methodological updates over 02_protein.R:
# 1. Exclude missing Genes and `;`-delimited multi-gene shared peptide groups.
# 2. Subset raw PXD051911 log2 intensities to the 46 MASLD participants FIRST,
#    apply <=50% missingness filter, run limma::normalizeBetweenArrays(..., method="quantile")
#    on those 46 MASLD samples alone, collapse by median, and z-score background genes
#    on those 46 MASLD samples alone.
# 3. Replace duplicate-row bootstrap resampling (replace=TRUE) with continuous
#    Dirichlet(1,...,1) participant weights + weighted midranks + lm.wfit across B=2000
#    draws so rank correlation estimates are free of ~36.8% tie artifacts.

suppressPackageStartupMessages({library(data.table); library(limma); library(ggplot2)})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
args <- commandArgs(TRUE)
out <- if (length(args) >= 1L && nzchar(args[1L])) {
  args[1L]
} else {
  file.path(BASE, "Analysis/Multimodal_Program_Projection/results/proteomics_20261118/trio_covariance")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20261118); setDTthreads(2); pdf.options(useDingbats = FALSE)
write <- function(x, name) fwrite(x, file.path(out, name), sep = "\t", na = "NA")

triad <- c("GNMT", "MAT1A", "CYP2C19")
background_sets <- c("HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_FATTY_ACID_METABOLISM")
writeLines(c(
  "seed=20261118", "bootstrap=2000", "bootstrap_type=Dirichlet_1_1_continuous_weights_tie_free",
  "normalization_scope=46_MASLD_only_quantile_norm_and_z",
  "biological_unit=participant", "primary_family=3_pairs",
  "adjustment_family=12_pair_by_model_tests",
  "background_family=3_pairs_sensitivity_only", "common_complete_cases_across_models=true",
  "histology_covariates=recorded_numeric_grades", "selection_conditioned_example=true",
  paste0("background_sets=", paste(background_sets, collapse = ";"))
), file.path(out, "analysis_definition.txt"))

raw_all <- fread(file.path(BASE, "data/PXD051911/liver_protein_quant.txt"))
md <- fread(file.path(BASE, "data/PXD051911/meta_data.txt"))
samples <- setdiff(names(raw_all), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% samples], by = "liver_proteomics_filename")
stopifnot(nrow(md) == 58L)
md[, `:=`(
  sample_id = liver_proteomics_filename,
  group = ifelse(saf_diagnosis == "No_MASLD", "Control", "MASLD"),
  batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")),
  age = as.numeric(alder), bmi_num = as.numeric(bmi), sex = factor(gender),
  Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)),
  Steatosis = as.numeric(steatosis_score), Ballooning = as.numeric(hepatocellular_ballooning_score),
  Inflammation = as.numeric(lobular_inflammation_score), NAS = as.numeric(nafld_activity_score)
)]
setorder(md, sample_id)
md_masld <- droplevels(md[group == "MASLD"])
stopifnot(nrow(md_masld) == 46L)

# Exclude missing and `;`-delimited shared PGs; subset to 46 MASLD FIRST
raw <- raw_all[!is.na(Genes) & trimws(Genes) != "" & !grepl(";", Genes, fixed = TRUE)]
x <- as.matrix(raw[, md_masld$sample_id, with = FALSE]); storage.mode(x) <- "double"
x[!is.finite(x) | x <= 0] <- NA_real_
x <- log2(x); rownames(x) <- trimws(raw$Genes)
x <- x[rowMeans(is.na(x)) <= 0.50, , drop = FALSE]
x <- normalizeBetweenArrays(x, method = "quantile")
indices <- split(seq_len(nrow(x)), rownames(x))
gene <- t(vapply(indices, function(i) if (length(i) == 1L) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x))))
colnames(gene) <- md_masld$sample_id; gene[!is.finite(gene)] <- NA_real_
stopifnot(all(triad %in% rownames(gene)))
for (g in triad) md_masld[, (g) := gene[g, sample_id]]

gmt <- strsplit(readLines(file.path(BASE, "Analysis/downstream_analysis/pathway_analysis/data/genesets/hallmark.gmt")), "\t")
bg <- setdiff(unique(unlist(lapply(gmt, function(v) if (v[1] %in% background_sets) v[-c(1, 2)] else character()))), triad)
bg <- intersect(bg, rownames(gene)); stopifnot(length(bg) >= 10L)
write(data.table(gene = bg, background = "three_fixed_Hallmark_metabolic_sets_46_masld"), "metabolic_background.tsv")
zz <- t(scale(t(gene[bg, , drop = FALSE])))
md_masld[, metabolic_background := apply(zz, 2, median, na.rm = TRUE)]

cols <- c(triad, "batch", "age", "bmi_num", "sex", "Fibrosis", "NAS", "Steatosis", "Ballooning", "Inflammation", "metabolic_background")
d <- droplevels(md_masld[complete.cases(md_masld[, ..cols])])
stopifnot(!anyDuplicated(d$sample_id), nrow(d) > 20L)
write(d[, c("sample_id", "group", cols), with = FALSE], "common_participants.tsv")

models <- list(
  baseline = c("batch", "age", "bmi_num", "sex"),
  fibrosis = c("batch", "age", "bmi_num", "sex", "Fibrosis"),
  nas = c("batch", "age", "bmi_num", "sex", "NAS"),
  components = c("batch", "age", "bmi_num", "sex", "Fibrosis", "Steatosis", "Ballooning", "Inflammation"),
  metabolic_background = c("batch", "age", "bmi_num", "sex", "metabolic_background")
)
pairs <- combn(triad, 2, simplify = FALSE)

w_midrank <- function(x_vec, w) {
  o <- order(x_vec)
  xs <- x_vec[o]; ws <- w[o]
  rx <- rle(xs)
  ends <- cumsum(rx$lengths)
  starts <- ends - rx$lengths + 1L
  ranks_sorted <- numeric(length(x_vec))
  prev_c <- 0
  for (k in seq_along(rx$lengths)) {
    idx_k <- starts[k]:ends[k]
    bw <- sum(ws[idx_k])
    ranks_sorted[idx_k] <- prev_c + 0.5 * bw
    prev_c <- prev_c + bw
  }
  out <- numeric(length(x_vec))
  out[o] <- ranks_sorted
  out
}

design_mats <- lapply(models, function(vars) model.matrix(reformulate(vars), data = d))
calc_w <- function(X_mat, pair, w = rep(1, nrow(d))) {
  r1 <- w_midrank(d[[pair[1L]]], w)
  r2 <- w_midrank(d[[pair[2L]]], w)
  e1 <- lm.wfit(X_mat, r1, w)$residuals
  e2 <- lm.wfit(X_mat, r2, w)$residuals
  v1 <- sum(w * e1^2); v2 <- sum(w * e2^2)
  if (!is.finite(v1) || !is.finite(v2) || v1 <= 1e-15 || v2 <= 1e-15) return(c(rho = NA_real_, p = NA_real_, df = NA_real_))
  r <- sum(w * e1 * e2) / sqrt(v1 * v2)
  df <- nrow(X_mat) - ncol(X_mat) - 1L
  c(rho = r, p = 2 * pt(-abs(r * sqrt(df / pmax(1e-15, 1 - r * r))), df), df = df)
}

observed <- rbindlist(lapply(names(models), function(model) rbindlist(lapply(pairs, function(pair) {
  v <- calc_w(design_mats[[model]], pair)
  data.table(pair = paste(pair, collapse = "–"), model = model, n = nrow(d), rho = v["rho"], p = v["p"], df = v["df"])
}))))

N_BOOT <- 2000L
draws <- array(NA_real_, dim = c(N_BOOT, length(models), 3L))
n_d <- nrow(d)
for (b in seq_len(N_BOOT)) {
  g <- rexp(n_d, rate = 1)
  w <- n_d * g / sum(g)
  for (j in seq_along(models)) {
    for (k in seq_along(pairs)) {
      draws[b, j, k] <- calc_w(design_mats[[j]], pairs[[k]], w = w)["rho"]
    }
  }
}

for (j in seq_along(models)) {
  for (k in seq_along(pairs)) {
    v <- draws[, j, k]
    delta <- v - draws[, 1L, k]
    pair_label <- paste(pairs[[k]], collapse = "–")
    idx <- which(observed$model == names(models)[j] & observed$pair == pair_label)
    observed[idx, `:=`(
      low = quantile(v, 0.025, na.rm = TRUE),
      high = quantile(v, 0.975, na.rm = TRUE),
      n_bootstrap = sum(is.finite(v)),
      delta_vs_baseline = rho - observed[model == "baseline" & pair == pair_label, rho],
      delta_low = quantile(delta, 0.025, na.rm = TRUE),
      delta_high = quantile(delta, 0.975, na.rm = TRUE)
    )]
  }
}
observed[, q_three_pairs := p.adjust(p, "BH", n = 3), by = model]
observed[model != "metabolic_background", q_twelve := p.adjust(p, "BH", n = 12)]
observed[, `:=`(seed = 20261118L, interval = "pointwise_dirichlet_weighted_participant_bootstrap", selection_conditioned = TRUE)]
write(observed, "protein_pair_covariance.tsv")
write(data.table(metric = c("all_liver_participants", "masld_participants_qnorm", "common_MASLD_participants", "background_genes"),
                 value = c(nrow(md), nrow(md_masld), nrow(d), length(bg))), "denominators.tsv")

p <- ggplot(observed[model != "metabolic_background"], aes(rho, pair, color = model, shape = model)) +
  geom_vline(xintercept = 0, color = "#9E9E9E", linewidth = 0.3) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0.15, position = position_dodge(width = 0.6), linewidth = 0.3) +
  geom_point(position = position_dodge(width = 0.6), size = 1) +
  labs(x = "Adjusted protein rank correlation (95% Dirichlet bootstrap interval)", y = NULL, color = NULL, shape = NULL) +
  theme_classic(base_size = 6) + theme(text = element_text(size = 6), legend.position = "bottom")
ggsave(file.path(out, "S3_metabolic_protein_covariance.pdf"), p, width = 4.6, height = 2.1, device = pdf, useDingbats = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines("Candidate 20261118: 46-MASLD-only qnorm/z + Dirichlet tie-free rank-OLS bootstrap.", file.path(out, "COMPLETE"))
cat("[02_protein_20261118] complete:\n")
print(observed[, .(pair, model, n, rho, low, high, p, q_twelve)])
