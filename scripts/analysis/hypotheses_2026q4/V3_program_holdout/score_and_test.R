#!/usr/bin/env Rscript
# V3 step 2: score the frozen Stromal ECM and Ductular injury programs in
# GSE256398 donor pseudobulk with the unchanged transport scoring function,
# fit the prespecified HC3 models, and build the expression-matched random-set
# null. Follows PRESPEC_V3.md.
# Usage: Rscript score_and_test.R <RESULT_DIR>

suppressPackageStartupMessages({library(data.table); library(edgeR)})
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
res <- args[[1]]
pb_dir <- file.path(res, "pseudobulk")
out <- file.path(res, "results")
if (dir.exists(out)) stop("results directory exists; refusing to overwrite")
if (!file.exists(file.path(res, "PRESPEC_V3.md"))) stop("PRESPEC_V3.md missing")
dir.create(out)

lib_dir <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/manuscript/program_observability_map"
source(file.path(lib_dir, "config.R"))
source(file.path(lib_dir, "analysis_lib.R"))

SEED <- 20261007L
N_RANDOM <- 1000L
N_BINS <- 20L
REGISTRY_SHA <- "b134b5f8e46a716e857af22271aa54e1294432b9f9ab735559d7f603cdd29fe7"
MEMBERSHIP_SHA <- "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b"
PROGRAMS <- c(stromal_ecm = "hotspot_hepatocytes_f05c535ae5bbc0b9",
              ductular_injury = "hotspot_hepatocytes_48f39dd4d817a10e")
CONTENT_LINEAGE <- c(stromal_ecm = "stellate_fibroblast", ductular_injury = "cholangiocyte")

if (sha256_file(program_registry) != REGISTRY_SHA) stop("registry hash changed")
if (sha256_file(program_membership) != MEMBERSHIP_SHA) stop("membership hash changed")
registry <- fread(program_registry)[program_uid %in% PROGRAMS]
membership <- fread(program_membership)[program_uid %in% PROGRAMS]
if (nrow(registry) != 2L) stop("program registry rows missing")
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))

meta <- fread(file.path(pb_dir, "donor_metadata.tsv"))
nuc <- fread(file.path(pb_dir, "nuclei_per_donor_lineage.tsv"))
setnames(nuc, 1, "donor")
meta[, sex := factor(sex, levels = c("female", "male"))]
meta[, dataset := "GSE256398"]

symbol_logcpm <- function(lineage, donors) {
  pb <- fread(file.path(pb_dir, sprintf("pseudobulk_%s.tsv.gz", lineage)))
  counts <- as.matrix(pb[, ..donors]); rownames(counts) <- pb$gene_id
  y <- calcNormFactors(DGEList(counts))
  lc <- cpm(y, log = TRUE, prior.count = 1)
  collapse_symbols(lc, rownames(lc), annotation)
}

fit_contrasts <- function(scores, sm, prefix_cols = list()) {
  rows <- list()
  for (feature in rownames(scores)) {
    d <- copy(sm); d[, score := scores[feature, donor]]
    specs <- list(
      C1_stage_ordinal_noncirrhotic = list(data = d[stage != "MASH_cirrhosis"],
                                           f = score ~ stage_ordinal + sex + age, term = "stage_ordinal",
                                           groups = "stage"),
      C2_disease_vs_healthy = list(data = d, f = score ~ disease + sex + age, term = "disease",
                                   groups = "disease"),
      C1_ols_nocov = list(data = d[stage != "MASH_cirrhosis"], f = score ~ stage_ordinal,
                          term = "stage_ordinal", groups = "stage"),
      C2_ols_nocov = list(data = d, f = score ~ disease, term = "disease", groups = "disease"))
    for (nm in names(specs)) {
      s <- specs[[nm]]
      grp <- table(s$data[[s$groups]])
      testable <- length(grp) >= 2L && all(grp >= 3L)
      fit <- if (testable) hc3_fit(NULL, s$data, s$f) else list(estimable = FALSE)
      if (isTRUE(fit$estimable)) {
        b <- fit$coefficients[[s$term]]; se <- fit$se[[s$term]]
        q <- qt(0.975, fit$df)
        ols <- summary(lm(s$f, data = s$data))$coefficients[s$term, ]
        rows[[length(rows) + 1L]] <- data.table(
          feature_id = feature, contrast = nm, n_donors = fit$n, df = fit$df,
          group_sizes = paste(names(grp), grp, sep = "=", collapse = ";"),
          beta = b, hc3_se = se, hc3_ci_low = b - q * se, hc3_ci_high = b + q * se,
          hc3_t = fit$statistic[[s$term]], hc3_p = fit$p_value[[s$term]],
          ols_se = ols[["Std. Error"]], ols_p = ols[["Pr(>|t|)"]],
          max_leverage = fit$max_leverage)
      } else {
        rows[[length(rows) + 1L]] <- data.table(
          feature_id = feature, contrast = nm, n_donors = nrow(s$data),
          group_sizes = paste(names(grp), grp, sep = "=", collapse = ";"),
          failure = if (!testable) "fewer_than_3_donors_in_a_group" else fit$failure_reason)
      }
    }
  }
  rbindlist(rows, fill = TRUE)
}

set.seed(SEED)

# ---------------------------------------------------------------- primary: hepatocytes
hep_ok <- nuc[hepatocyte >= 50L, donor]
sm <- meta[donor %in% hep_ok]
message("donors with >=50 hepatocyte nuclei: ", nrow(sm), " of ", nrow(meta))
S <- symbol_logcpm("hepatocyte", sm$donor)
scored <- score_programs(S, sm, membership, registry, program_weight_coverage)
fwrite(scored$coverage, file.path(out, "program_coverage_hepatocyte.tsv"), sep = "\t")
prim_scores <- scored$primary
fwrite(data.table(donor = colnames(prim_scores), t(prim_scores))[sm, on = "donor"],
       file.path(out, "donor_scores_hepatocyte.tsv"), sep = "\t")

res_all <- fit_contrasts(prim_scores, sm)
res_all[, program := names(PROGRAMS)[match(feature_id, PROGRAMS)]]
res_all[, lineage := "hepatocyte"]
fam <- res_all[contrast %in% c("C1_stage_ordinal_noncirrhotic", "C2_disease_vs_healthy")]
fam[, holm_p := p.adjust(hc3_p, method = "holm")]
fwrite(fam, file.path(out, "v3_primary_family.tsv"), sep = "\t")

# ------------------------------------------------- negative control: random gene sets
Z <- zscore_rows(S)
variable <- rownames(Z)[rowSums(is.finite(Z)) == ncol(Z)]
mean_expr <- rowMeans(S[variable, , drop = FALSE])
program_genes <- unique(toupper(membership$mapped_symbol))
pool <- setdiff(variable, program_genes)
brks <- unique(quantile(mean_expr, probs = seq(0, 1, length.out = N_BINS + 1L)))
bin_of <- setNames(cut(mean_expr, brks, include.lowest = TRUE, labels = FALSE), variable)

rand_membership <- list(); rand_registry <- list(); rand_parent <- character()
for (p in names(PROGRAMS)) {
  mm <- membership[program_uid == PROGRAMS[[p]] & mapped_symbol %in% variable,
                   .(weight = sum(original_l1_weight)), by = mapped_symbol]
  for (i in seq_len(N_RANDOM)) {
    pick <- character(nrow(mm))
    used <- character()
    for (g in seq_len(nrow(mm))) {
      cand <- setdiff(pool[bin_of[pool] == bin_of[[mm$mapped_symbol[g]]]], used)
      pick[g] <- cand[sample.int(length(cand), 1L)]
      used <- c(used, pick[g])
    }
    uid <- sprintf("random_%s_%04d", p, i)
    rand_membership[[length(rand_membership) + 1L]] <- data.table(
      program_uid = uid, mapped_symbol = pick, original_l1_weight = mm$weight)
    rand_registry[[length(rand_registry) + 1L]] <- data.table(program_uid = uid)
    rand_parent[uid] <- p
  }
}
rand_membership <- rbindlist(rand_membership); rand_registry <- rbindlist(rand_registry)
rs <- score_programs(S, sm, rand_membership, rand_registry, program_weight_coverage)$primary
rand_res <- fit_contrasts(rs, sm)[contrast %in% c("C1_stage_ordinal_noncirrhotic", "C2_disease_vs_healthy")]
rand_res[, program := rand_parent[feature_id]]
fwrite(rand_res, file.path(out, "v3_random_set_null.tsv.gz"), sep = "\t")
null_summary <- rand_res[, .(n_random_estimable = sum(is.finite(hc3_t))), by = .(program, contrast)]
null_summary <- merge(null_summary, fam[, .(program, contrast, observed_t = hc3_t, observed_beta = beta)],
                      by = c("program", "contrast"))
null_summary[, frac_random_t_ge_observed := mapply(function(p, c, t)
  mean(rand_res[program == p & contrast == c & is.finite(hc3_t), hc3_t] >= t), program, contrast, observed_t)]
null_summary[, frac_random_abs_t_ge_observed := mapply(function(p, c, t)
  mean(abs(rand_res[program == p & contrast == c & is.finite(hc3_t), hc3_t]) >= abs(t)), program, contrast, observed_t)]
null_summary[, random_t_median := mapply(function(p, c)
  median(rand_res[program == p & contrast == c, hc3_t], na.rm = TRUE), program, contrast)]
fwrite(null_summary, file.path(out, "v3_random_set_summary.tsv"), sep = "\t")

# ------------------------------------------- secondary: content lineage, sensitivity
sec <- list(res_all[contrast %in% c("C1_ols_nocov", "C2_ols_nocov")])
for (p in names(PROGRAMS)) {
  lin <- CONTENT_LINEAGE[[p]]
  ok <- nuc[get(lin) >= 20L, donor]
  sml <- meta[donor %in% ok]
  if (nrow(sml) < 6L) next
  Sl <- symbol_logcpm(lin, sml$donor)
  sl <- score_programs(Sl, sml, membership, registry[program_uid == PROGRAMS[[p]]],
                       program_weight_coverage)
  r <- fit_contrasts(sl$primary, sml)
  r[, `:=`(program = p, lineage = lin)]
  fwrite(sl$coverage, file.path(out, sprintf("program_coverage_%s.tsv", lin)), sep = "\t")
  sec[[length(sec) + 1L]] <- r
}
fwrite(rbindlist(sec, fill = TRUE), file.path(out, "v3_secondary_descriptive.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("done")
print(fam)
print(null_summary)
