# RV4: share of the 1,854 disease calls present in obese normal liver (GSE48452). Follows PRESPEC_RV4.md.
# Unit: participant (one biopsy each). Microarray: limma on deposited RMA values.
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
DISC <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
GEO <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z/downloads_untrusted/geo_microarray_components")
SEED <- 20261102
NPERM <- 1000L; NBOOT <- 1000L
wr <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t", na = "NA")

disc <- fread(DISC)[padj < 0.05 & abs(logFC) > 0.5]
stopifnot(nrow(disc) == 1854)

# ---------------- series matrix ----------------
lines <- readLines(gzfile(file.path(GEO, "GSE48452_series_matrix.txt.gz")))
field <- function(prefix) {
  hits <- grep(paste0("^", prefix, "\t"), lines, value = TRUE)
  lapply(hits, function(h) gsub('"', "", strsplit(h, "\t")[[1]][-1]))
}
gsm <- field("!Sample_geo_accession")[[1]]
ch <- field("!Sample_characteristics_ch1")
char <- function(key) {
  for (v in ch) if (all(startsWith(v, paste0(key, ": ")))) return(sub(paste0("^", key, ": "), "", v))
  stop("missing characteristic ", key)
}
meta <- data.table(gsm = gsm, group = char("group"), surgery = char("bariatric surgery"),
                   sex = char("Sex"), age = as.numeric(char("age")), bmi = as.numeric(char("bmi")))
b <- grep("^!series_matrix_table_begin", lines) + 1L; e <- grep("^!series_matrix_table_end", lines) - 1L
tab <- fread(text = lines[b:e])
expr <- as.matrix(tab[, -1]); rownames(expr) <- as.character(tab$ID_REF)
stopifnot(identical(colnames(expr), meta$gsm))

meta <- meta[surgery != "after surgery"]
meta[, grp := factor(c(Control = "lean", `Healthy obese` = "obese", Steatosis = "steatosis", Nash = "nash")[group],
                     levels = c("lean", "obese", "steatosis", "nash"))]
stopifnot(identical(as.integer(table(meta$grp)), c(12L, 16L, 9L, 17L)))
expr <- expr[, meta$gsm]

# ---------------- probe sets to genes ----------------
ann <- fread(file.path(GEO, "GPL11532_annot.txt"), skip = "ID\tGB_LIST", sep = "\t", quote = "", fill = TRUE)
ann <- ann[category == "main" & gene_assignment != "---"]
ann[, symbol := trimws(vapply(strsplit(gene_assignment, " // ", fixed = TRUE), `[`, "", 2))]
ann <- ann[!is.na(symbol) & symbol != "" & as.character(ID) %in% rownames(expr)]
ann[, mean_expr := rowMeans(expr[as.character(ID), , drop = FALSE])]
ann <- ann[order(-mean_expr)][!duplicated(symbol)]
G <- expr[as.character(ann$ID), ]; rownames(G) <- ann$symbol

calls <- disc[symbol %in% rownames(G), .(symbol, d = logFC)]

# ---------------- fits ----------------
contrasts_of <- c(obese = "grpobese - grplean", steatosis = "grpsteatosis - grplean", nash = "grpnash - grplean")
fit_all <- function(m, covars = "sex") {
  X <- model.matrix(as.formula(paste("~ 0 + grp +", paste(covars, collapse = " + "))), m)
  f <- lmFit(G[, m$gsm], X)
  cm <- makeContrasts(contrasts = contrasts_of, levels = X); colnames(cm) <- names(contrasts_of)
  eBayes(contrasts.fit(f, cm), trend = TRUE)
}
slope_of <- function(lfc) { r <- lfc[calls$symbol]; unname(coef(lm(r ~ calls$d))[2]) }

f0 <- fit_all(meta)
obs <- sapply(names(contrasts_of), function(k) slope_of(f0$coefficients[, k]))

# permutation null: shuffle labels among the samples of each contrast
set.seed(SEED)
null <- rbindlist(lapply(names(contrasts_of), function(k) {
  sub <- meta[grp %in% c("lean", k)]
  rbindlist(lapply(seq_len(NPERM), function(i) {
    m <- copy(meta); idx <- which(m$grp %in% c("lean", k))
    m$grp[idx] <- sample(m$grp[idx])
    data.table(contrast = k, perm = i, slope = slope_of(fit_all(m)$coefficients[, k]))
  }))
}))
wr(null, "perm_null.tsv.gz")

prim <- data.table(contrast = names(obs), n_calls_matched = nrow(calls), slope = obs,
                   null_median = null[, median(slope), by = contrast][match(names(obs), contrast), V1],
                   null_q95 = null[, quantile(slope, 0.95), by = contrast][match(names(obs), contrast), V1])
prim[, perm_p := mapply(function(k, o) (1 + null[contrast == k, sum(slope >= o)]) / (NPERM + 1), contrast, slope)]
prim[, holm_p := p.adjust(perm_p, "holm")]
prim[, sign_agree := sapply(contrast, function(k) mean(sign(f0$coefficients[calls$symbol, k]) == sign(calls$d)))]
prim[, sign_agree_up := sapply(contrast, function(k) calls[d > 0, mean(f0$coefficients[symbol, k] > 0)])]
prim[, sign_agree_down := sapply(contrast, function(k) calls[d < 0, mean(f0$coefficients[symbol, k] < 0)])]
wr(prim, "primary.tsv")

# bootstrap of the obese fraction, resampling participants within group
set.seed(SEED)
boot <- rbindlist(lapply(seq_len(NBOOT), function(i) {
  m <- meta[, .SD[sample.int(.N, .N, replace = TRUE)], by = grp]
  m[, gsm_orig := gsm]; m[, gsm := make.unique(gsm)]
  Gb <- G[, m$gsm_orig]; colnames(Gb) <- m$gsm
  X <- model.matrix(~ 0 + grp + sex, m)
  f <- lmFit(Gb, X)
  cm <- makeContrasts(contrasts = contrasts_of, levels = X); colnames(cm) <- names(contrasts_of)
  cf <- contrasts.fit(f, cm)$coefficients
  s <- sapply(names(contrasts_of), function(k) slope_of(cf[, k]))
  data.table(boot = i, obese = s[["obese"]], steatosis = s[["steatosis"]], nash = s[["nash"]])
}))
boot[, obese_fraction := obese / nash]
wr(boot, "bootstrap.tsv.gz")

# per-gene obese vs lean, one-sided in the discovery direction
tt <- f0$t[calls$symbol, "obese"]
p1 <- pt(ifelse(calls$d > 0, tt, -tt), df = f0$df.total[1], lower.tail = FALSE)
gt <- data.table(calls, obese_lfc = f0$coefficients[calls$symbol, "obese"], steatosis_lfc = f0$coefficients[calls$symbol, "steatosis"],
                 nash_lfc = f0$coefficients[calls$symbol, "nash"], obese_p_onesided = p1, obese_q_onesided = p.adjust(p1, "BH"))
wr(gt, "gene_table.tsv.gz")

fa <- fit_all(meta, c("sex", "age"))
sens <- data.table(analysis = "plus_age", contrast = names(contrasts_of),
                   slope = sapply(names(contrasts_of), function(k) slope_of(fa$coefficients[, k])))
wr(sens, "sensitivity.tsv")

sink(file.path(OUT, "summary.txt"))
cat(sprintf("Participants: %s\n", paste(names(table(meta$grp)), table(meta$grp), collapse = ", ")))
cat(sprintf("Array genes %d; discovery calls matched %d of 1854 (up %d, down %d)\n\n", nrow(G), nrow(calls), sum(calls$d > 0), sum(calls$d < 0)))
cat("Primary (Holm over 3)\n"); print(prim)
cat(sprintf("\nObese fraction slope(obese)/slope(nash) = %.3f; bootstrap 95%% interval %.3f to %.3f\n",
            obs[["obese"]] / obs[["nash"]], quantile(boot$obese_fraction, 0.025), quantile(boot$obese_fraction, 0.975)))
cat(sprintf("Steatosis fraction = %.3f\n", obs[["steatosis"]] / obs[["nash"]]))
cat(sprintf("Calls with obese-vs-lean one-sided BH q < 0.05: %d of %d\n", sum(gt$obese_q_onesided < 0.05), nrow(gt)))
cat("\nSensitivity\n"); print(sens)
sink()
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("done\n")
