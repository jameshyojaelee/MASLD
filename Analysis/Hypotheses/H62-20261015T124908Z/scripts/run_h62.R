#!/usr/bin/env Rscript
# H62 run script. Implements PRESPEC_H62.md (sha256 569d25c4fb4e6f0b89a90675fc3bd1c4e26d67d00baccba00105d8d184c63c35).
# Question: is the frozen adult disease-vs-control signature score higher in the 33 adolescent
# NAFLD participants than in the 34 No NAFLD participants of GSE66676 (sex, age, BMI adjusted)?
# Unit of inference: the participant (one biopsy each). Genes and probes are not replicates.
# Family: one confirmatory test (T1), alpha 0.05, no correction. Every other p is nominal.
#
# Modes:
#   default        : final run (real labels, B = 10,000 / 1,000 / 1,000 / 1,000)
#   --debug-null   : code check only. Everything up to and including PC-B and diagnostics 1-5 uses
#                    the real data; then the GSE66676 histology labels and the GSE185051 NAS values
#                    are permuted (seed 1, debug only) so the outcome code runs on null data with small B.
suppressMessages({ library(limma); library(edgeR); library(ggplot2) })
options(stringsAsFactors = FALSE, warn = 1)

## ---------------------------------------------------------------- arguments
a <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, default = NULL) { i <- match(k, a); if (is.na(i)) default else a[i + 1] }
OUT <- getarg("--out"); SEED <- as.integer(getarg("--seed", "20260923"))
DEBUG <- "--debug-null" %in% a
if (is.null(OUT) || !dir.exists(OUT)) stop("--out must be an existing directory")
if (SEED != 20260923L) stop("prespec seed is 20260923")
OUT <- normalizePath(OUT); RUN <- dirname(OUT)
B_FL <- 10000L; B_NC2 <- 1000L; B_NC2B <- 1000L; B_NC1 <- 1000L
if (DEBUG) {
  B_FL <- as.integer(getarg("--b-fl", "200")); B_NC2 <- as.integer(getarg("--b-small", "50"))
  B_NC2B <- B_NC2; B_NC1 <- B_NC2
}
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
PRESPEC_SHA <- "569d25c4fb4e6f0b89a90675fc3bd1c4e26d67d00baccba00105d8d184c63c35"
RESULT_FILE <- file.path(OUT, if (DEBUG) "RESULT_H62_DEBUG_NULL.md" else "RESULT_H62.md")

LOGF <- file.path(OUT, "run.log")
lg <- function(...) { s <- paste0(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = LOGF, append = TRUE) }
DIAG <- file.path(OUT, "diagnostics.tsv")
cat("item\tkey\tvalue\n", file = DIAG)
dg <- function(item, key, value) {
  value <- if (is.numeric(value)) vapply(value, function(v) format(v, digits = 10), "") else as.character(value)
  for (i in seq_along(value)) {
    cat(item, "\t", key[min(i, length(key))], "\t", value[i], "\n", sep = "", file = DIAG, append = TRUE)
  }
  lg("[diag ", item, "] ", paste(key, collapse = ","), " = ", paste(value, collapse = ", "))
}
lg("H62 run start ", format(Sys.time(), tz = "UTC", usetz = TRUE), " mode=", if (DEBUG) "DEBUG-NULL" else "FINAL",
   " B_FL=", B_FL, " B_NC2=", B_NC2, " B_NC2B=", B_NC2B, " B_NC1=", B_NC1)

STATE <- list(stage = "start")
write_stop <- function(rule, msg) {
  verdict <- if (rule == "S3") "stopped_by_rule S3 (PC-A failed): no number from the run is used" else
    paste0("untestable (", rule, "): ", msg)
  txt <- c("# RESULT H62", "",
           if (DEBUG) "DEBUG-NULL RUN: labels permuted; not a result." else NULL,
           paste0("Verdict: ", verdict), "",
           paste0("Stop rule ", rule, ": ", msg), "",
           "No score from this run enters the text.",
           paste0("Diagnostics written so far: ", DIAG))
  writeLines(txt, RESULT_FILE)
  if (!DEBUG) writeLines(txt, file.path(RUN, "RESULT_H62.md"))
  lg("STOP ", rule, ": ", msg)
  writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
  quit(save = "no", status = 0)
}
sha <- function(p) sub(" .*", "", system2("sha256sum", shQuote(p), stdout = TRUE))

## ---------------------------------------------------------------- prespec hash
for (p in c(file.path(RUN, "PRESPEC_H62.md"),
            file.path(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H62.md"))) {
  h <- sha(p); lg("prespec ", p, " sha256 ", h)
  if (h != PRESPEC_SHA) write_stop("S1", paste0("prespec sha256 mismatch at ", p, ": ", h))
}

## ---------------------------------------------------------------- diagnostic 1: input sha256
RF <- "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
R4 <- "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z"
inputs <- data.frame(
  role = c("GSE66676 series matrix", "GPL6244 annotation", "GSE185051 gene counts", "GSE185051 series matrix",
           "Adopted refit DEG table", "Refit unit table", "774-unit raw counts", "GENCODE v49 gene metadata",
           "Gate script", "Gate adult d", "Gate coverage", "Count script", "Count output", "Coverage script",
           "Coverage output", "Power script", "Power output", "Audit check script", "Audit check output",
           "Audit R checks", "Audit R output"),
  path = c("data/public_r4/GSE66676/GSE66676_series_matrix.txt.gz", "data/public_r4/GSE66676/GPL6244.annot.gz",
           "data/public_r4/GSE185051/GSE185051_Gene_count.csv.gz", "data/public_r4/GSE185051/GSE185051_series_matrix.txt.gz",
           file.path(RF, "deg_results.csv"), file.path(RF, "meta_matched.rds"),
           "figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds",
           "data/gencode_v49_gene_metadata.tsv.gz",
           file.path(R4, "gate/NG01/count_effect.R"), file.path(R4, "gate/NG01/adult_signature_effect_d.tsv"),
           file.path(R4, "gate/NG01/signature_coverage.tsv"), file.path(R4, "write/H62/count_n.py"),
           file.path(R4, "write/H62/count_n.out"), file.path(R4, "write/H62/coverage.py"),
           file.path(R4, "write/H62/coverage.out"), file.path(R4, "write/H62/mde.R"), file.path(R4, "write/H62/mde.out"),
           file.path(R4, "audit/H62/audit_checks.py"), file.path(R4, "audit/H62/audit_checks.out"),
           file.path(R4, "audit/H62/audit_checks.R"), file.path(R4, "audit/H62/audit_checks_R.out")),
  bytes = c(9304046, 7193437, 1963606, 4868, 5049936, 9082, 52980522, 1068672, rep(NA, 13)),
  sha256 = c("a4e255068ee8efdfd70b7538cd51fb27abbce29cc2e960a31a6a4cf7f33845d8",
             "99ee70f87e02c75ff1c8dc57b6c40c0482c7ef2760e004301315ac34344ac6dc",
             "2e9a97e17e4c41cb98d5cdd3405a97c47ca4644464b2398b655a06cd029e14f4",
             "34e5428b19d1ebc76865defbe73aeff4f15ac9b3d5c54f33f155cc0bdc2bb294",
             "0491ed10024388c3497dccb9f99bf669683edd5903c9d21c5308850d7c10bde6",
             "429e27376ad421f5f5fae19c75f5e760030e727fab65c5768ef562a9939edfc3",
             "f931e340546d679dc69dd0d6f712f4e73b1d805e80c741b22d3ed25498c68675",
             "49dcaae77323ffb20f5fd8399bb2899cd553e3634755a37b60a133d62549c189",
             "6e389f27236a9a7c0634cc22cc197d23a1f61db222b2b0bcf47f92eae3a410aa",
             "3e12988e7de57106b2d8cf11366587089d3f2d76c566377c7b94a35910ea5ec1",
             "c60157391d3a2dbf1daa6ebf0c25185d6d8e90cf5a0c1f67d4d0b1b6c6917fc5",
             "86fde1ea389bcae74f09aa175af20fc314fcce927f68a1c852697a244a4effd5",
             "3c148d7cb6155edec79a18982c8f73073fa797f3a04c14954d305f7047333ea3",
             "5c28ec0e20d180015adec9784fb12e45d7b9f2408fbd9c23162a386641832850",
             "69c1a123f1d51a2b79e9792906537f1f4770b535f909eace53c150eb3260cd1d",
             "f00c378bcd9b5a02cae553af332e7149ce89f278ac227453d990813d0d7a107e",
             "adbb8da8f4bd4174cd58ecf5edd682212d38d0ea45559411af6c94a20fb5d73f",
             "9b7bcb3a3d084306c972c4a5763dd825a93bd0ee4bf2af0327f18226a51213c4",
             "2916e5e5e2e1317f1191871c941c206d466e6d17213c2c42b8504ed0e013b91f",
             "064d73f6d9d8c76b797c4b3f358e2d49dd4ed817124e1c057bbed734ff572831",
             "3032793979e89ee41b3ddf0fba424ee366fa56a54c26f3dd2ebb06f4116d7fb6"))
P <- function(rel) file.path(ROOT, rel)
inputs$observed_sha256 <- vapply(inputs$path, function(p) if (file.exists(P(p))) sha(P(p)) else "MISSING", "")
inputs$observed_bytes <- vapply(inputs$path, function(p) if (file.exists(P(p))) file.size(P(p)) else NA_real_, 0)
inputs$match <- inputs$observed_sha256 == inputs$sha256 & (is.na(inputs$bytes) | inputs$bytes == inputs$observed_bytes)
write.table(inputs, file.path(OUT, "input_sha256.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (i in seq_len(nrow(inputs))) dg("1_sha256", inputs$path[i], paste(inputs$observed_sha256[i], if (inputs$match[i]) "match" else "MISMATCH"))
if (!all(inputs$match)) write_stop("S1", paste("sha256 or size mismatch:", paste(inputs$path[!inputs$match], collapse = "; ")))
lg("S1 passed: ", nrow(inputs), " of ", nrow(inputs), " input sha256 values match")

## ---------------------------------------------------------------- GSE66676 metadata (header rows)
read_series <- function(path) {
  L <- readLines(gzfile(path))
  b <- grep("^!series_matrix_table_begin", L); e <- grep("^!series_matrix_table_end", L)
  H <- L[seq_len(b - 1)]
  spl <- function(l) gsub('^"|"$', "", strsplit(l, "\t", fixed = TRUE)[[1]])
  meta <- list()
  for (l in H[startsWith(H, "!Sample_")]) {
    p <- spl(l); key <- p[1]; v <- p[-1]
    if (key == "!Sample_characteristics_ch1") {
      k <- sub(":.*", "", v[1]); meta[[paste0("ch:", k)]] <- ifelse(grepl(": ", v), sub("^[^:]*: ", "", v), v)
    } else if (is.null(meta[[key]])) meta[[key]] <- v else meta[[paste0(key, "#2")]] <- v
  }
  list(meta = meta, table_lines = if (length(b)) L[(b + 1):(e - 1)] else character())
}
g6 <- read_series(P(inputs$path[1]))
M6 <- g6$meta
hdr6 <- gsub('"', "", strsplit(g6$table_lines[1], "\t", fixed = TRUE)[[1]])
gsm <- M6[["!Sample_geo_accession"]]
if (!identical(hdr6[-1], gsm)) write_stop("S2", "GSE66676 matrix column order differs from !Sample_geo_accession")
hist_all <- M6[["ch:histology"]]; sex_all <- M6[["ch:gender"]]
age_all <- as.numeric(M6[["ch:age (yrs)"]]); bmi_all <- as.numeric(M6[["ch:bmi"]])
pid_all <- M6[["!Sample_description#2"]]
n0 <- sum(hist_all == "No NAFLD"); n1 <- sum(hist_all != "No NAFLD")
## diagnostic 2: n by histology, sex, age, BMI (must equal count_n.out)
fmt_grp <- function(g) {
  ix <- which(hist_all == g)
  sprintf("%s n %d F/M %d / %d age med %.2f range %.2f-%.2f bmi med %.2f range %.2f-%.2f", g, length(ix),
          sum(sex_all[ix] == "female"), sum(sex_all[ix] == "male"), median(age_all[ix]), min(age_all[ix]),
          max(age_all[ix]), median(bmi_all[ix]), min(bmi_all[ix]), max(bmi_all[ix]))
}
cn_out <- readLines(P(inputs$path[13]))
grp_lines <- vapply(sort(unique(hist_all)), fmt_grp, "")
for (s in grp_lines) dg("2_n_by_group", "GSE66676", s)
dg("2_n_by_group", c("No NAFLD", "NAFLD", "total", "female", "male", "missing_age", "missing_bmi", "distinct_participant_codes"),
   c(n0, n1, length(hist_all), sum(sex_all == "female"), sum(sex_all == "male"), sum(is.na(age_all)), sum(is.na(bmi_all)),
     length(unique(pid_all))))
if (n0 != 34 || n1 != 33 || length(hist_all) != 67) write_stop("S2", sprintf("n is %d / %d (total %d), not 34 / 33 (67)", n0, n1, length(hist_all)))
if (!all(grp_lines %in% cn_out)) write_stop("S2", paste("group summaries differ from count_n.out:", paste(grp_lines[!grp_lines %in% cn_out], collapse = " | ")))
if (any(is.na(age_all)) || any(is.na(bmi_all)) || length(unique(pid_all)) != 67) write_stop("S2", "missing covariate or duplicated participant code")
lg("S2 passed: 34 No NAFLD / 33 NAFLD, summaries equal count_n.out")

## ---------------------------------------------------------------- frozen signature
deg <- read.csv(P(inputs$path[5]))
sig <- deg[which(deg$padj < 0.05 & abs(deg$logFC) > 0.5), ]
sig$dir <- ifelse(sig$logFC > 0, "up", "down")
deg$id <- sub("\\..*", "", deg$gene); sig$id <- sub("\\..*", "", sig$gene)
up_id <- sig$id[sig$dir == "up"]; dn_id <- sig$id[sig$dir == "down"]
named <- !startsWith(sig$symbol, "ENSG")
lg("signature: up ", sum(sig$dir == "up"), " down ", sum(sig$dir == "down"), "; Ensembl-only symbols ", sum(!named),
   " (up ", sum(!named & sig$dir == "up"), ", down ", sum(!named & sig$dir == "down"), ")")
if (sum(sig$dir == "up") != 1333 || sum(sig$dir == "down") != 521 || sum(!named) != 364 || any(duplicated(sig$symbol)))
  write_stop("S5", "signature does not reproduce 1,333 up / 521 down / 364 Ensembl-only")
sym_up <- sig$symbol[sig$dir == "up" & named]; sym_dn <- sig$symbol[sig$dir == "down" & named]
sig_named <- sig$symbol[named]

## ---------------------------------------------------------------- mapping to GPL6244 matrix probes (before any value is read)
ann <- read.delim(gzfile(P(inputs$path[2])), skip = 27, quote = "", comment.char = "", colClasses = "character",
                  check.names = FALSE, na.strings = c("", "NA"))
ann <- ann[grepl("^[0-9]+$", ann$ID), c("ID", "Gene symbol")]
names(ann) <- c("ID", "symbol")
single <- ann[!is.na(ann$symbol) & !grepl("///", ann$symbol, fixed = TRUE), ]
mat_ids <- gsub('"', "", sub("\t.*", "", g6$table_lines[-1]))
if (length(mat_ids) != 28869 || anyDuplicated(mat_ids)) write_stop("S5", paste("matrix probe rows:", length(mat_ids)))
single_m <- single[single$ID %in% mat_ids, ]
probes_of <- split(single_m$ID, single_m$symbol)          # symbol -> matrix probes (single-symbol rows only)
dg("4_mapping", c("annot_single_symbol_rows", "annot_unique_single_symbols", "single_symbol_rows_in_matrix"),
   c(nrow(single), length(unique(single$symbol)), nrow(single_m)))
map_up <- probes_of[intersect(sym_up, names(probes_of))]; map_up <- map_up[order(names(map_up))]
map_dn <- probes_of[intersect(sym_dn, names(probes_of))]; map_dn <- map_dn[order(names(map_dn))]
n_up_sym <- length(map_up); n_dn_sym <- length(map_dn); n_up_pr <- sum(lengths(map_up)); n_dn_pr <- sum(lengths(map_dn))
mapping <- rbind(data.frame(direction = "up", symbol = rep(names(map_up), lengths(map_up)), probe = unlist(map_up, use.names = FALSE)),
                 data.frame(direction = "down", symbol = rep(names(map_dn), lengths(map_dn)), probe = unlist(map_dn, use.names = FALSE)))
mapping$probes_for_symbol <- ave(mapping$probe, mapping$symbol, FUN = length)
write.table(mapping, file.path(OUT, "signature_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
dg("4_mapping", c("mapped_up_symbols", "mapped_down_symbols", "mapped_up_probes", "mapped_down_probes"),
   c(n_up_sym, n_dn_sym, n_up_pr, n_dn_pr))
pps <- table(lengths(c(map_up, map_dn)))
dg("4_mapping", paste0("symbols_with_", names(pps), "_probes"), as.integer(pps))
multi <- lengths(c(map_up, map_dn)); multi <- multi[multi > 2]
dg("4_mapping", paste0("more_than_2_probes:", names(multi)), as.integer(multi))
pool <- setdiff(names(probes_of), sig_named)
dg("4_mapping", "NC1_pool_symbols", length(pool))
if (n_up_sym < 500 || n_dn_sym < 100) write_stop("S5", sprintf("coverage %d up / %d down below 500 / 100", n_up_sym, n_dn_sym))
if (n_up_sym != 881 || n_dn_sym != 183 || n_up_pr != 920 || n_dn_pr != 186)
  write_stop("S5", sprintf("recount %d / %d symbols, %d / %d probes differs from 881 / 183, 920 / 186", n_up_sym, n_dn_sym, n_up_pr, n_dn_pr))
if (length(pool) != 17928) write_stop("S5", paste("NC-1 pool is", length(pool), "not 17,928 (code or input mismatch)"))
lg("S5 passed: 881 / 183 symbols, 920 / 186 probes; NC-1 pool 17,928; mapping written before reading values")

## ---------------------------------------------------------------- shared score functions
zmat <- function(x, ref) {                 # same arithmetic as gate count_effect.R
  z <- (x - rowMeans(x[, ref, drop = FALSE])) / apply(x[, ref, drop = FALSE], 1, sd)
  z[!is.finite(z)] <- NA; z
}
sym_avg <- function(z, sets) {             # per-symbol mean of finite z over the symbol's features
  idx <- unlist(sets, use.names = FALSE); grp <- rep(seq_along(sets), lengths(sets))
  zz <- z[idx, , drop = FALSE]; ok <- !is.na(zz); zz[!ok] <- 0
  S <- rowsum(zz, grp, reorder = TRUE); C <- rowsum(ok * 1, grp, reorder = TRUE)
  A <- S / C; A[C == 0] <- NA; rownames(A) <- names(sets); A
}
score_ud <- function(A_up, A_dn) colMeans(A_up, na.rm = TRUE) - colMeans(A_dn, na.rm = TRUE)

## ---------------------------------------------------------------- PC-A: reproduce gate adult d (stop rule S3)
STATE$stage <- "PC-A"
m <- readRDS(P(inputs$path[6]))
cnt <- readRDS(P(inputs$path[7])); cnt <- as.matrix(cnt); rownames(cnt) <- sub("\\..*", "", rownames(cnt))
cnt <- cnt[, m$sample_id]
pca <- list()
for (ds in unique(m$dataset)) {
  mm <- m[m$dataset == ds, ]; y <- DGEList(cnt[, mm$sample_id]); keep <- filterByExpr(y, group = mm$group_binary)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE]); lc <- cpm(y, log = TRUE, prior.count = 2)
  ctl <- mm$group_binary == "Control"
  z <- zmat(lc, ctl)
  rows_up <- which(rownames(z) %in% up_id); rows_dn <- which(rownames(z) %in% dn_id)
  sc <- score_ud(sym_avg(z, as.list(rows_up)), sym_avg(z, as.list(rows_dn)))
  dfun <- function(a, b) (mean(a) - mean(b)) / sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  for (g in c("NAFL", "Borderline", "NASH", "AllDisease")) {
    sel <- if (g == "AllDisease") mm$group_binary == "Disease" else (mm$diagnosis_harmonized %in% g)
    if (sum(sel) < 3) next
    pca[[length(pca) + 1]] <- data.frame(dataset = ds, group = g, n_dis = sum(sel), n_ctl = sum(ctl), d = dfun(sc[sel], sc[ctl]))
  }
  sel <- mm$group_binary == "Disease" & !is.na(mm$fibrosis_stage) & mm$fibrosis_stage == 0
  if (sum(sel) >= 3) pca[[length(pca) + 1]] <- data.frame(dataset = ds, group = "F0_disease", n_dis = sum(sel), n_ctl = sum(ctl), d = dfun(sc[sel], sc[ctl]))
}
pca <- do.call(rbind, pca)
gate_d <- read.delim(P(inputs$path[10]))
cmp <- merge(gate_d, pca, by = c("dataset", "group"), suffixes = c("_gate", "_run"), all = TRUE, sort = FALSE)
cmp$abs_diff <- abs(cmp$d_gate - cmp$d_run)
write.table(cmp, file.path(OUT, "pcA_reproduction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (i in seq_len(nrow(cmp))) dg("3_PC-A", paste(cmp$dataset[i], cmp$group[i]), sprintf("gate %.9f run %.9f absdiff %.3g", cmp$d_gate[i], cmp$d_run[i], cmp$abs_diff[i]))
pcA_max <- max(cmp$abs_diff); dg("3_PC-A", c("rows", "max_abs_diff"), c(nrow(cmp), pcA_max))
pcA_pass <- nrow(cmp) == 18 && !any(is.na(cmp$abs_diff)) && all(cmp$n_dis_gate == cmp$n_dis_run) && pcA_max <= 1e-6
rm(cnt); invisible(gc())
if (!pcA_pass) write_stop("S3", sprintf("PC-A failed: %d rows, max |diff| %.3g", nrow(cmp), pcA_max))
lg("S3 passed: PC-A 18 of 18 rows reproduced, max |diff| ", signif(pcA_max, 3))

## ---------------------------------------------------------------- GSE66676 expression
STATE$stage <- "expression"
X <- read.delim(text = g6$table_lines, quote = "\"", check.names = FALSE,
                colClasses = c("character", rep("numeric", length(gsm))))
rownames(X) <- X[[1]]; X <- as.matrix(X[, -1]); stopifnot(identical(colnames(X), gsm), identical(rownames(X), mat_ids))
dg("expr", c("probes", "samples", "NA_values"), c(nrow(X), ncol(X), sum(is.na(X))))

## ---------------------------------------------------------------- PC-B: sex from Y probes (stop rule S4)
yprobes <- c(RPS4Y1 = "8176375", DDX3Y = "8176624", EIF1AY = "8176719", UTY = "8177137", KDM5D = "8177232")
ymean <- colMeans(X[yprobes, , drop = FALSE])
set.seed(SEED)
km <- kmeans(ymean, centers = 2, nstart = 25)
hi <- which.max(km$centers[, 1])
sex_call <- ifelse(km$cluster == hi, "male", "female")
disc <- which(sex_call != sex_all)
pcB <- data.frame(gsm = gsm, participant = pid_all, gender = sex_all, y_mean = ymean, cluster_call = sex_call, discordant = sex_call != sex_all)
write.table(pcB, file.path(OUT, "pcB_sex_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
dg("5_PC-B", c("cluster_centers_low_high", "n_called_male", "n_discordant", "discordant_ids"),
   c(paste(signif(sort(km$centers[, 1]), 5), collapse = ";"), sum(sex_call == "male"), length(disc),
     if (length(disc)) paste(gsm[disc], collapse = ";") else "none"))
dg("5_PC-B", "prediction_0_discordant", if (length(disc) == 0) "met" else "not met")
if (length(disc) > 3) write_stop("S4", sprintf("%d discordant sex calls (> 3): labels unreliable", length(disc)))
keep_s <- setdiff(seq_along(gsm), disc)
if (length(disc)) lg("S4: removed ", length(disc), " discordant participants: ", paste(gsm[disc], collapse = ","))
lg("S4 passed: ", length(disc), " discordant")

## ---------------------------------------------------------------- analysis set
X <- X[, keep_s, drop = FALSE]
hist <- hist_all[keep_s]; sex <- sex_all[keep_s]; age <- age_all[keep_s]; bmi <- bmi_all[keep_s]; ids <- gsm[keep_s]
if (DEBUG) {
  set.seed(1)   # DEBUG ONLY: null labels so outcome code never sees real labels while being debugged
  hist <- hist[sample.int(length(hist))]
  lg("DEBUG: histology labels permuted (seed 1); every number below is null data")
}
nafld <- as.integer(hist != "No NAFLD"); ref <- nafld == 0L
dat <- data.frame(nafld = nafld, sex = factor(sex, levels = c("female", "male")), age = age, bmi = bmi)
N <- nrow(dat); N1 <- sum(nafld); N0 <- sum(1 - nafld)
Xd <- model.matrix(~ nafld + sex + age + bmi, dat); jN <- which(colnames(Xd) == "nafld")
Zd <- model.matrix(~ sex + age + bmi, dat)
DF <- N - ncol(Xd)
make_tfun <- function(Xm, j) {
  q <- qr(Xm); if (q$rank < ncol(Xm)) stop("rank-deficient design")
  V <- chol2inv(qr.R(q)); dfr <- nrow(Xm) - ncol(Xm)
  function(Y) {                                # Y: n x k matrix; returns t of coefficient j per column
    Y <- as.matrix(Y); b <- qr.coef(q, Y); r <- qr.resid(q, Y)
    s2 <- colSums(r^2) / dfr; b[j, ] / sqrt(s2 * V[j, j])
  }
}
tfun <- make_tfun(Xd, jN)
pt2 <- function(t, df) 2 * pt(-abs(t), df)

## ---------------------------------------------------------------- scores
sig_sets_up <- lapply(map_up, function(p) match(p, rownames(X)))
sig_sets_dn <- lapply(map_dn, function(p) match(p, rownames(X)))
sig_rows <- sort(unique(c(unlist(sig_sets_up), unlist(sig_sets_dn))))
build_score <- function(Xs, refv, set_up, set_dn) {
  z <- zmat(Xs, refv); list(A_up = sym_avg(z, set_up), A_dn = sym_avg(z, set_dn), z = z)
}
S0 <- build_score(X, ref, sig_sets_up, sig_sets_dn)
score <- score_ud(S0$A_up, S0$A_dn)

## ---------------------------------------------------------------- diagnostics 6-8
sd_ref_all <- apply(X[, ref, drop = FALSE], 1, sd)
sdm <- sd_ref_all[sig_rows]
dg("6_refSD_mapped_probes", c("min", "p01", "median"), c(min(sdm), quantile(sdm, 0.01, names = FALSE), median(sdm)))
dg("6_refSD_mapped_probes", "probes_with_absz_gt10_any_participant", sum(apply(abs(S0$z[sig_rows, , drop = FALSE]) > 10, 1, any, na.rm = TRUE)))
dg("6_refSD_mapped_probes", "nonfinite_z_values_mapped_probes", sum(is.na(S0$z[sig_rows, ])))
C <- cor(X); diag(C) <- NA
medcor <- apply(C, 2, median, na.rm = TRUE)
q1 <- quantile(medcor, 0.25, names = FALSE); iqr <- IQR(medcor); thr <- q1 - 3 * iqr
dg("7_median_intersample_cor", c("min", "Q1", "median", "IQR", "threshold_Q1_minus_3IQR"), c(min(medcor), q1, median(medcor), iqr, thr))
dg("7_median_intersample_cor", "below_threshold", if (any(medcor < thr)) paste(names(medcor)[medcor < thr], collapse = ";") else "none")
write.table(data.frame(gsm = ids, median_cor = medcor), file.path(OUT, "diag7_median_cor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (g in 0:1) {
  ix <- nafld == g
  dg("8_covariate_balance", paste0(if (g) "NAFLD" else "No NAFLD", c(":n", ":female", ":male", ":age_median", ":age_mean", ":bmi_median", ":bmi_mean")),
     c(sum(ix), sum(sex[ix] == "female"), sum(sex[ix] == "male"), median(age[ix]), mean(age[ix]), median(bmi[ix]), mean(bmi[ix])))
}
dg("8_covariate_balance", "pearson_cor_bmi_nafld", cor(bmi, nafld))

## ---------------------------------------------------------------- diagnostic 9: NC-2 and NC-2b calibration
STATE$stage <- "NC-2"
fit_red <- qr.fitted(qr(Zd), score); res_red <- score - fit_red
set.seed(SEED)
P2 <- vapply(seq_len(B_NC2), function(b) sample.int(N), integer(N))
t_nc2 <- tfun(fit_red + matrix(res_red[P2], N))
nc2_frac <- mean(pt2(t_nc2, DF) < 0.05)
dg("9_NC-2", c("B", "fraction_param_p_lt_0.05", "prediction_0.035_to_0.065", "calibrated_0.02_to_0.08"),
   c(B_NC2, nc2_frac, if (nc2_frac >= 0.035 && nc2_frac <= 0.065) "met" else "not met", nc2_frac >= 0.02 && nc2_frac <= 0.08))
STATE$stage <- "NC-2b"
Xs_sig <- X[sig_rows, , drop = FALSE]
remap <- function(sets) lapply(sets, function(ix) match(ix, sig_rows))
su_r <- remap(sig_sets_up); sd_r <- remap(sig_sets_dn)
set.seed(SEED)
t_nc2b <- numeric(B_NC2B)
for (b in seq_len(B_NC2B)) {
  lab <- sample(nafld)
  Sb <- build_score(Xs_sig, lab == 0L, su_r, sd_r)
  Xb <- Xd; Xb[, jN] <- lab
  t_nc2b[b] <- make_tfun(Xb, jN)(score_ud(Sb$A_up, Sb$A_dn))
}
nc2b_frac <- mean(pt2(t_nc2b, DF) < 0.05)
dg("9_NC-2b", c("B", "fraction_param_p_lt_0.05", "prediction_0.035_to_0.065", "calibrated_0.02_to_0.08"),
   c(B_NC2B, nc2b_frac, if (nc2b_frac >= 0.035 && nc2b_frac <= 0.065) "met" else "not met", nc2b_frac >= 0.02 && nc2b_frac <= 0.08))

## ---------------------------------------------------------------- diagnostic 10: NC-1 random gene sets
STATE$stage <- "NC-1"
pool_sets <- lapply(probes_of[pool], function(p) match(p, rownames(X)))
sym_sd <- function(sets) vapply(sets, function(ix) mean(sd_ref_all[ix]), 0)
sd_pool <- sym_sd(pool_sets); sd_up <- sym_sd(sig_sets_up); sd_dn <- sym_sd(sig_sets_dn)
brk <- quantile(c(sd_pool, sd_up, sd_dn), probs = seq(0, 1, 0.1), names = FALSE)
if (anyDuplicated(brk)) stop("duplicated decile breakpoints")
dec <- function(v) cut(v, breaks = brk, include.lowest = TRUE, labels = FALSE)
dpool <- dec(sd_pool); dup_ <- dec(sd_up); ddn <- dec(sd_dn)
k_up <- tabulate(dup_, 10); k_dn <- tabulate(ddn, 10); n_pool_dec <- tabulate(dpool, 10)
dg("10_NC-1", paste0("decile", 1:10, ":pool/up/down"), paste(n_pool_dec, k_up, k_dn, sep = "/"))
if (any(k_up + k_dn > n_pool_dec)) stop("a decile has fewer pool symbols than the signature needs")
zpool <- zmat(X[unique(unlist(pool_sets)), , drop = FALSE], ref)
A_pool <- sym_avg(zpool, lapply(pool_sets, function(ix) match(rownames(X)[ix], rownames(zpool))))
pool_by_dec <- split(seq_along(pool), dpool)
set.seed(SEED)
nc1_t <- numeric(B_NC1)
for (b in seq_len(B_NC1)) {
  uu <- integer(0); dd <- integer(0)
  for (d in 1:10) {
    cand <- pool_by_dec[[as.character(d)]]; k <- k_up[d] + k_dn[d]
    if (k == 0) next
    pick <- cand[sample.int(length(cand), k)]
    uu <- c(uu, pick[seq_len(k_up[d])]); dd <- c(dd, pick[k_up[d] + seq_len(k_dn[d])])
  }
  nc1_t[b] <- tfun(score_ud(A_pool[uu, , drop = FALSE], A_pool[dd, , drop = FALSE]))
}
nc1_q <- quantile(nc1_t, c(0.05, 0.5, 0.95), names = FALSE)
dg("10_NC-1", c("B", "t_p05", "t_median", "t_p95", "prediction_median_in_-1.5_to_1.5"),
   c(B_NC1, nc1_q, if (abs(nc1_q[2]) <= 1.5) "met" else "not met"))

## ---------------------------------------------------------------- diagnostic 11: GSE185051 mapping
STATE$stage <- "GSE185051"
g18 <- read_series(P(inputs$path[4]))$meta
cnt18 <- read.csv(gzfile(P(inputs$path[3])), check.names = FALSE)
titles <- g18[["!Sample_title"]]
col_of <- ifelse(startsWith(titles, "LB-"), paste0("L", sub("^LB-", "", titles)), gsub("-", ".", titles, fixed = TRUE))
is_child <- g18[["ch:disease"]] == "NAFLD disease"
mapped <- col_of %in% colnames(cnt18)
dg("11_GSE185051", c("children", "children_mapped_to_columns", "controls", "controls_mapped_to_columns"),
   c(sum(is_child), sum(is_child & mapped), sum(!is_child), sum(!is_child & mapped)))
dupsym <- unique(cnt18$Gene[duplicated(cnt18$Gene)])
dg("11_GSE185051", c("duplicated_symbols_dropped", "rows_dropped"), c(paste(dupsym, collapse = ";"), sum(cnt18$Gene %in% dupsym)))
cnt18 <- cnt18[!cnt18$Gene %in% dupsym, ]
CM <- as.matrix(cnt18[, -1]); rownames(CM) <- cnt18$Gene
ch_ix <- which(is_child & mapped)
nas <- as.numeric(g18[["ch:nas"]][ch_ix]); nash_status <- g18[["ch:nash_status"]][ch_ix]
if (DEBUG) { set.seed(1); perm <- sample.int(length(ch_ix)); nas <- nas[perm]; nash_status <- nash_status[perm]
  lg("DEBUG: GSE185051 NAS and nash_status permuted (seed 1)") }
y18 <- DGEList(CM[, col_of[ch_ix]])
keep18 <- filterByExpr(y18)
y18 <- calcNormFactors(y18[keep18, , keep.lib.sizes = FALSE])
lc18 <- cpm(y18, log = TRUE, prior.count = 2)
z18 <- zmat(lc18, rep(TRUE, ncol(lc18)))
up18 <- intersect(sym_up, rownames(z18)); dn18 <- intersect(sym_dn, rownames(z18))
score18 <- score_ud(z18[up18, , drop = FALSE], z18[dn18, , drop = FALSE])
ls18 <- colSums(CM[, col_of[ch_ix]])
dg("11_GSE185051", c("mapped_up_before_filter", "mapped_down_before_filter", "genes_kept_filterByExpr", "mapped_up_after_filter",
                     "mapped_down_after_filter", "child_libsize_min", "child_libsize_median", "child_libsize_max"),
   c(length(intersect(sym_up, rownames(CM))), length(intersect(sym_dn, rownames(CM))), sum(keep18), length(up18), length(dn18),
     min(ls18), median(ls18), max(ls18)))
# adult controls: display only (genes kept in children; TMM over 57; z against the 52 children in that object)
all_ix <- which(mapped)
y57 <- DGEList(CM[keep18, col_of[all_ix]]); y57 <- calcNormFactors(y57)
lc57 <- cpm(y57, log = TRUE, prior.count = 2)
z57 <- zmat(lc57, is_child[all_ix])
score57 <- score_ud(z57[up18, , drop = FALSE], z57[dn18, , drop = FALSE])
ctrl_tab <- data.frame(column = col_of[all_ix][!is_child[all_ix]], score_display_only = score57[!is_child[all_ix]],
                       sex = g18[["ch:gender"]][all_ix][!is_child[all_ix]], age = g18[["ch:age"]][all_ix][!is_child[all_ix]])
write.table(ctrl_tab, file.path(OUT, "D6_adult_controls_display_only.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

lg("Diagnostics 1-11 written to ", DIAG, " before T1")

## ================================================================ T1 (primary)
STATE$stage <- "T1"
fit <- lm(score ~ nafld + sex + age + bmi, data = dat)
co <- summary(fit)$coefficients
b_T1 <- co["nafld", "Estimate"]; se_T1 <- co["nafld", "Std. Error"]; t_obs <- co["nafld", "t value"]; p_param <- co["nafld", "Pr(>|t|)"]
ci_T1 <- confint(fit)["nafld", ]
if (abs(tfun(score) - t_obs) > 1e-8) stop("fast t does not reproduce lm t")
sp <- sqrt(((N1 - 1) * var(score[nafld == 1]) + (N0 - 1) * var(score[nafld == 0])) / (N - 2))
d_adj <- b_T1 / sp; d_adj_ci <- ci_T1 / sp
# Freedman-Lane permutation, 10,000
set.seed(SEED)
PF <- vapply(seq_len(B_FL), function(b) sample.int(N), integer(N))
t_fl <- tfun(fit_red + matrix(res_red[PF], N))
p_fl <- (1 + sum(abs(t_fl) >= abs(t_obs))) / (B_FL + 1)
p_nc2b <- (1 + sum(abs(t_nc2b) >= abs(t_obs))) / (B_NC2B + 1)
# governing p (rules in T1, NC-2, NC-2b)
gov_rule <- "parametric p (parametric and Freedman-Lane p on the same side of 0.05; NC-2 and NC-2b calibrated)"
p_gov <- p_param
if ((p_param < 0.05) != (p_fl < 0.05)) { p_gov <- p_fl; gov_rule <- "Freedman-Lane p (parametric and permutation p on different sides of 0.05)" }
if (nc2_frac < 0.02 || nc2_frac > 0.08) { p_gov <- p_fl; gov_rule <- "Freedman-Lane p (NC-2 outside 0.02 to 0.08)" }
if (nc2b_frac < 0.02 || nc2b_frac > 0.08) { p_gov <- max(p_fl, p_nc2b); gov_rule <- "larger of Freedman-Lane and NC-2b permutation p (NC-2b outside 0.02 to 0.08)" }
nc1_rank <- mean(nc1_t >= t_obs); nc1_above95 <- t_obs > nc1_q[3]

## ---------------------------------------------------------------- PC-C (limma, T1 design, all probes)
STATE$stage <- "PC-C"
lf <- eBayes(lmFit(X, Xd))
tt <- topTable(lf, coef = "nafld", number = Inf, sort.by = "none")
pcc <- tt[c("8048205", "7992828"), c("logFC", "t", "P.Value")]; pcc$gene <- c("IGFBP2", "IL32")
pcc_dir_ok <- c(pcc$logFC[1] < 0, pcc$logFC[2] > 0)
pcc_met <- all(pcc_dir_ok) && any(pcc$P.Value < 0.05)
pcc_fail <- !any(pcc_dir_ok)
pcc_status <- if (pcc_met) "met" else if (pcc_fail) "FAILED (both genes opposite to the adult direction)" else "not met (partial; T1 reading unchanged)"

## ---------------------------------------------------------------- verdict (prespec Reading)
mde_fun <- function(na, nb, df) uniroot(function(d) { ncp <- d / sqrt(1 / na + 1 / nb); q <- qt(.975, df); 1 - pt(q, df, ncp) + pt(-q, df, ncp) - .8 }, c(.1, 2))$root
mde <- mde_fun(N0, N1, DF)
if (p_gov < 0.05 && b_T1 > 0) {
  verdict <- if (pcc_fail) "pass_pcc_failed" else "pass"
} else if (p_gov < 0.05 && b_T1 < 0) verdict <- "opposite" else verdict <- "null"
excl097 <- d_adj_ci[2] < 0.97
f2 <- function(x) formatC(x, format = "f", digits = 2)
fp <- function(p) if (p < 0.001) formatC(p, format = "e", digits = 1) else formatC(p, format = "g", digits = 2)
sentence <- switch(verdict,
  pass = paste0("The adult disease-state signature was also shifted in adolescent NAFLD: in ", N,
    " severely obese adolescents undergoing bariatric surgery (GSE66676), the frozen signature score was higher in ", N1,
    " with NAFLD than in ", N0, " without (adjusted difference ", f2(b_T1), ", d ", f2(d_adj), ", 95% CI ", f2(d_adj_ci[1]), " to ",
    f2(d_adj_ci[2]), ", P = ", fp(p_gov), "; sex, age and BMI adjusted).",
    if (nc1_above95) paste0(" The shift was larger than for 95% of random gene sets of the same size (competitive rank ", formatC(nc1_rank, format = "f", digits = 3), ").") else " Random gene sets of the same size shifted similarly.",
    " Most NAFLD participants had steatosis without NASH (", sum(hist == "NAFLD not NASH"), " of ", N1, ")."),
  pass_pcc_failed = paste0("In GSE66676 the frozen signature score was higher in ", N1, " adolescents with NAFLD than in ", N0,
    " without (adjusted difference ", f2(b_T1), ", d ", f2(d_adj), ", P = ", fp(p_gov), "), but IGFBP2 and IL32 did not move in the adult direction."),
  opposite = paste0("In GSE66676 the adult signature score was lower in ", N1, " adolescents with NAFLD than in ", N0,
    " without (d ", f2(d_adj), ", 95% CI ", f2(d_adj_ci[1]), " to ", f2(d_adj_ci[2]), ", P = ", fp(p_gov),
    "); controls were severely obese, the platform was an array and the score used symbol-mapped genes only."),
  null = paste0("In GSE66676, the adult signature score did not differ detectably between ", N1, " adolescents with NAFLD and ", N0,
    " without (d ", f2(d_adj), ", 95% CI ", f2(d_adj_ci[1]), " to ", f2(d_adj_ci[2]), "); the cohort had 80% power for d of ", f2(mde), ".",
    if (excl097) " The median adult NAFL-versus-control shift (d 0.97) is outside this interval." else ""))

## ---------------------------------------------------------------- descriptive outputs D1-D8 (nominal)
STATE$stage <- "descriptive"
lmT1 <- function(s, d = dat, form = s ~ nafld + sex + age + bmi) {
  d$s <- s; f <- lm(form, data = d); c(coef = coef(f)[["nafld"]], t = summary(f)$coefficients["nafld", "t value"],
                                       p = summary(f)$coefficients["nafld", "Pr(>|t|)"], n = nrow(d))
}
D <- list()
wt <- t.test(score[nafld == 1], score[nafld == 0])
cohen_d <- (mean(score[nafld == 1]) - mean(score[nafld == 0])) / sp
hedges_g <- cohen_d * (1 - 3 / (4 * N - 9))
D$D1 <- data.frame(output = "D1", item = c("welch_t", "welch_df", "welch_p_nominal", "mean_NAFLD", "mean_NoNAFLD", "cohen_d", "hedges_g", "pooled_sd"),
                   value = c(wt$statistic, wt$parameter, wt$p.value, mean(score[nafld == 1]), mean(score[nafld == 0]), cohen_d, hedges_g, sp))
s_up <- colMeans(S0$A_up, na.rm = TRUE); s_dn <- colMeans(S0$A_dn, na.rm = TRUE)
r_up <- lmT1(s_up); r_dn <- lmT1(s_dn)
D$D2 <- data.frame(output = "D2", item = paste0(rep(c("up_only_", "down_only_"), each = 4), names(r_up)), value = c(r_up, r_dn))
dev <- X - rowMeans(X[, ref, drop = FALSE])
s_unw <- score_ud(sym_avg(dev, sig_sets_up), sym_avg(dev, sig_sets_dn)); r_unw <- lmT1(s_unw)
D$D3 <- data.frame(output = "D3", item = paste0("unweighted_", names(r_unw)), value = r_unw)
hl <- c("No NAFLD", "NAFLD not NASH", "Borderline NASH", "Definite NASH")
d4 <- do.call(rbind, lapply(hl, function(g) { v <- score[hist == g]; data.frame(output = "D4", item = paste0(g, c(":n", ":median", ":Q1", ":Q3")),
  value = c(length(v), median(v), quantile(v, c(.25, .75), names = FALSE))) }))
sub_ix <- hist %in% c("No NAFLD", "NAFLD not NASH")
r_d4 <- lmT1(score[sub_ix], d = dat[sub_ix, ])
D$D4 <- rbind(d4, data.frame(output = "D4", item = paste0("NAFLDnotNASH_vs_NoNAFLD_", names(r_d4)), value = r_d4))
r_d5 <- lmT1(score, form = s ~ nafld + sex + age)
D$D5 <- data.frame(output = "D5", item = paste0("no_bmi_", names(r_d5)), value = r_d5)
sp6 <- suppressWarnings(cor.test(score18, nas, method = "spearman", exact = FALSE))
med_status <- tapply(score18, nash_status, median)
D$D6 <- data.frame(output = "D6", item = c("n_children", "spearman_rho_score_NAS", "spearman_p_nominal", paste0("median_score:", names(med_status)),
                                         paste0("adult_control_display_only:", ctrl_tab$column)),
                   value = c(length(score18), sp6$estimate, sp6$p.value, med_status, ctrl_tab$score_display_only))
lfc_sym_up <- vapply(map_up, function(p) mean(tt[p, "logFC"]), 0); lfc_sym_dn <- vapply(map_dn, function(p) mean(tt[p, "logFC"]), 0)
D$D7 <- data.frame(output = "D7", item = c("fraction_adult_sign_all", "fraction_adult_sign_up", "fraction_adult_sign_down", "n_symbols"),
                   value = c(mean(c(lfc_sym_up > 0, lfc_sym_dn < 0)), mean(lfc_sym_up > 0), mean(lfc_sym_dn < 0), length(lfc_sym_up) + length(lfc_sym_dn)))
best <- function(sets) lapply(sets, function(ix) ix[which.max(sd_ref_all[ix])])
s_d8 <- score_ud(sym_avg(S0$z, best(sig_sets_up)), sym_avg(S0$z, best(sig_sets_dn))); r_d8 <- lmT1(s_d8)
D$D8 <- data.frame(output = "D8", item = paste0("max_sd_probe_", names(r_d8)), value = r_d8)
Dall <- do.call(rbind, D)
write.table(Dall, file.path(OUT, "descriptive_D1_D8.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---------------------------------------------------------------- tables
t1 <- data.frame(item = c("n_NoNAFLD", "n_NAFLD", "coefficient", "se", "ci95_low", "ci95_high", "t", "df", "p_parametric", "p_freedman_lane",
                          "B_freedman_lane", "p_nc2b_regenerated", "p_governing", "governing_rule", "d_adj", "d_adj_ci_low", "d_adj_ci_high",
                          "hedges_g_unadjusted", "nc1_competitive_rank", "nc1_t_p95", "t_above_nc1_p95", "mde80_noncentral_t", "verdict"),
                 value = c(N0, N1, b_T1, se_T1, ci_T1, t_obs, DF, p_param, p_fl, B_FL, p_nc2b, p_gov, gov_rule, d_adj, d_adj_ci,
                           hedges_g, nc1_rank, nc1_q[3], nc1_above95, mde, verdict))
write.table(t1, file.path(OUT, "T1_result.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ctrls <- data.frame(control = c("PC-A", "PC-B", "PC-C IGFBP2", "PC-C IL32", "PC-C", "NC-1", "NC-2", "NC-2b"),
                    result = c(sprintf("18/18 rows, max |diff| %.2g", pcA_max), sprintf("%d discordant of %d", length(disc), length(gsm)),
                               sprintf("logFC %.3f, t %.2f, nominal p %.3g", pcc$logFC[1], pcc$t[1], pcc$P.Value[1]),
                               sprintf("logFC %.3f, t %.2f, nominal p %.3g", pcc$logFC[2], pcc$t[2], pcc$P.Value[2]), pcc_status,
                               sprintf("median t %.3f (p05 %.3f, p95 %.3f); competitive rank %.3f", nc1_q[2], nc1_q[1], nc1_q[3], nc1_rank),
                               sprintf("fraction p<0.05 = %.4f (B %d)", nc2_frac, B_NC2), sprintf("fraction p<0.05 = %.4f (B %d); regenerated-score p %.4g", nc2b_frac, B_NC2B, p_nc2b)),
                    prediction = c("exact reproduction", "0 discordant", "logFC < 0", "logFC > 0", "both in adult direction, >=1 nominal p<0.05",
                                   "median t in -1.5 to 1.5", "0.035 to 0.065", "0.035 to 0.065"),
                    met = c(pcA_pass, length(disc) == 0, pcc_dir_ok[1], pcc_dir_ok[2], pcc_met, abs(nc1_q[2]) <= 1.5,
                            nc2_frac >= 0.035 && nc2_frac <= 0.065, nc2b_frac >= 0.035 && nc2b_frac <= 0.065))
write.table(ctrls, file.path(OUT, "controls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gsm = ids, participant = pid_all[keep_s], histology = hist, nafld = nafld, sex = sex, age = age, bmi = bmi,
                       score = score, up_only = s_up, down_only = s_dn, unweighted = s_unw, max_sd_probe = s_d8),
            file.path(OUT, "scores_GSE66676.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(column = names(score18), nas = nas, nash_status = nash_status, score = score18),
            file.path(OUT, "scores_GSE185051_children.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(t_fl = t_fl, t_nc2 = t_nc2, t_nc2b = t_nc2b, nc1_t = nc1_t), file.path(OUT, "null_distributions.rds"))

## ---------------------------------------------------------------- supplementary panels (PDF, 6 pt, no bold)
STATE$stage <- "figures"
offs <- function(v) { r <- rank(v, ties.method = "first"); ((r - 1) %% 7 - 3) * 0.07 }
th <- theme_classic(base_size = 6) + theme(text = element_text(size = 6, colour = "black", face = "plain"),
  axis.text = element_text(size = 6, colour = "black"), axis.title = element_text(size = 6, colour = "black"),
  legend.position = "none", axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25))
fa <- data.frame(hist = factor(hist, levels = hl), score = score)
fa$x <- as.integer(fa$hist) + ave(fa$score, fa$hist, FUN = offs)
lab_a <- paste0(c("No\nNAFLD", "NAFLD\nnot NASH", "Borderline\nNASH", "Definite\nNASH"), "\n(n ", as.integer(table(fa$hist)), ")")
meds <- data.frame(x = seq_along(hl), m = vapply(hl, function(g) median(score[hist == g]), 0))
cols <- c("No NAFLD" = "#9E9E9E", "NAFLD not NASH" = "#F4A674", "Borderline NASH" = "#C9265E", "Definite NASH" = "#C9265E")
pA <- ggplot(fa, aes(x, score)) + geom_hline(yintercept = 0, linewidth = 0.2, colour = "#9E9E9E", linetype = "dashed") +
  geom_point(aes(colour = hist, shape = hist), size = 0.9) +
  geom_segment(data = meds, aes(x = x - 0.3, xend = x + 0.3, y = m, yend = m), linewidth = 0.4, colour = "black") +
  scale_colour_manual(values = cols) + scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_x_continuous(breaks = seq_along(hl), labels = lab_a) +
  labs(x = "GSE66676 histology (adolescents, bariatric surgery)", y = "Adult signature score (mean z up - mean z down)") + th
pdf.options(useDingbats = FALSE)
ggsave(file.path(OUT, "H62_panelA_GSE66676_score_by_histology.pdf"), pA, width = 2.8, height = 2.3, device = grDevices::pdf, useDingbats = FALSE)
fb <- data.frame(nas = nas, score = score18); fb$x <- fb$nas + ave(fb$score, fb$nas, FUN = offs) * 0.8
medb <- aggregate(score ~ nas, fb, median)
pB <- ggplot(fb, aes(x, score)) + geom_hline(yintercept = 0, linewidth = 0.2, colour = "#9E9E9E", linetype = "dashed") +
  geom_point(size = 0.9, colour = "#C9265E", shape = 16) +
  geom_segment(data = medb, aes(x = nas - 0.3, xend = nas + 0.3, y = score, yend = score), linewidth = 0.4, colour = "black") +
  scale_x_continuous(breaks = 0:7) +
  labs(x = "NAS (GSE185051 children with NAFLD)", y = "Adult signature score (mean z up - mean z down)") + th
ggsave(file.path(OUT, "H62_panelB_GSE185051_score_vs_NAS.pdf"), pB, width = 2.6, height = 2.3, device = grDevices::pdf, useDingbats = FALSE)
cap <- c(sprintf("Panel A. Per-participant adult signature score in GSE66676 (%d adolescents; points are participants; black bars are medians; z against the %d No NAFLD participants). NAFLD vs No NAFLD adjusted for sex, age and BMI: difference %.3f (95%% CI %.3f to %.3f), d %.2f, governing P %s.",
                 N, N0, b_T1, ci_T1[1], ci_T1[2], d_adj, fp(p_gov)),
         sprintf("Panel B. Score in %d GSE185051 children with NAFLD against NAS (z against the 52 children). Spearman rho %.2f, nominal P %s. The five adult controls are not shown and are not compared.",
                 length(score18), sp6$estimate, fp(sp6$p.value)))
message(paste(cap, collapse = "\n")); writeLines(cap, file.path(OUT, "panel_captions.txt"))

## ---------------------------------------------------------------- RESULT
dfmt <- function(x) trimws(formatC(x, format = "g", digits = 4))
res <- c("# RESULT H62", "",
  if (DEBUG) c("**DEBUG-NULL RUN: GSE66676 histology and GSE185051 NAS were permuted (seed 1); B reduced. Not a result.**", "") else NULL,
  paste0("Run: ", OUT), paste0("Prespec: ", file.path(RUN, "PRESPEC_H62.md"), " (sha256 ", PRESPEC_SHA, ", verified)"),
  paste0("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE)), "",
  paste0("## Verdict: ", verdict), "",
  paste0("Rule (prespec Reading): pass if governing p < 0.05 and coefficient > 0 with PC-A and PC-B passing and PC-C not failed; ",
         "opposite if governing p < 0.05 and coefficient < 0; otherwise null."), "",
  "## T1 (family of 1, alpha 0.05, no correction)", "",
  sprintf("- n: %d No NAFLD, %d NAFLD (participants, one biopsy each); model score ~ nafld + sex + age + bmi, df %d", N0, N1, DF),
  sprintf("- coefficient %.4f (SE %.4f), 95%% CI %.4f to %.4f, t %.3f", b_T1, se_T1, ci_T1[1], ci_T1[2], t_obs),
  sprintf("- d_adj %.3f (95%% CI %.3f to %.3f); Hedges' g unadjusted %.3f; pooled within-group SD %.4f", d_adj, d_adj_ci[1], d_adj_ci[2], hedges_g, sp),
  sprintf("- p parametric %s; Freedman-Lane p %s (B %d, floor %.4g); NC-2b regenerated-score p %s (B %d)", dfmt(p_param), dfmt(p_fl), B_FL, 1 / (B_FL + 1), dfmt(p_nc2b), B_NC2B),
  sprintf("- governing p %s: %s", dfmt(p_gov), gov_rule),
  sprintf("- 80%% power MDE for this n (noncentral t, df %d): d %.3f; median adult NAFL shift 0.97 %s by the d_adj CI", DF, mde, if (excl097) "excluded" else "not excluded"),
  sprintf("- NC-1 competitive rank %.3f; observed t %s the NC-1 95th percentile (%.3f)", nc1_rank, if (nc1_above95) "above" else "not above", nc1_q[3]), "",
  "## Controls", "", sprintf("- S4 (PC-B): %d discordant participant(s) removed before scoring: %s. Analysed n %d of %d.", length(disc), if (length(disc)) paste0(gsm[disc], " (", pid_all[disc], ", recorded ", sex_all[disc], ", ", hist_all[disc], ")", collapse = "; ") else "none", N, length(gsm)),
  paste0("- ", ctrls$control, ": ", ctrls$result, " | predicted ", ctrls$prediction, " | ", ifelse(ctrls$met, "met", "not met")), "",
  "## Descriptive outputs (nominal, outside the family; never called significant)", "",
  paste0("- ", Dall$output, " ", Dall$item, ": ", vapply(Dall$value, function(v) dfmt(as.numeric(v)), "")), "",
  "## Proposed text (returned as data; nothing installed)", "",
  "- Insert after `docs/manuscript/draft/fig3.md:40` (end of the leave-out paragraph), and one line in `docs/RESULTS.md` after line 47:",
  paste0("  \"", sentence, "\""),
  paste0("- GSE185051 clause (descriptive): \"in ", length(score18), " children the score correlated with NAS at rho ", f2(sp6$estimate), ", nominal P ", fp(sp6$p.value), "\""), "",
  "## Files", "", paste0("- ", list.files(OUT)))
writeLines(res, RESULT_FILE)
if (!DEBUG) { rt <- file.path(RUN, "RESULT_H62.md"); if (file.exists(rt)) stop("RUN/RESULT_H62.md exists; not overwriting"); file.copy(RESULT_FILE, rt) }
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
lg("H62 done: verdict ", verdict, " ; governing p ", dfmt(p_gov), " ; ", format(Sys.time(), tz = "UTC", usetz = TRUE))
