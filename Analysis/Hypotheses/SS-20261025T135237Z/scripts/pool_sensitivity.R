# SS POOL (PRESPEC_SS_addendum_POOL_sensitivity.md): pooled sensitivity of the Figure 3E adjacent steps.
# Usage: Rscript pool_sensitivity.R <new out dir>
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table); library(parallel); library(digest)})
SS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT <- commandArgs(trailingOnly = TRUE)[1]
if (file.exists(OUT)) stop("Output dir exists: ", OUT); dir.create(OUT, recursive = TRUE)
log_msg <- function(...) cat(format(Sys.time(), "%F %T"), ..., "\n")
p <- readLines(file.path(SS, "PRESPEC_POOL.sha256")); h <- sub(" .*", "", p)
stopifnot(identical(h, digest(file = file.path(SS, "PRESPEC_SS_addendum_POOL_sensitivity.md"), algo = "sha256")))
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

DISC <- file.path(ROOT, "figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions")
dA <- readRDS(file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/merged_dge.rds"))
man <- unique(fread(file.path(DISC, "stage_sample_manifest.tsv"))[, .(sample_id, dataset, inferred_sex, group = fib_group)])
stopifnot(!anyDuplicated(man$sample_id), nrow(man) == 387L)
man[, stage := as.integer(sub("F", "", group))][, source := "discovery"][, batch := dataset]

dB <- readRDS(file.path(SS, "out_s2_qc/dge_included.rds"))
stopifnot(identical(rownames(dA), rownames(dB)))
dB <- calcNormFactors(dB, method = "RLE")
u <- fread(file.path(SS, "out_s2_qc/units.tsv"))[include == TRUE]
l <- sub("\r$", "", readLines(file.path(ROOT, "data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv")))
st <- fread(text = paste(l, collapse = "\n"), colClasses = "character")[, .(sample_alias, stage = suppressWarnings(as.integer(nash_crn_kleiner_fibrosis_stage)))]
u <- merge(u, st, by = "sample_alias")[!is.na(stage) & !is.na(inferred_sex)]
stopifnot(nrow(u) == 493L)
ssu <- u[, .(sample_id = sample_alias, dataset = paste0("SS_", batch), inferred_sex, stage, source = "steatosite", batch)]
# SteatoSITE code 1 is F0-F1; it enters only the grouped step
ssu[, group := fifelse(stage == 1L, "F0F1", paste0("F", stage))]

counts <- cbind(dA$counts[, man$sample_id], dB$counts[, ssu$sample_id])
eff <- c(dA$samples[man$sample_id, "lib.size"], dB$samples[ssu$sample_id, "lib.size"])
nf <- c(dA$samples[man$sample_id, "norm.factors"], dB$samples[ssu$sample_id, "norm.factors"])
ALL <- DGEList(counts, lib.size = eff, norm.factors = nf)
U <- rbind(man[, .(sample_id, dataset, inferred_sex, stage, source)], ssu[, .(sample_id, dataset, inferred_sex, stage, source)])

# step membership: lo/hi stage sets per source
STEPS <- list(F0F1_to_F2 = list(lo = c(0L, 1L), hi = 2L), F2_to_F3 = list(lo = 2L, hi = 3L), F3_to_F4 = list(lo = 3L, hi = 4L))
units_for <- function(step, set, drop = NA_character_) {
  s <- STEPS[[step]]; x <- U[stage %in% c(s$lo, s$hi)]
  if (set == "D") x <- x[source == "discovery"]; if (set == "S") x <- x[source == "steatosite"]
  if (!is.na(drop)) x <- x[dataset != drop]
  x[, group := factor(fifelse(stage %in% s$hi, "hi", "lo"), levels = c("lo", "hi"))]
  setorder(x, sample_id); x }

fit_one <- function(step, set, drop = NA_character_, interaction = FALSE) {
  x <- units_for(step, set, drop)
  info <- data.frame(dataset = droplevels(factor(x$dataset)), inferred_sex = droplevels(factor(x$inferred_sex)),
                     group = x$group, ss = as.integer(x$source == "steatosite"), row.names = x$sample_id)
  des <- if (interaction) model.matrix(~ dataset + inferred_sex + group + group:ss, info) else
    model.matrix(~ dataset + inferred_sex + group, info)
  if (interaction) des <- des[, colnames(des) != "grouplo:ss", drop = FALSE]
  stopifnot(qr(des)$rank == ncol(des))
  v <- voomWithQualityWeights(ALL[, x$sample_id], des, plot = FALSE); f <- eBayes(lmFit(v, des))
  cf <- if (interaction) "grouphi:ss" else "grouphi"
  tt <- topTable(f, coef = cf, number = Inf, sort.by = "none")
  w <- v$weights; ssw <- if (any(info$ss == 1L) && any(info$ss == 0L)) rowSums(w[, info$ss == 1L, drop = FALSE]) / rowSums(w) else rep(NA_real_, nrow(w))
  list(key = paste(step, set, ifelse(is.na(drop), "all", paste0("minus_", drop)), ifelse(interaction, "interaction", "main"), sep = "|"),
       res = data.table(gene = rownames(tt), logFC = tt$logFC, P = tt$P.Value, q = p.adjust(tt$P.Value, "BH")),
       n = x[, .N, by = .(dataset_src = fifelse(source == "steatosite", "SteatoSITE", dataset), group)],
       ss_weight_share = median(ssw)) }

jobs <- list()
for (s in names(STEPS)) {
  for (set in c("D", "P", "S")) jobs[[length(jobs) + 1L]] <- list(s, set, NA_character_, FALSE)
  for (c in c("GSE130970", "GSE135251", "GSE162694")) for (set in c("D", "P")) jobs[[length(jobs) + 1L]] <- list(s, set, c, FALSE)
  jobs[[length(jobs) + 1L]] <- list(s, "P", NA_character_, TRUE) }
log_msg("fits:", length(jobs))
R <- mclapply(jobs, function(j) tryCatch(fit_one(j[[1]], j[[2]], j[[3]], j[[4]]), error = function(e) list(key = paste(unlist(j), collapse = "|"), error = conditionMessage(e))),
              mc.cores = NC, mc.preschedule = FALSE)
names(R) <- vapply(R, `[[`, "", "key")
err <- R[vapply(R, function(r) !is.null(r$error), TRUE)]
if (length(err)) { print(err); log_msg("fits with errors:", length(err)) }
R <- R[vapply(R, function(r) is.null(r$error), TRUE)]

# Gate: D reproduces frozen discovery for the exact steps
fr <- fread(file.path(DISC, "stage_extension_all_gene_results.tsv"), select = c("gene_id_versioned", "axis", "contrast", "logFC", "FDR"))[axis == "fibrosis_adjacent"]
gate <- rbindlist(lapply(c("F2_to_F3", "F3_to_F4"), function(s) {
  d <- R[[paste(s, "D|all|main", sep = "|")]]$res; z <- merge(d, fr[contrast == s], by.x = "gene", by.y = "gene_id_versioned")
  data.table(step = s, n = nrow(z), max_abs_diff = max(abs(z$logFC.x - z$logFC.y)), calls_refit = sum(z$q < 0.05), calls_frozen = sum(z$FDR < 0.05)) }))
fwrite(gate, file.path(OUT, "gate_reproduction.tsv"), sep = "\t"); print(gate)
if (any(gate$max_abs_diff >= 1e-6) || any(gate$calls_refit != gate$calls_frozen)) stop("GATE FAILED: discovery refit does not reproduce frozen results")

S <- rbindlist(lapply(names(R), function(k) { r <- R[[k]]; kk <- strsplit(k, "\\|")[[1]]
  data.table(step = kk[1], set = kk[2], omit = kk[3], term = kk[4], n_lo = sum(r$n[group == "lo", N]), n_hi = sum(r$n[group == "hi", N]),
             calls = sum(r$res$q < 0.05), up = sum(r$res$q < 0.05 & r$res$logFC > 0), down = sum(r$res$q < 0.05 & r$res$logFC < 0),
             ss_weight_share_median = r$ss_weight_share) }))
cmp <- rbindlist(lapply(names(STEPS), function(s) rbindlist(lapply(c("all", "minus_GSE130970", "minus_GSE135251", "minus_GSE162694"), function(o) {
  d <- R[[paste(s, "D", o, "main", sep = "|")]]$res; p <- R[[paste(s, "P", o, "main", sep = "|")]]$res
  if (is.null(d) || is.null(p)) return(NULL); z <- merge(d, p, by = "gene", suffixes = c(".D", ".P"))
  pc <- z[q.P < 0.05]; dc <- z[q.D < 0.05]
  data.table(step = s, omit = o, D_calls = nrow(dc), P_calls = nrow(pc), both = sum(z$q.D < 0.05 & z$q.P < 0.05),
             P_calls_with_D_sign = mean(sign(pc$logFC.P) == sign(pc$logFC.D)),
             D_calls_kept_in_P_same_sign = sum(dc$q.P < 0.05 & sign(dc$logFC.P) == sign(dc$logFC.D)),
             spearman_all_genes = cor(z$logFC.D, z$logFC.P, method = "spearman")) }))))
fwrite(S, file.path(OUT, "fit_summary.tsv"), sep = "\t"); fwrite(cmp, file.path(OUT, "pooled_vs_discovery.tsv"), sep = "\t")
fwrite(rbindlist(lapply(names(R), function(k) R[[k]]$n[, key := k])), file.path(OUT, "units_per_fit.tsv"), sep = "\t")
fwrite(rbindlist(lapply(names(R), function(k) R[[k]]$res[, key := k])), file.path(OUT, "all_gene_results.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(S); print(cmp); log_msg("done")
