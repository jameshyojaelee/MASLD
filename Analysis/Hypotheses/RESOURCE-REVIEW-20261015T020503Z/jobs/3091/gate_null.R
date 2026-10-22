#!/usr/bin/env Rscript
# Feasibility gate for candidate L03-03 (agent gate-L03-03, job 3091).
# Computes, per lineage, ONLY:
#   - donor n after the >=50 lineage-cell gate and sex availability
#   - how many of the 1,854 bulk calls pass filterByExpr
#   - the within-dataset stage-permutation NULL of the sign-agreement fraction
#     (end to end: voom + lmFit per permutation; filter is permutation-invariant
#     because within-dataset permutation keeps the per-stage donor counts)
#   - the positive control (COL1A1 slope in fibroblasts) and its permutation p
#   - per-lineage donor ambient load (decontX, N3 outputs) vs stage, nominal
# It does NOT compute the observed sign-agreement fraction (primary readout stays blind).
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table); library(parallel)})
set.seed(20261015)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT <- "out"
t0 <- Sys.time()
LINEAGES <- c("Hepatocytes", "Endothelial_cells", "Fibroblasts", "Cholangiocytes")
CT_NAME <- c(Hepatocytes = "Hepatocytes", Endothelial_cells = "Endothelial cells",
             Fibroblasts = "Fibroblasts", Cholangiocytes = "Cholangiocytes")
B_MAX <- 1000L
BUDGET_SEC <- 55  # per lineage permutation budget
NCORE <- 4L

# ---- roster, sex, run -> donor
ros <- fread(file.path(BASE, "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/current_source_diagnosis_donor_roster.tsv"))
ros <- ros[toupper(as.character(exclude_stage_analysis)) != "TRUE"]
stopifnot(nrow(ros) == 73)
ord <- c(Healthy = 0, Steatosis = 1, Steatohepatitis = 2, Cirrhosis = 3)
ros[, stage := ord[disease_stage_coarse]]
n3 <- fread(file.path(BASE, "Analysis/Hypotheses/N3-source-diagnosis-20261014/corrected/donor_meta.tsv"))
n3[, donor2 := sub("^(GSE174748|GSE189600)_", "", donor)]
ros <- merge(ros, n3[, .(donor = donor2, sex_call)], by = "donor", all.x = TRUE)

run2donor <- c()
for (ds in c("GSE244832", "GSE202379", "GSE185477")) {
  p <- fread(file.path(BASE, sprintf("data/%s/metadata/donor_pairing.csv", ds)), colClasses = "character")
  for (i in seq_len(nrow(p))) {
    s <- strsplit(p$rna_srrs[i], ";", fixed = TRUE)[[1]]
    run2donor <- c(run2donor, setNames(rep(paste0(ds, "_", p$donor_id[i]), length(s)), s))
  }
}
md <- fread(file.path(BASE, "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected.tsv"))
md[, donor := ifelse(sample %in% names(run2donor), run2donor[sample], sample)]
cellcols <- paste0("n_", LINEAGES)
cnt <- md[, lapply(.SD, sum, na.rm = TRUE), by = .(donor, dataset), .SDcols = cellcols]
ros <- merge(ros, cnt, by = c("donor", "dataset"), all.x = TRUE)

# ---- bulk calls (1,854) and gene map
co <- fread(file.path(BASE, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/coefficient_table.tsv.gz"))
co[, ens := sub("\\..*$", "", gene)]
calls <- co[padj < 0.05 & abs(logFC) > 0.5]
stopifnot(nrow(calls) == 1854)
gmap <- fread(cmd = paste("zcat", file.path(BASE, "Analysis/Hypotheses/N3-source-diagnosis-20261014/corrected/hep_donor_counts.tsv.gz"), "| cut -f1,2"))
setnames(gmap, c("symbol", "ensembl"))

# ---- per-lineage donor ambient load from N3 per-library decontX summaries
lq <- fread(file.path(BASE, "Analysis/Hypotheses/N3-ambient-20261013T060951Z/library_ambient_qc.tsv"))
amb <- rbindlist(lapply(lq$library, function(lib) {
  f <- file.path(BASE, "Analysis/Hypotheses/N3-ambient-20261013T060951Z/per_library", paste0(lib, ".rds"))
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f)$ct_summary
  data.table(library = lib, cell_type = x$cell_type, raw = x$raw_umis, dec = x$dec_umis, n = x$n_cells)
}))
amb <- merge(amb, lq[, .(library, donor)], by = "library")
amb[, donor := sub("^(GSE174748|GSE189600)_", "", donor)]
amb_d <- amb[, .(raw = sum(raw), dec = sum(dec)), by = .(donor, cell_type)]
amb_d[, contam := 1 - dec / raw]
fwrite(amb_d, file.path(OUT, "donor_lineage_ambient.tsv"), sep = "\t")
cat("ambient libraries read:", uniqueN(amb$library), " donors:", uniqueN(amb$donor), "\n")

summ <- list(); nulls <- list(); amb_res <- list(); pc <- list()
for (L in LINEAGES) {
  tl <- Sys.time()
  el <- ros[get(paste0("n_", L)) >= 50]
  n_gate <- nrow(el)
  el <- el[sex_call %in% c("Male", "Female")]
  # drop datasets with a single stage level among eligible donors? keep them (they add residual df only)
  pb <- fread(file.path(BASE, sprintf("Analysis/SingleCell/results_gpu_v2/pseudobulk/%s_pseudobulk.csv", L)))
  genes <- pb[[1]]; pb <- as.matrix(pb[, -1]); rownames(pb) <- genes
  stopifnot(identical(genes, gmap$symbol))
  runs <- colnames(pb)
  rdonor <- ifelse(runs %in% names(run2donor), run2donor[runs], runs)
  keep <- rdonor %in% el$donor
  M <- t(rowsum(t(pb[, keep]), group = rdonor[keep]))
  el <- el[match(colnames(M), el$donor)]
  stopifnot(identical(el$donor, colnames(M)))
  el[, dataset := factor(dataset)]; el[, sex := factor(sex_call)]
  design <- model.matrix(~ dataset + sex + stage, data = el)
  y <- DGEList(M)
  fk <- filterByExpr(y, group = el$disease_stage_coarse)
  y <- y[fk, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  ens <- gmap$ensembl[match(rownames(y), gmap$symbol)]
  idx <- which(ens %in% calls$ens)
  bsign <- sign(calls$logFC[match(ens[idx], calls$ens)])
  # within-dataset permutations (seeded), regenerated end to end
  perms <- lapply(seq_len(B_MAX), function(b) {
    s <- el$stage
    for (d in levels(el$dataset)) { w <- which(el$dataset == d); if (length(w) > 1) s[w] <- s[w][sample.int(length(w))] }
    s
  })
  fit_one <- function(s) {
    el2 <- copy(el); el2[, stage := s]
    des <- model.matrix(~ dataset + sex + stage, data = el2)
    v <- voom(y, des)
    f <- eBayes(lmFit(v, des))
    sl <- f$coefficients[, "stage"]
    c(agree = mean(sign(sl[idx]) == bsign),
      frac_p05 = mean(f$p.value[, "stage"] < 0.05),
      col1a1 = if ("COL1A1" %in% rownames(f)) unname(sl["COL1A1"]) else NA_real_)
  }
  done <- 0L; res <- list()
  while (done < B_MAX && as.numeric(difftime(Sys.time(), tl, units = "secs")) < BUDGET_SEC) {
    bb <- (done + 1L):min(done + 4L * NCORE, B_MAX)
    res <- c(res, mclapply(perms[bb], fit_one, mc.cores = NCORE))
    done <- max(bb)
  }
  R <- do.call(rbind, res)
  nulls[[L]] <- data.table(lineage = L, perm = seq_len(nrow(R)), R)
  # positive control (observed) only for fibroblasts; null slope for COL1A1 kept for all
  if (L == "Fibroblasts") {
    v <- voom(y, design); f <- eBayes(lmFit(v, design))
    tt <- topTable(f, coef = "stage", number = Inf, sort.by = "none")
    obs <- tt["COL1A1", ]
    pperm <- (1 + sum(R[, "col1a1"] >= obs$logFC)) / (nrow(R) + 1)
    pc[[L]] <- data.table(lineage = L, gene = "COL1A1", passes_filter = TRUE, logFC_per_step = obs$logFC,
                          t = obs$t, limma_p_two_sided = obs$P.Value, perm_p_one_sided = pperm,
                          n_perm = nrow(R), aveexpr = obs$AveExpr,
                          col1a1_bulk_logFC = co[symbol == "COL1A1", logFC])
  }
  # ambient load vs stage, nominal OLS (dataset adjusted)
  a <- merge(el[, .(donor, dataset, stage, sex)], amb_d[cell_type == CT_NAME[[L]], .(donor, contam)], by = "donor")
  a <- a[contam > 0 & contam < 1]
  if (nrow(a) > 5) {
    m <- lm(qlogis(contam) ~ dataset + stage, data = a)
    cf <- summary(m)$coefficients["stage", ]
    amb_res[[L]] <- data.table(lineage = L, n_donors_with_ambient = nrow(a), logit_contam_per_step = cf[1],
                               se = cf[2], p_nominal = cf[4], median_contam = median(a$contam))
  }
  tab <- el[, .N, by = .(dataset, disease_stage_coarse)]
  fwrite(tab, file.path(OUT, sprintf("n_%s.tsv", L)), sep = "\t")
  q <- quantile(R[, "agree"], c(0.025, 0.05, 0.5, 0.95, 0.975, 0.9875, 0.99))
  summ[[L]] <- data.table(lineage = L, n_gate50 = n_gate, n_fit = nrow(el),
                          n_datasets = uniqueN(el$dataset), genes_pass_filter = nrow(y),
                          calls_expressed = length(idx), calls_up = sum(bsign > 0), calls_down = sum(bsign < 0),
                          n_perm_done = nrow(R), null_mean = mean(R[, "agree"]), null_sd = sd(R[, "agree"]),
                          q025 = q[1], q05 = q[2], q50 = q[3], q95 = q[4], q975 = q[5], q9875 = q[6], q99 = q[7],
                          null_min = min(R[, "agree"]), null_max = max(R[, "agree"]),
                          mean_frac_nominal_p05 = mean(R[, "frac_p05"]),
                          secs = as.numeric(difftime(Sys.time(), tl, units = "secs")))
  cat(L, "done", nrow(R), "perms\n")
}
fwrite(rbindlist(summ), file.path(OUT, "null_summary.tsv"), sep = "\t")
fwrite(rbindlist(nulls), file.path(OUT, "null_draws.tsv"), sep = "\t")
fwrite(rbindlist(pc), file.path(OUT, "positive_control.tsv"), sep = "\t")
fwrite(rbindlist(amb_res), file.path(OUT, "ambient_vs_stage.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("total secs", as.numeric(difftime(Sys.time(), t0, units = "secs")), "\n")
