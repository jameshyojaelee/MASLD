# H51 shared definitions. Sourced by stage1_blinded.R and stage2_outcome.R.
suppressMessages({library(data.table); library(edgeR); library(limma); library(parallel)})

H <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H51-20261015T112617Z"
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
REL <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
RR <- file.path(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z")
DEC <- file.path(ROOT, "Analysis/Deconvolution/corrected_774_20261009T141331Z/results")
PRESPEC_SHA <- "9496f74d451834914e2aebf5053c8201994c8cc81cbc370edd30369b69e32165"
SEED <- 20260923L
DISC <- c("GSE130970", "GSE135251", "GSE162694")
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
DEBUG <- Sys.getenv("H51_DEBUG", "0") == "1"
B_PERM <- if (DEBUG) 8L else 1000L
B_BOOT <- if (DEBUG) 4L else 200L

INPUTS <- data.table(
  id = c("dge", "meta", "stage_ext", "sra", "pairs", "roster", "frac_GSE130970", "frac_GSE135251",
          "frac_GSE162694", "frac_GSE213621", "run3831", "sum3831", "n3831", "acc_bins"),
  path = c(file.path(REL, "merged_dge.rds"), file.path(REL, "meta_matched.rds"),
           file.path(ROOT, "figures/candidates/fig3e-s3-774-20261011T151616Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv"),
           file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/metadata/SraRunTable.csv"),
           file.path(ROOT, "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/full_filtered_related_pairs.tsv"),
           file.path(ROOT, "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/resource_identity_roster.tsv"),
           file.path(DEC, "GSE130970/GSE130970_bayesprism_proportions.tsv"),
           file.path(DEC, "GSE135251/GSE135251_bayesprism_proportions.tsv"),
           file.path(DEC, "GSE162694/GSE162694_bayesprism_proportions.tsv"),
           file.path(DEC, "GSE213621/GSE213621_bayesprism_proportions.tsv"),
           file.path(RR, "jobs/3831/run.R"), file.path(RR, "jobs/3831/out/discovery_grouped_summary.tsv"),
           file.path(RR, "jobs/3831/out/discovery_grouped_n.tsv"),
           file.path(RR, "scratch/validity-L10-03/accession_order_bins.tsv")),
  sha = c("73f3c53fcb9c847e98b0e4f0918c52a04a6b0b21744e5b448cb45e5cbcfcad3b",
          "429e27376ad421f5f5fae19c75f5e760030e727fab65c5768ef562a9939edfc3",
          "9e50bee20dab39f3e9e48b337769aab6af72d058d67bbda2cb071f053ddcd3b0",
          "025d98d6e90cc3194efbed92fa399cd9ab59ea18442885e9461a4922089a346d",
          "704697656622ef664f45aa19ed8779fc2135c242e5b3fad239988e96f5e412ea",
          "d6c244a5400c231f213b12564353ded438f56f40fb062ed07c86c3956bd0224f",
          "bbf57df3fbd076159934e8b8260d71674b41ef96b366180ab1d23b4f7f8d6f0f",
          "22c8f7aef4b47f39e979e1e7db2ba3a4fa2d2f585243505669aa0b089d2bb643",
          "99bcbc0305b892d203a30ea89f3974a383934ae795de2e2d33a5631d28c77b50",
          "cc91e212cfd7c874c3773ab5ab76e1fafd068b36d20196d7baa4b4cf04bd9e5a",
          "732a989ae771be4d0aa9186bdf8d14d62134182234c4d1ca321b55c19fdbf640",
          "71c8f1f0f496c1b892f3a9b692357649dd9fa0d03f275b22dc72196a51095c15",
          "290d9c465402a97c6c11995edb599e122a984664f9023104840051c7c98cac4d",
          "11cfd150a54c3286134830db6cd47f22751d2d02d24cac4a4e00ebb9f70a57c3"))
NONINPUT <- c(path = file.path(ROOT, "figures/candidates/continuum-774-panels-20261011T165944Z/stage_continuum_comparison/source_tables/transcript_stage_trend_cohort.tsv.gz"),
              sha = "a6d9db21465cd69321bee02cec5adc691aab049005e7cbb5bb0ad4af5ca6ba15")

sha256 <- function(f) sub(" .*", "", system2("sha256sum", shQuote(f), stdout = TRUE))

verify_hashes <- function(outdir, scripts) {
  ps <- sha256(file.path(H, "PRESPEC_H51.md"))
  if (ps != PRESPEC_SHA) stop("PRESPEC sha256 mismatch: ", ps)
  got <- vapply(INPUTS$path, sha256, "")
  chk <- data.table(INPUTS[, .(id, path, expected = sha)], observed = got)
  chk[, match := expected == observed]
  fwrite(chk, file.path(outdir, "input_sha256_check.tsv"), sep = "\t")
  if (!all(chk$match)) stop("Input sha256 mismatch: ", paste(chk[match == FALSE, id], collapse = ","))
  ni <- sha256(NONINPUT[["path"]])
  cat("Non-input overlap table sha256", if (ni == NONINPUT[["sha"]]) "matches" else "DIFFERS (recorded only)", "\n")
  sc <- data.table(script = scripts, sha256 = vapply(file.path(H, scripts), sha256, ""))
  fwrite(sc, file.path(outdir, "script_sha256.tsv"), sep = "\t")
  cat("HASH CHECK PASSED: prespec + ", nrow(chk), " inputs\n", sep = "")
}

bin_of_stage <- function(fs) fifelse(fs <= 1L, "F0F1", fifelse(fs == 2L, "F2", "F3F4"))
cond_to_bin <- c(Fibrosis_F0F1 = "F0F1", Fibrosis_F2 = "F2", Fibrosis_F3F4 = "F3F4")
BINS <- c("F0F1", "F2", "F3F4")
CF <- c(F3F4 = "binF3F4", F2 = "binF2")
SEXCF <- "factor(inferred_sex)M"

load_inputs <- function() {
  dge_all <- readRDS(INPUTS[id == "dge", path])
  meta <- as.data.table(readRDS(INPUTS[id == "meta", path]))
  stopifnot(identical(dim(dge_all), c(24115L, 774L)))
  ann <- unique(fread(INPUTS[id == "stage_ext", path], select = c("gene_id_versioned", "gene_name", "chromosome")))
  stopifnot(nrow(ann) == 24115L)
  ann <- ann[match(rownames(dge_all), gene_id_versioned)]
  stopifnot(identical(ann$gene_id_versioned, rownames(dge_all)))
  # Discovery set
  disc <- meta[dataset %in% DISC & group_binary == "Disease" & fibrosis_stage %in% 0:4 & !is.na(inferred_sex)]
  disc[, bin := factor(bin_of_stage(fibrosis_stage), levels = BINS)]
  setorder(disc, sample_id)
  # GSE213621 replication set; bin from `condition`
  g <- meta[dataset == "GSE213621" & group_binary == "Disease"]
  g[, bin := factor(unname(cond_to_bin[condition]), levels = BINS)]
  stopifnot(!anyNA(g$bin), !anyNA(g$inferred_sex))
  setorder(g, sample_id)
  g[, acc := as.numeric(sub("SRR", "", sample_id))]
  g[, acc_rank := frank(acc)]
  g[, acc_q := cut(acc_rank, quantile(acc_rank, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))]
  g[, stratum := paste(acc_q, inferred_sex, sep = "_")]
  g[, cell := paste(bin, inferred_sex, sep = "_")]
  disc[, cell := paste(bin, inferred_sex, sep = "_")]
  list(dge_all = dge_all, meta = meta, ann = ann, disc = disc, g = g)
}

subdge <- function(dge_all, ids, genes = NULL) {
  d <- dge_all[, match(ids, colnames(dge_all))]
  if (!is.null(genes)) d <- d[genes, , keep.lib.sizes = TRUE]
  d
}

# voomWithQualityWeights -> lmFit -> eBayes, limma defaults (as job 3831)
fit_qw <- function(d, des) {
  v <- voomWithQualityWeights(d, des, plot = FALSE)
  eBayes(lmFit(v, des))
}

ols_slope <- function(y, x) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2) }
conc <- function(y, x) mean(sign(y) == sign(x))

# Design builders
des_D0 <- function(info) model.matrix(~ factor(dataset) + factor(inferred_sex) + bin, info)
des_G0 <- function(info) model.matrix(~ factor(inferred_sex) + bin, info)

# 7-part lineage coding and ILR (hand-coded Helmert basis, -helmert(7).T)
LINEAGE <- list(hep_chol = c("Hepatocytes", "Cholangiocytes"), endothelial = "Endothelial cells",
                fibroblasts = "Fibroblasts", mac_mono = c("Macrophages", "Mono+mono derived cells"),
                other_myeloid = c("cDC1s", "cDC2s", "pDCs", "Neutrophils", "Basophils"),
                t_nk = c("T cells", "Circulating NK/NKT", "Resident NK"), b_plasma = c("B cells", "Plasma cells"))
ilr_basis <- function(D = 7L) {
  V <- matrix(0, D, D - 1L)
  for (k in seq_len(D - 1L)) { V[seq_len(k), k] <- -1 / sqrt(k * (k + 1)); V[k + 1L, k] <- k / sqrt(k * (k + 1)) }
  V
}
load_ilr <- function(ds, ids) {
  f <- read.delim(INPUTS[id == paste0("frac_", ds), path], check.names = FALSE)
  stopifnot(setequal(unlist(LINEAGE), colnames(f)), all(ids %in% rownames(f)))
  f <- as.matrix(f[ids, , drop = FALSE])
  P <- sapply(LINEAGE, function(cols) rowSums(f[, cols, drop = FALSE]))
  P <- P / rowSums(P)
  below <- mean(P < 1e-6)
  P <- pmax(P, 1e-6); P <- P / rowSums(P)
  L <- log(P); clr <- L - rowMeans(L)
  ilr <- clr %*% ilr_basis(7L)
  colnames(ilr) <- paste0("ilr", 1:6); rownames(ilr) <- ids
  list(ilr = ilr, share_below_floor = below)
}

z975 <- qnorm(0.975); z80 <- qnorm(0.80)
MDE_SLOPE_K <- z975 + z80                # S3 and gate MDE
MDE_GENE_K <- qnorm(1 - 0.005 / 2) + z80 # per-gene MDE_g = 3.65 x SE_g
# Table formatting for reports: logicals as text, whole numbers unrounded, others to `digits` significant digits
fmt <- function(dt, digits = 4) as.data.table(dt)[, lapply(.SD, function(v) if (is.logical(v)) as.character(v) else
  if (is.numeric(v) && any(abs(v - round(v)) > 0, na.rm = TRUE)) signif(v, digits) else v)]
log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n")
