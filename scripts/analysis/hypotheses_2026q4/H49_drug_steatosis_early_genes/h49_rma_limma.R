# H49 step B: per-compound RMA and limma for Open TG-GATEs human in vitro (PHH).
#
# Usage:
#   Rscript --vanilla h49_rma_limma.R <zip_dir> <cel_attribute_csv> <cdf_tarball> <out_dir> [compound ...]
#
# For each compound archive: validate zip member names, extract the CEL files
# to a session temp directory, RMA over all arrays of that compound experiment,
# then limma on log2 probe-set intensities. RMA is computed in R (affy's
# R-level convolution background, quantile normalization, log2, stats::medpolish
# per probe set) because preprocessCore's threaded C routines fail in this
# environment (pthread_create error 22). PM intensities are read with
# affy::read.affybatch and affy::pm using the CDF data object held in this R session;
# affy is restricted to that object so it cannot download or install anything.
# limma:
#   main model    ~0 + group (group = dose x time); contrasts treated - Control at the same time
#   vehicle model Control-24h arrays split into two one-array groups; contrast A - B
# eBayes(trend = TRUE, robust = TRUE) in both models.
# Writes one RDS per compound with probe-set t, logFC and df.total per contrast,
# plus per-probe-set mean intensity over the compound's Control arrays.
# Inputs are untrusted downloads; nothing is executed from them.

suppressPackageStartupMessages({
  library(affy)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: h49_rma_limma.R <zip_dir> <cel_attribute_csv> <cdf_tarball> <out_dir> [compound ...]")
zip_dir <- args[1]; attr_csv <- args[2]; cdf_tar <- args[3]; out_dir <- args[4]
only <- if (length(args) > 4) args[-(1:4)] else character(0)
dir.create(file.path(out_dir, "per_compound"), recursive = TRUE, showWarnings = FALSE)

# CDF probe layout: load the data object only (package R code is not run).
cdf_tmp <- file.path(tempdir(), "cdf")
dir.create(cdf_tmp)
untar(cdf_tar, files = "hgu133plus2cdf/data/hgu133plus2cdf.rda", exdir = cdf_tmp)
cdf_env <- new.env()
loaded <- load(file.path(cdf_tmp, "hgu133plus2cdf/data/hgu133plus2cdf.rda"), envir = cdf_env)
if (!"hgu133plus2cdf" %in% loaded || !is.environment(cdf_env$hgu133plus2cdf)) stop("CDF object not found in rda")
# affy::getCdfInfo looks up the raw CEL-header name ("HG-U133_Plus_2"); keep lookup to this session only
assign("HG-U133_Plus_2", cdf_env$hgu133plus2cdf, envir = .GlobalEnv)
bioc_opt <- getOption("BioC")
bioc_opt$affy$probesloc <- list(list(what = "environment", where = .GlobalEnv))
options(BioC = bioc_opt)
cat("CDF probe sets:", length(ls(cdf_env$hgu133plus2cdf)), "\n")

# RMA in R: convolution background per array (affy:::bg.adjust), quantile
# normalization (limma::normalizeQuantiles, ties = TRUE), log2, median polish
# per probe set (stats::medpolish, maxiter 10, eps 0.01).
rma_r <- function(cel_paths) {
  ab <- affy::read.affybatch(filenames = cel_paths)
  # PM rows and their probe-set names, built the way affy::pm() orders them
  psets <- as.list(affy::getCdfInfo(ab))
  psets <- psets[order(names(psets))]
  pm <- Biobase::exprs(ab)[unlist(lapply(psets, function(p) p[, 1]), use.names = FALSE), , drop = FALSE]
  ps <- rep(names(psets), vapply(psets, nrow, 1L))
  stopifnot(identical(unname(pm), unname(affy::pm(ab))), length(unique(ps)) == length(psets))
  rm(ab)
  bg <- apply(pm, 2, affy:::bg.adjust)
  lg <- log2(limma::normalizeQuantiles(bg, ties = TRUE))
  idx <- split(seq_len(nrow(lg)), ps)
  ex <- t(vapply(idx, function(i) {
    if (length(i) == 1) return(lg[i, ])
    mp <- stats::medpolish(lg[i, , drop = FALSE], trace.iter = FALSE, na.rm = FALSE)
    mp$overall + mp$col
  }, numeric(ncol(lg))))
  colnames(ex) <- colnames(pm)
  stopifnot(nrow(ex) == length(psets))
  ex
}

att <- read.csv(attr_csv, colClasses = "character", fileEncoding = "latin1")
att <- att[att$SPECIES == "Human" & att$EXP_TEST_TYPE == "in vitro", ]
att$time <- sub(" hr$", "h", att$SACRIFICE_PERIOD)
att$group <- make.names(paste(att$DOSE_LEVEL, att$time, sep = "_"))

zips <- sort(list.files(zip_dir, pattern = "\\.Human\\.in_vitro\\.Liver\\.zip$", full.names = TRUE))
cmp_of_zip <- sub("\\.Human\\.in_vitro\\.Liver\\.zip$", "", basename(zips))
cmp_key <- function(s) gsub("[^a-z0-9]", "", tolower(s))
att$key <- cmp_key(att$COMPOUND_NAME)
if (length(only)) {
  keep <- cmp_key(cmp_of_zip) %in% cmp_key(only)
  zips <- zips[keep]; cmp_of_zip <- cmp_of_zip[keep]
}
cat("archives to process:", length(zips), "\n")

safe_member <- function(m) !grepl("(^/)|(\\.\\.)|(\\\\)", m)

# One compound experiment per call; returns a status string. Compounds run in
# parallel worker processes (one compound per worker; no shared state).
process_one <- function(k) {
  zf <- zips[k]; cname <- cmp_of_zip[k]; key <- cmp_key(cname)
  out_rds <- file.path(out_dir, "per_compound", paste0(key, ".rds"))
  if (file.exists(out_rds)) stop("refusing to overwrite ", out_rds)
  a <- att[att$key == key, ]
  if (!nrow(a)) return("no attribute rows")
  mem <- unzip(zf, list = TRUE)$Name
  if (!all(safe_member(mem))) stop("unsafe member name in ", zf)
  cel_mem <- mem[grepl("/celfiles/[0-9]+\\.CEL$", mem)]
  bc <- sub("\\.CEL$", "", basename(cel_mem))
  if (!setequal(bc, a$BARCODE)) {
    return(sprintf("barcode mismatch: %d CEL vs %d attribute rows", length(bc), nrow(a)))
  }
  ex <- file.path(tempdir(), paste0("cel_", key))
  dir.create(ex)
  unzip(zf, files = cel_mem, exdir = ex, junkpaths = TRUE)
  files <- file.path(ex, paste0(a$BARCODE, ".CEL"))
  stopifnot(all(file.exists(files)))
  x <- rma_r(files)
  unlink(ex, recursive = TRUE)  # session temp copy only
  stopifnot(ncol(x) == nrow(a), !anyNA(x))
  colnames(x) <- a$BARCODE

  grp <- factor(a$group)
  design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)
  ctrl_levels <- grep("^Control_", levels(grp), value = TRUE)
  treat_levels <- setdiff(levels(grp), ctrl_levels)
  cons <- character(0)
  for (tl in treat_levels) {
    tm <- sub("^[^_]+_", "", tl)
    cl <- paste0("Control_", tm)
    if (cl %in% ctrl_levels) cons <- c(cons, paste0(tl, " - ", cl))
  }
  fit <- lmFit(x, design)
  cm <- makeContrasts(contrasts = cons, levels = design)
  colnames(cm) <- sub(" - .*$", "", cons)
  fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)

  # vehicle versus vehicle at 24 h
  veh <- NULL
  c24 <- which(a$group == "Control_24h")
  if (length(c24) == 2) {
    g2 <- as.character(a$group); g2[c24[1]] <- "VehA_24h"; g2[c24[2]] <- "VehB_24h"
    g2 <- factor(g2)
    d2 <- model.matrix(~0 + g2); colnames(d2) <- levels(g2)
    fv <- lmFit(x, d2)
    cmv <- makeContrasts(contrasts = "VehA_24h - VehB_24h", levels = d2)
    fv2 <- eBayes(contrasts.fit(fv, cmv), trend = TRUE, robust = TRUE)
    veh <- list(t = fv2$t[, 1], logFC = fv2$coefficients[, 1], df = fv2$df.total,
                barcodes = a$BARCODE[c24])
  }
  ctrl_cols <- grepl("^Control_", a$group)
  res <- list(
    compound = cname, key = key, exp_id = unique(a$EXP_ID),
    n_arrays = ncol(x), groups = table(a$group),
    doses = unique(a[, c("DOSE_LEVEL", "DOSE", "time")]),
    t = fit2$t, logFC = fit2$coefficients, df = fit2$df.total,
    ctrl_mean = rowMeans(x[, ctrl_cols, drop = FALSE]), n_ctrl = sum(ctrl_cols),
    veh = veh
  )
  saveRDS(res, out_rds, compress = "xz")
  s <- sprintf("ok: %d arrays, %d contrasts, vehicle=%s", ncol(x), length(cons), !is.null(veh))
  cat(k, cname, s, "\n")
  s
}

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
out <- parallel::mclapply(seq_along(zips), process_one, mc.cores = n_cores, mc.preschedule = FALSE)
z_status <- vapply(out, function(s) if (inherits(s, "try-error") || !is.character(s)) paste("error:", trimws(as.character(s)))[1] else s, "")
st <- data.frame(compound = cmp_of_zip, status = z_status, row.names = NULL)
write.table(st, file.path(out_dir, "rma_limma_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_rma_limma.txt"))
cat("done:", sum(grepl("^ok", st$status)), "ok of", nrow(st), "\n")
