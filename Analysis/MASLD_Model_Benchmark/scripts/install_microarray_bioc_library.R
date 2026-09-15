#!/usr/bin/env Rscript

# Build the portable Bioconductor 3.22 library overlay that activates the two
# blocked microarray transfer TaskSpecs (GSE49541 / GPL570, GSE83452 / GPL16686).
#
# The library tree is relocatable: it holds installed R packages only, so the
# staging directory can be frozen and moved. R itself comes from the
# R/4.5.2-mba module, which is a fixed system path and is never copied.
#
# Every required package is asserted at an exact version, and every source
# tarball this script downloads is checked against the exact Bioconductor 3.22
# repository MD5 before the manifest is written.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: install_microarray_bioc_library.R LIBRARY_DIR DOWNLOAD_DIR MANIFEST.tsv")
}
library_dir <- normalizePath(args[[1L]], mustWork = FALSE)
download_dir <- normalizePath(args[[2L]], mustWork = FALSE)
manifest_path <- args[[3L]]

BIOC_RELEASE <- "3.22"
R_VERSION <- "4.5.2"

# package -> exact required version.
REQUIRED <- c(
  affy = "1.88.0",
  affyio = "1.80.0",
  affxparser = "1.82.0",
  oligo = "1.74.0",
  frma = "1.62.0",
  SCAN.UPC = "2.52.0",
  AnnotationDbi = "1.72.0",
  org.Hs.eg.db = "3.22.0",
  hgu133plus2cdf = "2.18.0",
  hgu133plus2.db = "3.13.0",
  hgu133plus2frmavecs = "1.5.0",
  hugene20sttranscriptcluster.db = "8.8.0",
  pd.hugene.2.0.st = "3.14.1",
  pd.hg.u133.plus.2 = "3.12.0"
)

# Exact Bioconductor 3.22 repository MD5 for every source tarball this script is
# expected to fetch. AnnotationDbi 1.72.0 already ships with the R module, so no
# tarball is downloaded for it and only its version is asserted.
PINNED_MD5 <- c(
  "affy_1.88.0.tar.gz" = "85a03c51a040597094ab24e8067e5830",
  "affyio_1.80.0.tar.gz" = "0c5147bb9f0b683c9296156ee07a77d0",
  "affxparser_1.82.0.tar.gz" = "9f748ecae1c1ed5c10113e2dbced634d",
  "oligo_1.74.0.tar.gz" = "58b1b34f18141fb209591bc14e633a9e",
  "frma_1.62.0.tar.gz" = "cb2d731c4fc1edac3957c43ec4c55fcc",
  "SCAN.UPC_2.52.0.tar.gz" = "1365db0f1f73e73b39a85ccecf87a189",
  "org.Hs.eg.db_3.22.0.tar.gz" = "e80cac6ec018a95aea4f7530350e80a2",
  "hgu133plus2cdf_2.18.0.tar.gz" = "284fef2f0b777d7b53451538ddd53de3",
  "hgu133plus2.db_3.13.0.tar.gz" = "459fcf4880a9eaa25b373c5635fede3d",
  "hgu133plus2frmavecs_1.5.0.tar.gz" = "a4781cbcccc1ee17dfd16259f1c7bebc",
  "hugene20sttranscriptcluster.db_8.8.0.tar.gz" = "0b929a3a959662e8a7265f58b81b4e35",
  "pd.hugene.2.0.st_3.14.1.tar.gz" = "e484209aa0c2a839c3445d91c1a799ce",
  "pd.hg.u133.plus.2_3.12.0.tar.gz" = "8a87aa63c04e84266962bdde5226c06c"
)

if (!identical(as.character(getRversion()), R_VERSION)) {
  stop(sprintf("expected R %s, found %s", R_VERSION, getRversion()))
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop("BiocManager is absent from the R module")
}
if (!identical(as.character(BiocManager::version()), BIOC_RELEASE)) {
  stop(sprintf(
    "expected Bioconductor %s, found %s", BIOC_RELEASE, BiocManager::version()
  ))
}

dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))
system_library <- .libPaths()[-1L]

cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (is.na(cores) || cores < 1L) cores <- 1L
options(
  repos = BiocManager::repositories(version = BIOC_RELEASE),
  Ncpus = cores,
  timeout = 3600,
  warn = 1
)

# org.Hs.eg.db is installed first: both platform .db packages depend on it.
for (batch in list(
  c("affyio", "affxparser", "preprocessCore", "affy"),
  c("oligoClasses", "oligo"),
  c("frma", "SCAN.UPC"),
  c("org.Hs.eg.db"),
  c("hgu133plus2.db", "hugene20sttranscriptcluster.db"),
  c("hgu133plus2cdf", "hgu133plus2frmavecs"),
  c("pd.hugene.2.0.st", "pd.hg.u133.plus.2")
)) {
  BiocManager::install(
    batch,
    lib = library_dir,
    version = BIOC_RELEASE,
    update = FALSE,
    ask = FALSE,
    checkBuilt = FALSE,
    destdir = download_dir,
    type = "source",
    INSTALL_opts = "--no-multiarch",
    # preprocessCore built with threading is a known source of hangs inside
    # frma's quantile step on shared clusters; frma is the GPL570 method here.
    configure.args = c(preprocessCore = "--disable-threading")
  )
}

## ------------------------------------------------------------ exact versions
installed <- installed.packages()
missing <- setdiff(names(REQUIRED), rownames(installed))
if (length(missing)) {
  stop(sprintf("required package absent after install: %s", paste(missing, collapse = ", ")))
}
observed <- vapply(
  names(REQUIRED), function(p) unname(installed[p, "Version"]), character(1)
)
differing <- names(REQUIRED)[observed != REQUIRED]
if (length(differing)) {
  stop(sprintf(
    "required package version differs: %s",
    paste(sprintf("%s %s != %s", differing, observed[differing], REQUIRED[differing]),
          collapse = "; ")
  ))
}
# Every package the platform pipeline calls must resolve from the portable
# overlay or the fixed R module, never from a user library.
loaded_from <- vapply(
  names(REQUIRED), function(p) unname(installed[p, "LibPath"]), character(1)
)
foreign <- setdiff(loaded_from, c(library_dir, system_library))
if (length(foreign)) {
  stop(sprintf("required package resolved outside the runtime: %s",
               paste(foreign, collapse = "; ")))
}

## ------------------------------------------------------------ exact tarballs
tarballs <- list.files(download_dir, pattern = "\\.tar\\.gz$", full.names = TRUE)
names(tarballs) <- basename(tarballs)
pinned_missing <- setdiff(names(PINNED_MD5), names(tarballs))
if (length(pinned_missing)) {
  stop(sprintf("pinned source tarball was not downloaded: %s",
               paste(pinned_missing, collapse = ", ")))
}
observed_md5 <- unname(tools::md5sum(tarballs[names(PINNED_MD5)]))
if (!identical(observed_md5, unname(PINNED_MD5))) {
  stop("pinned Bioconductor source tarball MD5 differs")
}

## ----------------------------------------------------------------- manifest
overlay <- installed[installed[, "LibPath"] == library_dir, , drop = FALSE]
overlay <- overlay[order(rownames(overlay)), , drop = FALSE]
digest <- function(path) {
  if (is.na(path) || !nzchar(path)) return("not_downloaded")
  unname(tools::sha256sum(path))
}
records <- data.frame(
  package = rownames(overlay),
  version = unname(overlay[, "Version"]),
  library = "portable_overlay",
  stringsAsFactors = FALSE
)
records$source_tarball <- vapply(
  seq_len(nrow(records)),
  function(index) {
    name <- paste0(records$package[[index]], "_", records$version[[index]], ".tar.gz")
    if (name %in% names(tarballs)) name else "not_downloaded"
  },
  character(1)
)
records$source_sha256 <- vapply(
  records$source_tarball,
  function(name) if (identical(name, "not_downloaded")) "not_downloaded" else digest(tarballs[[name]]),
  character(1)
)
records$source_md5 <- vapply(
  records$source_tarball,
  function(name) {
    if (identical(name, "not_downloaded")) "not_downloaded" else unname(tools::md5sum(tarballs[[name]]))
  },
  character(1)
)
records$pinned <- ifelse(records$source_tarball %in% names(PINNED_MD5), "exact_md5_pinned", "resolved_dependency")

system_packages <- installed[installed[, "LibPath"] != library_dir, , drop = FALSE]
system_packages <- system_packages[order(rownames(system_packages)), , drop = FALSE]
system_records <- data.frame(
  package = rownames(system_packages),
  version = unname(system_packages[, "Version"]),
  library = "r_module_4.5.2_mba",
  source_tarball = "not_downloaded",
  source_sha256 = "not_downloaded",
  source_md5 = "not_downloaded",
  pinned = "module_provided",
  stringsAsFactors = FALSE
)
write.table(
  rbind(records, system_records), manifest_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(dirname(manifest_path), "install.sessionInfo.txt")
)
cat("LIBRARY_INSTALL\tpass\t", nrow(records), "\toverlay_packages\n", sep = "")
