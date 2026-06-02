options(repos = c(CRAN = "https://cloud.r-project.org"))

# Force installs into the active conda env library if provided
lib_sites <- c(Sys.getenv("R_LIBS_SITE"), Sys.getenv("R_LIBS"))
lib_sites <- lib_sites[lib_sites != ""]
if (length(lib_sites) > 0) {
  .libPaths(unique(lib_sites))
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Core Bioconductor packages (safety if not provided by conda)
BiocManager::install(
  c(
    "zellkonverter",
    "SingleCellExperiment",
    "SummarizedExperiment",
    "BiocParallel",
    "scater",
    "scran"
  ),
  ask = FALSE,
  update = FALSE
)

# MuSiC
BiocManager::install("TOAST", ask = FALSE, update = FALSE)
install.packages(c("MatrixModels", "quantreg", "MCMCpack", "nnls"), quiet = TRUE)
remotes::install_github("xuranw/MuSiC")
