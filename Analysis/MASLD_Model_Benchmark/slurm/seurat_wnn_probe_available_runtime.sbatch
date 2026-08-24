#!/usr/bin/env bash
#SBATCH --job-name=masld-seurat-wnn-runtime
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
export R_LIBS_USER=/nonexistent/masld_bench_r_user_library

module load R/4.5.2-mba
Rscript --vanilla - <<'RS'
required <- c(
  "Seurat", "SeuratObject", "Matrix", "Rcpp", "RcppAnnoy", "RcppHNSW",
  "igraph", "sctransform"
)
cat("r_version\t", R.version.string, "\n", sep = "")
cat("platform\t", R.version$platform, "\n", sep = "")
for (package in required) {
  available <- requireNamespace(package, quietly = TRUE)
  version <- if (available) as.character(packageVersion(package)) else "NOT_INSTALLED"
  cat(package, "\t", available, "\t", version, "\n", sep = "")
}
cat("project_data_read\tFALSE\n")
cat("outcomes_read\tFALSE\n")
cat("model_fit_executed\tFALSE\n")
sessionInfo()
RS
