#!/usr/bin/env bash
#SBATCH --job-name=MPRAnalyze
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/MPRAnalyze_env_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/MPRAnalyze_env_%j.err
set -euo pipefail
ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
export MASLD_PROJECT_ROOT=${ROOT}
ENV=${ROOT}/.mamba/seqfunc_mpra_v2
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/mpra_benchmark/v2
mkdir -p "${OUT}"
rm -f "${OUT}/environment_qc.json"
# jq must exist in the prefix before *any* Bioconductor data-package transaction:
# conda's GenomeInfoDbData post-link invokes yq, which shells out to jq.
if [[ -d ${ENV} ]]; then
  micromamba install -y -p "${ENV}" -c conda-forge jq
  micromamba env update -y -p "${ENV}" -f "${ROOT}/envs/seqfunc_mpra_v2.yml"
else
  micromamba create -y -p "${ENV}" -c conda-forge jq
  micromamba env update -y -p "${ENV}" -f "${ROOT}/envs/seqfunc_mpra_v2.yml"
fi
# The first build can leave GenomeInfoDbData registered in conda while its
# post-link payload is absent if jq was unavailable. Force its post-link again
# after jq is present, then validate actual R namespaces rather than package
# metadata alone.
# Keep these as two transactions. If they are force-reinstalled together,
# micromamba unlinks jq and can run the GenomeInfoDbData post-link before jq is
# relinked, recreating the exact partial-install state this repair addresses.
micromamba install -y --force-reinstall -p "${ENV}" -c conda-forge -c bioconda \
  bioconductor-genomeinfodbdata
micromamba run -p "${ENV}" Rscript -e 'pkgs <- c("GenomeInfoDbData","MPRAnalyze","DESeq2","BiocParallel","data.table","jsonlite"); ok <- vapply(pkgs, requireNamespace, logical(1), quietly=TRUE); print(ok); stopifnot(all(ok)); root <- Sys.getenv("MASLD_PROJECT_ROOT"); p <- file.path(root,"GWAS/finemapping/results/seqfunc/mpra_benchmark/v2/counts/HepG2.barcode_counts.tsv.gz"); z <- data.table::fread(cmd=paste("gzip -dc",shQuote(p)),nrows=2); stopifnot(nrow(z)==2L); sessionInfo()'
export OUT
micromamba run -p "${ENV}" Rscript -e 'jsonlite::write_json(list(environment_complete=TRUE,MPRAnalyze=as.character(packageVersion("MPRAnalyze")),DESeq2=as.character(packageVersion("DESeq2")),GenomeInfoDbData=as.character(packageVersion("GenomeInfoDbData"))), file.path(Sys.getenv("OUT"),"environment_qc.json"), auto_unbox=TRUE, pretty=TRUE)'
