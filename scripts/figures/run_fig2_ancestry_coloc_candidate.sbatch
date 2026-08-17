#!/usr/bin/env bash
#SBATCH --job-name=fig2_ancestry_coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fig2_ancestry_coloc_%j.out
#SBATCH --error=logs/fig2_ancestry_coloc_%j.err

set -euo pipefail

project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
python_bin="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"
rscript_bin="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq_r44/bin/Rscript"
stamp="$(date +%Y%m%d_%H%M%S)"
candidate_dir="${project_root}/figures/candidates/fig2_ancestry_coloc_${stamp}_${SLURM_JOB_ID}"

mkdir -p "${project_root}/logs" "${candidate_dir}"
export MASLD_PROJECT_ROOT="${project_root}"
export FIG2_CANDIDATE_DIR="${candidate_dir}"
export FIG2_PANEL="cascade"
export FIG2_RELEASE_STATE="provisional"
export PYTHONNOUSERSITE=1

cd "${project_root}"

{
  printf 'candidate_dir\t%s\n' "${candidate_dir}"
  printf 'release_state\t%s\n' "${FIG2_RELEASE_STATE}"
  printf 'coloc_input\t%s\n' "${FIG2_COLOC_INPUT:-GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv}"
  printf 'noneur_audit\t%s\n' "${FIG2_NONEUR_AUDIT:-RNA-seq/results/coloc_variant_classes/noneur_gws_audit.csv}"
  sha256sum "${FIG2_COLOC_INPUT:-GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv}"
  sha256sum "${FIG2_NONEUR_AUDIT:-RNA-seq/results/coloc_variant_classes/noneur_gws_audit.csv}"
} > "${candidate_dir}/input_manifest.tsv"

"${python_bin}" scripts/figures/select_fig2h_noneur_locus.py \
  > "${candidate_dir}/Fig2H_audit.log" 2>&1
"${python_bin}" scripts/figures/fig1_gwas_creative_options.py \
  > "${candidate_dir}/Fig2A_render.log" 2>&1
"${rscript_bin}" scripts/figures/fig2_resolution_plane.R \
  > "${candidate_dir}/Fig2B_render.log" 2>&1
GATE=gws "${python_bin}" scripts/figures/fig2_ancestry_unique_coloc_gated.py \
  > "${candidate_dir}/Fig2E_render.log" 2>&1
"${rscript_bin}" scripts/figures/fig2_tri_ancestry_coloc.R \
  > "${candidate_dir}/Fig2F_render.log" 2>&1

"${python_bin}" scripts/figures/validate_fig2_ancestry_coloc_candidate.py "${candidate_dir}" \
  > "${candidate_dir}/validation_stdout.log" 2>&1
python3 scripts/manuscript/validate_resource_scope.py \
  > "${candidate_dir}/resource_scope_validation.log" 2>&1

printf '%s\n' "${candidate_dir}"
