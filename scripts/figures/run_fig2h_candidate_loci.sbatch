#!/usr/bin/env bash
#SBATCH --job-name=fig2h_candidates
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/fig2h_candidates_%j.out
#SBATCH --error=logs/fig2h_candidates_%j.err

set -euo pipefail

project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
rscript_bin="/gpfs/commons/home/jameslee/micromamba/envs/plotgardener/bin/Rscript"
stamp="$(date +%Y%m%d_%H%M%S)"
candidate_dir="${project_root}/figures/candidates/fig2h_loci_${stamp}_${SLURM_JOB_ID}"
renderer="${project_root}/scripts/figures/figS09_locus_zoom_pg.R"

mkdir -p "${project_root}/logs" "${candidate_dir}"
cd "${project_root}"
export LZ_COMPACT_2H=1

LZ_WIN_START=121150000 LZ_WIN_END=121455000 \
LZ_EUR_EVIDENCE_LABEL=$'evaluated\nPP.H4 <= 0.5' \
LZ_EAS_EVIDENCE_LABEL="multi-signal COLOC" \
  "${rscript_bin}" "${renderer}" locus_GGT_chr12_121380544 GGT ACADS \
  "${candidate_dir}/Fig2H_candidate_ACADS_GGT_EAS_raw.pdf" \
  > "${candidate_dir}/ACADS_render.log" 2>&1

LZ_EUR_EVIDENCE_LABEL=$'evaluated\nPP.H4 <= 0.5' \
LZ_EAS_EVIDENCE_LABEL="multi-signal COLOC" \
"${rscript_bin}" "${renderer}" locus_AST_chr2_233520254 AST EFHD1 \
  "${candidate_dir}/Fig2H_candidate_EFHD1_AST_EAS_raw.pdf" \
  > "${candidate_dir}/EFHD1_render.log" 2>&1

for stem in Fig2H_candidate_ACADS_GGT_EAS Fig2H_candidate_EFHD1_AST_EAS; do
  gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dFirstPage=1 -dLastPage=1 \
    -sOutputFile="${candidate_dir}/${stem}.pdf" "${candidate_dir}/${stem}_raw.pdf"
done

for pdf in "${candidate_dir}"/Fig2H_candidate_ACADS_GGT_EAS.pdf \
           "${candidate_dir}"/Fig2H_candidate_EFHD1_AST_EAS.pdf; do
  pdfinfo "${pdf}" | awk '/Pages:|Page size:/' >> "${candidate_dir}/pdf_validation.txt"
  pdftotext "${pdf}" - >> "${candidate_dir}/artwork_text.txt"
done

if [[ "$(grep -c 'Pages:[[:space:]]*1' "${candidate_dir}/pdf_validation.txt")" -ne 2 ]]; then
  echo "Candidate PDF is not one page" >&2
  exit 1
fi

if grep -Eqi 'ABF-fallback|exploratory|causal gene|ancestry-specific' \
  "${candidate_dir}/artwork_text.txt"; then
  echo "Prohibited artwork label found" >&2
  exit 1
fi

{
  printf 'candidate\tpriority\ttrait\tgwas_ancestry\teqtl_ancestry\twindow_hg19\tstatus\n'
  printf 'ACADS\t1\tGGT\tEAS\tEUR\tchr12:121150000-121455000\tprovisional\n'
  printf 'EFHD1\t2\tAST\tEAS\tEUR\tchr2:233420254-233620254\tprovisional\n'
} > "${candidate_dir}/Fig2H_candidate_manifest.tsv"

printf '%s\n' "${candidate_dir}"
