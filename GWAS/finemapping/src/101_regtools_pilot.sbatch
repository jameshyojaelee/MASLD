#!/usr/bin/env bash
#SBATCH --job-name=regtools
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --array=0-19
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/regtools_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/regtools_%A_%a.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_leafcutter
MAN=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/pilot/pilot_manifest.tsv
mkdir -p "$(dirname "${MAN}")/junctions"
line=$(awk -F '\t' -v i="${SLURM_ARRAY_TASK_ID}" 'NR>1 && $1==i {print; exit}' "${MAN}")
line=${line//$'\r'/}
test -n "${line}"
sample=$(cut -f2 <<<"${line}"); bam=$(cut -f5 <<<"${line}"); sj=$(cut -f6 <<<"${line}"); out=$(cut -f7 <<<"${line}")
test -s "${bam}"; test -s "${sj}"
# All canonical STAR libraries are reverse-stranded; regtools 1.0 calls this RF.
local_bam="${out}.bam"
ln -sfn "${bam}" "${local_bam}"
"${ENV}/bin/samtools" index -@ "${SLURM_CPUS_PER_TASK}" -o "${local_bam}.bai.tmp" "${local_bam}"
mv "${local_bam}.bai.tmp" "${local_bam}.bai"
"${ENV}/bin/regtools" junctions extract -a 8 -m 50 -M 500000 -s RF "${local_bam}" -o "${out}.tmp"
test -s "${out}.tmp"; mv "${out}.tmp" "${out}"
awk -v s="${sample}" 'BEGIN{OFS="\t"} {n++; reads+=$5} END{print s,n,reads}' "${out}" > "${out}.qc.tsv"
