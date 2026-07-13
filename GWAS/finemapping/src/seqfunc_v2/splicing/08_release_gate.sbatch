#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/release_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/release_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2
test -s "${OUT}/primary_events_annotated.tsv.gz"
test -s "${OUT}/event_qc.json"
test -f "${OUT}/permutations/PASS"
printf 'PASS\t%s\n' "$(date -Is)" > "${OUT}/RELEASE_GATE.tsv"
