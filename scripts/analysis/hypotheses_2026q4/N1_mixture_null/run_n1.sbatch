#!/bin/bash
#SBATCH --job-name=N1_mixnull
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# N1 step one. Usage: sbatch --output=<out>/slurm-%j.log run_n1.sbatch <out_dir>
set -euo pipefail
out=$1
repo=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
py=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
ex=$repo/Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
test -s "$out/PRESPEC.md"
test ! -e "$out/summary.json"
$py -m pip freeze > "$out/pip_freeze.txt"
cp "$repo/scripts/analysis/hypotheses_2026q4/N1_mixture_null/n1_mixture_null.py" "$out/n1_mixture_null.as_run.py"
$py "$repo/scripts/analysis/hypotheses_2026q4/N1_mixture_null/n1_mixture_null.py" \
  --vcf-dir "$repo/RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/a1-full-20260924T025858Z/vcf" \
  --targets "$ex/targets/targets_annotated.tsv" \
  --identity "$ex/identity/sample_to_individual.tsv" \
  --kinship "$ex/identity/kinship_related_pairs.tsv" \
  --results-root "$repo/RNA-seq/Human/Patient_Cohorts/results" \
  --out "$out" --seed 20261007 --threads 8
echo "N1 done"
