#!/usr/bin/env bash
#SBATCH --job-name=repbench-embed-geneformer
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Geneformer V2-104M embeddings + random-init twin, 557 pool samples + 109 GSE268273.
# L40S is PINNED: torch 2.5.1+cu124 predates b6k sm_120 and cannot target it.
set -euo pipefail
umask 027
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONUNBUFFERED=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
PY=executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
D=executions/liver-repbench-embed-geneformer
"$PY" -s -E "$D/extract.py"
"$PY" -s -E "$D/score_neural.py"
