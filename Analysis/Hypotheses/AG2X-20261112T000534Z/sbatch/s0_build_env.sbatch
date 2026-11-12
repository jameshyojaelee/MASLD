#!/bin/bash
#SBATCH --job-name=ag2x_env
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/%x_%j.err
set -euo pipefail
R=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
O=$R/Analysis/Hypotheses/AG2X-20261112T000534Z
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba
if [ -d $MAMBA_ROOT_PREFIX/envs/alphagenome2 ]; then echo "env alphagenome2 already exists; not rebuilding"; exit 1; fi
$MM create -y -f $R/envs/alphagenome2.yml
PY=$MAMBA_ROOT_PREFIX/envs/alphagenome2/bin/python
$PY -c "import alphagenome, importlib.metadata as m; from alphagenome.models.v2 import dna_client, variant_scorers, model_condition; from alphagenome.models.v1 import dna_client as d1; print('alphagenome', m.version('alphagenome'), 'v2 import ok')"
$PY -m pip freeze > $O/stage0/pip_freeze_alphagenome2.txt
$PY --version > $O/stage0/python_version.txt
echo done
