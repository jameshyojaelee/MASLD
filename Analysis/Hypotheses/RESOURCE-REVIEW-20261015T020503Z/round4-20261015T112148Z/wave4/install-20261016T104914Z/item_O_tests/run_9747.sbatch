#!/bin/bash
#SBATCH --job-name=9747
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/wave4/install-20261016T104914Z/item_O_tests/9747_%j.out
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/src
/gpfs/commons/home/jameslee/micromamba/envs/borzoi/bin/python -m unittest discover -s tests/contract -t . -p test_registry_files.py 2>&1 | grep -E "^(FAIL|ERROR|Ran|OK|FAILED)|Error:" | head -60
