#!/bin/bash
#SBATCH -J 8120
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH --time 00:30:00
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -p cpu
#SBATCH -o /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG11/job_%j.out
G=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG11
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/public_r4/neale_30880
export TMPDIR=$G/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
B=https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs
( for s in female male; do wget -q -O $D/30880_irnt.gwas.imputed_v3.$s.varorder.tsv.bgz $B/30880_irnt.gwas.imputed_v3.$s.varorder.tsv.bgz & done; wait; ls -l $D > $G/download_sizes.txt; $PY $G/pc.py ) &
$PY $G/count.py
wait
$PY -m pip freeze > $G/pip_freeze.txt 2>/dev/null
