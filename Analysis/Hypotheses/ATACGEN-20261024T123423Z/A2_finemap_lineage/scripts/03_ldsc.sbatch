#!/bin/bash
#SBATCH --job-name=A2_ldsc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage/logs/ldsc_%j.out
# Ancestry-matched LD scores (1 Mb window) for every variant in the 1000G panels used by
# fine-mapping. Used only as a matching covariate for the Q1/Q2 matched expectation.
set -euo pipefail
A=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage
REF=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/data/ld_ref
PY=/gpfs/commons/home/jameslee/micromamba/envs/ldsc/bin/python
LDSC=/gpfs/commons/home/jameslee/micromamba/envs/ldsc/bin/ldsc.py
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p $A/work/ldsc
# largest jobs first for load balance
for pop in afr amr sas eur eas; do for c in $(seq 1 22); do
  n=$(wc -l < $REF/1kg_$pop/chr${c}_$pop.bim); echo "$n $pop $c"; done; done | sort -k1,1nr > $A/work/ldsc/tasks.txt
run_one(){ pop=$2; c=$3; mkdir -p $A/work/ldsc/$pop
  [[ -s $A/work/ldsc/$pop/chr$c.l2.ldscore.gz ]] && return 0
  $PY $LDSC --bfile $REF/1kg_$pop/chr${c}_$pop --l2 --ld-wind-kb 1000 --yes-really \
     --out $A/work/ldsc/$pop/chr$c > /dev/null 2>&1 || echo "FAIL $pop $c"; }
export -f run_one; export A REF PY LDSC
xargs -P 4 -L 1 bash -c 'run_one $0 $1 $2' < $A/work/ldsc/tasks.txt
ls $A/work/ldsc/*/chr*.l2.ldscore.gz | wc -l
