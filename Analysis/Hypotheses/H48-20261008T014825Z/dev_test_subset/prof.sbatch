#!/bin/bash
#SBATCH --job-name=H48_prof
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H48-20261008T014825Z/dev_test_subset/prof_%j.log
/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_gears/bin/python -I - <<'P'
import numpy as np, time
from scipy.stats import rankdata
rng=np.random.default_rng(0); M=rng.normal(size=(20000,400))
t=time.time(); rankdata(M,axis=1); print("rank 20000x400", time.time()-t)
mem=np.arange(400); t=time.time()
for r in range(20000): rng.choice(mem,size=20,replace=False)
print("20000 choices", time.time()-t)
P
