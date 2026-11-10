#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/verify_core_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/verify_core_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cp Analysis/SingleCell/scripts/diagnosis_correction/verify_results.py "$TASK/job_scripts/verify_core_${SLURM_JOB_ID}.py"
# Keep __file__ based project resolution in the canonical checker; Python
# parses its complete source before running the checks.
"$PY" Analysis/SingleCell/scripts/diagnosis_correction/verify_results.py --phase core
export SINGLECELL_BINARY_DE_OUT="$TASK/pseudobulk_de_verified_obs"
cp Analysis/SingleCell/scripts/pseudobulk_de.R "$TASK/job_scripts/binary_obs_regression_${SLURM_JOB_ID}.R"
"$RS" "$TASK/job_scripts/binary_obs_regression_${SLURM_JOB_ID}.R"
"$PY" - <<'PY'
from pathlib import Path
import numpy as np
import pandas as pd
base=Path('Analysis/SingleCell/candidates/source-diagnosis-correction-20261026')
a=sorted((base/'pseudobulk_de_v2').glob('*_de.csv'))
b=sorted((base/'pseudobulk_de_verified_obs').glob('*_de.csv'))
assert [p.name for p in a]==[p.name for p in b]
for old,new in zip(a,b):
 x=pd.read_csv(old).sort_values('gene');y=pd.read_csv(new).sort_values('gene')
 assert x.gene.tolist()==y.gene.tolist()
 assert np.allclose(x[['logFC','t_stat','pvalue','padj']],y[['logFC','t_stat','pvalue','padj']],atol=1e-10,rtol=1e-8,equal_nan=True)
(base/'GENE_OBS_REGRESSION_PASSED').touch()
print('Corrected atlas vs explicit metadata overlay: all lineage gene effects agree.')
PY
python3 scripts/manuscript/validate_resource_scope.py > "$TASK/resource_scope_check.txt"
touch "$TASK/CORE_CHECKS_COMPLETED"
