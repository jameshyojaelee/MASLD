#!/bin/bash
#SBATCH --job-name=masld_local_inference
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_local_inference_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_local_inference_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false HF_HUB_OFFLINE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
campaign_src=$campaign_root/scripts/analysis/alphagenome_campaign
campaign_base=$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
campaign_python=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
campaign_out=$campaign_base/inference-$SLURM_JOB_ID
mkdir -p "$campaign_out"
"$campaign_python" "$campaign_src/catalog_checks.py" > "$campaign_out/checks.txt" 2>&1
"$campaign_python" "$campaign_src/catalog_examples.py" --db "$campaign_base/catalog-21772728/regulatory.sqlite" --out "$campaign_out/examples"
"$campaign_python" - "$campaign_src" "$campaign_out" <<'PY'
import json,pathlib,subprocess,sys
src,out=map(pathlib.Path,sys.argv[1:])
request=json.loads((out/'examples/synthetic_haplotype_request.json').read_text())
for name,variants in [('variant',request['variants'][:1]),('haplotype',request['variants'])]:
    subprocess.run([sys.executable,str(src/'catalog_infer.py'),*variants,
        '--phase',request['phase'],'--phase-source',request['phase_source'],
        '--length','2048','--check-repeat','--out',str(out/(name+'.json'))],check=True)
PY
