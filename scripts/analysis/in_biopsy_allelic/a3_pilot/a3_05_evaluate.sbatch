#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
# A3 pilot gate G1 evaluation. Usage: sbatch --job-name=<n> --output=<exec>/logs/%j.out a3_05_evaluate.sbatch <exec_dir>
set -euo pipefail
ext=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
export PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/.local/bin/micromamba run -n allelic_qtl python \
  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/in_biopsy_allelic/a3_pilot/a3_05_evaluate.py \
  --exec-dir "$1" --gtex-egenes "$ext/allelic_refs/gtex_v8/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz" \
  --geuvadis-bam "$ext/geuvadis/bam"
