#!/bin/bash
#SBATCH --job-name=code-review
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu,io
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/diag_pnpla3_dedup-20261028T004839Z/%j.out
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/diag_pnpla3_dedup-20261028T004839Z
cat samples.txt | xargs -P 8 -I{} bash one.sh {} /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/diag_pnpla3_dedup-20261028T004839Z > pileup_long.tsv
echo EXIT_OK
