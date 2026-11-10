#!/bin/bash
#SBATCH --account=innovation --qos=innovation --job-name=code-review
#SBATCH --partition=cpu --prefer=v6 --cpus-per-task=1 --mem=16G --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/logs/ssx_b_gtf_%j.log
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z
python3 -I scripts/ssx_b_intron_gtf.py /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION/external_reference/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz release_genes.txt /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/out_b_gtf-$(date -u +%Y%m%dT%H%M%SZ)
echo EXIT_OK
