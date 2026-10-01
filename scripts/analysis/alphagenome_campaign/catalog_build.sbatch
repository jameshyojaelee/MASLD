#!/bin/bash
#SBATCH --job-name=masld_catalog_specimen
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_catalog_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_catalog_%j.err
set -euo pipefail
campaign_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
campaign_src=$campaign_root/scripts/analysis/alphagenome_campaign
campaign_out=$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/catalog-$SLURM_JOB_ID
mkdir -p "$campaign_out"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
python3 "$campaign_src/catalog_checks.py" > "$campaign_out/checks.txt" 2>&1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python "$campaign_src/review_prediction_checks.py" --out "$campaign_out/independent_prediction_checks"
python3 "$campaign_src/catalog_build.py" --out "$campaign_out"
python3 "$campaign_src/catalog_examples.py" --db "$campaign_out/regulatory.sqlite" --out "$campaign_out/examples"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python "$campaign_src/review_axes_catalog.py" --catalog "$campaign_out" --b2 "$campaign_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/review/b2-21771990" --out "$campaign_out/independent_axes_review"
python3 "$campaign_src/catalog_query.py" --db "$campaign_out/regulatory.sqlite" query variant GRCh38:chr15:43974532:G:A > "$campaign_out/example_variant.json"
python3 "$campaign_src/catalog_query.py" --db "$campaign_out/regulatory.sqlite" query gene ENSG00000223572 > "$campaign_out/example_gene.json"
python3 --version > "$campaign_out/environment.txt"
