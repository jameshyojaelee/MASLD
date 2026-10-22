#!/bin/bash
#SBATCH --job-name=9815
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1/slurm-9815-%j.out
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
F=GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv
O=Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1/positive_control_loci_rows.csv
head -1 $F > $O
grep -E ",(PNPLA3|SAMM50|PARVB|TM6SF2|SUGP1|MAU2|GATAD2A|HSD17B13|HSD17B11|GPAM|GCKR|MBOAT7|TMC4),ENSG" $F >> $O
wc -l $O
