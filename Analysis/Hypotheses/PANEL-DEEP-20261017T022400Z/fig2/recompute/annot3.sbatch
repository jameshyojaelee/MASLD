#!/bin/bash
#SBATCH --job-name=9823
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig2/recompute/annot3_%j.out
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig2/recompute
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
$R annotate_leads2.R out_ABCEF/C_union_leads_maxboth.tsv out_ABCEF/C_union_leads_maxboth_annot2.tsv
$R annotate_leads.R out_ABCEF/C_union_leads_maxboth.tsv out_ABCEF/C_union_leads_maxboth_annot.tsv
$R annotate_leads2.R out_ABCEF/B_shown_for_annot.tsv out_ABCEF/B_shown_annot2.tsv
echo done
