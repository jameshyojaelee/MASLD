#!/bin/bash
#SBATCH --job-name=DIG1gt67b
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.err
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
GTF=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/gencode.v49.primary_chr.gtf
sha256sum "$RUN/scripts/67b_genotype_scan_null.R" > "$RUN/script_sha256_job67b.txt"
GT=$RUN/scripts/67_gencode_v49_gene_type.tsv
test -e "$GT" || awk -F'\t' '$3=="gene"{match($9,/gene_id "[^"]+"/); id=substr($9,RSTART+9,RLENGTH-10); match($9,/gene_type "[^"]+"/); ty=substr($9,RSTART+11,RLENGTH-12); print id"\t"ty}' "$GTF" > "$GT"
wc -l "$GT"; head -2 "$GT"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$RUN/scripts/67b_genotype_scan_null.R" "$RUN/out_genotype_state" "$RUN/out_genotype_posthoc" "$GT" "$RUN/out_genotype_scan_null_67b"
