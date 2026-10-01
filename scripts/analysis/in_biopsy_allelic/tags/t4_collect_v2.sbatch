#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=12:00:00
# T4 for tag set v2: collect the v1 pileup at v2 SNVs and the v2-new pileup (collector skips
# indel records), then merge into one long table.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%j.out t4_collect_v2.sbatch <v2_exec_dir> <v2new_restricted_dir> <v1_restricted_dir> <merged_restricted_dir>
set -euo pipefail
py="env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"
t=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/in_biopsy_allelic/tags
v2=$1; nd=$2; v1d=$3; md=$4
mkdir -p "$md"; chmod 750 "$md"
$py $t/t4_collect.py --ad-dir "$v1d/ad" --tags "$v2/tags_v2_all_snvs.tsv" --out "$md/tag_allele_depths_v1pileup.tsv.gz"
$py $t/t4_collect.py --ad-dir "$nd/ad" --tags "$v2/tags_v2_new_snvs.tsv" --out "$md/tag_allele_depths_v2new.tsv.gz"
$py $t/t4_merge_v2.py --tags-dir "$v2/tags" --v1-ad "$md/tag_allele_depths_v1pileup.tsv.gz" \
  --v2new-ad "$md/tag_allele_depths_v2new.tsv.gz" \
  --out "$md/tag_allele_depths.tsv.gz" --depth-summary "$v2/tag_depth_summary_by_cohort.tsv"
sha256sum "$md/tag_allele_depths.tsv.gz" > "$md/tag_allele_depths.tsv.gz.sha256"
