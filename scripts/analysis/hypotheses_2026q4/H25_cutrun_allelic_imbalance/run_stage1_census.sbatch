#!/bin/bash
#SBATCH --job-name=h25_census
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
# H25 stage 1 census and stop rule. Usage: sbatch run_stage1_census.sbatch <hyp_dir> <data_dir>
set -euo pipefail
HD=$1
DATA=$2
S=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H25_cutrun_allelic_imbalance
ENV=/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl
OUT=$HD/stage1_census_$(date -u +%Y%m%dT%H%M%SZ)
$ENV/bin/python -I "$S/s02_stage1_census.py" --runs "$HD/inputs/runs.tsv" --targets "$HD/stage1_targets/targets.tsv" \
  --stage1-dir "$DATA/stage1" --out "$OUT"
$ENV/bin/python -I -m pip freeze > "$OUT/pip_freeze.txt"
module load Bowtie2/2.5.4-linux-x86_64 SAMtools/1.18-GCC-12.3.0 SRA-Toolkit/3.2.1 Java/17
set +o pipefail  # `| head` closes pipes early in the version block
{ bowtie2 --version | head -1; samtools --version | head -1; fasterq-dump --version 2>&1 | grep -m1 fasterq || true;
  /gpfs/commons/home/jameslee/bbmap/bbversion.sh; $ENV/bin/python --version; } > "$OUT/tool_versions.txt" 2>&1
# aggregate per-run read-processing statistics
{ echo -e "srr\tbait_pairs\twasp_tested\twasp_pass\tdedup_in\tdedup_out";
  for d in "$DATA"/stage1/SRR*; do s=$(basename "$d");
    bp=$(head -1 "$d/bowtie2.log" | awk '{print $1}');
    t=$($ENV/bin/python -I -c "import json,sys;j=json.load(open(sys.argv[1]));print(j['tested'],j['tested_pass'],sep='\t')" "$d/wasp.check_stats.json");
    dd=$($ENV/bin/python -I -c "import json,sys;j=json.load(open(sys.argv[1]));print(j['fragments_in'],j['fragments_after_dedup'],sep='\t')" "$d/counts.tsv.dedup.json");
    echo -e "$s\t$bp\t$t\t$dd"; done; } > "$OUT/read_processing_stats.tsv"
echo "$OUT"
