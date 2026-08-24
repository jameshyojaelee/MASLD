#!/usr/bin/env bash
#SBATCH --job-name=masld-gse105127-detail-probe
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export LC_ALL=C

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE="${ROOT}/executions/gse105127-rrbs-source-qc-21065724.staging/raw/GSM2819620__GSM2819620_CV_6610_207.MCSv3.20161129.hs37.cpg.filtered.CG.bed.gz"
gzip -dc "${SOURCE}" | awk '
  /^track / { next }
  $10 !~ /^0([.]0+)?$/ {
    seen++
    if (seen <= 20) print
    values[$10]++
  }
  END {
    for (value in values) print "detail_value", value, values[value] > "/dev/stderr"
  }
'
