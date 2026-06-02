#!/bin/bash
# Build unified driver TSV across 5 mega cohorts (Suppli/Govaere/Hoang/Bril/Chen).
# Two schemas exist (5-col vs 10-col); use header lookup to pull sample_id/layout/fastq_r{1,2}.
# Output columns: cohort sample_id layout fastq_r1 fastq_r2
set -euo pipefail

WT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom
OUT="$WT/RNA-seq/scripts/kallisto/sample_driver.tsv"

echo -e "cohort\tsample_id\tlayout\tfastq_r1\tfastq_r2" > "$OUT"
for COHORT in GSE126848 GSE135251 GSE130970 GSE162694 GSE213621; do
  SS="$BASE/$COHORT/metadata/samples_with_paths.tsv"
  if [[ ! -f "$SS" ]]; then echo "[WARN] missing $SS" >&2; continue; fi
  awk -F'\t' -v c="$COHORT" '
    NR==1 {
      for (i=1; i<=NF; i++) idx[$i] = i;
      next;
    }
    {
      sid = $(idx["sample_id"]);
      lay = $(idx["layout"]);
      r1  = $(idx["fastq_r1"]);
      r2  = $(idx["fastq_r2"]);
      print c"\t"sid"\t"lay"\t"r1"\t"r2;
    }' "$SS" >> "$OUT"
done

N=$(($(wc -l < "$OUT") - 1))
echo "Wrote $OUT with $N samples"
cut -f1 "$OUT" | tail -n +2 | sort | uniq -c
echo "---layouts---"
cut -f1,3 "$OUT" | tail -n +2 | sort | uniq -c
