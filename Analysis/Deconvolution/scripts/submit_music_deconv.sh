#!/usr/bin/env bash
set -euo pipefail

MICROMAMBA=${MICROMAMBA:-RNA-seq/in-house_MCD_RNAseq/micromamba}
ENV_PREFIX=${ENV_PREFIX:-RNA-seq/.mamba/music_deconv}
PARTITION=${PARTITION:-cpu}
CPUS_PER_TASK=${CPUS_PER_TASK:-8}
RUN_QC=${RUN_QC:-1}
RUN_SUMMARY=${RUN_SUMMARY:-1}

DATASETS=(
  "inhouse_MCD|mouse|RNA-seq/in-house_MCD_RNAseq/counts/featurecounts/gene_counts.txt|RNA-seq/in-house_MCD_RNAseq/metadata/samples.tsv"
  "GSE156918|mouse|RNA-seq/other_MCD_RNAseq/GSE156918/counts/featurecounts/gene_counts.txt|RNA-seq/other_MCD_RNAseq/GSE156918/metadata/samples.tsv"
  "GSE205974|mouse|RNA-seq/other_MCD_RNAseq/GSE205974/counts/featurecounts/gene_counts.txt|RNA-seq/other_MCD_RNAseq/GSE205974/metadata/samples.tsv"
  "GSE130970|human|RNA-seq/patient_RNAseq/results/GSE130970/counts/gene_counts_matrix.txt|RNA-seq/deconvolution/metadata/GSE130970_sra_metadata.tsv"
  "GSE135251|human|RNA-seq/patient_RNAseq/results/GSE135251/counts/gene_counts_matrix.txt|RNA-seq/deconvolution/metadata/GSE135251_sra_metadata.tsv"
)

get_sample_count() {
  local file=$1
  awk 'BEGIN{n=0} /^#/ {next} {if ($1=="Geneid") {print NF-6} else if ($1=="GeneID") {print NF-1} else {print NF-1} ; exit}' "$file"
}

mem_for_samples() {
  local samples=$1
  local mem
  mem=$(awk -v s="$samples" 'BEGIN{m=20 + (s*0.5); if (m < 40) m=40; if (m > 200) m=200; printf "%dG", (m+0.5)}')
  echo "$mem"
}

time_for_samples() {
  local samples=$1
  local hours
  hours=$(awk -v s="$samples" 'BEGIN{h=int((s+9)/10); if (h < 6) h=6; if (h > 48) h=48; print h}')
  if [[ "$hours" -ge 24 ]]; then
    local days=$((hours/24))
    local rem=$((hours%24))
    printf "%d-%02d:00:00" "$days" "$rem"
  else
    printf "%02d:00:00" "$hours"
  fi
}

job_ids=()

for entry in "${DATASETS[@]}"; do
  IFS='|' read -r name species counts_path metadata_path <<< "$entry"

  if [[ ! -f "$counts_path" ]]; then
    echo "Missing counts for ${name}: ${counts_path}" >&2
    continue
  fi

  if [[ ! -f "$metadata_path" ]]; then
    case "$name" in
      GSE130970)
        if [[ -x "${MICROMAMBA}" ]]; then
          "${MICROMAMBA}" run -p "${ENV_PREFIX}" python3 RNA-seq/deconvolution/scripts/00_convert_sra_metadata.py \
            --input RNA-seq/patient_RNAseq/data/metadata/GSE130970_SraRunTable.csv \
            --output "${metadata_path}"
        else
          python3 RNA-seq/deconvolution/scripts/00_convert_sra_metadata.py \
            --input RNA-seq/patient_RNAseq/data/metadata/GSE130970_SraRunTable.csv \
            --output "${metadata_path}"
        fi
        ;;
      GSE135251)
        if [[ -x "${MICROMAMBA}" ]]; then
          "${MICROMAMBA}" run -p "${ENV_PREFIX}" python3 RNA-seq/deconvolution/scripts/00_convert_sra_metadata.py \
            --input RNA-seq/patient_RNAseq/data/metadata/GSE135251_SraRunTable.csv \
            --output "${metadata_path}"
        else
          python3 RNA-seq/deconvolution/scripts/00_convert_sra_metadata.py \
            --input RNA-seq/patient_RNAseq/data/metadata/GSE135251_SraRunTable.csv \
            --output "${metadata_path}"
        fi
        ;;
      *)
        ;;
    esac
  fi

  if [[ ! -f "$metadata_path" ]]; then
    echo "Missing metadata for ${name}: ${metadata_path}" >&2
    continue
  fi

  samples=$(get_sample_count "$counts_path")
  if [[ -z "$samples" || "$samples" -le 0 ]]; then
    echo "Could not infer sample count for ${name}; defaulting to 50." >&2
    samples=50
  fi

  mem=$(mem_for_samples "$samples")
  time=$(time_for_samples "$samples")

  job_id=$(sbatch --parsable \
    --partition="$PARTITION" \
    --cpus-per-task="$CPUS_PER_TASK" \
    --mem="$mem" \
    --time="$time" \
    --job-name="music_${name}" \
    --export=ALL,DATASET_FILTER="$name",MICROMAMBA="$MICROMAMBA",ENV_PREFIX="$ENV_PREFIX",RUN_QC="$RUN_QC",RUN_SUMMARY=0 \
    RNA-seq/deconvolution/scripts/music_deconv_job.sh)

  echo "Submitted ${name} (samples=${samples}, mem=${mem}, time=${time}) -> ${job_id}"
  job_ids+=("$job_id")

done

if [[ "$RUN_SUMMARY" == "1" && ${#job_ids[@]} -gt 0 ]]; then
  dep=$(IFS=:; echo "${job_ids[*]}")
  summary_job=$(sbatch --parsable --dependency=afterok:${dep} RNA-seq/deconvolution/scripts/11_music_summary.sbatch)
  echo "Submitted summary job -> ${summary_job}"
fi
