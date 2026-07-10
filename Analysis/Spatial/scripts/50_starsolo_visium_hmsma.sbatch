#!/bin/bash
#SBATCH --job-name=STARsolo
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/starsolo_hmsma_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/starsolo_hmsma_%A_%a.out
set -euo pipefail
# HMSMA Visium (10x Visium v1: R1=16bp CB + 12bp UMI, R2=150bp cDNA) aligned with
# STARsolo — SpaceRanger requires an H&E image we don't have; STARsolo needs none
# (spot barcodes -> spot x gene). Ref = GENCODE v49 (matches the bulk RNA-seq that
# defines the prioritized universe). Spots are NOT tissue-filtered here; that is done
# downstream by UMI threshold before pseudobulk-per-donor DE.
module load STAR/2.7.10b-GCC-12.2.0
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SAMPLES=$BASE/Analysis/Spatial/metadata/HRA007511_samples.txt
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SAMPLES")
IN=$BASE/Analysis/Spatial/spaceranger_input/HRA007511/$SAMPLE
OUT=$BASE/Analysis/Spatial/results/starsolo_visium/HRA007511/$SAMPLE
WL=/nfs/sw/easybuild/software/SpaceRanger/2.0.0-GCC-11.2.0/lib/python/cellranger/barcodes/visium-v1.txt
IDX=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/star_index

# skip if already done
if [ -f "$OUT/Solo.out/Gene/raw/matrix.mtx" ] || [ -f "$OUT/Solo.out/Gene/raw/matrix.mtx.gz" ]; then
  echo "SKIP $SAMPLE (matrix exists)"; exit 0
fi
mkdir -p "$OUT"; cd "$OUT"
# gather all lanes (a sample can span >1 HRR run -> L001,L002,...); cDNA(R2) first, then CB(R1)
R2=$(ls "$IN"/${SAMPLE}_S1_L*_R2_001.fastq.gz | paste -sd,)
R1=$(ls "$IN"/${SAMPLE}_S1_L*_R1_001.fastq.gz | paste -sd,)
echo "SAMPLE=$SAMPLE"; echo "R2(cDNA)=$R2"; echo "R1(CB+UMI)=$R1"

STAR --runMode alignReads \
  --genomeDir "$IDX" \
  --readFilesIn "$R2" "$R1" \
  --readFilesCommand zcat \
  --runThreadN 16 \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist "$WL" \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloFeatures Gene \
  --soloCBmatchWLtype 1MM_multi \
  --soloUMIdedup 1MM_All \
  --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
  --outSAMtype None \
  --outFileNamePrefix "$OUT/"
gzip -f "$OUT"/Solo.out/Gene/raw/*.tsv "$OUT"/Solo.out/Gene/raw/*.mtx 2>/dev/null || true
echo "DONE $SAMPLE $(date)"
