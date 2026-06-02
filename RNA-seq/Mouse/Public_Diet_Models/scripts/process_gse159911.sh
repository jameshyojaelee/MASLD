#!/bin/bash
#SBATCH --job-name=extract_gse159911
#SBATCH --output=GSE159911/logs/extract_%j.log
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=6
#SBATCH --mem=12G
#SBATCH --time=4:00:00

# Setup
PROJECT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"
GSE_DIR="${PROJECT_DIR}/GSE159911"
SRA_CACHE="${GSE_DIR}/sra_cache"
FASTQ_DIR="${GSE_DIR}/fastq"
mkdir -p "$FASTQ_DIR" "$GSE_DIR/logs"

# Activate Env (hardcoded path as per existing usage)
eval "$(micromamba shell hook --shell bash)"
micromamba activate /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl

echo "Starting extraction for GSE159911..."

# Iterate over downloaded SRA directories/files
# prefetch usually creates sra_cache/SRRxxxx/SRRxxxx.sra
for sra_dir in ${SRA_CACHE}/SRR*; do
    if [ -d "$sra_dir" ]; then
        srr=$(basename "$sra_dir")
        sra_file="${sra_dir}/${srr}.sra"
        
        if [ -f "$sra_file" ]; then
            echo "Processing ${srr}..."
            
            # Extract
            # Use temp dir to avoid fastq clutter if valid
            fasterq-dump --split-files --threads 6 --outdir "$FASTQ_DIR" "$sra_file"
            
            # Compress (pigz)
            echo "Compressing ${srr}..."
            pigz -p 6 "$FASTQ_DIR/${srr}"*.fastq
            
            # Verify
            if ls "$FASTQ_DIR/${srr}"*_1.fastq.gz 1> /dev/null 2>&1; then
                echo "${srr} processed successfully. Removing .sra"
                rm "$sra_file"
                rmdir "$sra_dir"
            else
                echo "Error: Output fastq.gz not found for ${srr}"
            fi
        fi
    fi
done

echo "Extraction complete."

# Update Metadata
echo "Updating samples.tsv..."
cd "$PROJECT_DIR"
python scripts/generate_samples_sheet.py

echo "Done."
