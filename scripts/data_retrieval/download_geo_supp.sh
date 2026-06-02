#!/bin/bash
#SBATCH --job-name=download_geo_supp
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err
#SBATCH --partition=cpu,io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00

# Usage: sbatch download_geo_supp.sh "GSE163211 GSE202379 GSE244832 GSE289173 GSE292268" <OUTPUT_DIR>

GSE_LIST=$1
OUT_DIR=$2

if [ -z "$GSE_LIST" ] || [ -z "$OUT_DIR" ]; then
    echo "Error: Must provide GSE_LIST and OUT_DIR"
    exit 1
fi

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

for gse in $GSE_LIST; do
    echo "Downloading supplementary files for $gse..."
    
    mkdir -p "$gse"
    cd "$gse"
    
    # Extract the NNN format required for NCBI GEO FTP
    # e.g. GSE163211 -> GSE163nnn
    gse_prefix=$(echo "$gse" | sed -r 's/[0-9]{3}$/nnn/')
    ftp_url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${gse_prefix}/${gse}/suppl/"
    
    echo "Fetching from $ftp_url"
    
    # Use wget mirror to download all files in the directory
    wget -np -nd -r -A "*.tar*,*gz*,*tsv*,*csv*,*h5*,*txt*" "$ftp_url"
    
    if [ $? -ne 0 ]; then
        echo "Warning: wget encountered an error for $gse. It might not have supplementary files or the FTP timed out."
    else
        echo "Completed download for $gse"
    fi
    
    cd ..
done

echo "All supplementary downloads complete."
