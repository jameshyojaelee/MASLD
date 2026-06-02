#!/bin/bash
# Fetch SRR lists for given GSE accessions using NCBI Entrez Direct.

GSE_LIST=$1
OUTDIR=$2

if [ -z "$GSE_LIST" ] || [ -z "$OUTDIR" ]; then
    echo "Usage: $0 \"GSE1 GSE2 ...\" OUTDIR"
    exit 1
fi

mkdir -p "$OUTDIR"

for gse in $GSE_LIST; do
    echo "Fetching SRRs for $gse..."
    out_file="$OUTDIR/${gse}_SRR_list.txt"
    
    # esearch the SRA database for the GSE, efetch the runinfo CSV, extract the first column (Run ID), and filter for SRR/ERR
    esearch -db sra -query "$gse" | efetch -format runinfo | cut -d ',' -f 1 | grep -E "^[SE]RR" | sort -u > "$out_file"
    
    count=$(wc -l < "$out_file")
    if [ "$count" -eq 0 ]; then
        echo "Warning: No SRRs found for $gse."
    else
        echo "Saved $count SRRs to $out_file"
    fi
    
    sleep 2 # Honor NCBI rate limit policies
done
echo "All done."
