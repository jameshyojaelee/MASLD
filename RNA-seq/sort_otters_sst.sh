#!/bin/bash
# Sort and bgzip OTTERS format eQTL files for tabix indexing
# OTTERS training.py requires position-sorted, bgzipped input

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OTTERS_DIR="${BASE}/data/broadaway_eqtl/otters_format"

echo "=== Sorting and bgzipping OTTERS SST files ==="
echo "Start: $(date)"

for CHR in $(seq 1 22); do
    RAW="${OTTERS_DIR}/chr${CHR}_broadaway.txt"
    SORTED="${OTTERS_DIR}/chr${CHR}_broadaway_sorted.txt"

    if [[ ! -f "${RAW}" ]]; then
        echo "chr${CHR}: no file, skipping"
        continue
    fi

    echo -n "chr${CHR}: sorting..."
    # Keep header, sort rest by position (col 2, numeric)
    head -1 "${RAW}" > "${SORTED}"
    tail -n +2 "${RAW}" | sort -k2,2n >> "${SORTED}"

    # Replace original with sorted version
    mv "${SORTED}" "${RAW}"
    echo " done ($(wc -l < "${RAW}") lines)"
done

echo ""
echo "=== All files sorted: $(date) ==="
