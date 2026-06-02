#!/bin/bash
# 63_otters_twas.sh
# ---------------------------------------------------------------------------
# Phase 2A Step 3: Run OTTERS TWAS using trained weights + GWAS summary stats
#
# Input:
#   - data/broadaway_eqtl/otters_weights/chr{N}/  (from Script 62)
#   - GWAS/MR_Data/hg19/*_hg19.tsv.gz  (from Script 60)
#   - data/broadaway_eqtl/otters_format/gene_anno.txt
#   - data/1kg_eur/bed/chr{N}_eur.{bed,bim,fam}  (from Script 62 side effect)
#
# Output:
#   - RNA-seq/results/causal_inference/otters_broadaway/{GWAS}/
#
# Usage:
#   GWAS_NAME=ghodsian sbatch RNA-seq/run_otters_twas.sbatch
# ---------------------------------------------------------------------------

set -euo pipefail

# Configuration
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OTTERS_DIR="${BASE}/tools/OTTERS"
BROADAWAY_OTTERS="${BASE}/data/broadaway_eqtl/otters_format"
WEIGHT_DIR="${BASE}/data/broadaway_eqtl/otters_weights"
BED_DIR="${BASE}/data/1kg_eur/bed"
GWAS_HG19_DIR="${BASE}/GWAS/MR_Data/hg19"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

GWAS_NAME="${GWAS_NAME:-ghodsian}"
ANNO_FILE="${BROADAWAY_OTTERS}/gene_anno.txt"

echo "=== OTTERS TWAS: ${GWAS_NAME} ==="
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Start: $(date)"
echo ""

# Map GWAS name to file and column info
case "${GWAS_NAME}" in
    ghodsian)
        GWAS_FILE="${GWAS_HG19_DIR}/Ghodsian_2021_NAFLD_harmonised_hg19.tsv.gz"
        CHR_COL="hg19_chr"; POS_COL="hg19_pos"
        EA_COL="hm_effect_allele"; OA_COL="hm_other_allele"
        BETA_COL="hm_beta"; SE_COL="standard_error"
        ;;
    ukbb_alt)
        GWAS_FILE="${GWAS_HG19_DIR}/GCST90019492_UKBB_ALT_harmonised_hg19.tsv.gz"
        CHR_COL="hg19_chr"; POS_COL="hg19_pos"
        EA_COL="effect_allele"; OA_COL="other_allele"
        BETA_COL="beta"; SE_COL="standard_error"
        ;;
    ukbb_ast)
        GWAS_FILE="${GWAS_HG19_DIR}/GCST90019497_UKBB_AST_harmonised_hg19.tsv.gz"
        CHR_COL="hg19_chr"; POS_COL="hg19_pos"
        EA_COL="effect_allele"; OA_COL="other_allele"
        BETA_COL="beta"; SE_COL="standard_error"
        ;;
    ukbb_ggt)
        GWAS_FILE="${GWAS_HG19_DIR}/GCST90019507_UKBB_GGT_harmonised_hg19.tsv.gz"
        CHR_COL="hg19_chr"; POS_COL="hg19_pos"
        EA_COL="effect_allele"; OA_COL="other_allele"
        BETA_COL="beta"; SE_COL="standard_error"
        ;;
    pdff)
        GWAS_FILE="${GWAS_HG19_DIR}/GCST90267352_PDFF_Pazoki2022_hg19.tsv.gz"
        CHR_COL="hg19_chr"; POS_COL="hg19_pos"
        EA_COL="effect_allele"; OA_COL="other_allele"
        BETA_COL="beta"; SE_COL="standard_error"
        ;;
    *)
        echo "ERROR: Unknown GWAS_NAME: ${GWAS_NAME}"
        exit 1
        ;;
esac

if [[ ! -f "${GWAS_FILE}" ]]; then
    echo "ERROR: GWAS file not found: ${GWAS_FILE}"
    exit 1
fi

# Output directory
OUT_BASE="${BASE}/RNA-seq/results/causal_inference/otters_broadaway/${GWAS_NAME}"
mkdir -p "${OUT_BASE}"

# Ensure plink is in PATH (OTTERS calls 'plink')
PLINK19="/nfs/sw/easybuild/software/PLINK/1.9b_6.21-x86_64/plink"
TMPBIN=$(mktemp -d)
ln -s "${PLINK19}" "${TMPBIN}/plink"
export PATH="${TMPBIN}:${PATH}"

# Step 1: Format GWAS for OTTERS (per-chromosome)
echo "Formatting GWAS for OTTERS..."
GWAS_OTTERS_DIR="${OUT_BASE}/gwas_formatted"
mkdir -p "${GWAS_OTTERS_DIR}"

# Export variables for Python heredoc (single-quoted heredoc prevents bash expansion)
export GWAS_FILE CHR_COL POS_COL EA_COL OA_COL BETA_COL SE_COL GWAS_OTTERS_DIR

python3 - << 'PYEOF'
import pandas as pd
import numpy as np
import sys
import os
import gzip

gwas_file = os.environ['GWAS_FILE']
chr_col = os.environ['CHR_COL']
pos_col = os.environ['POS_COL']
ea_col = os.environ['EA_COL']
oa_col = os.environ['OA_COL']
beta_col = os.environ['BETA_COL']
se_col = os.environ['SE_COL']
out_dir = os.environ['GWAS_OTTERS_DIR']

print(f"  Reading {os.path.basename(gwas_file)}...")
df = pd.read_csv(gwas_file, sep='\t', dtype={chr_col: 'str', pos_col: 'Int64'})
print(f"  Total rows: {len(df):,}")

# Filter valid rows
df = df.dropna(subset=[chr_col, pos_col, beta_col, se_col])
df[chr_col] = pd.to_numeric(df[chr_col], errors='coerce')
df = df.dropna(subset=[chr_col])
df[chr_col] = df[chr_col].astype(int)
df[pos_col] = df[pos_col].astype(int)

# Compute Z-score
df['Z'] = df[beta_col].astype(float) / df[se_col].astype(float)
df = df[np.isfinite(df['Z'])]
print(f"  Valid variants: {len(df):,}")

# Write per-chromosome GWAS files
# OTTERS GWAS format: CHROM POS A1 A2 Z
for chrom in range(1, 23):
    chr_df = df[df[chr_col] == chrom]
    if len(chr_df) == 0:
        continue

    out = pd.DataFrame({
        'CHROM': chr_df[chr_col].values,
        'POS': chr_df[pos_col].values,
        'A1': chr_df[ea_col].str.upper().values,
        'A2': chr_df[oa_col].str.upper().values,
        'Z': chr_df['Z'].values
    })

    out_file = os.path.join(out_dir, f'chr{chrom}_gwas.txt')
    out.to_csv(out_file, sep='\t', index=False)
    print(f"  chr{chrom}: {len(out):,} variants")

print("  GWAS formatting complete.")
PYEOF

echo ""

# Step 2: Run OTTERS TWAS per chromosome
echo "Running OTTERS TWAS per chromosome..."

for CHROM in $(seq 1 22); do
    WEIGHT_CHR="${WEIGHT_DIR}/chr${CHROM}"
    BED_PREFIX="${BED_DIR}/chr${CHROM}_eur"
    GWAS_CHR="${GWAS_OTTERS_DIR}/chr${CHROM}_gwas.txt"
    CHR_OUT="${OUT_BASE}/chr${CHROM}"

    # Skip if no weights or GWAS for this chromosome
    if [[ ! -d "${WEIGHT_CHR}" ]]; then
        echo "  chr${CHROM}: no weights, skipping"
        continue
    fi
    if [[ ! -f "${GWAS_CHR}" ]]; then
        echo "  chr${CHROM}: no GWAS data, skipping"
        continue
    fi
    if [[ ! -f "${BED_PREFIX}.bed" ]]; then
        echo "  chr${CHROM}: no BED genotype, skipping"
        continue
    fi

    mkdir -p "${CHR_OUT}"

    echo "  chr${CHROM}: running TWAS..."
    python "${OTTERS_DIR}/testing.py" \
        --OTTERS_dir="${OTTERS_DIR}" \
        --anno_dir="${ANNO_FILE}" \
        --geno_dir="${BED_PREFIX}" \
        --weight_dir="${WEIGHT_CHR}" \
        --gwas_file="${GWAS_CHR}" \
        --out_dir="${CHR_OUT}" \
        --chrom="${CHROM}" \
        --models=P0.05,P0.001 \
        --window=1000000 \
        --thread="${THREADS}" 2>&1 | tail -5
done

# Step 3: Combine results across chromosomes with ACAT
echo ""
echo "Combining TWAS results with ACAT..."

export OUT_BASE GWAS_NAME ANNO_FILE

python3 - << 'PYEOF'
import pandas as pd
import numpy as np
from scipy import stats
import os
import glob

out_base = os.environ['OUT_BASE']
gwas_name = os.environ['GWAS_NAME']
anno_file = os.environ['ANNO_FILE']

# Load gene annotations for symbol mapping
anno = pd.read_csv(anno_file, sep='\t')
gene_map = dict(zip(anno['TargetID'], anno['GeneName']))

# Collect all per-chromosome results
all_results = []
for chrom in range(1, 23):
    chr_dir = os.path.join(out_base, f'chr{chrom}')
    # OTTERS outputs results in various formats; check for common patterns
    for pattern in ['*_ACAT.txt', '*_results.txt', '*_twas.txt', '*.results']:
        files = glob.glob(os.path.join(chr_dir, pattern))
        for f in files:
            try:
                df = pd.read_csv(f, sep='\t')
                if len(df) > 0:
                    df['chr'] = chrom
                    all_results.append(df)
            except Exception:
                pass

if not all_results:
    # Try to find individual gene results
    for chrom in range(1, 23):
        chr_dir = os.path.join(out_base, f'chr{chrom}')
        if not os.path.isdir(chr_dir):
            continue
        for f in glob.glob(os.path.join(chr_dir, '*.txt')):
            try:
                df = pd.read_csv(f, sep='\t')
                if len(df) > 0 and any(c in df.columns for c in ['Zscore', 'TWAS_Z', 'P', 'pvalue']):
                    df['chr'] = chrom
                    all_results.append(df)
            except Exception:
                pass

if not all_results:
    print("WARNING: No OTTERS TWAS results found. Check weight training output.")
    # Create empty output
    empty = pd.DataFrame(columns=['gene', 'gene_symbol', 'chr', 'otters_acat_z',
                                   'otters_acat_pval', 'fdr'])
    empty.to_csv(os.path.join(out_base, 'otters_twas_combined.csv'), index=False)
else:
    combined = pd.concat(all_results, ignore_index=True)
    print(f"Combined results: {len(combined):,} entries")

    # Map gene symbols
    if 'TargetID' in combined.columns:
        combined['gene_symbol'] = combined['TargetID'].map(gene_map)

    # Compute FDR
    if 'P' in combined.columns or 'pvalue' in combined.columns:
        pcol = 'P' if 'P' in combined.columns else 'pvalue'
        from statsmodels.stats.multitest import multipletests
        valid = combined[pcol].notna() & (combined[pcol] > 0)
        combined.loc[valid, 'fdr'] = multipletests(
            combined.loc[valid, pcol].values, method='fdr_bh')[1]

    out_file = os.path.join(out_base, 'otters_twas_combined.csv')
    combined.to_csv(out_file, index=False)
    print(f"Saved: {out_file}")

    # Summary
    if 'fdr' in combined.columns:
        n_sig = (combined['fdr'] < 0.05).sum()
        print(f"FDR < 0.05: {n_sig}")

print("ACAT combination complete.")
PYEOF

# Cleanup
rm -rf "${TMPBIN}"

echo ""
echo "=== OTTERS TWAS complete for ${GWAS_NAME} ==="
echo "Results: ${OUT_BASE}/"
ls -lh "${OUT_BASE}"/*.csv 2>/dev/null
echo "End: $(date)"
