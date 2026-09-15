#!/bin/bash
#SBATCH --job-name=sldsc_annot
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=GWAS/ldsc/program_heritability_20260830/logs/20_annot_%j.out
#SBATCH --error=GWAS/ldsc/program_heritability_20260830/logs/20_annot_%j.err
#
# Mark which reference panel SNPs fall inside each of the 123 gene windows
# (117 programs plus 6 backgrounds), one file per chromosome holding all 123
# columns. Then measure how correlated the 117 program annotations are, which
# sets the effective number of independent tests reported alongside the
# Benjamini-Hochberg correction.

set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"

BASE=GWAS/ldsc/program_heritability_20260830
PY=RNA-seq/.mamba/unified_env/bin/python
PLINK=data/gsmap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink

$PY $BASE/scripts/01_build_gene_windows.py

seq 1 22 | xargs -P 8 -I{} $PY $BASE/scripts/02_make_annot.py {} $PLINK $BASE/annot

echo ">> per-chromosome SNP counts written; combining"
$PY - <<'PYEOF'
import glob, pandas as pd
frames = [pd.read_csv(p, sep="\t") for p in
          sorted(glob.glob("GWAS/ldsc/program_heritability_20260830/annot/programs.*.snp_counts.tsv"))]
total = pd.concat(frames).groupby("annotation", as_index=False)[["n_snps", "n_snps_total"]].sum()
total["frac_reference_snps"] = total["n_snps"] / total["n_snps_total"]
total.sort_values("n_snps").to_csv(
    "GWAS/ldsc/program_heritability_20260830/annot/annotation_snp_counts.tsv",
    sep="\t", index=False)
print(total["n_snps"].describe())
PYEOF

$PY $BASE/scripts/04_effective_tests.py $BASE/annot $BASE/gene_sets/effective_tests.tsv

echo "ANNOTATIONS DONE"
