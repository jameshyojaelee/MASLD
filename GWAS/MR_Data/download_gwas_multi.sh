#!/bin/bash
# Download three MASLD/NAFLD GWAS summary statistics
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data

echo "=== Downloading GWAS summary statistics ==="

# 1. Ghodsian 2021 (N=778K, EHR-based NAFLD meta-analysis)
# GCST90091033 — harmonised version (GRCh37, GWAS-SSF format)
echo "Downloading Ghodsian 2021 (harmonised)..."
wget -c "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90091001-GCST90092000/GCST90091033/harmonised/34841290-GCST90091033-EFO_0003095.h.tsv.gz" \
  -O Ghodsian_2021_NAFLD_harmonised.tsv.gz

# Also download original for comparison
echo "Downloading Ghodsian 2021 (original)..."
wget -c "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90091001-GCST90092000/GCST90091033/GCST90091033_buildGRCh37.tsv.gz" \
  -O Ghodsian_2021_NAFLD_original.tsv.gz

# 2. Chen 2023 GOLDPlus (N=691K, Nature Genetics)
# GCST90271622 — harmonised version
# WARNING: betas/SEs = 1 (Z-score meta-analysis). Use for TWAS only, NOT MR.
echo "Downloading Chen 2023 GOLDPlus (harmonised)..."
wget -c "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90271001-GCST90272000/GCST90271622/harmonised/GCST90271622.h.tsv.gz" \
  -O Chen_2023_GOLDPlus_harmonised.tsv.gz

# 3. Namjou 2019 eMERGE (N=9,677)
# GCST008468
echo "Downloading Namjou 2019 eMERGE..."
wget -c "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008468/NamjouB_31311600_NAFLD.txt" \
  -O Namjou_2019_NAFLD_GWAS.txt

echo "=== Downloads complete ==="
ls -lh Ghodsian_2021_* Chen_2023_* Namjou_2019_*
