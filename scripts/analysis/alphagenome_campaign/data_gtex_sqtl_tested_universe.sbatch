#!/bin/bash
#SBATCH --job-name=masld_gtex_sqtl_universe
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_gtex_sqtl_universe_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_gtex_sqtl_universe_%j.err
set -euo pipefail
# The deposited splice work is bounded by an outcome-selected population: only
# GTEx v8 significant sQTL pairs are on disk, so no all-tested denominator exists
# and `new_effect_evaluation_ready_rows` is 0. The per-SNP all-associations files
# are not published in the public adult-gtex bucket (its
# sQTL-Tissue-Specific-All-SNP-Gene-Associations/ prefix holds only a zero-byte
# placeholder). The tested PHENOTYPE universe is published: the LeafCutter intron
# clusters actually quantified in Liver, plus the cluster-to-gene grouping. That
# is the denominator the splice family needs to stop being significance-selected.
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/gtex_v8_liver_sqtl
SRC=https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl
cd "$BASE"
for f in GTEx_Analysis_v8_sQTL_groups.tar.gz GTEx_Analysis_v8_sQTL_covariates.tar.gz GTEx_Analysis_v8_sQTL_phenotype_matrices.tar; do
  if [ ! -s "$f" ]; then
    echo "fetching $f"
    wget --no-verbose --tries=3 --timeout=120 --continue -O "$f" "$SRC/$f"
  else
    echo "already present: $f"
  fi
  sha256sum "$f" | tee -a fetch_receipt.sha256
done
mkdir -p liver_tested_universe
tar -xvf GTEx_Analysis_v8_sQTL_phenotype_matrices.tar -C liver_tested_universe --wildcards '*Liver*' || true
tar -xzvf GTEx_Analysis_v8_sQTL_groups.tar.gz -C liver_tested_universe --wildcards '*Liver*' || true
tar -xzvf GTEx_Analysis_v8_sQTL_covariates.tar.gz -C liver_tested_universe --wildcards '*Liver*' || true
ls -la liver_tested_universe
echo GTEX_SQTL_TESTED_UNIVERSE_DONE
