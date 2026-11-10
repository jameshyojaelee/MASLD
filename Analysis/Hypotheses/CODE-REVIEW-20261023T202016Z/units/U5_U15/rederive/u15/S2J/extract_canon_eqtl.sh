#!/bin/bash
# Pull Broadaway liver cis-eQTL gene-variant rows for the 21 canonical genes of Figure S2J
# (exact GeneSymbol match, column 13). Output: CHR POS NEA EA GeneSymbol ENSG.
set -uo pipefail
for ch in $(seq 1 22); do
  grep -wF -f /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/u15/S2J/canon_genes.txt /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/broadaway_eqtl/chr${ch}_marginal_summary_results.tsv |
    awk -F'\t' 'NR==FNR{g[$1]=1;next} ($13 in g){print $3"\t"$4"\t"$5"\t"$6"\t"$13"\t"$14}' /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/u15/S2J/canon_genes.txt -
done
exit 0
