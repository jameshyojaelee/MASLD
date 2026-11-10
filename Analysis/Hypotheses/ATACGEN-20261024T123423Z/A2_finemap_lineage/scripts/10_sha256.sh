#!/bin/bash
# sha256 of every input A2 reads (run on a compute node; large files).
set -euo pipefail
R=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
A=$R/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage
C=$R/Analysis/Multimodal_Program_Projection/candidates
SC=$R/Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk
{
  echo -e "sha256\tpath"
  for f in \
    $R/GWAS/finemapping/results/credible_sets.csv \
    $R/GWAS/finemapping/results/susie_all_results.csv \
    $R/GWAS/finemapping/config/gwas_trait_tier.tsv \
    $R/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_pip_architecture.tsv \
    $R/GWAS/finemapping/data/ld_ref/1kg_eur/approx_LD_blocks.txt \
    $C/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv \
    $C/atac-context-v3-candidate-2026-08-11-r1/genetics/context/variant_liftover_audit.tsv \
    $C/atac-context-v3-candidate-2026-08-11-r1/consensus_peak_manifest.tsv \
    $C/atac-context-v3-statistics-correction-20261022/da/da_peak_results.tsv.gz \
    $R/RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/annotation/GENCODE_v49/run_19749483/gene_identity.tsv \
    $R/data/broadaway_eqtl/hg19ToHg38.over.chain \
    $SC/sample_celltype_metadata.tsv \
    $SC/pseudobulk_*.tsv.gz \
    $R/data/GSE244832/metadata/donor_pairing.csv $R/data/GSE185477/metadata/donor_pairing.csv \
    $R/data/GSE202379/metadata/donor_pairing.csv $R/data/GSE136103/metadata/donor_pairing.csv; do
    [[ -e $f ]] && sha256sum "$f" | awk -v OFS='\t' '{print $1,$2}'
  done
  # panel files: one digest over the sorted per-file digests
  for d in panel_af 1kg_afr 1kg_amr 1kg_eas 1kg_eur 1kg_sas; do
    h=$(find $R/GWAS/finemapping/data/ld_ref/$d -type f \( -name '*.af.tsv.gz' -o -name '*.bed' -o -name '*.bim' -o -name '*.fam' \) | sort | xargs sha256sum | sha256sum | cut -d' ' -f1)
    echo -e "$h\t$R/GWAS/finemapping/data/ld_ref/$d (digest of sorted per-file sha256)"
  done
} > $A/inputs_sha256.tsv
