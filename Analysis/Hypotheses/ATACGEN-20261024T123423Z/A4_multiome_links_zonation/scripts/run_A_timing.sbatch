#!/bin/bash
#SBATCH --job-name=a4_timing
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation/logs/%x_%j.log
# A4 timing job: inputs (00), zonation pseudobulk (01), hepatocyte links with 2 bootstrap replicates (02 timing).
set -euo pipefail
A4=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A4_multiome_links_zonation
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$A4"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
EXE=$ROOT/Analysis/MASLD_Model_Benchmark/executions
H24=$ROOT/Analysis/Hypotheses/H24-20261008T211750Z/gse296875
sha256sum $H24/atac_peak_by_unit.f32 $H24/rna_gene_by_unit.f32 $H24/units.tsv $H24/peaks.tsv \
  $EXE/gse296875-donor-lineage-pseudobulk-21107249/pseudobulk/genes.tsv \
  $EXE/gse296875-fragment-membership-21063829/cell_membership.tsv.gz \
  $EXE/gse296875-crossmodal-recurrence-20260830T193000Z/rectangle/atac.data.i32 \
  $EXE/gse296875-crossmodal-recurrence-20260830T193000Z/rectangle/atac.indices.i32 \
  $EXE/gse296875-crossmodal-recurrence-20260830T193000Z/rectangle/atac.indptr.i32 \
  $EXE/gse296875-crossmodal-recurrence-20260830T193000Z/rectangle/nuclei.txt \
  $EXE/gse296875-phenotype-endpoints-20260825/endpoint_lock/donor_endpoints.tsv \
  $ROOT/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv \
  $ROOT/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz \
  $ROOT/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv \
  $ROOT/GWAS/finemapping/config/gwas_trait_tier.tsv \
  $ROOT/data/broadaway_eqtl/hg19ToHg38.over.chain \
  /gpfs/commons/home/jameslee/reference_genome/cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz \
  > input_sha256.txt
$ENV/bin/python -I -m pip freeze > pip_freeze.txt 2>/dev/null || true
/usr/bin/time -v $ENV/bin/Rscript scripts/00_inputs.R "$A4" 2> logs/time_00.txt
tail -n 25 logs/time_00.txt
/usr/bin/time -v $ENV/bin/python -I scripts/01_zonation_pseudobulk.py "$A4" 2> logs/time_01.txt
tail -n 25 logs/time_01.txt
/usr/bin/time -v $ENV/bin/Rscript scripts/02_links.R "$A4" timing hepatocyte 2 2 2> logs/time_02_timing.txt
tail -n 60 logs/time_02_timing.txt
