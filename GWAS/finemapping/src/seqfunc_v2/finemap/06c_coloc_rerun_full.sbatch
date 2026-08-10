#!/bin/bash
#SBATCH --job-name=susie_coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_rerun/coloc_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_rerun/coloc_%A_%a.err

# FULL-PORTFOLIO SuSiE-COLOC rerun: all 50 registry studies x 22 chromosomes.
#
# Supersedes 06b_coloc_gram_nonEUR.sbatch, which covered only the 29 non-EUR
# studies because the rebuilt LD panels were the sole change at the time. Two
# further changes now affect EVERY study, so the rerun must be portfolio-wide:
#
#   1. Rebuilt 1000G LD panels (non-EUR only) -- mean-imputed Gram, PSD by
#      construction. MEASURED IMPACT: none on COLOC. PIP correlation >=0.999998
#      and 0/39 lead-variant changes on the worst panels (min eigenvalue -0.62).
#   2. Regenerated eQTL SuSiE fits with LD_PANEL PINNED to polyfun and
#      ld_provenance stamped into every RDS. The previous fits used whatever
#      get_ld_base_dir() defaulted to at run time, unrecorded.
#   3. Palindromic variants RETAINED with an ancestry-matched frequency guard
#      instead of dropped outright (~12.6% of variants per locus, which excluded
#      PNPLA3 rs738409 from every locus). MEASURED IMPACT: ABF 0 threshold
#      crossings in 1,383 comparisons; SuSiE 1 in 109, that one a gene at
#      0.4965 vs 0.5006.
#
# PREDICTION ON RECORD: this should change almost nothing. A large shift is
# evidence of a problem with the rerun, not a discovery.
#
# KNOWN, ACCEPTED LOSS: 62 chr6 genes present in canonical are absent here.
# polyfun_eur deliberately omits the MHC -- exactly 5 blocks (chr6 112 defined,
# 107 on disk) spanning 28.0-32.7 Mb -- so those genes match 0 LD variants and
# are skipped. Canonical obtained them via an unrecorded fallback to a panel
# that includes the MHC, which is the very provenance defect the pin removes.
#
# SAFETY: COLOC_OUT_SUFFIX=_rerun writes to results/susie_coloc_rerun/<gwas>/.
# NOTHING canonical is touched -- gene_level_coloc.csv is firewall-guarded and
# feeds the atlas, so promotion is a separate, explicit decision.
#
# ENV: rnaseq, matching what produced the canonical results, so the before/after
# is not confounded by the data.table version difference against finemapping.
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
FM="${PROJECT_ROOT}/GWAS/finemapping"
SRC="${FM}/src"
LDREF="${FM}/data/ld_ref"
export MASLD_PROJECT_ROOT

# Rebuilt Gram panels: per-ancestry overrides beat the LD_PANEL dispatch in
# get_ld_base_dir(). EUR is unaffected (PolyFun, not rebuilt).
export AFR_LD_DIR="${LDREF}/1kg_afr_gram"
export AMR_LD_DIR="${LDREF}/1kg_amr_gram"
export EAS_LD_DIR="${LDREF}/1kg_eas_gram"
export SAS_LD_DIR="${LDREF}/1kg_sas_gram"
# Regenerated, LD-pinned eQTL SuSiE fits.
export EQTL_SUSIE_DIR="${FM}/results/eqtl_susie_polyfun"
export COLOC_OUT_SUFFIX="_rerun"

export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Study list comes from the registry itself rather than a hard-coded array, so it
# cannot drift out of sync the way the 29-study list in 06b did.
mapfile -t GWAS_LIST < <(awk -F'\t' 'NR>1 && $1!="" {print $1}' "${FM}/config/gwas_registry.tsv" | sort -u)
N_GWAS=${#GWAS_LIST[@]}

IDX=${SLURM_ARRAY_TASK_ID}
GI=$(( IDX / 22 ))
CHR=$(( IDX % 22 + 1 ))
if [ "$GI" -ge "$N_GWAS" ]; then echo "task ${IDX} beyond ${N_GWAS} studies -- nothing to do"; exit 0; fi
GWAS=${GWAS_LIST[$GI]}

echo "[$(date -u +%FT%TZ)] task ${IDX} -> ${GWAS} chr${CHR}  (gram LD + regenerated eQTL, suffix _rerun)"
cd "${SRC}"
GWAS_NAME="${GWAS}" CHR_FILTER="${CHR}" Rscript "${SRC}/06_susie_coloc.R"
echo "[$(date -u +%FT%TZ)] done ${GWAS} chr${CHR}"
