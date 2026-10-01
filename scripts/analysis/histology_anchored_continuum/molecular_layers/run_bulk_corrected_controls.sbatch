#!/usr/bin/env bash
#SBATCH --job-name=hac_ml_bulk_corr
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

# Continuum bulk transcript layer (10_bulk_transcript_deg.R) on the corrected-
# control F_five refit. The contract's F_five inputs and family sizes are
# replaced through HAC_CONTRACT_OVERRIDE; every other contract input is unchanged.
# Previous hac_ml_bulk runs used 2.4-2.7 GB and 15-27 min on one thread.
set -euo pipefail
umask 027
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
: "${HAC_ML_OUT_ROOT:?HAC_ML_OUT_ROOT is unset}"
: "${POOLED_ARM_ROOT:?POOLED_ARM_ROOT (validated refit arm) is unset}"
: "${STAGE_EXTENSION_ROOT:?STAGE_EXTENSION_ROOT is unset}"
test -f "${POOLED_ARM_ROOT}/VALIDATED.json"
test -f "${STAGE_EXTENSION_ROOT}/output_manifest.tsv"
test ! -e "${HAC_ML_OUT_ROOT}"
mkdir "${HAC_ML_OUT_ROOT}"

override="${HAC_ML_OUT_ROOT}/contract_override.json"
"${RNASEQ_ENV}/bin/python3" - "${POOLED_ARM_ROOT}" "${STAGE_EXTENSION_ROOT}" "${override}" <<'PY'
import csv, json, sys
from pathlib import Path
arm, stage, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
expected = json.loads((arm / "expected_counts.json").read_text())
with (arm / "model_input_manifest.tsv").open(newline="") as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))
assert len(manifest) == 1
override = {
    "dge_rds": manifest[0]["dge_path"],
    "canonical_deg_results": str(arm / "deg_results.csv"),
    "sample_manifest": str(stage / "five_cohort_sample_manifest.tsv"),
    "stage_results": str(stage / "stage_all_gene_results.tsv"),
    "stage_extension_results": str(stage / "stage_extension_all_gene_results.tsv"),
    "gene_family_size": expected["n_genes"],
    "canonical_deg_size": expected["canonical"]["n"],
}
out.write_text(json.dumps(override, indent=2) + "\n")
print(json.dumps(override, indent=2))
PY

export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export HAC_ML_OUT_ROOT
export HAC_CONTRACT_OVERRIDE="${override}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"${RNASEQ_ENV}/bin/Rscript" --vanilla "${SCRIPT_DIR}/10_bulk_transcript_deg.R"
echo "HAC_MOLECULAR_BULK_CORRECTED_COMPLETE ${HAC_ML_OUT_ROOT}"
