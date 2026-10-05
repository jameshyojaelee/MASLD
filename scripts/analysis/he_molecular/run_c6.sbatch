#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
# C6: Model C external test in TCGA-LIHC non-tumor liver (Amendment 6). Runs once.
# Downloads the Xena STAR TPM matrix, scores the targets (c3c), then runs c6.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%j.out run_c6.sbatch <c1b_exec_dir> <gtex_c2_exec_dir> <tcga_exec_dir>
set -euo pipefail
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mm=/gpfs/commons/home/jameslee/.local/bin/micromamba
he=$repo/scripts/analysis/he_molecular
ext=$repo/data/external/tcga_lihc
export PYTHONNOUSERSITE=1
mkdir -p "$ext" "$3"
tpm=$ext/TCGA-LIHC.star_tpm.tsv.gz
[ -s "$tpm" ] || curl -sS --fail -o "$tpm" https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.star_tpm.tsv.gz
sha256sum "$tpm" > "$3/star_tpm.sha256"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$he/c3c_tcga_targets.R" "$1" "$tpm" "$3/targets"
"$mm" run -n he_molecular python -m pip freeze > "$3/pip_freeze_c6.txt"
cd "$he"
"$mm" run -n he_molecular python c6_external_test.py --c2 "$2" \
  --tcga-targets "$3/targets/targets_tcga.tsv" --tcga-features "$3/embed/slide_features.tsv.gz" --out "$3/c6_external_test"
