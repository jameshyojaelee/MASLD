#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/joint_cluster_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/joint_cluster_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2
ENV=${ROOT}/.mamba/seqfunc_leafcutter
export PYTHONNOUSERSITE=1 MASLD_PROJECT_ROOT=${ROOT}
micromamba run -n rnaseq Rscript "${ROOT}/GWAS/finemapping/src/seqfunc_v2/splicing/00_prepare_manifest.R"
test "$(wc -l < "${OUT}/joint/juncfiles.txt")" -eq 846
if [[ ! -s "${OUT}/joint/joint846_perind_numers.counts.gz" ]]; then
  micromamba run -p "${ENV}" leafcutter-cluster -j "${OUT}/joint/juncfiles.txt" -r "${OUT}/joint" -o joint846 -l 500000 -m 30 -p 0.001
fi
test -s "${OUT}/joint/joint846_perind_numers.counts.gz"
# LeafCutter's numerized count header is the 846 sample IDs themselves (there
# is no leading cluster-name field). Validate the exact set without a
# zcat|head pipeline, which returns SIGPIPE=141 under `set -o pipefail`.
micromamba run -n rnaseq python - "${OUT}/joint/joint846_perind_numers.counts.gz" "${OUT}/canonical_846_manifest.tsv" <<'PY'
import csv, gzip, sys
counts, manifest = sys.argv[1:]
with gzip.open(counts, "rt") as handle:
    header = handle.readline().strip().split()
with open(manifest) as handle:
    expected = [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")]
if len(header) != 846 or len(set(header)) != 846 or set(header) != set(expected):
    raise SystemExit(f"joint header mismatch: columns={len(header)} unique={len(set(header))}")
print("PASS: joint LeafCutter header contains the exact canonical 846 sample IDs")
PY
