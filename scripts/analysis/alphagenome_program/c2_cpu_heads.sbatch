#!/usr/bin/env bash
#SBATCH --job-name=agp_c2-endogenous-head-cpu
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c2-cpu_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c2-cpu_%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT="${ROOT}/GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

"${PY}" "${ROOT}/scripts/analysis/alphagenome_program/c2_03_fit_heads.py" \
  --inputs "${OUT}/inputs" \
  --raw "${OUT}/raw" \
  --output "${OUT}" \
  --resamples 10000 \
  --bootstrap-seed 20260914

"${PY}" "${ROOT}/scripts/analysis/alphagenome_program/c2_04_label_precision.py" \
  --labels "${OUT}/inputs/currin_lead_labels.tsv.gz" \
  --output "${OUT}/tables/label_precision_posthoc.json"

"${PY}" -m pip freeze > "${OUT}/validation/pip_freeze_heads.txt"
"${PY}" - "${OUT}" <<'PY'
import hashlib, json, sys
from pathlib import Path
root = Path(sys.argv[1])
inputs = [
    Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1/liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz"),
    Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/fold_manifest.tsv"),
    Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/dna-lm-open-safetensors-21064556/checkpoints/hyenadna.model.safetensors"),
    Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/dna-lm-open-safetensors-21064556/checkpoints/caduceus.model.safetensors"),
    Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/scripts/fit_gse281364_dna_language_seeded_heads.py"),
    Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/scripts/gse281364_dna_lm_native_contract.py"),
    Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-train-605-21099008/evaluation/model_head_summary.tsv"),
    root / "prespec/PRESPECIFICATION.md",
    root / "inputs/currin_lead_labels.tsv.gz",
    root / "inputs/window_tokens.npy",
    root / "raw/hyenadna_allele_embeddings.npz",
    root / "raw/caduceus_allele_embeddings.npz",
]
for path in sorted(Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program").glob("c2_*")):
    inputs.append(path)
rows = ["path\tbytes\tsha256"]
for path in inputs:
    if not path.exists():
        rows.append(f"{path}\tMISSING\tMISSING")
        continue
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 << 20), b""):
            digest.update(block)
    rows.append(f"{path}\t{path.stat().st_size}\t{digest.hexdigest()}")
(root / "MANIFEST.tsv").write_text("\n".join(rows) + "\n")
print("\n".join(rows))
PY
printf 'agp c2 heads complete\n'
