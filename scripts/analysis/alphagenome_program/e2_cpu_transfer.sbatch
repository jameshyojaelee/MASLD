#!/usr/bin/env bash
#SBATCH --job-name=agp_e2-transfer-cpu
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e2-cpu_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e2-cpu_%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT="${ROOT}/GWAS/finemapping/results/alphagenome_program/e2-transfer-20260914T203752Z"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

"${PY}" "${ROOT}/scripts/analysis/alphagenome_program/e2_02_fit_transfer.py" \
  --inputs "${OUT}/inputs" \
  --raw "${OUT}/raw" \
  --output "${OUT}" \
  --resamples 10000 \
  --bootstrap-seed 20260914

"${PY}" "${ROOT}/scripts/analysis/alphagenome_program/e2_03_positive_control.py" \
  --inputs "${OUT}/inputs" \
  --raw "${OUT}/raw" \
  --output "${OUT}" \
  --resamples 2000 \
  --bootstrap-seed 20260914 \
  --noise-levels 1.0 2.0 4.0

"${PY}" -m pip freeze > "${OUT}/validation/pip_freeze_heads.txt"

"${PY}" - "${OUT}" <<'PY'
import hashlib, sys
from pathlib import Path
root = Path(sys.argv[1])
proj = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
inputs = [
    proj / "Analysis/MASLD_Model_Benchmark/executions/gse281364-outcome-blind-splits-21068885/split/elements.tsv",
    proj / "Analysis/MASLD_Model_Benchmark/executions/gse281364-validated-reconstruction-21066470/validated/replicate_outcomes.tsv.gz",
    proj / "Analysis/MASLD_Model_Benchmark/executions/gse281364-borzoi-native-fixture-21083008/fixture/source_group_map.tsv",
    proj / "GWAS/finemapping/results/alphagenome_program/e1-bridge-20260914T194145Z/bridge_variants.tsv",
    proj / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/fold_manifest.tsv",
    proj / "Analysis/MASLD_Model_Benchmark/executions/dna-lm-open-safetensors-21064556/checkpoints/hyenadna.model.safetensors",
    proj / "Analysis/MASLD_Model_Benchmark/executions/gse281364-hyenadna-window-ablation-20260905T192600Z/evaluation/arm_metrics.tsv",
    proj / "Analysis/MASLD_Model_Benchmark/scripts/fit_gse281364_dna_language_seeded_heads.py",
    proj / "Analysis/MASLD_Model_Benchmark/scripts/gse281364_dna_lm_native_contract.py",
    proj / "scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md",
    proj / "scripts/analysis/alphagenome_program/c2_02_extract_embeddings.py",
    proj / "scripts/analysis/alphagenome_program/c2_03_fit_heads.py",
    proj / "scripts/analysis/alphagenome_program/e2_01_build_inputs.py",
    proj / "scripts/analysis/alphagenome_program/e2_02_fit_transfer.py",
    proj / "scripts/analysis/alphagenome_program/e2_03_positive_control.py",
    proj / "scripts/analysis/alphagenome_program/e2_prep.sbatch",
    proj / "scripts/analysis/alphagenome_program/e2_gpu_embeddings.sbatch",
    proj / "scripts/analysis/alphagenome_program/e2_cpu_transfer.sbatch",
    root / "prespec/PRESPECIFICATION.md",
    root / "inputs/e2_units.tsv.gz",
    root / "inputs/window_tokens.npy",
    root / "inputs/window_meta.npz",
    root / "inputs/input_contract.json",
    root / "raw/hyenadna_allele_embeddings.npz",
    root / "raw/hyenadna_allele_embeddings.npz.receipt.json",
]
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

printf 'agp e2 transfer complete\n'
