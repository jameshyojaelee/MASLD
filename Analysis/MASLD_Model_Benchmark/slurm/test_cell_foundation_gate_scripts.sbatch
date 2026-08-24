#!/usr/bin/env bash
#SBATCH --job-name=masld-cellfm-gatetest
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONNOUSERSITE=1

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
micromamba_bin="/gpfs/commons/home/jameslee/.local/bin/micromamba"
python_env="/gpfs/commons/home/jameslee/micromamba/envs/decima"
target="${root}/executions/tests/cell-foundation-gate-scripts-${SLURM_JOB_ID}"
stage="${target}.staging"
if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite gate test: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"

bash -n \
  "${root}/slurm/cell_foundation_safe_checkpoint_inventory.sbatch" \
  "${root}/slurm/geneformer_safe_checkpoint_inventory.sbatch" \
  "${root}/slurm/cell_foundation_weights_only_inspection.sbatch" \
  "${root}/slurm/transcriptformer_safe_archive_inventory.sbatch" \
  "${root}/slurm/transcriptformer_safe_materialize.sbatch" \
  "${root}/slurm/cell_foundation_existing_runtime_audit.sbatch" \
  "${root}/slurm/build_cell_foundation_admission_ledger.sbatch"

python3 -m py_compile \
  "${root}/scripts/build_cell_foundation_admission_ledger.py" \
  "${root}/scripts/audit_existing_cell_foundation_runtime.py" \
  "${root}/scripts/safe_tar_materialize.py" \
  "${root}/scripts/torch_checkpoint_safe_inventory.py" \
  "${root}/scripts/safetensors_safe_inventory.py"

python3 "${root}/scripts/build_cell_foundation_admission_ledger.py" \
  --root "${root}" --output "${stage}/admission_ledger_test.json"
python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); assert len(d["records"]) == 16; assert d["sealed_or_test_outcomes_opened"] is False; assert d["champion_claims_allowed"] is False; assert all(r["promotion_status"] == "not_promoted" for r in d["records"])' \
  "${stage}/admission_ledger_test.json"

"${micromamba_bin}" run --prefix "${python_env}" python -c \
  'import sys,torch; torch.save({"state_dict":{"encoder.weight":torch.arange(12,dtype=torch.float32).reshape(3,4)}},sys.argv[1])' \
  "${stage}/synthetic_checkpoint.pt"
synthetic_sha256="$(sha256sum "${stage}/synthetic_checkpoint.pt" | cut -d' ' -f1)"
synthetic_size="$(stat -c '%s' "${stage}/synthetic_checkpoint.pt")"
"${micromamba_bin}" run --prefix "${python_env}" \
  python "${root}/scripts/torch_checkpoint_weights_only_inspect.py" \
  --checkpoint "${stage}/synthetic_checkpoint.pt" \
  --output "${stage}/synthetic_weights_only_inspection.json" \
  --expected-sha256 "${synthetic_sha256}" --expected-size "${synthetic_size}"
python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d["status"] == "pass"; assert d["weights_only"] is True; assert d["tensor_count"] == 1; assert d["tensor_numel"] == 12' \
  "${stage}/synthetic_weights_only_inspection.json"

python3 --version > "${stage}/python_version.txt" 2>&1
"${micromamba_bin}" run --prefix "${python_env}" python -m pip freeze \
  > "${stage}/pip_freeze.txt"
sha256sum \
  "${root}/slurm/test_cell_foundation_gate_scripts.sbatch" \
  "${root}/scripts/build_cell_foundation_admission_ledger.py" \
  "${root}/scripts/audit_existing_cell_foundation_runtime.py" \
  "${root}/scripts/torch_checkpoint_weights_only_inspect.py" \
  > "${stage}/source_manifest.sha256"
(
  cd "${stage}"
  sha256sum admission_ledger_test.json synthetic_checkpoint.pt \
    synthetic_weights_only_inspection.json python_version.txt pip_freeze.txt \
    > SHA256SUMS
)
chmod 0440 "${stage}"/*
mv "${stage}" "${target}"
printf 'Cell-foundation gate script tests: %s\n' "${target}"
