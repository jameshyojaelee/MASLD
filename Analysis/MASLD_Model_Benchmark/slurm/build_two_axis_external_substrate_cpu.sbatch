#!/usr/bin/env bash
#SBATCH --job-name=model-data-953
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 1, job 1 of 3: reproduce the Stage 0c memberships and build both arms.
#
# No outcome-versus-expression association is computed in either external arm.
# This job fixes what the assigned sets are and what survives the
# cross-namespace join, so the prespecification can be held back against a known
# substrate without any external answer existing yet.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/build_two_axis_external_substrate.py
ARM_A=${ROOT}/scripts/evaluate_gse267145_axis_count.py
ARM_B=${ROOT}/executions/model-data-880-21130257-gse135251-source/source
ACTIVATION=${ROOT}/config/evaluation/two_axis_activity_fibrosis_external_activation.json
EXTERNAL=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk

TARGET=${ROOT}/executions/model-data-953-${SLURM_JOB_ID}-two-axis-substrate
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the substrate build.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/build_two_axis_external_substrate_cpu.sbatch

# The gitignored activation record is pinned by digest, since git status cannot
# see it drift.
cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
3c5d3cc4e3f47b331c5f8ece04c213c304d49faa3ce3b974cef878c57e4f6ca4  config/evaluation/two_axis_activity_fibrosis_external_activation.json
5f26597ab6b6aeeaff867d54112554a33d4d0dc35be5f333e17e9010a7199fed  executions/model-check-1050-21182309-arm-b-axis-count/arm_b_axis_count/ARTIFACTS.json
9e7bdda3d4511a9bc0e520515f9ab12f4ca0729ca9a74a33f2c0491a4f253db7  executions/model-data-880-21130257-gse135251-source/source/outcomes/participant_endpoints.tsv
d6219f67fa676e888f47b0225a75a18036ddf83ba9bbd2f1991069e3781e39aa  executions/model-data-880-21130257-gse135251-source/source/molecular/rna_values.npy
3dcdd51055010215592f9ff7372c62824b7345565d1bda4428a29536b5e5c3a8  executions/model-data-880-21130257-gse135251-source/source/molecular/rna_feature_axis.tsv
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"
cat > "${STAGE}/external_inputs.sha256" <<'SUMS'
0d9aa05ce52b7c16501a947e12e25a9b0513a83a2a3a26751307de979a500334  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE130970/GSE130970_counts.tsv
a1b06183e9c9f4444db264afcf2fd18d0215a3d68a7163fd3bc0d3851acd3667  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE193066/GSE193066_counts.tsv
8e8e0661c03a431f7006cf06c1e9bfdf62e4030392c505a4fd93e5cda1a16cba  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE130970/GSE130970_metadata.tsv
667782f577e8ad7cc3d80d5b776d02a4e32c20e8fe54a964ad3102069e1ed7d9  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE193066/GSE193066_metadata.tsv
SUMS
sha256sum --check --strict "${STAGE}/external_inputs.sha256"

"${PYTHON}" -m unittest -v tests.unit.test_build_two_axis_external_substrate \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" "${SCRIPT}" \
  --analysis-script "${ARM_A}" \
  --arm-b-source "${ARM_B}" \
  --activation "${ACTIVATION}" \
  --external-root "${EXTERNAL}" \
  --output "${STAGE}/substrate" \
  --seed 20260827 \
  > "${STAGE}/build_stdout.json" 2> "${STAGE}/build_stderr.txt"

sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"
sha256sum --check --strict "${STAGE}/external_inputs.sha256"

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" "${ARM_A}" \
  slurm/build_two_axis_external_substrate_cpu.sbatch \
  tests/unit/test_build_two_axis_external_substrate.py \
  "${ACTIVATION}" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
payload = json.loads((root / "substrate" / "substrate.json").read_text(encoding="utf-8"))
if payload["no_outcome_association_was_computed_in_either_external_arm"] is not True:
    raise SystemExit("an external outcome association was computed before the seal")
for entry in payload["reproduction_of_stage_0c"]["invariants_reproduced"]:
    if entry["matches_stage_0c"] is not True:
        raise SystemExit("Stage 0c did not reproduce")
cells = payload["assignment"]["cell_sizes"]
marg = payload["assignment"]["marginals_reconstructed"]
if marg["activity_only_plus_both"] != 4388 or marg["fibrosis_only_plus_both"] != 1305:
    raise SystemExit("the assignment does not reconstruct Stage 0c's marginals")
if sum(cells.values()) != payload["reproduction_of_stage_0c"]["realised_universe"]:
    raise SystemExit("the four cells do not partition the universe")
if payload["activation_record"]["user_approved_scope_only"] is not False:
    raise SystemExit("the activation's user_approved flag was not carried as recorded")
verify_frozen_tree(root / "substrate")
freeze_tree(root, {
    "artifact_class": "two_axis_external_substrate_campaign",
    "stage": "stage_1", "job": "1_of_3",
    "stage_0c_reproduced_exactly": True,
    "outcome_association_computed_externally": False,
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({"cell_sizes": cells}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/substrate/ARTIFACTS.json"
