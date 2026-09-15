#!/usr/bin/env bash
#SBATCH --job-name=model-data-954
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 3, job 1 of 3: join composition to each arm and derive the eligible family.
# No outcome-versus-expression association is computed. The family is derived on
# the two kallisto arms only; GSE193066 is STAR and its detectability is
# reported separately rather than folded into the same rule.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/build_composition_substrate.py
EXTENSION=${ROOT}/scripts/multicovariate_partial.py
SEALED=${ROOT}/scripts/evaluate_gse267145_axis_count.py
PROPS=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/results
CROSSWALK=${ROOT}/executions/model-data-950-21181338-showcase-substrate/substrate/arm_b/participant_axis.tsv
ARM_B=${ROOT}/executions/model-data-880-21130257-gse135251-source/source
STAGE1=${ROOT}/executions/model-data-953-21182514-two-axis-substrate/substrate

TARGET=${ROOT}/executions/model-data-954-${SLURM_JOB_ID}-composition-substrate
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the composition substrate.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}" "${STAGE}/validation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/build_composition_substrate_cpu.sbatch

# The held-back instrument must be untouched: six held-back jobs hashed it.
cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
5f4e483dc22f1143f1c4447f60e308a5822dd86070540527abbd429e7278e559  scripts/evaluate_gse267145_axis_count.py
26268ea1d5e822a8d7824f2286aca04380f773aacea8e0a03608325cd48dee2f  executions/model-data-950-21181338-showcase-substrate/substrate/arm_b/ARTIFACTS.json
68837a8293d5386e6bd7b539dd9107f61d1cfddb783feec373d9853e57be84ab  executions/model-data-953-21182514-two-axis-substrate/substrate/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" -m unittest -v tests.unit.test_multicovariate_partial \
  tests.unit.test_build_composition_substrate \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" "${SCRIPT}" \
  --extension "${EXTENSION}" --sealed-instrument "${SEALED}" \
  --proportions-root "${PROPS}" --crosswalk "${CROSSWALK}" \
  --arm-b-source "${ARM_B}" --stage1-substrate "${STAGE1}" \
  --output "${STAGE}/substrate" \
  > "${STAGE}/build_stdout.json" 2> "${STAGE}/build_stderr.txt"

sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" "${EXTENSION}" "${SEALED}" \
  slurm/build_composition_substrate_cpu.sbatch \
  tests/unit/test_build_composition_substrate.py \
  tests/unit/test_multicovariate_partial.py > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
p = json.loads((root / "substrate" / "composition_substrate.json").read_text(encoding="utf-8"))
if p["no_outcome_association_was_computed"] is not True:
    raise SystemExit("an outcome association was computed before the seal")
fam = p["eligible_family"]["family"]
if len(fam) != 6:
    raise SystemExit(f"the family is {len(fam)}, not the approved six")
if p["eligible_family"]["derived_on"] != ["GSE135251", "GSE130970"]:
    raise SystemExit("the family was not derived on the two kallisto arms")
pooled = p["families_recorded_side_by_side"]["pooled_confounded_rule"]["family"]
if set(pooled) >= set(fam):
    raise SystemExit("the pooled rule did not lose lineages; the confound claim fails")
if p["quantification_is_not_uniform_across_the_arms"]["per_cohort"]["GSE193066"] != "STAR":
    raise SystemExit("GSE193066 was not recorded as STAR")
for cohort, arm in p["arms"].items():
    if arm["joined"] != arm["axis_members"]:
        raise SystemExit(f"{cohort} did not join every axis member")
    c = arm["family_conditioning_without_the_other_axis"]
    if c["is_rank_deficient"]:
        raise SystemExit(f"{cohort} family block is rank deficient")
if p["crosswalk"]["one_to_one_both_directions_asserted"] is not True:
    raise SystemExit("the crosswalk 1:1 property was not asserted")
verify_frozen_tree(root / "substrate")
freeze_tree(root, {"artifact_class": "composition_substrate_campaign",
                   "stage": "stage_3", "job": "1_of_3",
                   "family_size": len(fam),
                   "quantification_uniform_across_arms": False,
                   "status": "passed"})
verify_frozen_tree(root)
print(json.dumps({"family": fam, "pooled_would_give": pooled}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/substrate/ARTIFACTS.json"
