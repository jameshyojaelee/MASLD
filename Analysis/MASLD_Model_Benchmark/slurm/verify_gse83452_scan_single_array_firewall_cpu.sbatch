#!/usr/bin/env bash
#SBATCH --job-name=model-check-094
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --constraint=v4
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Prove the GSE83452 inductive separation on the GPL16686 code path, before any
# GPL16686 matrix exists.
#
# The GPL570 pass is not evidence here. oligo and SCAN.UPC are different
# packages with different summarization code, and oligo::rma summarizes by
# median polish over the probe-by-array matrix, which borrows across arrays by
# construction. The hazard is assumed present until measured absent, and the
# receipt records gpl570_evidence_carried_over: false so that assumption cannot
# be quietly relaxed.
#
# Pipeline contexts read one CEL per call and must be bit-identical. Multi-array
# estimator contexts are negative controls that must keep diverging. The 1e-12
# perturbation control and the distinct-target-digest check are identical to
# GPL570's.
#
# Raw CEL files are staged onto node-local scratch, checked against the frozen
# raw-source SHA-256 values, and removed on exit.
#
# --constraint=v4 and --cpus-per-task=4 are BOTH reproducibility constraints,
# not resource asks. Do not remove or loosen either one.
#
# This separation was run twice on SHA-256-identical CEL bytes and every digest
# changed. Diagnosed in
# executions/diagnostic-node-numerical-reproducibility-2026-08-25/: the SCAN
# output digest is a pure function of
#   (whether the CPU exposes AVX-512) x (BLAS thread count).
# Nothing else moves it. Node identity does not, and neither does the OS: two
# nodes of the same class agree bitwise, and Ubuntu 22.04 and 24.04 agree bitwise
# within each class. Four cells, four digests:
#   no AVX-512 + 4 threads  4e4a136015...    no AVX-512 + 1 thread  d824f0857b...
#   AVX-512    + 4 threads  e885118fc8...    AVX-512    + 1 thread  bcb4cc4331...
#
# It is NOT a vendor split, and --constraint=intel would NOT fix it. Intel v4
# (Xeon E5-2690 v4; avx, avx2, fma; no AVX-512) gives 4e4a136015, while Intel v5
# (Xeon Gold 6226R; avx512f, avx512dq) gives e885118fc8 -- the same digest as AMD
# v6 EPYC 9254, which also has AVX-512. Intel v5 and AMD v6 agree with each other
# and disagree with Intel v4. The node feature that must be pinned is therefore
# the microarchitecture generation (v3/v4/v5/v6), which is what selects the
# OpenBLAS kernel. Vendor is not the axis.
#
# The cause is OpenBLAS, by two mechanisms with two different fixes:
#   1. Kernel dispatch. OpenBLAS selects a micro-kernel per microarchitecture, so
#      crossprod, solve, chol, qr and svd differ across dispatch classes even
#      single-threaded. Threads CANNOT fix this; only pinning the generation does.
#   2. Threaded reduction order. A multithreaded reduction combines partial sums
#      in a thread-count-dependent order, so gemm and gemv differ by thread count
#      as well. Pinning threads fixes this one.
# SCAN hits both: mybeta and vbeta compute solve(t(Xw)%*%Xw) %*% t(Xw) %*% yw with
# the reduction running over the full probe count.
#
# So neither knob alone is sufficient. Threads alone leaves kernel dispatch; the
# generation pin alone leaves thread count. Setting OPENBLAS_NUM_THREADS=1 does
# NOT make this reproducible: single-threaded AVX-512 and single-threaded
# non-AVX-512 still differ on 840 of 53617 features. Thread count alone, holding
# the node fixed, moves 774 features. cpus-per-task=4 is what pins the thread
# count, because the exports below derive OPENBLAS_NUM_THREADS from
# SLURM_CPUS_PER_TASK. Trimming it as an unnecessary resource ask silently
# reopens this.
#
# Measured, 15 cells: the map (AVX-512 present, thread count) -> digest is exactly
# one-to-one with no exceptions, across Intel E5-2697 v3, Intel E5-2690 v4, Intel
# Gold 6226R v5 and AMD EPYC 9254 v6, and across both Ubuntu generations.
#   no AVX-512: v3, v4     AVX-512: v5, v6
#
# v4 specifically: run 1 of this separation (job 21112166, ne1dc4-074) was v4, and
# 5 of the 8 shards that produced the frozen GSE83452 matrix were v4, so pinning
# v4 preserves continuity rather than choosing arbitrarily. v4 is also the largest
# Intel pool in the cpu partition, so the pin does not starve scheduling. v3
# happens to give the same digest as v4, so --constraint="v3|v4" would also work;
# v4 alone is deliberately the tighter choice, since "v3 has no AVX-512" is a fact
# about today's hardware rather than a guarantee about the feature label.
#
# Scope: this affects the SCAN separation only. The oligo path that actually built
# the frozen 231-array matrix, and the fRMA path behind GSE49541, were both tested
# across dispatch class, glibc and thread count and are bitwise node-stable, so
# they need no such pin -- which is also why those 8 shards could span v4 and v5
# and still agree. The magnitude here is small either way (max abs 1.0e-8 over
# 2.58% of features, r = 1 to machine precision); the pin exists so the separation's
# bit-identity requirements mean what it says, not because the numbers are at risk.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
export R_LIBS_USER=/nonexistent/masld_bench_r_user_library

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark

RUNTIME=${ROOT}/executions/environments/microarray-r-bioc-21104442
RUNTIME_SHA256=8115e0966804dfc4b7bfac3f95016896d099d82fc51f3428fab1ec34182a7d07
R_MODULE=R/4.5.2-mba
PYTHON_MODULE=python/3.11.5-GCCcore-13.2.0

RAW_AUDIT=${ROOT}/executions/model-data-065-21079931
RAW_AUDIT_SHA256=8b30a2265c421f2917ed886d99cfaa08ee8d874abd40e68c8c72d48c726781d1
CEL_TAR=${RAW_AUDIT}/raw/GSE83452_RAW.tar
CEL_MANIFEST=${RAW_AUDIT}/GSE83452_cel_manifest.tsv
CEL_MANIFEST_SHA256=a7c20245ed131b27c56b256c450e9337d5f32d51a54e1b12e14f5e78339988f1

PLATFORM_AUDIT=${ROOT}/executions/model-data-067-21080193-platform
PLATFORM_AUDIT_SHA256=0ce6af4196d2fb0905d6ca3d07fbbeb3bf1ce5b37776a63a726cc1d21794f4bf

PREPROCESSING=${ROOT}/config/evaluation/microarray_transfer_preprocessing.toml
PREPROCESSING_SHA256=11162c381017bd16736433f478c0c82d97726ae2d896957d27efa5aee899c27e
LABEL_SEMANTICS=${ROOT}/config/evaluation/microarray_transfer_label_semantics.toml
LABEL_SEMANTICS_SHA256=b66c20ef5c6bc2ce7597affe7816358bd7fd4ca65df8e3f4be1096dbddaf0513

TARGET=${ROOT}/executions/model-check-094-${SLURM_JOB_ID}-firewall
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
CEL_WORK=${TMPDIR:-/tmp}/masld_gse83452_${SLURM_JOB_ID}

cleanup() {
  rm -rf "${CEL_WORK}"
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed GSE83452 firewall stage: %s\n' "${FAILED}" >&2
  fi
}
trap cleanup EXIT
for path in "${TARGET}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite GSE83452 firewall artifact: %s\n' "${path}" >&2
    exit 1
  fi
done
install -d -m 0750 "${STAGE}" "${STAGE}/validation" "${STAGE}/firewall"
install -d -m 0700 "${CEL_WORK}"

printf '%s  %s\n' \
  "${RUNTIME_SHA256}"         "${RUNTIME}/ARTIFACTS.json" \
  "${RAW_AUDIT_SHA256}"       "${RAW_AUDIT}/ARTIFACTS.json" \
  "${CEL_MANIFEST_SHA256}"    "${CEL_MANIFEST}" \
  "${PLATFORM_AUDIT_SHA256}"  "${PLATFORM_AUDIT}/ARTIFACTS.json" \
  "${PREPROCESSING_SHA256}"   "${PREPROCESSING}" \
  "${LABEL_SEMANTICS_SHA256}" "${LABEL_SEMANTICS}" \
  | sha256sum --check --strict

cd "${ROOT}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

module purge
module load "${PYTHON_MODULE}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
python3 -m unittest -v \
  tests.unit.test_stage_microarray_cel_arrays \
  tests.contract.test_microarray_transfer_contract \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

python3 scripts/stage_microarray_cel_arrays.py \
  --tar "${CEL_TAR}" \
  --manifest "${CEL_MANIFEST}" \
  --output "${CEL_WORK}" \
  --series GSE83452 \
  --expected-records 231 \
  --platform-id GPL16686 \
  --receipt "${STAGE}/firewall/cel_staging_receipt.json" \
  > "${STAGE}/validation/cel_staging_stdout.json" \
  2> "${STAGE}/validation/cel_staging_stderr.txt"

module purge
module load "${R_MODULE}"
R --version > "${STAGE}/validation/R.version.txt"
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/verify_gse83452_scan_single_array_firewall.R" \
  "${CEL_WORK}" \
  "${CEL_MANIFEST}" \
  "${STAGE}/firewall" \
  > "${STAGE}/validation/scan_firewall_stdout.txt" \
  2> "${STAGE}/validation/scan_firewall_stderr.txt"
cp "${STAGE}/firewall/R_sessionInfo.txt" "${STAGE}/validation/R_sessionInfo.txt"

rm -rf "${CEL_WORK}"

# This block records the source identity of the run. Until 2026-08-25 the fourth
# line hashed slurm/verify_gse83452_oligo_single_array_firewall_cpu.sbatch -- the
# OLIGO wrapper -- rather than this file, almost certainly because this wrapper
# was copied from that one and the line was not updated. The four sibling
# microarray wrappers all hash themselves correctly; this was the only one.
#
# The effect was that the SCAN separation's recorded provenance never covered the
# script that actually ran, so any edit to this wrapper was invisible in it.
# Two such edits were made on 2026-08-25 and neither appears in the source.sha256
# of the runs that preceded them: the module purge / PYTHON_MODULE reload fix, and
# the --constraint=v4 reproducibility pin. Runs before this fix
# (model-check-094-21112166-firewall.failed and -21113244-firewall) are preserved
# as they are and are not disturbed; from here on the block hashes this file.
sha256sum \
  scripts/stage_microarray_cel_arrays.py \
  scripts/verify_gse83452_scan_single_array_firewall.R \
  tests/unit/test_stage_microarray_cel_arrays.py \
  slurm/verify_gse83452_scan_single_array_firewall_cpu.sbatch \
  config/evaluation/microarray_transfer_preprocessing.toml \
  config/evaluation/microarray_transfer_label_semantics.toml \
  "${RUNTIME}/ARTIFACTS.json" \
  "${RAW_AUDIT}/ARTIFACTS.json" \
  "${CEL_MANIFEST}" \
  > "${STAGE}/source.sha256"

if find "${STAGE}" -type l -print -quit | grep -q .; then
  printf 'GSE83452 firewall stage contains a symlink\n' >&2
  exit 1
fi
if find "${STAGE}" -type f -name '*.CEL' -print -quit | grep -q .; then
  printf 'GSE83452 firewall stage contains a raw CEL file\n' >&2
  exit 1
fi

module purge
module load "${PYTHON_MODULE}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
stage = Path(sys.argv[1]).resolve(strict=True)
firewall = stage / "firewall"
r = json.loads((firewall / "scan_single_array_firewall.json").read_text(encoding="utf-8"))
if (
    r["status"] != "pass_single_array_isolation_on_real_CEL"
    or r["series"] != "GSE83452" or r["platform_id"] != "GPL16686"
    or r["method"] != "SCAN_with_pd.hugene.2.0.st"
    or r["every_target_pipeline_bitwise_identical"] is not True
    or r["every_perturbation_control_detected"] is not True
    or r["distinct_targets_have_distinct_digests"] is not True
    or r["gpl570_evidence_carried_over"] is not False
    or r["labels_read"] is not False
    or r["model_training_activated"] is not False
):
    raise SystemExit("SCAN firewall receipt differs")
for rec in r["per_target"]:
    if any(int(v) != 0 for v in rec["pipeline_differing_elements"].values()):
        raise SystemExit("a SCAN target differs across pipeline contexts")
metadata = {
    "artifact_class": "gse83452_gpl16686_scan_single_array_firewall",
    "series": "GSE83452", "platform_id": "GPL16686",
    "method": "SCAN_with_pd.hugene.2.0.st",
    "scan_features": r["scan_features"],
    "targets_verified": len(r["per_target"]),
    "held_array_summary_invariant_to_other_arrays_in_the_run": True,
    "perturbation_control_detected_for_every_target": True,
    "multi_array_call_is_possible": r["multi_array_call_is_possible"],
    "scan_pools_across_arrays_when_batched": r["scan_pools_across_arrays_when_batched"],
    "foreach_sequential_backend_registered": True,
    "gpl570_evidence_carried_over": False,
    "labels_read": False, "model_training_activated": False,
    "sealed_outcomes_read": False, "status": "passed",
}
freeze_tree(firewall, metadata)
verify_frozen_tree(firewall)
freeze_tree(stage, metadata | {"artifact_class": "gse83452_scan_firewall_campaign"})
verify_frozen_tree(stage)
PY

trap - EXIT
rm -rf "${CEL_WORK}"
mv "${STAGE}" "${TARGET}"
printf 'GSE83452 GPL16686 single-array firewall: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/firewall/ARTIFACTS.json"
