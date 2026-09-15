#!/usr/bin/env bash
#SBATCH --job-name=model-probe-021
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# How far do TranscriptFormer's own outputs move between torch 2.5.1+cu124 and
# torch 2.8.0+cu128 when the GPU pool is held fixed?
#
# L40S deliberately. cu124 cannot initialise on B6K at all, so an L40S node is
# the only place both runtimes can run and therefore the only place the pool can
# be held fixed while the runtime varies.
#
# Thresholds were frozen BEFORE this job was written, in
# config/evaluation/transcriptformer_runtime_version_effect_prespecification_20260830.json
# sha256 f1798a2529949a7f18d714281198d4cf1e4fa4b3aa41025c367bbd9c27215780
# and that digest is checked below, so the job cannot run against moved goalposts.
#
# This reuses the extraction the supported run already used, unchanged and bound
# by its recorded digest, rather than a fresh invocation. The runtime is the only
# thing that differs between the two arms: same script, same checkpoint, same
# fixture, same offline guard, same deterministic cuBLAS setting, same node, same
# GPU, one job.
#
# ⛔ PROBE ONLY. The embeddings this writes are NOT the production re-extraction
# and may not enter a ranking. The tournament re-run would have to be on B6K;
# this is on L40S precisely because it answers a different question.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark

SUPPORTED_RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
SUPPORTED_RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PROBE_RUNTIME=${ROOT}/executions/environments/transcriptformer-cu128-UNSUPPORTED-PROBE-21247291
PROBE_RUNTIME_SHA256=b87ddc4339a252a3f4be3c7c05436dc7e4403099c59e0cc8c9ff32728c457682

CHECKPOINT=${ROOT}/executions/transcriptformer-materialized-21066178/bundle/tf_sapiens
FIXTURE=${ROOT}/executions/model-data-051-21077382
FIXTURE_ARTIFACTS_SHA256=ba40ce9555bb6d7bc278be2a80015b993f9f3650fa06f59505ad279c88cc117a
EXTRACT_SCRIPT_SHA256=cdcc828a24c7176e08b002a59bb80fde1fe89e29d7ba6fadc00cc58a1e3b1a48
OFFLINE_GUARD_SHA256=ed4bd726db8210d11c034b806c54fa541afce62664ee0846912022036e414635
PRESPEC=${ROOT}/config/evaluation/transcriptformer_runtime_version_effect_prespecification_20260830.json
PRESPEC_SHA256=f1798a2529949a7f18d714281198d4cf1e4fa4b3aa41025c367bbd9c27215780
REFERENCE_EMBEDDINGS=${ROOT}/executions/model-training-500-21080428/embeddings/common_embeddings.npz

TARGET=${ROOT}/executions/transcriptformer-runtime-version-effect-PROBE-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/validation" \
  "${STAGE}/cu124/cache" \
  "${STAGE}/cu128/cache"

cd "${ROOT}"
sha256sum --check --strict <<EOF
${SUPPORTED_RUNTIME_SHA256}  ${SUPPORTED_RUNTIME}/ARTIFACTS.json
${PROBE_RUNTIME_SHA256}  ${PROBE_RUNTIME}/SHA256SUMS
${FIXTURE_ARTIFACTS_SHA256}  ${FIXTURE}/ARTIFACTS.json
${EXTRACT_SCRIPT_SHA256}  ${ROOT}/scripts/transcriptformer_extract_unlabeled_embeddings.py
${OFFLINE_GUARD_SHA256}  ${ROOT}/scripts/offline_guard/sitecustomize.py
${PRESPEC_SHA256}  ${PRESPEC}
EOF

nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv \
  > "${STAGE}/validation/nvidia-smi.csv"

export PYTHONPATH="${ROOT}/scripts/offline_guard:${ROOT}:${ROOT}/src"

# One extraction, run twice. The only argument that differs is the interpreter
# and the output directory; the caches are separated so neither arm can read a
# compiled kernel the other produced.
extract () {   # extract <arm> <python>
  local arm=$1 python_bin=$2
  export XDG_CACHE_HOME="${STAGE}/${arm}/cache"
  export MPLCONFIGDIR="${STAGE}/${arm}/cache/matplotlib"
  export TORCHINDUCTOR_CACHE_DIR="${STAGE}/${arm}/cache/torchinductor"
  export TRITON_CACHE_DIR="${STAGE}/${arm}/cache/triton"
  "${python_bin}" -c "import json,torch;print(json.dumps({'arm':'${arm}','torch':torch.__version__,'cuda':torch.version.cuda,'device':torch.cuda.get_device_name(0),'capability':list(torch.cuda.get_device_capability(0))}))" \
    > "${STAGE}/validation/${arm}-runtime.json"
  "${python_bin}" scripts/transcriptformer_extract_unlabeled_embeddings.py \
    --checkpoint-root "${CHECKPOINT}" \
    --fixture-root "${FIXTURE}/fixture" \
    --output "${STAGE}/${arm}/embeddings" \
    --expected-rows 50000 \
    > "${STAGE}/validation/${arm}-embedding-stdout.txt" \
    2> "${STAGE}/validation/${arm}-embedding-stderr.txt"
  "${python_bin}" -m pip freeze --all > "${STAGE}/validation/${arm}-pip-freeze.txt"
}

extract cu124 "${SUPPORTED_RUNTIME}/env/bin/python"
extract cu128 "${PROBE_RUNTIME}/env/bin/python"

# Compare. Thresholds come from the frozen record rather than from this file, so
# there is one place they can be read and no way to restate them differently here.
"${SUPPORTED_RUNTIME}/env/bin/python" - \
  "${STAGE}" "${PRESPEC}" "${REFERENCE_EMBEDDINGS}" \
  <<'PY' > "${STAGE}/validation/embedding-comparison-stdout.json"
import json
import sys
from pathlib import Path

import numpy as np

stage, prespec_path, reference_path = (Path(a) for a in sys.argv[1:4])
prespec = json.loads(prespec_path.read_text(encoding="utf-8"))
secondary = prespec["secondary_criterion"]
rms = float(secondary["reference_scale_measured_first"]["rms"])

a = np.load(stage / "cu124/embeddings/common_embeddings.npz")["embeddings"]
b = np.load(stage / "cu128/embeddings/common_embeddings.npz")["embeddings"]
if a.shape != b.shape:
    raise SystemExit(f"embedding shapes differ: {a.shape} vs {b.shape}")

a64 = a.astype(np.float64)
b64 = b.astype(np.float64)
delta = np.abs(a64 - b64)
max_abs = float(delta.max())
bitwise_identical = bool(max_abs == 0.0)

norm_a = np.linalg.norm(a64, axis=1)
norm_b = np.linalg.norm(b64, axis=1)
denom = np.maximum(norm_a * norm_b, 1e-300)
cosine = (a64 * b64).sum(axis=1) / denom

# Also compare the cu124 arm against the already-published reference embeddings.
# Those were produced by the same runtime on an L40S, so any difference here is
# run-to-run variation rather than a version effect, and it is the floor against
# which the version effect has to be read.
reference = np.load(reference_path)["embeddings"].astype(np.float64)
run_to_run_max_abs = (
    float(np.abs(a64 - reference).max()) if reference.shape == a64.shape else None
)

behaving = max_abs <= 1e-3 * rms and float(np.median(cosine)) >= 0.9999
fragile = max_abs > 1e-2 * rms or float(np.median(cosine)) < 0.999

result = {
    "arms": {"a": "torch 2.5.1+cu124 (supported)", "b": "torch 2.8.0+cu128 (unsupported probe)"},
    "pool_held_fixed": "one L40S node, one GPU, one job",
    "rows": int(a.shape[0]),
    "width": int(a.shape[1]),
    "reference_rms_from_prespecification": rms,
    "bitwise_identical": bitwise_identical,
    "max_absolute_difference": max_abs,
    "max_absolute_difference_as_fraction_of_rms": max_abs / rms,
    "mean_absolute_difference": float(delta.mean()),
    "per_row_cosine_min": float(cosine.min()),
    "per_row_cosine_median": float(np.median(cosine)),
    "per_row_cosine_p01": float(np.percentile(cosine, 1)),
    "same_runtime_run_to_run_max_absolute_difference": run_to_run_max_abs,
    "secondary_threshold_behaving": 1e-3 * rms,
    "secondary_threshold_fragile": 1e-2 * rms,
    "secondary_verdict": (
        "behaving" if behaving else "fragile" if fragile else "inconclusive"
    ),
    "secondary_is_descriptive_only": True,
    "primary_verdict": "not_available_until_heads_are_fit",
    "primary_note": (
        "A refit head can absorb a rotation, a global scale change, or a "
        "coordinate permutation, so these embedding numbers do not decide the "
        "question. The donor-class-balanced macro-F1 comparison does."
    ),
    "prespecification_sha256_verified": True,
}
(stage / "validation/embedding_comparison.json").write_text(
    json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(result, indent=2, sort_keys=True))
PY

cat > "${STAGE}/PROBE_ONLY.json" <<EOF
{
  "schema_version": "masld-bench-runtime-version-effect-probe-v1",
  "status": "runtime_version_effect_probe",
  "production_use_permitted": false,
  "may_produce_benchmark_metrics": false,
  "may_enter_a_ranking": false,
  "purpose": "Measure the torch 2.5.1 to 2.8.0 effect on TranscriptFormer with the GPU pool held fixed on L40S.",
  "not_the_production_re_extraction": "The tournament re-run would have to be on B6K. This is on L40S deliberately, to hold the pool fixed while the runtime varies.",
  "prespecification": "config/evaluation/transcriptformer_runtime_version_effect_prespecification_20260830.json",
  "prespecification_sha256": "${PRESPEC_SHA256}",
  "unsupported_environment": "${PROBE_RUNTIME}",
  "authorised_by": "team lead, 2026-08-30, probe only"
}
EOF

# The extraction script writes its own ARTIFACTS.json into the output directory
# (scripts/transcriptformer_extract_unlabeled_embeddings.py:370), so the tree is
# already sealed by the producer. Freezing it again here would collide with that
# record rather than add to it. Record the digests instead: these are exactly the
# EMBEDDINGS_ARTIFACTS_SHA256 values the head bundle would consume if the metric
# stage is authorised.
for arm in cu124 cu128; do
  sha256sum "${STAGE}/${arm}/embeddings/ARTIFACTS.json"
done > "${STAGE}/validation/embeddings-artifacts.sha256"
cat "${STAGE}/validation/embeddings-artifacts.sha256"

sha256sum \
  scripts/transcriptformer_extract_unlabeled_embeddings.py \
  scripts/offline_guard/sitecustomize.py \
  slurm/probe_transcriptformer_runtime_version_effect_l40s.sbatch \
  "${PRESPEC}" \
  "${SUPPORTED_RUNTIME}/ARTIFACTS.json" \
  "${PROBE_RUNTIME}/SHA256SUMS" \
  "${FIXTURE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'runtime version-effect probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/cu124/embeddings/ARTIFACTS.json" "${TARGET}/cu128/embeddings/ARTIFACTS.json"
cat "${TARGET}/validation/embedding_comparison.json"
