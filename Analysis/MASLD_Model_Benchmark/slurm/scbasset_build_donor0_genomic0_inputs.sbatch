#!/usr/bin/env bash
#SBATCH --job-name=masld-scbasset-inputs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=5
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONHASHSEED=20260824
: "${OUTER_FOLD:=0}"
case "${OUTER_FOLD}" in 0|1|2|3|4) ;; *) printf 'Invalid OUTER_FOLD\n' >&2; exit 2 ;; esac
GENOMIC_VALID_FOLD=$(( (OUTER_FOLD + 1) % 5 ))
SPLIT_ID="donor${OUTER_FOLD}_genomic${OUTER_FOLD}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PSEUDOBULK="${ROOT}/executions/gse296875-training-fold-pseudobulk-21063993"
MEMBERSHIP="${ROOT}/executions/gse296875-fragment-membership-21063829"
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
REFERENCE="${ROOT}/executions/sequence-reference-contract-21063670"
FASTA_ROOT="${ROOT}/executions/sequence-fasta-21063808"
FASTA="${FASTA_ROOT}/GRCh38.p14.sequence_model.fa"
BLACKLIST_ROOT="${ROOT}/executions/encode-grch38-blacklist-21063857"
BLACKLIST_GZ="${BLACKLIST_ROOT}/ENCFF356LFX.bed.gz"
MACS_PROBE="${ROOT}/executions/macs2-runtime-probe-21063873"
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
TARGET="${ROOT}/executions/scbasset-${SPLIT_ID}-inputs-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PSEUDOBULK_ARTIFACTS_SHA256=3190aa76cd5872f50d7b0c71bcb92520f1f2428110ce2fb4c26856994bc17e5b
SPLIT_ARTIFACTS_SHA256=00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f
PYTHONPATH_CORE="${ROOT}:${CORE}/site-packages"

preserve_failure() {
  status=$?
  trap - EXIT
  if (( status != 0 )) && [[ -d "${STAGE}" ]]; then
    failed="${TARGET}.failed"
    if [[ -e "${failed}" ]]; then
      printf 'Refusing to overwrite preserved failed staging: %s\n' "${failed}" >&2
    else
      mv "${STAGE}" "${failed}"
      printf 'Preserved failed staging: %s\n' "${failed}" >&2
    fi
  fi
  exit "${status}"
}
trap preserve_failure EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite scBasset input campaign: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/tagalign" \
  "${STAGE}/macs2" \
  "${STAGE}/regions" \
  "${STAGE}/validation"

printf '%s  %s\n' "${PSEUDOBULK_ARTIFACTS_SHA256}" "${PSEUDOBULK}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' 582387b78ec25215d45820c24aace12f130b826e91755bf19ac4581730611e64 \
  "${MEMBERSHIP}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${SPLIT_ARTIFACTS_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' 7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d \
  "${REFERENCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 \
  "${FASTA_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' ed20dceeb315cffa2bc3de050d846eceefc552451fef741f0dd95b0149d10fce \
  "${FASTA}" | sha256sum --check --strict
printf '%s  %s\n' fdbc0d2fbfb23781349a6ba880e29cf4e921cb30857075d5a91395f164c5a540 \
  "${BLACKLIST_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 5c5fdf646e72c13f099380d92348c150850083b33bc40659e0f92a299d211760 \
  "${MACS_PROBE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${PSEUDOBULK}" "${MEMBERSHIP}" "${SPLIT}" "${REFERENCE}" \
  "${FASTA_ROOT}" "${BLACKLIST_ROOT}" "${MACS_PROBE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("scBasset input frozen-tree verification passed")
PY

module purge
module load singularity/4.2.2
cd "${ROOT}"
lineages=(cholangiocyte fibroblast hepatocyte macrophage t_cell)
pids=()
for lineage in "${lineages[@]}"; do
  singularity exec --cleanenv \
    --bind "${ROOT}:${ROOT}:rw" \
    --env "PYTHONPATH=${PYTHONPATH_CORE}" \
    --env PYTHONNOUSERSITE=1 \
    --env PYTHONDONTWRITEBYTECODE=1 \
    "${IMAGE}" python scripts/emit_training_tagalign.py \
      --pseudobulk "${PSEUDOBULK}" \
      --pseudobulk-artifacts-sha256 "${PSEUDOBULK_ARTIFACTS_SHA256}" \
      --split-contract "${SPLIT}" \
      --split-artifacts-sha256 "${SPLIT_ARTIFACTS_SHA256}" \
      --donor-test-fold "${OUTER_FOLD}" \
      --genomic-test-fold "${OUTER_FOLD}" \
      --genomic-valid-fold "${GENOMIC_VALID_FOLD}" \
      --lineage "${lineage}" \
      --output-directory "${STAGE}/tagalign/${lineage}" \
      --summary "${STAGE}/tagalign/${lineage}.summary.json" \
      > "${STAGE}/validation/tagalign.${lineage}.stdout.txt" &
  pids+=("$!")
done
tagalign_status=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then
    tagalign_status=1
  fi
done
if (( tagalign_status != 0 )); then
  printf 'At least one scBasset tagAlign emitter failed\n' >&2
  exit 1
fi
gzip -t "${STAGE}"/tagalign/*/genomic_fold_*.tagAlign.gz

module purge
module load MACS2/2.2.9.1-foss-2023b
macs2 --version > "${STAGE}/validation/macs2.version.txt"
python - "${SPLIT}/genomic_folds.tsv" "${STAGE}/validation/fold_gsizes.tsv" \
  "${OUTER_FOLD}" "${GENOMIC_VALID_FOLD}" <<'PY'
from pathlib import Path
import csv
import sys

with Path(sys.argv[1]).open(encoding="utf-8", newline="") as handle:
    sizes = {int(row["genomic_fold"]): int(row["total_bp"]) for row in csv.DictReader(handle, delimiter="\t")}
held = {int(sys.argv[3]), int(sys.argv[4])}
with Path(sys.argv[2]).open("x", encoding="utf-8") as handle:
    for fold in sorted(set(sizes).difference(held)):
        handle.write(f"{fold}\t{sizes[fold]}\n")
PY
pids=()
while IFS=$'\t' read -r genomic_fold genome_size; do
  treatments=()
  for lineage in "${lineages[@]}"; do
    treatments+=("${STAGE}/tagalign/${lineage}/genomic_fold_${genomic_fold}.tagAlign.gz")
  done
  macs2 callpeak \
    --treatment "${treatments[@]}" \
    --format BED \
    --gsize "${genome_size}" \
    --name "genomic_fold_${genomic_fold}" \
    --outdir "${STAGE}/macs2" \
    --nomodel \
    --shift -75 \
    --extsize 150 \
    --keep-dup all \
    --qvalue 0.01 \
    --call-summits \
    > "${STAGE}/validation/macs2.genomic_fold_${genomic_fold}.stdout.txt" \
    2> "${STAGE}/validation/macs2.genomic_fold_${genomic_fold}.stderr.txt" &
  pids+=("$!")
done < "${STAGE}/validation/fold_gsizes.tsv"
macs_status=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then
    macs_status=1
  fi
done
if (( macs_status != 0 )); then
  printf 'At least one pooled-five-lineage MACS2 call failed\n' >&2
  exit 1
fi
python - "${STAGE}/macs2" "${STAGE}/regions/macs2.training.narrowPeak" \
  "${STAGE}/validation/fold_gsizes.tsv" <<'PY'
from pathlib import Path
import sys

source = Path(sys.argv[1])
target = Path(sys.argv[2])
with Path(sys.argv[3]).open(encoding="utf-8") as handle:
    folds = [int(line.split("\t", 1)[0]) for line in handle if line.strip()]
rows = 0
with target.open("x", encoding="utf-8") as writer:
    for fold in folds:
        path = source / f"genomic_fold_{fold}_peaks.narrowPeak"
        with path.open(encoding="utf-8") as reader:
            for line in reader:
                writer.write(line)
                rows += 1
if rows == 0:
    raise SystemExit("pooled-five-lineage MACS2 output is empty")
print(rows)
PY

module purge
module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_CORE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python "${ROOT}/scripts/filter_training_peaks.py" \
    --narrowpeak "${STAGE}/regions/macs2.training.narrowPeak" \
    --blacklist "${BLACKLIST_GZ}" \
    --fasta "${FASTA}" \
    --split-contract "${SPLIT}" \
    --genomic-test-fold "${OUTER_FOLD}" \
    --genomic-valid-fold "${GENOMIC_VALID_FOLD}" \
    --output "${STAGE}/regions/training.filtered.narrowPeak" \
    --summary "${STAGE}/regions/training.filtered.summary.json" \
    > "${STAGE}/validation/peak-filter.stdout.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_CORE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python "${ROOT}/scripts/assemble_chrombpnet_fold_regions.py" \
    --training-peaks "${STAGE}/regions/training.filtered.narrowPeak" \
    --split-contract "${SPLIT}" \
    --split-id "${SPLIT_ID}" \
    --output-peaks "${STAGE}/regions/fold.regions.bed" \
    --output-fold "${STAGE}/regions/fold.json" \
    --output-summary "${STAGE}/regions/fold.summary.json" \
    > "${STAGE}/validation/region-assembly.stdout.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "${IMAGE}" python "${ROOT}/scripts/scbasset_build_fold_inputs.py" \
    --pseudobulk "${PSEUDOBULK}" \
    --membership "${MEMBERSHIP}/cell_membership.tsv.gz" \
    --donor-folds "${MEMBERSHIP}/donor_folds.tsv" \
    --regions-bed "${STAGE}/regions/fold.regions.bed" \
    --fold-json "${STAGE}/regions/fold.json" \
    --fasta "${FASTA}" \
    --donor-test-fold "${OUTER_FOLD}" \
    --donor-valid-fold "${GENOMIC_VALID_FOLD}" \
    --minimum-cell-fraction 0.05 \
    --workers 5 \
    --output "${STAGE}/inputs" \
    > "${STAGE}/validation/fold-inputs.stdout.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:ro" \
  --env "PYTHONPATH=${ROOT}/src:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"

sha256sum \
  "${ROOT}/scripts/emit_training_tagalign.py" \
  "${ROOT}/scripts/filter_training_peaks.py" \
  "${ROOT}/scripts/assemble_chrombpnet_fold_regions.py" \
  "${ROOT}/scripts/scbasset_build_fold_inputs.py" \
  "${ROOT}/slurm/scbasset_build_donor0_genomic0_inputs.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${PSEUDOBULK_ARTIFACTS_SHA256}" \
  "${OUTER_FOLD}" "${GENOMIC_VALID_FOLD}" "${SPLIT_ID}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "inputs" / "summary.json").read_text())
outer_fold = int(sys.argv[3])
valid_fold = int(sys.argv[4])
split_id = sys.argv[5]
train_folds = sorted(set(range(5)).difference({outer_fold, valid_fold}))
if (
    summary["status"] != "pass"
    or summary["split_id"] != split_id
    or summary["donor_train_folds"] != train_folds
    or summary["donor_valid_fold"] != valid_fold
    or summary["donor_test_fold"] != outer_fold
    or summary["training_cells"] < 1000
    or summary["validation_ccres"] != 16000
    or summary["test_ccres"] != 16000
    or summary["held_donor_atac_read"]
    or not summary["genomic_test_atac_parsed_only_to_exclude"]
    or summary["genomic_test_atac_used_for_matrix"]
    or summary["genomic_test_atac_exported"]
    or summary["model_training_input_contains_genomic_test_atac"]
    or not summary["intermediate_lineage_matrices_removed"]
    or summary["missing_evidence_encoded_as_zero"]
    or any(
        row["genomic_test_fragments_parsed_and_excluded"] <= 0
        for row in summary["worker_summaries"]
    )
    or (root / "inputs" / "lineage_matrices").exists()
):
    raise SystemExit("scBasset production input summary differs")
freeze_tree(
    root,
    {
        "artifact_class": "scbasset_training_inputs",
        "model_id": "scbasset",
        "dataset_id": "gse296875",
        "split_id": split_id,
        "donor_test_fold": outer_fold,
        "donor_valid_fold": valid_fold,
        "biological_outer_unit": "donor",
        "pseudobulk_artifacts_sha256": sys.argv[2],
        "membership_artifacts_sha256": "582387b78ec25215d45820c24aace12f130b826e91755bf19ac4581730611e64",
        "split_artifacts_sha256": "00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f",
        "reference_artifacts_sha256": "7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d",
        "training_initialization": "from_scratch",
        "public_tutorial_weights_used": False,
        "held_donor_atac_used": False,
        "genomic_test_atac_parsed_only_to_exclude": True,
        "genomic_test_atac_used": False,
        "genomic_test_atac_exported": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scBasset %s inputs: %s\n' "${SPLIT_ID}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
