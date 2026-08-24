#!/usr/bin/env bash
#SBATCH --job-name=masld-chrombpnet-inputs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
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
: "${LINEAGE:=hepatocyte}"

case "${OUTER_FOLD}" in 0|1|2|3|4) ;; *) printf 'Invalid OUTER_FOLD\n' >&2; exit 1 ;; esac
case "${LINEAGE}" in cholangiocyte|fibroblast|hepatocyte|macrophage|t_cell) ;;
  *) printf 'Invalid LINEAGE\n' >&2; exit 1 ;;
esac
SPLIT_ID="donor${OUTER_FOLD}_genomic${OUTER_FOLD}"
GENOMIC_VALID_FOLD=$(( (OUTER_FOLD + 1) % 5 ))

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PSEUDOBULK="${ROOT}/executions/gse296875-training-fold-pseudobulk-21063993"
PSEUDOBULK_ARTIFACTS_SHA256=3190aa76cd5872f50d7b0c71bcb92520f1f2428110ce2fb4c26856994bc17e5b
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
REFERENCE="${ROOT}/executions/sequence-reference-contract-21063670"
CHROM_SIZES="${REFERENCE}/GRCh38.p14.sequence_model.chrom.sizes"
FASTA_ROOT="${ROOT}/executions/sequence-fasta-21063808"
FASTA="${FASTA_ROOT}/GRCh38.p14.sequence_model.fa"
BLACKLIST_ROOT="${ROOT}/executions/encode-grch38-blacklist-21063857"
BLACKLIST_GZ="${BLACKLIST_ROOT}/ENCFF356LFX.bed.gz"
MACS_PROBE="${ROOT}/executions/macs2-runtime-probe-21063873"
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
PREP="${ROOT}/executions/environments/chrombpnet-native-preprocessing-overlay-21063949"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/chrombpnet-v1.0.1-21063022"
SOURCE="${ACQUISITION}/source/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c"
BEDTOOLS_RUNTIME="${ROOT}/executions/environments/bedtools-bioconda-2.31.1-21063932"
BEDTOOLS_BIN="${BEDTOOLS_RUNTIME}/runtime/bin"
TARGET="${ROOT}/executions/chrombpnet-${LINEAGE}-${SPLIT_ID}-inputs-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
PYTHONPATH_CORE="${ROOT}:${CORE}/site-packages"
PYTHONPATH_PREP="${PREP}/site-packages:${CORE}/site-packages:${SOURCE}"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed ChromBPNet input stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite ChromBPNet input campaign: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/macs2" \
  "${STAGE}/regions" \
  "${STAGE}/nonpeaks" \
  "${STAGE}/validation" \
  "${STAGE}/mplconfig"

printf '%s  %s\n' "${PSEUDOBULK_ARTIFACTS_SHA256}" \
  "${PSEUDOBULK}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
  "${SPLIT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d \
  "${REFERENCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 1901777f65ef618569eabb62967c5ace55a63b4bef8795af7a8838fc7042858b \
  "${CHROM_SIZES}" | sha256sum --check --strict
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
printf '%s  %s\n' f1fc72d551449ec0f917d5453183870c8486720f72ff27d963adad5685720f73 \
  "${PREP}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 04bd59c825874f2d2d7dbf38b741633ccacad4735e95997a7477288e9c9656c4 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 81ed99952bb5caf92cb129e7ea87713468f1eebb143ba0efe6542bfd1a00bc68 \
  "${BEDTOOLS_RUNTIME}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${PSEUDOBULK}" "${SPLIT}" "${REFERENCE}" "${FASTA_ROOT}" \
  "${BLACKLIST_ROOT}" "${MACS_PROBE}" "${BEDTOOLS_RUNTIME}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("input frozen-tree verification passed")
PY
python - "${BLACKLIST_GZ}" "${STAGE}/regions/blacklist.bed" <<'PY'
from pathlib import Path
import gzip
import shutil
import sys

source = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2])
with gzip.open(source, "rb") as reader, target.open("xb") as writer:
    shutil.copyfileobj(reader, writer)
PY

module purge
module load singularity/4.2.2
cd "${ROOT}"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_CORE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python scripts/emit_training_tagalign.py \
    --pseudobulk "${PSEUDOBULK}" \
    --pseudobulk-artifacts-sha256 "${PSEUDOBULK_ARTIFACTS_SHA256}" \
    --split-contract "${SPLIT}" \
    --split-artifacts-sha256 \
      00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
    --donor-test-fold "${OUTER_FOLD}" \
    --genomic-test-fold "${OUTER_FOLD}" \
    --genomic-valid-fold "${GENOMIC_VALID_FOLD}" \
    --lineage "${LINEAGE}" \
    --output-directory "${STAGE}/tagalign" \
    --summary "${STAGE}/tagalign/summary.json" \
    > "${STAGE}/validation/tagalign.stdout.txt"
gzip -t "${STAGE}"/tagalign/genomic_fold_*.tagAlign.gz

module purge
module load MACS2/2.2.9.1-foss-2023b
macs2 --version > "${STAGE}/validation/macs2.version.txt"
python - "${SPLIT}/genomic_folds.tsv" "${STAGE}/tagalign/summary.json" \
  "${STAGE}/validation/fold_gsizes.tsv" <<'PY'
from pathlib import Path
import csv
import json
import sys

sizes = {}
with Path(sys.argv[1]).open(encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        sizes[int(row["genomic_fold"])] = int(row["total_bp"])
summary = json.loads(Path(sys.argv[2]).read_text())
with Path(sys.argv[3]).open("x", encoding="utf-8") as handle:
    for fold in summary["allowed_genomic_folds"]:
        handle.write(f"{fold}\t{sizes[fold]}\n")
PY
pids=()
while IFS=$'\t' read -r genomic_fold genome_size; do
  macs2 callpeak \
    --treatment "${STAGE}/tagalign/genomic_fold_${genomic_fold}.tagAlign.gz" \
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
  printf 'At least one training-fold MACS2 call failed\n' >&2
  exit 1
fi
python - "${STAGE}/validation/fold_gsizes.tsv" "${STAGE}/macs2" \
  "${STAGE}/regions/macs2.training.narrowPeak" <<'PY'
from pathlib import Path
import sys

with Path(sys.argv[1]).open(encoding="utf-8") as handle:
    folds = [int(line.split("\t", 1)[0]) for line in handle if line.strip()]
if len(folds) != 3 or len(set(folds)) != 3 or any(fold not in range(5) for fold in folds):
    raise SystemExit("training genomic fold roster differs")
source = Path(sys.argv[2]).resolve(strict=True)
target = Path(sys.argv[3])
rows = 0
with target.open("x", encoding="utf-8") as writer:
    for fold in folds:
        path = source / f"genomic_fold_{fold}_peaks.narrowPeak"
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"MACS2 fold {fold} narrowPeak is empty")
        with path.open(encoding="utf-8") as reader:
            for line in reader:
                writer.write(line)
                rows += 1
if rows == 0:
    raise SystemExit("combined MACS2 narrowPeak is empty")
print(rows)
PY
python -m pip freeze > "${STAGE}/validation/macs2.pip-freeze.txt"

module purge
module load singularity/4.2.2
cd "${ROOT}"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_CORE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python scripts/filter_training_peaks.py \
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
  "${IMAGE}" python scripts/assemble_chrombpnet_fold_regions.py \
    --training-peaks "${STAGE}/regions/training.filtered.narrowPeak" \
    --split-contract "${SPLIT}" \
    --split-id "${SPLIT_ID}" \
    --output-peaks "${STAGE}/regions/chrombpnet.regions.bed" \
    --output-fold "${STAGE}/regions/chrombpnet.fold.json" \
    --output-summary "${STAGE}/regions/chrombpnet.regions.summary.json" \
    > "${STAGE}/validation/region-assembly.stdout.txt"

export SINGULARITYENV_PREPEND_PATH="${BEDTOOLS_BIN}"
"${BEDTOOLS_BIN}/bedtools" --version > "${STAGE}/validation/bedtools.version.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_PREP}" \
  --env "MPLCONFIGDIR=${STAGE}/mplconfig" \
  --env MPLBACKEND=Agg \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m chrombpnet.CHROMBPNET prep nonpeaks \
    --genome "${FASTA}" \
    --output-prefix "${STAGE}/nonpeaks/chrombpnet" \
    --peaks "${STAGE}/regions/chrombpnet.regions.bed" \
    --chrom-sizes "${CHROM_SIZES}" \
    --chr-fold-path "${STAGE}/regions/chrombpnet.fold.json" \
    --inputlen 2114 \
    --stride 1000 \
    --neg-to-pos-ratio-train 2 \
    --blacklist-regions "${STAGE}/regions/blacklist.bed" \
    --seed 20260824 \
    > "${STAGE}/validation/nonpeaks.stdout.txt" \
    2> "${STAGE}/validation/nonpeaks.stderr.txt"

singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_CORE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${PSEUDOBULK}" "${STAGE}" \
    "${OUTER_FOLD}" "${LINEAGE}" "${SPLIT_ID}" <<'PY'
from pathlib import Path
import csv
import json
import sys

pseudobulk = Path(sys.argv[1]).resolve(strict=True)
root = Path(sys.argv[2]).resolve(strict=True)
outer_fold = int(sys.argv[3])
lineage = sys.argv[4]
split_id = sys.argv[5]
with (pseudobulk / "training_pseudobulk_manifest.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    source_rows = [
        row
        for row in csv.DictReader(handle, delimiter="\t")
        if int(row["donor_test_fold"]) == outer_fold and row["lineage_id"] == lineage
    ]
if len(source_rows) != 1:
    raise SystemExit("source pseudobulk row differs")
source = source_rows[0]
tagalign = json.loads((root / "tagalign" / "summary.json").read_text())
if (
    tagalign["included_unique_fragments"]
    + tagalign["excluded_genomic_valid_fragments"]
    + tagalign["excluded_genomic_test_fragments"]
    != int(source["unique_fragments"])
):
    raise SystemExit("tagAlign fragment accounting differs")
if tagalign["included_tn5_events"] != 2 * tagalign["included_unique_fragments"]:
    raise SystemExit("tagAlign Tn5 event accounting differs")
peak_filter = json.loads(
    (root / "regions" / "training.filtered.summary.json").read_text()
)
if (
    peak_filter["held_genomic_test"] != 0
    or peak_filter["held_genomic_valid"] != 0
    or peak_filter["retained"] < 1000
):
    raise SystemExit("training-only peak contract differs")
regions = json.loads(
    (root / "regions" / "chrombpnet.regions.summary.json").read_text()
)
if (
    regions["split_id"] != split_id
    or tagalign["allowed_genomic_folds"] != regions["genomic_train_folds"]
    or peak_filter["allowed_genomic_folds"] != regions["genomic_train_folds"]
    or regions["donor_test_fold"] != outer_fold
    or regions["role_counts"]["train"] != peak_filter["retained"]
    or regions["role_counts"]["valid"] != 16000
    or regions["role_counts"]["test"] != 16000
    or regions["held_atac_used_for_region_selection"]
    or regions["model_training_may_use_test_role"]
):
    raise SystemExit("assembled region contract differs")
fold = json.loads((root / "regions" / "chrombpnet.fold.json").read_text())
role_by_contig = {
    contig: role for role, contigs in fold.items() for contig in contigs
}
nonpeak_counts = {"train": 0, "valid": 0, "test": 0}
with (root / "nonpeaks" / "chrombpnet_negatives.bed").open(
    encoding="utf-8"
) as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 10 or fields[0] not in role_by_contig:
            raise SystemExit("ChromBPNet nonpeak row differs")
        nonpeak_counts[role_by_contig[fields[0]]] += 1
if any(value == 0 for value in nonpeak_counts.values()):
    raise SystemExit("ChromBPNet nonpeak role is empty")
validation = {
    "schema_version": "masld-bench-chrombpnet-training-inputs-v1",
    "status": "pass",
    "dataset_id": "gse296875",
    "lineage_id": lineage,
    "split_id": split_id,
    "biological_outer_unit": "donor",
    "source_donors": int(source["donors"]),
    "source_nuclei": int(source["nuclei"]),
    "source_unique_fragments": int(source["unique_fragments"]),
    "included_training_unique_fragments": tagalign["included_unique_fragments"],
    "macs2_calls": 3,
    "training_peaks": peak_filter["retained"],
    "fixed_validation_regions": 16000,
    "fixed_test_regions": 16000,
    "nonpeak_role_counts": nonpeak_counts,
    "held_donor_validation_used_for_peak_discovery": False,
    "held_donor_test_used_for_peak_discovery": False,
    "held_genomic_validation_used_for_peak_discovery": False,
    "held_genomic_test_used_for_peak_discovery": False,
    "test_outcomes_available_to_training": False,
}
(root / "validation" / "summary.json").write_text(
    json.dumps(validation, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(validation, sort_keys=True))
PY
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_PREP}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/chrombpnet-preprocessing.pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/emit_training_tagalign.py" \
  "${ROOT}/scripts/filter_training_peaks.py" \
  "${ROOT}/scripts/assemble_chrombpnet_fold_regions.py" \
  "${ROOT}/slurm/build_chrombpnet_hepatocyte_donor0_genomic0_inputs.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${PSEUDOBULK}" "${SPLIT}" "${REFERENCE}" "${FASTA_ROOT}" \
  "${BLACKLIST_ROOT}" "${MACS_PROBE}" "${BEDTOOLS_RUNTIME}" \
  "${STAGE}/validation/input-verification.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

payload = {
    Path(value).name: verify_frozen_tree(Path(value)) for value in sys.argv[1:-1]
}
Path(sys.argv[-1]).write_text(
    json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
PY
PYTHONPATH="${ROOT}/src" python - "${STAGE}" \
  "${PSEUDOBULK_ARTIFACTS_SHA256}" "${LINEAGE}" "${SPLIT_ID}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "chrombpnet_training_inputs",
            "dataset_id": "gse296875",
            "lineage_id": sys.argv[3],
            "split_id": sys.argv[4],
            "biological_outer_unit": "donor",
            "pseudobulk_artifacts_sha256": sys.argv[2],
            "split_artifacts_sha256": "00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f",
            "reference_artifacts_sha256": "7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d",
            "fasta_artifacts_sha256": "7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39",
            "blacklist_artifacts_sha256": "fdbc0d2fbfb23781349a6ba880e29cf4e921cb30857075d5a91395f164c5a540",
            "macs2_version": "2.2.9.1",
            "macs2_mode": "full_depth_per_training_genomic_fold",
            "held_atac_used_for_peak_discovery": False,
            "test_outcomes_available_to_training": False,
            "seed": 20260824,
            "status": "passed",
        },
    )
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet training inputs: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
