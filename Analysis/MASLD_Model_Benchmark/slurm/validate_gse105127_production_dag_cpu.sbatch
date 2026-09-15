#!/usr/bin/env bash
#SBATCH --job-name=data-check-200
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-check-200-%j.out
#SBATCH --error=executions/logs/data-check-200-%j.err

set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
FINAL=${ROOT}/executions/gse105127-production-validation-${SLURM_JOB_ID}
STAGE=${FINAL}.staging
FAILED=${FINAL}.failed

cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONPATH="${ROOT}/src:${ROOT}"
if [[ -e ${FINAL} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to replace validation attempt\n' >&2
  exit 1
fi
mkdir -p "${STAGE}/runtime-fixture"
trap 'status=$?; if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi; exit ${status}' EXIT
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
"${PYTHON}" -m pip freeze > "${STAGE}/runtime-fixture/pip-freeze.txt"
scontrol show job "${SLURM_JOB_ID}" -o > "${STAGE}/runtime-fixture/slurm-job.txt"

for script in slurm/gse105127_*.sbatch slurm/validate_gse105127_production_dag_cpu.sbatch slurm/gse105127_production_common.sh slurm/submit_gse105127_production_dag.sh; do
  bash -n "${script}"
done
"${PYTHON}" - <<'PY'
import ast
from pathlib import Path

for path in sorted(Path("scripts").glob("*gse105127*.py")):
    ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
PY
"${PYTHON}" -m unittest \
  tests.unit.test_build_gse105127_assay_native_activation \
  tests.unit.test_gse105127_production_dag

module purge
module load SAMtools/1.18-GCC-12.3.0 BEDTools/2.31.0-GCC-12.3.0 crossmap/0.7.0
samtools --version | head -1 > "${STAGE}/runtime-fixture/samtools.version.txt"
bedtools --version > "${STAGE}/runtime-fixture/bedtools.version.txt"
CrossMap bed -h > "${STAGE}/runtime-fixture/crossmap.help.txt" 2>&1
rg --fixed-strings -- '--chromid' "${STAGE}/runtime-fixture/crossmap.help.txt"
rg --fixed-strings -- '--unmap-file' "${STAGE}/runtime-fixture/crossmap.help.txt"

printf '>chr1\nACGTACGT\n' > "${STAGE}/runtime-fixture/tiny.fa"
samtools faidx "${STAGE}/runtime-fixture/tiny.fa"
printf 'chr1\t1\t3\tcpg_000000001\n' > "${STAGE}/runtime-fixture/tiny.bed"
bedtools getfasta -fi "${STAGE}/runtime-fixture/tiny.fa" \
  -bed "${STAGE}/runtime-fixture/tiny.bed" -nameOnly -tab \
  > "${STAGE}/runtime-fixture/tiny.sequence.tsv"
printf 'cpg_000000001\tCG\n' | cmp - "${STAGE}/runtime-fixture/tiny.sequence.tsv"

curl --fail --location --retry 3 \
  --output "${STAGE}/runtime-fixture/hg19ToHg38.over.chain.gz" \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
curl --fail --location --retry 3 \
  --output "${STAGE}/runtime-fixture/hg38ToHg19.over.chain.gz" \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
printf '%s  %s\n' 35887f73fe5e2231656504d1f6430900 "${STAGE}/runtime-fixture/hg19ToHg38.over.chain.gz" | md5sum --check --strict
printf '%s  %s\n' ff3031d93792f4cbb86af44055efd903 "${STAGE}/runtime-fixture/hg38ToHg19.over.chain.gz" | md5sum --check --strict
printf 'chr1\t1000000\t1000002\tcpg_000000001\n' > "${STAGE}/runtime-fixture/source.bed"
CrossMap bed --chromid l --unmap-file "${STAGE}/runtime-fixture/forward.unmap.bed" \
  "${STAGE}/runtime-fixture/hg19ToHg38.over.chain.gz" \
  "${STAGE}/runtime-fixture/source.bed" "${STAGE}/runtime-fixture/forward.bed"
CrossMap bed --chromid l --unmap-file "${STAGE}/runtime-fixture/reverse.unmap.bed" \
  "${STAGE}/runtime-fixture/hg38ToHg19.over.chain.gz" \
  "${STAGE}/runtime-fixture/forward.bed" "${STAGE}/runtime-fixture/reverse.bed"
cmp "${STAGE}/runtime-fixture/source.bed" "${STAGE}/runtime-fixture/reverse.bed"

module purge
module load RSEM/1.3.3-foss-2022a STAR/2.7.10b-GCC-11.3.0
rsem-calculate-expression --version > "${STAGE}/runtime-fixture/rsem.version.txt" 2>&1
STAR --version > "${STAGE}/runtime-fixture/star.version.txt"
rg --fixed-strings 'RSEM v1.3.1' "${STAGE}/runtime-fixture/rsem.version.txt"
printf '2.7.10b\n' | cmp - "${STAGE}/runtime-fixture/star.version.txt"
rsem-calculate-expression --help > "${STAGE}/runtime-fixture/rsem-calculate.help.txt" 2>&1 || true
rsem-prepare-reference --help > "${STAGE}/runtime-fixture/rsem-prepare.help.txt" 2>&1 || true
[[ -s "${STAGE}/runtime-fixture/rsem-calculate.help.txt" ]]
[[ -s "${STAGE}/runtime-fixture/rsem-prepare.help.txt" ]]
for option in --star-gzipped-read-file --no-bam-output --strandedness --seed; do
  rg --fixed-strings -- "${option}" "${STAGE}/runtime-fixture/rsem-calculate.help.txt"
done
rg --fixed-strings -- '--star-sjdboverhang' "${STAGE}/runtime-fixture/rsem-prepare.help.txt"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1])
(root / "receipt.json").write_text(
    json.dumps(
        {
            "schema_version": "masld-bench-gse105127-production-validation-v1",
            "status": "passed",
            "labels_accessed": False,
            "fit_or_score_performed": False,
            "runtime_fixture": "CrossMap_roundtrip_bedtools_name_RSEM_STAR",
        },
        indent=2,
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
freeze_tree(root, {"artifact_class": "gse105127_production_validation", "status": "passed"})
verify_frozen_tree(root)
PY
mv "${STAGE}" "${FINAL}"
"${PYTHON}" - "${FINAL}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
PY
trap - EXIT
printf '%s\n' "${FINAL}"
