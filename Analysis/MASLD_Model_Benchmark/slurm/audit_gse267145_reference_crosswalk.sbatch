#!/usr/bin/env bash
#SBATCH --job-name=masld-gse267145-ref-xwalk
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
MATRIX="${ROOT}/executions/gse267145-matrix-audit-21064333"
SEQUENCE="${ROOT}/executions/sequence-fasta-21063808"
GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
TARGET="${ROOT}/executions/gse267145-reference-crosswalk-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite GSE267145 reference audit: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse267145_reference_crosswalk
PYTHONPATH="src" python - "${MATRIX}" "${SEQUENCE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("input frozen-tree verification passed")
PY
gzip -t "${GTF}"
sha256sum --check --strict <<'EOF'
73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9  /gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
EOF
PYTHONPATH="src:." python scripts/audit_gse267145_reference_crosswalk.py \
  --matrix-audit "${MATRIX}" \
  --matrix-artifacts-sha256 4e0d17f392530f4e6f0b54346d5d9333e450b7606d9698fbf7253a4f2a0e0417 \
  --sequence-reference "${SEQUENCE}" \
  --sequence-artifacts-sha256 7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 \
  --gencode49-gtf "${GTF}" \
  --gencode49-gtf-sha256 73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9 \
  --ensembl98-url https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz \
  --output "${STAGE}"
sha256sum \
  scripts/audit_gse267145_reference_crosswalk.py \
  tests/unit/test_audit_gse267145_reference_crosswalk.py \
  slurm/audit_gse267145_reference_crosswalk.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gse267145_reference_crosswalk",
            "cohort_family_ids": ["gse267145_znf469_human_liver"],
            "rna_gene_crosswalk_frozen": True,
            "h3_primary_contig_bounds_pass": True,
            "h3_coordinate_base_resolved": False,
            "sequence_extraction_allowed": False,
            "model_training_activated": False,
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'GSE267145 reference crosswalk: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
