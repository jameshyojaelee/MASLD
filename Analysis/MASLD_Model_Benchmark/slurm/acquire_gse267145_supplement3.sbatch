#!/usr/bin/env bash
#SBATCH --job-name=gse267145-s3
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
JOIN="${ROOT}/executions/gse267145-authoritative-join-21064930"
SOURCE_URL=https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11071482/supplementaryFiles
TARGET="${ROOT}/executions/gse267145-supplement3-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE267145 supplement artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/raw" "${STAGE}/validation"
printf '%s  %s\n' \
  e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580 \
  "${JOIN}/ARTIFACTS.json" | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${JOIN}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]).resolve(strict=True))
PY
curl --proto '=https' --tlsv1.2 --fail --location --retry 5 \
  --output "${STAGE}/raw/supplementaryFiles.zip" "${SOURCE_URL}"
python - "${STAGE}" "${SOURCE_URL}" <<'PY'
from pathlib import Path
import json
import sys
from zipfile import ZipFile, BadZipFile

root = Path(sys.argv[1]).resolve(strict=True)
path = root / "raw/media-3.xlsx"
with ZipFile(root / "raw/supplementaryFiles.zip") as archive:
    matches = [name for name in archive.namelist() if Path(name).name == "media-3.xlsx"]
    if len(matches) != 1:
        raise SystemExit("Europe PMC supplement package does not contain one media-3.xlsx")
    path.write_bytes(archive.read(matches[0]))
try:
    with ZipFile(path) as archive:
        names = set(archive.namelist())
        required = {"[Content_Types].xml", "xl/workbook.xml", "xl/sharedStrings.xml"}
        if not required.issubset(names) or archive.testzip() is not None:
            raise SystemExit("Supplementary workbook structure differs")
        workbook_members = sorted(name for name in names if name.startswith("xl/worksheets/sheet"))
except BadZipFile as error:
    raise SystemExit("Supplementary workbook is not a valid XLSX") from error
if not workbook_members:
    raise SystemExit("Supplementary workbook has no worksheets")
receipt = {
    "schema_version": "masld-bench-gse267145-supplement3-acquisition-v1",
    "status": "pass_primary_source_acquisition",
    "source_url": sys.argv[2],
    "source_pmc": "PMC11071482",
    "source_doi": "10.1101/2024.04.25.591188",
    "supplement_label": "Supplementary Table 3: Differential peaks for all CUT&RUN datasets",
    "worksheet_members": workbook_members,
    "coordinate_semantics_resolved": False,
    "model_fitted_or_selected": False,
    "development_source_owned_outcomes_present": True,
    "sealed_or_test_outcomes_read": False,
}
(root / "validation/receipt.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum slurm/acquire_gse267145_supplement3.sbatch > "${STAGE}/source.sha256"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree
freeze_tree(
    Path(sys.argv[1]).resolve(strict=True),
    {
        "artifact_class": "gse267145_primary_source_supplement3",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "coordinate_semantics_resolved": False,
        "model_fitted_or_selected": False,
        "development_source_owned_outcomes_present": True,
        "sealed_or_test_outcomes_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE267145 Supplementary Table 3: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
