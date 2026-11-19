#!/usr/bin/env bash
#SBATCH --job-name=final_compile_verify_261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/reports/final_compile_verify_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/reports/final_compile_verify_20261118_%j.err

set -euo pipefail

ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
TECTONIC="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/tectonic"
PY_SPATIAL="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
PYTEST_SPATIAL="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/pytest"
PYTEST_RNASEQ="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/pytest"
PY_BENCH="${ROOT}/Analysis/MASLD_Model_Benchmark/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"

cd "${ROOT}"
echo "=== [0/6] Updating Fig 4A schematic (650K unique cells) ==="
"${PY_SPATIAL}" figures/candidates/fig4-20261118/make_fig4a_20261118.py

echo "=== [1/6] Compiling main.tex and supplement.tex with Tectonic ==="
cd "${ROOT}/docs/manuscript/latex"
"${TECTONIC}" -X compile --outdir build --keep-logs main.tex > build/main-console.txt 2>&1
"${TECTONIC}" -X compile --outdir build --keep-logs supplement.tex > build/supplement-console.txt 2>&1
ls -lh build/main.pdf build/supplement.pdf
pdfinfo build/main.pdf | grep -E "Pages|File size" || true
pdfinfo build/supplement.pdf | grep -E "Pages|File size" || true
echo "main.tex error lines:" $(grep -c "^error:" build/main-console.txt || true)
echo "supplement.tex error lines:" $(grep -c "^error:" build/supplement-console.txt || true)

echo "=== [2/6] Rendering compiled main figure plate pages (150 DPI) ==="
PREV_DIR="${ROOT}/reports/compiled_plate_previews_20261118"
mkdir -p "${PREV_DIR}"
"${PY_SPATIAL}" - <<'PYEOF'
import subprocess, re
from pathlib import Path

root = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
main_pdf = root / "docs/manuscript/latex/build/main.pdf"
out_dir = root / "reports/compiled_plate_previews_20261118"
info = subprocess.check_output(["pdfinfo", str(main_pdf)], text=True)
n_pages = int(re.search(r"Pages:\s+(\d+)", info).group(1))
print(f"main.pdf has {n_pages} pages")
fig_pages = {}
for p in range(1, n_pages + 1):
    txt = subprocess.check_output(["pdftotext", "-f", str(p), "-l", str(p), str(main_pdf), "-"], text=True)
    for k in range(1, 7):
        if re.search(rf"(?m)^Figure {k}\s*$", txt) and k not in fig_pages:
            fig_pages[k] = p
print("Detected main figure plate pages:", fig_pages)
for k, p in sorted(fig_pages.items()):
    prefix = out_dir / f"main_plate_fig{k}_p{p}"
    subprocess.check_call(["pdftoppm", "-png", "-r", "150", "-f", str(p), "-l", str(p), "-singlefile", str(main_pdf), str(prefix)])
    print(f"Rendered {prefix}.png")
PYEOF

cd "${ROOT}"
echo "=== [3/6] Running 5 unit test suites ==="
echo "--- Suite 1: scripts/portal ---"
"${PYTEST_SPATIAL}" -q scripts/portal/test_evidence_passport_bundle.py scripts/portal/test_spatial_resource_semantics.py scripts/portal/test_gene_catalog_observability.py scripts/portal/test_gene_catalog_v2.py
echo "--- Suite 2: scripts/figures ---"
"${PYTEST_SPATIAL}" -q scripts/figures/tests/test_fig2_cluster_inference.py scripts/figures/test_fig2_locus_coordinates.py
echo "--- Suite 3: scripts/analysis/masld_two_models ---"
"${PYTEST_SPATIAL}" -q scripts/analysis/masld_two_models/test_scientific_invariants.py
echo "--- Suite 4: scripts/analysis/alt_contig_recount (rnaseq env w/ pysam) ---"
"${PYTEST_RNASEQ}" -q scripts/analysis/alt_contig_recount/test_annotation_correspondence.py scripts/analysis/alt_contig_recount/test_count_execution_evidence.py scripts/analysis/alt_contig_recount/test_native_library_subsample.py scripts/analysis/alt_contig_recount/test_verified_fragments.py
echo "--- Suite 5: Analysis/MASLD_Model_Benchmark/tests ---"
(cd "${ROOT}/Analysis/MASLD_Model_Benchmark" && PYTHONPATH=src "${PY_BENCH}" -m pytest -q tests)

echo "=== [4/6] Running 3 .1 model release test_release.sh scripts ==="
bash Analysis/MASLD_Model_Benchmark/release/masld-severity-v1.2.1/test_release.sh
bash Analysis/MASLD_Model_Benchmark/release/masld-panel-v1.1/test_release.sh
bash Analysis/MASLD_Model_Benchmark/release/masld-liver-chromatin-state-v1.2.1/test_release.sh

echo "=== [5/6] Running validate_resource_scope.py --strict-latex ==="
python3 scripts/manuscript/validate_resource_scope.py --strict-latex

echo "=== [6/6] Running verify_remediation_20261118.py (N01-N23) ==="
python3 reports/verify_remediation_20261118.py

echo "=== ALL 6 STAGES COMPLETED WITH EXIT 0 ==="
