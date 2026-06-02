#!/usr/bin/env bash
#SBATCH --job-name=setup_gsmap_env
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/setup_gsmap_env_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/setup_gsmap_env_%j.err
##############################################################################
# Create gsmap conda env with python=3.11, install gsMap, and download the
# gsMap resource bundle (~15GB) including:
#   - LD reference (1000G EUR Phase 3)
#   - HapMap3 SNPs
#   - LD weights
#   - Gene annotation (gencode.v39lift37, hg19)
#   - Pre-computed SNP-gene weights
##############################################################################
set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RESOURCE_DIR="${PROJECT_ROOT}/data/gsmap_resource"

eval "$(micromamba shell hook -s bash)"

# ── Step 1: Create gsmap env ────────────────────────────────────────────────
echo "=== Creating gsmap environment (python=3.11, numpy<2.0) ==="
micromamba env remove -n gsmap -y 2>/dev/null || true
micromamba create -n gsmap python=3.11 "numpy<2.0" -y -c conda-forge 2>&1 | tail -5

micromamba activate gsmap

echo ""
echo "Python: $(which python)"
echo "Env:    ${CONDA_PREFIX}"

# ── Step 2: Install gsMap via pip ───────────────────────────────────────────
echo ""
echo "=== Installing gsMap ==="
pip install gsMap pyyaml 2>&1 | tail -10

# ── Step 3: Download resource bundle ───────────────────────────────────────
echo ""
echo "=== Downloading gsMap resource bundle to ${RESOURCE_DIR} ==="
mkdir -p "${RESOURCE_DIR}"

# Download from official gsMap resource server
cd "${RESOURCE_DIR}"

echo "  Downloading gsMap_resource.tar.gz (~15GB)..."
wget -q --show-progress https://yanglab.westlake.edu.cn/data/gsMap/gsMap_resource.tar.gz

echo ""
echo "  Extracting..."
tar -xzf gsMap_resource.tar.gz
rm -f gsMap_resource.tar.gz

# Move contents up if nested in a subdirectory
if [ -d "${RESOURCE_DIR}/gsMap_resource" ]; then
    mv "${RESOURCE_DIR}/gsMap_resource"/* "${RESOURCE_DIR}/"
    rmdir "${RESOURCE_DIR}/gsMap_resource"
fi

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo ""
echo "=== Resource directory contents ==="
du -sh "${RESOURCE_DIR}"
ls -lh "${RESOURCE_DIR}/"

# ── Step 4: Verify ─────────────────────────────────────────────────────────
echo ""
echo "=== Verification ==="
python << 'PYEOF'
packages = [
    ("gsMap", "gsmap"),
    ("numpy", "numpy"),
    ("pandas", "pandas"),
    ("scipy", "scipy"),
]

results = {}
for name, imp in packages:
    try:
        mod = __import__(imp)
        ver = getattr(mod, "__version__", "OK")
        results[name] = ver
        print(f"  {name}: {ver}")
    except Exception as e:
        results[name] = "FAILED"
        print(f"  {name}: FAILED ({str(e)[:120]})")

all_ok = all("FAILED" not in str(v) for v in results.values())
print(f"\n{'ALL PACKAGES OK!' if all_ok else 'Some packages had issues — see above'}")
PYEOF

echo ""
echo "=== Complete at $(date) ==="
