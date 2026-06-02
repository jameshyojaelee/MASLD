#!/usr/bin/env bash
#SBATCH --job-name=install_commot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/install_commot_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/install_commot_%j.err
##############################################################################
# Install COMMOT (COMMunication analysis by Optimal Transport) into the
# existing spatial environment.
##############################################################################
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

echo "Python: $(which python)"
echo "Env:    ${CONDA_PREFIX}"

# ── Step 1: Install commot ─────────────────────────────────────────────────
echo ""
echo "=== Installing COMMOT ==="
pip install commot 2>&1 | tail -10

# ── Step 1b: Patch numpy 2.x incompatibility (np.Inf -> np.inf) ──────────
echo ""
echo "=== Patching COMMOT for numpy 2.x compatibility ==="
COMMOT_DIR="${CONDA_PREFIX}/lib/python3.12/site-packages/commot"
if [ -d "$COMMOT_DIR" ]; then
    grep -rl 'np\.Inf' "$COMMOT_DIR" 2>/dev/null | while read f; do
        echo "  Patching: $f"
        sed -i 's/np\.Inf/np.inf/g' "$f"
    done
    echo "  Patch complete"
else
    echo "  WARNING: COMMOT directory not found at $COMMOT_DIR"
fi

# ── Step 2: Verify ─────────────────────────────────────────────────────────
echo ""
echo "=== Verification ==="
python << 'PYEOF'
packages = [
    ("commot", "commot"),
    ("squidpy", "squidpy"),
    ("scanpy", "scanpy"),
    ("cell2location", "cell2location"),
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
