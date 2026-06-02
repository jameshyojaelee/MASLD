#!/usr/bin/env bash
#SBATCH --job-name=spatial_multiomics_env
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/spatial_multiomics_env_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/spatial_multiomics_env_%j.err
##############################################################################
# Create spatial_multiomics environment with:
#   - Core spatial packages (scanpy, squidpy, anndata, cell2location)
#   - torch 2.6.0+cu124 (pinned) + torch-geometric
#   - MultiSP (+ episcanpy dep), ONTraC
#
# Key: Install ONTraC with --no-deps to prevent it pulling torch 2.7.x
#      and downgrading numpy. Then manually install its real deps.
##############################################################################
set -euo pipefail

eval "$(micromamba shell hook -s bash)"

# ── Step 1: Create fresh env ──────────────────────────────────────────────
echo "=== Creating spatial_multiomics env ==="
micromamba env remove -n spatial_multiomics -y 2>/dev/null || true
micromamba create -n spatial_multiomics python=3.12 -y -c conda-forge 2>&1 | tail -5

micromamba activate spatial_multiomics

echo ""
echo "Python: $(which python)"
echo "Env:    ${CONDA_PREFIX}"

# ── Step 2: Install PyTorch 2.6.0+cu124 (PINNED) ─────────────────────────
echo ""
echo "=== Installing torch 2.6.0+cu124 (pinned) ==="
pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124 2>&1 | tail -5

# ── Step 3: Core scientific stack ─────────────────────────────────────────
echo ""
echo "=== Installing core packages ==="
pip install "numpy>=2.0" "pandas>=2.2" scipy scikit-learn matplotlib seaborn 2>&1 | tail -5
pip install scanpy==1.12 squidpy==1.8.1 anndata 2>&1 | tail -5
pip install scvi-tools cell2location 2>&1 | tail -5
pip install pyyaml statsmodels 2>&1 | tail -5

# ── Step 4: torch-geometric ──────────────────────────────────────────────
echo ""
echo "=== Installing torch-geometric ==="
pip install torch-geometric 2>&1 | tail -5
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv \
    -f https://data.pyg.org/whl/torch-2.6.0+cu124.html 2>&1 | tail -5

# ── Step 5: MultiSP + episcanpy ──────────────────────────────────────────
echo ""
echo "=== Installing MultiSP + episcanpy ==="
pip install episcanpy 2>&1 | tail -5
pip install multisp 2>&1 | tail -5

# ── Step 6: ONTraC (--no-deps to avoid torch/numpy conflicts) ────────────
echo ""
echo "=== Installing ONTraC (--no-deps) ==="
pip install ONTraC --no-deps 2>&1 | tail -5

# Install ONTraC's actual runtime deps that aren't already installed
pip install harmonypy igraph leidenalg session-info 2>&1 | tail -5

# ── Step 7: snapatac2 ───────────────────────────────────────────────────
echo ""
echo "=== Installing snapatac2 ==="
pip install snapatac2 2>&1 | tail -5

# ── Step 8: Verify ───────────────────────────────────────────────────────
echo ""
echo "=== Verification ==="
python << 'PYEOF'
import sys

required = [
    ("torch", "torch"),
    ("scanpy", "scanpy"),
    ("squidpy", "squidpy"),
    ("anndata", "anndata"),
    ("cell2location", "cell2location"),
    ("numpy", "numpy"),
    ("ONTraC", "ONTraC"),
]

optional = [
    ("torch_geometric", "torch_geometric"),
    ("torch_scatter", "torch_scatter"),
    ("torch_sparse", "torch_sparse"),
    ("multisp", "multisp"),
    ("episcanpy", "episcanpy"),
    ("snapatac2", "snapatac2"),
]

failed_required = []
for name, imp in required + optional:
    try:
        mod = __import__(imp)
        ver = getattr(mod, "__version__", "OK")
        is_optional = (name, imp) in optional
        print(f"  {'[opt]' if is_optional else '[req]'} {name}: {ver}")
    except Exception as e:
        is_optional = (name, imp) in optional
        label = "OPTIONAL" if is_optional else "REQUIRED"
        print(f"  [{label}] {name}: FAILED ({str(e)[:80]})")
        if not is_optional:
            failed_required.append(name)

# Critical checks
import torch
print(f"\n  torch.version.cuda: {torch.version.cuda}")
print(f"  torch version: {torch.__version__}")

import numpy
print(f"  numpy version: {numpy.__version__}")
assert numpy.__version__.startswith("2"), f"numpy must be >=2.0, got {numpy.__version__}"
assert "2.6.0" in torch.__version__, f"torch must be 2.6.0, got {torch.__version__}"

if failed_required:
    print(f"\nFAILED REQUIRED: {failed_required}")
    sys.exit(1)
else:
    print("\nALL REQUIRED PACKAGES OK!")
PYEOF

echo ""
echo "=== Complete at $(date) ==="
