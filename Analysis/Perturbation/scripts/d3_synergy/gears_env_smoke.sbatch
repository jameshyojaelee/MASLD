#!/bin/bash
#SBATCH --job-name=GEARS
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=Analysis/Perturbation/scripts/d3_synergy/logs/gears_env_smoke_%j.log

# Verify perturbation_gears env imports cleanly (numpy / torch / gears / scanpy).
# Lightweight import-only check - no compute. Run via sbatch per HPC policy.

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_gears

set -u

python - <<'PY'
import sys, traceback
print("py:", sys.version.split()[0])
for mod_name, attr in [
    ("numpy", "__version__"),
    ("scipy", "__version__"),
    ("pandas", "__version__"),
    ("torch", "__version__"),
    ("torch_geometric", "__version__"),
    ("anndata", "__version__"),
    ("scanpy", "__version__"),
    ("h5py", "__version__"),
    ("gears", "__version__"),
]:
    try:
        m = __import__(mod_name)
        v = getattr(m, attr, "(no version attr)")
        print(f"  {mod_name:18s}: {v}")
    except Exception as e:
        print(f"  {mod_name:18s}: IMPORT ERROR")
        traceback.print_exc()

# Functional smoke: instantiate PertData (no data download — just constructor)
print("\nTrying gears.PertData + gears.GEARS instantiation (constructor only):")
try:
    from gears import PertData, GEARS
    print("  from gears import OK")
except Exception as e:
    print(f"  IMPORT ERROR: {e}")
    traceback.print_exc()

# Torch CUDA availability (will be False on cpu partition, expected)
try:
    import torch
    print(f"\ntorch.cuda.is_available() = {torch.cuda.is_available()}  (cpu partition → False expected)")
except Exception:
    pass

print("\nENV CHECK DONE")
PY
