#!/bin/bash
#SBATCH --job-name=mcp_loo_build_LiverAtlas
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --mem=200G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/470c_build_LiverAtlas_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/470c_build_LiverAtlas_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"
source "$HOME/.bashrc"
micromamba activate rapids_singlecell

python - <<'PY'
from pathlib import Path
import anndata as ad

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
inp = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/atlas_cnmf_global.h5ad"
out_dir = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs"

print(f"[build_loo_LiverAtlas] reading {inp}")
atlas = ad.read_h5ad(inp)
print(f"[build_loo_LiverAtlas] atlas: {atlas.shape}")
sub = atlas[atlas.obs["dataset"] != "Liver_Atlas"].copy()
out = out_dir / "atlas_cnmf_global_loo_Liver_Atlas.h5ad"
print(f"[build_loo_LiverAtlas] writing {out.name} ({sub.shape[0]:,} cells)")
sub.write_h5ad(out, compression="gzip")
print("[build_loo_LiverAtlas] DONE.")
PY
