#!/usr/bin/env bash
#SBATCH --job-name=build_test_v1_1_releases
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/release/build_test_v1_1_releases_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/release/build_test_v1_1_releases_%j.err
set -euo pipefail

BENCH="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
REL="${BENCH}/release"
PY_TF="${BENCH}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
PY_RN="${HOME}/micromamba/envs/rnaseq/bin/python"
export PYTHONNOUSERSITE=1

echo "=== 1. Build masld-panel-v1.1 weights & expected.json ==="
"${PY_TF}" -u "${REL}/masld-panel-v1.1/build_panel_v1_1.py"

echo "=== 2. Build masld-liver-chromatin-state-v1.2.1 weights ==="
"${PY_RN}" -u "${REL}/masld-liver-chromatin-state-v1.2.1/build_weights_v1_2_1.py"

echo "=== 3. Generate MANIFEST.sha256.json for all three .1 bundles ==="
"${PY_RN}" - "${REL}" <<'PYEOF'
import hashlib, json, os, sys
from pathlib import Path

rel_root = Path(sys.argv[1])
bundles = [
    ("masld-severity-v1.2.1", "masld-severity-v1.2 (score.py and weights bitwise identical; kleiner_sensitivity.py replaces kleiner.py)"),
    ("masld-panel-v1.1", "masld-panel-v1 (label-free abundance/biotype refit + sample-independent NaN handling)"),
    ("masld-liver-chromatin-state-v1.2.1", "masld-liver-chromatin-state-v1.2 (v1.1 states bitwise; form-specific prof_reliable_<form> masks)"),
]

for name, note in bundles:
    bdir = rel_root / name
    files = {}
    for p in sorted(bdir.rglob("*")):
        if p.is_dir():
            continue
        rel = p.relative_to(bdir).as_posix()
        if rel == "MANIFEST.sha256.json" or "__pycache__" in rel or rel.endswith(".pyc") or os.path.basename(rel).startswith("."):
            continue
        raw = p.read_bytes()
        files[rel] = {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
    manifest = {"release": name, "_supersedes": note, "files": files}
    (bdir / "MANIFEST.sha256.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {name}/MANIFEST.sha256.json ({len(files)} files)")
PYEOF

chmod +x "${REL}/masld-severity-v1.2.1/test_release.sh" \
         "${REL}/masld-panel-v1.1/test_release.sh" \
         "${REL}/masld-liver-chromatin-state-v1.2.1/test_release.sh"

echo "=== 4. Run test_release.sh for masld-severity-v1.2.1 ==="
bash "${REL}/masld-severity-v1.2.1/test_release.sh"

echo "=== 5. Run test_release.sh for masld-panel-v1.1 ==="
bash "${REL}/masld-panel-v1.1/test_release.sh"

echo "=== 6. Run test_release.sh for masld-liver-chromatin-state-v1.2.1 ==="
bash "${REL}/masld-liver-chromatin-state-v1.2.1/test_release.sh"

echo "=== ALL THREE .1 RELEASES BUILT AND VERIFIED ==="
