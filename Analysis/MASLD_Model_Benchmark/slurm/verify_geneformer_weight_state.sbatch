#!/usr/bin/env bash
#SBATCH --job-name=model-check-911
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/geneformer-weight-state-verification-${SLURM_JOB_ID}"
STAGING="${TARGET}.staging"
cd "${ROOT}"; module purge; module load h5py/3.9.0-foss-2023a
mkdir -p "${STAGING}/audit" "${STAGING}/validation"

# READ-ONLY. Re-derive every acquired body from the bytes and compare to the
# bundle. Nothing is corrected here; this establishes what is true.
PYTHONPATH=${ROOT}:${ROOT}/src python - <<'PY' | tee "${STAGING}/audit/geneformer_weight_verification.txt"
import hashlib, json
from pathlib import Path
A=Path("executions/geneformer-acquisition-21014827")
d=json.loads(Path("config/artifacts/models/geneformer/checkpoints.json").read_text())

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as fh:
        for b in iter(lambda: fh.read(64*1024*1024), b""): h.update(b)
    return h.hexdigest()

print("=== 3 weight bodies vs bundle artifacts ===")
tot=0; bad=[]
for a in d["artifacts"]:
    p=A/"source"/a["path"]
    assert p.is_file(), f"MISSING {p}"
    s=p.stat().st_size; g=sha(p); tot+=s
    ok = (s==a["size_bytes"] and g==a["sha256"])
    print(f"  {'ok ' if ok else 'BAD'} {s:>13,}  {g[:16]}...  {a['path']}")
    if not ok: bad.append(a["path"])
print(f"  total weight bytes: {tot:,}")

print("\n=== EXPECTED_SHA256SUMS entries re-derived ===")
n=0; miss=0
for line in (A/"EXPECTED_SHA256SUMS").read_text().splitlines():
    if not line.strip(): continue
    exp,rel=line.split(None,1); rel=rel.strip()
    p=A/"source"/rel
    if not p.is_file(): print(f"  MISSING {rel}"); miss+=1; continue
    g=sha(p); n+=1
    if g!=exp: bad.append(rel); print(f"  BAD  {rel}")
print(f"  {n} entries re-derived, {miss} missing, {len(bad)} mismatches")
if bad: raise SystemExit("digest mismatches: "+", ".join(bad))

print("\n=== licence position (the point in question) ===")
print(f"  bundle 'license'        = {d.get('license')}")
print(f"  bundle 'weight_license' = {d.get('weight_license','<KEY ABSENT>')}")
print(f"  bundle 'code_license'   = {d.get('code_license','<KEY ABSENT>')}")
print(f"  weight_content_downloaded = {d['weight_content_downloaded']}  (bodies present -> stale)")

print("\n=== duplicate copies in abandoned staging dirs ===")
for s in ("geneformer-acquisition-21014704.staging","geneformer-acquisition-21014802.staging"):
    sp=Path("executions")/s
    if sp.is_dir():
        w=[q for q in sp.rglob("*.safetensors")]
        print(f"  {s}: {len(w)} weight bodies, {sum(q.stat().st_size for q in w):,} bytes")
print("\nNOTHING WAS MODIFIED OR DELETED BY THIS JOB")
PY

pip freeze --all > "${STAGING}/validation/pip-freeze.txt" 2>&1 || true
sha256sum "${ROOT}/config/artifacts/models/geneformer/checkpoints.json" \
          "${ROOT}/config/artifacts/models/geneformer/native_task_disposition_20260825.json" \
          "${ROOT}/slurm/verify_geneformer_weight_state.sbatch" > "${STAGING}/source.sha256"
PYTHONPATH=${ROOT}/src python - <<PY
from pathlib import Path
from masld_bench.artifacts import freeze_tree
print("ARTIFACTS.json sha256:", freeze_tree(Path("${STAGING}"), metadata={
    "artifact_class":"geneformer_weight_state_verification_read_only",
    "files_modified":0,"weights_deleted":False,"correction_applied":False,
    "sealed_outcomes_read":False,"models_fit":False,"models_scored":False}))
PY
mv "${STAGING}" "${TARGET}"; echo "FROZEN: ${TARGET}"
