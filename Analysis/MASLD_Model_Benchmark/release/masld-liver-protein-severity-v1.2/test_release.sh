#!/usr/bin/env bash
# Smoke test for masld-liver-protein-severity-v1, through the released code path only.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# A default that WORKS. The first draft defaulted to python3, which has no pandas on
# the build machine, so a stranger's first run died on ModuleNotFoundError. Override
# with PY= or PYTHON=; if neither resolves to a usable interpreter this FAILS LOUDLY
# with an actionable message instead of a bare traceback.
DEFAULT_PY="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
PY="${PY:-${PYTHON:-}}"
if [ -z "$PY" ]; then
  if [ -x "$DEFAULT_PY" ]; then PY="$DEFAULT_PY"; else PY="python3"; fi
fi
if ! "$PY" -c 'import numpy, scipy, pandas' >/dev/null 2>&1; then
  echo "ERROR: interpreter '$PY' is missing numpy, scipy or pandas." >&2
  echo "       Re-run as:  PY=/path/to/python bash test_release.sh" >&2
  exit 1
fi
echo "== interpreter: $PY"
DEPOSIT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/PXD051911"
echo "== building the training matrix from the deposit"
"$PY" - "$DEPOSIT" "$HERE" <<'EOF'
import sys, pandas as pd
dep, here = sys.argv[1], sys.argv[2]
L = pd.read_csv(f"{dep}/liver_protein_quant.txt", sep="	")
M = pd.read_csv(f"{dep}/meta_data.txt", sep="	", dtype=str)
LIV = M[M.liver_proteomics_filename.notna() & (M.liver_proteomics_filename != "NA")]
X = L.set_index("ProteinAccessions").iloc[:, 2:][list(LIV.liver_proteomics_filename)]
X.columns = list(LIV.patient_name)
X.to_csv(f"{here}/.smoke_matrix.tsv", sep="	", index_label="protein")
print("  wrote", X.shape[0], "proteins x", X.shape[1], "samples")
EOF
echo "== scoring through score.py"
"$PY" "$HERE/score.py" --abundance "$HERE/.smoke_matrix.tsv"       --out "$HERE/.smoke_scores.tsv" --json "$HERE/.smoke_report.json"
echo "== asserting the released path reproduces the fitted model, and that it refuses"
"$PY" - "$HERE" <<'EOF'
import hashlib, importlib.util, json, sys
sys.dont_write_bytecode = True
import numpy as np, pandas as pd
here = sys.argv[1]
exp = json.load(open(f"{here}/expected.json"))
sc = pd.read_csv(f"{here}/.smoke_scores.tsv", sep="	")
assert len(sc) == exp["n_samples"], f'{len(sc)} scores, expected {exp["n_samples"]}'
assert list(sc.sample_id.astype(str)) == exp["sample_order"], "sample order changed"
d = float(np.abs(sc.severity.to_numpy(float) - np.array(exp["in_sample_scores"])).max())
print(f"  max |released - fitted| = {d:.3e} (tolerance {exp['tolerance']:.0e})")
assert d < exp["tolerance"], f"released path differs from the fitted model by {d:.3e}"
rep = json.load(open(f"{here}/.smoke_report.json"))
assert rep["n_matched"] == exp["n_axis"], f'matched {rep["n_matched"]} != {exp["n_axis"]}'
assert abs(rep["axis_coverage"] - 1.0) < 1e-9, rep["axis_coverage"]
spec = importlib.util.spec_from_file_location("psc", f"{here}/score.py")
psc = importlib.util.module_from_spec(spec); spec.loader.exec_module(psc)
mat, ids, samples = psc.read_matrix(f"{here}/.smoke_matrix.tsv")
n_ref = 0
cases = [
    ("zero join", lambda: psc.score(mat, np.array([f"XX{i}" for i in range(len(ids))]))),
    ("low coverage", lambda: psc.score(mat[:200], ids[:200])),
    ("already logged", lambda: psc.score(np.log2(np.clip(mat, 1e-6, None)), ids)),
    ("negative values", lambda: psc.score(-np.abs(mat), ids)),
    ("min-coverage under the floor", lambda: psc.score(mat, ids, min_coverage=0.001)),
]
for name, fn in cases:
    try:
        fn(); print(f"  FAILED TO REFUSE: {name}")
    except psc.ScoreError:
        n_ref += 1; print(f"  refused: {name}")
assert n_ref == exp["n_refusals_expected"], f"{n_ref} refusals, expected {exp['n_refusals_expected']}"
print("")
print(f"OK  {len(sc)} samples scored, max diff {d:.3e}, {n_ref}/{len(cases)} refusals fired")
print(f"    leave-one-out Spearman (from the build, not recomputed here): {exp['loo_spearman']:.4f}")
EOF
rm -f "$HERE/.smoke_matrix.tsv" "$HERE/.smoke_scores.tsv" "$HERE/.smoke_report.json"
rm -rf "$HERE/__pycache__"
