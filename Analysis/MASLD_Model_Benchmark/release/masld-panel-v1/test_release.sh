#!/usr/bin/env bash
# End-to-end check. Reproduces the number in expected.json from the frozen artifact.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH="$(cd "$HERE/../.." && pwd)"
PY="${MASLD_PY:-$BENCH/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python}"
export PYTHONNOUSERSITE=1
M="$BENCH/executions/standardized-metric-panel-20260902T130923Z/results/GSE268273_axis_counts.tsv"
P="$BENCH/executions/gse268273-evaluator-phenotype-20260830T205701Z/gse268273_participant_phenotype.tsv"
[ -f "$M" ] || { echo "SKIP: evaluation matrix not present at $M"; exit 0; }

echo "== 1. score 109 samples from counts =="
"$PY" -s -E "$HERE/score.py" --counts "$M" --out /tmp/_panel_scores.tsv

echo "== 2. assert the frozen figure, the Cq identity, and the refusal contract =="
"$PY" -s -E - "$HERE" "$M" "$P" <<'PYEOF'
import sys, json, numpy as np, pandas as pd
here, mpath, ppath = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, here)
sys.path.insert(0, str(__import__("pathlib").Path(here).parents[1] /
                       "executions/standardized-metric-panel-20260902T130923Z"))
from score import load_panel, read_table, score, ScoreError
from panel_metrics import spearman

exp = json.load(open(f"{here}/expected.json"))
panel = load_panel(f"{here}/weights/masld_panel_v1.npz")
genes, samples, vals = read_table(mpath)
s, rep = score(vals, genes, panel)
assert rep["pair_coverage"] == 1.0, rep

ph = pd.read_csv(ppath, sep="\t")
m = pd.DataFrame({"sample_id": samples, "s": s}).merge(ph, left_on="sample_id", right_on="row_id")
y = m["fibrosis_stage"].astype(str).str.extract(r"(\d)")[0].astype(float).values
r = spearman(m["s"].values, y)
d = abs(r - exp["spearman_vs_kleiner_fibrosis"])
print(f"  Spearman {r:.4f} vs expected {exp['spearman_vs_kleiner_fibrosis']:.4f}  |diff| {d:.6f}")
assert d <= exp["tolerance"], f"FROZEN FIGURE MOVED by {d}"

need = sorted(set(panel["a"]) | set(panel["b"]))
idx = {g: i for i, g in enumerate(genes)}
sub = np.array([vals[idx[g]] for g in need])
cq = -np.log2(np.maximum(sub / sub.sum(0) * 1e6, 1e-12)) + 22.0
s2, _ = score(cq, need, panel, is_cq=True)
dd = np.abs(s - s2).max()
print(f"  counts vs synthetic-Cq: max |diff| {dd:.3e}")
assert dd < 1e-9, "delta-Cq identity BROKEN"

s3, _ = score(cq + 7.5, need, panel, is_cq=True)
do = np.abs(s2 - s3).max()
print(f"  global +7.5 cycle offset: max |diff| {do:.3e}")
assert do < 1e-9, "NOT invariant to a global Cq offset"

keep = [g for g in genes if g != "ENSG00000164318"]      # EGFLAM, 10 of 20 denominators
kv = np.array([vals[idx[g]] for g in keep])
try:
    score(kv, keep, panel)
    raise SystemExit("REFUSAL CONTRACT DID NOT FIRE when EGFLAM was removed")
except ScoreError as e:
    print(f"  refusal fires without EGFLAM: {str(e)[:70]}...")

print("ALL CHECKS PASSED")
PYEOF
