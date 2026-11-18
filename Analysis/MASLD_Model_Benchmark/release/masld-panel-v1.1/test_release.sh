#!/usr/bin/env bash
# End-to-end check for masld-panel-v1.1 (OOD Spearman, Cq identity, batch-NaN invariance, refusal contract, manifest).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH="$(cd "$HERE/../.." && pwd)"
PY="${MASLD_PY:-$BENCH/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python}"
export PYTHONNOUSERSITE=1
M="$BENCH/executions/standardized-metric-panel-20260902T130923Z/results/GSE268273_axis_counts.tsv"
P="$BENCH/executions/gse268273-evaluator-phenotype-20260830T205701Z/gse268273_participant_phenotype.tsv"
TMP="$(mktemp -d)"
trap 'rm -rf "${TMP}"' EXIT
[ -f "$M" ] || { echo "SKIP: evaluation matrix not present at $M"; exit 0; }

echo "== 1. score 109 samples from counts =="
"$PY" -s -E "$HERE/score.py" --counts "$M" --out "${TMP}/panel_scores.tsv"

echo "== 2. assert OOD Spearman, Cq identity, batch-NaN invariance, refusal, and manifest =="
"$PY" -s -E - "$HERE" "$M" "$P" <<'PYEOF'
import hashlib, json, os, sys
import numpy as np, pandas as pd

here, mpath, ppath = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, here)
sys.path.insert(0, str(__import__("pathlib").Path(here).parents[1] /
                       "executions/standardized-metric-panel-20260902T130923Z"))
from score import load_panel, read_table, score, ScoreError
from panel_metrics import spearman

exp = json.load(open(f"{here}/expected.json"))
panel = load_panel(f"{here}/weights/masld_panel_v1_1.npz")
genes, samples, vals = read_table(mpath)
s, rep = score(vals, genes, panel)
assert rep["pair_coverage_min"] == 1.0, rep

# Check no OR2I1 / olfactory / low-abundance genes in weights metadata
fit_meta = json.load(open(f"{here}/weights/fit_metadata.json"))
assert fit_meta["min_panel_median_cpm_pool_v2"] >= 5.0
assert fit_meta["max_panel_frac_zero_pool_v2"] <= 0.05
assert "OR2I1" not in set(pd.read_csv(f"{here}/weights/panel_pairs.tsv", sep="\t")[["symbol_a", "symbol_b"]].values.ravel())

ph = pd.read_csv(ppath, sep="\t")
m = pd.DataFrame({"sample_id": samples, "s": s}).merge(ph, left_on="sample_id", right_on="row_id")
y = m["fibrosis_stage"].astype(str).str.extract(r"(\d)")[0].astype(float).values
r = spearman(m["s"].values, y)
d = abs(r - exp["spearman_vs_kleiner_fibrosis"])
print(f"  GSE268273 Spearman {r:.6f} vs expected {exp['spearman_vs_kleiner_fibrosis']:.6f}  |diff| {d:.6f}")
assert d <= exp["tolerance"], f"FROZEN FIGURE MOVED by {d}"

need = sorted(set(panel["a"]) | set(panel["b"]))
idx = {g: i for i, g in enumerate(genes)}
sub = np.array([vals[idx[g]] for g in need])
cq = -np.log2(np.maximum(sub / sub.sum(0) * 1e6, 1e-12)) + 22.0
s2, _ = score(cq, need, panel, is_cq=True)
dd = float(np.abs(s - s2).max())
print(f"  counts vs synthetic-Cq: max |diff| {dd:.3e}")
assert dd < 1e-9, "delta-Cq identity BROKEN"

s3, _ = score(cq + 7.5, need, panel, is_cq=True)
do = float(np.abs(s2 - s3).max())
print(f"  global +7.5 cycle offset: max |diff| {do:.3e}")
assert do < 1e-9, "NOT invariant to a global Cq offset"

# Batch-NaN sample invariance test (N07): Patient A solo vs batched with Patient B/C having NaNs
a_cq_solo = cq[:, [0]]
s_a_solo, _ = score(a_cq_solo, need, panel, is_cq=True)

cq_batch = cq[:, :5].copy()
cq_batch[0, 1] = np.nan
cq_batch[1, 2] = np.nan
cq_batch[2, 3] = np.nan
s_cq_batch, rep_cq_batch = score(cq_batch, need, panel, is_cq=True)
drift_cq = float(abs(s_cq_batch[0] - s_a_solo[0]))
print(f"  batch-NaN Cq invariance (Patient A solo vs batched with NaN patients): |diff| {drift_cq:.3e}")
assert drift_cq < 1e-12, f"Patient A Cq score drifted by {drift_cq}"

# Also verify 100-patient Cq batch with sporadic NaNs across different genes in different samples
cq_100 = np.tile(cq[:, :1], (1, 100))
cq_100[0, 1] = np.nan
cq_100[1, 2] = np.nan
cq_100[2, 3] = np.nan
s_cq_100, _ = score(cq_100, need, panel, is_cq=True)
assert np.max(np.abs(s_cq_100[[0] + list(range(4, 100))] - s_a_solo[0])) < 1e-12

# Verify counts-mode batch-NaN invariance
a_cnt_solo = vals[:, [0]]
s_cnt_solo, _ = score(a_cnt_solo, genes, panel, is_cq=False)
cnt_batch = vals[:, :3].copy()
g_pos = genes.index(need[0])
cnt_batch[g_pos, 1] = np.nan
s_cnt_batch, _ = score(cnt_batch, genes, panel, is_cq=False)
drift_cnt = float(abs(s_cnt_batch[0] - s_cnt_solo[0]))
print(f"  batch-NaN counts invariance: |diff| {drift_cnt:.3e}")
assert drift_cnt < 1e-12, f"Patient A counts score drifted by {drift_cnt}"

# Verify refusal below 0.70 coverage on a single patient
cq_drop = cq[:, [0]].copy()
# Drop > 30% of pairs by setting enough needed genes to NaN
n_drop = 0
for g_i in range(len(need)):
    cq_drop[g_i, 0] = np.nan
    n_drop += 1
    try:
        s_try, rep_try = score(cq_drop, need, panel, is_cq=True)
    except ScoreError as e:
        print(f"  refusal fires after dropping {n_drop} gene(s) below 0.70 coverage: {str(e)[:80]}...")
        break
else:
    raise SystemExit("REFUSAL CONTRACT DID NOT FIRE when pair coverage dropped below 0.70")

# Verify MANIFEST.sha256.json
m_path = f"{here}/MANIFEST.sha256.json"
assert os.path.exists(m_path), "missing MANIFEST.sha256.json"
mf = json.load(open(m_path))
for rel, v in mf["files"].items():
    p = os.path.join(here, rel)
    assert hashlib.sha256(open(p, "rb").read()).hexdigest() == v["sha256"], f"sha256 mismatch: {rel}"
    assert os.path.getsize(p) == v["bytes"], f"bytes mismatch: {rel}"
print(f"  manifest verified: {len(mf['files'])} files")
print("ALL CHECKS PASSED")
PYEOF
