#!/usr/bin/env bash
# End-to-end test for masld-liver-chromatin-state-v1.
#
# Builds a caller-shaped counts TSV from the deposited fixture (VERSIONED ids on purpose), runs
# score.py through the released path only, and asserts:
#   1. it scores every sample at full coverage and strips version suffixes,
#   2. the shipped heads reproduce the release's OWN all-99 target in sample -- a self-consistency
#      guard: an 81-PC ridge that cannot fit its own target is broken,
#   3. the guards actually refuse: zero join, thin coverage, already-normalised, already-logged,
#   4. fit_report.json reports the OUT-OF-FOLD accuracy, not the in-sample fit.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
FIX="${BENCH}/executions/model-data-064-21079902/fixture"
PY="${PY:-$HOME/micromamba/envs/rnaseq/bin/python}"
TMP="$(mktemp -d)"; trap 'rm -rf "${TMP}"' EXIT

echo "== masld-liver-chromatin-state-v1 test"
for f in "${HERE}/score.py" "${HERE}/weights/chromatin_state_v1.npz" "${HERE}/weights/fit_report.json"; do
  [ -e "$f" ] || { echo "FAIL: missing $f"; exit 1; }
done

echo "-- building a caller-shaped counts TSV (versioned ids on purpose)"
"${PY}" - "$FIX" "$TMP" <<'PYEOF'
import sys, numpy as np, pandas as pd
fix, tmp = sys.argv[1], sys.argv[2]
r = np.load(f"{fix}/molecular/rna_values.npy")
g = pd.read_csv(f"{fix}/molecular/rna_feature_axis.tsv", sep="\t")["stable_gene_id"].astype(str)
p = pd.read_csv(f"{fix}/molecular/participant_axis.tsv", sep="\t")["participant_id"].astype(str)
d = pd.DataFrame(r.T, index=[x + ".9" for x in g], columns=p); d.index.name = "gene_id"
d.to_csv(f"{tmp}/counts.tsv", sep="\t")
print(f"   {d.shape[0]} genes x {d.shape[1]} samples")
PYEOF

echo "-- running score.py"
"${PY}" "${HERE}/score.py" --counts "${TMP}/counts.tsv" --out "${TMP}/states.tsv" --json "${TMP}/rep.json" \
  || { echo "FAIL: score.py errored"; exit 1; }

echo "-- asserting output, coverage, self-consistency, and that the reported accuracy is out-of-fold"
"${PY}" - "${TMP}" "${HERE}" "${BENCH}" <<'PYEOF'
import sys, json, numpy as np, pandas as pd
from scipy.stats import spearmanr
tmp, here, bench = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, f"{bench}/executions/chromatin-beyond-histology-20260907T161521Z")
import common as C0
s = pd.read_csv(f"{tmp}/states.tsv", sep="\t"); rep = json.load(open(f"{tmp}/rep.json"))
fr = json.load(open(f"{here}/weights/fit_report.json"))
W = np.load(f"{here}/weights/chromatin_state_v1.npz", allow_pickle=True)
fail = []
if s.shape[0] != 99: fail.append(f"scored {s.shape[0]} samples, expected 99")
if abs(rep["axis_coverage"] - 1.0) > 1e-9: fail.append(f"axis coverage {rep['axis_coverage']}")
if not all(f"chromatin_state_k{k+1}" in s.columns for k in range(10)): fail.append("missing state columns")
if not np.isfinite(s.filter(like="chromatin_state").to_numpy()).all(): fail.append("non-finite states")
# rebuild the release's own all-99 measured target and check the shipped heads fit it
S = C0.load_substrate(); Craw = C0.logcpm(S["h3_raw"]); n = len(S["ids"]); raw = S["h3_raw"]
lib = raw.sum(1); ncol = raw.shape[1]; Pm = raw / lib[:, None]
with np.errstate(divide="ignore", invalid="ignore"):
    ent = -np.nansum(np.where(Pm > 0, Pm * np.log(Pm), 0.0), axis=1)
sa = np.sort(raw, axis=1); gini = ((2*np.arange(1, ncol+1)-ncol-1)*sa).sum(1)/(ncol*sa.sum(1))
iqr = np.array([float(np.subtract(*np.percentile(np.log2(r[r > 0]), [75, 25]))) for r in raw])
kt = int(round(0.05*ncol)); top5 = (-np.sort(-raw, axis=1))[:, :kt].sum(1)/lib
CONC = np.column_stack([gini, ent, iqr, top5])
Zi = np.column_stack([np.ones(n), CONC]); b, *_ = np.linalg.lstsq(Zi, Craw, rcond=None); Cres = Craw - Zi @ b
V = np.asarray(W["pls_V"], float); muc = np.asarray(W["chrom_mu"], float); sdc = np.asarray(W["chrom_sd"], float)
sm = np.asarray(W["score_mean"], float); ss = np.asarray(W["score_sd"], float)
TARGET = (((Cres - muc)/sdc) @ V - sm)/ss
rhos = [float(spearmanr(s[f"chromatin_state_k{k+1}"], TARGET[:, k]).correlation) for k in range(10)]
low = [(k+1, round(r, 3)) for k, r in enumerate(rhos) if r < 0.80]
if low: fail.append(f"heads do not fit their own target in sample: {low}")
oof = [fr["component_accuracy_out_of_fold"][f"k{k+1}"]["oof_partial_given_histology_sex_rnadesc"] for k in range(10)]
if max(oof) > 0.90:
    fail.append(f"fit_report's out-of-fold max is {max(oof)} -- implausibly high, is it the in-sample fit?")
if not all(fr["component_accuracy_out_of_fold"][f"k{k+1}"]["lower_bound_above_zero"] for k in range(10)):
    fail.append("a released component lacks an out-of-fold lower bound above zero")
print(f"   samples {s.shape[0]}, coverage {rep['axis_coverage']:.4f}, versioned ids handled")
print(f"   in-sample fit to own target: min {min(rhos):.3f}, max {max(rhos):.3f} (guard 0.80)")
print(f"   out-of-fold partials: max {max(oof):+.3f}, min {min(oof):+.3f}; all ten lower bounds > 0")
if fail:
    print("FAIL:\n  " + "\n  ".join(fail)); sys.exit(1)
print("   PASS")
PYEOF
[ $? -eq 0 ] || exit 1

echo "-- asserting the guards refuse"
"${PY}" - "${TMP}" <<'PYEOF'
import sys, numpy as np, pandas as pd
tmp = sys.argv[1]
d = pd.read_csv(f"{tmp}/counts.tsv", sep="\t", index_col=0)
d.head(500).rename(index={i: f"SYMBOL{n}" for n, i in enumerate(d.index[:500])}).to_csv(f"{tmp}/nojoin.tsv", sep="\t")
d.iloc[:int(0.30*len(d))].to_csv(f"{tmp}/thin.tsv", sep="\t")
(d / d.sum(axis=0) * 1e6).to_csv(f"{tmp}/cpm.tsv", sep="\t")
np.log2(d + 1).to_csv(f"{tmp}/logged.tsv", sep="\t")
PYEOF
fails=0
for case in nojoin thin cpm logged; do
  if "${PY}" "${HERE}/score.py" --counts "${TMP}/${case}.tsv" --out "${TMP}/x.tsv" >/dev/null 2>"${TMP}/${case}.err"; then
    echo "   FAIL: ${case} was ACCEPTED"; fails=$((fails+1))
  else
    echo "   refused ${case}: $(head -c 100 "${TMP}/${case}.err" | tr '\n' ' ')"
  fi
done
[ "${fails}" -eq 0 ] || { echo "TEST FAILED: ${fails} guard(s) did not refuse"; exit 1; }

echo
echo "== ALL CHECKS PASSED"
