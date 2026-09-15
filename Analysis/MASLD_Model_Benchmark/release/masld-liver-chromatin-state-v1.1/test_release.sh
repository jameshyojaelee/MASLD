#!/usr/bin/env bash
# End-to-end test for masld-liver-chromatin-state-v1.1. Asserts: (1) it scores every fixture sample at full coverage with ten
# state columns plus the steatosis head; (2) the ten state columns are BITWISE what v1's score.py returns on the same input;
# (3) every v1 weight array is bitwise preserved; (4) fit_report numbers trace to the exact-recipe lane JSON and are NOT the
# in-sample fit; (5) the steatosis head fits its own all-99 target; (6) the guards refuse.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; BENCH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
V1="${BENCH}/release/masld-liver-chromatin-state-v1"; FIX="${BENCH}/executions/model-data-064-21079902/fixture"; LANEA="${BENCH}/executions/chromatin-state-exact-recipe-eval-20260907T230940Z"
PY="${PY:-$HOME/micromamba/envs/rnaseq/bin/python}"; export PYTHONNOUSERSITE=1; TMP="$(mktemp -d)"; trap 'rm -rf "${TMP}"' EXIT
echo "== masld-liver-chromatin-state-v1.1 test"
for f in "${HERE}/score.py" "${HERE}/weights/chromatin_state_v1_1.npz" "${HERE}/weights/fit_report.json" "${V1}/weights/chromatin_state_v1.npz"; do [ -e "$f" ] || { echo "FAIL: missing $f"; exit 1; }; done
"${PY}" - "$FIX" "$TMP" <<'PYEOF'
import sys, numpy as np, pandas as pd
fix, tmp = sys.argv[1], sys.argv[2]; r = np.load(f"{fix}/molecular/rna_values.npy")
g = pd.read_csv(f"{fix}/molecular/rna_feature_axis.tsv", sep="\t")["stable_gene_id"].astype(str); p = pd.read_csv(f"{fix}/molecular/participant_axis.tsv", sep="\t")["participant_id"].astype(str)
d = pd.DataFrame(r.T, index=[x + ".9" for x in g], columns=p); d.index.name = "gene_id"; d.to_csv(f"{tmp}/counts.tsv", sep="\t"); print(f"   {d.shape[0]} genes x {d.shape[1]} samples (versioned ids)")
PYEOF
echo "-- running v1.1 and v1 score.py on the same input"
"${PY}" "${HERE}/score.py" --counts "${TMP}/counts.tsv" --out "${TMP}/s11.tsv" --json "${TMP}/rep.json" || { echo "FAIL: v1.1 score.py errored"; exit 1; }
"${PY}" "${V1}/score.py" --counts "${TMP}/counts.tsv" --out "${TMP}/s1.tsv" 2>/dev/null || { echo "FAIL: v1 score.py errored"; exit 1; }
"${PY}" - "${TMP}" "${HERE}" "${V1}" "${LANEA}" "${BENCH}" <<'PYEOF'
import sys, json, numpy as np, pandas as pd
from scipy.stats import spearmanr
tmp, here, v1, lanea, bench = sys.argv[1:6]; fail = []
s11 = pd.read_csv(f"{tmp}/s11.tsv", sep="\t"); s1 = pd.read_csv(f"{tmp}/s1.tsv", sep="\t"); rep = json.load(open(f"{tmp}/rep.json")); fr = json.load(open(f"{here}/weights/fit_report.json")); la = json.load(open(f"{lanea}/out/exact_recipe_results.json"))
if s11.shape[0] != 99: fail.append(f"scored {s11.shape[0]} samples")
if abs(rep["axis_coverage"] - 1.0) > 1e-9: fail.append("coverage != 1")
cols = [f"chromatin_state_k{k+1}" for k in range(10)]
if list(s11.columns) != ["sample_id"] + cols + ["steatosis_chromatin_axis", "axis_coverage"]: fail.append(f"columns {list(s11.columns)}")
if not (s11[cols].to_numpy() == s1[cols].to_numpy()).all(): fail.append("ten state columns differ from v1 output")
W1 = np.load(f"{v1}/weights/chromatin_state_v1.npz", allow_pickle=True); W2 = np.load(f"{here}/weights/chromatin_state_v1_1.npz", allow_pickle=True)
bad = [k for k in W1.files if not np.array_equal(W1[k], W2[k])]
if bad: fail.append(f"v1 arrays changed: {bad}")
for k in range(10):
    a = fr["component_accuracy_out_of_fold_exact_recipe"][f"k{k+1}"]["oof_partial_given_histology_sex_rnadesc"]; b = la["exact_recipe"][f"k{k+1}"]["partial"]
    if abs(a - b) > 1e-12: fail.append(f"fit_report k{k+1} {a} does not trace to lane JSON {b}")
if max(fr["component_accuracy_out_of_fold_exact_recipe"][f"k{k+1}"]["oof_partial_given_histology_sex_rnadesc"] for k in range(10)) > 0.90: fail.append("reported OOF max > 0.90: is it the in-sample fit?")
if fr["oof_residual_sd_exact"] != la["oof_residual_sd_exact"]: fail.append("oof_residual_sd does not trace to lane JSON")
if abs(fr["steatosis_head"]["oof_partial_given_histology_sex_rnadesc"] - la["steatosis_head"]["partial_hist_sex_rnadesc"]) > 1e-12: fail.append("steatosis head number does not trace")
if fr["steatosis_head"]["in_sample_fit_means_nothing"] < 0.80: fail.append("steatosis head does not fit its own target")
# steatosis column reproduces the head applied to the fixture PC scores
sys.path.insert(0, f"{bench}/executions/chromatin-beyond-histology-20260907T161521Z"); import common as C0
S = C0.load_substrate(); X = C0.logcpm(S["rna_raw"]); Z = ((X - W2["pca_mean"]) @ W2["pca_components"]) / W2["pca_score_sd"]; ste = ((Z - W2["ste_mu"]) / W2["ste_sd"]) @ W2["ste_W"] + W2["ste_b"]
if np.abs(ste[:, 0] - s11["steatosis_chromatin_axis"].to_numpy()).max() > 1e-5: fail.append("steatosis column does not reproduce from weights")
print(f"   samples {s11.shape[0]}, coverage {rep['axis_coverage']:.4f}, 10 state columns bitwise = v1, {len(W1.files)} v1 arrays preserved")
print(f"   fit_report exact-recipe partials: k1 {fr['component_accuracy_out_of_fold_exact_recipe']['k1']['oof_partial_given_histology_sex_rnadesc']:+.3f} k2 {fr['component_accuracy_out_of_fold_exact_recipe']['k2']['oof_partial_given_histology_sex_rnadesc']:+.3f}; steatosis head {fr['steatosis_head']['oof_partial_given_histology_sex_rnadesc']:+.3f}; external verdict: {fr['external_rna_side_evaluation']['verdict']}")
if fail: print("FAIL:\n  " + "\n  ".join(fail)); sys.exit(1)
print("   PASS")
PYEOF
[ $? -eq 0 ] || exit 1
echo "-- asserting the guards refuse"
"${PY}" - "${TMP}" <<'PYEOF'
import sys, numpy as np, pandas as pd
tmp = sys.argv[1]; d = pd.read_csv(f"{tmp}/counts.tsv", sep="\t", index_col=0)
d.head(500).rename(index={i: f"SYMBOL{n}" for n, i in enumerate(d.index[:500])}).to_csv(f"{tmp}/nojoin.tsv", sep="\t"); d.iloc[:int(0.30*len(d))].to_csv(f"{tmp}/thin.tsv", sep="\t")
(d / d.sum(axis=0) * 1e6).to_csv(f"{tmp}/cpm.tsv", sep="\t"); np.log2(d + 1).to_csv(f"{tmp}/logged.tsv", sep="\t")
PYEOF
fails=0
for case in nojoin thin cpm logged; do
  if "${PY}" "${HERE}/score.py" --counts "${TMP}/${case}.tsv" --out "${TMP}/x.tsv" >/dev/null 2>"${TMP}/${case}.err"; then echo "   FAIL: ${case} was ACCEPTED"; fails=$((fails+1)); else echo "   refused ${case}: $(head -c 90 "${TMP}/${case}.err" | tr '\n' ' ')"; fi
done
[ "${fails}" -eq 0 ] || { echo "TEST FAILED: ${fails} guard(s) did not refuse"; exit 1; }
echo; echo "== ALL CHECKS PASSED"
