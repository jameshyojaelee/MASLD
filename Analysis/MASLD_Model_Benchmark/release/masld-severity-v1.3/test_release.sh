#!/usr/bin/env bash
# Release test for the masld-pairwise-conformal-v1 CANDIDATE package (B1, model-program-20261015T223347Z).
#
# Scores GSE213621 (confirmatory cohort) and GSE281797 from deposited counts THROUGH score.py ONLY,
# forms within-cohort pairs, joins the deposited stage labels, and asserts that coverage, call rate
# and call accuracy at the sealed 0.90 threshold equal the B1 result (RESULT.md, results.json):
# GSE213621 to 1e-6 (wave-3 ruling), GSE281797 to 1e-9 (PRESPEC section 11). Integer pair counts must
# be equal. Nothing here imports B1, RC or SUB code.
#
#   bash test_release.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
BENCH=${ROOT}/Analysis/MASLD_Model_Benchmark
G213_COUNTS=${BENCH}/executions/external-validation-expansion-20260901T190000Z/GSE213621_axis_counts.tsv
G213_META=${ROOT}/Analysis/Deconvolution/bulk/GSE213621/GSE213621_metadata.tsv
G281_COUNTS=${BENCH}/executions/nas-endpoint-20260902T234500Z/matrices/GSE281797_axis_counts.tsv
G281_PH=${ROOT}/data/GSE281797/metadata/sample_phenotypes_resolved.tsv
G225_COUNTS=${BENCH}/executions/nas-endpoint-20260902T234500Z/matrices/GSE225740_axis_counts.tsv
PY="${PY:-${BENCH}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python}"
# Working files go to TEST_WORKDIR (default: a new mktemp directory). They are left in place.
TMP="${TEST_WORKDIR:-$(mktemp -d)}"
mkdir -p "${TMP}"
echo "   workdir: ${TMP}"

echo "== masld-pairwise-conformal-v1 candidate release test"
echo "   package: ${HERE}"
echo "   python:  ${PY}"
echo "   host:    $(hostname)  job ${SLURM_JOB_ID:-none}  $(date -u +%FT%TZ)"

for f in "${HERE}/score.py" "${HERE}/weights/b1_pairwise_refit_weights.npz" "${HERE}/thresholds.json" \
         "${HERE}/expected.json" "${HERE}/MANIFEST.sha256.json" \
         "${HERE}/test_fixtures/GSE213621_eval_units.tsv" "${HERE}/test_fixtures/GSE281797_eval_units.tsv" \
         "${HERE}/test_fixtures/reference_scores.tsv" \
         "${G213_COUNTS}" "${G213_META}" "${G281_COUNTS}" "${G281_PH}" "${G225_COUNTS}"; do
  [ -e "${f}" ] || { echo "FAIL: missing ${f}"; exit 1; }
done

# ---------------------------------------------------------------- 0. manifest
echo "-- verifying MANIFEST.sha256.json"
"${PY}" - "${HERE}" <<'PYEOF'
import hashlib, json, sys
here = sys.argv[1]
m = json.load(open(f"{here}/MANIFEST.sha256.json"))
bad = []
for rel, rec in m["files"].items():
    h = hashlib.sha256(open(f"{here}/{rel}", "rb").read()).hexdigest()
    if h != rec["sha256"]:
        bad.append(rel)
if bad:
    print(f"   FAIL: sha256 differs for {bad}"); sys.exit(1)
print(f"   PASS: {len(m['files'])} files match")
PYEOF

# ---------------------------------------------------------------- 1. caller-shaped inputs
echo "-- writing caller-shaped counts TSVs with VERSIONED gene ids"
"${PY}" - "${G213_COUNTS}" "${G281_COUNTS}" "${G225_COUNTS}" "${TMP}" <<'PYEOF'
import sys, pandas as pd
g213, g281, g225, tmp = sys.argv[1:5]
for name, path in (("GSE213621", g213), ("GSE281797", g281), ("GSE225740", g225)):
    d = pd.read_csv(path, sep="\t", index_col=0)
    d.index = [f"{str(i)}.5" for i in d.index]
    d.index.name = "gene_id"
    d.to_csv(f"{tmp}/{name}_counts.tsv", sep="\t")
    print(f"   {name}: {d.shape[0]} genes x {d.shape[1]} samples")
PYEOF

# ---------------------------------------------------------------- 2. the released path
for c in GSE213621 GSE281797 GSE225740; do
  echo "-- score.py on ${c}"
  "${PY}" "${HERE}/score.py" --counts "${TMP}/${c}_counts.tsv" --out "${TMP}/${c}_scores.tsv" \
    --pairs "${TMP}/${c}_pairs.tsv" --json "${TMP}/${c}_report.json" 2>&1 | sed 's/^/   /'
done

# ---------------------------------------------------------------- 3. assertions
echo "-- asserting coverage, call rate and accuracy against expected.json"
"${PY}" - "${HERE}" "${TMP}" "${G213_META}" "${G281_PH}" <<'PYEOF'
import json, sys, numpy as np, pandas as pd
here, tmp, meta213, ph281 = sys.argv[1:5]
exp = json.load(open(f"{here}/expected.json"))
fail = []

def evaluate(cohort, units, label):
    pr = pd.read_csv(f"{tmp}/{cohort}_pairs.tsv", sep="\t", dtype={"sample_i": str, "sample_j": str})
    u = set(units)
    pr = pr[pr.sample_i.isin(u) & pr.sample_j.isin(u)]
    n_units = len(u)
    if len(pr) != n_units * (n_units - 1) // 2:
        fail.append(f"{cohort}: {len(pr)} unit pairs, expected {n_units*(n_units-1)//2}")
    dy = pr.sample_i.map(label).values - pr.sample_j.map(label).values
    call = pr.call.values
    disc = dy != 0
    right = ((call == "i_above_j") & (dy > 0)) | ((call == "j_above_i") & (dy < 0))
    called = call != "abstain"
    n_disc = int(disc.sum()); n_call = int((called & disc).sum())
    n_cor = int((right & disc).sum()); n_cov = int(((~called | right) & disc).sum())
    return dict(n_units=n_units, n_discordant_pairs=n_disc, n_covered=n_cov, n_calls=n_call,
                n_correct=n_cor, coverage=n_cov / n_disc, call_rate=n_call / n_disc,
                call_accuracy=n_cor / n_call)

# GSE213621: genetic-individual roster (fixture), coarse classes, Control excluded
md = pd.read_csv(meta213, sep="\t", dtype=str, keep_default_na=False)
st = md.set_index("sample_id")["fibrotic_stage"]
U = pd.read_csv(f"{here}/test_fixtures/GSE213621_eval_units.tsv", sep="\t", dtype=str)
cls = {"F0F1": 0.0, "F2": 1.0, "F3F4": 2.0}
u213 = [r for r in U.run if st[r] in cls]
lab213 = {r: cls[st[r]] for r in u213}
# GSE281797: deposited participant key, sex_concordant == False excluded (fixture), fibrosis 0-4
ph = pd.read_csv(ph281, sep="\t", dtype=str, keep_default_na=False).set_index("sample_id")
u281 = pd.read_csv(f"{here}/test_fixtures/GSE281797_eval_units.tsv", sep="\t", dtype=str).sample_id.tolist()
lab281 = {s: float(ph.fibrosis[s]) for s in u281}
for cohort, units, lab in (("GSE213621", u213, lab213), ("GSE281797", u281, lab281)):
    e = exp[cohort]
    r = evaluate(cohort, units, lab)
    rep = json.load(open(f"{tmp}/{cohort}_report.json"))
    if rep["plane"] != e["plane"]:
        fail.append(f"{cohort}: plane {rep['plane']} != {e['plane']}")
    if abs(rep["q"] - e["q"]) > 0:
        fail.append(f"{cohort}: q {rep['q']!r} != sealed {e['q']!r}")
    if rep["n_versioned_ids_stripped"] != rep["n_input_genes"]:
        fail.append(f"{cohort}: version stripping did not fire on every id")
    for k in ("n_units", "n_discordant_pairs", "n_covered", "n_calls", "n_correct"):
        if r[k] != e[k]:
            fail.append(f"{cohort}: {k} {r[k]} != {e[k]}")
    for k in ("coverage", "call_rate", "call_accuracy"):
        if abs(r[k] - e[k]) > e["tolerance"]:
            fail.append(f"{cohort}: {k} {r[k]:.10f} != {e[k]:.10f} (tolerance {e['tolerance']})")
    print(f"   {cohort} plane {rep['plane']} ({rep['n_plane_genes']} genes) q {rep['q']:.6f}: "
          f"{r['n_units']} units, {r['n_discordant_pairs']} discordant pairs")
    print(f"      coverage      {r['coverage']:.10f}  expected {e['coverage']:.10f}")
    print(f"      call rate     {r['call_rate']:.10f}  expected {e['call_rate']:.10f}")
    print(f"      call accuracy {r['call_accuracy']:.10f}  expected {e['call_accuracy']:.10f}")
    print(f"      |diff| max    {max(abs(r[k]-e[k]) for k in ('coverage','call_rate','call_accuracy')):.3e}"
          f"  (tolerance {e['tolerance']}); integer counts {r['n_covered']}/{r['n_calls']}/{r['n_correct']}")

# per-sample scores against the refit's own evaluation scores (dual path in the export job)
ref = pd.read_csv(f"{here}/test_fixtures/reference_scores.tsv", sep="\t", dtype={"sample_id": str})
for cohort in ("GSE213621", "GSE281797", "GSE225740"):
    s = pd.read_csv(f"{tmp}/{cohort}_scores.tsv", sep="\t", dtype={"sample_id": str})
    m = s.merge(ref[ref.cohort == cohort], on="sample_id")
    dmax = float(np.max(np.abs(m.refit_score_x - m.refit_score_y)))
    if len(m) != len(s) or dmax > exp["reference_score_tolerance"]:
        fail.append(f"{cohort}: per-sample score max |diff| {dmax:.3e} over {len(m)}/{len(s)} samples")
    if not (m.plane_x == m.plane_y).all():
        fail.append(f"{cohort}: plane differs from the reference")
    print(f"   {cohort}: per-sample refit score vs reference max |diff| {dmax:.3e} ({len(m)} samples)")
if fail:
    print("FAIL:\n  " + "\n  ".join(fail)); sys.exit(1)
print("   PASS")
PYEOF

# ---------------------------------------------------------------- 4. one sample at a time
echo "-- asserting scores do not depend on the batch (subset of 5 samples scores identically)"
"${PY}" - "${TMP}" <<'PYEOF'
import sys, pandas as pd
tmp = sys.argv[1]
d = pd.read_csv(f"{tmp}/GSE281797_counts.tsv", sep="\t", index_col=0)
d.iloc[:, :5].to_csv(f"{tmp}/sub5.tsv", sep="\t")
PYEOF
"${PY}" "${HERE}/score.py" --counts "${TMP}/sub5.tsv" --out "${TMP}/sub5_scores.tsv" 2>/dev/null
"${PY}" - "${TMP}" <<'PYEOF'
import sys, numpy as np, pandas as pd
tmp = sys.argv[1]
a = pd.read_csv(f"{tmp}/GSE281797_scores.tsv", sep="\t").iloc[:5]
b = pd.read_csv(f"{tmp}/sub5_scores.tsv", sep="\t")
if not np.array_equal(a.refit_score.values, b.refit_score.values):
    print("   FAIL: scores changed with the batch"); sys.exit(1)
print("   PASS: bitwise identical")
PYEOF

# ---------------------------------------------------------------- 5. refusals
echo "-- asserting the guards refuse"
"${PY}" - "${TMP}" <<'PYEOF'
import sys, numpy as np, pandas as pd
tmp = sys.argv[1]
d = pd.read_csv(f"{tmp}/GSE281797_counts.tsv", sep="\t", index_col=0)
d.iloc[10:].to_csv(f"{tmp}/noplane.tsv", sep="\t")                     # 26,619 genes: no plane
d.head(500).rename(index={i: f"SYM{n}" for n, i in enumerate(d.index[:500])}).to_csv(
    f"{tmp}/nojoin.tsv", sep="\t")
(d / d.sum(axis=0) * 1e6).to_csv(f"{tmp}/cpm.tsv", sep="\t")
np.log2(d / d.sum(axis=0) * 1e6 + 1).to_csv(f"{tmp}/log.tsv", sep="\t")
PYEOF
fails=0
for case in noplane nojoin cpm log; do
  if "${PY}" "${HERE}/score.py" --counts "${TMP}/${case}.tsv" --out "${TMP}/x.tsv" >/dev/null 2>"${TMP}/${case}.err"; then
    echo "   FAIL: ${case} was accepted"; fails=$((fails + 1))
  else
    echo "   refused ${case}: $(head -c 120 "${TMP}/${case}.err" | tr '\n' ' ')"
  fi
done
if "${PY}" "${HERE}/score.py" --counts "${TMP}/GSE281797_counts.tsv" --out "${TMP}/x.tsv" \
     --accession GSE193084 >/dev/null 2>"${TMP}/deny.err"; then
  echo "   FAIL: denylist did not fire"; fails=$((fails + 1))
else
  echo "   refused GSE193084: $(head -c 100 "${TMP}/deny.err" | tr '\n' ' ')"
fi
[ "${fails}" -eq 0 ] || { echo "RELEASE TEST FAILED: ${fails} guard(s) did not refuse"; exit 1; }

echo
echo "== ALL CHECKS PASSED"
