#!/usr/bin/env bash
# Smoke test for encoder-benchmark-v1.
#
# Runs the released evaluator over the released reference, THROUGH THE RELEASED
# PATH ONLY, and asserts three things:
#
#   1. it reproduces the published Panel 7F paired deltas;
#   2. the four encoders with no clean held-out study are WITHHELD, and asking
#      for a number from them RAISES rather than returning one;
#   3. leaderboard.tsv contains none of them, and not_comparable.tsv has no
#      numeric column at all.
#
#   bash test_release.sh
#
# Expected values are read from expected_benchmark.json, written by the build.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PY="${PY:-${BENCH}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python}"
TMP="$(mktemp -d)"
trap 'rm -rf "${TMP}"' EXIT

echo "== encoder-benchmark-v1 smoke test"
echo "   python:  ${PY}"

for f in "${HERE}/evaluate.py" "${HERE}/exposure.py" "${HERE}/metrics.py" \
         "${HERE}/expected_benchmark.json" "${HERE}/reference/predictions.npz" \
         "${HERE}/reference/row_contract.tsv" "${HERE}/reference/exposure_declaration.json"; do
  [ -e "${f}" ] || { echo "FAIL: missing ${f}"; exit 1; }
done

# ------------------------------------------------------------ 1. released path
echo "-- running evaluate.py against the liver-specific encoder, both heads"
for head in linear two_layer_mlp; do
  "${PY}" "${HERE}/evaluate.py" \
    --predictions "${HERE}/reference/predictions.npz" \
    --row-contract "${HERE}/reference/row_contract.tsv" \
    --exposure "${HERE}/reference/exposure_declaration.json" \
    --reference "scvi_liver_latent::${head}" \
    --out "${TMP}/${head}" > "${TMP}/${head}.log"
  echo "   ${head}: $(grep -c . "${TMP}/${head}/leaderboard.tsv") leaderboard lines"
done

# ------------------------------------------------------------- 2. the numbers
echo "-- asserting the published deltas"
"${PY}" - "${TMP}" "${HERE}/expected_benchmark.json" <<'PYEOF'
import json, sys, csv
tmp, expected_path = sys.argv[1], sys.argv[2]
exp = json.load(open(expected_path))
got = {}
for head in ("linear", "two_layer_mlp"):
    with open(f"{tmp}/{head}/leaderboard.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            got[row["model_id"]] = row
fail = []
for model, e in exp["expected_deltas"].items():
    if model not in got:
        fail.append(f"{model} is absent from the leaderboard")
        continue
    d = float(got[model]["delta"])
    if abs(d - e["delta"]) > exp["tolerance"]:
        fail.append(f"{model} delta {d:.9f} != {e['delta']:.9f}")
    if got[model]["verdict"] != e["verdict"]:
        fail.append(f"{model} verdict {got[model]['verdict']} != {e['verdict']}")
    print(f"   {model:<44} {d:+.6f}  {got[model]['verdict']}")
if fail:
    print("FAIL:\n  " + "\n  ".join(fail)); sys.exit(1)
print("   PASS")
PYEOF

# ------------------------------------------- 3. the withheld models are withheld
echo "-- asserting the exposure refusal"
"${PY}" - "${TMP}" "${HERE}/expected_benchmark.json" <<'PYEOF'
import json, sys, csv
tmp, expected_path = sys.argv[1], sys.argv[2]
exp = json.load(open(expected_path))
must_be_withheld = set(exp["withheld_for_exposure"])
fail = []
for head in ("linear", "two_layer_mlp"):
    lb = {r["model_id"] for r in csv.DictReader(open(f"{tmp}/{head}/leaderboard.tsv"),
                                                delimiter="\t")}
    with open(f"{tmp}/{head}/not_comparable.tsv") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        nc = {r["model_id"] for r in rdr}
        cols = rdr.fieldnames
    skipped = {r["model_id"] for r in csv.DictReader(
        open(f"{tmp}/{head}/skipped_for_head.tsv"), delimiter="\t")}
    # Under the default same-head policy this run evaluates only models carrying
    # `head`; the other head's variant is SKIPPED FOR HEAD, not withheld for
    # exposure. Those are different reasons and the release keeps them apart.
    this_head = {m for m in must_be_withheld if m.endswith(f"::{head}")}
    other_head = must_be_withheld - this_head
    leaked = must_be_withheld & lb
    if leaked:
        fail.append(f"{head}: {sorted(leaked)} appear in the LEADERBOARD")
    absent = this_head - nc
    if absent:
        fail.append(f"{head}: {sorted(absent)} are not in not_comparable.tsv")
    misfiled = other_head - skipped
    if misfiled:
        fail.append(f"{head}: {sorted(misfiled)} are neither evaluated nor skipped for head")
    numeric = [c for c in cols if c in ("delta", "lower", "upper", "median", "estimate")]
    if numeric:
        fail.append(f"{head}: not_comparable.tsv carries numeric columns {numeric}")
    print(f"   {head}: withheld for exposure {sorted(this_head)}; "
          f"{len(skipped)} skipped for a differing head")
print(f"   not_comparable.tsv columns: {cols}")
if fail:
    print("FAIL:\n  " + "\n  ".join(fail)); sys.exit(1)
print("   PASS")
PYEOF

# ------------------------------- 4. asking a withheld model for a number must RAISE
echo "-- asserting that a withheld model REFUSES to produce a number"
"${PY}" - "${HERE}" <<'PYEOF'
import sys
sys.path.insert(0, sys.argv[1])
from exposure import (ComparableDelta, ExposureLedger, NotComparable,
                      NotComparableError, NoCleanStudyError, leaderboard)
nc = NotComparable(model_id="geneformer_v2_316m::two_layer_mlp", reference_id="ref",
                   head_id="two_layer_mlp", metric="donor_class_balanced_macro_f1",
                   reason="no_clean_held_out_study_exists_for_this_encoder",
                   confounded_studies=("GSE136103", "GSE185477", "Liver_Atlas"),
                   unresolved_studies=("GSE174748", "GSE189600", "GSE202379", "GSE244832"))
fail = []
for label, fn, want in (
    ("reading .delta", lambda: nc.delta, NotComparableError),
    ("reading .lower", lambda: nc.lower, NotComparableError),
    ("reading .estimate", lambda: nc.estimate, NotComparableError),
    ("ranking it", lambda: leaderboard([nc]), NotComparableError),
    ("building a delta with no clean study",
     lambda: ComparableDelta(model_id="x", reference_id="r", head_id="h", metric="m",
                             clean_studies=(), confounded_studies=(), unresolved_studies=(),
                             delta=0.9, lower=0.8, upper=1.0, median=0.9,
                             probability_model_better=1.0, n_resamples=1), NoCleanStudyError),
    ("a ledger that drops a study",
     lambda: ExposureLedger(model_id="x", clean=("A",), confounded=(), unresolved=(),
                            roster=("A", "B")), Exception),
):
    try:
        fn()
        fail.append(f"{label} did NOT raise")
    except want:
        print(f"   refused: {label}")
    except Exception as e:
        fail.append(f"{label} raised {type(e).__name__}, expected {want.__name__}")
# and a record of a withheld model must carry no numeric key
rec = nc.to_record()
bad = [k for k in rec if k in ("delta", "lower", "upper", "median", "estimate")]
if bad:
    fail.append(f"to_record() carries numeric keys {bad}")
print(f"   record keys: {sorted(rec)}")
if fail:
    print("FAIL:\n  " + "\n  ".join(fail)); sys.exit(1)
print("   PASS")
PYEOF

echo
echo "== ALL CHECKS PASSED"
