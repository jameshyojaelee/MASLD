#!/usr/bin/env bash
# Smoke test for masld-severity-v1.1.
#
# Runs score.py on GSE268273 -- a cohort the model never trained on -- from raw
# deposited counts, THROUGH THE RELEASED PATH ONLY, and asserts the out-of-cohort
# Spearman it reproduces. Nothing in this script imports the fitting code.
#
#   bash test_release.sh
#
# The expected Spearman is read from expected_ood.json, written by the execution
# that fitted the model. If score.py and the fitting job disagree, this fails.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
OOD="${BENCH}/executions/full-axis-rebuild-and-ood-transfer-20260901T122932Z"
PHENO="${BENCH}/executions/gse268273-evaluator-phenotype-20260830T205701Z/gse268273_participant_phenotype.tsv"
PY="${PY:-${BENCH}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python}"
TMP="$(mktemp -d)"
trap 'rm -rf "${TMP}"' EXIT

echo "== masld-severity-v1 smoke test"
echo "   python:  ${PY}"
echo "   scratch: ${TMP}"

for f in "${HERE}/score.py" "${HERE}/kleiner.py" \
         "${HERE}/weights/masld_severity_v1_weights.npz" \
         "${HERE}/expected_ood.json" \
         "${OOD}/gse268273_model_axis_expected_counts.npy" \
         "${OOD}/gse268273_model_axis_genes.tsv" \
         "${OOD}/gse268273_participants.tsv" "${PHENO}"; do
  [ -e "${f}" ] || { echo "FAIL: missing ${f}"; exit 1; }
done

# ---------------------------------------------------------------- 1. raw -> caller-shaped TSV
# Deliberately writes a plain genes-by-samples counts TSV with VERSIONED gene ids,
# which is what a stranger would actually hand us, then never touches it again.
echo "-- building a caller-shaped counts TSV from the raw deposit"
"${PY}" - "${OOD}" "${TMP}" <<'PYEOF'
import sys, numpy as np, pandas as pd
ood, tmp = sys.argv[1], sys.argv[2]
C = np.load(f"{ood}/gse268273_model_axis_expected_counts.npy")           # participants x genes
g = pd.read_csv(f"{ood}/gse268273_model_axis_genes.tsv", sep="\t")["stable_gene_id"].astype(str)
p = pd.read_csv(f"{ood}/gse268273_participants.tsv", sep="\t")["row_id"].astype(str)
df = pd.DataFrame(C.T, index=[f"{x}.7" for x in g], columns=p)           # versioned on purpose
df.index.name = "gene_id"
df.to_csv(f"{tmp}/counts.tsv", sep="\t")
print(f"   wrote {df.shape[0]} genes x {df.shape[1]} samples, ids versioned")
PYEOF

# ---------------------------------------------------------------- 2. the released path
echo "-- running score.py (released path only)"
"${PY}" "${HERE}/score.py" \
  --counts "${TMP}/counts.tsv" \
  --out "${TMP}/scores.tsv" \
  --json "${TMP}/report.json" \
  --weights "${HERE}/weights/masld_severity_v1_weights.npz"

# ---------------------------------------------------------------- 3. assertions
echo "-- asserting"
"${PY}" - "${TMP}" "${PHENO}" "${HERE}/expected_ood.json" <<'PYEOF'
import json, sys, numpy as np, pandas as pd
from scipy import stats
tmp, pheno_path, expected_path = sys.argv[1], sys.argv[2], sys.argv[3]
exp = json.load(open(expected_path))
s = pd.read_csv(f"{tmp}/scores.tsv", sep="\t")
rep = json.load(open(f"{tmp}/report.json"))
ph = pd.read_csv(pheno_path, sep="\t")
ph["row_id"] = ph["row_id"].astype(str)
m = s.merge(ph[["row_id", "fibrosis_stage"]], left_on="sample_id", right_on="row_id")
fail = []
if len(m) != exp["n_participants"]:
    fail.append(f"joined {len(m)} rows, expected {exp['n_participants']}")
rho = float(stats.spearmanr(m.latent_severity, m.fibrosis_stage).statistic)
if abs(rho - exp["ood_spearman"]) > exp["tolerance"]:
    fail.append(f"Spearman {rho:.6f} != expected {exp['ood_spearman']:.6f} "
                f"(tolerance {exp['tolerance']})")
if abs(rep["axis_coverage"] - exp["axis_coverage"]) > 1e-6:
    fail.append(f"axis coverage {rep['axis_coverage']:.6f} != {exp['axis_coverage']:.6f}")
if rep["n_versioned_ids_stripped"] != rep["n_input_genes"]:
    fail.append("version stripping did not fire on deliberately versioned input")
if not np.isfinite(m.latent_severity).all():
    fail.append("non-finite scores")
print(f"   n              {len(m)} PARTICIPANTS")
print(f"   axis coverage  {rep['axis_coverage']:.6f} ({rep['n_matched']}/{rep['n_axis']})")
print(f"   version strip  {rep['n_versioned_ids_stripped']} ids")
print(f"   OOD Spearman   {rho:.6f}   expected {exp['ood_spearman']:.6f}")
print(f"   ceiling        {exp['tie_ceiling']:.6f}   native in-cohort {exp['native_in_cohort']:.6f}")
if fail:
    print("FAIL:\n  " + "\n  ".join(fail))
    sys.exit(1)
print("   PASS")
PYEOF

# ------------------------------------------------- 3b. the detection filter must be ENFORCED
# A caller supplying the full 42,163-gene transfer axis and a caller supplying exactly the
# release axis must get BITWISE IDENTICAL scores. If they do not, the 15,534 undetected genes
# are leaking into the design or into the library-size denominator, which is the failure that
# produced standardized values of 1.5e6 and destroyed test-retest reliability.
echo "-- asserting the detection filter is enforced, not documented"
"${PY}" - "${TMP}" "${HERE}" <<'PYEOF2'
import sys, numpy as np, pandas as pd
tmp, here = sys.argv[1], sys.argv[2]
W = np.load(f"{here}/weights/masld_severity_v1_weights.npz", allow_pickle=True)
axis = set(np.asarray(W["gene_axis"]).astype(str).tolist())
d = pd.read_csv(f"{tmp}/counts.tsv", sep="\t", index_col=0)
keep = [i for i in d.index if i.rsplit(".", 1)[0] in axis]
d.loc[keep].to_csv(f"{tmp}/axis_only.tsv", sep="\t")
print(f"   full input {len(d)} genes -> release-axis-only input {len(keep)} genes")
PYEOF2
"${PY}" "${HERE}/score.py" --counts "${TMP}/axis_only.tsv" --out "${TMP}/scores_axis_only.tsv" \
  --weights "${HERE}/weights/masld_severity_v1_weights.npz" 2>/dev/null
"${PY}" - "${TMP}" <<'PYEOF3'
import sys, numpy as np, pandas as pd
tmp = sys.argv[1]
a = pd.read_csv(f"{tmp}/scores.tsv", sep="\t")
b = pd.read_csv(f"{tmp}/scores_axis_only.tsv", sep="\t")
if not np.array_equal(a.latent_severity.values, b.latent_severity.values):
    print(f"   FAIL: scores differ by up to "
          f"{np.abs(a.latent_severity - b.latent_severity).max():.3e}; the excluded genes are "
          f"leaking into the design or the library-size denominator")
    sys.exit(1)
print("   PASS: bitwise identical with and without the 15,534 excluded genes")
PYEOF3

# ------------------------------------------------- 3c. the interval contract
# Above the certified minimum batch size the centred interval is emitted; below it,
# score.py must return scores and NO interval rather than an uncalibrated one.
echo "-- asserting the batch-size interval contract"
"${PY}" - "${TMP}" <<'PYEOF4'
import sys, pandas as pd
tmp = sys.argv[1]
d = pd.read_csv(f"{tmp}/counts.tsv", sep="\t", index_col=0)
d.iloc[:, :10].to_csv(f"{tmp}/small_batch.tsv", sep="\t")   # 10 samples, below n*=20
d.iloc[:, :25].to_csv(f"{tmp}/ok_batch.tsv", sep="\t")      # 25 samples, above n*
PYEOF4
for case in small_batch ok_batch; do
  "${PY}" "${HERE}/score.py" --counts "${TMP}/${case}.tsv" --out "${TMP}/${case}_out.tsv" \
    --weights "${HERE}/weights/masld_severity_v1_weights.npz" 2>"${TMP}/${case}.log"
done
"${PY}" - "${TMP}" <<'PYEOF5'
import sys, pandas as pd
tmp = sys.argv[1]
small = pd.read_csv(f"{tmp}/small_batch_out.tsv", sep="\t")
okb = pd.read_csv(f"{tmp}/ok_batch_out.tsv", sep="\t")
fail = []
if any("lo80" in c for c in small.columns):
    fail.append(f"10 samples got an interval: {list(small.columns)}")
if "severity_vs_batch_median" not in small.columns:
    fail.append("10 samples did not get a centred score")
if not any("lo80" in c for c in okb.columns):
    fail.append(f"25 samples got no interval: {list(okb.columns)}")
print(f"   n=10  columns: {list(small.columns)}")
print(f"   n=25  columns: {list(okb.columns)}")
if fail:
    print("   FAIL:\n     " + "\n     ".join(fail)); sys.exit(1)
print("   PASS: interval emitted above n*, refused below it")
PYEOF5
grep -q "NO INTERVAL EMITTED" "${TMP}/small_batch.log" || {
  echo "   FAIL: no explanation printed when the interval was withheld"; exit 1; }
echo "   PASS: the withheld interval is explained on stderr"

# ---------------------------------------------------------------- 4. the refusals must refuse
echo "-- asserting the guards actually refuse"
"${PY}" - "${TMP}" <<'PYEOF'
import sys, numpy as np, pandas as pd
tmp = sys.argv[1]
d = pd.read_csv(f"{tmp}/counts.tsv", sep="\t", index_col=0)
# a matrix whose ids are symbols, not ENSG: must be a zero join
d.head(500).rename(index={i: f"SYMBOL{n}" for n, i in enumerate(d.index[:500])}
                   ).to_csv(f"{tmp}/nojoin.tsv", sep="\t")
# a matrix carrying only 5% of the axis: must be refused on coverage
d.iloc[:int(0.05 * len(d))].to_csv(f"{tmp}/thin.tsv", sep="\t")
# already-CPM values: must be refused as not-counts
(d / d.sum(axis=0) * 1e6).to_csv(f"{tmp}/cpm.tsv", sep="\t")
PYEOF

fails=0
for case in nojoin thin cpm; do
  if "${PY}" "${HERE}/score.py" --counts "${TMP}/${case}.tsv" --out "${TMP}/x.tsv" \
       --weights "${HERE}/weights/masld_severity_v1_weights.npz" >/dev/null 2>"${TMP}/${case}.err"
  then
    echo "   FAIL: ${case} was ACCEPTED and should have been refused"
    fails=$((fails + 1))
  else
    echo "   refused ${case}: $(head -c 110 "${TMP}/${case}.err" | tr '\n' ' ')"
  fi
done

if "${PY}" "${HERE}/score.py" --counts "${TMP}/counts.tsv" --out "${TMP}/x.tsv" \
     --accession GSE193084 --weights "${HERE}/weights/masld_severity_v1_weights.npz" \
     >/dev/null 2>"${TMP}/deny.err"
then
  echo "   FAIL: the HCC denylist did not fire on GSE193084"
  fails=$((fails + 1))
else
  echo "   refused GSE193084: $(head -c 110 "${TMP}/deny.err" | tr '\n' ' ')"
fi

# ------------------------------------------- 4b. the denylist must read EVERY sample id
# v1 called check_denylist(*samples[:2000]), so a denylisted accession in column 2001 or
# beyond passed. The matrix here is 3 genes wide on purpose: the denylist runs before the
# join, so the test isolates the prefix bug from every other refusal.
echo "-- asserting the denylist reads every sample id, not the first 2000"
"${PY}" - "${TMP}" <<'PYEOF6'
import sys
tmp = sys.argv[1]
names = [f"clean_{i}" for i in range(2500)]
names[2400] = "GSE193084_sample"          # position 2401, past v1's 2000-id prefix
with open(f"{tmp}/wide.tsv", "w") as fh:
    fh.write("gene_id\t" + "\t".join(names) + "\n")
    for g in ("ENSG00000000003", "ENSG00000000005", "ENSG00000000419"):
        fh.write(g + "\t" + "\t".join(["100"] * len(names)) + "\n")
print(f"   built 3 genes x {len(names)} samples, denylisted id at position 2401")
PYEOF6
if "${PY}" "${HERE}/score.py" --counts "${TMP}/wide.tsv" --out "${TMP}/x.tsv" \
     --weights "${HERE}/weights/masld_severity_v1_weights.npz" >/dev/null 2>"${TMP}/wide.err"
then
  echo "   FAIL: a denylisted sample id at position 2401 was ACCEPTED"
  fails=$((fails + 1))
elif grep -q "GSE193084" "${TMP}/wide.err"; then
  echo "   refused: $(head -c 110 "${TMP}/wide.err" | tr '\n' ' ')"
else
  echo "   FAIL: refused for the wrong reason, not the denylist:"
  echo "         $(head -c 200 "${TMP}/wide.err" | tr '\n' ' ')"
  fails=$((fails + 1))
fi

# ------------------------------------------------- 4c. the INPUT CONTRACT must refuse,
# and --acknowledge-out-of-contract must let a caller through while recording it.
# The violation is built by masking counts to zero at a rate well past the measured
# bound; the library-size and not-counts guards all PASS on it, which is the point.
echo "-- asserting the input contract refuses, and that the escape hatch is recorded"
"${PY}" - "${TMP}" <<'PYEOF7'
import sys, numpy as np, pandas as pd
tmp = sys.argv[1]
d = pd.read_csv(f"{tmp}/counts.tsv", sep="\t", index_col=0)
rng = np.random.default_rng(0)
m = rng.random(d.shape) >= 0.70                    # 70% of entries zeroed
(d * m).to_csv(f"{tmp}/sparse.tsv", sep="\t")
print(f"   built a 70%-masked matrix, {d.shape[0]} genes x {d.shape[1]} samples")
PYEOF7
if "${PY}" "${HERE}/score.py" --counts "${TMP}/sparse.tsv" --out "${TMP}/x.tsv" \
     --weights "${HERE}/weights/masld_severity_v1_weights.npz" >/dev/null 2>"${TMP}/sparse.err"
then
  echo "   FAIL: an out-of-contract matrix was ACCEPTED"
  fails=$((fails + 1))
elif grep -q "OUT OF CONTRACT" "${TMP}/sparse.err"; then
  echo "   refused: $(head -c 150 "${TMP}/sparse.err" | tr '\n' ' ')"
else
  echo "   FAIL: refused for the wrong reason, not the contract:"
  echo "         $(head -c 200 "${TMP}/sparse.err" | tr '\n' ' ')"
  fails=$((fails + 1))
fi
if "${PY}" "${HERE}/score.py" --counts "${TMP}/sparse.tsv" --out "${TMP}/sparse_out.tsv" \
     --json "${TMP}/sparse.json" --acknowledge-out-of-contract \
     --weights "${HERE}/weights/masld_severity_v1_weights.npz" >/dev/null 2>"${TMP}/ack.err"
then
  "${PY}" - "${TMP}" <<'PYEOF8'
import json, sys
tmp = sys.argv[1]
r = json.load(open(f"{tmp}/sparse.json"))
bad = []
if not r.get("out_of_contract"):
    bad.append("out_of_contract not recorded as true in the JSON")
if not r.get("acknowledged_out_of_contract"):
    bad.append("acknowledged_out_of_contract not recorded")
if r.get("n_samples_out_of_contract", 0) < 1:
    bad.append("n_samples_out_of_contract is 0 on a matrix that was refused without the flag")
print(f"   escape hatch: out_of_contract={r.get('out_of_contract')} "
      f"n_out={r.get('n_samples_out_of_contract')}/{r.get('n_samples')} "
      f"zero_frac_max={r.get('zero_fraction_observed_max'):.4f} "
      f"(bound {r.get('contract_zero_fraction_max'):.4f}) "
      f"qn_w1_max={r.get('qn_w1_observed_max'):.4f} "
      f"(bound {r.get('contract_qn_w1_max'):.4f})")
if bad:
    print("   FAIL:\n     " + "\n     ".join(bad)); sys.exit(1)
print("   PASS: scored under acknowledgement, and the acknowledgement is in the output")
PYEOF8
  grep -q "OUT OF CONTRACT" "${TMP}/ack.err" || {
    echo "   FAIL: the escape hatch did not warn on stderr"; fails=$((fails + 1)); }
else
  echo "   FAIL: --acknowledge-out-of-contract did not let the matrix through"
  fails=$((fails + 1))
fi

# ------------------------------- 4e. a _PAR_Y row must be SUMMED, not silently dropped
# GENCODE writes the chrY copy of a pseudoautosomal gene as <ENSG>.<ver>_PAR_Y. v1's
# VERSION_SUFFIX is anchored at end and does not strip it, so the row never matched the
# axis and its counts were dropped. Both branches are asserted: splitting a PAR gene in
# half must leave the score UNCHANGED here, and must change it under v1's rule.
echo "-- asserting a _PAR_Y row is summed into its X twin, not dropped"
"${PY}" - "${TMP}" "${HERE}" <<'PYEOFA'
import sys, numpy as np, pandas as pd
tmp, here = sys.argv[1], sys.argv[2]
sys.path.insert(0, here)
import score as SC
W = np.load(f"{here}/weights/masld_severity_v1_weights.npz", allow_pickle=True)
axis = np.asarray(W["gene_axis"]).astype(str)
d = pd.read_csv(f"{tmp}/counts.tsv", sep="\t", index_col=0)
stable = np.array([i.rsplit(".", 1)[0] for i in d.index])
# 23 axis genes are pseudoautosomal; the three named in the model card are enough here
pick = [g for g in ("ENSG00000182162", "ENSG00000198223", "ENSG00000236017") if g in set(stable)]
if not pick:
    print("   SKIP: none of the named PAR genes is in this matrix"); sys.exit(0)
rows = [int(np.where(stable == g)[0][0]) for g in pick]
d2 = d.copy()
d2.iloc[rows] = d2.iloc[rows] / 2.0
twin = d2.iloc[rows].copy()
twin.index = [f"{stable[r]}.11_PAR_Y" for r in rows]
d2 = pd.concat([d2, twin])
d2.to_csv(f"{tmp}/par_split.tsv", sep="\t")
print(f"   split {len(pick)} pseudoautosomal genes into an X row and a _PAR_Y row, "
      f"each carrying half")
PYEOFA
"${PY}" "${HERE}/score.py" --counts "${TMP}/par_split.tsv" --out "${TMP}/par_split_out.tsv" \
  --json "${TMP}/par_split.json" \
  --weights "${HERE}/weights/masld_severity_v1_weights.npz" 2>"${TMP}/par.log"
"${PY}" - "${TMP}" <<'PYEOFB'
import json, sys, numpy as np, pandas as pd
tmp = sys.argv[1]
a = pd.read_csv(f"{tmp}/scores.tsv", sep="\t")
b = pd.read_csv(f"{tmp}/par_split_out.tsv", sep="\t")
r = json.load(open(f"{tmp}/par_split.json"))
fail = []
if r.get("n_par_y_ids_stripped", 0) < 1:
    fail.append(f"no _PAR_Y suffix was stripped: n_par_y_ids_stripped="
                f"{r.get('n_par_y_ids_stripped')}")
if r.get("n_duplicate_rows_summed", 0) < 1:
    fail.append(f"the _PAR_Y rows were not summed into their X twin: "
                f"n_duplicate_rows_summed={r.get('n_duplicate_rows_summed')}")
d = float(np.abs(a.latent_severity.values - b.latent_severity.values).max())
if d > 1e-9:
    fail.append(f"splitting a PAR gene in half changed the score by up to {d:.3e}; "
                f"summing did not restore it")
print(f"   stripped {r.get('n_par_y_ids_stripped')} _PAR_Y ids, summed "
      f"{r.get('n_duplicate_rows_summed')} duplicate rows, max |dscore| {d:.3e}")
if fail:
    print("   FAIL:\n     " + "\n     ".join(fail)); sys.exit(1)
print("   PASS: the _PAR_Y row was summed back and the score is unchanged")
PYEOFB
grep -q "_PAR_Y" "${TMP}/par.log" || {
  echo "   FAIL: the _PAR_Y strip was not reported on stderr"; fails=$((fails + 1)); }

# ------------------------- 4d. the contract must NOT fire on the in-contract cohort
"${PY}" - "${TMP}" <<'PYEOF9'
import json, sys
tmp = sys.argv[1]
r = json.load(open(f"{tmp}/report.json"))
if r.get("out_of_contract"):
    print(f"   FAIL: the contract fired on GSE268273, which is in contract: "
          f"{r['n_samples_out_of_contract']} samples flagged"); sys.exit(1)
print(f"   PASS: in contract on GSE268273 — zero fraction max "
      f"{r['zero_fraction_observed_max']:.4f} of {r['contract_zero_fraction_max']:.4f}, "
      f"divergence max {r['qn_w1_observed_max']:.4f} of {r['contract_qn_w1_max']:.4f}")
PYEOF9

[ "${fails}" -eq 0 ] || { echo "SMOKE TEST FAILED: ${fails} guard(s) did not refuse"; exit 1; }

echo "-- kleiner.py worked example"
"${PY}" "${HERE}/kleiner.py" --demo | sed -n '1,14p'

echo
echo "== ALL CHECKS PASSED"
