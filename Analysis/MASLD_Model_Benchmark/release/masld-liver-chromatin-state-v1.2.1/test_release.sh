#!/usr/bin/env bash
# End-to-end test for masld-liver-chromatin-state-v1.2.1.
# Asserts:
#  (1) eleven state columns are BITWISE v1.1's on the fixture;
#  (2) every v1.1 weight array is bitwise preserved;
#  (3) ALL 5 profile forms return their form-specific prof_reliable_<form> mask (states_cis=4847, rrr_offset_cis=4364, and distinct masks for rrr/rrr_cis/rrr_offset_states_cis);
#  (4) --profile on the 39 GSE296875 donors reproduces deposited transfer correlations to 1e-6;
#  (5) guards refuse invalid inputs;
#  (6) MANIFEST.sha256.json verifies.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
V11="${BENCH}/release/masld-liver-chromatin-state-v1.1"
FIX="${BENCH}/executions/model-data-064-21079902/fixture"
TL="${BENCH}/executions/chromatin-profile-transfer-gse296875-20260908T161209Z"
STB="$(ls -d ${BENCH}/executions/chromatin-stable-rrr-forms-* | tail -1)"
DON="${BENCH}/executions/gse296875-crossmodal-recurrence-20260830T193000Z/donor"
PY="${PY:-$HOME/micromamba/envs/rnaseq/bin/python}"
export PYTHONNOUSERSITE=1
TMP="$(mktemp -d)"
trap 'rm -rf "${TMP}"' EXIT

echo "== masld-liver-chromatin-state-v1.2.1 test"
for f in "${HERE}/score.py" "${HERE}/weights/chromatin_state_v1_2_1.npz" "${HERE}/weights/fit_report.json" \
         "${HERE}/weights/predictable_regions.bed" "${HERE}/MANIFEST.sha256.json" \
         "${V11}/weights/chromatin_state_v1_1.npz"; do
  [ -e "$f" ] || { echo "FAIL: missing $f"; exit 1; }
done

"${PY}" - "$FIX" "$TMP" "$DON" <<'PYEOF'
import sys, numpy as np, pandas as pd
fix, tmp, don = sys.argv[1:4]
r = np.load(f"{fix}/molecular/rna_values.npy")
g = pd.read_csv(f"{fix}/molecular/rna_feature_axis.tsv", sep="\t")["stable_gene_id"].astype(str)
p = pd.read_csv(f"{fix}/molecular/participant_axis.tsv", sep="\t")["participant_id"].astype(str)
d = pd.DataFrame(r.T, index=[x + ".9" for x in g], columns=p)
d.index.name = "gene_id"
d.to_csv(f"{tmp}/counts.tsv", sep="\t")
R = np.load(f"{don}/donor_rna_all.npy")
genes = pd.read_csv(f"{don}/rna_genes.tsv", sep="\t")
ax = pd.read_csv(f"{don}/donor_axis.tsv", sep="\t")
dd = pd.DataFrame(R.T, index=genes.ensembl_id.astype(str), columns=[f"donor_{x}" for x in ax.donor_id])
dd.index.name = "gene_id"
dd.to_csv(f"{tmp}/donors.tsv", sep="\t")
PYEOF

echo "-- states: v1.2.1 vs v1.1 on fixture + default rrr_offset_cis profile"
"${PY}" "${HERE}/score.py" --counts "${TMP}/counts.tsv" --out "${TMP}/s121.tsv" --json "${TMP}/rep121.json" \
  --profile "${TMP}/prof_fix.npz" --transport raw || { echo "FAIL: v1.2.1 score.py errored"; exit 1; }
"${PY}" "${V11}/score.py" --counts "${TMP}/counts.tsv" --out "${TMP}/s11.tsv" 2>/dev/null || { echo "FAIL: v1.1 score.py errored"; exit 1; }

echo "-- profile: all 5 forms on the 39 donors (marginal transport)"
for form in states_cis rrr_cis rrr rrr_offset_cis rrr_offset_states_cis; do
  "${PY}" "${HERE}/score.py" --counts "${TMP}/donors.tsv" --out "${TMP}/sd.tsv" \
    --profile "${TMP}/prof_${form}.npz" --profile-form "${form}" 2>/dev/null || { echo "FAIL: profile ${form} errored"; exit 1; }
done

"${PY}" - "${TMP}" "${HERE}" "${V11}" "${TL}" "${DON}" "${BENCH}" "${STB}" <<'PYEOF'
import os, sys, json, hashlib
import numpy as np, pandas as pd
from scipy.stats import rankdata

tmp, here, v11, tl, don, bench, stb = sys.argv[1:8]
fail = []
s121 = pd.read_csv(f"{tmp}/s121.tsv", sep="\t")
s11 = pd.read_csv(f"{tmp}/s11.tsv", sep="\t")
cols = [f"chromatin_state_k{k+1}" for k in range(10)] + ["steatosis_chromatin_axis"]
if not (s121[cols].to_numpy() == s11[cols].to_numpy()).all():
    fail.append("state/steatosis columns differ from v1.1 output")

W1 = np.load(f"{v11}/weights/chromatin_state_v1_1.npz", allow_pickle=True)
W21 = np.load(f"{here}/weights/chromatin_state_v1_2_1.npz", allow_pickle=True)
bad = [k for k in W1.files if k != "version" and not np.array_equal(W1[k], W21[k])]
if bad:
    fail.append(f"v1.1 arrays changed: {bad}")

# Assert all 5 forms emit their own prof_reliable_<form> mask
pf_fix = np.load(f"{tmp}/prof_fix.npz", allow_pickle=True)
if int(pf_fix["reliable"].sum()) != 4364:
    fail.append(f"default rrr_offset_cis reliable mask has {int(pf_fix['reliable'].sum())} regions, expected 4364")

seen_masks = {}
for form in ("states_cis", "rrr", "rrr_cis", "rrr_offset_cis", "rrr_offset_states_cis"):
    pf = np.load(f"{tmp}/prof_{form}.npz", allow_pickle=True)
    m = np.asarray(pf["reliable"], bool)
    expected = np.asarray(W21[f"prof_reliable_{form}"], bool)
    if not np.array_equal(m, expected):
        fail.append(f"{form}: emitted reliable mask does not match W['prof_reliable_{form}']")
    seen_masks[form] = m
    print(f"   {form:22s}: reliable={int(m.sum())} regions (matches prof_reliable_{form})")

if int(seen_masks["states_cis"].sum()) != 4847:
    fail.append("states_cis reliable count != 4847")
if np.array_equal(seen_masks["states_cis"], seen_masks["rrr_offset_cis"]):
    fail.append("states_cis and rrr_offset_cis masks must differ")

# Transfer correlation reproduction
per = np.load(f"{tl}/out/transfer_per_region.npz", allow_pickle=True)
regs, pk = per["region_index"], per["peak_index"]
A = np.load(f"{don}/donor_atac_all.npy")
lib = A.sum(1, keepdims=True)
Alog = np.log2(A / lib * 1e6 + 1.0)
def conc_block(raw):
    lb = raw.sum(1); ncol = raw.shape[1]; Pm = raw / lb[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -np.nansum(np.where(Pm > 0, Pm * np.log(Pm), 0.0), axis=1)
    sa = np.sort(raw, axis=1)
    gini = ((2 * np.arange(1, ncol + 1) - ncol - 1) * sa).sum(1) / (ncol * sa.sum(1))
    iqr = np.array([float(np.subtract(*np.percentile(np.log2(r[r > 0]), [75, 25]))) for r in raw])
    kt = int(round(0.05 * ncol)); top5 = (-np.sort(-raw, axis=1))[:, :kt].sum(1) / lb
    return np.column_stack([gini, ent, iqr, top5])
Zi = np.column_stack([np.ones(39), conc_block(A)])
b, *_ = np.linalg.lstsq(Zi, Alog, rcond=None)
ZY_raw = (Alog - Zi @ b)[:, pk]
def zrank(M):
    Rk = rankdata(M, axis=0); Rk = Rk - Rk.mean(0); sd = np.sqrt((Rk ** 2).sum(0)); return Rk / np.where(sd < 1e-12, 1, sd)
ZY = zrank(ZY_raw)
stab = np.load(f"{stb}/out/stable_transfer_rho.npz" if os.path.exists(f"{stb}/out/stable_transfer_rho.npz") else f"{stb}/out_finish/stable_transfer_rho.npz", allow_pickle=True)
for form, src, key in (
    ("states_cis", per, "rho_marginal_F1_states10_cis"),
    ("rrr_cis", stab, "rho_marginal_F3_rrr_cis_stable"),
    ("rrr", stab, "rho_marginal_F2_rrr_stable"),
    ("rrr_offset_cis", stab, "rho_marginal_rrr_offset_cis"),
    ("rrr_offset_states_cis", stab, "rho_marginal_rrr_offset_union"),
):
    P = np.load(f"{tmp}/prof_{form}.npz", allow_pickle=True)["profile"].astype(np.float64)[:, regs]
    rho = (zrank(P) * ZY).sum(0)
    d = float(np.abs(rho - src[key]).max())
    if d > 1e-6:
        fail.append(f"{form}: transfer rho max|diff| {d:.2e}")

if fail:
    print("FAIL:\n  " + "\n  ".join(fail))
    sys.exit(1)
print("   PASS: state bitwise v1.1, all 5 form-specific prof_reliable_<form> masks & transfer rho verified")
PYEOF
[ $? -eq 0 ] || exit 1

echo "-- asserting guards refuse"
"${PY}" - "${TMP}" <<'PYEOF'
import sys, numpy as np, pandas as pd
tmp = sys.argv[1]; d = pd.read_csv(f"{tmp}/counts.tsv", sep="\t", index_col=0)
d.head(500).rename(index={i: f"SYMBOL{n}" for n, i in enumerate(d.index[:500])}).to_csv(f"{tmp}/nojoin.tsv", sep="\t")
d.iloc[:int(0.30 * len(d))].to_csv(f"{tmp}/thin.tsv", sep="\t")
(d / d.sum(axis=0) * 1e6).to_csv(f"{tmp}/cpm.tsv", sep="\t")
np.log2(d + 1).to_csv(f"{tmp}/logged.tsv", sep="\t")
d.iloc[:, :3].to_csv(f"{tmp}/three.tsv", sep="\t")
PYEOF
fails=0
for case in nojoin thin cpm logged; do
  if "${PY}" "${HERE}/score.py" --counts "${TMP}/${case}.tsv" --out "${TMP}/x.tsv" >/dev/null 2>"${TMP}/${case}.err"; then
    echo "   FAIL: ${case} accepted"; fails=$((fails+1))
  fi
done
if "${PY}" "${HERE}/score.py" --counts "${TMP}/three.tsv" --out "${TMP}/x.tsv" --profile "${TMP}/x.npz" >/dev/null 2>"${TMP}/three.err"; then
  echo "   FAIL: 3-sample marginal accepted"; fails=$((fails+1))
fi
[ "${fails}" -eq 0 ] || { echo "TEST FAILED"; exit 1; }

echo "-- verifying MANIFEST.sha256.json"
"${PY}" - "${HERE}" <<'PYEOF'
import sys, json, hashlib, os
here = sys.argv[1]; m = json.load(open(f"{here}/MANIFEST.sha256.json"))
for f, v in m["files"].items():
    p = f"{here}/{f}"
    assert hashlib.sha256(open(p, "rb").read()).hexdigest() == v["sha256"], f"sha256 mismatch: {f}"
    assert os.path.getsize(p) == v["bytes"], f"bytes mismatch: {f}"
print(f"   {len(m['files'])} manifest files verify")
PYEOF
echo "== ALL CHECKS PASSED"
