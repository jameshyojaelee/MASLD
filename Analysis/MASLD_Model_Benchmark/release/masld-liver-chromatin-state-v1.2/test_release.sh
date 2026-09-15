#!/usr/bin/env bash
# End-to-end test for masld-liver-chromatin-state-v1.2. Asserts: (1) the eleven state columns are BITWISE v1.1's on the fixture;
# (2) every v1.1 weight array is bitwise preserved; (3) --profile on the 39 GSE296875 donors reproduces the transfer lane's deposited
# per-region correlations to 1e-6 for the primary, rrr_cis and rrr forms (the card's transfer numbers, end to end from the release code);
# (4) fit_report sources verify by sha256, the card skill traces to the exact-recipe lane, and no in-cohort number exceeds 0.30;
# (5) the guards refuse (v1.1's four, plus a 3-sample batch under marginal transport); (6) the manifest verifies.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; BENCH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
V11="${BENCH}/release/masld-liver-chromatin-state-v1.1"; FIX="${BENCH}/executions/model-data-064-21079902/fixture"
TL="${BENCH}/executions/chromatin-profile-transfer-gse296875-20260908T161209Z"; STB="$(ls -d ${BENCH}/executions/chromatin-stable-rrr-forms-* | tail -1)"; DON="${BENCH}/executions/gse296875-crossmodal-recurrence-20260830T193000Z/donor"
PY="${PY:-$HOME/micromamba/envs/rnaseq/bin/python}"; export PYTHONNOUSERSITE=1; TMP="$(mktemp -d)"; trap 'rm -rf "${TMP}"' EXIT
echo "== masld-liver-chromatin-state-v1.2 test"
for f in "${HERE}/score.py" "${HERE}/weights/chromatin_state_v1_2.npz" "${HERE}/weights/fit_report.json" "${HERE}/weights/predictable_regions.bed" "${HERE}/MANIFEST.sha256.json" "${V11}/weights/chromatin_state_v1_1.npz"; do [ -e "$f" ] || { echo "FAIL: missing $f"; exit 1; }; done
"${PY}" - "$FIX" "$TMP" "$DON" <<'PYEOF'
import sys, numpy as np, pandas as pd
fix, tmp, don = sys.argv[1:4]; r = np.load(f"{fix}/molecular/rna_values.npy")
g = pd.read_csv(f"{fix}/molecular/rna_feature_axis.tsv", sep="\t")["stable_gene_id"].astype(str); p = pd.read_csv(f"{fix}/molecular/participant_axis.tsv", sep="\t")["participant_id"].astype(str)
d = pd.DataFrame(r.T, index=[x + ".9" for x in g], columns=p); d.index.name = "gene_id"; d.to_csv(f"{tmp}/counts.tsv", sep="\t")
R = np.load(f"{don}/donor_rna_all.npy"); genes = pd.read_csv(f"{don}/rna_genes.tsv", sep="\t"); ax = pd.read_csv(f"{don}/donor_axis.tsv", sep="\t")
dd = pd.DataFrame(R.T, index=genes.ensembl_id.astype(str), columns=[f"donor_{x}" for x in ax.donor_id]); dd.index.name = "gene_id"; dd.to_csv(f"{tmp}/donors.tsv", sep="\t")
print(f"   fixture {d.shape[0]} genes x {d.shape[1]} samples; donors {dd.shape[0]} genes x {dd.shape[1]} samples")
PYEOF
echo "-- states: v1.2 vs v1.1 on the fixture (+ a raw-transport profile)"
"${PY}" "${HERE}/score.py" --counts "${TMP}/counts.tsv" --out "${TMP}/s12.tsv" --json "${TMP}/rep12.json" --profile "${TMP}/prof_fix.npz" --transport raw || { echo "FAIL: v1.2 score.py errored"; exit 1; }
"${PY}" "${V11}/score.py" --counts "${TMP}/counts.tsv" --out "${TMP}/s11.tsv" 2>/dev/null || { echo "FAIL: v1.1 score.py errored"; exit 1; }
echo "-- profile: the three forms on the 39 donors (marginal transport)"
for form in states_cis rrr_cis rrr rrr_offset_cis rrr_offset_states_cis; do
  "${PY}" "${HERE}/score.py" --counts "${TMP}/donors.tsv" --out "${TMP}/sd.tsv" --profile "${TMP}/prof_${form}.npz" --profile-form ${form} 2>/dev/null || { echo "FAIL: profile ${form} errored"; exit 1; }
done
"${PY}" - "${TMP}" "${HERE}" "${V11}" "${TL}" "${DON}" "${BENCH}" "${STB}" <<'PYEOF'
import sys, json, hashlib, numpy as np, pandas as pd
from scipy.stats import rankdata
tmp, here, v11, tl, don, bench, stb = sys.argv[1:8]; fail = []
s12 = pd.read_csv(f"{tmp}/s12.tsv", sep="\t"); s11 = pd.read_csv(f"{tmp}/s11.tsv", sep="\t")
cols = [f"chromatin_state_k{k+1}" for k in range(10)] + ["steatosis_chromatin_axis"]
if s12.shape[0] != 99: fail.append(f"scored {s12.shape[0]} fixture samples")
if not (s12[cols].to_numpy() == s11[cols].to_numpy()).all(): fail.append("state/steatosis columns differ from v1.1 output")
W1 = np.load(f"{v11}/weights/chromatin_state_v1_1.npz", allow_pickle=True); W2 = np.load(f"{here}/weights/chromatin_state_v1_2.npz", allow_pickle=True)
bad = [k for k in W1.files if k != "version" and not np.array_equal(W1[k], W2[k])]
if bad: fail.append(f"v1.1 arrays changed: {bad}")
pf = np.load(f"{tmp}/prof_fix.npz", allow_pickle=True)
if pf["profile"].shape != (99, 96460) or pf["profile"].dtype != np.float64: fail.append(f"fixture profile shape/dtype {pf['profile'].shape} {pf['profile'].dtype}")
if int(pf["reliable"].sum()) != 4847: fail.append(f"reliable mask has {int(pf['reliable'].sum())} regions, not 4847")
per = np.load(f"{tl}/out/transfer_per_region.npz", allow_pickle=True); regs = per["region_index"]; pk = per["peak_index"]
A = np.load(f"{don}/donor_atac_all.npy"); lib = A.sum(1, keepdims=True); Alog = np.log2(A / lib * 1e6 + 1.0)
def conc_block(raw):
    lib = raw.sum(1); ncol = raw.shape[1]; Pm = raw / lib[:, None]
    with np.errstate(divide="ignore", invalid="ignore"): ent = -np.nansum(np.where(Pm > 0, Pm * np.log(Pm), 0.0), axis=1)
    sa = np.sort(raw, axis=1); gini = ((2 * np.arange(1, ncol + 1) - ncol - 1) * sa).sum(1) / (ncol * sa.sum(1))
    iqr = np.array([float(np.subtract(*np.percentile(np.log2(r[r > 0]), [75, 25]))) for r in raw]); kt = int(round(0.05 * ncol)); top5 = (-np.sort(-raw, axis=1))[:, :kt].sum(1) / lib
    return np.column_stack([gini, ent, iqr, top5])
Zi = np.column_stack([np.ones(39), conc_block(A)]); b, *_ = np.linalg.lstsq(Zi, Alog, rcond=None); Y = (Alog - Zi @ b)[:, pk]
def zrank(M):
    Rk = rankdata(M, axis=0); Rk = Rk - Rk.mean(0); sd = np.sqrt((Rk ** 2).sum(0)); return Rk / np.where(sd < 1e-12, 1, sd)
ZY = zrank(Y)
import os
stab = np.load(f"{stb}/out/stable_transfer_rho.npz" if os.path.exists(f"{stb}/out/stable_transfer_rho.npz") else f"{stb}/out_finish/stable_transfer_rho.npz", allow_pickle=True)
for form, src, key in (("states_cis", per, "rho_marginal_F1_states10_cis"), ("rrr_cis", stab, "rho_marginal_F3_rrr_cis_stable"), ("rrr", stab, "rho_marginal_F2_rrr_stable"),
                       ("rrr_offset_cis", stab, "rho_marginal_rrr_offset_cis"), ("rrr_offset_states_cis", stab, "rho_marginal_rrr_offset_union")):
    P = np.load(f"{tmp}/prof_{form}.npz", allow_pickle=True)["profile"].astype(np.float64)[:, regs]; rho = (zrank(P) * ZY).sum(0); d = float(np.abs(rho - src[key]).max())
    print(f"   {form}: per-region transfer rho reproduced from the release code, max|diff| {d:.2e} (mean rho {rho.mean():+.4f})")
    if d > 1e-6: fail.append(f"{form}: transfer rho not reproduced (max|diff| {d:.2e})")
fr = json.load(open(f"{here}/weights/fit_report.json"))
ex = fr["profile_head"]["in_cohort_exact_shipped_recipe"]
if ex["states_cis_fixed_100kb_10_states"]["skill"] > 0.30: fail.append("in-cohort skill > 0.30: is it the in-sample fit?")
if abs(ex["states_cis_fixed_100kb_10_states"]["skill"] - json.load(open(f"{bench}/{fr['sources']['exact_f1']['path']}"))["skill"]["fixed"]["skill"]) > 1e-12: fail.append("card skill does not trace to the exact-recipe lane")
card = open(f"{here}/MODEL_CARD.md").read()
for num in (f"{ex['states_cis_fixed_100kb_10_states']['skill']:+.4f}", f"{fr['profile_head']['external_transfer_gse296875']['marginal']['F1_states10_cis']['summaries']['mean_reliable']:+.4f}"):
    if num not in card: fail.append(f"card lacks the number {num} printed in fit_report")
bed = pd.read_csv(f"{here}/weights/predictable_regions.bed", sep="\t")
W2r = np.asarray(W2["prof_reliable"], bool); W2o = np.asarray(W2["prof_reliable_offset_form"], bool)
if len(bed) != int((W2r | W2o).sum()) or int(bed["reliable_states_cis"].sum()) != 4847 or (bed.loc[bed["reliable_states_cis"], [c for c in bed.columns if c.startswith("skill_fold")]].to_numpy() <= 0.2).any(): fail.append("predictable_regions.bed does not match the two reliable masks")
print(f"   states bitwise v1.1; {len(W1.files)} v1.1 arrays preserved; fixture profile 99 x 96460 float64; reliable 4847")
if fail: print("FAIL:\n  " + "\n  ".join(fail)); sys.exit(1)
print("   PASS")
PYEOF
[ $? -eq 0 ] || exit 1
echo "-- asserting the guards refuse"
"${PY}" - "${TMP}" <<'PYEOF'
import sys, numpy as np, pandas as pd
tmp = sys.argv[1]; d = pd.read_csv(f"{tmp}/counts.tsv", sep="\t", index_col=0)
d.head(500).rename(index={i: f"SYMBOL{n}" for n, i in enumerate(d.index[:500])}).to_csv(f"{tmp}/nojoin.tsv", sep="\t"); d.iloc[:int(0.30*len(d))].to_csv(f"{tmp}/thin.tsv", sep="\t")
(d / d.sum(axis=0) * 1e6).to_csv(f"{tmp}/cpm.tsv", sep="\t"); np.log2(d + 1).to_csv(f"{tmp}/logged.tsv", sep="\t"); d.iloc[:, :3].to_csv(f"{tmp}/three.tsv", sep="\t")
PYEOF
fails=0
for case in nojoin thin cpm logged; do
  if "${PY}" "${HERE}/score.py" --counts "${TMP}/${case}.tsv" --out "${TMP}/x.tsv" >/dev/null 2>"${TMP}/${case}.err"; then echo "   FAIL: ${case} was ACCEPTED"; fails=$((fails+1)); else echo "   refused ${case}: $(head -c 90 "${TMP}/${case}.err" | tr '\n' ' ')"; fi
done
if "${PY}" "${HERE}/score.py" --counts "${TMP}/three.tsv" --out "${TMP}/x.tsv" --profile "${TMP}/x.npz" >/dev/null 2>"${TMP}/three.err"; then echo "   FAIL: 3-sample batch under marginal transport was ACCEPTED"; fails=$((fails+1)); else echo "   refused 3-sample marginal profile: $(head -c 90 "${TMP}/three.err" | tr '\n' ' ')"; fi
"${PY}" "${HERE}/score.py" --counts "${TMP}/three.tsv" --out "${TMP}/x.tsv" --profile "${TMP}/x.npz" --transport raw >/dev/null 2>&1 || { echo "   FAIL: 3-sample batch under RAW transport was refused"; fails=$((fails+1)); }
[ "${fails}" -eq 0 ] || { echo "TEST FAILED: ${fails} guard(s)"; exit 1; }
echo "-- manifest"
"${PY}" - "${HERE}" <<'PYEOF'
import sys, json, hashlib, os
here = sys.argv[1]; m = json.load(open(f"{here}/MANIFEST.sha256.json")); bad = []
for f, v in m["files"].items():
    p = f"{here}/{f}"; h = hashlib.sha256(open(p, "rb").read()).hexdigest()
    if h != v["sha256"] or os.path.getsize(p) != v["bytes"]: bad.append(f)
if bad: print("FAIL: manifest mismatch " + ", ".join(bad)); sys.exit(1)
print(f"   {len(m['files'])} files verify")
PYEOF
[ $? -eq 0 ] || exit 1
echo; echo "== ALL CHECKS PASSED"
