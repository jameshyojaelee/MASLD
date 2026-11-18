#!/usr/bin/env bash
# End-to-end verification for masld-severity-v1.2.1 (v1.2 scorer bitwise + kleiner_sensitivity guards + manifest).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
V12="${BENCH}/release/masld-severity-v1.2"
V11="${BENCH}/release/masld-severity-v1.1"
OOD="${BENCH}/executions/full-axis-rebuild-and-ood-transfer-20260901T122932Z"
PHENO="${BENCH}/executions/gse268273-evaluator-phenotype-20260830T205701Z/gse268273_participant_phenotype.tsv"
PY="${PY:-${BENCH}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python}"
TMP="$(mktemp -d)"
trap 'rm -rf "${TMP}"' EXIT

echo "== masld-severity-v1.2.1 smoke test"
echo "   python:  ${PY}"
echo "   scratch: ${TMP}"

for f in "${HERE}/score.py" "${HERE}/kleiner_sensitivity.py" \
         "${HERE}/weights/masld_severity_v1_weights.npz" \
         "${HERE}/expected_ood.json" \
         "${HERE}/MANIFEST.sha256.json" \
         "${OOD}/gse268273_model_axis_expected_counts.npy" \
         "${OOD}/gse268273_model_axis_genes.tsv" \
         "${OOD}/gse268273_participants.tsv" "${PHENO}"; do
  [ -e "${f}" ] || { echo "FAIL: missing ${f}"; exit 1; }
done

echo "-- asserting score.py & weights are bitwise identical to v1.2"
cmp -s "${HERE}/score.py" "${V12}/score.py" || { echo "FAIL: score.py differs from v1.2"; exit 1; }
for w in "${V12}/weights/"*; do
  bw="$(basename "${w}")"
  cmp -s "${w}" "${HERE}/weights/${bw}" || { echo "FAIL: weights/${bw} differs from v1.2"; exit 1; }
done

echo "-- building GSE268273 counts TSV and running score.py"
"${PY}" - "${OOD}" "${TMP}" <<'PYEOF'
import sys, numpy as np, pandas as pd
ood, tmp = sys.argv[1], sys.argv[2]
C = np.load(f"{ood}/gse268273_model_axis_expected_counts.npy")
g = pd.read_csv(f"{ood}/gse268273_model_axis_genes.tsv", sep="\t")["stable_gene_id"].astype(str)
p = pd.read_csv(f"{ood}/gse268273_participants.tsv", sep="\t")["row_id"].astype(str)
df = pd.DataFrame(C.T, index=[f"{x}.7" for x in g], columns=p)
df.index.name = "gene_id"
df.to_csv(f"{tmp}/counts.tsv", sep="\t")
PYEOF

"${PY}" "${HERE}/score.py" \
  --counts "${TMP}/counts.tsv" \
  --out "${TMP}/scores.tsv" \
  --json "${TMP}/report.json" \
  --weights "${HERE}/weights/masld_severity_v1_weights.npz"

"${PY}" "${V11}/score.py" --counts "${TMP}/counts.tsv" --out "${TMP}/scores_v11.tsv" \
  --weights "${V11}/weights/masld_severity_v1_weights.npz" 2>/dev/null

"${PY}" - "${TMP}" "${PHENO}" "${HERE}/expected_ood.json" <<'PYEOF'
import json, sys, numpy as np, pandas as pd
from scipy import stats
tmp, pheno_path, expected_path = sys.argv[1], sys.argv[2], sys.argv[3]
exp = json.load(open(expected_path))
s = pd.read_csv(f"{tmp}/scores.tsv", sep="\t")
s11 = pd.read_csv(f"{tmp}/scores_v11.tsv", sep="\t")
ph = pd.read_csv(pheno_path, sep="\t")
ph["row_id"] = ph["row_id"].astype(str)
m = s.merge(ph[["row_id", "fibrosis_stage"]], left_on="sample_id", right_on="row_id")
rho = float(stats.spearmanr(m.latent_severity, m.fibrosis_stage).statistic)
assert len(m) == exp["n_participants"]
assert abs(rho - exp["ood_spearman"]) <= exp["tolerance"], f"{rho} != {exp['ood_spearman']}"
assert np.array_equal(s11.latent_severity.values, s.latent_severity.values)
print(f"   PASS: GSE268273 OOD Spearman {rho:.6f} (bitwise identical to v1.1/v1.2)")
PYEOF

echo "-- testing kleiner_sensitivity.py feasibility and predictor-vs-response guards"
"${PY}" "${HERE}/kleiner_sensitivity.py" --stage-var 1.44 --sigma-e 0.4655 \
  --assumption "hypothetical sensitivity check" --r 0.62 --stage-predictor-beta 0.41 > "${TMP}/kl_ok.json"
"${PY}" - "${TMP}/kl_ok.json" <<'PYEOF'
import json, sys
r = json.load(open(sys.argv[1]))
assert r["feasible"] is True and abs(r["conditional_latent_correlation"]) <= 1.0
print(f"   PASS: valid run r_latent={r['conditional_latent_correlation']:.4f}")
PYEOF

if "${PY}" "${HERE}/kleiner_sensitivity.py" --stage-var 0.5091 --sigma-e 0.6433 \
     --assumption "gse193066 infeasible check" --r 0.65 >/dev/null 2>"${TMP}/kl_inf.err"; then
  echo "FAIL: infeasible r_latent > 1.0 was ACCEPTED"; exit 1
else
  echo "   PASS: refused |r_latent| > 1.0: $(cat "${TMP}/kl_inf.err" | tr '\n' ' ' | head -c 100)"
fi

if "${PY}" "${HERE}/kleiner_sensitivity.py" --stage-var 1.44 --sigma-e 0.4655 \
     --assumption "invalid response slope" --stage-role response --stage-predictor-beta 0.41 >/dev/null 2>"${TMP}/kl_resp.err"; then
  echo "FAIL: --stage-role response was ACCEPTED"; exit 1
else
  echo "   PASS: refused --stage-role response: $(cat "${TMP}/kl_resp.err" | tr '\n' ' ' | head -c 100)"
fi

echo "-- verifying MANIFEST.sha256.json"
"${PY}" - "${HERE}" <<'PYEOF'
import hashlib, json, os, sys
here = sys.argv[1]
m = json.load(open(f"{here}/MANIFEST.sha256.json"))
for rel, v in m["files"].items():
    p = os.path.join(here, rel)
    assert hashlib.sha256(open(p, "rb").read()).hexdigest() == v["sha256"], f"sha256 mismatch: {rel}"
    assert os.path.getsize(p) == v["bytes"], f"size mismatch: {rel}"
print(f"   PASS: {len(m['files'])} manifest entries verify")
PYEOF

echo "== ALL CHECKS PASSED"
