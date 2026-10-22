"""Independent recomputation of the GSE281364 MPRA rows (Figure 7B/7C).

Inputs (analysis outputs, not figure tables):
  executions/gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz
  executions/model-check-220-21088696/contract/row_universe.tsv   (element -> long-range block)
  executions/model-cpu-train-605-21099008/fit/oof_predictions.tsv.gz (DNA LMs + controls)
  executions/model-cpu-train-604-21097189/fit/oof_predictions.tsv.gz (Sei, Enformer + controls)

Metric (as recorded in the run receipts): per-context Spearman of the five-seed mean
prediction against the mean of four replicate log2(alt RNA/DNA) - log2(ref RNA/DNA)
(pseudocount 0.5), Fisher-z macro over HepG2_control and HepG2_PAOA. Gain is
candidate minus allele_identity_ridge. Interval: paired bootstrap over the 239
long-range blocks, 10,000 draws, seed 20260825 (same draw for every candidate).
Spearman is computed here with pandas ranks, not with the producer's function.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
OUT = Path(sys.argv[1])
CONTEXTS = ["HepG2_control", "HepG2_PAOA"]
PSEUDO = 0.5
B = 10000
SEED = 20260825

rows = pd.read_csv(ROOT / "model-check-220-21088696/contract/row_universe.tsv", sep="\t", dtype=str)
elements = sorted(rows["element_id"].unique())
blk = rows.drop_duplicates("element_id").set_index("element_id")["long_range_block_id"]
assert len(elements) == 1033 and blk.nunique() == 239

oc = pd.read_csv(ROOT / "gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz", sep="\t",
                 dtype={"element_id": str})
oc = oc[oc["element_id"].isin(elements) & oc["context_id"].isin(CONTEXTS) & (oc["assay_state"] == "observed")]
oc["act"] = np.log2((oc["RNA"] + PSEUDO) / (oc["DNA"] + PSEUDO))
wide = oc.pivot_table(index=["element_id", "context_id", "experimental_replicate"], columns="allele", values="act")
wide["delta"] = wide["alt"] - wide["ref"]
nrep = wide.groupby(level=[0, 1])["delta"].count()
assert (nrep == 4).all() and len(nrep) == 2066, (nrep.value_counts(), len(nrep))
obs = wide.groupby(level=[0, 1])["delta"].mean().unstack("context_id").loc[elements]

preds = {}
for run in ["model-cpu-train-605-21099008", "model-cpu-train-604-21097189"]:
    p = pd.read_csv(ROOT / run / "fit/oof_predictions.tsv.gz", sep="\t",
                    usecols=["seed", "assay_context_id", "element_id", "model_id", "head_id", "prediction"],
                    dtype={"element_id": str})
    nseed = p.groupby(["model_id", "head_id"])["seed"].nunique()
    assert (nseed == 5).all()
    ens = p.groupby(["model_id", "head_id", "assay_context_id", "element_id"])["prediction"].mean()
    for (m, h), sub in ens.groupby(level=[0, 1]):
        key = f"{m}::{h}"
        mat = sub.droplevel([0, 1]).unstack("assay_context_id").loc[elements, CONTEXTS]
        if key in preds:
            assert np.allclose(preds[key].values, mat.values), key  # controls shared by both runs
            continue
        preds[key] = mat

def rank_spearman(a, b):
    ra = pd.Series(a).rank().to_numpy()
    rb = pd.Series(b).rank().to_numpy()
    if np.all(rb == rb[0]):
        return np.nan
    ra = ra - ra.mean(); rb = rb - rb.mean()
    return float((ra * rb).sum() / np.sqrt((ra ** 2).sum() * (rb ** 2).sum()))

def macro(vals):
    if any(not np.isfinite(v) for v in vals):
        return np.nan
    v = np.clip(np.asarray(vals), -0.999999, 0.999999)
    return float(np.tanh(np.mean(np.arctanh(v))))

point = {k: macro([rank_spearman(obs[c].to_numpy(), v[c].to_numpy()) for c in CONTEXTS]) for k, v in preds.items()}
ref = "available_simple_controls::allele_identity_ridge"

blocks = blk.loc[elements].to_numpy()
uniq = sorted(set(blocks.tolist()))
by_block = {b: np.flatnonzero(blocks == b) for b in uniq}
rng = np.random.default_rng(SEED)
cands = [k for k in preds if np.isfinite(point[k])]
O = {c: obs[c].to_numpy() for c in CONTEXTS}
P = {k: {c: preds[k][c].to_numpy() for c in CONTEXTS} for k in cands}
draws = {k: np.full(B, np.nan) for k in cands}
from scipy.stats import rankdata
def fsp(a, b):
    ra = rankdata(a); rb = rankdata(b)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return np.nan if d == 0 else float((ra * rb).sum() / d)
for it in range(B):
    s = rng.integers(0, len(uniq), size=len(uniq))
    idx = np.concatenate([by_block[uniq[i]] for i in s])
    for k in cands:
        draws[k][it] = macro([fsp(O[c][idx], P[k][c][idx]) for c in CONTEXTS])

np.savez_compressed(OUT / "mpra_recomputed_draws.npz", **{k: v for k, v in draws.items()})
displayed = {
    "hyenadna": None, "sei": "sei::native_40class_ridge", "nucleotide_transformer": None, "caduceus": None,
    "dnabert2": None, "Enformer SAR": "enformer_crested_restricted_port::rc_ensemble_sar_ridge",
    "Enformer SAD": "enformer_crested_restricted_port::rc_ensemble_sad_ridge",
}
recs = []
alpha = 0.05
m = 7
for k in cands:
    g = draws[k] - draws[ref]
    ok = np.isfinite(g)
    recs.append({
        "candidate": k, "primary_score": point[k], "gain_vs_allele_ridge": point[k] - point[ref],
        "valid_draws": int(ok.sum()),
        "gain_ci95_low": np.quantile(g[ok], 0.025), "gain_ci95_high": np.quantile(g[ok], 0.975),
        "bonf7_two_sided_low": np.quantile(g[ok], alpha / m / 2), "bonf7_two_sided_high": np.quantile(g[ok], 1 - alpha / m / 2),
        "bonf7_one_sided_low": np.quantile(g[ok], alpha / m),
        "frac_draws_gain_le_0": float((g[ok] <= 0).mean()),
        "boot_se": float(np.std(g[ok], ddof=1)),
    })
res = pd.DataFrame(recs).sort_values("gain_vs_allele_ridge", ascending=False)
res.to_csv(OUT / "mpra_recomputed.tsv", sep="\t", index=False, float_format="%.6f")
print(res.to_string())
print("allele ridge point", point[ref])
json.dump({"allele_identity_ridge_point": point[ref], "n_elements": 1033, "n_blocks": 239, "B": B, "seed": SEED},
          open(OUT / "mpra_recomputed_meta.json", "w"), indent=1)
