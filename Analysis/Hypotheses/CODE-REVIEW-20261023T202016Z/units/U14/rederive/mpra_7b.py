"""U14 independent re-derivation of the Figure 7B reporter-assay (GSE281364) gains and intervals.

Written from METHODS.md "Reporter-assay rows" and the task inputs, not from the review producer
(02_mpra_bonferroni.py) or the benchmark evaluator.

Outcome per element and HepG2 context: mean over the 4 experimental replicates of
log2((alt_RNA+0.5)/(alt_DNA+0.5)) - log2((ref_RNA+0.5)/(ref_DNA+0.5)) (task definition, pseudocount 0.5).
Prediction per element and context: mean of the 5 fixed-seed out-of-fold predictions.
Score: tanh(mean over the two contexts of arctanh(Spearman)). Gain = row score - allele-identity ridge score.
Bootstrap: resample the 239 long-range blocks with replacement (both contexts and all rows on the same draw),
recompute every score on the resampled elements, percentile intervals (95% and Bonferroni-7 two-sided).

usage: mpra_7b.py OUTDIR B SEED [B SEED ...]
"""
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mpra_fastboot import check_against_rankdata, draw_block_multiplicities, fisher_scores  # noqa: E402

EX = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
OUTCOMES = EX / "gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz"
FITS = [EX / "model-cpu-train-605-21099008/fit/oof_predictions.tsv.gz",
        EX / "model-cpu-train-604-21097189/fit/oof_predictions.tsv.gz"]
CONTEXTS = ["HepG2_control", "HepG2_PAOA"]
REF = ("available_simple_controls", "allele_identity_ridge")
ROWS = [("hyenadna", "delta_ridge"), ("sei", "native_40class_ridge"),
        ("nucleotide_transformer", "delta_ridge"), ("caduceus", "delta_ridge"),
        ("dnabert2", "delta_ridge"), ("enformer_crested_restricted_port", "rc_ensemble_sar_ridge"),
        ("enformer_crested_restricted_port", "rc_ensemble_sad_ridge")]
EXTRA = [("hyenadna", "full_ridge"), ("nucleotide_transformer", "full_ridge"),
         ("caduceus", "full_ridge"), ("dnabert2", "full_ridge")]   # the 4 heads that make the family of 11


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
runs = [(int(sys.argv[i]), int(sys.argv[i + 1])) for i in range(2, len(sys.argv), 2)]

# predictions: element x context x candidate, mean over 5 seeds
pred_frames = []
for f in FITS:
    d = pd.read_csv(f, sep="\t", usecols=["seed", "assay_context_id", "element_id", "long_range_block_id",
                                           "model_id", "head_id", "prediction"])
    pred_frames.append(d)
pred = pd.concat(pred_frames)
n_before = len(pred)
pred = pred.drop_duplicates(["seed", "assay_context_id", "element_id", "model_id", "head_id"])
# the allele ridge is present in both fit files; check the two copies agree, keep one
print("prediction rows", n_before, "after dedup on keys", len(pred))
cands = [REF] + ROWS + EXTRA
seed_counts = pred.groupby(["model_id", "head_id", "element_id", "assay_context_id"])["seed"].nunique()
assert (seed_counts == 5).all(), "every candidate x element x context must have 5 seeds"
mean_pred = pred.groupby(["model_id", "head_id", "element_id", "assay_context_id"])["prediction"].mean()
elements = sorted(pred["element_id"].unique())
assert len(elements) == 1033, len(elements)
blk = pred.drop_duplicates("element_id").set_index("element_id")["long_range_block_id"].loc[elements]
assert pred.groupby("element_id")["long_range_block_id"].nunique().max() == 1
blocks = sorted(blk.unique())
assert len(blocks) == 239, len(blocks)

# check that duplicated ridge rows in the two fit files are identical
r = pred_frames[0].query("model_id == @REF[0] and head_id == @REF[1]").set_index(
    ["seed", "assay_context_id", "element_id"])["prediction"].sort_index()
r2 = pred_frames[1].query("model_id == @REF[0] and head_id == @REF[1]").set_index(
    ["seed", "assay_context_id", "element_id"])["prediction"].sort_index()
ridge_maxdiff = float(np.abs(r - r2.loc[r.index]).max())
print("allele ridge 605 vs 604 max abs diff", ridge_maxdiff)

P = np.empty((len(cands), 2, len(elements)))
for i, (m, h) in enumerate(cands):
    for j, c in enumerate(CONTEXTS):
        P[i, j] = mean_pred.loc[(m, h)].xs(c, level="assay_context_id").loc[elements].to_numpy()
# The allele-identity ridge takes about 60 distinct values per context; averaging 5 seeds leaves
# 1e-17 differences between values that are equal, which break ties at random (ridge_ties_check.py).
# Rounding to 12 decimals restores the ties; it does not create ties in the embedding heads.
n_distinct_raw = [[len(np.unique(P[i, j])) for j in range(2)] for i in range(len(cands))]
P = np.round(P, 12)
n_distinct_round = [[len(np.unique(P[i, j])) for j in range(2)] for i in range(len(cands))]
print("distinct prediction values per candidate x context, raw vs rounded:",
      [(f"{m}::{h}", a, b) for (m, h), a, b in zip(cands, n_distinct_raw, n_distinct_round)])

# outcomes
o = pd.read_csv(OUTCOMES, sep="\t")
o = o[o["context_id"].isin(CONTEXTS) & o["element_id"].isin(elements)]
print("outcome rows for 1033 elements x 2 contexts:", len(o), "assay_state counts:", o["assay_state"].value_counts().to_dict())
o = o[o["assay_state"] == "observed"]
w = o.pivot_table(index=["element_id", "context_id", "experimental_replicate"], columns="allele",
                  values=["DNA", "RNA"], aggfunc="first")
assert w.notna().all().all(), "every element x context x replicate needs both alleles"
lr = (np.log2((w[("RNA", "alt")] + 0.5) / (w[("DNA", "alt")] + 0.5))
      - np.log2((w[("RNA", "ref")] + 0.5) / (w[("DNA", "ref")] + 0.5)))
nrep = lr.groupby(level=[0, 1]).size()
assert (nrep == 4).all() and len(nrep) == 2066, (nrep.value_counts(), len(nrep))
y_mean = lr.groupby(level=[0, 1]).mean()
Y = np.vstack([y_mean.xs(c, level="context_id").loc[elements].to_numpy() for c in CONTEXTS])


def scores(idx):
    """Fisher-mean Spearman for every candidate on element index array idx."""
    zs = np.zeros(len(cands))
    for j in range(2):
        ry = rankdata(Y[j, idx])
        rp = rankdata(P[:, j, idx], axis=1)
        ry = ry - ry.mean()
        rp = rp - rp.mean(axis=1, keepdims=True)
        rho = (rp @ ry) / np.sqrt((rp ** 2).sum(1) * (ry ** 2).sum())
        zs += np.arctanh(rho)
    return np.tanh(zs / 2)


all_idx = np.arange(len(elements))
point = scores(all_idx)
point_rows = {f"{m}::{h}": float(point[i]) for i, (m, h) in enumerate(cands)}
print("point scores", json.dumps(point_rows, indent=1))

elem_by_block = {b: np.flatnonzero(blk.to_numpy() == b) for b in blocks}
block_list = [elem_by_block[b] for b in blocks]

rows = []
meta = {"inputs": {str(p): sha256(p) for p in FITS + [OUTCOMES]}, "runs": [],
        "unit": "long_range_block_id (239 blocks; elements of a drawn block enter together, both contexts)",
        "ridge_605_vs_604_maxdiff": ridge_maxdiff, "generator": "numpy.random.default_rng(seed).integers"}
elem_block = np.searchsorted(blocks, blk.to_numpy())
for B, seed in runs:
    rng = np.random.default_rng(seed)
    counts = draw_block_multiplicities(rng, len(blocks), B)
    M = counts[:, elem_block]
    worst = check_against_rankdata(Y, P, M[:20], None)
    assert worst < 1e-10, worst
    G = np.vstack([fisher_scores(Y, P, M[s:s + 1000]) for s in range(0, B, 1000)])
    meta.setdefault("fast_vs_rankdata_maxdiff_first20", {})[f"{B}_{seed}"] = worst
    np.save(out / f"mpra_boot_scores_B{B}_seed{seed}.npy", G)
    gains = G[:, 1:] - G[:, [0]]
    for i, (m, h) in enumerate(cands[1:]):
        g = gains[:, i]
        fam = "F7" if (m, h) in ROWS else "F11_extra"
        rows.append({"B": B, "seed": seed, "model": m, "head": h, "family": fam,
                     "gain": point[i + 1] - point[0],
                     "lo95": np.percentile(g, 2.5), "hi95": np.percentile(g, 97.5),
                     "bonf7_lo": np.percentile(g, 100 * 0.025 / 7), "bonf7_hi": np.percentile(g, 100 * (1 - 0.025 / 7)),
                     "bonf11_lo": np.percentile(g, 100 * 0.025 / 11), "bonf11_hi": np.percentile(g, 100 * (1 - 0.025 / 11)),
                     "frac_le0": float(np.mean(g <= 0))})
    meta["runs"].append({"B": B, "seed": seed})
    print("done B", B, "seed", seed, flush=True)

t = pd.DataFrame(rows)
t.to_csv(out / "mpra_7b_rederived.tsv", sep="\t", index=False, float_format="%.6f")
meta["point_scores"] = point_rows
json.dump(meta, open(out / "mpra_7b_meta.json", "w"), indent=1)
print(t.to_string())
