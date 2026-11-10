"""U14: HyenaDNA reporter-assay gain over the sequence-enriched ridge rungs (FIG7_HANDOFF item 1).

Same outcome, score and paired long-range-block bootstrap as mpra_7b.py; references are the ladder rungs
B0 (allele identity), B1 (CpG), B2 (composition), B3 (flank and dinucleotide; the PRESPEC decision rung).

usage: mpra_ladder.py OUTDIR B SEED [B SEED ...]
"""
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
HYENA = EX / "model-cpu-train-605-21099008/fit/oof_predictions.tsv.gz"
LADDER = EX / "mpra-enriched-baseline-20260905T192000Z-21454084/ladder/oof_predictions.tsv.gz"
CONTEXTS = ["HepG2_control", "HepG2_PAOA"]
RUNGS = ["allele_identity_ridge_b0", "enriched_b1_cpg", "enriched_b2_composition", "enriched_b3_flank_dinucleotide"]

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
runs = [(int(sys.argv[i]), int(sys.argv[i + 1])) for i in range(2, len(sys.argv), 2)]

h = pd.read_csv(HYENA, sep="\t", usecols=["seed", "assay_context_id", "element_id", "long_range_block_id",
                                          "model_id", "head_id", "prediction"])
h = h[(h.model_id == "hyenadna") & (h.head_id == "delta_ridge")]
elements = sorted(h.element_id.unique())
assert len(elements) == 1033
blk = h.drop_duplicates("element_id").set_index("element_id")["long_range_block_id"].loc[elements].to_numpy()
lad = pd.read_csv(LADDER, sep="\t")
assert set(lad.element_id) == set(elements)

cands = ["hyenadna"] + RUNGS
P = np.empty((len(cands), 2, len(elements)))
for j, c in enumerate(CONTEXTS):
    P[0, j] = h[h.assay_context_id == c].groupby("element_id")["prediction"].mean().loc[elements].to_numpy()
    for i, r in enumerate(RUNGS):
        x = lad[(lad.rung == r) & (lad.assay_context_id == c)]
        assert x.groupby("element_id")["seed"].nunique().eq(5).all()
        P[i + 1, j] = x.groupby("element_id")["prediction"].mean().loc[elements].to_numpy()
P = np.round(P, 12)   # restore ties broken by 1e-17 averaging noise (see mpra_7b.py)

o = pd.read_csv(OUTCOMES, sep="\t")
o = o[o.context_id.isin(CONTEXTS) & o.element_id.isin(elements) & (o.assay_state == "observed")]
w = o.pivot_table(index=["element_id", "context_id", "experimental_replicate"], columns="allele",
                  values=["DNA", "RNA"], aggfunc="first")
lr = (np.log2((w[("RNA", "alt")] + 0.5) / (w[("DNA", "alt")] + 0.5))
      - np.log2((w[("RNA", "ref")] + 0.5) / (w[("DNA", "ref")] + 0.5)))
y = lr.groupby(level=[0, 1]).mean()
Y = np.vstack([y.xs(c, level="context_id").loc[elements].to_numpy() for c in CONTEXTS])


def scores(idx):
    zs = np.zeros(len(cands))
    for j in range(2):
        ry = rankdata(Y[j, idx]); ry = ry - ry.mean()
        rp = rankdata(P[:, j, idx], axis=1); rp = rp - rp.mean(axis=1, keepdims=True)
        zs += np.arctanh((rp @ ry) / np.sqrt((rp ** 2).sum(1) * (ry ** 2).sum()))
    return np.tanh(zs / 2)


point = scores(np.arange(len(elements)))
print("point", dict(zip(cands, point)))
blocks = sorted(set(blk))
block_list = [np.flatnonzero(blk == b) for b in blocks]
rows = []
for B, seed in runs:
    rng = np.random.default_rng(seed)
    M = draw_block_multiplicities(rng, len(blocks), B)[:, np.searchsorted(blocks, blk)]
    worst = check_against_rankdata(Y, P, M[:20], None)
    assert worst < 1e-10, worst
    G = np.vstack([fisher_scores(Y, P, M[s:s + 1000]) for s in range(0, B, 1000)])
    for i, r in enumerate(RUNGS):
        g = G[:, 0] - G[:, i + 1]
        rows.append({"reference": r, "B": B, "seed": seed, "ref_score": point[i + 1], "hyenadna_score": point[0],
                     "gain": point[0] - point[i + 1], "lo95": np.percentile(g, 2.5), "hi95": np.percentile(g, 97.5),
                     "bonf7_lo": np.percentile(g, 100 * 0.025 / 7), "bonf7_hi": np.percentile(g, 100 * (1 - 0.025 / 7)),
                     "frac_le0": float(np.mean(g <= 0))})
    print("done", B, seed, flush=True)
t = pd.DataFrame(rows)
t.to_csv(out / "mpra_ladder_rederived.tsv", sep="\t", index=False, float_format="%.6f")
print(t.to_string())
