"""Verifier (F7-02): HyenaDNA gain over allele ridge (B0) and sequence-feature ridges (B1-B3), recomputed from raw
reporter counts and OOF predictions. Own block bootstrap (seed 777, 4,000 draws), percentile 95%."""
import gzip, numpy as np, pandas as pd
from scipy.stats import rankdata
X = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/"
CTX = ["HepG2_control", "HepG2_PAOA"]
o = pd.read_csv(X + "gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz", sep="\t")
o = o[o.context_id.isin(CTX) & (o.assay_state == "observed")]
w = o.pivot_table(index=["element_id", "context_id", "experimental_replicate"], columns="allele", values=["DNA", "RNA"])
pc = 0.5
d = np.log2((w[("RNA", "alt")] + pc) / (w[("DNA", "alt")] + pc)) - np.log2((w[("RNA", "ref")] + pc) / (w[("DNA", "ref")] + pc))
y = d.groupby(level=[0, 1]).mean()
h = pd.read_csv(X + "model-cpu-train-605-21099008/fit/oof_predictions.tsv.gz", sep="\t")
h["m"] = h.model_id + "|" + h.head_id
P = {}
for key, name in [("hyenadna|delta_ridge", "HyenaDNA"), ("available_simple_controls|allele_identity_ridge", "B0_605")]:
    P[name] = h[h.m == key].groupby(["element_id", "assay_context_id"]).prediction.mean()
lad = pd.read_csv(X + "mpra-enriched-baseline-20260905T192000Z-21454084/ladder/oof_predictions.tsv.gz", sep="\t")
for r, g in lad.groupby("rung"):
    P[r] = g.groupby(["element_id", "assay_context_id"]).prediction.mean()
blocks = h.drop_duplicates("element_id").set_index("element_id").long_range_block_id
els = sorted(set(P["HyenaDNA"].index.get_level_values(0)))
print("elements", len(els), "blocks", blocks.loc[els].nunique())
Y = {c: y.xs(c, level=1).loc[els].to_numpy() for c in CTX}
M = {k: {c: v.xs(c, level=1).loc[els].to_numpy() for c in CTX} for k, v in P.items()}
def sp(a, b):
    ra, rb = rankdata(a) - (len(a) + 1) / 2, rankdata(b) - (len(b) + 1) / 2
    return (ra * rb).sum() / np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
def macro(k, idx=None):
    v = [sp(Y[c] if idx is None else Y[c][idx], M[k][c] if idx is None else M[k][c][idx]) for c in CTX]
    return np.tanh(np.mean(np.arctanh(np.clip(v, -0.999999, 0.999999))))
pt = {k: macro(k) for k in M}
print({k: round(v, 5) for k, v in pt.items()})
blk = blocks.loc[els].to_numpy(); ub = np.unique(blk); members = {b: np.where(blk == b)[0] for b in ub}
rng = np.random.default_rng(777); B = 4000
refs = ["allele_identity_ridge_b0", "enriched_b1_cpg", "enriched_b2_composition", "enriched_b3_flank_dinucleotide"]
draws = {r: [] for r in refs}
for _ in range(B):
    idx = np.concatenate([members[b] for b in rng.choice(ub, len(ub), replace=True)])
    hy = macro("HyenaDNA", idx)
    for r in refs: draws[r].append(hy - macro(r, idx))
for r in refs:
    a = np.array(draws[r]); print(r, "gain", round(pt["HyenaDNA"] - pt[r], 5), "95%", np.round(np.percentile(a, [2.5, 97.5]), 5),
          "Bonf7", np.round(np.percentile(a, [100 * 0.025 / 7, 100 * (1 - 0.025 / 7)]), 5), "P(gain<=0)", round((a <= 0).mean(), 4))
