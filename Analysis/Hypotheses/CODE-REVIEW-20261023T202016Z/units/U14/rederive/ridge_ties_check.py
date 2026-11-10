"""Check whether the allele-identity ridge score depends on floating-point tie breaking."""
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, rankdata
EX = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
d = pd.read_csv(EX / "model-cpu-train-605-21099008/fit/oof_predictions.tsv.gz", sep="\t",
                usecols=["seed", "assay_context_id", "element_id", "model_id", "head_id", "prediction"])
d = d[(d.model_id == "available_simple_controls") & (d.head_id == "allele_identity_ridge")]
elements = sorted(d.element_id.unique())
o = pd.read_csv(EX / "gse281364-replicate-outcomes-21066307/outcomes/replicate_outcomes.tsv.gz", sep="\t")
o = o[o.context_id.isin(["HepG2_control", "HepG2_PAOA"]) & o.element_id.isin(elements) & (o.assay_state == "observed")]
w = o.pivot_table(index=["element_id", "context_id", "experimental_replicate"], columns="allele", values=["DNA", "RNA"], aggfunc="first")
lr = np.log2((w[("RNA","alt")]+0.5)/(w[("DNA","alt")]+0.5)) - np.log2((w[("RNA","ref")]+0.5)/(w[("DNA","ref")]+0.5))
y = lr.groupby(level=[0,1]).mean()
zs = {k: [] for k in ["groupby_mean", "vstack_mean", "rounded12"]}
for c in ["HepG2_control", "HepG2_PAOA"]:
    yc = y.xs(c, level="context_id").loc[elements].to_numpy()
    dc = d[d.assay_context_id == c]
    a = dc.groupby("element_id")["prediction"].mean().loc[elements].to_numpy()
    seeds = sorted(dc.seed.unique())
    mat = np.vstack([dc[dc.seed == s].set_index("element_id")["prediction"].loc[elements].to_numpy() for s in seeds])
    b = np.mean(mat, axis=0)
    print(c, "distinct per-seed values", [len(np.unique(m)) for m in mat], "distinct mean (groupby)", len(np.unique(a)),
          "distinct mean (vstack)", len(np.unique(b)), "distinct rounded 1e-12", len(np.unique(np.round(b, 12))),
          "max |a-b|", np.abs(a - b).max())
    zs["groupby_mean"].append(np.arctanh(spearmanr(yc, a)[0]))
    zs["vstack_mean"].append(np.arctanh(spearmanr(yc, b)[0]))
    zs["rounded12"].append(np.arctanh(spearmanr(yc, np.round(b, 12))[0]))
for k, v in zs.items():
    print(k, "ridge score", np.tanh(np.mean(v)))
