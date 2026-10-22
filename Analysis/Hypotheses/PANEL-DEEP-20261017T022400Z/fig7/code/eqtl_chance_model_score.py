"""Within-block label-permutation chance for the 7D signed-eQTL strata, using each stratum's best head's own
held-out score (keeps the score's between-block structure), 1,000 permutations, seed 9874."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata
EX = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
OUT = Path(sys.argv[1])
A = pd.read_parquet(EX / "dnalm_heads-20260831T183647Z/tables/analysis_rows.parquet", columns=["abs_z", "beta_alt", "ld_fold", "block"])
SO = {M: np.load(EX / f"eqtl_sign_reliability-20260831T195702Z/tables/oof_{M}.npz") for M in ["hyenadna", "nucleotide_transformer"]}
st = pd.read_csv(OUT / "eqtl_strata_recomputed.tsv", sep="\t").set_index("stratum")
rng = np.random.default_rng(9874)
res = {}
for s in st.index:
    M, fs = st.loc[s, "best_head"].split("::")
    pos = SO[M][f"__index__{s}"]; tl = np.nonzero(A.ld_fold.to_numpy()[pos] == 4)[0]; g = pos[tl]
    sc = SO[M][f"{s}__{fs}__beta_alt"][tl]; lb = (A.beta_alt.to_numpy()[g] > 0).astype(int); bk = A.block.to_numpy()[g]
    r = rankdata(sc); n1 = lb.sum(); n0 = len(lb) - n1
    order = np.argsort(bk, kind="stable"); segs = np.split(order, np.nonzero(np.diff(bk[order]))[0] + 1)
    vals = []
    for _ in range(1000):
        yp = lb.copy()
        for sg in segs:
            if len(sg) > 1:
                yp[sg] = lb[rng.permutation(sg)]
        vals.append((r[yp == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
    vals = np.array(vals)
    res[s] = {"head": st.loc[s, "best_head"], "chance_median": float(np.median(vals)), "lo": float(np.quantile(vals, .025)), "hi": float(np.quantile(vals, .975))}
    print(s, res[s], flush=True)
json.dump(res, open(OUT / "eqtl_chance_model_score.json", "w"), indent=1)
