"""Recompute p, BH q and q<0.05 counts from the stored AN4 family null draws; MC spread of the count."""
import glob, sys
import numpy as np, pandas as pd
W = sys.argv[1]
fam = pd.read_csv(f"{W}/combined/family_all_programs.tsv", sep="\t")
fam = fam[fam.testable.astype(bool)]
def bh(p):
    n=len(p); o=np.argsort(p, kind="mergesort"); a=p[o]*n/np.arange(1,n+1); a=np.minimum.accumulate(a[::-1])[::-1]; q=np.empty(n); q[o]=np.minimum(a,1); return q
rng = np.random.default_rng(20261023)
for ds in ["GSE192741", "Vu_et_al_2025"]:
    nd = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"{W}/famr2_*_s*/family_null_draws_{ds}.tsv.gz")])
    g = fam[fam.dataset == ds].set_index("program_id")
    k = {pid: int((d.residual_moran_i.to_numpy() >= g.loc[pid, "residual_moran_i"]).sum()) for pid, d in nd.groupby("program_id")}
    B = nd.groupby("program_id").size()
    p = np.array([(1 + k[pid]) / (B[pid] + 1) for pid in g.index])
    q = bh(p)
    print(ds, "programs with stored null", len(k), "B unique", B.unique(), "p equal recorded", np.isclose(p, g.residual_pvalue.to_numpy()).sum(),
          "q<0.05", (q < 0.05).sum(), "recorded", (g.residual_padj_family < 0.05).sum())
    # parametric MC spread of the count across hypothetical seeds: k ~ Binomial(B, p_hat)
    ph = np.clip(np.array([k[pid] for pid in g.index]) / 9999, 0, 1)
    counts = []
    for _ in range(2000):
        kk = rng.binomial(9999, ph)
        counts.append((bh((kk + 1) / 10000) < 0.05).sum())
    counts = np.array(counts)
    print("  MC count over 2000 pseudo-seeds: mean %.2f sd %.2f min %d max %d; P(count==recorded) %.3f; 2.5-97.5%% %s"
          % (counts.mean(), counts.std(), counts.min(), counts.max(), np.mean(counts == (q < 0.05).sum()), np.percentile(counts, [2.5, 97.5])))
