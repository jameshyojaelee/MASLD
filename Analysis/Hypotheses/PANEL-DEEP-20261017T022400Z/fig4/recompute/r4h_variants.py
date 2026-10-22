"""Fig 4H: grid of reasonable estimand variants on the 367 disease units, to locate which one the panel
plots and how rho / same-direction counts move across them. System-level only (light)."""
import numpy as np, pandas as pd
from scipy import stats
root = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out = f"{root}/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
HAC = f"{root}/RNA-seq/results/histology_anchored_continuum"
SYS = f"{HAC}/molecular_systems_inference/hac-system-inference-corrected774-20261009T171017Z/tables"
CORE = f"{HAC}/candidates/hac-continuum-corrected774-20261009T171017Z"
EVAL = ["GSE162694", "GSE213621"]


def ols(X, Y):
    XtX = np.linalg.pinv(X.T @ X); B = XtX @ X.T @ Y; R = Y - X @ B
    df = X.shape[0] - np.linalg.matrix_rank(X); s2 = (R ** 2).sum(0) / df
    return B, np.sqrt(np.outer(np.diag(XtX), s2))


zs = lambda x: (x - x.mean()) / x.std(ddof=1)
fp = pd.read_csv(f"{CORE}/projection/participant_scores.tsv", sep="\t")
fp = fp[fp.axis_id == "fixed_projection"][["sample_id", "dataset", "axis_raw", "axis_percentile", "inferred_sex", "fibrosis_stage"]]
sp = pd.read_csv(f"{CORE}/unsupervised/participant_scores.tsv", sep="\t")
sp = sp[sp.axis_id == "signature_pc1"][["sample_id", "dataset", "axis_raw"]].rename(columns={"axis_raw": "axis_sigpc1"})
cov = fp.rename(columns={"axis_raw": "axis_fp"}).merge(sp, on=["sample_id", "dataset"])
cov = cov[cov.dataset.isin(EVAL)]
cov["disease"] = cov.fibrosis_stage.notna() & cov.inferred_sex.notna() & ~((cov.dataset == "GSE213621") & (cov.fibrosis_stage == 0))
ss = pd.read_csv(f"{SYS}/participant_system_scores.tsv.gz", sep="\t")
ss = ss[ss.dataset.isin(EVAL)]
rows = []
for agg in ["collection_balanced", "feature_balanced_median"]:
    W = ss[ss.aggregation_method == agg].pivot_table(index="sample_id", columns="community_id", values="system_score_raw")
    for ymode in ["zall", "zdis"]:
        for axis in ["axis_fp", "axis_sigpc1", "axis_percentile"]:
            for stage_model in ["stage_z_plus_sex", "stage_z_only"]:
                per = {}
                for ds in EVAL:
                    call = cov[cov.dataset == ds].set_index("sample_id")
                    Yall = W.loc[call.index].values
                    Yz = (Yall - Yall.mean(0)) / Yall.std(0, ddof=1)
                    c = call[call.disease]
                    idx = call.index.get_indexer(c.index)
                    Y = Yz[idx] if ymode == "zall" else zs(Yall[idx])
                    sex = (c.inferred_sex == "M").astype(float).values
                    st = zs(c.fibrosis_stage.astype(float).values)
                    X1 = np.column_stack([np.ones(len(c)), st, sex]) if stage_model == "stage_z_plus_sex" else np.column_stack([np.ones(len(c)), st])
                    B1, S1 = ols(X1, Y)
                    lev = sorted(c.fibrosis_stage.unique())
                    D = np.column_stack([(c.fibrosis_stage.values == L).astype(float) for L in lev[1:]])
                    X2 = np.column_stack([np.ones(len(c)), zs(c[axis].values.astype(float)), D, sex])
                    B2, S2 = ols(X2, Y)
                    per[ds] = (B1[1], S1[1], B2[1], S2[1])
                for meta in ["fixed", "dl_random"]:
                    def m(b, s):
                        w = 1 / s ** 2; bf = (w * b).sum(0) / w.sum(0)
                        if meta == "fixed":
                            return bf
                        Q = (w * (b - bf) ** 2).sum(0); k = b.shape[0]
                        tau2 = np.maximum(0, (Q - (k - 1)) / (w.sum(0) - (w ** 2).sum(0) / w.sum(0)))
                        ws = 1 / (s ** 2 + tau2); return (ws * b).sum(0) / ws.sum(0)
                    sb = m(np.vstack([per[d][0] for d in EVAL]), np.vstack([per[d][1] for d in EVAL]))
                    cb = m(np.vstack([per[d][2] for d in EVAL]), np.vstack([per[d][3] for d in EVAL]))
                    rho = stats.spearmanr(sb, cb).correlation
                    rows.append(dict(agg=agg, ymode=ymode, axis=axis, stage_model=stage_model, meta=meta, rho=round(rho, 4),
                                     same=int((np.sign(sb) == np.sign(cb)).sum()), both_pos=int(((sb > 0) & (cb > 0)).sum()),
                                     both_neg=int(((sb < 0) & (cb < 0)).sum()), n_near0_x=int((np.abs(sb) < 0.01).sum()),
                                     n_near0_y=int((np.abs(cb) < 0.01).sum()), sys=list(W.columns), sb=sb, cb=cb))
R = pd.DataFrame(rows)
print(R.drop(columns=["sys", "sb", "cb"]).to_string())
R.to_pickle(f"{out}/r4h_variants.pkl")
