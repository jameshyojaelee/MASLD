"""Blind check of R08-08 (TranscriptFormer clean-subset weightings) and R08-09 (Fig 7F per-study ranges).
Sources: encoder-field-score-21314196/evaluation/{study_metrics.tsv, paired_deltas_by_held_out_study.tsv,
bootstrap_distributions.npz}. Donor counts read from study_metrics.tsv, not the Fig 7F grid.
No new resampling: intervals are reaggregated from the saved 10,000 paired draws."""
import json
import numpy as np, pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
EV = f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions/encoder-field-score-21314196/evaluation"
OUT = f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F7b"
res = {}
sm = pd.read_csv(f"{EV}/study_metrics.tsv", sep="\t")
don = sm.drop_duplicates("held_out_study").set_index("held_out_study").donors.astype(int)
rows_ = sm.drop_duplicates("held_out_study").set_index("held_out_study").rows.astype(int)
res["donors"] = don.to_dict(); res["rows"] = rows_.to_dict()
f1 = sm.set_index(["model_id", "held_out_study"]).donor_class_balanced_macro_f1
clean = ["GSE136103", "GSE174748", "GSE189600", "GSE202379", "GSE244832"]
z = np.load(f"{EV}/bootstrap_distributions.npz")
tf = {}
for head in ["linear", "two_layer_mlp"]:
    m, r = f"transcriptformer_tf_sapiens::{head}", f"scvi_liver_latent::{head}"
    per = {s: float(f1[m, s] - f1[r, s]) for s in clean}
    draws = {s: z[f"{m}::study_macro_f1::{s}"] - z[f"{r}::study_macro_f1::{s}"] for s in clean}
    out = {"per_study": per}
    for name, ss, w in [("equal5", clean, {s: 1 for s in clean}),
                        ("donor5", clean, {s: don[s] for s in clean}),
                        ("equal4_noGSE189600", [s for s in clean if s != "GSE189600"], {s: 1 for s in clean})]:
        W = np.array([w[s] for s in ss], float); W /= W.sum()
        pt = float(sum(W[i] * per[s] for i, s in enumerate(ss)))
        dv = sum(W[i] * draws[s] for i, s in enumerate(ss))
        lo, hi = np.percentile(dv, [2.5, 97.5])
        out[name] = dict(point=pt, lo=float(lo), hi=float(hi), excludes0=bool(lo > 0 or hi < 0), donors=int(sum(don[s] for s in ss)))
    out["GSE189600_share_of_equal5_mean"] = per["GSE189600"] / 5
    tf[head] = out
res["R08-08"] = tf

# ---------- R08-09 from the per-study paired-delta source ----------
pd_ = pd.read_csv(f"{EV}/paired_deltas_by_held_out_study.tsv", sep="\t")
pd_ = pd_[(pd_.metric == "donor_class_balanced_macro_f1")]
blocks = [b for b in pd_.model_id.unique() if "::" in b and b.split("::")[0] != "scvi_liver_latent"]
rows = []
for b in blocks:
    d = pd_[pd_.model_id == b].set_index("held_out_study")
    assert len(d) == 7 and (d.reference_id == f"scvi_liver_latent::{b.split('::')[1]}").all(), b
    x = d.oriented_delta
    o = x.drop("GSE189600")
    eq = x.mean(); dw = np.average(x, weights=don.loc[x.index])
    rows.append(dict(block=b, g189600=x["GSE189600"], o_min=o.min(), o_max=o.max(), o_pos=int((o > 0).sum()),
                     eq=eq, dw=dw, flip=bool(np.sign(eq) != np.sign(dw))))
t = pd.DataFrame(rows)
big = t.g189600 < -0.4
res["R08-09"] = dict(
    n_blocks=len(t),
    big_blocks=t[big].block.tolist(), big_range=[t[big].g189600.min(), t[big].g189600.max()],
    small_range=[t[~big].g189600.min(), t[~big].g189600.max()],
    other6_range=[t.o_min.min(), t.o_max.max()], other6_pos=int(t.o_pos.sum()), other6_total=6 * len(t),
    sign_flips=t[t.flip].block.tolist(),
    table=t.to_dict(orient="records"))
json.dump(res, open(f"{OUT}/chk_r0808_r0809.json", "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in res.items()}, indent=1, default=float))
