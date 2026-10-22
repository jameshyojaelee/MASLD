"""Independent re-derivation of the Fig 4C primary stage model (review R04).

Inputs (read-only):
  v3 donor scores (64-donor era) + corrected donor metadata (2026-10-14).
Model: score ~ stage_ordinal + factor(dataset), stages Healthy=0/Steatosis=1/Steatohepatitis=2,
OLS and HC3 SE, BH over 117 programs.
Also: per-dataset design, informative donors, residual SD by dataset, within-dataset slopes,
and a within-dataset stage-label permutation null (seeded) for all 117 programs jointly.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
V3 = ROOT + "/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot"
META = ROOT + "/Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
META_OLD = ROOT + "/Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/originals/donor_metadata_stage_corrected_donor.tsv"
PUB = ROOT + "/figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv"
OUT = sys.argv[1]
NPERM = int(sys.argv[2]) if len(sys.argv) > 2 else 0
STAGES = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2}

scores = pd.read_csv(V3 + "/donor_program_scores_primary.tsv", sep="\t",
                     usecols=["cell_type", "module", "donor", "score"])


def bh(p, n):
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out


def load(meta_path):
    m = pd.read_csv(meta_path, sep="\t", keep_default_na=False)
    m = m[m.exclude_stage_analysis.astype(str).str.upper().isin(["FALSE", "F", "0"])]
    m = m[m.disease_stage_coarse.isin(STAGES)]
    m = m.assign(stage=m.disease_stage_coarse.map(STAGES))[["sample", "dataset", "stage"]]
    return scores.merge(m.rename(columns={"sample": "donor"}), on="donor")


def fit(d):
    ds = pd.get_dummies(d.dataset, drop_first=True, dtype=float).to_numpy()
    x = np.column_stack([np.ones(len(d)), d.stage.to_numpy(float), ds])
    y = d.score.to_numpy(float)
    inv = np.linalg.inv(x.T @ x)
    b = inv @ x.T @ y
    r = y - x @ b
    df = len(y) - np.linalg.matrix_rank(x)
    se = np.sqrt(r @ r / df * inv[1, 1])
    h = np.einsum("ij,jk,ik->i", x, inv, x)
    hc = inv @ (x.T @ ((r / (1 - h))[:, None] ** 2 * x)) @ inv
    # informative donors: those whose dataset has >1 stage level
    nlev = d.groupby("dataset").stage.transform("nunique")
    return dict(beta=b[1], se=se, p=2 * stats.t.sf(abs(b[1] / se), df),
                hc3_p=2 * stats.t.sf(abs(b[1] / np.sqrt(hc[1, 1])), df), df=df,
                n=len(y), n_ds=d.dataset.nunique(), n_informative=int((nlev > 1).sum()),
                resid=pd.Series(r, index=d.index))


rows = {}
for label, path in [("corrected", META), ("original", META_OLD)]:
    dd = load(path)
    res = []
    for (ct, mod), d in dd.groupby(["cell_type", "module"]):
        f = fit(d)
        f.pop("resid")
        res.append(dict(cell_type=ct, module=mod, **f))
    res = pd.DataFrame(res)
    assert len(res) == 117, len(res)
    res["q"] = bh(res.p, 117)
    res["hc3_q"] = bh(res.hc3_p, 117)
    rows[label] = res
    print(f"[{label}] donors max={dd.donor.nunique()} programs={len(res)} "
          f"q<0.05={int((res.q < .05).sum())} hc3_q<0.05={int((res.hc3_q < .05).sum())}")

cor = rows["corrected"]
pub = pd.read_csv(PUB, sep="\t")
cmp_ = cor.merge(pub[["cell_type", "module", "module_name", "primary_beta", "primary_qvalue",
                      "primary_hc3_qvalue", "primary_n_donors"]], on=["cell_type", "module"])
cmp_["abs_dq"] = (cmp_.q - cmp_.primary_qvalue).abs()
cmp_["abs_dhc3q"] = (cmp_.hc3_q - cmp_.primary_hc3_qvalue).abs()
print("max |q - published q| =", cmp_.abs_dq.max(), " max |hc3 q diff| =", cmp_.abs_dhc3q.max())
print("n_donors match:", (cmp_.n == cmp_.primary_n_donors).all())
focal = cmp_[(cmp_.cell_type == "hepatocytes") & cmp_.module.isin([8, 20])]
print(focal[["module", "module_name", "beta", "q", "hc3_q", "n", "n_ds", "n_informative",
             "primary_qvalue", "primary_hc3_qvalue"]].to_string(index=False))
orig = rows["original"]
of = orig[(orig.cell_type == "hepatocytes") & orig.module.isin([8, 20])]
print("original-metadata focal:\n", of[["module", "beta", "q", "hc3_q", "n"]].to_string(index=False))
print("informative donors by lineage (corrected):")
print(cor.groupby("cell_type")[["n", "n_informative", "n_ds"]].agg(["min", "max"]).to_string())
cmp_.to_csv(OUT + "/rederived_primary_effects.tsv", sep="\t", index=False)

# Design table and per-dataset residual SD / within-dataset slopes for focal programs.
dd = load(META)
design = dd.drop_duplicates(["cell_type", "donor"]).pivot_table(
    index=["cell_type", "dataset"], columns="stage", values="donor", aggfunc="count", fill_value=0)
print("stage x dataset donors per lineage (0=Healthy,1=Steatosis,2=Steatohepatitis):")
print(design.to_string())
design.to_csv(OUT + "/stage_by_dataset_by_lineage.tsv", sep="\t")
for mod in [20, 8]:
    d = dd[(dd.cell_type == "hepatocytes") & (dd.module == mod)].copy()
    f = fit(d)
    d["resid"] = f["resid"]
    print(f"hep module {mod}: residual SD by dataset")
    print(d.groupby("dataset").agg(n=("resid", "size"), resid_sd=("resid", "std")).to_string())
    for ds, g in d.groupby("dataset"):
        if g.stage.nunique() > 1 and len(g) > 2:
            sl = stats.linregress(g.stage, g.score)
            print(f"  within {ds}: n={len(g)} slope={sl.slope:.4f} se={sl.stderr:.4f} p={sl.pvalue:.4g}")
    # Levene / Brown-Forsythe across datasets with >=3 donors
    groups = [g.resid.to_numpy() for _, g in d.groupby("dataset") if len(g) >= 3]
    print("  Brown-Forsythe across datasets (n>=3):", stats.levene(*groups, center="median"))
    # leave-one-dataset-out
    for ds in sorted(d.dataset.unique()):
        g = d[d.dataset != ds]
        if g.groupby("dataset").stage.nunique().max() > 1:
            ff = fit(g)
            print(f"  LODO -{ds}: beta={ff['beta']:.4f} p={ff['p']:.4g} hc3_p={ff['hc3_p']:.4g} n={ff['n']}")

if NPERM:
    rng = np.random.default_rng(20261015)
    donors = dd.drop_duplicates("donor")[["donor", "dataset", "stage"]].reset_index(drop=True)
    progs = list(dd.groupby(["cell_type", "module"]))
    obs_t = []
    prep = []
    for (ct, mod), d in progs:
        d = d.set_index("donor").loc[:, ["dataset", "score"]]
        y = d.score - d.groupby("dataset").score.transform("mean")
        prep.append((d.index.to_numpy(), d.dataset.to_numpy(), y.to_numpy(), d.score.to_numpy()))
    donor_idx = {k: i for i, k in enumerate(donors.donor)}
    base_stage = donors.stage.to_numpy(float)
    ds_arr = donors.dataset.to_numpy()

    def tstats(stage_vec):
        out = np.empty(len(prep))
        for i, (dn, dsv, yc, _) in enumerate(prep):
            s = stage_vec[[donor_idx[x] for x in dn]]
            sc = s - pd.Series(s).groupby(dsv).transform("mean").to_numpy()
            sxx = sc @ sc
            b = (sc @ yc) / sxx
            r = yc - b * sc
            k = len(np.unique(dsv))
            df = len(yc) - k - 1
            out[i] = b / np.sqrt(r @ r / df / sxx)
        return out

    obs = tstats(base_stage)
    groups = {ds: np.where(ds_arr == ds)[0] for ds in np.unique(ds_arr)}
    exceed = np.zeros(len(prep))
    maxexceed = np.zeros(len(prep))
    for _ in range(NPERM):
        perm = base_stage.copy()
        for ds, idx in groups.items():
            perm[idx] = base_stage[rng.permutation(idx)]
        t = tstats(perm)
        exceed += np.abs(t) >= np.abs(obs) - 1e-12
        maxexceed += np.abs(t).max() >= np.abs(obs) - 1e-12
    pperm = (exceed + 1) / (NPERM + 1)
    pfwer = (maxexceed + 1) / (NPERM + 1)
    keys = pd.DataFrame([k for k, _ in progs], columns=["cell_type", "module"])
    keys["t_obs"] = obs
    keys["perm_p"] = pperm
    keys["perm_q_bh117"] = bh(pperm, 117)
    keys["maxT_fwer_p"] = pfwer
    keys.to_csv(OUT + "/within_dataset_permutation.tsv", sep="\t", index=False)
    print(f"permutation null: {NPERM} within-dataset stage permutations, seed 20261015, floor {1/(NPERM+1):.2e}")
    print(keys.sort_values("perm_p").head(6).to_string(index=False))
