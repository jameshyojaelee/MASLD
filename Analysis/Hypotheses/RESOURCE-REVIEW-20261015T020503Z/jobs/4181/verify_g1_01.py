"""verify2-G1-01: independent check of the Fig 4C Ductular injury (hepatocytes__20) call.

Read-only inputs. Writes only into argv[1].
Questions:
  1. Published scores and q reproduce from per-cell scores + corrected roster?
  2. Lineage cell counts (and total donor cells) for the two GSE189600 donors.
  3. hepatocytes__20 q under: lineage 50/20-cell gate (all programs), donor-total 50-cell gate,
     drop GSE189600 for this program only, leave-one-dataset-out, leave-one-donor-out.
  4. v2-QC-passing hepatocyte counts for GSE189600 donors.
"""
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import stats

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1])
PARQ = BASE / "Analysis/SingleCell/results_gpu_v2/hotspot_modules/cell_scores_all.parquet"
V3 = BASE / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot"
META = BASE / "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
REG = BASE / "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv"
PRIM = BASE / "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/primary_effects.tsv"
V2DIR = BASE / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/per_dataset"
STAGE = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2}
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
pd.set_option("display.max_rows", 200)


def sample_to_donor():
    pairing = {ds: BASE / f"data/{ds}/metadata/donor_pairing.csv"
               for ds in ["GSE244832", "GSE185477", "GSE202379", "GSE136103"]}
    m = {}
    for ds, fp in pairing.items():
        dp = pd.read_csv(fp, dtype=str)
        for _, r in dp.iterrows():
            for s in str(r["rna_srrs"]).split(";"):
                if s.strip():
                    m[s.strip()] = f"{ds}_{r['donor_id']}"
    runs = pd.read_csv(BASE / "Liver_Atlas/metadata/SraRunTable.csv", dtype=str, usecols=["Run", "shortfilename"])
    smp = pd.read_csv(BASE / "Liver_Atlas/metadata/GSE192740_sampleInfo_scRNAseq.tsv", sep="\t", dtype=str,
                      usecols=["characteristics: shortFileName", "title"])
    smp = smp.rename(columns={"characteristics: shortFileName": "shortfilename"})
    lv = runs.merge(smp, on="shortfilename")
    lv["d"] = lv["title"].str.extract(r"(H[0-9]+)", expand=False)
    lv = lv.dropna(subset=["d"])
    m.update(dict(zip(lv["Run"], "Liver_Atlas_" + lv["d"])))
    return m


s2d = sample_to_donor()

# v2 QC-passing cell ids (obs index only)
v2 = []
for ds in ["GSE136103", "GSE174748", "GSE185477", "GSE189600", "GSE202379", "GSE244832", "Liver_Atlas"]:
    with h5py.File(V2DIR / f"{ds}_v2.h5ad", "r") as h:
        obs = h["obs"]
        v2.append(obs[obs.attrs.get("_index", "_index")].asstr()[:])
v2 = pd.Index(np.concatenate(v2))

pf = pq.ParquetFile(PARQ)
parts = []
cellparts = []
for rg in range(pf.metadata.num_row_groups):
    t = pf.read_row_group(rg, columns=["cell_id", "module", "score", "cell_type"]).to_pandas()
    t["sample"] = t["cell_id"].str.split("_", n=1).str[0]
    t["donor"] = t["sample"].map(s2d).fillna(t["sample"])
    t["mod"] = t["module"].str.rsplit("__", n=1).str[1].astype(int)
    t["qc"] = t["cell_id"].isin(v2)
    parts.append(t.groupby(["cell_type", "mod", "donor", "qc"]).score.agg(["sum", "count"]))
    # distinct cells per lineage/donor counted from module 1 rows only
    one = t[t["mod"] == 1]
    cellparts.append(one.groupby(["cell_type", "donor", "qc"]).cell_id.nunique())
    del t
A = pd.concat(parts).groupby(level=[0, 1, 2, 3]).sum().reset_index()
CC = pd.concat(cellparts).groupby(level=[0, 1, 2]).sum().rename("n").reset_index()

allc = A.groupby(["cell_type", "mod", "donor"]).agg(s=("sum", "sum"), n=("count", "sum")).reset_index()
allc["score"] = allc.s / allc.n
qcc = A[A.qc].groupby(["cell_type", "mod", "donor"]).agg(s=("sum", "sum"), n=("count", "sum")).reset_index()
qcc["score"] = qcc.s / qcc.n

meta = pd.read_csv(META, sep="\t", keep_default_na=False, dtype=str)
meta = meta[(meta.exclude_stage_analysis == "FALSE") & meta.disease_stage_coarse.isin(STAGE)]
meta = meta.assign(stage=meta.disease_stage_coarse.map(STAGE)).rename(columns={"sample": "donor"})[["donor", "dataset", "stage"]]
print("[roster] donors", meta.donor.nunique(), meta.groupby("dataset").stage.value_counts().to_dict())

# donor total cells across the 5 scored lineages
tot = CC.groupby("donor").n.sum().rename("n_total_5lineages")

# reproduce published scores
pub = pd.read_csv(V3 / "donor_program_scores_primary.tsv", sep="\t", usecols=["cell_type", "module", "donor", "score"])
pub = pub.rename(columns={"module": "mod"})
chk = pub.merge(allc, on=["cell_type", "mod", "donor"], how="left", suffixes=("_pub", "_me"))
print(f"[regen scores] rows={len(pub)} matched={chk.score_me.notna().sum()} max|diff|={np.nanmax(np.abs(chk.score_pub - chk.score_me)):.3e}")


def ols(d):
    d = d[d.groupby("dataset").donor.transform("size") > 1]  # singleton dataset: own dummy, zero leverage-free info
    X = pd.get_dummies(d.dataset, drop_first=True, dtype=float)
    X.insert(0, "stage", d.stage.astype(float).values)
    X.insert(0, "const", 1.0)
    X = X.to_numpy()
    y = d.score.to_numpy(float)
    if np.linalg.matrix_rank(X) < X.shape[1] or len(y) <= X.shape[1]:
        return None
    XtXi = np.linalg.inv(X.T @ X)
    b = XtXi @ X.T @ y
    e = y - X @ b
    df = len(y) - X.shape[1]
    se = np.sqrt(e @ e / df * XtXi[1, 1])
    h = np.sum((X @ XtXi) * X, axis=1)
    V = XtXi @ (X.T * ((e / (1 - h)) ** 2)) @ X @ XtXi
    return dict(beta=b[1], se=se, p=2 * stats.t.sf(abs(b[1] / se), df),
                hc3_p=2 * stats.t.sf(abs(b[1] / np.sqrt(V[1, 1])), df), n=len(y), n_ds=d.dataset.nunique())


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def family(tab, filt=None):
    rows = []
    for (ct, mod), d in tab.groupby(["cell_type", "mod"]):
        dd = filt(ct, mod, d) if filt else d
        f = ols(dd)
        rows.append(dict(cell_type=ct, mod=mod, **(f or {})))
    r = pd.DataFrame(rows).dropna(subset=["p"])
    r["q"] = bh(r.p)
    r["hc3_q"] = bh(r.hc3_p)
    return r


R_all = allc.merge(meta, on="donor")
R_all = R_all.merge(tot, left_on="donor", right_index=True, how="left")
R_qc = qcc.merge(meta, on="donor")

base = family(R_all)
reg = pd.read_csv(REG, sep="\t")
cmp_ = base.merge(reg[["cell_type", "module", "primary_qvalue", "primary_hc3_qvalue", "primary_n_donors"]],
                  left_on=["cell_type", "mod"], right_on=["cell_type", "module"])
print(f"[regen q] programs={len(base)} matched={len(cmp_)} max|q diff|={np.abs(cmp_.q - cmp_.primary_qvalue).max():.3e} "
      f"max|hc3 q diff|={np.abs(cmp_.hc3_q - cmp_.primary_hc3_qvalue).max():.3e} n equal={(cmp_.n == cmp_.primary_n_donors).all()}")
print("[published] programs with registry primary_qvalue<0.05:")
print(reg.loc[reg.primary_qvalue < 0.05, ["cell_type", "module", "module_name", "primary_beta", "primary_qvalue", "primary_hc3_qvalue", "primary_n_donors"]].to_string(index=False))

foc = lambda r: r[(r.cell_type == "hepatocytes") & (r["mod"] == 20)][["beta", "se", "p", "q", "hc3_q", "n", "n_ds"]]
res = []


def rec(name, r):
    f = foc(r).iloc[0].to_dict()
    f.update(arm=name, n_calls_q05=int((r.q < 0.05).sum()), calls=";".join(f"{a}__{b}" for a, b in zip(r[r.q < 0.05].cell_type, r[r.q < 0.05]["mod"])))
    res.append(f)


rec("published_no_gate", base)
rec("lineage_gate50_all_programs", family(R_all, lambda ct, m, d: d[d.n >= 50]))
rec("lineage_gate20_all_programs", family(R_all, lambda ct, m, d: d[d.n >= 20]))
rec("donor_total_gate50_all_programs", family(R_all, lambda ct, m, d: d[d.n_total_5lineages >= 50]))
rec("drop_GSE189600_hep20_only", family(R_all, lambda ct, m, d: d[d.dataset != "GSE189600"] if (ct, m) == ("hepatocytes", 20) else d))
rec("drop_GSE189600_all_programs", family(R_all, lambda ct, m, d: d[d.dataset != "GSE189600"]))
for ds in ["GSE202379", "GSE244832", "GSE174748", "GSE185477"]:
    rec(f"drop_{ds}_hep20_only", family(R_all, lambda ct, m, d, ds=ds: d[d.dataset != ds] if (ct, m) == ("hepatocytes", 20) else d))
rec("v2qc_no_gate", family(R_qc))
rec("v2qc_lineage_gate50", family(R_qc, lambda ct, m, d: d[d.n >= 50]))
arms = pd.DataFrame(res)[["arm", "beta", "se", "p", "q", "hc3_q", "n", "n_ds", "n_calls_q05", "calls"]]
print("[arms] hepatocytes__20 under each arm; BH over the estimable 117-program family")
print(arms.to_string(index=False))
arms.to_csv(OUT / "arms_hep20.tsv", sep="\t", index=False)

# leave-one-donor-out for hepatocytes__20 (other 116 unchanged)
hep_donors = R_all[(R_all.cell_type == "hepatocytes") & (R_all["mod"] == 20)].donor.tolist()
lodo = []
h20 = R_all[(R_all.cell_type == "hepatocytes") & (R_all["mod"] == 20)]
is20 = ((base.cell_type == "hepatocytes") & (base["mod"] == 20)).to_numpy()
for dn in hep_donors:
    f = ols(h20[h20.donor != dn])
    p = base.p.to_numpy().copy()
    p[is20] = f["p"]
    lodo.append(dict(donor=dn, beta=f["beta"], q=bh(p)[is20][0]))
lodo = pd.DataFrame(lodo).merge(meta, on="donor")
lodo.to_csv(OUT / "lodo_hep20.tsv", sep="\t", index=False)
print(f"[LOO donor] donors={len(lodo)}; q>=0.05 after dropping: {int((lodo.q >= 0.05).sum())}")
print(lodo.sort_values("q", ascending=False).head(8).to_string(index=False))

# GSE189600 donor detail
g = R_all[(R_all.dataset == "GSE189600") & (R_all["mod"] == 20) & (R_all.cell_type == "hepatocytes")][["donor", "stage", "n", "score", "n_total_5lineages"]]
gq = R_qc[(R_qc.dataset == "GSE189600") & (R_qc["mod"] == 20) & (R_qc.cell_type == "hepatocytes")][["donor", "n", "score"]].rename(columns={"n": "n_v2qc", "score": "score_v2qc"})
print("[GSE189600] hepatocytes__20 donor detail (n = hepatocyte cells; n_total = cells in 5 scored lineages):")
print(g.merge(gq, on="donor", how="left").to_string(index=False))
print("[GSE189600] cells per lineage:")
print(CC[CC.donor.isin(["SRR22677587", "SRR22677595"])].pivot_table(index=["donor"], columns=["cell_type", "qc"], values="n", fill_value=0).to_string())

# hepatocyte cell counts across roster donors: how many below 50 / 20
h = R_all[(R_all.cell_type == "hepatocytes") & (R_all["mod"] == 20)]
print(f"[hep counts] roster donors with hepatocytes: {len(h)}; <50: {(h.n < 50).sum()} {h[h.n < 50].donor.tolist()}; <20: {(h.n < 20).sum()}")
# within-dataset stage contrast for GSE189600 and slope from the other informative datasets
for ds in ["GSE189600", "GSE202379", "GSE244832"]:
    d = h[h.dataset == ds]
    print(f"[within {ds}] n={len(d)} stage means:", d.groupby("stage").score.mean().round(4).to_dict())
