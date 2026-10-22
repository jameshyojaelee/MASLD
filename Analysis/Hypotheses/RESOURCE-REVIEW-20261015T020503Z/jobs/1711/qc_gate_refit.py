"""Review G1: does the Fig 4C primary call survive a per-donor cell gate or the harmonized v2 cell QC?

Read-only inputs:
  cell_scores_all.parquet (per-cell Hotspot scores; frozen weights; not recomputed)
  v2 per-dataset h5ad obs names (cells passing MT<20%, min_genes 200, min_counts 500, scrublet 0.25)
  V3 donor_program_scores_primary.tsv (published pooled-cell donor scores, used only to verify regeneration)
  corrected donor metadata 2026-10-14 (the published primary roster)
  published program registry (Fig 4C q values)
Model (copied from RUN/scratch/review-R04/rederive_primary.py, which matched published q to 2e-14):
  score ~ stage_ordinal + dataset, OLS t-test and HC3, BH over the 117 programs.
Run-to-donor map copied from Analysis/SingleCell/scripts/hotspot_modules/511_pooled_cell_donor_scores.py.
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
V2DIR = BASE / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/per_dataset"
V3 = BASE / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot"
META = BASE / "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
PUB = BASE / "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv"
STAGES = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2}
DATASETS = ["GSE136103", "GSE174748", "GSE185477", "GSE189600", "GSE202379", "GSE244832", "Liver_Atlas"]
pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 30)

# ---- run -> donor map (verbatim logic of 511_pooled_cell_donor_scores.py) ----
DONOR_PAIRING = {
    "GSE244832": "data/GSE244832/metadata/donor_pairing.csv",
    "GSE185477": "data/GSE185477/metadata/donor_pairing.csv",
    "GSE202379": "data/GSE202379/metadata/donor_pairing.csv",
    "GSE136103": "data/GSE136103/metadata/donor_pairing.csv",
}


def build_sample_to_donor():
    m = {}
    for ds, rel in DONOR_PAIRING.items():
        dp = pd.read_csv(BASE / rel, dtype=str)
        for _, row in dp.iterrows():
            donor = f"{ds}_{row['donor_id']}"
            for srr in str(row["rna_srrs"]).split(";"):
                srr = srr.strip()
                if srr:
                    m[srr] = donor
    runs = pd.read_csv(BASE / "Liver_Atlas/metadata/SraRunTable.csv", dtype=str, usecols=["Run", "shortfilename"])
    samples = pd.read_csv(BASE / "Liver_Atlas/metadata/GSE192740_sampleInfo_scRNAseq.tsv", sep="\t", dtype=str,
                          usecols=["characteristics: shortFileName", "title"]).rename(
        columns={"characteristics: shortFileName": "shortfilename"})
    liver = runs.merge(samples, on="shortfilename", how="inner")
    liver["donor_short"] = liver["title"].map(lambda v: (mm.group(0) if (mm := re.search(r"H[0-9]+", v)) else None))
    liver = liver.dropna(subset=["donor_short"])
    assert liver["Run"].nunique() == 48 and liver["donor_short"].nunique() == 19
    lm = dict(zip(liver["Run"], "Liver_Atlas_" + liver["donor_short"]))
    assert not set(m).intersection(lm)
    m.update(lm)
    return m


s2d = build_sample_to_donor()
print(f"[map] {len(s2d)} run ids mapped to donors")

# ---- v2 QC-passing cell ids (obs names only; X never read) ----
v2_names = []
v2_sample_ds = {}
v2_n = {}
for ds in DATASETS:
    with h5py.File(V2DIR / f"{ds}_v2.h5ad", "r") as h:
        obs = h["obs"]
        idx = obs.attrs.get("_index", "_index")
        names = obs[idx].asstr()[:]
        v2_names.append(names)
        v2_n[ds] = len(names)
        cats = obs["sample"]["categories"].asstr()[:]
        for s in cats:
            v2_sample_ds[s] = ds
v2_index = pd.Index(np.concatenate(v2_names))
assert v2_index.is_unique, "v2 obs names not unique across datasets"
print("[v2] cells per dataset:", v2_n, "total", sum(v2_n.values()))

# ---- scan per-cell scores ----
pf = pq.ParquetFile(PARQ)
aggs = []
cells = []
for rg in range(pf.metadata.num_row_groups):
    t = pf.read_row_group(rg, columns=["cell_id", "module", "score", "cell_type"]).to_pandas()
    t["sample"] = t["cell_id"].str.split("_", n=1).str[0]
    t["donor"] = t["sample"].map(s2d).fillna(t["sample"])
    t["qc"] = t["cell_id"].isin(v2_index)
    aggs.append(t.groupby(["cell_type", "module", "donor", "qc"]).score.agg(["sum", "count"]))
    one = t[t["module"].str.endswith("__1")]
    cells.append(one.groupby(["cell_type", "sample", "qc"]).size())
    del t
A = pd.concat(aggs).groupby(level=[0, 1, 2, 3]).sum().reset_index()
A["module"] = A["module"].str.split("__").str[1].astype(int)
C = pd.concat(cells).groupby(level=[0, 1, 2]).sum().rename("n").reset_index()

allc = A.groupby(["cell_type", "module", "donor"]).agg(s=("sum", "sum"), n=("count", "sum")).reset_index()
allc["score"] = allc.s / allc.n
qcc = A[A.qc].groupby(["cell_type", "module", "donor"]).agg(s=("sum", "sum"), n=("count", "sum")).reset_index()
qcc["score"] = qcc.s / qcc.n

# Cell counts must be identical across modules within a lineage.
spread = allc.groupby(["cell_type", "donor"]).n.agg(lambda x: x.max() - x.min()).max()
spread_qc = qcc.groupby(["cell_type", "donor"]).n.agg(lambda x: x.max() - x.min()).max()
print(f"[check] max within-lineage spread of per-module cell counts: all={spread} qc={spread_qc}")

# ---- regeneration check against the published pooled-cell scores ----
pub_scores = pd.read_csv(V3 / "donor_program_scores_primary.tsv", sep="\t", usecols=["cell_type", "module", "donor", "score"])
chk = pub_scores.merge(allc[["cell_type", "module", "donor", "score"]], on=["cell_type", "module", "donor"],
                       how="left", suffixes=("_pub", "_regen"))
print(f"[regen] published rows={len(pub_scores)} matched={chk.score_regen.notna().sum()} "
      f"max|diff|={np.nanmax(np.abs(chk.score_pub - chk.score_regen)):.3e}")

# ---- overlap of v1 atlas cells with v2 QC-passing cells, by dataset ----
C["dataset"] = C["sample"].map(v2_sample_ds).fillna("not_in_v2_sample_list")
ov = C.pivot_table(index="dataset", columns="qc", values="n", aggfunc="sum", fill_value=0)
ov.columns = ["v1_fail_v2" if not c else "v1_pass_v2" for c in ov.columns]
ov["v1_cells_in_5_lineages"] = ov.sum(axis=1)
ov["frac_pass"] = ov.v1_pass_v2 / ov.v1_cells_in_5_lineages
print("[overlap] v1 atlas cells in the 5 scored lineages vs v2 QC-passing set:")
print(ov.to_string())
ov.to_csv(OUT / "v1_cells_passing_v2_by_dataset.tsv", sep="\t")


# ---- model (from review-R04 template) ----
def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out


def fit(d):
    # A dataset with one donor is absorbed by its own dummy (zero residual, leverage 1);
    # dropping it leaves beta, SE and df unchanged and avoids 0/0 in HC3.
    keep = d.groupby("dataset").donor.transform("size") > 1
    d = d[keep]
    ds = pd.get_dummies(d.dataset, drop_first=True, dtype=float).to_numpy()
    x = np.column_stack([np.ones(len(d)), d.stage.to_numpy(float), ds])
    y = d.score.to_numpy(float)
    if np.linalg.matrix_rank(x) < x.shape[1] or len(y) - x.shape[1] < 1:
        return None
    inv = np.linalg.inv(x.T @ x)
    b = inv @ x.T @ y
    r = y - x @ b
    df = len(y) - x.shape[1]
    se = np.sqrt(r @ r / df * inv[1, 1])
    h = np.einsum("ij,jk,ik->i", x, inv, x)
    hc = inv @ (x.T @ ((r / (1 - h))[:, None] ** 2 * x)) @ inv
    nlev = d.groupby("dataset").stage.transform("nunique")
    return dict(beta=b[1], se=se, p=2 * stats.t.sf(abs(b[1] / se), df),
                hc3_p=2 * stats.t.sf(abs(b[1] / np.sqrt(hc[1, 1])), df), df=df, n=len(y),
                n_ds=d.dataset.nunique(), n_informative=int((nlev > 1).sum()))


meta = pd.read_csv(META, sep="\t", keep_default_na=False)
meta = meta[meta.exclude_stage_analysis.astype(str).str.upper().isin(["FALSE", "F", "0"])]
meta = meta[meta.disease_stage_coarse.isin(STAGES)]
meta = meta.assign(stage=meta.disease_stage_coarse.map(STAGES))[["sample", "dataset", "stage"]].rename(
    columns={"sample": "donor"})
print(f"[meta] primary roster donors: {meta.donor.nunique()}")

roster_all = allc.merge(meta, on="donor")
roster_qc = qcc.merge(meta, on="donor")

# ---- per-lineage gate failures among roster donors ----
rows = []
for label, tab in [("all_cells", roster_all), ("v2_qc_cells", roster_qc)]:
    per = tab.drop_duplicates(["cell_type", "donor"])
    for ct, g in per.groupby("cell_type"):
        base = roster_all.drop_duplicates(["cell_type", "donor"])
        n_roster = int((base.cell_type == ct).sum())
        rows.append(dict(cells=label, cell_type=ct, roster_donors=n_roster, donors_with_cells=len(g),
                         donors_zero=n_roster - len(g), lt20=int((g.n < 20).sum()) + (n_roster - len(g)),
                         lt50=int((g.n < 50).sum()) + (n_roster - len(g)), min_cells=int(g.n.min()),
                         median_cells=float(g.n.median())))
gates = pd.DataFrame(rows)
print("[gate] roster donor x lineage cell support (lt20/lt50 include donors with zero cells):")
print(gates.to_string(index=False))
gates.to_csv(OUT / "gate_failures_by_lineage.tsv", sep="\t", index=False)
low = roster_all.drop_duplicates(["cell_type", "donor"])
low = low[low.n < 50][["cell_type", "donor", "dataset", "stage", "n"]].merge(
    roster_qc.drop_duplicates(["cell_type", "donor"])[["cell_type", "donor", "n"]].rename(columns={"n": "n_qc"}),
    on=["cell_type", "donor"], how="left")
print("[gate] roster donor x lineage rows with <50 cells (all cells), with v2-QC-passing count:")
print(low.sort_values(["cell_type", "n"]).to_string(index=False))
qlow = roster_qc.drop_duplicates(["cell_type", "donor"])
print("[gate] hepatocyte roster donors under v2 QC with <50 passing cells:")
hb = roster_all[(roster_all.cell_type == "hepatocytes")].drop_duplicates("donor")[["donor", "dataset", "stage", "n"]]
hb = hb.merge(qlow[qlow.cell_type == "hepatocytes"][["donor", "n"]].rename(columns={"n": "n_qc"}), on="donor", how="left")
hb["n_qc"] = hb.n_qc.fillna(0).astype(int)
print(hb[hb.n_qc < 50].sort_values("n_qc").to_string(index=False))


# ---- arms ----
def run_arm(name, tab, min_cells, drop_ds=None):
    res = []
    for (ct, mod), d in tab.groupby(["cell_type", "module"]):
        d = d[d.n >= min_cells]
        if drop_ds:
            d = d[d.dataset != drop_ds]
        f = fit(d)
        if f is None:
            res.append(dict(cell_type=ct, module=mod, estimable=False))
        else:
            res.append(dict(cell_type=ct, module=mod, estimable=True, **f))
    res = pd.DataFrame(res)
    est = res.estimable
    res.loc[est, "q"] = bh(res.loc[est, "p"])
    res.loc[est, "hc3_q"] = bh(res.loc[est, "hc3_p"])
    res["arm"] = name
    return res


arms = [
    ("A0_published_all_cells_no_gate", roster_all, 1, None),
    ("A1_all_cells_gate50", roster_all, 50, None),
    ("A2_all_cells_gate20", roster_all, 20, None),
    ("A3_v2qc_cells_no_gate", roster_qc, 1, None),
    ("A4_v2qc_cells_gate50", roster_qc, 50, None),
    ("A5_v2qc_cells_gate20", roster_qc, 20, None),
    ("R_all_cells_drop_GSE189600", roster_all, 1, "GSE189600"),
]
out = []
summ = []
for name, tab, mc, dd in arms:
    r = run_arm(name, tab, mc, dd)
    out.append(r)
    calls = r[r.q < 0.05]
    summ.append(dict(arm=name, estimable=int(r.estimable.sum()), bh117_calls=len(calls),
                     hc3_bh117_calls=int((r.hc3_q < 0.05).sum()),
                     called=";".join(f"{a}__{b}" for a, b in zip(calls.cell_type, calls.module)),
                     hep_n_max=int(r[r.cell_type == "hepatocytes"].n.max())))
eff = pd.concat(out)
eff.to_csv(OUT / "arm_effects_117.tsv", sep="\t", index=False)
summ = pd.DataFrame(summ)
print("[arms] BH over 117 programs (OLS) and BH over 117 (HC3):")
print(summ.to_string(index=False))
summ.to_csv(OUT / "arm_summary.tsv", sep="\t", index=False)
foc = eff[(eff.cell_type == "hepatocytes") & eff.module.isin([8, 20])][
    ["arm", "module", "beta", "se", "p", "q", "hc3_p", "hc3_q", "n", "n_ds", "n_informative"]]
print("[arms] focal programs (hepatocytes__20 Ductular injury BICC1; hepatocytes__8 Stromal ECM IGFBP7):")
print(foc.to_string(index=False))
foc.to_csv(OUT / "arm_focal.tsv", sep="\t", index=False)

pub = pd.read_csv(PUB, sep="\t")
a0 = eff[eff.arm == "A0_published_all_cells_no_gate"].merge(
    pub[["cell_type", "module", "primary_qvalue", "primary_hc3_qvalue", "primary_n_donors"]], on=["cell_type", "module"])
print(f"[regen] A0 vs published registry: max|q diff|={np.abs(a0.q - a0.primary_qvalue).max():.3e} "
      f"max|hc3 q diff|={np.abs(a0.hc3_q - a0.primary_hc3_qvalue).max():.3e} "
      f"n match={(a0.n == a0.primary_n_donors).all()}")

# Donor-level hepatocytes__20 scores under all vs QC cells, for the record.
h20 = roster_all[(roster_all.cell_type == "hepatocytes") & (roster_all.module == 20)][
    ["donor", "dataset", "stage", "n", "score"]].merge(
    roster_qc[(roster_qc.cell_type == "hepatocytes") & (roster_qc.module == 20)][["donor", "n", "score"]].rename(
        columns={"n": "n_qc", "score": "score_qc"}), on="donor", how="left")
h20.to_csv(OUT / "hep20_donor_scores_all_vs_qc.tsv", sep="\t", index=False)
print("[h20] Spearman(all-cell score, QC-cell score) over roster donors with QC cells:",
      stats.spearmanr(h20.dropna().score, h20.dropna().score_qc))
print("[h20] by dataset median |score_qc - score|:")
print(h20.dropna().assign(d=lambda x: (x.score_qc - x.score).abs()).groupby("dataset").d.median().to_string())
