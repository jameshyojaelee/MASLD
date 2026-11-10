"""U8a: donor means of stored per-cell Hotspot scores over QC-passing cells only (AN7 S1 input).

QC-passing = mitochondrial fraction < 20%, >= 200 genes, >= 500 counts, Scrublet < 0.25 (METHODS:343),
applied here to the obs columns of the harmonized v2 per-dataset h5ad files; also reports whether the
v2 obs set is already filtered.  Compares with AN7 out/donor_scores_v2qc_cells.tsv.
"""
import sys
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
V2 = ROOT / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/per_dataset"
CELL = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules/cell_scores_all.parquet"
AN7 = ROOT / "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN7/20261021T124200Z/out/donor_scores_v2qc_cells.tsv"
PAIR = {d: ROOT / f"data/{d}/metadata/donor_pairing.csv" for d in ["GSE244832", "GSE185477", "GSE202379", "GSE136103"]}


def col(o, name):
    x = o[name]
    if isinstance(x, h5py.Group):  # categorical
        cats = x["categories"][:]
        return cats[x["codes"][:]]
    return x[:]


ids_all, ids_pass = [], []
for f in sorted(V2.glob("*_v2.h5ad")):
    with h5py.File(f, "r") as h:
        o = h["obs"]
        idx = np.array([s.decode() if isinstance(s, bytes) else s for s in o[o.attrs.get("_index", "_index")][:]])
        mt = col(o, "pct_counts_mt").astype(float)
        ng = col(o, "n_genes_by_counts").astype(float)
        tc = col(o, "total_counts").astype(float)
        sc = col(o, "scrublet_score").astype(float)
        ok = (mt < 20) & (ng >= 200) & (tc >= 500) & (sc < 0.25)
        print(f.name, len(idx), "pass", int(ok.sum()), "max mt %.2f min genes %d min counts %d max scrublet %.3f"
              % (np.nanmax(mt), np.nanmin(ng), np.nanmin(tc), np.nanmax(sc)), flush=True)
        ids_all.append(idx); ids_pass.append(idx[ok])
ids_all = set(np.concatenate(ids_all)); ids_pass = set(np.concatenate(ids_pass))
print("v2 cells", len(ids_all), "passing own thresholds", len(ids_pass))

run2donor = {}
for ds, fp in PAIR.items():
    for _, r in pd.read_csv(fp, dtype=str).iterrows():
        for s in str(r["rna_srrs"]).split(";"):
            if s.strip():
                run2donor[s.strip()] = f"{ds}_{r['donor_id']}"

tab = pq.read_table(CELL, columns=["cell_id", "module", "score", "cell_type"])
df = pd.DataFrame({"cell_id": tab.column("cell_id").to_pandas().astype("category"),
                   "module": tab.column("module").to_pandas().astype("category"),
                   "score": tab.column("score").to_numpy(),
                   "cell_type": tab.column("cell_type").to_pandas().astype("category")})
del tab
cats = pd.Series(df["cell_id"].cat.categories)
in_all = cats.isin(ids_all).values; in_pass = cats.isin(ids_pass).values
print("scored cells", len(cats), "in v2 obs", int(in_all.sum()), "passing", int(in_pass.sum()))
codes = df["cell_id"].cat.codes.values
for label, keep in [("v2obs", in_all), ("own_thresholds", in_pass)]:
    sub = df.loc[keep[codes], ["cell_id", "module", "score", "cell_type"]].copy()
    smp = cats.str.split("_", n=1).str[0].values[sub["cell_id"].cat.codes.values]
    sub["donor"] = pd.Series(smp).map(run2donor).fillna(pd.Series(smp)).values
    g = sub.groupby(["cell_type", "module", "donor"], observed=True)["score"].agg(["mean", "size"]).reset_index()
    g["module"] = g["module"].astype(str).str.split("__").str[1].astype(int)
    g.to_csv(OUT / f"donor_scores_qc_{label}.tsv", sep="\t", index=False)
    a = pd.read_csv(AN7, sep="\t")
    m = a.merge(g, on=["cell_type", "module", "donor"], how="outer", indicator=True)
    both = m["_merge"] == "both"
    print(label, "AN7 rows", len(a), "own rows", len(g), "matched", int(both.sum()),
          "max |score diff|", float((m.loc[both, "score"] - m.loc[both, "mean"]).abs().max()),
          "n_cells mismatches", int((m.loc[both, "n_cells"] != m.loc[both, "size"]).sum()), flush=True)
