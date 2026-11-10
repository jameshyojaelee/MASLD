"""U8a: recompute pooled-cell donor program scores from the stored per-cell Hotspot scores.

Own run->donor map built from the donor_pairing files; runs without a pairing entry stay
their own donor (GSE174748, GSE189600 are one run per listed donor in the 2026-10-14
metadata).  Compares to donor_program_scores_primary.tsv (v3 candidate) for every row.
Also writes cells per donor and lineage.
"""
import sys
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq


ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
CELL = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules/cell_scores_all.parquet"
PRIM = ROOT / ("Analysis/Multimodal_Program_Projection/candidates/"
               "program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv")
META = ROOT / ("Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/"
               "metadata/donor_metadata_stage_corrected_donor.tsv")
PAIR = {d: ROOT / f"data/{d}/metadata/donor_pairing.csv" for d in ["GSE244832", "GSE185477", "GSE202379", "GSE136103"]}

run2donor = {}
for ds, fp in PAIR.items():
    dp = pd.read_csv(fp, dtype=str)
    for _, r in dp.iterrows():
        for srr in str(r["rna_srrs"]).split(";"):
            srr = srr.strip()
            if srr:
                assert srr not in run2donor, srr
                run2donor[srr] = f"{ds}_{r['donor_id']}"
print("pairing runs", len(run2donor), "donors", len(set(run2donor.values())))

tab = pq.read_table(CELL, columns=["cell_id", "module", "score", "cell_type"])
print("cell score rows", tab.num_rows)
df = pd.DataFrame({
    "cell_id": tab.column("cell_id").to_pandas().astype("category"),
    "module": tab.column("module").to_pandas().astype("category"),
    "score": tab.column("score").to_numpy(),
    "cell_type": tab.column("cell_type").to_pandas().astype("category"),
})
del tab
# sample = token before the first "_"; check against the token before the last "_"
cats = pd.Series(df["cell_id"].cat.categories)
samp_first = cats.str.split("_", n=1).str[0]
samp_last = cats.str.rsplit("_", n=1).str[0]
n_disagree = int((samp_first != samp_last).sum())
print("cell ids", len(cats), "first-vs-last-underscore sample disagreements", n_disagree)
if n_disagree:
    print(cats[samp_first != samp_last].head(10).tolist())
df["sample"] = pd.Categorical(samp_first.values[df["cell_id"].cat.codes.values])
df["sample"] = df["sample"].astype(str)
df["donor"] = df["sample"].map(run2donor).fillna(df["sample"])
df["program"] = df["module"].astype(str).str.replace("__", "_", regex=False)

meta = pd.read_csv(META, sep="\t", dtype=str)
keep = set(meta["sample"])
print("samples in cell scores", df["sample"].nunique(), "donors", df["donor"].nunique(),
      "donors with metadata", len(set(df["donor"].unique()) & keep))
runs_per = (df.loc[df["cell_type"] == "hepatocytes", ["donor", "sample"]].drop_duplicates()
            .groupby("donor").size())
small = meta.loc[meta["dataset"].isin(["GSE174748", "GSE189600"]), "sample"]
print("runs per donor for GSE174748/GSE189600 donors (hepatocytes):", runs_per.reindex(small).to_dict())
print("GSE189600/GSE174748 runs present:", sorted(s for s in df["sample"].unique() if s in set(small)))

g = df.groupby(["cell_type", "program", "donor"], observed=True)["score"].agg(["mean", "size"]).reset_index()
g.columns = ["cell_type", "program", "donor", "score_recomputed", "n_cells"]
g.to_csv(OUT / "donor_scores_recomputed.tsv.gz", sep="\t", index=False)
cells = g.groupby(["cell_type", "donor"])["n_cells"].first().reset_index()
cells.to_csv(OUT / "cells_per_donor_lineage.tsv", sep="\t", index=False)

prim = pd.read_csv(PRIM, sep="\t", usecols=["cell_type", "module", "donor", "score"])
prim["program"] = prim["cell_type"] + "_" + prim["module"].astype(str)
m = prim.merge(g, on=["cell_type", "program", "donor"], how="left", indicator=True)
print("primary rows", len(prim), "matched", int((m["_merge"] == "both").sum()))
d = (m["score"] - m["score_recomputed"]).abs()
print("max |primary - recomputed| =", float(d.max()), "; rows > 1e-9:", int((d > 1e-9).sum()))
m.to_csv(OUT / "primary_vs_recomputed.tsv.gz", sep="\t", index=False)

# T-cell donors with < 20 T cells among the 53 fitted T-cell donors (GSE202379, GSE244832)
fit_donors = set(meta.loc[(meta["dataset"].isin(["GSE202379", "GSE244832"])) &
                          (meta["exclude_stage_analysis"].str.upper() == "FALSE") &
                          (meta["disease_stage_coarse"].isin(["Healthy", "Steatosis", "Steatohepatitis"])), "sample"])
tc = cells[(cells["cell_type"] == "tcells") & cells["donor"].isin(fit_donors)]
print("fitted T-cell donors", len(tc), "with < 20 T cells", int((tc["n_cells"] < 20).sum()))
for f in [CELL, PRIM, META]:
    h = hashlib.sha256(Path(f).read_bytes()).hexdigest()
    print("sha256", h, f)
