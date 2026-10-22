"""AN7 step 1: donor x program Hotspot scores from stored per-cell scores, all cells and v2 QC-passing cells.

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

# ---- AN7 additions: write donor score tables (score = mean of stored per-cell scores; n = cells averaged) ----
assert np.nanmax(np.abs(chk.score_pub - chk.score_regen)) < 1e-10, "all-cell regeneration does not match published"
cols = ["cell_type", "module", "donor", "score", "n"]
allc[cols].rename(columns={"n": "n_cells"}).to_csv(OUT / "donor_scores_all_cells.tsv", sep="\t", index=False, float_format="%.17g")
qcc[cols].rename(columns={"n": "n_cells"}).to_csv(OUT / "donor_scores_v2qc_cells.tsv", sep="\t", index=False, float_format="%.17g")
print(f"[write] all-cell rows {len(allc)}, QC-cell rows {len(qcc)}")
