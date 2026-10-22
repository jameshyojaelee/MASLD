"""F4abcd G1-03: cells per donor per scored lineage behind the 117-program donor scores.

Counts atlas cells per (lineage, donor) from cell_scores_all.parquet (module 1 of each lineage;
counts are identical across modules within a lineage), maps runs to donors with the
511_pooled_cell_donor_scores.py logic (as copied in review job 1711), and summarizes donors
below 20 and 50 cells in (a) the 62-donor OLS roster and (b) the H59 fitted universe
(GSE202379 and GSE244832 donors of the 60-donor universe).
"""
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1])
PARQ = BASE / "Analysis/SingleCell/results_gpu_v2/hotspot_modules/cell_scores_all.parquet"
META = BASE / "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
SCORES = BASE / "figures/candidates/source-diagnosis-programs-20261014T161500Z_22237556/hotspot/donor_program_scores_primary.tsv"
STAGES = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2}
DONOR_PAIRING = {"GSE244832": "data/GSE244832/metadata/donor_pairing.csv",
                 "GSE185477": "data/GSE185477/metadata/donor_pairing.csv",
                 "GSE202379": "data/GSE202379/metadata/donor_pairing.csv",
                 "GSE136103": "data/GSE136103/metadata/donor_pairing.csv"}


def build_sample_to_donor():
    m = {}
    for ds, rel in DONOR_PAIRING.items():
        dp = pd.read_csv(BASE / rel, dtype=str)
        for _, row in dp.iterrows():
            for srr in str(row["rna_srrs"]).split(";"):
                srr = srr.strip()
                if srr:
                    m[srr] = f"{ds}_{row['donor_id']}"
    runs = pd.read_csv(BASE / "Liver_Atlas/metadata/SraRunTable.csv", dtype=str, usecols=["Run", "shortfilename"])
    samples = pd.read_csv(BASE / "Liver_Atlas/metadata/GSE192740_sampleInfo_scRNAseq.tsv", sep="\t", dtype=str,
                          usecols=["characteristics: shortFileName", "title"]).rename(
        columns={"characteristics: shortFileName": "shortfilename"})
    liver = runs.merge(samples, on="shortfilename", how="inner")
    liver["donor_short"] = liver["title"].map(lambda v: (mm.group(0) if (mm := re.search(r"H[0-9]+", v)) else None))
    liver = liver.dropna(subset=["donor_short"])
    m.update(dict(zip(liver["Run"], "Liver_Atlas_" + liver["donor_short"])))
    return m


s2d = build_sample_to_donor()
pf = pq.ParquetFile(PARQ)
parts = []
for rg in range(pf.metadata.num_row_groups):
    t = pf.read_row_group(rg, columns=["cell_id", "module", "cell_type"]).to_pandas()
    t = t[t["module"].str.endswith("__1")]
    t["sample"] = t["cell_id"].str.split("_", n=1).str[0]
    t["donor"] = t["sample"].map(s2d).fillna(t["sample"])
    parts.append(t.groupby(["cell_type", "donor"]).size())
    del t
n = pd.concat(parts).groupby(level=[0, 1]).sum().rename("n_cells").reset_index()

# donors that carry a score for each lineage in the published 62-donor primary table
sc = pd.read_csv(SCORES, sep="\t", usecols=["cell_type", "module", "donor", "dataset", "stage_ordinal",
                                            "disease_stage_coarse", "exclude_stage_analysis"])
sc = sc[(sc.exclude_stage_analysis.astype(str).str.upper() == "FALSE") & sc.disease_stage_coarse.isin(STAGES)]
roster = sc.drop_duplicates(["cell_type", "donor"])[["cell_type", "donor", "dataset", "stage_ordinal"]]
tab = roster.merge(n, on=["cell_type", "donor"], how="left")
tab["n_cells"] = tab["n_cells"].fillna(0).astype(int)
tab["in_h59_fit"] = tab.dataset.isin(["GSE202379", "GSE244832"])
tab.sort_values(["cell_type", "n_cells"]).to_csv(OUT / "cells_per_donor_lineage_62roster.tsv", sep="\t", index=False)
rows = []
for label, sub in [("ols62_roster", tab), ("h59_fit_GSE202379_GSE244832", tab[tab.in_h59_fit])]:
    for ct, g in sub.groupby("cell_type"):
        rows.append(dict(universe=label, cell_type=ct, donors=len(g), lt20=int((g.n_cells < 20).sum()),
                         lt50=int((g.n_cells < 50).sum()), min_cells=int(g.n_cells.min()),
                         median_cells=float(g.n_cells.median())))
S = pd.DataFrame(rows)
S.to_csv(OUT / "cells_per_donor_lineage_summary.tsv", sep="\t", index=False)
print(S.to_string(index=False))
print(tab[tab.n_cells < 5].sort_values(["cell_type", "n_cells"]).to_string(index=False))
