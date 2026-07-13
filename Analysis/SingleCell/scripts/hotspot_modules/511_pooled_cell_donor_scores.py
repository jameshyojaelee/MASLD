"""Pooled-cell (cell-count-weighted) donor-level Hotspot module scores.

SENSITIVITY arm for the run-vs-donor pseudoreplication fix. The canonical fix
(511_donorlevel_disease_stage.R) collapses run-level module scores to one per
biological donor by the UNWEIGHTED mean of a donor's per-run scores. This script
produces the alternative POOLED-CELL donor score: the mean module score over ALL
of a donor's cells (equivalent to treating the donor as one sample from the
start), which is the exact donor-level statistic. Comparing the two confirms the
weighting choice does not change the disease-stage conclusions.

Module DEFINITIONS are NOT recomputed - per-cell scores in cell_scores_all.parquet
are reused verbatim; only the cell -> donor aggregation changes (sample -> donor).

Output (NEW path, does not touch canonical donor_scores_all.tsv):
  results_gpu_v2/hotspot_modules/donor_collapse/donor_scores_all_weighted.tsv
  columns: sample (= dataset-prefixed donor id), module, score, cell_type
"""
from __future__ import annotations
import os
from pathlib import Path
import pandas as pd

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
RES = BASE / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
OUT = RES / "donor_collapse"
OUT.mkdir(parents=True, exist_ok=True)

# Mirror lib_donor_collapse.R::build_srr_to_donor_map in python. Keys = raw atlas
# sample ids (SRR run id, or GSM for GSE136103); values = "{dataset}_{donor_id}".
DONOR_PAIRING = {
    "GSE244832": "data/GSE244832/metadata/donor_pairing.csv",
    "GSE185477": "data/GSE185477/metadata/donor_pairing.csv",
    "GSE202379": "data/GSE202379/metadata/donor_pairing.csv",
    "GSE136103": "data/GSE136103/metadata/donor_pairing.csv",
}


def build_sample_to_donor() -> dict[str, str]:
    m: dict[str, str] = {}
    for ds, rel in DONOR_PAIRING.items():
        fp = BASE / rel
        dp = pd.read_csv(fp, dtype=str)
        for _, row in dp.iterrows():
            donor = f"{ds}_{row['donor_id']}"
            for srr in str(row["rna_srrs"]).split(";"):
                srr = srr.strip()
                if srr:
                    m[srr] = donor
    return m


def main() -> None:
    s2d = build_sample_to_donor()
    print(f"[511py] sample->donor map: {len(s2d)} run ids")

    cs = pd.read_parquet(RES / "cell_scores_all.parquet",
                         columns=["cell_id", "module", "score", "cell_type"])
    print(f"[511py] cell_scores_all: {len(cs):,} rows")

    # cell_id == "{sample}_{barcode}"; the leading token before the first "_" is the
    # atlas sample id (validated: the derived set is an exact bijection with the 223
    # donor_scores_all.tsv samples).
    cs["sample"] = cs["cell_id"].str.split("_", n=1).str[0]
    cs["donor"] = cs["sample"].map(s2d).fillna(cs["sample"])

    n_runs = cs["sample"].nunique()
    n_donors = cs["donor"].nunique()
    print(f"[511py] cells map to {n_runs} runs -> {n_donors} true donors")

    # Pooled-cell donor score = mean over ALL of a donor's cells (cell-count weighted).
    donor = (cs.groupby(["cell_type", "module", "donor"], as_index=False)["score"]
               .mean()
               .rename(columns={"donor": "sample"}))
    donor.to_csv(OUT / "donor_scores_all_weighted.tsv", sep="\t", index=False)
    print(f"[511py] wrote {OUT/'donor_scores_all_weighted.tsv'} "
          f"({len(donor):,} rows, {donor['sample'].nunique()} donors)")


if __name__ == "__main__":
    main()
