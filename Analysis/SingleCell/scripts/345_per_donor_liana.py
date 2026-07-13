#!/usr/bin/env python
"""
345_per_donor_liana.py

Per-TRUE-DONOR LIANA rank_aggregate on the 9-lineage subset of the MASLD scRNA
atlas. Runs one biological donor (or a batch of donors) per SLURM array task.
Outputs a tidy parquet per donor.

PSEUDOREPLICATION FIX (2026-07-12):
  obs["sample"] in the integrated atlas is a SEQUENCING RUN (SRR/GSM), not a
  biological donor. Four datasets have multiple runs per donor:
    GSE244832  117 runs -> 18 donors
    GSE185477   21 runs ->  3 donors
    GSE202379   67 runs -> 46 donors
    GSE136103   20 GSMs -> 10 donors  (CD45-sort-fraction multiplicity, GSM-keyed)
  Previously this script ran LIANA once per RUN, so those datasets contributed
  2-8x pseudo-independent points downstream. We now POOL the cells across all
  runs of a true donor into one AnnData subset and run LIANA ONCE per true
  donor (the single-cell analog of the R utility's raw-count pooling in
  lib_donor_collapse.R -- for single cells we take the UNION of the donor's
  cells rather than summing counts, since each cell is its own observation).
  All other datasets are 1 run = 1 donor and pass through unchanged.
  Donor IDs are dataset-prefixed (e.g. "GSE202379_P3") via
  data/{GSE}/metadata/donor_pairing.csv.

9 lineages (CellTypist labels):
  Hepatocytes, Macrophages, Fibroblasts, Endothelial cells, Cholangiocytes,
  T cells, B cells, Resident NK, Plasma cells

Per-donor gating (applied AFTER pooling a donor's runs):
  - drop lineages with < MIN_CELLS_PER_LINEAGE (default 30)
  - skip donor if < 2 surviving lineages
  - skip donor if total cells < MIN_TOTAL_CELLS (default 200)

Output columns (per parquet):
  sample (= biological donor ID), biological_donor, dataset,
  source, target, ligand_complex, receptor_complex,
  magnitude_rank, specificity_rank, n_source_cells, n_target_cells

Usage:
  python 345_per_donor_liana.py --donor GSE202379_P3
  python 345_per_donor_liana.py --donor-list /path/to/donor_batch.txt
  python 345_per_donor_liana.py --array-task-id 3 --array-batch-size 4

SLURM array launch (see run_345_346_donorcollapse.sbatch).
"""

from __future__ import annotations
import argparse
import gc
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

H5AD_PATH  = PROJECT_ROOT / "Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
DONOR_META = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
# NEW donor-collapsed output locations (do NOT clobber the run-level per_donor_lr/).
OUT_DIR    = Path(os.environ.get("PER_DONOR_LR_DIR", str(
    PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/per_donor_lr_donorcollapsed")))
COUNTS_TSV = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/lineage_cell_counts_donorcollapsed.tsv"

# Datasets whose obs["sample"] is a sequencing RUN (or, for GSE136103, a
# CD45-sort-fraction GSM), not a biological donor.
DONOR_PAIRING_FILES = {
    "GSE244832": PROJECT_ROOT / "data/GSE244832/metadata/donor_pairing.csv",
    "GSE185477": PROJECT_ROOT / "data/GSE185477/metadata/donor_pairing.csv",
    "GSE202379": PROJECT_ROOT / "data/GSE202379/metadata/donor_pairing.csv",
    "GSE136103": PROJECT_ROOT / "data/GSE136103/metadata/donor_pairing.csv",
}

OUT_DIR.mkdir(parents=True, exist_ok=True)

KEEP_LINEAGES = {
    "Hepatocytes", "Macrophages", "Fibroblasts", "Endothelial cells",
    "Cholangiocytes", "T cells", "B cells", "Resident NK", "Plasma cells",
}
MIN_CELLS_PER_LINEAGE = 30
MIN_TOTAL_CELLS       = 200
EXPR_PROP             = 0.10
N_PERMS               = 100  # LIANA permutations; 100 fast, 1000 slow


def build_srr_to_donor() -> dict:
    """Map each SRR run ID -> dataset-prefixed biological donor ID.

    Reads data/{GSE}/metadata/donor_pairing.csv (column `rna_srrs` is a
    semicolon-joined list of the runs belonging to one biological donor).
    Runs not covered pass through as their own donor upstream (see
    build_donor_roster).
    """
    m = {}
    for ds, fp in DONOR_PAIRING_FILES.items():
        if not fp.exists():
            print(f"[srr-map] WARNING donor_pairing.csv missing for {ds}: {fp}",
                  file=sys.stderr)
            continue
        dp = pd.read_csv(fp, dtype=str)
        for _, r in dp.iterrows():
            donor_id = f"{ds}_{r['donor_id']}"
            for srr in str(r["rna_srrs"]).split(";"):
                srr = srr.strip()
                if srr:
                    m[srr] = donor_id
    return m


def build_donor_roster() -> list:
    """Return an ordered list of (donor_id, dataset, [member_samples]).

    The analysis universe is donor_metadata.tsv (run-level `sample`). Runs
    belonging to the same biological donor (per donor_pairing.csv) are grouped;
    all other samples become their own donor. Order is deterministic (sorted by
    donor_id) so SLURM array slicing is stable across tasks.
    """
    dm = pd.read_csv(DONOR_META, sep="\t", dtype={"sample": str})
    srr2donor = build_srr_to_donor()
    dm["biological_donor"] = dm["sample"].map(lambda s: srr2donor.get(s, s))
    # dataset is invariant within a donor; take the first
    roster = []
    for donor_id, grp in dm.groupby("biological_donor"):
        member_samples = grp["sample"].tolist()
        dataset = str(grp["dataset"].iloc[0])
        roster.append((donor_id, dataset, member_samples))
    roster.sort(key=lambda x: x[0])
    return roster


def run_donor(adata_full: ad.AnnData, donor_id: str, dataset: str,
              member_samples: list) -> dict:
    """Run LIANA once for one TRUE donor, pooling cells across its runs."""
    t0 = time.time()
    member_set = set(member_samples)
    sub = adata_full[adata_full.obs["sample"].isin(member_set)].to_memory()
    sub = sub[sub.obs["cell_type"].isin(KEEP_LINEAGES)].copy()
    sample = donor_id  # output key = biological donor

    if sub.n_obs < MIN_TOTAL_CELLS:
        return {"sample": sample, "status": "skipped_total_cells",
                "n_runs": len(member_samples),
                "n_cells": int(sub.n_obs), "elapsed_s": time.time() - t0}

    # Cell-type counts; drop lineages with too few cells (set their cells aside)
    counts = sub.obs["cell_type"].value_counts()
    keep_cts = counts[counts >= MIN_CELLS_PER_LINEAGE].index.tolist()
    if len(keep_cts) < 2:
        return {"sample": sample, "status": "skipped_lineages",
                "n_cells": int(sub.n_obs),
                "n_lineages": len(keep_cts), "elapsed_s": time.time() - t0}
    sub = sub[sub.obs["cell_type"].isin(keep_cts)].copy()
    sub.obs["cell_type"] = sub.obs["cell_type"].astype(str)

    # Run LIANA on this donor. Use logged + raw counts; depends on what's in .X
    import liana as li
    try:
        li.mt.rank_aggregate(
            sub,
            groupby="cell_type",
            resource_name="consensus",
            expr_prop=EXPR_PROP,
            min_cells=5,
            n_perms=N_PERMS,
            use_raw=False,   # .X is normalized log1p (CellTypist-annotated h5ad)
            verbose=False,
            seed=42,
        )
    except Exception as e:
        return {"sample": sample, "status": f"liana_failed:{type(e).__name__}:{e}",
                "n_cells": int(sub.n_obs), "elapsed_s": time.time() - t0}

    res = sub.uns.get("liana_res")
    if res is None or len(res) == 0:
        return {"sample": sample, "status": "empty_result",
                "n_cells": int(sub.n_obs), "elapsed_s": time.time() - t0}

    # Annotate with donor + cell counts per cell-type
    res = res.copy()
    res["sample"] = sample                 # biological donor ID (join key)
    res["biological_donor"] = donor_id
    res["dataset"] = dataset
    ct_counts = counts.to_dict()           # pooled cell counts per lineage
    res["n_source_cells"] = res["source"].map(ct_counts)
    res["n_target_cells"] = res["target"].map(ct_counts)
    keep_cols = ["sample", "biological_donor", "dataset",
                 "source", "target", "ligand_complex", "receptor_complex",
                 "magnitude_rank", "specificity_rank",
                 "n_source_cells", "n_target_cells"]
    res = res[[c for c in keep_cols if c in res.columns]]

    out_parquet = OUT_DIR / f"{sample}_lr_scores.parquet"
    res.to_parquet(out_parquet, index=False)

    return {"sample": sample, "status": "ok",
            "n_runs": len(member_samples),
            "n_cells": int(sub.n_obs),
            "n_lineages": len(keep_cts),
            "n_lr": int(len(res)),
            "out": str(out_parquet),
            "elapsed_s": round(time.time() - t0, 1)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", help="single donor (sample) ID")
    ap.add_argument("--donor-list", help="file with one donor per line")
    ap.add_argument("--array-task-id", type=int,
                    help="SLURM array task id (1-based)")
    ap.add_argument("--array-batch-size", type=int, default=4,
                    help="donors per array task")
    ap.add_argument("--write-counts", action="store_true",
                    help="also emit lineage_cell_counts sidecar (donor x lineage)")
    args = ap.parse_args()

    # Resolve TRUE-donor roster: [(donor_id, dataset, [member_samples]), ...]
    roster = build_donor_roster()
    roster_by_id = {donor_id: (dataset, members)
                    for donor_id, dataset, members in roster}
    print(f"[main] donor roster: {len(roster)} true donors "
          f"(from {sum(len(m) for _, _, m in roster)} run-level samples)")

    if args.donor:
        donor_ids = [args.donor]
    elif args.donor_list:
        donor_ids = [L.strip() for L in open(args.donor_list) if L.strip()]
    elif args.array_task_id is not None:
        i = args.array_task_id - 1
        donor_ids = [d for d, _, _ in
                     roster[i * args.array_batch_size:(i + 1) * args.array_batch_size]]
    else:
        donor_ids = [d for d, _, _ in roster]

    if not donor_ids:
        print("[main] no donors to process; exiting", file=sys.stderr)
        sys.exit(0)

    print(f"[main] processing {len(donor_ids)} donors")
    print(f"[main] h5ad: {H5AD_PATH}")
    print(f"[main] out dir: {OUT_DIR}")
    adata_full = ad.read_h5ad(H5AD_PATH, backed="r")
    print(f"[main] full atlas shape: {adata_full.shape}")

    if args.write_counts:
        # Donor x lineage cell counts (pooled across each donor's runs).
        srr2donor = build_srr_to_donor()
        obs = adata_full.obs[["sample", "cell_type"]].copy()
        obs["sample"] = obs["sample"].astype(str)
        obs["biological_donor"] = obs["sample"].map(lambda s: srr2donor.get(s, s))
        obs = obs[obs["cell_type"].isin(KEEP_LINEAGES)]
        wide = (obs.groupby(["biological_donor", "cell_type"]).size()
                  .unstack(fill_value=0).reset_index())
        wide = wide.rename(columns={"biological_donor": "sample"})
        wide.columns = [f"n_{c.replace(' ', '_')}"
                        if c != "sample" else c for c in wide.columns]
        wide.to_csv(COUNTS_TSV, sep="\t", index=False)
        print(f"[counts] wrote {COUNTS_TSV} ({wide.shape})")

    log_rows = []
    for d in donor_ids:
        dataset, members = roster_by_id.get(d, (None, [d]))
        try:
            row = run_donor(adata_full, d, dataset, members)
        except Exception as e:
            row = {"sample": d, "status": f"unhandled:{type(e).__name__}:{e}"}
        print(f"[donor] {d}: {row.get('status')} "
              f"({row.get('n_runs', '?')} runs, "
              f"{row.get('n_cells', '?')} cells, "
              f"{row.get('n_lr', '?')} LR pairs, "
              f"{row.get('elapsed_s', '?')} s)",
              flush=True)
        log_rows.append(row)
        gc.collect()

    log = pd.DataFrame(log_rows)
    log_csv = OUT_DIR / f"log_{donor_ids[0]}_{donor_ids[-1]}.csv"
    log.to_csv(log_csv, index=False)
    print(f"[main] log written to {log_csv}")

if __name__ == "__main__":
    main()
