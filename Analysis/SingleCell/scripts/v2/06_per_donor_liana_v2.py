#!/usr/bin/env python
"""
06_per_donor_liana_v2.py

Phase 0.5 v2 equivalent of Script 345.

Per-donor LIANA rank_aggregate on the v2 atlas (all cell types, not just
hepatocytes). LIANA needs hepatocyte + non-hepatocyte cells in the same
AnnData object so source/target lineage pairs are tested. We read Agent 2's
full v2 atlas:
  Analysis/SingleCell/results_gpu_v2_phase05/atlas/scalesc_human_annotated_celltypist_v2.h5ad

A pure hepatocyte atlas (from Agent 3) is NOT sufficient for CCC.

PSEUDOREPLICATION FIX (2026-07-12) — DONOR COLLAPSE
  obs["sample"] in the v2 atlas is a SEQUENCING RUN (SRR / GSM), not a
  biological donor. Five* datasets carry multiple runs per donor; the four
  that survive LIANA gating here are:
    GSE244832  117 runs ->  18 donors  (SRR-keyed; 82 present)
    GSE202379   67 runs ->  46 donors  (SRR-keyed; 60 present)
    GSE185477   21 runs ->   3 donors  (SRR-keyed; 21 present)
    GSE136103   20 GSMs ->  10 donors  (GSM-keyed; CD45 sort fractions)
  Previously this script ran LIANA once per RUN, so those datasets contributed
  2-8x pseudo-independent points to every downstream per-donor statistic (the
  stage LMM, the Fig3H L-R trajectory heatmap). We now POOL the cells across
  all runs of a true donor into one AnnData subset and run LIANA ONCE per true
  donor -- the single-cell analog of the raw-count pooling in the bulk
  lib_donor_collapse: for single cells we take the UNION of the donor's cells
  (each cell is its own observation) rather than summing counts. All other
  datasets are 1 run = 1 donor and pass through unchanged. Donor IDs are
  dataset-prefixed (e.g. "GSE202379_P3") via data/{GSE}/metadata/
  donor_pairing.csv. This mirrors the v1 fix in 345_per_donor_liana.py.

  Outputs go to a NEW dir (per_donor_lr_v2_dc) so the run-level "before"
  artifacts (per_donor_lr_v2/, all_donor_lr_scores_v2.tsv.gz) are preserved
  untouched for before/after comparison. The team lead promotes the _dc
  outputs after review.

Outputs:
  Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2/per_donor_lr_v2_dc/{donor}_lr_scores.parquet
  + a lineage_cell_counts_v2_dc.tsv sidecar (donor x lineage) written once
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

PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)

# v2 inputs
H5AD_V2_FULL = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/scalesc_human_annotated_celltypist_v2.h5ad"
)
H5AD_V2_HEP = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2_annotated.h5ad"
)

# v1 fallback (for donor roster only; never used for cells)
DONOR_META_V1 = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
)
DONOR_META_V2 = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv"
)

V2_STAGE_DIR = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2"
)
# Donor-collapsed outputs (do NOT clobber the run-level per_donor_lr_v2/).
OUT_DIR = Path(
    os.environ.get("PER_DONOR_LR_DIR", str(V2_STAGE_DIR / "per_donor_lr_v2_dc"))
)
COUNTS_TSV = V2_STAGE_DIR / "lineage_cell_counts_v2_dc.tsv"

# Datasets whose obs["sample"] is a sequencing RUN, not a biological donor.
DONOR_PAIRING_FILES = {
    "GSE244832": PROJECT_ROOT / "data/GSE244832/metadata/donor_pairing.csv",
    "GSE202379": PROJECT_ROOT / "data/GSE202379/metadata/donor_pairing.csv",
    "GSE185477": PROJECT_ROOT / "data/GSE185477/metadata/donor_pairing.csv",
    "GSE136103": PROJECT_ROOT / "data/GSE136103/metadata/donor_pairing.csv",
}

KEEP_LINEAGES = {
    "Hepatocytes",
    "Macrophages",
    "Fibroblasts",
    "Endothelial cells",
    "Cholangiocytes",
    "T cells",
    "B cells",
    "Resident NK",
    "Plasma cells",
}
MIN_CELLS_PER_LINEAGE = 30
MIN_TOTAL_CELLS = 200
EXPR_PROP = 0.10
N_PERMS = 100


def resolve_atlas() -> Path:
    """Pick the full v2 atlas; halt if absent."""
    if H5AD_V2_FULL.exists():
        return H5AD_V2_FULL
    sys.exit(
        f"FATAL: v2 full atlas not found at {H5AD_V2_FULL}. "
        "LIANA requires all cell types; halting."
    )


def build_srr_to_donor() -> dict:
    """Map each run ID (SRR or GSM) -> dataset-prefixed biological donor ID.

    Reads data/{GSE}/metadata/donor_pairing.csv (column `rna_srrs` is a
    semicolon-joined list of the runs -- SRR or GSM -- belonging to one
    biological donor). Runs not covered pass through as their own donor
    upstream (see build_donor_roster).
    """
    m = {}
    for ds, fp in DONOR_PAIRING_FILES.items():
        if not fp.exists():
            print(
                f"[srr-map] WARNING donor_pairing.csv missing for {ds}: {fp}",
                file=sys.stderr,
            )
            continue
        dp = pd.read_csv(fp, dtype=str)
        for _, r in dp.iterrows():
            donor_id = f"{ds}_{r['donor_id']}"
            for srr in str(r["rna_srrs"]).split(";"):
                srr = srr.strip()
                if srr:
                    m[srr] = donor_id
    return m


def get_run_roster() -> pd.DataFrame:
    """Prefer v2 donor metadata; fall back to v1 if not yet built.

    Returns the RUN-level roster (obs["sample"] == a sequencing run).
    """
    if DONOR_META_V2.exists():
        roster = pd.read_csv(DONOR_META_V2, sep="\t", dtype={"sample": str})
        print(f"[roster] using v2 donor metadata: {len(roster)} run-level samples")
    else:
        roster = pd.read_csv(DONOR_META_V1, sep="\t", dtype={"sample": str})
        print(
            f"[roster] v2 missing; falling back to v1: {len(roster)} run-level samples"
        )
    return roster


def build_donor_roster() -> list:
    """Return an ordered list of (donor_id, dataset, [member_samples]).

    The analysis universe is the run-level donor_metadata (obs["sample"]).
    Runs belonging to the same biological donor (per donor_pairing.csv) are
    grouped; all other samples become their own donor. Order is deterministic
    (sorted by donor_id) so SLURM array slicing is stable across tasks.
    """
    dm = get_run_roster()
    srr2donor = build_srr_to_donor()
    dm["biological_donor"] = dm["sample"].map(lambda s: srr2donor.get(s, s))
    roster = []
    for donor_id, grp in dm.groupby("biological_donor"):
        member_samples = grp["sample"].tolist()
        dataset = str(grp["dataset"].iloc[0]) if "dataset" in grp.columns else "NA"
        roster.append((donor_id, dataset, member_samples))
    roster.sort(key=lambda x: x[0])
    return roster


def run_donor(
    adata_full: ad.AnnData, donor_id: str, dataset: str, member_samples: list
) -> dict:
    """Run LIANA once for one TRUE donor, pooling cells across its runs."""
    t0 = time.time()
    member_set = set(member_samples)
    sub = adata_full[adata_full.obs["sample"].isin(member_set)].to_memory()
    sample = donor_id  # output key = biological donor
    if "cell_type" not in sub.obs.columns:
        return {
            "sample": sample,
            "status": "missing_cell_type_col",
            "n_runs": len(member_samples),
            "n_cells": int(sub.n_obs),
            "elapsed_s": time.time() - t0,
        }
    sub = sub[sub.obs["cell_type"].isin(KEEP_LINEAGES)].copy()
    if sub.n_obs < MIN_TOTAL_CELLS:
        return {
            "sample": sample,
            "status": "skipped_total_cells",
            "n_runs": len(member_samples),
            "n_cells": int(sub.n_obs),
            "elapsed_s": time.time() - t0,
        }
    counts = sub.obs["cell_type"].value_counts()
    keep_cts = counts[counts >= MIN_CELLS_PER_LINEAGE].index.tolist()
    if len(keep_cts) < 2:
        return {
            "sample": sample,
            "status": "skipped_lineages",
            "n_runs": len(member_samples),
            "n_cells": int(sub.n_obs),
            "n_lineages": len(keep_cts),
            "elapsed_s": time.time() - t0,
        }
    sub = sub[sub.obs["cell_type"].isin(keep_cts)].copy()
    sub.obs["cell_type"] = sub.obs["cell_type"].astype(str)

    import liana as li

    try:
        li.mt.rank_aggregate(
            sub,
            groupby="cell_type",
            resource_name="consensus",
            expr_prop=EXPR_PROP,
            min_cells=5,
            n_perms=N_PERMS,
            use_raw=False,
            verbose=False,
            seed=42,
        )
    except Exception as e:
        return {
            "sample": sample,
            "status": f"liana_failed:{type(e).__name__}:{e}",
            "n_runs": len(member_samples),
            "n_cells": int(sub.n_obs),
            "elapsed_s": time.time() - t0,
        }

    res = sub.uns.get("liana_res")
    if res is None or len(res) == 0:
        return {
            "sample": sample,
            "status": "empty_result",
            "n_runs": len(member_samples),
            "n_cells": int(sub.n_obs),
            "elapsed_s": time.time() - t0,
        }
    res = res.copy()
    res["sample"] = sample  # biological donor ID (join key)
    res["biological_donor"] = donor_id
    res["dataset"] = dataset
    ct_counts = counts.to_dict()
    res["n_source_cells"] = res["source"].map(ct_counts)
    res["n_target_cells"] = res["target"].map(ct_counts)
    keep_cols = [
        "sample",
        "biological_donor",
        "dataset",
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
        "magnitude_rank",
        "specificity_rank",
        "n_source_cells",
        "n_target_cells",
    ]
    res = res[[c for c in keep_cols if c in res.columns]]
    out_parquet = OUT_DIR / f"{sample}_lr_scores.parquet"
    res.to_parquet(out_parquet, index=False)
    return {
        "sample": sample,
        "status": "ok",
        "n_runs": len(member_samples),
        "n_cells": int(sub.n_obs),
        "n_lineages": len(keep_cts),
        "n_lr": int(len(res)),
        "out": str(out_parquet),
        "elapsed_s": round(time.time() - t0, 1),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor")
    ap.add_argument("--donor-list")
    ap.add_argument("--array-task-id", type=int)
    ap.add_argument("--array-batch-size", type=int, default=8)
    ap.add_argument("--write-counts", action="store_true")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    V2_STAGE_DIR.mkdir(parents=True, exist_ok=True)

    h5ad_path = resolve_atlas()
    print(f"[atlas] {h5ad_path}")
    adata_full = ad.read_h5ad(h5ad_path, backed="r")
    print(
        f"[atlas] shape={adata_full.shape}; "
        f"obs cols include cell_type={'cell_type' in adata_full.obs.columns}"
    )

    roster = build_donor_roster()
    roster_by_id = {
        donor_id: (dataset, members) for donor_id, dataset, members in roster
    }
    all_donors = [d for d, _, _ in roster]
    n_multi = sum(1 for _, _, m in roster if len(m) > 1)
    print(
        f"[roster] {len(roster)} TRUE donors "
        f"(from {sum(len(m) for _, _, m in roster)} run-level samples); "
        f"{n_multi} donors pool >1 run"
    )

    if args.donor:
        donors = [args.donor]
    elif args.donor_list:
        donors = [L.strip() for L in open(args.donor_list) if L.strip()]
    elif args.array_task_id is not None:
        i = args.array_task_id - 1
        donors = all_donors[
            i * args.array_batch_size : (i + 1) * args.array_batch_size
        ]
    else:
        donors = all_donors

    if not donors:
        print("[main] no donors to process; exiting", file=sys.stderr)
        sys.exit(0)
    print(f"[main] processing {len(donors)} donors")

    if args.write_counts and "cell_type" in adata_full.obs.columns:
        srr2donor = build_srr_to_donor()
        obs = adata_full.obs[["sample", "cell_type"]].copy()
        obs["sample"] = obs["sample"].astype(str)
        obs["biological_donor"] = obs["sample"].map(lambda s: srr2donor.get(s, s))
        obs = obs[obs["cell_type"].isin(KEEP_LINEAGES)]
        wide = (
            obs.groupby(["biological_donor", "cell_type"], observed=True)
            .size()
            .unstack(fill_value=0)
            .reset_index()
            .rename(columns={"biological_donor": "sample"})
        )
        wide.columns = [
            f"n_{c.replace(' ', '_')}" if c != "sample" else c
            for c in wide.columns
        ]
        wide.to_csv(COUNTS_TSV, sep="\t", index=False)
        print(f"[counts] wrote {COUNTS_TSV} ({wide.shape})")

    log_rows = []
    for d in donors:
        dataset, members = roster_by_id.get(d, ("NA", [d]))
        try:
            row = run_donor(adata_full, d, dataset, members)
        except Exception as e:
            row = {"sample": d, "status": f"unhandled:{type(e).__name__}:{e}"}
        print(
            f"[donor] {d}: {row.get('status')} "
            f"({row.get('n_runs', '?')} runs, "
            f"{row.get('n_cells', '?')} cells, "
            f"{row.get('n_lr', '?')} LR, {row.get('elapsed_s', '?')} s)",
            flush=True,
        )
        log_rows.append(row)
        gc.collect()

    log = pd.DataFrame(log_rows)
    log_csv = OUT_DIR / f"log_{donors[0]}_{donors[-1]}.csv"
    log.to_csv(log_csv, index=False)
    print(f"[main] log written to {log_csv}")


if __name__ == "__main__":
    main()
