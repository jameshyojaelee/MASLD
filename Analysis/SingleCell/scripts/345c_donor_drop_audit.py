"""
Audit donors dropped by Script 345 (per-donor LIANA).

Read-only consolidation of per-donor LIANA logs + donor metadata + lineage cell
counts. Outputs:
  - donor_drop_audit.tsv          (per-donor, 269 rows)
  - donor_drop_audit_summary.txt  (summary stats)

Gating rules in Script 345:
  - <200 total cells across 9 KEEP_LINEAGES         -> status 'skipped_total_cells'
  - <2 lineages with >=30 cells (MIN_CELLS_PER_LINEAGE) -> status 'skipped_lineages'

The 9 KEEP_LINEAGES are:
  Hepatocytes, Macrophages, Fibroblasts, Endothelial cells, Cholangiocytes,
  T cells, B cells, Resident NK, Plasma cells.
"""
from __future__ import annotations
import glob
import pandas as pd
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ST   = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
LOG_DIR    = ST / "per_donor_lr"
META_TSV   = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
LINEAGE_TSV = ST / "lineage_cell_counts.tsv"

OUT_TSV  = ST / "donor_drop_audit.tsv"
OUT_TXT  = ST / "donor_drop_audit_summary.txt"

MIN_CELLS_PER_LINEAGE = 30
MIN_TOTAL_CELLS       = 200
KEEP_LINEAGES_COLS = [
    "n_Hepatocytes", "n_Macrophages", "n_Fibroblasts", "n_Endothelial_cells",
    "n_Cholangiocytes", "n_T_cells", "n_B_cells", "n_Resident_NK", "n_Plasma_cells",
]

# ----------------- 1. consolidate logs -----------------
log_files = sorted(glob.glob(str(LOG_DIR / "log_*.csv")))
print(f"[audit] found {len(log_files)} log CSVs")
log = pd.concat([pd.read_csv(f) for f in log_files], ignore_index=True)
# de-dupe in case a sample appears in >1 log
log = log.drop_duplicates(subset="sample", keep="last")
print(f"[audit] consolidated unique samples in logs: {len(log)}")

# normalize schema
for col in ("n_cells", "n_lineages", "n_lr", "elapsed_s"):
    if col not in log.columns:
        log[col] = pd.NA

# ----------------- 2. join metadata -----------------
meta = pd.read_csv(META_TSV, sep="\t")
print(f"[audit] donor_metadata rows: {len(meta)}")

# keep canonical 5 fractions (the cols in donor_metadata)
frac_cols = [c for c in meta.columns if c.startswith("frac_")]
meta_keep = meta[["sample", "dataset", "condition", "disease_stage_coarse",
                  "n_cells"] + frac_cols].rename(columns={"n_cells": "meta_n_cells"})

audit = meta_keep.merge(log, on="sample", how="left")
# samples in logs but not in meta (should be 0)
orphan = set(log["sample"]) - set(meta["sample"])
if orphan:
    print(f"[audit] WARNING: {len(orphan)} samples in logs not in donor_metadata: "
          f"{sorted(orphan)[:5]}...")

# samples with no log entry -> status missing
audit["status"] = audit["status"].fillna("no_log_entry")

# ----------------- 3. join lineage counts -----------------
lin = pd.read_csv(LINEAGE_TSV, sep="\t")
keep_lin = lin[["sample"] + [c for c in KEEP_LINEAGES_COLS if c in lin.columns]].copy()

# Per-donor, compute total cells in the 9-lineage subset and which lineages are below floor
keep_lin["n_total_9lin"] = keep_lin[KEEP_LINEAGES_COLS].sum(axis=1)
keep_lin["n_lin_ge30"]   = (keep_lin[KEEP_LINEAGES_COLS] >= MIN_CELLS_PER_LINEAGE).sum(axis=1)
keep_lin["n_lin_ge20"]   = (keep_lin[KEEP_LINEAGES_COLS] >= 20).sum(axis=1)

def lineages_below_floor(row, floor=MIN_CELLS_PER_LINEAGE):
    below = []
    for col in KEEP_LINEAGES_COLS:
        v = row[col]
        if pd.notna(v) and v < floor:
            below.append(f"{col[2:]}({int(v)})")
    return ";".join(below)

def lineages_at_or_above(row, floor=MIN_CELLS_PER_LINEAGE):
    above = []
    for col in KEEP_LINEAGES_COLS:
        v = row[col]
        if pd.notna(v) and v >= floor:
            above.append(f"{col[2:]}({int(v)})")
    return ";".join(above)

keep_lin["lineages_below_30"] = keep_lin.apply(lineages_below_floor, axis=1)
keep_lin["lineages_ge30"]     = keep_lin.apply(lineages_at_or_above, axis=1)

audit = audit.merge(
    keep_lin[["sample", "n_total_9lin", "n_lin_ge30", "n_lin_ge20",
              "lineages_below_30", "lineages_ge30"]],
    on="sample", how="left",
)

# Rescue flag: would relaxing MIN_CELLS_PER_LINEAGE to 20 save this donor?
def rescue_status(row):
    if row["status"] != "skipped_lineages":
        return ""
    if pd.notna(row["n_lin_ge20"]) and row["n_lin_ge20"] >= 2 \
       and pd.notna(row["n_total_9lin"]) and row["n_total_9lin"] >= MIN_TOTAL_CELLS:
        return "rescuable_at_20"
    return "still_skipped_at_20"

audit["rescue_if_floor20"] = audit.apply(rescue_status, axis=1)

# Reorder
front = ["sample", "dataset", "condition", "disease_stage_coarse",
         "status", "n_cells", "n_lineages", "n_lr", "elapsed_s",
         "n_total_9lin", "n_lin_ge30", "n_lin_ge20",
         "lineages_below_30", "lineages_ge30", "rescue_if_floor20"]
audit = audit[front + [c for c in audit.columns if c not in front]]

audit.to_csv(OUT_TSV, sep="\t", index=False)
print(f"[audit] wrote {OUT_TSV} ({len(audit)} rows)")

# ----------------- 4. summary -----------------
lines = []
def w(s=""):
    lines.append(s)
    print(s)

w("=" * 78)
w("Donor drop audit -- Script 345 per-donor LIANA")
w("=" * 78)
w(f"donor_metadata rows           : {len(meta)}")
w(f"unique samples in logs        : {len(log)}")
w(f"samples missing from logs     : {(audit['status']=='no_log_entry').sum()}")
w("")

status_ct = audit["status"].value_counts(dropna=False)
w("Status counts (audit table):")
for s, n in status_ct.items():
    w(f"  {s:30s} {n}")
w("")

dropped = audit[audit["status"].isin(
    ["skipped_total_cells", "skipped_lineages", "no_log_entry"])].copy()
w(f"Total dropped donors: {len(dropped)}")
w("")

# Per dataset
w("Per dataset (processed = ok, dropped = skipped_* or no_log_entry):")
g = audit.groupby("dataset")["status"].apply(
    lambda s: pd.Series({
        "processed_ok": (s == "ok").sum(),
        "skipped_total_cells": (s == "skipped_total_cells").sum(),
        "skipped_lineages":    (s == "skipped_lineages").sum(),
        "no_log_entry":        (s == "no_log_entry").sum(),
        "total": len(s),
    })
).unstack().fillna(0).astype(int)
w(g.to_string())
w("")

# Per disease_stage_coarse
w("Per disease_stage_coarse:")
g2 = audit.groupby("disease_stage_coarse")["status"].apply(
    lambda s: pd.Series({
        "processed_ok": (s == "ok").sum(),
        "skipped_total_cells": (s == "skipped_total_cells").sum(),
        "skipped_lineages":    (s == "skipped_lineages").sum(),
        "no_log_entry":        (s == "no_log_entry").sum(),
        "total": len(s),
    })
).unstack().fillna(0).astype(int)
g2["pct_dropped"] = ((g2["skipped_total_cells"] + g2["skipped_lineages"]
                     + g2["no_log_entry"]) / g2["total"] * 100).round(1)
w(g2.to_string())
w("")

# Per-donor detail for each skipped donor
w("Per-donor detail for each dropped donor:")
w("-" * 78)
disp_cols = ["sample", "dataset", "condition", "disease_stage_coarse",
             "status", "n_cells", "n_total_9lin", "n_lin_ge30", "n_lin_ge20",
             "rescue_if_floor20", "lineages_ge30", "lineages_below_30"]
for _, r in dropped.sort_values(["status", "disease_stage_coarse", "sample"]).iterrows():
    w(f"  sample={r['sample']}  dataset={r['dataset']}  "
      f"stage={r['disease_stage_coarse']}  cond={r['condition']}")
    w(f"    status={r['status']}  log_n_cells={r['n_cells']}  "
      f"n_total_9lin={r['n_total_9lin']}  "
      f"n_lin_ge30={r['n_lin_ge30']}  n_lin_ge20={r['n_lin_ge20']}  "
      f"rescue_if_floor20={r['rescue_if_floor20']}")
    w(f"    lineages_ge30      : {r['lineages_ge30'] or '(none)'}")
    w(f"    lineages_below_30  : {r['lineages_below_30'] or '(none)'}")
w("")

# Rescue summary
resc = (audit["rescue_if_floor20"] == "rescuable_at_20").sum()
w(f"Donors rescuable at MIN_CELLS_PER_LINEAGE=20 (still need >=200 total cells): {resc}")

# Stage-bias chi-square-ish quick view (just print proportions)
w("")
w("Drop-rate by stage (headline bias check):")
for stg, row in g2.iterrows():
    dropped_n = int(row["skipped_total_cells"] + row["skipped_lineages"] + row["no_log_entry"])
    tot = int(row["total"])
    w(f"  {str(stg):24s}  {dropped_n:>3d}/{tot:<3d}  ({row['pct_dropped']}%)")

with open(OUT_TXT, "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"[audit] wrote {OUT_TXT}")
