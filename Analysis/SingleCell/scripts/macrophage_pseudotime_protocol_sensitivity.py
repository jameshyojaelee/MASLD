#!/usr/bin/env python
# =============================================================================
# macrophage_pseudotime_protocol_sensitivity.py
#
# Sensitivity check: recompute per-donor mean macrophage consensus pseudotime
# vs F_stage_augmented_clean Spearman rho, both INCLUDING and EXCLUDING the
# GSE136103 / Liver_Atlas donors flagged by exclude_stage_analysis == TRUE.
#
# Methodology mirrors 344_build_donor_metadata_extended.R lines 164-183:
#   - Filter consensus_pseudotime_all.csv to cell_type == "Macrophages"
#   - Extract sample from cell_id via re.sub(r"_[^_]+$", "", cell_id)
#   - Donor mean of consensus_pseudotime by sample
#   - Spearman vs F_stage_augmented_clean (fallback: F_stage_augmented)
#
# Output: macrophage_pseudotime_protocol_sensitivity.tsv
# =============================================================================
import os
import re
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
PT_FILE = os.path.join(
    BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime/consensus_pseudotime_all.csv"
)
META_FILE = os.path.join(
    BASE,
    "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv",
)
OUT_FILE = os.path.join(
    BASE,
    "Analysis/SingleCell/results_gpu_v2/pseudotime/macrophage_pseudotime_protocol_sensitivity.tsv",
)

print(f"[load] pseudotime: {PT_FILE}")
pt = pd.read_csv(PT_FILE)
# Header: index column is unnamed (cell_id), then cell_type, ...
pt = pt.rename(columns={pt.columns[0]: "cell_id"})
print(f"[load] pseudotime rows: {len(pt):,}  cell_types: {pt['cell_type'].unique().tolist()}")

# Macrophages only
mac = pt.loc[pt["cell_type"] == "Macrophages"].copy()
print(f"[filter] macrophage cells: {len(mac):,}")

# Extract sample id (mirror 344.R: sub('_[^_]+$', '', cell_id))
mac["sample"] = mac["cell_id"].str.replace(r"_[^_]+$", "", regex=True)

# Per-donor mean
donor_mac = (
    mac.groupby("sample", as_index=False)
    .agg(
        macrophage_pseudotime_mean=("consensus_pseudotime", "mean"),
        n_macrophages_with_pt=("consensus_pseudotime", "size"),
    )
)
print(f"[aggregate] donors with macrophage pseudotime: {len(donor_mac):,}")

# Load donor metadata
print(f"[load] donor metadata: {META_FILE}")
meta = pd.read_csv(META_FILE, sep="\t", low_memory=False)
print(f"[load] donor rows: {len(meta):,}; columns include {[c for c in meta.columns if 'F_stage' in c or 'exclude' in c]}")

# Merge
donor = meta.merge(
    donor_mac.rename(columns={
        "macrophage_pseudotime_mean": "macrophage_pseudotime_mean_recomputed",
        "n_macrophages_with_pt": "n_macrophages_with_pt_recomputed",
    }),
    on="sample",
    how="left",
)

# Sanity check that 344.R's stored values match what we recompute
both_present = donor.dropna(subset=["macrophage_pseudotime_mean", "macrophage_pseudotime_mean_recomputed"])
if len(both_present):
    diff = (both_present["macrophage_pseudotime_mean"] - both_present["macrophage_pseudotime_mean_recomputed"]).abs()
    print(f"[sanity] vs 344-stored: n={len(both_present)} max|diff|={diff.max():.2e} mean|diff|={diff.mean():.2e}")

# Choose F-stage axis. NOTE: F_stage_augmented_clean is set to NA for all
# exclude_stage_analysis==TRUE donors by design, which would make the
# full-vs-clean contrast vacuous. Use F_stage_augmented (retains values for
# all donors) and apply the exclusion ourselves to compare matched contrasts.
F_AXIS = "F_stage_augmented"
if F_AXIS not in donor.columns or donor[F_AXIS].notna().sum() < 30:
    F_AXIS = "F_stage_augmented_clean"
print(f"[axis] F-stage column = {F_AXIS}")

# Coerce F_stage to numeric
donor["F_stage_num"] = pd.to_numeric(donor[F_AXIS], errors="coerce")
print(
    f"[axis] non-NA F_stage = {donor['F_stage_num'].notna().sum()} "
    f"(full set, before macrophage merge filter)"
)

# Identify GSE136103 / Liver_Atlas exclusion
EXCL_COL = "exclude_stage_analysis"
if EXCL_COL in donor.columns:
    excl = donor[EXCL_COL].astype(str).str.upper() == "TRUE"
else:
    # Fallback: exclude by dataset name
    excl = donor["dataset"].isin(["GSE136103", "Liver_Atlas"])
print(f"[exclude] flagged donors (exclude_stage_analysis==TRUE): {excl.sum()}")
print(f"[exclude] by dataset: {donor.loc[excl, 'dataset'].value_counts().to_dict()}")

# --- Full set: all donors with mac pseudotime + F-stage ---
full_mask = donor["macrophage_pseudotime_mean_recomputed"].notna() & donor["F_stage_num"].notna()
full = donor.loc[full_mask].copy()
rho_full, p_full = spearmanr(full["macrophage_pseudotime_mean_recomputed"], full["F_stage_num"])
n_full = len(full)

# Per-stage donor counts & per-stage mean for delta
def stage_summary(df):
    out = (
        df.groupby("F_stage_num")
        .agg(n=("sample", "size"), mean_pt=("macrophage_pseudotime_mean_recomputed", "mean"))
        .reset_index()
        .sort_values("F_stage_num")
    )
    return out

full_stage = stage_summary(full)
delta_full = float(full_stage["mean_pt"].max() - full_stage["mean_pt"].min())
print("\n[FULL set, including GSE136103/Liver_Atlas]")
print(full_stage.to_string(index=False))
print(f"  rho = {rho_full:.4f}  p = {p_full:.3e}  n_donors = {n_full}  delta = {delta_full:.4f}")

# --- Clean set: exclude flagged donors ---
clean_mask = full_mask & (~excl)
clean = donor.loc[clean_mask].copy()
rho_clean, p_clean = spearmanr(
    clean["macrophage_pseudotime_mean_recomputed"], clean["F_stage_num"]
)
n_clean = len(clean)
clean_stage = stage_summary(clean)
delta_clean = float(clean_stage["mean_pt"].max() - clean_stage["mean_pt"].min())

print("\n[CLEAN set, excluding GSE136103/Liver_Atlas]")
print(clean_stage.to_string(index=False))
print(f"  rho = {rho_clean:.4f}  p = {p_clean:.3e}  n_donors = {n_clean}  delta = {delta_clean:.4f}")

# --- Excluded subset: are GSE136103 macrophages systematically higher/lower? ---
excl_mask = full_mask & excl
excl_df = donor.loc[excl_mask].copy()
print(f"\n[EXCLUDED donors with pseudotime + F-stage] n={len(excl_df)}")
if len(excl_df) > 0:
    mean_excl = excl_df["macrophage_pseudotime_mean_recomputed"].mean()
    mean_kept = clean["macrophage_pseudotime_mean_recomputed"].mean()
    print(f"  mean pseudotime EXCLUDED = {mean_excl:.4f}")
    print(f"  mean pseudotime KEPT    = {mean_kept:.4f}")
    print(f"  delta (excl - kept)     = {mean_excl - mean_kept:+.4f}")
    print(excl_df.groupby("dataset").agg(
        n=("sample", "size"),
        mean_pt=("macrophage_pseudotime_mean_recomputed", "mean"),
    ).to_string())

# --- Also: cell-level check (whole tissue macrophage pseudotime distribution) ---
# Does GSE136103 macrophage pseudotime distribution differ from the rest?
mac_with_dataset = mac.merge(meta[["sample", "dataset"]].drop_duplicates(), on="sample", how="left")
excl_dataset = mac_with_dataset["dataset"].isin(["GSE136103", "Liver_Atlas"])
print(f"\n[cell-level] GSE136103/Liver_Atlas macrophage cells: {excl_dataset.sum():,} / {len(mac_with_dataset):,}")
print(
    f"  pseudotime mean (excluded datasets) = "
    f"{mac_with_dataset.loc[excl_dataset, 'consensus_pseudotime'].mean():.4f}"
)
print(
    f"  pseudotime mean (kept datasets)     = "
    f"{mac_with_dataset.loc[~excl_dataset, 'consensus_pseudotime'].mean():.4f}"
)

# --- Cross-check: also run on F_stage_augmented_clean (NA for excluded donors)
clean_axis_mask = (
    donor["macrophage_pseudotime_mean_recomputed"].notna()
    & donor.get("F_stage_augmented_clean", pd.Series([np.nan] * len(donor))).notna()
)
clean_axis = donor.loc[clean_axis_mask].copy()
clean_axis["F_stage_clean_num"] = pd.to_numeric(clean_axis["F_stage_augmented_clean"], errors="coerce")
if len(clean_axis) > 3:
    rho_axis_clean, p_axis_clean = spearmanr(
        clean_axis["macrophage_pseudotime_mean_recomputed"], clean_axis["F_stage_clean_num"]
    )
    n_axis_clean = len(clean_axis)
    print(
        f"\n[cross-check] F_stage_augmented_clean axis (donors with non-NA _clean): "
        f"rho={rho_axis_clean:+.4f}  n={n_axis_clean}"
    )
else:
    rho_axis_clean, p_axis_clean, n_axis_clean = np.nan, np.nan, 0

# --- Verdict ---
delta_rho = rho_full - rho_clean
verdict = "robust" if abs(delta_rho) < 0.05 else "update_manuscript"

print("\n" + "=" * 70)
print(f"FINAL: rho_full = {rho_full:+.4f}  rho_clean = {rho_clean:+.4f}  "
      f"delta_rho = {delta_rho:+.4f}  verdict = {verdict}")
print("=" * 70)

# Write output TSV
out = pd.DataFrame({
    "metric": [
        "spearman_rho_full",
        "spearman_rho_clean",
        "delta_rho",
        "delta_full",
        "delta_clean",
        "n_donors_full",
        "n_donors_clean",
        "n_donors_excluded",
        "p_full",
        "p_clean",
        "F_stage_axis",
        "spearman_rho_axis_clean_xcheck",
        "n_donors_axis_clean_xcheck",
        "verdict",
    ],
    "value": [
        f"{rho_full:.6f}",
        f"{rho_clean:.6f}",
        f"{delta_rho:.6f}",
        f"{delta_full:.6f}",
        f"{delta_clean:.6f}",
        n_full,
        n_clean,
        n_full - n_clean,
        f"{p_full:.3e}",
        f"{p_clean:.3e}",
        F_AXIS,
        f"{rho_axis_clean:.6f}" if not np.isnan(rho_axis_clean) else "NA",
        n_axis_clean,
        verdict,
    ],
})
out.to_csv(OUT_FILE, sep="\t", index=False)
print(f"\n[write] {OUT_FILE}")

# Also write per-stage breakdown for reference
stage_out = OUT_FILE.replace(".tsv", "_per_stage.tsv")
full_stage["set"] = "full"
clean_stage["set"] = "clean"
pd.concat([full_stage, clean_stage], ignore_index=True).to_csv(stage_out, sep="\t", index=False)
print(f"[write] {stage_out}")
