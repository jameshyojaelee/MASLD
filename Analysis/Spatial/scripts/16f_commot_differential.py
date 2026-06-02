#!/usr/bin/env python3
"""
16f_commot_differential.py — Differential COMMOT communication analysis.

Compares Healthy vs Steatotic conditions: per-pathway signaling changes,
spatial range shifts, sender/receiver identity by cell type, zone-specific
signaling, and concordance with existing CellPhoneDB results.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, C2L_PREFIX, load_config, load_zonation_scores,
    strip_c2l_prefix, save_csv, print_header, print_step,
)
from spatial_stats import get_commot_pathways

# Pathway-level obsp matrices are keyed 'commot-CellChat-<PATHWAY>' with no 2nd
# '-' (LR pairs like 'commot-CellChat-HGF-MET' have one); the global aggregate
# is 'commot-CellChat-total-total'. Used to restrict per-pathway extraction.
COMMOT_DB = "CellChat"
COMMOT_PREFIX = f"commot-{COMMOT_DB}-"


def load_commot_adata(condition):
    """Load COMMOT results h5ad for a condition."""
    cond_label = condition.replace(" ", "_").replace("/", "_")
    path = RESULTS_DIR / "commot" / f"adata_commot_{cond_label}.h5ad"
    if not path.exists():
        print(f"  WARNING: {path} not found")
        return None
    adata = sc.read_h5ad(path)
    print(f"  Loaded {condition}: {adata.n_obs} spots")
    return adata


def _pathway_spot_scores(adata, pathway):
    """Per-spot total communication for one pathway (sum over the obsp matrix).

    The COMMOT pathway-level obsp matrix is a spot x spot signaling matrix; the
    per-spot score is its row-sum (total signaling involving that spot). Returns
    a 1-D array aligned to adata.obs_names.
    """
    key = f"{COMMOT_PREFIX}{pathway}"
    mat = adata.obsp[key]
    if hasattr(mat, "toarray"):
        return np.asarray(mat.sum(axis=1)).ravel()
    return np.asarray(mat).sum(axis=1).ravel()


def get_sum_pathway_scores(adata, direction):
    """Per-spot per-pathway sender/receiver scores from the obsm 'sum' table.

    Fixed (F078): COMMOT stores the per-spot pathway sums as a DataFrame in
    ``adata.obsm['commot-CellChat-sum-sender']`` (and ``-receiver``), whose
    COLUMNS are named ``s-<PATHWAY>`` / ``r-<PATHWAY>`` (e.g. 's-TGFb') for real
    pathways and ``s-<LIG>-<REC>`` for LR pairs plus ``s-total-total``. The old
    code read the obsm *key* and parsed 'sum' as the pathway; the per-pathway
    signal is actually in the column names. This returns ``{pathway: 1-D array}``
    restricted to genuine pathway-level columns (no 2nd '-', excludes
    'total-total').

    ``direction`` is 'sender' or 'receiver'.
    """
    key = f"{COMMOT_PREFIX}sum-{direction}"
    if key not in adata.obsm:
        return {}
    df = adata.obsm[key]
    if not hasattr(df, "columns"):
        return {}
    prefix = "s-" if direction == "sender" else "r-"
    out = {}
    for col in df.columns:
        if not col.startswith(prefix):
            continue
        rest = col[len(prefix):]
        if "-" in rest:            # LR pair (lig-rec) -> not a pathway
            continue
        if rest.lower() in ("total", "sum"):
            continue
        out[rest] = np.asarray(df[col].values, dtype=float)
    return out


def extract_pathway_scores(adata, condition, donor_col="sample_id"):
    """Extract per-pathway signaling from genuine pathway-level COMMOT outputs.

    Fixed (F078): only the real pathway-level obsp matrices
    ``commot-CellChat-<PATHWAY>`` are used (LR-pair matrices like
    ``commot-CellChat-HGF-MET`` and the global ``total``/``sum`` aggregates in
    obsp/obsm are excluded). The old code split keys on '-' and took parts[2],
    turning every ligand gene into a fake pathway and ingesting 'total'/'sum'.

    Returns ``(pathway_scores, donor_long)``:
      * ``pathway_scores`` — list of dicts (one per real pathway) with the
        condition-pooled ``total_communication`` and ``mean_communication``
        (kept for back-compat; pooled totals remain pseudoreplicated — see
        F077). ``mean_communication`` is the mean per-spot score.
      * ``donor_long`` — long DataFrame [pathway, condition, donor, mean_score]
        giving the per-donor mean per-spot score, the unit of evidence for the
        donor-level differential test (F077).
    """
    pathways = get_commot_pathways(adata, db=COMMOT_DB)

    donors = (adata.obs[donor_col].astype(str).values
              if donor_col in adata.obs.columns else None)

    pathway_scores = []
    donor_rows = []
    for pw in pathways:
        spot_scores = _pathway_spot_scores(adata, pw)
        pathway_scores.append({
            "pathway": pw,
            "condition": condition,
            "total_communication": float(np.nansum(spot_scores)),
            "mean_communication": float(np.nanmean(spot_scores)),
        })
        if donors is not None:
            s = pd.Series(spot_scores).groupby(donors).mean()
            for donor, val in s.items():
                donor_rows.append({"pathway": pw, "condition": condition,
                                   "donor": donor, "mean_score": float(val)})

    donor_long = pd.DataFrame(donor_rows)
    return pathway_scores, donor_long


def compute_pathway_differential(scores_h, scores_m, donor_h=None, donor_m=None):
    """Per-pathway differential between conditions with a donor-level test.

    Fixed (F077): the headline fold change is now computed on the per-donor mean
    per-spot score (mean over donors), NOT on the spot-pooled total — a single
    donor's slice size / depth no longer drives the FC. We still emit the pooled
    ``total_communication_*`` columns (back-compat for 16h), but the canonical
    ``total_communication_fc`` is the donor-mean ratio and a donor-level
    Mann-Whitney p-value is added (n=2 vs n=2 here — descriptive, see caveat in
    main()). Column names are preserved; ``total_communication_fc`` now means a
    donor-mean ratio rather than a spot-pooled-total ratio.
    """
    df_h = pd.DataFrame(scores_h).set_index("pathway")
    df_m = pd.DataFrame(scores_m).set_index("pathway")

    # Per-donor mean per-spot score, indexed [pathway -> {donor: mean_score}].
    donor_h = donor_h if donor_h is not None else pd.DataFrame()
    donor_m = donor_m if donor_m is not None else pd.DataFrame()

    all_pathways = sorted(set(df_h.index) | set(df_m.index))
    records = []

    for pw in all_pathways:
        record = {"pathway": pw}

        # Pooled totals / spot-means (kept for back-compat; pseudoreplicated).
        for col in ["total_communication", "mean_communication"]:
            val_h = df_h.loc[pw, col] if pw in df_h.index and col in df_h.columns else 0.0
            val_m = df_m.loc[pw, col] if pw in df_m.index and col in df_m.columns else 0.0
            record[f"{col}_healthy"] = val_h
            record[f"{col}_steatotic"] = val_m
            record[f"{col}_delta"] = val_m - val_h

        # Donor-level evidence (F077): per-donor mean per-spot score.
        dh = donor_h.loc[donor_h["pathway"] == pw, "mean_score"].values \
            if len(donor_h) else np.array([])
        dm = donor_m.loc[donor_m["pathway"] == pw, "mean_score"].values \
            if len(donor_m) else np.array([])
        donor_mean_h = float(np.nanmean(dh)) if len(dh) else np.nan
        donor_mean_m = float(np.nanmean(dm)) if len(dm) else np.nan
        record["donor_mean_healthy"] = donor_mean_h
        record["donor_mean_steatotic"] = donor_mean_m
        record["n_donors_healthy"] = int(len(dh))
        record["n_donors_steatotic"] = int(len(dm))

        # Canonical FC = donor-mean ratio (falls back to pooled mean if no donor
        # info). NOTE: same column name as before; statistical meaning is now a
        # donor-mean ratio, not a spot-pooled total ratio.
        if np.isfinite(donor_mean_h) and np.isfinite(donor_mean_m):
            record["total_communication_fc"] = donor_mean_m / max(donor_mean_h, 1e-10)
        else:
            mh = record.get("mean_communication_healthy", 0.0)
            mm = record.get("mean_communication_steatotic", 0.0)
            record["total_communication_fc"] = mm / max(mh, 1e-10)

        # Donor-level Mann-Whitney (descriptive at n=2 vs n=2; padj added below).
        pval = np.nan
        if len(dh) >= 1 and len(dm) >= 1:
            try:
                from scipy.stats import mannwhitneyu
                pval = float(mannwhitneyu(dh, dm, alternative="two-sided").pvalue)
            except ValueError:
                pval = np.nan
        record["donor_mwu_pval"] = pval

        records.append(record)

    diff = pd.DataFrame(records)
    # BH-FDR over donor-level p-values (present even if all-nan).
    diff["donor_mwu_padj"] = np.nan
    if "donor_mwu_pval" in diff and diff["donor_mwu_pval"].notna().any():
        from statsmodels.stats.multitest import multipletests
        mask = diff["donor_mwu_pval"].notna()
        diff.loc[mask, "donor_mwu_padj"] = multipletests(
            diff.loc[mask, "donor_mwu_pval"], method="fdr_bh")[1]
    return diff


def compute_sender_receiver_by_celltype(adata, condition):
    """Aggregate sender/receiver scores by dominant cell type, per real pathway.

    Fixed (F078): reads the per-pathway columns of the obsm 'sum-sender' /
    'sum-receiver' tables (via get_sum_pathway_scores) instead of parsing the
    obsm key as pathway='sum'. Output is now one row per
    (cell_type, real pathway, direction) rather than collapsing to 'sum'.
    """
    if "cell_type_dominant" not in adata.obs.columns:
        return pd.DataFrame()

    cell_types = adata.obs["cell_type_dominant"].copy()
    cell_types = cell_types.apply(strip_c2l_prefix).astype("category")

    records = []
    for direction in ("sender", "receiver"):
        pw_scores = get_sum_pathway_scores(adata, direction)
        for pathway, spot_scores in pw_scores.items():
            for ct in cell_types.unique():
                mask = (cell_types == ct).values
                if mask.sum() == 0:
                    continue
                records.append({
                    "condition": condition,
                    "pathway": pathway,
                    "direction": direction,
                    "cell_type": ct,
                    "mean_score": float(np.nanmean(spot_scores[mask])),
                    "total_score": float(np.nansum(spot_scores[mask])),
                    "n_spots": int(mask.sum()),
                })

    return pd.DataFrame(records)


def compute_zonation_signaling(adata, condition, n_bins=5):
    """Aggregate signaling scores by zonation bins."""
    zon_scores = load_zonation_scores()
    if len(zon_scores) == 0:
        # Try loading from adata.obs
        if "zonation_score" not in adata.obs.columns:
            return pd.DataFrame()
        zon_col = adata.obs["zonation_score"]
    else:
        shared = adata.obs.index.intersection(zon_scores.index)
        if len(shared) == 0:
            return pd.DataFrame()
        zon_col = zon_scores.loc[shared, "zonation_score"] \
            if "zonation_score" in zon_scores.columns else pd.Series(dtype=float)

    if len(zon_col) == 0:
        return pd.DataFrame()

    # Bin zonation scores (1=PP, n_bins=PC)
    zon_bins = pd.qcut(zon_col, q=n_bins, labels=False, duplicates="drop") + 1
    # Per-spot bin aligned to adata order (NaN where no zonation score).
    bin_per_spot = pd.Series(zon_bins).reindex(adata.obs.index)

    records = []
    # Fixed (F078): use per-pathway columns of the obsm 'sum' tables, not the
    # obsm key parsed as pathway='sum'.
    for direction in ("sender", "receiver"):
        pw_scores = get_sum_pathway_scores(adata, direction)
        for pathway, spot_scores in pw_scores.items():
            for b in range(1, n_bins + 1):
                valid = (bin_per_spot == b).values
                if valid.sum() == 0:
                    continue
                records.append({
                    "condition": condition,
                    "pathway": pathway,
                    "direction": direction,
                    "zonation_bin": b,
                    "mean_score": float(np.nanmean(spot_scores[valid])),
                    "n_spots": int(valid.sum()),
                })

    return pd.DataFrame(records)


def compare_with_cellphonedb(commot_diff, cpdb_dir):
    """Compare COMMOT pathway-level results with CellPhoneDB from 05b."""
    concordance = []

    # Load CellPhoneDB differential results
    cpdb_diff_path = cpdb_dir / "differential_lr_pairs.csv"
    if not cpdb_diff_path.exists():
        print("  WARNING: CellPhoneDB differential results not found")
        return pd.DataFrame()

    cpdb_diff = pd.read_csv(cpdb_diff_path, index_col=0)

    # Load CellPhoneDB per-condition results
    cpdb_conditions = {}
    for csv_path in cpdb_dir.glob("ligrec_*.csv"):
        cond = csv_path.stem.replace("ligrec_", "")
        cpdb_conditions[cond] = pd.read_csv(csv_path, index_col=0)

    # Compare top LR pairs
    cpdb_lr_set = set()
    if "lr_pair" in cpdb_diff.columns:
        cpdb_lr_set = set(cpdb_diff["lr_pair"].dropna().astype(str))

    # COMMOT pathways that are also in CellPhoneDB
    commot_pathways = set(commot_diff["pathway"]) if len(commot_diff) > 0 else set()

    concordance.append({
        "metric": "n_cpdb_differential_lr",
        "value": len(cpdb_lr_set),
    })
    concordance.append({
        "metric": "n_commot_pathways",
        "value": len(commot_pathways),
    })

    # Compare pathway-level direction of change
    if len(commot_diff) > 0 and len(cpdb_diff) > 0:
        # Count disease-emergent in CellPhoneDB
        n_emergent = (cpdb_diff["category"] == "disease_emergent").sum() \
            if "category" in cpdb_diff.columns else 0
        n_lost = (cpdb_diff["category"] == "disease_lost").sum() \
            if "category" in cpdb_diff.columns else 0
        concordance.append({
            "metric": "cpdb_disease_emergent",
            "value": n_emergent,
        })
        concordance.append({
            "metric": "cpdb_disease_lost",
            "value": n_lost,
        })

        # COMMOT upregulated pathways
        if "total_communication_fc" in commot_diff.columns:
            n_up = (commot_diff["total_communication_fc"] > 1.5).sum()
            n_down = (commot_diff["total_communication_fc"] < 0.67).sum()
            concordance.append({
                "metric": "commot_upregulated_pathways",
                "value": int(n_up),
            })
            concordance.append({
                "metric": "commot_downregulated_pathways",
                "value": int(n_down),
            })

    return pd.DataFrame(concordance)


def main():
    print_header("16f: Differential COMMOT Communication")

    config = load_config()
    commot_config = config["commot"]
    output_dir = RESULTS_DIR / "commot"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Load COMMOT results ──
    print_step("Loading COMMOT results", 1, 5)
    adata_h = load_commot_adata("Healthy")
    adata_m = load_commot_adata("Steatotic")

    if adata_h is None or adata_m is None:
        print("  ERROR: Need both Healthy and Steatotic COMMOT results. "
              "Run 16c_run_commot.py first.")
        sys.exit(1)

    # ── Per-pathway total signaling change ──
    print_step("Computing per-pathway signaling changes", 2, 5)
    # F077: aggregate COMMOT scores per donor BEFORE the fold change, and report
    # the honest n. GSE192741 is n=2 Healthy (JBO018/JBO022) vs n=2 Steatotic
    # donors; spots within a slice are NOT independent (autocorrelation is the
    # whole point of COMMOT), so the donor is the unit of evidence. The FC is now
    # a donor-mean ratio and donor_mwu_pval is a donor-level test — DESCRIPTIVE
    # at n=2 vs n=2, not a significance claim. Document this limitation.
    scores_h, donor_h = extract_pathway_scores(adata_h, "Healthy",
                                               donor_col="sample_id")
    scores_m, donor_m = extract_pathway_scores(adata_m, "Steatotic",
                                               donor_col="sample_id")
    n_donor_h = donor_h["donor"].nunique() if len(donor_h) else 0
    n_donor_m = donor_m["donor"].nunique() if len(donor_m) else 0
    print(f"  Healthy pathways: {len(scores_h)} (n_donors={n_donor_h})")
    print(f"  Steatotic pathways: {len(scores_m)} (n_donors={n_donor_m})")
    if n_donor_h <= 2 or n_donor_m <= 2:
        print(f"  CAVEAT (F077): n={n_donor_h} vs n={n_donor_m} donors — "
              f"donor-level p-values are DESCRIPTIVE only; no valid "
              f"disease-vs-control significance is possible at this n.")

    diff_df = compute_pathway_differential(scores_h, scores_m, donor_h, donor_m)
    save_csv(diff_df, "differential_communication.csv", subdir="commot")
    # Also persist the per-donor long table so downstream / reviewers can see
    # the n=2 vs n=2 distributions behind the FC.
    donor_all = pd.concat([donor_h, donor_m], ignore_index=True) \
        if (len(donor_h) or len(donor_m)) else pd.DataFrame()
    if len(donor_all):
        save_csv(donor_all, "differential_communication_per_donor.csv",
                 subdir="commot")

    # Report top changes
    if "total_communication_fc" in diff_df.columns:
        diff_sorted = diff_df.sort_values("total_communication_fc", ascending=False)
        print(f"\n  Top upregulated pathways (Steatotic vs Healthy, donor-mean FC):")
        for _, row in diff_sorted.head(5).iterrows():
            print(f"    {row['pathway']}: FC={row['total_communication_fc']:.2f}")
        print(f"\n  Top downregulated pathways:")
        for _, row in diff_sorted.tail(5).iterrows():
            print(f"    {row['pathway']}: FC={row['total_communication_fc']:.2f}")

    # ── Sender/receiver by cell type ──
    print_step("Sender/receiver identity by cell type", 3, 5)
    sr_h = compute_sender_receiver_by_celltype(adata_h, "Healthy")
    sr_m = compute_sender_receiver_by_celltype(adata_m, "Steatotic")

    if len(sr_h) > 0 or len(sr_m) > 0:
        sr_combined = pd.concat([sr_h, sr_m], ignore_index=True)
        save_csv(sr_combined, "sender_receiver_shift.csv", subdir="commot")
        print(f"  Cell type x pathway records: {len(sr_combined)}")

        # Identify biggest shifts
        if len(sr_h) > 0 and len(sr_m) > 0:
            merged = sr_h.merge(sr_m, on=["pathway", "direction", "cell_type"],
                                suffixes=("_h", "_m"), how="outer")
            merged["delta"] = merged["mean_score_m"].fillna(0) - \
                merged["mean_score_h"].fillna(0)
            top_shifts = merged.sort_values("delta", ascending=False).head(5)
            print(f"\n  Top cell type signaling shifts (Steatotic - Healthy):")
            for _, row in top_shifts.iterrows():
                print(f"    {row['cell_type']} / {row['pathway']} ({row['direction']}): "
                      f"delta={row['delta']:.4f}")

    # ── Zone-specific signaling ──
    print_step("Zone-specific signaling", 4, 5)
    zon_h = compute_zonation_signaling(adata_h, "Healthy")
    zon_m = compute_zonation_signaling(adata_m, "Steatotic")

    if len(zon_h) > 0 or len(zon_m) > 0:
        zon_combined = pd.concat([zon_h, zon_m], ignore_index=True)
        save_csv(zon_combined, "zonation_signaling_matrix.csv", subdir="commot")
        print(f"  Zonation x pathway records: {len(zon_combined)}")

    # ── CellPhoneDB concordance ──
    print_step("CellPhoneDB concordance", 5, 5)
    cpdb_dir = RESULTS_DIR / "communication"
    concordance = compare_with_cellphonedb(diff_df, cpdb_dir)
    if len(concordance) > 0:
        save_csv(concordance, "cellphonedb_concordance.csv", subdir="commot")
        for _, row in concordance.iterrows():
            print(f"    {row['metric']}: {row['value']}")

    # ── Summary ──
    print(f"\n  Summary:")
    print(f"    Pathways analyzed: {len(diff_df)}")
    if "total_communication_fc" in diff_df.columns:
        n_up = (diff_df["total_communication_fc"] > 1.5).sum()
        n_down = (diff_df["total_communication_fc"] < 0.67).sum()
        print(f"    Upregulated (FC>1.5): {n_up}")
        print(f"    Downregulated (FC<0.67): {n_down}")

    print_header("16f: Complete")


if __name__ == "__main__":
    main()
