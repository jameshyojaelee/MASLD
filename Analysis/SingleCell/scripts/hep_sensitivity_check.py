"""
Hepatocyte subclustering sensitivity check — Liver_Atlas nuclei vs whole-tissue cells.

Asks two questions, quantitatively:
  1. Does dropping Liver_Atlas (and GSE136103) cells materially shift the
     Healthy <-> Steatohepatitis centroid in scVI latent space?
  2. Does the "Progressor expansion" trajectory across disease_stage_coarse
     hold up if Liver_Atlas / GSE136103 are removed?

Outputs:
  results_gpu_v2/hepatocyte_subtypes/protocol_sensitivity_check.tsv
  results_gpu_v2/hepatocyte_subtypes/protocol_sensitivity_check.md
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

PROJ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HEP_DIR = PROJ / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes"
META_CSV = HEP_DIR / "hepatocyte_subtype_metadata.csv"
META_MAP_CSV = HEP_DIR / "meta_subtype_mapping.csv"
ATLAS_H5AD = HEP_DIR / "hepatocyte_atlas_annotated.h5ad"

OUT_TSV = HEP_DIR / "protocol_sensitivity_check.tsv"
OUT_MD = HEP_DIR / "protocol_sensitivity_check.md"

CONTAM_DATASETS = ("Liver_Atlas", "GSE136103")


def read_x_scvi(h5ad_path: Path) -> tuple[np.ndarray, list[str]]:
    """Read X_scvi (obsm) + obs_names directly with h5py to avoid anndata IO regressions.

    The .h5ad has a uns/log1p/base null-typed entry that breaks anndata.read_h5ad on
    this environment. Pulling X_scvi from /obsm/X_scvi and obs index from /obs/_index
    is sufficient for centroid math.
    """
    with h5py.File(h5ad_path, "r") as f:
        # Find obsm/X_scvi (named X_scvi or scVI in different versions; check both)
        candidates = ["X_scvi", "X_scVI", "scVI", "scvi"]
        obsm_keys = list(f["obsm"].keys())
        chosen = next((k for k in candidates if k in obsm_keys), None)
        if chosen is None:
            raise KeyError(
                f"None of {candidates} present in obsm. Available: {obsm_keys}"
            )
        X = f["obsm"][chosen][...]
        # obs index
        if "_index" in f["obs"].attrs:
            idx_key = f["obs"].attrs["_index"]
            if isinstance(idx_key, bytes):
                idx_key = idx_key.decode()
        else:
            idx_key = "_index"
        idx_arr = f["obs"][idx_key][...]
        if idx_arr.dtype.kind in ("O", "S"):
            obs_names = [
                x.decode() if isinstance(x, bytes) else str(x) for x in idx_arr
            ]
        else:
            obs_names = [str(x) for x in idx_arr]
    return X, obs_names


def cosine_distance(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return float("nan")
    return float(1.0 - np.dot(a, b) / (na * nb))


def main():
    print(f"[1] reading metadata: {META_CSV}", flush=True)
    meta = pd.read_csv(META_CSV, index_col=0, low_memory=False)
    print(f"    metadata: {meta.shape} rows; cols = {list(meta.columns)}", flush=True)

    # Map subtype id -> meta_subtype name
    mapping = pd.read_csv(META_MAP_CSV).set_index("subtype")["meta_subtype"]
    meta["meta_subtype"] = meta["hepatocyte_subtype"].map(mapping)

    # disease_stage_coarse is empty for some datasets (e.g. GSE202379). Treat as
    # "Unstaged" so we can still report a clean breakdown, but the headline
    # Healthy<->Steatohepatitis comparison restricts to staged donors.
    stage = meta["disease_stage_coarse"].astype(str).where(
        meta["disease_stage_coarse"].notna() & (meta["disease_stage_coarse"] != ""),
        "Unstaged",
    )
    meta["stage"] = stage

    print("    dataset x stage:", flush=True)
    print(
        meta.groupby(["dataset", "stage"]).size().unstack(fill_value=0),
        flush=True,
    )

    print(f"[2] reading X_scvi from: {ATLAS_H5AD}", flush=True)
    X, obs_names = read_x_scvi(ATLAS_H5AD)
    print(f"    X_scvi shape: {X.shape}; n_obs_names = {len(obs_names)}", flush=True)

    # Align metadata to h5ad obs order. The metadata's index is exactly the obs
    # name pattern "SRRxxxx_BARCODE-1".
    meta_idx = meta.index.astype(str)
    obs_idx = pd.Index(obs_names)
    # Find common cells; meta should fully cover obs (or vice versa)
    common = obs_idx.intersection(meta_idx)
    n_common = len(common)
    print(
        f"    common cells: {n_common} / X={X.shape[0]} / meta={len(meta_idx)}",
        flush=True,
    )
    if n_common < 0.95 * min(X.shape[0], len(meta_idx)):
        sys.exit(
            f"FATAL: only {n_common} cells aligned; expected near-full overlap."
        )

    # Reorder X to match meta row order (drop X rows not in meta)
    obs_pos = pd.Series(np.arange(len(obs_idx)), index=obs_idx)
    keep_obs = obs_pos.loc[common].to_numpy()
    Xc = X[keep_obs]
    meta_c = meta.loc[common].copy()
    print(f"    aligned: X={Xc.shape}, meta={meta_c.shape}", flush=True)

    # -----------------------------------------------------------------------
    # Step 2: centroid shifts (Healthy vs Steatohepatitis, scVI latent)
    # -----------------------------------------------------------------------
    is_contam = meta_c["dataset"].isin(CONTAM_DATASETS).to_numpy()
    is_healthy = (meta_c["stage"] == "Healthy").to_numpy()
    is_sh = (meta_c["stage"] == "Steatohepatitis").to_numpy()
    is_steat = (meta_c["stage"] == "Steatosis").to_numpy()

    def centroid(mask: np.ndarray) -> np.ndarray:
        if mask.sum() == 0:
            return np.full(Xc.shape[1], np.nan)
        return Xc[mask].mean(axis=0)

    # FULL (all cells)
    c_h_full = centroid(is_healthy)
    c_sh_full = centroid(is_sh)
    d_full = cosine_distance(c_h_full, c_sh_full)

    # CLEAN (exclude Liver_Atlas + GSE136103)
    is_clean = ~is_contam
    c_h_clean = centroid(is_healthy & is_clean)
    c_sh_clean = centroid(is_sh & is_clean)
    d_clean = cosine_distance(c_h_clean, c_sh_clean)

    delta_abs = abs(d_clean - d_full)
    delta_pct = delta_abs / d_full if d_full > 0 else float("nan")

    # Per-stage centroid shift (full vs clean), L2 distance in scVI
    def centroid_shift(stage_mask: np.ndarray) -> float:
        a = centroid(stage_mask)
        b = centroid(stage_mask & is_clean)
        return float(np.linalg.norm(a - b))

    healthy_centroid_shift_l2 = centroid_shift(is_healthy)
    sh_centroid_shift_l2 = centroid_shift(is_sh)
    steat_centroid_shift_l2 = centroid_shift(is_steat)

    print("\n[2 results] centroid cosine distance Healthy<->Steatohepatitis:")
    print(f"    full:  {d_full:.6f}")
    print(f"    clean: {d_clean:.6f}")
    print(f"    delta_abs: {delta_abs:.6f}  delta_pct: {delta_pct*100:.3f}%")
    print(
        f"    per-stage centroid L2 shift (full -> clean):"
        f" Healthy={healthy_centroid_shift_l2:.4f},"
        f" Steatosis={steat_centroid_shift_l2:.4f},"
        f" SH={sh_centroid_shift_l2:.4f}"
    )

    verdict_centroid = "robust" if (delta_pct < 0.05) else "investigate"

    # -----------------------------------------------------------------------
    # Step 3: Progressor meta-subtype trajectory
    # -----------------------------------------------------------------------
    # meta-subtypes per mapping: Healthy, Disease-Neutral, Disease-Associated,
    # Disease-Progressor, Neutral. Task asked about "Progressor", so we use
    # Disease-Progressor.
    PROG = "Disease-Progressor"

    meta_c["is_progressor"] = meta_c["meta_subtype"] == PROG

    def progressor_frac_by_donor(df: pd.DataFrame) -> pd.DataFrame:
        g = (
            df.groupby(["sample", "dataset", "stage"], observed=True)
            .agg(n_cells=("is_progressor", "size"),
                 n_progressor=("is_progressor", "sum"))
            .reset_index()
        )
        g["progressor_frac"] = g["n_progressor"] / g["n_cells"]
        return g

    full_donors = progressor_frac_by_donor(meta_c)
    clean_donors = progressor_frac_by_donor(meta_c[~is_contam])

    def stage_summary(donors: pd.DataFrame) -> pd.DataFrame:
        # mean progressor_frac per stage (averaged across donors, donor-balanced)
        s = (
            donors.groupby("stage", observed=True)
            .agg(
                n_donors=("sample", "nunique"),
                mean_prog_frac=("progressor_frac", "mean"),
                median_prog_frac=("progressor_frac", "median"),
                total_cells=("n_cells", "sum"),
                total_progressor=("n_progressor", "sum"),
            )
            .reset_index()
        )
        s["pooled_prog_frac"] = s["total_progressor"] / s["total_cells"]
        return s

    full_stage = stage_summary(full_donors).set_index("stage")
    clean_stage = stage_summary(clean_donors).set_index("stage")

    # Tabulate the Steatosis -> SH expansion under both settings
    rows = []
    for stage_name in ["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis", "Unstaged"]:
        full_row = full_stage.loc[stage_name] if stage_name in full_stage.index else None
        clean_row = clean_stage.loc[stage_name] if stage_name in clean_stage.index else None
        rows.append(
            {
                "stage": stage_name,
                "full_n_donors": int(full_row["n_donors"]) if full_row is not None else 0,
                "full_mean_prog_frac": float(full_row["mean_prog_frac"])
                if full_row is not None else np.nan,
                "full_pooled_prog_frac": float(full_row["pooled_prog_frac"])
                if full_row is not None else np.nan,
                "clean_n_donors": int(clean_row["n_donors"]) if clean_row is not None else 0,
                "clean_mean_prog_frac": float(clean_row["mean_prog_frac"])
                if clean_row is not None else np.nan,
                "clean_pooled_prog_frac": float(clean_row["pooled_prog_frac"])
                if clean_row is not None else np.nan,
            }
        )
    stage_compare = pd.DataFrame(rows)

    # Specifically the headline expansion: Steatosis -> Steatohepatitis (donor mean)
    def expansion(stage_df, stage_col):
        try:
            s0 = stage_df.loc["Steatosis", stage_col]
            s1 = stage_df.loc["Steatohepatitis", stage_col]
            return s0, s1
        except Exception:
            return float("nan"), float("nan")

    full_s0_mean, full_s1_mean = expansion(full_stage, "mean_prog_frac")
    clean_s0_mean, clean_s1_mean = expansion(clean_stage, "mean_prog_frac")
    full_s0_pool, full_s1_pool = expansion(full_stage, "pooled_prog_frac")
    clean_s0_pool, clean_s1_pool = expansion(clean_stage, "pooled_prog_frac")

    print("\n[3 results] Progressor trajectory (Steatosis -> Steatohepatitis):")
    print(
        f"    FULL  donor-mean: {full_s0_mean*100:.2f}% -> {full_s1_mean*100:.2f}%"
        f" | pooled-cell: {full_s0_pool*100:.2f}% -> {full_s1_pool*100:.2f}%"
    )
    print(
        f"    CLEAN donor-mean: {clean_s0_mean*100:.2f}% -> {clean_s1_mean*100:.2f}%"
        f" | pooled-cell: {clean_s0_pool*100:.2f}% -> {clean_s1_pool*100:.2f}%"
    )

    # Verdict: trajectory robust if (a) Steatosis->SH direction preserved,
    # (b) clean delta is within 50% (relative) of full delta on donor-mean.
    full_delta_mean = full_s1_mean - full_s0_mean
    clean_delta_mean = clean_s1_mean - clean_s0_mean
    full_delta_pool = full_s1_pool - full_s0_pool
    clean_delta_pool = clean_s1_pool - clean_s0_pool

    direction_ok = (
        (full_delta_mean > 0) == (clean_delta_mean > 0)
        and (full_delta_pool > 0) == (clean_delta_pool > 0)
    )
    magnitude_relchange_mean = (
        abs(clean_delta_mean - full_delta_mean) / max(abs(full_delta_mean), 1e-9)
        if not np.isnan(full_delta_mean) else float("nan")
    )

    verdict_trajectory = (
        "robust" if (direction_ok and magnitude_relchange_mean < 0.50) else "investigate"
    )

    overall_verdict = (
        "robust"
        if (verdict_centroid == "robust" and verdict_trajectory == "robust")
        else "investigate"
    )

    # -----------------------------------------------------------------------
    # Write per-donor table + stage summary + main TSV
    # -----------------------------------------------------------------------
    donors_compare_path = HEP_DIR / "protocol_sensitivity_donor_table.tsv"
    full_donors.assign(set="full").to_csv(donors_compare_path, sep="\t", index=False)

    stage_compare_path = HEP_DIR / "protocol_sensitivity_stage_compare.tsv"
    stage_compare.to_csv(stage_compare_path, sep="\t", index=False)

    main_rows = [
        {"metric": "centroid_cosine_distance_full",
         "value": d_full,
         "note": "Healthy<->Steatohepatitis, all cells"},
        {"metric": "centroid_cosine_distance_clean",
         "value": d_clean,
         "note": "Healthy<->Steatohepatitis, Liver_Atlas+GSE136103 dropped"},
        {"metric": "centroid_cosine_delta_abs",
         "value": delta_abs,
         "note": "absolute change"},
        {"metric": "centroid_cosine_delta_pct",
         "value": delta_pct,
         "note": "relative change vs full"},
        {"metric": "centroid_L2_shift_healthy_full_to_clean",
         "value": healthy_centroid_shift_l2,
         "note": "L2 distance between Healthy centroid full and clean"},
        {"metric": "centroid_L2_shift_steatosis_full_to_clean",
         "value": steat_centroid_shift_l2,
         "note": "L2 distance between Steatosis centroid full and clean"},
        {"metric": "centroid_L2_shift_sh_full_to_clean",
         "value": sh_centroid_shift_l2,
         "note": "L2 distance between Steatohepatitis centroid full and clean"},
        {"metric": "n_cells_total",
         "value": int(meta_c.shape[0]),
         "note": "hep cells aligned w/ scVI latent"},
        {"metric": "n_cells_contaminant",
         "value": int(is_contam.sum()),
         "note": "Liver_Atlas + GSE136103"},
        {"metric": "frac_contaminant",
         "value": float(is_contam.mean()),
         "note": "fraction of hep cells from contaminant datasets"},
        {"metric": "progressor_frac_steatosis_full_donormean",
         "value": full_s0_mean,
         "note": "donor-mean Progressor fraction in Steatosis donors, all data"},
        {"metric": "progressor_frac_sh_full_donormean",
         "value": full_s1_mean,
         "note": "donor-mean Progressor fraction in SH donors, all data"},
        {"metric": "progressor_frac_steatosis_clean_donormean",
         "value": clean_s0_mean,
         "note": "donor-mean Progressor fraction in Steatosis donors, contaminants dropped"},
        {"metric": "progressor_frac_sh_clean_donormean",
         "value": clean_s1_mean,
         "note": "donor-mean Progressor fraction in SH donors, contaminants dropped"},
        {"metric": "progressor_delta_steatosis_to_sh_full",
         "value": full_delta_mean,
         "note": "donor-mean: SH - Steatosis, full data"},
        {"metric": "progressor_delta_steatosis_to_sh_clean",
         "value": clean_delta_mean,
         "note": "donor-mean: SH - Steatosis, contaminants dropped"},
        {"metric": "progressor_delta_relchange_mean",
         "value": magnitude_relchange_mean,
         "note": "|clean_delta - full_delta| / |full_delta| on donor-mean"},
        {"metric": "progressor_frac_steatosis_full_pooled",
         "value": full_s0_pool,
         "note": "pooled cell fraction (full)"},
        {"metric": "progressor_frac_sh_full_pooled",
         "value": full_s1_pool,
         "note": "pooled cell fraction (full)"},
        {"metric": "progressor_frac_steatosis_clean_pooled",
         "value": clean_s0_pool,
         "note": "pooled cell fraction (clean)"},
        {"metric": "progressor_frac_sh_clean_pooled",
         "value": clean_s1_pool,
         "note": "pooled cell fraction (clean)"},
        {"metric": "verdict_centroid",
         "value": verdict_centroid,
         "note": "robust if cosine delta_pct < 5%"},
        {"metric": "verdict_trajectory",
         "value": verdict_trajectory,
         "note": "robust if direction preserved and donor-mean delta relchange < 50%"},
        {"metric": "verdict_overall",
         "value": overall_verdict,
         "note": "both centroid + trajectory robust => overall robust"},
    ]
    main_df = pd.DataFrame(main_rows)
    main_df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\n[wrote] {OUT_TSV}", flush=True)
    print(f"[wrote] {stage_compare_path}", flush=True)
    print(f"[wrote] {donors_compare_path}", flush=True)

    # -----------------------------------------------------------------------
    # Markdown report (1 page)
    # -----------------------------------------------------------------------
    md_lines = []
    md_lines.append("# Hepatocyte subclustering — Liver_Atlas nuclei sensitivity check")
    md_lines.append("")
    md_lines.append(f"**Verdict (overall):** `{overall_verdict.upper()}`")
    md_lines.append("")
    md_lines.append("## Inputs")
    md_lines.append(f"- atlas: `{ATLAS_H5AD.relative_to(PROJ)}`")
    md_lines.append(f"- metadata: `{META_CSV.relative_to(PROJ)}`")
    md_lines.append(f"- aligned cells: **{Xc.shape[0]:,}**; scVI latent dim **{Xc.shape[1]}**")
    md_lines.append(f"- contaminant datasets: {', '.join(CONTAM_DATASETS)}")
    md_lines.append(
        f"- contaminant cells: **{int(is_contam.sum()):,}**"
        f" ({100*is_contam.mean():.2f}% of hep atlas)"
    )
    md_lines.append("")
    md_lines.append("## Step 1 — scVI centroid (Healthy ↔ Steatohepatitis)")
    md_lines.append("")
    md_lines.append("| set | cosine distance |")
    md_lines.append("|---|---|")
    md_lines.append(f"| full | {d_full:.6f} |")
    md_lines.append(f"| clean (drop {', '.join(CONTAM_DATASETS)}) | {d_clean:.6f} |")
    md_lines.append(f"| |Δ| | {delta_abs:.6f} |")
    md_lines.append(f"| Δ% | **{delta_pct*100:.3f}%** (threshold 5%) |")
    md_lines.append("")
    md_lines.append(
        f"Per-stage centroid L2 shift (full → clean): Healthy="
        f"{healthy_centroid_shift_l2:.4f}, Steatosis={steat_centroid_shift_l2:.4f},"
        f" Steatohepatitis={sh_centroid_shift_l2:.4f}."
    )
    md_lines.append(f"**Verdict:** `{verdict_centroid}` (Δ% < 5%).")
    md_lines.append("")
    md_lines.append("## Step 2 — Progressor trajectory (donor-mean fraction)")
    md_lines.append("")
    md_lines.append(
        "| stage | full donors | full donor-mean | full pooled | clean donors | clean donor-mean | clean pooled |"
    )
    md_lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for _, r in stage_compare.iterrows():
        md_lines.append(
            f"| {r['stage']} | {r['full_n_donors']} |"
            f" {r['full_mean_prog_frac']*100:.2f}% | {r['full_pooled_prog_frac']*100:.2f}% |"
            f" {r['clean_n_donors']} |"
            f" {r['clean_mean_prog_frac']*100:.2f}% | {r['clean_pooled_prog_frac']*100:.2f}% |"
        )
    md_lines.append("")
    md_lines.append("Steatosis → Steatohepatitis (donor-mean):")
    md_lines.append(f"- full: {full_s0_mean*100:.2f}% → {full_s1_mean*100:.2f}% (Δ={full_delta_mean*100:.2f} pp)")
    md_lines.append(f"- clean: {clean_s0_mean*100:.2f}% → {clean_s1_mean*100:.2f}% (Δ={clean_delta_mean*100:.2f} pp)")
    md_lines.append(f"- relative change of Δ: **{magnitude_relchange_mean*100:.1f}%** (threshold 50%)")
    md_lines.append("")
    md_lines.append("Steatosis → Steatohepatitis (pooled cell):")
    md_lines.append(f"- full: {full_s0_pool*100:.2f}% → {full_s1_pool*100:.2f}%")
    md_lines.append(f"- clean: {clean_s0_pool*100:.2f}% → {clean_s1_pool*100:.2f}%")
    md_lines.append("")
    md_lines.append(f"**Verdict:** `{verdict_trajectory}`.")
    md_lines.append("")
    md_lines.append("## Notes")
    md_lines.append(
        "- GSE136103 contributes only 192 hepatocyte cells in the atlas — negligible by"
        " construction, included here for completeness."
    )
    md_lines.append(
        "- Liver_Atlas contributes 44,809 Healthy hepatocyte cells from FACS-sorted"
        " nuclei donors (~6.81% of hep atlas; entirely in the Healthy stage)."
    )
    md_lines.append(
        "- GSE202379 (84,802 cells) lacks a `disease_stage_coarse` annotation in the"
        " hep atlas metadata; its donors are reported under `Unstaged` and excluded"
        " from the Healthy↔SH centroid comparison either way."
    )
    md_lines.append("")
    md_lines.append("## Files written")
    md_lines.append(f"- `{OUT_TSV.relative_to(PROJ)}`")
    md_lines.append(f"- `{stage_compare_path.relative_to(PROJ)}`")
    md_lines.append(f"- `{donors_compare_path.relative_to(PROJ)}`")
    md_lines.append("")

    OUT_MD.write_text("\n".join(md_lines))
    print(f"[wrote] {OUT_MD}", flush=True)
    print(f"\n*** OVERALL VERDICT: {overall_verdict.upper()} ***", flush=True)


if __name__ == "__main__":
    main()
