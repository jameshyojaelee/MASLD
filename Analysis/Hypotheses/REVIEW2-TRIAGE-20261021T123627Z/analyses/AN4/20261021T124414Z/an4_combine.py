#!/usr/bin/env python3
"""Combine AN4 family slices, rank focal programs, build the Catalog spatial sidecar.

Light pandas only (login-node safe). Reads the slice outputs written by
an4_spatial.py and the r2 all-16 sensitivity tables; writes to combined/.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
R2 = PROJECT / "Analysis/Multimodal_Program_Projection/candidates/spatial-full-composition-sensitivity-candidate-2026-08-11-r2"
CATALOG_PC = PROJECT / "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_program_context.parquet"
R9_SIDECAR = PROJECT / ("Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/"
                        "catalog_r9d-20261016T105550Z/build_out/sidecar_program_context_r9.tsv")
FOCAL = {"hotspot_hepatocytes_f05c535ae5bbc0b9": "Stromal ECM (IGFBP7)",
         "hotspot_hepatocytes_48f39dd4d817a10e": "Ductular injury (BICC1)"}
OUT = HERE / "combined"
OUT.mkdir(exist_ok=True)


def bh(p: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=p.index)
    v = p.notna()
    x = p[v].to_numpy(float)
    o = np.argsort(x)
    adj = x[o] * len(x) / np.arange(1, len(x) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    vals = np.empty_like(adj)
    vals[o] = np.minimum(adj, 1)
    out[v] = vals
    return out


def main() -> int:
    slices = sorted(HERE.glob("famr2_*_s*"))
    frames, hashes = [], []
    for d in slices:
        for f in d.glob("family_results_*.tsv"):
            if f.name.endswith(".partial.tsv"):
                continue
            frames.append(pd.read_csv(f, sep="\t").assign(slice_dir=d.name))
        for f in d.glob("matched_set_hashes_*.tsv"):
            hashes.append(pd.read_csv(f, sep="\t").assign(slice_dir=d.name))
    fam = pd.concat(frames, ignore_index=True)
    hs = pd.concat(hashes, ignore_index=True)
    # every slice job must have drawn identical matched sets for every program
    nuniq = hs.groupby(["dataset", "program_id"])["matched_set_sha256"].nunique()
    hash_consistent = bool((nuniq == 1).all())
    fam = fam.drop_duplicates(["dataset", "program_id"])
    status = []
    rows = []
    for ds, part in fam.groupby("dataset"):
        part = part.copy()
        t = part["testable"].astype(bool) & part["residual_pvalue"].notna()
        part.loc[t, "residual_padj_family"] = bh(part.loc[t, "residual_pvalue"]).values
        part.loc[t, "rank_excess"] = part.loc[t, "excess"].rank(ascending=False, method="min")
        part.loc[t, "rank_z"] = part.loc[t, "residual_moran_z"].rank(ascending=False, method="min")
        n = int(t.sum())
        ex = part.loc[t, "excess"].to_numpy(float)
        zz = part.loc[t, "residual_moran_z"].to_numpy(float)
        part.loc[t, "pct_excess_le"] = [100 * np.mean(ex <= v) for v in part.loc[t, "excess"]]
        part.loc[t, "pct_z_le"] = [100 * np.mean(zz <= v) for v in part.loc[t, "residual_moran_z"]]
        part["robust_family"] = t & (part["residual_padj_family"] < 0.05) & part["sensitivity_sign_agree"].fillna(False).astype(bool)
        rows.append(part)
        status.append({
            "dataset": ds,
            "n_programs": len(part),
            "n_testable": n,
            "n_with_full_null": int(part["residual_pvalue"].notna().sum()),
            "n_excess_gt0": int((part.loc[t, "excess"] > 0).sum()),
            "n_bh_q_lt_0.05": int((part.loc[t, "residual_padj_family"] < 0.05).sum()),
            "n_robust_q_and_sign": int(part["robust_family"].sum()),
            "n_p_at_floor": int((part.loc[t, "n_null_ge_observed"] == 0).sum()),
            "excess_median": float(np.median(ex)),
            "excess_q90": float(np.quantile(ex, 0.9)),
            "matched_sets_identical_across_slices": hash_consistent,
        })
    fam = pd.concat(rows, ignore_index=True)
    fam.to_csv(OUT / "family_all_programs.tsv", sep="\t", index=False)
    st = pd.DataFrame(status)
    st.to_csv(OUT / "family_summary_by_source.tsv", sep="\t", index=False)
    # by lineage
    t = fam["testable"].astype(bool)
    bylin = fam[t].groupby(["dataset", "cell_type"]).agg(
        n=("program_id", "size"),
        n_excess_gt0=("excess", lambda s: int((s > 0).sum())),
        n_q_lt_005=("residual_padj_family", lambda s: int((s < 0.05).sum())),
        median_excess=("excess", "median"),
    ).reset_index()
    bylin.to_csv(OUT / "family_summary_by_lineage.tsv", sep="\t", index=False)
    foc = fam[fam["program_id"].isin(FOCAL)][[
        "dataset", "program_id", "legacy_program_id", "program_name", "n_measured", "retained_l1_weight",
        "residual_moran_i", "residual_null_mean", "residual_null_sd", "excess", "residual_moran_z",
        "residual_pvalue", "residual_padj_family", "rank_excess", "pct_excess_le", "rank_z", "pct_z_le",
        "sensitivity_sign_agree", "robust_family"]].merge(st[["dataset", "n_testable"]], on="dataset")
    foc.to_csv(OUT / "focal_family_rank.tsv", sep="\t", index=False)

    # r2 primary (two-program) reference and 71%/66% check
    cmp = pd.read_csv(R2 / "comparison_to_primary.tsv", sep="\t")
    cmp["fraction_removed"] = 1 - cmp["full_composition_centered_moran_i"] / cmp["primary_centered_moran_i"]
    cmp.to_csv(OUT / "r2_reduction_check.tsv", sep="\t", index=False)

    # Catalog sidecar: replacement values for the four Visium spatial program-context rows
    pc = pd.read_parquet(CATALOG_PC)
    pc = pc[pc["assay"] == "Visium_spatial_transcriptomics"].copy()
    r9 = pd.read_csv(R9_SIDECAR, sep="\t")
    r2res = pd.read_csv(R2 / "full_composition_program_results.tsv", sep="\t")
    lin = None
    lf = HERE / "focal_lineage" / "focal_lineage_results.tsv"
    if lf.exists():
        lin = pd.read_csv(lf, sep="\t")
    side = []
    for r in pc.itertuples(index=False):
        a = r2res[(r2res.program_id == r.program_uid) & (r2res.dataset == r.dataset_id)].iloc[0]
        f = foc[(foc.program_id == r.program_uid) & (foc.dataset == r.dataset_id)]
        r9row = r9[r9.program_context_id == r.program_context_id].iloc[0]
        rec = {
            "program_context_id": r.program_context_id,
            "program_uid": r.program_uid,
            "program_label": r.program_label,
            "dataset_id": r.dataset_id,
            "assay": r.assay,
            "effect_unit": r.effect_unit,
            "estimate_r8_hepatocyte_only": float(r.estimate),
            "p_value_r8": float(r.p_value),
            "q_value_r8": float(r.q_value),
            "call_state_r8": r.call_state,
            "call_state_r9d": r9row["call_state_r9_candidate"],
            "residualization_new": "all_16_cell2location_q05_plus_log1p_total_counts_and_log1p_detected_genes",
            "estimate_all16": float(a["residual_moran_i"] - a["residual_null_mean"]),
            "residual_moran_i_all16": float(a["residual_moran_i"]),
            "null_mean_all16": float(a["residual_null_mean"]),
            "null_sd_all16": float(a["residual_null_sd"]),
            "p_value_all16": float(a["residual_pvalue"]),
            "q_value_all16_two_program_bh": float(a["residual_padj"]),
            "n_null": int(a["n_null"]),
            "n_measured_genes": int(a["n_measured"]),
            "retained_l1_weight": float(a["retained_l1_weight"]),
            "within_source_support_all16": bool(a["robust"]),
            "evidence_state_r2_table": a["evidence_state"],
            # r9d convention: unsupported -> indeterminate; supported in Vu (no donor key) -> source_dependent
            "call_state_all16": ("indeterminate" if not bool(a["robust"]) else
                                 ("source_dependent" if r.dataset_id == "Vu_et_al_2025" else "supported")),
            "matched_set_sha256": a["matched_set_sha256"],
            "source_file": str(R2 / "full_composition_program_results.tsv"),
            "family_rank_excess": float(f["rank_excess"].iloc[0]) if len(f) else np.nan,
            "family_n_testable": int(f["n_testable"].iloc[0]) if len(f) else np.nan,
            "family_pct_excess_le": float(f["pct_excess_le"].iloc[0]) if len(f) else np.nan,
            "family_q_bh_all_testable": float(f["residual_padj_family"].iloc[0]) if len(f) else np.nan,
        }
        if lin is not None:
            l = lin[(lin.scheme == "content_lineage") & (lin.program_id == r.program_uid) & (lin.dataset == r.dataset_id)]
            if len(l):
                rec.update({
                    "content_lineage_matched": l["match_lineage"].iloc[0],
                    "excess_content_lineage_null": float(l["excess"].iloc[0]),
                    "p_content_lineage_null": float(l["residual_pvalue"].iloc[0]),
                    "q_content_lineage_null_two_program_bh": float(l["residual_padj_two_program"].iloc[0]),
                })
        if not rec["within_source_support_all16"]:
            rec["proposed_state_wording"] = "indeterminate: within-section autocorrelation not above matched genes; not negative"
        elif rec["call_state_all16"] == "source_dependent":
            rec["proposed_state_wording"] = ("source_dependent: within-array autocorrelation above matched genes; "
                                             "arrays have no donor key; not donor-level")
        else:
            rec["proposed_state_wording"] = "within-section autocorrelation above matched genes; not donor-level"
        side.append(rec)
    pd.DataFrame(side).to_csv(OUT / "catalog_spatial_rows_all16_sidecar.tsv", sep="\t", index=False)
    print(st.to_string())
    print(foc.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
