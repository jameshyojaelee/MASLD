#!/usr/bin/env python3
"""Deduplicated 16-lineage Fig 4C program and gene specificity on unique_counts/pseudobulk."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
UC_DIR = ROOT / "Analysis/Hypotheses/NUCLEUS-FIX-20261107/unique_counts"
PB_DIR = UC_DIR / "pseudobulk"
OUT_DIR = ROOT / "Analysis/Hypotheses/NUCLEUS-FIX-20261107/candidates/sc_20261118"
HOTSPOT_DIR = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"

SOURCE_TO_LINEAGE = {
    "hepatocytes": "Hepatocytes",
    "cholangiocytes": "Cholangiocytes",
    "fibroblasts": "Fibroblasts",
    "macrophages": "Macrophages",
    "endothelial": "Endothelial_cells",
    "tcells": "T_cells",
}
SIX_LINEAGES = list(SOURCE_TO_LINEAGE.values())

MARKERS_AND_FIG4_GENES = [
    "ALB", "APOB", "HAL", "GLUL", "CPS1", "ASS1", "PCK1", "CYP3A4",
    "KRT19", "KRT7", "EPCAM", "SOX9", "SPP1", "BICC1", "CFTR", "TACSTD2",
    "PDGFRB", "COL1A1", "COL1A2", "COL3A1", "DCN", "ACTA2", "IGFBP7", "TIMP1", "LUM",
    "CD163", "MARCO", "TREM2", "CD68", "LYZ", "C1QA", "VSIG4", "CSF1R",
    "PECAM1", "STAB2", "VWF", "CDH5", "CLEC4G", "CLEC4M", "OIT3",
    "CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "GZMB", "KLRB1", "PRF1",
    "MS4A1", "CD79A", "BANK1", "JCHAIN", "MZB1", "SDC1",
    "FCER1A", "CD1C", "CLEC9A", "XCR1", "LILRA4", "TCF4", "S100A8", "S100A9", "CSF3R",
    "CLC", "MS4A2", "CPA3", "TPSAB1", "CD14", "FCGR3A", "CCR2",
    "CCL20", "CCL19", "CCL2", "CXCL8", "CXCL10",
    "PNPLA3", "HSD17B13", "SERPINA1", "TM6SF2", "MBOAT7", "GCKR", "MARC1",
    "MAT1A", "BHMT", "MTHFR", "MTR", "MTAP", "AMD1", "GNMT", "CTH",
]


def to_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def norm_key(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9+]", "_", s)


def load_pseudobulks(donor_meta: pd.DataFrame, lib_meta: pd.DataFrame):
    n_cols = [c for c in lib_meta.columns if c.startswith("n_") and c not in {"n_cells", "n_source_records", "n_donor_cells"}]
    col_map = {norm_key(c[2:]): c for c in n_cols}
    lib_idx = lib_meta.set_index("sample")

    pb_paths = sorted(PB_DIR.glob("*_pseudobulk.csv"))
    assert len(pb_paths) == 16, f"Expected 16 pseudobulk CSVs, found {len(pb_paths)}"

    genes_ref = None
    lineage_donor_logcpm = {}
    lineage_donor_info = {}

    for p in pb_paths:
        ct = p.name.replace("_pseudobulk.csv", "")
        df = pd.read_csv(p)
        genes = df.iloc[:, 0].astype(str)
        mat = df.iloc[:, 1:].to_numpy(dtype=np.float64)
        if genes.duplicated().any():
            mat_df = pd.DataFrame(mat, index=genes, columns=df.columns[1:]).groupby(level=0, sort=False).sum()
            genes = mat_df.index.astype(str)
            mat = mat_df.to_numpy(dtype=np.float64)
        if genes_ref is None:
            genes_ref = pd.Index(genes)
        else:
            assert (pd.Index(genes) == genes_ref).all(), f"Gene order mismatch in {p.name}"

        caps = list(df.columns[1:])
        lm = lib_idx.loc[caps].copy()
        ncell_col = col_map[norm_key(ct)]
        lm["ct_cells"] = lm[ncell_col].astype(int)

        # Store separate collapsed donor log2(CPM+1) matrices for prep_eligible vs all_preps
        for prep_mode in ("prep_eligible", "all_preps"):
            if prep_mode == "prep_eligible":
                mask = to_bool(lm["prep_eligible"]).to_numpy()
            else:
                mask = np.ones(len(lm), dtype=bool)
            lm_sub = lm.iloc[mask]
            mat_sub = mat[:, mask]
            donors_u = pd.Index(lm_sub["donor"].astype(str).unique())
            donor_counts = np.zeros((len(genes_ref), len(donors_u)), dtype=np.float64)
            donor_cells = np.zeros(len(donors_u), dtype=np.int64)
            donor_caps = np.zeros(len(donors_u), dtype=np.int64)
            d_arr = lm_sub["donor"].astype(str).to_numpy()
            c_arr = lm_sub["ct_cells"].to_numpy(dtype=np.int64)
            for j, d in enumerate(donors_u):
                idx = np.where(d_arr == d)[0]
                donor_counts[:, j] = mat_sub[:, idx].sum(axis=1)
                donor_cells[j] = int(c_arr[idx].sum())
                donor_caps[j] = len(idx)

            totals = donor_counts.sum(axis=0)
            valid = totals > 0
            cpm = np.zeros_like(donor_counts)
            cpm[:, valid] = 1e6 * donor_counts[:, valid] / totals[valid][None, :]
            logcpm = np.log2(cpm + 1.0)  # shape (n_genes, n_donors)

            info = pd.DataFrame({
                "donor": donors_u,
                "lineage": ct,
                "n_lineage_cells": donor_cells,
                "n_captures": donor_caps,
                "lib_size": totals,
            }).merge(donor_meta, on="donor", how="left")
            lineage_donor_logcpm[(ct, prep_mode)] = (logcpm.T, donors_u)  # (n_donors, n_genes)
            lineage_donor_info[(ct, prep_mode)] = info

    return genes_ref, lineage_donor_logcpm, lineage_donor_info


def score_programs(
    genes: pd.Index,
    lineage_donor_logcpm: dict,
    lineage_donor_info: dict,
    reg: pd.DataFrame,
    mem: pd.DataFrame,
    arm_mask_fn,
    prep_mode: str,
    min_cells: int,
    lineages_subset: list[str] | None = None,
    min_donors_per_lineage: int = 2,
):
    cts_all = sorted({ct for ct, pm in lineage_donor_logcpm if pm == prep_mode})
    if lineages_subset is not None:
        cts_all = [c for c in cts_all if c in set(lineages_subset)]

    means = {}
    donor_counts = {}
    for ct in cts_all:
        mat, donors_u = lineage_donor_logcpm[(ct, prep_mode)]
        info = lineage_donor_info[(ct, prep_mode)]
        keep = arm_mask_fn(info) & (info["n_lineage_cells"] >= min_cells) & (info["lib_size"] > 0)
        n_k = int(keep.sum())
        donor_counts[ct] = n_k
        if n_k >= min_donors_per_lineage:
            means[ct] = mat[keep.to_numpy()].mean(axis=0)

    active_cts = sorted(means.keys())
    mean_df = pd.DataFrame(means, index=genes).T  # rows = active_cts, cols = genes
    sd = mean_df.std(axis=0, ddof=1)
    good = sd > 0
    z_df = (mean_df.loc[:, good] - mean_df.loc[:, good].mean(axis=0)) / sd[good]
    z_cols = z_df.columns
    col_idx = {g: i for i, g in enumerate(z_cols)}
    z_mat = z_df.to_numpy(dtype=np.float64)  # (n_ct, n_good_genes)

    scores = np.full((len(reg), len(active_cts)), np.nan, dtype=np.float64)
    cov_w = np.zeros(len(reg), dtype=np.float64)
    for i, r in enumerate(reg.itertuples(index=False)):
        sub_m = mem[mem["program_uid"] == r.program_uid]
        idx_list = []
        w_list = []
        for g, w in zip(sub_m["mapped_symbol"].astype(str), sub_m["original_l1_weight"].astype(float)):
            j = col_idx.get(g)
            if j is not None:
                idx_list.append(j)
                w_list.append(w)
        if not w_list:
            continue
        w_arr = np.asarray(w_list, dtype=np.float64)
        w_arr = w_arr / w_arr.sum()
        cov_w[i] = sum(w_list) / float(sub_m["original_l1_weight"].sum())
        scores[i, :] = z_mat[:, idx_list] @ w_arr

    score_df = pd.DataFrame(scores, index=reg["program_uid"], columns=active_cts)
    return score_df, mean_df, z_df, donor_counts, cov_w


def bootstrap_specificity(
    genes: pd.Index,
    lineage_donor_logcpm: dict,
    lineage_donor_info: dict,
    donor_meta: pd.DataFrame,
    reg: pd.DataFrame,
    mem: pd.DataFrame,
    arm_mask_fn,
    prep_mode: str,
    min_cells: int,
    n_boot: int = 1000,
    seed: int = 20261118,
):
    rng = np.random.default_rng(seed)
    dm_elig = donor_meta[arm_mask_fn(donor_meta)].copy()
    strata = dm_elig.groupby(["dataset", "condition_harmonized"]).groups
    strata_arrays = {k: dm_elig.loc[idx, "donor"].astype(str).to_numpy() for k, idx in strata.items()}
    all_elig_donors = set(dm_elig["donor"].astype(str))

    # Pre-index matrices for fast bootstrap across all 16 lineages & 6 source lineages
    cts_16 = sorted({ct for ct, pm in lineage_donor_logcpm if pm == prep_mode})
    donor_lookup = {}
    for ct in cts_16:
        mat, donors_u = lineage_donor_logcpm[(ct, prep_mode)]
        info = lineage_donor_info[(ct, prep_mode)]
        keep = (info["donor"].isin(all_elig_donors)) & (info["n_lineage_cells"] >= min_cells) & (info["lib_size"] > 0)
        d_map = {}
        for row_idx in np.where(keep.to_numpy())[0]:
            d_map[str(donors_u[row_idx])] = mat[row_idx]
        donor_lookup[ct] = d_map

    # Precompute program gene weights indexed over `genes`
    g_map = {g: i for i, g in enumerate(genes)}
    prog_gene_idx = []
    prog_gene_w = []
    prog_src = []
    for r in reg.itertuples(index=False):
        sub_m = mem[mem["program_uid"] == r.program_uid]
        idxs, ws = [], []
        for g, w in zip(sub_m["mapped_symbol"].astype(str), sub_m["original_l1_weight"].astype(float)):
            j = g_map.get(g)
            if j is not None:
                idxs.append(j)
                ws.append(w)
        prog_gene_idx.append(np.asarray(idxs, dtype=np.int64))
        ws_arr = np.asarray(ws, dtype=np.float64)
        prog_gene_w.append(ws_arr / ws_arr.sum())
        prog_src.append(SOURCE_TO_LINEAGE[r.cell_type])

    n_prog = len(reg)
    frac_6 = np.zeros(n_boot, dtype=np.float64)
    frac_16 = np.zeros(n_boot, dtype=np.float64)
    prog_hits_6 = np.zeros(n_prog, dtype=np.int64)
    prog_hits_16 = np.zeros(n_prog, dtype=np.int64)
    valid_boot_6 = 0
    valid_boot_16 = 0

    for b in range(n_boot):
        boot_donors = []
        for arr in strata_arrays.values():
            sampled = rng.choice(arr, size=len(arr), replace=True)
            boot_donors.extend(sampled.tolist())
        # Count multiplicity of each donor in bootstrap draw
        d_counts = {}
        for d in boot_donors:
            d_counts[d] = d_counts.get(d, 0) + 1

        ct_means = {}
        for ct in cts_16:
            d_map = donor_lookup[ct]
            hits = [(d_map[d], cnt) for d, cnt in d_counts.items() if d in d_map]
            tot_cnt = sum(cnt for _, cnt in hits)
            if tot_cnt >= 2:
                acc = np.zeros(len(genes), dtype=np.float64)
                for vec, cnt in hits:
                    acc += cnt * vec
                ct_means[ct] = acc / float(tot_cnt)

        # Evaluate across 6 source lineages
        if all(ct in ct_means for ct in SIX_LINEAGES):
            M6 = np.vstack([ct_means[ct] for ct in SIX_LINEAGES])
            sd6 = M6.std(axis=0, ddof=1)
            good6 = sd6 > 0
            Z6 = np.zeros_like(M6)
            Z6[:, good6] = (M6[:, good6] - M6[:, good6].mean(axis=0)) / sd6[good6]
            hits = 0
            for i in range(n_prog):
                idx = prog_gene_idx[i]
                g_ok = good6[idx]
                if not g_ok.any():
                    continue
                w = prog_gene_w[i][g_ok]
                w = w / w.sum()
                sc = Z6[:, idx[g_ok]] @ w
                top_ct = SIX_LINEAGES[int(np.argmax(sc))]
                if top_ct == prog_src[i]:
                    hits += 1
                    prog_hits_6[i] += 1
            frac_6[valid_boot_6] = hits / float(n_prog)
            valid_boot_6 += 1

        # Evaluate across all active lineages (of 16)
        active_16 = sorted(ct_means.keys())
        if len(active_16) >= 6:
            M16 = np.vstack([ct_means[ct] for ct in active_16])
            sd16 = M16.std(axis=0, ddof=1)
            good16 = sd16 > 0
            Z16 = np.zeros_like(M16)
            Z16[:, good16] = (M16[:, good16] - M16[:, good16].mean(axis=0)) / sd16[good16]
            hits16 = 0
            for i in range(n_prog):
                idx = prog_gene_idx[i]
                g_ok = good16[idx]
                if not g_ok.any():
                    continue
                w = prog_gene_w[i][g_ok]
                w = w / w.sum()
                sc = Z16[:, idx[g_ok]] @ w
                top_ct = active_16[int(np.argmax(sc))]
                if top_ct == prog_src[i]:
                    hits16 += 1
                    prog_hits_16[i] += 1
            frac_16[valid_boot_16] = hits16 / float(n_prog)
            valid_boot_16 += 1

    f6 = frac_6[:valid_boot_6]
    f16 = frac_16[:valid_boot_16]
    return {
        "n_boot_valid_6": valid_boot_6,
        "boot_med_6": float(np.median(f6)) if len(f6) else np.nan,
        "boot_ci_low_6": float(np.quantile(f6, 0.025)) if len(f6) else np.nan,
        "boot_ci_high_6": float(np.quantile(f6, 0.975)) if len(f6) else np.nan,
        "n_boot_valid_16": valid_boot_16,
        "boot_med_16": float(np.median(f16)) if len(f16) else np.nan,
        "boot_ci_low_16": float(np.quantile(f16, 0.025)) if len(f16) else np.nan,
        "boot_ci_high_16": float(np.quantile(f16, 0.975)) if len(f16) else np.nan,
        "prog_boot_rate_6": prog_hits_6 / float(max(valid_boot_6, 1)),
        "prog_boot_rate_16": prog_hits_16 / float(max(valid_boot_16, 1)),
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    donor_meta = pd.read_csv(UC_DIR / "metadata/donor_metadata.tsv", sep="\t")
    lib_meta = pd.read_csv(UC_DIR / "metadata/library_metadata.tsv", sep="\t")

    donor_meta["exclude_flag"] = to_bool(donor_meta["exclude"])
    donor_meta["binary_elig_flag"] = to_bool(donor_meta["binary_eligible"])

    reg = pd.read_csv(HOTSPOT_DIR / "program_registry_v2.tsv", sep="\t")
    mem = pd.read_csv(HOTSPOT_DIR / "program_membership_v2.tsv", sep="\t")
    assert len(reg) == 117

    genes, lineage_donor_logcpm, lineage_donor_info = load_pseudobulks(donor_meta, lib_meta)

    arms = [
        (
            "excl_GSE136103_hm_min10",
            lambda d: (~to_bool(d["exclude"])) & to_bool(d["binary_eligible"]) & d["condition_harmonized"].isin(["Healthy", "MASLD"]),
            "prep_eligible",
            10,
            True,
        ),
        (
            "excl_GSE136103_hm_min30",
            lambda d: (~to_bool(d["exclude"])) & to_bool(d["binary_eligible"]) & d["condition_harmonized"].isin(["Healthy", "MASLD"]),
            "prep_eligible",
            30,
            False,
        ),
        (
            "incl_GSE136103_hm_min10",
            lambda d: (d["dataset"] != "Liver_Atlas") & to_bool(d["binary_eligible"]) & d["condition_harmonized"].isin(["Healthy", "MASLD"]),
            "prep_eligible",
            10,
            True,
        ),
        (
            "incl_GSE136103_hm_min30",
            lambda d: (d["dataset"] != "Liver_Atlas") & to_bool(d["binary_eligible"]) & d["condition_harmonized"].isin(["Healthy", "MASLD"]),
            "prep_eligible",
            30,
            False,
        ),
        (
            "all_cohorts_hm_allpreps_min10",
            lambda d: d["condition_harmonized"].isin(["Healthy", "MASLD", "Other_liver_disease", "Unresolved"]),
            "all_preps",
            10,
            False,
        ),
    ]

    match_rows = []
    prog_summary_tables = {}
    long_16_rows = []
    primary_mean16 = None
    primary_z16 = None
    primary_dcounts = None

    for arm_name, mask_fn, prep_mode, min_cells, do_boot in arms:
        sc16, mean16, z16, dcounts16, cov16 = score_programs(
            genes, lineage_donor_logcpm, lineage_donor_info, reg, mem,
            mask_fn, prep_mode, min_cells, lineages_subset=None, min_donors_per_lineage=2,
        )
        sc6, mean6, z6, dcounts6, cov6 = score_programs(
            genes, lineage_donor_logcpm, lineage_donor_info, reg, mem,
            mask_fn, prep_mode, min_cells, lineages_subset=SIX_LINEAGES, min_donors_per_lineage=2,
        )

        if arm_name == "excl_GSE136103_hm_min10":
            primary_mean16 = mean16
            primary_z16 = z16
            primary_dcounts = dcounts16

        boot = None
        if do_boot:
            boot = bootstrap_specificity(
                genes, lineage_donor_logcpm, lineage_donor_info, donor_meta, reg, mem,
                mask_fn, prep_mode, min_cells, n_boot=1000, seed=20261118,
            )

        rows = []
        for idx_p, r in enumerate(reg.itertuples(index=False)):
            src_lin = SOURCE_TO_LINEAGE[r.cell_type]
            row16 = sc16.loc[r.program_uid].dropna()
            row6 = sc6.loc[r.program_uid].dropna()
            top1_16 = str(row16.idxmax()) if len(row16) else ""
            top1_6 = str(row6.idxmax()) if len(row6) else ""
            src_sc16 = float(row16.get(src_lin, np.nan))
            src_sc6 = float(row6.get(src_lin, np.nan))
            others6 = row6.drop(labels=[src_lin], errors="ignore")
            margin6 = src_sc6 - float(others6.max()) if len(others6) else np.nan
            others16 = row16.drop(labels=[src_lin], errors="ignore")
            margin16 = src_sc16 - float(others16.max()) if len(others16) else np.nan
            rank16 = int((row16 > src_sc16).sum() + 1) if np.isfinite(src_sc16) else np.nan

            rows.append({
                "arm": arm_name,
                "program_uid": r.program_uid,
                "source_cell_type": r.cell_type,
                "source_lineage": src_lin,
                "module": r.module,
                "module_name": r.module_name,
                "primary_selected": r.primary_selected,
                "robust_display": r.robust_display,
                "n_active_lineages_16": len(row16),
                "top1_lineage_16": top1_16,
                "diagonal_match_16": top1_16 == src_lin,
                "source_score_z16": src_sc16,
                "margin_over_next_16": margin16,
                "rank_of_source_in_16": rank16,
                "top1_lineage_6": top1_6,
                "diagonal_match_6": top1_6 == src_lin,
                "source_score_z6": src_sc6,
                "margin_over_next_6": margin6,
                "boot_match_rate_6": float(boot["prog_boot_rate_6"][idx_p]) if boot is not None else np.nan,
                "boot_match_rate_16": float(boot["prog_boot_rate_16"][idx_p]) if boot is not None else np.nan,
            })

            if arm_name in {"excl_GSE136103_hm_min10", "incl_GSE136103_hm_min10"}:
                for tgt_lin, val in row16.items():
                    long_16_rows.append({
                        "arm": arm_name,
                        "program_uid": r.program_uid,
                        "source_cell_type": r.cell_type,
                        "source_lineage": src_lin,
                        "module": r.module,
                        "module_name": r.module_name,
                        "target_lineage": tgt_lin,
                        "n_donors_target": dcounts16[tgt_lin],
                        "score_z16": float(val),
                        "is_source_lineage": tgt_lin == src_lin,
                    })

        prog_df = pd.DataFrame(rows)
        prog_summary_tables[arm_name] = prog_df
        n_match_6 = int(prog_df["diagonal_match_6"].sum())
        n_match_16 = int(prog_df["diagonal_match_16"].sum())
        n_act16 = int(prog_df["n_active_lineages_16"].iloc[0])
        match_rows.append({
            "arm": arm_name,
            "prep_mode": prep_mode,
            "min_cells": min_cells,
            "n_active_lineages_16": n_act16,
            "active_lineages_16": ";".join(f"{k}({dcounts16[k]})" for k in sorted(sc16.columns)),
            "dropped_lineages_16": ";".join(f"{k}({dcounts16[k]})" for k in sorted(dcounts16) if k not in set(sc16.columns)),
            "n_programs": len(prog_df),
            "diagonal_match_6": n_match_6,
            "diagonal_frac_6": n_match_6 / float(len(prog_df)),
            "diagonal_match_16": n_match_16,
            "diagonal_frac_16": n_match_16 / float(len(prog_df)),
            "boot_med_frac_6": boot["boot_med_6"] if boot else np.nan,
            "boot_ci_low_6": boot["boot_ci_low_6"] if boot else np.nan,
            "boot_ci_high_6": boot["boot_ci_high_6"] if boot else np.nan,
            "boot_med_frac_16": boot["boot_med_16"] if boot else np.nan,
            "boot_ci_low_16": boot["boot_ci_low_16"] if boot else np.nan,
            "boot_ci_high_16": boot["boot_ci_high_16"] if boot else np.nan,
        })

    # Write Fig 4C program specificity outputs
    pd.DataFrame(match_rows).to_csv(OUT_DIR / "fig4c_program_match_counts_bootstrap.tsv", sep="\t", index=False)
    pd.concat(prog_summary_tables.values(), ignore_index=True).to_csv(
        OUT_DIR / "fig4c_program_specificity_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(long_16_rows).to_csv(OUT_DIR / "fig4c_program_specificity_16lineages.tsv", sep="\t", index=False)

    # Compute gene-level 16-lineage specificity table for primary excl_GSE136103_hm_min10
    assert primary_mean16 is not None and primary_z16 is not None
    active_cts = list(primary_mean16.columns if False else primary_mean16.index)
    # Build (genes x lineages) mean log2(CPM+1) and across-lineage Z
    m_mat = primary_mean16.T  # rows = genes, cols = active_cts
    z_mat = pd.DataFrame(np.nan, index=m_mat.index, columns=m_mat.columns)
    z_sub = primary_z16.T
    z_mat.loc[z_sub.index, z_sub.columns] = z_sub

    arr = m_mat.to_numpy(dtype=np.float64)
    row_max = arr.max(axis=1)
    tau = np.zeros(len(m_mat), dtype=np.float64)
    pos = row_max > 0
    tau[pos] = np.sum(1.0 - (arr[pos] / row_max[pos, None]), axis=1) / float(arr.shape[1] - 1)

    z_arr = z_mat.to_numpy(dtype=np.float64)
    all_nan = np.all(~np.isfinite(z_arr), axis=1)
    safe_z = np.where(np.isfinite(z_arr), z_arr, -np.inf)
    order = np.argsort(-safe_z, axis=1)
    top1_idx = order[:, 0]
    top2_idx = order[:, 1]
    top1_lin = [active_cts[j] if not all_nan[i] else "" for i, j in enumerate(top1_idx)]
    top2_lin = [active_cts[j] if not all_nan[i] else "" for i, j in enumerate(top2_idx)]
    top1_z = [float(z_arr[i, j]) if not all_nan[i] else np.nan for i, j in enumerate(top1_idx)]
    top2_z = [float(z_arr[i, j]) if not all_nan[i] else np.nan for i, j in enumerate(top2_idx)]

    gene_spec = pd.DataFrame({
        "gene": m_mat.index.astype(str),
        "max_mean_log2cpm1": row_max,
        "tau_16lineages": tau,
        "top1_lineage": top1_lin,
        "top1_z": top1_z,
        "top2_lineage": top2_lin,
        "top2_z": top2_z,
        "margin_z": np.asarray(top1_z) - np.asarray(top2_z),
    })
    for ct in active_cts:
        gene_spec[f"mean_{ct}"] = m_mat[ct].to_numpy()
        gene_spec[f"z_{ct}"] = z_mat[ct].to_numpy()

    gene_spec.to_csv(OUT_DIR / "fig4c_gene_specificity_16lineages.tsv.gz", sep="\t", index=False, compression="gzip")

    # Save curated Fig 4 / marker + top-hub subset in plain TSV for direct inspection
    top_hubs = mem.sort_values("original_l1_weight", ascending=False).groupby("program_uid").head(1)["mapped_symbol"].astype(str).tolist()
    focus = set(MARKERS_AND_FIG4_GENES) | set(top_hubs)
    gene_spec[gene_spec["gene"].isin(focus)].to_csv(
        OUT_DIR / "fig4c_marker_and_fig4_genes_16lineages.tsv", sep="\t", index=False
    )

    # Also save lineage donor counts census across all 16 lineages
    census_rows = []
    for arm_label, mfn, pmode in [
        ("excl_GSE136103_prep_eligible", lambda d: (~to_bool(d["exclude"])) & to_bool(d["binary_eligible"]) & d["condition_harmonized"].isin(["Healthy", "MASLD"]), "prep_eligible"),
        ("incl_GSE136103_prep_eligible", lambda d: (d["dataset"] != "Liver_Atlas") & to_bool(d["binary_eligible"]) & d["condition_harmonized"].isin(["Healthy", "MASLD"]), "prep_eligible"),
        ("excl_GSE136103_all_preps", lambda d: (~to_bool(d["exclude"])) & to_bool(d["binary_eligible"]) & d["condition_harmonized"].isin(["Healthy", "MASLD"]), "all_preps"),
    ]:
        for ct in sorted({c for c, _ in lineage_donor_info}):
            info = lineage_donor_info[(ct, pmode)]
            sub = info[mfn(info)]
            census_rows.append({
                "arm": arm_label,
                "lineage": ct,
                "donors_ge10": int((sub["n_lineage_cells"] >= 10).sum()),
                "donors_ge30": int((sub["n_lineage_cells"] >= 30).sum()),
                "donors_ge50": int((sub["n_lineage_cells"] >= 50).sum()),
                "healthy_ge50": int(((sub["n_lineage_cells"] >= 50) & (sub["condition_harmonized"] == "Healthy")).sum()),
                "masld_ge50": int(((sub["n_lineage_cells"] >= 50) & (sub["condition_harmonized"] == "MASLD")).sum()),
                "total_cells_ge10": int(sub.loc[sub["n_lineage_cells"] >= 10, "n_lineage_cells"].sum()),
            })
    pd.DataFrame(census_rows).to_csv(OUT_DIR / "fig4c_lineage_donor_census_16lineages.tsv", sep="\t", index=False)
    print("DONE Fig 4C 16-lineage specificity", flush=True)


if __name__ == "__main__":
    main()
