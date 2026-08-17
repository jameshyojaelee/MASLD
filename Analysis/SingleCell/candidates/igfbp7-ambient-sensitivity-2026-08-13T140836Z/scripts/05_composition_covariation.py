#!/usr/bin/env python
"""Correction-independent lineage-covariation test for the IGFBP7 SENSITIVITY ANALYSIS.

SCOPE (binding): sensitivity analysis only; the frozen programs are read-only.
Nothing is rediscovered, refit, reweighted, renamed or re-selected.

Ambient contamination of hepatocyte barcodes by another lineage's transcripts is
proportional to how much of THAT LINEAGE'S material was in the same droplet
suspension. So an ambient-driven hepatocyte program should track its OWN putative
source lineage's donor fraction.

The discriminating comparison is therefore paired: each program is correlated
with its own source-lineage fraction, and the frozen registry's two unambiguous
non-hepatocyte contamination programs supply the benchmark magnitude:
  hep module 1  Leukocyte / immune (DOCK2)     -> leukocyte fraction
  hep module 2  Sinusoidal endothelial (STAB2) -> endothelial fraction
The question is whether IGFBP7-vs-fibroblast reaches the magnitude that known
contamination programs reach against their own sources. Correlating every
program against fibroblasts alone does not discriminate and is not used here.

Every program is additionally correlated against ALL major lineage fractions, so
that a program tracking one lineage specifically can be told apart from a
program tracking generic non-hepatocyte content.

Uses STORED frozen donor scores only; no ambient correction enters this script.
The fibroblast-fraction-adjusted stage coefficient is DESCRIPTIVE ONLY:
fibroblast expansion is part of MASH biology, so conditioning on it is adjusting
for a mediator, not producing a corrected stage estimate.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(os.environ["MASLD_PROJECT_ROOT"])
CAND = Path(os.environ["CAND_ROOT"])
WORK = CAND / "work"
RES = CAND / "results"

HS_ROOT = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
DONOR_SCORES = HS_ROOT / "donor_collapse/donor_scores_all_weighted.tsv"
PROGRAM_ROOT = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / "program-context-v2-candidate-2026-08-07/hotspot"
)
META = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
PAIRINGS = {
    "GSE244832": ROOT / "data/GSE244832/metadata/donor_pairing.csv",
    "GSE185477": ROOT / "data/GSE185477/metadata/donor_pairing.csv",
    "GSE202379": ROOT / "data/GSE202379/metadata/donor_pairing.csv",
    "GSE136103": ROOT / "data/GSE136103/metadata/donor_pairing.csv",
}
PRIMARY_STAGES = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2}

LEUKOCYTE = [
    "T cells", "Macrophages", "Mono+mono derived cells", "Circulating NK/NKT",
    "Resident NK", "Neutrophils", "B cells", "cDC1s", "cDC2s", "pDCs",
    "Plasma cells", "Basophils",
]
# lineage aggregates measured against every program
LINEAGES = {
    "hepatocyte": ["Hepatocytes"],
    "fibroblast": ["Fibroblasts"],
    "cholangiocyte": ["Cholangiocytes"],
    "endothelial": ["Endothelial cells"],
    "leukocyte": LEUKOCYTE,
    "macrophage": ["Macrophages"],
    "tcell": ["T cells"],
    "non_hepatocyte": None,  # 1 - hepatocyte fraction
}
# frozen-registry programs whose putative source lineage is unambiguous
SOURCE_LINEAGE = {
    1: ("leukocyte", "known non-hepatocyte contamination benchmark"),
    2: ("endothelial", "known non-hepatocyte contamination benchmark"),
    8: ("fibroblast", "program under test"),
    20: ("cholangiocyte", "paired comparator"),
    12: ("hepatocyte", "hepatocyte-intrinsic control"),
    21: ("hepatocyte", "hepatocyte-intrinsic control"),
    16: ("hepatocyte", "hepatocyte-intrinsic control"),
}


def donor_map() -> dict[str, str]:
    out: dict[str, str] = {}
    for dataset, path in PAIRINGS.items():
        table = pd.read_csv(path, dtype=str)
        for row in table.itertuples(index=False):
            donor = f"{dataset}_{row.donor_id}"
            for s in str(row.rna_srrs).split(";"):
                if s.strip():
                    out[s.strip()] = donor
    return out


def donor_metadata(mapping: dict[str, str]) -> pd.DataFrame:
    t = pd.read_csv(META, sep="\t", dtype=str)[
        ["sample", "dataset", "disease_stage_coarse", "exclude_stage_analysis"]
    ]
    t["donor"] = t["sample"].map(mapping).fillna(t["sample"])
    t["exclude"] = t["exclude_stage_analysis"].str.lower().isin({"true", "t", "1", "yes"})
    rows = []
    for donor, g in t.groupby("donor", sort=False):
        ds = sorted(set(g["dataset"].dropna()) - {""})
        st = sorted(set(g["disease_stage_coarse"].dropna()) - {""})
        rows.append(
            {
                "donor": donor,
                "dataset": ds[0] if ds else None,
                "stage": st[0] if st else None,
                "exclude": bool(g["exclude"].any()),
                "n_sequencing_samples": int(len(g)),
            }
        )
    return pd.DataFrame(rows)


def ols_hc3(y: np.ndarray, X: np.ndarray) -> dict:
    XtXi = np.linalg.inv(X.T @ X)
    b = XtXi @ X.T @ y
    e = y - X @ b
    n, k = X.shape
    df = n - k
    se = np.sqrt(float(e @ e) / df * np.diag(XtXi))
    h = np.clip(np.diag(X @ XtXi @ X.T), 0, 1 - 1e-12)
    V = XtXi @ (X.T @ (((e / (1 - h)) ** 2)[:, None] * X)) @ XtXi
    se3 = np.sqrt(np.diag(V))
    return {"beta": b, "se": se, "hc3_se": se3,
            "p": 2 * stats.t.sf(np.abs(b / se), df),
            "hc3_p": 2 * stats.t.sf(np.abs(b / se3), df), "df": df, "n": n}


def main() -> None:
    manifest = json.loads((WORK / "manifest.json").read_text(encoding="utf-8"))
    frames = []
    for dataset in sorted(manifest["datasets"]):
        m = pd.read_csv(WORK / dataset / "meta.csv.gz", usecols=["cell_type", "sample"])
        m["dataset"] = dataset
        frames.append(m)
    cells = pd.concat(frames, ignore_index=True)

    mapping = donor_map()
    cells["donor"] = cells["sample"].astype(str).map(mapping).fillna(cells["sample"].astype(str))
    comp = cells.groupby(["donor", "cell_type"]).size().unstack(fill_value=0).astype(float)
    total = comp.sum(axis=1)

    lin = pd.DataFrame(index=comp.index)
    for name, members in LINEAGES.items():
        if name == "non_hepatocyte":
            lin[name] = 1 - comp.get("Hepatocytes", 0) / total
        else:
            present = [c for c in members if c in comp.columns]
            lin[name] = comp[present].sum(axis=1) / total if present else np.nan
    lin["n_cells_total"] = total
    lin = lin.reset_index()
    lin.to_csv(RES / "05_donor_lineage_fractions.tsv", sep="\t", index=False)

    scores = pd.read_csv(DONOR_SCORES, sep="\t").rename(columns={"sample": "donor"})
    registry = pd.read_csv(PROGRAM_ROOT / "program_registry_v2.tsv", sep="\t")
    meta = donor_metadata(mapping)

    # ---- donor census for the frozen primary model ----------------------
    fm = meta[(~meta["exclude"]) & meta["stage"].isin(PRIMARY_STAGES)]
    census = pd.crosstab(fm["dataset"], fm["stage"]).reindex(
        columns=list(PRIMARY_STAGES), fill_value=0
    )
    census["total_donors"] = census.sum(axis=1)
    census.loc["ALL"] = census.sum(axis=0)
    census.to_csv(RES / "05_donor_census_primary_model.tsv", sep="\t")
    print("=== donor census entering the frozen primary stage model ===", flush=True)
    print(census.to_string(), flush=True)
    print(
        f"\nsequencing samples in donor metadata: {int(meta['n_sequencing_samples'].sum())}; "
        f"biological donors after collapse: {len(meta)}; "
        f"donors flagged exclude_stage_analysis: {int(meta['exclude'].sum())}; "
        f"donors in the primary model: {len(fm)}\n",
        flush=True,
    )

    long_rows, paired_rows = [], []
    hep_modules = registry[registry["cell_type"] == "hepatocytes"]
    for r in hep_modules.itertuples(index=False):
        s = scores[scores["module"] == f"hepatocytes__{r.module}"][["donor", "score"]]
        f = (
            s.merge(meta, on="donor", how="inner", validate="one_to_one")
            .merge(lin, on="donor", how="inner", validate="one_to_one")
        )
        f = f[(~f["exclude"]) & f["stage"].isin(PRIMARY_STAGES)].copy()
        if len(f) < 10:
            continue
        f["stage_ordinal"] = f["stage"].map(PRIMARY_STAGES).astype(float)
        y = f["score"].to_numpy(float)
        D = pd.get_dummies(f["dataset"], drop_first=True, dtype=float).to_numpy(float)
        ones = np.ones(len(f))
        # residualise on dataset so associations are within-dataset
        Xd = np.column_stack([ones, D])
        Pd = Xd @ np.linalg.pinv(Xd)
        y_r = y - Pd @ y

        per_lineage = {}
        for name in LINEAGES:
            v = f[name].to_numpy(float)
            v_r = v - Pd @ v
            per_lineage[name] = {
                "pearson_within_dataset": float(np.corrcoef(y_r, v_r)[0, 1]),
                "spearman_within_dataset": float(stats.spearmanr(y_r, v_r).statistic),
                "pearson_marginal": float(np.corrcoef(y, v)[0, 1]),
            }
            long_rows.append(
                {
                    "module": r.module, "module_name": r.module_name,
                    "robust_display": r.robust_display, "lineage": name,
                    "n_donors": len(f), **per_lineage[name],
                }
            )

        base = ols_hc3(y, np.column_stack([ones, f["stage_ordinal"].to_numpy(float), D]))
        src, role = SOURCE_LINEAGE.get(r.module, (None, ""))
        row = {
            "module": r.module, "module_name": r.module_name,
            "robust_display": r.robust_display, "n_donors": len(f),
            "source_lineage": src, "calibrator_role": role,
            "stage_beta_unadjusted": float(base["beta"][1]),
            "stage_hc3_p_unadjusted": float(base["hc3_p"][1]),
        }
        if src is not None:
            row["spearman_vs_own_source"] = per_lineage[src]["spearman_within_dataset"]
            row["pearson_vs_own_source"] = per_lineage[src]["pearson_within_dataset"]
            v = f[src].to_numpy(float)
            adj = ols_hc3(
                y, np.column_stack([ones, f["stage_ordinal"].to_numpy(float), v, D])
            )
            row["stage_beta_source_adjusted"] = float(adj["beta"][1])
            row["stage_hc3_p_source_adjusted"] = float(adj["hc3_p"][1])
            row["pct_stage_beta_retained"] = float(100 * adj["beta"][1] / base["beta"][1])
        # strongest lineage association regardless of assignment
        best = max(
            LINEAGES, key=lambda n: abs(per_lineage[n]["spearman_within_dataset"])
        )
        row["strongest_lineage"] = best
        row["strongest_lineage_spearman"] = per_lineage[best]["spearman_within_dataset"]
        row["spearman_vs_fibroblast"] = per_lineage["fibroblast"]["spearman_within_dataset"]
        row["spearman_vs_non_hepatocyte"] = per_lineage["non_hepatocyte"]["spearman_within_dataset"]
        row["fibroblast_minus_non_hepatocyte"] = (
            row["spearman_vs_fibroblast"] - row["spearman_vs_non_hepatocyte"]
        )
        paired_rows.append(row)

    long = pd.DataFrame(long_rows)
    long.to_csv(RES / "05_program_by_lineage_matrix.tsv", sep="\t", index=False)
    paired = pd.DataFrame(paired_rows)
    paired.to_csv(RES / "05_own_source_calibration.tsv", sep="\t", index=False)

    print("=== each program vs ITS OWN putative source lineage (the valid calibration) ===", flush=True)
    cal = paired[paired["source_lineage"].notna()].copy()
    cal = cal.sort_values("spearman_vs_own_source", ascending=False)
    print(
        cal[
            ["module", "module_name", "source_lineage", "calibrator_role",
             "spearman_vs_own_source", "stage_beta_unadjusted",
             "stage_beta_source_adjusted", "pct_stage_beta_retained",
             "stage_hc3_p_source_adjusted"]
        ].to_string(index=False),
        flush=True,
    )

    print("\n=== program x lineage Spearman matrix (within-dataset) ===", flush=True)
    mat = long.pivot_table(
        index=["module", "module_name"], columns="lineage",
        values="spearman_within_dataset",
    )
    mat = mat[list(LINEAGES)]
    print(mat.round(3).to_string(), flush=True)

    print("\n=== is IGFBP7 fibroblast-specific or just non-hepatocyte-generic? ===", flush=True)
    print(
        paired.sort_values("fibroblast_minus_non_hepatocyte", ascending=False)[
            ["module", "module_name", "spearman_vs_fibroblast",
             "spearman_vs_non_hepatocyte", "fibroblast_minus_non_hepatocyte",
             "strongest_lineage", "strongest_lineage_spearman"]
        ].head(10).to_string(index=False),
        flush=True,
    )
    print("[composition] DONE", flush=True)


if __name__ == "__main__":
    main()
