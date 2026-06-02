#!/usr/bin/env python
"""S4-1: Compositional retest of cell-type proportion shifts with scCODA.

Inputs (read from project, NOT worktree):
  - Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv
        per-cell hepatocyte subtype assignments (incl meta-subtype label).
  - Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/subtype_disease_enrichment.csv
        original Fisher per-cluster p-values (43 raw subtypes).
  - Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/meta_subtype_mapping.csv
        raw-subtype -> meta-subtype (Healthy/Progressor/Moderate/Stable/Neutral/Disease-Neutral)
  - Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv

For each proportion claim we contrast Healthy vs MASLD donors (binary
covariate) with scCODA's Bayesian Dirichlet-multinomial. Donor IS the
replicate; cells per donor are the multinomial draw. scCODA does not
fit a donor random effect explicitly (the Dirichlet-multinomial
likelihood already captures donor overdispersion) but the per-donor row
structure preserves donor-level variance.

Output:
  Analysis/SingleCell/results_gpu_v2/proportion_compositional/sccoda_results.tsv
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SC = PROJECT_ROOT / "Analysis/SingleCell"
HEP = SC / "results_gpu_v2/hepatocyte_subtypes"
DONOR_META_PATH = SC / "results_gpu_v2/mcp/inputs/donor_metadata.tsv"

OUT = SC / "results_gpu_v2/proportion_compositional"
OUT.mkdir(parents=True, exist_ok=True)


def per_donor_counts(meta_cell: pd.DataFrame, label_col: str) -> pd.DataFrame:
    """Wide donor x label cell-count matrix."""
    counts = (
        meta_cell.groupby(["sample", label_col])
        .size()
        .unstack(fill_value=0)
    )
    return counts


def run_sccoda(counts: pd.DataFrame, donor_cov: pd.DataFrame,
               formula: str, reference_cell_type: str | None = None,
               n_iter: int = 20000):
    """Fit scCODA for one binary covariate.

    counts:    donor x label cell counts (rows = donors).
    donor_cov: DataFrame indexed by sample with covariate columns.
    formula:   patsy formula, e.g. "C(condition, Treatment('Healthy'))".
    Returns per-label DataFrame with effect / inclusion-probability / Final-Parameter.
    """
    from sccoda.util import comp_ana as mod
    import anndata as ad

    # ensure aligned
    shared = counts.index.intersection(donor_cov.index)
    counts = counts.loc[shared]
    donor_cov = donor_cov.loc[shared]

    # drop labels with zero variance (all-zero across donors)
    counts = counts.loc[:, counts.sum(axis=0) > 0]

    if reference_cell_type is None:
        # auto-pick reference = most abundant label observed in >=80% of donors
        nz = (counts > 0).mean(axis=0)
        candidates = counts.loc[:, nz >= 0.8]
        if candidates.shape[1] == 0:
            reference_cell_type = "automatic"
        else:
            reference_cell_type = candidates.sum(axis=0).idxmax()

    adata = ad.AnnData(
        X=counts.values.astype(float),
        obs=donor_cov.copy(),
        var=pd.DataFrame(index=counts.columns.astype(str)),
    )
    adata.obs.index = counts.index.astype(str)

    print(f"   scCODA: {adata.n_obs} donors x {adata.n_vars} labels; reference = {reference_cell_type!r}; formula = {formula}")
    model = mod.CompositionalAnalysis(
        adata,
        formula=formula,
        reference_cell_type=reference_cell_type,
    )
    try:
        result = model.sample_hmc(num_results=n_iter)
    except Exception as exc:  # pragma: no cover
        print(f"      HMC failed: {exc}; retrying NUTS")
        result = model.sample_nuts(num_results=n_iter // 2)
    eff = result.effect_df.copy().reset_index()
    eff.columns = [c if isinstance(c, str) else "_".join(map(str, c)) for c in eff.columns]
    return eff, reference_cell_type


def main():
    # 1. Load per-cell metadata (43 raw subtypes + meta-subtype label)
    meta_path = HEP / "hepatocyte_subtype_metadata.csv"
    cell = pd.read_csv(meta_path, index_col=0)
    print(f"[load] hepatocyte cells: {len(cell):,}  donors: {cell['sample'].nunique()}")

    mapping = pd.read_csv(HEP / "meta_subtype_mapping.csv")
    raw_to_meta = dict(zip(mapping["subtype"].astype(str), mapping["meta_subtype"].astype(str)))
    cell["meta_subtype"] = cell["hepatocyte_subtype"].astype(str).map(raw_to_meta).fillna("Other")

    # 2. Donor covariate
    donor = pd.read_csv(DONOR_META_PATH, sep="\t")
    cond_col = "condition_harmonized" if "condition_harmonized" in donor.columns else "condition"
    donor["condition_binary"] = donor[cond_col].map(
        lambda x: "Healthy" if str(x).strip().lower() in {"healthy", "normal", "control"} else "MASLD"
    )
    donor_cov = donor.set_index("sample")[["condition_binary"]]

    # 3. Original Fisher p-values
    fisher_orig = pd.read_csv(HEP / "subtype_disease_enrichment.csv")
    fisher_orig["subtype"] = fisher_orig["subtype"].astype(str)

    results: list[dict] = []

    # ---- Claim A: META-SUBTYPE proportions (Progressor 1.3% -> 63.1%) ----
    print("\n=== Claim A: meta-subtype proportions ===")
    counts_meta = per_donor_counts(cell, "meta_subtype")
    counts_meta.to_csv(OUT / "counts_meta_subtype_per_donor.tsv", sep="\t")
    sc_meta, ref_meta = run_sccoda(
        counts_meta, donor_cov,
        formula="C(condition_binary, Treatment('Healthy'))",
        reference_cell_type=None,
        n_iter=20000,
    )
    sc_meta.to_csv(OUT / "sccoda_meta_subtype.tsv", sep="\t", index=False)
    for _, row in sc_meta.iterrows():
        results.append({
            "claim": "meta_subtype_Healthy_vs_MASLD",
            "label": str(row.get("Cell Type", row.get("index", "?"))),
            "p_original_min": float(fisher_orig["fisher_pval"].min()),
            "sccoda_final_param": float(row.get("Final Parameter", np.nan)),
            "sccoda_inclusion_prob": float(row.get("Inclusion probability", np.nan)),
            "sccoda_credible": bool(row.get("Inclusion probability", 0) >= 0.95),
            "reference_label": ref_meta,
        })

    # ---- Claim B: RAW 43-subtype proportions (Fisher per-cluster) ----
    print("\n=== Claim B: raw 43-subtype proportions ===")
    counts_raw = per_donor_counts(cell.assign(subtype=cell["hepatocyte_subtype"].astype(str)), "subtype")
    counts_raw.to_csv(OUT / "counts_raw_subtype_per_donor.tsv", sep="\t")
    sc_raw, ref_raw = run_sccoda(
        counts_raw, donor_cov,
        formula="C(condition_binary, Treatment('Healthy'))",
        reference_cell_type=None,
        n_iter=20000,
    )
    sc_raw.to_csv(OUT / "sccoda_raw_subtype.tsv", sep="\t", index=False)
    fisher_map = dict(zip(fisher_orig["subtype"], fisher_orig["fisher_pval"]))
    for _, row in sc_raw.iterrows():
        lab = str(row.get("Cell Type", row.get("index", "?")))
        results.append({
            "claim": f"raw_subtype_{lab}_Healthy_vs_MASLD",
            "label": lab,
            "p_original_min": float(fisher_map.get(lab, np.nan)),
            "sccoda_final_param": float(row.get("Final Parameter", np.nan)),
            "sccoda_inclusion_prob": float(row.get("Inclusion probability", np.nan)),
            "sccoda_credible": bool(row.get("Inclusion probability", 0) >= 0.95),
            "reference_label": ref_raw,
        })

    out_df = pd.DataFrame(results)
    out_df.to_csv(OUT / "sccoda_results.tsv", sep="\t", index=False)
    print(f"\n[done] wrote {OUT/'sccoda_results.tsv'}  ({len(out_df)} rows)")


if __name__ == "__main__":
    main()
