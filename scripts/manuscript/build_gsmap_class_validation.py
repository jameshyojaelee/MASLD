#!/usr/bin/env python3
"""Test frozen evidence classes with gsMap's native continuous PCC statistic.

The primary comparison is convergent versus genetic-only genes, matched on
bulk expression, gene biotype, and median spatial gene-specificity score. The
legacy thresholded risk-gene derivation is not used.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
RELEASE_ID = os.environ.get("MANUSCRIPT_RELEASE_ID", "2026-07-15-r2")
OUT = ROOT / "RNA-seq/results/manuscript_release" / RELEASE_ID
PCC_FILE = ROOT / "Analysis/Spatial/results/gsmap/gsmap_pcc_by_trait.csv"
N_PERM = 10000
SEED = 42
rng = np.random.default_rng(SEED)

ENZYME_TRAITS = {"ukbb_alt", "ukbb_ast", "ukbb_ggt", "mvp_alt", "mvp_ast"}


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def matched_pool(
    data: pd.DataFrame,
    comparator: str,
    controls_per_case: int = 5,
) -> pd.DataFrame:
    cases = data[data["primary_evidence_class"].eq("convergent")].copy()
    controls = data[data["primary_evidence_class"].eq(comparator)].copy()
    selected_controls: list[int] = []
    selected_cases: list[int] = []
    used: set[int] = set()
    for idx, case in cases.sort_values("gene").iterrows():
        pool = controls[
            (controls["protein_coding"] == case["protein_coding"])
            & (controls["expr_bin"] == case["expr_bin"])
            & (controls["gss_bin"] == case["gss_bin"])
            & (~controls.index.isin(used))
        ].copy()
        if len(pool) < controls_per_case:
            pool = controls[
                (controls["protein_coding"] == case["protein_coding"])
                & (~controls.index.isin(used))
            ].copy()
        if pool.empty:
            continue
        pool["distance"] = (
            (pool["bulk_AveExpr"] - case["bulk_AveExpr"]).abs()
            + (pool["median_gss_mean"] - case["median_gss_mean"]).abs()
        )
        chosen = pool.sort_values(["distance", "gene"]).head(controls_per_case)
        selected_cases.append(idx)
        selected_controls.extend(chosen.index.tolist())
        used.update(chosen.index.tolist())
    return data.loc[selected_cases + selected_controls].copy()


def permutation_test(pool: pd.DataFrame, comparator: str) -> tuple[float, float, int, int]:
    label = pool["primary_evidence_class"].eq("convergent").to_numpy()
    values = pool["pcc_mean"].to_numpy(float)
    n_case = int(label.sum())
    n_control = int((~label).sum())
    if n_case < 3 or n_control < 3:
        return np.nan, np.nan, n_case, n_control
    observed = float(values[label].mean() - values[~label].mean())
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        case_idx = rng.choice(len(values), size=n_case, replace=False)
        perm_label = np.zeros(len(values), dtype=bool)
        perm_label[case_idx] = True
        null[i] = values[perm_label].mean() - values[~perm_label].mean()
    p_value = float((np.sum(null >= observed) + 1) / (N_PERM + 1))
    return observed, p_value, n_case, n_control


def main() -> None:
    pcc = pd.read_csv(PCC_FILE, low_memory=False)
    required = {
        "trait", "cohort", "gene", "pcc_mean", "median_gss_mean",
        "primary_evidence_class", "bulk_AveExpr", "gene_biotype",
    }
    missing = required.difference(pcc.columns)
    if missing:
        raise RuntimeError(
            f"PCC table lacks frozen class columns {sorted(missing)}; "
            "rerun Analysis/Spatial/scripts/15l_gsmap_pcc_narrowing.py"
        )

    pcc = pcc[pcc["primary_evidence_class"].isin(
        ["convergent", "genetic_only", "disease_state_only"]
    )].copy()
    pcc = pcc[
        np.isfinite(pcc["pcc_mean"])
        & np.isfinite(pcc["bulk_AveExpr"])
        & np.isfinite(pcc["median_gss_mean"])
    ].copy()
    pcc["protein_coding"] = pcc["gene_biotype"].eq("protein_coding")
    pcc["expr_bin"] = pcc.groupby(["trait", "cohort"])["bulk_AveExpr"].transform(
        lambda x: pd.qcut(x.rank(method="first"), 5, labels=False, duplicates="drop")
    )
    pcc["gss_bin"] = pcc.groupby(["trait", "cohort"])["median_gss_mean"].transform(
        lambda x: pd.qcut(x.rank(method="first"), 5, labels=False, duplicates="drop")
    )
    pcc["trait_scope"] = np.where(pcc["trait"].isin(ENZYME_TRAITS), "enzyme", "direct_disease")

    rows = []
    matched_frames = []
    for (trait, cohort), group in pcc.groupby(["trait", "cohort"], sort=True):
        for comparator in ["genetic_only", "disease_state_only"]:
            pool = matched_pool(group, comparator)
            effect, p_value, n_case, n_control = permutation_test(pool, comparator)
            rows.append({
                "analysis_release_id": RELEASE_ID,
                "trait": trait,
                "trait_scope": "enzyme" if trait in ENZYME_TRAITS else "direct_disease",
                "cohort": cohort,
                "comparison": f"convergent_vs_{comparator}",
                "n_convergent": n_case,
                "n_comparator": n_control,
                "mean_pcc_difference": effect,
                "permutation_p": p_value,
            })
            if not pool.empty:
                pool = pool.copy()
                pool["comparison"] = f"convergent_vs_{comparator}"
                matched_frames.append(pool)

    results = pd.DataFrame(rows)
    finite = results["permutation_p"].notna()
    results["q_value"] = np.nan
    results.loc[finite, "q_value"] = bh(results.loc[finite, "permutation_p"])
    results["positive_significant"] = (
        results["mean_pcc_difference"].gt(0) & results["q_value"].lt(0.05)
    )
    results.to_csv(OUT / "gsmap_class_validation.tsv", sep="\t", index=False)
    if matched_frames:
        pd.concat(matched_frames, ignore_index=True).to_csv(
            OUT / "gsmap_class_matched_genes.tsv", sep="\t", index=False
        )

    primary = results[
        results["trait_scope"].eq("direct_disease")
        & results["comparison"].eq("convergent_vs_genetic_only")
    ]
    replication = (
        primary.groupby("trait")
        .agg(
            n_cohorts=("cohort", "nunique"),
            n_positive_significant=("positive_significant", "sum"),
            min_effect=("mean_pcc_difference", "min"),
            max_q=("q_value", "max"),
        )
        .reset_index()
    )
    replication["replicated_both_cohorts"] = (
        replication["n_cohorts"].ge(2)
        & replication["n_positive_significant"].eq(replication["n_cohorts"])
    )
    replication.to_csv(OUT / "gsmap_trait_replication.tsv", sep="\t", index=False)

    n_replicated = int(replication["replicated_both_cohorts"].sum())
    specificity_controls_available = False
    headline_pass = bool(n_replicated >= 2 and specificity_controls_available)
    gates_path = OUT / "acceptance_gates.tsv"
    gates = pd.read_csv(gates_path, sep="\t")
    gates = gates[gates["gate"] != "gsmap_direct_disease_article"]
    gates = pd.concat([
        gates,
        pd.DataFrame([{
            "analysis_release_id": RELEASE_ID,
            "gate": "gsmap_direct_disease_article",
            "passed": headline_pass,
            "criterion": (
                ">=2 direct-disease traits show positive matched PCC effects at q<0.05 "
                "in both cohorts, with prespecified non-liver specificity controls"
            ),
            "detail": (
                f"replicated_direct_traits={n_replicated}; "
                f"specificity_controls_available={specificity_controls_available}"
            ),
        }]),
    ], ignore_index=True)
    gates.to_csv(gates_path, sep="\t", index=False)

    ledger_path = ROOT / "docs/manuscript/release/claim_ledger.tsv"
    ledger = pd.read_csv(ledger_path, sep="\t")
    hit = ledger["claim_id"].eq("C_GSMAP_RETIRED")
    ledger.loc[hit, "status"] = "failed_acceptance_gate" if not headline_pass else "supported"
    ledger.loc[hit, "allowed_wording"] = (
        "Native gsMap PCC was tested against frozen evidence classes; the prespecified "
        f"headline gate did not pass ({n_replicated} direct traits replicated in both cohorts; "
        "non-liver specificity controls unavailable)."
        if not headline_pass else
        f"Native gsMap PCC replicated class enrichment for {n_replicated} direct traits in both cohorts."
    )
    ledger.loc[hit, "source_output"] = "gsmap_class_validation.tsv; gsmap_trait_replication.tsv"
    ledger.to_csv(ledger_path, sep="\t", index=False)

    print(results.to_string(index=False))
    print("\nReplication:\n", replication.to_string(index=False))
    print(f"\nheadline_pass={headline_pass}")


if __name__ == "__main__":
    main()
