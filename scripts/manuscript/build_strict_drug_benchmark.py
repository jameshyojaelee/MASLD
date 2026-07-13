#!/usr/bin/env python3
"""Build an auditable MASLD clinical-target benchmark for the manuscript.

Only curated approvals, curated active programs, and MASLD-specific
ClinicalTrials.gov records define positive clinical evidence. Open Targets
animal-model and literature evidence is retained as ``association_only`` and is
never treated as a drug-development stage.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
RELEASE_ID = os.environ.get("MANUSCRIPT_RELEASE_ID", "2026-07-10-r1")
OUT = ROOT / "RNA-seq/results/manuscript_release" / RELEASE_ID
OUT.mkdir(parents=True, exist_ok=True)
SEED = 42
N_BOOT = 5000
rng = np.random.default_rng(SEED)


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.upper().isin({"TRUE", "T", "1", "YES"})


def auc(score: np.ndarray, label: np.ndarray) -> float:
    ok = np.isfinite(score)
    score = score[ok]
    label = label[ok].astype(bool)
    n_pos = int(label.sum())
    n_neg = int((~label).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = stats.rankdata(score)
    u = ranks[label].sum() - n_pos * (n_pos + 1) / 2
    return float(u / (n_pos * n_neg))


def stratified_boot_diff(
    score_a: np.ndarray,
    score_b: np.ndarray,
    label: np.ndarray,
) -> tuple[float, float, float]:
    ok = np.isfinite(score_a) & np.isfinite(score_b)
    a, b, y = score_a[ok], score_b[ok], label[ok].astype(bool)
    pos = np.flatnonzero(y)
    neg = np.flatnonzero(~y)
    observed = auc(a, y) - auc(b, y)
    boot = np.empty(N_BOOT)
    for i in range(N_BOOT):
        idx = np.concatenate([
            rng.choice(pos, size=len(pos), replace=True),
            rng.choice(neg, size=len(neg), replace=True),
        ])
        boot[i] = auc(a[idx], y[idx]) - auc(b[idx], y[idx])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return float(observed), float(lo), float(hi)


def make_matches(df: pd.DataFrame, controls_per_positive: int = 5) -> pd.DataFrame:
    """Greedy exact-stratum matching with expression-distance tie breaking."""
    positives = df[df["strict_positive"]].copy()
    controls = df[df["eligible_control"]].copy()
    used: set[int] = set()
    selected: list[int] = []
    matched_positive: list[int] = []

    for pidx, p in positives.sort_values("symbol").iterrows():
        pool = controls[
            (controls["tdl_group"] == p["tdl_group"])
            & (controls["protein_coding"] == p["protein_coding"])
            & (controls["expr_bin"] == p["expr_bin"])
            & (~controls.index.isin(used))
        ].copy()
        if len(pool) < controls_per_positive:
            pool = controls[
                (controls["tdl_group"] == p["tdl_group"])
                & (controls["protein_coding"] == p["protein_coding"])
                & (~controls.index.isin(used))
            ].copy()
        if pool.empty:
            continue
        pool["distance"] = (pool["bulk_AveExpr"] - p["bulk_AveExpr"]).abs()
        chosen = pool.sort_values(["distance", "symbol"]).head(controls_per_positive)
        matched_positive.append(pidx)
        selected.extend(chosen.index.tolist())
        used.update(chosen.index.tolist())

    out = df.loc[matched_positive + selected].copy()
    out["matched_role"] = np.where(out["strict_positive"], "positive", "control")
    return out


def main() -> None:
    conv = pd.read_csv(
        ROOT / "RNA-seq/results/multi_evidence/convergence_evidence.csv",
        low_memory=False,
    ).rename(columns={"human_symbol": "symbol"})
    drug = pd.read_csv(
        ROOT / "data/external/drug_targets/drug_target_classification.tsv",
        sep="\t",
        low_memory=False,
    )
    release = pd.read_csv(OUT / "evidence_class_table.tsv", sep="\t", low_memory=False)

    curated = drug["curated_status"].fillna("").str.strip().str.lower()
    ct_trial = as_bool(drug["ct_has_masld_trial"])
    ct_phase = pd.to_numeric(drug["ct_masld_trial_max_phase"], errors="coerce")
    approved = curated.eq("approved")
    discontinued = curated.eq("discontinued")
    clinical = curated.eq("active") | (ct_trial & ct_phase.ge(1))
    association = (
        as_bool(drug["masld_preclinical_evidence"])
        | pd.to_numeric(drug["ot_masld_animal_model"], errors="coerce").fillna(0).gt(0)
        | pd.to_numeric(drug["ot_masld_literature"], errors="coerce").fillna(0).gt(0.2)
    )

    drug["clinical_evidence_status"] = np.select(
        [approved, discontinued, clinical, association],
        ["approved", "discontinued", "clinical", "association_only"],
        default="none",
    )
    drug["strict_label_source"] = np.select(
        [approved, discontinued, curated.eq("active"), ct_trial & ct_phase.ge(1), association],
        ["curated_approved", "curated_discontinued", "curated_active",
         "clinicaltrials_masld", "opentargets_association"],
        default="none",
    )
    drug["strict_positive"] = drug["clinical_evidence_status"].isin(["approved", "clinical"])

    label_cols = [
        "symbol", "in_atlas", "clinical_evidence_status", "strict_label_source",
        "strict_positive", "curated_status", "ct_has_masld_trial",
        "ct_masld_trial_max_phase", "pharos_tdl", "max_phase_any",
        "dgidb_n_drugs", "ot_masld_literature", "ot_masld_animal_model",
    ]
    drug[label_cols].sort_values("symbol").to_csv(
        OUT / "strict_drug_labels.tsv", sep="\t", index=False
    )

    score_cols = {
        "convergence_score": "convergence_score",
        "modality_count": "n_modalities_active",
        "genetics_pp4": "coloc_genetic_pp4",
        "bulk_log_bf": "log_BF_S1",
    }
    keep_conv = ["symbol", "excluded_from_ranking", *score_cols.values()]
    df = conv[keep_conv].merge(drug[label_cols], on="symbol", how="left")
    df = df.merge(
        release[["symbol", "bulk_AveExpr", "gene_biotype", "analysis_release_id"]],
        on="symbol",
        how="left",
    )
    df["strict_positive"] = df["strict_positive"].fillna(False).astype(bool)
    df["clinical_evidence_status"] = df["clinical_evidence_status"].fillna("none")
    df["protein_coding"] = df["gene_biotype"].eq("protein_coding")
    df["tdl_group"] = df["pharos_tdl"].fillna("none").replace({
        "Tclin": "clinical_chemistry", "Tchem": "clinical_chemistry",
        "Tbio": "biology_dark", "Tdark": "biology_dark",
    })
    df["expr_bin"] = pd.qcut(
        df["bulk_AveExpr"].rank(method="first"), 5,
        labels=False, duplicates="drop",
    )
    druggable = (
        df["pharos_tdl"].notna()
        | pd.to_numeric(df["max_phase_any"], errors="coerce").notna()
        | pd.to_numeric(df["dgidb_n_drugs"], errors="coerce").fillna(0).gt(0)
    )
    df["eligible_control"] = (~df["strict_positive"]) & druggable
    ranking = df[
        ~df["excluded_from_ranking"].fillna(False).astype(bool)
        & df["bulk_AveExpr"].notna()
    ].copy()
    matched = make_matches(ranking)
    matched.to_csv(OUT / "strict_drug_matched_universe.tsv", sep="\t", index=False)

    coverage = (
        df.groupby("clinical_evidence_status", dropna=False)
        .agg(
            n_all=("symbol", "size"),
            n_in_ranking=("excluded_from_ranking", lambda x: int((~x.fillna(False).astype(bool)).sum())),
        )
        .reset_index()
    )
    coverage["analysis_release_id"] = RELEASE_ID
    coverage.to_csv(OUT / "strict_drug_coverage.tsv", sep="\t", index=False)

    perf_rows = []
    comparison_rows = []
    y = matched["strict_positive"].to_numpy(bool)
    for label, col in score_cols.items():
        value = pd.to_numeric(matched[col], errors="coerce").to_numpy(float)
        perf_rows.append({
            "analysis_release_id": RELEASE_ID,
            "score": label,
            "n_positive": int(y.sum()),
            "n_control": int((~y).sum()),
            "auc": auc(value, y),
        })

    conv_score = pd.to_numeric(matched["convergence_score"], errors="coerce").to_numpy(float)
    for label, col in score_cols.items():
        if label == "convergence_score":
            continue
        baseline = pd.to_numeric(matched[col], errors="coerce").to_numpy(float)
        delta, lo, hi = stratified_boot_diff(conv_score, baseline, y)
        comparison_rows.append({
            "analysis_release_id": RELEASE_ID,
            "comparison": f"convergence_score_minus_{label}",
            "auc_difference": delta,
            "ci_low": lo,
            "ci_high": hi,
            "material_superiority": bool(delta >= 0.03 and lo > 0),
        })

    perf = pd.DataFrame(perf_rows)
    comparisons = pd.DataFrame(comparison_rows)
    perf.to_csv(OUT / "strict_drug_benchmark.tsv", sep="\t", index=False)
    comparisons.to_csv(OUT / "strict_drug_comparisons.tsv", sep="\t", index=False)

    n_approved_ranked = int((
        ranking["clinical_evidence_status"].eq("approved")
    ).sum())
    clinical_gate = bool(
        len(comparisons) > 0
        and comparisons["material_superiority"].all()
        and n_approved_ranked >= 2
    )
    verdict = {
        "analysis_release_id": RELEASE_ID,
        "seed": SEED,
        "n_boot": N_BOOT,
        "n_matched_positive": int(y.sum()),
        "n_matched_control": int((~y).sum()),
        "n_approved_in_ranking_universe": n_approved_ranked,
        "strict_preclinical_intervention_labels_available": False,
        "ordered_gradient_testable": False,
        "clinical_target_superiority_gate": clinical_gate,
        "headline_pass": False,
        "decision": (
            "FAIL: strict labels do not support an ordered drug-development gradient, "
            "and the atlas has insufficient approved-target coverage. Use clinical examples "
            "as mechanism-class boundaries; do not claim recovery of both approved targets, "
            "a full gradient, or prospective validation."
        ),
    }
    with open(OUT / "strict_drug_verdict.json", "w", encoding="utf-8") as handle:
        json.dump(verdict, handle, indent=2)

    gates_path = OUT / "acceptance_gates.tsv"
    if gates_path.exists():
        gates = pd.read_csv(gates_path, sep="\t")
        gates = gates[gates["gate"] != "drug_calibration_article"]
        gates = pd.concat([
            gates,
            pd.DataFrame([{
                "analysis_release_id": RELEASE_ID,
                "gate": "drug_calibration_article",
                "passed": False,
                "criterion": (
                    "ordered strict stages testable; >=2 approved targets in ranking; "
                    "convergence exceeds modality count and best channel by >=0.03 with CI>0"
                ),
                "detail": (
                    f"approved_in_ranking={n_approved_ranked}; matched_positive={int(y.sum())}; "
                    "strict_preclinical_intervention_labels=0"
                ),
            }]),
        ], ignore_index=True)
        gates.to_csv(gates_path, sep="\t", index=False)

    ledger_path = ROOT / "docs/manuscript/release/claim_ledger.tsv"
    if ledger_path.exists():
        ledger = pd.read_csv(ledger_path, sep="\t")
        hit = ledger["claim_id"].eq("C_DRUG_RETIRED")
        ledger.loc[hit, "status"] = "failed_acceptance_gate"
        ledger.loc[hit, "allowed_wording"] = (
            f"In the strict matched benchmark ({int(y.sum())} clinical/approved targets; "
            f"{int((~y).sum())} controls), convergence AUC was "
            f"{perf.loc[perf.score.eq('convergence_score'), 'auc'].iloc[0]:.3f}; "
            "the analysis did not support an ordered drug-development gradient."
        )
        ledger.loc[hit, "source_output"] = "strict_drug_benchmark.tsv; strict_drug_verdict.json"
        ledger.to_csv(ledger_path, sep="\t", index=False)

    print(json.dumps(verdict, indent=2))
    print(perf.to_string(index=False))
    print(comparisons.to_string(index=False))


if __name__ == "__main__":
    main()
