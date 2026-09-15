#!/usr/bin/env python3
"""Print the C2/C3 numbers needed for RESULTS.md in one pass."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    out = arguments.output
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 50)

    summary = pd.read_csv(out / "tables/head_summary.tsv", sep="\t")
    print("=== head_summary (macro) ===")
    print(
        summary[
            [
                "stratum",
                "arm",
                "variants",
                "blocks",
                "macro_spearman",
                "macro_ci_low",
                "macro_ci_high",
                "gain_vs_allele_identity",
                "gain_ci_low",
                "gain_ci_high",
                "positive_gain_seeds",
                "comparable_seeds",
            ]
        ].round(5).to_string(index=False)
    )
    print()
    print("=== head_summary (pooled OOF) ===")
    print(
        summary[
            [
                "stratum",
                "arm",
                "pooled_spearman",
                "pooled_ci_low",
                "pooled_ci_high",
                "pooled_gain_vs_allele_identity",
                "pooled_gain_ci_low",
                "pooled_gain_ci_high",
            ]
        ].round(5).to_string(index=False)
    )

    per_seed = pd.read_csv(out / "tables/per_seed_metrics.tsv", sep="\t")
    print()
    print("=== per-seed spread (primary stratum) ===")
    spread = (
        per_seed[per_seed["stratum"] == "primary_all_variants"]
        .groupby("arm")["macro_spearman"]
        .agg(["min", "median", "max"])
        .round(5)
    )
    print(spread.to_string())

    per_fold = pd.read_csv(out / "tables/per_fold_metrics.tsv", sep="\t")
    print()
    print("=== per-fold Spearman, seed ensemble, primary stratum ===")
    print(
        per_fold[per_fold["stratum"] == "primary_all_variants"]
        .pivot(index="arm", columns="held_out_fold", values="spearman_seed_ensemble")
        .round(4)
        .to_string()
    )

    audit = pd.read_csv(out / "tables/selection_audit.tsv", sep="\t")
    print()
    print("=== selected alpha / feature count, primary campaign ===")
    print(
        audit.groupby("arm")
        .agg(
            fits=("selected_alpha", "size"),
            constant_mean_fits=("selected_alpha", lambda s: int((s.astype(str) == "constant_training_mean").sum())),
            median_k=("selected_feature_count", "median"),
            median_outer_train=("outer_training_variants", "median"),
        )
        .to_string()
    )

    curve_path = out / "tables/learning_curve.tsv"
    if curve_path.exists():
        curve = pd.read_csv(curve_path, sep="\t")
        print()
        print("=== C3 learning curve (macro Spearman, mean and sd over 5 seeds) ===")
        table = (
            curve.groupby(["model", "requested_training_size", "mean_training_variants_per_fold"])["macro_spearman"]
            .agg(["mean", "std", "min", "max"])
            .round(5)
            .reset_index()
            .sort_values(["model", "mean_training_variants_per_fold"])
        )
        print(table.to_string(index=False))
        print()
        print("=== C3 top-end increment (largest minus second largest size) ===")
        for model in sorted(curve["model"].unique()):
            sub = curve[curve["model"] == model]
            sizes = sorted(sub["mean_training_variants_per_fold"].unique())
            top, second = sizes[-1], sizes[-2]
            paired = []
            for seed in sorted(sub["seed"].unique()):
                a = sub[(sub["mean_training_variants_per_fold"] == top) & (sub["seed"] == seed)]["macro_spearman"].iloc[0]
                b = sub[(sub["mean_training_variants_per_fold"] == second) & (sub["seed"] == seed)]["macro_spearman"].iloc[0]
                paired.append(a - b)
            print(
                f"{model}: n_train {second:.0f} -> {top:.0f}, per-seed increments "
                f"{[round(v, 5) for v in paired]}, mean {np.mean(paired):.5f}"
            )

    precision = out / "tables/label_precision_posthoc.json"
    if precision.exists():
        print()
        print("=== post-hoc label precision ===")
        print(json.dumps(json.loads(precision.read_text()), indent=2, sort_keys=True))

    for name in ("hyenadna", "caduceus"):
        receipt = out / f"raw/{name}_allele_embeddings.npz.receipt.json"
        if receipt.exists():
            data = json.loads(receipt.read_text())
            print()
            print(
                f"=== {name} embedding receipt === variants={data['variants']} "
                f"reproducibility_max_abs_diff={data['reproducibility_max_abs_diff']:.3e} "
                f"fwd_vs_rc_max_abs_diff={data['forward_vs_reverse_complement_max_abs_diff']:.3e} "
                f"elapsed_s={data['elapsed_seconds']} checkpoint_sha256={data['checkpoint_sha256'][:16]}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
