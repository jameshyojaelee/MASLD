#!/usr/bin/env python3
"""Paired, five-fold sensitivity for the fixed 20% acquisition budget."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main(source, out):
    table = pd.read_csv(source, sep="\t")
    primary = table.loc[np.isclose(table.budget_fraction, .2)]
    if set(primary.policy) != {"random", "histology", "rna_diversity", "global_local_disagreement"}:
        raise ValueError("Acquisition policy census differs")
    random = primary.loc[primary.policy.eq("random")].groupby("outer_fold").agg(
        capture=("captured_pool_error_fraction", "mean"),
        refit=("evaluation_MSE_reduction", "mean"))
    if len(random) != 5 or not np.all(primary.loc[primary.policy.eq("random")].groupby("outer_fold").size() == 3):
        raise ValueError("Expected five folds and three random orderings per fold")
    rows = []
    rng = np.random.default_rng(20260922)
    for policy in ("histology", "rna_diversity", "global_local_disagreement"):
        item = primary.loc[primary.policy.eq(policy)].set_index("outer_fold").sort_index()
        if len(item) != 5 or not item.index.equals(random.index):
            raise ValueError("Different evaluation folds")
        for metric, column in (("captured_error_fraction", "captured_pool_error_fraction"),
                               ("held_participant_MSE_reduction", "evaluation_MSE_reduction")):
            diff = item[column].to_numpy()-random["capture" if metric == "captured_error_fraction" else "refit"].to_numpy()
            draws = diff[rng.integers(0, 5, size=(10000, 5))].mean(1)
            null = np.array([np.mean(diff*np.array(signs)) for signs in itertools.product((-1, 1), repeat=5)])
            p = float(np.mean(np.abs(null) >= abs(diff.mean())-1e-15))
            rows.append({"policy": policy, "metric": metric, "n_disjoint_evaluation_folds": 5,
                         "difference_vs_random": float(diff.mean()),
                         "fold_bootstrap_low95": float(np.quantile(draws, .025)),
                         "fold_bootstrap_high95": float(np.quantile(draws, .975)),
                         "exact_fold_signflip_p_two_sided": p,
                         "fold_differences": json.dumps(diff.tolist())})
    output = pd.DataFrame(rows)
    for metric in output.metric.unique():
        index = output.index[output.metric.eq(metric)].to_numpy()
        order = index[np.argsort(output.loc[index, "exact_fold_signflip_p_two_sided"].to_numpy())]
        previous = 0.0
        for rank, row in enumerate(order):
            previous = max(previous, min(1.0, float(output.at[row, "exact_fold_signflip_p_two_sided"])*(len(order)-rank)))
            output.loc[row, "holm_p_within_metric"] = previous
    output.to_csv(out, sep="\t", index=False)
    print(output.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    main(a.source, a.out)
