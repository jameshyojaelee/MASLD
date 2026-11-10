#!/usr/bin/env python3
"""Independently reconstruct the six wider-ridge donor comparisons."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import socket

import numpy as np
import pandas as pd


METHODS = ("distance", "cicero", "directnet_adapted")
ALL_METHODS = ("promoter_only", *METHODS)
GRIDS = {"original_replay": [.1, 1., 10., 100.],
         "wide": [.1, 1., 10., 100., 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]}


def seed_for(method):
    value = str(("rna", method, "promoter_only"))
    return (42 + int(hashlib.sha256(value.encode()).hexdigest()[:8], 16)) % 2**32


def inference(differences, method):
    values = np.asarray(differences, float)
    rng = np.random.default_rng(seed_for(method))
    resampled = np.array([rng.choice(values, size=len(values), replace=True).mean() for _ in range(2000)])
    estimate = values.mean()
    if len(values) <= 15:
        statistics = []
        for pattern in range(2 ** len(values)):
            signs = 2 * ((pattern >> np.arange(len(values))) & 1) - 1
            statistics.append(np.mean(values * signs))
        probability = np.mean(np.asarray(statistics) >= estimate - 1e-15)
    else:
        exceeds = 0
        for _ in range(10000):
            exceeds += np.mean(values * rng.choice([-1, 1], size=len(values))) >= estimate - 1e-15
        probability = (exceeds + 1) / 10001
    lower, upper = np.quantile(resampled, [.025, .975])
    return dict(improvement=float(estimate), ci_low=float(lower), ci_high=float(upper), pvalue=float(probability))


def holm(probabilities):
    ordered = sorted(enumerate(probabilities), key=lambda value: (value[1], value[0]))
    result, previous = [None] * len(ordered), 0.
    for rank, (index, value) in enumerate(ordered):
        previous = max(previous, min(1., value * (len(ordered) - rank)))
        result[index] = previous
    return result


def blocks(path):
    """Pandas chunks preserve complete gene groups using one trailing carry."""
    carry = None
    for chunk in pd.read_csv(path, sep="\t", chunksize=100000, dtype={"donor_id": str}, float_precision="round_trip"):
        if carry is not None:
            chunk = pd.concat([carry, chunk], ignore_index=True)
        last = chunk.gene_id.iloc[-1]
        carry = chunk[chunk.gene_id.eq(last)].copy()
        complete = chunk[~chunk.gene_id.eq(last)]
        if len(complete):
            yield complete
    if carry is not None and len(carry):
        yield carry


def check(root, summary, output):
    if output.exists():
        raise FileExistsError(f"Refusing overwrite: {output}")
    if not os.environ.get("SLURM_JOB_ID") or "login" in socket.gethostname():
        raise RuntimeError("Full prediction rederivation requires a compute allocation")
    summary_state = json.loads((summary / "completion.json").read_text())
    if summary_state.get("status") != "completed" or summary_state.get("n_folds") != 10:
        raise ValueError("Completed ten-fold summary required")
    totals = defaultdict(lambda: np.zeros(4, dtype=float))
    matched_genes = defaultdict(set)
    alpha_counts, model_counts, prediction_counts = Counter(), Counter(), Counter()
    nrecords, fold_receipts = 0, []
    for lineage, expected_n in [("hepatocyte", 39), ("stellate", 18)]:
        roster = pd.read_csv(root / lineage / "prepared/fold_assignments.tsv", sep="\t", dtype={"donor_id": str})
        if roster.donor_id.nunique() != expected_n or roster.groupby("donor_id").outer_fold.nunique().max() != 1:
            raise ValueError("Invalid donor roster")
        for fold in range(5):
            directory = root / lineage / f"fold{fold}"
            state = json.loads((directory / "completion.json").read_text())
            if state.get("status") != "completed" or state.get("technical_check") is not False:
                raise ValueError("Incomplete or bounded fold")
            held = roster[roster.outer_fold.eq(fold)].groupby("donor_id").size()
            seen_genes, fold_rows = set(), 0
            for frame in blocks(directory / "predictions.tsv"):
                fold_rows += len(frame)
                nrecords += len(frame)
                if not frame.outer_fold.eq(fold).all() or not frame.endpoint.eq("rna").all() or not frame.lineage.eq(lineage).all():
                    raise ValueError("Wrong endpoint, lineage or outer fold")
                if not frame.target_unit.eq("log1p_RNA_counts_per_10000").all() or not frame.input_representation.eq("binary").all():
                    raise ValueError("Wrong RNA effect units or ATAC representation")
                if set(frame.method) - set(ALL_METHODS) or set(frame.grid) != set(GRIDS):
                    raise ValueError("Unexpected method or grid")
                current_genes = set(frame.gene_id)
                if current_genes & seen_genes:
                    raise ValueError("A gene appears in disjoint prediction blocks")
                seen_genes.update(current_genes)
                keys = ["gene_id", "method", "donor_id"]
                if frame.duplicated(keys + ["grid"]).any():
                    raise ValueError("Duplicate prediction")
                if not frame.donor_id.isin(held.index).all() or not np.array_equal(frame.n_metacells, frame.donor_id.map(held)):
                    raise ValueError("Wrong held donor or aggregate count")
                numeric = frame[["ridge_alpha", "observed", "predicted", "squared_error", "standardized_squared_error"]].to_numpy(float)
                if not np.isfinite(numeric).all() or (frame[["squared_error", "standardized_squared_error"]] < 0).any().any():
                    raise ValueError("Invalid prediction values")
                old = frame[frame.grid.eq("original_replay")].set_index(keys).sort_index()
                wide = frame[frame.grid.eq("wide")].set_index(keys).sort_index()
                if not old.index.equals(wide.index):
                    raise ValueError("Original/wide coverage differs")
                for field in ["n_metacells", "observed", "target_unit", "input_representation"]:
                    if not old[field].equals(wide[field]):
                        raise ValueError(f"Original/wide {field} differs")
                for grid, selected in frame.groupby("grid", sort=False):
                    grouped = selected.groupby(["gene_id", "method"])
                    if not grouped.ridge_alpha.nunique().eq(1).all() or not grouped.donor_id.nunique().eq(len(held)).all():
                        raise ValueError("Alpha or donor coverage varies within a held model")
                    unique_models = grouped.ridge_alpha.first().reset_index()
                    if not unique_models.ridge_alpha.isin(GRIDS[grid]).all():
                        raise ValueError("Selected alpha outside grid")
                    for (method, alpha), count in unique_models.groupby(["method", "ridge_alpha"]).size().items():
                        alpha_counts[lineage, grid, method, float(alpha)] += int(count)
                        model_counts[lineage, grid, method] += int(count)
                    for method, count in selected.method.value_counts().items():
                        prediction_counts[lineage, grid, method] += int(count)
                scored = frame[frame.grid.eq("wide")]
                promoter = scored[scored.method.eq("promoter_only")]
                for method in METHODS:
                    candidate = scored[scored.method.eq(method)]
                    joined = candidate.merge(promoter, on=["gene_id", "donor_id"], suffixes=("_candidate", "_baseline"), validate="one_to_one")
                    if len(joined) != len(candidate):
                        raise ValueError("Missing promoter baseline on candidate genes")
                    if not joined.observed_candidate.equals(joined.observed_baseline):
                        raise ValueError("Candidate and baseline targets differ")
                    joined["difference"] = joined.standardized_squared_error_baseline - joined.standardized_squared_error_candidate
                    matched_genes[lineage, method].update(joined.gene_id)
                    sums = joined.groupby("donor_id").agg(
                        difference=("difference", "sum"), count=("difference", "size"),
                        candidate=("standardized_squared_error_candidate", "sum"),
                        baseline=("standardized_squared_error_baseline", "sum"))
                    for donor, values in sums.iterrows():
                        totals[lineage, method, donor] += values.to_numpy(float)
            if fold_rows != state["n_prediction_rows"] or fold_rows != 2 * state["n_original_predictions_verified"]:
                raise ValueError("Fold prediction count differs from completion")
            fold_receipts.append(dict(lineage=lineage, fold=fold, prediction_records=fold_rows))
    comparisons, donor_rows = [], []
    for lineage, expected_n in [("hepatocyte", 39), ("stellate", 18)]:
        for method in METHODS:
            donors = sorted(donor for (population, candidate, donor) in totals if (population, candidate) == (lineage, method))
            if len(donors) != expected_n:
                raise ValueError("A comparison lacks original donors")
            values = np.array([totals[lineage, method, donor] for donor in donors])
            donor_difference = values[:, 0] / values[:, 1]
            candidate_mse = np.mean(values[:, 2] / values[:, 1])
            baseline_mse = np.mean(values[:, 3] / values[:, 1])
            result = dict(lineage=lineage, method=method, **inference(donor_difference, method),
                          n_donors=len(donors), n_common_predictions=int(values[:, 1].sum()),
                          n_common_genes=len(matched_genes[lineage, method]),
                          candidate_standardized_mse=float(candidate_mse), baseline_standardized_mse=float(baseline_mse))
            comparisons.append(result)
            for donor, value in zip(donors, values):
                donor_rows.append(dict(lineage=lineage, method=method, donor_id=donor,
                    n_common_genes=int(value[1]), mean_loss_difference=float(value[0]/value[1]),
                    candidate_standardized_mse=float(value[2]/value[1]), baseline_standardized_mse=float(value[3]/value[1])))
    adjusted = holm([row["pvalue"] for row in comparisons])
    for row, value in zip(comparisons, adjusted):
        row["holm_pvalue"] = value
    expected = pd.read_csv(summary / "comparisons.tsv", sep="\t", float_precision="round_trip")
    if len(expected) != 6:
        raise ValueError("Summary comparison family must contain exactly six rows")
    actual = pd.DataFrame(comparisons)
    maximum = {}
    for field in ["improvement", "ci_low", "ci_high", "pvalue", "holm_pvalue", "n_donors", "n_common_predictions", "n_common_genes", "candidate_standardized_mse", "baseline_standardized_mse"]:
        merged = actual.merge(expected, on=["lineage", "method"], suffixes=("_actual", "_expected"), validate="one_to_one")
        if len(merged) != 6:
            raise ValueError("Summary does not contain exactly six matching comparisons")
        delta = float(np.max(np.abs(merged[field+"_actual"] - merged[field+"_expected"])))
        maximum[field] = delta
        if delta > 1e-10:
            raise ValueError(f"Independent comparison mismatch: {field}: {delta}")
    saved_donors = pd.read_csv(summary / "donor_scores.tsv", sep="\t", dtype={"donor_id": str}, float_precision="round_trip")
    saved_donors = saved_donors[saved_donors.grid.eq("wide")]
    donor_comparison = pd.DataFrame(donor_rows).merge(saved_donors, on=["lineage", "method", "donor_id"], suffixes=("_actual", "_expected"), validate="one_to_one")
    if len(donor_comparison) != len(donor_rows) or len(saved_donors) != len(donor_rows):
        raise ValueError("Donor comparison coverage differs")
    donor_max = {}
    for field in ["n_common_genes", "mean_loss_difference", "candidate_standardized_mse", "baseline_standardized_mse"]:
        delta = float(np.max(np.abs(donor_comparison[field+"_actual"] - donor_comparison[field+"_expected"])))
        donor_max[field] = delta
        if delta > 1e-10:
            raise ValueError(f"Independent donor mismatch: {field}: {delta}")
    alpha_table = pd.read_csv(summary / "ridge_alpha_distribution.tsv", sep="\t", float_precision="round_trip")
    expected_keys = {(lineage, grid, method, alpha) for lineage in ["hepatocyte", "stellate"] for grid in GRIDS for method in ALL_METHODS for alpha in GRIDS[grid]}
    actual_keys = set()
    alpha_fraction_max = 0.
    for row in alpha_table.itertuples():
        key = row.lineage, row.grid, row.method, float(row.ridge_alpha)
        if key in actual_keys:
            raise ValueError("Duplicate alpha occupancy")
        actual_keys.add(key)
        total = model_counts[key[:3]]
        if int(row.n_models) != alpha_counts[key] or int(row.total_models) != total:
            raise ValueError("Independent alpha occupancy count mismatch")
        delta = abs(float(row.fraction_models) - alpha_counts[key] / total)
        alpha_fraction_max = max(alpha_fraction_max, delta)
        if delta > 1e-14:
            raise ValueError("Independent alpha occupancy fraction mismatch")
    if actual_keys != expected_keys:
        raise ValueError("Alpha occupancy family incomplete")
    coverage = pd.read_csv(summary / "coverage.tsv", sep="\t")
    if len(coverage) != 16 or coverage.duplicated(["lineage", "grid", "method"]).any():
        raise ValueError("Coverage family incomplete or duplicated")
    for row in coverage.itertuples():
        key = row.lineage, row.grid, row.method
        if int(row.n_predictions) != prediction_counts[key] or int(row.n_models) != model_counts[key]:
            raise ValueError("Independent prediction/model coverage mismatch")
    receipt = dict(status="completed", source_root=str(root.resolve()), summary=str(summary.resolve()),
                   independent_implementation=True, n_folds=len(fold_receipts), n_prediction_records=nrecords,
                   n_comparisons=6, n_donor_contrast_scores=len(donor_rows),
                   n_alpha_occupancy_rows=len(actual_keys), alpha_fraction_max_difference=alpha_fraction_max,
                   comparison_max_absolute_differences=maximum, donor_max_absolute_differences=donor_max,
                   bootstrap_reps=2000, signflip_reps=10000, seed=42, family_size=6,
                   original_wide_coverage_equal=True, selected_alpha_constant_within_held_model=True,
                   folds=fold_receipts, comparisons=comparisons)
    with output.open("x") as handle:
        json.dump(receipt, handle, indent=2)
        handle.write("\n")
    print(json.dumps({key: receipt[key] for key in ["status", "n_prediction_records", "n_comparisons", "n_donor_contrast_scores", "comparison_max_absolute_differences"]}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ["root", "summary", "output"]:
        parser.add_argument(f"--{option}", required=True, type=Path)
    args = parser.parse_args()
    check(args.root, args.summary, args.output)
