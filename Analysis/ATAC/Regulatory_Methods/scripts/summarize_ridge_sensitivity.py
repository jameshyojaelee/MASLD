#!/usr/bin/env python3
"""Reduce a separate wider-ridge RNA sensitivity to paired donor comparisons.

The primary sensitivity family is three distal methods by two populations versus
the promoter-only model tuned on the same wider grid. Gene observations are first
averaged within donor. The original 48-comparison family is retained unchanged.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import network_compare as core


METHODS = ("distance", "cicero", "directnet_adapted")
ALL_METHODS = ("promoter_only", *METHODS)
LINEAGES = ("hepatocyte", "stellate")
GRIDS = {"original_replay": (0.1, 1., 10., 100.),
         "wide": (0.1, 1., 10., 100., 1000., 10000., 100000., 1000000., 10000000., 100000000.)}


def read(path):
    with Path(path).open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write(path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    with Path(path).open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def completed(path):
    root = Path(path)
    if any(root.glob("INVALIDATED*.json")):
        raise ValueError(f"Invalidated sensitivity component: {root}")
    receipt = json.loads((root / "completion.json").read_text())
    if receipt.get("status") != "completed":
        raise ValueError(f"Incomplete sensitivity component: {root}")
    if receipt.get("technical_check") is not False or receipt.get("max_genes") != 0:
        raise ValueError("Bounded technical checks cannot enter the full sensitivity")
    if receipt.get("endpoint") != "rna" or receipt.get("original_grid") != list(GRIDS["original_replay"]) or receipt.get("wide_grid") != list(GRIDS["wide"]):
        raise ValueError("Completed fold endpoint or ridge grid differs")
    return receipt


def roster(path, expected_donors):
    result = {}
    counts = Counter()
    for row in read(path):
        donor, fold = row["donor_id"], int(row["outer_fold"])
        if donor in result and result[donor][0] != fold:
            raise ValueError("Donor occurs in more than one outer fold")
        counts[donor] += 1
        result[donor] = (fold, 0)
    if len(result) != expected_donors or {x[0] for x in result.values()} != set(range(5)):
        raise ValueError("Original five-fold donor roster is incomplete")
    return {donor: (fold, counts[donor]) for donor, (fold, _) in result.items()}


def number(row, field, nonnegative=False):
    value = float(row[field])
    if not math.isfinite(value) or (nonnegative and value < 0):
        raise ValueError(f"Invalid {field} in prediction")
    return value


def add(summary, contrast, donor, gene, candidate, baseline):
    item = summary[contrast][donor]
    item[0] += baseline - candidate
    item[1] += 1
    item[2] += candidate
    item[3] += baseline


def reduce_fold(path, lineage, fold, donors, totals, genes, alpha, coverage, model_genes):
    """Stream one fold; retain at most one gene's donor/method/grid rows."""
    held = {donor for donor, (assigned, _) in donors.items() if assigned == fold}
    previous = None
    nrows = 0
    for gene, records in itertools.groupby(read(path), key=lambda row: row["gene_id"]):
        if previous is not None and gene <= previous:
            raise ValueError("Gene groups must be strictly ordered within each fold")
        previous = gene
        index = {}
        for row in records:
            nrows += 1
            grid, method, donor = row["grid"], row["method"], row["donor_id"]
            if grid not in GRIDS or method not in ALL_METHODS:
                raise ValueError("Unexpected grid or method")
            if row["endpoint"] != "rna" or row.get("cohort", "GSE296875") != "GSE296875":
                raise ValueError("Sensitivity must contain genuine GSE296875 RNA predictions")
            if row.get("lineage", lineage) != lineage or int(row["outer_fold"]) != fold:
                raise ValueError("Prediction context or fold mismatch")
            if donor not in held or int(row["n_metacells"]) != donors[donor][1]:
                raise ValueError("Prediction donor fold or aggregate coverage differs from original roster")
            key = grid, method, donor
            if key in index:
                raise ValueError("Duplicated grid/gene/method/donor prediction")
            number(row, "predicted")
            number(row, "observed")
            number(row, "squared_error", True)
            number(row, "standardized_squared_error", True)
            if number(row, "ridge_alpha", True) not in GRIDS[grid]:
                raise ValueError("Selected alpha is outside the recorded grid")
            index[key] = row
            coverage[grid, method] += 1
        old = {(method, donor) for grid, method, donor in index if grid == "original_replay"}
        new = {(method, donor) for grid, method, donor in index if grid == "wide"}
        if not old or old != new:
            raise ValueError("Original/wide prediction coverage differs")
        available = {method for method, donor in old}
        if "promoter_only" not in available:
            raise ValueError("Promoter-only baseline is missing")
        for method in available:
            if {donor for candidate, donor in old if candidate == method} != held:
                raise ValueError("An evaluable model omits held donors")
            for grid in GRIDS:
                models = [index[grid, method, donor] for donor in sorted(held)]
                alphas = {float(row["ridge_alpha"]) for row in models}
                if len(alphas) != 1:
                    raise ValueError("Ridge tuning varies across held donors")
                alpha[grid, method, alphas.pop()] += 1
                model_genes[grid, method].add(gene)
            for donor in sorted(held):
                a, b = index["wide", method, donor], index["original_replay", method, donor]
                if (a["observed"], a.get("target_unit"), a["n_metacells"]) != (b["observed"], b.get("target_unit"), b["n_metacells"]):
                    raise ValueError("Grid change altered held target or aggregate identity")
                contrast = "wide_vs_original", method
                add(totals, contrast, donor, gene, float(a["standardized_squared_error"]), float(b["standardized_squared_error"]))
                genes[contrast].add(gene)
        for grid in GRIDS:
            for method in METHODS:
                if method not in available:
                    continue
                contrast = grid, method
                for donor in sorted(held):
                    candidate, base = index[grid, method, donor], index[grid, "promoter_only", donor]
                    if candidate["observed"] != base["observed"] or candidate.get("target_unit") != base.get("target_unit"):
                        raise ValueError("Candidate and promoter-only target differ")
                    add(totals, contrast, donor, gene, float(candidate["standardized_squared_error"]), float(base["standardized_squared_error"]))
                genes[contrast].add(gene)
    if not nrows:
        raise ValueError("Empty completed fold")
    return nrows


def aggregate_row(lineage, contrast, values, genes, seed, reps):
    grid, method = contrast
    ordered = sorted(values.items())
    differences = [item[0] / item[1] for _, item in ordered]
    if grid == "wide_vs_original":
        rng_seed = core.stable_seed(seed, ("rna", method, "original_grid"))
    else:
        # Preserve the original inference seed for direct reference reconstruction.
        rng_seed = core.stable_seed(seed, ("rna", method, "promoter_only"))
    statistics = core.paired_inference(differences, rng_seed, reps)
    donor_count = len(ordered)
    candidate = sum(item[2] / item[1] for _, item in ordered) / donor_count if donor_count else float("nan")
    baseline = sum(item[3] / item[1] for _, item in ordered) / donor_count if donor_count else float("nan")
    result = dict(cohort="GSE296875", lineage=lineage,
                  source_population="Hepatocyte" if lineage == "hepatocyte" else "Mesenchymal",
                  grid=grid, endpoint="rna", method=method,
                  baseline="same_method_original_grid" if grid == "wide_vs_original" else "promoter_only",
                  n_common_predictions=sum(item[1] for _, item in ordered), n_common_genes=len(genes),
                  **statistics, candidate_standardized_mse=candidate, baseline_standardized_mse=baseline,
                  relative_mse_reduction=statistics["improvement"] / baseline if baseline > 0 else float("nan"),
                  fraction_positive_donor_gains=sum(value > 0 for value in differences) / donor_count if donor_count else float("nan"))
    donor_rows = [dict(lineage=lineage, grid=grid, method=method,
                       baseline=result["baseline"], donor_id=donor,
                       n_common_genes=item[1], mean_loss_difference=item[0] / item[1],
                       candidate_standardized_mse=item[2] / item[1], baseline_standardized_mse=item[3] / item[1])
                  for donor, item in ordered]
    return result, donor_rows


def check_reference(rows, reference):
    original = {}
    for row in read(reference):
        if row["cohort"] == "GSE296875" and row["endpoint"] == "rna" and row["baseline"] == "promoter_only":
            key = row["lineage"], row["method"]
            if key in original:
                raise ValueError("Duplicate original reference comparison")
            original[key] = row
    expected = {(lineage, method) for lineage in LINEAGES for method in METHODS}
    if set(original) != expected:
        raise ValueError("Original reference does not contain all six RNA comparisons")
    for row in rows:
        old = original[row["lineage"], row["method"]]
        for field in ("n_donors", "n_common_predictions", "n_common_genes", "improvement", "ci_low", "ci_high", "pvalue"):
            if not math.isclose(float(row[field]), float(old[field]), rel_tol=1e-9, abs_tol=1e-11):
                raise ValueError(f"Original-grid replay fails reference comparison: {row['lineage']}/{row['method']}/{field}")
        row.update(holm_pvalue=float(old["holm_pvalue"]), family_size=48,
                   family=old["family"], predictive_improvement=old["predictive_improvement"] in {"True", "true"},
                   interpretation="original_grid_reference_with_original_48_comparison_multiplicity")


def adjusted(rows, family):
    values = core.holm([row["pvalue"] for row in rows])
    for row, value in zip(rows, values):
        row.update(holm_pvalue=float(value), family_size=len(rows), family=family,
                   predictive_improvement=bool(row["ci_low"] > 0 and value < .05),
                   interpretation="secondary_sensitivity_conditional_on_fixed_donor_cv_predictions")


def summarize(root, reference, output, seed=42, reps=2000):
    if output.exists():
        raise FileExistsError(f"Refusing overwrite: {output}")
    summaries, old_rows, grid_rows, donor_rows, alpha_rows, coverage_rows, receipts = [], [], [], [], [], [], []
    for lineage in LINEAGES:
        prepared = root / lineage / "prepared"
        donors = roster(prepared / "fold_assignments.tsv", 39 if lineage == "hepatocyte" else 18)
        totals = defaultdict(lambda: defaultdict(lambda: [0., 0, 0., 0.]))
        genes, model_genes = defaultdict(set), defaultdict(set)
        alpha, coverage = Counter(), Counter()
        for fold in range(5):
            directory = root / lineage / f"fold{fold}"
            receipt = completed(directory)
            if int(receipt.get("outer_fold", fold)) != fold:
                raise ValueError("Fold completion identity differs")
            if receipt.get("seed") != seed or receipt.get("lineage") != lineage or receipt.get("group") != f"GSE296875_{lineage}":
                raise ValueError("Fold completion seed or population differs")
            if "prepared_root" in receipt and Path(receipt["prepared_root"]).resolve() != prepared.resolve():
                raise ValueError("Fold uses a different prepared input")
            nrows = reduce_fold(directory / "predictions.tsv", lineage, fold, donors, totals, genes, alpha, coverage, model_genes)
            if receipt.get("n_prediction_rows") != nrows or receipt.get("n_original_predictions_verified", -1) * 2 != nrows:
                raise ValueError("Completed prediction or original-replay count differs from observed rows")
            print(f"Reduced {lineage} fold={fold}: {nrows:,} paired-grid prediction rows", flush=True)
            receipts.append(dict(lineage=lineage, outer_fold=fold, rows=nrows,
                                 source=str(directory.resolve()), completion_sha256=digest(directory / "completion.json"),
                                 prediction_bytes=(directory / "predictions.tsv").stat().st_size))
        for contrast in [(grid, method) for grid in GRIDS for method in METHODS] + [("wide_vs_original", method) for method in ALL_METHODS]:
            if set(totals[contrast]) != set(donors):
                raise ValueError(f"Comparison lacks original donors: {lineage}/{contrast}")
            row, donor = aggregate_row(lineage, contrast, totals[contrast], genes[contrast], seed, reps)
            (summaries if contrast[0] == "wide" else old_rows if contrast[0] == "original_replay" else grid_rows).append(row)
            donor_rows.extend(donor)
        for grid in GRIDS:
            for method in ALL_METHODS:
                nmodels = sum(count for (g, m, a), count in alpha.items() if (g, m) == (grid, method))
                if not nmodels:
                    raise ValueError("A planned method/grid has no selected ridge models")
                for value in GRIDS[grid]:
                    count = alpha[grid, method, value]
                    alpha_rows.append(dict(lineage=lineage, grid=grid, method=method, ridge_alpha=value,
                                           n_models=count, total_models=nmodels, fraction_models=count / nmodels,
                                           lower_boundary=value == min(GRIDS[grid]), upper_boundary=value == max(GRIDS[grid])))
                coverage_rows.append(dict(lineage=lineage, grid=grid, method=method,
                                          n_predictions=coverage[grid, method], n_models=nmodels,
                                          n_genes=len(model_genes[grid, method]), n_donors=len(donors), n_outer_folds=5))
    if len(summaries) != 6 or len(grid_rows) != 8 or len(receipts) != 10:
        raise ValueError("The fixed sensitivity family is incomplete")
    check_reference(old_rows, reference)
    adjusted(summaries, "secondary_wide_grid_6_distal_method_by_population_RNA_vs_promoter_only")
    adjusted(grid_rows, "descriptive_8_method_by_population_wide_vs_original_grid")
    output.mkdir(parents=True, exist_ok=False)
    for filename, rows in [("comparisons.tsv", summaries), ("original_grid_reference.tsv", old_rows),
                           ("grid_effects.tsv", grid_rows), ("donor_scores.tsv", donor_rows),
                           ("ridge_alpha_distribution.tsv", alpha_rows), ("coverage.tsv", coverage_rows)]:
        write(output / filename, rows)
    success = sum(row["predictive_improvement"] for row in summaries)
    lines = ["Separate wider-ridge RNA-prediction sensitivity", "",
             f"Completed original donor folds: {len(receipts)}/10; paired populations: 39 Hepatocyte donors and 18 source-Mesenchymal donors.",
             "Original grid: " + ", ".join(map(str, GRIDS["original_replay"])),
             "Wider grid: " + ", ".join(map(str, GRIDS["wide"])),
             f"Distal-vs-promoter-only support: {success}/6 comparisons (positive 95% donor-bootstrap interval and Holm p<0.05 across six).", ""]
    for row in summaries:
        lines.append(f"{row['source_population']} {row['method']}: gain {row['improvement']:.10g} [{row['ci_low']:.10g}, {row['ci_high']:.10g}], Holm p={row['holm_pvalue']:.10g}; n={row['n_donors']} donors, {row['n_common_genes']} matched genes.")
    lines.extend(["", "Upper-bound alpha frequencies (unique fold/gene/model fits):"])
    for lineage in LINEAGES:
        for method in ALL_METHODS:
            boundary = {row["grid"]: row for row in alpha_rows
                        if row["lineage"] == lineage and row["method"] == method and row["upper_boundary"]}
            lines.append(f"{lineage} {method}: original alpha=100 {boundary['original_replay']['fraction_models']:.4%}; wider alpha=1e8 {boundary['wide']['fraction_models']:.4%}.")
    lines.extend(["", "Positive gain means reduced held-donor RNA MSE beyond promoter accessibility alone, standardized by outer-training target variance.",
                  "Every comparison uses contrast-matched genes, then averages squared prediction errors within donor before averaging donors.",
                  "Selected-alpha frequencies count unique outer-fold/gene/method fits; donor prediction rows are not independent ridge fits.",
                  "Original-grid replay reproduces the six original RNA comparisons; their original 48-contrast Holm probabilities are retained unchanged.",
                  "Grid-effect comparisons describe how retuning changes each method, including promoter-only, and have a separate eight-contrast Holm family.",
                  "This is a secondary sensitivity prompted by original upper-bound tuning. It does not replace the original comparison or provide a newly untouched validation cohort.",
                  "Intervals use 2,000 donor bootstrap resamples conditional on fitted cross-validation predictions; sign-flip p-values assume independent symmetric donor differences.",
                  "Shared training folds mean uncertainty from refitting models is not included. Mesenchymal cells are source-defined, not purified stellate cells.",
                  "No causal regulatory claim or MASH replication follows from these predictions."])
    (output / "REPORT.txt").write_text("\n".join(lines) + "\n")
    (output / "completion.json").write_text(json.dumps(dict(
        status="completed", n_folds=len(receipts), n_primary_comparisons=len(summaries),
        n_supported_primary_comparisons=success, primary_holm_family_size=6,
        n_original_reference_comparisons=len(old_rows), n_descriptive_grid_comparisons=len(grid_rows),
        seed=seed, bootstrap_reps=reps, original_reference=str(reference.resolve()),
        original_reference_sha256=digest(reference), source_root=str(root.resolve()),
        comparison_unit="donor", conditional_on_fitted_cv_predictions=True,
        wider_grid_secondary_sensitivity=True, scientific_adoption=False, folds=receipts), indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--original-comparisons", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--bootstrap-reps", type=int, default=2000)
    args = parser.parse_args()
    core.require_compute()
    if args.bootstrap_reps != 2000:
        raise ValueError("The specified sensitivity uses 2,000 donor bootstrap resamples")
    summarize(args.root, args.original_comparisons, args.output, args.seed, args.bootstrap_reps)


if __name__ == "__main__":
    main()
