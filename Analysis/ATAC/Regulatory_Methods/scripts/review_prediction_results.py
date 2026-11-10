#!/usr/bin/env python3
"""Independently reduce completed ordered prediction tables to donor comparisons."""
import argparse
import csv
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


def read(path):
    with Path(path).open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def unique_rows(path, fields):
    result = {}
    for row in read(path):
        key = tuple(row[field] for field in fields)
        if key in result:
            raise ValueError(f"Duplicate recorded key in {path}: {key}")
        result[key] = row
    return result


def contrasts(endpoint):
    base = "promoter_only" if endpoint == "rna" else "intercept"
    return [("distance", base), ("cicero", base), ("cicero", "distance"),
            ("directnet_adapted", base), ("directnet_adapted", "distance"),
            ("directnet_adapted", "cicero")]


def inference(values, seed, reps=2000):
    """Donor bootstrap and one-sided symmetric sign null, fixed predictions."""
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite donor differences")
    if len(values) < 5:
        return dict(n_donors=len(values), improvement=float("nan"), ci_low=float("nan"),
                    ci_high=float("nan"), pvalue=1., status="untestable")
    generator = np.random.default_rng(seed)
    means = [np.mean(generator.choice(values, size=len(values), replace=True)) for _ in range(reps)]
    observed = float(np.mean(values))
    if len(values) <= 15:
        exceed = 0
        for bits in range(1 << len(values)):
            signed = [value * (1 if bits & (1 << i) else -1) for i, value in enumerate(values)]
            exceed += np.mean(signed) >= observed - 1e-15
        probability = exceed / (1 << len(values))
    else:
        exceed = sum(np.mean(values * generator.choice([-1, 1], size=len(values))) >= observed - 1e-15
                     for _ in range(10000))
        probability = (exceed + 1) / 10001
    low, high = np.quantile(means, [.025, .975])
    return dict(n_donors=len(values), improvement=observed, ci_low=float(low),
                ci_high=float(high), pvalue=float(probability), status="evaluated")


def score_influence(values):
    """Summarize donor score influence without refitting any model."""
    values = list(values)
    omitted = [math.fsum(values[:i] + values[i+1:]) / (len(values)-1)
               for i in range(len(values))] if len(values) > 1 else []
    return dict(fraction_positive_donor_gains=sum(v > 0 for v in values)/len(values) if values else float("nan"),
                min_donor_gain=min(values, default=float("nan")),
                max_donor_gain=max(values, default=float("nan")),
                min_leave_one_donor_out_mean_gain=min(omitted, default=float("nan")),
                max_leave_one_donor_out_mean_gain=max(omitted, default=float("nan")))


def planned():
    for cohort in ("GSE244832", "GSE281367", "GSE296875"):
        for lineage in ("hepatocyte", "stellate"):
            for endpoint in (["promoter_accessibility", "rna"] if cohort == "GSE296875" else ["promoter_accessibility"]):
                for method, baseline in contrasts(endpoint):
                    yield cohort, lineage, endpoint, method, baseline


def holm(probabilities):
    adjusted = [1.] * len(probabilities)
    previous = 0.
    for rank, index in enumerate(sorted(range(len(probabilities)), key=probabilities.__getitem__)):
        previous = max(previous, min(1., probabilities[index] * (len(probabilities) - rank)))
        adjusted[index] = previous
    return adjusted


def check_number(actual, expected, label):
    a, b = float(actual), float(expected)
    if math.isnan(a) and math.isnan(b):
        return
    if not math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-13):
        raise ValueError(f"Mismatch {label}: independently {a}, recorded {b}")


def reduce_predictions(path, roster, cohort, lineage, endpoints):
    totals = defaultdict(lambda: defaultdict(lambda: [0., 0, 0., 0.]))
    genes = defaultdict(set)
    coverage = Counter()
    previous = None
    rows_seen = 0
    for group_key, group in itertools.groupby(read(path), key=lambda row: (int(row["outer_fold"]), row["gene_id"])):
        if previous is not None and group_key <= previous:
            raise ValueError("Prediction fold/gene groups are not strictly ordered")
        previous = group_key
        indexed = {}
        for row in group:
            rows_seen += 1
            donor = row["donor_id"]
            if (row["cohort"], row["lineage"]) != (cohort, lineage):
                raise ValueError("Prediction context mismatch")
            if donor not in roster or roster[donor] != (group_key[0], int(row["n_metacells"])):
                raise ValueError("Held donor fold or aggregate count mismatch")
            key = row["endpoint"], row["method"], donor
            if key in indexed:
                raise ValueError("Duplicate prediction")
            if key[0] not in endpoints or key[1] not in {"distance", "cicero", "directnet_adapted", "intercept", "promoter_only"}:
                raise ValueError("Unexpected endpoint or method")
            loss = float(row["standardized_squared_error"])
            if not math.isfinite(loss) or loss < 0:
                raise ValueError("Invalid standardized loss")
            indexed[key] = loss
            coverage[key[:2]] += 1
        for endpoint in endpoints:
            for method, baseline in contrasts(endpoint):
                contrast = endpoint, method, baseline
                for donor in sorted(roster):
                    a, b = (endpoint, method, donor), (endpoint, baseline, donor)
                    if a not in indexed or b not in indexed:
                        continue
                    item = totals[contrast][donor]
                    item[0] += indexed[b] - indexed[a]
                    item[1] += 1
                    item[2] += indexed[a]
                    item[3] += indexed[b]
                    genes[contrast].add(group_key[1])
    return totals, genes, coverage, rows_seen


def review(root):
    complete = json.loads((root / "completion.json").read_text())
    if complete.get("status") != "completed" or any(root.glob("INVALIDATED*.json")):
        raise ValueError("Incomplete or invalidated analysis")
    settings = json.loads((root / "run_settings.json").read_text())
    if settings.get("pilot", True) or settings.get("source_pilot", False) or settings.get("analysis_scope") != "genome_wide":
        raise ValueError("Not a genome-wide completed comparison")
    folders = list(root.glob("*/predictions.tsv"))
    if len(folders) != 1:
        raise ValueError("Expected one prediction context")
    folder = folders[0].parent
    cohort, lineage = folder.name.split("_", 1)
    roster, metacells = {}, set()
    counts = Counter()
    for row in read(folder / "fold_assignments.tsv"):
        donor, fold = row["donor_id"], int(row["outer_fold"])
        if row["metacell_id"] in metacells or (donor in roster and roster[donor] != fold):
            raise ValueError("Duplicate aggregate or donor leakage")
        if (row["cohort"], row["lineage"]) != (cohort, lineage):
            raise ValueError("Fold context mismatch")
        metacells.add(row["metacell_id"])
        roster[donor] = fold
        counts[donor] += 1
    roster = {donor: (fold, counts[donor]) for donor, fold in roster.items()}
    eligible = {row["donor_id"] for row in read(root / f"{folder.name}_eligibility.tsv")
                if row["eligible"].lower() == "true"}
    if set(roster) != eligible:
        raise ValueError("Eligible donor roster differs from evaluated roster")
    endpoints = ["promoter_accessibility", "rna"] if cohort == "GSE296875" else ["promoter_accessibility"]
    totals, genes, coverage, nrows = reduce_predictions(folders[0], roster, cohort, lineage, endpoints)
    recorded_donors = unique_rows(folder / "donor_errors.tsv", ("endpoint", "method", "baseline", "donor_id"))
    recorded = unique_rows(folder / "comparisons.tsv", ("endpoint", "method", "baseline"))
    results, donor_rows = [], []
    expected_donor_keys = set()
    for endpoint in endpoints:
        for method, baseline in contrasts(endpoint):
            key = endpoint, method, baseline
            donor_values = []
            matched_candidate, matched_baseline = [], []
            for donor, (total, count, candidate, reference) in sorted(totals[key].items()):
                donor_key = (*key, donor)
                expected_donor_keys.add(donor_key)
                old = recorded_donors[donor_key]
                check_number(total / count, old["mean_loss_difference"], str(donor_key))
                check_number(count, old["n_common_genes"], str(donor_key))
                donor_values.append(total / count)
                matched_candidate.append(candidate / count)
                matched_baseline.append(reference / count)
                donor_rows.append(dict(cohort=cohort, lineage=lineage, endpoint=endpoint, method=method,
                                       baseline=baseline, donor_id=donor, mean_loss_difference=total/count,
                                       n_common_genes=count))
            seed = (int(settings["seed"]) + int(hashlib.sha256(str(key).encode()).hexdigest()[:8], 16)) % 2**32
            stats = inference(donor_values, seed, int(settings["bootstrap_reps"]))
            stats.update(n_common_predictions=sum(v[1] for v in totals[key].values()), n_common_genes=len(genes[key]))
            for field, value in stats.items():
                if field == "status":
                    if value != recorded[key][field]:
                        raise ValueError("Inference status mismatch")
                else:
                    check_number(value, recorded[key][field], str((key, field)))
            results.append(dict(cohort=cohort, lineage=lineage, endpoint=endpoint, method=method, baseline=baseline,
                                matched_candidate_loss=float(np.mean(matched_candidate)) if matched_candidate else float("nan"),
                                matched_baseline_loss=float(np.mean(matched_baseline)) if matched_baseline else float("nan"),
                                min_common_genes_per_donor=min((v[1] for v in totals[key].values()), default=0),
                                max_common_genes_per_donor=max((v[1] for v in totals[key].values()), default=0),
                                **score_influence(donor_values), **stats))
    if set(recorded_donors) != expected_donor_keys or set(recorded) != {tuple(r[k] for k in ("endpoint", "method", "baseline")) for r in results}:
        raise ValueError("Unexpected recorded donor/contrast keys")
    return results, donor_rows, [dict(cohort=cohort, lineage=lineage, endpoint=e, method=m, prediction_rows=n) for (e,m),n in sorted(coverage.items())], dict(root=str(root), prediction_rows=nrows, donors=len(roster), metacells=len(metacells))


def write(path, rows):
    if not rows:
        return
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pooled", type=Path)
    parser.add_argument("--require-full-family", action="store_true")
    args = parser.parse_args()
    comparisons, donors, coverage, contexts = [], [], [], []
    for root in args.network:
        result, d, c, context = review(root)
        comparisons.extend(result); donors.extend(d); coverage.extend(c); contexts.append(context)
    fields = ("cohort", "lineage", "endpoint", "method", "baseline")
    indexed = {tuple(row[f] for f in fields): row for row in comparisons}
    if len(indexed) != len(comparisons) or set(indexed) - set(planned()):
        raise ValueError("Duplicate or unplanned comparison")
    if args.require_full_family and set(indexed) != set(planned()):
        raise ValueError("Full family not complete")
    keys = list(planned())
    for key, adjusted in zip(keys, holm([indexed[k]["pvalue"] if k in indexed else 1. for k in keys])):
        if key in indexed:
            indexed[key]["holm_pvalue"] = adjusted
            indexed[key]["family_size"] = 48
            indexed[key]["predictive_improvement"] = adjusted < .05 and indexed[key]["ci_low"] > 0
    if args.pooled:
        pooled = {tuple(row[f] for f in fields): row for row in read(args.pooled)}
        for key, row in indexed.items():
            check_number(row["holm_pvalue"], pooled[key]["holm_pvalue"], str(key))
            if str(row["predictive_improvement"]).lower() != pooled[key]["predictive_improvement"].lower():
                raise ValueError("Pooled improvement mismatch")
    args.output.mkdir(parents=True, exist_ok=False)
    write(args.output / "independent_comparisons.tsv", comparisons)
    write(args.output / "independent_donor_differences.tsv", donors)
    write(args.output / "method_coverage.tsv", coverage)
    (args.output / "verification.json").write_text(json.dumps(dict(status="passed", contexts=contexts,
        comparisons=len(comparisons), full_family=len(indexed)==48, numpy=np.__version__), indent=2)+"\n")
    (args.output / "REPORT.txt").write_text(
        "Prediction table verification passed for the listed completed contexts.\n"
        "Losses are donor-mean metacell MSE divided by training target variance; positive differences favor the candidate.\n"
        "Genes contribute to each donor mean and are not independent inferential replicates. Coverage is contrast matched.\n"
        "Bootstrap uncertainty conditions on the fitted cross-validation predictions, not model retraining. Sign flips require independent donors and symmetric donor differences.\n"
        "Leave-one-donor-out means omit fixed donor scores only; they do not refit selection, tuning or prediction models and do not alter success criteria.\n"
        "The fixed Holm family contains 48 comparisons; absent contexts have P=1 and provisional adjustments.\n"
        "GSE296875 stellate denotes source-dependent Mesenchymal context, not a purified stellate population.\n"
        "Small positive values alone do not establish material predictive gains or causal regulation.\n\n" +
        "\n".join(
            f"{r['cohort']} {r['lineage']} {r['endpoint']} {r['method']} vs {r['baseline']}: "
            f"donors={r['n_donors']}, matched_genes={r['n_common_genes']}, "
            f"donor_gene_range={r['min_common_genes_per_donor']}..{r['max_common_genes_per_donor']}, "
            f"standardized_MSE_gain={r['improvement']:.8g} "
            f"CI=[{r['ci_low']:.8g},{r['ci_high']:.8g}], "
            f"matched_baseline_MSE={r['matched_baseline_loss']:.8g}, "
            f"fraction_positive_donors={r['fraction_positive_donor_gains']:.6g}, "
            f"donor_gain_range=[{r['min_donor_gain']:.8g},{r['max_donor_gain']:.8g}], "
            f"fixed_score_LODO_mean_range=[{r['min_leave_one_donor_out_mean_gain']:.8g},"
            f"{r['max_leave_one_donor_out_mean_gain']:.8g}], "
            f"Holm48_P={r['holm_pvalue']:.8g}, supported={r['predictive_improvement']}"
            for r in comparisons) + "\n")


if __name__ == "__main__":
    main()
