#!/usr/bin/env python3
"""Combine completed candidate comparisons without selecting a method or changing claims."""
import argparse
import csv
import json
import math
from pathlib import Path


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def holm(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    result, previous = [1.0] * len(values), 0.0
    for rank, index in enumerate(order):
        previous = max(previous, min(1.0, values[index] * (len(values) - rank)))
        result[index] = previous
    return result


def planned_comparisons():
    for cohort in ("GSE244832", "GSE281367", "GSE296875"):
        endpoints = ("promoter_accessibility", "rna") if cohort == "GSE296875" else ("promoter_accessibility",)
        for lineage in ("hepatocyte", "stellate"):
            for endpoint in endpoints:
                base = "promoter_only" if endpoint == "rna" else "intercept"
                for method, baseline in (("distance", base), ("cicero", base),
                                         ("cicero", "distance"), ("directnet_adapted", base),
                                         ("directnet_adapted", "distance"), ("directnet_adapted", "cicero")):
                    yield (cohort, lineage, endpoint, method, baseline)


def source_is_technical(root):
    settings_path = root / "run_settings.json"
    completion = json.loads((root / "completion.json").read_text())
    if not settings_path.exists():
        return True
    settings = json.loads(settings_path.read_text())
    return bool(settings.get("pilot", True) or settings.get("source_pilot", False)
                or settings.get("analysis_scope", "") == "bounded_gene_development"
                or completion.get("technical_pilot", False)
                or "technical_pilot" in completion.get("interpretation", ""))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--da", required=True, type=Path)
    parser.add_argument("--network", action="append", default=[], type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--technical-pilot", action="store_true")
    args = parser.parse_args()
    if not (args.da / "run_summary.json").exists():
        raise ValueError("DA analysis has not completed")
    by_key = {}
    keys = ["cohort", "lineage", "endpoint", "method", "baseline"]
    for root in args.network:
        if not (root / "completion.json").exists():
            raise ValueError(f"Network analysis has not completed: {root}")
        completion = json.loads((root / "completion.json").read_text())
        if completion.get("status") != "completed" or any(root.glob("INVALIDATED*.json")):
            raise ValueError(f"Network analysis is incomplete or invalidated: {root}")
        technical = source_is_technical(root)
        path = root / "comparisons.tsv"
        if not path.exists():
            continue
        for row in read_tsv(path):
            key = tuple(row[field] for field in keys)
            if key in by_key:
                raise ValueError(f"Multiple results for the same prespecified comparison: {key}")
            row["source"] = str(path.resolve())
            row["source_technical_pilot"] = technical
            by_key[key] = row
    planned = list(planned_comparisons())
    if set(by_key) - set(planned):
        raise ValueError("Unexpected comparison outside the prespecified family")
    rows, pvalues = [], []
    for key in planned:
        row = {**dict(zip(keys, key)), **by_key.get(key, {})}
        value = row.get("pvalue", "")
        p = float(value) if value not in ("", "nan", "NA") else 1.0
        if not math.isfinite(p) or not 0 <= p <= 1:
            raise ValueError("Invalid probability")
        row["pvalue"] = p
        row.setdefault("status", "not_run")
        row.setdefault("source", "")
        row["lineage_state"] = "source_dependent_mesenchymal_context" if key[:2] == ("GSE296875", "stellate") else "source_label"
        rows.append(row)
        pvalues.append(p)
    family_complete = all(key in by_key and by_key[key].get("status") == "evaluated"
                          and by_key[key].get("pvalue", "") not in ("", "nan", "NA")
                          for key in planned)
    for row, adjusted in zip(rows, holm(pvalues)):
        row["source_family"] = row.get("family", "")
        row["source_family_complete"] = row.get("family_complete", "")
        row["family"] = "prespecified_48_cohort_lineage_endpoint_method_baseline_comparisons"
        row["family_complete"] = family_complete
        row["holm_pvalue"] = adjusted
        row["family_size"] = len(planned)
        row["technical_pilot"] = bool(args.technical_pilot or row.get("source_technical_pilot", True))
        lower = row.get("ci_low", "")
        lower = float(lower) if lower not in ("", "nan", "NA") else -math.inf
        row["predictive_improvement"] = bool(not row["technical_pilot"] and adjusted < 0.05 and lower > 0)
        row["interpretation"] = "conditional_on_fixed_donor_cross_validation_predictions_not_causal"
    args.output.mkdir(parents=True, exist_ok=False)
    fields = keys + sorted(set().union(*(set(row) for row in rows)) - set(keys))
    with (args.output / "network_comparisons.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    da_summary = read_tsv(args.da / "cohort_summary.tsv")
    reproduction = read_tsv(args.da / "baseline_reproduction.tsv")
    payload = {"candidate_only": True, "technical_pilot": bool(args.technical_pilot or any(
                   row["technical_pilot"] for row in rows if row["source"])),
               "da_root": str(args.da.resolve()), "da_summary_rows": len(da_summary),
               "baseline_checks": len(reproduction), "network_family_size": len(planned),
               "network_family_complete": family_complete,
               "network_comparisons_evaluated": sum(row["status"] == "evaluated" for row in rows),
               "network_comparisons_not_run_or_untestable": sum(row["status"] != "evaluated" for row in rows),
               "scientific_adoption": False, "canonical_manuscript_changed": False}
    (args.output / "summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload))


if __name__ == "__main__":
    main()
