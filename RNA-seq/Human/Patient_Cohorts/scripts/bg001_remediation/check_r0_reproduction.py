#!/usr/bin/env python3
"""Hard gate: the candidate runner must reproduce the frozen canonical R0 fit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


FIELDS = ("logFC", "SE", "P.Value", "treat_p")
ABS_TOL = 1e-10


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve(strict=True)
    reference_path = root / "frozen_sets/canonical_deg_results.csv"
    candidate_path = root / "arms/R0/results/integration/deg_results.csv"
    reference = load_csv(reference_path)
    candidate = load_csv(candidate_path)
    failures: list[str] = []
    metrics: dict[str, object] = {}
    frozen_hashes = root / "manifests/baseline_frozen_files.sha256"
    expected_reference_hash = None
    for line in frozen_hashes.read_text().splitlines():
        value, path = line.split(maxsplit=1)
        if Path(path).resolve() == reference_path.resolve():
            expected_reference_hash = value
            break
    observed_reference_hash = sha256(reference_path)
    metrics["reference_sha256"] = observed_reference_hash
    if expected_reference_hash is None or observed_reference_hash != expected_reference_hash:
        failures.append("frozen canonical DEG hash differs from baseline manifest")
    ref_genes = [row["gene"] for row in reference]
    cand_genes = [row["gene"] for row in candidate]
    with (root / "frozen_sets/canonical_genes.tsv").open(newline="") as handle:
        frozen_genes = [row["gene"] for row in csv.DictReader(handle, delimiter="\t")]
    if ref_genes != frozen_genes:
        failures.append("frozen canonical DEG and DGE gene order differ")
    if ref_genes != cand_genes:
        failures.append("gene identity/order differs")
    metrics["reference_genes"] = len(ref_genes)
    metrics["candidate_genes"] = len(cand_genes)
    max_delta = {field: 0.0 for field in FIELDS}
    if len(reference) == len(candidate) and ref_genes == cand_genes:
        for old, new in zip(reference, candidate):
            for field in FIELDS:
                a = float(old[field])
                b = float(new[field])
                if math.isnan(a) and math.isnan(b):
                    continue
                if not math.isfinite(a) or not math.isfinite(b):
                    if a != b:
                        failures.append(f"{field} nonfinite mask differs at {old['gene']}")
                    continue
                delta = abs(a - b)
                max_delta[field] = max(max_delta[field], delta)
                if delta > ABS_TOL and len(failures) < 20:
                    failures.append(f"{field} delta {delta:.3g} at {old['gene']}")
    metrics["max_abs_delta"] = max_delta

    def treat(rows: list[dict[str, str]]) -> dict[str, int]:
        result: dict[str, int] = {}
        for row in rows:
            fdr = float(row["treat_fdr"])
            if math.isfinite(fdr) and fdr < 0.05:
                result[row["gene"]] = 1 if float(row["logFC"]) > 0 else -1
        return result

    old_treat = treat(reference)
    new_treat = treat(candidate)
    metrics["reference_treat"] = len(old_treat)
    metrics["candidate_treat"] = len(new_treat)
    metrics["candidate_up"] = sum(value > 0 for value in new_treat.values())
    metrics["candidate_down"] = sum(value < 0 for value in new_treat.values())
    if old_treat != new_treat:
        failures.append("TREAT membership or direction differs")
    if (len(new_treat), metrics["candidate_up"], metrics["candidate_down"]) != (1918, 1419, 499):
        failures.append("Canonical 1918/1419/499 TREAT contract differs")

    with (root / "frozen_sets/canonical_samples.tsv").open(newline="") as handle:
        frozen_samples = [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")]
    with (root / "arms/R0/results/integration/model_design.tsv").open(newline="") as handle:
        design_reader = csv.DictReader(handle, delimiter="\t")
        design_columns = design_reader.fieldnames or []
        candidate_design = list(design_reader)
    with (root / "frozen_sets/canonical_model_design.tsv").open(newline="") as handle:
        frozen_design_reader = csv.DictReader(handle, delimiter="\t")
        frozen_design_columns = frozen_design_reader.fieldnames or []
        frozen_design = list(frozen_design_reader)
    model_samples = [row["sample_id"] for row in candidate_design]
    metrics["model_samples"] = len(model_samples)
    metrics["design_columns"] = design_columns
    if model_samples != frozen_samples:
        failures.append("Canonical model sample identity/order differs")
    if "group_binaryDisease" not in design_columns:
        failures.append("Disease coefficient absent from design")
    if design_columns != frozen_design_columns:
        failures.append("Model design columns/order differ from frozen baseline")
    elif candidate_design != frozen_design:
        failures.append("Model design values differ from frozen baseline")
    if len(model_samples) != 846 or len(candidate) != 27638:
        failures.append("Canonical R0 dimensions differ")

    report = {"status": "PASS" if not failures else "FAIL", "metrics": metrics, "failures": failures}
    output = root / "comparisons/R0_reproduction.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))
    if failures:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
