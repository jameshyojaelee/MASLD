#!/usr/bin/env python3
"""Independently rederive the candidate biotype counts and logistic model."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def standardize(values: np.ndarray) -> np.ndarray:
    return (values - values.mean()) / values.std(ddof=1)


def fit_logistic(
    design: np.ndarray, outcome: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    beta = np.zeros(design.shape[1], dtype=float)
    for _ in range(100):
        eta = np.clip(design @ beta, -40.0, 40.0)
        probability = 1.0 / (1.0 + np.exp(-eta))
        weights = probability * (1.0 - probability)
        information = design.T @ (weights[:, None] * design)
        step = np.linalg.solve(information, design.T @ (outcome - probability))
        beta += step
        if np.max(np.abs(step)) < 1e-12:
            break
    else:
        raise RuntimeError("independent logistic fit did not converge")
    eta = np.clip(design @ beta, -40.0, 40.0)
    probability = 1.0 / (1.0 + np.exp(-eta))
    weights = probability * (1.0 - probability)
    covariance = np.linalg.inv(design.T @ (weights[:, None] * design))
    return beta, np.sqrt(np.diag(covariance))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(args.candidate.is_dir(), "candidate directory is missing")
    require(not args.output.exists(), "review output already exists")

    manifest = read_tsv(args.candidate / "output_manifest.tsv")
    for row in manifest:
        path = args.candidate / row["relative_path"]
        require(
            path.stat().st_size == int(row["size_bytes"]), f"size drift: {path.name}"
        )
        require(sha256(path) == row["sha256"], f"checksum drift: {path.name}")

    genes = read_tsv(args.candidate / "biotype_association_gene_covariates.tsv")
    require(len(genes) == 23_370, "gene universe drift")
    require(
        len({row["gene_id_versioned"] for row in genes}) == len(genes), "duplicate gene"
    )

    expected: dict[str, dict[str, int]] = {}
    for row in genes:
        biotype = row["display_biotype"]
        values = expected.setdefault(
            biotype,
            {"n_tested": 0, "n_treat_positive": 0, "n_treat_up": 0, "n_treat_down": 0},
        )
        values["n_tested"] += 1
        positive = row["treat_positive"] == "TRUE"
        values["n_treat_positive"] += int(positive)
        values["n_treat_up"] += int(positive and row["direction"] == "up")
        values["n_treat_down"] += int(positive and row["direction"] == "down")
    observed_counts = read_tsv(args.candidate / "bulk_deg_counts_by_biotype.tsv")
    for row in observed_counts:
        require(
            expected[row["display_biotype"]]
            == {field: int(row[field]) for field in expected[row["display_biotype"]]},
            f"biotype count drift: {row['display_biotype']}",
        )

    model_rows = [
        row for row in genes if row["gene_type"] in {"protein_coding", "lncRNA"}
    ]
    y = np.array([row["treat_positive"] == "TRUE" for row in model_rows], dtype=float)
    is_lncrna = np.array(
        [row["gene_type"] == "lncRNA" for row in model_rows], dtype=float
    )
    abundance = standardize(np.array([float(row["AveExpr"]) for row in model_rows]))
    variability = standardize(
        np.log1p(np.array([float(row["expression_variability"]) for row in model_rows]))
    )
    length = standardize(
        np.log(np.array([float(row["gene_length_bp"]) for row in model_rows]))
    )
    detection = standardize(
        np.array([float(row["cohort_detection_count"]) for row in model_rows])
    )
    design = np.column_stack(
        (np.ones(len(model_rows)), is_lncrna, abundance, variability, length, detection)
    )
    beta, standard_error = fit_logistic(design, y)
    reported = read_tsv(args.candidate / "biotype_deg_association.tsv")
    require(len(reported) == 1, "association result must have one row")
    require(
        math.isclose(beta[1], float(reported[0]["log_odds_ratio"]), abs_tol=1e-10),
        "coefficient drift",
    )
    require(
        math.isclose(
            standard_error[1], float(reported[0]["standard_error"]), abs_tol=1e-8
        ),
        "standard-error drift: "
        f"independent={standard_error[1]:.17g}; "
        f"reported={float(reported[0]['standard_error']):.17g}",
    )

    args.output.mkdir(parents=True)
    payload = {
        "status": "pass",
        "candidate": str(args.candidate.resolve()),
        "n_genes": len(genes),
        "n_model_genes": len(model_rows),
        "independent_log_odds_ratio": beta[1],
        "independent_standard_error": standard_error[1],
        "interpretation": "descriptive_gene_level_association_not_biological_replication",
    }
    (args.output / "VALIDATED.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
