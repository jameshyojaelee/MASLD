#!/usr/bin/env python3
"""Source-only probe: does a non-degenerate converged elastic-net C exist?

Five of six registered grid points for the PCA elastic net failed to converge on
GSE267145 and the sixth collapsed to zero coefficients.  With four
advanced-fibrosis participants in a five-dimensional projection the classes are
completely separable, which is exactly the regime where a logistic fit either
diverges or is shrunk away entirely.

This probe walks a dense regularisation path and records, per configuration,
whether it converged and whether it kept a coefficient.  It reads GSE267145
only.  It fits nothing on GSE49541, reads no external expression, and reads no
outcome from either study beyond the source fibrosis label the source fit
already uses.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.fit_gse267145_fibrosis_transfer_source_models import (
    build_labels,
    build_representation,
    fit_pipeline,
    log2_cpm_complete_axis,
    raw_decision,
    read_tsv,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    root = arguments.benchmark_root
    contract = json.loads(arguments.contract.read_text(encoding="utf-8"))
    inputs = contract["source_inputs"]
    molecular = root / inputs["molecular_path"]
    outcomes = root / inputs["outcomes_path"]
    activation = root / contract["activation"]["path"]

    _, source_axis = read_tsv(molecular / "rna_feature_axis.tsv")
    _, endpoint_rows = read_tsv(outcomes / "participant_endpoints.tsv")
    _, axis_rows = read_tsv(activation / "common_stable_gene_axis.tsv")
    shared_ids = [row["stable_gene_id"] for row in axis_rows]
    source_ids = [row["stable_gene_id"] for row in source_axis]
    lookup = {value: index for index, value in enumerate(source_ids)}
    labels, keep, _ = build_labels(
        endpoint_rows=endpoint_rows, mapping=contract["endpoint_mapping"]
    )
    raw = np.load(molecular / "rna_values.npy", mmap_mode="r", allow_pickle=False)
    columns = [lookup[value] for value in shared_ids]
    log2_shared = log2_cpm_complete_axis(raw)[np.ix_(keep, columns)]
    eligibility = np.asarray(raw, dtype=np.float64)[np.ix_(keep, columns)]
    representation = build_representation(name="gene_median", log2_shared=log2_shared)
    full = np.arange(len(labels), dtype=np.int64)

    path = [
        0.0003, 0.001, 0.003, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05,
        0.08, 0.1, 0.2, 0.3, 0.5, 1.0, 3.0, 10.0,
    ]
    rows = []
    for reducer, feature_count, pca_components in (("pca", 1000, 5),):
        config = {
            "representation": "gene_median",
            "reducer": reducer,
            "classifier": "elastic_net",
        }
        for l1_ratio in (0.0, 0.5, 1.0):
            for c_value in path:
                state = fit_pipeline(
                    model_id=f"gene_median_{reducer}_elastic_net",
                    config=config,
                    representation=representation,
                    eligibility=eligibility,
                    gene_ids=shared_ids,
                    fitting=full,
                    labels=labels,
                    feature_count=feature_count,
                    pca_components=pca_components,
                    hyperparameters={"c": c_value, "l1_ratio": l1_ratio},
                    seed=1701,
                )
                converged = bool(state.pop("_converged"))
                state.pop("_design_columns", None)
                nonzero = int(state["nonzero_coefficients"][0])
                scores = raw_decision(
                    config=config, state=state, representation=representation
                )
                rows.append(
                    {
                        "reducer": reducer,
                        "l1_ratio": l1_ratio,
                        "c": c_value,
                        "converged": converged,
                        "nonzero_coefficients": nonzero,
                        "distinct_full_fit_scores": int(
                            len(np.unique(np.round(scores, 12)))
                        ),
                        "max_abs_coefficient": float(
                            np.max(np.abs(state.get("coef", np.zeros(1))))
                        ),
                        "usable": bool(converged and nonzero > 0),
                    }
                )
    summary = {
        "schema_version": "masld-bench-gse267145-fibrosis-regularization-probe-v1",
        "source_participants": int(len(labels)),
        "advanced_f3_f4": int(np.sum(labels == 1)),
        "external_expression_values_read": False,
        "external_labels_read": False,
        "usable_configurations": {
            reducer: sorted(
                {
                    (row["l1_ratio"], row["c"])
                    for row in rows
                    if row["reducer"] == reducer and row["usable"]
                }
            )
            for reducer in ("pca",)
        },
        "path": rows,
    }
    with arguments.output.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(summary, indent=2, sort_keys=True, default=str) + "\n")
    print(
        json.dumps(
            {
                "pca_usable": len(summary["usable_configurations"]["pca"]),
                "usable_configurations": summary["usable_configurations"],
            },
            indent=2,
            sort_keys=True,
            default=str,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
