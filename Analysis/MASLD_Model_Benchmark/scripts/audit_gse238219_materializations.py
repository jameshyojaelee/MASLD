#!/usr/bin/env python3
"""Aggregate frozen GSE238219 deposited-sample materialization receipts."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

from masld_bench.artifacts import verify_frozen_tree


SAMPLES = ("GSM7660623", "GSM7660624", "GSM7660625", "GSM7660626", "GSM7660627")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution-root", type=Path, required=True)
    parser.add_argument("--array-job-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise RuntimeError("refusing to overwrite cross-sample audit")
    summaries = []
    for index, sample in enumerate(SAMPLES):
        source = args.execution_root / f"gse238219-pseudobulk-{sample}-{args.array_job_id}_{index}"
        verify_frozen_tree(source)
        summaries.append(
            json.loads((source / "materialization_summary.json").read_text(encoding="utf-8"))
        )
    source_hashes = {summary["source_sha256"] for summary in summaries}
    gene_hashes = {summary["gene_axis_sha256"] for summary in summaries}
    if len(source_hashes) != 1 or len(gene_hashes) != 1:
        raise RuntimeError("source or gene axis differs across deposited samples")
    totals: Counter[str] = Counter()
    for summary in summaries:
        totals.update(summary["cell_assignment_counts"])
    result = {
        "schema_version": "masld-bench-gse238219-cross-sample-gate-v1",
        "source_sha256": next(iter(source_hashes)),
        "gene_axis_sha256": next(iter(gene_hashes)),
        "gene_axis_identical_across_samples": True,
        "guide_reference_identical_across_samples": len(
            {summary["guide_reference_source_sha256"] for summary in summaries}
        ) == 1,
        "deposited_sample_count": len(summaries),
        "deposited_series_labels": [summary["deposited_series_label"] for summary in summaries],
        "per_sample": summaries,
        "cell_assignment_totals": dict(sorted(totals.items())),
        "independent_biological_unit_count": "UNRESOLVED",
        "reason": "GEO exposes five deposited sample objects, but the source paper reports differential expression over two Perturb-seq replicates. Deposit, culture, transduction, capture, and sequencing ancestry remain unresolved.",
        "allowed_use": "raw-count preprocessing and descriptive topology only",
        "prohibited_use": "five-sample biological confidence intervals, five-fold biological replication claims, or model ranking before ancestry resolution",
        "source_perturbation_effect_tables_read": False,
        "gse313774_accessed": False,
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "cross_sample_gate.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
