#!/usr/bin/env python3
"""Independently validate the five-cohort lncRNA robustness candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path


COHORTS = {"GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"}
N_ALL = 23_370
N_LNCRNA = 6_249


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_table(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    require(path.is_file() and not path.is_symlink(), f"missing table: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    require(rows and all(None not in row for row in rows), f"malformed table: {path}")
    return rows


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bh(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [1.0] * len(values)
    running = 1.0
    total = len(values)
    for reversed_rank, index in enumerate(reversed(order), start=1):
        rank = total - reversed_rank + 1
        running = min(running, values[index] * total / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def direction(value: float) -> int:
    return int(value > 0) - int(value < 0)


def as_bool(value: str) -> bool:
    require(value in {"TRUE", "FALSE", "true", "false"}, f"invalid boolean: {value}")
    return value.lower() == "true"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    require(
        args.candidate.is_dir() and not args.candidate.is_symlink(), "invalid candidate"
    )
    require(not args.output.exists(), "review output already exists")

    manifest = read_table(args.candidate / "output_manifest.tsv")
    require(len(manifest) == 8, "output manifest cardinality drift")
    for row in manifest:
        path = args.candidate / row["relative_path"]
        require(
            path.stat().st_size == int(row["size_bytes"]), f"size drift: {path.name}"
        )
        require(sha256(path) == row["sha256"], f"checksum drift: {path.name}")

    source = read_table(args.candidate / "source_manifest.tsv")
    source_by_role = {row["role"]: row for row in source}
    require(len(source_by_role) == len(source), "duplicate source roles")
    for row in source:
        path = Path(row["path"])
        require(path.is_file() and not path.is_symlink(), f"invalid source: {path}")
        require(
            path.stat().st_size == int(row["size_bytes"]), f"source size drift: {path}"
        )
        require(sha256(path) == row["sha256"], f"source checksum drift: {path}")

    primary_path = Path(source_by_role["validated_primary_all_gene_result"]["path"])
    primary = read_table(primary_path, delimiter=",")
    require(len(primary) == N_ALL, "all-gene family drift")
    require(len({row["gene"] for row in primary}) == N_ALL, "duplicate all-gene IDs")
    recomputed_bh = bh([float(row["treat_p"]) for row in primary])
    require(
        max(
            abs(value - float(row["treat_fdr"]))
            for value, row in zip(recomputed_bh, primary)
        )
        < 1e-10,
        "all-gene TREAT BH rederivation failed",
    )
    require(
        sum(float(row["treat_fdr"]) < 0.05 for row in primary) == 1_616,
        "TREAT count drift",
    )

    results = read_table(args.candidate / "bulk_lncrna_results.tsv")
    require(len(results) == N_LNCRNA, "lncRNA universe drift")
    require(
        len({row["gene_id_versioned"] for row in results}) == N_LNCRNA,
        "duplicate lncRNA IDs",
    )
    result_by_gene = {row["gene_id_versioned"]: row for row in results}
    primary_by_gene = {row["gene"]: row for row in primary}
    require(
        set(result_by_gene) <= set(primary_by_gene), "lncRNA ID not in all-gene family"
    )
    for gene, row in result_by_gene.items():
        source_row = primary_by_gene[gene]
        for field in ("logFC", "treat_p", "treat_fdr", "t"):
            require(
                math.isclose(
                    float(row[field]),
                    float(source_row[field]),
                    rel_tol=0,
                    abs_tol=1e-12,
                ),
                f"primary statistic drift: {gene} {field}",
            )
    require(
        sum(float(row["treat_fdr"]) < 0.05 for row in results) == 434,
        "lncRNA positive drift",
    )

    cohort = read_table(args.candidate / "bulk_lncrna_cohort_effects.tsv")
    loco = read_table(args.candidate / "bulk_lncrna_leave_one_cohort_out.tsv")
    require(len(cohort) == N_LNCRNA * 5, "cohort family drift")
    require(len(loco) == N_LNCRNA * 5, "leave-one-cohort family drift")
    cohort_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    loco_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in cohort:
        cohort_by_gene[row["gene_id_versioned"]].append(row)
    for row in loco:
        loco_by_gene[row["gene_id_versioned"]].append(row)

    high_confidence: list[str] = []
    for gene, row in result_by_gene.items():
        pooled_direction = direction(float(row["logFC"]))
        cohort_rows = cohort_by_gene[gene]
        loco_rows = loco_by_gene[gene]
        require(
            {item["cohort"] for item in cohort_rows} == COHORTS,
            f"cohort coverage drift: {gene}",
        )
        require(
            {item["excluded_cohort"] for item in loco_rows} == COHORTS,
            f"LOCO coverage drift: {gene}",
        )
        cohort_agree = sum(
            direction(float(item["logFC"])) == pooled_direction for item in cohort_rows
        )
        materially_opposite = sum(
            (pooled_direction > 0 and float(item["CI_high"]) < 0)
            or (pooled_direction < 0 and float(item["CI_low"]) > 0)
            for item in cohort_rows
        )
        loco_agree = sum(
            direction(float(item["logFC"])) == pooled_direction for item in loco_rows
        )
        loco_positive = sum(float(item["treat_fdr"]) < 0.05 for item in loco_rows)
        eligible = (
            row["is_canonical_chromosome"] == "true"
            and row["mapping_status"] == "mapping_unambiguous"
            and row["main_text_eligible"] == "true"
        )
        expected = (
            float(row["treat_fdr"]) < 0.05
            and eligible
            and cohort_agree >= 4
            and materially_opposite == 0
            and loco_agree == 5
            and loco_positive >= 4
            and as_bool(row["equal_weight_direction_agrees"])
            and float(row["equal_weight_treat_fdr"]) < 0.05
        )
        require(
            as_bool(row["high_confidence"]) == expected,
            f"high-confidence drift: {gene}",
        )
        if expected:
            high_confidence.append(gene)

    verdict = {
        row["metric"]: row["value"]
        for row in read_table(args.candidate / "bulk_lncrna_verdict.tsv")
    }
    require(verdict["source_gate"] == "pass", "source gate verdict drift")
    require(int(verdict["tested_all_genes"]) == N_ALL, "all-gene verdict drift")
    require(int(verdict["tested_lncrna"]) == N_LNCRNA, "lncRNA verdict drift")
    require(
        int(verdict["primary_treat_all_genes"]) == 1_616,
        "all-gene positive verdict drift",
    )
    require(
        int(verdict["primary_treat_lncrna"]) == 434, "lncRNA positive verdict drift"
    )
    require(
        int(verdict["high_confidence_lncrna"]) == len(high_confidence),
        "high-confidence verdict drift",
    )

    example = (
        min(
            high_confidence,
            key=lambda gene: (-abs(float(result_by_gene[gene]["t"])), gene),
        )
        if high_confidence
        else None
    )
    args.output.mkdir(parents=True)
    payload = {
        "status": "pass",
        "candidate": str(args.candidate.resolve()),
        "biological_unit": "human_sample",
        "multiple_testing": "TREAT_BH_across_complete_23370_gene_family",
        "n_all_genes": N_ALL,
        "n_lncrna": N_LNCRNA,
        "n_all_treat_positive": 1_616,
        "n_lncrna_treat_positive": 434,
        "n_high_confidence_lncrna": len(high_confidence),
        "deterministic_disease_state_example": example,
        "example_rule": "largest_absolute_pooled_t_then_versioned_Ensembl_ID",
    }
    (args.output / "VALIDATED.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
