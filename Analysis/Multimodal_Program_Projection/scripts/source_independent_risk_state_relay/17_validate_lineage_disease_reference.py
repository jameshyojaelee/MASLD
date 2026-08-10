#!/usr/bin/env python3
"""Independently validate the Plan 45 lineage disease-reference release."""

from __future__ import annotations

import csv
import json
import math
import os
from collections import defaultdict
from pathlib import Path

import numpy as np

from relay_common import PROJECT_ROOT, sha256_file


CANDIDATE_ID = os.environ.get(
    "PLAN45_DISEASE_REFERENCE_CANDIDATE_ID",
    "source-independent-risk-state-relay-lineage-disease-reference-2026-08-10",
).strip()
CANDIDATE_ROOT = (
    PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates" / CANDIDATE_ID
)
PRIMARY_LINEAGES = {
    "Hepatocytes",
    "Macrophages",
    "Fibroblasts",
    "Cholangiocytes",
    "Endothelial_cells",
}
PAIRING_FILES = {
    "GSE244832": "data/GSE244832/metadata/donor_pairing.csv",
    "GSE185477": "data/GSE185477/metadata/donor_pairing.csv",
    "GSE202379": "data/GSE202379/metadata/donor_pairing.csv",
    "GSE136103": "data/GSE136103/metadata/donor_pairing.csv",
}
CELL_COLUMNS = {
    "Hepatocytes": "n_Hepatocytes",
    "Macrophages": "n_Macrophages",
    "Fibroblasts": "n_Fibroblasts",
    "Cholangiocytes": "n_Cholangiocytes",
    "Endothelial_cells": "n_Endothelial_cells",
}


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def comma_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle)


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "t", "1", "yes", "y"}


def parse_float(value: str) -> float:
    return float(value) if value.strip() else math.nan


def bh(values: list[float]) -> list[float]:
    out = [math.nan] * len(values)
    valid = [(index, value) for index, value in enumerate(values) if math.isfinite(value)]
    if not valid:
        return out
    ordered = sorted(valid, key=lambda item: item[1])
    m = len(ordered)
    running = 1.0
    for rank_index in range(m - 1, -1, -1):
        original_index, value = ordered[rank_index]
        rank = rank_index + 1
        running = min(running, value * m / rank)
        out[original_index] = min(1.0, running)
    return out


def build_run_to_donor() -> dict[str, str]:
    mapping: dict[str, str] = {}
    for dataset, relative in PAIRING_FILES.items():
        for row in comma_rows(PROJECT_ROOT / relative):
            donor = f"{dataset}_{row['donor_id']}"
            for run in row["rna_srrs"].split(";"):
                run = run.strip()
                if not run:
                    continue
                previous = mapping.setdefault(run, donor)
                if previous != donor:
                    raise RuntimeError(f"Run maps to two donors: {run}")
    return mapping


def validate_exact_lineage_counts(design_path: Path) -> None:
    metadata_path = (
        PROJECT_ROOT
        / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
    )
    run_to_donor = build_run_to_donor()
    expected: dict[tuple[str, str], float] = defaultdict(float)
    with metadata_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            donor = run_to_donor.get(row["sample"], row["sample"])
            for lineage, column in CELL_COLUMNS.items():
                text = row[column].strip()
                expected[(donor, lineage)] += float(text) if text else 0.0
    observed_keys = set()
    for row in rows(design_path):
        key = (row["donor"], row["cell_type"])
        observed_keys.add(key)
        observed = float(row["n_lineage_cells"])
        if not math.isclose(observed, expected[key], rel_tol=0, abs_tol=1e-8):
            raise RuntimeError(f"Exact lineage-cell count mismatch: {key}")
        expected_primary = (
            row["disease_stage_coarse"] in {"Healthy", "Steatosis", "Steatohepatitis"}
            and not parse_bool(row["exclude_stage_analysis"])
            and row["dataset"] != "GSE189600"
            and observed >= 50
            and parse_bool(row["count_matrix_present"])
            and parse_float(row["library_size"]) > 0
        )
        if parse_bool(row["primary_eligible"]) != expected_primary:
            raise RuntimeError(f"Primary design gate mismatch: {key}")
    if not set(expected).issubset(observed_keys):
        raise RuntimeError("Donor-lineage design omits metadata donors")


def validate_bh(coefficients_path: Path) -> tuple[int, int]:
    families: dict[tuple[str, str], list[tuple[str, float, float, float, float]]] = defaultdict(list)
    row_count = 0
    for row in rows(coefficients_path):
        row_count += 1
        families[(row["cell_type"], row["model_id"])].append(
            (
                row["gene_symbol"],
                parse_float(row["pvalue"]),
                parse_float(row["qvalue"]),
                parse_float(row["treat_pvalue"]),
                parse_float(row["treat_qvalue"]),
            )
        )
    for family, records in families.items():
        ordinary = bh([record[1] for record in records])
        interval = bh([record[3] for record in records])
        for record, expected_q, expected_treat_q in zip(records, ordinary, interval):
            gene, _, observed_q, _, observed_treat_q = record
            if math.isfinite(expected_q) and not math.isclose(
                observed_q, expected_q, rel_tol=2e-7, abs_tol=2e-10
            ):
                raise RuntimeError(f"BH mismatch: {family} / {gene}")
            if math.isfinite(expected_treat_q) and not math.isclose(
                observed_treat_q, expected_treat_q, rel_tol=2e-7, abs_tol=2e-10
            ):
                raise RuntimeError(f"TREAT BH mismatch: {family} / {gene}")
    return row_count, len(families)


def validate_simple_bh(path: Path, family_columns: tuple[str, ...]) -> tuple[int, int]:
    families: dict[tuple[str, ...], list[tuple[str, float, float]]] = defaultdict(list)
    row_count = 0
    for row in rows(path):
        row_count += 1
        key = tuple(row[column] for column in family_columns)
        families[key].append(
            (row["gene_symbol"], parse_float(row["pvalue"]), parse_float(row["qvalue"]))
        )
    for family, records in families.items():
        expected = bh([record[1] for record in records])
        for record, expected_q in zip(records, expected):
            gene, _, observed_q = record
            if math.isfinite(expected_q) and not math.isclose(
                observed_q, expected_q, rel_tol=2e-7, abs_tol=2e-10
            ):
                raise RuntimeError(f"BH mismatch: {family} / {gene}")
    return row_count, len(families)


def validate_reference(
    reference_path: Path, lodo_path: Path
) -> tuple[int, int, dict[str, int], dict[str, int]]:
    lodo: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows(lodo_path):
        if parse_bool(row["estimable"]):
            value = parse_float(row["logFC"])
            if math.isfinite(value):
                lodo[(row["cell_type"], row["gene_symbol"])].append(value)
    rows_seen = 0
    eligible = 0
    eligible_by_lineage: dict[str, int] = defaultdict(int)
    strict_by_lineage: dict[str, int] = defaultdict(int)
    ordinary_fdr_by_lineage: dict[str, int] = defaultdict(int)
    within_direction_by_lineage: dict[str, int] = defaultdict(int)
    fdr_and_within_by_lineage: dict[str, int] = defaultdict(int)
    lineages = set()
    for row in rows(reference_path):
        rows_seen += 1
        lineages.add(row["cell_type"])
        key = (row["cell_type"], row["gene_symbol"])
        values = lodo.get(key, [])
        beta = parse_float(row["logFC"])
        n_estimable = int(row["n_lodo_estimable"])
        if n_estimable != len(values):
            raise RuntimeError(f"LODO count mismatch: {key}")
        if values and math.isfinite(beta) and beta != 0:
            same = sum(value * beta > 0 for value in values)
            agreement = same / len(values)
            observed = parse_float(row["lodo_direction_agreement_fraction"])
            if not math.isclose(observed, agreement, rel_tol=0, abs_tol=1e-12):
                raise RuntimeError(f"LODO direction mismatch: {key}")
        nuisance = (
            parse_bool(row["static_mitochondrial"])
            or parse_bool(row["static_ribosomal"])
            or parse_bool(row["static_annotation_identity"])
        )
        if parse_bool(row["static_nuisance"]) != nuisance:
            raise RuntimeError(f"Nuisance union mismatch: {key}")
        agreement = parse_float(row["lodo_direction_agreement_fraction"])
        expected_eligible = (
            parse_bool(row["estimable"])
            and math.isfinite(beta)
            and beta != 0
            and not nuisance
            and n_estimable >= 3
            and math.isfinite(agreement)
            and agreement >= 0.75
        )
        if parse_bool(row["primary_reference_eligible"]) != expected_eligible:
            raise RuntimeError(f"Reference eligibility mismatch: {key}")
        loading = parse_float(row["primary_raw_loading"])
        if expected_eligible:
            eligible += 1
            eligible_by_lineage[row["cell_type"]] += 1
            if not math.isclose(loading, beta, rel_tol=0, abs_tol=1e-12):
                raise RuntimeError(f"Primary loading mismatch: {key}")
        elif loading != 0:
            raise RuntimeError(f"Ineligible reference has nonzero loading: {key}")
        if parse_bool(row["strict_treat_reference_eligible"]):
            strict_by_lineage[row["cell_type"]] += 1
        ordinary_fdr = (
            expected_eligible
            and math.isfinite(parse_float(row["qvalue"]))
            and parse_float(row["qvalue"]) < 0.05
        )
        within_direction = (
            expected_eligible
            and parse_bool(row["within_dataset_direction_sensitivity_pass"])
        )
        if ordinary_fdr:
            ordinary_fdr_by_lineage[row["cell_type"]] += 1
        if within_direction:
            within_direction_by_lineage[row["cell_type"]] += 1
        if ordinary_fdr and within_direction:
            fdr_and_within_by_lineage[row["cell_type"]] += 1
    if lineages != PRIMARY_LINEAGES:
        raise RuntimeError("Reference lineage universe drift")
    if eligible == 0:
        raise RuntimeError("Reference contains no eligible loadings")
    return (
        rows_seen,
        eligible,
        dict(eligible_by_lineage),
        dict(strict_by_lineage),
        dict(ordinary_fdr_by_lineage),
        dict(within_direction_by_lineage),
        dict(fdr_and_within_by_lineage),
    )


def validate_design_audit(gate_path: Path, design_path: Path) -> None:
    design_rows = list(rows(design_path))
    for gate in rows(gate_path):
        lineage = gate["cell_type"]
        selected = [
            row for row in design_rows
            if row["cell_type"] == lineage and parse_bool(row["primary_eligible"])
        ]
        if len(selected) != int(gate["n_donors"]):
            raise RuntimeError(f"Primary donor count mismatch: {lineage}")
        if not parse_bool(gate["source_gate_pass"]):
            raise RuntimeError(f"Frozen source gate did not pass: {lineage}")
        datasets = sorted({row["dataset"] for row in selected})
        matrix = []
        for row in selected:
            matrix.append(
                [1.0]
                + [float(row["dataset"] == value) for value in datasets[1:]]
                + [float(row["disease_binary"])]
            )
        x = np.asarray(matrix, dtype=float)
        rank = int(np.linalg.matrix_rank(x))
        condition = float(np.linalg.cond(x))
        if rank != int(gate["design_rank"]) or x.shape[1] != int(gate["n_model_columns"]):
            raise RuntimeError(f"Independent design-rank mismatch: {lineage}")
        if not math.isclose(
            condition, float(gate["condition_number"]), rel_tol=1e-8, abs_tol=1e-8
        ):
            raise RuntimeError(f"Independent condition-number mismatch: {lineage}")


def main() -> None:
    seal_path = CANDIDATE_ROOT / "LINEAGE_DISEASE_REFERENCE_SEALED.json"
    required_names = [
        "donor_lineage_design",
        "lineage_disease_coefficients",
        "lineage_design_audit",
        "lineage_disease_lodo",
        "lineage_within_dataset_effects",
        "frozen_lineage_disease_reference",
        "static_nuisance_gene_registry",
        "source_gate_status",
        "donor_metadata_collapse_audit",
        "lineage_disease_input_manifest",
    ]
    required = {name: CANDIDATE_ROOT / f"{name}.tsv" for name in required_names}
    if not seal_path.is_file():
        raise RuntimeError(f"Missing disease-reference seal: {seal_path}")
    for path in required.values():
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing disease-reference artifact: {path}")
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    forbidden = {
        "failed_plan43_or_plan44_score_reused": False,
        "inferred_f_stage_used": False,
        "run_level_rows_used_as_replicates": False,
        "experimental_outcomes_inspected": False,
        "experimental_targets_frozen": False,
    }
    for key, expected in forbidden.items():
        if seal.get(key) is not expected:
            raise RuntimeError(f"Disease-reference firewall violation: {key}")
    for name, path in required.items():
        expected = seal["output_sha256"].get(name)
        if expected != sha256_file(path):
            raise RuntimeError(f"Disease-reference hash mismatch: {name}")

    validate_exact_lineage_counts(required["donor_lineage_design"])
    validate_design_audit(
        required["source_gate_status"], required["donor_lineage_design"]
    )
    coefficient_rows, bh_families = validate_bh(
        required["lineage_disease_coefficients"]
    )
    lodo_rows, lodo_bh_families = validate_simple_bh(
        required["lineage_disease_lodo"], ("cell_type", "held_out_dataset")
    )
    within_rows, within_bh_families = validate_simple_bh(
        required["lineage_within_dataset_effects"], ("cell_type", "dataset")
    )
    (
        reference_rows,
        eligible,
        eligible_by_lineage,
        strict_by_lineage,
        ordinary_fdr_by_lineage,
        within_direction_by_lineage,
        fdr_and_within_by_lineage,
    ) = validate_reference(
        required["frozen_lineage_disease_reference"],
        required["lineage_disease_lodo"],
    )
    print(
        "Plan 45 lineage-disease reference validation passed: "
        f"lineages={len(PRIMARY_LINEAGES)}; coefficient_rows={coefficient_rows}; "
        f"primary_BH_families={bh_families}; lodo_rows={lodo_rows}; "
        f"lodo_BH_families={lodo_bh_families}; within_rows={within_rows}; "
        f"within_BH_families={within_bh_families}; reference_rows={reference_rows}; "
        f"eligible_loadings={eligible}; eligible_by_lineage={eligible_by_lineage}; "
        f"ordinary_FDR_by_lineage={ordinary_fdr_by_lineage}; "
        f"within_direction_by_lineage={within_direction_by_lineage}; "
        f"FDR_and_within_by_lineage={fdr_and_within_by_lineage}; "
        f"strict_TREAT_by_lineage={strict_by_lineage}; targets not frozen"
    )


if __name__ == "__main__":
    main()
