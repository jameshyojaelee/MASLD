#!/usr/bin/env python3
"""Validate exact tests and gate arithmetic for Plan 45 lineage CV."""

from __future__ import annotations

import csv
import itertools
import json
import math
import os
from collections import defaultdict
from pathlib import Path

from relay_common import PROJECT_ROOT, sha256_file


CANDIDATE_ID = os.environ.get(
    "PLAN45_DISEASE_CV_CANDIDATE_ID",
    "source-independent-risk-state-relay-lineage-disease-reference-cv-2026-08-10",
).strip()
ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates" / CANDIDATE_ID
LINEAGES = {
    "Hepatocytes", "Macrophages", "Fibroblasts", "Cholangiocytes",
    "Endothelial_cells",
}
DATASETS = {"GSE174748", "GSE244832"}
RULES = {
    "continuous_all", "sign_all", "continuous_top50_abs", "continuous_top25_abs"
}


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def flag(value: str) -> bool:
    return value.strip().lower() == "true"


def exact(scores: list[float], labels: list[int]) -> tuple[float, float, float, int]:
    disease_n = sum(labels)
    observed = (
        sum(score for score, label in zip(scores, labels) if label) / disease_n
        - sum(score for score, label in zip(scores, labels) if not label)
        / (len(labels) - disease_n)
    )
    null = []
    for disease_indices in itertools.combinations(range(len(labels)), disease_n):
        chosen = set(disease_indices)
        null.append(
            sum(scores[index] for index in chosen) / disease_n
            - sum(scores[index] for index in range(len(labels)) if index not in chosen)
            / (len(labels) - disease_n)
        )
    one = sum(value >= observed - 1e-14 for value in null) / len(null)
    two = sum(abs(value) >= abs(observed) - 1e-14 for value in null) / len(null)
    return observed, one, two, len(null)


def main() -> None:
    names = [
        "cv_fold_design", "cv_gene_loadings", "cv_donor_scores",
        "cv_fold_effects", "cv_gate_status", "cv_input_manifest",
    ]
    paths = {name: ROOT / f"{name}.tsv" for name in names}
    seal_path = ROOT / "LINEAGE_DISEASE_CV_SEALED.json"
    if not seal_path.is_file():
        raise RuntimeError("CV seal absent")
    seal = json.loads(seal_path.read_text())
    for key in (
        "experimental_outcomes_inspected", "experimental_targets_frozen",
        "failed_plan43_or_plan44_score_reused",
    ):
        if seal.get(key) is not False:
            raise RuntimeError(f"CV firewall violation: {key}")
    for name, path in paths.items():
        if not path.is_file() or sha256_file(path) != seal["output_sha256"][name]:
            raise RuntimeError(f"CV artifact/hash failure: {name}")

    score_groups: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows(paths["cv_donor_scores"]):
        score_groups[(row["cell_type"], row["heldout_dataset"], row["weighting_rule"])].append(row)
    expected_keys = set(itertools.product(LINEAGES, DATASETS, RULES))
    if set(score_groups) != expected_keys:
        raise RuntimeError("CV donor-score family drift")

    effects = {
        (row["cell_type"], row["heldout_dataset"], row["weighting_rule"]): row
        for row in rows(paths["cv_fold_effects"])
    }
    if set(effects) != expected_keys:
        raise RuntimeError("CV effect family drift")
    for key, group in score_groups.items():
        scores = [float(row["score"]) for row in group]
        labels = [int(row["disease_binary"]) for row in group]
        observed, one, two, n_null = exact(scores, labels)
        row = effects[key]
        checks = (
            (observed, float(row["effect_disease_minus_control"])),
            (one, float(row["p_one_sided_positive"])),
            (two, float(row["p_two_sided"])),
        )
        if any(not math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12) for a, b in checks):
            raise RuntimeError(f"Exact-test mismatch: {key}")
        if n_null != int(row["n_exact_permutations"]):
            raise RuntimeError(f"Exact-null census mismatch: {key}")

    gates = list(rows(paths["cv_gate_status"]))
    if {row["cell_type"] for row in gates} != LINEAGES:
        raise RuntimeError("CV gate lineage drift")
    passing = 0
    for gate in gates:
        lineage = gate["cell_type"]
        lineage_effects = [row for key, row in effects.items() if key[0] == lineage]
        expected = (
            float(gate["gse244832_effect"]) > 0
            and float(gate["gse244832_p_one_sided"]) < 0.05
            and float(gate["gse174748_effect"]) > 0
            and flag(gate["all_primary_leave_one_positive"])
            and flag(gate["all_weighting_sensitivities_positive"])
            and int(gate["minimum_nonzero_loadings"]) >= 500
        )
        if flag(gate["cv_gate_pass"]) != expected:
            raise RuntimeError(f"CV gate arithmetic mismatch: {lineage}")
        passing += int(expected)
        if len(lineage_effects) != 8:
            raise RuntimeError(f"CV lineage effect count mismatch: {lineage}")
    if passing != int(seal["n_lineages_passing"]):
        raise RuntimeError("CV seal pass count mismatch")
    print(
        "Plan 45 lineage-disease cross-validation passed: "
        f"families={len(effects)}; exact tests independently rederived; "
        f"lineages_passing={passing}/{len(LINEAGES)}; targets not frozen"
    )


if __name__ == "__main__":
    main()
