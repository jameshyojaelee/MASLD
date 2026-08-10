#!/usr/bin/env python3
"""Validate the non-overwriting Plan 45 lineage-CV BH correction."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path

from relay_common import PROJECT_ROOT, sha256_file


CANDIDATE_ID = os.environ.get(
    "PLAN45_DISEASE_CV_CORRECTION_ID",
    "source-independent-risk-state-relay-lineage-disease-reference-cv-fdr-correction-2026-08-10",
).strip()
ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates" / CANDIDATE_ID


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def bh(pvalues: list[float]) -> list[float]:
    order = sorted(range(len(pvalues)), key=pvalues.__getitem__)
    out = [1.0] * len(pvalues)
    running = 1.0
    m = len(pvalues)
    for reverse_index in range(m - 1, -1, -1):
        original = order[reverse_index]
        running = min(running, pvalues[original] * m / (reverse_index + 1))
        out[original] = min(1.0, running)
    return out


def main() -> None:
    seal_path = ROOT / "LINEAGE_DISEASE_CV_FDR_CORRECTION_SEALED.json"
    effect_path = ROOT / "cv_fold_effects_multiplicity_corrected.tsv"
    gate_path = ROOT / "cv_gate_status_multiplicity_corrected.tsv"
    manifest_path = ROOT / "cv_multiplicity_correction_manifest.tsv"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    for name, path in {
        "cv_fold_effects_multiplicity_corrected": effect_path,
        "cv_gate_status_multiplicity_corrected": gate_path,
        "cv_multiplicity_correction_manifest": manifest_path,
    }.items():
        if sha256_file(path) != seal["output_sha256"][name]:
            raise RuntimeError(f"CV correction hash mismatch: {name}")
    primary = [
        row for row in rows(effect_path)
        if row["heldout_dataset"] == "GSE244832"
        and row["weighting_rule"] == "continuous_all"
    ]
    if len(primary) != 5:
        raise RuntimeError("Corrected primary family is not five lineages")
    expected = bh([float(row["p_one_sided_positive"]) for row in primary])
    expected_by_lineage = {
        row["cell_type"]: value for row, value in zip(primary, expected)
    }
    for row in primary:
        observed = float(row["primary_gse244832_lineage_family_q"])
        if not math.isclose(observed, expected_by_lineage[row["cell_type"]], abs_tol=1e-15):
            raise RuntimeError(f"Corrected BH mismatch: {row['cell_type']}")
    gates = list(rows(gate_path))
    passing = 0
    for row in gates:
        qvalue = expected_by_lineage[row["cell_type"]]
        expected_pass = (
            float(row["gse244832_effect"]) > 0
            and qvalue < 0.05
            and float(row["gse174748_effect"]) > 0
            and row["all_primary_leave_one_positive"].lower() == "true"
            and row["all_weighting_sensitivities_positive"].lower() == "true"
            and int(row["minimum_nonzero_loadings"]) >= 500
        )
        observed_pass = row["multiplicity_corrected_cv_gate_pass"].lower() == "true"
        if expected_pass != observed_pass:
            raise RuntimeError(f"Corrected gate mismatch: {row['cell_type']}")
        passing += int(expected_pass)
    if passing != seal["n_lineages_passing"]:
        raise RuntimeError("Corrected seal pass count mismatch")
    print(
        "Plan 45 lineage-CV multiplicity validation passed: "
        f"BH_family=5; lineages_passing={passing}/5; targets not frozen"
    )


if __name__ == "__main__":
    main()
