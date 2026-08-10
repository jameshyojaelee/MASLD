#!/usr/bin/env python3
"""Run exact biological-unit label/sign permutations for the claim-bearing contrasts."""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd

from public_functional_common import CANDIDATE_ROOT, read_tsv, require_sealed, write_tsv


FIELDS = [
    "dataset_id", "contrast_id", "program_uid", "lineage", "permutation_unit",
    "permutation_scheme", "n_biological", "n_permutations", "observed", "null_mean",
    "null_sd", "empirical_p", "null_center_tolerance", "null_centered", "specification_sha256",
]


def exact_sign_flip(values: np.ndarray) -> tuple[np.ndarray, float]:
    signs = np.asarray(list(itertools.product([-1.0, 1.0], repeat=len(values))))
    null = (signs * values).mean(axis=1)
    observed = float(values.mean())
    p = float((np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1))
    return null, p


def pcls(primary_uids: set[str], specification: str) -> list[dict[str, object]]:
    scores = pd.read_csv(CANDIDATE_ROOT / "analyses/GSE200418/per_sample_program_scores.tsv", sep="\t")
    scores = scores[(scores.scoring_scheme == "weighted") & scores.program_uid.isin(primary_uids)]
    metadata = pd.DataFrame([row for row in read_tsv(CANDIDATE_ROOT / "sample_manifest.tsv") if row["dataset_id"] == "GSE200418"])
    metadata["condition_time"] = metadata.condition + "_" + metadata.timepoint
    merged = scores.merge(metadata[["sample_id", "biological_unit_id", "condition_time"]], on="sample_id", validate="many_to_one")
    rows = []
    for uid, frame in merged.groupby("program_uid"):
        differences = []
        for donor, donor_frame in frame.groupby("biological_unit_id"):
            numerator = donor_frame.loc[donor_frame.condition_time == "GFIPO_48h", "program_score"]
            denominator = donor_frame.loc[donor_frame.condition_time == "GFI_48h", "program_score"]
            if len(numerator) and len(denominator):
                differences.append(float(numerator.mean() - denominator.mean()))
        values = np.asarray(differences)
        null, p = exact_sign_flip(values)
        tolerance = max(1e-12, 0.01 * float(null.std(ddof=1)))
        rows.append(record("GSE200418", "PCLS_GFIPO_GFI_48h", uid, "", "donor", "exact_within_donor_condition_sign_flip", values, null, p, tolerance, specification))
    return rows


def biopsy(primary_uids: set[str], specification: str) -> list[dict[str, object]]:
    path = CANDIDATE_ROOT / "analyses/GSE106737/per_sample_program_scores.tsv"
    if not path.is_file():
        return []
    scores = pd.read_csv(path, sep="\t")
    scores = scores[(scores.scoring_scheme == "weighted") & scores.program_uid.isin(primary_uids)]
    metadata = pd.DataFrame([row for row in read_tsv(CANDIDATE_ROOT / "sample_manifest.tsv") if row["dataset_id"] == "GSE106737"])
    merged = scores.merge(metadata[["sample_id", "biological_unit_id", "condition", "timepoint"]], on="sample_id", validate="many_to_one")
    rows = []
    for uid, frame in merged.groupby("program_uid"):
        deltas = []
        labels = []
        for (participant, condition), participant_frame in frame.groupby(["biological_unit_id", "condition"]):
            if condition not in {"lifestyle_responder", "lifestyle_nonresponder"}:
                continue
            before = participant_frame.loc[participant_frame.timepoint == "baseline", "program_score"]
            after = participant_frame.loc[participant_frame.timepoint == "followup", "program_score"]
            if len(before) == 1 and len(after) == 1:
                deltas.append(float(after.iloc[0] - before.iloc[0]))
                labels.append(condition)
        values = np.asarray(deltas)
        n_responder = labels.count("lifestyle_responder")
        observed_mask = np.asarray([label == "lifestyle_responder" for label in labels])
        observed = float(values[observed_mask].mean() - values[~observed_mask].mean())
        null = []
        for selected in itertools.combinations(range(len(values)), n_responder):
            mask = np.zeros(len(values), dtype=bool)
            mask[list(selected)] = True
            null.append(values[mask].mean() - values[~mask].mean())
        null_array = np.asarray(null)
        p = float((np.sum(np.abs(null_array) >= abs(observed)) + 1) / (len(null_array) + 1))
        tolerance = max(1e-12, 0.01 * float(null_array.std(ddof=1)))
        rows.append(record("GSE106737", "LSI_RESPONSE_DIFF", uid, "", "participant", "exact_10_of_20_response_label_permutation", values, null_array, p, tolerance, specification, observed=observed))
    return rows


def schlo(primary_uids: set[str], specification: str) -> list[dict[str, object]]:
    scores = pd.read_csv(CANDIDATE_ROOT / "analyses/GSE207889/per_sample_program_scores.tsv", sep="\t")
    scores = scores[(scores.scoring_scheme == "weighted") & scores.program_uid.isin(primary_uids)]
    metadata = pd.read_csv(CANDIDATE_ROOT / "preprocessed/GSE207889/pseudobulk_manifest.tsv", sep="\t")
    merged = scores.merge(metadata[["pseudobulk_id", "condition", "replicate", "lineage"]], on=["pseudobulk_id", "lineage"], validate="many_to_one")
    definitions = {"SCHLO_PA_CONTROL": ("PA", "CONTROL_PA"), "SCHLO_OA_CONTROL": ("OA", "CONTROL_OA"), "SCHLO_TGFB_CONTROL": ("TGFB", "CONTROL_TGFB")}
    rows = []
    for (uid, lineage), frame in merged.groupby(["program_uid", "lineage"]):
        for contrast_id, (numerator, denominator) in definitions.items():
            differences = []
            for replicate, rep_frame in frame.groupby("replicate"):
                left = rep_frame.loc[rep_frame.condition == numerator, "program_score"]
                right = rep_frame.loc[rep_frame.condition == denominator, "program_score"]
                if len(left) == 1 and len(right) == 1:
                    differences.append(float(left.iloc[0] - right.iloc[0]))
            values = np.asarray(differences)
            if len(values) != 2:
                raise RuntimeError(f"Expected two scHLO differences for {uid}/{lineage}/{contrast_id}")
            null, p = exact_sign_flip(values)
            tolerance = max(1e-12, 0.01 * float(null.std(ddof=1)))
            rows.append(record("GSE207889", contrast_id, uid, lineage, "source_replicate", "exact_within_replicate_treatment_sign_flip", values, null, p, tolerance, specification))
    return rows


def record(dataset: str, contrast: str, uid: str, lineage: str, unit: str, scheme: str, values: np.ndarray, null: np.ndarray, p: float, tolerance: float, specification: str, observed: float | None = None) -> dict[str, object]:
    null_mean = float(null.mean())
    return {
        "dataset_id": dataset, "contrast_id": contrast, "program_uid": uid, "lineage": lineage,
        "permutation_unit": unit, "permutation_scheme": scheme, "n_biological": len(values),
        "n_permutations": len(null), "observed": float(values.mean()) if observed is None else observed,
        "null_mean": null_mean, "null_sd": float(null.std(ddof=1)), "empirical_p": p,
        "null_center_tolerance": tolerance, "null_centered": str(abs(null_mean) <= tolerance).lower(),
        "specification_sha256": specification,
    }


def main() -> None:
    seal = require_sealed()
    primary_uids = {row["program_uid"] for row in read_tsv(CANDIDATE_ROOT / "frozen_inputs/external_test_programs.tsv")}
    rows = pcls(primary_uids, str(seal["specification_sha256"])) + schlo(primary_uids, str(seal["specification_sha256"])) + biopsy(primary_uids, str(seal["specification_sha256"]))
    if any(row["null_centered"] != "true" for row in rows):
        raise RuntimeError("At least one biological-unit null is not centered")
    write_tsv(CANDIDATE_ROOT / "null_permutation_audit.tsv", rows, FIELDS)
    print(f"NULL_PERMUTATIONS_COMPLETE\t{len(rows)} rows")


if __name__ == "__main__":
    main()
