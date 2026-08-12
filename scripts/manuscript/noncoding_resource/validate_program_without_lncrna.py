#!/usr/bin/env python3
"""Validate the prespecified leave-all-lncRNA-out program sensitivity.

This module does not refit programs. It consumes biological-donor scores and
association coefficients produced in the frozen programs' original
single-cell scoring framework, then checks the frozen 5% lncRNA-weight trigger,
score concordance, and association direction. A bulk fibrosis-stage model is
not a substitute for the original assay and cannot satisfy this contract.
"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
import pandas as pd


TRIGGER_WEIGHT = 0.05
MIN_CORRELATION = 0.90


class ProgramSensitivityError(ValueError):
    """Raised when a future sensitivity input violates the frozen contract."""


def _require_columns(frame: pd.DataFrame, fields: Iterable[str], label: str) -> None:
    missing = sorted(set(fields) - set(frame.columns))
    if missing:
        raise ProgramSensitivityError(f"{label} missing columns: {missing}")


def validate_program_sensitivity(
    program_content: pd.DataFrame,
    sample_scores: pd.DataFrame,
    stage_effects: pd.DataFrame,
) -> pd.DataFrame:
    """Return one complete sensitivity row for each triggered frozen program."""

    _require_columns(
        program_content,
        (
            "program_uid",
            "lncrna_l1_fraction",
            "requires_leave_all_lncrna_out_sensitivity",
        ),
        "program_content",
    )
    _require_columns(
        sample_scores,
        ("program_uid", "sample_id", "original_score", "without_lncrna_score"),
        "sample_scores",
    )
    _require_columns(
        stage_effects,
        ("program_uid", "original_stage_beta", "without_lncrna_stage_beta"),
        "stage_effects",
    )
    if program_content["program_uid"].duplicated().any():
        raise ProgramSensitivityError("program_content has duplicate program_uid")
    if sample_scores.duplicated(["program_uid", "sample_id"]).any():
        raise ProgramSensitivityError("sample_scores has duplicate program/sample rows")
    if stage_effects["program_uid"].duplicated().any():
        raise ProgramSensitivityError("stage_effects has duplicate program_uid")

    content = program_content.copy()
    content["lncrna_l1_fraction"] = pd.to_numeric(
        content["lncrna_l1_fraction"], errors="raise"
    )
    declared = content["requires_leave_all_lncrna_out_sensitivity"].map(
        {True: True, False: False, "true": True, "false": False}
    )
    if declared.isna().any():
        raise ProgramSensitivityError("invalid sensitivity-trigger flag")
    expected = content["lncrna_l1_fraction"] >= TRIGGER_WEIGHT
    if not np.array_equal(declared.to_numpy(dtype=bool), expected.to_numpy(dtype=bool)):
        raise ProgramSensitivityError("declared trigger disagrees with frozen 5% rule")
    triggered = content.loc[expected, ["program_uid", "lncrna_l1_fraction"]]

    scores = sample_scores.merge(
        triggered[["program_uid"]],
        on="program_uid",
        how="inner",
        validate="many_to_one",
    )
    effects = stage_effects.merge(
        triggered[["program_uid"]], on="program_uid", how="inner", validate="one_to_one"
    )
    if set(scores["program_uid"]) != set(triggered["program_uid"]):
        raise ProgramSensitivityError(
            "sample scores do not cover every triggered program"
        )
    if set(effects["program_uid"]) != set(triggered["program_uid"]):
        raise ProgramSensitivityError(
            "stage effects do not cover every triggered program"
        )

    numeric = ("original_score", "without_lncrna_score")
    for field in numeric:
        scores[field] = pd.to_numeric(scores[field], errors="raise")
        if not np.isfinite(scores[field]).all():
            raise ProgramSensitivityError(f"non-finite {field}")
    for field in ("original_stage_beta", "without_lncrna_stage_beta"):
        effects[field] = pd.to_numeric(effects[field], errors="raise")
        if not np.isfinite(effects[field]).all():
            raise ProgramSensitivityError(f"non-finite {field}")

    rows: list[dict[str, object]] = []
    indexed_effects = effects.set_index("program_uid")
    indexed_content = triggered.set_index("program_uid")
    for uid, group in scores.groupby("program_uid", sort=True):
        if len(group) < 3:
            raise ProgramSensitivityError(
                f"{uid} has fewer than three biological samples"
            )
        pearson = float(
            group["original_score"].corr(
                group["without_lncrna_score"], method="pearson"
            )
        )
        spearman = float(
            group["original_score"].corr(
                group["without_lncrna_score"], method="spearman"
            )
        )
        if not math.isfinite(pearson) or not math.isfinite(spearman):
            raise ProgramSensitivityError(f"{uid} has degenerate score variation")
        original_beta = float(indexed_effects.loc[uid, "original_stage_beta"])
        without_beta = float(indexed_effects.loc[uid, "without_lncrna_stage_beta"])
        direction_preserved = (
            original_beta != 0
            and without_beta != 0
            and math.copysign(1.0, original_beta) == math.copysign(1.0, without_beta)
        )
        passed = (
            pearson >= MIN_CORRELATION
            and spearman >= MIN_CORRELATION
            and direction_preserved
        )
        rows.append(
            {
                "program_uid": uid,
                "lncrna_l1_fraction": float(
                    indexed_content.loc[uid, "lncrna_l1_fraction"]
                ),
                "n_biological_samples": len(group),
                "pearson_correlation": pearson,
                "spearman_correlation": spearman,
                "original_stage_beta": original_beta,
                "without_lncrna_stage_beta": without_beta,
                "stage_direction_preserved": direction_preserved,
                "sensitivity_passed": passed,
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values("program_uid", kind="mergesort")
        .reset_index(drop=True)
    )
