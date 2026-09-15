"""Answer the frozen external-support question, prespecification proved first.

Stage 1, job 3 of 3.

**Both external metadata files carry a ragged N+1 header.** GSE130970 has 35
header fields against 36 data fields, GSE193066 has 87 against 88, and in both
every data column sits one position to the right of the name that labels it. A
name-based parser reads GSE193066's `fibrosis stage` column and gets *sex*, and
reads `nafld activity score` and gets *fibrosis*. Nothing downstream would look
wrong. So this script resolves every column by position through a measured
offset, refuses any raggedness other than exactly one, and then validates the
contents of each resolved column against what it is supposed to hold. Measure
column contents, not column names.

The criteria are read out of the frozen prespecification rather than restated.
S1 asks whether an assigned set is elevated on its own axis against an
expression-matched background. S2 asks whether it is not elevated on the other
axis, and because S2 accepts a null it is met only when the arm's set-level
minimum detectable effect shows it could have detected an elevation of the size
seen on the own axis; otherwise it is indeterminate rather than met.

Every statistic is imported from the Stage 0b analysis so training and external
arms share one instrument.
"""

from __future__ import annotations

import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path
from typing import Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


class ExternalError(RuntimeError):
    """Raised when the substrate or the frozen requirements do not hold."""


EXCLUSIVE_CELLS = ("activity_only", "fibrosis_only")
OWN_AXIS = {"activity_only": "activity", "fibrosis_only": "fibrosis"}
OTHER_AXIS = {"activity_only": "fibrosis", "fibrosis_only": "activity"}

FROZEN_THRESHOLDS = {
    "s1_the_assigned_set_is_elevated_on_its_own_axis": (
        "the set's median absolute partial association on its own axis exceeds "
        "the 95th percentile of a null that permutes assignment among "
        "expression-matched genes"
    ),
    "s2_the_assigned_set_is_not_elevated_on_the_other_axis": (
        "the set's median absolute partial association on the other axis does "
        "NOT exceed the 95th percentile of the same matched-permutation null, "
        "AND the arm's set-level minimum detectable effect is at or below the "
        "elevation observed on the set's own axis; if the arm could not have "
        "detected an elevation that size, the result is indeterminate rather "
        "than met"
    ),
}

ARM_AXES = {"GSE130970": ("activity", "fibrosis"), "GSE193066": ("activity",)}


def _import(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ExternalError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_ragged_table(path: Path) -> tuple[list[dict[str, str]], dict[str, object]]:
    """Read a table whose data rows carry one field more than its header.

    The offset is measured, required to be exactly one, and required to be the
    same on every row. Anything else is refused rather than guessed at, because
    a wrong offset produces a table that looks entirely normal and is silently
    wrong in every column.
    """

    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    widths = Counter(len(row) for row in rows)
    if len(widths) != 1:
        raise ExternalError(f"{path.name}: inconsistent row widths {dict(widths)}")
    width = next(iter(widths))
    offset = width - len(header)
    if offset not in (0, 1):
        raise ExternalError(
            f"{path.name}: header has {len(header)} fields against {width} data "
            f"fields, an offset of {offset}. Only 0 or 1 is handled; refusing "
            "to guess which column is which."
        )
    resolved = [
        {name: row[index + offset] for index, name in enumerate(header)}
        for row in rows
    ]
    return resolved, {
        "header_fields": len(header),
        "data_fields": width,
        "measured_offset": offset,
        "is_ragged": offset == 1,
        "note": (
            "every data column sits one position right of the name labelling it; "
            "resolved by measured position, never by name"
            if offset == 1 else "header and data align"
        ),
    }


def integer_column(rows: Sequence[dict[str, str]], name: str, low: int, high: int
                   ) -> np.ndarray:
    """Resolve a column and prove its contents are the ordinal it claims to be."""

    values = [row[name].strip() for row in rows]
    if not all(v.isdigit() for v in values):
        raise ExternalError(
            f"column {name!r} does not hold integers; it holds "
            f"{sorted(set(values))[:6]}. This is the ragged-header failure: the "
            "name does not label the data under it."
        )
    array = np.asarray([int(v) for v in values], dtype=float)
    if array.min() < low or array.max() > high:
        raise ExternalError(
            f"column {name!r} ranges [{array.min()}, {array.max()}], outside "
            f"its declared [{low}, {high}]")
    return array


def decile_bins(mean_expression: np.ndarray, n_bins: int = 10) -> np.ndarray:
    edges = np.quantile(mean_expression, np.linspace(0, 1, n_bins + 1)[1:-1])
    return np.searchsorted(edges, mean_expression, side="right")


def matched_background(bins, target, pool, generator):
    """Draw a background matching the target's bin histogram exactly."""

    drawn: list[int] = []
    for bin_id, needed in Counter(bins[target].tolist()).items():
        available = pool[bins[pool] == bin_id]
        if available.size < needed:
            return None
        drawn.extend(generator.choice(available, size=needed, replace=False).tolist())
    return np.asarray(drawn, dtype=int)


def matched_null(statistic, bins, target, pool, *, n_draws, generator):
    """Median statistic when assignment is permuted among matched genes."""

    draws = np.empty(n_draws, dtype=float)
    skipped = 0
    filled = 0
    while filled < n_draws:
        background = matched_background(bins, target, pool, generator)
        if background is None:
            skipped += 1
            if skipped > n_draws:
                raise ExternalError("the matched background pool is exhausted")
            continue
        draws[filled] = float(np.median(statistic[background]))
        filled += 1
    return {
        "null": "assignment_permuted_among_expression_matched_genes",
        "n_draws": n_draws,
        "null_median": float(np.median(draws)),
        "null_percentile_95": float(np.percentile(draws, 95)),
        "null_percentile_99": float(np.percentile(draws, 99)),
        "draws_skipped_for_an_exhausted_bin": skipped,
    }


def set_level_mde(statistic, bins, target, pool, floor, *, generator, n_draws, grid):
    """Smallest injected median shift detected at 80% power against the floor.

    Derived for the median-shift statistic actually used. A single-correlation
    MDE is a different estimand with different power and is reported only as
    context for the per-gene claim this stage does not make.
    """

    curve = []
    detectable = None
    for shift in grid:
        hits = 0
        used = 0
        for _ in range(n_draws):
            background = matched_background(bins, target, pool, generator)
            if background is None:
                continue
            used += 1
            hits += float(np.median(statistic[background] + shift)) > floor
        power = hits / used if used else 0.0
        curve.append({"injected_median_shift": float(shift), "power": power})
        if detectable is None and power >= 0.80:
            detectable = float(shift)
    return {
        "definition": (
            "the smallest median shift in absolute partial association, "
            "injected into a random expression-matched set of the same size, "
            "that clears the matched-permutation 95th percentile at 80% power"),
        "minimum_detectable_median_shift": detectable,
        "power_curve": curve,
        "is_a_gate": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prespecification", type=Path, required=True)
    parser.add_argument("--analysis-script", type=Path, required=True)
    parser.add_argument("--substrate", type=Path, required=True)
    parser.add_argument("--external-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--null-draws", type=int, default=2000)
    parser.add_argument("--mde-draws", type=int, default=400)
    parser.add_argument("--seed", type=int, default=20260827)
    arguments = parser.parse_args()

    if arguments.output.exists():
        raise ExternalError("refusing to overwrite an external-support result")
    analysis = _import(arguments.analysis_script, "arm_a_analysis")

    verify_frozen_tree(arguments.prespecification)
    prespec = json.loads(
        (arguments.prespecification / "two_axis_external_prespecification.json"
         ).read_text(encoding="utf-8"))
    for name, expected in FROZEN_THRESHOLDS.items():
        if prespec["criteria"][name]["threshold"] != expected:
            raise ExternalError(f"{name} threshold drifted since the freeze")
    if not prespec["the_deliverable_and_its_honesty_constraint"][
            "external_support_is_set_level_not_per_gene"]:
        raise ExternalError("the frozen set-level constraint is missing")
    print("verified the frozen prespecification", flush=True)

    verify_frozen_tree(arguments.substrate)
    substrate = json.loads(
        (arguments.substrate / "substrate.json").read_text(encoding="utf-8"))
    assignment = {
        line.split("\t")[0]: line.split("\t")[1]
        for line in (arguments.substrate / "gene_assignment.tsv"
                     ).read_text(encoding="utf-8").splitlines()[1:]
    }
    generator = np.random.default_rng(arguments.seed)

    results, header_audit = {}, {}
    for cohort, axes_served in ARM_AXES.items():
        meta, audit = read_ragged_table(
            arguments.external_root / cohort / f"{cohort}_metadata.tsv")
        header_audit[cohort] = audit
        stable_ids = [
            line.split("\t")[0]
            for line in (arguments.substrate / "arms" / f"{cohort}_gene_axis.tsv"
                         ).read_text(encoding="utf-8").splitlines()[1:]
        ]
        matrix = np.load(arguments.substrate / "arms" / f"{cohort}_cpm.npy")
        samples = (arguments.substrate / "arms" / f"{cohort}_sample_axis.tsv"
                   ).read_text(encoding="utf-8").splitlines()[1:]
        if matrix.shape != (len(samples), len(stable_ids)):
            raise ExternalError(f"{cohort}: matrix does not match its axes")
        by_run = {row["Run"]: row for row in meta}
        ordered = [by_run[s] for s in samples]

        if cohort == "GSE130970":
            activity = integer_column(ordered, "nafld_activity_score", 0, 8)
            fibrosis = integer_column(ordered, "fibrosis_stage", 0, 4)
            components = sum(
                integer_column(ordered, name, 0, 3)
                for name in ("steatosis_grade", "cytological_ballooning_grade",
                             "lobular_inflammation_grade"))
            if not np.array_equal(components, activity):
                raise ExternalError(
                    "GSE130970: the deposited NAS is not steatosis + ballooning "
                    "+ lobular inflammation on every row")
            unit = {"unit": "bulk_rna_sample", "n": len(ordered),
                    "distinct_people_assertable": False,
                    "nas_identity_verified_on_every_row": True}
        else:
            biopsy = [row["biopsy"].strip() for row in ordered]
            titles = [row["!Sample_title"].strip() for row in ordered]
            keep = np.asarray([i for i, b in enumerate(biopsy) if b == "1st biopsy"],
                              dtype=int)
            chosen = {titles[i] for i in keep}
            if keep.size != 106 or len(chosen) != 106:
                raise ExternalError(
                    f"GSE193066: first-biopsy selection gave {keep.size} rows and "
                    f"{len(chosen)} distinct participants; the frozen contract "
                    "requires exactly 106 of each")
            ordered = [ordered[i] for i in keep]
            matrix = matrix[keep, :]
            activity = integer_column(ordered, "nafld activity score", 0, 8)
            fibrosis = integer_column(ordered, "fibrosis stage", 0, 4)
            if int(fibrosis.max()) == 4:
                raise ExternalError(
                    "GSE193066 cross-sectional fibrosis realises stage 4; the "
                    "omitted-level ruling that scopes this arm depends on it not")
            unit = {"unit": "participant", "n": 106,
                    "rows_in_the_deposited_matrix": len(samples),
                    "distinct_participants_asserted": len(chosen),
                    "note": prespec["unit_of_inference"]["GSE193066"][
                        "REQUIRED_NOTE_ON_THE_SUBSTRATE_ARTIFACT"]}

        library = matrix.sum(axis=1)
        outcomes = {"activity": activity, "fibrosis": fibrosis}
        usable = matrix.min(axis=0) != matrix.max(axis=0)
        ranks = analysis.average_ranks(matrix[:, usable])
        ranks -= ranks.mean(axis=0, keepdims=True)
        usable_ids = [stable_ids[i] for i in np.flatnonzero(usable)]
        partial = {}
        for exposure, covariate in (("activity", "fibrosis"), ("fibrosis", "activity")):
            unit_covariate = analysis.unit_ranks(outcomes[covariate])
            scaled, norms = analysis.unit_columns(
                analysis.residualize(ranks, unit_covariate))
            residual = analysis.residualize(
                analysis.centred_ranks(outcomes[exposure]), unit_covariate)
            vector = scaled.T @ (residual / np.sqrt((residual**2).sum()))
            vector[norms == 0.0] = 0.0
            partial[exposure] = np.abs(vector)
        print(f"{cohort}: partials on {int(usable.sum())} genes", flush=True)

        cells = np.asarray([assignment.get(g, "unmapped") for g in usable_ids])
        bins = decile_bins(matrix[:, usable].mean(axis=0))
        pool = np.flatnonzero(cells == "neither")
        arm = {
            "unit_of_inference": unit,
            "axes_this_arm_serves": list(axes_served),
            "genes_entering": int(usable.sum()),
            "library_size_ratio": float(library.max() / library.min()),
            "library_size_orthogonality_to_each_outcome": {
                name: analysis.spearman_correlation(library.tolist(), v.tolist())
                for name, v in outcomes.items()},
            "background_pool_size": int(pool.size),
            "cells": {},
        }
        for cell in EXCLUSIVE_CELLS:
            target = np.flatnonzero(cells == cell)
            entry = {"joined_genes": int(target.size),
                     "assigned_genes": substrate["assignment"]["cell_sizes"][cell],
                     "applicable": bool(target.size and pool.size)}
            if entry["applicable"]:
                for role, axis in (("own", OWN_AXIS[cell]), ("other", OTHER_AXIS[cell])):
                    if axis not in axes_served:
                        entry[role] = {"axis": axis, "served_by_this_arm": False}
                        continue
                    statistic = partial[axis]
                    observed = float(np.median(statistic[target]))
                    null = matched_null(statistic, bins, target, pool,
                                        n_draws=arguments.null_draws,
                                        generator=generator)
                    block = {
                        "axis": axis, "served_by_this_arm": True,
                        "observed_median_abs_partial": observed,
                        "matched_null": null,
                        "exceeds_the_matched_null_p95": bool(
                            observed > null["null_percentile_95"]),
                        "elevation_over_the_null_median": observed - null["null_median"],
                    }
                    if role == "other":
                        elevation = entry.get("own", {}).get(
                            "elevation_over_the_null_median")
                        mde = set_level_mde(
                            statistic, bins, target, pool, null["null_percentile_95"],
                            generator=generator, n_draws=arguments.mde_draws,
                            grid=np.linspace(0.0, 0.12, 25))
                        detectable = mde["minimum_detectable_median_shift"]
                        block["set_level_mde"] = mde
                        block["own_axis_elevation_for_comparison"] = elevation
                        block["arm_could_have_detected_an_elevation_that_size"] = bool(
                            elevation is not None and detectable is not None
                            and detectable <= elevation)
                    entry[role] = block
            arm["cells"][cell] = entry
        results[cohort] = arm
        print(f"{cohort}: cells done", flush=True)

    def conditions(role: str) -> list[dict[str, object]]:
        out = []
        for cohort, arm in results.items():
            for cell, entry in arm["cells"].items():
                block = entry.get(role, {})
                if not entry.get("applicable") or not block.get("served_by_this_arm"):
                    out.append({"condition": f"{cohort}:{cell}", "applicable": False,
                                "met": False,
                                "why_not_applicable": "this arm does not serve that axis"})
                    continue
                if role == "own":
                    out.append({"condition": f"{cohort}:{cell}", "applicable": True,
                                "met": bool(block["exceeds_the_matched_null_p95"]),
                                "observed": block["observed_median_abs_partial"],
                                "null_p95": block["matched_null"]["null_percentile_95"]})
                else:
                    powered = block["arm_could_have_detected_an_elevation_that_size"]
                    quiet = not block["exceeds_the_matched_null_p95"]
                    out.append({
                        "condition": f"{cohort}:{cell}", "applicable": bool(powered),
                        "met": bool(quiet and powered),
                        "evidence_state": ("tested_negative" if (quiet and powered)
                                           else "indeterminate" if quiet else "elevated"),
                        "why_not_applicable": (
                            None if powered else
                            "the set-level MDE exceeds the own-axis elevation, so an "
                            "absence here is indeterminate rather than tested_negative"),
                        "observed": block["observed_median_abs_partial"],
                        "null_p95": block["matched_null"]["null_percentile_95"],
                        "set_level_mde": block["set_level_mde"][
                            "minimum_detectable_median_shift"]})
        return out

    s1 = analysis.evaluate_gate("s1_the_assigned_set_is_elevated_on_its_own_axis",
                                conditions("own"))
    s2 = analysis.evaluate_gate(
        "s2_the_assigned_set_is_not_elevated_on_the_other_axis", conditions("other"))

    if "NO_APPLICABLE_CONDITIONS" in (s1["verdict"], s2["verdict"]):
        outcome = "INDETERMINATE"
    elif not s1["passed"]:
        outcome = "NOT_EXTERNALLY_SUPPORTED"
    elif not s2["passed"]:
        outcome = "TRANSFERS_WITHOUT_SPECIFICITY"
    else:
        second_arm = any(e["condition"].startswith("GSE193066") and e["applicable"]
                         and e["met"] for e in s1["conditions"])
        outcome = "EXTERNALLY_SUPPORTED" if second_arm else "SUPPORTED_SINGLE_ARM"
    print(f"outcome: {outcome}", flush=True)

    payload = {
        "schema_version": "masld-bench-two-axis-external-result-v1",
        "prespec_id": prespec["prespec_id"],
        "outcome": outcome,
        "external_support_is_set_level_not_per_gene": True,
        "honest_provenance": prespec["honest_provenance"],
        "ragged_header_audit": {
            "what_was_found": (
                "both external metadata files carry one more data field than "
                "header field, so every data column sits one position right of "
                "the name labelling it. Read by name, GSE193066's 'fibrosis "
                "stage' returns sex and 'nafld activity score' returns fibrosis, "
                "and nothing downstream would look wrong."),
            "how_it_is_handled": (
                "every column is resolved by measured position, the offset is "
                "required to be exactly one and identical on every row, and each "
                "resolved column's contents are validated against the ordinal it "
                "claims to be. Measure column contents, not names."),
            "per_cohort": header_audit,
            "discovered_after_the_freeze_and_changes_no_criterion": True,
        },
        "arms": results,
        "criteria": {
            "s1_the_assigned_set_is_elevated_on_its_own_axis": {
                "threshold": FROZEN_THRESHOLDS[
                    "s1_the_assigned_set_is_elevated_on_its_own_axis"], "gate": s1},
            "s2_the_assigned_set_is_not_elevated_on_the_other_axis": {
                "threshold": FROZEN_THRESHOLDS[
                    "s2_the_assigned_set_is_not_elevated_on_the_other_axis"], "gate": s2},
        },
        "multiplicity_determination": {
            "recorded_at_analysis_time_not_pre_registered": True,
            "what_the_family_would_be": (
                "six gate tests: S1 and S2 on activity_only and fibrosis_only in "
                "GSE130970, and on activity_only in GSE193066, each read at a "
                "95th percentile"),
            "why_no_correction_is_applied": (
                "the decision rule is conjunctive, not disjunctive. Every "
                "positive outcome requires several tests to pass together, and a "
                "conjunction makes a false positive rarer rather than commoner: "
                "under the null its probability is bounded by that of its least "
                "likely member."),
            "a_correction_to_that_reasoning": (
                "The familywise control comes almost entirely from S1, not from "
                "the conjunction as a whole. S2 is met when a test FAILS to clear "
                "its p95, so under the null S2 is met with probability about 0.95 "
                "and contributes essentially nothing to false-positive control. "
                "Its role is to prevent a different error, claiming a specificity "
                "that is not there, and that error is governed by power. The "
                "set-level MDE guard rather than a multiplicity correction is the "
                "right instrument for it; a correction applied to S2 would make a "
                "spurious S2 pass more likely, not less."),
            "is_a_gate": False,
        },
        "claim_boundary": prespec["claim_boundary"],
        "controls": {
            "seed": arguments.seed, "null_draws": arguments.null_draws,
            "mde_draws": arguments.mde_draws,
            "statistics_imported_from_stage_0b": True, "no_model_was_fitted": True},
    }
    arguments.output.mkdir(mode=0o750, parents=True)
    (arguments.output / "two_axis_external_support.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    freeze_tree(arguments.output, {
        "artifact_class": "two_axis_external_support_result",
        "prespec_id": prespec["prespec_id"], "outcome": outcome,
        "external_support_is_set_level": True, "model_fitted": False,
        "status": "passed"})
    verify_frozen_tree(arguments.output)
    print(json.dumps({
        "outcome": outcome,
        "s1": {"verdict": s1["verdict"],
               "met_over_applicable":
                   f"{s1['n_conditions_met']}/{s1['n_conditions_applicable']}"},
        "s2": {"verdict": s2["verdict"],
               "met_over_applicable":
                   f"{s2['n_conditions_met']}/{s2['n_conditions_applicable']}"},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
