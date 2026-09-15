"""Answer the frozen composition-dependence question, prespecification first.

Stage 3, job 3 of 3.

A composition-dependent gene is not an artifact. Fibrosis genuinely involves
cell-population change, so composition is mediation. T1 asks whether a
composition-INDEPENDENT component exists; its failure is named
``SET_IS_COMPOSITION_MEDIATED`` and is a finding about mechanism, never a
failure of the two-axis result. T2 is the deliverable: a per-gene label.

Both criteria can be satisfied by an absence -- T1 by a set that stops clearing
its null, T2 by a gene that stops clearing BH -- so both are power guarded. A
drop the arm could not have detected is ``indeterminate``, never mediated.

Six lineages is a stricter adjustment than five, which makes a mediated verdict
more likely for a partly technical reason. The set-level MDE is what separates
that from a real mediation finding.

Every statistic comes from the held-back Stage 0b instrument, extended for several
covariates by ``multicovariate_partial`` and proved to reduce to it exactly at
one covariate.
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


class CompositionEvalError(RuntimeError):
    """Raised when the frozen requirements do not hold."""


EXCLUSIVE_CELLS = ("activity_only", "fibrosis_only")
OWN_AXIS = {"activity_only": "activity", "fibrosis_only": "fibrosis"}
ARM_AXES = {"GSE135251": ("activity", "fibrosis"),
            "GSE130970": ("activity", "fibrosis"),
            "GSE193066": ("activity",)}

FROZEN_THRESHOLDS = {
    "t1_a_composition_independent_component_exists": (
        "after adjusting for the six lineage proportions in addition to the "
        "other axis, the assigned set's median absolute partial association on "
        "its own axis still exceeds the 95th percentile of a null that permutes "
        "assignment among expression-matched genes"
    ),
    "t2_per_gene_composition_dependence": (
        "a gene is composition_independent if it retains BH 0.05 significance on "
        "its own axis after the composition adjustment, composition_mediated if "
        "it loses it while the arm's set-level MDE shows the loss was "
        "detectable, and indeterminate otherwise"
    ),
}


def _import(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise CompositionEvalError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def decile_bins(mean_expression: np.ndarray, n_bins: int = 10) -> np.ndarray:
    edges = np.quantile(mean_expression, np.linspace(0, 1, n_bins + 1)[1:-1])
    return np.searchsorted(edges, mean_expression, side="right")


def matched_background(bins, target, pool, generator):
    drawn: list[int] = []
    for bin_id, needed in Counter(bins[target].tolist()).items():
        available = pool[bins[pool] == bin_id]
        if available.size < needed:
            return None
        drawn.extend(generator.choice(available, size=needed, replace=False).tolist())
    return np.asarray(drawn, dtype=int)


def matched_null(statistic, bins, target, pool, *, n_draws, generator):
    draws, skipped, filled = np.empty(n_draws), 0, 0
    while filled < n_draws:
        background = matched_background(bins, target, pool, generator)
        if background is None:
            skipped += 1
            if skipped > n_draws:
                raise CompositionEvalError("the matched background pool is exhausted")
            continue
        draws[filled] = float(np.median(statistic[background]))
        filled += 1
    return {"null": "assignment_permuted_among_expression_matched_genes",
            "n_draws": n_draws,
            "null_median": float(np.median(draws)),
            "null_percentile_95": float(np.percentile(draws, 95)),
            "draws_skipped_for_an_exhausted_bin": skipped}


def set_level_mde(statistic, bins, target, pool, floor, *, generator, n_draws, grid):
    curve, detectable = [], None
    for shift in grid:
        hits = used = 0
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
    return {"minimum_detectable_median_shift": detectable, "power_curve": curve,
            "is_a_gate": False}


def bh_mask(analysis, associations: np.ndarray, residual_df: int) -> np.ndarray:
    critical = analysis.bh_critical_abs_r(int(associations.size), residual_df)
    count = analysis.bh_count_from_abs_r(np.abs(associations), critical)
    if count == 0:
        return np.zeros(associations.size, dtype=bool)
    threshold = np.sort(np.abs(associations))[::-1][count - 1]
    return np.abs(associations) >= threshold


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prespecification", type=Path, required=True)
    parser.add_argument("--sealed-instrument", type=Path, required=True)
    parser.add_argument("--extension", type=Path, required=True)
    parser.add_argument("--stage1-substrate", type=Path, required=True)
    parser.add_argument("--composition-substrate", type=Path, required=True)
    parser.add_argument("--arm-b-source", type=Path, required=True)
    parser.add_argument("--external-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--null-draws", type=int, default=2000)
    parser.add_argument("--mde-draws", type=int, default=400)
    parser.add_argument("--seed", type=int, default=20260827)
    arguments = parser.parse_args()

    if arguments.output.exists():
        raise CompositionEvalError("refusing to overwrite a composition result")
    mc = _import(arguments.extension, "multicovariate")
    analysis = mc.load_sealed(arguments.sealed_instrument)

    verify_frozen_tree(arguments.prespecification)
    prespec = json.loads(
        (arguments.prespecification / "composition_prespecification.json"
         ).read_text(encoding="utf-8"))
    for name, expected in FROZEN_THRESHOLDS.items():
        if prespec["criteria"][name]["threshold"] != expected:
            raise CompositionEvalError(f"{name} threshold drifted since the freeze")
    if prespec["criteria"]["t1_a_composition_independent_component_exists"][
            "failure_is_named"] != "SET_IS_COMPOSITION_MEDIATED":
        raise CompositionEvalError("T1's failure name drifted")
    family = prespec["substrate"]["eligible_family"]
    print(f"verified the frozen prespecification; family {family}", flush=True)

    verify_frozen_tree(arguments.composition_substrate)
    verify_frozen_tree(arguments.stage1_substrate)
    assignment = {
        l.split("\t")[0]: l.split("\t")[1]
        for l in (arguments.stage1_substrate / "gene_assignment.tsv"
                  ).read_text(encoding="utf-8").splitlines()[1:]}
    generator = np.random.default_rng(arguments.seed)

    results = {}
    for cohort, axes_served in ARM_AXES.items():
        comp = np.load(arguments.composition_substrate / "arms"
                       / f"{cohort}_family_proportions.npy")
        axis_rows = [l.split("\t") for l in (
            arguments.composition_substrate / "arms" / f"{cohort}_sample_axis.tsv"
        ).read_text(encoding="utf-8").splitlines()[1:]]

        if cohort == "GSE135251":
            matrix = np.load(arguments.arm_b_source / "molecular" / "rna_values.npy")
            gene_ids = [l.split("\t")[0] for l in (
                arguments.arm_b_source / "molecular" / "rna_feature_axis.tsv"
            ).read_text(encoding="utf-8").splitlines()[1:]]
            rows = analysis.read_table(
                arguments.arm_b_source / "outcomes" / "participant_endpoints.tsv")
            outcomes = {
                "activity": np.asarray([int(r["nas_score"]) for r in rows], float),
                "fibrosis": np.asarray([int(r["fibrosis_stage"]) for r in rows], float)}
            keep = np.arange(matrix.shape[0])
        else:
            matrix = np.load(arguments.stage1_substrate / "arms" / f"{cohort}_cpm.npy")
            gene_ids = [l.split("\t")[0] for l in (
                arguments.stage1_substrate / "arms" / f"{cohort}_gene_axis.tsv"
            ).read_text(encoding="utf-8").splitlines()[1:]]
            samples = (arguments.stage1_substrate / "arms" / f"{cohort}_sample_axis.tsv"
                       ).read_text(encoding="utf-8").splitlines()[1:]
            external = _import(
                Path("scripts/evaluate_two_axis_external_support.py"), f"ext_{cohort}")
            meta, _ = external.read_ragged_table(
                arguments.external_root / cohort / f"{cohort}_metadata.tsv")
            by_run = {r["Run"]: r for r in meta}
            ordered = [by_run[s] for s in samples]
            if cohort == "GSE193066":
                keep = np.asarray(
                    [i for i, r in enumerate(ordered)
                     if r["biopsy"].strip() == "1st biopsy"], dtype=int)
                titles = {ordered[i]["!Sample_title"].strip() for i in keep}
                if keep.size != 106 or len(titles) != 106:
                    raise CompositionEvalError(
                        f"GSE193066 gave {keep.size} rows / {len(titles)} participants, "
                        "not the frozen 106")
                ordered = [ordered[i] for i in keep]
                act, fib = "nafld activity score", "fibrosis stage"
            else:
                keep = np.arange(len(ordered))
                act, fib = "nafld_activity_score", "fibrosis_stage"
            outcomes = {
                "activity": external.integer_column(ordered, act, 0, 8),
                "fibrosis": external.integer_column(ordered, fib, 0, 4)}
            matrix = matrix[keep, :]
        comp = comp[keep, :]
        if comp.shape[0] != matrix.shape[0]:
            raise CompositionEvalError(f"{cohort}: composition and expression disagree")

        usable = matrix.min(axis=0) != matrix.max(axis=0)
        ranks = analysis.average_ranks(matrix[:, usable])
        ranks -= ranks.mean(axis=0, keepdims=True)
        usable_ids = [gene_ids[i] for i in np.flatnonzero(usable)]
        cells = np.asarray([assignment.get(g, "unmapped") for g in usable_ids])
        bins = decile_bins(matrix[:, usable].mean(axis=0))
        pool = np.flatnonzero(cells == "neither")
        lineage_covs = [comp[:, j] for j in range(comp.shape[1])]

        unadjusted, adjusted, conditioning = {}, {}, {}
        for axis in axes_served:
            other = "fibrosis" if axis == "activity" else "activity"
            base, base_report = mc.partial_associations(
                analysis, ranks, outcomes[axis], [outcomes[other]])
            full, full_report = mc.partial_associations(
                analysis, ranks, outcomes[axis], [outcomes[other]] + lineage_covs)
            unadjusted[axis] = (np.abs(base), base_report)
            adjusted[axis] = (np.abs(full), full_report)
            conditioning[axis] = full_report
        print(f"{cohort}: adjusted associations done "
              f"(rank {conditioning[axes_served[0]]['numerical_rank']}/"
              f"{conditioning[axes_served[0]]['n_covariates']}, "
              f"df {conditioning[axes_served[0]]['residual_df']})", flush=True)

        arm = {"axes_served": list(axes_served), "n": int(matrix.shape[0]),
               "genes_entering": int(usable.sum()),
               "conditioning": {a: {k: conditioning[a][k] for k in
                                    ("numerical_rank", "n_covariates",
                                     "condition_number", "residual_df")}
                                for a in axes_served},
               "cells": {}}
        for cell in EXCLUSIVE_CELLS:
            axis = OWN_AXIS[cell]
            target = np.flatnonzero(cells == cell)
            entry = {"joined_genes": int(target.size), "own_axis": axis}
            if axis not in axes_served or target.size == 0 or pool.size == 0:
                entry["applicable"] = False
                entry["why_not_applicable"] = (
                    "this arm does not serve the cell's own axis"
                    if axis not in axes_served else "the cell or pool is empty")
                arm["cells"][cell] = entry
                continue
            entry["applicable"] = True
            statistic = adjusted[axis][0]
            observed = float(np.median(statistic[target]))
            null = matched_null(statistic, bins, target, pool,
                                n_draws=arguments.null_draws, generator=generator)
            before = float(np.median(unadjusted[axis][0][target]))
            drop = before - observed
            mde = set_level_mde(statistic, bins, target, pool,
                                null["null_percentile_95"], generator=generator,
                                n_draws=arguments.mde_draws,
                                grid=np.linspace(0.0, 0.12, 25))
            detectable = mde["minimum_detectable_median_shift"]
            entry.update({
                "median_abs_partial_before_adjustment": before,
                "median_abs_partial_after_adjustment": observed,
                "drop_from_the_composition_adjustment": drop,
                "matched_null": null,
                "exceeds_the_matched_null_p95": bool(
                    observed > null["null_percentile_95"]),
                "set_level_mde": mde,
                "the_drop_was_detectable": bool(
                    detectable is not None and drop >= detectable),
            })
            # ---- T2, per gene ----
            before_mask = bh_mask(analysis, unadjusted[axis][0],
                                  unadjusted[axis][1]["residual_df"])
            after_mask = bh_mask(analysis, statistic,
                                 adjusted[axis][1]["residual_df"])
            was = before_mask[target]
            still = after_mask[target]
            lost = was & ~still
            classes = np.where(
                still, "composition_independent",
                np.where(lost & entry["the_drop_was_detectable"],
                         "composition_mediated", "indeterminate"))
            entry["t2_per_gene"] = {
                "n_genes_in_the_cell": int(target.size),
                "significant_before_adjustment": int(was.sum()),
                "composition_independent": int((classes == "composition_independent").sum()),
                "composition_mediated": int((classes == "composition_mediated").sum()),
                "indeterminate": int((classes == "indeterminate").sum()),
                "fraction_independent_of_those_significant_before": (
                    float(still.sum() / was.sum()) if was.sum() else None),
                "a_mediated_gene_is_not_an_artifact": True,
            }
            arm["cells"][cell] = entry
        results[cohort] = arm
        print(f"{cohort}: cells done", flush=True)

    conditions = []
    for cohort, arm in results.items():
        for cell, entry in arm["cells"].items():
            if not entry.get("applicable"):
                conditions.append({"condition": f"{cohort}:{cell}", "applicable": False,
                                   "met": False,
                                   "why_not_applicable": entry.get("why_not_applicable")})
                continue
            held = entry["exceeds_the_matched_null_p95"]
            detectable = entry["the_drop_was_detectable"]
            conditions.append({
                "condition": f"{cohort}:{cell}",
                "applicable": bool(held or detectable),
                "met": bool(held),
                "verdict": ("composition_independent_component_exists" if held
                            else "set_is_composition_mediated" if detectable
                            else "indeterminate"),
                "why_not_applicable": (
                    None if (held or detectable) else
                    "the set stopped clearing its null but the drop was below the "
                    "arm's set-level MDE, so mediation is indeterminate rather "
                    "than demonstrated"),
                "observed_after": entry["median_abs_partial_after_adjustment"],
                "null_p95": entry["matched_null"]["null_percentile_95"]})
    t1 = analysis.evaluate_gate(
        "t1_a_composition_independent_component_exists", conditions)

    applicable = [c for c in t1["conditions"] if c["applicable"]]
    if t1["verdict"] == "NO_APPLICABLE_CONDITIONS":
        outcome = "INDETERMINATE"
    elif t1["passed"]:
        outcome = "COMPOSITION_INDEPENDENT_COMPONENT_EXISTS"
    elif any(c["met"] for c in applicable):
        outcome = "COMPOSITION_INDEPENDENT_IN_ONE_ARM_ONLY"
    elif all(c.get("verdict") == "set_is_composition_mediated" for c in applicable):
        outcome = "SET_IS_COMPOSITION_MEDIATED"
    else:
        outcome = "INDETERMINATE"
    print(f"outcome: {outcome}", flush=True)

    payload = {
        "schema_version": "masld-bench-composition-result-v1",
        "prespec_id": prespec["prespec_id"],
        "outcome": outcome,
        "the_framing_that_governs_every_reading": prespec[
            "the_framing_that_governs_every_reading"],
        "correction_to_the_frozen_substrate_artifact": prespec[
            "correction_to_the_frozen_substrate_artifact"],
        "honest_provenance": prespec["honest_provenance"],
        "eligible_family": family,
        "arms": results,
        "criteria": {
            "t1_a_composition_independent_component_exists": {
                "threshold": FROZEN_THRESHOLDS[
                    "t1_a_composition_independent_component_exists"],
                "gate": t1,
                "failure_is_named": "SET_IS_COMPOSITION_MEDIATED",
                "six_lineages_is_a_stricter_test_than_five": prespec["criteria"][
                    "t1_a_composition_independent_component_exists"][
                    "six_lineages_is_a_stricter_test_than_five"]},
            "t2_per_gene_composition_dependence": {
                "threshold": FROZEN_THRESHOLDS["t2_per_gene_composition_dependence"],
                "is_the_deliverable": True,
                "per_arm_per_cell": {
                    cohort: {cell: e.get("t2_per_gene")
                             for cell, e in arm["cells"].items() if e.get("applicable")}
                    for cohort, arm in results.items()}},
        },
        "closure": prespec["closure"],
        "controls": {**prespec["controls"], "seed": arguments.seed,
                     "null_draws": arguments.null_draws,
                     "mde_draws": arguments.mde_draws},
        "claim_boundary": prespec["claim_boundary"],
    }
    arguments.output.mkdir(mode=0o750, parents=True)
    (arguments.output / "composition_dependence.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    freeze_tree(arguments.output, {
        "artifact_class": "composition_dependence_result",
        "prespec_id": prespec["prespec_id"], "outcome": outcome,
        "model_fitted": False, "status": "passed"})
    verify_frozen_tree(arguments.output)
    print(json.dumps({"outcome": outcome, "t1_verdict": t1["verdict"],
                      "met_over_applicable":
                          f"{t1['n_conditions_met']}/{t1['n_conditions_applicable']}"},
                     indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
