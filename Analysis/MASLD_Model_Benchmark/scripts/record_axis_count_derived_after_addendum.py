"""Two measurements the Stage 0b check did not make, recorded as derived-after.

Stage 0b, job 3 of 3. This is an additive overlay. It alters no frozen verdict,
registers no check, and could not have changed the outcome: the outcome was
already computed and held-back when these were commissioned.

Both were requested after the Stage 0b analysis job had completed. Recording
them as though they came first would be the precise failure the two-job pattern
exists to prevent, so the ordering is written down as it happened.

**Split-half reliability per axis.**  Stage 0's addendum established that a
less reliably estimated axis correlates lower with everything, which is a
measurement artefact rather than biology. The same confound reaches a *count*:
a noisier axis yields fewer BH-surviving genes for a measurement reason. Stage 0
measured this for the three NAS components and for fibrosis and necrosis, but
not for either composite, so a count comparison across the Stage 0b axes had no
reliability to sit beside it. This supplies it.

It also tests a prediction made before the measurement: summing two positively
correlated 0-2 grades should give a 0-4 grade with higher reliability than
either component, so A = ballooning + lobular_inflammation should exceed both.
The prediction is recorded here verbatim and the measurement is reported
whichever way it falls, because a failed prediction about the instrument is
more informative than a confirmed one.

**The A-vs-F partial at n=99.**  Stage 0b's D1 tested (S+A) against F. It never
tested A against F on its own. That contrast is computed here through the
identical code path -- the functions are imported from the analysis script
rather than reimplemented -- so it is comparable to the four families in the
frozen result without being one of them.

Neither number is a check. Neither may be cited as pre-registered.
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


class AddendumError(RuntimeError):
    """Raised when the overlay cannot be built against the frozen result."""


def _load_analysis(path: Path):
    """Import the frozen analysis module so the code path is provably identical."""

    spec = importlib.util.spec_from_file_location("axis_count_analysis", path)
    if spec is None or spec.loader is None:
        raise AddendumError("cannot import the analysis module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


#: Recorded before the split-half numbers were computed.
RELIABILITY_PREDICTION = (
    "Summing two positively correlated 0-2 grades should give a 0-4 grade with "
    "higher reliability than either component, so A = ballooning + "
    "lobular_inflammation should exceed both ballooning (0.479) and lobular "
    "inflammation (0.518) as Stage 0 measured them. Recorded as a prediction "
    "about the instrument, not as a result, and reported whichever way it falls."
)


def split_half_reliability(
    analysis,
    raw_genes: np.ndarray,
    outcomes: dict[str, np.ndarray],
    *,
    n_splits: int,
    seed: int,
) -> dict[str, dict[str, float]]:
    """How well an axis's own gene ordering reproduces on independent halves.

    Ranks are recomputed inside each half, because a Spearman estimated on 49 or
    50 participants is the thing whose reproducibility is in question. A gene
    that varies across all 99 can still be constant inside a half; those are
    dropped for that split only, on the intersection of the two halves, and the
    surviving count is reported rather than assumed to be the whole universe.
    This is the procedure Stage 0 used, so the numbers are comparable to its.
    """

    generator = np.random.default_rng(seed)
    n = raw_genes.shape[0]
    collected: dict[str, list[float]] = {axis: [] for axis in outcomes}
    surviving: list[int] = []
    for _ in range(n_splits):
        order = generator.permutation(n)
        halves = (order[: n // 2], order[n // 2 :])
        keep = np.ones(raw_genes.shape[1], dtype=bool)
        blocks = []
        for indices in halves:
            values = raw_genes[indices, :]
            keep &= values.min(axis=0) != values.max(axis=0)
            blocks.append((indices, values))
        if not keep.any():
            continue
        surviving.append(int(keep.sum()))
        standardized = []
        for indices, values in blocks:
            ranked = analysis.average_ranks(values[:, keep])
            ranked -= ranked.mean(axis=0, keepdims=True)
            scaled, norms = analysis.unit_columns(ranked)
            if np.any(norms == 0.0):
                raise AddendumError("a constant column survived the split filter")
            standardized.append((indices, scaled))
        for axis, outcome in outcomes.items():
            vectors = []
            usable = True
            for indices, block in standardized:
                values = outcome[indices]
                if values.min() == values.max():
                    usable = False
                    break
                vectors.append(block.T @ analysis.unit_ranks(values))
            if not usable:
                continue
            collected[axis].append(
                _spearman_of_vectors(analysis, vectors[0], vectors[1])
            )
    summary: dict[str, dict[str, float]] = {}
    median_surviving = float(np.median(surviving)) if surviving else float("nan")
    for axis, values in collected.items():
        array = np.asarray(values, dtype=float)
        if array.size == 0:
            summary[axis] = {"n_usable_splits": 0}
            continue
        half = float(np.median(array))
        summary[axis] = {
            "n_usable_splits": int(array.size),
            "median_genes_non_constant_in_both_halves": median_surviving,
            "median_half_length_reliability": half,
            "percentile_5": float(np.percentile(array, 5)),
            "percentile_95": float(np.percentile(array, 95)),
            # Spearman-Brown from half length back to the full cohort.
            "spearman_brown_full_length_reliability": (
                2.0 * half / (1.0 + half) if half > -1.0 else float("nan")
            ),
        }
    return summary


def _spearman_of_vectors(analysis, left: np.ndarray, right: np.ndarray) -> float:
    from scipy.stats import rankdata

    left_ranked = rankdata(left, method="average")
    right_ranked = rankdata(right, method="average")
    left_centred = left_ranked - left_ranked.mean()
    right_centred = right_ranked - right_ranked.mean()
    denominator = float(
        np.sqrt((left_centred**2).sum() * (right_centred**2).sum())
    )
    if denominator == 0.0:
        raise AddendumError("a constant association vector reached the comparison")
    return float(left_centred @ right_centred / denominator)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-script", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--endpoints", type=Path, required=True)
    parser.add_argument("--molecular", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--count-null-draws", type=int, default=20000)
    parser.add_argument("--null-block", type=int, default=250)
    parser.add_argument("--splits", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260827)
    arguments = parser.parse_args()

    if arguments.output.exists():
        raise AddendumError("refusing to overwrite an addendum")

    analysis = _load_analysis(arguments.analysis_script)
    verify_frozen_tree(arguments.result)
    frozen = json.loads(
        (arguments.result / "axis_count.json").read_text(encoding="utf-8")
    )
    if frozen["prespec_id"] != "gse267145_axis_count_v1":
        raise AddendumError("this is not the Stage 0b result")
    outcome = frozen["outcome"]

    rows = analysis.read_table(arguments.endpoints)
    participants = [row["participant_id"] for row in rows]
    n = len(rows)
    components = {
        name: np.asarray([int(row[name]) for row in rows], dtype=float)
        for name in analysis.COMPONENTS
    }
    fibrosis = np.asarray([int(row["fibrosis"]) for row in rows], dtype=float)
    saf_activity = sum(components[name] for name in analysis.SAF_ACTIVITY)
    nas_activity = sum(components[name] for name in analysis.NAS_ACTIVITY)
    axes = {
        **components,
        "fibrosis": fibrosis,
        "saf_activity_sum": saf_activity,
        "nas_activity_sum": nas_activity,
    }
    sexes = [row["recorded_sex"] for row in rows]

    axis_rows = analysis.read_table(arguments.molecular / "participant_axis.tsv")
    if [row["participant_id"] for row in axis_rows] != participants:
        raise AddendumError("the molecular participant axis is out of order")
    features = analysis.read_table(arguments.molecular / "rna_feature_axis.tsv")
    values = np.load(arguments.molecular / "rna_values.npy")
    if values.shape != (n, len(features)):
        raise AddendumError("the expression matrix does not match its axes")

    ordered = np.sort(values, axis=0)
    realised = ((ordered[1:] != ordered[:-1]).sum(axis=0) + 1) >= 2
    del ordered
    universe = np.flatnonzero(realised)
    realised_values = values[:, universe]
    if int(universe.size) != frozen["criteria"][
        "d1_activity_and_fibrosis_are_independent"
    ]["directions"][0]["realised_universe"]["non_constant_across_all_participants"]:
        raise AddendumError("the realised universe differs from the frozen result")
    print(f"realised gene universe: {universe.size}", flush=True)

    gene_ranks = analysis.average_ranks(realised_values)
    gene_ranks -= gene_ranks.mean(axis=0, keepdims=True)

    # ---- the A-vs-F partial, identical code path, never a check ----
    families = {}
    for offset, (exposure, covariate) in enumerate(
        (("saf_activity_sum", "fibrosis"), ("fibrosis", "saf_activity_sum"))
    ):
        unit_covariate = analysis.unit_ranks(axes[covariate])
        residual_genes = analysis.residualize(gene_ranks, unit_covariate)
        scaled, norms = analysis.unit_columns(residual_genes)
        keep = norms > 0.0
        result = analysis.partial_family(
            scaled[:, keep],
            analysis.residualize(
                analysis.centred_ranks(axes[exposure]), unit_covariate
            ),
            unit_covariate,
            exposure=exposure,
            covariate=covariate,
            n_participants=n,
            n_draws=arguments.count_null_draws,
            block=arguments.null_block,
            seed=arguments.seed + 70 + offset,
            strata=sexes,
            null_name="within_recorded_sex_residual_permutation",
        )
        result["is_a_gate"] = False
        result["dropped_as_collinear_with_the_covariate"] = int((~keep).sum())
        families[f"{exposure}|{covariate}"] = result
        print(f"derived-after family done: {exposure}|{covariate}", flush=True)

    # ---- split-half reliability per axis ----
    reliability = split_half_reliability(
        analysis,
        realised_values,
        axes,
        n_splits=arguments.splits,
        seed=arguments.seed + 1,
    )
    print("split-half reliability done", flush=True)

    stage0 = {
        "ballooning": 0.4790336956496051,
        "lobular_inflammation": 0.5184070284950335,
        "steatosis": 0.5529578702179281,
        "fibrosis": 0.35742181194650596,
    }
    measured_a = reliability["saf_activity_sum"][
        "spearman_brown_full_length_reliability"
    ]
    prediction_held = bool(
        measured_a > stage0["ballooning"] and measured_a > stage0["lobular_inflammation"]
    )

    payload = {
        "schema_version": "masld-bench-axis-count-addendum-v1",
        "prespec_id": frozen["prespec_id"],
        "record_is_additive_overlay": True,
        "frozen_verdict_altered": False,
        "in_force_for_the_frozen_decision": False,
        "is_a_gate": False,
        "frozen_artifacts_this_record_covers": {
            "outcome_in_the_frozen_result": outcome,
            "result_artifacts_sha256": json.loads(
                (arguments.result / "ARTIFACTS.json").read_text(encoding="utf-8")
            )["artifacts"][0]["sha256"],
        },
        "honesty_about_ordering": {
            "frozen_before_the_analysis": False,
            "recorded_after_the_analysis": True,
            "could_it_have_changed_the_outcome": False,
            "why_it_could_not": (
                "The outcome was computed and sealed before either measurement "
                "here was commissioned. Neither is a gate and neither enters "
                "the decision rule, so there is nothing for them to move."
            ),
            "why_it_is_not_backdated": (
                "The request for both arrived after the analysis job had "
                "started. Recording them as pre-registered would be the precise "
                "failure this pattern exists to prevent. Derivable-later is not "
                "recorded-before."
            ),
            "may_not_be_cited_as_pre_registered": True,
        },
        "split_half_reliability": {
            "what_it_is": (
                "how well one axis's own gene ordering reproduces on two "
                "independent halves of the cohort, Spearman-Brown corrected "
                "back to full length. The same procedure Stage 0 used, so the "
                "numbers are comparable to its."
            ),
            "why_it_belongs_beside_a_count": (
                "a less reliably estimated axis yields fewer BH-surviving genes "
                "for a measurement reason rather than a biological one, so a "
                "count comparison across axes is not interpretable without it"
            ),
            "per_axis": reliability,
            "stage0_reference_values": stage0,
            "prediction_recorded_before_the_measurement": RELIABILITY_PREDICTION,
            "prediction_held": prediction_held,
            "measured_saf_activity_full_length_reliability": measured_a,
        },
        "saf_activity_versus_fibrosis_partial": {
            "why_it_is_here": (
                "Stage 0b's D1 tested the full NAS sum against fibrosis. It "
                "never tested activity{ballooning + lobular_inflammation} "
                "against fibrosis on its own. This supplies that contrast "
                "through the identical imported code path."
            ),
            "is_a_gate": False,
            "derived_after_the_frozen_outcome": True,
            "families": families,
        },
        "claim_boundary": (
            "This record supplies two measurements the frozen gate did not "
            "make. It authorises no aspect-specific attribution, changes no "
            "outcome, and neither number may be presented as pre-registered."
        ),
    }

    arguments.output.mkdir(mode=0o750, parents=True)
    (arguments.output / "derived_after_addendum.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(arguments.output, {
        "artifact_class": "gse267145_axis_count_derived_after_addendum",
        "prespec_id": frozen["prespec_id"],
        "frozen_verdict_altered": False,
        "is_a_gate": False,
        "status": "passed",
    })
    verify_frozen_tree(arguments.output)

    print(json.dumps({
        "outcome_unchanged": outcome,
        "reliability_prediction_held": prediction_held,
        "split_half_full_length": {
            axis: round(entry.get("spearman_brown_full_length_reliability", float("nan")), 4)
            for axis, entry in reliability.items()
        },
        "a_vs_f_counts": {
            key: entry["genes_bh_below_0_05"] for key, entry in families.items()
        },
        "a_vs_f_floors": {
            key: entry["count_null"]["null_percentile_95_count"]
            for key, entry in families.items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
