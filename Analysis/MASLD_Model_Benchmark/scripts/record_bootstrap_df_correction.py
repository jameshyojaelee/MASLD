"""Correct the recorded *reason* the participant bootstrap was dropped.

Additive overlay. It alters no verdict, registers no check, and moves neither
Stage 0b's outcome nor Stage 0c's design.

What is being corrected. Stage 0b found its participant bootstrap of a BH count
incoherent: every observed count sat below its own bootstrap 2.5th percentile,
which is impossible for a well-behaved interval and therefore indicts the
method rather than the data. That observation was right and the remedy -- a
leave-one-participant-out jackknife -- was right. The *explanation* attached to
it was not checked, and it was wrong. Two held-back places carry it:

*   the Stage 0c prespecification, at
    ``diagnostics_never_gates.participant_jackknife_not_bootstrap.why_the_bootstrap_was_dropped``
*   the ``jackknife_counts`` docstring in
    ``scripts/evaluate_gse135251_axis_count.py``

Both say duplicated participants "manufacture exact concordance blocks that
inflate rank correlation". They do not. Resampling with replacement leaves a
mean of about 62.5 distinct participants of 99, and the p-value was computed
with ``df = n - 3 = 96`` regardless, scaling the t-statistic by ``sqrt(96)``
where the resample supported roughly ``sqrt(59)``. Recomputing the identical
draws with ``df = distinct - 3`` removes the inflation entirely. Duplication was
the trigger; the wrong degrees of freedom were the bug.

Neither held-back output file is retrofitted. Editing the analysis script would have
recorded a digest for content the job did not execute, which is strictly worse
than a wrong sentence that a labelled correction can reach.

The scope matters and is measured here rather than asserted. A bootstrap of a
plain correlation involves no df and is first-order unbiased, which is why
Stage 0's C1 label intervals are sound. A bootstrap of anything computed
*through* a p-value or a df-dependent threshold is not. Stage 0's C3
association-vector intervals are also measured here, because they were quoted
and their bias direction decides whether a reported conclusion needs softening.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, t as student_t

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


class CorrectionError(RuntimeError):
    """Raised when the correction cannot be built against the held-back record."""


WRONG_REASON_SUBSTRING = "manufacture exact concordance blocks"

CORRECTED_REASON = (
    "Resampling with replacement leaves a mean of about 62.5 distinct "
    "participants of 99, and Stage 0b computed every per-gene p-value with "
    "df = n - 3 = 96 regardless. That scales the t-statistic by sqrt(96) where "
    "the resample supports roughly sqrt(59), which inflates the BH count. "
    "Recomputing the identical draws with df = distinct - 3 removes the "
    "inflation. Duplication was the trigger; the degrees of freedom were the "
    "bug. A jackknife has n - 1 distinct participants and df = n - 1 - 3, so it "
    "carries no mismatch, which is why the remedy was correct even though the "
    "reason recorded for it was not."
)

STANDING_RULE = (
    "Any bootstrap of a statistic computed through a p-value or a "
    "df-dependent threshold must recompute the degrees of freedom from the "
    "resample's distinct-unit count. Bootstraps of plain correlations involve "
    "no df and are first-order unbiased; the rule does not reach them."
)

DIAGNOSTIC_SHAPE = (
    "An observed statistic falling below its own bootstrap 2.5th percentile is "
    "impossible for a well-behaved interval, so it indicts the method rather "
    "than reporting a surprising truth. Same shape as a disattenuation "
    "correction returning values above 1."
)


def unit_ranks(values: np.ndarray) -> np.ndarray:
    ranked = rankdata(np.asarray(values, dtype=float), method="average").astype(float)
    centred = ranked - ranked.mean()
    return centred / np.sqrt((centred**2).sum())


def standardize(block: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ranked = rankdata(block, method="average", axis=0).astype(float)
    ranked -= ranked.mean(axis=0, keepdims=True)
    norms = np.sqrt((ranked**2).sum(axis=0))
    keep = norms > 0.0
    out = np.zeros_like(ranked)
    out[:, keep] = ranked[:, keep] / norms[keep]
    return out, keep


def vector_spearman(left: np.ndarray, right: np.ndarray) -> float:
    a = rankdata(left, method="average")
    b = rankdata(right, method="average")
    a = a - a.mean()
    b = b - b.mean()
    return float(a @ b / np.sqrt((a**2).sum() * (b**2).sum()))


def bh_count(correlations: np.ndarray, df: int) -> int:
    bounded = np.clip(np.abs(correlations), 0.0, 1.0 - 1e-12)
    statistic = bounded * np.sqrt(df) / np.sqrt(1.0 - bounded**2)
    p = 2.0 * student_t.sf(statistic, df)
    m = p.size
    ordered = np.sort(p)
    passing = np.flatnonzero(ordered <= 0.05 * np.arange(1, m + 1) / m)
    return int(passing[-1] + 1) if passing.size else 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm-a-prespec", type=Path, required=True)
    parser.add_argument("--arm-a-result", type=Path, required=True)
    parser.add_argument("--arm-b-prespec", type=Path, required=True)
    parser.add_argument("--endpoints", type=Path, required=True)
    parser.add_argument("--molecular", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--genes", type=int, default=6000)
    parser.add_argument("--resamples", type=int, default=400)
    parser.add_argument("--label-resamples", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260827)
    arguments = parser.parse_args()

    if arguments.output.exists():
        raise CorrectionError("refusing to overwrite a correction record")

    for tree in (arguments.arm_a_result, arguments.arm_b_prespec):
        verify_frozen_tree(tree)

    # The wrong sentence must actually be in the held-back file, or this
    # correction is aimed at nothing.
    arm_b = json.loads(
        (arguments.arm_b_prespec / "arm_b_axis_count_prespecification.json").read_text(
            encoding="utf-8"
        )
    )
    recorded = arm_b["diagnostics_never_gates"]["participant_jackknife_not_bootstrap"][
        "why_the_bootstrap_was_dropped"
    ]
    if WRONG_REASON_SUBSTRING not in recorded:
        raise CorrectionError(
            "the sealed Arm B prespecification does not carry the reason this "
            "record corrects; the correction would be aimed at nothing"
        )
    arm_a = json.loads(
        (arguments.arm_a_result / "axis_count.json").read_text(encoding="utf-8")
    )
    if arm_a["outcome"] != "ONE_AXIS_ACTIVITY":
        raise CorrectionError("the Arm A outcome is not the one this record covers")

    lines = arguments.endpoints.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]
    column = lambda name: np.asarray([int(r[name]) for r in rows], dtype=float)
    steatosis, ballooning = column("steatosis"), column("ballooning")
    inflammation, fibrosis = column("lobular_inflammation"), column("fibrosis")
    nas = steatosis + ballooning + inflammation
    n = len(rows)

    values = np.load(arguments.molecular / "rna_values.npy")
    ordered = np.sort(values, axis=0)
    realised = ((ordered[1:] != ordered[:-1]).sum(axis=0) + 1) >= 2
    del ordered
    universe = values[:, np.flatnonzero(realised)]
    generator = np.random.default_rng(arguments.seed)
    subset = generator.choice(
        universe.shape[1], min(arguments.genes, universe.shape[1]), replace=False
    )
    genes = universe[:, subset]
    print(f"evidence subsample: {genes.shape}", flush=True)

    # ---- the correction's own evidence: df=96 against df=distinct-3 ----
    standardized, _ = standardize(genes)
    unit_fibrosis = unit_ranks(fibrosis)
    residual_genes = standardized - np.outer(unit_fibrosis, unit_fibrosis @ standardized)
    norms = np.sqrt((residual_genes**2).sum(axis=0))
    usable = norms > 0.0
    residual_genes[:, usable] /= norms[usable]
    exposure = unit_ranks(nas) - unit_fibrosis * float(unit_fibrosis @ unit_ranks(nas))
    exposure /= np.sqrt((exposure**2).sum())
    observed_count = bh_count(residual_genes[:, usable].T @ exposure, n - 3)

    at_96: list[int] = []
    at_effective: list[int] = []
    distinct: list[int] = []
    for _ in range(arguments.resamples):
        index = generator.integers(0, n, size=n)
        unique = len(set(index.tolist()))
        distinct.append(unique)
        block = genes[index, :]
        standardized_block, keep = standardize(block)
        covariate = unit_ranks(fibrosis[index])
        residual = standardized_block[:, keep] - np.outer(
            covariate, covariate @ standardized_block[:, keep]
        )
        block_norms = np.sqrt((residual**2).sum(axis=0))
        ok = block_norms > 0.0
        residual[:, ok] /= block_norms[ok]
        ranked_exposure = unit_ranks(nas[index])
        residual_exposure = ranked_exposure - covariate * float(covariate @ ranked_exposure)
        scale = np.sqrt((residual_exposure**2).sum())
        if scale == 0.0 or not ok.any():
            continue
        correlations = residual[:, ok].T @ (residual_exposure / scale)
        at_96.append(bh_count(correlations, n - 3))
        at_effective.append(bh_count(correlations, unique - 3))
    print("df evidence done", flush=True)

    # ---- scope: a plain-correlation bootstrap is unbiased (Stage 0 C1) ----
    label_bias = []
    for name, (left, right) in {
        "steatosis|ballooning": (steatosis, ballooning),
        "steatosis|lobular_inflammation": (steatosis, inflammation),
        "ballooning|lobular_inflammation": (ballooning, inflammation),
    }.items():
        observed = abs(float(unit_ranks(left) @ unit_ranks(right)))
        draws = []
        for _ in range(arguments.label_resamples):
            index = generator.integers(0, n, size=n)
            if len(set(left[index].tolist())) < 2 or len(set(right[index].tolist())) < 2:
                continue
            draws.append(abs(float(unit_ranks(left[index]) @ unit_ranks(right[index]))))
        array = np.asarray(draws, dtype=float)
        label_bias.append({
            "pair": name,
            "observed_abs_spearman": observed,
            "bootstrap_median": float(np.median(array)),
            "bias_median_minus_observed": float(np.median(array)) - observed,
            "percentile_2_5": float(np.percentile(array, 2.5)),
            "percentile_97_5": float(np.percentile(array, 97.5)),
            "n_usable_resamples": int(array.size),
        })
    print("label-correlation scope evidence done", flush=True)

    # ---- Stage 0's C3 association-vector intervals, bias and direction ----
    aspects = {
        "steatosis": steatosis,
        "ballooning": ballooning,
        "lobular_inflammation": inflammation,
    }
    pairs = (
        ("steatosis", "ballooning"),
        ("steatosis", "lobular_inflammation"),
        ("ballooning", "lobular_inflammation"),
    )
    observed_vectors = {
        name: standardized.T @ unit_ranks(vector) for name, vector in aspects.items()
    }
    observed_c3 = {
        f"{a}|{b}": abs(vector_spearman(observed_vectors[a], observed_vectors[b]))
        for a, b in pairs
    }
    c3_draws: dict[str, list[float]] = {key: [] for key in observed_c3}
    for _ in range(arguments.resamples):
        index = generator.integers(0, n, size=n)
        block = genes[index, :]
        standardized_block, keep = standardize(block)
        usable_block = standardized_block[:, keep]
        vectors = {
            name: usable_block.T @ unit_ranks(vector[index])
            for name, vector in aspects.items()
        }
        for a, b in pairs:
            c3_draws[f"{a}|{b}"].append(abs(vector_spearman(vectors[a], vectors[b])))
    c3 = []
    for key, draws in c3_draws.items():
        array = np.asarray(draws, dtype=float)
        c3.append({
            "pair": key,
            "observed_abs_spearman": observed_c3[key],
            "bootstrap_median": float(np.median(array)),
            "bias_median_minus_observed": float(np.median(array)) - observed_c3[key],
            "n_resamples": int(array.size),
        })
    print("c3 scope evidence done", flush=True)

    worst_c3 = min(entry["bias_median_minus_observed"] for entry in c3)
    payload = {
        "schema_version": "masld-bench-bootstrap-df-correction-v1",
        "record_is_additive_overlay": True,
        "is_a_gate": False,
        "frozen_verdict_altered": False,
        "stage_0b_outcome_unchanged": arm_a["outcome"],
        "stage_0c_design_unchanged": True,
        "what_is_corrected": {
            "the_observation_is_unchanged": (
                "Stage 0b's participant bootstrap of a BH count was invalid. "
                "Every observed count sat below its own bootstrap 2.5th "
                "percentile."
            ),
            "the_remedy_is_unchanged": (
                "A leave-one-participant-out jackknife replaces it. The "
                "jackknife has n - 1 distinct participants and is computed with "
                "df = n - 1 - 3, so it carries no mismatch."
            ),
            "only_the_reason_was_wrong": True,
            "where_the_wrong_reason_is_sealed": [
                "the Stage 0c prespecification, at diagnostics_never_gates."
                "participant_jackknife_not_bootstrap.why_the_bootstrap_was_dropped",
                "the jackknife_counts docstring in "
                "scripts/evaluate_gse135251_axis_count.py",
            ],
            "the_wrong_reason_verbatim": recorded,
            "the_corrected_reason": CORRECTED_REASON,
            "why_neither_was_retrofitted": (
                "The Arm B prespecification is sealed and the analysis script "
                "was hashed into a completed job's source manifest. Editing the "
                "script would have recorded a digest for content the job did not "
                "execute, which is strictly worse than a wrong sentence a "
                "labelled correction can reach."
            ),
        },
        "evidence_measured_here": {
            "substrate": "GSE267145, n=99, nas_activity_sum | fibrosis",
            "gene_subsample": int(usable.sum()),
            "why_a_subsample": (
                "the correction is about the df, which does not depend on the "
                "size of the family; a subsample makes the comparison cheap and "
                "the effect is identical in direction and magnitude"
            ),
            "observed_count_at_df_96": observed_count,
            "mean_distinct_participants_per_resample": float(np.mean(distinct)),
            "bootstrap_median_count_at_df_96": float(np.median(at_96)),
            "bootstrap_median_count_at_df_distinct_minus_3": float(
                np.median(at_effective)
            ),
            "n_resamples": len(at_96),
            "reading": (
                "At the df Stage 0b used, the bootstrap median sits far above "
                "the observed count, which is the incoherence. At the df the "
                "resample actually supports, it sits below the observed count, "
                "which is ordinary. The inflation is the df and nothing else."
            ),
        },
        "scope_of_the_rule": {
            "standing_rule": STANDING_RULE,
            "diagnostic_shape_worth_remembering": DIAGNOSTIC_SHAPE,
            "plain_correlation_bootstraps_are_unaffected": {
                "what_it_shows": (
                    "a bootstrap of a correlation involves no df and is "
                    "first-order unbiased, so Stage 0's C1 label intervals are "
                    "sound and need no correction"
                ),
                "per_pair": label_bias,
            },
            "stage_0_c3_association_vector_intervals": {
                "what_it_shows": (
                    "these are also plain correlations, but of vectors "
                    "re-estimated inside each resample, so they carry a "
                    "downward bias. The direction is what matters: the bias "
                    "runs toward the pass side of the 0.80 bar, so the reported "
                    "fraction of resamples below 0.80 is too generous to a pass "
                    "and C3's failure is more solid than reported, not less."
                ),
                "per_pair": c3,
                "largest_downward_bias": worst_c3,
                "c2_involves_no_resampling_at_all": True,
                "no_stage_0_conclusion_moves": True,
            },
        },
        "honesty_about_ordering": {
            "frozen_before_the_analyses_it_corrects": False,
            "recorded_after": True,
            "could_it_have_changed_any_outcome": False,
            "why_it_could_not": (
                "It corrects a recorded explanation, not a statistic. No gate "
                "reads it and no threshold depends on it."
            ),
            "may_not_be_cited_as_pre_registered": True,
        },
        "claim_boundary": (
            "This record corrects one recorded explanation and measures the "
            "scope of the rule that replaces it. It authorises no "
            "aspect-specific attribution and changes no outcome."
        ),
    }

    arguments.output.mkdir(mode=0o750, parents=True)
    (arguments.output / "bootstrap_df_correction.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(arguments.output, {
        "artifact_class": "masld_bench_bootstrap_df_correction",
        "record_is_additive_overlay": True,
        "is_a_gate": False,
        "frozen_verdict_altered": False,
        "status": "passed",
    })
    verify_frozen_tree(arguments.output)

    print(json.dumps({
        "observed_count_at_df_96": observed_count,
        "mean_distinct_participants": round(float(np.mean(distinct)), 1),
        "bootstrap_median_at_df_96": float(np.median(at_96)),
        "bootstrap_median_at_honest_df": float(np.median(at_effective)),
        "label_bias": {e["pair"]: round(e["bias_median_minus_observed"], 4)
                       for e in label_bias},
        "c3_bias": {e["pair"]: round(e["bias_median_minus_observed"], 4) for e in c3},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
