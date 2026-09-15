"""Record the diagnostics' decision role, and say plainly that it came late.

Stage 0 of the MASLD showcase model, job 3 of 3. This is an additive overlay in
the shape of the W2 interpretive record: it names the frozen output files it
covers by digest, alters no frozen verdict, and registers no check.

The instruction to freeze this addendum *before* the analysis arrived after the
analysis job had already completed. Backdating it would be the exact failure
this campaign has been guarding against all week, so the record states the
ordering it actually had:

* the addendum was recorded AFTER analysis job 21181334 finished
* it was NOT in force for that decision
* it could not have changed that decision, because the declared role is
  asymmetric and may only downgrade a GO, and the decision was STOP
* it is prospective for any future aspect-separability run

The declared role keeps the asymmetry exactly as specified: a reliability
diagnostic may downgrade a GO to ``indeterminate`` and may never upgrade a STOP
to a GO. Conservative in one direction only.

Where this record departs from the proposed threshold, it says so and gives the
evidence. The proposal was that a C3 pair counts as genuinely different gene
orderings only when the cross-aspect ``|rho|`` sits materially below the
geometric mean of the two aspects' own split-half reliabilities. That is the
classical attenuation ceiling, and it assumes the two vectors' sampling errors
are independent. They are not: both come from the same participants and the
same expression matrix. The same assumption already failed loudly once in this
lane, returning disattenuated correlations of 1.50, 1.59 and 1.68. Applied as a
check it would also reject the one contrast this lane agrees is real -- fibrosis
-- which is measured here rather than asserted.

The second half of the record is a post-hoc observation about which pairs
breached which threshold. It was derived by looking at the answer. It is
recorded, marked post-hoc, and authorises nothing.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from masld_bench.hashing import sha256_file


class AddendumError(RuntimeError):
    """Raised when the addendum would assert more than it should."""


def load_analysis(root: Path):
    path = root / "scripts" / "evaluate_gse267145_aspect_separability.py"
    spec = importlib.util.spec_from_file_location("aspect_analysis_for_addendum", path)
    if spec is None or spec.loader is None:
        raise AddendumError("cannot load the analysis module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def contrast_axis_reliability(
    analysis, molecular: Path, endpoints: Path, *, n_splits: int, seed: int
) -> dict[str, dict[str, float]]:
    """Split-half reliability for the contrast axes, on the identical pipeline.

    This is the evidence for departing from the proposed threshold. Fibrosis is
    the axis this lane treats as genuinely different from the three aspects, so
    if the proposed rule would reject fibrosis too, the rule is measuring the
    wrong thing. Computed here rather than assumed.
    """

    rows = analysis.read_table(endpoints)
    values = np.load(molecular / "rna_values.npy")
    ordered = np.sort(values, axis=0)
    realised = ((ordered[1:] != ordered[:-1]).sum(axis=0) + 1) >= 2
    del ordered
    outcomes = {
        axis: np.asarray([int(row[axis]) for row in rows], dtype=float)
        for axis in ("fibrosis", "lobular_necrosis")
    }
    return analysis.split_half_reliability(
        values[:, realised], outcomes, n_splits=n_splits, seed=seed
    )


def build(
    *,
    prespec_digest: str,
    result_digest: str,
    result: dict,
    contrast_reliability: dict,
    analysis_job_id: int,
    analysis_end: str,
) -> tuple[dict, dict]:
    criteria = result["criteria"]
    aspect_reliability = result["diagnostics_not_criteria"]["split_half_reliability"][
        "per_aspect"
    ]
    geometric_ceiling = {}
    for pair in (("steatosis", "ballooning"), ("steatosis", "lobular_inflammation"),
                 ("ballooning", "lobular_inflammation")):
        product = (
            aspect_reliability[pair[0]]["spearman_brown_full_length_reliability"]
            * aspect_reliability[pair[1]]["spearman_brown_full_length_reliability"]
        )
        geometric_ceiling["|".join(pair)] = product**0.5

    fibrosis_full = contrast_reliability["fibrosis"][
        "spearman_brown_full_length_reliability"
    ]
    steatosis_full = aspect_reliability["steatosis"][
        "spearman_brown_full_length_reliability"
    ]
    fibrosis_steatosis_ceiling = (fibrosis_full * steatosis_full) ** 0.5
    fibrosis_steatosis_observed = next(
        entry["spearman_of_association_vectors"]
        for entry in result["diagnostics_not_criteria"][
            "contrast_axes_through_the_same_pipeline"
        ]["pairs"]
        if entry["contrast_axis"] == "fibrosis" and entry["aspect"] == "steatosis"
    )

    addendum = {
        "schema_version": "masld-bench-aspect-separability-addendum-v1",
        "record_is_additive_overlay": True,
        "frozen_verdict_altered": False,
        "promotion_gate_registered": False,
        "frozen_artifacts_this_record_covers": {
            "prespecification_artifacts_sha256": prespec_digest,
            "result_artifacts_sha256": result_digest,
            "decision_in_the_frozen_result": result["decision"],
        },
        "honesty_about_ordering": {
            "frozen_before_the_analysis": False,
            "recorded_after_the_analysis": True,
            "analysis_job_id": analysis_job_id,
            "analysis_end_local": analysis_end,
            "why_it_is_not_backdated": (
                "The instruction to freeze this addendum before the analysis "
                "arrived after the analysis job had already completed. "
                "Recording it as though it came first would be the precise "
                "failure it exists to prevent. Derivable-later is not "
                "recorded-before."
            ),
            "in_force_for_this_decision": False,
            "could_it_have_changed_this_decision": False,
            "why_it_could_not": (
                "The declared role is asymmetric: a reliability diagnostic may "
                "downgrade a GO to indeterminate and may never upgrade a STOP "
                "to a GO. The frozen decision is STOP, so there is no GO for it "
                "to act on. Every frozen criterion failed on its own frozen "
                "threshold."
            ),
            "prospective_from": "any future aspect-separability run in this lane",
        },
        "declared_decision_role": {
            "may_downgrade_a_go_to_indeterminate": True,
            "may_never_upgrade_a_stop_to_a_go": True,
            "rationale": (
                "A noisy association vector correlates poorly with everything, "
                "including a second measurement of itself, so a low "
                "cross-aspect rho is consistent both with different biology and "
                "with being too noisy to tell. A reliability diagnostic "
                "separates those two readings. Nothing rescues a failed frozen "
                "criterion, so the instrument is allowed to act in one "
                "direction only."
            ),
            "threshold_is_a_judgment_call": True,
        },
        "reliability_gate": {
            "proposed_form_recorded_verbatim": (
                "a C3 pair only counts as genuinely different orderings if the "
                "cross-aspect |rho| is materially below the geometric mean of "
                "the two aspects' own split-half reliabilities"
            ),
            "adopted_form": (
                "Per-aspect signal check. A C3 pair counts as genuinely "
                "different gene orderings only if BOTH aspects' own association "
                "vectors are reliably estimated, judged by split-half "
                "reliability sitting clear of zero on this substrate. If either "
                "aspect's vector is not reliably estimated, a low cross-aspect "
                "rho is uninformative rather than passing, and a GO is "
                "downgraded to indeterminate."
            ),
            "why_the_proposed_form_was_not_adopted": (
                "It is the classical attenuation ceiling, which assumes the two "
                "vectors' sampling errors are independent. Here both vectors "
                "come from the same 99 participants and the same expression "
                "matrix, so their errors are shared. That assumption already "
                "failed loudly in this lane: dividing the observed cross-aspect "
                "values by this ceiling returned 1.50, 1.59 and 1.68, which are "
                "impossible for a correlation. A ceiling that the data can "
                "legitimately exceed cannot be used as a gate."
            ),
            "measured_counterexample": {
                "what_it_shows": (
                    "Applied as a gate, the proposed form would also reject the "
                    "one contrast this lane treats as genuinely different."
                ),
                "fibrosis_split_half_full_length_reliability": fibrosis_full,
                "steatosis_split_half_full_length_reliability": steatosis_full,
                "geometric_mean_ceiling_for_fibrosis_vs_steatosis": (
                    fibrosis_steatosis_ceiling
                ),
                "observed_fibrosis_vs_steatosis": fibrosis_steatosis_observed,
                "observed_exceeds_the_proposed_ceiling": (
                    abs(fibrosis_steatosis_observed) > fibrosis_steatosis_ceiling
                ),
                "reading": (
                    "Fibrosis is the axis whose separation makes this a finding "
                    "about the three aspects rather than a resolution limit. A "
                    "rule that classifies fibrosis as uninformative is "
                    "measuring the wrong thing."
                ),
            },
            "geometric_mean_ceilings_for_the_three_aspect_pairs": geometric_ceiling,
            "what_the_correct_pair_level_reference_is": (
                "the measured shared-sample floor already in the frozen result: "
                "permute both outcomes independently through the identical "
                "pipeline, which reproduces the shared error structure by "
                "construction instead of assuming it away"
            ),
            "threshold_is_a_judgment_call": True,
            "not_in_force_for_the_frozen_decision": True,
        },
        "severity_residualized_sensitivity_role": {
            "reported": True,
            "gated_on": False,
            "question_it_answers": (
                "whether any apparent aspect specificity is severity wearing "
                "three labels"
            ),
            "caveat_carried_from_the_frozen_result": (
                "the component sum is the sum of these three aspects, so "
                "residualizing the parts on their own total forces the "
                "residuals to be negatively dependent; negative values there "
                "are arithmetic, not opposing biology"
            ),
        },
        "claim_boundary": (
            "This addendum declares how a diagnostic may act on a future "
            "decision. It alters no frozen verdict, registers no promotion "
            "gate, and authorises no aspect-specific attribution."
        ),
    }

    observations = {
        "schema_version": "masld-bench-aspect-separability-post-hoc-v1",
        "status": "recorded_post_hoc_observation_only",
        "explicitly_post_hoc": True,
        "derived_by_looking_at_the_answer": True,
        "authorises_nothing": True,
        "no_work_was_built_toward_it": True,
        "frozen_verdict_altered": False,
        "observation": {
            "what_was_noticed": (
                "Of the three C1 pairs, only ballooning|lobular_inflammation "
                "breaches the 0.80 threshold, at 0.8181, while both steatosis "
                "pairs sit lower at 0.6800 and 0.7284. Separately, fibrosis "
                "separates from all three aspects at gene level."
            ),
            "c1_pairwise_spearman": {
                "|".join(entry["pair"]): entry["spearman"]
                for entry in criteria["c1_labels_are_not_redundant"]["pairs"]
            },
            "why_it_is_not_a_plan": (
                "This pattern was read off the answer to the frozen question. "
                "Any successor question it suggests needs its own "
                "prespecification, written before anything further is examined, "
                "or it is post hoc dressed as a plan. Nothing here has been "
                "tuned in its direction and no successor analysis has been "
                "built."
            ),
            "who_decides_whether_to_open_it": "the campaign lead, not this record",
        },
        "claim_boundary": (
            "A recorded observation is not a hypothesis under test and not a "
            "result. It may not be cited as evidence for any successor model."
        ),
    }
    return addendum, observations


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--prespecification", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--molecular", type=Path, required=True)
    parser.add_argument("--endpoints", type=Path, required=True)
    parser.add_argument("--analysis-job-id", type=int, required=True)
    parser.add_argument("--analysis-end", type=str, required=True)
    parser.add_argument("--splits", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260828)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    if arguments.output.exists():
        raise AddendumError("refusing to overwrite an addendum")

    verify_frozen_tree(arguments.prespecification)
    verify_frozen_tree(arguments.result)
    result = json.loads(
        (arguments.result / "aspect_separability.json").read_text(encoding="utf-8")
    )
    if result["decision"] != "STOP":
        raise AddendumError(
            "this addendum records a STOP; the frozen decision is not a STOP"
        )

    analysis = load_analysis(arguments.root)
    print("computing contrast-axis reliability for the threshold argument",
          flush=True)
    contrast = contrast_axis_reliability(
        analysis,
        arguments.molecular,
        arguments.endpoints,
        n_splits=arguments.splits,
        seed=arguments.seed,
    )

    addendum, observations = build(
        prespec_digest=sha256_file(arguments.prespecification / "ARTIFACTS.json"),
        result_digest=sha256_file(arguments.result / "ARTIFACTS.json"),
        result=result,
        contrast_reliability=contrast,
        analysis_job_id=arguments.analysis_job_id,
        analysis_end=arguments.analysis_end,
    )
    addendum["contrast_axis_split_half_reliability"] = contrast

    if addendum["honesty_about_ordering"]["frozen_before_the_analysis"] is not False:
        raise AddendumError("the addendum must not claim it came first")
    if addendum["frozen_verdict_altered"] is not False:
        raise AddendumError("the addendum must not alter a frozen verdict")

    arguments.output.mkdir(mode=0o750, parents=True)
    (arguments.output / "diagnostic_role_addendum.json").write_text(
        json.dumps(addendum, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (arguments.output / "post_hoc_observations.json").write_text(
        json.dumps(observations, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(arguments.output, {
        "artifact_class": "gse267145_aspect_separability_diagnostic_addendum",
        "record_is_additive_overlay": True,
        "frozen_verdict_altered": False,
        "promotion_gate_registered": False,
        "in_force_for_the_frozen_decision": False,
        "covers_decision": result["decision"],
        "status": "passed",
    })
    verify_frozen_tree(arguments.output)
    print(json.dumps({
        "covers_decision": result["decision"],
        "in_force_for_this_decision": False,
        "fibrosis_reliability": contrast["fibrosis"][
            "spearman_brown_full_length_reliability"],
        "proposed_ceiling_would_reject_fibrosis": addendum["reliability_gate"][
            "measured_counterexample"]["observed_exceeds_the_proposed_ceiling"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
