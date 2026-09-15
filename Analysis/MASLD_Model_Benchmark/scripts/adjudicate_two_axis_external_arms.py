"""Adjudicate GSE130970 and GSE193066 and write the two-axis activation record.

Records, requirements and label marginals only.  Nothing here reads an expression
value, fits anything, or cross-tabulates the two axes against each other.

Three things are decided here:

1.  A per-cohort verdict for GSE130970 and GSE193066 against the fibrosis label
    standard already registered in
    ``config/evaluation/bulk_existing_source_pool_composition_activation.json``.
    Both were in the family roster but neither was ever assessed, so they are
    unassessed rather than rejected and each needs its own verdict.

2.  A new activation for a two-axis activity-plus-fibrosis task.  The existing
    activation is scoped to ``bulk_lineage_composition_fibrosis_transfer`` and
    its own scope clause says silence means the question was out of scope and
    remains unsettled, so a different task may not inherit it.

3.  An explicit ruling on whether the frozen ``[source_training]`` NASH label
    requirements reaches a graded two-axis endpoint at all.

The detectable-effect figures are computed from the OUTCOME MARGINAL alone.  The
attainable ceiling is a property of the outcome's tie structure and is fixed
before any model exists; reporting it now is the difference between knowing an
arm's limit in advance and discovering it from a null result.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from math import sqrt, tanh
from pathlib import Path
from statistics import NormalDist
from typing import Any, Mapping, Sequence

from scripts.build_showcase_aspect_lineage_substrate import (
    SubstrateError,
    sha256_file,
    write_json,
)

POWER = 0.80
NOMINAL_ALPHA = 0.05


def read_rowname_offset(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read an R ``write.table(row.names=TRUE)`` table, asserting the N+1 shape."""

    lines = [l for l in path.read_text(encoding="utf-8", errors="replace").split("\n") if l != ""]
    if len(lines) < 2:
        raise SubstrateError(f"table has no data rows: {path}")
    header = lines[0].split("\t")
    rows: list[dict[str, str]] = []
    for index, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(header) + 1:
            raise SubstrateError(
                f"{path} line {index} carries {len(fields)} fields, not "
                f"{len(header) + 1}; the R rowname offset does not hold"
            )
        record = dict(zip(header, fields[1:]))
        record["__rowname"] = fields[0]
        rows.append(record)
    return header, rows


def attainable_ceiling(group_sizes: Sequence[int]) -> float:
    """Largest |Spearman rho| any untied continuous predictor can reach.

    Spearman is Pearson on ranks.  A tied ordinal outcome is scored by midranks,
    so its rank vector has less variance than 1..n and no predictor can correlate
    with it perfectly.  The maximum is attained by a predictor that orders the
    outcome exactly, and depends only on the outcome's tie structure.  It is
    therefore knowable before any model is fit.
    """

    if any(g <= 0 for g in group_sizes) or len(group_sizes) < 2:
        raise SubstrateError(f"degenerate outcome marginal: {list(group_sizes)}")
    outcome_ranks: list[float] = []
    start = 1
    for size in group_sizes:
        outcome_ranks += [(start + start + size - 1) / 2.0] * size
        start += size
    n = len(outcome_ranks)
    predictor_ranks = [float(i) for i in range(1, n + 1)]
    mx = sum(predictor_ranks) / n
    my = sum(outcome_ranks) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(predictor_ranks, outcome_ranks))
    vx = sum((a - mx) ** 2 for a in predictor_ranks)
    vy = sum((b - my) ** 2 for b in outcome_ranks)
    if vy == 0.0:
        raise SubstrateError("outcome has zero rank variance")
    return cov / sqrt(vx * vy)


def minimum_detectable_rho(n: int, alpha: float, power: float = POWER) -> float:
    """Smallest |Spearman rho| detectable at the given alpha and power (Fisher z)."""

    if n <= 4:
        raise SubstrateError(f"n={n} is too small for a Fisher-z power statement")
    z_alpha = NormalDist().inv_cdf(1 - alpha / 2)
    z_power = NormalDist().inv_cdf(power)
    return tanh((z_alpha + z_power) / sqrt(n - 3))


def detectable_effect_block(
    *, label: str, group_sizes: Sequence[int], family_sizes: Sequence[int]
) -> dict[str, Any]:
    n = sum(group_sizes)
    ceiling = attainable_ceiling(group_sizes)
    mde = {
        "nominal_alpha_0.05": minimum_detectable_rho(n, NOMINAL_ALPHA),
        **{
            f"bh_single_true_family_{k}": minimum_detectable_rho(n, NOMINAL_ALPHA / k)
            for k in family_sizes
        },
    }
    worst = max(mde.values())
    return {
        "endpoint": label,
        "n": n,
        "marginal": list(group_sizes),
        "levels_realised": len(group_sizes),
        "smallest_level_n": min(group_sizes),
        "attainable_ceiling_rho": ceiling,
        "ceiling_is_a_property_of": "the outcome tie structure alone, fixed before any model exists",
        "minimum_detectable_rho_at_80_percent_power": mde,
        "headroom_at_the_strictest_family": ceiling - worst,
        "arm_can_detect_an_effect_below_its_ceiling": bool(ceiling > worst),
        "ceiling_caveat": (
            "This assumes an untied continuous predictor. Cross-fitted model "
            "predictions carry their own ties, which can only lower the realised "
            "ceiling, never raise it."
        ),
        "family_size_caveat": (
            "The BH family sizes here are illustrative arithmetic, not a "
            "prescription. The confirmatory family is Stage 1's to fix, and the "
            "formula is recorded so it can be recomputed for whatever family is "
            "chosen."
        ),
    }


def adjudicate_gse130970(metadata: Path, counts: Path) -> dict[str, Any]:
    header, rows = read_rowname_offset(metadata)
    fibrosis = Counter(r["fibrosis_stage"] for r in rows)
    activity = Counter(r["nafld_activity_score"] for r in rows)
    mismatches = [
        r["sample_id"]
        for r in rows
        if int(r["steatosis_grade"])
        + int(r["cytological_ballooning_grade"])
        + int(r["lobular_inflammation_grade"])
        != int(r["nafld_activity_score"])
    ]
    control_fields = [
        f for f in header
        if f.lower() in {"disease", "condition", "disease_state", "group", "status",
                         "disease_subtype", "group_in_paper", "diagnosis"}
    ]
    with counts.open(encoding="utf-8") as handle:
        counts_columns = handle.readline().rstrip("\n").split("\t")[1:]
    keys = {r["sample_id"] for r in rows}
    levels = sorted(int(k) for k in fibrosis)
    all_five = levels == [0, 1, 2, 3, 4]
    return {
        "cohort": "GSE130970",
        "prior_state": "in the bulk_existing_source_pool roster, never assessed",
        "fields": {
            "activity": "nafld_activity_score",
            "fibrosis": "fibrosis_stage",
            "components": [
                "steatosis_grade",
                "cytological_ballooning_grade",
                "lobular_inflammation_grade",
            ],
            "covariates_deposited": ["age_at_biopsy", "sex"],
        },
        "unit_census": {
            "metadata_rows": len(rows),
            "distinct_sample_id": len({r["sample_id"] for r in rows}),
            "distinct_biosample": len({r["BioSample"] for r in rows}),
            "distinct_geo_accession": len({r["GEO_Accession (exp)"] for r in rows}),
            "donor_key_exists": False,
            "arithmetic_unit": "bulk_rna_sample",
            "distinct_people_assertable": False,
            "why_not": (
                "No donor key is deposited and no repeat-specimen indicator "
                "exists. There is no evidence of repeat sampling, but absence of "
                "an indicator is not proof of one specimen per person, so the "
                "sample count may not be restated as a participant count."
            ),
        },
        "expression_coverage": {
            "counts_columns": len(counts_columns),
            "metadata_rows": len(rows),
            "columns_without_metadata": len(set(counts_columns) - keys),
            "metadata_without_column": len(keys - set(counts_columns)),
        },
        "fibrosis_standard_assessment": {
            "field": "fibrosis_stage",
            "scale": "0 to 4",
            "levels_present": levels,
            "all_five_levels_distinguishable": all_five,
            "graded_samples": len(rows),
            "distribution": [fibrosis[str(i)] for i in range(5)],
            "control_arm_excluded": 0,
            "control_arm_identifiable": bool(control_fields),
            "verdict": (
                "admissible_for_the_graded_primary"
                if all_five
                else "not_admissible_for_the_graded_primary"
            ),
        },
        "activity_assessment": {
            "field": "nafld_activity_score",
            "scale": "0 to 6 realised of a 0 to 8 possible range",
            "distribution": dict(sorted(activity.items(), key=lambda kv: int(kv[0]))),
            "internal_consistency": (
                "the deposited score equals steatosis + ballooning + lobular "
                "inflammation on every row"
            ),
            "component_sum_mismatches": len(mismatches),
            "saf_activity_A_derivable": True,
            "verdict": "admissible_as_a_graded_activity_axis",
        },
        "both_axes_on_every_row": all(
            r["nafld_activity_score"] != "" and r["fibrosis_stage"] != "" for r in rows
        ),
        "control_field_adjudication": {
            "control_or_disease_field_present": bool(control_fields),
            "fields_found": control_fields,
            "finding": (
                "The metadata carries no disease, condition or control column of "
                "any kind. Every row is a graded liver biopsy."
            ),
            "permits": [
                "All 78 samples enter the graded fibrosis endpoint. The standard's "
                "exclusion rule targets an identifiable healthy-control arm, and "
                "there is none to identify here.",
                "The stage-0 samples are admissible as graded values. Their "
                "fibrosis_stage of 0 is a deposited biopsy grade, not a control "
                "arm mapped onto the low end, which is the substitution the "
                "standard actually prohibits.",
            ],
            "forbids": [
                "any control-versus-disease contrast from this cohort",
                "describing its low end as healthy liver",
                "using it to establish a disease/non-disease boundary",
                "asserting that it contains no control samples",
            ],
            "precise_statement": (
                "No control arm is IDENTIFIABLE in this cohort. That is not the "
                "same as no control arm existing. The field is absent, so the "
                "standard's control-exclusion rule can be neither satisfied nor "
                "violated by evidence, and the honest record is that the question "
                "is unanswerable from the deposit rather than answered."
            ),
        },
    }


def adjudicate_gse193066(metadata: Path, counts: Path) -> dict[str, Any]:
    import re

    header, rows = read_rowname_offset(metadata)
    title = re.compile(r"^HUnafld(\d+)(?:_(\d+))?$")
    unmatched = [r["!Sample_title"] for r in rows if not title.fullmatch(r["!Sample_title"])]
    if unmatched:
        raise SubstrateError(f"unparsable GSE193066 titles: {unmatched[:3]}")
    for r in rows:
        r["__participant"] = title.fullmatch(r["!Sample_title"]).group(1)
    participants = {r["__participant"] for r in rows}
    first = [r for r in rows if r["biopsy"] == "1st biopsy"]
    second = [r for r in rows if r["biopsy"] == "2nd biopsy"]
    first_p = {r["__participant"] for r in first}
    second_p = {r["__participant"] for r in second}
    per = Counter(r["__participant"] for r in rows)
    fib_first = Counter(r["fibrosis stage"] for r in first)
    fib_all = Counter(r["fibrosis stage"] for r in rows)
    act_first = Counter(r["nafld activity score"] for r in first)
    levels_first = sorted(int(k) for k in fib_first)
    with counts.open(encoding="utf-8") as handle:
        counts_columns = handle.readline().rstrip("\n").split("\t")[1:]
    keys = {r["sample_id"] for r in rows}
    suffix_forms = Counter(
        (r["biopsy"], title.fullmatch(r["!Sample_title"]).group(2) or "none")
        for r in rows
    )
    return {
        "cohort": "GSE193066",
        "prior_state": "in the bulk_existing_source_pool roster, never assessed",
        "fields": {
            "activity": "nafld activity score",
            "fibrosis": "fibrosis stage",
            "participant_key": "!Sample_title",
            "timepoint": "biopsy",
            "covariates_deposited": ["age", "Sex_y"],
        },
        "unit_census": {
            "metadata_rows": len(rows),
            "distinct_participants": len(participants),
            "participants_with_two_specimens": sum(1 for v in per.values() if v == 2),
            "participants_with_one_specimen": sum(1 for v in per.values() if v == 1),
            "first_biopsy_rows": len(first),
            "second_biopsy_rows": len(second),
            "second_biopsy_participants_with_a_first": len(second_p & first_p),
            "orphan_second_biopsies": sorted(second_p - first_p),
            "arithmetic_unit": "participant",
            "rows_are_not_units": (
                f"{len(rows)} rows are {len(participants)} participants; "
                f"{len(second)} of the rows are repeat specimens of people already "
                f"present, not independent units"
            ),
            "donor_key_exists": True,
            "donor_key_caveat": (
                "The participant key is derivable from !Sample_title but the "
                "column is heterogeneous: paired participants carry _1 and _2, "
                "unpaired participants carry no suffix at all. Stripping _2 pairs "
                "zero of 58 and silently reports 164 participants. The suffix "
                "forms actually present are recorded here so no reader has to "
                "rediscover this."
            ),
            "suffix_forms_present": {f"{k[0]}|suffix={k[1]}": v for k, v in sorted(suffix_forms.items())},
        },
        "expression_coverage": {
            "counts_columns": len(counts_columns),
            "metadata_rows": len(rows),
            "columns_without_metadata": len(set(counts_columns) - keys),
            "metadata_without_column": len(keys - set(counts_columns)),
        },
        "fibrosis_standard_assessment": {
            "field": "fibrosis stage",
            "cross_sectional_set": "first biopsy only, one row per participant",
            "cross_sectional_n": len(first),
            "levels_present_cross_sectional": levels_first,
            "all_five_levels_distinguishable": levels_first == [0, 1, 2, 3, 4],
            "distribution_cross_sectional": [fib_first.get(str(i), 0) for i in range(5)],
            "distribution_all_rows": [fib_all.get(str(i), 0) for i in range(5)],
            "omitted_level": 4,
            "verdict": "not_admissible_for_the_graded_primary_admissible_as_a_separate_omitted_level_arm",
            "reason": (
                "The cross-sectional set omits stage 4 entirely. The registered "
                "standard states that a scale which pools stages OR OMITS A LEVEL "
                "is a different measurement instrument, reported as a separate "
                "arm and never meta-analysed against a graded scale. Stage 4 "
                "appears in the deposit only among second biopsies, so recovering "
                "it means using repeat specimens of participants already counted, "
                "which the frozen roster rule forbids treating as independent."
            ),
        },
        "activity_assessment": {
            "field": "nafld activity score",
            "scale": "1 to 8 realised, cross-sectional set",
            "distribution": dict(sorted(act_first.items(), key=lambda kv: int(kv[0]))),
            "components_deposited": False,
            "internal_consistency_checkable": False,
            "why_not": "the three NAS components are not deposited, so the sum cannot be re-derived",
            "verdict": "admissible_as_a_graded_activity_axis",
        },
        "both_axes_on_every_row": all(
            r["nafld activity score"] != "" and r["fibrosis stage"] != "" for r in rows
        ),
    }


def source_training_scope_ruling(contract: Path) -> dict[str, Any]:
    text = contract.read_text(encoding="utf-8")
    return {
        "question": (
            "Does the frozen [source_training] NASH label contract reach a graded "
            "two-axis activity-plus-fibrosis endpoint?"
        ),
        "ruling": "out_of_scope",
        "contract_path": str(contract),
        "contract_sha256": sha256_file(contract),
        "reasoning": [
            "[source_training].positive_required defines the admission unit as an "
            "explicit participant-level histology-derived steatohepatitis / NASH / "
            "MASH label requiring steatosis plus ballooning plus lobular "
            "inflammation. It governs a binary positive/negative endpoint.",
            "prohibited_positive_substitutions, which contains fibrosis_stage_only, "
            "is a list of things that may not STAND IN FOR that NASH positive. It "
            "is a rule about substitution, not a rule about which endpoints may "
            "exist.",
            "Neither two-axis endpoint is a NASH label and neither is offered as a "
            "substitute for one. A graded activity sum on a named ordinal scale is "
            "a measured quantity, not a positive-class assignment; a graded "
            "fibrosis stage is a different histological axis entirely.",
            "The contract therefore does not reach these endpoints, and its "
            "prohibitions are not violated by them because they never engage.",
        ],
        "what_still_applies": [
            "admission_unit = participant_baseline_liver_biopsy, as the unit "
            "discipline for any endpoint",
            "required_metadata, as the provenance floor",
            "missing_required: ambiguous or unavailable histology is excluded and "
            "retained in the audit rather than imputed",
        ],
        "what_does_not_apply": [
            "positive_required and negative_required, which have no counterpart in "
            "a graded endpoint",
            "prohibited_positive_substitutions, including fibrosis_stage_only",
            "prohibited_negative_substitutions",
            "eligible_source_roster = [] and transfer_training_allowed_now = false, "
            "which gate NASH transfer training specifically",
        ],
        "precedent": (
            "The 2026-08-26 activation reached the same conclusion for a fibrosis "
            "endpoint: it activated GSE162694 and GSE240729 for a graded fibrosis "
            "task while recording that both fail the NASH label contract."
        ),
        "why_recorded_as_a_ruling": (
            "The registered scope clause states that silence means the question "
            "was out of scope and remains unsettled, and that a later record may "
            "not treat the absence of a clause as permission. This conclusion is "
            "therefore written down with its reasoning rather than assumed from "
            "the contract's silence."
        ),
        "authority": "technical adjudication by the team lead and this producer; NOT user-approved",
        "contract_quotes_verified_present": all(
            token in text
            for token in ("prohibited_positive_substitutions", "fibrosis_stage_only",
                          "positive_required", "admission_unit")
        ),
    }


def build_activation_record(
    *,
    gse130970: Mapping[str, Any],
    gse193066: Mapping[str, Any],
    scope_ruling: Mapping[str, Any],
    prior_activation: Path,
    family_sizes: Sequence[int],
) -> dict[str, Any]:
    prior = json.loads(prior_activation.read_text(encoding="utf-8"))
    standard = prior["fibrosis_label_standard_proposed"]

    graded_arm = detectable_effect_block(
        label="GSE130970 fibrosis_stage 0-4",
        group_sizes=gse130970["fibrosis_standard_assessment"]["distribution"],
        family_sizes=family_sizes,
    )
    graded_activity = detectable_effect_block(
        label="GSE130970 nafld_activity_score",
        group_sizes=[
            v for _, v in sorted(
                gse130970["activity_assessment"]["distribution"].items(),
                key=lambda kv: int(kv[0]),
            )
        ],
        family_sizes=family_sizes,
    )
    omitted_arm = detectable_effect_block(
        label="GSE193066 fibrosis stage 0-3, cross-sectional",
        group_sizes=[
            g for g in gse193066["fibrosis_standard_assessment"][
                "distribution_cross_sectional"] if g > 0
        ],
        family_sizes=family_sizes,
    )
    omitted_activity = detectable_effect_block(
        label="GSE193066 nafld activity score, cross-sectional",
        group_sizes=[
            v for _, v in sorted(
                gse193066["activity_assessment"]["distribution"].items(),
                key=lambda kv: int(kv[0]),
            )
        ],
        family_sizes=family_sizes,
    )

    return {
        "schema_version": "masld-bench-source-roster-admission-record-v1",
        "record_id": "two_axis_activity_fibrosis_external_activation_v1",
        "status": "registered_active",
        "record_is_additive_overlay": True,
        "frozen_bases_edited": False,
        "frozen_bases_read_only": [
            "config/evaluation/cross_cohort_expansion.toml",
            "config/evaluation/microarray_transfer_label_semantics.toml",
            "config/evaluation/bulk_existing_source_pool_composition_activation.json",
        ],
        "cohort_family_id": "bulk_existing_source_pool",
        "activation_requested_for": "two_axis_activity_fibrosis_external_check",
        "activation_is_task_specific": True,
        "why_a_new_record_was_needed": (
            "The 2026-08-26 activation is scoped to "
            "bulk_lineage_composition_fibrosis_transfer, and its own scope clause "
            "states that silence means a question was out of scope and remains "
            "unsettled. A two-axis activity-plus-fibrosis task is a different task "
            "and may not inherit that activation."
        ),
        "inherited_standard": {
            "from": str(prior_activation),
            "sha256": sha256_file(prior_activation),
            "record_id": prior["record_id"],
            "clauses_lifted_verbatim": {
                "positive_required": standard["positive_required"],
                "graded_scale_required": standard["graded_scale_required"],
                "control_arm_handling": standard["control_arm_handling"],
                "pooled_scale_handling": standard["pooled_scale_handling"],
            },
            "lift_is_permitted_by": (
                "the prior record's scope clause, which states a later contract "
                "may lift these clauses verbatim"
            ),
        },
        "arms": {
            "graded_primary": {
                "arm_id": "gse130970_two_axis_graded",
                "cohort": "GSE130970",
                "samples": gse130970["unit_census"]["metadata_rows"],
                "unit": "bulk_rna_sample",
                "donor_key_exists": False,
                "fibrosis_scale": "0 to 4, all five levels",
                "activity_scale": "NAS 0 to 6 realised",
                "verdict": gse130970["fibrosis_standard_assessment"]["verdict"],
            },
            "omitted_level_arm": {
                "arm_id": "gse193066_two_axis_omitted_level",
                "cohort": "GSE193066",
                "rows": gse193066["unit_census"]["metadata_rows"],
                "participants": gse193066["unit_census"]["distinct_participants"],
                "cross_sectional_n": gse193066["fibrosis_standard_assessment"][
                    "cross_sectional_n"],
                "unit": "participant",
                "donor_key_exists": True,
                "fibrosis_scale": "0 to 3, level 4 omitted",
                "activity_scale": "NAS 1 to 8 realised",
                "verdict": gse193066["fibrosis_standard_assessment"]["verdict"],
            },
        },
        "separate_arm_rule": {
            "arms_are_never_pooled": True,
            "arms_are_never_meta_analysed_against_each_other": True,
            "reason": (
                "GSE193066's cross-sectional fibrosis scale omits stage 4. Under "
                "the lifted standard an omitted level makes a different "
                "measurement instrument, and two different instruments may not be "
                "combined into one estimate. Each arm is reported on its own scale "
                "with its own interval."
            ),
            "applies_to_the_fibrosis_axis": True,
            "applies_to_the_activity_axis": True,
            "why_the_activity_axis_too": (
                "The arms are also not pooled on activity. GSE130970 realises NAS "
                "0-6 and GSE193066 realises 1-8, and the two are additionally "
                "denominated in different units, samples against participants. "
                "Pooling one axis while separating the other would produce two "
                "estimates over different denominators in the same record."
            ),
        },
        "detectable_effect": {
            "statistic": "Spearman rho, consistent with the prior activation's rank primary",
            "power": POWER,
            "computed_from": "the outcome marginal alone; no expression value was read",
            "graded_primary_fibrosis": graded_arm,
            "graded_primary_activity": graded_activity,
            "omitted_level_arm_fibrosis": omitted_arm,
            "omitted_level_arm_activity": omitted_activity,
            "interpretation": (
                "Both arms clear their own ceilings comfortably, so the tie "
                "structure of the outcome is not the binding constraint in either "
                "one; n is. The fibrosis ceiling is high in both arms because "
                "neither has a single dominant stage. This is the substantive "
                "contrast with GSE267145, whose 71-of-99 stage-0 mass drives its "
                "ceiling down to roughly 0.79 and leaves far less headroom."
            ),
            "what_this_does_not_establish": (
                "A ceiling above an MDE means an effect of that size COULD be "
                "detected. It does not predict that any effect exists, and it is "
                "not a substitute for the observed in-cohort effect size, which "
                "belongs to Stage 0c and is not reproduced here."
            ),
        },
        "source_training_scope_ruling": scope_ruling,
        "adjudications": {"GSE130970": gse130970, "GSE193066": gse193066},
        "forbidden_under_this_activation": [
            "pooling or meta-analysing the two arms against each other on either axis",
            "treating a sample as a donor, or letting a sample count stand in for a "
            "donor count in a power statement",
            "restating GSE130970's 78 samples as 78 participants; no donor key exists",
            "treating GSE193066's 164 rows as 164 units; they are 106 participants "
            "with 58 repeat specimens",
            "deriving GSE193066's participant key by stripping the _2 suffix alone; "
            "unpaired participants carry no suffix and the strip silently fails",
            "mapping any healthy-control arm to fibrosis stage 0",
            "describing GSE130970's stage-0 samples as controls, or asserting the "
            "cohort contains no controls",
            "recovering GSE193066 stage 4 by adding second biopsies to a "
            "cross-sectional estimate",
            "reading the prior activation's fibrosis-task approval as covering this task",
            "claiming user approval for any per-cohort verdict in this record",
        ],
        "claim_boundary": {
            "external_development_only": True,
            "champion_eligible": False,
            "diagnostic_or_prognostic_claim_allowed": False,
            "clinical_claim_allowed": False,
            "in_cohort_two_axis_result_remains_in_cohort_until_an_arm_reports": True,
        },
        "approval": {
            "scope_decision": {
                "decision": "both cohorts admitted as separate arms",
                "approved_by": "user",
                "relayed_by": "team lead",
                "approval_date": "2026-08-27",
            },
            "technical_adjudication": {
                "content": (
                    "every per-cohort verdict, the fibrosis-standard assessment, "
                    "the source_training scope ruling, the control-field "
                    "adjudication and the detectable-effect figures"
                ),
                "authored_by": "team lead and this producer",
                "user_approved": False,
                "note": (
                    "The user approved the scope only. No per-cohort verdict in "
                    "this record has been seen by the user, and none may be "
                    "presented as user-approved."
                ),
            },
            "prior_activation_approval_does_not_extend": {
                "prior_record": prior["record_id"],
                "prior_approval_date": prior["approval"]["approval_date"],
                "prior_scope": prior["approval"]["note"],
                "extends_to_this_record": False,
            },
        },
        "status_note": "prepared and frozen; Stage 1 is not designed in this record",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bulk-root", type=Path, required=True)
    parser.add_argument("--prior-activation", type=Path, required=True)
    parser.add_argument("--label-contract", type=Path, required=True)
    parser.add_argument("--family-sizes", type=int, nargs="+", default=[4, 12])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    a = adjudicate_gse130970(
        args.bulk_root / "GSE130970" / "GSE130970_metadata.tsv",
        args.bulk_root / "GSE130970" / "GSE130970_counts.tsv",
    )
    b = adjudicate_gse193066(
        args.bulk_root / "GSE193066" / "GSE193066_metadata.tsv",
        args.bulk_root / "GSE193066" / "GSE193066_counts.tsv",
    )
    ruling = source_training_scope_ruling(args.label_contract)
    record = build_activation_record(
        gse130970=a,
        gse193066=b,
        scope_ruling=ruling,
        prior_activation=args.prior_activation,
        family_sizes=args.family_sizes,
    )
    write_json(args.output / "two_axis_activity_fibrosis_external_activation.json", record)

    de = record["detectable_effect"]
    print(json.dumps({
        "GSE130970_verdict": a["fibrosis_standard_assessment"]["verdict"],
        "GSE130970_n": a["unit_census"]["metadata_rows"],
        "GSE193066_verdict": b["fibrosis_standard_assessment"]["verdict"],
        "GSE193066_participants": b["unit_census"]["distinct_participants"],
        "source_training_ruling": ruling["ruling"],
        "fibrosis_ceiling_graded": de["graded_primary_fibrosis"]["attainable_ceiling_rho"],
        "fibrosis_mde_graded_nominal": de["graded_primary_fibrosis"][
            "minimum_detectable_rho_at_80_percent_power"]["nominal_alpha_0.05"],
        "fibrosis_ceiling_omitted": de["omitted_level_arm_fibrosis"]["attainable_ceiling_rho"],
        "fibrosis_mde_omitted_nominal": de["omitted_level_arm_fibrosis"][
            "minimum_detectable_rho_at_80_percent_power"]["nominal_alpha_0.05"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
