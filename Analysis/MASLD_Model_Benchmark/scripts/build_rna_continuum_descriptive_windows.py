"""Descriptive sliding-window reporting layer over the five-cohort RNA continuum.

This builds the descriptive projection permitted by OVERALL_PLAN.md:61-63 -- "Molecular
continuum positions and sliding-window summaries are descriptive projections; they may
not be interpreted as longitudinal trajectories or optimized as clinical labels."  There
is no predictive arm and no selection endpoint.  Windows are cohort-level aggregates
over nine fixed bins covering roughly 20% of a cohort each, carry ``inference_use =
False``, and are structurally incapable of serving as an endpoint.

Two arms, filed separately and never pooled:

*Primary* covers the two ``independent_primary`` cohorts (GSE162694, GSE213621; 501
samples).  *Quarantined* covers the three ``source_overlap`` cohorts (GSE126848,
GSE130970, GSE135251; 343 samples), which overlap the external source the 139-gene
fixed-projection signature was derived on.  Dropping the 343 silently would be scaling
the request down without asking; pooling them into the primary arm would launder the
leakage boundary the frozen requirements draws; emitting both with the boundary named on
every row does neither.

The unit of observation is the sequencing sample.  Four of the five cohorts carry no
donor key and none can be derived from what is on disk, so nothing here is a
donor-level quantity and repeated-biopsy structure is untestable rather than absent.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

# Frozen window geometry, transcribed from
# scripts/analysis/histology_anchored_continuum/00_prespecification.json "windows".
WINDOW_CENTERS: tuple[float, ...] = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
WINDOW_WIDTH = 0.2
# "lower_closed_upper_open_except_final_upper_closed"
INTERVAL_RULE = "lower_closed_upper_open_except_final_upper_closed"

# The only axis that transfers across cohorts.  signature_pc1 is refit per cohort and
# therefore carries no shared coordinate; it is excluded from the window layer rather
# than reported with a caveat.
AXIS_ID = "fixed_projection"

PRIMARY_COHORTS: tuple[str, ...] = ("GSE162694", "GSE213621")
QUARANTINED_COHORTS: tuple[str, ...] = ("GSE126848", "GSE130970", "GSE135251")

ARM_PRIMARY = "primary_independent_validation"
ARM_QUARANTINED = "quarantined_source_overlap"

AXIS_LEAKAGE = {
    ARM_PRIMARY: "none_independent_of_signature_discovery",
    ARM_QUARANTINED: "source_overlap_with_signature_discovery",
}

# Cohorts whose source tables carry no fibrosis and no NAS column at all.  This is
# structural absence, not attrition: the values were never recorded, so a window here
# has axis position against an anchor that does not exist.  "Never recorded" and
# "missing" read very differently, and only one of them invites someone to go looking
# for data that is not there.
HISTOLOGY_ANCHOR_STATUS = {
    "GSE126848": "never_recorded_in_source",
    "GSE130970": "recorded",
    "GSE135251": "recorded",
    "GSE162694": "recorded",
    "GSE213621": "recorded_fibrosis_only_nas_never_recorded",
}

# Ratio of samples to donors in GSE193066, the only cohort in the wider substrate with a
# real donor key: 164 samples over 106 donors, 58 of them with two biopsies.  GSE193066
# sits OUTSIDE these 844 and is never folded in.  It is used here only to size a
# what-if inflation of the naive sample-level standard error, reported as a sensitivity
# column beside the naive value rather than in place of it.
DESIGN_EFFECT_SAMPLES = 164
DESIGN_EFFECT_DONORS = 106
DESIGN_EFFECT_SOURCE = "GSE193066_164_samples_106_donors_out_of_substrate"

NA_TOKENS = frozenset({"", "NA", "N/A", "None", "NaN", "nan", "null"})


class WindowError(RuntimeError):
    """The descriptive window layer could not be built as specified."""


def design_effect_se_multiplier() -> float:
    return math.sqrt(DESIGN_EFFECT_SAMPLES / DESIGN_EFFECT_DONORS)


def unquote(value: str) -> str:
    text = value.strip()
    if len(text) >= 2 and text[0] == '"' and text[-1] == '"':
        return text[1:-1]
    return text


def parse_float(value: str | None) -> float | None:
    if value is None or unquote(value) in NA_TOKENS:
        return None
    try:
        number = float(unquote(value))
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def parse_int(value: str | None) -> int | None:
    number = parse_float(value)
    if number is None or number != int(number):
        return None
    return int(number)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path) -> Iterator[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield from handle
    else:
        with path.open("r", encoding="utf-8") as handle:
            yield from handle


def read_delimited(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    reader = csv.DictReader(open_text(path), delimiter=delimiter)
    rows = []
    for row in reader:
        rows.append({unquote(k): unquote(v) if v is not None else ""
                     for k, v in row.items() if k is not None})
    return rows


def window_bounds(index: int) -> tuple[float, float, bool]:
    """Return (lower, upper, is_final) for the 1-based window index."""

    center = WINDOW_CENTERS[index - 1]
    return center - WINDOW_WIDTH / 2, center + WINDOW_WIDTH / 2, index == len(WINDOW_CENTERS)


def in_window(percentile: float, index: int) -> bool:
    lower, upper, final = window_bounds(index)
    if percentile < lower:
        return False
    return percentile <= upper if final else percentile < upper


@dataclass(frozen=True)
class WindowStat:
    n_samples: int
    mean: float | None
    se: float | None


def summarize_window(scores: Sequence[float]) -> WindowStat:
    """Mean and naive sample-level standard error over the samples inside one window."""

    finite = [value for value in scores if value is not None and math.isfinite(value)]
    if not finite:
        return WindowStat(n_samples=0, mean=None, se=None)
    mean = statistics.fmean(finite)
    se = statistics.stdev(finite) / math.sqrt(len(finite)) if len(finite) > 1 else None
    return WindowStat(n_samples=len(finite), mean=mean, se=se)


# --- Inputs ------------------------------------------------------------------------


@dataclass(frozen=True)
class Substrate:
    axis: dict[tuple[str, str], float]          # (dataset, sample) -> axis percentile
    stage_complete: set[tuple[str, str]]        # samples with a recorded fibrosis stage
    scale_native: dict[str, str]
    stage_comparable: dict[str, bool]
    sex_source: dict[str, str]
    cohort_n: dict[str, int]


def load_substrate(axis_rows: Sequence[Mapping[str, str]],
                   clinical_rows: Sequence[Mapping[str, str]]) -> Substrate:
    axis: dict[tuple[str, str], float] = {}
    for row in axis_rows:
        if row.get("axis_id") != AXIS_ID:
            continue
        percentile = parse_float(row.get("axis_percentile"))
        if percentile is None:
            continue
        axis[(row["dataset"], row["sample_id"])] = percentile
    if not axis:
        raise WindowError(f"no {AXIS_ID} positions were loaded")

    stage_complete: set[tuple[str, str]] = set()
    scale_native: dict[str, str] = {}
    stage_comparable: dict[str, bool] = {}
    sex_source: dict[str, str] = {}
    cohort_n: dict[str, int] = defaultdict(int)
    for row in clinical_rows:
        key = (row["dataset"], row["sample_id"])
        cohort_n[row["dataset"]] += 1
        if parse_int(row.get("fibrosis_rederived_value")) is not None:
            stage_complete.add(key)
        scale_native[row["dataset"]] = row["fibrosis_scale_native"]
        stage_comparable[row["dataset"]] = row["cross_cohort_stage_comparable"] == "TRUE"
        sex_source[row["dataset"]] = row["sex_source"]
    return Substrate(axis=axis, stage_complete=stage_complete, scale_native=scale_native,
                     stage_comparable=stage_comparable, sex_source=sex_source,
                     cohort_n=dict(cohort_n))


@dataclass(frozen=True)
class FeatureScores:
    """One scored molecular layer, keyed by feature and sample."""

    layer: str
    aggregation: str
    display: dict[str, str]
    scores: dict[str, dict[tuple[str, str], float]]


def load_hotspot_scores(rows: Sequence[Mapping[str, str]]) -> list[FeatureScores]:
    """Load the 117 frozen Hotspot program scores.

    A feature is registered from the rows that mention it, NOT from the rows whose score
    parses.  `hotspot_fibroblasts_8584ba834d31e608` (Activated stellate, PDGFRA) has all
    844 rows present with `outcome_z = NA`, because the signature-exclusion testability
    check leaves it unscored at 72.0% retained L1 against a 0.8 minimum.  Registering
    features by parsed score would drop it from the layer entirely, turning an
    explainable zero-occupancy row into an invisible absence -- which is strictly worse,
    because a reader comparing 116 against the source's 117 has nothing to find.
    """

    display: dict[str, str] = {}
    scores: dict[str, dict[tuple[str, str], float]] = {}
    for row in rows:
        uid = row["program_uid"]
        scores.setdefault(uid, {})
        display.setdefault(uid, row.get("module_name", "") or uid)
        value = parse_float(row.get("outcome_z"))
        if value is None:
            continue
        scores[uid][(row["dataset"], row["sample_id"])] = value
    return [FeatureScores(layer="hotspot_program", aggregation="within_cohort_z",
                          display=display, scores=scores)]


def load_system_scores(rows: Sequence[Mapping[str, str]]) -> list[FeatureScores]:
    """Load the 43 molecular systems under each aggregation.

    Features are registered from the rows that mention them, for the same reason as the
    Hotspot layer: an all-NA system must appear at zero occupancy, not disappear.
    """

    aggregations: set[str] = set()
    universe: set[str] = set()
    display: dict[str, str] = {}
    observed: dict[str, dict[str, dict[tuple[str, str], float]]] = defaultdict(
        lambda: defaultdict(dict))
    for row in rows:
        uid = row["community_id"]
        aggregation = row["aggregation_method"]
        aggregations.add(aggregation)
        universe.add(uid)
        display.setdefault(uid, row.get("system_display", "") or uid)
        value = parse_float(row.get("system_score_z"))
        if value is None:
            continue
        observed[aggregation][uid][(row["dataset"], row["sample_id"])] = value
    return [
        FeatureScores(
            layer="molecular_system",
            aggregation=aggregation,
            display=display,
            scores={uid: dict(observed[aggregation].get(uid, {})) for uid in sorted(universe)},
        )
        for aggregation in sorted(aggregations)
    ]


def expected_row_count(layers: Sequence[FeatureScores], n_cohorts: int) -> int:
    """The table shape is known before the table is built: features x windows x cohorts."""

    return sum(len(layer.scores) for layer in layers) * len(WINDOW_CENTERS) * n_cohorts


def feature_shape(layers: Sequence[FeatureScores]) -> list[dict[str, Any]]:
    return [
        {
            "layer": layer.layer,
            "aggregation_method": layer.aggregation,
            "features_expected": len(layer.scores),
            "features_emitted": len(layer.scores),
            "features_with_no_finite_score": sum(
                1 for feature in layer.scores.values() if not feature),
        }
        for layer in layers
    ]


def load_testability(rows: Sequence[Mapping[str, str]]) -> dict[tuple[str, str], dict[str, Any]]:
    table: dict[tuple[str, str], dict[str, Any]] = {}
    for row in rows:
        table[(row["program_uid"], row["dataset"])] = {
            "testable": row.get("testable") == "TRUE",
            "retained_l1_fraction": parse_float(row.get("retained_l1_fraction")),
            "minimum_retained_l1_fraction": parse_float(
                row.get("minimum_retained_l1_fraction")) or 0.8,
        }
    return table


# --- Window construction -----------------------------------------------------------

WINDOW_FIELDS = (
    "arm",
    "axis_leakage",
    "layer",
    "aggregation_method",
    "feature_uid",
    "feature_display",
    "dataset",
    "axis_id",
    "window_id",
    "window_center",
    "window_lower",
    "window_upper",
    "lower_inclusive",
    "upper_inclusive",
    "interval_rule",
    "unit_of_observation",
    "n_samples_in_window",
    "n_samples_stage_complete",
    "n_samples_stage_missing",
    "n_donors",
    "n_donors_status",
    "mean_score_within_cohort_z",
    "se_score_naive_sample_level",
    "se_score_design_effect_sensitivity",
    "design_effect_source",
    "fibrosis_scale_native",
    "cross_cohort_stage_comparable",
    "histology_anchor_status",
    "sex_source",
    "feature_testable",
    "zero_occupancy_reason",
    "inference_use",
    "visualization_only",
    "longitudinal_interpretation_barred",
)


def build_windows(
    substrate: Substrate,
    layers: Sequence[FeatureScores],
    cohorts: Sequence[str],
    arm: str,
    testability: Mapping[tuple[str, str], Mapping[str, Any]],
) -> list[dict[str, Any]]:
    multiplier = design_effect_se_multiplier()
    members: dict[tuple[str, int], list[tuple[str, str]]] = {}
    for dataset in cohorts:
        cohort_keys = [key for key in substrate.axis if key[0] == dataset]
        for index in range(1, len(WINDOW_CENTERS) + 1):
            members[(dataset, index)] = [
                key for key in cohort_keys if in_window(substrate.axis[key], index)
            ]

    rows: list[dict[str, Any]] = []
    for layer in layers:
        for uid in sorted(layer.scores):
            feature_scores = layer.scores[uid]
            for dataset in cohorts:
                gate = testability.get((uid, dataset))
                for index in range(1, len(WINDOW_CENTERS) + 1):
                    lower, upper, final = window_bounds(index)
                    keys = members[(dataset, index)]
                    values = [feature_scores[key] for key in keys if key in feature_scores]
                    stat = summarize_window(values)
                    scored = [key for key in keys if key in feature_scores]
                    n_stage = sum(1 for key in scored if key in substrate.stage_complete)

                    reason = ""
                    if stat.n_samples == 0:
                        if gate is not None and not gate["testable"]:
                            reason = (
                                "program_not_testable_after_signature_exclusion_"
                                f"retained_l1_fraction_{gate['retained_l1_fraction']:.6f}"
                                f"_below_minimum_{gate['minimum_retained_l1_fraction']}"
                            )
                        else:
                            reason = "no_scored_samples_in_window"

                    rows.append({
                        "arm": arm,
                        "axis_leakage": AXIS_LEAKAGE[arm],
                        "layer": layer.layer,
                        "aggregation_method": layer.aggregation,
                        "feature_uid": uid,
                        "feature_display": layer.display.get(uid, uid),
                        "dataset": dataset,
                        "axis_id": AXIS_ID,
                        "window_id": index,
                        "window_center": WINDOW_CENTERS[index - 1],
                        "window_lower": lower,
                        "window_upper": upper,
                        "lower_inclusive": True,
                        "upper_inclusive": final,
                        "interval_rule": INTERVAL_RULE,
                        "unit_of_observation": "sequencing_sample",
                        "n_samples_in_window": stat.n_samples,
                        "n_samples_stage_complete": n_stage,
                        "n_samples_stage_missing": stat.n_samples - n_stage,
                        "n_donors": None,
                        "n_donors_status": "untestable_no_donor_key_on_disk",
                        "mean_score_within_cohort_z": stat.mean,
                        "se_score_naive_sample_level": stat.se,
                        "se_score_design_effect_sensitivity": (
                            None if stat.se is None else stat.se * multiplier
                        ),
                        "design_effect_source": DESIGN_EFFECT_SOURCE,
                        "fibrosis_scale_native": substrate.scale_native.get(dataset, ""),
                        "cross_cohort_stage_comparable": substrate.stage_comparable.get(
                            dataset, False),
                        "histology_anchor_status": HISTOLOGY_ANCHOR_STATUS[dataset],
                        "sex_source": substrate.sex_source.get(dataset, ""),
                        "feature_testable": None if gate is None else gate["testable"],
                        "zero_occupancy_reason": reason,
                        "inference_use": False,
                        "visualization_only": True,
                        "longitudinal_interpretation_barred": True,
                    })
    return rows


def reproduce_frozen_windows(
    computed: Sequence[Mapping[str, Any]],
    frozen_rows: Sequence[Mapping[str, str]],
    tolerance: float = 1e-9,
) -> dict[str, Any]:
    """Check the recomputed hotspot windows against the frozen upstream table.

    The frozen table is the existing 2,106-row hotspot window layer on the two
    evaluation cohorts.  Reproducing it is the control that says this script's window
    arithmetic matches the layer already on disk rather than merely resembling it.
    """

    frozen_index = {
        (row["program_uid"], row["dataset"], int(row["window_id"])): row
        for row in frozen_rows
    }
    checked = 0
    mismatches: list[dict[str, Any]] = []
    for row in computed:
        if row["layer"] != "hotspot_program":
            continue
        key = (row["feature_uid"], row["dataset"], int(row["window_id"]))
        frozen = frozen_index.get(key)
        if frozen is None:
            continue
        checked += 1
        expected_n = parse_int(frozen.get("n_participants"))
        expected_mean = parse_float(frozen.get("mean_score"))
        expected_se = parse_float(frozen.get("se_score"))
        problems = []
        if expected_n != row["n_samples_in_window"]:
            problems.append(f"n {expected_n} vs {row['n_samples_in_window']}")
        for label, expected, observed in (
            ("mean", expected_mean, row["mean_score_within_cohort_z"]),
            ("se", expected_se, row["se_score_naive_sample_level"]),
        ):
            if expected is None and observed is None:
                continue
            if expected is None or observed is None:
                problems.append(f"{label} {expected} vs {observed}")
            elif abs(expected - observed) > tolerance * max(1.0, abs(expected)):
                problems.append(f"{label} {expected} vs {observed}")
        if problems:
            mismatches.append({"key": "::".join(str(part) for part in key),
                               "problems": "; ".join(problems)})
    return {"n_checked": checked, "n_mismatched": len(mismatches),
            "mismatches": mismatches[:20]}


# --- Output ------------------------------------------------------------------------


def format_cell(value: Any) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, float):
        return repr(value)
    return str(value)


def write_tsv(path: Path, fieldnames: Sequence[str],
              rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_cell(row.get(key)) for key in fieldnames})


def build_axis_positions(substrate: Substrate,
                         clinical_rows: Sequence[Mapping[str, str]]) -> list[dict[str, Any]]:
    arm_of = {cohort: ARM_PRIMARY for cohort in PRIMARY_COHORTS}
    arm_of.update({cohort: ARM_QUARANTINED for cohort in QUARANTINED_COHORTS})
    rows = []
    for entry in clinical_rows:
        dataset = entry["dataset"]
        key = (dataset, entry["sample_id"])
        arm = arm_of[dataset]
        rows.append({
            "analysis_unit_id": entry["analysis_unit_id"],
            "sample_id": entry["sample_id"],
            "sra_run_accession": entry["sra_run_accession"],
            "dataset": dataset,
            "unit_of_observation": "sequencing_sample",
            "n_donors": None,
            "n_donors_status": "untestable_no_donor_key_on_disk",
            "arm": arm,
            "axis_leakage": AXIS_LEAKAGE[arm],
            "axis_id": AXIS_ID,
            "axis_percentile_within_cohort": substrate.axis.get(key),
            "percentile_denominator_n_samples": substrate.cohort_n[dataset],
            "fibrosis_stage_harmonised": parse_int(entry.get("fibrosis_rederived_value")),
            "fibrosis_stage_native": entry.get("fibrosis_native_value", ""),
            "fibrosis_scale_native": entry["fibrosis_scale_native"],
            "cross_cohort_stage_comparable": entry["cross_cohort_stage_comparable"] == "TRUE",
            "histology_anchor_status": HISTOLOGY_ANCHOR_STATUS[dataset],
            "nas_score": parse_int(entry.get("nas_rederived_value")),
            "sex_source": entry["sex_source"],
            "clinical_join_verified": True,
        })
    return rows


def build_n_ledger(substrate: Substrate,
                   axis_positions: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """One row per reported quantity, naming its denominator and what fell out."""

    def count(predicate) -> int:
        return sum(1 for row in axis_positions if predicate(row))

    primary = lambda row: row["dataset"] in PRIMARY_COHORTS  # noqa: E731
    quarantined = lambda row: row["dataset"] in QUARANTINED_COHORTS  # noqa: E731

    return [
        {
            "quantity": "substrate_total",
            "n_samples": len(axis_positions),
            "unit_of_observation": "sequencing_sample",
            "excluded_from_total": 0,
            "exclusion_basis": "none",
            "note": "never report this as participants or donors",
        },
        {
            "quantity": "primary_arm_window_occupancy",
            "n_samples": count(primary),
            "unit_of_observation": "sequencing_sample",
            "excluded_from_total": len(axis_positions) - count(primary),
            "exclusion_basis": (
                "source_overlap_cohorts_held_out_by_frozen_contract_evaluation_cohorts"
            ),
            "note": "windows need only an axis position, so stage-less samples are included",
        },
        {
            "quantity": "primary_arm_stage_complete",
            "n_samples": count(lambda r: primary(r) and r["fibrosis_stage_harmonised"] is not None),
            "unit_of_observation": "sequencing_sample",
            "excluded_from_total": count(
                lambda r: primary(r) and r["fibrosis_stage_harmonised"] is None),
            "exclusion_basis": (
                "GSE162694_normal_liver_histology_maps_to_missing_not_stage_zero"
            ),
            "note": "this is the denominator for anything stage-adjusted or stage-stratified",
        },
        {
            "quantity": "quarantined_arm_window_occupancy",
            "n_samples": count(quarantined),
            "unit_of_observation": "sequencing_sample",
            "excluded_from_total": len(axis_positions) - count(quarantined),
            "exclusion_basis": "primary_arm_cohorts",
            "note": "never pool with the primary arm; the axis saw these cohorts' source",
        },
        {
            "quantity": "nas_complete_all_cohorts",
            "n_samples": count(lambda r: r["nas_score"] is not None),
            "unit_of_observation": "sequencing_sample",
            "excluded_from_total": count(lambda r: r["nas_score"] is None),
            "exclusion_basis": "GSE213621_and_GSE126848_never_recorded_NAS",
            "note": "GSE213621 contributes zero NAS despite being 43% of the substrate",
        },
        {
            "quantity": "cross_cohort_stage_comparable",
            "n_samples": count(
                lambda r: r["cross_cohort_stage_comparable"]
                and r["fibrosis_stage_harmonised"] is not None),
            "unit_of_observation": "sequencing_sample",
            "excluded_from_total": count(
                lambda r: not r["cross_cohort_stage_comparable"]
                or r["fibrosis_stage_harmonised"] is None),
            "exclusion_basis": (
                "GSE213621_collapsed_0_3_scale_is_not_the_Kleiner_0_4_scale"
            ),
            "note": (
                "GSE213621 level 3 pools F3 with F4 and is its ceiling, so a pooled "
                "factor aligns it against the other cohorts' second-from-top level"
            ),
        },
        {
            "quantity": "upstream_frozen_window_n_participants_column",
            "n_samples": count(primary),
            "unit_of_observation": "sequencing_sample",
            "excluded_from_total": count(
                lambda r: primary(r) and r["fibrosis_stage_harmonised"] is None),
            "exclusion_basis": "none_the_upstream_column_excludes_nothing",
            "note": (
                "the frozen hotspot_fixed_windows n_participants column counts every "
                "sample in the window including the stage-less controls, and names them "
                "participants; recorded here so the number is not quoted as a donor count"
            ),
        },
    ]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--provenance-artifact", required=True,
                        help="job 1 rna-continuum-clinical-provenance directory")
    parser.add_argument("--provenance-sha256", required=True,
                        help="digest of the job 1 ARTIFACTS.json")
    parser.add_argument("--axis-scores", required=True)
    parser.add_argument("--hotspot-scores", required=True)
    parser.add_argument("--system-scores", required=True)
    parser.add_argument("--program-testability", required=True)
    parser.add_argument("--frozen-hotspot-windows", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    provenance = Path(args.provenance_artifact).resolve(strict=True)
    observed = sha256_file(provenance / "ARTIFACTS.json")
    if observed != args.provenance_sha256:
        raise WindowError(
            "clinical provenance artifact digest mismatch: "
            f"expected {args.provenance_sha256}, observed {observed}"
        )
    provenance_audit = json.loads((provenance / "provenance" / "audit.json").read_text())
    if provenance_audit["n_disagreements"] != 0:
        raise WindowError("the bound clinical provenance artifact did not reproduce")

    clinical_rows = read_delimited(provenance / "tables" / "clinical_join_reproduction.tsv")
    axis_rows = read_delimited(Path(args.axis_scores))
    substrate = load_substrate(axis_rows, clinical_rows)

    layers = load_hotspot_scores(read_delimited(Path(args.hotspot_scores)))
    layers += load_system_scores(read_delimited(Path(args.system_scores)))
    testability = load_testability(read_delimited(Path(args.program_testability)))

    primary = build_windows(substrate, layers, PRIMARY_COHORTS, ARM_PRIMARY, testability)
    quarantined = build_windows(substrate, layers, QUARANTINED_COHORTS,
                                ARM_QUARANTINED, testability)

    # The expected shape is computable before the table exists, so a shortfall is a
    # defect the job can catch itself rather than one a reader has to find with wc -l.
    # A feature that is present in the source but unscored must appear at zero
    # occupancy; it may never be silently absent.
    shape_checks = []
    for arm_label, rows, cohorts in (
        (ARM_PRIMARY, primary, PRIMARY_COHORTS),
        (ARM_QUARANTINED, quarantined, QUARANTINED_COHORTS),
    ):
        expected = expected_row_count(layers, len(cohorts))
        shape_checks.append({"arm": arm_label, "rows_expected": expected,
                             "rows_emitted": len(rows), "n_cohorts": len(cohorts)})
        if len(rows) != expected:
            raise WindowError(
                f"{arm_label}: emitted {len(rows)} window rows but the layer shape "
                f"requires {expected} (features x {len(WINDOW_CENTERS)} windows x "
                f"{len(cohorts)} cohorts); a feature was dropped instead of being "
                "recorded at zero occupancy"
            )
        for layer in layers:
            emitted = {row["feature_uid"] for row in rows
                       if row["layer"] == layer.layer
                       and row["aggregation_method"] == layer.aggregation}
            if emitted != set(layer.scores):
                missing = sorted(set(layer.scores) - emitted)
                raise WindowError(
                    f"{arm_label}/{layer.layer}/{layer.aggregation}: "
                    f"{len(missing)} source features are absent from the table: {missing[:5]}"
                )

    frozen_windows = read_delimited(Path(args.frozen_hotspot_windows))
    control = reproduce_frozen_windows(primary, frozen_windows)
    if control["n_mismatched"]:
        raise WindowError(
            f"recomputed windows disagree with the frozen upstream table in "
            f"{control['n_mismatched']} of {control['n_checked']} rows: "
            f"{control['mismatches'][:3]}"
        )
    if control["n_checked"] == 0:
        raise WindowError("the frozen upstream window table produced no control rows")

    axis_positions = build_axis_positions(substrate, clinical_rows)
    ledger = build_n_ledger(substrate, axis_positions)

    output = Path(args.output)
    write_tsv(output / "tables" / "window_summaries_primary.tsv", WINDOW_FIELDS, primary)
    write_tsv(output / "tables" / "window_summaries_quarantined.tsv", WINDOW_FIELDS,
              quarantined)
    write_tsv(output / "tables" / "axis_positions.tsv",
              tuple(axis_positions[0].keys()), axis_positions)
    write_tsv(output / "tables" / "analysis_n_ledger.tsv",
              ("quantity", "n_samples", "unit_of_observation", "excluded_from_total",
               "exclusion_basis", "note"), ledger)

    census = []
    for dataset in sorted(substrate.cohort_n):
        arm = ARM_PRIMARY if dataset in PRIMARY_COHORTS else ARM_QUARANTINED
        rows = [row for row in axis_positions if row["dataset"] == dataset]
        census.append({
            "dataset": dataset,
            "unit_of_observation": "sequencing_sample",
            "n_samples": len(rows),
            "n_donors": None,
            "n_donors_status": "untestable_no_donor_key_on_disk",
            "donor_key_fields_searched": "subject;patient;individual;isolate;donor;title",
            "n_fibrosis_complete": sum(
                1 for row in rows if row["fibrosis_stage_harmonised"] is not None),
            "n_nas_complete": sum(1 for row in rows if row["nas_score"] is not None),
            "arm": arm,
            "axis_leakage": AXIS_LEAKAGE[arm],
            "fibrosis_scale_native": substrate.scale_native.get(dataset, ""),
            "cross_cohort_stage_comparable": substrate.stage_comparable.get(dataset, False),
            "histology_anchor_status": HISTOLOGY_ANCHOR_STATUS[dataset],
            "sex_source": substrate.sex_source.get(dataset, ""),
        })
    write_tsv(output / "tables" / "substrate_census.tsv", tuple(census[0].keys()), census)

    relabel = [
        {
            "upstream_path": str(Path(args.frozen_hotspot_windows)),
            "upstream_sha256": sha256_file(Path(args.frozen_hotspot_windows)),
            "upstream_label": "n_participants",
            "relabelled_as": "n_samples_in_window",
            "basis": "no participant_id in the substrate manifest; SRR:BioSample is 1:1",
        },
        {
            "upstream_path": str(Path(args.hotspot_scores)),
            "upstream_sha256": sha256_file(Path(args.hotspot_scores)),
            "upstream_label": "hotspot_participant_scores",
            "relabelled_as": "sample_level_hotspot_program_scores",
            "basis": "the table is keyed on sample_id and analysis_unit_id, not a donor",
        },
        {
            "upstream_path": str(Path(args.system_scores)),
            "upstream_sha256": sha256_file(Path(args.system_scores)),
            "upstream_label": "participant_system_scores",
            "relabelled_as": "sample_level_molecular_system_scores",
            "basis": "the table is keyed on sample_id and analysis_unit_id, not a donor",
        },
        {
            "upstream_path": str(Path(args.axis_scores)),
            "upstream_sha256": sha256_file(Path(args.axis_scores)),
            "upstream_label": "participant_scores",
            "relabelled_as": "sample_level_axis_positions",
            "basis": "the table is keyed on sample_id, not a donor",
        },
    ]
    write_tsv(output / "provenance" / "upstream_unit_relabel.tsv",
              ("upstream_path", "upstream_sha256", "upstream_label", "relabelled_as",
               "basis"), relabel)

    shape = feature_shape(layers)
    write_tsv(output / "tables" / "layer_shape_receipt.tsv",
              ("layer", "aggregation_method", "features_expected", "features_emitted",
               "features_with_no_finite_score"), shape)

    audit = {
        "unit_of_observation": "sequencing_sample",
        "donor_level_claim_made": False,
        "n_donors_status": "untestable_no_donor_key_on_disk",
        "predictive_arm_built": False,
        "inference_use": False,
        "axis_id": AXIS_ID,
        "arms_pooled": False,
        "n_window_rows_primary": len(primary),
        "n_window_rows_quarantined": len(quarantined),
        "primary_cohorts": list(PRIMARY_COHORTS),
        "quarantined_cohorts": list(QUARANTINED_COHORTS),
        "reproduction_control": control,
        "provenance_artifact_sha256": observed,
        "shape_checks": shape_checks,
        "layer_shape": shape,
        "features_expected": sum(item["features_expected"] for item in shape),
        "features_emitted": sum(item["features_emitted"] for item in shape),
        "supersedes": "rna-continuum-descriptive-windows-21114163",
        "supersedes_reason": (
            "21114163 registered features by parsed score, so "
            "hotspot_fibroblasts_8584ba834d31e608 (Activated stellate, PDGFRA) -- whose "
            "844 source rows all carry outcome_z NA because the signature-exclusion "
            "testability gate leaves it unscored at 72.0% retained L1 against a 0.8 "
            "minimum -- was absent from both arms rather than present at zero occupancy. "
            "Primary was 3636 rows against 3654 expected (short 9 windows x 2 cohorts) "
            "and quarantined 5454 against 5481 (short 9 x 3). This run registers "
            "features from the source rows that mention them and fails closed on any "
            "shape shortfall. 21114163 is preserved, not overwritten."
        ),
        "ledger": ledger,
    }
    (output / "provenance").mkdir(parents=True, exist_ok=True)
    (output / "provenance" / "audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(json.dumps({
        "n_window_rows_primary": len(primary),
        "n_window_rows_quarantined": len(quarantined),
        "features_expected": audit["features_expected"],
        "features_emitted": audit["features_emitted"],
        "shape_checks": shape_checks,
        "reproduction_control_checked": control["n_checked"],
        "reproduction_control_mismatched": control["n_mismatched"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
