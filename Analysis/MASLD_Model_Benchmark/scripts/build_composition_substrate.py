"""Join lineage composition to each axis arm and derive the eligible family.

Stage 3, job 1 of 3. No outcome-versus-expression association is computed here.

**The three arms are not on one quantification.** GSE135251 and GSE130970 are
BayesPrism on kallisto counts (77,078 features, first key ``TSPAN6``);
**GSE193066 is STAR** (37,606 features, first key ``A1BG``). It carries no
``_star_backup`` because it was never in the kallisto driver's ``DATASETS``,
which is exactly the five cohorts that do have one. So its canonical proportions
sit on the same footing as the other cohorts' backups. The within-arm pairing --
proportions and expression from the same quantification -- holds in every arm,
which is what the covariate validity depends on; the uniformity does not.

**That makes the detection floor partly a pipeline property.** STAR outputs
contain no value below 1e-8 at all while kallisto reaches 1e-22 here, so a
lineage can fall under a 1e-4 floor for a quantification reason rather than a
biological one. The eligible family is therefore derived on the **two kallisto
arms only**, and GSE193066's detectability is reported separately and per
quantification. Evaluating one rule across incomparable substrates is the
asymmetric-denominator shape this project has a standing rule about.

Three families are recorded side by side so the correction is visible rather
than absorbed: the family the fixed rule gives on the kallisto arms, the family
the confounded pooled rule would have given, and the earlier five that were
derived from the ``star_backup`` table before the variant pairing was resolved.
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


class CompositionError(RuntimeError):
    """Raised when the composition substrate does not hold."""


KALLISTO_ARMS = ("GSE135251", "GSE130970")
ALL_ARMS = ("GSE135251", "GSE130970", "GSE193066")
QUANTIFICATION = {
    "GSE135251": "kallisto", "GSE130970": "kallisto", "GSE193066": "STAR",
}
DETECTION_FLOOR = 1e-4
SAMPLE_FRACTION = 0.90
ROW_SUM_TOLERANCE = 1e-9

#: Recorded so the correction is inspectable, not to be used.
SUPERSEDED_FIVE = (
    "Hepatocytes", "Macrophages", "pDCs", "Plasma cells", "cDC1s",
)


def _import(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise CompositionError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_ragged(path: Path) -> tuple[list[str], list[str], np.ndarray, dict]:
    """Read a proportions table whose rows carry one field more than the header.

    The leading unnamed field is the sample key. The offset is measured and
    required to be exactly one; anything else is refused rather than guessed.
    """

    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    widths = Counter(len(r) for r in rows)
    if len(widths) != 1:
        raise CompositionError(f"{path.name}: inconsistent row widths {dict(widths)}")
    width = next(iter(widths))
    if width - len(header) != 1:
        raise CompositionError(
            f"{path.name}: header {len(header)} against data {width}; only an "
            "offset of exactly one is handled"
        )
    keys = [r[0] for r in rows]
    values = np.asarray([[float(x) for x in r[1:]] for r in rows], dtype=float)
    sums = values.sum(axis=1)
    if np.max(np.abs(sums - 1.0)) > ROW_SUM_TOLERANCE:
        raise CompositionError(f"{path.name}: rows do not sum to 1 within tolerance")
    return header, keys, values, {
        "header_fields": len(header),
        "data_fields": width,
        "measured_offset": 1,
        "rows_summing_to_exactly_1.0": int((sums == 1.0).sum()),
        "max_abs_deviation_from_1": float(np.max(np.abs(sums - 1.0))),
        "why_tolerance_not_equality": (
            "only a minority of rows are bitwise 1.0; an equality assertion "
            "fails on most of them while the composition is sound"
        ),
    }


def detectability(header: Sequence[str], values: np.ndarray) -> dict[str, float]:
    return {
        header[i]: float((values[:, i] >= DETECTION_FLOOR).mean())
        for i in range(len(header))
    }


def family_from(tables: dict[str, dict[str, float]], arms: Sequence[str]) -> list[str]:
    names = list(next(iter(tables.values())))
    return sorted(
        n for n in names if all(tables[a][n] >= SAMPLE_FRACTION for a in arms)
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extension", type=Path, required=True)
    parser.add_argument("--sealed-instrument", type=Path, required=True)
    parser.add_argument("--proportions-root", type=Path, required=True)
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--arm-b-source", type=Path, required=True)
    parser.add_argument("--stage1-substrate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    if arguments.output.exists():
        raise CompositionError("refusing to overwrite a composition substrate")
    mc = _import(arguments.extension, "multicovariate")
    analysis = mc.load_sealed(arguments.sealed_instrument)

    # ---- crosswalk, bound and asserted ----
    lines = arguments.crosswalk.read_text(encoding="utf-8").splitlines()
    head = lines[0].split("\t")
    rows = [dict(zip(head, l.split("\t"), strict=True)) for l in lines[1:]]
    gsm = [r["participant_id"] for r in rows]
    srr = [r["run_accession"] for r in rows]
    if not (len(gsm) == len(set(gsm)) == len(srr) == len(set(srr)) == 180):
        raise CompositionError("the bound crosswalk is not 1:1 over 180 rows")
    to_srr = dict(zip(gsm, srr))
    print(f"crosswalk bound and asserted 1:1 over {len(gsm)} rows", flush=True)

    # ---- proportions, canonical only ----
    tables, ragged, keyed = {}, {}, {}
    for cohort in ALL_ARMS:
        header, keys, values, audit = read_ragged(
            arguments.proportions_root / cohort / f"{cohort}_bayesprism_proportions.tsv")
        tables[cohort] = detectability(header, values)
        ragged[cohort] = audit
        keyed[cohort] = (header, keys, values)
        print(f"{cohort} ({QUANTIFICATION[cohort]}): {values.shape}", flush=True)

    lineages = list(keyed["GSE135251"][0])
    for cohort in ALL_ARMS[1:]:
        if list(keyed[cohort][0]) != lineages:
            raise CompositionError(f"{cohort} carries a different lineage roster")

    family = family_from(tables, KALLISTO_ARMS)
    pooled = family_from(tables, ALL_ARMS)
    print(f"family on the kallisto arms: {len(family)} -> {family}", flush=True)
    print(f"pooled (confounded) rule would give: {len(pooled)}", flush=True)

    # ---- variant sensitivity on the family itself ----
    variant = {}
    for cohort in KALLISTO_ARMS:
        backup = (arguments.proportions_root / cohort
                  / f"{cohort}_bayesprism_proportions_star_backup.tsv")
        if not backup.exists():
            continue
        h, _, v, _ = read_ragged(backup)
        variant[cohort] = detectability(h, v)
    backup_family = family_from(variant, tuple(variant)) if variant else []

    small_value = {}
    for cohort in ALL_ARMS:
        v = keyed[cohort][2]
        positive = v[v > 0]
        small_value[cohort] = {
            "quantification": QUANTIFICATION[cohort],
            "min_positive": float(positive.min()),
            "values_in_0_to_1e-8": int(((v > 0) & (v < 1e-8)).sum()),
        }

    # ---- join each arm and measure conditioning ----
    arms = {}
    (arguments.output / "arms").mkdir(mode=0o750, parents=True, exist_ok=True)
    for cohort in ALL_ARMS:
        header, keys, values = keyed[cohort]
        index = {k: i for i, k in enumerate(keys)}
        if cohort == "GSE135251":
            axis = [l.split("\t")[0] for l in (
                arguments.arm_b_source / "molecular" / "participant_axis.tsv"
            ).read_text(encoding="utf-8").splitlines()[1:]]
            wanted = [to_srr[g] for g in axis]
            join_note = "GSM axis joined through the bound crosswalk"
        else:
            wanted = (arguments.stage1_substrate / "arms"
                      / f"{cohort}_sample_axis.tsv").read_text(
                          encoding="utf-8").splitlines()[1:]
            axis = wanted
            join_note = "SRR axis joined directly; no derivation"
        missing = [w for w in wanted if w not in index]
        if missing:
            raise CompositionError(
                f"{cohort}: {len(missing)} of {len(wanted)} axis members have no "
                "proportion row")
        block = values[[index[w] for w in wanted], :]
        columns = [lineages.index(n) for n in family]
        proportions = block[:, columns]
        _, conditioning = mc.covariate_basis(
            analysis, [proportions[:, j] for j in range(proportions.shape[1])])
        arms[cohort] = {
            "quantification": QUANTIFICATION[cohort],
            "join": join_note,
            "axis_members": len(wanted),
            "joined": len(wanted) - len(missing),
            "join_fraction": 1.0,
            "family_conditioning_without_the_other_axis": conditioning,
            "ragged_audit": ragged[cohort],
            "detectability_of_the_family": {n: tables[cohort][n] for n in family},
        }
        np.save(arguments.output / "arms" / f"{cohort}_family_proportions.npy",
                proportions)
        (arguments.output / "arms" / f"{cohort}_sample_axis.tsv").write_text(
            "axis_member\tproportion_key\n"
            + "".join(f"{a}\t{w}\n" for a, w in zip(axis, wanted)), encoding="utf-8")
        print(f"{cohort}: joined {len(wanted)}/{len(wanted)}, "
              f"rank {conditioning['numerical_rank']}/{conditioning['n_covariates']}, "
              f"cond {conditioning['condition_number']:.2f}", flush=True)

    payload = {
        "schema_version": "masld-bench-composition-substrate-v1",
        "no_outcome_association_was_computed": True,
        "quantification_is_not_uniform_across_the_arms": {
            "per_cohort": QUANTIFICATION,
            "why_gse193066_has_no_backup": (
                "it was never in the kallisto driver's DATASETS, which is "
                "exactly the five cohorts that do have a _star_backup. Its "
                "canonical proportions sit on the same footing as the other "
                "cohorts' backups."
            ),
            "what_still_holds": (
                "the within-arm pairing: proportions and expression come from "
                "the same quantification in every arm, which is what the "
                "covariate validity depends on"
            ),
            "what_does_not": (
                "uniformity across arms. 'Canonical across all three' gives "
                "kallisto, kallisto, STAR, and must not be described as one "
                "quantification."
            ),
            "second_independent_ground_for_the_separate_arm_rule": (
                "the two external arms are quantified differently as well as "
                "graded differently"
            ),
        },
        "small_value_structure_is_a_pipeline_property": {
            "per_cohort": small_value,
            "reading": (
                "STAR outputs carry no value below 1e-8; kallisto reaches 1e-22 "
                "here. A lineage can fall under a 1e-4 floor for a "
                "quantification reason rather than a biological one."
            ),
        },
        "eligible_family": {
            "rule": (
                f"detectable at >= {DETECTION_FLOOR} in at least "
                f"{SAMPLE_FRACTION:.0%} of samples in every kallisto arm"
            ),
            "derived_on": list(KALLISTO_ARMS),
            "why_not_pooled": (
                "one arm is on a different quantification, so a pooled rule is "
                "confounded with the pipeline. Evaluating a rule across "
                "incomparable substrates is the asymmetric-denominator shape."
            ),
            "family": family,
            "size": len(family),
            "residual_df_with_the_other_axis": f"n - 2 - {len(family) + 1}",
        },
        "families_recorded_side_by_side": {
            "why": (
                "the family changed twice during design. Recording all three "
                "makes the correction inspectable rather than absorbed."
            ),
            "adopted_kallisto_rule": family,
            "pooled_confounded_rule": {
                "family": pooled,
                "size": len(pooled),
                "what_it_would_have_cost": sorted(set(family) - set(pooled)),
                "reading": (
                    "the confound was material, not cosmetic: it would have lost "
                    f"{len(set(family) - set(pooled))} of {len(family)} lineages "
                    "to STAR's detection floor rather than to biology"
                ),
            },
            "superseded_five_from_the_star_backup_table": {
                "family": list(SUPERSEDED_FIVE),
                "why_it_was_wrong": (
                    "derived from the GSE135251 _star_backup table before the "
                    "variant pairing was resolved. That variant pairs with no "
                    "frozen axis in this benchmark."
                ),
                "entered_but_absent_from_it": sorted(set(family) - set(SUPERSEDED_FIVE)),
                "present_in_it_but_dropped": sorted(set(SUPERSEDED_FIVE) - set(family)),
            },
        },
        "variant_sensitivity_on_the_family": {
            "is_a_gate": False,
            "what_it_shows": (
                "the eligible family is itself variant-dependent. Only "
                "Hepatocytes and Macrophages are stable across variants; the "
                "reassurance that the confirmatory family was not at risk was "
                "drawn from that stable minority."
            ),
            "family_under_the_backup_variant": backup_family,
            "family_under_canonical": family,
            "stable_across_variants": sorted(set(family) & set(backup_family)),
        },
        "arms": arms,
        "crosswalk": {
            "bound_not_derived": True,
            "rows": len(gsm),
            "one_to_one_both_directions_asserted": True,
        },
        "ragged_headers_are_systemic_in_this_tree": (
            "five file families now carry the same R write.table(row.names=TRUE) "
            "signature: bulk metadata, external metadata, the source metadata, "
            "the counts tables and these proportions. It is a property of the "
            "tree, not a per-file quirk."
        ),
        "a_manifest_must_record_what_it_excluded": (
            "a glob that silently drops a sibling variant produces a substrate "
            "that looks complete and is not"
        ),
    }
    (arguments.output / "composition_substrate.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    freeze_tree(arguments.output, {
        "artifact_class": "composition_substrate",
        "family_size": len(family),
        "quantification_uniform_across_arms": False,
        "outcome_association_computed": False,
        "status": "passed"})
    verify_frozen_tree(arguments.output)
    print(json.dumps({
        "family": family,
        "pooled_would_give": pooled,
        "backup_variant_family": backup_family,
        "conditioning": {c: round(a["family_conditioning_without_the_other_axis"][
            "condition_number"], 2) for c, a in arms.items()},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
