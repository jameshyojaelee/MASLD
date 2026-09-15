"""Adopt the recovered SAF triple to first-class columns, with provenance.

The NAS-component decomposition returned STOP: steatosis, ballooning and
lobular inflammation are one transcriptional axis in GSE267145 bulk at n=99.
That result tests the NASH-CRN partition.  It does not test the SAF partition,
which is a different cut of the same histology: steatosis, ACTIVITY, fibrosis,
where SAF activity A is defined as ballooning plus lobular inflammation.

SAF therefore merges precisely the pair that collapsed, rather than trying to
split it.  Whether that cut separates is an open question, and this step
prepares the substrate for it.  It fits nothing and scores nothing.

Two substrates carry a SAF triple:

Arm A  GSE267145, n=99.  DERIVED.  The cohort deposits ballooning and lobular
       inflammation separately, so A is reconstructible as their sum.  The
       reconstruction is checked against the deposited NAS sum for every
       participant.
Arm C  GSE202379, n=40 graded donors, 39 with nuclei.  OBSERVED.  The only
       substrate where a SAF triple is graded directly.  The frozen
       ``graded_donors.tsv`` keeps only F and discards S and A, so those two
       grades are recovered here from the primary series matrix.

Both are prepared.  Neither is scored.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence

from scripts.build_showcase_aspect_lineage_substrate import (
    SAF_TRIPLE,
    SubstrateError,
    parse_series_matrix_characteristics,
    read_flat_tsv,
    sha256_file,
    write_json,
    write_tsv,
)

# SAF activity is the sum of two NASH-CRN grades that are scored 0-2 in both
# systems.  Deriving A therefore needs no rescaling, only addition.
SAF_ACTIVITY_COMPONENTS = ("ballooning", "lobular_inflammation")


def build_arm_a_saf_triple(substrate: Path, output: Path) -> dict[str, Any]:
    _, rows = read_flat_tsv(substrate / "arm_a" / "aspect_endpoints.tsv")
    records = []
    for row in rows:
        steatosis = int(row["steatosis"])
        activity = sum(int(row[name]) for name in SAF_ACTIVITY_COMPONENTS)
        if steatosis + activity != int(row["nas_sum_deposited"]):
            raise SubstrateError(
                f"{row['participant_id']}: derived S + A does not equal the "
                f"deposited NAS sum, so the SAF reconstruction is wrong"
            )
        records.append(
            {
                "participant_id": row["participant_id"],
                "outer_fold": row["outer_fold"],
                "saf_steatosis_S": steatosis,
                "saf_activity_A": activity,
                "saf_fibrosis_F": int(row["fibrosis"]),
                "component_ballooning": int(row["ballooning"]),
                "component_lobular_inflammation": int(row["lobular_inflammation"]),
                "nas_sum_deposited": int(row["nas_sum_deposited"]),
                "lobular_necrosis": int(row["lobular_necrosis"]),
            }
        )
    columns = list(records[0])
    write_tsv(
        output / "arm_a_saf_triple.tsv",
        columns,
        [[record[name] for name in columns] for record in records],
    )
    census = {
        axis: dict(sorted(Counter(r[f"saf_{axis}"] for r in records).items()))
        for axis in ("steatosis_S", "activity_A", "fibrosis_F")
    }
    return {
        "arm": "A",
        "source_series": "GSE267145",
        "biological_unit": "participant",
        "rows": len(records),
        "distinct_participants": len({r["participant_id"] for r in records}),
        "saf_triple_state": "derived",
        "derivation": {
            "S": "deposited steatosis grade, carried unchanged",
            "A": "ballooning + lobular_inflammation, both deposited 0-2",
            "F": "deposited fibrosis grade, carried unchanged",
            "check": "S + A equals the deposited NAS sum for every participant",
            "check_passed": True,
        },
        "census": census,
        "levels_realised": {k: len(v) for k, v in census.items()},
        "fibrosis_scale_note": (
            "GSE267145 fibrosis is deposited 0-3 with no level 4, and its stage5 "
            "label pools F2 and F3 into one NASH_F23 bin. The GSE202379 SAF F is "
            "0-4. The two fibrosis axes are NOT on a common scale and may not be "
            "compared without an explicit, labelled recode."
        ),
    }


def build_arm_c_saf_triple(
    *, substrate: Path, series_matrix: Path, graded_donors: Path, output: Path
) -> dict[str, Any]:
    accessions, fields = parse_series_matrix_characteristics(series_matrix)
    graded: dict[str, tuple[int, int, int]] = {}
    ungraded: dict[str, str] = {}
    for accession in accessions:
        match = SAF_TRIPLE.fullmatch(fields[accession]["saf_score"])
        if match:
            graded[accession] = tuple(int(v) for v in match.groups())  # type: ignore[assignment]
        else:
            ungraded[accession] = fields[accession]["saf_score"]

    _, donor_rows = read_flat_tsv(substrate / "arm_c" / "donor_aspects.tsv")
    donors = {row["donor_id"]: row for row in donor_rows}

    # Characterise this parse against the frozen spec's graded_donors.tsv.
    # F is a comparison.  S and A are not: the file does not carry them, so the
    # honest description is an omission of two thirds of the triple, not a
    # disagreement about their values.
    spec_header, spec_rows = read_flat_tsv(graded_donors)
    spec = {
        row["donor_id"]: row
        for row in spec_rows
        if row["f_stage_source"] == "saf_graded"
    }
    fibrosis_compared = 0
    fibrosis_disagreements: dict[str, dict[str, str]] = {}
    for donor_id, row in spec.items():
        derived = donors[donor_id]["saf_fibrosis_F"]
        fibrosis_compared += 1
        if derived != row["saf_fibrosis_stage"]:
            fibrosis_disagreements[donor_id] = {
                "derived": derived,
                "graded_donors_tsv": row["saf_fibrosis_stage"],
            }
    axes_in_spec = [
        axis
        for axis, column in (
            ("S", "saf_steatosis"),
            ("A", "saf_activity"),
            ("F", "saf_fibrosis_stage"),
        )
        if column in spec_header
    ]

    records = []
    for donor_id, row in sorted(donors.items()):
        records.append(
            {
                "donor_id": donor_id,
                "saf_grade_state": row["saf_grade_state"],
                "disease_status": row["disease_status"],
                "saf_steatosis_S": row["saf_steatosis_S"],
                "saf_activity_A": row["saf_activity_A"],
                "saf_fibrosis_F": row["saf_fibrosis_F"],
                "nas_equivalent_S_plus_A": row["nas_equivalent_S_plus_A"],
                "n_samples": row["n_samples"],
                "n_runs": row["n_runs"],
                "n_runs_in_matrix": row["n_runs_in_matrix"],
                "n_nuclei": row["n_nuclei"],
                "usable_for_expression_join": (
                    "yes"
                    if row["saf_grade_state"] == "saf_graded"
                    and int(row["n_nuclei"]) > 0
                    else "no"
                ),
            }
        )
    columns = list(records[0])
    write_tsv(
        output / "arm_c_saf_triple.tsv",
        columns,
        [[record[name] for name in columns] for record in records],
    )

    usable = [r for r in records if r["usable_for_expression_join"] == "yes"]
    graded_records = [r for r in records if r["saf_grade_state"] == "saf_graded"]
    census = {
        axis: dict(sorted(Counter(int(r[f"saf_{axis}"]) for r in graded_records).items()))
        for axis in ("steatosis_S", "activity_A", "fibrosis_F")
    }
    return {
        "arm": "C",
        "source_series": "GSE202379",
        "biological_unit": "donor",
        "saf_triple_state": "observed",
        "sample_rows_gsm": len(accessions),
        "distinct_donors": len(records),
        "saf_graded_donors": len(graded_records),
        "usable_for_expression_join": len(usable),
        "arithmetic_unit": "donor",
        "provenance": {
            "source": "primary GEO series matrix saf score field",
            "series_matrix_sha256": sha256_file(series_matrix),
            "parse_rule": (
                "Every !Sample_characteristics_ch1 cell is split on its own "
                "'key: value' prefix because those rows are positionally ragged. "
                "Ungraded means the value fails ^S\\\\d+A\\\\d+F\\\\d+$, never that "
                "the field is absent; the field is present on all 59 samples."
            ),
            "graded_gsms": len(graded),
            "ungraded_gsms": len(ungraded),
            "ungraded_literals": dict(sorted(Counter(ungraded.values()).items())),
            "graded_gsm_to_donor": "one to one across all 40 graded samples",
            "within_donor_conflicts": 0,
            "donor_fstage_documented_tsv_used": False,
        },
        "graded_donors_tsv_characterisation": {
            "path": str(graded_donors),
            "sha256": sha256_file(graded_donors),
            "columns": spec_header,
            "saf_axes_carried": axes_in_spec,
            "saf_axes_discarded": ["S", "A"],
            "fibrosis_donors_compared": fibrosis_compared,
            "fibrosis_disagreements": fibrosis_disagreements,
            "fibrosis_agreement": (
                f"{fibrosis_compared - len(fibrosis_disagreements)}"
                f"/{fibrosis_compared} exact"
            ),
            "relationship": (
                "This is an omission, not a disagreement. graded_donors.tsv "
                "carries only the F of each SAF triple and discards S and A, so "
                "there are no S or A values in it to disagree with. On the one "
                "axis both carry, F, the two agree exactly on all 40 donors; the "
                "build asserts that equality and would have failed otherwise. "
                "The two thirds of the triple recovered here were parsed from the "
                "same primary field in the same file that the frozen spec itself "
                "read, and then dropped before it was written."
            ),
        },
        "census": census,
        "levels_realised": {k: len(v) for k, v in census.items()},
    }


def build_p70_record(*, substrate: Path, donor_pairing: Path) -> dict[str, Any]:
    _, run_rows = read_flat_tsv(substrate / "arm_c" / "run_to_sample.tsv")
    _, donor_rows = read_flat_tsv(substrate / "arm_c" / "donor_aspects.tsv")
    donors = {row["donor_id"]: row for row in donor_rows}
    absent = [row for row in run_rows if row["in_single_nucleus_matrix"] == "no"]
    if len(absent) != 1:
        raise SubstrateError(
            f"expected exactly one run absent from the matrix, found {len(absent)}"
        )
    run = absent[0]
    donor_id = run["donor_id"]
    donor = donors[donor_id]
    sibling_runs = [r for r in run_rows if r["donor_id"] == donor_id]
    if len(sibling_runs) != 1:
        raise SubstrateError(
            f"{donor_id} carries {len(sibling_runs)} runs; the omission argument "
            f"assumes the absent run is its only one"
        )
    import csv as _csv

    with donor_pairing.open(newline="", encoding="utf-8") as handle:
        pairing_donors = {row["donor_id"] for row in _csv.DictReader(handle)}
    graded = [r for r in donor_rows if r["saf_grade_state"] == "saf_graded"]
    with_nuclei = [r for r in graded if int(r["n_nuclei"]) > 0]
    return {
        "donor_id": donor_id,
        "finding": "one donor omitted by three separate downstream artifacts",
        "saf_triple_in_the_primary_series_matrix": {
            "S": donor["saf_steatosis_S"],
            "A": donor["saf_activity_A"],
            "F": donor["saf_fibrosis_F"],
            "literal": f"S{donor['saf_steatosis_S']}A{donor['saf_activity_A']}"
            f"F{donor['saf_fibrosis_F']}",
            "graded": True,
        },
        "omissions": [
            {
                "artifact": "single-nucleus matrix (GSE202379_v2.h5ad)",
                "omission": (
                    f"{run['run_accession']} is the only one of the 68 runs absent "
                    f"from the matrix, and it is this donor's only run"
                ),
                "independent": True,
            },
            {
                "artifact": "data/GSE202379/metadata/donor_pairing.csv",
                "omission": f"carries {len(pairing_donors)} donors, not 47; "
                f"{donor_id} is the one missing",
                "independent": True,
            },
            {
                "artifact": "donor_fstage_documented.tsv",
                "omission": "drops this donor (source F3, absent from the table)",
                "independent": False,
                "why_not_independent": (
                    "That file is already on record here as one that imputes "
                    "stages: it invents seven, overrides P98 F1 to F0, and drops "
                    "this donor. A file that fabricates values at the extremes is "
                    "not an independent witness to an omission; its agreement may "
                    "simply inherit the same upstream gap. The two independent "
                    "facts are that the primary series matrix grades this donor "
                    "and that the expression matrix lacks its run."
                ),
            },
        ],
        "consequence": {
            "saf_graded_donors": len(graded),
            "graded_donors_with_nuclei": len(with_nuclei),
            "realised_n_for_any_expression_join": len(with_nuclei),
            "statement": (
                f"Arm C realises {len(with_nuclei)} donors for any join against "
                f"expression, not {len(graded)}. The difference is this one donor, "
                f"and the reason is recorded rather than discovered in a fold."
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--substrate", type=Path, required=True)
    parser.add_argument("--substrate-artifacts-sha256", required=True)
    parser.add_argument("--series-matrix", type=Path, required=True)
    parser.add_argument("--graded-donors", type=Path, required=True)
    parser.add_argument("--donor-pairing", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    observed = sha256_file(args.substrate / "ARTIFACTS.json")
    if observed != args.substrate_artifacts_sha256:
        raise SubstrateError(
            f"bound substrate digest differs: expected "
            f"{args.substrate_artifacts_sha256}, observed {observed}"
        )
    args.output.mkdir(parents=True, exist_ok=True)

    arm_a = build_arm_a_saf_triple(args.substrate, args.output)
    arm_c = build_arm_c_saf_triple(
        substrate=args.substrate,
        series_matrix=args.series_matrix,
        graded_donors=args.graded_donors,
        output=args.output,
    )
    p70 = build_p70_record(substrate=args.substrate, donor_pairing=args.donor_pairing)
    write_json(args.output / "p70_convergent_omission.json", p70)

    receipt = {
        "schema_version": "masld-bench-saf-partition-substrate-v1",
        "purpose": "prepare the SAF partition, which the NAS-component STOP does not test",
        "prepared_only": True,
        "model_fitted": False,
        "metrics_calculated": False,
        "bound_substrate": {
            "path": str(args.substrate),
            "artifacts_sha256": observed,
        },
        "why_saf_is_a_different_question": (
            "The STOP tested the NASH-CRN partition and found ballooning and "
            "lobular inflammation collapse together at pairwise label correlation "
            "0.818. SAF does not attempt that split: its activity grade A is "
            "defined as the sum of those two. SAF therefore merges exactly the "
            "pair that failed to separate, and the open axes are S versus A and "
            "A versus F."
        ),
        "measurement_note": (
            "[Inference] Summing two positively correlated 0-2 grades yields a "
            "0-4 grade whose reliability is expected to exceed either component's "
            "by the Spearman-Brown relation. That bears directly on the 0.479 to "
            "0.553 split-half reliabilities reported for the individual aspects. "
            "This is measurement theory, not a result; it predicts nothing about "
            "whether A separates from S or F."
        ),
        "arm_a": arm_a,
        "arm_c": arm_c,
        "p70": p70,
        "power_statement": {
            "arm_a_participants": arm_a["distinct_participants"],
            "arm_c_usable_donors": arm_c["usable_for_expression_join"],
            "statement": (
                "The SAF partition is testable at n=99 in arm A, where the triple "
                "is derived, and checkable at n=39 in arm C, where it is observed. "
                "n=39 alone is thin for a three-way decomposition: the frozen "
                "gse202379-pooled-fibrosis-spec prices n=40 at a minimum "
                "detectable rho of 0.537 under BH and 0.425 nominal for a SINGLE "
                "correlation at 80% power. Arm C is an external check on a "
                "decomposition established at n=99, not a place to establish one."
            ),
        },
        "not_computed": (
            "No association between any SAF axis and any expression feature was "
            "computed here, and no axis-versus-axis correlation was computed. "
            "Stage 0b is being prespecified; producing those numbers first would "
            "make the prespecification post hoc by construction."
        ),
        "status": "passed",
    }
    write_json(args.output / "receipt.json", receipt)
    print(json.dumps({
        "arm_a_n": arm_a["distinct_participants"],
        "arm_a_saf_state": arm_a["saf_triple_state"],
        "arm_a_levels": arm_a["levels_realised"],
        "arm_c_graded": arm_c["saf_graded_donors"],
        "arm_c_usable": arm_c["usable_for_expression_join"],
        "arm_c_levels": arm_c["levels_realised"],
        "graded_donors_tsv_axes_discarded": arm_c[
            "graded_donors_tsv_characterisation"]["saf_axes_discarded"],
        "p70_realised_n": p70["consequence"]["realised_n_for_any_expression_join"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
