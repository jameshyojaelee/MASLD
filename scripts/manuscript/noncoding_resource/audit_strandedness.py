#!/usr/bin/env python3
"""Audit the frozen five-cohort count contract and RNA-seq strandedness.

This is a source gate. A failed cohort is recorded as a failure; the script does
not drop samples or relax the preregistered fourfold threshold.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import statistics
from pathlib import Path


COHORTS = ("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
STAR_COHORTS = tuple(x for x in COHORTS if x != "GSE162694")
MIN_RATIO = 4.0
AMENDED_POLICY = "cohort_reverse_orientation_v2"
STRICT_POLICY = "all_samples_fourfold_v1"
MIN_FOURFOLD_FRACTION = 0.90
MIN_ANY_REVERSE_RATIO = 1.0
MAX_GROUP_FAILURE_FRACTION = 0.10
MAX_GROUP_FAILURE_SPREAD = 0.05


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_star(path: Path) -> tuple[int, int, int]:
    totals = [0, 0, 0]
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not fields or fields[0].startswith("N_"):
                continue
            if len(fields) != 4:
                raise ValueError(f"Unexpected STAR gene-count row in {path}")
            for index in range(3):
                totals[index] += int(fields[index + 1])
    return tuple(totals)  # type: ignore[return-value]


def sample_id_from_star(path: Path) -> str:
    return path.name.split(".ReadsPerGene", 1)[0].split("_ReadsPerGene", 1)[0]


def parse_featurecounts_assigned(path: Path) -> dict[str, int]:
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        assigned = next(reader)
    if not assigned or assigned[0] != "Assigned":
        raise ValueError(
            f"Expected Assigned as first featureCounts summary row: {path}"
        )
    return {
        Path(column).name.split(".Aligned", 1)[0]: int(value)
        for column, value in zip(header[1:], assigned[1:])
    }


def count_source_rows(manifest: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with manifest.open() as handle:
        selected = {
            row["dataset"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["dataset"] in COHORTS
        }
    for cohort in COHORTS:
        source = selected.get(cohort)
        if source is None:
            rows.append(
                {
                    "dataset": cohort,
                    "count_path": "",
                    "sha256_matches": False,
                    "featurecounts_v2_1_1": False,
                    "reverse_stranded_s2": False,
                    "gencode_v49": False,
                    "status": "fail_missing_manifest_row",
                }
            )
            continue
        path = Path(source["count_path"])
        header = path.open().readline().rstrip("\n") if path.is_file() else ""
        checks = {
            "sha256_matches": path.is_file() and sha256(path) == source["sha256"],
            "featurecounts_v2_1_1": "featureCounts v2.1.1" in header,
            "reverse_stranded_s2": '"-s" "2"' in header,
            "gencode_v49": "gencode.v49.chr_patch_hapl_scaff.annotation.gtf" in header,
        }
        rows.append(
            {
                "dataset": cohort,
                "count_path": str(path),
                **checks,
                "status": "pass" if all(checks.values()) else "fail_count_contract",
            }
        )
    return rows


def amended_cohort_pass(row: dict[str, object]) -> bool:
    """Apply the approved cohort-level orientation rule without dropping samples."""

    n = int(row["n_audited"])
    if n <= 0:
        return False
    below = int(row["n_below_fourfold"])
    return (
        float(row["median_reverse_same_ratio"]) >= MIN_RATIO
        and (n - below) / n >= MIN_FOURFOLD_FRACTION
        and float(row["minimum_ratio"]) > MIN_ANY_REVERSE_RATIO
    )


def stage_concentration_audit(
    sample_rows: list[dict[str, object]], metadata_path: Path
) -> list[dict[str, object]]:
    """Check that lower-ratio GSE213621 libraries are not stage concentrated."""

    with metadata_path.open(newline="") as handle:
        metadata = {
            row["sample_id"]: row
            for row in csv.DictReader(handle)
            if row["dataset"] == "GSE213621"
        }
    target = [row for row in sample_rows if row["dataset"] == "GSE213621"]
    missing = sorted(
        str(row["sample_id"]) for row in target if row["sample_id"] not in metadata
    )
    if missing:
        raise ValueError(f"GSE213621 strand samples missing metadata: {missing[:5]}")

    groups: dict[str, list[bool]] = {}
    for row in target:
        condition = metadata[str(row["sample_id"])]["condition"]
        if not condition:
            raise ValueError(f"GSE213621 sample lacks condition: {row['sample_id']}")
        groups.setdefault(condition, []).append(not bool(row["pass_fourfold"]))

    rows: list[dict[str, object]] = []
    fractions: list[float] = []
    for condition in sorted(groups):
        flags = groups[condition]
        fraction = sum(flags) / len(flags)
        fractions.append(fraction)
        rows.append(
            {
                "dataset": "GSE213621",
                "condition": condition,
                "n_audited": len(flags),
                "n_below_fourfold": sum(flags),
                "below_fourfold_fraction": fraction,
                "group_fraction_gate_passed": fraction <= MAX_GROUP_FAILURE_FRACTION,
            }
        )
    spread = max(fractions) - min(fractions)
    overall_pass = (
        len(rows) == 4
        and all(bool(row["group_fraction_gate_passed"]) for row in rows)
        and spread <= MAX_GROUP_FAILURE_SPREAD
    )
    for row in rows:
        row["maximum_group_failure_fraction"] = MAX_GROUP_FAILURE_FRACTION
        row["observed_group_failure_spread"] = spread
        row["maximum_group_failure_spread"] = MAX_GROUP_FAILURE_SPREAD
        row["stage_concentration_gate_passed"] = overall_pass
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--featurecounts-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--policy",
        choices=(STRICT_POLICY, AMENDED_POLICY),
        default=STRICT_POLICY,
    )
    parser.add_argument("--metadata", type=Path)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"Refusing to overwrite {args.output_dir}")
    args.output_dir.mkdir(parents=True)

    f_five = (
        args.project_root
        / "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION/arms/F_five"
    )
    count_rows = count_source_rows(f_five / "provenance/effective_count_sources.tsv")
    write_tsv(
        args.output_dir / "count_source_audit.tsv", count_rows, list(count_rows[0])
    )

    sample_rows: list[dict[str, object]] = []
    cohort_rows: list[dict[str, object]] = []
    current_results = args.project_root / "RNA-seq/Human/Patient_Cohorts/results"
    for cohort in STAR_COHORTS:
        paths = sorted(
            (current_results / cohort / "alignments/star").glob(
                "*/*.ReadsPerGene.out.tab"
            )
        )
        if not paths:
            cohort_rows.append(
                {
                    "dataset": cohort,
                    "audit_source": "STAR_ReadsPerGene",
                    "n_audited": 0,
                    "median_reverse_same_ratio": "",
                    "minimum_ratio": "",
                    "n_below_fourfold": "",
                    "all_samples_pass": False,
                    "status": "fail_missing_star_summaries",
                }
            )
            continue
        ratios: list[float] = []
        for path in paths:
            unstranded, same, reverse = parse_star(path)
            ratio = reverse / same if same else math.inf
            ratios.append(ratio)
            sample_rows.append(
                {
                    "dataset": cohort,
                    "sample_id": sample_id_from_star(path),
                    "audit_source": "STAR_ReadsPerGene",
                    "unstranded_assigned": unstranded,
                    "same_stranded_assigned": same,
                    "reverse_stranded_assigned": reverse,
                    "reverse_same_ratio": ratio,
                    "pass_fourfold": ratio >= MIN_RATIO,
                    "source_path": str(path),
                }
            )
        cohort_rows.append(
            {
                "dataset": cohort,
                "audit_source": "STAR_ReadsPerGene",
                "n_audited": len(ratios),
                "median_reverse_same_ratio": statistics.median(ratios),
                "minimum_ratio": min(ratios),
                "n_below_fourfold": sum(x < MIN_RATIO for x in ratios),
                "all_samples_pass": all(x >= MIN_RATIO for x in ratios),
                "status": "pass"
                if all(x >= MIN_RATIO for x in ratios)
                else "fail_sample_threshold",
            }
        )

    summaries = {
        strand: parse_featurecounts_assigned(
            args.featurecounts_run / f"gene_counts_s{strand}.txt.summary"
        )
        for strand in (0, 1, 2)
    }
    sample_sets = [set(x) for x in summaries.values()]
    if not sample_sets[0] or any(x != sample_sets[0] for x in sample_sets[1:]):
        raise SystemExit(
            "GSE162694 featureCounts sample identities differ across strand modes"
        )
    ratios = []
    for sample in sorted(sample_sets[0]):
        same, reverse = summaries[1][sample], summaries[2][sample]
        ratio = reverse / same if same else math.inf
        ratios.append(ratio)
        sample_rows.append(
            {
                "dataset": "GSE162694",
                "sample_id": sample,
                "audit_source": "featureCounts_s0_s1_s2",
                "unstranded_assigned": summaries[0][sample],
                "same_stranded_assigned": same,
                "reverse_stranded_assigned": reverse,
                "reverse_same_ratio": ratio,
                "pass_fourfold": ratio >= MIN_RATIO,
                "source_path": str(args.featurecounts_run),
            }
        )
    cohort_rows.append(
        {
            "dataset": "GSE162694",
            "audit_source": "featureCounts_s0_s1_s2",
            "n_audited": len(ratios),
            "median_reverse_same_ratio": statistics.median(ratios),
            "minimum_ratio": min(ratios),
            "n_below_fourfold": sum(x < MIN_RATIO for x in ratios),
            "all_samples_pass": all(x >= MIN_RATIO for x in ratios),
            "status": "pass"
            if all(x >= MIN_RATIO for x in ratios)
            else "fail_sample_threshold",
        }
    )

    sample_rows.sort(key=lambda x: (str(x["dataset"]), str(x["sample_id"])))
    cohort_rows.sort(key=lambda x: str(x["dataset"]))
    for row in cohort_rows:
        n_audited = int(row["n_audited"])
        if n_audited:
            n_below = int(row["n_below_fourfold"])
            row["fourfold_fraction"] = (n_audited - n_below) / n_audited
            row["amended_policy_cohort_pass"] = amended_cohort_pass(row)
        else:
            row["fourfold_fraction"] = ""
            row["amended_policy_cohort_pass"] = False
    write_tsv(
        args.output_dir / "strandedness_audit.tsv", sample_rows, list(sample_rows[0])
    )
    write_tsv(
        args.output_dir / "strandedness_cohort_summary.tsv",
        cohort_rows,
        list(cohort_rows[0]),
    )

    count_pass = all(row["status"] == "pass" for row in count_rows)
    stage_rows: list[dict[str, object]] = []
    if args.policy == STRICT_POLICY:
        strand_pass = len(cohort_rows) == len(COHORTS) and all(
            row["status"] == "pass" for row in cohort_rows
        )
        policy_detail = "Every audited sample must have reverse/same >=4."
    else:
        if args.metadata is None:
            raise SystemExit("--metadata is required for the amended policy")
        stage_rows = stage_concentration_audit(sample_rows, args.metadata)
        stage_pass = bool(stage_rows) and all(
            bool(row["stage_concentration_gate_passed"]) for row in stage_rows
        )
        strand_pass = (
            len(cohort_rows) == len(COHORTS)
            and all(amended_cohort_pass(row) for row in cohort_rows)
            and stage_pass
        )
        policy_detail = (
            "Each cohort: median reverse/same >=4, >=90% of samples >=4, and "
            "every sample >1; GSE213621 condition groups must each have <=10% "
            "below fourfold and a maximum-minus-minimum failure fraction <=5 percentage points."
        )
        write_tsv(
            args.output_dir / "stage_concentration_audit.tsv",
            stage_rows,
            list(stage_rows[0]),
        )
    status = "pass" if count_pass and strand_pass else "fail"
    failed = sorted(
        str(row["dataset"]) for row in cohort_rows if row["status"] != "pass"
    )
    gate_rows = [
        {
            "gate": "count_contract",
            "status": "pass" if count_pass else "fail",
            "detail": "All five count sources use featureCounts v2.1.1, GENCODE v49, and -s 2.",
        },
        {
            "gate": "strandedness",
            "status": "pass" if strand_pass else "fail",
            "detail": policy_detail
            if strand_pass
            else f"Policy {args.policy} failed; strict-rule cohorts below fourfold: {','.join(failed)}",
        },
        {
            "gate": "lncrna_source_gate",
            "status": status,
            "detail": "Eligible for lncRNA inference."
            if status == "pass"
            else "Stop lncRNA inference; reopen the bulk counting substrate or amend the gate explicitly.",
        },
    ]
    write_tsv(args.output_dir / "source_gate_status.tsv", gate_rows, list(gate_rows[0]))
    policy_rows = [
        {"field": "policy_id", "value": args.policy},
        {"field": "decision_date", "value": "2026-08-11"},
        {
            "field": "decision_basis",
            "value": "explicit_user_approval_after_strict_gate_failure",
        },
        {"field": "sample_exclusions", "value": "0"},
        {"field": "policy_detail", "value": policy_detail},
    ]
    write_tsv(args.output_dir / "gate_policy.tsv", policy_rows, list(policy_rows[0]))


if __name__ == "__main__":
    main()
