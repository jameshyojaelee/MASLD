#!/usr/bin/env python3
"""Source-gate paired human cohorts without reading expression outcomes."""

from __future__ import annotations

import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

from bridge_common import CANDIDATE_ROOT, PROJECT_ROOT, atomic_write_text, require_validated_seal, write_tsv


PLAN41_ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09"


def matrix_metadata(path: Path) -> tuple[dict[str, list[str]], dict[str, list[str]], dict[str, str]]:
    sample: dict[str, list[str]] = {}
    characteristics: list[list[str]] = []
    series: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Series_"):
                row = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                series.setdefault(row[0], row[1] if len(row) > 1 else "")
            elif line.startswith("!Sample_"):
                row = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                if row[0] == "!Sample_characteristics_ch1":
                    characteristics.append(row[1:])
                else:
                    sample.setdefault(row[0], row[1:])
    fields: dict[str, list[str]] = {}
    for values in characteristics:
        keys = {value.split(":", 1)[0].strip().lower() for value in values if ":" in value}
        if len(keys) == 1:
            fields[next(iter(keys))] = [value.split(":", 1)[1].strip() for value in values]
    return sample, fields, series


def gate_83452() -> tuple[list[dict[str, object]], dict[str, object]]:
    path = CANDIDATE_ROOT / "sources/GSE83452/GSE83452_series_matrix.txt.gz"
    if not path.is_file():
        return [], gate("GSE83452", "skipped_source_unavailable", False, 0, 0, "matrix missing")
    sample, fields, series = matrix_metadata(path)
    titles = sample.get("!Sample_title", [])
    accessions = sample.get("!Sample_geo_accession", [])
    n = len(accessions)
    required = {"sample name", "liver status", "type of intervention", "time"}
    if n != 231 or not required.issubset(fields):
        return [], gate("GSE83452", "rejected_qc", False, n, 0, f"metadata fields={sorted(fields)}")
    baseline_by_number: dict[str, int] = {}
    for index, title in enumerate(titles):
        match = re.fullmatch(r"liver biopsy (\d+)", title)
        if match and fields["time"][index] == "baseline":
            baseline_by_number[match.group(1)] = index
    pairs: list[dict[str, object]] = []
    for follow_index, title in enumerate(titles):
        if fields["time"][follow_index] != "follow-up":
            continue
        match = re.fullmatch(r"liver biopsy (\d+) \((\d+)\)", title)
        if not match or match.group(2) not in baseline_by_number:
            continue
        baseline_index = baseline_by_number[match.group(2)]
        intervention = fields["type of intervention"][follow_index]
        baseline_status = fields["liver status"][baseline_index]
        followup_status = fields["liver status"][follow_index]
        if fields["type of intervention"][baseline_index] != intervention:
            raise RuntimeError(f"GSE83452 intervention mismatch for {title}")
        primary = intervention == "Diet" and baseline_status == "NASH" and followup_status in {"NASH", "no NASH"}
        direction = intervention == "BS" and baseline_status == "NASH" and followup_status == "no NASH"
        classification = "resolver" if followup_status == "no NASH" else "persistent_nash" if followup_status == "NASH" else "undefined"
        pairs.append(pair_row(
            "GSE83452", f"GSE83452_{match.group(2)}", accessions[baseline_index], accessions[follow_index],
            intervention, baseline_status, followup_status, "", "", classification,
            primary, direction, False,
            "GEO parenthetical follow-up title linked to exact baseline title; deposited intervention/time/liver-status fields",
        ))
    primary_counts = Counter(row["response_class"] for row in pairs if row["primary_eligible"] == "true")
    passed = (
        len(pairs) == 60
        and primary_counts == Counter({"persistent_nash": 16, "resolver": 6})
        and sum(row["direction_check_eligible"] == "true" for row in pairs) == 14
        and "152 patients at baseline" in series.get("!Series_overall_design", "")
    )
    detail = f"231 arrays; 60 title-linked pairs; Diet baseline-NASH resolver/persistent=6/16; BS resolver direction-check=14"
    return pairs, gate("GSE83452", "pass" if passed else "rejected_qc", passed, n, len(pairs), detail)


def gate_48452() -> tuple[list[dict[str, object]], dict[str, object]]:
    path = CANDIDATE_ROOT / "sources/GSE48452/GSE48452_series_matrix.txt.gz"
    if not path.is_file():
        return [], gate("GSE48452", "skipped_source_unavailable", False, 0, 0, "matrix missing")
    sample, fields, _ = matrix_metadata(path)
    titles = sample.get("!Sample_title", [])
    accessions = sample.get("!Sample_geo_accession", [])
    required = {"other_id", "bariatric surgery", "nas"}
    if len(accessions) != 73 or not required.issubset(fields):
        return [], gate("GSE48452", "rejected_qc", False, len(accessions), 0, f"metadata fields={sorted(fields)}")
    file_ids = [f"{title.rsplit(', ', 1)[-1]}.CEL" for title in titles]
    by_file = {value: index for index, value in enumerate(file_ids)}
    pairs: list[dict[str, object]] = []
    seen: set[tuple[int, int]] = set()
    for index, other in enumerate(fields["other_id"]):
        if other == "NA" or other not in by_file:
            continue
        partner = by_file[other]
        key = tuple(sorted((index, partner)))
        if key in seen:
            continue
        seen.add(key)
        symmetric = fields["other_id"][partner] == file_ids[index]
        timing = {fields["bariatric surgery"][index], fields["bariatric surgery"][partner]}
        explicit = timing == {"before surgery", "after surgery"}
        before = index if fields["bariatric surgery"][index] == "before surgery" else partner if fields["bariatric surgery"][partner] == "before surgery" else None
        after = index if fields["bariatric surgery"][index] == "after surgery" else partner if fields["bariatric surgery"][partner] == "after surgery" else None
        ambiguous_baseline = (not symmetric) and after is not None and before is None and fields["bariatric surgery"][partner if after == index else index] == "NA"
        if ambiguous_baseline:
            before = partner if after == index else index
        if before is None or after is None:
            continue
        primary = symmetric and explicit
        sensitivity = ambiguous_baseline
        delta_nas = float(fields["nas"][after]) - float(fields["nas"][before])
        pairs.append(pair_row(
            "GSE48452", f"GSE48452_{file_ids[before].removesuffix('.CEL')}",
            accessions[before], accessions[after], "bariatric_surgery", "", "",
            fields["nas"][before], fields["nas"][after], "continuous_delta_nas",
            primary, False, sensitivity,
            "GEO symmetric other_id and explicit before/after fields" if primary else "one-way other_id from an after-surgery sample to a counterpart with ambiguous deposited timing; sensitivity only",
            delta_nas,
        ))
    n_primary = sum(row["primary_eligible"] == "true" for row in pairs)
    n_sensitivity = sum(row["sensitivity_only"] == "true" for row in pairs)
    passed = n_primary == 13 and n_primary >= 10 and n_sensitivity == 3
    detail = f"73 arrays; 13 symmetric explicit pairs; 3 one-way links with ambiguous counterpart timing retained as sensitivity only"
    return pairs, gate("GSE48452", "pass" if passed else "rejected_qc", passed, 73, n_primary, detail)


def gate_106737() -> tuple[list[dict[str, object]], dict[str, object]]:
    manifest = PLAN41_ROOT / "sample_manifest.tsv"
    if not manifest.is_file():
        return [], gate("GSE106737", "skipped_source_unavailable", False, 0, 0, "terminal Plan 41 manifest missing")
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["dataset_id"] == "GSE106737" and row["include_in_inference"] == "true"]
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["biological_unit_id"], []).append(row)
    pairs: list[dict[str, object]] = []
    for participant, records in sorted(grouped.items()):
        by_time = {row["timepoint"]: row for row in records}
        if set(by_time) != {"baseline", "followup"}:
            continue
        condition = records[0]["condition"]
        pairs.append(pair_row(
            "GSE106737", f"GSE106737_{participant}", by_time["baseline"]["sample_id"],
            by_time["followup"]["sample_id"], "RYGB" if condition == "RYGB_responder" else "lifestyle",
            "NASH" if condition != "lifestyle_nonresponder" else "source_group3",
            "histologic_responder" if "responder" in condition and "nonresponder" not in condition else "nonresponder",
            "", "", condition, condition.startswith("lifestyle_"), condition == "RYGB_responder", False,
            "terminal validated Plan 41 participant pairing",
        ))
    counts = Counter(row["response_class"] for row in pairs)
    passed = counts == Counter({"RYGB_responder": 21, "lifestyle_responder": 10, "lifestyle_nonresponder": 10})
    return pairs, gate("GSE106737", "pass_reuse_discovery" if passed else "rejected_qc", passed, 111, len(pairs), f"validated Plan 41 pairs={dict(counts)}")


def pair_row(dataset: str, pair_id: str, baseline: str, followup: str, intervention: str,
             baseline_status: str, followup_status: str, baseline_nas: str, followup_nas: str,
             response: str, primary: bool, direction: bool, sensitivity: bool, basis: str,
             delta_nas: float | str = "") -> dict[str, object]:
    return {
        "dataset_id": dataset, "pair_id": pair_id, "participant_id": pair_id,
        "baseline_sample_id": baseline, "followup_sample_id": followup,
        "intervention": intervention, "baseline_histology": baseline_status,
        "followup_histology": followup_status, "baseline_nas": baseline_nas,
        "followup_nas": followup_nas, "delta_nas": delta_nas,
        "response_class": response, "primary_eligible": str(primary).lower(),
        "direction_check_eligible": str(direction).lower(),
        "sensitivity_only": str(sensitivity).lower(), "source_basis": basis,
    }


def gate(dataset: str, status: str, authorized: bool, samples: int, complete: int, detail: str) -> dict[str, object]:
    return {
        "dataset_id": dataset, "gate_id": "PAIRED_HUMAN_SOURCE",
        "status": status, "inference_authorized": str(authorized).lower(),
        "n_source_samples": samples, "n_primary_complete_units": complete,
        "detail": detail,
    }


def main() -> None:
    require_validated_seal()
    all_pairs: list[dict[str, object]] = []
    gates: list[dict[str, object]] = []
    for function in (gate_106737, gate_83452, gate_48452):
        pairs, gate_row = function()
        all_pairs.extend(pairs); gates.append(gate_row)
    pair_fields = list(all_pairs[0])
    suffix = ""
    if (CANDIDATE_ROOT / "paired_sample_manifest.tsv").exists():
        version = 2
        while (CANDIDATE_ROOT / f"paired_sample_manifest_v{version}.tsv").exists():
            version += 1
        suffix = f"_v{version}"
    pair_path = CANDIDATE_ROOT / f"paired_sample_manifest{suffix}.tsv"
    gate_path = CANDIDATE_ROOT / f"source_gates/human_pair_gates{suffix}.tsv"
    write_tsv(pair_path, all_pairs, pair_fields)
    write_tsv(gate_path, gates, list(gates[0]))
    fingerprint_path = CANDIDATE_ROOT / "source_gates/sample_fingerprint_matches.tsv"
    duplicate_count = 0
    if fingerprint_path.is_file():
        with fingerprint_path.open(newline="", encoding="utf-8") as handle:
            duplicate_count = sum(row["duplicate_threshold_pass"] == "TRUE" for row in csv.DictReader(handle, delimiter="\t"))
    overlap_status = "confirmed_overlap" if duplicate_count > 0 else "unresolved_probable_overlap"
    overlap = [
        {
            "cohort_a": "GSE106737", "cohort_b": "GSE83452",
            "publication_a": "GEO series without linked PubMed ID; Inserm/Antwerp investigators",
            "publication_b": "PMID:28679947; Antwerp University Hospital",
            "investigator_overlap": "Francque;Lefebvre;Staels",
            "sample_namespace_overlap": "not_authoritatively_resolved",
            "expression_fingerprint_status": f"{duplicate_count}_reciprocal_matches_at_r_ge_0.9999" if duplicate_count else "pending",
            "overlap_verdict": overlap_status,
            "meta_analysis_contribution": "collapse_to_one_antwerp_inserm_contribution",
            "reason": "shared investigators/recruitment context and compatible intervention counts; independence is not publicly established",
        },
        {
            "cohort_a": "GSE106737/GSE83452", "cohort_b": "GSE48452",
            "publication_a": "Antwerp/Inserm intervention cohorts",
            "publication_b": "PMID:23931760; Universitaetsklinikum Regensburg",
            "investigator_overlap": "none_in_deposited_series",
            "sample_namespace_overlap": "none",
            "expression_fingerprint_status": "not_required_different_platform_and_namespace",
            "overlap_verdict": "independent_by_public_provenance",
            "meta_analysis_contribution": "separate_contribution",
            "reason": "distinct publication, recruitment institution, investigators, and sample namespace",
        },
    ]
    overlap_path = CANDIDATE_ROOT / f"source_gates/cohort_overlap_audit{suffix}.tsv"
    write_tsv(overlap_path, overlap, list(overlap[0]))
    effective = [
        {"dataset_id": "GSE106737", "effective_cohort_id": "ANTWERP_INSERM_SHARED", "independent_meta_contribution": "false", "reason": overlap_status},
        {"dataset_id": "GSE83452", "effective_cohort_id": "ANTWERP_INSERM_SHARED", "independent_meta_contribution": "true", "reason": overlap_status},
        {"dataset_id": "GSE48452", "effective_cohort_id": "REGENSBURG_PAIRED", "independent_meta_contribution": "true", "reason": "independent_by_public_provenance"},
    ]
    effective_path = CANDIDATE_ROOT / f"source_gates/effective_cohort_registry{suffix}.tsv"
    write_tsv(effective_path, effective, list(effective[0]))
    active = {
        "active_pair_manifest": str(pair_path.relative_to(CANDIDATE_ROOT)),
        "active_human_gate": str(gate_path.relative_to(CANDIDATE_ROOT)),
        "active_overlap_audit": str(overlap_path.relative_to(CANDIDATE_ROOT)),
        "active_effective_cohort_registry": str(effective_path.relative_to(CANDIDATE_ROOT)),
        "correction_reason": "GSE48452 one-way links have ambiguous counterpart timing and are sensitivity-only; primary 13 symmetric pairs unchanged" if suffix else "",
        "n_effective_independent_cohorts": 2,
        "human_three_independent_cohort_gate_possible": False,
    }
    atomic_write_text(CANDIDATE_ROOT / "source_gates/HUMAN_GATE_ACTIVE.json", json.dumps(active, indent=2, sort_keys=True) + "\n")
    print("HUMAN_PAIR_GATES_COMPLETE")


if __name__ == "__main__":
    main()
