#!/usr/bin/env python3
"""Audit GSE267145 coordinates and freeze a coordinate-independent task requirement."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
from zipfile import ZipFile

from masld_bench.registry import load_task_spec


NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REGION = re.compile(r"chr(?:[0-9]+|X|Y):[0-9]+-[0-9]+")
CELL_COLUMN = re.compile(r"([A-Z]+)[0-9]+")


class CoordinateActivationError(RuntimeError):
    """Raised when primary evidence or the task-specific requirement differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_artifact(root: Path, expected_sha256: str) -> None:
    observed = sha256_file(root / "ARTIFACTS.json")
    if observed != expected_sha256:
        raise CoordinateActivationError(
            f"input ARTIFACTS SHA-256 differs for {root}: {observed}"
        )


def column_index(reference: str) -> int:
    match = CELL_COLUMN.fullmatch(reference)
    if match is None:
        raise CoordinateActivationError(f"invalid XLSX cell reference: {reference}")
    value = 0
    for character in match.group(1):
        value = value * 26 + ord(character) - ord("A") + 1
    return value - 1


def shared_strings(archive: ZipFile) -> list[str]:
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return ["".join(item.itertext()) for item in root.findall(f"{NS}si")]


def worksheet_rows(
    archive: ZipFile, member: str, strings: list[str]
) -> list[dict[str, str]]:
    root = ET.fromstring(archive.read(member))
    observed: list[list[str]] = []
    for row in root.findall(f".//{NS}sheetData/{NS}row"):
        values: dict[int, str] = {}
        for cell in row.findall(f"{NS}c"):
            reference = cell.attrib.get("r", "")
            value_node = cell.find(f"{NS}v")
            if value_node is None or value_node.text is None:
                value = ""
            elif cell.attrib.get("t") == "s":
                try:
                    value = strings[int(value_node.text)]
                except (IndexError, ValueError) as error:
                    raise CoordinateActivationError(
                        f"invalid shared-string reference: {reference}"
                    ) from error
            else:
                value = value_node.text
            values[column_index(reference)] = value
        width = max(values, default=-1) + 1
        observed.append([values.get(index, "") for index in range(width)])
    if not observed:
        raise CoordinateActivationError(f"worksheet has no rows: {member}")
    header = observed[0]
    if len(header) != len(set(header)):
        raise CoordinateActivationError("worksheet header is duplicated")
    return [
        {name: row[index] if index < len(row) else "" for index, name in enumerate(header)}
        for row in observed[1:]
    ]


def read_matrix_regions(path: Path) -> tuple[list[str], set[str]]:
    ordered: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        if not handle.readline():
            raise CoordinateActivationError("H3K27ac matrix is empty")
        for line_number, line in enumerate(handle, start=2):
            token = line.split(maxsplit=1)[0].strip('"')
            if REGION.fullmatch(token) is None:
                raise CoordinateActivationError(
                    f"H3K27ac region label differs at line {line_number}"
                )
            ordered.append(token)
    if len(ordered) != 96_460 or len(set(ordered)) != len(ordered):
        raise CoordinateActivationError("H3K27ac matrix region axis differs")
    return ordered, set(ordered)


def coordinate_claims(text: str) -> dict[str, list[str]]:
    normalized = " ".join(text.lower().split())
    patterns = {
        "zero_based_half_open": (
            r"\b0[- ]based(?:,?\s+half[- ]open)?\b",
            r"\bzero[- ]based(?:,?\s+half[- ]open)?\b",
            r"\bbed(?:3|6|12)?\b.{0,100}\bhalf[- ]open\b",
        ),
        "one_based_inclusive": (
            r"\b1[- ]based(?:,?\s+inclusive)?\b",
            r"\bone[- ]based(?:,?\s+inclusive)?\b",
            r"\binclusive\s+(?:end(?:point)?|coordinates?)\b",
        ),
    }
    return {
        state: sorted({match.group(0) for pattern in members for match in re.finditer(pattern, normalized)})
        for state, members in patterns.items()
    }


def xml_text(payload: bytes) -> str:
    try:
        root = ET.fromstring(payload)
    except ET.ParseError:
        return payload.decode("utf-8", errors="strict")
    values = list(root.itertext())
    values.extend(value for element in root.iter() for value in element.attrib.values())
    return "\n".join(values)


def audit_supplement(
    workbook: Path, matrix: Path, article_xml: Path, geo_soft: Path
) -> dict[str, object]:
    ordered_matrix, matrix_regions = read_matrix_regions(matrix)
    with ZipFile(workbook) as archive:
        strings = shared_strings(archive)
        rows = worksheet_rows(archive, "xl/worksheets/sheet1.xml", strings)
        workbook_text = "\n".join(
            xml_text(archive.read(name))
            for name in sorted(archive.namelist())
            if name.endswith(".xml")
        )
    required = {"seqnames", "start", "end", "peaks"}
    if not rows or not required <= set(rows[0]):
        raise CoordinateActivationError("Supplementary Table 3 liver sheet differs")
    supplement_regions: list[str] = []
    peak_key_matches = 0
    for row_number, row in enumerate(rows, start=2):
        contig = row["seqnames"]
        try:
            start = int(row["start"])
            end = int(row["end"])
        except ValueError as error:
            raise CoordinateActivationError(
                f"noninteger supplement interval at row {row_number}"
            ) from error
        region = f"{contig}:{start}-{end}"
        if REGION.fullmatch(region) is None or start < 0 or end <= start:
            raise CoordinateActivationError(
                f"invalid supplement interval at row {row_number}"
            )
        if region not in matrix_regions:
            raise CoordinateActivationError(
                f"supplement region is absent from count matrix at row {row_number}"
            )
        supplement_regions.append(region)
        if row["peaks"] == region.replace(":", "_").replace("-", "_"):
            peak_key_matches += 1
    if len(rows) != 14_348 or len(set(supplement_regions)) != len(supplement_regions):
        raise CoordinateActivationError("Supplementary Table 3 liver region axis differs")
    if peak_key_matches != len(rows):
        raise CoordinateActivationError("supplement peaks column does not match coordinates")
    article_text = xml_text(article_xml.read_bytes())
    with gzip.open(geo_soft, "rt", encoding="utf-8", errors="strict") as handle:
        geo_text = handle.read()
    claims = {
        "workbook": coordinate_claims(workbook_text),
        "article": coordinate_claims(article_text),
        "geo_subseries_metadata": coordinate_claims(geo_text),
    }
    direct_zero = any(claims[source]["zero_based_half_open"] for source in claims)
    direct_one = any(claims[source]["one_based_inclusive"] for source in claims)
    if direct_zero and direct_one:
        raise CoordinateActivationError("primary sources contain conflicting coordinate claims")
    resolved = direct_zero ^ direct_one
    convention = (
        "zero_based_half_open"
        if direct_zero
        else "one_based_inclusive"
        if direct_one
        else "unresolved"
    )
    ordered_axis_sha = sha256(
        ("\n".join(ordered_matrix) + "\n").encode("utf-8")
    ).hexdigest()
    supplement_axis_sha = sha256(
        ("\n".join(supplement_regions) + "\n").encode("utf-8")
    ).hexdigest()
    return {
        "schema_version": "masld-bench-gse267145-coordinate-audit-v1",
        "status": "pass_resolved" if resolved else "pass_unresolved_fail_closed",
        "official_supplement": "Supplementary Table 3, liver_h3k27ac_cutrun_NORvsNASH",
        "deposited_h3k27ac_regions": len(ordered_matrix),
        "supplement_liver_regions": len(supplement_regions),
        "supplement_regions_are_exact_matrix_subset": True,
        "supplement_peak_keys_match_seqnames_start_end": True,
        "matrix_region_axis_sha256": ordered_axis_sha,
        "supplement_region_axis_sha256": supplement_axis_sha,
        "direct_coordinate_claims": claims,
        "coordinate_semantics_resolved": resolved,
        "coordinate_convention": convention,
        "bounds_compatibility_used_as_semantic_evidence": False,
        "tool_default_or_file_extension_used_as_semantic_evidence": False,
        "interpretation": (
            "The official differential-region table proves that its seqnames/start/end triples reuse a subset of the deposited count-matrix feature identifiers. Neither the workbook, source article, nor official GEO subseries metadata explicitly states the endpoint convention. File labels, in-bounds coordinates, ChIPseeker/GRanges conventions, and mention of BEDtools are insufficient to distinguish 0-based half-open from 1-based inclusive coordinates."
            if not resolved
            else "A direct primary-source statement resolves the coordinate convention."
        ),
        "sequence_extraction_allowed": resolved,
    }


def read_participant_contract(join_path: Path, qc_path: Path) -> dict[str, object]:
    with join_path.open(encoding="utf-8", newline="") as handle:
        participants = list(csv.DictReader(handle, delimiter="\t"))
    with qc_path.open(encoding="utf-8", newline="") as handle:
        qc = list(csv.DictReader(handle, delimiter="\t"))
    if len(participants) != 99 or len(qc) != 99:
        raise CoordinateActivationError("paired participant or QC axis differs")
    ids = [row["participant_id"] for row in participants]
    if len(set(ids)) != 99 or set(ids) != {row["participant_id"] for row in qc}:
        raise CoordinateActivationError("participant join and QC axes differ")
    if any(row["pairing"] != "same_sample_different_aliquot" for row in participants):
        raise CoordinateActivationError("paired topology differs")
    if any(row["hard_qc_state"] != "observed" for row in qc):
        raise CoordinateActivationError("a paired participant fails hard QC")
    folds = Counter(int(row["outer_fold"]) for row in participants)
    if set(folds) != {0, 1, 2, 3, 4}:
        raise CoordinateActivationError("outer participant folds differ")
    return {
        "participants": 99,
        "biological_unit": "participant",
        "pairing": "same_sample_different_aliquot",
        "outer_fold_counts": {str(key): folds[key] for key in sorted(folds)},
        "hard_qc_passed_participants": 99,
        "cells_or_assay_records_are_replicates": False,
    }


def validate_task_contract(task_path: Path, gate_path: Path) -> dict[str, object]:
    task = load_task_spec(task_path)
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    if task.task_id != "paired_bulk_rna_h3k27ac":
        raise CoordinateActivationError("task_id differs")
    if task.unit_of_inference != "participant":
        raise CoordinateActivationError("task unit is not participant")
    if task.required_pairing_levels[0].value != "same_sample_different_aliquot":
        raise CoordinateActivationError("task pairing differs")
    if task.promotion_gate_config_sha256 != sha256_file(gate_path):
        raise CoordinateActivationError("promotion gate SHA-256 differs")
    if gate.get("champion_eligible") is not False or gate.get("claim_mode") != "development_only":
        raise CoordinateActivationError("single-cohort promotion gate is unsafe")
    text = "\n".join((task.endpoint, *task.admission_gates, task.claim_gate)).lower()
    required = (
        "opaque",
        "sequence extraction",
        "outer training",
        "single-cell",
        "derivative weights",
        "source-owned",
    )
    if any(token not in text for token in required):
        raise CoordinateActivationError("task admission gates omit a required guard")
    return {
        "task_id": task.task_id,
        "task_status": task.status.value,
        "task_spec_sha256": sha256_file(task_path),
        "promotion_gate_id": task.promotion_gate_id,
        "promotion_gate_sha256": sha256_file(gate_path),
        "champion_eligible": False,
        "external_claim_eligible": False,
    }


def activation_contract(
    coordinate: dict[str, object], participants: dict[str, object], task: dict[str, object]
) -> dict[str, object]:
    return {
        "schema_version": "masld-bench-gse267145-task-activation-v1",
        "status": "coordinate_independent_paired_bulk_active_sequence_blocked",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "task_id": task["task_id"],
        "biological_unit": participants["biological_unit"],
        "participants": participants["participants"],
        "pairing": participants["pairing"],
        "task_specific_model_fitting_allowed": True,
        "dataset_wide_or_sequence_activation": False,
        "h3_feature_identity": "opaque_source_region_string",
        "coordinate_semantics_resolved": coordinate["coordinate_semantics_resolved"],
        "sequence_extraction_allowed": coordinate["sequence_extraction_allowed"],
        "allowed_lanes": {
            "rna_conditioned_h3_profile": {
                "query_modalities": ["bulk_rna"],
                "held_h3_access": "evaluator_only",
                "primary_endpoint": "participant_macro_h3_profile_multinomial_deviance_skill",
            },
            "observed_pair_multiview": {
                "query_modalities": ["bulk_rna", "h3k27ac_cutrun"],
                "held_h3_access": "observed_by_task_definition",
                "secondary_endpoints": [
                    "paired_participant_retrieval_mrr",
                    "paired_participant_retrieval_top1",
                ],
                "cannot_support_rna_conditioned_claim": True,
            },
        },
        "eligible_bulk_model_families": [
            "pca_ridge",
            "reduced_rank_regression",
            "pls2",
            "sparse_cca",
            "mofa_plus",
            "generalized_low_rank_multiview",
            "regularized_paired_bulk_autoencoder_screening_only",
        ],
        "mandatory_baselines": [
            "training_mean_h3_profile",
            "pca_ridge",
            "reduced_rank_regression",
            "pls2",
            "sparse_cca",
        ],
        "single_cell_model_families_allowed": False,
        "preprocessing": {
            "rna_measurement": "nonnegative_fractional_gene_expression_estimates",
            "rna_allowed_transforms": ["log1p_library_scaled", "rank_or_quantile_mapper"],
            "rna_integer_count_likelihood_allowed": False,
            "rna_allowed_genes": 42_163,
            "rna_masked_retired_or_absent_genes": 1_122,
            "h3_measurement": "nonnegative_integer_region_counts",
            "h3_profile_target": "count_composition_conditioned_on_observed_evaluator_total",
            "all_learned_transforms_and_feature_selection_fit": "outer_training_participants_only",
            "missing_assay_encoded_as_zero": False,
        },
        "source_outcome_firewall": {
            "fields": [
                "stage",
                "fibrosis",
                "steatosis",
                "ballooning",
                "lobular_inflammation",
                "lobular_necrosis",
                "sex",
                "other_clinical_or_pathology_metadata",
            ],
            "used_only_in_already_frozen_fold_stratification": True,
            "model_input_allowed": False,
            "fitting_target_allowed": False,
            "model_selection_or_error_selector_allowed": False,
            "separate_task_spec_required": True,
        },
        "blocked_actions": [
            "sequence_extraction",
            "coordinate_conversion",
            "interval_overlap",
            "genomic_distance",
            "motif_analysis",
            "coordinate_derived_gene_or_peak_annotation",
            "sequence_or_context_histone_head_training",
            "single_cell_method_ported_by_convenience",
            "source_outcome_supervision_or_selection_under_this_task",
            "task_champion_or_external_transfer_claim",
            "raw_or_processed_source_data_redistribution",
            "derivative_weight_release_without_model_specific_terms_and_legal_review",
        ],
        "rights": {
            "internal_nonclinical_research": "allowed_without_new_DUA",
            "source_data_redistribution": "prohibited_by_project_policy",
            "derivative_weight_release": "requires_model_specific_terms_and_legal_review",
        },
        "coordinate_independent_due_to_unresolved_coordinates": not bool(
            coordinate["coordinate_semantics_resolved"]
        ),
        "champion_eligible": task["champion_eligible"],
        "external_claim_eligible": task["external_claim_eligible"],
        "model_fitted_or_selected_by_this_audit": False,
    }


# ---------------------------------------------------------------------------
# Structural evidence for the coordinate convention.
#
# coordinate_claims() above accepts only a literal primary-source statement and
# is deliberately left untouched. The two tests below are a SEPARATE, sibling
# class of evidence: prespecified, two-sided and falsifiable. They can freeze a
# convention by structural inference, and they record that they did so by
# inference and never as a primary-source statement. A later primary-source
# sentence that contradicts them voids the freeze.
# ---------------------------------------------------------------------------

STRUCTURAL_MIN_OFFSET_MATCHES = 1_000
STRUCTURAL_MIN_WIDTH_RATIO = 100.0


def parse_region(label: str) -> tuple[str, int, int]:
    contig, span = label.split(":", 1)
    start, end = span.split("-", 1)
    return contig, int(start), int(end)


def load_chromhmm_breakpoints(bed_path: Path) -> dict[str, set[int]]:
    """Segment boundaries of a 0-based half-open ChromHMM BED, per contig."""
    breakpoints: dict[str, set[int]] = {}
    with gzip.open(bed_path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("track", "#")):
                continue
            fields = line.rstrip("\n").split("\t")
            contig, start, end = fields[0], int(fields[1]), int(fields[2])
            bucket = breakpoints.setdefault(contig, set())
            bucket.add(start)
            bucket.add(end)
    if not breakpoints:
        raise CoordinateActivationError("ChromHMM BED is empty")
    return breakpoints


def chromhmm_grid_and_offset_test(
    rows: list[dict[str, str]], breakpoints: dict[str, set[int]]
) -> dict[str, object]:
    """Test A. Supplementary Table 3 carries ChromHMM-derived start_position and
    end_position columns. ChromHMM emits 0-based half-open 200-bp bins, so the
    external file's breakpoints are known to be 0-based. Under a 1-based
    inclusive reading, (end - start + 1) is a multiple of 200 and
    start_position - 1 lands on a breakpoint; under a 0-based half-open reading
    the same counts appear at (end - start) and at start_position itself.
    The two hypotheses swap the four counts wholesale, so the test can fail."""
    tested = inclusive_mult = half_open_mult = 0
    hit_minus_one = hit_zero = hit_plus_one = 0
    for row in rows:
        start_text = str(row.get("start_position", "")).strip()
        end_text = str(row.get("end_position", "")).strip()
        if not (start_text.lstrip("-").isdigit() and end_text.lstrip("-").isdigit()):
            continue
        start, end = int(start_text), int(end_text)
        contig = str(row.get("seqnames", ""))
        tested += 1
        if (end - start + 1) % 200 == 0:
            inclusive_mult += 1
        if (end - start) % 200 == 0:
            half_open_mult += 1
        bucket = breakpoints.get(contig, set())
        if (start - 1) in bucket:
            hit_minus_one += 1
        if start in bucket:
            hit_zero += 1
        if (start + 1) in bucket:
            hit_plus_one += 1
    one_based = (
        inclusive_mult >= STRUCTURAL_MIN_WIDTH_RATIO * max(half_open_mult, 1)
        and hit_minus_one >= STRUCTURAL_MIN_OFFSET_MATCHES
        and hit_zero == 0
    )
    zero_based = (
        half_open_mult >= STRUCTURAL_MIN_WIDTH_RATIO * max(inclusive_mult, 1)
        and hit_zero >= STRUCTURAL_MIN_OFFSET_MATCHES
        and hit_minus_one == 0
    )
    if one_based and zero_based:
        raise CoordinateActivationError("Test A returned both conventions at once")
    verdict = (
        "one_based_inclusive" if one_based
        else "zero_based_half_open" if zero_based
        else "indeterminate"
    )
    return {
        "test": "A_chromhmm_grid_and_offset",
        "rows_tested": tested,
        "inclusive_width_multiple_of_200": inclusive_mult,
        "half_open_width_multiple_of_200": half_open_mult,
        "offset_minus_one_breakpoint_matches": hit_minus_one,
        "offset_zero_breakpoint_matches": hit_zero,
        "offset_plus_one_breakpoint_matches": hit_plus_one,
        "min_offset_matches_required": STRUCTURAL_MIN_OFFSET_MATCHES,
        "min_width_ratio_required": STRUCTURAL_MIN_WIDTH_RATIO,
        "what_this_prints_under_zero_based_half_open": (
            "widths become a multiple of 200 under end-start rather than end-start+1, and the "
            "breakpoint hits move from offset -1 to offset 0"
        ),
        "what_this_test_touches": (
            "the supplement's annotation columns, not the deposited region axis; it reaches "
            "the deposited keys only through the inference that one R session wrote both"
        ),
        "verdict": verdict,
    }


def reduce_postcondition_gap_test(regions: list[str]) -> dict[str, object]:
    """Test B. Sort the deposited regions within each contig and look at the gap
    next_start - prev_end between adjacent regions. GenomicRanges::reduce under
    1-based inclusive arithmetic merges overlapping AND abutting ranges, so no
    adjacent pair can have a gap of 0 or 1 and the floor is 2. bedtools merge on
    0-based half-open intervals merges only overlapping intervals, so gap-1 pairs
    survive at roughly the rate of gap-2 pairs. This touches the deposited axis
    directly and needs no external file."""
    by_contig: dict[str, list[tuple[int, int]]] = {}
    for label in regions:
        contig, start, end = parse_region(label)
        by_contig.setdefault(contig, []).append((start, end))
    pairs = 0
    gaps_le_zero = gaps_eq_one = gaps_eq_two = gaps_eq_three = 0
    min_gap: int | None = None
    for spans in by_contig.values():
        spans.sort()
        for (_, prev_end), (next_start, _) in zip(spans, spans[1:]):
            gap = next_start - prev_end
            pairs += 1
            min_gap = gap if min_gap is None else min(min_gap, gap)
            if gap <= 0:
                gaps_le_zero += 1
            elif gap == 1:
                gaps_eq_one += 1
            elif gap == 2:
                gaps_eq_two += 1
            elif gap == 3:
                gaps_eq_three += 1
    one_based = gaps_le_zero == 0 and gaps_eq_one == 0 and gaps_eq_two >= 1
    zero_based = gaps_le_zero == 0 and gaps_eq_one >= 1
    verdict = (
        "one_based_inclusive" if one_based
        else "zero_based_half_open" if zero_based
        else "indeterminate"
    )
    return {
        "test": "B_reduce_postcondition_gap",
        "same_contig_adjacent_pairs": pairs,
        "min_gap": min_gap,
        "n_gap_le_zero": gaps_le_zero,
        "n_gap_eq_one": gaps_eq_one,
        "n_gap_eq_two": gaps_eq_two,
        "n_gap_eq_three": gaps_eq_three,
        "what_this_prints_under_zero_based_half_open": (
            "at least one adjacent pair with gap exactly 1, since bedtools merge does not "
            "merge book-ended 0-based intervals; observing zero such pairs among ~96,000 "
            "when ~7 are expected is p of order 1e-3"
        ),
        "what_this_test_cannot_exclude": (
            "that the narrowPeak BED entered R without the +1 that rtracklayer::import "
            "applies; the merge and count would still be 1-based on the integers as given, "
            "and the labels would be BED numbers mis-provenanced but counted 1-based"
        ),
        "verdict": verdict,
    }


def coordinate_structural_evidence(
    test_a: dict[str, object],
    test_b: dict[str, object],
    direct_claims: dict[str, dict[str, list[str]]] | None = None,
) -> dict[str, object]:
    """Agreement rule. Freeze only if both tests are non-indeterminate and agree,
    and no direct primary-source statement contradicts. A direct statement always
    wins: if one contradicts, the structural freeze is voided."""
    va, vb = str(test_a["verdict"]), str(test_b["verdict"])
    agree = va == vb and va != "indeterminate"
    direct_zero = direct_one = False
    if direct_claims:
        direct_zero = any(v.get("zero_based_half_open") for v in direct_claims.values())
        direct_one = any(v.get("one_based_inclusive") for v in direct_claims.values())
    primary_found = direct_zero or direct_one
    contradicted = agree and (
        (va == "one_based_inclusive" and direct_zero)
        or (va == "zero_based_half_open" and direct_one)
    )
    if agree and not contradicted:
        convention = va
        status = "frozen_by_structural_inference"
    elif contradicted:
        convention = "unresolved"
        status = "voided_by_contradicting_primary_source_statement"
    elif va == "indeterminate" or vb == "indeterminate":
        convention = "unresolved"
        status = "structural_evidence_indeterminate_fail_closed"
    else:
        convention = "unresolved"
        status = "structural_evidence_conflicts_fail_closed"
    return {
        "schema_version": "masld-bench-gse267145-coordinate-structural-evidence-v1",
        "test_a": test_a,
        "test_b": test_b,
        "tests_agree": agree,
        "coordinate_convention_structural": convention,
        "status": status,
        "resolution_class": (
            "structural_inference_two_independent_prespecified_tests"
            if status == "frozen_by_structural_inference" else "none"
        ),
        "primary_source_statement_found": primary_found,
        "coordinate_semantics_resolved_by_primary_source": False,
        "bounds_compatibility_used_as_semantic_evidence": False,
        "tool_default_or_file_extension_used_as_semantic_evidence": False,
        "voids_if_primary_statement_contradicts": True,
        "voided_by_primary_source": contradicted,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--supplement", type=Path, required=True)
    parser.add_argument("--supplement-artifacts-sha256", required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--metadata-artifacts-sha256", required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--matrix-artifacts-sha256", required=True)
    parser.add_argument("--join", type=Path, required=True)
    parser.add_argument("--join-artifacts-sha256", required=True)
    parser.add_argument("--measurement", type=Path, required=True)
    parser.add_argument("--measurement-artifacts-sha256", required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--reference-artifacts-sha256", required=True)
    parser.add_argument("--qc-rights", type=Path, required=True)
    parser.add_argument("--qc-rights-artifacts-sha256", required=True)
    parser.add_argument("--task-spec", type=Path, required=True)
    parser.add_argument("--task-spec-sha256", required=True)
    parser.add_argument("--promotion-gate", type=Path, required=True)
    parser.add_argument("--promotion-gate-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for root, expected in (
        (args.supplement, args.supplement_artifacts_sha256),
        (args.metadata, args.metadata_artifacts_sha256),
        (args.matrix, args.matrix_artifacts_sha256),
        (args.join, args.join_artifacts_sha256),
        (args.measurement, args.measurement_artifacts_sha256),
        (args.reference, args.reference_artifacts_sha256),
        (args.qc_rights, args.qc_rights_artifacts_sha256),
    ):
        verify_artifact(root, expected)
    if sha256_file(args.task_spec) != args.task_spec_sha256:
        raise CoordinateActivationError("task spec SHA-256 differs")
    if sha256_file(args.promotion_gate) != args.promotion_gate_sha256:
        raise CoordinateActivationError("promotion gate SHA-256 differs")
    args.output.mkdir(parents=True, exist_ok=False)
    coordinate = audit_supplement(
        args.supplement / "raw" / "media-3.xlsx",
        args.matrix / "raw" / "GSE267119_H3K27ac.txt.gz",
        args.join / "raw" / "PMC11071482.xml",
        args.metadata / "raw" / "GSE267119_family.soft.gz",
    )
    participants = read_participant_contract(
        args.join / "participant_join.tsv", args.qc_rights / "participant_qc.tsv"
    )
    task = validate_task_contract(args.task_spec, args.promotion_gate)
    contract = activation_contract(coordinate, participants, task)
    outputs = {
        "coordinate_audit.json": coordinate,
        "participant_split_contract.json": participants,
        "task_spec_receipt.json": task,
        "activation_contract.json": contract,
    }
    for filename, payload in outputs.items():
        (args.output / filename).write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(json.dumps(contract, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
