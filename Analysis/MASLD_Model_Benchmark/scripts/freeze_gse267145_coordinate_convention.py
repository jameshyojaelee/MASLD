#!/usr/bin/env python3
"""Freeze the GSE267145 H3K27ac coordinate convention by two prespecified structural tests.

This does NOT set `coordinate_semantics_resolved`; that flag means "a primary source stated it"
and stays False. The freeze is recorded under its own keys and is voided by any later
primary-source statement that contradicts it.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from zipfile import ZipFile

from scripts.audit_gse267145_coordinate_activation import (
    chromhmm_grid_and_offset_test,
    coordinate_claims,
    coordinate_structural_evidence,
    load_chromhmm_breakpoints,
    read_matrix_regions,
    reduce_postcondition_gap_test,
    shared_strings,
    sha256_file,
    worksheet_rows,
    xml_text,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--chromhmm-bed", type=Path, required=True)
    parser.add_argument("--article-xml", type=Path, required=True)
    parser.add_argument("--geo-soft", type=Path, required=True)
    parser.add_argument("--prior-audit", type=Path, required=True,
                        help="coordinate_audit.json of the frozen fail-closed activation")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)

    with ZipFile(args.workbook) as archive:
        strings = shared_strings(archive)
        rows = worksheet_rows(archive, "xl/worksheets/sheet1.xml", strings)
        workbook_text = "\n".join(
            xml_text(archive.read(name)) for name in sorted(archive.namelist())
            if name.endswith(".xml"))
    ordered, _ = read_matrix_regions(args.matrix)
    breakpoints = load_chromhmm_breakpoints(args.chromhmm_bed)
    with gzip.open(args.geo_soft, "rt", encoding="utf-8", errors="strict") as handle:
        geo_text = handle.read()
    direct = {
        "workbook": coordinate_claims(workbook_text),
        "article": coordinate_claims(xml_text(args.article_xml.read_bytes())),
        "geo_subseries_metadata": coordinate_claims(geo_text),
    }
    test_a = chromhmm_grid_and_offset_test(rows, breakpoints)
    test_b = reduce_postcondition_gap_test(ordered)
    evidence = coordinate_structural_evidence(test_a, test_b, direct_claims=direct)
    evidence["direct_coordinate_claims"] = direct
    evidence["inputs_sha256"] = {
        "workbook": sha256_file(args.workbook), "matrix": sha256_file(args.matrix),
        "chromhmm_bed": sha256_file(args.chromhmm_bed),
        "article_xml": sha256_file(args.article_xml), "geo_soft": sha256_file(args.geo_soft),
        "prior_coordinate_audit": sha256_file(args.prior_audit)}
    prior = json.loads(args.prior_audit.read_text(encoding="utf-8"))
    if prior.get("coordinate_semantics_resolved") is not False:
        raise SystemExit("prior audit must be the fail-closed activation; refusing")
    (args.output / "structural_evidence.json").write_text(
        json.dumps(evidence, indent=2), encoding="utf-8")

    frozen = evidence["status"] == "frozen_by_structural_inference"
    convention = evidence["coordinate_convention_structural"]
    freeze = {
        "schema_version": "masld-bench-gse267145-coordinate-convention-v1",
        "coordinate_convention_frozen": convention if frozen else "unresolved",
        "frozen": frozen,
        "resolution_class": evidence["resolution_class"],
        "primary_source_statement_found": evidence["primary_source_statement_found"],
        "coordinate_semantics_resolved_by_primary_source": False,
        "coordinate_semantics_resolved_flag_in_prior_audit": prior["coordinate_semantics_resolved"],
        "bounds_compatibility_used_as_semantic_evidence": False,
        "tool_default_or_file_extension_used_as_semantic_evidence": False,
        "tier1_gene_and_interval_operations_allowed": frozen,
        "tier2_base_resolution_sequence_allowed_under_counted_interval_semantics": frozen,
        "tier2_called_peak_base_resolution_allowed": False,
        "residual_uncertainty": (
            "Test B proves the merge and count steps ran under 1-based arithmetic on the "
            "integers as given. It cannot exclude that the IDR narrowPeak entered R without "
            "the +1 that rtracklayer::import applies, in which case the labels are BED numbers "
            "mis-provenanced but counted 1-based. Invisible at gene scale (7 of 96,460 regions "
            "change block at 25 kb). For sequence, fetch(chrom, start-1, end) returns the "
            "COUNTED interval under either scenario; base-resolution correspondence to the "
            "CALLED peak stays blocked."
        ) if frozen else "not frozen",
        "depositor_query_sent": None,
        "depositor_reply": None,
        "voids_if_primary_statement_contradicts": True,
        "voided_by_primary_source": evidence["voided_by_primary_source"],
        "test_a_summary": {k: test_a[k] for k in (
            "rows_tested", "inclusive_width_multiple_of_200", "half_open_width_multiple_of_200",
            "offset_minus_one_breakpoint_matches", "offset_zero_breakpoint_matches",
            "offset_plus_one_breakpoint_matches", "verdict")},
        "test_b_summary": {k: test_b[k] for k in (
            "same_contig_adjacent_pairs", "min_gap", "n_gap_le_zero", "n_gap_eq_one",
            "n_gap_eq_two", "n_gap_eq_three", "verdict")},
    }
    (args.output / "convention_freeze.json").write_text(
        json.dumps(freeze, indent=2), encoding="utf-8")
    print(json.dumps({"status": evidence["status"], "convention": freeze["coordinate_convention_frozen"],
                      "test_a": freeze["test_a_summary"], "test_b": freeze["test_b_summary"]},
                     indent=1))
    return 0 if frozen else 2


if __name__ == "__main__":
    raise SystemExit(main())
