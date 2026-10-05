#!/usr/bin/env python3
"""Recover qualified BAM inputs and outcome-free pilots without reading BAMs."""
import argparse
import csv
import hashlib
import json
import os
import re
import shlex
from collections import Counter
from pathlib import Path


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path, rows):
    with path.open("x") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sample_id(path):
    matches = set(re.findall(r"(?:SRR|ERR|DRR)\d+", path))
    if len(matches) != 1:
        raise ValueError("BAM column lacks one unambiguous run ID: " + path)
    return next(iter(matches))


def count_header(path):
    with path.open() as handle:
        command_header, columns_header = handle.readline(), handle.readline()
    columns = columns_header.rstrip("\n").split("\t")
    if columns[:6] != ["Geneid", "Chr", "Start", "End", "Strand", "Length"]:
        raise ValueError("Count header schema differs: " + str(path))
    bam_columns = columns[6:]
    if not bam_columns or any(not Path(p).is_absolute() or not p.endswith(".bam") for p in bam_columns):
        raise ValueError("Count columns are not absolute BAM paths")
    command = shlex.split(command_header.split("Command:", 1)[1].split("; chunks=", 1)[0])
    command_bams = [token for token in command if token.endswith(".bam")]
    if command_bams and command_bams != bam_columns:
        raise ValueError("Command BAM order differs from count columns")
    if "-a" not in command or "-s" not in command:
        raise ValueError("Native annotation or explicit strandedness missing")
    gtf = command[command.index("-a") + 1]
    strand = int(command[command.index("-s") + 1])
    if strand not in {0, 1, 2} or not Path(gtf).is_file():
        raise ValueError("Native annotation/strandedness is invalid")
    paired = "-p" in command
    if paired != ("--countReadPairs" in command):
        raise ValueError("Paired layout and fragment counting disagree")
    flags = [token for token in command if token not in command_bams]
    return {
        "bams": {sample_id(p): p for p in bam_columns},
        "n_bam_columns": len(bam_columns),
        "layout": "paired" if paired else "single", "strandedness": strand,
        "native_gtf": gtf, "original_counting_flags": json.dumps(flags),
        "header_sha256": hashlib.sha256((command_header + columns_header).encode()).hexdigest(),
        "bam_path_evidence": "count_columns_and_command" if command_bams else "merged_count_columns",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rosters", type=Path, required=True)
    parser.add_argument("--count-sources", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    union_path = args.rosters / "count_union_776_sample_ids.tsv"
    union = read_tsv(union_path)
    if len(union) != 776 or len({r["sample_id"] for r in union}) != 776:
        raise ValueError("Qualified union is not 776 distinct libraries")
    primary_path = args.rosters / "primary_bulk_metadata.tsv"
    strict_path = args.rosters / "strict_depth_bulk_metadata.tsv"
    primary, strict = read_tsv(primary_path), read_tsv(strict_path)
    primary_ids = {r["sample_id"] for r in primary}
    strict_ids = {r["sample_id"] for r in strict}
    if primary_ids | strict_ids != {r["sample_id"] for r in union}:
        raise ValueError("Union differs from primary plus strict-depth libraries")
    metadata = {r["sample_id"]: r for r in strict + primary}
    sources = {r["dataset"]: r for r in read_tsv(args.count_sources)}
    cohorts = sorted({r["dataset"] for r in union})
    if len(cohorts) != 5:
        raise ValueError("Exactly five cohorts expected")
    headers = {cohort: count_header(Path(sources[cohort]["count_path"])) for cohort in cohorts}
    for cohort, header in headers.items():
        if len(header["bams"]) != header["n_bam_columns"]:
            raise ValueError("Repeated sample columns: " + cohort)
    rows = []
    for record in sorted(union, key=lambda r: (r["dataset"], r["sample_id"])):
        sid, cohort = record["sample_id"], record["dataset"]
        header = headers[cohort]
        meta = metadata[sid]
        if meta["dataset"] != cohort:
            raise ValueError("Roster dataset disagreement: " + sid)
        bam = Path(header["bams"][sid])
        info = bam.stat()
        if not bam.is_file() or info.st_size <= 0 or not os.access(bam, os.R_OK):
            raise ValueError("BAM missing, empty or unreadable: " + str(bam))
        rows.append({"sample_id": sid, "dataset": cohort, "bam_path": str(bam),
            "layout": header["layout"], "strandedness": header["strandedness"],
            "native_gtf": header["native_gtf"], "bam_size_bytes": info.st_size,
            "bam_mtime_ns": info.st_mtime_ns, "original_counting_flags": header["original_counting_flags"],
            "bam_path_evidence": header["bam_path_evidence"], "count_header_sha256": header["header_sha256"],
            "source_count_path": sources[cohort]["count_path"],
            "source_count_recorded_sha256": sources[cohort]["sha256"],
            "full_called": int(meta["full_called"]), "analysis_unit_id": meta["analysis_unit_id"],
            "unit_evidence": meta["unit_evidence"], "primary_selected": sid in primary_ids,
            "strict_depth_selected": sid in strict_ids})
    if len({r["bam_path"] for r in rows}) != 776:
        raise ValueError("Distinct qualified libraries point to repeated BAMs")
    pilots = []
    for cohort in cohorts:
        eligible = sorted((r for r in rows if r["dataset"] == cohort and r["primary_selected"]),
                          key=lambda r: (r["full_called"], r["sample_id"]))
        pilot = dict(eligible[(len(eligible) - 1) // 2])
        pilot["pilot_selection"] = "primary_roster_lower_median_full_called_tie_SRR"
        pilot["primary_cohort_n"] = len(eligible)
        pilot["coverage_order_one_based"] = (len(eligible) - 1) // 2 + 1
        pilots.append(pilot)
    write_tsv(args.out / "library_manifest_776.tsv", rows)
    write_tsv(args.out / "pilot_libraries.tsv", pilots)
    summary = {"status": "headers_and_776_BAM_stats_checked_no_BAM_scans", "libraries": len(rows),
        "primary_libraries": len(primary_ids), "strict_libraries": len(strict_ids),
        "cohort_counts": dict(Counter(r["dataset"] for r in rows)),
        "total_bam_bytes": sum(r["bam_size_bytes"] for r in rows),
        "total_bam_TB_decimal": sum(r["bam_size_bytes"] for r in rows) / 1e12,
        "total_bam_TiB": sum(r["bam_size_bytes"] for r in rows) / 2**40,
        "pilots": [{k: r[k] for k in ("sample_id", "dataset", "full_called", "bam_size_bytes", "layout", "strandedness")} for r in pilots],
        "selection_uses_disease_effects": False, "BAM_content_reads": 0,
        "input_sha256": {str(p): sha(p) for p in (union_path, primary_path, strict_path, args.count_sources)},
        "count_file_full_hashes_rechecked": False,
        "count_file_evidence": "first two headers read; full hashes retained from effective-count-source ledger",
        "script_sha256": sha(Path(__file__))}
    with (args.out / "preparation.json").open("x") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
