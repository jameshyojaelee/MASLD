#!/usr/bin/env python3
"""Build the GSE105127 CpG crosswalk with the exact legacy-reference requirements."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sqlite3

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from scripts.build_gse105127_cpg_crosswalk import (
    GSE105127CrosswalkError,
    PENDING_CHAIN,
    PENDING_TARGET_CHECK,
    TARGET_DECOMPRESSED_SHA256,
    _run,
    classify_forward,
    classify_reverse,
    collect_collapsed,
    create_database,
    export_crosswalk,
    getfasta_states,
    load_mapping,
    merge_union,
    sha256_file,
)
from scripts.gse105127_exact_reference_contract import (
    REFERENCE_RECEIPT_STATUS,
    admits_exact_reference,
)


def build_crosswalk(
    *,
    plan_root: Path,
    reference_root: Path,
    collapsed_root: Path,
    target_fasta: Path,
    output: Path,
    work: Path,
) -> dict[str, object]:
    if output.exists() or work.exists():
        raise GSE105127CrosswalkError(
            "refusing to overwrite crosswalk output or work"
        )
    verify_frozen_tree(plan_root)
    reference_manifest = verify_frozen_tree(reference_root)
    reference = json.loads(
        (reference_root / "receipt.json").read_text(encoding="utf-8")
    )
    if not admits_exact_reference(reference_manifest, reference):
        raise GSE105127CrosswalkError("exact legacy reference bundle differs")
    if sha256_file(target_fasta) != TARGET_DECOMPRESSED_SHA256:
        raise GSE105127CrosswalkError(
            "decompressed GRCh38.p14 FASTA identity differs"
        )
    for executable in ("bedtools", "CrossMap"):
        if shutil.which(executable) is None:
            raise GSE105127CrosswalkError(
                f"required executable is absent: {executable}"
            )
    paths = collect_collapsed(plan_root, collapsed_root)
    output.mkdir(parents=True)
    work.mkdir(parents=True)
    union = work / "union.tsv"
    merge_union(paths, union)
    database = work / "crosswalk.sqlite"
    source_check = work / "source.check.bed"
    unique_sites = create_database(union, database, source_check)
    source_sequence = work / "source.sequence.tsv"
    getfasta_states(
        database,
        reference_root / "human_g1k_v37.fasta",
        source_check,
        source_sequence,
        False,
    )
    valid_source = work / "source.valid.hg19.bed"
    connection = sqlite3.connect(database)
    with valid_source.open("x", encoding="utf-8") as handle:
        for identifier, chrom, start, end in connection.execute(
            "SELECT id,source_chrom,source_start,source_end FROM sites "
            "WHERE state=? ORDER BY id",
            (PENDING_CHAIN,),
        ):
            handle.write(f"{chrom}\t{start}\t{end}\tcpg_{identifier:09d}\n")
    connection.close()
    forward = work / "forward.hg38.bed"
    forward_unmap = work / "forward.unmap.bed"
    _run(
        [
            "CrossMap",
            "bed",
            "--chromid",
            "l",
            "--unmap-file",
            str(forward_unmap),
            str(reference_root / "hg19ToHg38.over.chain.gz"),
            str(valid_source),
            str(forward),
        ]
    )
    load_mapping(database, forward, "forward")
    target_bed = work / "target.check.bed"
    classify_forward(database, target_bed)
    target_sequence = work / "target.sequence.tsv"
    getfasta_states(database, target_fasta, target_bed, target_sequence, True)
    reverse_input = work / "reverse.input.hg38.bed"
    connection = sqlite3.connect(database)
    with reverse_input.open("x", encoding="utf-8") as handle:
        for identifier, chrom, start, end in connection.execute(
            "SELECT id,target_chrom,target_start,target_end FROM sites "
            "WHERE state=? ORDER BY id",
            (PENDING_TARGET_CHECK,),
        ):
            handle.write(f"{chrom}\t{start}\t{end}\tcpg_{identifier:09d}\n")
    connection.close()
    reverse = work / "reverse.hg19.bed"
    reverse_unmap = work / "reverse.unmap.bed"
    _run(
        [
            "CrossMap",
            "bed",
            "--chromid",
            "l",
            "--unmap-file",
            str(reverse_unmap),
            str(reference_root / "hg38ToHg19.over.chain.gz"),
            str(reverse_input),
            str(reverse),
        ]
    )
    load_mapping(database, reverse, "reverse")
    classify_reverse(database)
    crosswalk = output / "cpg_crosswalk.tsv.gz"
    counts = export_crosswalk(database, crosswalk)
    receipt = {
        "schema_version": "masld-bench-gse105127-cpg-crosswalk-v1",
        "status": "passed_failure_aware_roundtrip",
        "source_build": "1000_Genomes_GRCh37",
        "source_reference_status": REFERENCE_RECEIPT_STATUS,
        "target_build": "GRCh38.p14",
        "unique_source_cpg_intervals": unique_sites,
        "mapping_states": dict(sorted(counts.items())),
        "accepted_state": "mapped_unique_cpg",
        "all_source_intervals_retained": sum(counts.values()) == unique_sites,
        "coordinate_convention": "zero_based_half_open_two_base_CpG",
        "source_and_target_CG_verified": True,
        "roundtrip_required": True,
        "missing_or_failed_as_zero": False,
        "labels_accessed": False,
        "fit_or_score_performed": False,
        "collapsed_sample_files": 57,
        "crosswalk_sha256": sha256_file(crosswalk),
    }
    if receipt["all_source_intervals_retained"] is not True:
        raise GSE105127CrosswalkError(
            "crosswalk silently dropped a source CpG"
        )
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(
        output,
        {
            "artifact_class": "gse105127_failure_aware_cpg_crosswalk",
            "status": "passed",
            "source_reference_status": REFERENCE_RECEIPT_STATUS,
        },
    )
    verify_frozen_tree(output)
    shutil.rmtree(work)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--collapsed-root", type=Path, required=True)
    parser.add_argument("--target-fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    arguments = parser.parse_args()
    print(json.dumps(build_crosswalk(**vars(arguments)), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
