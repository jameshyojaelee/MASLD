#!/usr/bin/env python3
"""Build a failure-aware round-trip GRCh37-to-GRCh38.p14 CpG crosswalk."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re
import shlex
import shutil
import sqlite3
import subprocess

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


TARGET_DECOMPRESSED_SHA256 = "3ed0c28ded22eac00112e47331c3e146f5c0a50b9dbe2d15dac818ce2a8103df"
SAFE_PATH = re.compile(r"^[A-Za-z0-9_./-]+$")
PRIMARY_TARGET = {f"chr{value}" for value in range(1, 23)} | {"chrX", "chrY", "chrM"}
SOURCE_ALIAS = {f"chr{value}": str(value) for value in range(1, 23)} | {
    "chrX": "X", "chrY": "Y", "chrM": "MT",
}
FINAL_STATES = {
    "source_contig_unrecognized", "source_not_cpg", "strand_pair_inconsistent",
    "unmapped_chain", "mapped_nonunique", "mapped_length_changed",
    "mapped_nonprimary", "target_not_cpg", "roundtrip_failed", "mapped_unique_cpg",
}
PENDING_CHAIN = "pending_chain"
PENDING_TARGET_CHECK = "pending_target_check"


class GSE105127CrosswalkError(RuntimeError):
    """Raised when a failure-aware CpG crosswalk invariant differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _run(command: list[str], *, stdout: Path | None = None) -> None:
    if stdout is None:
        result = subprocess.run(command, check=False)
    else:
        with stdout.open("xb") as handle:
            result = subprocess.run(command, stdout=handle, check=False)
    if result.returncode != 0:
        raise GSE105127CrosswalkError(f"command failed: {command}")


def collect_collapsed(plan_root: Path, collapsed_root: Path) -> list[Path]:
    with (plan_root / "rrbs_rows.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 57:
        raise GSE105127CrosswalkError("RRBS plan row census differs")
    paths = []
    for row in rows:
        source = collapsed_root / "participants" / row["row_id"]
        verify_frozen_tree(source)
        receipt = json.loads((source / "receipt.json").read_text(encoding="utf-8"))
        path = source / "cpg_counts.hg19.tsv.gz"
        if (
            receipt.get("status") != "passed"
            or receipt.get("row_id") != row["row_id"]
            or receipt.get("labels_accessed") is not False
            or receipt.get("fit_or_score_performed") is not False
            or not path.is_file()
        ):
            raise GSE105127CrosswalkError("collapsed RRBS member differs")
        paths.append(path.resolve())
    return paths


def merge_union(paths: list[Path], output: Path) -> None:
    if any(SAFE_PATH.fullmatch(str(path)) is None for path in paths):
        raise GSE105127CrosswalkError("collapsed path is unsafe for merge shell")
    streams = " ".join(
        f"<(gzip -dc {shlex.quote(str(path))} | tail -n +2 | cut -f1-4)"
        for path in paths
    )
    command = (
        "set -euo pipefail; export LC_ALL=C; "
        "sort -m -u -t $'\\t' -k1,1n -k3,3n -k4,4n " + streams
    )
    with output.open("xb") as handle:
        result = subprocess.run(["bash", "-c", command], stdout=handle, check=False)
    if result.returncode != 0 or output.stat().st_size == 0:
        raise GSE105127CrosswalkError("sorted RRBS CpG union failed")


def create_database(union: Path, database: Path, source_check: Path) -> int:
    connection = sqlite3.connect(database)
    connection.executescript(
        "PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA temp_store=FILE;"
        "CREATE TABLE sites (id INTEGER PRIMARY KEY, source_chrom TEXT NOT NULL, source_start INTEGER NOT NULL, source_end INTEGER NOT NULL, state TEXT NOT NULL, target_chrom TEXT, target_start INTEGER, target_end INTEGER);"
        "CREATE TABLE forward (id INTEGER PRIMARY KEY, mappings INTEGER NOT NULL, chrom TEXT, start INTEGER, end INTEGER);"
        "CREATE TABLE reverse (id INTEGER PRIMARY KEY, mappings INTEGER NOT NULL, chrom TEXT, start INTEGER, end INTEGER);"
    )
    count = 0
    previous = (-1, -1)
    batch = []
    with union.open(encoding="utf-8") as handle, source_check.open("x", encoding="utf-8") as check:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise GSE105127CrosswalkError("union BED width differs")
            rank, chrom, start, end = int(fields[0]), fields[1], int(fields[2]), int(fields[3])
            if (rank, start) <= previous or end != start + 2:
                raise GSE105127CrosswalkError("union BED order or interval length differs")
            previous = (rank, start)
            count += 1
            identifier = f"cpg_{count:09d}"
            state = PENDING_CHAIN if chrom in SOURCE_ALIAS else "source_contig_unrecognized"
            batch.append((count, chrom, start, end, state))
            if chrom in SOURCE_ALIAS:
                check.write(f"{SOURCE_ALIAS[chrom]}\t{start}\t{end}\t{identifier}\n")
            if len(batch) == 100_000:
                connection.executemany("INSERT INTO sites VALUES (?,?,?,?,?,NULL,NULL,NULL)", batch)
                connection.commit(); batch.clear()
    if batch:
        connection.executemany("INSERT INTO sites VALUES (?,?,?,?,?,NULL,NULL,NULL)", batch)
        connection.commit()
    connection.close()
    if count == 0:
        raise GSE105127CrosswalkError("union CpG census is empty")
    return count


def getfasta_states(database: Path, fasta: Path, bed: Path, output: Path, target: bool) -> None:
    _run(["bedtools", "getfasta", "-fi", str(fasta), "-bed", str(bed), "-nameOnly", "-tab"], stdout=output)
    connection = sqlite3.connect(database)
    expected = connection.execute(
        "SELECT COUNT(*) FROM sites WHERE state=? AND target_chrom IS NOT NULL"
        if target
        else "SELECT COUNT(*) FROM sites WHERE state=?",
        (PENDING_TARGET_CHECK if target else PENDING_CHAIN,),
    ).fetchone()[0]
    updates = []
    observed = 0
    with output.open(encoding="utf-8") as handle:
        for line in handle:
            observed += 1
            name, sequence = line.rstrip("\n").split("\t", 1)
            match = re.match(r"cpg_([0-9]+)", name)
            if match is None:
                raise GSE105127CrosswalkError("bedtools CpG identifier differs")
            identifier = int(match.group(1))
            if sequence.upper() != "CG":
                updates.append(("target_not_cpg" if target else "source_not_cpg", identifier))
            if len(updates) == 100_000:
                connection.executemany("UPDATE sites SET state=? WHERE id=?", updates)
                connection.commit(); updates.clear()
    if updates:
        connection.executemany("UPDATE sites SET state=? WHERE id=?", updates)
        connection.commit()
    if observed != expected:
        raise GSE105127CrosswalkError("bedtools sequence-output census differs")
    connection.close()


def load_mapping(database: Path, path: Path, table: str) -> None:
    connection = sqlite3.connect(database)
    batch = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise GSE105127CrosswalkError("CrossMap BED width differs")
            match = re.fullmatch(r"cpg_([0-9]+)", fields[3])
            if match is None:
                raise GSE105127CrosswalkError("CrossMap CpG identifier differs")
            batch.append((int(match.group(1)), fields[0], int(fields[1]), int(fields[2])))
            if len(batch) == 100_000:
                _upsert_mapping(connection, table, batch); batch.clear()
    if batch:
        _upsert_mapping(connection, table, batch)
    connection.commit(); connection.close()


def _upsert_mapping(connection: sqlite3.Connection, table: str, rows: list[tuple[int, str, int, int]]) -> None:
    if table not in {"forward", "reverse"}:
        raise GSE105127CrosswalkError("invalid mapping table")
    connection.executemany(
        f"INSERT INTO {table}(id,mappings,chrom,start,end) VALUES (?,1,?,?,?) "
        "ON CONFLICT(id) DO UPDATE SET mappings=mappings+1",
        rows,
    )


def classify_forward(database: Path, target_bed: Path) -> None:
    connection = sqlite3.connect(database)
    connection.execute(
        "UPDATE sites SET state='unmapped_chain' WHERE state=? AND id NOT IN (SELECT id FROM forward)",
        (PENDING_CHAIN,),
    )
    connection.execute(
        "UPDATE sites SET state='mapped_nonunique' WHERE state=? AND id IN (SELECT id FROM forward WHERE mappings<>1)",
        (PENDING_CHAIN,),
    )
    connection.execute(
        "UPDATE sites SET state='mapped_length_changed' WHERE state=? AND id IN (SELECT id FROM forward WHERE mappings=1 AND end-start<>2)",
        (PENDING_CHAIN,),
    )
    placeholders = ",".join("?" for _ in PRIMARY_TARGET)
    connection.execute(
        f"UPDATE sites SET state='mapped_nonprimary' WHERE state=? AND id IN (SELECT id FROM forward WHERE mappings=1 AND end-start=2 AND chrom NOT IN ({placeholders}))",
        (PENDING_CHAIN, *sorted(PRIMARY_TARGET)),
    )
    connection.execute(
        "UPDATE sites SET state=?, target_chrom=(SELECT chrom FROM forward WHERE forward.id=sites.id), target_start=(SELECT start FROM forward WHERE forward.id=sites.id), target_end=(SELECT end FROM forward WHERE forward.id=sites.id) WHERE state=? AND id IN (SELECT id FROM forward WHERE mappings=1)",
        (PENDING_TARGET_CHECK, PENDING_CHAIN),
    )
    connection.commit()
    with target_bed.open("x", encoding="utf-8") as handle:
        for identifier, chrom, start, end in connection.execute(
            "SELECT id,target_chrom,target_start,target_end FROM sites WHERE state=? AND target_chrom IS NOT NULL ORDER BY id",
            (PENDING_TARGET_CHECK,),
        ):
            handle.write(f"{chrom}\t{start}\t{end}\tcpg_{identifier:09d}\n")
    connection.close()


def classify_reverse(database: Path) -> None:
    connection = sqlite3.connect(database)
    connection.execute(
        "UPDATE sites SET state=CASE WHEN id IN ("
        "SELECT s.id FROM sites s JOIN reverse r ON s.id=r.id "
        "WHERE r.mappings=1 AND s.source_chrom=r.chrom AND s.source_start=r.start AND s.source_end=r.end"
        ") THEN 'mapped_unique_cpg' ELSE 'roundtrip_failed' END WHERE state=?",
        (PENDING_TARGET_CHECK,),
    )
    connection.commit(); connection.close()


def export_crosswalk(database: Path, output: Path) -> Counter[str]:
    connection = sqlite3.connect(database)
    counts: Counter[str] = Counter({state: 0 for state in FINAL_STATES})
    with output.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as binary:
            header = "cpg_id\tsource_chrom\tsource_start0\tsource_end0\ttarget_chrom\ttarget_start0\ttarget_end0\tmapping_state\n"
            binary.write(header.encode())
            for row in connection.execute(
                "SELECT id,source_chrom,source_start,source_end,COALESCE(target_chrom,''),COALESCE(target_start,''),COALESCE(target_end,''),state FROM sites ORDER BY id"
            ):
                counts[str(row[7])] += 1
                binary.write((f"cpg_{row[0]:09d}\t" + "\t".join(str(value) for value in row[1:]) + "\n").encode())
    connection.close()
    if sum(counts.values()) == 0 or not set(counts) <= FINAL_STATES:
        raise GSE105127CrosswalkError("final crosswalk state census differs")
    return counts


def build_crosswalk(
    *, plan_root: Path, reference_root: Path, collapsed_root: Path,
    target_fasta: Path, output: Path, work: Path,
) -> dict[str, object]:
    if output.exists() or work.exists():
        raise GSE105127CrosswalkError("refusing to overwrite crosswalk output or work")
    verify_frozen_tree(plan_root); verify_frozen_tree(reference_root)
    reference = json.loads((reference_root / "receipt.json").read_text(encoding="utf-8"))
    if (
        reference.get("status") != "passed"
        or reference.get("labels_accessed") is not False
        or reference.get("fit_or_score_performed") is not False
    ):
        raise GSE105127CrosswalkError("reference bundle differs")
    if sha256_file(target_fasta) != TARGET_DECOMPRESSED_SHA256:
        raise GSE105127CrosswalkError("decompressed GRCh38.p14 FASTA identity differs")
    for executable in ("bedtools", "CrossMap"):
        if shutil.which(executable) is None:
            raise GSE105127CrosswalkError(f"required executable is absent: {executable}")
    paths = collect_collapsed(plan_root, collapsed_root)
    output.mkdir(parents=True); work.mkdir(parents=True)
    union = work / "union.tsv"
    merge_union(paths, union)
    database = work / "crosswalk.sqlite"
    source_check = work / "source.check.bed"
    unique_sites = create_database(union, database, source_check)
    source_sequence = work / "source.sequence.tsv"
    getfasta_states(database, reference_root / "human_g1k_v37.fasta", source_check, source_sequence, False)
    valid_source = work / "source.valid.hg19.bed"
    connection = sqlite3.connect(database)
    with valid_source.open("x", encoding="utf-8") as handle:
        for identifier, chrom, start, end in connection.execute(
            "SELECT id,source_chrom,source_start,source_end FROM sites WHERE state=? ORDER BY id",
            (PENDING_CHAIN,),
        ):
            handle.write(f"{chrom}\t{start}\t{end}\tcpg_{identifier:09d}\n")
    connection.close()
    forward = work / "forward.hg38.bed"
    forward_unmap = work / "forward.unmap.bed"
    _run(["CrossMap", "bed", "--chromid", "l", "--unmap-file", str(forward_unmap), str(reference_root / "hg19ToHg38.over.chain.gz"), str(valid_source), str(forward)])
    load_mapping(database, forward, "forward")
    target_bed = work / "target.check.bed"
    classify_forward(database, target_bed)
    target_sequence = work / "target.sequence.tsv"
    getfasta_states(database, target_fasta, target_bed, target_sequence, True)
    reverse_input = work / "reverse.input.hg38.bed"
    connection = sqlite3.connect(database)
    with reverse_input.open("x", encoding="utf-8") as handle:
        for identifier, chrom, start, end in connection.execute(
            "SELECT id,target_chrom,target_start,target_end FROM sites WHERE state=? ORDER BY id",
            (PENDING_TARGET_CHECK,),
        ):
            handle.write(f"{chrom}\t{start}\t{end}\tcpg_{identifier:09d}\n")
    connection.close()
    reverse = work / "reverse.hg19.bed"
    reverse_unmap = work / "reverse.unmap.bed"
    _run(["CrossMap", "bed", "--chromid", "l", "--unmap-file", str(reverse_unmap), str(reference_root / "hg38ToHg19.over.chain.gz"), str(reverse_input), str(reverse)])
    load_mapping(database, reverse, "reverse")
    classify_reverse(database)
    crosswalk = output / "cpg_crosswalk.tsv.gz"
    counts = export_crosswalk(database, crosswalk)
    receipt = {
        "schema_version": "masld-bench-gse105127-cpg-crosswalk-v1",
        "status": "passed_failure_aware_roundtrip",
        "source_build": "1000_Genomes_GRCh37",
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
        raise GSE105127CrosswalkError("crosswalk silently dropped a source CpG")
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    freeze_tree(output, {"artifact_class": "gse105127_failure_aware_cpg_crosswalk", "status": "passed"})
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
    args = parser.parse_args()
    print(json.dumps(build_crosswalk(
        plan_root=args.plan_root, reference_root=args.reference_root,
        collapsed_root=args.collapsed_root, target_fasta=args.target_fasta,
        output=args.output, work=args.work,
    ), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
