#!/usr/bin/env python3
"""Freeze GSE105127 references while admitting one exact legacy gzip stream."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import shutil
import subprocess

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from scripts.freeze_gse105127_reference_bundle import (
    FORWARD_MD5,
    GSE105127ReferenceError,
    PRIMARY_SOURCE,
    PRIMARY_TARGET,
    REVERSE_MD5,
    SOURCE_BYTES,
    TARGET_FASTA_SHA256,
    TARGET_GTF_SHA256,
    digest_file,
)


SOURCE_SHA256_COMPRESSED = "8b6c538abf0dd92d3f3020f36cc1dd67ce004ffa421c2781205f1eb690bdb442"
SOURCE_BYTES_UNCOMPRESSED = 3_153_506_519
SOURCE_MD5_UNCOMPRESSED = "0ce84c872fc0072a885926823dcd0338"
SOURCE_SHA256_UNCOMPRESSED = "2f9cd9e853a9284c53884e6a551b1c7284795dd053f255d630aeeb114d1fa81f"
LEGACY_DIAGNOSTIC_ARTIFACTS_SHA256 = "0b69e0d7a941b4b86f50b62b66534df503163c0eb2e46d219db6964047f550fb"
ALPHABET_DIAGNOSTIC_ARTIFACTS_SHA256 = "2d2407dc694c20890f5ffd21a05bffe2fa3d671300c92cd725ca32d173e1dd0b"
SOURCE_CACHE_ARTIFACTS_SHA256 = "171eadb68760252213acb936f5e59ba85263e2e40041d98705d59233cf3ce0a9"
LEGACY_ALPHABET_COUNTS = {
    "A": 845_903_867,
    "C": 585_709_826,
    "G": 586_026_532,
    "M": 1,
    "N": 237_019_516,
    "R": 2,
    "T": 847_144_995,
}
SOURCE_CONTIGS_TOTAL = 84


def validate_legacy_diagnostic(receipt: dict[str, object]) -> None:
    expected = {
        "schema_version": "masld-bench-gse105127-legacy-reference-stream-diagnostic-v1",
        "status": "pass_exact_legacy_stream",
        "compressed_bytes": SOURCE_BYTES,
        "compressed_sha256": SOURCE_SHA256_COMPRESSED,
        "gzip_exit_status": 2,
        "legacy_warning": "trailing garbage ignored",
        "uncompressed_bytes": SOURCE_BYTES_UNCOMPRESSED,
        "uncompressed_md5": SOURCE_MD5_UNCOMPRESSED,
        "uncompressed_sha256": SOURCE_SHA256_UNCOMPRESSED,
        "fatal_gzip_diagnostic": False,
        "fit_or_score_performed": False,
        "labels_accessed": False,
    }
    if receipt != expected:
        raise GSE105127ReferenceError("legacy reference diagnostic differs")


def validate_alphabet_diagnostic(receipt: dict[str, object]) -> None:
    expected = {
        "schema_version": "masld-bench-gse105127-reference-alphabet-diagnostic-v1",
        "status": "passed_outcome_free_census",
        "source_bytes": SOURCE_BYTES_UNCOMPRESSED,
        "source_sha256": SOURCE_SHA256_UNCOMPRESSED,
        "header_rows": 84,
        "sequence_rows": 51_696_747,
        "alphabet_counts": LEGACY_ALPHABET_COUNTS,
        "alphabet": sorted(LEGACY_ALPHABET_COUNTS),
        "labels_accessed": False,
        "fit_or_score_performed": False,
    }
    if receipt != expected:
        raise GSE105127ReferenceError("legacy FASTA alphabet diagnostic differs")


def fasta_inventory_legacy_exact(path: Path) -> list[dict[str, object]]:
    """Inventory the exact legacy FASTA without broadening its observed alphabet."""
    rows = []
    counts: Counter[int] = Counter()
    name: str | None = None
    length = 0
    allowed = frozenset(b"ACGMRNT")
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    rows.append({"contig": name, "length": length})
                name = line[1:].split(None, 1)[0].decode("ascii")
                length = 0
                continue
            sequence = line.strip()
            counts.update(sequence)
            if any(value not in allowed for value in sequence):
                raise GSE105127ReferenceError("legacy FASTA alphabet differs")
            length += len(sequence)
    if name is not None:
        rows.append({"contig": name, "length": length})
    if not rows or len({str(row["contig"]) for row in rows}) != len(rows):
        raise GSE105127ReferenceError("FASTA contig inventory differs")
    observed = {chr(value): count for value, count in sorted(counts.items())}
    if observed != LEGACY_ALPHABET_COUNTS:
        raise GSE105127ReferenceError("legacy FASTA alphabet census differs")
    return rows


def validate_source_inventory_names(inventory: list[dict[str, object]]) -> list[str]:
    names = [str(row["contig"]) for row in inventory]
    if len(names) != SOURCE_CONTIGS_TOTAL:
        raise GSE105127ReferenceError("1000 Genomes contig count differs")
    if tuple(names[: len(PRIMARY_SOURCE)]) != PRIMARY_SOURCE:
        raise GSE105127ReferenceError("1000 Genomes primary contig prefix differs")
    return names


def validate_source_cache(root: Path) -> None:
    verify_frozen_tree(root)
    if digest_file(root / "ARTIFACTS.json", "sha256") != SOURCE_CACHE_ARTIFACTS_SHA256:
        raise GSE105127ReferenceError("reference source cache ARTIFACTS SHA-256 differs")
    receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
    expected_files = [
        {
            "algorithm": "sha256",
            "bytes": SOURCE_BYTES,
            "digest": SOURCE_SHA256_COMPRESSED,
            "name": "human_g1k_v37.fasta.gz",
        },
        {
            "algorithm": "md5",
            "bytes": 227_698,
            "digest": FORWARD_MD5,
            "name": "hg19ToHg38.over.chain.gz",
        },
        {
            "algorithm": "md5",
            "bytes": 1_246_411,
            "digest": REVERSE_MD5,
            "name": "hg38ToHg19.over.chain.gz",
        },
    ]
    expected = {
        "schema_version": "masld-bench-gse105127-reference-source-cache-v1",
        "status": "passed_exact_source_bytes",
        "files": expected_files,
        "labels_accessed": False,
        "fit_or_score_performed": False,
    }
    if receipt != expected:
        raise GSE105127ReferenceError("reference source cache receipt differs")


def copy_cached_source(source: Path, target: Path) -> dict[str, object]:
    if target.exists():
        raise GSE105127ReferenceError(f"refusing to overwrite cached source: {target}")
    with source.open("rb") as input_handle, target.open("xb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle, length=8 * 1024 * 1024)
    return {
        "source_cache_path": str(source),
        "bytes": target.stat().st_size,
        "sha256": digest_file(target, "sha256"),
    }


def decompress_exact_legacy_source(source: Path, target: Path, stderr_path: Path) -> dict[str, object]:
    if source.stat().st_size != SOURCE_BYTES:
        raise GSE105127ReferenceError("1000 Genomes compressed FASTA size differs")
    if digest_file(source, "sha256") != SOURCE_SHA256_COMPRESSED:
        raise GSE105127ReferenceError("1000 Genomes compressed FASTA SHA-256 differs")
    if target.exists() or stderr_path.exists():
        raise GSE105127ReferenceError("refusing to overwrite decompressed reference")
    with target.open("xb") as output:
        result = subprocess.run(
            ["gzip", "-cd", str(source)], stdout=output, stderr=subprocess.PIPE, check=False
        )
    stderr = result.stderr.decode("utf-8", errors="strict")
    stderr_path.write_text(stderr, encoding="utf-8")
    normalized = stderr.strip()
    if result.returncode != 2 or not normalized.endswith(
        ": decompression OK, trailing garbage ignored"
    ):
        raise GSE105127ReferenceError("legacy gzip disposition differs")
    if normalized.count("\n") != 0:
        raise GSE105127ReferenceError("legacy gzip emitted an additional diagnostic")
    observed = {
        "bytes": target.stat().st_size,
        "md5": digest_file(target, "md5"),
        "sha256": digest_file(target, "sha256"),
    }
    if observed != {
        "bytes": SOURCE_BYTES_UNCOMPRESSED,
        "md5": SOURCE_MD5_UNCOMPRESSED,
        "sha256": SOURCE_SHA256_UNCOMPRESSED,
    }:
        raise GSE105127ReferenceError("1000 Genomes uncompressed FASTA identity differs")
    return observed


def build_reference_bundle(
    *,
    plan_root: Path,
    target_fasta: Path,
    target_gtf: Path,
    legacy_diagnostic_root: Path,
    alphabet_diagnostic_root: Path,
    source_cache_root: Path,
    implementation_sha256: str,
    output: Path,
) -> dict[str, object]:
    if output.exists():
        raise GSE105127ReferenceError(f"refusing to overwrite reference bundle: {output}")
    if digest_file(Path(__file__), "sha256") != implementation_sha256:
        raise GSE105127ReferenceError("exact reference implementation SHA-256 differs")
    verify_frozen_tree(plan_root)
    verify_frozen_tree(legacy_diagnostic_root)
    verify_frozen_tree(alphabet_diagnostic_root)
    validate_source_cache(source_cache_root)
    if digest_file(legacy_diagnostic_root / "ARTIFACTS.json", "sha256") != LEGACY_DIAGNOSTIC_ARTIFACTS_SHA256:
        raise GSE105127ReferenceError("legacy diagnostic ARTIFACTS SHA-256 differs")
    validate_legacy_diagnostic(json.loads(
        (legacy_diagnostic_root / "diagnostic_receipt.json").read_text(encoding="utf-8")
    ))
    if digest_file(alphabet_diagnostic_root / "ARTIFACTS.json", "sha256") != ALPHABET_DIAGNOSTIC_ARTIFACTS_SHA256:
        raise GSE105127ReferenceError("alphabet diagnostic ARTIFACTS SHA-256 differs")
    validate_alphabet_diagnostic(json.loads(
        (alphabet_diagnostic_root / "alphabet_receipt.json").read_text(encoding="utf-8")
    ))
    plan = json.loads((plan_root / "references.json").read_text(encoding="utf-8"))
    if digest_file(target_fasta, "sha256") != TARGET_FASTA_SHA256:
        raise GSE105127ReferenceError("GRCh38.p14 FASTA SHA-256 differs")
    if digest_file(target_gtf, "sha256") != TARGET_GTF_SHA256:
        raise GSE105127ReferenceError("GENCODE v49 GTF SHA-256 differs")
    for path in (target_fasta, target_gtf):
        if subprocess.run(["gzip", "-t", str(path)], check=False).returncode != 0:
            raise GSE105127ReferenceError(f"gzip integrity failed: {path}")
    if shutil.which("samtools") is None:
        raise GSE105127ReferenceError("samtools is unavailable")

    output.mkdir(parents=True)
    source_gz = output / "human_g1k_v37.fasta.gz"
    forward = output / "hg19ToHg38.over.chain.gz"
    reverse = output / "hg38ToHg19.over.chain.gz"
    receipts = [
        copy_cached_source(source_cache_root / source_gz.name, source_gz),
        copy_cached_source(source_cache_root / forward.name, forward),
        copy_cached_source(source_cache_root / reverse.name, reverse),
    ]
    for receipt, url in zip(receipts, (
        plan["source_fasta_url"],
        plan["forward_chain_url"],
        plan["reverse_chain_url"],
    )):
        receipt["provenance_url"] = url
    if source_gz.stat().st_size != SOURCE_BYTES or digest_file(source_gz, "sha256") != SOURCE_SHA256_COMPRESSED:
        raise GSE105127ReferenceError("cached 1000 Genomes compressed FASTA differs")
    if digest_file(forward, "md5") != FORWARD_MD5 or digest_file(reverse, "md5") != REVERSE_MD5:
        raise GSE105127ReferenceError("UCSC lift-chain MD5 differs")
    for path in (forward, reverse):
        if subprocess.run(["gzip", "-t", str(path)], check=False).returncode != 0:
            raise GSE105127ReferenceError(f"downloaded gzip integrity failed: {path}")
    source_fa = output / "human_g1k_v37.fasta"
    legacy_stream = decompress_exact_legacy_source(
        source_gz, source_fa, output / "human_g1k_v37.fasta.decompression.stderr.txt"
    )
    inventory = fasta_inventory_legacy_exact(source_fa)
    source_names = validate_source_inventory_names(inventory)
    subprocess.run(["samtools", "faidx", str(source_fa)], check=True)
    fai = list(csv.reader((output / "human_g1k_v37.fasta.fai").open(), delimiter="\t"))
    if [row[0] for row in fai] != source_names:
        raise GSE105127ReferenceError("rebuilt source FASTA index differs")
    target_fai = target_fasta.parent.parent / "fasta/genome.fa.fai"
    if not target_fai.is_file():
        raise GSE105127ReferenceError("project target FASTA index is absent")
    target_names = [line.split("\t", 1)[0] for line in target_fai.read_text().splitlines()]
    if not set(PRIMARY_TARGET) <= set(target_names):
        raise GSE105127ReferenceError("GRCh38.p14 primary target contigs are incomplete")

    receipt = {
        "schema_version": "masld-bench-gse105127-reference-bundle-v2",
        "status": "passed_exact_legacy_stream",
        "source_assembly": "1000_Genomes_GRCh37_human_g1k_v37",
        "source_compressed_sha256": SOURCE_SHA256_COMPRESSED,
        "source_uncompressed_sha256": SOURCE_SHA256_UNCOMPRESSED,
        "source_uncompressed_md5": SOURCE_MD5_UNCOMPRESSED,
        "source_contigs_total": len(inventory),
        "source_primary_contigs": len(PRIMARY_SOURCE),
        "source_non_primary_contigs": len(inventory) - len(PRIMARY_SOURCE),
        "forward_chain_sha256": digest_file(forward, "sha256"),
        "reverse_chain_sha256": digest_file(reverse, "sha256"),
        "target_fasta_path": str(target_fasta),
        "target_fasta_sha256": TARGET_FASTA_SHA256,
        "target_gtf_path": str(target_gtf),
        "target_gtf_sha256": TARGET_GTF_SHA256,
        "implementation_sha256": implementation_sha256,
        "legacy_diagnostic_root": str(legacy_diagnostic_root),
        "legacy_diagnostic_artifacts_sha256": LEGACY_DIAGNOSTIC_ARTIFACTS_SHA256,
        "alphabet_diagnostic_root": str(alphabet_diagnostic_root),
        "alphabet_diagnostic_artifacts_sha256": ALPHABET_DIAGNOSTIC_ARTIFACTS_SHA256,
        "source_cache_root": str(source_cache_root),
        "source_cache_artifacts_sha256": SOURCE_CACHE_ARTIFACTS_SHA256,
        "source_alphabet_counts": LEGACY_ALPHABET_COUNTS,
        "legacy_stream_admission": legacy_stream,
        "labels_accessed": False,
        "fit_or_score_performed": False,
        "download_receipts": receipts,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(output, {
        "artifact_class": "gse105127_reference_bundle",
        "status": "passed_exact_legacy_stream",
    })
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--target-fasta", type=Path, required=True)
    parser.add_argument("--target-gtf", type=Path, required=True)
    parser.add_argument("--legacy-diagnostic-root", type=Path, required=True)
    parser.add_argument("--alphabet-diagnostic-root", type=Path, required=True)
    parser.add_argument("--source-cache-root", type=Path, required=True)
    parser.add_argument("--implementation-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build_reference_bundle(
        plan_root=args.plan_root,
        target_fasta=args.target_fasta,
        target_gtf=args.target_gtf,
        legacy_diagnostic_root=args.legacy_diagnostic_root,
        alphabet_diagnostic_root=args.alphabet_diagnostic_root,
        source_cache_root=args.source_cache_root,
        implementation_sha256=args.implementation_sha256,
        output=args.output,
    ), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
