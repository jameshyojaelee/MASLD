#!/usr/bin/env python3
"""Freeze the exact GRCh37 source, lift chains, and GRCh38.p14 target requirements."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import md5, sha256
import json
from pathlib import Path
import shutil
import subprocess
import urllib.request

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


SOURCE_BYTES = 892_331_003
SOURCE_MD5_UNCOMPRESSED = "0ce84c872fc0072a885926823dcd0338"
FORWARD_MD5 = "35887f73fe5e2231656504d1f6430900"
REVERSE_MD5 = "ff3031d93792f4cbb86af44055efd903"
TARGET_FASTA_SHA256 = "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215"
TARGET_GTF_SHA256 = "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
PRIMARY_SOURCE = tuple(str(value) for value in range(1, 23)) + ("X", "Y", "MT")
PRIMARY_TARGET = tuple(f"chr{value}" for value in range(1, 23)) + ("chrX", "chrY", "chrM")


class GSE105127ReferenceError(RuntimeError):
    """Raised when a reference asset differs from its frozen requirement."""


def digest_file(path: Path, algorithm: str) -> str:
    digest = {"md5": md5, "sha256": sha256}[algorithm]()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path, max_bytes: int) -> dict[str, object]:
    temporary = path.with_suffix(path.suffix + ".part")
    if path.exists() or temporary.exists():
        raise GSE105127ReferenceError(f"refusing to overwrite download: {path}")
    request = urllib.request.Request(url, headers={"User-Agent": "masld-bench-reference/1.0"})
    digest = sha256()
    size = 0
    with urllib.request.urlopen(request, timeout=300) as response, temporary.open("xb") as handle:
        declared = response.headers.get("Content-Length")
        while True:
            block = response.read(8 * 1024 * 1024)
            if not block:
                break
            size += len(block)
            if size > max_bytes:
                raise GSE105127ReferenceError(f"reference download exceeds cap: {url}")
            digest.update(block)
            handle.write(block)
        if declared is not None and int(declared) != size:
            raise GSE105127ReferenceError("reference Content-Length differs")
    temporary.rename(path)
    return {"url": url, "bytes": size, "sha256": digest.hexdigest()}


def fasta_inventory(path: Path) -> list[dict[str, object]]:
    rows = []
    name: str | None = None
    length = 0
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    rows.append({"contig": name, "length": length})
                name = line[1:].split(None, 1)[0].decode("ascii")
                length = 0
            else:
                sequence = line.strip()
                if any(value not in b"ACGTNacgtn" for value in sequence):
                    raise GSE105127ReferenceError("FASTA has a non-IUPAC primary sequence base")
                length += len(sequence)
    if name is not None:
        rows.append({"contig": name, "length": length})
    if not rows or len({str(row["contig"]) for row in rows}) != len(rows):
        raise GSE105127ReferenceError("FASTA contig inventory differs")
    return rows


def build_reference_bundle(
    *, plan_root: Path,
    target_fasta: Path,
    target_gtf: Path,
    output: Path,
) -> dict[str, object]:
    if output.exists():
        raise GSE105127ReferenceError(f"refusing to overwrite reference bundle: {output}")
    verify_frozen_tree(plan_root)
    plan = json.loads((plan_root / "references.json").read_text(encoding="utf-8"))
    if digest_file(target_fasta, "sha256") != TARGET_FASTA_SHA256:
        raise GSE105127ReferenceError("GRCh38.p14 FASTA SHA-256 differs")
    if digest_file(target_gtf, "sha256") != TARGET_GTF_SHA256:
        raise GSE105127ReferenceError("GENCODE v49 GTF SHA-256 differs")
    for path in (target_fasta, target_gtf):
        result = subprocess.run(["gzip", "-t", str(path)], check=False)
        if result.returncode != 0:
            raise GSE105127ReferenceError(f"gzip integrity failed: {path}")
    if shutil.which("samtools") is None:
        raise GSE105127ReferenceError("samtools is unavailable")
    output.mkdir(parents=True)
    source_gz = output / "human_g1k_v37.fasta.gz"
    forward = output / "hg19ToHg38.over.chain.gz"
    reverse = output / "hg38ToHg19.over.chain.gz"
    receipts = [
        download(plan["source_fasta_url"], source_gz, SOURCE_BYTES),
        download(plan["forward_chain_url"], forward, 16 * 1024 * 1024),
        download(plan["reverse_chain_url"], reverse, 16 * 1024 * 1024),
    ]
    if source_gz.stat().st_size != SOURCE_BYTES:
        raise GSE105127ReferenceError("1000 Genomes compressed FASTA size differs")
    if digest_file(forward, "md5") != FORWARD_MD5 or digest_file(reverse, "md5") != REVERSE_MD5:
        raise GSE105127ReferenceError("UCSC lift-chain MD5 differs")
    for path in (source_gz, forward, reverse):
        if subprocess.run(["gzip", "-t", str(path)], check=False).returncode != 0:
            raise GSE105127ReferenceError(f"downloaded gzip integrity failed: {path}")
    source_fa = output / "human_g1k_v37.fasta"
    with gzip.open(source_gz, "rb") as source, source_fa.open("xb") as target:
        shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
    if digest_file(source_fa, "md5") != SOURCE_MD5_UNCOMPRESSED:
        raise GSE105127ReferenceError("1000 Genomes uncompressed FASTA MD5 differs")
    inventory = fasta_inventory(source_fa)
    if tuple(str(row["contig"]) for row in inventory) != PRIMARY_SOURCE:
        raise GSE105127ReferenceError("1000 Genomes primary contig order differs")
    subprocess.run(["samtools", "faidx", str(source_fa)], check=True)
    fai = list(csv.reader((output / "human_g1k_v37.fasta.fai").open(), delimiter="\t"))
    if [row[0] for row in fai] != list(PRIMARY_SOURCE):
        raise GSE105127ReferenceError("rebuilt source FASTA index differs")
    target_fai = target_fasta.parent.parent / "fasta/genome.fa.fai"
    if not target_fai.is_file():
        raise GSE105127ReferenceError("project target FASTA index is absent")
    target_names = [line.split("\t", 1)[0] for line in target_fai.read_text().splitlines()]
    if not set(PRIMARY_TARGET) <= set(target_names):
        raise GSE105127ReferenceError("GRCh38.p14 primary target contigs are incomplete")
    receipt = {
        "schema_version": "masld-bench-gse105127-reference-bundle-v1",
        "status": "passed",
        "source_assembly": "1000_Genomes_GRCh37_human_g1k_v37",
        "source_compressed_sha256": digest_file(source_gz, "sha256"),
        "source_uncompressed_sha256": digest_file(source_fa, "sha256"),
        "source_uncompressed_md5": SOURCE_MD5_UNCOMPRESSED,
        "source_primary_contigs": len(inventory),
        "forward_chain_sha256": digest_file(forward, "sha256"),
        "reverse_chain_sha256": digest_file(reverse, "sha256"),
        "target_fasta_path": str(target_fasta),
        "target_fasta_sha256": TARGET_FASTA_SHA256,
        "target_gtf_path": str(target_gtf),
        "target_gtf_sha256": TARGET_GTF_SHA256,
        "labels_accessed": False,
        "fit_or_score_performed": False,
        "download_receipts": receipts,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(output, {"artifact_class": "gse105127_reference_bundle", "status": "passed"})
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--target-fasta", type=Path, required=True)
    parser.add_argument("--target-gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build_reference_bundle(
        plan_root=args.plan_root, target_fasta=args.target_fasta,
        target_gtf=args.target_gtf, output=args.output,
    ), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
