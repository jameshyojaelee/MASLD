#!/usr/bin/env python3
"""Build the exact GRCh38.p14/GENCODE v49 RSEM+STAR reference."""

from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess

from masld_bench.artifacts import freeze_tree, verify_frozen_tree


FASTA_SHA256 = "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215"
GTF_SHA256 = "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"


class GSE105127RSEMReferenceError(RuntimeError):
    """Raised when the RNA quantification reference differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_output(command: list[str]) -> str:
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise GSE105127RSEMReferenceError(f"version probe failed: {command}")
    return (result.stdout + result.stderr).strip()


def build_reference(
    *, plan_root: Path,
    target_fasta: Path,
    target_gtf: Path,
    output: Path,
    threads: int,
) -> dict[str, object]:
    if output.exists():
        raise GSE105127RSEMReferenceError(f"refusing to overwrite RSEM reference: {output}")
    verify_frozen_tree(plan_root)
    if sha256_file(target_fasta) != FASTA_SHA256 or sha256_file(target_gtf) != GTF_SHA256:
        raise GSE105127RSEMReferenceError("target FASTA or GTF identity differs")
    if shutil.which("rsem-prepare-reference") is None or shutil.which("STAR") is None:
        raise GSE105127RSEMReferenceError("RSEM or STAR is unavailable")
    output.mkdir(parents=True)
    inputs = output / "inputs"
    inputs.mkdir()
    fasta = inputs / "GRCh38.p14.genome.fa"
    gtf = inputs / "gencode.v49.annotation.gtf"
    for source, target in ((target_fasta, fasta), (target_gtf, gtf)):
        with gzip.open(source, "rb") as read_handle, target.open("xb") as write_handle:
            shutil.copyfileobj(read_handle, write_handle, length=8 * 1024 * 1024)
    prefix_root = output / "reference"
    prefix_root.mkdir()
    prefix = prefix_root / "gse105127_gencode_v49"
    stdout = output / "rsem_prepare.stdout.txt"
    stderr = output / "rsem_prepare.stderr.txt"
    command = [
        "rsem-prepare-reference", "--gtf", str(gtf), "--star",
        "--star-path", str(Path(shutil.which("STAR") or "").parent),
        "--star-sjdboverhang", "75",
        "--num-threads", str(threads), str(fasta), str(prefix),
    ]
    with stdout.open("x", encoding="utf-8") as out, stderr.open("x", encoding="utf-8") as err:
        result = subprocess.run(command, stdout=out, stderr=err, text=True, check=False)
    if result.returncode != 0:
        raise GSE105127RSEMReferenceError("rsem-prepare-reference failed")
    required_suffixes = (".grp", ".ti", ".seq", ".chrlist", ".transcripts.fa")
    if any(
        not Path(str(prefix) + suffix).is_file()
        or Path(str(prefix) + suffix).stat().st_size == 0
        for suffix in required_suffixes
    ):
        raise GSE105127RSEMReferenceError("RSEM reference members are incomplete")
    if any(
        not (prefix_root / name).is_file() or (prefix_root / name).stat().st_size == 0
        for name in ("Genome", "SA", "SAindex")
    ):
        raise GSE105127RSEMReferenceError("STAR reference members are incomplete")
    versions = {
        "rsem": command_output(["rsem-calculate-expression", "--version"]),
        "star": command_output(["STAR", "--version"]),
    }
    if "RSEM v1.3.1" not in versions["rsem"] or versions["star"] != "2.7.10b":
        raise GSE105127RSEMReferenceError("RSEM/STAR versions differ")
    fasta.unlink()
    gtf.unlink()
    inputs.rmdir()
    receipt = {
        "schema_version": "masld-bench-gse105127-rsem-reference-v1",
        "status": "passed",
        "assembly": "GRCh38.p14",
        "annotation": "GENCODE_v49",
        "fasta_sha256": FASTA_SHA256,
        "gtf_sha256": GTF_SHA256,
        "rsem_version": versions["rsem"],
        "star_version": versions["star"],
        "star_sjdboverhang": 75,
        "reference_prefix": "reference/gse105127_gencode_v49",
        "labels_accessed": False,
        "fit_or_score_performed": False,
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    freeze_tree(output, {"artifact_class": "gse105127_rsem_star_reference", "status": "passed"})
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--target-fasta", type=Path, required=True)
    parser.add_argument("--target-gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    args = parser.parse_args()
    print(json.dumps(build_reference(
        plan_root=args.plan_root, target_fasta=args.target_fasta,
        target_gtf=args.target_gtf, output=args.output, threads=args.threads,
    ), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
