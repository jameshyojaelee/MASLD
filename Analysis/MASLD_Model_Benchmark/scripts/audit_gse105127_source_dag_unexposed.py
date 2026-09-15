#!/usr/bin/env python3
"""Independently audit the source-only GSE105127 production DAG.

The audit reads plans, manifests, receipts, reference bases, and compressed
source bytes.  It does not parse RNA quantification outputs, methylation values,
labels, outcomes, predictions, or evaluator assets, and it never fits or scores.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
from typing import Any, Iterable, Mapping, Sequence
import zlib

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


CONTROL_SHA256 = "4cdd6acecfcbef95ece3c45b7319a609f2555039b74916f8f8db0552a9521a1b"
PLAN_SHA256 = "c4f91d2a3fedc7b3b1802a0d6cf3d1b1803ce1bea70cc343f6e18f744783b973"
REFERENCE_SHA256 = "9575aa24dc39f855b961eecf53e29243b4a45fca6770c7337e05901bc7d1ec89"
RSEM_REFERENCE_SHA256 = "137f103ec6676681deb3354f2e0942476d6a24ebfaddaa46ecf9c979c0946dfa"
LEGACY_DIAGNOSTIC_SHA256 = "0b69e0d7a941b4b86f50b62b66534df503163c0eb2e46d219db6964047f550fb"
ALPHABET_DIAGNOSTIC_SHA256 = "2d2407dc694c20890f5ffd21a05bffe2fa3d671300c92cd725ca32d173e1dd0b"
SOURCE_CACHE_SHA256 = "171eadb68760252213acb936f5e59ba85263e2e40041d98705d59233cf3ce0a9"
SOURCE_COMPRESSED_BYTES = 892_331_003
SOURCE_COMPRESSED_SHA256 = "8b6c538abf0dd92d3f3020f36cc1dd67ce004ffa421c2781205f1eb690bdb442"
SOURCE_UNCOMPRESSED_BYTES = 3_153_506_519
SOURCE_UNCOMPRESSED_MD5 = "0ce84c872fc0072a885926823dcd0338"
SOURCE_UNCOMPRESSED_SHA256 = "2f9cd9e853a9284c53884e6a551b1c7284795dd053f255d630aeeb114d1fa81f"
TARGET_FASTA_SHA256 = "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215"
TARGET_GTF_SHA256 = "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
FORWARD_CHAIN_MD5 = "35887f73fe5e2231656504d1f6430900"
REVERSE_CHAIN_MD5 = "ff3031d93792f4cbb86af44055efd903"
EXPECTED_PARTICIPANTS = 19
EXPECTED_ROWS = 57
EXPECTED_BUNDLES = 8
EXPECTED_RNA_BYTES = 79_576_393_854
EXPECTED_RNA_READS = 1_140_779_295
EXPECTED_RNA_BASES = 86_699_226_420
EXPECTED_RRBS_BYTES = 8_192_644_943
EXPECTED_ALPHABET_COUNTS = {
    "A": 845_903_867,
    "C": 585_709_826,
    "G": 586_026_532,
    "M": 1,
    "N": 237_019_516,
    "R": 2,
    "T": 847_144_995,
}
PRIMARY_SOURCE = tuple([str(value) for value in range(1, 23)] + ["X", "Y", "MT"])
RNA_FIELDS = (
    "bundle_id", "row_id", "participant_group_id", "zone", "pairing_topology",
    "sample_accession", "run_accession", "library_layout", "read_length",
    "read_count", "base_count", "fastq_url", "fastq_bytes", "fastq_md5",
    "relative_fastq_path",
)
RRBS_FIELDS = (
    "bundle_id", "row_id", "participant_group_id", "zone", "pairing_topology",
    "sample_accession", "source_bed_path", "source_bed_bytes", "source_bed_sha256",
    "source_cytosine_rows", "relative_collapsed_path",
)
BUNDLE_FIELDS = ("bundle_id", "participants", "participant_zone_rows", "combined_source_bytes")
FORBIDDEN_DATA_SURFACES = (
    "evaluator_only", "outcomes.tsv", "predictions", "config/evaluation",
    "phenotype", "fibrosis", "nas_score", "clinical_covariates",
)
FORBIDDEN_CALLS = frozenset(
    {"fit", "fit_predict", "fit_transform", "predict", "predict_proba", "score", "train"}
)
SOURCE_CODE_HASHES = {
    "scripts/plan_gse105127_production_dag.py": "73aa9ac393273aeb214783baf04dd2085799748807504b655692cfdaed2b9cdb",
    "scripts/freeze_gse105127_reference_bundle_legacy_exact.py": "e3610b025be3c4c9238fb111c4225bdde47f69eb22cc9ebe83309eff3b8b462e",
    "scripts/download_gse105127_rna_bundle_resumable.py": "2b93452ecbefbb38a64721b0a736b90fd1f126e5bc5d04b6c977fd563c5416fe",
    "scripts/collapse_gse105127_rrbs_bundle.py": "ee830a4b51d01df1b1a8beced6e7484ec92448e0d4f5304198776b6911e7d7e6",
    "scripts/build_gse105127_rsem_reference.py": "d360d7b4204bf9603c3796f9569b609ba12a3b1bf9d1091ca6fb1fe714234938",
    "scripts/quantify_gse105127_rna_bundle.py": "43e20d33c2308fce131de129bdd4063bfd1e6c184ee4ce0a89b3aec3befce18f",
    "scripts/consolidate_gse105127_rna_quantification.py": "cdb77fcdaded3f8ac73b93a7c98047f98c4a5a41abd11ff961e4e97c66e87d0f",
    "scripts/build_gse105127_cpg_crosswalk.py": "2676dbeb4db530a90f0816d34df480649d66b50755b0978922d95f418d37736e",
    "scripts/finalize_gse105127_production_activation.py": "fb891b529c545972c14f120c341256c993e3d54cbf13f8032499036d7b9d9ac3",
}
SHELL_CODE_HASHES = {
    "slurm/gse105127_production_common.sh": "1f1fc22613d5320b88b94e1c8bf1e54a75776369cf91d758978136e18627c07f",
    "slurm/gse105127_stage_reference_legacy_exact_io.sbatch": "0ccf56856b16e09629cd8402c7692cebcc9ea806b05f9c61fa0110b16015b273",
    "slurm/gse105127_stage_rna_download_resumable_io.sbatch": "bad698411349f9874d99ef91355c6c69a83b01ded5c95691d210cf33e441cf7b",
    "slurm/gse105127_stage_rsem_reference_cpu.sbatch": "229f7a42da864f3cb8b7b230b862c61109152ab67d5d52ae845e08b055804ac7",
    "slurm/gse105127_stage_rna_quant_cpu.sbatch": "f6b2d2783673fcf9933600b9f7d298b8749d07279ac600c07741a1e490047af6",
    "slurm/gse105127_stage_consolidate_cpu.sbatch": "325ec1abc472f5d274f00b8f99ab96a4664c05a4f1478d2d7be43b6c29dada6e",
    "slurm/gse105127_stage_crosswalk_io.sbatch": "41fa1227b0f9927e8e11011147ddd4129e5e12924e623b1d87c06c124b2270dd",
    "slurm/gse105127_stage_finalize_cpu.sbatch": "e893cbb3aeef1948278ca4d0ad4c1d7a6d7555195b7fe87c4280798f74ec5af2",
}


class GSE105127IndependentAuditError(ValueError):
    """Raised when the source-only GSE105127 DAG differs."""


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE105127IndependentAuditError(f"TSV lacks a header: {path}")
        return tuple(reader.fieldnames), list(reader)


def digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_pinned_tree(root: Path, expected_sha256: str, label: str) -> dict[str, Any]:
    if sha256_file(root / "ARTIFACTS.json") != expected_sha256:
        raise GSE105127IndependentAuditError(f"{label} ARTIFACTS hash differs")
    return verify_frozen_tree(root)


def gzip_test(path: Path) -> tuple[int, str]:
    result = subprocess.run(
        ["gzip", "-t", str(path)], check=False, capture_output=True, text=True
    )
    return result.returncode, result.stderr.strip()


def inspect_legacy_gzip_member(path: Path) -> dict[str, Any]:
    """Authenticate one gzip member and inventory exact bytes trailing it."""
    compressed_sha = hashlib.sha256()
    uncompressed_sha = hashlib.sha256()
    uncompressed_md5 = hashlib.md5(usedforsecurity=False)
    decompressor = zlib.decompressobj(wbits=31)
    total_input = 0
    total_output = 0
    tail_sha = hashlib.sha256()
    tail_bytes = 0
    tail_prefix = bytearray()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            compressed_sha.update(block)
            total_input += len(block)
            if decompressor.eof:
                tail_sha.update(block)
                tail_bytes += len(block)
                tail_prefix.extend(block[: max(0, 16 - len(tail_prefix))])
                continue
            output = decompressor.decompress(block)
            uncompressed_sha.update(output)
            uncompressed_md5.update(output)
            total_output += len(output)
            if decompressor.eof and decompressor.unused_data:
                unused = decompressor.unused_data
                tail_sha.update(unused)
                tail_bytes += len(unused)
                tail_prefix.extend(unused[:16])
    output = decompressor.flush()
    uncompressed_sha.update(output)
    uncompressed_md5.update(output)
    total_output += len(output)
    if not decompressor.eof or decompressor.unconsumed_tail:
        raise GSE105127IndependentAuditError("legacy gzip member is incomplete")
    return {
        "compressed_bytes": total_input,
        "compressed_sha256": compressed_sha.hexdigest(),
        "gzip_member_bytes": total_input - tail_bytes,
        "trailing_bytes": tail_bytes,
        "trailing_sha256": tail_sha.hexdigest(),
        "trailing_prefix_hex": bytes(tail_prefix[:16]).hex(),
        "trailing_starts_with_gzip_magic": bytes(tail_prefix[:2]) == b"\x1f\x8b",
        "uncompressed_bytes": total_output,
        "uncompressed_md5": uncompressed_md5.hexdigest(),
        "uncompressed_sha256": uncompressed_sha.hexdigest(),
    }


def validate_legacy_reference_rule(source_gz: Path, source_fa: Path) -> dict[str, Any]:
    if source_gz.stat().st_size != SOURCE_COMPRESSED_BYTES:
        raise GSE105127IndependentAuditError("legacy compressed source byte count differs")
    member = inspect_legacy_gzip_member(source_gz)
    expected = {
        "compressed_bytes": SOURCE_COMPRESSED_BYTES,
        "compressed_sha256": SOURCE_COMPRESSED_SHA256,
        "uncompressed_bytes": SOURCE_UNCOMPRESSED_BYTES,
        "uncompressed_md5": SOURCE_UNCOMPRESSED_MD5,
        "uncompressed_sha256": SOURCE_UNCOMPRESSED_SHA256,
    }
    for key, value in expected.items():
        if member.get(key) != value:
            raise GSE105127IndependentAuditError(f"legacy gzip {key} differs")
    if member["trailing_bytes"] <= 0 or member["trailing_starts_with_gzip_magic"]:
        raise GSE105127IndependentAuditError("legacy gzip trailing-byte disposition differs")
    status, stderr = gzip_test(source_gz)
    if (
        status != 2
        or not stderr.endswith(": decompression OK, trailing garbage ignored")
        or stderr.count("\n") != 0
    ):
        raise GSE105127IndependentAuditError("legacy compressed source gzip-test disposition differs")
    if (
        source_fa.stat().st_size != SOURCE_UNCOMPRESSED_BYTES
        or digest_file(source_fa, "md5") != SOURCE_UNCOMPRESSED_MD5
        or digest_file(source_fa, "sha256") != SOURCE_UNCOMPRESSED_SHA256
    ):
        raise GSE105127IndependentAuditError("admitted decompressed reference identity differs")
    return {
        **member,
        "gzip_test_exit_status": status,
        "gzip_test_passed": False,
        "compressed_source_admission": "provenance_only_not_sequence_input",
        "decompressed_reference_admission": "eligible_only_by_exact_derived_reference_rule",
        "campaign_wide_gzip_t_conflict": True,
    }


def fasta_inventory(path: Path) -> tuple[list[tuple[str, int]], dict[str, int]]:
    inventory: list[tuple[str, int]] = []
    alphabet: Counter[str] = Counter()
    name: str | None = None
    length = 0
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    inventory.append((name, length))
                name = line[1:].split(None, 1)[0].decode("ascii")
                length = 0
                continue
            sequence = line.strip()
            try:
                text = sequence.decode("ascii")
            except UnicodeDecodeError as error:
                raise GSE105127IndependentAuditError("reference FASTA is not ASCII") from error
            alphabet.update(text)
            length += len(text)
    if name is not None:
        inventory.append((name, length))
    return inventory, dict(sorted(alphabet.items()))


def validate_fasta_inventory(source_fa: Path, source_fai: Path) -> dict[str, Any]:
    inventory, alphabet = fasta_inventory(source_fa)
    fields, _ = (), []
    fai_rows = []
    with source_fai.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 5:
                raise GSE105127IndependentAuditError("source FASTA index row differs")
            fai_rows.append((row[0], int(row[1])))
    names = [name for name, _ in inventory]
    if (
        len(inventory) != 84
        or len(set(names)) != 84
        or tuple(names[:25]) != PRIMARY_SOURCE
        or inventory != fai_rows
        or alphabet != EXPECTED_ALPHABET_COUNTS
    ):
        raise GSE105127IndependentAuditError("admitted decompressed FASTA census differs")
    return {"contigs": 84, "primary_contigs": 25, "non_primary_contigs": 59, "alphabet": alphabet}


def validate_topology(
    rna_fields: tuple[str, ...],
    rna_rows: Sequence[Mapping[str, str]],
    rrbs_fields: tuple[str, ...],
    rrbs_rows: Sequence[Mapping[str, str]],
    bundle_fields: tuple[str, ...],
    bundle_rows: Sequence[Mapping[str, str]],
) -> dict[str, Any]:
    if rna_fields != RNA_FIELDS or rrbs_fields != RRBS_FIELDS or bundle_fields != BUNDLE_FIELDS:
        raise GSE105127IndependentAuditError("source plan schema differs")
    if len(rna_rows) != EXPECTED_ROWS or len(rrbs_rows) != EXPECTED_ROWS or len(bundle_rows) != EXPECTED_BUNDLES:
        raise GSE105127IndependentAuditError("source plan census differs")
    rna_by_row = {row["row_id"]: row for row in rna_rows}
    rrbs_by_row = {row["row_id"]: row for row in rrbs_rows}
    if len(rna_by_row) != EXPECTED_ROWS or len(rrbs_by_row) != EXPECTED_ROWS or set(rna_by_row) != set(rrbs_by_row):
        raise GSE105127IndependentAuditError("RNA/RRBS row identity differs")
    participant_zones: defaultdict[str, set[str]] = defaultdict(set)
    participant_bundles: defaultdict[str, set[str]] = defaultdict(set)
    bundle_summary: defaultdict[str, dict[str, int]] = defaultdict(
        lambda: {"participants": 0, "participant_zone_rows": 0, "combined_source_bytes": 0}
    )
    participants_by_bundle: defaultdict[str, set[str]] = defaultdict(set)
    for row_id in sorted(rna_by_row):
        rna = rna_by_row[row_id]
        rrbs = rrbs_by_row[row_id]
        identity = ("participant_group_id", "zone", "bundle_id", "pairing_topology")
        if any(rna[key] != rrbs[key] for key in identity):
            raise GSE105127IndependentAuditError("RNA/RRBS participant-zone topology differs")
        if (
            rna["pairing_topology"] != "adjacent_section"
            or rna["zone"] not in {"CV", "IZ", "PP"}
            or rna["library_layout"] != "SINGLE"
            or int(rna["read_length"]) != 76
            or int(rna["base_count"]) != 76 * int(rna["read_count"])
            or re.fullmatch(r"[0-9a-f]{32}", rna["fastq_md5"]) is None
            or int(rna["fastq_bytes"]) <= 0
            or int(rrbs["source_bed_bytes"]) <= 0
            or int(rrbs["source_cytosine_rows"]) <= 0
            or re.fullmatch(r"[0-9a-f]{64}", rrbs["source_bed_sha256"]) is None
        ):
            raise GSE105127IndependentAuditError("participant-zone source contract differs")
        participant = rna["participant_group_id"]
        bundle = rna["bundle_id"]
        participant_zones[participant].add(rna["zone"])
        participant_bundles[participant].add(bundle)
        participants_by_bundle[bundle].add(participant)
        bundle_summary[bundle]["participant_zone_rows"] += 1
        bundle_summary[bundle]["combined_source_bytes"] += int(rna["fastq_bytes"]) + int(rrbs["source_bed_bytes"])
    if (
        len(participant_zones) != EXPECTED_PARTICIPANTS
        or any(zones != {"CV", "IZ", "PP"} for zones in participant_zones.values())
        or any(len(bundles) != 1 for bundles in participant_bundles.values())
    ):
        raise GSE105127IndependentAuditError("participant containment or zone completeness differs")
    for bundle, participants in participants_by_bundle.items():
        bundle_summary[bundle]["participants"] = len(participants)
    expected_summary = {
        row["bundle_id"]: {
            "participants": int(row["participants"]),
            "participant_zone_rows": int(row["participant_zone_rows"]),
            "combined_source_bytes": int(row["combined_source_bytes"]),
        }
        for row in bundle_rows
    }
    if dict(bundle_summary) != expected_summary or set(expected_summary) != {str(value) for value in range(8)}:
        raise GSE105127IndependentAuditError("bundle accounting differs")
    if (
        sum(int(row["fastq_bytes"]) for row in rna_rows) != EXPECTED_RNA_BYTES
        or sum(int(row["read_count"]) for row in rna_rows) != EXPECTED_RNA_READS
        or sum(int(row["base_count"]) for row in rna_rows) != EXPECTED_RNA_BASES
        or sum(int(row["source_bed_bytes"]) for row in rrbs_rows) != EXPECTED_RRBS_BYTES
    ):
        raise GSE105127IndependentAuditError("source plan totals differ")
    return {"participants": 19, "rows": 57, "zones_per_participant": 3, "pairing_topology": "adjacent_section"}


def audit_python_surface(path: Path, expected_sha256: str) -> None:
    if sha256_file(path) != expected_sha256:
        raise GSE105127IndependentAuditError(f"source code hash differs: {path}")
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=path.as_posix())
    calls = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Attribute):
                calls.add(node.func.attr)
            elif isinstance(node.func, ast.Name):
                calls.add(node.func.id)
    if calls & FORBIDDEN_CALLS:
        raise GSE105127IndependentAuditError(f"modeling call entered source code: {path}")
    lowered = source.lower()
    for forbidden in ("evaluator_only", "outcomes.tsv", "predictions"):
        if forbidden in lowered:
            raise GSE105127IndependentAuditError(f"forbidden data path entered source code: {path}")


def audit_shell_surface(path: Path, expected_sha256: str) -> None:
    if sha256_file(path) != expected_sha256:
        raise GSE105127IndependentAuditError(f"wrapper hash differs: {path}")
    lowered = path.read_text(encoding="utf-8").lower()
    for forbidden in ("evaluator_only", "outcomes.tsv", "predictions", "config/evaluation"):
        if forbidden in lowered:
            raise GSE105127IndependentAuditError(f"forbidden data path entered wrapper: {path}")


def validate_submission_record(path: Path) -> str:
    fields, rows = read_tsv(path)
    expected_fields = ("JobID", "JobName", "State", "ExitCode", "SubmitLine")
    if fields != expected_fields:
        raise GSE105127IndependentAuditError("SLURM record schema differs")
    by_id = {row["JobID"].split("_", 1)[0]: row for row in rows}
    required = {
        "21083518": ("afterok:21083517", "gse105127_stage_plan_cpu.sbatch"),
        "21083520": ("afterok:21083518", "gse105127_stage_rrbs_io.sbatch"),
        "21083522": ("afterok:21083518", "gse105127_stage_rsem_reference_cpu.sbatch"),
        "21087896": ("afterany:21083521", "gse105127_stage_rna_download_resumable_io.sbatch"),
        "21087897": ("afterok:21087896", "gse105127_stage_rna_quant_cpu.sbatch"),
        "21087898": ("afterok:21087897", "gse105127_stage_consolidate_cpu.sbatch"),
        "21088014": ("GSE105127_CAMPAIGN_ID=9675cd48993468d4", "gse105127_stage_reference_legacy_exact_io.sbatch"),
        "21088066": ("afterok:21088014", "gse105127_stage_crosswalk_io.sbatch"),
        "21088068": ("afterok:21088066:21087898", "gse105127_stage_finalize_cpu.sbatch"),
    }
    for job_id, substrings in required.items():
        row = by_id.get(job_id)
        if row is None:
            raise GSE105127IndependentAuditError(f"SLURM record lacks job {job_id}")
        if any(value not in row["SubmitLine"] for value in substrings):
            raise GSE105127IndependentAuditError(f"SLURM dependency differs for job {job_id}")
        if row["JobName"] == "" or row["ExitCode"] == "":
            raise GSE105127IndependentAuditError(f"SLURM accounting differs for job {job_id}")
    return sha256_file(path)


def validate_source_stage_receipts(
    plan_rna: Sequence[Mapping[str, str]],
    plan_rrbs: Sequence[Mapping[str, str]],
    rna_raw_root: Path,
    rrbs_root: Path,
) -> dict[str, Any]:
    rna_by_row = {row["row_id"]: row for row in plan_rna}
    rrbs_by_row = {row["row_id"]: row for row in plan_rrbs}
    rna_receipts: dict[str, dict[str, Any]] = {}
    rrbs_receipts: dict[str, dict[str, Any]] = {}
    for bundle_id in range(8):
        for root, expected_class in (
            (rna_raw_root / f"bundle_{bundle_id:02d}", "gse105127_rna_download_bundle"),
            (rrbs_root / f"bundle_{bundle_id:02d}", "gse105127_rrbs_collapse_bundle"),
        ):
            manifest = verify_frozen_tree(root)
            receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
            if (
                manifest.get("metadata", {}).get("status") != "passed"
                or receipt.get("status") != "complete"
                or receipt.get("labels_accessed") is not False
                or receipt.get("fit_or_score_performed") is not False
                or receipt.get("bundle_id") != bundle_id
            ):
                raise GSE105127IndependentAuditError("source bundle receipt differs")
    for row_id, plan in rna_by_row.items():
        participant = rna_raw_root / "participants" / row_id
        verify_frozen_tree(participant)
        receipt = json.loads((participant / "receipt.json").read_text(encoding="utf-8"))
        fastq = participant / "reads.fastq.gz"
        if (
            receipt.get("status") != "complete_verified"
            or receipt.get("row_id") != row_id
            or receipt.get("participant_group_id") != plan["participant_group_id"]
            or receipt.get("zone") != plan["zone"]
            or receipt.get("pairing_topology") != "adjacent_section"
            or receipt.get("read_length") != 76
            or receipt.get("read_count") != int(plan["read_count"])
            or receipt.get("base_count") != int(plan["base_count"])
            or receipt.get("bytes") != int(plan["fastq_bytes"])
            or receipt.get("md5") != plan["fastq_md5"]
            or receipt.get("gzip_integrity") is not True
            or receipt.get("labels_accessed") is not False
            or receipt.get("fit_or_score_performed") is not False
            or gzip_test(fastq)[0] != 0
        ):
            raise GSE105127IndependentAuditError("RNA source participant receipt differs")
        rna_receipts[row_id] = receipt
    for row_id, plan in rrbs_by_row.items():
        participant = rrbs_root / "participants" / row_id
        verify_frozen_tree(participant)
        receipt = json.loads((participant / "receipt.json").read_text(encoding="utf-8"))
        collapsed = participant / "cpg_counts.hg19.tsv.gz"
        if (
            receipt.get("status") != "passed"
            or receipt.get("row_id") != row_id
            or receipt.get("participant_group_id") != plan["participant_group_id"]
            or receipt.get("zone") != plan["zone"]
            or receipt.get("pairing_topology") != "adjacent_section"
            or receipt.get("source_bed_sha256") != plan["source_bed_sha256"]
            or receipt.get("labels_accessed") is not False
            or receipt.get("fit_or_score_performed") is not False
            or gzip_test(collapsed)[0] != 0
        ):
            raise GSE105127IndependentAuditError("RRBS source participant receipt differs")
        rrbs_receipts[row_id] = receipt
    if set(rna_receipts) != set(rna_by_row) or set(rrbs_receipts) != set(rrbs_by_row):
        raise GSE105127IndependentAuditError("source participant receipt census differs")
    return {"rna_participants": 57, "rrbs_participants": 57, "rna_bundles": 8, "rrbs_bundles": 8}


def run(
    *,
    project_root: Path,
    campaign_root: Path,
    legacy_diagnostic_root: Path,
    alphabet_diagnostic_root: Path,
    source_cache_root: Path,
    target_fasta: Path,
    target_gtf: Path,
    submission_record: Path,
    output: Path,
) -> None:
    if output.exists():
        raise GSE105127IndependentAuditError("refusing to overwrite independent audit")
    control = campaign_root / "control"
    plan = campaign_root / "plan"
    reference = campaign_root / "reference"
    rsem_reference = campaign_root / "rsem_reference"
    rna_raw = campaign_root / "rna_raw"
    rrbs = campaign_root / "rrbs_collapsed"
    require_pinned_tree(control, CONTROL_SHA256, "campaign control")
    require_pinned_tree(plan, PLAN_SHA256, "production plan")
    require_pinned_tree(reference, REFERENCE_SHA256, "legacy reference")
    require_pinned_tree(rsem_reference, RSEM_REFERENCE_SHA256, "RSEM reference")
    require_pinned_tree(legacy_diagnostic_root, LEGACY_DIAGNOSTIC_SHA256, "legacy diagnostic")
    require_pinned_tree(alphabet_diagnostic_root, ALPHABET_DIAGNOSTIC_SHA256, "alphabet diagnostic")
    require_pinned_tree(source_cache_root, SOURCE_CACHE_SHA256, "reference source cache")

    plan_receipt = json.loads((plan / "plan.json").read_text(encoding="utf-8"))
    if (
        plan_receipt.get("labels_accessed") is not False
        or plan_receipt.get("fit_or_score_performed") is not False
        or plan_receipt.get("historical_outcome_qc_table_opened") is not False
        or plan_receipt.get("participants") != 19
        or plan_receipt.get("participant_zone_rows") != 57
    ):
        raise GSE105127IndependentAuditError("production plan receipt differs")
    rna_fields, rna_rows = read_tsv(plan / "rna_rows.tsv")
    rrbs_fields, rrbs_rows = read_tsv(plan / "rrbs_rows.tsv")
    bundle_fields, bundle_rows = read_tsv(plan / "bundles.tsv")
    topology = validate_topology(rna_fields, rna_rows, rrbs_fields, rrbs_rows, bundle_fields, bundle_rows)

    legacy_rule = validate_legacy_reference_rule(
        reference / "human_g1k_v37.fasta.gz", reference / "human_g1k_v37.fasta"
    )
    fasta_census = validate_fasta_inventory(
        reference / "human_g1k_v37.fasta", reference / "human_g1k_v37.fasta.fai"
    )
    reference_receipt = json.loads((reference / "receipt.json").read_text(encoding="utf-8"))
    if (
        reference_receipt.get("source_compressed_sha256") != SOURCE_COMPRESSED_SHA256
        or reference_receipt.get("source_uncompressed_sha256") != SOURCE_UNCOMPRESSED_SHA256
        or reference_receipt.get("source_contigs_total") != 84
        or reference_receipt.get("labels_accessed") is not False
        or reference_receipt.get("fit_or_score_performed") is not False
    ):
        raise GSE105127IndependentAuditError("legacy reference receipt differs")
    for chain, expected_md5 in (
        (reference / "hg19ToHg38.over.chain.gz", FORWARD_CHAIN_MD5),
        (reference / "hg38ToHg19.over.chain.gz", REVERSE_CHAIN_MD5),
    ):
        if gzip_test(chain)[0] != 0 or digest_file(chain, "md5") != expected_md5:
            raise GSE105127IndependentAuditError("lift chain integrity differs")
    for compressed, expected_sha in ((target_fasta, TARGET_FASTA_SHA256), (target_gtf, TARGET_GTF_SHA256)):
        if gzip_test(compressed)[0] != 0 or sha256_file(compressed) != expected_sha:
            raise GSE105127IndependentAuditError("target reference integrity differs")
    rsem_receipt = json.loads((rsem_reference / "receipt.json").read_text(encoding="utf-8"))
    if (
        rsem_receipt.get("status") != "passed"
        or rsem_receipt.get("fasta_sha256") != TARGET_FASTA_SHA256
        or rsem_receipt.get("gtf_sha256") != TARGET_GTF_SHA256
        or rsem_receipt.get("labels_accessed") is not False
        or rsem_receipt.get("fit_or_score_performed") is not False
    ):
        raise GSE105127IndependentAuditError("RSEM reference receipt differs")

    source_stages = validate_source_stage_receipts(rna_rows, rrbs_rows, rna_raw, rrbs)
    for relative, expected in SOURCE_CODE_HASHES.items():
        audit_python_surface(project_root / relative, expected)
    for relative, expected in SHELL_CODE_HASHES.items():
        audit_shell_surface(project_root / relative, expected)
    submission_sha = validate_submission_record(submission_record)

    output.mkdir(mode=0o750)
    receipt = {
        "schema_version": "masld-bench-gse105127-independent-source-dag-audit-v1",
        "status": "passed_audit_source_mechanics_verified_promotion_blocked",
        "source_mechanics_independently_rederived": True,
        "participants": 19,
        "participant_zone_rows": 57,
        "topology": topology,
        "source_stages": source_stages,
        "reference_fasta_census": fasta_census,
        "legacy_compressed_source": legacy_rule,
        "gzip_t_rule_resolution": (
            "compressed legacy blob fails gzip -t and is provenance-only; only the exact "
            "decompressed FASTA is admitted after compressed identity, one-member CRC/ISIZE, "
            "uncompressed identity, alphabet, contig, and faidx checks"
        ),
        "strict_every_sequence_job_begins_with_gzip_t_contract_passed": False,
        "strict_gzip_t_procedural_blockers": [
            "legacy reference materialization intentionally consumes a compressed blob whose gzip -t exits 2",
            "RSEM reference construction checks hashes but does not run gzip -t before gzip.open",
            "RRBS collapse uses gzip -dc without a preceding gzip -t",
            "RNA quantification consumes frozen FASTQ gzip files without a job-local preceding gzip -t",
            "CpG crosswalk consumes collapsed gzip files and a compressed lift chain without a job-local preceding gzip -t",
        ],
        "required_campaign_revision": (
            "add explicit fail-closed gzip -t preflights before every valid compressed input; "
            "keep the exact legacy compressed FASTA provenance-only and allow only its independently "
            "authenticated decompressed derivative downstream"
        ),
        "current_dag_promotion_eligible": False,
        "compressed_legacy_source_usable_by_sequence_jobs": False,
        "decompressed_legacy_reference_eligible_for_lift_mechanics": True,
        "rna_or_methylation_values_parsed": False,
        "labels_read": False,
        "outcomes_read": False,
        "predictions_read": False,
        "fit_or_score_performed": False,
        "modeling_or_scoring_authorized_by_this_receipt": False,
        "sealed_evaluation_authorized_by_this_receipt": False,
        "clean_or_sealed_champion_eligible_by_this_receipt": False,
        "rna_quantification_values_audited": False,
        "methylation_values_audited": False,
        "submission_record_sha256": submission_sha,
        "control_artifacts_sha256": CONTROL_SHA256,
        "plan_artifacts_sha256": PLAN_SHA256,
        "reference_artifacts_sha256": REFERENCE_SHA256,
        "rsem_reference_artifacts_sha256": RSEM_REFERENCE_SHA256,
    }
    write_json_exclusive(output / "independent_source_dag_audit.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse105127_independent_source_dag_audit",
            "source_mechanics_independently_rederived": True,
            "compressed_legacy_source_usable_by_sequence_jobs": False,
            "decompressed_legacy_reference_eligible_for_lift_mechanics": True,
            "strict_every_sequence_job_begins_with_gzip_t_contract_passed": False,
            "current_dag_promotion_eligible": False,
            "modeling_or_scoring_authorized": False,
            "clean_or_sealed_champion_eligible": False,
            "status": "passed_audit_source_mechanics_verified_promotion_blocked",
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--campaign-root", required=True, type=Path)
    parser.add_argument("--legacy-diagnostic-root", required=True, type=Path)
    parser.add_argument("--alphabet-diagnostic-root", required=True, type=Path)
    parser.add_argument("--source-cache-root", required=True, type=Path)
    parser.add_argument("--target-fasta", required=True, type=Path)
    parser.add_argument("--target-gtf", required=True, type=Path)
    parser.add_argument("--submission-record", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(**vars(args))
    print(json.dumps({"output": args.output.as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
