#!/usr/bin/env python3
"""Run one outcome-blind LS-GKM production-command scale probe."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import subprocess
import time
from typing import Any, Iterable, Mapping, Sequence

from scripts.lsgkm_safe_adapter import (
    canonical_kmer,
    canonical_score_table,
    deltasvm_scores,
    direct_gkmsvm_alt_minus_ref,
    file_sha256,
    parse_predictions,
    validate_fasta,
    validate_project_model,
)


SCHEMA = "masld-bench-lsgkm-gse281364-dinucleotide-scale-probe-input-v1"
SCORE_FIELDS = (
    "split_id",
    "model_seed",
    "model_state_id",
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
    "outer_fold",
    "contig",
    "variant_pos0",
    "ref",
    "alt",
    "direct_gkmsvm_alt_minus_ref",
    "deltasvm_alt_minus_ref",
)


class LSGKMScaleProbeError(RuntimeError):
    """Raised when the one-fit scale-probe requirements are not met."""


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise LSGKMScaleProbeError(f"JSON object differs: {path}")
    return value


def resolve_file(root: Path, relative_text: str, expected_sha256: str) -> Path:
    relative = Path(relative_text)
    if relative.is_absolute() or ".." in relative.parts:
        raise LSGKMScaleProbeError("unsafe project-relative file")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    if path.is_symlink() or not path.is_file() or file_sha256(path) != expected_sha256:
        raise LSGKMScaleProbeError(f"file authority differs: {relative_text}")
    return path


def verify_manifest_members(
    root: Path,
    artifacts_path: str,
    artifacts_sha256: str,
    bindings: Iterable[tuple[str, str]],
) -> None:
    artifacts = resolve_file(root, artifacts_path, artifacts_sha256)
    manifest = load_json(artifacts)
    roster = {
        str(member["path"]): str(member["sha256"])
        for member in manifest.get("artifacts", [])
        if isinstance(member, dict)
    }
    for relative_text, expected_sha256 in bindings:
        member = resolve_file(root, relative_text, expected_sha256)
        relative = member.relative_to(artifacts.parent).as_posix()
        if roster.get(relative) != expected_sha256:
            raise LSGKMScaleProbeError(f"member absent from authority: {relative}")


def read_tsv_gzip(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    if not rows:
        raise LSGKMScaleProbeError(f"TSV is empty: {path}")
    return rows


def read_fasta_gzip(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    sequence: list[str] = []
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(sequence)))
                identifier, sequence = line[1:], []
            elif identifier is None or not line:
                raise LSGKMScaleProbeError("gzip FASTA structure differs")
            else:
                sequence.append(line)
    if identifier is not None:
        records.append((identifier, "".join(sequence)))
    if not records:
        raise LSGKMScaleProbeError("gzip FASTA is empty")
    return records


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> int:
    if path.exists():
        raise LSGKMScaleProbeError(f"refusing to overwrite: {path}")
    count = 0
    with path.open("x", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n{sequence}\n")
            count += 1
    return count


def write_tsv_gzip(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> int:
    if path.exists():
        raise LSGKMScaleProbeError(f"refusing to overwrite: {path}")
    count = 0
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                for row in rows:
                    writer.writerow({field: row[field] for field in fields})
                    count += 1
    return count


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != "one_fit_outcome_blind_scale_probe_authorized"
        or config.get("design_id")
        != "lsgkm_hepatocyte_atac_peak_exact_dinucleotide_null_v1"
    ):
        raise LSGKMScaleProbeError("scale probe config identity differs")
    runtime = config["runtime"]
    if (
        runtime.get("gkmtrain_argv")
        != ["-t", "2", "-l", "11", "-k", "7", "-d", "3", "-c", "1", "-e", "0.001", "-w", "1", "-m", "4096", "-T", "1", "-z"]
        or runtime.get("gkmpredict_argv") != ["-T", "1"]
        or runtime.get("native_fit_seed_available") is not False
        or runtime.get("resume_supported") is not False
    ):
        raise LSGKMScaleProbeError("scale probe runtime differs")
    fit = config["active_fit"]
    if (
        fit.get("split_id") != "donor0_genomic0"
        or fit.get("model_seed") != 1103
        or fit.get("positive_count") != 10000
        or fit.get("negative_count") != 10000
        or fit.get("sequence_length_bp") != 300
        or fit.get("positive_class_weight") != 1
        or fit.get("model_state_id") != "hepatocyte"
    ):
        raise LSGKMScaleProbeError("active fit differs")
    firewall = config["action_firewall"]
    if (
        firewall.get("one_fit_production_training_authorized") is not True
        or firewall.get("held_sequence_scoring_authorized") is not True
        or firewall.get("remaining_24_fits_authorized") is not False
        or firewall.get("benchmark_metric_calculation_authorized") is not False
        or any(
            firewall.get(field) is not False
            for field in (
                "outcome_access_authorized",
                "reporter_count_access_authorized",
                "sealed_asset_access_authorized",
                "control_prediction_value_access_authorized",
            )
        )
    ):
        raise LSGKMScaleProbeError("scale probe action firewall differs")


def run_timed(
    argv: Sequence[str], *, stdout_path: Path, stderr_path: Path, time_path: Path
) -> dict[str, Any]:
    started = time.monotonic()
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        completed = subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(time_path), *argv],
            stdout=stdout,
            stderr=stderr,
            check=False,
        )
    elapsed = time.monotonic() - started
    if completed.returncode != 0:
        raise LSGKMScaleProbeError(f"native command failed: {argv[0]}")
    time_text = time_path.read_text(encoding="utf-8")
    match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", time_text)
    if match is None:
        raise LSGKMScaleProbeError("GNU time MaxRSS receipt is missing")
    return {
        "argv": list(argv),
        "returncode": completed.returncode,
        "wall_seconds_monotonic": elapsed,
        "maximum_resident_set_kbytes": int(match.group(1)),
        "stdout_sha256": file_sha256(stdout_path),
        "stderr_sha256": file_sha256(stderr_path),
        "time_receipt_sha256": file_sha256(time_path),
    }


def run_probe(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    config = load_json(config_path)
    validate_config(config)
    if output.exists():
        raise LSGKMScaleProbeError("scale probe output exists")
    for directory in ("inputs", "models", "predictions", "logs"):
        (output / directory).mkdir(parents=True, exist_ok=True, mode=0o750)

    runtime = config["runtime"]
    verify_manifest_members(
        root,
        runtime["artifacts_path"],
        runtime["artifacts_sha256"],
        (
            (runtime["gkmtrain_path"], runtime["gkmtrain_sha256"]),
            (runtime["gkmpredict_path"], runtime["gkmpredict_sha256"]),
        ),
    )
    inputs = config["inputs"]
    input_bindings = (
        (inputs["independent_audit_path"], inputs["independent_audit_sha256"]),
        (inputs["positive_path"], inputs["positive_sha256"]),
        (inputs["negative_path"], inputs["negative_sha256"]),
        (inputs["pair_manifest_path"], inputs["pair_manifest_sha256"]),
        (inputs["scoring_manifest_path"], inputs["scoring_manifest_sha256"]),
        (inputs["scoring_fasta_path"], inputs["scoring_fasta_sha256"]),
    )
    verify_manifest_members(
        root, inputs["artifacts_path"], inputs["artifacts_sha256"], input_bindings
    )
    audit = load_json(root / inputs["independent_audit_path"])
    if (
        audit.get("status") != "pass_independent_emitted_fasta_invariant_audit"
        or audit.get("fit_count") != 25
        or audit.get("audited_pair_rows") != 250000
        or audit.get("outcomes_read") is not False
    ):
        raise LSGKMScaleProbeError("input audit differs")

    positive_records = read_fasta_gzip(root / inputs["positive_path"])
    negative_records = read_fasta_gzip(root / inputs["negative_path"])
    pairs = read_tsv_gzip(root / inputs["pair_manifest_path"])
    if len(positive_records) != 10000 or len(negative_records) != 10000 or len(pairs) != 10000:
        raise LSGKMScaleProbeError("fit input denominator differs")
    positive_path = output / "inputs/positive.fa"
    negative_path = output / "inputs/negative.fa"
    write_fasta(positive_path, positive_records)
    write_fasta(negative_path, negative_records)
    validated_positive = validate_fasta(positive_path, expected_length=300)
    validated_negative = validate_fasta(negative_path, expected_length=300)
    if len(validated_positive) != 10000 or len(validated_negative) != 10000:
        raise LSGKMScaleProbeError("safe-adapter fit input differs")

    fit = config["active_fit"]
    model_prefix = output / "models/fit"
    train_argv = [
        str(root / runtime["gkmtrain_path"]),
        *runtime["gkmtrain_argv"],
        str(positive_path),
        str(negative_path),
        str(model_prefix),
    ]
    train_receipt = run_timed(
        train_argv,
        stdout_path=output / "logs/gkmtrain.stdout.txt",
        stderr_path=output / "logs/gkmtrain.stderr.txt",
        time_path=output / "logs/gkmtrain.time.txt",
    )
    model_path = output / "models/fit.model.txt.gz"
    model_receipt = validate_project_model(model_path, support_vector_length=300)

    scoring_rows = read_tsv_gzip(root / inputs["scoring_manifest_path"])
    scoring_records = dict(read_fasta_gzip(root / inputs["scoring_fasta_path"]))
    held_rows = [row for row in scoring_rows if row["scoring_split_id"] == fit["split_id"]]
    if not held_rows or any(row["scoring_split_id"] != fit["split_id"] for row in held_rows):
        raise LSGKMScaleProbeError("held scoring roster differs")
    held_fasta_records: list[tuple[str, str]] = []
    for row in held_rows:
        element = row["element_id"]
        held_fasta_records.extend(
            (
                (f"{element}|REF", scoring_records[f"{element}|REF"]),
                (f"{element}|ALT", scoring_records[f"{element}|ALT"]),
            )
        )
    held_fasta_path = output / "inputs/held_alleles.fa"
    write_fasta(held_fasta_path, held_fasta_records)
    validate_fasta(held_fasta_path, expected_length=300)
    direct_path = output / "predictions/held_allele_native.tsv"
    direct_receipt = run_timed(
        [
            str(root / runtime["gkmpredict_path"]),
            *runtime["gkmpredict_argv"],
            str(held_fasta_path),
            str(model_path),
            str(direct_path),
        ],
        stdout_path=output / "logs/gkmpredict_alleles.stdout.txt",
        stderr_path=output / "logs/gkmpredict_alleles.stderr.txt",
        time_path=output / "logs/gkmpredict_alleles.time.txt",
    )
    direct = parse_predictions(
        direct_path, expected_identifiers=[identifier for identifier, _sequence in held_fasta_records]
    )

    required_kmers: set[str] = set()
    for row in held_rows:
        element = row["element_id"]
        for allele in ("REF", "ALT"):
            sequence = scoring_records[f"{element}|{allele}"][140:161]
            for start in range(11):
                required_kmers.add(canonical_kmer(sequence[start : start + 11]))
    kmer_records = [
        (f"kmer_{hashlib.sha256(sequence.encode('ascii')).hexdigest()}", sequence)
        for sequence in sorted(required_kmers)
    ]
    kmer_path = output / "inputs/required_canonical_11mers.fa"
    write_fasta(kmer_path, kmer_records)
    validate_fasta(kmer_path, expected_length=11)
    kmer_score_path = output / "predictions/required_canonical_11mers.native.tsv"
    kmer_receipt = run_timed(
        [
            str(root / runtime["gkmpredict_path"]),
            *runtime["gkmpredict_argv"],
            str(kmer_path),
            str(model_path),
            str(kmer_score_path),
        ],
        stdout_path=output / "logs/gkmpredict_11mers.stdout.txt",
        stderr_path=output / "logs/gkmpredict_11mers.stderr.txt",
        time_path=output / "logs/gkmpredict_11mers.time.txt",
    )
    kmer_native = parse_predictions(
        kmer_score_path,
        expected_identifiers=[identifier for identifier, _sequence in kmer_records],
    )
    table = canonical_score_table(
        (sequence, kmer_native[identifier]) for identifier, sequence in kmer_records
    )

    score_rows: list[dict[str, object]] = []
    for row in held_rows:
        element = row["element_id"]
        reference = scoring_records[f"{element}|REF"]
        alternative = scoring_records[f"{element}|ALT"]
        direct_score = direct_gkmsvm_alt_minus_ref(
            reference_score=direct[f"{element}|REF"],
            alternative_score=direct[f"{element}|ALT"],
        )
        delta = deltasvm_scores(reference[140:161], alternative[140:161], table)[
            "canonical_alt_minus_ref"
        ]
        if not math.isfinite(direct_score) or not math.isfinite(delta):
            raise LSGKMScaleProbeError("non-finite scale-probe score")
        score_rows.append(
            {
                "split_id": fit["split_id"],
                "model_seed": fit["model_seed"],
                "model_state_id": fit["model_state_id"],
                "element_id": element,
                "source_locus_group_id": row["source_locus_group_id"],
                "long_range_block_id": row["long_range_block_id"],
                "outer_fold": row["outer_fold"],
                "contig": row["contig"],
                "variant_pos0": row["variant_pos0"],
                "ref": row["ref"],
                "alt": row["alt"],
                "direct_gkmsvm_alt_minus_ref": f"{direct_score:.17g}",
                "deltasvm_alt_minus_ref": f"{delta:.17g}",
            }
        )
    score_path = output / "predictions/held_element_scores.tsv.gz"
    write_tsv_gzip(score_path, SCORE_FIELDS, score_rows)
    command_receipts = {
        "gkmtrain": train_receipt,
        "gkmpredict_held_alleles": direct_receipt,
        "gkmpredict_required_11mers": kmer_receipt,
    }
    receipt = {
        "schema_version": "masld-bench-lsgkm-gse281364-dinucleotide-scale-probe-v1",
        "status": "pass_one_fit_production_command_scale_probe",
        "design_id": config["design_id"],
        "config_sha256": file_sha256(config_path),
        "runtime_artifacts_sha256": runtime["artifacts_sha256"],
        "input_artifacts_sha256": inputs["artifacts_sha256"],
        "active_fit": fit,
        "native_fit_seed_available": False,
        "model_seed_role": "negative_generation_identity_only",
        "positive_records": len(positive_records),
        "negative_records": len(negative_records),
        "held_elements_scored": len(score_rows),
        "held_REF_ALT_sequences_scored": len(held_fasta_records),
        "required_canonical_11mers_scored": len(kmer_records),
        "model_receipt": model_receipt,
        "command_receipts": command_receipts,
        "maximum_resident_set_kbytes": max(
            value["maximum_resident_set_kbytes"] for value in command_receipts.values()
        ),
        "total_native_wall_seconds": sum(
            value["wall_seconds_monotonic"] for value in command_receipts.values()
        ),
        "held_element_scores_sha256": file_sha256(score_path),
        "readouts": ["direct_gkmsvm", "deltasvm"],
        "readouts_are_independent_model_families": False,
        "production_fits_executed": 1,
        "production_predictions_generated": True,
        "benchmark_metrics_calculated": False,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "control_prediction_values_read": False,
        "sealed_assets_read": False,
        "claim_limits": config["claim_limits"],
    }
    (output / "scale_probe_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    run_probe(
        arguments.project_root.resolve(strict=True),
        arguments.config.resolve(strict=True),
        arguments.output,
    )


if __name__ == "__main__":
    main()
