#!/usr/bin/env python3
"""Build and audit a synthetic, outcome-blind LS-GKM runtime fixture."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import tomllib
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import verify_frozen_tree
from scripts.lsgkm_safe_adapter import (
    LSGKMSafetyError,
    canonical_kmer,
    deltasvm_scores,
    direct_gkmsvm_alt_minus_ref,
    file_sha256,
    parse_predictions,
    reverse_complement,
    validate_fasta,
    validate_project_model,
)


SCHEMA = "masld-bench-lsgkm-runtime-fixture-v1"


class RuntimeFixtureError(RuntimeError):
    """Raised when source, runtime, fixture, or activation requirement differs."""


def load_config(path: Path) -> dict[str, Any]:
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, tomllib.TOMLDecodeError) as error:
        raise RuntimeFixtureError(f"invalid runtime config: {error}") from error
    if not isinstance(value, dict) or value.get("schema_version") != SCHEMA:
        raise RuntimeFixtureError("runtime config schema differs")
    return value


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise RuntimeFixtureError(f"invalid JSON authority: {path}") from error
    if not isinstance(value, dict):
        raise RuntimeFixtureError(f"JSON authority is not an object: {path}")
    return value


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_fasta(path: Path, records: Sequence[tuple[str, str]]) -> None:
    with path.open("x", encoding="ascii", newline="") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n{sequence}\n")


def random_sequence(generator: random.Random, length: int) -> str:
    return "".join(generator.choice("ACGT") for _ in range(length))


def insert_motif(sequence: str, motif: str, positions: Sequence[int]) -> str:
    values = list(sequence)
    for start in positions:
        if start < 0 or start + len(motif) > len(values):
            raise RuntimeFixtureError("fixture motif position is out of bounds")
        values[start : start + len(motif)] = motif
    return "".join(values)


def command_receipt(
    command: Sequence[str], *, stdout_path: Path, stderr_path: Path
) -> dict[str, Any]:
    with stdout_path.open("x", encoding="utf-8") as stdout, stderr_path.open(
        "x", encoding="utf-8"
    ) as stderr:
        completed = subprocess.run(
            list(command),
            stdin=subprocess.DEVNULL,
            stdout=stdout,
            stderr=stderr,
            check=False,
            shell=False,
        )
    if completed.returncode != 0:
        raise RuntimeFixtureError(f"native command failed: {' '.join(command)}")
    return {
        "argv": list(command),
        "returncode": completed.returncode,
        "stdout_sha256": file_sha256(stdout_path),
        "stderr_sha256": file_sha256(stderr_path),
    }


def validate_authorities(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    tree_contracts: dict[str, dict[str, Any]] = {}
    file_contracts: dict[str, dict[str, Any]] = {}
    for name, authority in config["authorities"].items():
        relative = Path(authority["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise RuntimeFixtureError("unsafe authority path")
        path = (root / relative).resolve(strict=True)
        path.relative_to(root)
        if path.is_dir():
            expected = authority["artifacts_sha256"]
            observed = file_sha256(path / "ARTIFACTS.json")
            if observed != expected:
                raise RuntimeFixtureError(f"frozen tree differs: {name}")
            verify_frozen_tree(path)
            tree_contracts[name] = {
                "path": relative.as_posix(),
                "artifacts_sha256": observed,
            }
        else:
            observed = file_sha256(path)
            if observed != authority["sha256"]:
                raise RuntimeFixtureError(f"file authority differs: {name}")
            file_contracts[name] = {
                "path": relative.as_posix(),
                "sha256": observed,
            }

    pseudobulk = root / config["authorities"]["gse296875_training_pseudobulk"]["path"]
    pseudobulk_contract = load_json(pseudobulk / "contract.json")
    observed_roster = pseudobulk_contract.get("lineages")
    split_tree = root / config["authorities"]["sequence_split"]["path"]
    split_contract = load_json(split_tree / "contract.json")
    with (split_tree / "crossed_outer_splits.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        split_rows = list(csv.DictReader(handle, delimiter="\t"))
    activation = config["activation"]
    if (
        observed_roster != activation["observed_frozen_gse296875_roster"]
        or observed_roster == activation["registered_state_roster"]
        or len(split_rows) != activation["crossed_outer_splits"]
        or split_contract.get("crossed_outer_splits") != len(split_rows)
        or split_contract.get("donor_folds") != activation["donor_folds"]
        or split_contract.get("genomic_folds") != activation["genomic_folds"]
    ):
        raise RuntimeFixtureError("roster or crossed-split activation census differs")
    expected_fits = (
        len(split_rows)
        * len(observed_roster)
        * len(activation["fixed_seeds"])
    )
    if expected_fits != activation[
        "provisional_shared_fit_count_if_observed_roster_is_adopted"
    ]:
        raise RuntimeFixtureError("provisional shared-fit count differs")
    return {
        "file_authorities": file_contracts,
        "tree_authorities": tree_contracts,
        "registered_state_roster": activation["registered_state_roster"],
        "observed_frozen_gse296875_roster": observed_roster,
        "roster_match": False,
        "crossed_outer_splits": len(split_rows),
        "donor_folds": split_contract["donor_folds"],
        "genomic_folds": split_contract["genomic_folds"],
        "provisional_shared_fit_count": expected_fits,
    }


def verify_source(source_root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    receipt = load_json(source_root / "source_receipt.json")
    records = receipt.get("files", [])
    expected = {record["path"]: record["sha256"] for record in config["source"]["files"]}
    observed = {record.get("path"): record.get("sha256") for record in records}
    discovered = {
        path.relative_to(source_root).as_posix()
        for path in source_root.rglob("*")
        if path.is_file() and path.name not in {"gkmtrain", "gkmpredict", "libsvm.o", "libsvm_gkm.o"}
    }
    if (
        receipt.get("status") != "pass_exact_source_only_no_historical_weights"
        or receipt.get("revision") != config["source_revision"]
        or receipt.get("historical_weights_acquired") is not False
        or receipt.get("historical_implementation_acquired") is not False
        or receipt.get("repository_cloned") is not False
        or receipt.get("tests_directory_acquired") is not False
        or observed != expected
        or discovered != set(expected) | {"source_receipt.json"}
    ):
        raise RuntimeFixtureError("pinned source receipt or file surface differs")
    for relative, digest in expected.items():
        if file_sha256(source_root / relative) != digest:
            raise RuntimeFixtureError(f"source bytes changed after acquisition: {relative}")
    return receipt


def build_fixture(arguments: argparse.Namespace) -> dict[str, Any]:
    if arguments.output.exists():
        raise RuntimeFixtureError("fixture output exists")
    root = arguments.project_root.resolve(strict=True)
    config = load_config(arguments.config)
    if (
        config.get("status")
        != "prespecified_outcome_blind_cpu_runtime_and_synthetic_fixture"
        or config.get("outcome_access_authorized") is not False
        or config.get("biological_data_access_authorized") is not False
        or config.get("production_training_authorized") is not False
        or config.get("production_prediction_authorized") is not False
        or config.get("historical_weights_allowed") is not False
        or config.get("native_resume_support") != "unsupported_retry_from_scratch_only"
    ):
        raise RuntimeFixtureError("action firewall differs")
    source_root = arguments.source_root.resolve(strict=True)
    source_receipt = verify_source(source_root, config)
    binaries = {
        name: (source_root / "src" / name).resolve(strict=True)
        for name in config["runtime"]["targets"]
    }
    if any(path.is_symlink() or not path.is_file() for path in binaries.values()):
        raise RuntimeFixtureError("native binary surface differs")
    authorities = validate_authorities(root, config)

    output = arguments.output
    for relative in ("fixture", "models", "predictions", "logs", "contract"):
        (output / relative).mkdir(parents=True, exist_ok=True, mode=0o750)
    fixture = config["fixture"]
    generator = random.Random(fixture["seed"])
    positive: list[tuple[str, str]] = []
    negative: list[tuple[str, str]] = []
    for index in range(fixture["positive_records"]):
        background = random_sequence(generator, fixture["sequence_length_bp"])
        positive.append(
            (
                f"positive_{index:03d}",
                insert_motif(
                    background, fixture["positive_motif"], fixture["motif_positions"]
                ),
            )
        )
        negative.append(
            (
                f"negative_{index:03d}",
                insert_motif(
                    background, fixture["negative_motif"], fixture["motif_positions"]
                ),
            )
        )
    write_fasta(output / "fixture/positive.fa", positive)
    write_fasta(output / "fixture/negative.fa", negative)
    validate_fasta(
        output / "fixture/positive.fa",
        expected_length=fixture["sequence_length_bp"],
    )
    validate_fasta(
        output / "fixture/negative.fa",
        expected_length=fixture["sequence_length_bp"],
    )

    probe_background = random_sequence(generator, fixture["sequence_length_bp"])
    probe_records = [
        (
            "positive_probe",
            insert_motif(
                probe_background, fixture["positive_motif"], fixture["motif_positions"]
            ),
        ),
        (
            "negative_probe",
            insert_motif(
                probe_background, fixture["negative_motif"], fixture["motif_positions"]
            ),
        ),
    ]
    write_fasta(output / "fixture/probe.fa", probe_records)
    write_fasta(
        output / "fixture/probe_rc.fa",
        [(f"{identifier}_rc", reverse_complement(sequence)) for identifier, sequence in probe_records],
    )
    validate_fasta(output / "fixture/probe.fa", expected_length=300)
    validate_fasta(output / "fixture/probe_rc.fa", expected_length=300)

    variant_background = random_sequence(generator, 300)
    alternative = insert_motif(variant_background, fixture["positive_motif"], [147])
    reference_values = list(alternative)
    if reference_values[150] != fixture["positive_motif"][3]:
        raise RuntimeFixtureError("variant fixture center does not map to motif")
    reference_values[150] = "A" if reference_values[150] != "A" else "T"
    reference = "".join(reference_values)
    if [index for index, pair in enumerate(zip(reference, alternative, strict=True)) if pair[0] != pair[1]] != [150]:
        raise RuntimeFixtureError("full-window allele fixture is not a center SNP")
    allele_records = [("variant_REF", reference), ("variant_ALT", alternative)]
    write_fasta(output / "fixture/alleles.fa", allele_records)
    validate_fasta(output / "fixture/alleles.fa", expected_length=300)

    commands: list[dict[str, Any]] = []
    model_receipts: list[dict[str, Any]] = []
    for repeat in range(1, fixture["train_repeats"] + 1):
        prefix = output / f"models/synthetic_repeat_{repeat}"
        command = [
            binaries["gkmtrain"].as_posix(),
            "-t", "2", "-l", "11", "-k", "7", "-d", "3",
            "-c", "1", "-e", "0.001", "-w", "1", "-m", "4096",
            "-T", "1", "-z",
            (output / "fixture/positive.fa").as_posix(),
            (output / "fixture/negative.fa").as_posix(),
            prefix.as_posix(),
        ]
        commands.append(
            command_receipt(
                command,
                stdout_path=output / f"logs/train_{repeat}.stdout.txt",
                stderr_path=output / f"logs/train_{repeat}.stderr.txt",
            )
        )
        model_receipts.append(
            validate_project_model(prefix.with_suffix(".model.txt.gz"), support_vector_length=300)
        )
    model_hashes = [receipt["sha256"] for receipt in model_receipts]
    if len(set(model_hashes)) != 1:
        raise RuntimeFixtureError("repeat LS-GKM training is not byte deterministic")

    model = output / "models/synthetic_repeat_1.model.txt.gz"
    prediction_jobs = (
        ("probe_repeat_1", output / "fixture/probe.fa"),
        ("probe_repeat_2", output / "fixture/probe.fa"),
        ("probe_rc", output / "fixture/probe_rc.fa"),
        ("alleles_repeat_1", output / "fixture/alleles.fa"),
        ("alleles_repeat_2", output / "fixture/alleles.fa"),
    )
    parsed_predictions: dict[str, dict[str, float]] = {}
    for label, fasta in prediction_jobs:
        prediction = output / f"predictions/{label}.tsv"
        command = [
            binaries["gkmpredict"].as_posix(),
            "-T", "1", fasta.as_posix(), model.as_posix(), prediction.as_posix(),
        ]
        commands.append(
            command_receipt(
                command,
                stdout_path=output / f"logs/{label}.stdout.txt",
                stderr_path=output / f"logs/{label}.stderr.txt",
            )
        )
        identifiers = [identifier for identifier, _ in validate_fasta(fasta, expected_length=300)]
        parsed_predictions[label] = parse_predictions(
            prediction, expected_identifiers=identifiers
        )
    if (
        file_sha256(output / "predictions/probe_repeat_1.tsv")
        != file_sha256(output / "predictions/probe_repeat_2.tsv")
        or file_sha256(output / "predictions/alleles_repeat_1.tsv")
        != file_sha256(output / "predictions/alleles_repeat_2.tsv")
    ):
        raise RuntimeFixtureError("repeat native prediction bytes differ")
    probe = parsed_predictions["probe_repeat_1"]
    probe_rc = parsed_predictions["probe_rc"]
    tolerance = fixture["reverse_complement_absolute_tolerance"]
    if (
        probe["positive_probe"] <= 0
        or probe["negative_probe"] >= 0
        or probe["positive_probe"] <= probe["negative_probe"]
        or abs(probe["positive_probe"] - probe_rc["positive_probe_rc"]) > tolerance
        or abs(probe["negative_probe"] - probe_rc["negative_probe_rc"]) > tolerance
    ):
        raise RuntimeFixtureError("synthetic direction or reverse-complement fixture failed")
    alleles = parsed_predictions["alleles_repeat_1"]
    direct_delta = direct_gkmsvm_alt_minus_ref(
        reference_score=alleles["variant_REF"],
        alternative_score=alleles["variant_ALT"],
    )
    if not math.isfinite(direct_delta) or direct_delta <= 0:
        raise RuntimeFixtureError("full-window ALT-minus-REF direction fixture failed")

    nrkmers_path = output / "fixture/nrkmers_3.fa"
    commands.append(
        command_receipt(
            [
                sys.executable,
                (source_root / "scripts/nrkmers.py").as_posix(),
                str(fixture["nrkmers_probe_length"]),
                nrkmers_path.as_posix(),
            ],
            stdout_path=output / "logs/nrkmers.stdout.txt",
            stderr_path=output / "logs/nrkmers.stderr.txt",
        )
    )
    nrkmers = validate_fasta(
        nrkmers_path, expected_length=fixture["nrkmers_probe_length"]
    )
    if (
        len(nrkmers) != fixture["nrkmers_expected_records"]
        or len({canonical_kmer(sequence) for _, sequence in nrkmers}) != len(nrkmers)
    ):
        raise RuntimeFixtureError("nonredundant k-mer generator fixture differs")

    reference_context = "AACCGGTTAACCGGTTAACCG"
    alternative_context = reference_context[:10] + "G" + reference_context[11:]
    keys = {
        canonical_kmer(context[start : start + 11])
        for context in (reference_context, alternative_context)
        for start in range(11)
    }
    weights = {
        key: (int(sha256(key.encode("ascii")).hexdigest()[:8], 16) / 2**32) - 0.5
        for key in keys
    }
    delta = deltasvm_scores(reference_context, alternative_context, weights)
    manual_reference = math.fsum(
        weights[canonical_kmer(reference_context[start : start + 11])]
        for start in range(11)
    )
    manual_alternative = math.fsum(
        weights[canonical_kmer(alternative_context[start : start + 11])]
        for start in range(11)
    )
    reverse_delta = deltasvm_scores(
        reverse_complement(reference_context),
        reverse_complement(alternative_context),
        weights,
    )
    if (
        delta["reference_11mer_sum"] != manual_reference
        or delta["alternative_11mer_sum"] != manual_alternative
        or delta["native_ref_minus_alt"] != manual_reference - manual_alternative
        or delta["canonical_alt_minus_ref"] != -delta["native_ref_minus_alt"]
        or reverse_delta != delta
    ):
        raise RuntimeFixtureError("clean-room deltaSVM sign or strand fixture failed")

    activation_plan = {
        "schema_version": "masld-bench-lsgkm-variant-activation-plan-v1",
        "status": config["activation"]["status"],
        "family_id": config["family_id"],
        "shared_fit_readouts": config["algorithm"]["shared_fit_readouts"],
        "independent_model_families": 1,
        "authorities": authorities,
        "fit_contract": config["activation"]["fit_contract"],
        "claim_contract": config["activation"]["claim_contract"],
        "missing_activation_authorities": config["activation"]["missing_activation_authorities"],
        "existing_chrombpnet_negatives_admissible": False,
        "variant_scoring_source": config["activation"]["variant_scoring_source"],
        "full_11mer_table_records_per_shared_fit": config["activation"]["full_11mer_table_records_per_shared_fit"],
        "provisional_shared_fit_count_if_observed_roster_is_adopted": authorities["provisional_shared_fit_count"],
        "readout_count_per_shared_fit": 2,
        "independent_evidence_count_per_shared_fit": 1,
        "outcomes_read": False,
        "production_training_executed": False,
        "production_predictions_generated": False,
    }
    write_json(output / "contract/variant_scoring_activation_plan.json", activation_plan)
    runtime_receipt = {
        "schema_version": "masld-bench-lsgkm-runtime-fixture-receipt-v1",
        "status": "pass_cpu_runtime_and_synthetic_fixture",
        "family_id": config["family_id"],
        "source_revision": config["source_revision"],
        "source_tree_sha256": config["source_tree_sha256"],
        "source_receipt_sha256": file_sha256(source_root / "source_receipt.json"),
        "source_file_count": source_receipt["admitted_file_count"],
        "binary_sha256": {name: file_sha256(path) for name, path in binaries.items()},
        "models": model_receipts,
        "repeat_model_sha256_equal": True,
        "repeat_prediction_bytes_equal": True,
        "positive_probe_score": probe["positive_probe"],
        "negative_probe_score": probe["negative_probe"],
        "reverse_complement_tolerance": tolerance,
        "reverse_complement_fixture_pass": True,
        "direct_gkmsvm_alt_minus_ref": direct_delta,
        "direct_gkmsvm_direction_fixture_pass": True,
        "deltasvm_hand_calculated_fixture": delta,
        "deltasvm_reverse_complement_fixture_pass": True,
        "deltasvm_sign_conversion_once": True,
        "nrkmers_probe_length": fixture["nrkmers_probe_length"],
        "nrkmers_probe_records": len(nrkmers),
        "native_resume_support": "unsupported_retry_from_scratch_only",
        "native_commands": commands,
        "shared_fit_readouts": ["direct_gkmsvm", "deltasvm"],
        "independent_model_families": 1,
        "historical_weights_read": False,
        "historical_implementation_read": False,
        "outcomes_read": False,
        "biological_data_read_for_training": False,
        "benchmark_metrics_calculated": False,
        "synthetic_model_training_executed": True,
        "production_model_training_executed": False,
        "production_predictions_generated": False,
    }
    write_json(output / "runtime_receipt.json", runtime_receipt)
    return runtime_receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    try:
        result = build_fixture(parse_args())
    except LSGKMSafetyError as error:
        raise SystemExit(f"LS-GKM safety failure: {error}") from error
    print(json.dumps(result, sort_keys=True))
