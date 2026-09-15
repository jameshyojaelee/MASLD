#!/usr/bin/env python3
"""Run one authorized outcome-blind LS-GKM full-campaign fit and held scoring pass."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

from scripts.freeze_lsgkm_gse281364_full_campaign_authorization import (
    expected_run_fits,
    validate_config as validate_campaign_config,
)
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
from scripts.run_lsgkm_gse281364_scale_probe import (
    SCORE_FIELDS,
    load_json,
    read_fasta_gzip,
    read_tsv_gzip,
    resolve_file,
    run_timed,
    verify_manifest_members,
    write_fasta,
    write_tsv_gzip,
)


class LSGKMFitError(RuntimeError):
    """Raised when a full-campaign fit or authority differs."""


def manifest_roster(path: Path) -> dict[str, str]:
    manifest = load_json(path)
    members = manifest.get("artifacts")
    if not isinstance(members, list):
        raise LSGKMFitError("artifact manifest differs")
    roster: dict[str, str] = {}
    for member in members:
        if not isinstance(member, dict) or not isinstance(member.get("path"), str):
            raise LSGKMFitError("artifact member differs")
        relative = member["path"]
        if relative in roster or not isinstance(member.get("sha256"), str):
            raise LSGKMFitError("artifact member roster differs")
        roster[relative] = member["sha256"]
    return roster


def member_binding(
    root: Path, authority_path: str, roster: Mapping[str, str], relative_from_authority: str
) -> tuple[str, str]:
    project_relative = (Path(authority_path).parent / relative_from_authority).as_posix()
    expected = roster.get(relative_from_authority)
    if expected is None:
        raise LSGKMFitError(f"required input member is absent: {relative_from_authority}")
    resolve_file(root, project_relative, expected)
    return project_relative, expected


def validate_fit_request(
    config: Mapping[str, Any], authorization: Mapping[str, Any], split_id: str, model_seed: int
) -> None:
    validate_campaign_config(config)
    requested = (split_id, model_seed)
    if requested not in set(expected_run_fits(config)):
        raise LSGKMFitError("fit is not in the remaining authorized rectangle")
    if authorization.get("status") != "authorized_remaining_24_outcome_blind_fits":
        raise LSGKMFitError("full-campaign authorization status differs")
    authorized = {
        (row.get("split_id"), row.get("model_seed"))
        for row in authorization.get("authorized_run_fits", [])
        if isinstance(row, dict)
    }
    if authorized != set(expected_run_fits(config)) or requested not in authorized:
        raise LSGKMFitError("authorization fit roster differs")
    firewall = authorization.get("action_firewall", {})
    if (
        firewall.get("remaining_24_fits_authorized") is not True
        or firewall.get("held_sequence_scoring_authorized") is not True
        or any(
            firewall.get(field) is not False
            for field in (
                "probe_fit_retraining_authorized",
                "outcome_access_authorized",
                "reporter_count_access_authorized",
                "sealed_asset_access_authorized",
                "control_prediction_value_access_authorized",
                "benchmark_metric_calculation_authorized",
                "aggregation_or_evaluation_authorized",
            )
        )
    ):
        raise LSGKMFitError("authorization firewall differs")


def bind_authorization(
    root: Path, artifacts_path: Path, artifacts_sha256: str
) -> dict[str, Any]:
    relative_authority = artifacts_path.resolve(strict=True).relative_to(root).as_posix()
    authority = resolve_file(root, relative_authority, artifacts_sha256)
    roster = manifest_roster(authority)
    receipt_relative = "authorization/full_campaign_authorization.json"
    receipt_sha = roster.get(receipt_relative)
    if receipt_sha is None:
        raise LSGKMFitError("authorization receipt is absent from authority")
    receipt_path = authority.parent / receipt_relative
    if file_sha256(receipt_path) != receipt_sha:
        raise LSGKMFitError("authorization receipt authority differs")
    return load_json(receipt_path)


def _read_scoring_rows(path: Path, split_id: str) -> list[dict[str, str]]:
    rows = read_tsv_gzip(path)
    held = [row for row in rows if row["scoring_split_id"] == split_id]
    if not held or len({row["element_id"] for row in held}) != len(held):
        raise LSGKMFitError("held scoring roster differs")
    return held


def run_fit(
    root: Path,
    config_path: Path,
    authorization_artifacts: Path,
    authorization_artifacts_sha256: str,
    split_id: str,
    model_seed: int,
    output: Path,
) -> dict[str, Any]:
    config = load_json(config_path)
    authorization = bind_authorization(
        root, authorization_artifacts, authorization_artifacts_sha256
    )
    validate_fit_request(config, authorization, split_id, model_seed)
    if output.exists():
        raise LSGKMFitError("fit output exists")
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
    input_authority = resolve_file(
        root, inputs["artifacts_path"], inputs["artifacts_sha256"]
    )
    roster = manifest_roster(input_authority)
    base = "materialization"
    positive_binding = member_binding(
        root,
        inputs["artifacts_path"],
        roster,
        f"{base}/inputs/{split_id}.positive.fa.gz",
    )
    negative_binding = member_binding(
        root,
        inputs["artifacts_path"],
        roster,
        f"{base}/inputs/{split_id}.seed{model_seed}.negative.fa.gz",
    )
    pair_binding = member_binding(
        root,
        inputs["artifacts_path"],
        roster,
        f"{base}/inputs/{split_id}.seed{model_seed}.pairs.tsv.gz",
    )
    fixed_bindings = (
        (inputs["independent_audit_path"], inputs["independent_audit_sha256"]),
        (inputs["materialization_contract_path"], inputs["materialization_contract_sha256"]),
        (inputs["scoring_manifest_path"], inputs["scoring_manifest_sha256"]),
        (inputs["scoring_fasta_path"], inputs["scoring_fasta_sha256"]),
    )
    verify_manifest_members(
        root,
        inputs["artifacts_path"],
        inputs["artifacts_sha256"],
        (*fixed_bindings, positive_binding, negative_binding, pair_binding),
    )

    audit = load_json(root / inputs["independent_audit_path"])
    if (
        audit.get("status") != "pass_independent_emitted_fasta_invariant_audit"
        or audit.get("fit_count") != 25
        or audit.get("audited_pair_rows") != 250000
        or audit.get("outcomes_read") is not False
    ):
        raise LSGKMFitError("input audit differs")
    positive_records = read_fasta_gzip(root / positive_binding[0])
    negative_records = read_fasta_gzip(root / negative_binding[0])
    pairs = read_tsv_gzip(root / pair_binding[0])
    if (
        len(positive_records) != 10000
        or len(negative_records) != 10000
        or len(pairs) != 10000
        or any(int(row["seed"]) != model_seed for row in pairs)
        or any(row["split_id"] != split_id for row in pairs)
    ):
        raise LSGKMFitError("fit input denominator or identity differs")
    positive_path = output / "inputs/positive.fa"
    negative_path = output / "inputs/negative.fa"
    write_fasta(positive_path, positive_records)
    write_fasta(negative_path, negative_records)
    if (
        len(validate_fasta(positive_path, expected_length=300)) != 10000
        or len(validate_fasta(negative_path, expected_length=300)) != 10000
    ):
        raise LSGKMFitError("safe-adapter fit input differs")

    model_prefix = output / "models/fit"
    train_receipt = run_timed(
        [
            str(root / runtime["gkmtrain_path"]),
            *runtime["gkmtrain_argv"],
            str(positive_path),
            str(negative_path),
            str(model_prefix),
        ],
        stdout_path=output / "logs/gkmtrain.stdout.txt",
        stderr_path=output / "logs/gkmtrain.stderr.txt",
        time_path=output / "logs/gkmtrain.time.txt",
    )
    model_path = output / "models/fit.model.txt.gz"
    model_receipt = validate_project_model(model_path, support_vector_length=300)

    held_rows = _read_scoring_rows(root / inputs["scoring_manifest_path"], split_id)
    expected_held = config["fit_grid"]["splits"][split_id]["expected_held_elements"]
    if len(held_rows) != expected_held:
        raise LSGKMFitError("held scoring denominator differs")
    scoring_records = dict(read_fasta_gzip(root / inputs["scoring_fasta_path"]))
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
        direct_path,
        expected_identifiers=[identifier for identifier, _ in held_fasta_records],
    )

    required_kmers: set[str] = set()
    for row in held_rows:
        element = row["element_id"]
        for allele in ("REF", "ALT"):
            context = scoring_records[f"{element}|{allele}"][140:161]
            for start in range(11):
                required_kmers.add(canonical_kmer(context[start : start + 11]))
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
        expected_identifiers=[identifier for identifier, _ in kmer_records],
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
            raise LSGKMFitError("non-finite held score")
        score_rows.append(
            {
                "split_id": split_id,
                "model_seed": model_seed,
                "model_state_id": config["fit_grid"]["model_state_id"],
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
        "schema_version": "masld-bench-lsgkm-gse281364-dinucleotide-fit-v1",
        "status": "pass_one_authorized_full_campaign_fit",
        "design_id": config["design_id"],
        "config_sha256": file_sha256(config_path),
        "authorization_artifacts_sha256": authorization_artifacts_sha256,
        "runtime_artifacts_sha256": runtime["artifacts_sha256"],
        "input_artifacts_sha256": inputs["artifacts_sha256"],
        "split_id": split_id,
        "model_seed": model_seed,
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
    (output / "fit_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--authorization-artifacts", type=Path, required=True)
    parser.add_argument("--authorization-artifacts-sha256", required=True)
    parser.add_argument("--split-id", required=True)
    parser.add_argument("--model-seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    run_fit(
        arguments.project_root.resolve(strict=True),
        arguments.config.resolve(strict=True),
        arguments.authorization_artifacts,
        arguments.authorization_artifacts_sha256,
        arguments.split_id,
        arguments.model_seed,
        arguments.output,
    )


if __name__ == "__main__":
    main()
