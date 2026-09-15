#!/usr/bin/env python3
"""Freeze outcome-blind donor and whole-chromosome masks for GSE296875."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import importlib.metadata
import json
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.adapters.rna_atac_classical import fold_index
from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-observed-multiome-mask-plan-v1"
PLAN_ID = "gse296875_observed_multiome_mask_plan_20260825"


class ObservedMultiomeMaskPlanError(ValueError):
    """Raised when an identifier-only mask plan does not meet its frozen requirements."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canonical_digest(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ObservedMultiomeMaskPlanError(f"JSON object required: {path}")
    return value


def _resolve_bound_file(
    root: Path, record: Mapping[str, Any], *, label: str
) -> Path:
    path = (root / str(record.get("path", ""))).resolve(strict=True)
    path.relative_to(root)
    if path.is_symlink() or not path.is_file():
        raise ObservedMultiomeMaskPlanError(f"{label} is not a regular file")
    if _digest(path) != record.get("sha256"):
        raise ObservedMultiomeMaskPlanError(f"{label} SHA-256 drifted")
    return path


def _resolve_bound_tree(
    root: Path, record: Mapping[str, Any], *, label: str
) -> Path:
    path = (root / str(record.get("path", ""))).resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise ObservedMultiomeMaskPlanError(f"{label} ARTIFACTS.json drifted")
    return path


def validate_config(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("plan_id") != PLAN_ID
        or config.get("dataset_id") != "gse296875"
        or config.get("dataset_view_id") != "gse296875_rna_atac_smoke_1000_v1"
        or config.get("stage") != "smoke"
    ):
        raise ObservedMultiomeMaskPlanError("mask-plan identity differs")
    parents = config.get("parents")
    expected_parents = {
        "dataset_view_tree",
        "dataset_view_config",
        "target_axis_contract",
        "factorized_fixture",
        "donor_fold_implementation",
    }
    if not isinstance(parents, dict) or set(parents) != expected_parents:
        raise ObservedMultiomeMaskPlanError("parent artifact roster differs")
    resolved = {
        "dataset_view_tree": _resolve_bound_tree(
            root, parents["dataset_view_tree"], label="dataset view tree"
        ),
        "dataset_view_config": _resolve_bound_file(
            root, parents["dataset_view_config"], label="dataset view config"
        ),
        "target_axis_contract": _resolve_bound_tree(
            root, parents["target_axis_contract"], label="target-axis contract"
        ),
        "factorized_fixture": _resolve_bound_tree(
            root, parents["factorized_fixture"], label="factorized fixture"
        ),
        "donor_fold_implementation": _resolve_bound_file(
            root,
            parents["donor_fold_implementation"],
            label="donor-fold implementation",
        ),
    }
    h5_record = config.get("input_h5")
    if not isinstance(h5_record, dict):
        raise ObservedMultiomeMaskPlanError("input HDF5 contract is absent")
    resolved["input_h5"] = _resolve_bound_file(root, h5_record, label="input HDF5")
    expected_lineages = [
        "cholangiocyte",
        "fibroblast",
        "hepatocyte",
        "macrophage",
        "t_cell",
    ]
    if (
        h5_record.get("expected_nuclei") != 1000
        or h5_record.get("expected_donors") != 39
        or h5_record.get("expected_lineages") != expected_lineages
        or h5_record.get("expected_peaks") != 306706
        or h5_record.get("pairing_level") != "same_nucleus"
        or h5_record.get("coordinate_system") != "0_based_half_open"
    ):
        raise ObservedMultiomeMaskPlanError("input HDF5 census differs")

    split = config.get("split")
    canonical = [f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"]
    if (
        not isinstance(split, dict)
        or split.get("donor_folds") != 5
        or split.get("donor_seed") != 20260821
        or split.get("genomic_folds") != 5
        or split.get("genomic_seed") != 20260825
        or split.get("canonical_contigs") != canonical
        or split.get("component_policy")
        != "whole_chromosome_indivisible_greedy_peak_count_balance"
        or split.get("largest_receptive_field_buffer_bp") != 524288
        or split.get("targets_per_genomic_fold") != 1000
        or split.get("observed_atac_input_peaks_per_held_fold") != 2000
        or split.get("target_selection_policy")
        != "identifier_hash_only_within_genomic_fold"
        or split.get("input_selection_policy")
        != "identifier_hash_only_excluding_all_held_fold_chromosomes"
    ):
        raise ObservedMultiomeMaskPlanError("split or mask contract differs")
    if not split.get("donor_hash_namespace") or not split.get("peak_hash_namespace"):
        raise ObservedMultiomeMaskPlanError("hash namespaces are absent")

    firewall = config.get("firewall")
    allowed = {
        "obs/donor_id",
        "obs/broad_label",
        "atac/peak_id",
        "atac/chromosome",
        "atac/bed_start_0based",
        "atac/bed_end_half_open",
    }
    if not isinstance(firewall, dict) or set(firewall.get("allowed_h5_datasets", [])) != allowed:
        raise ObservedMultiomeMaskPlanError("HDF5 metadata allowlist differs")
    false_fields = {
        "rna_count_values_read",
        "atac_count_values_read",
        "development_outcomes_read",
        "sealed_outcomes_read",
        "benchmark_metrics_calculated",
        "raw_donor_ids_exported",
        "partial_rectangle_ranking_allowed",
    }
    if any(firewall.get(field) is not False for field in false_fields):
        raise ObservedMultiomeMaskPlanError("outcome firewall is open")
    disposition = config.get("disposition")
    if (
        not isinstance(disposition, dict)
        or disposition.get("biological_identifier_plan") is not True
        or disposition.get("biological_matrix_execution_authorized") is not False
        or disposition.get("production_ranking_authorized") is not False
        or disposition.get("champion_claim_allowed") is not False
        or disposition.get("next_gate")
        != "separate_model_input_and_evaluator_outcome_materialization"
    ):
        raise ObservedMultiomeMaskPlanError("plan disposition differs")
    return resolved


def _decode_strings(values: Iterable[Any], *, label: str) -> list[str]:
    result: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            value = value.decode("utf-8")
        text = str(value)
        if not text:
            raise ObservedMultiomeMaskPlanError(f"empty {label}")
        result.append(text)
    return result


def read_identifier_axes(path: Path) -> dict[str, Any]:
    """Read the six allowed identifier/coordinate datasets and no matrix values."""

    import h5py

    with h5py.File(path, "r") as handle:
        if (
            handle.attrs.get("schema_version") != "masld-bench-multimodal-h5-v1"
            or handle.attrs.get("view_id") != "gse296875_rna_atac_smoke_1000_v1"
            or handle.attrs.get("pairing_level") != "same_nucleus"
            or handle.attrs.get("bed_coordinate_system") != "0_based_half_open"
        ):
            raise ObservedMultiomeMaskPlanError("multimodal HDF5 attributes differ")
        donors = _decode_strings(handle["obs/donor_id"][:], label="donor ID")
        lineages = _decode_strings(handle["obs/broad_label"][:], label="lineage")
        peak_ids = _decode_strings(handle["atac/peak_id"][:], label="peak ID")
        chromosomes = _decode_strings(
            handle["atac/chromosome"][:], label="chromosome"
        )
        starts = [int(value) for value in handle["atac/bed_start_0based"][:]]
        ends = [int(value) for value in handle["atac/bed_end_half_open"][:]]
    return {
        "donors": donors,
        "lineages": lineages,
        "peak_ids": peak_ids,
        "chromosomes": chromosomes,
        "starts": starts,
        "ends": ends,
    }


def _salted_hash(namespace: str, *parts: object) -> str:
    payload = "\0".join((namespace, *(str(part) for part in parts)))
    return sha256(payload.encode("utf-8")).hexdigest()


def _priority(seed: int, purpose: str, identifier: str) -> tuple[bytes, str]:
    digest = sha256(f"{seed}\0{purpose}\0{identifier}".encode("utf-8")).digest()
    return digest, identifier


def assign_chromosome_folds(
    chromosomes: Sequence[str], *, canonical_contigs: Sequence[str], folds: int, seed: int
) -> tuple[dict[str, int], dict[str, int]]:
    counts = Counter(chromosomes)
    if not counts or not set(counts).issubset(set(canonical_contigs)):
        raise ObservedMultiomeMaskPlanError("a peak is outside the canonical chromosome roster")
    order = sorted(
        canonical_contigs,
        key=lambda contig: (-counts[contig], _priority(seed, "chromosome", contig), contig),
    )
    totals = [0] * folds
    assignment: dict[str, int] = {}
    for contig in order:
        fold = min(range(folds), key=lambda value: (totals[value], value))
        assignment[contig] = fold
        totals[fold] += counts[contig]
    if set(assignment.values()) != set(range(folds)):
        raise ObservedMultiomeMaskPlanError("one or more genomic folds are empty")
    return assignment, {contig: counts[contig] for contig in canonical_contigs}


def build_plan(
    axes: Mapping[str, Sequence[Any]],
    *,
    donor_folds: int,
    donor_seed: int,
    genomic_folds: int,
    genomic_seed: int,
    canonical_contigs: Sequence[str],
    targets_per_fold: int,
    input_peaks_per_fold: int,
    donor_namespace: str,
    peak_namespace: str,
) -> dict[str, Any]:
    donors = list(axes["donors"])
    lineages = list(axes["lineages"])
    peak_ids = list(axes["peak_ids"])
    chromosomes = list(axes["chromosomes"])
    starts = [int(value) for value in axes["starts"]]
    ends = [int(value) for value in axes["ends"]]
    if len(donors) != len(lineages) or not donors:
        raise ObservedMultiomeMaskPlanError("biological row axes differ")
    if not (
        len(peak_ids) == len(chromosomes) == len(starts) == len(ends)
        and len(set(peak_ids)) == len(peak_ids)
    ):
        raise ObservedMultiomeMaskPlanError("peak axes differ or IDs are duplicated")
    if any(start < 0 or end <= start for start, end in zip(starts, ends, strict=True)):
        raise ObservedMultiomeMaskPlanError("BED coordinates are invalid")

    donor_ids = sorted(set(donors))
    donor_assignment = {
        donor: fold_index(donor, seed=donor_seed, outer_folds=donor_folds)
        for donor in donor_ids
    }
    if set(donor_assignment.values()) != set(range(donor_folds)):
        raise ObservedMultiomeMaskPlanError("one or more donor folds are empty")
    donor_lineage_counts = Counter(zip(donors, lineages, strict=True))
    donor_rows = [
        {
            "donor_hash": _salted_hash(donor_namespace, donor),
            "lineage": lineage,
            "donor_fold": donor_assignment[donor],
            "nuclei": count,
        }
        for (donor, lineage), count in sorted(donor_lineage_counts.items())
    ]

    chromosome_assignment, chromosome_counts = assign_chromosome_folds(
        chromosomes,
        canonical_contigs=canonical_contigs,
        folds=genomic_folds,
        seed=genomic_seed,
    )
    genomic_rows = [
        {
            "chromosome": contig,
            "genomic_fold": chromosome_assignment[contig],
            "eligible_peaks": chromosome_counts[contig],
        }
        for contig in canonical_contigs
    ]
    peaks = [
        {
            "peak_id": peak_id,
            "chromosome": chromosome,
            "bed_start_0based": start,
            "bed_end_half_open": end,
            "genomic_fold": chromosome_assignment[chromosome],
        }
        for peak_id, chromosome, start, end in zip(
            peak_ids, chromosomes, starts, ends, strict=True
        )
    ]

    target_rows: list[dict[str, Any]] = []
    input_rows: list[dict[str, Any]] = []
    for held_fold in range(genomic_folds):
        held = [row for row in peaks if row["genomic_fold"] == held_fold]
        allowed = [row for row in peaks if row["genomic_fold"] != held_fold]
        if len(held) < targets_per_fold or len(allowed) < input_peaks_per_fold:
            raise ObservedMultiomeMaskPlanError("peak budget exceeds a genomic fold")
        held.sort(
            key=lambda row: _priority(genomic_seed, f"target:{held_fold}", row["peak_id"])
        )
        allowed.sort(
            key=lambda row: _priority(genomic_seed, f"input:{held_fold}", row["peak_id"])
        )
        held_chromosomes = {
            contig for contig, fold in chromosome_assignment.items() if fold == held_fold
        }
        for row in held[:targets_per_fold]:
            target_rows.append(
                {
                    "target_hash": _salted_hash(peak_namespace, "target", row["peak_id"]),
                    "chromosome": row["chromosome"],
                    "bed_start_0based": row["bed_start_0based"],
                    "bed_end_half_open": row["bed_end_half_open"],
                    "genomic_fold": held_fold,
                }
            )
        for row in allowed[:input_peaks_per_fold]:
            if row["chromosome"] in held_chromosomes:
                raise ObservedMultiomeMaskPlanError("held chromosome entered ATAC input")
            input_rows.append(
                {
                    "input_hash": _salted_hash(
                        peak_namespace, "input", held_fold, row["peak_id"]
                    ),
                    "chromosome": row["chromosome"],
                    "bed_start_0based": row["bed_start_0based"],
                    "bed_end_half_open": row["bed_end_half_open"],
                    "source_genomic_fold": row["genomic_fold"],
                    "held_genomic_fold": held_fold,
                }
            )
    if len({row["target_hash"] for row in target_rows}) != len(target_rows):
        raise ObservedMultiomeMaskPlanError("target hashes are duplicated")
    if len({row["input_hash"] for row in input_rows}) != len(input_rows):
        raise ObservedMultiomeMaskPlanError("input hashes are duplicated")
    return {
        "donor_rows": donor_rows,
        "genomic_rows": genomic_rows,
        "target_rows": target_rows,
        "input_rows": input_rows,
        "donor_fold_counts": Counter(donor_assignment.values()),
        "genomic_fold_peak_counts": Counter(
            chromosome_assignment[chromosome] for chromosome in chromosomes
        ),
        "lineages": sorted(set(lineages)),
        "raw_donor_ids": donor_ids,
    }


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_plan(output: Path, plan: Mapping[str, Any]) -> None:
    _write_tsv(
        output / "donor_lineage_units.tsv",
        ("donor_hash", "lineage", "donor_fold", "nuclei"),
        plan["donor_rows"],
    )
    _write_tsv(
        output / "genomic_folds.tsv",
        ("chromosome", "genomic_fold", "eligible_peaks"),
        plan["genomic_rows"],
    )
    _write_tsv(
        output / "target_peaks.tsv",
        (
            "target_hash",
            "chromosome",
            "bed_start_0based",
            "bed_end_half_open",
            "genomic_fold",
        ),
        plan["target_rows"],
    )
    _write_tsv(
        output / "observed_atac_input_peaks.tsv",
        (
            "input_hash",
            "chromosome",
            "bed_start_0based",
            "bed_end_half_open",
            "source_genomic_fold",
            "held_genomic_fold",
        ),
        plan["input_rows"],
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    config = _load_json(config_path)
    resolved = validate_config(root, config)
    axes = read_identifier_axes(resolved["input_h5"])
    h5_contract = config["input_h5"]
    if (
        len(axes["donors"]) != h5_contract["expected_nuclei"]
        or len(set(axes["donors"])) != h5_contract["expected_donors"]
        or sorted(set(axes["lineages"])) != h5_contract["expected_lineages"]
        or len(axes["peak_ids"]) != h5_contract["expected_peaks"]
    ):
        raise ObservedMultiomeMaskPlanError("identifier-axis census differs")
    split = config["split"]
    plan = build_plan(
        axes,
        donor_folds=split["donor_folds"],
        donor_seed=split["donor_seed"],
        genomic_folds=split["genomic_folds"],
        genomic_seed=split["genomic_seed"],
        canonical_contigs=split["canonical_contigs"],
        targets_per_fold=split["targets_per_genomic_fold"],
        input_peaks_per_fold=split["observed_atac_input_peaks_per_held_fold"],
        donor_namespace=split["donor_hash_namespace"],
        peak_namespace=split["peak_hash_namespace"],
    )
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    write_plan(output, plan)
    exported_tokens = {
        token
        for path in sorted(output.glob("*.tsv"))
        for line in path.read_text(encoding="utf-8").splitlines()
        for token in line.split("\t")
    }
    if set(plan["raw_donor_ids"]).intersection(exported_tokens):
        raise ObservedMultiomeMaskPlanError("a raw donor ID entered plan outputs")
    receipt = {
        "schema_version": "masld-bench-observed-multiome-mask-plan-receipt-v1",
        "plan_id": PLAN_ID,
        "dataset_id": "gse296875",
        "dataset_view_id": "gse296875_rna_atac_smoke_1000_v1",
        "stage": "smoke",
        "pairing_level": "same_nucleus",
        "biological_unit": "donor_by_lineage",
        "nuclei": len(axes["donors"]),
        "donors": len(set(axes["donors"])),
        "lineages": plan["lineages"],
        "donor_lineage_units": len(plan["donor_rows"]),
        "donor_fold_counts": {
            str(fold): plan["donor_fold_counts"][fold]
            for fold in range(split["donor_folds"])
        },
        "canonical_chromosomes": len(split["canonical_contigs"]),
        "genomic_fold_peak_counts": {
            str(fold): plan["genomic_fold_peak_counts"][fold]
            for fold in range(split["genomic_folds"])
        },
        "target_peaks": len(plan["target_rows"]),
        "observed_atac_input_peak_views": len(plan["input_rows"]),
        "whole_chromosome_components": True,
        "largest_receptive_field_buffer_bp": split[
            "largest_receptive_field_buffer_bp"
        ],
        "held_chromosome_in_observed_atac_input": False,
        "metadata_fields_read": config["firewall"]["allowed_h5_datasets"],
        "rna_count_values_read": False,
        "atac_count_values_read": False,
        "development_outcomes_read": False,
        "sealed_outcomes_read": False,
        "benchmark_metrics_calculated": False,
        "raw_donor_ids_exported": False,
        "biological_matrix_execution_authorized": False,
        "production_ranking_authorized": False,
        "champion_claim_allowed": False,
        "next_gate": config["disposition"]["next_gate"],
        "plan_identity_sha256": _canonical_digest(
            {
                "donor_rows": plan["donor_rows"],
                "genomic_rows": plan["genomic_rows"],
                "target_rows": plan["target_rows"],
                "input_rows": plan["input_rows"],
            }
        ),
        "runtime": {
            "python": sys.version.split()[0],
            "h5py": importlib.metadata.version("h5py"),
        },
    }
    write_json_exclusive(output / "receipt.json", receipt)
    source_paths = [
        config_path,
        Path(__file__).resolve(strict=True),
        root / "tests/unit/test_gse296875_observed_multiome_mask_plan.py",
        resolved["dataset_view_config"],
        resolved["donor_fold_implementation"],
    ]
    (output / "source.sha256").write_text(
        "".join(f"{_digest(path)}  {path}\n" for path in source_paths),
        encoding="utf-8",
    )
    digest = freeze_tree(
        output,
        metadata={
            "artifact_class": "gse296875_observed_multiome_mask_plan",
            "plan_id": PLAN_ID,
            "stage": "smoke",
            "matrix_values_read": False,
            "sealed_outcomes_accessed": False,
            "biological_matrix_execution_authorized": False,
        },
    )
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
