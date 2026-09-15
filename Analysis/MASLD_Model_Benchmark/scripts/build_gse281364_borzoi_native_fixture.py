#!/usr/bin/env python3
"""Build an outcome-blind, long-range-safe Borzoi MPRA screen fixture."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import gzip
from hashlib import sha256
import json
from pathlib import Path
from typing import Iterable, Mapping

from scripts.alphagenome_sei_build_fixture import IndexedFasta


COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
FIXTURE_FIELDS = (
    "fixture_id",
    "element_id",
    "source_locus_group_id",
    "borzoi_long_range_group_id",
    "outer_fold",
    "contig",
    "variant_pos0",
    "variant_pos1",
    "input_start0",
    "input_end0",
    "variant_index0",
    "ref",
    "alt",
    "reference_sequence_sha256",
    "alternative_sequence_sha256",
    "reference_reverse_complement_sha256",
    "alternative_reverse_complement_sha256",
    "allele_effect_sign",
)
SOURCE_GROUP_FIELDS = (
    "source_locus_group_id",
    "borzoi_long_range_group_id",
    "outer_fold",
    "contig",
    "minimum_variant_pos0",
    "maximum_variant_pos0",
    "source_elements",
)
COMPONENT_FIELDS = (
    "borzoi_long_range_group_id",
    "outer_fold",
    "contig",
    "minimum_variant_pos0",
    "maximum_variant_pos0",
    "source_locus_groups",
    "source_elements",
)
TRACK_FIELDS = (
    "index",
    "identifier",
    "file",
    "clip",
    "clip_soft",
    "scale",
    "sum_stat",
    "strand_pair",
    "description",
    "assay",
    "native_screen_role",
)


class BorzoiFixtureError(ValueError):
    """Raised when a frozen Borzoi fixture invariant differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise BorzoiFixtureError(f"TSV has no header: {path.name}")
        return [dict(row) for row in reader]


def write_tsv(
    path: Path, fields: Iterable[str], rows: Iterable[Mapping[str, object]]
) -> None:
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


def load_contract(path: Path) -> dict[str, object]:
    if path.is_symlink() or not path.is_file():
        raise BorzoiFixtureError("screen contract is not a regular file")
    contract = json.loads(path.read_text(encoding="utf-8"))
    if (
        contract.get("schema_version")
        != "masld-bench-gse281364-borzoi-native-screen-v1"
        or contract.get("dataset_id") != "gse281364"
        or contract.get("model_id") != "borzoi_ensemble"
    ):
        raise BorzoiFixtureError("screen contract identity differs")
    sequence = contract.get("sequence_contract", {})
    firewall = contract.get("leakage_firewall", {})
    gate = contract.get("execution_gate", {})
    if (
        sequence.get("input_length_bp") != 524_288
        or sequence.get("variant_index0") != 262_144
        or sequence.get("allele_effect_sign") != "ALT_minus_REF"
        or firewall.get("borzoi_receptive_field_buffer_bp") != 524_288
        or firewall.get("outer_folds") != 5
        or firewall.get("selection_or_fold_fit_on_outcomes") is not False
        or any(
            gate.get(field) is not False
            for field in (
                "checkpoint_download_allowed",
                "checkpoint_deserialization_allowed",
                "model_forward_allowed",
                "model_metrics_allowed_in_fixture_or_adapter",
                "open_champion_eligible",
            )
        )
        or gate.get("weight_terms") != "UNDECLARED_BLOCKED"
    ):
        raise BorzoiFixtureError("geometry, leakage, or execution firewall differs")
    return contract


def validate_authorities(project_root: Path, contract: Mapping[str, object]) -> dict[str, Path]:
    roots: dict[str, Path] = {}
    authorities = contract["authorities"]
    for key in (
        "outcome_blind_split",
        "common_locus_selection",
        "reference",
        "admission_sources",
    ):
        record = authorities[key]
        root = (project_root / record["root"]).resolve(strict=True)
        if project_root.resolve() not in root.parents or root.is_symlink():
            raise BorzoiFixtureError(f"authority escapes project root: {key}")
        if digest(root / "ARTIFACTS.json") != record["artifacts_sha256"]:
            raise BorzoiFixtureError(f"authority artifact changed: {key}")
        roots[key] = root
    for key in ("checkpoint_contract", "exposure_audit", "development_crosswalk"):
        record = authorities[key]
        path = (project_root / record["path"]).resolve(strict=True)
        if project_root.resolve() not in path.parents or path.is_symlink():
            raise BorzoiFixtureError(f"authority file escapes project root: {key}")
        if digest(path) != record["sha256"]:
            raise BorzoiFixtureError(f"authority file changed: {key}")

    split_receipt = json.loads((roots["outcome_blind_split"] / "split/receipt.json").read_text())
    selection_receipt = json.loads(
        (roots["common_locus_selection"] / "fixture/receipt.json").read_text()
    )
    reference_artifacts = json.loads((roots["reference"] / "ARTIFACTS.json").read_text())
    source_receipt = json.loads(
        (roots["admission_sources"] / "review/source_admission_receipt.json").read_text()
    )
    checkpoints = json.loads(
        (project_root / authorities["checkpoint_contract"]["path"]).read_text()
    )
    exposure = json.loads(
        (project_root / authorities["exposure_audit"]["path"]).read_text()
    )
    if (
        split_receipt.get("status") != "pass_outcome_blind_split_contract"
        or split_receipt.get("outcomes_read") is not False
        or split_receipt.get("reporter_counts_read") is not False
        or selection_receipt.get("status") != "pass_outcome_blind_sei_fixture"
        or selection_receipt.get("outcomes_read") is not False
        or reference_artifacts.get("metadata", {}).get("build") != "GRCh38.p14"
        or reference_artifacts.get("metadata", {}).get("indexed") is not True
        or source_receipt.get("checkpoint_bytes_downloaded") is not False
        or source_receipt.get("checkpoint_bytes_deserialized") is not False
        or source_receipt.get("model_forward_executed") is not False
        or checkpoints.get("weight_license") != "UNDECLARED"
        or checkpoints.get("canonical_ensemble", {}).get("ensemble_size") != 4
        or "same genomic split"
        not in checkpoints.get("canonical_ensemble", {}).get("member_relationship", "")
        or exposure.get("checkpoint_findings", {})
        .get("borzoi_ensemble", {})
        .get("exposure_state")
        != "target_label_unexposed"
    ):
        raise BorzoiFixtureError("outcome, reference, ensemble, terms, or exposure authority differs")
    expected_members = list(range(4))
    members = checkpoints["canonical_ensemble"]["members"]
    if [member.get("replicate") for member in members] != expected_members:
        raise BorzoiFixtureError("official ensemble member order differs")
    reference_path = roots["reference"] / authorities["reference"]["fasta"]
    if digest(reference_path) != authorities["reference"]["fasta_sha256"]:
        raise BorzoiFixtureError("reference FASTA bytes changed")
    target_path = roots["admission_sources"] / authorities["admission_sources"][
        "target_manifest"
    ]
    if digest(target_path) != authorities["admission_sources"]["target_manifest_sha256"]:
        raise BorzoiFixtureError("native target manifest changed")
    roots["reference_fasta"] = reference_path
    roots["target_manifest"] = target_path
    return roots


def load_native_tracks(
    path: Path, contract: Mapping[str, object]
) -> list[dict[str, object]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        raw_fields = list(reader.fieldnames or [])
        fields = ["index" if field == "" else field for field in raw_fields]
        rows = [
            {("index" if key == "" else str(key)): str(value) for key, value in row.items()}
            for row in reader
        ]
    output = contract["native_output_contract"]
    if fields != output["manifest_fields"] or len(rows) != output["human_tracks"]:
        raise BorzoiFixtureError("native target schema or row count differs")
    try:
        indexes = [int(row["index"]) for row in rows]
        strand_pairs = [int(row["strand_pair"]) for row in rows]
    except ValueError as exc:
        raise BorzoiFixtureError("native target numeric field differs") from exc
    if indexes != list(range(len(rows))) or any(
        pair < 0
        or pair >= len(rows)
        or strand_pairs[pair] != index
        for index, pair in enumerate(strand_pairs)
    ):
        raise BorzoiFixtureError("target order or strand-pair involution differs")
    assay_counts: dict[str, int] = defaultdict(int)
    primary_roles = output["primary_track_roles"]
    primary_identifiers = output["primary_track_identifiers"]
    rendered: list[dict[str, object]] = []
    for row in rows:
        assay = row["description"].split(":", 1)[0].upper()
        assay_counts[assay] += 1
        index = row["index"]
        role = primary_roles.get(index, "unselected_native_track")
        if index in primary_identifiers and row["identifier"] != primary_identifiers[index]:
            raise BorzoiFixtureError(f"primary native track identity differs: {index}")
        rendered.append({**row, "assay": assay, "native_screen_role": role})
    if dict(assay_counts) != output["assay_track_counts"]:
        raise BorzoiFixtureError("native assay track counts differ")
    selected = [int(row["index"]) for row in rendered if row["native_screen_role"] != "unselected_native_track"]
    if selected != output["primary_track_indices"]:
        raise BorzoiFixtureError("primary native track order differs")
    if any(
        strand_pairs[index] not in output["primary_track_indices"]
        for index in output["primary_track_indices"]
    ):
        raise BorzoiFixtureError("primary track set is not strand-pair closed")
    return rendered


def build_long_range_groups(
    source_groups: list[dict[str, str]], *, buffer_bp: int, folds: int, seed: str
) -> tuple[list[dict[str, object]], dict[str, tuple[str, int]]]:
    by_contig: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in source_groups:
        by_contig[row["contig"]].append(row)
    components: list[dict[str, object]] = []
    for contig in sorted(by_contig):
        ordered = sorted(
            by_contig[contig], key=lambda row: int(row["minimum_variant_pos0"])
        )
        current: list[dict[str, str]] = []
        current_end = -1
        for row in ordered:
            start = int(row["minimum_variant_pos0"])
            end = int(row["maximum_variant_pos0"])
            if current and start - current_end >= buffer_bp:
                components.append(component_record(contig, current))
                current = []
                current_end = -1
            current.append(row)
            current_end = max(current_end, end)
        if current:
            components.append(component_record(contig, current))

    component_totals = [0] * folds
    contig_totals: dict[str, list[int]] = defaultdict(lambda: [0] * folds)
    ordered_components = sorted(
        components,
        key=lambda row: (
            -int(row["source_locus_groups"]),
            -int(row["source_elements"]),
            sha256(f"{seed}\0{row['borzoi_long_range_group_id']}".encode()).hexdigest(),
        ),
    )
    assignment: dict[str, int] = {}
    for row in ordered_components:
        contig = str(row["contig"])
        weight = int(row["source_locus_groups"])
        fold = min(
            range(folds),
            key=lambda candidate: (
                contig_totals[contig][candidate],
                component_totals[candidate],
                sha256(
                    f"{seed}\0{row['borzoi_long_range_group_id']}\0{candidate}".encode()
                ).hexdigest(),
            ),
        )
        group_id = str(row["borzoi_long_range_group_id"])
        assignment[group_id] = fold
        contig_totals[contig][fold] += weight
        component_totals[fold] += weight
        row["outer_fold"] = fold
    if set(assignment.values()) != set(range(folds)):
        raise BorzoiFixtureError("a long-range outer fold is empty")

    source_map: dict[str, tuple[str, int]] = {}
    for component in components:
        component_id = str(component["borzoi_long_range_group_id"])
        for source_group in component.pop("_member_ids"):
            if source_group in source_map:
                raise BorzoiFixtureError("source group maps to multiple long-range groups")
            source_map[source_group] = (component_id, assignment[component_id])
    verify_long_range_firewall(components, buffer_bp)
    return sorted(components, key=lambda row: str(row["borzoi_long_range_group_id"])), source_map


def component_record(contig: str, rows: list[dict[str, str]]) -> dict[str, object]:
    member_ids = sorted(row["outer_locus_sequence_group_id"] for row in rows)
    identifier = "borzoi_lr_" + sha256("\n".join(member_ids).encode()).hexdigest()[:20]
    return {
        "borzoi_long_range_group_id": identifier,
        "outer_fold": -1,
        "contig": contig,
        "minimum_variant_pos0": min(int(row["minimum_variant_pos0"]) for row in rows),
        "maximum_variant_pos0": max(int(row["maximum_variant_pos0"]) for row in rows),
        "source_locus_groups": len(rows),
        "source_elements": sum(int(row["elements"]) for row in rows),
        "_member_ids": member_ids,
    }


def verify_long_range_firewall(components: list[dict[str, object]], buffer_bp: int) -> None:
    by_contig: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in components:
        by_contig[str(row["contig"])].append(row)
    for rows in by_contig.values():
        ordered = sorted(rows, key=lambda row: int(row["minimum_variant_pos0"]))
        for left, right in zip(ordered, ordered[1:]):
            distance = int(right["minimum_variant_pos0"]) - int(left["maximum_variant_pos0"])
            if distance < buffer_bp:
                raise BorzoiFixtureError("long-range components violate receptive-field buffer")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def allele_window(
    reference: IndexedFasta,
    row: Mapping[str, str],
    *,
    window: int,
    center: int,
    allowed_contigs: set[str],
) -> tuple[str, str, int, int]:
    contig = row["contig"]
    if contig not in allowed_contigs:
        raise BorzoiFixtureError("variant is not on an admitted primary contig")
    position0 = int(row["variant_pos0"])
    if int(row["variant_pos1"]) != position0 + 1:
        raise BorzoiFixtureError("0-based and 1-based coordinates differ")
    start, end = position0 - center, position0 - center + window
    if start < 0 or end > reference.contig_length(contig):
        raise BorzoiFixtureError("Borzoi window crosses a contig boundary")
    sequence = reference.fetch(contig, start, end)
    ref, alt = row["genomic_ref"], row["genomic_alt"]
    if (
        len(sequence) != window
        or set(sequence) - set("ACGTN")
        or sequence[center] != ref
        or ref not in "ACGT"
        or alt not in "ACGT"
        or ref == alt
    ):
        raise BorzoiFixtureError("reference sequence, REF, or ALT differs")
    alternative = sequence[:center] + alt + sequence[center + 1 :]
    if sum(left != right for left, right in zip(sequence, alternative)) != 1:
        raise BorzoiFixtureError("alternative sequence differs at more than one base")
    return sequence, alternative, start, end


def build(*, project_root: Path, contract_path: Path, output: Path) -> dict[str, object]:
    if output.exists():
        raise BorzoiFixtureError("output exists")
    project_root = project_root.resolve(strict=True)
    contract = load_contract(contract_path)
    roots = validate_authorities(project_root, contract)
    sequence_contract = contract["sequence_contract"]
    firewall = contract["leakage_firewall"]
    split_root = roots["outcome_blind_split"]
    selection_root = roots["common_locus_selection"]

    source_groups = read_tsv(split_root / "split/groups.tsv")
    split_rows = {
        row["element_id"]: row
        for row in read_tsv(split_root / "split/elements.tsv")
    }
    selection_rows = read_tsv(selection_root / "fixture/manifest.tsv")
    if (
        len(source_groups) != firewall["expected_source_locus_groups"]
        or len(split_rows) != 4_359
        or len(selection_rows) != firewall["expected_source_locus_groups"]
        or len({row["outer_locus_sequence_group_id"] for row in selection_rows})
        != firewall["expected_source_locus_groups"]
    ):
        raise BorzoiFixtureError("outcome-blind source census differs")
    selected: list[dict[str, str]] = []
    for selection_row in selection_rows:
        row = split_rows.get(selection_row["element_id"])
        if row is None or any(
            row[field] != selection_row[field]
            for field in (
                "outer_locus_sequence_group_id",
                "contig",
                "variant_pos0",
            )
        ):
            raise BorzoiFixtureError("common outcome-blind locus selection differs")
        selected.append(row)

    components, source_map = build_long_range_groups(
        source_groups,
        buffer_bp=int(firewall["borzoi_receptive_field_buffer_bp"]),
        folds=int(firewall["outer_folds"]),
        seed=str(firewall["fold_seed"]),
    )
    if (
        len(components) != firewall["expected_long_range_components"]
        or set(source_map)
        != {row["outer_locus_sequence_group_id"] for row in source_groups}
    ):
        raise BorzoiFixtureError("long-range component census differs")

    reference = IndexedFasta(
        roots["reference_fasta"], Path(f"{roots['reference_fasta']}.fai")
    )
    fixture_rows: list[dict[str, object]] = []
    for row in selected:
        sequence, alternative, start, end = allele_window(
            reference,
            row,
            window=int(sequence_contract["input_length_bp"]),
            center=int(sequence_contract["variant_index0"]),
            allowed_contigs=set(sequence_contract["allowed_contigs"]),
        )
        source_group = row["outer_locus_sequence_group_id"]
        long_group, fold = source_map[source_group]
        fixture_rows.append(
            {
                "fixture_id": "borzoi_" + sha256(row["element_id"].encode()).hexdigest()[:24],
                "element_id": row["element_id"],
                "source_locus_group_id": source_group,
                "borzoi_long_range_group_id": long_group,
                "outer_fold": fold,
                "contig": row["contig"],
                "variant_pos0": row["variant_pos0"],
                "variant_pos1": row["variant_pos1"],
                "input_start0": start,
                "input_end0": end,
                "variant_index0": sequence_contract["variant_index0"],
                "ref": row["genomic_ref"],
                "alt": row["genomic_alt"],
                "reference_sequence_sha256": sha256(sequence.encode()).hexdigest(),
                "alternative_sequence_sha256": sha256(alternative.encode()).hexdigest(),
                "reference_reverse_complement_sha256": sha256(
                    reverse_complement(sequence).encode()
                ).hexdigest(),
                "alternative_reverse_complement_sha256": sha256(
                    reverse_complement(alternative).encode()
                ).hexdigest(),
                "allele_effect_sign": sequence_contract["allele_effect_sign"],
            }
        )

    track_rows = load_native_tracks(roots["target_manifest"], contract)
    source_by_id = {row["outer_locus_sequence_group_id"]: row for row in source_groups}
    source_group_rows = []
    for source_group in sorted(source_map):
        row = source_by_id[source_group]
        long_group, fold = source_map[source_group]
        source_group_rows.append(
            {
                "source_locus_group_id": source_group,
                "borzoi_long_range_group_id": long_group,
                "outer_fold": fold,
                "contig": row["contig"],
                "minimum_variant_pos0": row["minimum_variant_pos0"],
                "maximum_variant_pos0": row["maximum_variant_pos0"],
                "source_elements": row["elements"],
            }
        )

    output.mkdir(parents=True)
    write_tsv(output / "manifest.tsv", FIXTURE_FIELDS, fixture_rows)
    write_tsv(output / "source_group_map.tsv", SOURCE_GROUP_FIELDS, source_group_rows)
    write_tsv(output / "long_range_components.tsv", COMPONENT_FIELDS, components)
    write_tsv(output / "native_tracks.tsv", TRACK_FIELDS, track_rows)
    fold_loci = {
        str(fold): sum(int(row["outer_fold"]) == fold for row in fixture_rows)
        for fold in range(int(firewall["outer_folds"]))
    }
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-gse281364-borzoi-native-fixture-v1",
        "status": "pass_fixture_checkpoint_execution_blocked",
        "dataset_id": "gse281364",
        "model_id": "borzoi_ensemble",
        "fixture_elements": len(fixture_rows),
        "source_locus_groups": len(source_groups),
        "borzoi_long_range_groups": len(components),
        "outer_folds": firewall["outer_folds"],
        "fold_fixture_elements": fold_loci,
        "receptive_field_buffer_bp": firewall["borzoi_receptive_field_buffer_bp"],
        "input_length_bp": sequence_contract["input_length_bp"],
        "variant_index0": sequence_contract["variant_index0"],
        "native_human_tracks": len(track_rows),
        "primary_native_tracks": len(contract["native_output_contract"]["primary_track_indices"]),
        "official_ensemble_members": 4,
        "orientation_predictions_per_allele": 8,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_outcomes_read": False,
        "checkpoint_bytes_loaded": False,
        "model_forward_executed": False,
        "model_metrics_calculated": False,
        "weight_terms": "UNDECLARED_BLOCKED",
        "champion_eligible": False,
        "claim_role": "exposed_development_MPRA_sequence_native_diagnostic_only",
        "mandatory_baselines": contract["mandatory_baselines"],
        "contract_sha256": digest(contract_path),
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--contract", dest="contract_path", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    print(json.dumps(build(**vars(parser.parse_args(argv))), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
