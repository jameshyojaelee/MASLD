#!/usr/bin/env python3
"""Materialize outcome-blind LS-GKM exact-dinucleotide-null inputs."""

from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import defaultdict
import csv
from dataclasses import dataclass
import gzip
import hashlib
import io
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from scripts.freeze_lsgkm_gse281364_dinucleotide_null_readiness import (
    DNA,
    SEEDS,
    canonical_sequence,
    exact_dinucleotide_shuffle,
    file_sha256,
)


SCHEMA = "masld-bench-lsgkm-gse281364-dinucleotide-materialization-input-v1"
CONTEXTS = {"HepG2_control", "HepG2_PAOA"}
PAIR_FIELDS = (
    "split_id",
    "seed",
    "pair_index",
    "positive_id",
    "negative_id",
    "contig",
    "start0",
    "end0",
    "positive_sequence_sha256",
    "negative_sequence_sha256",
)
SCORING_FIELDS = (
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
    "outer_fold",
    "scoring_split_id",
    "genomic_test_fold",
    "contig",
    "variant_pos0",
    "variant_index0",
    "ref",
    "alt",
    "reference_sequence_sha256",
    "alternative_sequence_sha256",
    "allele_effect_sign",
)


class DinucleotideMaterializationError(RuntimeError):
    """Raised when source, split, sequence, or output invariants differ."""


@dataclass(frozen=True)
class PeakWindow:
    identifier: str
    contig: str
    start: int
    end: int
    sequence: str

    @property
    def sequence_sha256(self) -> str:
        return hashlib.sha256(self.sequence.encode("ascii")).hexdigest()

    @property
    def canonical(self) -> str:
        return canonical_sequence(self.sequence)


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise DinucleotideMaterializationError(f"JSON object differs: {path}")
    return value


def resolve_file(root: Path, relative_text: str, expected_sha256: str) -> Path:
    relative = Path(relative_text)
    if relative.is_absolute() or ".." in relative.parts:
        raise DinucleotideMaterializationError("unsafe project-relative file")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    if not path.is_file() or path.is_symlink() or file_sha256(path) != expected_sha256:
        raise DinucleotideMaterializationError(f"file authority differs: {relative_text}")
    return path


def manifest_member_map(manifest: Mapping[str, Any]) -> dict[str, str]:
    members = manifest.get("artifacts")
    if not isinstance(members, list):
        raise DinucleotideMaterializationError("artifact manifest differs")
    result: dict[str, str] = {}
    for member in members:
        if not isinstance(member, dict) or set(member) != {"path", "sha256", "size_bytes"}:
            raise DinucleotideMaterializationError("artifact member differs")
        result[str(member["path"])] = str(member["sha256"])
    return result


def verify_authority(
    root: Path,
    binding: Mapping[str, Any],
    member_keys: Sequence[tuple[str, str]],
) -> dict[str, Any]:
    artifacts = resolve_file(root, str(binding["artifacts_path"]), str(binding["artifacts_sha256"]))
    manifest = load_json(artifacts)
    if manifest.get("schema_version") != "masld-bench-artifacts-v1":
        raise DinucleotideMaterializationError("artifact schema differs")
    roster = manifest_member_map(manifest)
    authority_root = artifacts.parent
    for path_key, sha_key in member_keys:
        member = resolve_file(root, str(binding[path_key]), str(binding[sha_key]))
        try:
            relative = member.relative_to(authority_root).as_posix()
        except ValueError as error:
            raise DinucleotideMaterializationError("member outside artifact authority") from error
        if roster.get(relative) != binding[sha_key]:
            raise DinucleotideMaterializationError(f"member absent from authority: {relative}")
    return manifest


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    if not rows:
        raise DinucleotideMaterializationError(f"TSV is empty: {path}")
    return rows


def write_tsv_gzip(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> int:
    if path.exists():
        raise DinucleotideMaterializationError(f"refusing to overwrite: {path}")
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


def write_fasta_gzip(path: Path, records: Iterable[tuple[str, str]]) -> int:
    if path.exists():
        raise DinucleotideMaterializationError(f"refusing to overwrite: {path}")
    count = 0
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
            for identifier, sequence in records:
                compressed.write(f">{identifier}\n{sequence}\n".encode("ascii"))
                count += 1
    return count


def load_interval_index(path: Path) -> dict[str, tuple[list[int], list[int]]]:
    grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise DinucleotideMaterializationError("blacklist interval differs")
            grouped[fields[0]].append((start, end))
    result: dict[str, tuple[list[int], list[int]]] = {}
    for contig, intervals in grouped.items():
        merged: list[list[int]] = []
        for start, end in sorted(intervals):
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        result[contig] = ([row[0] for row in merged], [row[1] for row in merged])
    return result


def overlaps(
    index: Mapping[str, tuple[list[int], list[int]]], contig: str, start: int, end: int
) -> bool:
    starts, ends = index.get(contig, ([], []))
    position = bisect_right(ends, start)
    return position < len(starts) and starts[position] < end


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status")
        != "outcome_blind_input_materialization_authorized_training_blocked"
        or config.get("design_id")
        != "lsgkm_hepatocyte_atac_peak_exact_dinucleotide_null_v1"
    ):
        raise DinucleotideMaterializationError("materialization config identity differs")
    firewall = config["action_firewall"]
    if firewall.get("production_input_materialization_authorized") is not True or any(
        firewall.get(field) is not False
        for field in (
            "outcome_access_authorized",
            "reporter_count_access_authorized",
            "sealed_asset_access_authorized",
            "prediction_value_access_authorized",
            "production_training_authorized",
            "production_prediction_authorized",
        )
    ):
        raise DinucleotideMaterializationError("materialization action firewall differs")
    design = config["design"]
    if (
        design.get("window_bp") != 300
        or design.get("variant_index0") != 150
        or design.get("positive_count_per_split") != 10000
        or design.get("negative_count_per_split_seed") != 10000
        or design.get("fixed_seeds") != SEEDS
        or design.get("max_shuffle_attempts_per_positive_seed") != 256
        or design.get("shared_positive_roster_across_seeds") is not True
        or design.get("exact_monomer_and_dinucleotide_preservation") is not True
        or design.get("canonical_exact_and_reverse_complement_collision_rejection") is not True
        or design.get("full_test_and_validation_chromosome_exclusion") is not True
        or design.get("emitted_FASTA_independent_invariant_audit_required") is not True
        or design.get("shared_fit_count") != 25
        or design.get("independent_model_families") != 1
    ):
        raise DinucleotideMaterializationError("materialization design differs")
    terminal = config["terminal_gate"]
    if (
        terminal.get("minimum_selected_positive_windows_every_split") != 10000
        or terminal.get("exact_pairs_every_split_seed") != 10000
        or terminal.get("all_25_fit_inputs_required") is not True
    ):
        raise DinucleotideMaterializationError("materialization terminal gate differs")


def fold_contract(
    genomic_path: Path, crossed_path: Path
) -> tuple[dict[str, int], dict[int, dict[str, str]]]:
    fold_by_contig: dict[str, int] = {}
    for row in read_tsv(genomic_path):
        fold = int(row["genomic_fold"])
        for contig in row["contigs"].split(","):
            if contig in fold_by_contig:
                raise DinucleotideMaterializationError("duplicate genomic-fold contig")
            fold_by_contig[contig] = fold
    diagonal: dict[int, dict[str, str]] = {}
    for row in read_tsv(crossed_path):
        if (
            row["donor_test_fold"] == row["genomic_test_fold"]
            and row["split_id"]
            == f"donor{row['donor_test_fold']}_genomic{row['genomic_test_fold']}"
        ):
            diagonal[int(row["genomic_test_fold"])] = row
    if len(fold_by_contig) != 24 or set(diagonal) != set(range(5)):
        raise DinucleotideMaterializationError("genomic split contract differs")
    return fold_by_contig, diagonal


def build_scoring_sequences(
    *,
    fasta: Any,
    row_path: Path,
    allele_path: Path,
    fold_by_contig: Mapping[str, int],
    diagonal: Mapping[int, Mapping[str, str]],
) -> tuple[list[dict[str, object]], list[tuple[str, str]], set[str]]:
    row_groups: dict[str, tuple[str, str, str]] = {}
    seeds: set[int] = set()
    contexts: set[str] = set()
    for row in read_tsv(row_path):
        if row["study_id"] != "gse281364" or row["stratum"] != "all":
            raise DinucleotideMaterializationError("scoring row identity differs")
        seeds.add(int(row["seed"]))
        contexts.add(row["assay_context_id"])
        identity = (
            row["source_locus_group_id"],
            row["long_range_block_id"],
            row["outer_fold"],
        )
        prior = row_groups.setdefault(row["element_id"], identity)
        if prior != identity:
            raise DinucleotideMaterializationError("scoring row group differs")
    if len(row_groups) != 1033 or seeds != set(SEEDS) or contexts != CONTEXTS:
        raise DinucleotideMaterializationError("scoring row universe differs")

    allele_rows = read_tsv(allele_path)
    if len(allele_rows) != 1033 or len({row["element_id"] for row in allele_rows}) != 1033:
        raise DinucleotideMaterializationError("scoring allele universe differs")
    output_rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str]] = []
    canonical: set[str] = set()
    for source in sorted(allele_rows, key=lambda row: row["element_id"]):
        element = source["element_id"]
        row_group = row_groups.get(element)
        if row_group is None:
            raise DinucleotideMaterializationError("scoring allele missing from rows")
        contig = source["contig"]
        position = int(source["variant_pos0"])
        ref, alt = source["ref"], source["alt"]
        start, end = position - 150, position + 150
        reference = str(fasta[contig][start:end]).upper()
        if (
            len(reference) != 300
            or set(reference) - DNA
            or len(ref) != 1
            or len(alt) != 1
            or reference[150] != ref
        ):
            raise DinucleotideMaterializationError(f"scoring REF/ALT geometry differs: {element}")
        alternative = reference[:150] + alt + reference[151:]
        test_fold = fold_by_contig[contig]
        if row_group != (
            source["source_locus_group_id"],
            source["borzoi_long_range_group_id"],
            f"fold-{source['outer_fold']}",
        ):
            raise DinucleotideMaterializationError("scoring long-range fold join differs")
        split_id = diagonal[test_fold]["split_id"]
        ref_hash = hashlib.sha256(reference.encode("ascii")).hexdigest()
        alt_hash = hashlib.sha256(alternative.encode("ascii")).hexdigest()
        fasta_records.extend(((f"{element}|REF", reference), (f"{element}|ALT", alternative)))
        canonical.update((canonical_sequence(reference), canonical_sequence(alternative)))
        output_rows.append(
            {
                "element_id": element,
                "source_locus_group_id": row_group[0],
                "long_range_block_id": row_group[1],
                "outer_fold": row_group[2],
                "scoring_split_id": split_id,
                "genomic_test_fold": test_fold,
                "contig": contig,
                "variant_pos0": position,
                "variant_index0": 150,
                "ref": ref,
                "alt": alt,
                "reference_sequence_sha256": ref_hash,
                "alternative_sequence_sha256": alt_hash,
                "allele_effect_sign": "ALT_minus_REF",
            }
        )
    return output_rows, fasta_records, canonical


def read_peak_candidates(
    *,
    fasta: Any,
    peak_path: Path,
    train_contigs: set[str],
    blacklist: Mapping[str, tuple[list[int], list[int]]],
    scoring_canonical: set[str],
) -> tuple[list[PeakWindow], dict[str, int]]:
    raw = 0
    excluded_nontraining = 0
    excluded_blacklist = 0
    excluded_sequence = 0
    excluded_scoring_collision = 0
    deduplicated = 0
    by_canonical: dict[str, PeakWindow] = {}
    with peak_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            raw += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 10:
                raise DinucleotideMaterializationError("training narrowPeak row differs")
            contig = fields[0]
            center = int(fields[1]) + int(fields[9])
            start, end = center - 150, center + 150
            if contig not in train_contigs:
                excluded_nontraining += 1
                continue
            if start < 0 or overlaps(blacklist, contig, start, end):
                excluded_blacklist += 1
                continue
            sequence = str(fasta[contig][start:end]).upper()
            if len(sequence) != 300 or set(sequence) - DNA:
                excluded_sequence += 1
                continue
            canonical = canonical_sequence(sequence)
            if canonical in scoring_canonical:
                excluded_scoring_collision += 1
                continue
            candidate = PeakWindow(
                identifier=f"peak_{line_number:06d}",
                contig=contig,
                start=start,
                end=end,
                sequence=sequence,
            )
            if canonical in by_canonical:
                deduplicated += 1
            else:
                by_canonical[canonical] = candidate
    candidates = sorted(
        by_canonical.values(), key=lambda peak: (peak.contig, peak.start, peak.identifier)
    )
    return candidates, {
        "raw_training_peaks": raw,
        "eligible_unique_peak_windows": len(candidates),
        "excluded_nontraining_contig": excluded_nontraining,
        "excluded_blacklist_or_edge": excluded_blacklist,
        "excluded_non_ACGT_or_incomplete": excluded_sequence,
        "excluded_scoring_allele_collision": excluded_scoring_collision,
        "deduplicated_exact_or_reverse_complement": deduplicated,
    }


def candidate_order_key(design_id: str, split_id: str, peak: PeakWindow) -> tuple[str, str]:
    material = (
        f"{design_id}|{split_id}|{peak.identifier}|{peak.sequence_sha256}"
    ).encode("utf-8")
    return hashlib.sha256(material).hexdigest(), peak.identifier


def select_rectangle(
    *,
    candidates: Sequence[PeakWindow],
    scoring_canonical: set[str],
    design_id: str,
    split_id: str,
    selected_count: int,
    max_attempts: int,
) -> tuple[list[PeakWindow], dict[int, list[str]], int]:
    eligible_candidates = [
        peak for peak in candidates if peak.canonical not in scoring_canonical
    ]
    all_positive_canonical = {peak.canonical for peak in eligible_candidates}
    used_by_seed = {
        seed: set(all_positive_canonical).union(scoring_canonical) for seed in SEEDS
    }
    selected: list[PeakWindow] = []
    negatives: dict[int, list[str]] = {seed: [] for seed in SEEDS}
    rejected_unshufflable = 0
    for peak in sorted(
        eligible_candidates,
        key=lambda value: candidate_order_key(design_id, split_id, value),
    ):
        local: set[str] = set()
        proposals: dict[int, str] = {}
        for seed in SEEDS:
            for attempt in range(max_attempts):
                negative = exact_dinucleotide_shuffle(
                    peak.sequence,
                    design_id=design_id,
                    split_id=split_id,
                    positive_id=peak.identifier,
                    model_seed=seed,
                    attempt=attempt,
                )
                canonical = canonical_sequence(negative)
                if canonical not in used_by_seed[seed] and canonical not in local:
                    proposals[seed] = negative
                    local.add(canonical)
                    break
            else:
                proposals = {}
                break
        if not proposals:
            rejected_unshufflable += 1
            continue
        selected.append(peak)
        for seed, negative in proposals.items():
            negatives[seed].append(negative)
            used_by_seed[seed].add(canonical_sequence(negative))
        if len(selected) == selected_count:
            break
    if len(selected) != selected_count or any(
        len(negatives[seed]) != selected_count for seed in SEEDS
    ):
        raise DinucleotideMaterializationError(
            f"terminal exact-dinucleotide rectangle gate failed: {split_id} selected={len(selected)}"
        )
    return selected, negatives, rejected_unshufflable


def build(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.project_root.resolve(strict=True)
    config = load_json(arguments.config.resolve(strict=True))
    validate_config(config)
    output = arguments.output
    if output.exists():
        raise DinucleotideMaterializationError("materialization output exists")
    (output / "inputs").mkdir(parents=True, mode=0o750)
    (output / "scoring").mkdir(mode=0o750)
    (output / "contract").mkdir(mode=0o750)

    readiness_artifacts = resolve_file(
        root,
        config["readiness"]["artifacts_path"],
        config["readiness"]["artifacts_sha256"],
    )
    readiness = load_json(
        resolve_file(
            root,
            config["readiness"]["receipt_path"],
            config["readiness"]["receipt_sha256"],
        )
    )
    if (
        readiness.get("status")
        != "ready_for_outcome_blind_input_materialization_and_scale_probe_only"
        or readiness.get("design_id") != config["design_id"]
        or readiness.get("outcomes_read") is not False
        or readiness.get("prediction_values_read") is not False
    ):
        raise DinucleotideMaterializationError("readiness receipt differs")

    sources = config["source_authorities"]
    verify_authority(
        root,
        sources["sequence_fasta"],
        (("member_path", "member_sha256"), ("index_path", "index_sha256")),
    )
    verify_authority(
        root,
        sources["sequence_split"],
        (
            ("genomic_folds_path", "genomic_folds_sha256"),
            ("crossed_splits_path", "crossed_splits_sha256"),
        ),
    )
    verify_authority(root, sources["blacklist"], (("member_path", "member_sha256"),))
    verify_authority(root, sources["scoring_rows"], (("member_path", "member_sha256"),))
    verify_authority(root, sources["scoring_alleles"], (("member_path", "member_sha256"),))
    for training in config["training_sources"]:
        verify_authority(
            root,
            training,
            (("peak_path", "peak_sha256"), ("fold_path", "fold_sha256")),
        )

    fasta_path = root / sources["sequence_fasta"]["member_path"]
    fold_by_contig, diagonal = fold_contract(
        root / sources["sequence_split"]["genomic_folds_path"],
        root / sources["sequence_split"]["crossed_splits_path"],
    )
    blacklist = load_interval_index(root / sources["blacklist"]["member_path"])

    from pyfaidx import Fasta

    fasta = Fasta(str(fasta_path), as_raw=True, sequence_always_upper=True, rebuild=False)
    try:
        scoring_rows, scoring_records, scoring_canonical = build_scoring_sequences(
            fasta=fasta,
            row_path=root / sources["scoring_rows"]["member_path"],
            allele_path=root / sources["scoring_alleles"]["member_path"],
            fold_by_contig=fold_by_contig,
            diagonal=diagonal,
        )
        write_tsv_gzip(
            output / "contract/scoring_elements.tsv.gz", SCORING_FIELDS, scoring_rows
        )
        write_fasta_gzip(output / "scoring/alleles.fa.gz", scoring_records)

        split_summaries: list[dict[str, Any]] = []
        for training in config["training_sources"]:
            split_id = training["split_id"]
            test_fold = int(split_id[-1])
            split = diagonal[test_fold]
            fold_json = load_json(root / training["fold_path"])
            train_contigs = set(fold_json["train"])
            valid_contigs = set(fold_json["valid"])
            test_contigs = set(fold_json["test"])
            expected_train = {
                contig
                for contig, fold in fold_by_contig.items()
                if str(fold) in split["genomic_train_folds"].split(",")
            }
            expected_valid = {
                contig
                for contig, fold in fold_by_contig.items()
                if fold == int(split["genomic_valid_fold"])
            }
            expected_test = {
                contig for contig, fold in fold_by_contig.items() if fold == test_fold
            }
            if (
                train_contigs != expected_train
                or valid_contigs != expected_valid
                or test_contigs != expected_test
                or train_contigs & (valid_contigs | test_contigs)
            ):
                raise DinucleotideMaterializationError(f"outer-safe fold differs: {split_id}")
            candidates, census = read_peak_candidates(
                fasta=fasta,
                peak_path=root / training["peak_path"],
                train_contigs=train_contigs,
                blacklist=blacklist,
                scoring_canonical=scoring_canonical,
            )
            selected, negatives, rejected = select_rectangle(
                candidates=candidates,
                scoring_canonical=scoring_canonical,
                design_id=config["design_id"],
                split_id=split_id,
                selected_count=config["design"]["positive_count_per_split"],
                max_attempts=config["design"]["max_shuffle_attempts_per_positive_seed"],
            )
            positive_path = output / "inputs" / f"{split_id}.positive.fa.gz"
            write_fasta_gzip(
                positive_path, ((peak.identifier, peak.sequence) for peak in selected)
            )
            seed_summaries: dict[str, dict[str, Any]] = {}
            for seed in SEEDS:
                negative_path = output / "inputs" / f"{split_id}.seed{seed}.negative.fa.gz"
                pair_path = output / "inputs" / f"{split_id}.seed{seed}.pairs.tsv.gz"
                write_fasta_gzip(
                    negative_path,
                    (
                        (f"{selected[index].identifier}|shuffle_seed{seed}", sequence)
                        for index, sequence in enumerate(negatives[seed])
                    ),
                )
                pair_rows = []
                for index, (peak, negative) in enumerate(zip(selected, negatives[seed])):
                    pair_rows.append(
                        {
                            "split_id": split_id,
                            "seed": seed,
                            "pair_index": index,
                            "positive_id": peak.identifier,
                            "negative_id": f"{peak.identifier}|shuffle_seed{seed}",
                            "contig": peak.contig,
                            "start0": peak.start,
                            "end0": peak.end,
                            "positive_sequence_sha256": peak.sequence_sha256,
                            "negative_sequence_sha256": hashlib.sha256(
                                negative.encode("ascii")
                            ).hexdigest(),
                        }
                    )
                write_tsv_gzip(pair_path, PAIR_FIELDS, pair_rows)
                seed_summaries[str(seed)] = {
                    "pair_count": len(pair_rows),
                    "negative_fasta_sha256": file_sha256(negative_path),
                    "pair_manifest_sha256": file_sha256(pair_path),
                }
            split_summaries.append(
                {
                    "split_id": split_id,
                    "genomic_test_fold": test_fold,
                    "donor_test_fold": int(split["donor_test_fold"]),
                    "train_contigs": sorted(train_contigs),
                    "held_validation_contigs": sorted(valid_contigs),
                    "held_test_contigs": sorted(test_contigs),
                    "source_census": census,
                    "candidates_examined_through_selection": len(selected) + rejected,
                    "rejected_unshufflable": rejected,
                    "selected_positive_windows": len(selected),
                    "positive_fasta_sha256": file_sha256(positive_path),
                    "seeds": seed_summaries,
                }
            )
    finally:
        fasta.close()

    contract = {
        "schema_version": "masld-bench-lsgkm-gse281364-dinucleotide-materialization-v1",
        "status": "inputs_materialized_pending_independent_emitted_fasta_audit",
        "design_id": config["design_id"],
        "readiness_artifacts_sha256": file_sha256(readiness_artifacts),
        "config_sha256": file_sha256(arguments.config),
        "scientific_estimand": readiness["scientific_estimand"],
        "not_an_estimand": readiness["not_an_estimand"],
        "split_count": 5,
        "shared_fit_count": 25,
        "positive_FASTA_count": 5,
        "negative_FASTA_count": 25,
        "pair_manifest_count": 25,
        "positive_roster_shared_across_seeds_within_split": True,
        "scoring_universe": {
            "elements": len(scoring_rows),
            "long_range_blocks": len({row["long_range_block_id"] for row in scoring_rows}),
            "contexts": sorted(CONTEXTS),
            "seeds": SEEDS,
            "standardized_rows_per_readout": 10330,
            "scoring_element_manifest_sha256": file_sha256(
                output / "contract/scoring_elements.tsv.gz"
            ),
            "scoring_allele_fasta_sha256": file_sha256(output / "scoring/alleles.fa.gz"),
            "allele_effect_sign": "ALT_minus_REF",
        },
        "split_summaries": split_summaries,
        "runtime_fit_contract": readiness["runtime_capabilities"],
        "readouts": ["direct_gkmsvm", "deltasvm"],
        "readouts_are_independent_families": False,
        "claim_limits": config["claim_limits"],
        "independent_emitted_FASTA_audit_passed": False,
        "production_training_authorized": False,
        "production_fits_executed": 0,
        "production_predictions_generated": 0,
        "benchmark_metrics_calculated": False,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "prediction_values_read": False,
        "sealed_assets_read": False,
    }
    (output / "contract/materialization_contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(contract, sort_keys=True))
    return contract


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    build(arguments)


if __name__ == "__main__":
    main()
