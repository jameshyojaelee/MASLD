#!/usr/bin/env python3
"""Validate reconstructed GSE281364 oligos and assemble replicate-safe MPRA data."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import re
from typing import Iterable


INTERVAL = re.compile(r"^(chr(?:[0-9]+|X|Y)):(\d+)-(\d+)$")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class ReconstructionValidationError(ValueError):
    """Raised when reconstructed sequences or replicate topology differ."""


class UnionFind:
    def __init__(self, values: Iterable[str]) -> None:
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if left_root < right_root:
            self.parent[right_root] = left_root
        else:
            self.parent[left_root] = right_root


class IndexedFasta:
    """Minimal read-only index for an uncompressed, fixed-line-width FASTA."""

    def __init__(self, path: Path) -> None:
        self.handle = path.open("rb")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        name: str | None = None
        sequence_offset = 0
        line_bases = 0
        line_bytes = 0
        length = 0
        saw_short_line = False

        def commit() -> None:
            if name is None:
                return
            if name in self.index or line_bases == 0 or length == 0:
                raise ReconstructionValidationError("reference FASTA index differs")
            self.index[name] = (sequence_offset, line_bases, line_bytes, length)

        while True:
            line_start = self.handle.tell()
            line = self.handle.readline()
            if not line:
                break
            if line.startswith(b">"):
                commit()
                name = line[1:].split(None, 1)[0].decode("ascii")
                if not name:
                    raise ReconstructionValidationError("empty reference FASTA identifier")
                sequence_offset = self.handle.tell()
                line_bases = 0
                line_bytes = 0
                length = 0
                saw_short_line = False
                continue
            if name is None:
                raise ReconstructionValidationError("reference sequence precedes identifier")
            bases = line.rstrip(b"\r\n")
            if not bases:
                raise ReconstructionValidationError("blank reference FASTA sequence line")
            if line_bases == 0:
                line_bases = len(bases)
                line_bytes = len(line)
            elif saw_short_line or len(bases) > line_bases or (
                len(bases) == line_bases and len(line) != line_bytes
            ):
                raise ReconstructionValidationError(
                    f"reference FASTA line width differs near byte {line_start}"
                )
            elif len(bases) < line_bases:
                saw_short_line = True
            length += len(bases)
        commit()
        if not self.index:
            raise ReconstructionValidationError("reference FASTA is empty")

    def fetch(self, contig: str, start: int, end: int) -> str:
        if contig not in self.index:
            raise ReconstructionValidationError(f"reference contig is absent: {contig}")
        offset, line_bases, line_bytes, length = self.index[contig]
        if start < 0 or end < start or end > length:
            raise ReconstructionValidationError("reference fetch coordinates differ")
        seek = offset + (start // line_bases) * line_bytes + start % line_bases
        self.handle.seek(seek)
        remaining = end - start
        pieces: list[bytes] = []
        while remaining:
            piece = self.handle.readline().rstrip(b"\r\n")
            if not piece:
                raise ReconstructionValidationError("reference FASTA fetch truncated")
            take = min(remaining, len(piece))
            pieces.append(piece[:take])
            remaining -= take
        return b"".join(pieces).decode("ascii")

    def close(self) -> None:
        self.handle.close()


def _sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _sha256_text(value: str) -> str:
    return sha256(value.encode("ascii")).hexdigest()


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def _canonical_sequence(sequence: str) -> str:
    reverse = _reverse_complement(sequence)
    return min(sequence, reverse)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else path.open
    with opener(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else opener(
        "r", encoding="utf-8", newline=""
    ) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ReconstructionValidationError(f"cannot write empty table: {path.name}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _write_gzip_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ReconstructionValidationError(f"cannot write empty table: {path.name}")
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(
                    text, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)


def _unique_reference_alignment(
    fasta: IndexedFasta,
    contig: str,
    interval_start: int,
    interval_end: int,
    reference_sequence: str,
) -> tuple[str, int, int]:
    search_start = max(0, interval_start - 2)
    search_end = interval_end + 2
    genomic = fasta.fetch(contig, search_start, search_end).upper()
    candidates: list[tuple[str, int]] = []
    for orientation, query in (
        ("forward", reference_sequence),
        ("reverse", _reverse_complement(reference_sequence)),
    ):
        offset = genomic.find(query)
        while offset >= 0:
            candidates.append((orientation, search_start + offset))
            offset = genomic.find(query, offset + 1)
    candidates = sorted(set(candidates))
    if len(candidates) != 1:
        raise ReconstructionValidationError(
            f"reference alignment is not unique for {contig}:{interval_start}-{interval_end}: "
            f"{candidates!r}"
        )
    orientation, start0 = candidates[0]
    return orientation, start0, start0 + len(reference_sequence)


def _below_qc_pair_element(
    element_id: str,
    reference: dict[str, str],
    alternative: dict[str, str],
    reason: str,
) -> dict[str, object]:
    reference_sequence = reference["consensus_107bp"]
    alternative_sequence = alternative["consensus_107bp"]
    return {
        "element_id": element_id,
        "pair_state": "paired_constructs_below_sequence_qc",
        "contig": "join_unresolved",
        "interval_start_as_published": "join_unresolved",
        "interval_end_as_published": "join_unresolved",
        "reference_match_orientation": "join_unresolved",
        "reference_match_start0": "join_unresolved",
        "reference_match_end0": "join_unresolved",
        "variant_pos0": "join_unresolved",
        "variant_pos1": "join_unresolved",
        "genomic_ref": "join_unresolved",
        "genomic_alt": "join_unresolved",
        "oligo_difference_index0": "join_unresolved",
        "ref_sequence_107bp": reference_sequence,
        "alt_sequence_107bp": alternative_sequence,
        "ref_sequence_sha256": _sha256_text(reference_sequence),
        "alt_sequence_sha256": _sha256_text(alternative_sequence),
        "canonical_ref_sequence_sha256": _sha256_text(
            _canonical_sequence(reference_sequence)
        ),
        "canonical_alt_sequence_sha256": _sha256_text(
            _canonical_sequence(alternative_sequence)
        ),
        "ref_construct_barcodes": int(reference["barcode_count"]),
        "alt_construct_barcodes": int(alternative["barcode_count"]),
        "sequence_state": "below_qc",
        "sequence_qc_reason": reason,
    }


def _group_elements(elements: list[dict[str, object]]) -> None:
    paired = [row for row in elements if row["pair_state"] == "paired_snv"]
    union = UnionFind(str(row["element_id"]) for row in paired)
    by_contig: dict[str, list[dict[str, object]]] = {}
    for row in paired:
        by_contig.setdefault(str(row["contig"]), []).append(row)
    for rows in by_contig.values():
        ordered = sorted(rows, key=lambda row: int(row["variant_pos0"]))
        for left, right in zip(ordered, ordered[1:]):
            if int(right["variant_pos0"]) - int(left["variant_pos0"]) < 6000:
                union.union(str(left["element_id"]), str(right["element_id"]))
    sequence_owner: dict[str, str] = {}
    for row in paired:
        element_id = str(row["element_id"])
        for field in ("canonical_ref_sequence_sha256", "canonical_alt_sequence_sha256"):
            sequence_hash = str(row[field])
            owner = sequence_owner.setdefault(sequence_hash, element_id)
            union.union(element_id, owner)
    groups: dict[str, list[str]] = {}
    for row in paired:
        element_id = str(row["element_id"])
        groups.setdefault(union.find(element_id), []).append(element_id)
    group_ids = {
        element_id: "outer_" + _sha256_text("\n".join(sorted(members)))[:20]
        for members in groups.values()
        for element_id in members
    }
    for row in elements:
        element_id = str(row["element_id"])
        row["outer_locus_sequence_group_id"] = group_ids.get(
            element_id, "unpaired_" + _sha256_text(element_id)[:20]
        )


def _load_replicate_outcomes(
    path: Path,
    elements: dict[str, dict[str, object]],
) -> tuple[list[dict[str, object]], dict[str, list[str]]]:
    output: list[dict[str, object]] = []
    sample_contexts: dict[str, list[str]] = {}
    observed_keys: set[tuple[str, str, str]] = set()
    expected_conditions = {
        "HepG2": {"control", "PAOA"},
        "LX2": {"control", "TGFb"},
    }
    for row in _read_tsv(path):
        cell_line = row["cell_line"]
        if cell_line not in expected_conditions:
            raise ReconstructionValidationError("MPRA cell line differs")
        element_id = row["element_id"]
        if element_id not in elements:
            raise ReconstructionValidationError(
                f"replicate outcome has an unknown element: {element_id}"
            )
        condition = row["condition"]
        if condition not in expected_conditions[cell_line]:
            raise ReconstructionValidationError("MPRA condition differs")
        allele = row["allele"]
        element = elements[element_id]
        expected_alleles = (
            {"ref", "alt"}
            if element["pair_state"] in {
                "paired_snv",
                "paired_constructs_below_sequence_qc",
            }
            else {"ref"}
            if element["pair_state"] == "unpaired_reference_construct"
            else {"alt"}
        )
        if allele not in expected_alleles:
            raise ReconstructionValidationError("MPRA allele topology differs")
        replicate = int(row["experimental_replicate"])
        if replicate not in {1, 2, 3, 4}:
            raise ReconstructionValidationError("MPRA replicate index differs")
        sample_id = row["sample_id"]
        context_id = row["context_id"]
        if context_id != f"{cell_line}_{condition}":
            raise ReconstructionValidationError("MPRA context identifier differs")
        key = (element_id, allele, sample_id)
        if key in observed_keys:
            raise ReconstructionValidationError("duplicate replicate outcome")
        observed_keys.add(key)
        sample_contexts.setdefault(context_id, []).append(sample_id)
        assay_state = row["assay_state"]
        if assay_state == "observed":
            dna: int | str = int(row["DNA"])
            rna: int | str = int(row["RNA"])
            barcodes: int | str = int(row["n_barcodes"])
            if min(dna, rna, barcodes) < 0:
                raise ReconstructionValidationError("negative MPRA count")
            missing_reason = "not_applicable"
        elif assay_state == "below_qc":
            if any(row[field] != "" for field in ("DNA", "RNA", "n_barcodes")):
                raise ReconstructionValidationError("missing MPRA evidence encoded as a value")
            dna = rna = barcodes = ""
            missing_reason = row["missing_reason"]
            if missing_reason != "absent_from_source_reproduction_aggregate":
                raise ReconstructionValidationError("MPRA missing reason differs")
        else:
            raise ReconstructionValidationError("MPRA assay state differs")
        output.append(
            {
                "element_id": element_id,
                "allele": allele,
                "context_id": context_id,
                "cell_line": cell_line,
                "condition": condition,
                "experimental_replicate": replicate,
                "sample_id": sample_id,
                "DNA": dna,
                "RNA": rna,
                "n_barcodes": barcodes,
                "assay_state": assay_state,
                "missing_reason": missing_reason,
                "pairing": "same_sample_different_aliquot",
                "biological_unit": "experimental_replicate",
                "donor_id": "not_applicable",
            }
        )
    normalized_contexts = {
        context: sorted(set(samples)) for context, samples in sample_contexts.items()
    }
    if set(normalized_contexts) != {
        "HepG2_control",
        "HepG2_PAOA",
        "LX2_control",
        "LX2_TGFb",
    } or any(len(samples) != 4 for samples in normalized_contexts.values()):
        raise ReconstructionValidationError("four-context replicate topology differs")
    expected_rows = sum(
        32
        if row["pair_state"] in {
            "paired_snv",
            "paired_constructs_below_sequence_qc",
        }
        else 16
        for row in elements.values()
    )
    if len(output) != expected_rows:
        raise ReconstructionValidationError(
            f"replicate outcome census differs: {len(output)} != {expected_rows}"
        )
    return sorted(
        output,
        key=lambda row: (
            str(row["element_id"]),
            str(row["allele"]),
            str(row["context_id"]),
            int(row["experimental_replicate"]),
        ),
    ), normalized_contexts


def validate(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists() or any(
        path.is_symlink()
        for path in (
            arguments.constructs,
            arguments.reconstruction_stats,
            arguments.reference_fasta,
            arguments.replicate_outcomes,
        )
    ):
        raise ReconstructionValidationError("reconstruction validation request differs")
    arguments.output.mkdir(mode=0o750)
    stats = json.loads(arguments.reconstruction_stats.read_text())
    if (
        stats.get("status") != "pass"
        or stats.get("constructs") != 10797
        or stats.get("first_pass_paired_reads") != 71_771_807
        or stats.get("first_pass_paired_reads") != stats.get("second_pass_paired_reads")
        or stats.get("observed_count_columns_loaded")
        or stats.get("sealed_outcomes_loaded")
    ):
        raise ReconstructionValidationError("reconstruction statistics differ")
    construct_rows = _read_tsv(arguments.constructs)
    if len(construct_rows) != 10797 or len({row["construct"] for row in construct_rows}) != 10797:
        raise ReconstructionValidationError("reconstructed construct census differs")
    constructs = {row["construct"]: row for row in construct_rows}
    alternative_names = sorted(name for name in constructs if name.endswith("_Mut"))
    if len(alternative_names) != 5361:
        raise ReconstructionValidationError("reconstructed alternative census differs")

    fasta = IndexedFasta(arguments.reference_fasta)
    elements: list[dict[str, object]] = []
    paired_names: set[str] = set()
    alignment_offsets: dict[str, int] = {}
    paired_alternative_names = [
        name for name in alternative_names if name[:-4] in constructs
    ]
    alternative_only_names = sorted(set(alternative_names) - set(paired_alternative_names))
    if len(paired_alternative_names) != 5_355 or len(alternative_only_names) != 6:
        raise ReconstructionValidationError("reconstructed ref/alt join topology differs")
    for alternative_name in paired_alternative_names:
        element_id = alternative_name[:-4]
        reference = constructs[element_id]
        alternative = constructs[alternative_name]
        paired_names.update({element_id, alternative_name})
        reference_sequence = reference["consensus_107bp"]
        alternative_sequence = alternative["consensus_107bp"]
        minimum_coverage = min(
            int(reference["minimum_position_coverage"]),
            int(alternative["minimum_position_coverage"]),
        )
        minimum_fraction = min(
            float(reference["minimum_consensus_fraction"]),
            float(alternative["minimum_consensus_fraction"]),
        )
        sequences_resolved = (
            len(reference_sequence) == 107
            and len(alternative_sequence) == 107
            and not (set(reference_sequence + alternative_sequence) - set("ACGT"))
        )
        differences = (
            [
                index
                for index, (ref, alt) in enumerate(
                    zip(reference_sequence, alternative_sequence, strict=True)
                )
                if ref != alt
            ]
            if sequences_resolved
            else []
        )
        qc_reason: str | None = None
        if minimum_coverage < 20:
            qc_reason = "below_20x_minimum_position_coverage"
        elif minimum_fraction <= 0.50:
            qc_reason = "no_strict_per_base_consensus_majority"
        elif not sequences_resolved:
            qc_reason = "consensus_sequence_unresolved"
        elif len(differences) != 1:
            qc_reason = "paired_consensus_not_exactly_one_snv"
        if qc_reason is not None:
            elements.append(
                _below_qc_pair_element(
                    element_id, reference, alternative, qc_reason
                )
            )
            continue
        match = INTERVAL.fullmatch(element_id)
        if match is None:
            elements.append(
                _below_qc_pair_element(
                    element_id,
                    reference,
                    alternative,
                    "paired_construct_identifier_not_a_locus",
                )
            )
            continue
        contig, start_text, end_text = match.groups()
        interval_start = int(start_text)
        interval_end = int(end_text)
        if interval_end - interval_start != 126:
            elements.append(
                _below_qc_pair_element(
                    element_id,
                    reference,
                    alternative,
                    "published_interval_width_differs",
                )
            )
            continue
        try:
            orientation, match_start0, match_end0 = _unique_reference_alignment(
                fasta, contig, interval_start, interval_end, reference_sequence
            )
        except ReconstructionValidationError:
            elements.append(
                _below_qc_pair_element(
                    element_id,
                    reference,
                    alternative,
                    "grch38_reference_alignment_not_unique",
                )
            )
            continue
        difference_index0 = differences[0]
        if orientation == "forward":
            variant_pos0 = match_start0 + difference_index0
            genomic_ref = reference_sequence[difference_index0]
            genomic_alt = alternative_sequence[difference_index0]
        else:
            variant_pos0 = match_start0 + 106 - difference_index0
            genomic_ref = _reverse_complement(reference_sequence[difference_index0])
            genomic_alt = _reverse_complement(alternative_sequence[difference_index0])
        observed_reference = fasta.fetch(contig, variant_pos0, variant_pos0 + 1).upper()
        if observed_reference != genomic_ref:
            elements.append(
                _below_qc_pair_element(
                    element_id,
                    reference,
                    alternative,
                    "reconstructed_ref_allele_differs_from_grch38",
                )
            )
            continue
        offset_key = f"{orientation}:{match_start0 - interval_start}"
        alignment_offsets[offset_key] = alignment_offsets.get(offset_key, 0) + 1
        canonical_ref = _canonical_sequence(reference_sequence)
        canonical_alt = _canonical_sequence(alternative_sequence)
        elements.append(
            {
                "element_id": element_id,
                "pair_state": "paired_snv",
                "contig": contig,
                "interval_start_as_published": interval_start,
                "interval_end_as_published": interval_end,
                "reference_match_orientation": orientation,
                "reference_match_start0": match_start0,
                "reference_match_end0": match_end0,
                "variant_pos0": variant_pos0,
                "variant_pos1": variant_pos0 + 1,
                "genomic_ref": genomic_ref,
                "genomic_alt": genomic_alt,
                "oligo_difference_index0": difference_index0,
                "ref_sequence_107bp": reference_sequence,
                "alt_sequence_107bp": alternative_sequence,
                "ref_sequence_sha256": _sha256_text(reference_sequence),
                "alt_sequence_sha256": _sha256_text(alternative_sequence),
                "canonical_ref_sequence_sha256": _sha256_text(canonical_ref),
                "canonical_alt_sequence_sha256": _sha256_text(canonical_alt),
                "ref_construct_barcodes": int(reference["barcode_count"]),
                "alt_construct_barcodes": int(alternative["barcode_count"]),
                "sequence_state": "observed_reconstructed_from_raw_plasmid_reads",
                "sequence_qc_reason": "not_applicable",
            }
        )
    fasta.close()

    reference_only_names = sorted(
        name
        for name in constructs
        if not name.endswith("_Mut") and name not in paired_names
    )
    if len(reference_only_names) != 81:
        raise ReconstructionValidationError("unpaired-reference construct census differs")
    for element_id in reference_only_names:
        row = constructs[element_id]
        sequence = row["consensus_107bp"]
        sequence_pass = (
            int(row["minimum_position_coverage"]) >= 20
            and float(row["minimum_consensus_fraction"]) > 0.50
            and len(sequence) == 107
            and not (set(sequence) - set("ACGT"))
        )
        elements.append(
            {
                "element_id": element_id,
                "pair_state": "unpaired_reference_construct",
                "contig": "not_applicable",
                "interval_start_as_published": "not_applicable",
                "interval_end_as_published": "not_applicable",
                "reference_match_orientation": "not_applicable",
                "reference_match_start0": "not_applicable",
                "reference_match_end0": "not_applicable",
                "variant_pos0": "not_applicable",
                "variant_pos1": "not_applicable",
                "genomic_ref": "not_applicable",
                "genomic_alt": "not_applicable",
                "oligo_difference_index0": "not_applicable",
                "ref_sequence_107bp": sequence,
                "alt_sequence_107bp": "join_unresolved",
                "ref_sequence_sha256": _sha256_text(sequence),
                "alt_sequence_sha256": "join_unresolved",
                "canonical_ref_sequence_sha256": _sha256_text(
                    _canonical_sequence(sequence)
                ),
                "canonical_alt_sequence_sha256": "join_unresolved",
                "ref_construct_barcodes": int(row["barcode_count"]),
                "alt_construct_barcodes": "join_unresolved",
                "sequence_state": (
                    "observed_reconstructed_from_raw_plasmid_reads"
                    if sequence_pass
                    else "below_qc"
                ),
                "sequence_qc_reason": (
                    "not_applicable"
                    if sequence_pass
                    else "unpaired_construct_consensus_below_qc"
                ),
            }
        )
    for alternative_name in alternative_only_names:
        element_id = alternative_name[:-4]
        row = constructs[alternative_name]
        sequence = row["consensus_107bp"]
        sequence_pass = (
            int(row["minimum_position_coverage"]) >= 20
            and float(row["minimum_consensus_fraction"]) > 0.50
            and len(sequence) == 107
            and not (set(sequence) - set("ACGT"))
        )
        elements.append(
            {
                "element_id": element_id,
                "pair_state": "unpaired_alternative_construct",
                "contig": "join_unresolved",
                "interval_start_as_published": "join_unresolved",
                "interval_end_as_published": "join_unresolved",
                "reference_match_orientation": "join_unresolved",
                "reference_match_start0": "join_unresolved",
                "reference_match_end0": "join_unresolved",
                "variant_pos0": "join_unresolved",
                "variant_pos1": "join_unresolved",
                "genomic_ref": "join_unresolved",
                "genomic_alt": "join_unresolved",
                "oligo_difference_index0": "join_unresolved",
                "ref_sequence_107bp": "join_unresolved",
                "alt_sequence_107bp": sequence,
                "ref_sequence_sha256": "join_unresolved",
                "alt_sequence_sha256": _sha256_text(sequence),
                "canonical_ref_sequence_sha256": "join_unresolved",
                "canonical_alt_sequence_sha256": _sha256_text(
                    _canonical_sequence(sequence)
                ),
                "ref_construct_barcodes": "join_unresolved",
                "alt_construct_barcodes": int(row["barcode_count"]),
                "sequence_state": (
                    "observed_alt_only_ref_join_unresolved"
                    if sequence_pass
                    else "below_qc"
                ),
                "sequence_qc_reason": (
                    "pair_join_unresolved"
                    if sequence_pass
                    else "unpaired_construct_consensus_below_qc"
                ),
            }
        )
    elements = sorted(elements, key=lambda row: str(row["element_id"]))
    _group_elements(elements)
    element_map = {str(row["element_id"]): row for row in elements}
    outcomes, sample_contexts = _load_replicate_outcomes(
        arguments.replicate_outcomes, element_map
    )
    elements_path = arguments.output / "elements.tsv"
    outcomes_path = arguments.output / "replicate_outcomes.tsv.gz"
    constructs_path = arguments.output / "construct_reconstruction.tsv"
    _write_tsv(elements_path, elements)
    _write_gzip_tsv(outcomes_path, outcomes)
    standardized_constructs: list[dict[str, object]] = []
    for row in sorted(construct_rows, key=lambda value: value["construct"]):
        sequence = row["consensus_107bp"]
        standardized_constructs.append(
            {
                **row,
                "sequence_sha256": _sha256_text(sequence),
                "canonical_sequence_sha256": _sha256_text(
                    _canonical_sequence(sequence)
                ),
            }
        )
    _write_tsv(constructs_path, standardized_constructs)
    group_count = len(
        {row["outer_locus_sequence_group_id"] for row in elements if row["pair_state"] == "paired_snv"}
    )
    receipt = {
        "schema_version": "masld-bench-gse281364-mpra-reconstruction-v1",
        "status": "pass",
        "constructs": len(construct_rows),
        "paired_construct_elements": len(paired_alternative_names),
        "qualified_paired_snv_elements": sum(
            row["pair_state"] == "paired_snv" for row in elements
        ),
        "paired_construct_elements_below_sequence_qc": sum(
            row["pair_state"] == "paired_constructs_below_sequence_qc"
            for row in elements
        ),
        "reference_only_construct_elements": len(reference_only_names),
        "alternative_only_construct_elements": len(alternative_only_names),
        "exactly_one_base_difference_pairs": sum(
            row["pair_state"] == "paired_snv" for row in elements
        ),
        "grch38_ref_matched_pairs": sum(
            row["pair_state"] == "paired_snv" for row in elements
        ),
        "sequence_qc_rule": "both_constructs_minimum_position_coverage_at_least_20_strict_per_base_majority_exactly_one_snv_unique_grch38_reference_alignment",
        "sequence_qc_rule_selection_timing": "defined_after_outcome_free_raw_reconstruction_diagnostic_before_any_MPRA_model_evaluation",
        "raw_reconstruction_original_threshold_passing_constructs": stats.get(
            "passing_constructs"
        ),
        "raw_reconstruction_original_status_used_for_admission": False,
        "published_interval_width_bp": 126,
        "reconstructed_observed_oligo_bp": 107,
        "alignment_offset_orientation_census": alignment_offsets,
        "outer_locus_sequence_groups": group_count,
        "outer_group_contract": "merge_6000bp_overlapping_loci_and_exact_or_reverse_complement_sequence_identities",
        "replicate_outcome_rows": len(outcomes),
        "observed_replicate_outcome_rows": sum(
            row["assay_state"] == "observed" for row in outcomes
        ),
        "below_qc_replicate_outcome_rows": sum(
            row["assay_state"] == "below_qc" for row in outcomes
        ),
        "missing_evidence_encoded_as_zero": False,
        "contexts": sample_contexts,
        "experimental_replicates": 16,
        "replicates_are_independent_donors": False,
        "donor_count": 0,
        "pairing": "same_sample_different_aliquot_DNA_RNA",
        "elements_sha256": _sha256_file(elements_path),
        "constructs_sha256": _sha256_file(constructs_path),
        "replicate_outcomes_sha256": _sha256_file(outcomes_path),
        "reference_build": "GRCh38.p14",
        "reference_annotation": "GENCODE_v49_not_used_for_MPRA_sequence_reconstruction",
        "outcome_role": "exposed_development_MPRA_only",
        "eQTL_or_ieQTL_supervision": False,
        "sealed_outcomes_loaded": False,
        "champion_eligible": False,
        "release_state": "internal_source_terms_apply_no_sequence_or_raw_redistribution_without_rights_review",
    }
    (arguments.output / "reconstruction_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--constructs", type=Path, required=True)
    parser.add_argument("--reconstruction-stats", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--replicate-outcomes", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    validate(parser.parse_args())


if __name__ == "__main__":
    main()
