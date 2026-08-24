#!/usr/bin/env python3
"""Freeze whole-chromosome folds and outcome-independent cCRE windows."""

from __future__ import annotations

import argparse
from bisect import bisect_left
import csv
import gzip
from hashlib import sha256
import heapq
import json
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence


CCRE_CLASSES = (
    "PLS",
    "pELS",
    "dELS",
    "CA-H3K4me3",
    "CA-CTCF",
    "CA-TF",
    "CA",
    "TF",
)
SEED = 20260824
FOLDS = 5
INPUT_LENGTH = 2114
OUTPUT_LENGTH = 1000
CCRES_PER_CLASS_PER_FOLD = 2000
CANDIDATE_MULTIPLIER = 5


class SequenceSplitError(ValueError):
    """Raised when a sequence split or evaluation window violates its contract."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_chrom_sizes(path: Path) -> list[tuple[str, int]]:
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            contig, size_text = line.rstrip("\n").split("\t")
            rows.append((contig, int(size_text)))
    expected = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])
    if tuple(contig for contig, _ in rows) != expected:
        raise SequenceSplitError("sequence-model chromosome order differs")
    if any(size <= INPUT_LENGTH for _, size in rows):
        raise SequenceSplitError("sequence-model chromosome is too short")
    return rows


def assign_balanced_folds(chrom_sizes: Sequence[tuple[str, int]]) -> dict[str, int]:
    bins: list[tuple[list[str], int]] = [([], 0) for _ in range(FOLDS)]
    for contig, size in sorted(chrom_sizes, key=lambda row: (-row[1], row[0])):
        fold = min(range(FOLDS), key=lambda value: (bins[value][1], value))
        bins[fold][0].append(contig)
        bins[fold] = (bins[fold][0], bins[fold][1] + size)
    assignment = {
        contig: fold for fold, (contigs, _size) in enumerate(bins) for contig in contigs
    }
    if set(assignment) != {contig for contig, _ in chrom_sizes}:
        raise SequenceSplitError("genomic fold assignment is incomplete")
    return assignment


def digest_integer(*values: object, seed: int = SEED) -> int:
    payload = "\0".join((str(seed), *(str(value) for value in values)))
    return int.from_bytes(sha256(payload.encode("utf-8")).digest(), "big")


def _candidate_from_fields(
    fields: Sequence[str],
    *,
    fold_by_contig: Mapping[str, int],
    chrom_size_by_contig: Mapping[str, int],
    seed: int,
) -> tuple[int, str, str, int, int, str, int, int, int, int, int] | None:
    if len(fields) < 10:
        raise SequenceSplitError("ENCODE cCRE row has fewer than ten fields")
    contig, start_text, end_text, identifier = fields[:4]
    ccre_class = fields[9]
    if contig not in fold_by_contig or ccre_class not in CCRE_CLASSES:
        return None
    start, end = int(start_text), int(end_text)
    center = (start + end) // 2
    output_start = center - OUTPUT_LENGTH // 2
    output_end = output_start + OUTPUT_LENGTH
    input_start = center - INPUT_LENGTH // 2
    input_end = input_start + INPUT_LENGTH
    if input_start < 0 or input_end > chrom_size_by_contig[contig]:
        return None
    digest = digest_integer(ccre_class, identifier, seed=seed)
    return (
        digest,
        identifier,
        contig,
        start,
        end,
        ccre_class,
        fold_by_contig[contig],
        output_start,
        output_end,
        input_start,
        input_end,
    )


def candidate_heaps(
    *,
    ccre_bed: Path,
    fold_by_contig: Mapping[str, int],
    chrom_size_by_contig: Mapping[str, int],
    per_class_per_fold: int,
    candidate_multiplier: int,
    seed: int,
) -> tuple[
    dict[tuple[int, str], list[tuple[int, str, tuple[object, ...]]]],
    dict[tuple[int, str], int],
]:
    capacity = per_class_per_fold * candidate_multiplier
    if capacity < per_class_per_fold:
        raise SequenceSplitError("candidate capacity is smaller than selection quota")
    heaps: dict[tuple[int, str], list[tuple[int, str, tuple[object, ...]]]] = {
        (fold, ccre_class): []
        for fold in range(FOLDS)
        for ccre_class in CCRE_CLASSES
    }
    eligible = {key: 0 for key in heaps}
    with gzip.open(ccre_bed, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            candidate = _candidate_from_fields(
                fields,
                fold_by_contig=fold_by_contig,
                chrom_size_by_contig=chrom_size_by_contig,
                seed=seed,
            )
            if candidate is None:
                continue
            digest, identifier, *_rest = candidate
            key = (int(candidate[6]), str(candidate[5]))
            eligible[key] += 1
            item = (-digest, str(identifier), candidate)
            heap = heaps[key]
            if len(heap) < capacity:
                heapq.heappush(heap, item)
            elif digest < -heap[0][0]:
                heapq.heapreplace(heap, item)
    return heaps, eligible


def filter_candidate_heaps_by_sequence(
    heaps: Mapping[
        tuple[int, str], Sequence[tuple[int, str, tuple[object, ...]]]
    ],
    *,
    fetch_sequence: Callable[[str, int, int], str],
) -> tuple[
    dict[tuple[int, str], list[tuple[int, str, tuple[object, ...]]]],
    dict[tuple[int, str], int],
]:
    """Exclude ambiguous input windows before quota and overlap selection."""

    filtered = {key: [] for key in heaps}
    excluded = {key: 0 for key in heaps}
    ordered = sorted(
        (
            (str(item[2][2]), int(item[2][9]), key, item)
            for key, values in heaps.items()
            for item in values
        ),
        key=lambda row: (row[0], row[1], row[2]),
    )
    for contig, start, key, item in ordered:
        end = int(item[2][10])
        sequence = fetch_sequence(contig, start, end).upper()
        if len(sequence) != INPUT_LENGTH or set(sequence).difference("ACGT"):
            excluded[key] += 1
            continue
        filtered[key].append(item)
    return filtered, excluded


def _add_if_nonoverlapping(
    occupied: dict[str, list[tuple[int, int]]], candidate: tuple[object, ...]
) -> bool:
    contig = str(candidate[2])
    start, end = int(candidate[9]), int(candidate[10])
    intervals = occupied.setdefault(contig, [])
    index = bisect_left(intervals, (start, end))
    if index > 0 and intervals[index - 1][1] > start:
        return False
    if index < len(intervals) and end > intervals[index][0]:
        return False
    intervals.insert(index, (start, end))
    return True


def select_windows(
    heaps: Mapping[tuple[int, str], Sequence[tuple[int, str, tuple[object, ...]]]],
    *,
    per_class_per_fold: int,
) -> list[tuple[object, ...]]:
    selected = []
    for fold in range(FOLDS):
        candidates = {
            ccre_class: sorted(
                (item[2] for item in heaps[(fold, ccre_class)]),
                key=lambda row: (int(row[0]), str(row[1])),
            )
            for ccre_class in CCRE_CLASSES
        }
        indices = {ccre_class: 0 for ccre_class in CCRE_CLASSES}
        counts = {ccre_class: 0 for ccre_class in CCRE_CLASSES}
        occupied: dict[str, list[tuple[int, int]]] = {}
        while any(counts[value] < per_class_per_fold for value in CCRE_CLASSES):
            progress = False
            for ccre_class in CCRE_CLASSES:
                if counts[ccre_class] >= per_class_per_fold:
                    continue
                values = candidates[ccre_class]
                while indices[ccre_class] < len(values):
                    candidate = values[indices[ccre_class]]
                    indices[ccre_class] += 1
                    if _add_if_nonoverlapping(occupied, candidate):
                        selected.append(candidate)
                        counts[ccre_class] += 1
                        progress = True
                        break
            if not progress:
                raise SequenceSplitError(
                    f"fold {fold} cannot satisfy nonoverlapping cCRE quotas: {counts}"
                )
    selected.sort(key=lambda row: (int(row[6]), str(row[2]), int(row[7]), str(row[1])))
    return selected


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_contract(
    *,
    ccre_acquisition: Path,
    ccre_artifacts_sha256: str,
    ccre_bed_sha256: str,
    reference_contract: Path,
    reference_artifacts_sha256: str,
    chrom_sizes_sha256: str,
    sequence_fasta_contract: Path,
    sequence_fasta_artifacts_sha256: str,
    fasta_sha256: str,
    output: Path,
    per_class_per_fold: int = CCRES_PER_CLASS_PER_FOLD,
    candidate_multiplier: int = CANDIDATE_MULTIPLIER,
    seed: int = SEED,
) -> None:
    ccre_acquisition = ccre_acquisition.resolve(strict=True)
    reference_contract = reference_contract.resolve(strict=True)
    sequence_fasta_contract = sequence_fasta_contract.resolve(strict=True)
    if output.exists():
        raise SequenceSplitError("refusing to overwrite sequence split contract")
    if sha256_file(ccre_acquisition / "ARTIFACTS.json") != ccre_artifacts_sha256:
        raise SequenceSplitError("cCRE acquisition ARTIFACTS SHA-256 differs")
    if sha256_file(reference_contract / "ARTIFACTS.json") != reference_artifacts_sha256:
        raise SequenceSplitError("reference contract ARTIFACTS SHA-256 differs")
    if (
        sha256_file(sequence_fasta_contract / "ARTIFACTS.json")
        != sequence_fasta_artifacts_sha256
    ):
        raise SequenceSplitError("sequence FASTA ARTIFACTS SHA-256 differs")
    ccre_bed = ccre_acquisition / "ENCFF420VPZ.bed.gz"
    chrom_sizes_path = reference_contract / "GRCh38.p14.sequence_model.chrom.sizes"
    fasta = sequence_fasta_contract / "GRCh38.p14.sequence_model.fa"
    if sha256_file(ccre_bed) != ccre_bed_sha256:
        raise SequenceSplitError("cCRE BED SHA-256 differs")
    if sha256_file(chrom_sizes_path) != chrom_sizes_sha256:
        raise SequenceSplitError("sequence chromosome sizes SHA-256 differs")
    if sha256_file(fasta) != fasta_sha256:
        raise SequenceSplitError("sequence FASTA SHA-256 differs")
    chrom_sizes = read_chrom_sizes(chrom_sizes_path)
    size_by_contig = dict(chrom_sizes)
    fold_by_contig = assign_balanced_folds(chrom_sizes)
    heaps, eligible = candidate_heaps(
        ccre_bed=ccre_bed,
        fold_by_contig=fold_by_contig,
        chrom_size_by_contig=size_by_contig,
        per_class_per_fold=per_class_per_fold,
        candidate_multiplier=candidate_multiplier,
        seed=seed,
    )
    if any(value < per_class_per_fold for value in eligible.values()):
        raise SequenceSplitError("a genomic-fold/cCRE-class stratum is under quota")
    import pyfaidx

    reference = pyfaidx.Fasta(
        str(fasta),
        as_raw=True,
        sequence_always_upper=True,
        rebuild=False,
        read_ahead=1_000_000,
    )
    try:
        heaps, ambiguous_excluded = filter_candidate_heaps_by_sequence(
            heaps,
            fetch_sequence=lambda contig, start, end: str(reference[contig][start:end]),
        )
    finally:
        reference.close()
    if any(len(values) < per_class_per_fold for values in heaps.values()):
        raise SequenceSplitError(
            "an unambiguous genomic-fold/cCRE-class candidate pool is under quota"
        )
    selected = select_windows(heaps, per_class_per_fold=per_class_per_fold)

    output.mkdir(mode=0o750)
    fold_rows = []
    for fold in range(FOLDS):
        contigs = [contig for contig, _ in chrom_sizes if fold_by_contig[contig] == fold]
        fold_rows.append(
            {
                "genomic_fold": fold,
                "contigs": ",".join(contigs),
                "total_bp": sum(size_by_contig[contig] for contig in contigs),
            }
        )
    write_tsv(output / "genomic_folds.tsv", ("genomic_fold", "contigs", "total_bp"), fold_rows)

    crossed_rows = []
    for donor_test in range(FOLDS):
        donor_valid = (donor_test + 1) % FOLDS
        donor_train = [value for value in range(FOLDS) if value not in {donor_test, donor_valid}]
        for genomic_test in range(FOLDS):
            genomic_valid = (genomic_test + 1) % FOLDS
            genomic_train = [value for value in range(FOLDS) if value not in {genomic_test, genomic_valid}]
            crossed_rows.append(
                {
                    "split_id": f"donor{donor_test}_genomic{genomic_test}",
                    "donor_train_folds": ",".join(map(str, donor_train)),
                    "donor_valid_fold": donor_valid,
                    "donor_test_fold": donor_test,
                    "genomic_train_folds": ",".join(map(str, genomic_train)),
                    "genomic_valid_fold": genomic_valid,
                    "genomic_test_fold": genomic_test,
                }
            )
    write_tsv(
        output / "crossed_outer_splits.tsv",
        (
            "split_id",
            "donor_train_folds",
            "donor_valid_fold",
            "donor_test_fold",
            "genomic_train_folds",
            "genomic_valid_fold",
            "genomic_test_fold",
        ),
        crossed_rows,
    )

    window_rows = []
    selected_counts = {(fold, value): 0 for fold in range(FOLDS) for value in CCRE_CLASSES}
    for index, candidate in enumerate(selected, start=1):
        (
            digest,
            identifier,
            contig,
            ccre_start,
            ccre_end,
            ccre_class,
            fold,
            output_start,
            output_end,
            input_start,
            input_end,
        ) = candidate
        selected_counts[(int(fold), str(ccre_class))] += 1
        window_rows.append(
            {
                "contig": contig,
                "output_start": output_start,
                "output_end": output_end,
                "window_id": f"ccre_{index:06d}",
                "genomic_fold": fold,
                "window_class": "encode_ccre",
                "ccre_class": ccre_class,
                "ccre_id": identifier,
                "ccre_start": ccre_start,
                "ccre_end": ccre_end,
                "input_start": input_start,
                "input_end": input_end,
                "selection_hash": f"{int(digest):064x}",
            }
        )
    write_tsv(
        output / "ccre_evaluation_windows.tsv",
        (
            "contig",
            "output_start",
            "output_end",
            "window_id",
            "genomic_fold",
            "window_class",
            "ccre_class",
            "ccre_id",
            "ccre_start",
            "ccre_end",
            "input_start",
            "input_end",
            "selection_hash",
        ),
        window_rows,
    )
    census_rows = []
    for fold in range(FOLDS):
        for ccre_class in CCRE_CLASSES:
            census_rows.append(
                {
                    "genomic_fold": fold,
                    "ccre_class": ccre_class,
                    "eligible_ccres": eligible[(fold, ccre_class)],
                    "candidate_pool_ambiguous_excluded": ambiguous_excluded[
                        (fold, ccre_class)
                    ],
                    "selected_windows": selected_counts[(fold, ccre_class)],
                    "registry_weight": eligible[(fold, ccre_class)] / per_class_per_fold,
                }
            )
    write_tsv(
        output / "ccre_selection_census.tsv",
        (
            "genomic_fold",
            "ccre_class",
            "eligible_ccres",
            "candidate_pool_ambiguous_excluded",
            "selected_windows",
            "registry_weight",
        ),
        census_rows,
    )
    contract = {
        "schema_version": "masld-bench-sequence-split-contract-v2",
        "dataset_id": "gse296875",
        "biological_outer_unit": "donor",
        "genomic_outer_unit": "whole_chromosome_group",
        "donor_folds": FOLDS,
        "genomic_folds": FOLDS,
        "crossed_outer_splits": FOLDS * FOLDS,
        "genomic_fold_algorithm": "largest_chromosome_first_greedy_balance",
        "boundary_buffer_required": False,
        "boundary_buffer_reason": "whole_contigs_do_not_share_linear_boundaries",
        "input_length": INPUT_LENGTH,
        "output_length": OUTPUT_LENGTH,
        "evaluation_source": "ENCODE_v4_ENCFF420VPZ",
        "evaluation_source_is_outcome_independent": True,
        "held_atac_used_for_window_selection": False,
        "sequence_ambiguity_policy": "exclude_non_ACGT_input_windows_before_selection",
        "selection_seed": seed,
        "ccre_classes": list(CCRE_CLASSES),
        "ccres_per_class_per_fold": per_class_per_fold,
        "ccre_windows": len(selected),
        "ccre_input_windows_nonoverlapping_within_chromosome": True,
        "primary_endpoint_chromosomes": "autosomes",
        "secondary_endpoint_chromosomes": "chrX_chrY",
        "ccre_artifacts_sha256": ccre_artifacts_sha256,
        "reference_artifacts_sha256": reference_artifacts_sha256,
        "sequence_fasta_artifacts_sha256": sequence_fasta_artifacts_sha256,
        "sequence_fasta_sha256": fasta_sha256,
    }
    (output / "contract.json").write_text(
        json.dumps(contract, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ccre-acquisition", type=Path, required=True)
    parser.add_argument("--ccre-artifacts-sha256", required=True)
    parser.add_argument("--ccre-bed-sha256", required=True)
    parser.add_argument("--reference-contract", type=Path, required=True)
    parser.add_argument("--reference-artifacts-sha256", required=True)
    parser.add_argument("--chrom-sizes-sha256", required=True)
    parser.add_argument("--sequence-fasta-contract", type=Path, required=True)
    parser.add_argument("--sequence-fasta-artifacts-sha256", required=True)
    parser.add_argument("--fasta-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    build_contract(
        ccre_acquisition=arguments.ccre_acquisition,
        ccre_artifacts_sha256=arguments.ccre_artifacts_sha256,
        ccre_bed_sha256=arguments.ccre_bed_sha256,
        reference_contract=arguments.reference_contract,
        reference_artifacts_sha256=arguments.reference_artifacts_sha256,
        chrom_sizes_sha256=arguments.chrom_sizes_sha256,
        sequence_fasta_contract=arguments.sequence_fasta_contract,
        sequence_fasta_artifacts_sha256=arguments.sequence_fasta_artifacts_sha256,
        fasta_sha256=arguments.fasta_sha256,
        output=arguments.output,
    )
    print(arguments.output.resolve(strict=True))


if __name__ == "__main__":
    main()
