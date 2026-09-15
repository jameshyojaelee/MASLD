#!/usr/bin/env python3
"""Build outcome-blind, 300-bp LS-GKM matched training inputs."""

from __future__ import annotations

import argparse
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from scripts.lsgkm_safe_adapter import reverse_complement


SCHEMA = "masld-bench-lsgkm-gse281364-input-activation-v1"
PRIMARY_CONTIGS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
ROW_SEEDS = (1103, 2909, 4721, 6673, 8111)
LEGACY_SEEDS = (11, 29, 47, 71, 101)
TSS_BIN_EDGES = (1_000, 10_000, 100_000, 1_000_000)
MATCH_FIELDS = (
    "seed",
    "pair_id",
    "split_id",
    "genomic_test_fold",
    "positive_id",
    "negative_id",
    "positive_contig",
    "negative_contig",
    "positive_start0",
    "positive_end0",
    "negative_start0",
    "negative_end0",
    "positive_atac_decile",
    "negative_atac_decile",
    "positive_tss_distance_bin",
    "negative_tss_distance_bin",
    "positive_gc_fraction",
    "negative_gc_fraction",
    "positive_mappability",
    "negative_mappability",
    "positive_log1p_atac",
    "negative_log1p_atac",
    "positive_tss_distance",
    "negative_tss_distance",
    "positive_repeat_fraction",
    "negative_repeat_fraction",
    "standardized_covariate_distance",
    "positive_sequence_sha256",
    "negative_sequence_sha256",
)
SCORE_FIELDS = (
    "seed",
    "study_id",
    "assay_context_id",
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
    "outer_fold",
    "contig",
    "variant_pos0",
    "ref",
    "alt",
    "scoring_split_id",
    "genomic_test_fold",
    "model_state_id",
    "context_specific_prediction",
)
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class MatchedInputError(RuntimeError):
    """Raised when an authority, split, interval, or match differs."""


@dataclass(frozen=True)
class Window:
    identifier: str
    contig: str
    start: int
    end: int
    sequence: str
    gc: float
    mappability: float
    log1p_atac: float
    tss_distance: int
    tss_bin: int
    repeat_fraction: float
    atac_decile: int = -1

    @property
    def canonical_sequence(self) -> str:
        reverse = reverse_complement(self.sequence)
        return min(self.sequence, reverse)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_frozen_tree(root: Path) -> dict[str, Any]:
    """Verify the read-only tree protocol without importing the Python 3.11 control plane."""

    if not root.is_dir() or root.is_symlink():
        raise MatchedInputError(f"invalid frozen-tree root: {root}")
    manifest_path = root / "ARTIFACTS.json"
    complete_path = root / "COMPLETE"
    if (
        not manifest_path.is_file()
        or manifest_path.is_symlink()
        or not complete_path.is_file()
        or complete_path.is_symlink()
    ):
        raise MatchedInputError(f"incomplete frozen tree: {root}")
    manifest = load_json(manifest_path)
    complete = load_json(complete_path)
    raw_artifacts = manifest.get("artifacts")
    if (
        set(manifest) != {"schema_version", "metadata", "artifacts"}
        or manifest.get("schema_version") != "masld-bench-artifacts-v1"
        or not isinstance(manifest.get("metadata"), dict)
        or not isinstance(raw_artifacts, list)
        or complete
        != {
            "artifact_count": len(raw_artifacts),
            "manifest_sha256": file_sha256(manifest_path),
            "schema_version": "masld-bench-complete-v1",
        }
    ):
        raise MatchedInputError(f"frozen-tree control file differs: {root}")
    expected: dict[str, tuple[str, int]] = {}
    for item in raw_artifacts:
        if not isinstance(item, dict) or set(item) != {"path", "sha256", "size_bytes"}:
            raise MatchedInputError("artifact manifest member differs")
        relative = Path(str(item["path"]))
        relative_text = relative.as_posix()
        checksum = str(item["sha256"])
        size = item["size_bytes"]
        if (
            relative.is_absolute()
            or not relative.parts
            or ".." in relative.parts
            or relative_text in expected
            or not SHA256.fullmatch(checksum)
            or isinstance(size, bool)
            or not isinstance(size, int)
            or size < 0
        ):
            raise MatchedInputError("unsafe artifact manifest member")
        expected[relative_text] = (checksum, size)
    observed: set[str] = set()
    for path in root.rglob("*"):
        if path.is_symlink():
            raise MatchedInputError(f"symlink in frozen tree: {path}")
        relative = path.relative_to(root).as_posix()
        if path.is_file() and relative not in {"ARTIFACTS.json", "COMPLETE"}:
            observed.add(relative)
    if observed != set(expected):
        raise MatchedInputError("frozen-tree file inventory differs")
    for relative, (expected_hash, expected_size) in expected.items():
        path = root / relative
        if path.stat().st_size != expected_size or file_sha256(path) != expected_hash:
            raise MatchedInputError(f"frozen artifact differs: {relative}")
    return manifest


def hash_id(*parts: object) -> str:
    return hashlib.sha256("\0".join(map(str, parts)).encode("utf-8")).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise MatchedInputError(f"invalid JSON: {path}") from error
    if not isinstance(value, dict):
        raise MatchedInputError(f"JSON authority is not an object: {path}")
    return value


def load_config(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise MatchedInputError("invalid LS-GKM activation config") from error
    if not isinstance(value, dict) or value.get("schema_version") != SCHEMA:
        raise MatchedInputError("LS-GKM activation schema differs")
    return value


def read_tsv(path: Path, fields: Sequence[str] | None = None) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if fields is not None and tuple(reader.fieldnames or ()) != tuple(fields):
            raise MatchedInputError(f"TSV header differs: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise MatchedInputError(f"TSV is empty: {path}")
    return rows


def write_tsv_gzip(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    if path.exists():
        raise MatchedInputError(f"refusing to overwrite: {path}")
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for row in rows:
                    writer.writerow({field: row[field] for field in fields})


def write_fasta_gzip(path: Path, records: Iterable[tuple[str, str]]) -> int:
    if path.exists():
        raise MatchedInputError(f"refusing to overwrite: {path}")
    count = 0
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
            for identifier, sequence in records:
                compressed.write(f">{identifier}\n{sequence}\n".encode("ascii"))
                count += 1
    return count


def validate_tree(root: Path, binding: Mapping[str, object], label: str) -> Path:
    relative = Path(str(binding.get("path", "")))
    if relative.is_absolute() or ".." in relative.parts:
        raise MatchedInputError(f"unsafe tree authority: {label}")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    expected = str(binding.get("artifacts_sha256", ""))
    if not SHA256.fullmatch(expected) or file_sha256(path / "ARTIFACTS.json") != expected:
        raise MatchedInputError(f"tree authority differs: {label}")
    verify_frozen_tree(path)
    return path


def validate_file(root: Path, binding: Mapping[str, object], label: str) -> Path:
    configured = Path(str(binding.get("path", "")))
    path = configured if configured.is_absolute() else root / configured
    if ".." in configured.parts or path.is_symlink() or not path.is_file():
        raise MatchedInputError(f"unsafe file authority: {label}")
    expected = str(binding.get("sha256", ""))
    if not SHA256.fullmatch(expected) or file_sha256(path) != expected:
        raise MatchedInputError(f"file authority differs: {label}")
    return path.resolve(strict=True)


def merged_intervals(rows: Iterable[tuple[str, int, int]]) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for contig, start, end in rows:
        if contig in PRIMARY_CONTIGS and start >= 0 and end > start:
            grouped[contig].append((start, end))
    output: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for contig in PRIMARY_CONTIGS:
        merged: list[list[int]] = []
        for start, end in sorted(grouped.get(contig, ())):
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            elif end > merged[-1][1]:
                merged[-1][1] = end
        output[contig] = (
            np.asarray([item[0] for item in merged], dtype=np.int64),
            np.asarray([item[1] for item in merged], dtype=np.int64),
        )
    return output


def interval_overlap_bases(index: Mapping[str, tuple[np.ndarray, np.ndarray]], contig: str, start: int, end: int) -> int:
    starts, ends = index.get(contig, (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)))
    position = int(np.searchsorted(ends, start, side="right"))
    total = 0
    while position < len(starts) and int(starts[position]) < end:
        total += max(0, min(end, int(ends[position])) - max(start, int(starts[position])))
        position += 1
    return total


def overlap_exists(index: Mapping[str, tuple[np.ndarray, np.ndarray]], contig: str, start: int, end: int) -> bool:
    return interval_overlap_bases(index, contig, start, end) > 0


def load_repeat_index(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    def rows() -> Iterable[tuple[str, int, int]]:
        with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
            for line in handle:
                fields = line.split()
                if not fields or not fields[0].isdigit():
                    continue
                if len(fields) not in {14, 15}:
                    raise MatchedInputError("RepeatMasker row width differs")
                yield fields[4], int(fields[5]) - 1, int(fields[6])

    return merged_intervals(rows())


def load_bed_index(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    def rows() -> Iterable[tuple[str, int, int]]:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                yield fields[0], int(fields[1]), int(fields[2])

    return merged_intervals(rows())


def load_tss(path: Path) -> dict[str, np.ndarray]:
    positions: dict[str, set[int]] = defaultdict(set)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene" or fields[0] not in PRIMARY_CONTIGS:
                continue
            start1, end1 = int(fields[3]), int(fields[4])
            tss0 = start1 - 1 if fields[6] == "+" else end1 - 1
            positions[fields[0]].add(tss0)
    if any(not positions.get(contig) for contig in PRIMARY_CONTIGS):
        raise MatchedInputError("GENCODE v49 primary-contig TSS roster is incomplete")
    return {contig: np.asarray(sorted(positions[contig]), dtype=np.int64) for contig in PRIMARY_CONTIGS}


def nearest_tss_distance(index: Mapping[str, np.ndarray], contig: str, center: int) -> int:
    values = index[contig]
    position = int(np.searchsorted(values, center))
    candidates = []
    if position < len(values):
        candidates.append(abs(int(values[position]) - center))
    if position:
        candidates.append(abs(int(values[position - 1]) - center))
    if not candidates:
        raise MatchedInputError(f"no TSS on {contig}")
    return min(candidates)


def tss_bin(distance: int) -> int:
    return bisect_left(TSS_BIN_EDGES, distance + 1)


def mean_with_missing_zero(total: float | None, width: int) -> float:
    if width <= 0:
        raise MatchedInputError("interval width differs")
    numeric = 0.0 if total is None else float(total)
    value = numeric / width
    if not math.isfinite(value) or value < 0.0 or value > 1.0:
        raise MatchedInputError("bounded interval mean differs")
    return value


def read_peak_windows(path: Path) -> list[tuple[str, str, int, int]]:
    rows: list[tuple[str, str, int, int]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 10:
                raise MatchedInputError("positive narrowPeak width differs")
            center = int(fields[1]) + int(fields[9])
            start, end = center - 150, center + 150
            if start < 0 or end - start != 300:
                raise MatchedInputError("positive 300-bp geometry differs")
            rows.append((f"peak_{line_number:06d}", fields[0], start, end))
    return rows


def read_candidate_windows(path: Path, allowed_contigs: set[str]) -> list[tuple[str, str, int, int]]:
    rows: list[tuple[str, str, int, int]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 10 or fields[0] not in allowed_contigs:
                continue
            center = int(fields[1]) + int(fields[9])
            start, end = center - 150, center + 150
            if start < 0 or end - start != 300:
                continue
            rows.append((f"candidate_{line_number:06d}", fields[0], start, end))
    if not rows:
        raise MatchedInputError("training candidate universe is empty")
    return rows


def decile_boundaries(values: Sequence[float]) -> tuple[float, ...]:
    if not values or any(not math.isfinite(value) for value in values):
        raise MatchedInputError("ATAC decile input differs")
    return tuple(float(value) for value in np.quantile(np.asarray(values), np.arange(1, 10) / 10.0))


def assign_decile(value: float, boundaries: Sequence[float]) -> int:
    return int(np.searchsorted(np.asarray(boundaries), value, side="right"))


def annotate_window(
    raw: tuple[str, str, int, int],
    *,
    fasta: Any,
    umap: Any,
    atac: Any,
    repeats: Mapping[str, tuple[np.ndarray, np.ndarray]],
    tsses: Mapping[str, np.ndarray],
) -> Window | None:
    identifier, contig, start, end = raw
    sequence = str(fasta[contig][start:end]).upper()
    if len(sequence) != 300 or any(base not in "ACGT" for base in sequence):
        return None
    umap_sum = umap.stats(contig, start, end, type="sum", exact=True)[0]
    mappability = mean_with_missing_zero(umap_sum, end - start)
    if mappability < 0.80:
        return None
    atac_value = atac.stats(contig, start, end, type="sum", exact=True)[0]
    atac_sum = 0.0 if atac_value is None else float(atac_value)
    if not math.isfinite(atac_sum) or atac_sum < 0:
        raise MatchedInputError("ATAC window sum differs")
    center = start + 150
    distance = nearest_tss_distance(tsses, contig, center)
    return Window(
        identifier=identifier,
        contig=contig,
        start=start,
        end=end,
        sequence=sequence,
        gc=(sequence.count("G") + sequence.count("C")) / 300.0,
        mappability=mappability,
        log1p_atac=math.log1p(atac_sum),
        tss_distance=distance,
        tss_bin=tss_bin(distance),
        repeat_fraction=interval_overlap_bases(repeats, contig, start, end) / 300.0,
    )


def with_decile(window: Window, boundaries: Sequence[float]) -> Window:
    return Window(
        identifier=window.identifier,
        contig=window.contig,
        start=window.start,
        end=window.end,
        sequence=window.sequence,
        gc=window.gc,
        mappability=window.mappability,
        log1p_atac=window.log1p_atac,
        tss_distance=window.tss_distance,
        tss_bin=window.tss_bin,
        repeat_fraction=window.repeat_fraction,
        atac_decile=assign_decile(window.log1p_atac, boundaries),
    )


def deduplicate_windows(windows: Iterable[Window]) -> list[Window]:
    retained: dict[str, Window] = {}
    for window in sorted(windows, key=lambda value: (value.contig, value.start, value.identifier)):
        retained.setdefault(window.canonical_sequence, window)
    return sorted(retained.values(), key=lambda value: (value.contig, value.start, value.identifier))


def exclude_cross_class_sequence_duplicates(
    positives: Sequence[Window], negatives: Sequence[Window]
) -> list[Window]:
    positive_sequences = {window.canonical_sequence for window in positives}
    return [
        window for window in negatives if window.canonical_sequence not in positive_sequences
    ]


def match_windows(positives: Sequence[Window], negatives: Sequence[Window], seed: int) -> list[tuple[Window, Window, float]]:
    from scipy.spatial import cKDTree

    pools: dict[tuple[str, int, int], list[Window]] = defaultdict(list)
    for candidate in negatives:
        pools[(candidate.contig, candidate.atac_decile, candidate.tss_bin)].append(candidate)
    trees: dict[tuple[str, int, int], Any] = {}
    for key, values in pools.items():
        values.sort(key=lambda value: (value.gc, value.mappability, value.repeat_fraction, value.start))
        trees[key] = cKDTree(
            np.asarray(
                [
                    (value.gc / 0.02, value.mappability / 0.05, value.repeat_fraction / 0.05)
                    for value in values
                ],
                dtype=np.float64,
            )
        )
    ordered = sorted(positives, key=lambda value: hash_id("positive-order", seed, value.identifier, value.contig, value.start))
    used: set[str] = set()
    matches: list[tuple[Window, Window, float]] = []
    for positive in ordered:
        key = (positive.contig, positive.atac_decile, positive.tss_bin)
        pool = pools.get(key, ())
        tree = trees.get(key)
        if tree is None:
            continue
        point = np.asarray(
            (positive.gc / 0.02, positive.mappability / 0.05, positive.repeat_fraction / 0.05),
            dtype=np.float64,
        )
        eligible: list[tuple[float, str, Window]] = []
        for index in tree.query_ball_point(point, r=3.0 + 1e-12, p=1.0):
            candidate = pool[int(index)]
            if candidate.identifier in used:
                continue
            gc_delta = abs(candidate.gc - positive.gc)
            map_delta = abs(candidate.mappability - positive.mappability)
            repeat_delta = abs(candidate.repeat_fraction - positive.repeat_fraction)
            if gc_delta > 0.02 + 1e-12 or map_delta > 0.05 + 1e-12 or repeat_delta > 0.05 + 1e-12:
                continue
            distance = gc_delta / 0.02 + map_delta / 0.05 + repeat_delta / 0.05
            eligible.append((distance, hash_id("candidate-tie", seed, positive.identifier, candidate.identifier), candidate))
        if eligible:
            distance, _tie, candidate = min(eligible, key=lambda item: (item[0], item[1]))
            used.add(candidate.identifier)
            matches.append((positive, candidate, distance))
    return sorted(matches, key=lambda item: (item[0].contig, item[0].start, item[0].identifier))


def outer_fold_join_matches(row_universe_fold: str, borzoi_fold: str) -> bool:
    match = re.fullmatch(r"fold-([0-4])", row_universe_fold)
    return match is not None and borzoi_fold == match.group(1)


def build_scoring_map(row_universe: Path, borzoi_manifest: Path, fold_by_contig: Mapping[str, int], split_by_test_fold: Mapping[int, str]) -> list[dict[str, object]]:
    rows = read_tsv(row_universe)
    if {row["assay_context_id"] for row in rows} != set(CONTEXTS):
        raise MatchedInputError("GSE281364 context roster differs")
    if {int(row["seed"]) for row in rows} != set(ROW_SEEDS):
        raise MatchedInputError("GSE281364 seed roster differs")
    if {row["stratum"] for row in rows} != {"all"} or len(rows) != 10_330:
        raise MatchedInputError("GSE281364 row universe differs")
    manifest = {row["element_id"]: row for row in read_tsv(borzoi_manifest)}
    if len(manifest) != 1_033:
        raise MatchedInputError("Borzoi allele manifest denominator differs")
    output: list[dict[str, object]] = []
    for row in rows:
        source = manifest.get(row["element_id"])
        if source is None or source["source_locus_group_id"] != row["source_locus_group_id"] or source["borzoi_long_range_group_id"] != row["long_range_block_id"] or not outer_fold_join_matches(row["outer_fold"], source["outer_fold"]):
            raise MatchedInputError("row-universe allele join differs")
        contig = source["contig"]
        test_fold = fold_by_contig[contig]
        output.append({
            "seed": int(row["seed"]),
            "study_id": "gse281364",
            "assay_context_id": row["assay_context_id"],
            "element_id": row["element_id"],
            "source_locus_group_id": row["source_locus_group_id"],
            "long_range_block_id": row["long_range_block_id"],
            "outer_fold": row["outer_fold"],
            "contig": contig,
            "variant_pos0": int(source["variant_pos0"]),
            "ref": source["ref"],
            "alt": source["alt"],
            "scoring_split_id": split_by_test_fold[test_fold],
            "genomic_test_fold": test_fold,
            "model_state_id": "hepatocyte",
            "context_specific_prediction": "false",
        })
    return output


def match_row(seed: int, split_id: str, test_fold: int, positive: Window, negative: Window, distance: float) -> dict[str, object]:
    pair_id = hash_id("lsgkm-pair-v1", seed, split_id, positive.identifier, negative.identifier)
    return {
        "seed": seed,
        "pair_id": pair_id,
        "split_id": split_id,
        "genomic_test_fold": test_fold,
        "positive_id": positive.identifier,
        "negative_id": negative.identifier,
        "positive_contig": positive.contig,
        "negative_contig": negative.contig,
        "positive_start0": positive.start,
        "positive_end0": positive.end,
        "negative_start0": negative.start,
        "negative_end0": negative.end,
        "positive_atac_decile": positive.atac_decile,
        "negative_atac_decile": negative.atac_decile,
        "positive_tss_distance_bin": positive.tss_bin,
        "negative_tss_distance_bin": negative.tss_bin,
        "positive_gc_fraction": f"{positive.gc:.9f}",
        "negative_gc_fraction": f"{negative.gc:.9f}",
        "positive_mappability": f"{positive.mappability:.9f}",
        "negative_mappability": f"{negative.mappability:.9f}",
        "positive_log1p_atac": f"{positive.log1p_atac:.9f}",
        "negative_log1p_atac": f"{negative.log1p_atac:.9f}",
        "positive_tss_distance": positive.tss_distance,
        "negative_tss_distance": negative.tss_distance,
        "positive_repeat_fraction": f"{positive.repeat_fraction:.9f}",
        "negative_repeat_fraction": f"{negative.repeat_fraction:.9f}",
        "standardized_covariate_distance": f"{distance:.12f}",
        "positive_sequence_sha256": hashlib.sha256(positive.sequence.encode("ascii")).hexdigest(),
        "negative_sequence_sha256": hashlib.sha256(negative.sequence.encode("ascii")).hexdigest(),
    }


def validate_config_boundary(config: Mapping[str, Any]) -> None:
    if (
        config.get("status") != "prespecified_outcome_blind_input_materialization"
        or config.get("outcome_access_authorized") is not False
        or config.get("sealed_asset_access_authorized") is not False
        or config.get("production_training_authorized") is not False
        or config.get("production_prediction_authorized") is not False
    ):
        raise MatchedInputError("action firewall differs")
    design = config.get("design", {})
    if (
        design.get("active_state_roster") != ["hepatocyte"]
        or tuple(design.get("active_assay_contexts", ())) != CONTEXTS
        or design.get("context_specific_prediction") is not False
        or tuple(design.get("fixed_seeds", ())) != ROW_SEEDS
        or tuple(design.get("superseded_legacy_seeds", ())) != LEGACY_SEEDS
        or design.get("diagonal_split_count") != 5
        or design.get("shared_fit_count") != 25
        or design.get("readouts_per_shared_fit") != 2
        or design.get("independent_model_families") != 1
    ):
        raise MatchedInputError("prospective minimal-fit design differs")


def matching_gate_result(
    *,
    matched_pairs: int,
    eligible_positive_windows: int,
    minimum_matched_pairs: int,
    minimum_positive_coverage: float,
) -> dict[str, object]:
    coverage = (
        matched_pairs / eligible_positive_windows
        if eligible_positive_windows
        else 0.0
    )
    return {
        "matched_pairs": matched_pairs,
        "eligible_positive_windows": eligible_positive_windows,
        "positive_coverage": coverage,
        "minimum_matched_pairs": minimum_matched_pairs,
        "minimum_positive_coverage": minimum_positive_coverage,
        "pair_count_gate_passed": matched_pairs >= minimum_matched_pairs,
        "positive_coverage_gate_passed": coverage >= minimum_positive_coverage,
        "passed": (
            matched_pairs >= minimum_matched_pairs
            and coverage >= minimum_positive_coverage
        ),
    }


def build(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.project_root.resolve(strict=True)
    config = load_config(arguments.config)
    validate_config_boundary(config)
    if arguments.output.exists():
        raise MatchedInputError("matched-input output exists")
    arguments.output.mkdir(parents=True, mode=0o750)
    (arguments.output / "inputs").mkdir(mode=0o750)
    (arguments.output / "contract").mkdir(mode=0o750)

    trees = {name: validate_tree(root, binding, name) for name, binding in config["tree_authorities"].items()}
    files = {name: validate_file(root, binding, name) for name, binding in config["file_authorities"].items()}
    track_receipt = load_json(trees["matching_tracks"] / "tracks/receipt.json")
    if track_receipt.get("outcomes_read") is not False or track_receipt.get("sealed_assets_read") is not False:
        raise MatchedInputError("matching-track firewall differs")
    umap_path = trees["matching_tracks"] / "tracks/source/k100.Umap.MultiTrackMappability.bw"
    rmsk_path = trees["matching_tracks"] / "tracks/source/hg38.fa.out.gz"
    fasta_path = trees["sequence_fasta"] / "GRCh38.p14.sequence_model.fa"
    row_universe = trees["row_universe"] / "contract/row_universe.tsv"
    borzoi_manifest = trees["borzoi_fixture"] / "fixture/manifest.tsv"

    split_rows = read_tsv(trees["sequence_split"] / "crossed_outer_splits.tsv")
    genomic_rows = read_tsv(trees["sequence_split"] / "genomic_folds.tsv")
    fold_by_contig = {
        contig: int(row["genomic_fold"])
        for row in genomic_rows
        for contig in row["contigs"].split(",")
    }
    diagonal = {
        int(row["genomic_test_fold"]): row
        for row in split_rows
        if row["split_id"] == f"donor{row['donor_test_fold']}_genomic{row['genomic_test_fold']}"
        and row["donor_test_fold"] == row["genomic_test_fold"]
    }
    if len(diagonal) != 5 or set(fold_by_contig) != set(PRIMARY_CONTIGS):
        raise MatchedInputError("diagonal split or contig roster differs")
    split_by_test_fold = {fold: row["split_id"] for fold, row in diagonal.items()}
    scoring = build_scoring_map(row_universe, borzoi_manifest, fold_by_contig, split_by_test_fold)
    write_tsv_gzip(arguments.output / "contract/scoring_row_map.tsv.gz", SCORE_FIELDS, scoring)

    repeats = load_repeat_index(rmsk_path)
    tsses = load_tss(files["gencode_v49_gtf"])
    blacklist = load_bed_index(trees["encode_blacklist"] / "ENCFF356LFX.bed.gz")

    from pyfaidx import Fasta
    import pyBigWig

    fasta = Fasta(str(fasta_path), as_raw=True, sequence_always_upper=True, rebuild=False)
    umap = pyBigWig.open(str(umap_path))
    pseudobulk_manifest = read_tsv(trees["training_pseudobulk"] / "training_pseudobulk_manifest.tsv")
    pseudobulk_index = {
        (int(row["donor_test_fold"]), row["lineage_id"]): row for row in pseudobulk_manifest
    }
    source_trees = config["training_input_trees"]
    summaries: list[dict[str, Any]] = []
    try:
        for test_fold in range(5):
            split = diagonal[test_fold]
            split_id = split["split_id"]
            source = validate_tree(root, source_trees[split_id], f"training_inputs_{split_id}")
            source_meta = load_json(source / "ARTIFACTS.json")["metadata"]
            if source_meta.get("lineage_id") != "hepatocyte" or source_meta.get("split_id") != split_id or source_meta.get("held_atac_used_for_peak_discovery") is not False:
                raise MatchedInputError("training-input tree topology differs")
            allowed_contigs = {
                contig for contig, fold in fold_by_contig.items() if str(fold) in split["genomic_train_folds"].split(",")
            }
            fold_json = load_json(source / "regions/chrombpnet.fold.json")
            if set(fold_json["train"]) != allowed_contigs or set(fold_json["test"]) != {contig for contig, fold in fold_by_contig.items() if fold == test_fold}:
                raise MatchedInputError("training or held-contig roster differs")
            all_peaks = load_bed_index(source / "regions/macs2.training.narrowPeak")
            raw_positives = read_peak_windows(source / "regions/training.filtered.narrowPeak")
            raw_negatives = read_candidate_windows(source / "nonpeaks/chrombpnet_negatives.bed", allowed_contigs)
            pseudo = pseudobulk_index[(test_fold, "hepatocyte")]
            atac_path = trees["training_pseudobulk"] / pseudo["bigwig_path"]
            if file_sha256(atac_path) != pseudo["bigwig_sha256"]:
                raise MatchedInputError("training-only ATAC bigWig differs")
            atac = pyBigWig.open(str(atac_path))
            try:
                positive_windows: list[Window] = []
                for raw in raw_positives:
                    if raw[1] not in allowed_contigs or overlap_exists(blacklist, raw[1], raw[2], raw[3]):
                        continue
                    window = annotate_window(raw, fasta=fasta, umap=umap, atac=atac, repeats=repeats, tsses=tsses)
                    if window is not None:
                        positive_windows.append(window)
                negative_windows: list[Window] = []
                for raw in raw_negatives:
                    if overlap_exists(blacklist, raw[1], raw[2], raw[3]) or overlap_exists(all_peaks, raw[1], raw[2], raw[3]):
                        continue
                    window = annotate_window(raw, fasta=fasta, umap=umap, atac=atac, repeats=repeats, tsses=tsses)
                    if window is not None:
                        negative_windows.append(window)
            finally:
                atac.close()
            positive_windows = deduplicate_windows(positive_windows)
            negative_windows = deduplicate_windows(negative_windows)
            negative_windows = exclude_cross_class_sequence_duplicates(
                positive_windows, negative_windows
            )
            boundaries = decile_boundaries([window.log1p_atac for window in positive_windows + negative_windows])
            positive_windows = [with_decile(window, boundaries) for window in positive_windows]
            negative_windows = [with_decile(window, boundaries) for window in negative_windows]
            split_summary: dict[str, Any] = {
                "split_id": split_id,
                "genomic_test_fold": test_fold,
                "donor_test_fold": int(split["donor_test_fold"]),
                "donor_valid_fold": int(split["donor_valid_fold"]),
                "donor_train_folds": [int(value) for value in split["donor_train_folds"].split(",")],
                "genomic_train_folds": [int(value) for value in split["genomic_train_folds"].split(",")],
                "held_test_contigs": sorted(fold_json["test"]),
                "held_valid_contigs": sorted(fold_json["valid"]),
                "raw_positive_windows": len(raw_positives),
                "eligible_unique_positive_windows": len(positive_windows),
                "raw_training_candidate_windows": len(raw_negatives),
                "eligible_unique_training_candidate_windows": len(negative_windows),
                "atac_decile_boundaries": list(boundaries),
                "seeds": {},
            }
            for seed in ROW_SEEDS:
                matched = match_windows(positive_windows, negative_windows, seed)
                minimum = int(config["matching"]["minimum_matched_pairs_per_fit"])
                gate = matching_gate_result(
                    matched_pairs=len(matched),
                    eligible_positive_windows=len(positive_windows),
                    minimum_matched_pairs=minimum,
                    minimum_positive_coverage=float(
                        config["matching"]["minimum_positive_coverage"]
                    ),
                )
                coverage = float(gate["positive_coverage"])
                if gate["passed"] is not True:
                    failure = {
                        "schema_version": "masld-bench-lsgkm-matching-gate-failure-v1",
                        "status": "terminally_blocked_prespecified_matching_gate",
                        "split_id": split_id,
                        "seed": seed,
                        "genomic_test_fold": test_fold,
                        "raw_positive_windows": len(raw_positives),
                        "eligible_unique_positive_windows": len(positive_windows),
                        "raw_training_candidate_windows": len(raw_negatives),
                        "eligible_unique_training_candidate_windows": len(negative_windows),
                        "gate": gate,
                        "input_fastas_written": 0,
                        "production_fits_executed": 0,
                        "production_predictions_generated": 0,
                        "outcomes_read": False,
                        "reporter_counts_read": False,
                        "sealed_assets_read": False,
                    }
                    (arguments.output / "contract/matching_gate_failure.json").write_text(
                        json.dumps(failure, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8",
                    )
                    raise MatchedInputError(f"matched-input coverage failed: {split_id} seed {seed}")
                stem = arguments.output / "inputs" / f"{split_id}.seed{seed}"
                rows = [match_row(seed, split_id, test_fold, positive, negative, distance) for positive, negative, distance in matched]
                write_tsv_gzip(Path(f"{stem}.matches.tsv.gz"), MATCH_FIELDS, rows)
                positive_count = write_fasta_gzip(
                    Path(f"{stem}.positive.fa.gz"),
                    (
                        (rows[index]["pair_id"] + "|positive", matched[index][0].sequence)
                        for index in range(len(matched))
                    ),
                )
                negative_count = write_fasta_gzip(
                    Path(f"{stem}.negative.fa.gz"),
                    (
                        (rows[index]["pair_id"] + "|negative", matched[index][1].sequence)
                        for index in range(len(matched))
                    ),
                )
                if positive_count != len(matched) or negative_count != len(matched):
                    raise MatchedInputError("FASTA pair count differs")
                split_summary["seeds"][str(seed)] = {
                    "matched_pairs": len(matched),
                    "positive_coverage": coverage,
                    "match_manifest_sha256": file_sha256(Path(f"{stem}.matches.tsv.gz")),
                    "positive_fasta_sha256": file_sha256(Path(f"{stem}.positive.fa.gz")),
                    "negative_fasta_sha256": file_sha256(Path(f"{stem}.negative.fa.gz")),
                    "maximum_standardized_covariate_distance": max(distance for _positive, _negative, distance in matched),
                }
            summaries.append(split_summary)
    finally:
        umap.close()
        fasta.close()

    matched_counts = [seed["matched_pairs"] for summary in summaries for seed in summary["seeds"].values()]
    contract = {
        "schema_version": "masld-bench-lsgkm-gse281364-production-activation-contract-v1",
        "status": "inputs_frozen_production_fit_not_executed",
        "family_id": "shared_lsgkm_gkmsvm_deltasvm",
        "activation_config_sha256": file_sha256(arguments.config),
        "tree_authority_artifacts_sha256": {
            name: binding["artifacts_sha256"]
            for name, binding in sorted(config["tree_authorities"].items())
        },
        "training_input_artifacts_sha256": {
            name: binding["artifacts_sha256"]
            for name, binding in sorted(config["training_input_trees"].items())
        },
        "file_authority_sha256": {
            name: binding["sha256"]
            for name, binding in sorted(config["file_authorities"].items())
        },
        "active_state_roster": ["hepatocyte"],
        "excluded_state_roster": {
            "cholangiocyte": "no_matched_GSE281364_reporter_context",
            "fibroblast": "LX2_rows_absent_from_the_frozen_common_row_universe_and_fibroblast_is_not_silently_relabeled_stellate",
            "macrophage": "no_matched_GSE281364_reporter_context",
            "t_cell": "no_matched_GSE281364_reporter_context",
            "endothelial_cell": "not_in_the_frozen_primary_GSE296875_training_pseudobulk_roster_and_no_matched_reporter_context",
            "stellate_cell": "not_an_observed_GSE296875_label;Mesenchymal_is_frozen_as_fibroblast",
        },
        "assay_contexts": list(CONTEXTS),
        "context_specific_prediction": False,
        "context_duplication_rule": "one_static_hepatocyte_sequence_score_is_joined_to_both_HepG2_context_rows_without_claiming_PAOA_modulation",
        "row_universe": {
            "elements": 1_033,
            "contexts": 2,
            "seeds": list(ROW_SEEDS),
            "standardized_rows_per_readout": 10_330,
            "scoring_map_sha256": file_sha256(arguments.output / "contract/scoring_row_map.tsv.gz"),
        },
        "seed_adjudication": {
            "active": list(ROW_SEEDS),
            "superseded_untrained_registry_seeds": list(LEGACY_SEEDS),
            "reason": "exact_join_to_the_frozen_five_seed_GSE281364_row_universe",
        },
        "split_design": {
            "design": "five_diagonal_donor_safe_whole_chromosome_held_partitions",
            "shared_fits": 25,
            "formula": "5_scoring_partitions_x_1_hepatocyte_state_x_5_matching_seeds",
            "reason_not_625": "GSE281364_has_zero_donors_and_one_HepG2_aligned_state;one_whole_chromosome_held_model_per_target_contig_and_seed_is_sufficient_for_locus_safe_static_scoring",
            "factorial_donor_by_genomic_robustness_supported": False,
            "all_25_crossed_partitions_required_for_this_target": False,
            "whole_held_chromosome_is_stricter_than_524kb_target_locus_buffer": True,
        },
        "matched_input_design": {
            "sequence_window_bp": 300,
            "positive_center": "training_only_MACS2_summit",
            "negative_candidate_source": "training_partition_ChromBPNet_candidate_universe_recentered_from_2114bp_to_exact_300bp",
            "negative_candidate_2114bp_intervals_are_final_training_sequences": False,
            "negative_exclusions": [
                "any_training_MACS2_peak_overlap",
                "ENCODE_GRCh38_blacklist_overlap",
                "non_ACGT_sequence",
                "within_class_exact_or_reverse_complement_duplicate",
                "positive_class_exact_or_reverse_complement_duplicate",
                "mappability_below_0.80",
            ],
            "exact_match_strata": [
                "chromosome",
                "training_only_ATAC_decile",
                "GENCODE_v49_nearest_TSS_distance_bin",
            ],
            "continuous_match_tolerances": {
                "absolute_GC_fraction": 0.02,
                "absolute_Umap_k100_multi_read_mappability": 0.05,
                "absolute_RepeatMasker_union_fraction": 0.05,
            },
            "assignment": "deterministic_seeded_greedy_nearest_without_replacement",
            "mappability_window_mean": "exact_bigWig_sum_divided_by_300_with_uncovered_bases_zero",
            "Umap_source_absent_primary_contigs": ["chrY"],
            "source_absent_chrY_windows_pass_0.80_gate": False,
            "atac_covariate": "log1p_training_donor_only_hepatocyte_Tn5_insertion_sum",
            "repeat_covariate": "union_fraction_after_RepeatMasker_one_based_inclusive_to_zero_based_half_open_conversion",
            "minimum_matched_pairs_per_fit": int(config["matching"]["minimum_matched_pairs_per_fit"]),
            "minimum_positive_coverage": float(config["matching"]["minimum_positive_coverage"]),
        },
        "input_summary": summaries,
        "matched_pairs_per_fit_minimum": min(matched_counts),
        "matched_pairs_per_fit_maximum": max(matched_counts),
        "readouts": ["direct_gkmsvm_ALT_minus_REF", "deltasvm_canonical_ALT_minus_REF"],
        "readouts_are_independent_models": False,
        "independent_model_families": 1,
        "production_output_counts": {
            "model_checkpoints": 25,
            "canonical_11mer_tables": 25,
            "unique_element_seed_scores_per_readout": 5_165,
            "standardized_prediction_rows_per_readout": 10_330,
            "standardized_prediction_rows_both_readouts": 20_660,
        },
        "recommended_job_bundle": {
            "jobs": 5,
            "one_job_per_split": True,
            "parallel_fits_per_job": 5,
            "cpus_per_job": 5,
            "memory_gb_per_job": 32,
            "walltime_hours_per_job": 48,
            "total_requested_cpu_hours": 1_200,
            "total_requested_memory_gb_concurrent": 160,
            "resource_estimate_status": "conservative_preproduction_estimate_requires_one_split_scale_timing_probe_before_submission",
        },
        "claim_limits": {
            "development_MPRA_only": True,
            "champion_eligible": False,
            "external_claim_eligible": False,
            "HepG2_is_primary_hepatocyte_equivalent": False,
            "PAOA_condition_specific_effect_supported": False,
            "stellate_or_LX2_effect_supported": False,
            "donor_conditioning_supported": False,
            "target_gene_supported": False,
            "signed_expression_effect_supported": False,
            "eqtl_or_ieqtl_effect_supported": False,
            "causal_gene_supported": False,
            "task_role": "secondary_sequence_accessibility_or_reporter_direction_diagnostic",
        },
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_assets_read": False,
        "production_training_executed": False,
        "production_predictions_generated": False,
    }
    (arguments.output / "contract/production_activation_contract.json").write_text(
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
