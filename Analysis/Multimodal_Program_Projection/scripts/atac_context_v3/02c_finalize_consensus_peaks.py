#!/usr/bin/env python3
"""Finalize consensus peaks from completed, audited cohort-native peak calls."""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl
import snapatac2 as snap

from atac_context_v3_lib import (
    COHORTS,
    LINEAGES,
    PEAK_MANIFEST_COLUMNS,
    RELEASE_ID,
    STANDARD_CHROMS,
    ContractError,
    candidate_root,
    default_candidate_root,
    require_new_path,
    write_tsv,
)


BLACKLIST = Path("/gpfs/commons/home/jameslee/reference_genome/blacklists/hg38-blacklist.v2.bed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, default=default_candidate_root())
    return parser.parse_args()


def normalize_peak_frame(frame: pl.DataFrame) -> pl.DataFrame:
    rename = {}
    for old, new in {"chromosome": "chrom", "start": "start0"}.items():
        if old in frame.columns and new not in frame.columns:
            rename[old] = new
    frame = frame.rename(rename)
    required = {"chrom", "start0", "end"}
    if not required.issubset(frame.columns):
        raise ContractError(f"peak columns missing: {required - set(frame.columns)}")
    return frame.filter(pl.col("chrom").is_in(STANDARD_CHROMS)).sort(
        ["chrom", "start0", "end"]
    )


def merge_input_frame(frame: pl.DataFrame) -> pl.DataFrame:
    normalized = normalize_peak_frame(frame)
    if "start" in normalized.columns:
        raise ContractError("ambiguous merge schema contains both start and start0")
    return normalized.rename({"start0": "start"}).with_columns(
        pl.col("start").cast(pl.UInt64),
        pl.col("end").cast(pl.UInt64),
        pl.col("score").cast(pl.UInt16),
        pl.col("signal_value").cast(pl.Float64),
        pl.col("p_value").cast(pl.Float64),
        pl.col("q_value").cast(pl.Float64),
        pl.col("peak").cast(pl.UInt64),
    )


def normalize_merged_peak_frame(frame: pl.DataFrame) -> pl.DataFrame:
    """Convert SnapATAC2's range strings to the v3 500-bp BED schema."""
    if "Peaks" not in frame.columns:
        raise ContractError("SnapATAC2 merged output is missing Peaks")
    parsed = frame.with_columns(
        pl.col("Peaks").str.extract(r"^([^:]+):(\d+)-(\d+)$", 1).alias("chrom"),
        pl.col("Peaks").str.extract(r"^([^:]+):(\d+)-(\d+)$", 2).cast(pl.UInt64).alias("start0"),
        pl.col("Peaks").str.extract(r"^([^:]+):(\d+)-(\d+)$", 3).cast(pl.UInt64).alias("raw_end"),
    )
    if parsed.select(
        pl.any_horizontal(pl.col(["chrom", "start0", "raw_end"]).is_null()).any()
    ).item():
        raise ContractError("failed to parse SnapATAC2 merged peak coordinates")
    parsed = parsed.with_columns((pl.col("start0") + 500).alias("end"))
    if parsed.filter(pl.col("raw_end") - pl.col("start0") != 501).height:
        raise ContractError("SnapATAC2 returned a clipped or non-501-bp merged peak")
    return parsed.select(["chrom", "start0", "end"]).sort(["chrom", "start0", "end"])


def interval_index(frame: pl.DataFrame) -> dict[str, tuple[list[int], list[int]]]:
    result = {}
    for chrom in STANDARD_CHROMS:
        ordered = sorted(
            (int(start), int(end))
            for start, end in frame.filter(pl.col("chrom") == chrom)
            .select(["start0", "end"])
            .iter_rows()
        )
        starts = [value[0] for value in ordered]
        prefix_max = []
        running = -1
        for _, end in ordered:
            running = max(running, end)
            prefix_max.append(running)
        result[chrom] = (starts, prefix_max)
    return result


def overlap_flags(query: pl.DataFrame, subject: pl.DataFrame) -> list[bool]:
    from bisect import bisect_left

    index = interval_index(subject)
    flags = []
    for chrom, start, end in query.select(["chrom", "start0", "end"]).iter_rows():
        starts, prefix_max = index.get(chrom, ([], []))
        position = bisect_left(starts, int(end)) - 1
        flags.append(position >= 0 and prefix_max[position] > int(start))
    return flags


def write_bed(path: Path, frame: pl.DataFrame) -> None:
    require_new_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        for chrom, start, end in frame.select(["chrom", "start0", "end"]).iter_rows():
            handle.write(f"{chrom}\t{start}\t{end}\n")


def main() -> None:
    args = parse_args()
    out = candidate_root(args.candidate_root)
    manifest_path = out / "consensus_peak_manifest.tsv"
    require_new_path(manifest_path)
    support_path = out / "peaks/native_peak_support.tsv"
    if not support_path.is_file():
        raise ContractError("completed native support audit is required")

    native = {}
    for cohort in COHORTS:
        for lineage in LINEAGES:
            path = out / "peaks/native" / f"{cohort}.{lineage}.tsv"
            if not path.is_file():
                raise ContractError(f"missing completed native peak table: {path}")
            native[(cohort, lineage)] = normalize_peak_frame(
                pl.read_csv(path, separator="\t")
            )

    blacklist = pl.read_csv(
        BLACKLIST,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start0", "end", "name"],
        truncate_ragged_lines=True,
    ).select(["chrom", "start0", "end"]).filter(pl.col("chrom").is_in(STANDARD_CHROMS))

    rows = []
    for lineage in LINEAGES:
        merged = normalize_merged_peak_frame(
            snap.tl.merge_peaks(
                {
                    cohort: merge_input_frame(native[(cohort, lineage)])
                    for cohort in COHORTS
                },
                snap.genome.hg38,
                half_width=250,
            )
        )
        merged = merged.filter(~pl.Series("blacklisted", overlap_flags(merged, blacklist)))
        if merged.select(["chrom", "start0", "end"]).is_duplicated().any():
            raise ContractError(f"duplicated consensus coordinates: {lineage}")
        flags244 = overlap_flags(merged, native[("GSE244832", lineage)])
        flags281 = overlap_flags(merged, native[("GSE281367", lineage)])
        write_bed(out / "peaks/consensus" / f"{lineage}.bed", merged)
        for index, ((chrom, start, end), f244, f281) in enumerate(
            zip(merged.select(["chrom", "start0", "end"]).iter_rows(), flags244, flags281),
            start=1,
        ):
            width = int(end) - int(start)
            if width != 500:
                raise ContractError(f"consensus peak width drift: {lineage} {chrom}:{start}-{end}")
            rows.append(
                {
                    "release_id": RELEASE_ID,
                    "lineage": lineage,
                    "peak_id": f"{lineage}_peak_{index:07d}",
                    "chrom": chrom,
                    "start0": start,
                    "end": end,
                    "width": width,
                    "supported_gse244832": str(bool(f244)).upper(),
                    "supported_gse281367": str(bool(f281)).upper(),
                    "supported_both": str(bool(f244 and f281)).upper(),
                    "blacklist_overlap": "FALSE",
                }
            )
    write_tsv(manifest_path, PEAK_MANIFEST_COLUMNS, rows)
    print(f"Finalized consensus peaks from completed native calls: {manifest_path}")


if __name__ == "__main__":
    main()
