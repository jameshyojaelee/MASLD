#!/usr/bin/env python3
"""Call condition-blind donor-replicated peaks and build lineage consensus peaks."""

from __future__ import annotations

import argparse
import bisect
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
from snap_helpers import build_dataset, close_dataset


BLACKLIST = Path("/gpfs/commons/home/jameslee/reference_genome/blacklists/hg38-blacklist.v2.bed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, default=default_candidate_root())
    parser.add_argument("--n-jobs", type=int, default=16)
    parser.add_argument("--tempdir", type=Path, required=True)
    return parser.parse_args()


def normalize_peak_frame(frame: pl.DataFrame) -> pl.DataFrame:
    rename = {}
    for old, new in {
        "chromosome": "chrom",
        "start": "start0",
        "end": "end",
    }.items():
        if old in frame.columns and new not in frame.columns:
            rename[old] = new
    frame = frame.rename(rename)
    required = {"chrom", "start0", "end"}
    if not required.issubset(frame.columns):
        raise ContractError(f"MACS3 peak columns missing: {required - set(frame.columns)}")
    return frame.filter(pl.col("chrom").is_in(STANDARD_CHROMS)).sort(["chrom", "start0", "end"])


def write_native(path: Path, frame: pl.DataFrame) -> None:
    require_new_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = frame.columns
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in frame.iter_rows():
            handle.write("\t".join(str(value) for value in row) + "\n")


def write_bed(path: Path, frame: pl.DataFrame) -> None:
    require_new_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in frame.select(["chrom", "start0", "end"]).iter_rows():
            handle.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")


def overlap_flags(consensus: pl.DataFrame, native: pl.DataFrame) -> list[bool]:
    by_chrom: dict[str, tuple[list[int], list[int]]] = {}
    for chrom in STANDARD_CHROMS:
        rows = native.filter(pl.col("chrom") == chrom).select(["start0", "end"]).iter_rows()
        ordered = sorted((int(start), int(end)) for start, end in rows)
        starts = [value[0] for value in ordered]
        prefix_max = []
        running = -1
        for _, end in ordered:
            running = max(running, end)
            prefix_max.append(running)
        by_chrom[chrom] = (starts, prefix_max)
    result = []
    for chrom, start, end in consensus.select(["chrom", "start0", "end"]).iter_rows():
        starts, prefix_max = by_chrom.get(str(chrom), ([], []))
        upper = bisect.bisect_left(starts, int(end))
        result.append(upper > 0 and prefix_max[upper - 1] > int(start))
    return result


def load_blacklist() -> pl.DataFrame:
    frame = pl.read_csv(
        BLACKLIST,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start0", "end", "name"],
        truncate_ragged_lines=True,
    )
    return frame.select(["chrom", "start0", "end"]).filter(
        pl.col("chrom").is_in(STANDARD_CHROMS)
    )


def main() -> None:
    args = parse_args()
    out = candidate_root(args.candidate_root)
    if not (out / "input_manifest.tsv").is_file():
        raise ContractError("input freeze must complete before peak calling")
    target = out / "peaks"
    require_new_path(target)
    target.mkdir(parents=True)
    args.tempdir.mkdir(parents=True, exist_ok=True)

    native: dict[tuple[str, str], pl.DataFrame] = {}
    for cohort in COHORTS:
        dataset_file = args.tempdir / f"{cohort}.peaks.h5ads"
        dataset, adatas = build_dataset(cohort, dataset_file)
        try:
            result = snap.tl.macs3(
                dataset,
                groupby="lineage_v3",
                replicate="donor_id",
                selections=set(LINEAGES),
                qvalue=0.05,
                replicate_qvalue=0.05,
                call_broad_peaks=False,
                blacklist=BLACKLIST,
                shift=-100,
                extsize=200,
                tempdir=args.tempdir,
                inplace=False,
                n_jobs=args.n_jobs,
            )
            if result is None:
                raise ContractError(f"MACS3 returned no peak dictionary for {cohort}")
            for lineage in LINEAGES:
                if lineage not in result:
                    raise ContractError(f"MACS3 omitted {cohort}/{lineage}")
                frame = normalize_peak_frame(result[lineage])
                native[(cohort, lineage)] = frame
                write_native(target / "native" / f"{cohort}.{lineage}.tsv", frame)
                write_bed(target / "native" / f"{cohort}.{lineage}.bed", frame)
        finally:
            close_dataset(dataset, adatas)
            if dataset_file.exists():
                dataset_file.unlink()

    blacklist = load_blacklist()
    manifest_rows = []
    for lineage in LINEAGES:
        merged = snap.tl.merge_peaks(
            {
                "GSE244832": native[("GSE244832", lineage)],
                "GSE281367": native[("GSE281367", lineage)],
            },
            snap.genome.hg38,
            half_width=250,
        )
        merged = normalize_peak_frame(merged)
        blacklisted = overlap_flags(merged, blacklist)
        merged = merged.filter(~pl.Series("blacklisted", blacklisted))
        flags244 = overlap_flags(merged, native[("GSE244832", lineage)])
        flags281 = overlap_flags(merged, native[("GSE281367", lineage)])
        write_bed(target / "consensus" / f"{lineage}.bed", merged)
        for index, ((chrom, start, end), f244, f281) in enumerate(
            zip(merged.select(["chrom", "start0", "end"]).iter_rows(), flags244, flags281),
            start=1,
        ):
            width = int(end) - int(start)
            if width != 500:
                raise ContractError(f"consensus peak width drift: {lineage} {chrom}:{start}-{end}")
            manifest_rows.append(
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
    write_tsv(out / "consensus_peak_manifest.tsv", PEAK_MANIFEST_COLUMNS, manifest_rows)
    print(f"Wrote condition-blind consensus peaks: {target}")


if __name__ == "__main__":
    main()
