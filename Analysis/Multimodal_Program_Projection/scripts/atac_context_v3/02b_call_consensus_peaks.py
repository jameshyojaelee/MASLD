#!/usr/bin/env python3
"""Call condition-blind peaks with explicit two-donor reproducibility support.

SnapATAC2 2.8.0 requires a pooled peak to overlap a peak in every replicate
when ``replicate=`` is supplied.  That behavior is inappropriate when many
donors have shallow lineage pseudobulks.  This producer therefore calls pooled
lineage peaks and eligible donor-lineage peaks separately, then retains pooled
peaks observed in at least two donors.
"""

from __future__ import annotations

import argparse
import bisect
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl
import snapatac2 as snap

from atac_context_v3_lib import (
    COHORTS,
    LINEAGES,
    MIN_CELLS,
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
MIN_SUPPORTING_DONORS = 2
POOLED_QVALUE = 0.05
DONOR_QVALUE = 0.10

SUPPORT_COLUMNS = (
    "release_id",
    "cohort",
    "lineage",
    "chrom",
    "start0",
    "end",
    "pooled_qvalue",
    "eligible_donors",
    "min_supporting_donors",
    "n_supporting_donors",
    "supporting_donors",
    "retained",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, default=default_candidate_root())
    parser.add_argument("--n-jobs", type=int, default=16)
    parser.add_argument("--tempdir", type=Path, required=True)
    return parser.parse_args()


def normalize_peak_frame(frame: pl.DataFrame) -> pl.DataFrame:
    rename = {}
    for old, new in {"chromosome": "chrom", "start": "start0"}.items():
        if old in frame.columns and new not in frame.columns:
            rename[old] = new
    frame = frame.rename(rename)
    required = {"chrom", "start0", "end"}
    if not required.issubset(frame.columns):
        raise ContractError(f"MACS3 peak columns missing: {required - set(frame.columns)}")
    return frame.filter(pl.col("chrom").is_in(STANDARD_CHROMS)).sort(
        ["chrom", "start0", "end"]
    )


def write_native(path: Path, frame: pl.DataFrame) -> None:
    require_new_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        handle.write("\t".join(frame.columns) + "\n")
        for row in frame.iter_rows():
            handle.write("\t".join(str(value) for value in row) + "\n")


def write_bed(path: Path, frame: pl.DataFrame) -> None:
    require_new_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        for chrom, start, end in frame.select(["chrom", "start0", "end"]).iter_rows():
            handle.write(f"{chrom}\t{start}\t{end}\n")


def interval_index(frame: pl.DataFrame) -> dict[str, tuple[list[int], list[int]]]:
    result: dict[str, tuple[list[int], list[int]]] = {}
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


def overlap_flags(
    query: pl.DataFrame, subject: pl.DataFrame | dict[str, tuple[list[int], list[int]]]
) -> list[bool]:
    indexed = interval_index(subject) if isinstance(subject, pl.DataFrame) else subject
    flags = []
    for chrom, start, end in query.select(["chrom", "start0", "end"]).iter_rows():
        starts, prefix_max = indexed.get(str(chrom), ([], []))
        upper = bisect.bisect_left(starts, int(end))
        flags.append(upper > 0 and prefix_max[upper - 1] > int(start))
    return flags


def eligible_donor_groups(dataset) -> dict[str, list[str]]:
    donors = np.asarray(dataset.obs["donor_id"]).astype(str).tolist()
    lineages = np.asarray(dataset.obs["lineage_v3"]).astype(str).tolist()
    counts = Counter(zip(donors, lineages))
    groups = {lineage: [] for lineage in LINEAGES}
    group_values = []
    for donor, lineage in zip(donors, lineages):
        if lineage in groups and counts[(donor, lineage)] >= MIN_CELLS:
            key = f"{lineage}__{donor}"
            groups[lineage].append(key)
            group_values.append(key)
        else:
            group_values.append("excluded")
    dataset.obs["donor_lineage_v3"] = np.asarray(group_values, dtype=str)
    return {lineage: sorted(set(values)) for lineage, values in groups.items()}


def call_supported_native_peaks(dataset, tempdir: Path, n_jobs: int):
    groups = eligible_donor_groups(dataset)
    for lineage, values in groups.items():
        if len(values) < MIN_SUPPORTING_DONORS:
            raise ContractError(f"fewer than two eligible peak-support donors: {lineage}")
    pooled = snap.tl.macs3(
        dataset,
        groupby="lineage_v3",
        selections=set(LINEAGES),
        qvalue=POOLED_QVALUE,
        call_broad_peaks=False,
        blacklist=BLACKLIST,
        shift=-100,
        extsize=200,
        tempdir=tempdir,
        inplace=False,
        n_jobs=n_jobs,
    )
    donor = snap.tl.macs3(
        dataset,
        groupby="donor_lineage_v3",
        selections={value for values in groups.values() for value in values},
        qvalue=DONOR_QVALUE,
        call_broad_peaks=False,
        blacklist=BLACKLIST,
        shift=-100,
        extsize=200,
        tempdir=tempdir,
        inplace=False,
        n_jobs=n_jobs,
    )
    if pooled is None or donor is None:
        raise ContractError("MACS3 returned no peak dictionary")

    retained: dict[str, pl.DataFrame] = {}
    audit_rows: list[dict[str, object]] = []
    for lineage in LINEAGES:
        if lineage not in pooled:
            raise ContractError(f"MACS3 omitted pooled lineage: {lineage}")
        frame = normalize_peak_frame(pooled[lineage])
        support_flags: dict[str, list[bool]] = {}
        for group in groups[lineage]:
            if group not in donor:
                raise ContractError(f"MACS3 omitted donor-lineage group: {group}")
            support_flags[group] = overlap_flags(frame, normalize_peak_frame(donor[group]))
        keep = []
        for index, (chrom, start, end) in enumerate(
            frame.select(["chrom", "start0", "end"]).iter_rows()
        ):
            supporters = [
                group.rsplit("__", 1)[1]
                for group in groups[lineage]
                if support_flags[group][index]
            ]
            retain = len(supporters) >= MIN_SUPPORTING_DONORS
            keep.append(retain)
            pooled_qvalue = frame[index, "q_value"] if "q_value" in frame.columns else "NA"
            audit_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "cohort": "",
                    "lineage": lineage,
                    "chrom": chrom,
                    "start0": start,
                    "end": end,
                    "pooled_qvalue": pooled_qvalue,
                    "eligible_donors": len(groups[lineage]),
                    "min_supporting_donors": MIN_SUPPORTING_DONORS,
                    "n_supporting_donors": len(supporters),
                    "supporting_donors": ",".join(supporters),
                    "retained": str(retain).upper(),
                }
            )
        retained[lineage] = frame.filter(pl.Series("retained", keep))
        if retained[lineage].height == 0:
            raise ContractError(f"no two-donor-supported peaks retained: {lineage}")
    return retained, audit_rows


def load_blacklist() -> pl.DataFrame:
    return pl.read_csv(
        BLACKLIST,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start0", "end", "name"],
        truncate_ragged_lines=True,
    ).select(["chrom", "start0", "end"]).filter(pl.col("chrom").is_in(STANDARD_CHROMS))


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
    support_rows = []
    for cohort in COHORTS:
        dataset_file = args.tempdir / f"{cohort}.peaks.h5ads"
        dataset, adatas = build_dataset(cohort, dataset_file)
        try:
            cohort_native, cohort_support = call_supported_native_peaks(
                dataset, args.tempdir, args.n_jobs
            )
            for row in cohort_support:
                row["cohort"] = cohort
            support_rows.extend(cohort_support)
            for lineage, frame in cohort_native.items():
                native[(cohort, lineage)] = frame
                write_native(target / "native" / f"{cohort}.{lineage}.tsv", frame)
                write_bed(target / "native" / f"{cohort}.{lineage}.bed", frame)
        finally:
            close_dataset(dataset, adatas)
            if dataset_file.exists():
                dataset_file.unlink()

    write_tsv(target / "native_peak_support.tsv", SUPPORT_COLUMNS, support_rows)
    blacklist = load_blacklist()
    manifest_rows = []
    for lineage in LINEAGES:
        merged = normalize_peak_frame(
            snap.tl.merge_peaks(
                {
                    "GSE244832": native[("GSE244832", lineage)],
                    "GSE281367": native[("GSE281367", lineage)],
                },
                snap.genome.hg38,
                half_width=250,
            )
        )
        merged = merged.filter(~pl.Series("blacklisted", overlap_flags(merged, blacklist)))
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
    print(f"Wrote condition-blind two-donor-supported consensus peaks: {target}")


if __name__ == "__main__":
    main()
