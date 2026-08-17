#!/usr/bin/env python3
"""Independent consensus-only validation for the peak recovery boundary."""

from __future__ import annotations

from pathlib import Path

import polars as pl


ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CANDIDATE = ROOT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
BLACKLIST = Path("/gpfs/commons/home/jameslee/reference_genome/blacklists/hg38-blacklist.v2.bed")
LINEAGES = ("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk")
STANDARD = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def has_overlap(query: pl.DataFrame, subject: pl.DataFrame) -> bool:
    for chrom in STANDARD:
        left = query.filter(pl.col("chrom") == chrom).sort("start0")
        right = subject.filter(pl.col("chrom") == chrom).sort("start0")
        right_rows = list(right.select(["start0", "end"]).iter_rows())
        j = 0
        for start, end in left.select(["start0", "end"]).iter_rows():
            while j < len(right_rows) and right_rows[j][1] <= start:
                j += 1
            if j < len(right_rows) and right_rows[j][0] < end:
                return True
    return False


manifest = pl.read_csv(CANDIDATE / "consensus_peak_manifest.tsv", separator="\t")
assert manifest.height > 0
assert set(manifest["lineage"].unique()) == set(LINEAGES)
assert set(manifest["chrom"].unique()).issubset(STANDARD)
assert manifest.filter(pl.col("width") != 500).is_empty()
assert manifest.filter(pl.col("end") - pl.col("start0") != 500).is_empty()
assert manifest.filter(pl.col("blacklist_overlap") != False).is_empty()
assert manifest.filter(~(pl.col("supported_gse244832") | pl.col("supported_gse281367"))).is_empty()
assert manifest.filter(
    pl.col("supported_both") != (pl.col("supported_gse244832") & pl.col("supported_gse281367"))
).is_empty()
assert not manifest.select(["lineage", "chrom", "start0", "end"]).is_duplicated().any()

blacklist = pl.read_csv(
    BLACKLIST,
    separator="\t",
    has_header=False,
    new_columns=["chrom", "start0", "end", "name"],
    truncate_ragged_lines=True,
).select(["chrom", "start0", "end"]).filter(pl.col("chrom").is_in(STANDARD))

summary = []
for lineage in LINEAGES:
    observed = manifest.filter(pl.col("lineage") == lineage).sort(["chrom", "start0", "end"])
    prior = observed.with_columns(pl.col("end").shift(1).over("chrom").alias("prior_end"))
    assert prior.filter(pl.col("prior_end").is_not_null() & (pl.col("start0") < pl.col("prior_end"))).is_empty()
    assert not has_overlap(observed.select(["chrom", "start0", "end"]), blacklist)
    bed = pl.read_csv(
        CANDIDATE / "peaks/consensus" / f"{lineage}.bed",
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start0", "end"],
    )
    expected = observed.select(["chrom", "start0", "end"])
    assert bed.equals(expected), lineage
    summary.append(
        {
            "lineage": lineage,
            "n_peaks": observed.height,
            "n_shared": observed.filter(pl.col("supported_both")).height,
        }
    )

print(pl.DataFrame(summary))
print(f"Independent consensus recovery validation passed: {manifest.height} peaks")
