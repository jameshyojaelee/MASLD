#!/usr/bin/env python3
"""Read-only Currin split/identity audit; no effect values or predictions used.

Coordinates are zero-based half-open variant-centered input windows, using the
same formula as i1_common.window_bounds. Sequence identity between different
chromosomes and donor-level LD are not measured by interval overlap.
"""
from __future__ import annotations

import argparse
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915"
LABELS = ROOT / "GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz"
FOLDS = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/fold_manifest.tsv"
FAI = Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa.fai")
FIELDS = ("key", "chr", "pos_hg38", "ref", "alt", "peak_id", "peak_start_hg38", "peak_stop_hg38", "heldout_fold", "block_1mb")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        result = []
        for row in csv.DictReader(handle, delimiter="\t"):
            row["key"] = row.get("key", row.get("lead_variant_id", "").removeprefix("chr"))
            # Deliberately discard every effect/p-value column from C2 labels.
            row = {k: row[k] for k in FIELDS}
            for k in ("pos_hg38", "peak_start_hg38", "peak_stop_hg38", "heldout_fold"):
                row[k] = int(row[k])
            result.append(row)
    return result


def write_table(path, rows, columns=None):
    rows = list(rows)
    columns = columns or (list(rows[0]) if rows else [])
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def identity_summary(name, rows, authority):
    results = []
    for field, make_key in (
        ("variant", lambda r: r["key"]),
        ("variant_target", lambda r: (r["key"], r["peak_id"])),
        ("peak_target", lambda r: (r["chr"], r["peak_id"])),
        ("chromosome", lambda r: r["chr"]),
        ("historical_bootstrap_1mb_bin", lambda r: r["block_1mb"]),
    ):
        counts, foldsets = Counter(), defaultdict(set)
        for row in rows:
            key = make_key(row)
            counts[key] += 1
            foldsets[key].add(row["heldout_fold"])
        results.append(dict(population=name, identity=field, rows=len(rows), distinct=len(counts),
            repeated_identities=sum(v > 1 for v in counts.values()),
            rows_in_repeated_identities=sum(v for v in counts.values() if v > 1),
            cross_fold_identities=sum(len(x) > 1 for x in foldsets.values()),
            fold_authority_disagreements=sum(r["heldout_fold"] != authority.get(r["chr"]) for r in rows)))
    return results


def interval_degrees(rows, length, chrom_lengths):
    """O(n log n): half-open overlap degrees, independently of fold grouping."""
    groups = defaultdict(list)
    intervals = []
    for i, row in enumerate(rows):
        raw_start = row["pos_hg38"] - 1 - length // 2
        raw_end = raw_start + length
        # Out-of-bounds requests cannot create real genomic overlap in padding.
        start, end = max(0, raw_start), min(chrom_lengths[row["chr"]], raw_end)
        intervals.append((start, end, raw_start, raw_end))
        groups[row["chr"]].append(i)
    degrees = []
    index = {}
    for chrom, indices in groups.items():
        for field in ("heldout_fold", "block_1mb"):
            subgroups = defaultdict(list)
            for i in indices:
                subgroups[rows[i][field]].append(i)
            for value, members in subgroups.items():
                index[(chrom, field, value)] = (
                    sorted(intervals[i][0] for i in members), sorted(intervals[i][1] for i in members))
        index[(chrom, "all", "all")] = (
            sorted(intervals[i][0] for i in indices), sorted(intervals[i][1] for i in indices))

    def count(chrom, field, value, start, end):
        starts, ends = index[(chrom, field, value)]
        return bisect_left(starts, end) - bisect_right(ends, start)

    for i, row in enumerate(rows):
        start, end, raw_start, raw_end = intervals[i]
        total = count(row["chr"], "all", "all", start, end) - 1
        same = count(row["chr"], "heldout_fold", row["heldout_fold"], start, end) - 1
        same_bin = count(row["chr"], "block_1mb", row["block_1mb"], start, end) - 1
        degrees.append(dict(key=row["key"], heldout_fold=row["heldout_fold"], chromosome=row["chr"],
            window_start0=raw_start, window_end0=raw_end,
            out_of_bounds=int(raw_start < 0 or raw_end > chrom_lengths[row["chr"]]),
            overlapping_rows=total, same_fold_overlapping_rows=same,
            cross_fold_overlapping_rows=total - same, cross_bootstrap_bin_overlapping_rows=total - same_bin))
    return degrees


def self_tests():
    rows = [dict(key=str(i), chr="chr1", pos_hg38=p, heldout_fold=f, block_1mb=b)
            for i, (p, f, b) in enumerate(((6, 0, "a"), (10, 1, "a"), (16, 0, "b"), (17, 1, "b")))]
    result = interval_degrees(rows, 10, {"chr1": 100})
    # [0,10) and [10,20) just touch: no overlap; [4,14) overlaps both.
    expected = [[], [], [], []]
    for i in range(len(rows)):
        for j in range(len(rows)):
            if i != j and max(result[i]["window_start0"], result[j]["window_start0"]) < min(result[i]["window_end0"], result[j]["window_end0"]):
                expected[i].append(j)
    for i, rec in enumerate(result):
        assert rec["overlapping_rows"] == len(expected[i])
        assert rec["cross_fold_overlapping_rows"] == sum(rows[i]["heldout_fold"] != rows[j]["heldout_fold"] for j in expected[i])
        assert rec["cross_bootstrap_bin_overlapping_rows"] == sum(rows[i]["block_1mb"] != rows[j]["block_1mb"] for j in expected[i])
    rows[-1]["chr"] = "chr2"
    assert interval_degrees(rows, 10, {"chr1": 100, "chr2": 100})[-1]["overlapping_rows"] == 0
    return "half_open_boundaries_cross_fold_cross_bin_and_chromosome_identity_passed"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run executable scientific checks on a compute node through the authorized launcher")
    args.out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, args.out / "executed_data_overlap_audit.py")
    checks = [self_tests()]
    authority = {}
    with FOLDS.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for chrom in row["test_chromosomes"].split(","):
                if chrom in authority:
                    raise ValueError("Repeated chromosome in fold authority")
                authority[chrom] = int(row["fold"])
    with FAI.open() as handle:
        chrom_lengths = {parts[0]: int(parts[1]) for parts in (line.split("\t") for line in handle)}
    paths = {"original_C2": LABELS,
        "frozen_2048_manifest": BASE / "model/frozen/2048/manifest.tsv",
        "frozen_16384_manifest": BASE / "model/frozen/16384/manifest.tsv",
        "native_manifest": BASE / "model/native/sequence_manifest.tsv"}
    populations = {name: read_rows(path) for name, path in paths.items()}
    reference = populations["original_C2"]
    agreements = []
    for name, rows in populations.items():
        agreements.append(dict(population=name, rows=len(rows), expected_rows=32322,
            same_identity_and_order_as_C2=rows == reference,
            unique_variant_keys=len({r["key"] for r in rows})))
    write_table(args.out / "manifest_agreement.tsv", agreements)
    identity = []
    for name, rows in populations.items():
        identity.extend(identity_summary(name, rows, authority))
    if not all(row["same_identity_and_order_as_C2"] and row["rows"] == 32322 for row in agreements):
        raise ValueError("Manifest identities differ; written agreement table identifies the mismatch")
    coverage_rows, coverage_paths = [], []
    matched = np.ones(len(reference), dtype=bool)
    for length in (2048, 16384):
        path = BASE / f"model/frozen/{length}/coverage.npz"
        coverage_paths.append(path)
        with np.load(path, allow_pickle=False) as archive:
            for pool in ("variant", "symmetric", "target"):
                mask = archive[pool]
                if mask.shape != matched.shape or mask.dtype != np.bool_:
                    raise ValueError("Coverage shape or type differs")
                matched &= mask
                coverage_rows.append(dict(length=length, pooling=pool, admitted_rows=int(mask.sum()), requested_rows=len(reference)))
    populations["all_six_frozen_views_common"] = [r for r, keep in zip(reference, matched) if keep]
    identity.extend(identity_summary("all_six_frozen_views_common", populations["all_six_frozen_views_common"], authority))
    write_table(args.out / "identity_and_fold_audit.tsv", identity)
    write_table(args.out / "length_pooling_coverage.tsv", coverage_rows)
    summaries, by_fold = [], []
    for name in ("original_C2", "all_six_frozen_views_common"):
        rows = populations[name]
        for length in (2048, 16384, 1048576):
            degree = interval_degrees(rows, length, chrom_lengths)
            write_table(args.out / f"{name}_{length}_overlap_degrees.tsv.gz", degree)
            n = len(rows)
            total_pairs = sum(r["overlapping_rows"] for r in degree) // 2
            cross_pairs = sum(r["cross_fold_overlapping_rows"] for r in degree) // 2
            cross_rows = sum(r["cross_fold_overlapping_rows"] > 0 for r in degree)
            summaries.append(dict(population=name, input_length_bp=length, rows=n,
                overlapping_pairs=total_pairs, cross_fold_overlapping_pairs=cross_pairs,
                cross_fold_pair_fraction_among_overlapping=cross_pairs / total_pairs if total_pairs else 0.,
                rows_with_any_overlap=sum(r["overlapping_rows"] > 0 for r in degree),
                rows_with_cross_fold_overlap=cross_rows, cross_fold_row_fraction=cross_rows / n if n else 0.,
                cross_bootstrap_bin_overlapping_pairs=sum(r["cross_bootstrap_bin_overlapping_rows"] for r in degree) // 2,
                rows_with_cross_bootstrap_bin_overlap=sum(r["cross_bootstrap_bin_overlapping_rows"] > 0 for r in degree),
                out_of_bounds_requested_windows=sum(r["out_of_bounds"] for r in degree)))
            for fold in range(5):
                held = [r for r in degree if r["heldout_fold"] == fold]
                by_fold.append(dict(population=name, length_bp=length, held_fold=fold,
                    held_rows=len(held), held_rows_overlapping_other_fold=sum(r["cross_fold_overlapping_rows"] > 0 for r in held),
                    held_vs_other_fold_overlap_pairs=sum(r["cross_fold_overlapping_rows"] for r in held)))
    write_table(args.out / "window_overlap_summary.tsv", summaries)
    write_table(args.out / "window_overlap_by_held_fold.tsv", by_fold)
    exposure = [
        dict(scope="Currin_population", observation="32322 historically significance-selected and q-value-deduplicated source leads; same publication and donor cohort throughout folds",
             consequence="Development comparison; chromosome holdout does not create independent donors or independent studies", source="scripts/analysis/alphagenome_program/c2_01_build_variants.py:118-151"),
        dict(scope="fold_authority", observation="Established hepatocyte_5fold_v2 chromosome test groups are the C2 fold authority",
             consequence="Any cis-LD relationship on a chromosome remains in one fold; pairwise donor LD is not reconstructed", source=str(FOLDS.relative_to(ROOT))),
        dict(scope="LD_and_bootstrap", observation="Historical block_1mb is chr:floor(position1/1000000); no measured-LD-block partition for these32322 variants is supplied by C2",
             consequence="Bootstrap bins are not validated independent LD units; adjacent-bin windows may overlap and long-range LD may cross bins", source="scripts/analysis/alphagenome_program/c2_01_build_variants.py:160-164"),
        dict(scope="frozen_native", observation="Local AlphaGenome all-folds checkpoint has encoder reference-locus exposure possible or seen; declared ENCODE/GTEx/FANTOM5/4DN track corpus",
             consequence="Coordinate-separated evaluation is held for downstream label fitting, not established held genomic sequence during pretraining", source="Analysis/MASLD_Model_Benchmark/config/artifacts/models/alphagenome/exposure_audit.json:16-25;docs/technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md:112"),
        dict(scope="Currin_foundation_exposure", observation="No checkpoint-bound donor/accession training manifest establishes Currin donor or exact target-label nonexposure",
             consequence="Currin pretraining source overlap remains unresolved; do not relabel this exposed development source as protected", source="GWAS/finemapping/results/alphagenome_program/a2-rights-exposure-20260914T193626Z/tables/rights_exposure.tsv:163"),
        dict(scope="frozen_representations", observation="Extracted allele features are fixed before fitting; scaling and scalar heads use training folds in the inspected fit code",
             consequence="All-folds encoder exposure is separate from downstream label leakage; the audit does not verify every fitted checkpoint gradient", source="scripts/analysis/alphagenome_campaign/model_frozen_fit.py:62-76"),
        dict(scope="adapter_representations", observation="Trainable adapters use whole-chromosome held folds with representations recomputed during training",
             consequence="Current coordinate overlap audit covers direct genomic input overlap; it does not prove absence of homologous sequence or source-donor overlap", source="scripts/analysis/alphagenome_campaign/model_train.py:50-56"),
    ]
    write_table(args.out / "source_and_pretraining_exposure.tsv", exposure)
    receipt = dict(status="completed_identity_interval_and_documented_exposure_audit", job_id=os.environ["SLURM_JOB_ID"],
        python=sys.version, numpy=np.__version__, node=platform.node(), checks=checks,
        input_sha256={str(p): digest(p) for p in [*paths.values(), *coverage_paths, FOLDS, FAI, Path(__file__)]},
        requested_rows=len(reference), common_frozen_rows=int(matched.sum()),
        no_effect_values_used=True, no_predictions_used=True, no_protected_outcomes_read=True,
        native_intervals="1048576bp requested centered windows; actual native score coverage is a separate result",
        boundary_convention="[pos1-1-length/2,pos1-1+length/2); overlap only of physical genomic interval",
        limitations=["Homologous sequence across chromosomes untested", "Donor LD not measured;1Mb bins are not LD blocks",
                     "Foundation source and reference-locus exposure not resolved by downstream folds",
                     "Intervals do not prove independent biological replication or generalization to another study"])
    (args.out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(dict(receipt=receipt, overlap=summaries), indent=2), flush=True)


if __name__ == "__main__":
    main()
