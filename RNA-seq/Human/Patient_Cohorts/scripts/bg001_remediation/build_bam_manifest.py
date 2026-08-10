#!/usr/bin/env python3
"""Build the frozen BG-001 BAM manifest from the five affected cohorts."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


EXPECTED = {
    "GSE130970": 78,
    "GSE135251": 216,
    "GSE174478": 93,
    "GSE213621": 367,
    "GSE240729": 66,
}
DOCUMENTED_MISSING = {
    "GSE174478": ("SRR14551000", "documented_failed_run"),
    "GSE240729": ("SRR25630203", "documented_failed_run"),
}
SUFFIX = ".Aligned.sortedByCoord.out.bam"
FIELDS = (
    "dataset",
    "sample_id",
    "bam_path",
    "layout",
    "strandedness",
    "expected_status",
    "exclusion_reason",
    "bam_size",
    "bam_mtime",
    "bam_sha256",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def sample_from_bam(path: Path) -> str:
    if not path.name.endswith(SUFFIX):
        raise ValueError(f"Unexpected BAM name: {path}")
    return path.name[: -len(SUFFIX)]


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve(strict=True)
    patient_root = project_root / "RNA-seq/Human/Patient_Cohorts"
    rows: list[dict[str, str | int]] = []

    for dataset, expected_count in EXPECTED.items():
        bam_root = patient_root / "results" / dataset / "alignments" / "star"
        if not bam_root.is_dir():
            raise SystemExit(f"Missing BAM root for {dataset}: {bam_root}")
        discovered = sorted(bam_root.glob(f"*/*{SUFFIX}"), key=lambda p: sample_from_bam(p))
        if len(discovered) != expected_count:
            raise SystemExit(
                f"{dataset}: expected {expected_count} BAMs, found {len(discovered)}"
            )
        by_sample = {sample_from_bam(path): path for path in discovered}
        if len(by_sample) != len(discovered):
            raise SystemExit(f"{dataset}: discovered duplicate BAM sample IDs")
        current_counts = patient_root / "results" / dataset / "counts/featurecounts/gene_counts.txt"
        if not current_counts.is_file():
            raise SystemExit(f"Missing current matrix used to freeze BAM order: {current_counts}")
        with current_counts.open() as handle:
            handle.readline()
            header = handle.readline().rstrip("\n").split("\t")
        ordered_samples = [sample_from_bam(Path(value)) for value in header[6:]]
        if len(ordered_samples) != expected_count or set(ordered_samples) != set(by_sample):
            raise SystemExit(f"{dataset}: current matrix header and discovered BAM set differ")
        bams = [by_sample[sample] for sample in ordered_samples]
        seen: set[str] = set()
        for bam in bams:
            sample = sample_from_bam(bam)
            if sample in seen:
                raise SystemExit(f"{dataset}: duplicate sample ID {sample}")
            seen.add(sample)
            stat = bam.stat()
            rows.append(
                {
                    "dataset": dataset,
                    "sample_id": sample,
                    "bam_path": str(bam.resolve(strict=True)),
                    "layout": "paired",
                    "strandedness": "2",
                    "expected_status": "included",
                    "exclusion_reason": "",
                    "bam_size": stat.st_size,
                    "bam_mtime": str(int(stat.st_mtime)),
                    "bam_sha256": "",
                }
            )

    for dataset, (sample, reason) in DOCUMENTED_MISSING.items():
        rows.append(
            {
                "dataset": dataset,
                "sample_id": sample,
                "bam_path": "",
                "layout": "paired",
                "strandedness": "2",
                "expected_status": "documented_unavailable",
                "exclusion_reason": reason,
                "bam_size": "",
                "bam_mtime": "",
                "bam_sha256": "",
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise SystemExit(f"Refusing to overwrite manifest: {args.output}")
    tmp = args.output.with_suffix(args.output.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, args.output)
    print(f"Wrote {args.output} with {len(rows)} rows ({sum(EXPECTED.values())} included)")


if __name__ == "__main__":
    main()
