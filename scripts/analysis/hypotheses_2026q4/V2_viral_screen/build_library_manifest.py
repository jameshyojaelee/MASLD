#!/usr/bin/env python3
"""V2 target screen: build the array manifest of libraries to screen.

A library is included only when it appears both in the cohort's SRA run table
and as FASTQ on disk. The manifest is ordered by descending compressed FASTQ
size so that the longest jobs start first.

Usage:
  build_library_manifest.py --out <manifest.tsv>
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys

COHORTS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621",
           "GSE167523", "GSE174478", "GSE193066", "GSE240729", "PRJNA512027"]

POOLED_FIVE = {"GSE126848", "GSE130970", "GSE135251", "GSE162694",
               "GSE213621"}

ROOT = ("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
        "RNA-seq/Human/Patient_Cohorts")
RAW = os.path.join(ROOT, "data", "raw")

RUN_TABLES = [
    os.path.join(ROOT, "pipelines", "custom", "%s", "metadata",
                 "SraRunTable.csv"),
    os.path.join(ROOT, "archive", "old_data_metadata", "metadata",
                 "%s_SraRunTable.csv"),
]

ACC_RE = re.compile(r"^([A-Z]+[0-9]+)")


def run_table_runs(cohort: str) -> set[str]:
    for tpl in RUN_TABLES:
        path = tpl % cohort
        if os.path.exists(path):
            with open(path, newline="") as fh:
                rows = list(csv.DictReader(fh))
            for key in ("Run", "run", "Run_accession"):
                if rows and key in rows[0]:
                    return {r[key].strip() for r in rows if r[key].strip()}
    return set()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    records = []
    summary = []
    for cohort in COHORTS:
        table_runs = run_table_runs(cohort)
        fq_dir = os.path.join(RAW, cohort, "fastq")
        by_lib: dict[str, list[str]] = {}
        for path in sorted(glob.glob(os.path.join(fq_dir, "*.fastq.gz"))):
            m = ACC_RE.match(os.path.basename(path))
            if not m:
                continue
            by_lib.setdefault(m.group(1), []).append(path)

        on_disk = set(by_lib)
        included = sorted(on_disk & table_runs) if table_runs else sorted(on_disk)
        summary.append({
            "cohort": cohort, "in_run_table": len(table_runs),
            "on_disk": len(on_disk), "included": len(included),
            "disk_not_in_table": len(on_disk - table_runs) if table_runs else 0,
            "table_not_on_disk": len(table_runs - on_disk) if table_runs else 0,
        })

        for srr in included:
            files = sorted(by_lib[srr])
            fq1 = fq2 = ""
            for p in files:
                base = os.path.basename(p)
                if base == "%s_1.fastq.gz" % srr:
                    fq1 = p
                elif base == "%s_2.fastq.gz" % srr:
                    fq2 = p
                elif base == "%s.fastq.gz" % srr and not fq1:
                    fq1 = p
            if not fq1:
                fq1 = files[0]
            size = sum(os.path.getsize(p) for p in (fq1, fq2) if p)
            records.append({
                "cohort": cohort, "library": srr,
                "layout": "paired" if fq2 else "single",
                "fastq1": fq1, "fastq2": fq2, "bytes_gz": size,
                "in_pooled_five": int(cohort in POOLED_FIVE),
            })

    records.sort(key=lambda r: -r["bytes_gz"])
    cols = ["array_index", "cohort", "library", "layout", "in_pooled_five",
            "bytes_gz", "fastq1", "fastq2"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for i, r in enumerate(records, start=1):
            r["array_index"] = i
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    total = sum(r["bytes_gz"] for r in records)
    sys.stderr.write("manifest: %d libraries, %.2f TB compressed FASTQ\n"
                     % (len(records), total / 1e12))
    for s in summary:
        sys.stderr.write(
            "  %-12s run_table=%-4d on_disk=%-4d included=%-4d "
            "disk_only=%-3d table_only=%d\n"
            % (s["cohort"], s["in_run_table"], s["on_disk"], s["included"],
               s["disk_not_in_table"], s["table_not_on_disk"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
