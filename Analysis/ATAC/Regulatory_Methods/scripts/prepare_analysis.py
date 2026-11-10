#!/usr/bin/env python3
"""Prepare fixed gene coordinates and record the candidate's execution context."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess


def gene_records(gtf):
    seen = set()
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) != 9 or row[2] != "gene" or row[0] not in {
                *(f"chr{x}" for x in range(1, 23)), "chrX", "chrY"
            }:
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', row[8]))
            gene_id = attrs["gene_id"].split(".")[0]
            # PAR copies and other repeated stable IDs are not silently merged.
            if gene_id in seen:
                raise ValueError(f"Duplicate primary-assembly stable gene ID: {gene_id}")
            seen.add(gene_id)
            start0, end = int(row[3]) - 1, int(row[4])
            strand = row[6]
            if strand not in ("+", "-") or end <= start0:
                raise ValueError("Invalid gene coordinates")
            tss0 = start0 if strand == "+" else end - 1
            yield dict(gene_id=gene_id, gene_name=attrs.get("gene_name", gene_id),
                       chrom=row[0], start0=start0, end=end, tss0=tss0, strand=strand,
                       promoter_start0=max(0, tss0 - 2000), promoter_end=tss0 + 2001,
                       gene_type=attrs.get("gene_type", ""))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    genes = list(gene_records(config["gtf"]))
    with (args.output / "genes.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(genes[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(genes)
    tracked = subprocess.run(["git", "status", "--short"], cwd=config["project_root"],
                             capture_output=True, text=True, check=True)
    (args.output / "working_tree_at_start.txt").write_text(tracked.stdout)
    summary = dict(candidate_only=True, seed=config["seed"], n_genes=len(genes),
                   tss_definition="gene boundary, zero based; minus strand end minus one",
                   promoter_definition="[max(0,tss0-2000), tss0+2001)",
                   host=platform.node(), slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                   config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest(),
                   gtf=config["gtf"], scientific_adoption=False)
    (args.output / "preparation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
