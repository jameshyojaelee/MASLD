#!/usr/bin/env python3
"""Audit a guide CSV using exact sequence-derived genomic target intervals."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from sequence_overlap_patch import SequenceReference, overlap_pairs

DEFAULT_FASTA = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
                     "Cas13_Library_Design/data/guides/cache/nt_screen/transcriptome.fa")
DEFAULT_GTF = Path("/gpfs/commons/home/jameslee/reference_genome/"
                   "refdata-cellranger-GRCm39-vM38/genes/genes.gtf.gz")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("guide_csv", type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    args = parser.parse_args()

    with args.guide_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows = [row for row in rows
            if row.get("control_class", "target") != "non_targeting"
            and (row.get("gene_id_mouse") or row.get("gene_id_base"))]
    for row in rows:
        row["gene_id_base"] = (row.get("gene_id_mouse") or row.get("gene_id_base")).split(".")[0]

    by_gene: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_gene[row["gene_id_base"]].append(row)
    reference = SequenceReference(args.fasta, args.gtf, by_gene)
    pairs: list[dict] = []
    for gene_rows in by_gene.values():
        pairs.extend(overlap_pairs(gene_rows, reference))

    fields = ["gene_id_mouse", "gene_symbol_mouse", "guide_id_1", "guide_id_2",
              "shared_genomic_nt", "overlap_intervals"]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(pairs)

    print(f"guides audited: {len(rows):,}")
    print(f"genes audited: {len(by_gene):,}")
    print(f"candidate rows located: {reference.candidate_rows_located:,}/"
          f"{reference.candidate_rows:,}")
    print(f"physical-overlap pairs: {len(pairs):,}")
    print(f"genes with overlap: {len({p['gene_id_mouse'] for p in pairs}):,}")
    print(f"report: {args.report}")


if __name__ == "__main__":
    main()
