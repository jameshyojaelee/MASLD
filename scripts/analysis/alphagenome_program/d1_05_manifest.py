#!/usr/bin/env python3
"""D1: input manifest (path, size, sha256 for small inputs) and the environment record."""
from __future__ import annotations

import hashlib
import pathlib
import sys

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SMALL = 200 * 1024 * 1024  # sha256 only below this; the BAMs are 6-40 GB

INPUTS = [
    "GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables/atac_targets_in_peaks.tsv",
    "GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables/allelic_sites.tsv",
    "GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables/gse244832_allelic_sites.tsv",
    "GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/73_allelic_prespec.json",
    "scripts/analysis/alphagenome_atlas/72_atac_allelic_counts.sh",
    "scripts/analysis/alphagenome_atlas/73_atac_allelic_analysis.py",
    "Analysis/ATAC/Human_External/results/coembed/umap_export_joint.tsv.gz",
    "Analysis/ATAC/Human_Multiome/results/snapatac2/cell_donor_condition_map.tsv.gz",
    "scripts/analysis/alphagenome_program/d1_prespec.json",
    "scripts/analysis/alphagenome_program/d1_01_lineage_barcodes.py",
    "scripts/analysis/alphagenome_program/d1_02_extract_and_count.sh",
    "scripts/analysis/alphagenome_program/d1_03_assemble_counts.py",
    "scripts/analysis/alphagenome_program/d1_04_mapping_bias.py",
]
EXTRA_GLOBS = ["Analysis/ATAC/Human_External/cellranger/Z*/outs/possorted_bam.bam",
               "Analysis/ATAC/Human_Multiome/results/alignment/D*.filtered.bam"]
REFS = ["/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa",
        "/gpfs/commons/home/jameslee/reference_genome/bowtie2/GRCh38/GRCh38.1.bt2"]


def row(p: pathlib.Path) -> str:
    size = p.stat().st_size
    if size <= SMALL:
        h = hashlib.sha256(p.read_bytes()).hexdigest()
    else:
        h = "not_hashed_over_200MB"
    return f"{p}\t{size}\t{h}"


def main() -> None:
    out = pathlib.Path(sys.argv[1]).resolve()
    lines = ["path\tsize_bytes\tsha256"]
    for rel in INPUTS:
        lines.append(row(PROJECT / rel))
    for g in EXTRA_GLOBS:
        for p in sorted(PROJECT.glob(g)):
            lines.append(row(p))
    for r in REFS:
        lines.append(row(pathlib.Path(r)))
    (out / "MANIFEST.tsv").write_text("\n".join(lines) + "\n")
    print(f"[d1_05] manifest rows {len(lines) - 1}")


if __name__ == "__main__":
    main()
