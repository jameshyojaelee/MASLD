#!/usr/bin/env python3
"""F2/F3 step 4: MANIFEST.tsv of every input this package read, with sha256 for the small ones."""

from __future__ import annotations

import hashlib
import os
import pathlib

OUT = pathlib.Path(os.environ["F2_OUT_ROOT"])
PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BIG = 200 * 1024 * 1024

INPUTS = [
    PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/tables/haplotype_pairs.tsv",
    PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/tables/haplotype_null_effects.tsv",
    PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/tables/haplotype_summary.json",
    PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/50_haplotype_prespec.json",
    PROJECT / "scripts/analysis/alphagenome_atlas/50_haplotype_additivity.py",
    PROJECT / "scripts/analysis/alphagenome_atlas/atlas_query.py",
    PROJECT / "scripts/analysis/alphagenome_atlas/lib_atlas.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2_recipe.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2_01_build_null_panel.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2_01_liftover.R",
    PROJECT / "scripts/analysis/alphagenome_program/f2_02_score.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2_03_analyse.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2_00_prespec.json",
    PROJECT / "scripts/analysis/alphagenome_program/f2_00_amendment_01.json",
    PROJECT / "data/broadaway_eqtl/hg19ToHg38.over.chain",
    PROJECT / "Analysis/MASLD_Model_Benchmark/executions/reference-build-21062075/gencode.v49.primary_analysis.annotation.gtf.gz",
    pathlib.Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"),
    pathlib.Path("/gpfs/commons/home/jameslee/reference_genome/dbsnp/GCF_000001405.40.gz"),
]
INPUTS += sorted((PROJECT / "data/1kg_eur").glob("chr*_eur.bim"))
INPUTS += sorted((OUT / "tables" / "null_panel").glob("chr*.tsv.gz"))


def sha256(path: pathlib.Path) -> str:
    d = hashlib.sha256()
    with open(path, "rb") as h:
        for b in iter(lambda: h.read(1 << 20), b""):
            d.update(b)
    return d.hexdigest()


def main() -> None:
    with open(OUT / "MANIFEST.tsv", "w") as out:
        out.write("path\tbytes\tsha256\trole\n")
        for p in INPUTS:
            if not p.exists():
                out.write(f"{p}\tMISSING\t\tinput\n")
                continue
            n = p.stat().st_size
            digest = sha256(p) if n <= BIG else "not_hashed_over_200MB"
            out.write(f"{p}\t{n}\t{digest}\tinput\n")
    print(f"MANIFEST.tsv written: {len(INPUTS)} entries")


if __name__ == "__main__":
    main()
