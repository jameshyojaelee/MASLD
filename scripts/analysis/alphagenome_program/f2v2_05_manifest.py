#!/usr/bin/env python3
"""f2-haplotype-v2 step 5: MANIFEST.tsv of every input this package read, with sha256 for the small ones."""

from __future__ import annotations

import hashlib
import os
import pathlib

OUT = pathlib.Path(os.environ["F2V2_OUT_ROOT"])
PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
V1 = PROJECT / "GWAS/finemapping/results/alphagenome_program/f2-haplotype-20260914T230433Z"
BIG = 200 * 1024 * 1024

INPUTS = [
    # the superseded package, read only
    V1 / "raw/scored_rows.jsonl",
    V1 / "tables/f3_observed_effects_n2500.tsv",
    V1 / "tables/f3_null_effects_n2500.tsv",
    V1 / "tables/f3_null_effects_n1000.tsv",
    V1 / "tables/f3_summary_n2500.json",
    V1 / "tables/f3_summary_n1000.json",
    V1 / "tables/f2_gnmt_verdict_n2500.json",
    V1 / "tables/null_panel/chr6.tsv.gz",
    V1 / "prespec/f2_00_prespec.json",
    V1 / "prespec/f2_00_prespec.sha256",
    V1 / "RESULTS.md",
    # the recipe and the API wrapper, read only
    PROJECT / "scripts/analysis/alphagenome_program/f2_recipe.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2_02_score.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2_03_analyse.py",
    PROJECT / "scripts/analysis/alphagenome_atlas/atlas_query.py",
    PROJECT / "scripts/analysis/alphagenome_atlas/lib_atlas.py",
    # this package's own code and prespecification
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_00_prespec.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_00_prespec.json",
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_00_amendment_01.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_00_amendment_01.json",
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_01_draw_genespan_null.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_02_score_genespan_null.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_02b_loader_equivalence.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_02b_probe.sbatch",
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_03_recompute.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_07_mark_superseded.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_04_write_results.py",
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_02_score.sbatch",
    PROJECT / "scripts/analysis/alphagenome_program/f2v2_06_final.sbatch",
    PROJECT / "scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md",
    # reference
    PROJECT / "Analysis/MASLD_Model_Benchmark/executions/reference-build-21062075/gencode.v49.primary_analysis.annotation.gtf.gz",
    pathlib.Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"),
    PROJECT / "GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/tables/haplotype_pairs.tsv",
]
OUTPUTS = [
    OUT / "tables/gnmt_genespan_null_pairs.tsv",
    OUT / "tables/gnmt_genespan_null_design.json",
    OUT / "raw/scored_genespan_null.jsonl",
    OUT / "tables/f2v2_gnmt_verdict.json",
    OUT / "tables/f2v2_both_active_stratum.json",
    OUT / "tables/f2v2_readout_geometry.json",
    OUT / "tables/f2v2_accounting.json",
    OUT / "tables/f2v2_summary.json",
    OUT / "tables/f2v2_loader_equivalence.json",
    OUT / "tables/f2v2_genespan_null_effects.tsv",
    OUT / "prespec/f2v2_00_prespec.json",
    OUT / "prespec/f2v2_00_amendment_01.json",
    OUT / "RESULTS.md",
]


def sha256(path: pathlib.Path) -> str:
    d = hashlib.sha256()
    with open(path, "rb") as h:
        for b in iter(lambda: h.read(1 << 20), b""):
            d.update(b)
    return d.hexdigest()


def main() -> None:
    n = 0
    with open(OUT / "MANIFEST.tsv", "w") as out:
        out.write("path\tbytes\tsha256\trole\n")
        for role, paths in (("input", INPUTS), ("output", OUTPUTS)):
            for p in paths:
                n += 1
                if not p.exists():
                    out.write(f"{p}\tMISSING\t\t{role}\n")
                    continue
                size = p.stat().st_size
                out.write(f"{p}\t{size}\t{sha256(p) if size <= BIG else 'not_hashed_over_200MB'}\t{role}\n")
    print(f"MANIFEST.tsv written: {n} entries")


if __name__ == "__main__":
    main()
