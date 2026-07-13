#!/usr/bin/env python3
"""Preflight paths plus focused tests of prioritization and physical non-overlap."""
from __future__ import annotations

import sys
from pathlib import Path

RUN = Path(__file__).resolve().parent
ORIGINAL_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GUIDE_SCRIPTS = ORIGINAL_ROOT / "Cas13_Library_Design/scripts/guides"


def run_unit_tests() -> None:
    import test_sequence_overlap_patch as tests
    names = sorted(name for name in dir(tests) if name.startswith("test_"))
    for name in names:
        getattr(tests, name)()
        print(f"[preflight] PASS {name}")


def run_prioritization_tests() -> None:
    sys.path.insert(0, str(GUIDE_SCRIPTS))
    import guide_selection as gs
    from sequence_overlap_patch import sequence_far_enough

    old_far_enough = gs._far_enough
    gs._far_enough = sequence_far_enough
    try:
        def candidate(name, isoforms, score, start, region="CDS"):
            return {
                "name": name,
                "n_isoforms_targeted": isoforms,
                "combined_score": score,
                "position": start,
                "region": region,
                "_actual_intervals": (("chr1", start, start + 22),),
            }

        # High transcript coverage must outrank a higher-scoring low-coverage guide.
        rows = [candidate("constitutive_a", 5, 0.80, 100),
                candidate("constitutive_b", 5, 0.81, 200),
                candidate("isoform_specific_high_score", 1, 0.999, 300)]
        selected = gs._select_for_gene(rows, 2)
        assert {row["name"] for row in selected} == {"constitutive_a", "constitutive_b"}
        print("[preflight] PASS constitutive-first transcript-coverage priority")

        # CDS remains preferred over UTR even when UTR efficacy scores are higher.
        rows = [candidate("cds_a", 2, 0.80, 100),
                candidate("cds_b", 2, 0.81, 200),
                candidate("utr_high_score", 9, 0.999, 300, "3'UTR")]
        selected = gs._select_for_gene(rows, 2)
        assert {row["name"] for row in selected} == {"cds_a", "cds_b"}
        print("[preflight] PASS CDS-first priority")

        # Even a high-coverage/high-score candidate is rejected if it physically overlaps.
        rows = [candidate("first", 5, 0.99, 100),
                candidate("overlapping", 5, 0.98, 110),
                candidate("nonoverlapping", 5, 0.80, 200)]
        selected = gs._select_for_gene(rows, 2)
        assert len(selected) == 2
        assert not ({"first", "overlapping"} <= {row["name"] for row in selected})
        assert sequence_far_enough(selected[1], selected[:1], 23)
        print("[preflight] PASS physical non-overlap overrides score")
    finally:
        gs._far_enough = old_far_enough


def check_paths() -> None:
    paths = [
        Path("/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"),
        ORIGINAL_ROOT / "Cas13_Library_Design/data/guides/cache/guides_vM38.parquet",
        ORIGINAL_ROOT / "Cas13_Library_Design/data/guides/cache/nt_screen/transcriptome.fa",
        ORIGINAL_ROOT / "Cas13_Library_Design/data/depmap_liver_essentiality.csv",
        Path("/gpfs/commons/home/jameslee/reference_genome/"
             "refdata-cellranger-GRCm39-vM38/genes/genes.gtf.gz"),
        RUN / "frozen_cas13_library.csv",
    ]
    for path in paths:
        assert path.is_file(), f"missing required input: {path}"
        print(f"[preflight] READABLE {path}")
    output = RUN / "output"
    assert not output.exists() or not any(path.is_file() for path in output.rglob("*")), \
        f"output directory contains files; refusing overwrite: {output}"
    print(f"[preflight] OUTPUT READY {output}")


if __name__ == "__main__":
    check_paths()
    run_unit_tests()
    run_prioritization_tests()
    print("[preflight] READY TO SUBMIT")
