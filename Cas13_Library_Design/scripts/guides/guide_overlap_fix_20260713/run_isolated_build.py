#!/usr/bin/env python3
"""Rebuild target/control guides with sequence-derived overlap filtering.

All writable paths are redirected beneath this script's directory. Production
scripts, indexes, reference files, and outputs are read-only inputs.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

RUN = Path(__file__).resolve().parent
OUT = RUN / "output"
ORIGINAL_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ORIGINAL_GUIDES = ORIGINAL_ROOT / "Cas13_Library_Design/data/guides"
GUIDE_SCRIPTS = ORIGINAL_ROOT / "Cas13_Library_Design/scripts/guides"
FASTA = ORIGINAL_GUIDES / "cache/nt_screen/transcriptome.fa"
GTF = Path("/gpfs/commons/home/jameslee/reference_genome/"
           "refdata-cellranger-GRCm39-vM38/genes/genes.gtf.gz")
ROSTER = RUN / "frozen_cas13_library.csv"

ROSTER_COLUMNS = [
    "gene_id_mouse", "gene_symbol_mouse", "biotype", "tier", "library_arm",
    "gene_symbol_human", "has_human_evidence", "has_coloc", "hep_substrate",
    "mouse_hep_substrate", "mouse_untestable", "sc_disease_celltype",
    "direction_conflict", "ortholog_ambiguous", "is_positive_control",
    "pos_control_direction", "control_benchmark_eligible",
]


def create_frozen_roster() -> None:
    source = ORIGINAL_GUIDES / "cas13_library_guides_vM38.csv"
    genes: dict[str, dict[str, str]] = {}
    with source.open(newline="") as handle:
        for row in csv.DictReader(handle):
            genes.setdefault(row["gene_id_mouse"],
                             {column: row.get(column, "") for column in ROSTER_COLUMNS})
    with ROSTER.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROSTER_COLUMNS)
        writer.writeheader()
        writer.writerows(genes.values())
    print(f"[isolated] frozen roster: {len(genes):,} genes -> {ROSTER}", flush=True)


def main() -> None:
    if OUT.exists() and any(path.is_file() for path in OUT.rglob("*")):
        raise SystemExit(f"Refusing to overwrite non-empty output directory: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "cache").mkdir(exist_ok=True)
    if not ROSTER.exists():
        create_frozen_roster()

    sys.path.insert(0, str(GUIDE_SCRIPTS))
    sys.path.insert(0, str(RUN))
    import guides_config as cfg

    # Writable locations are isolated. All large/reference inputs remain read-only.
    cfg.GUIDES_DIR = OUT
    cfg.CACHE_DIR = OUT / "cache"
    cfg.LIBRARY_CSV = ROSTER
    cfg.index_path = lambda release: ORIGINAL_GUIDES / "cache" / f"guides_{release}.parquet"
    cfg.LIVER_ESSENTIALITY_CACHE = (
        ORIGINAL_ROOT / "Cas13_Library_Design/data/depmap_liver_essentiality.csv"
    )

    import guide_selection as gs
    import pandas as pd
    from sequence_overlap_patch import SequenceReference, install

    roster = pd.read_csv(ROSTER)
    target_ids = set(roster["gene_id_mouse"].dropna().astype(str).str.split(".").str[0])
    essential_map = gs.human_to_mouse_ids(cfg.ESSENTIAL_GENES_HUMAN)
    target_ids.update(essential_map["gene_id_base"].dropna().astype(str).str.split(".").str[0])

    print(f"[isolated] loading sequence reference for {len(target_ids):,} genes", flush=True)
    reference = SequenceReference(FASTA, GTF, target_ids)
    print(f"[isolated] reference: {len(reference.transcripts):,} transcripts; "
          f"{sum(bool(t.sequence) for t in reference.transcripts.values()):,} with sequence",
          flush=True)
    install(gs, reference)

    import build_library_guides
    sys.argv = ["build_library_guides.py", "--release", "vM38", "--n", "4"]
    build_library_guides.main()

    metadata = {
        "mode": "isolated_sequence_derived_physical_nonoverlap",
        "production_scripts_modified": False,
        "guides_per_gene": 4,
        "reference_fasta": str(FASTA),
        "reference_gtf": str(GTF),
        "source_index": str(cfg.index_path("vM38")),
        "source_roster_reconstruction": str(ORIGINAL_GUIDES / "cas13_library_guides_vM38.csv"),
        "candidate_rows": reference.candidate_rows,
        "candidate_rows_located": reference.candidate_rows_located,
        "candidate_rows_unlocated": reference.candidate_rows_unlocated,
        "overlap_rule": "reject if exact target_seq-derived GTF genomic intervals share >=1 nt",
    }
    with (OUT / "SEQUENCE_OVERLAP_BUILD.json").open("x") as handle:
        json.dump(metadata, handle, indent=2)
        handle.write("\n")
    print(f"[isolated] sequence-location stats: {metadata}", flush=True)


if __name__ == "__main__":
    main()
