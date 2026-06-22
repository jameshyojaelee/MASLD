"""
guides_config.py — shared paths + constants for Cas13 sgRNA guide generation.

Single source of truth for: the upstream (sfriedman) genome-wide guide pools
(READ-ONLY), the cached query index, the library roster, the control rosters,
and the guide-design parameters. Imported by build_guide_index.py,
guide_selection.py, pull_guides.py, and build_library_guides.py.
"""
from __future__ import annotations

import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Project root (override with MASLD_PROJECT_ROOT)
# ---------------------------------------------------------------------------
PROJECT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

CAS13LIB = PROJECT / "Cas13_Library_Design"

# ---------------------------------------------------------------------------
# Upstream guide pools (sfriedman / track-cas13) — READ ONLY, never written to.
# Each release file is a post-basic-criteria pool (no homopolymers, target
# unique per isoform, TIGER>=0.75 OR Cas13Design>=0.75, off-target filtered).
# Covers protein-coding + lncRNA only (miRNA biotype is excluded upstream).
# ---------------------------------------------------------------------------
SF_BASE = Path(
    "/gpfs/commons/groups/sanjana_lab/sfriedman/track-cas13/outputs"
)

RELEASES = {
    "vM38": {
        "genome": SF_BASE / "mouse_GENCODEvM38_v1/sgrna/gencode.vM38.5.genome.csv",
        "posthoc": SF_BASE / "mouse_GENCODEvM38_v1/summary/gencode.vM38.posthoc.guide_counts.tsv",
    },
    "vM37": {
        "genome": SF_BASE / "mouse_GENCODEvM37_v1/sgrna/gencode.vM37.5.genome.csv",
        "posthoc": SF_BASE / "mouse_GENCODEvM37_v1/summary/gencode.vM37.posthoc.guide_counts.tsv",
    },
}
DEFAULT_RELEASE = "vM38"  # GENCODE vM38 / GRCm39 is the project-canonical mouse reference

# ---------------------------------------------------------------------------
# Outputs / cache (under Cas13_Library_Design/data/guides/)
# ---------------------------------------------------------------------------
GUIDES_DIR = CAS13LIB / "data" / "guides"
CACHE_DIR = GUIDES_DIR / "cache"          # git-ignored derived parquet index


def index_path(release: str) -> Path:
    """Cached parquet index for a release (built by build_guide_index.py)."""
    return CACHE_DIR / f"guides_{release}.parquet"


# ---------------------------------------------------------------------------
# Library roster + annotations
# ---------------------------------------------------------------------------
LIBRARY_CSV = CAS13LIB / "data" / "cas13_library_v3.0.csv"   # v8 target list
LIBRARY_VERSION = "v8"

# Library annotation columns merged onto the guide table (keyed by gene_id_mouse).
# v7: positive controls are folded into the target roster, so the guide table
# carries is_positive_control / pos_control_direction instead of a separate pull.
# v8: adds the mouse-hepatocyte gate axis (mouse_hep_substrate), the mouse_untestable
# flag (hep-failing genes kept via control/COLOC exemption), and library_arm
# (primary vs exploratory_lncrna).
LIBRARY_ANNOT_COLS = [
    "tier", "library_arm", "gene_symbol_human", "has_human_evidence", "has_coloc",
    "hep_substrate", "mouse_hep_substrate", "mouse_untestable", "sc_disease_celltype",
    "direction_conflict", "ortholog_ambiguous", "is_positive_control",
    "pos_control_direction", "control_benchmark_eligible",
]

# Biotypes we can build guides for from the upstream pool (miRNA deferred).
TARGET_BIOTYPES = ("protein_coding", "lncRNA")

# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------
# Positive controls: 65 human-symbol genes with a steatosis-direction label.
POS_CONTROL_CSV = PROJECT / "results" / "library" / "positive_control.csv"
POS_CONTROL_SYMBOL_COL = "Gene symbol"
POS_CONTROL_DIR_COL = "Steatosis_Change_upon_KD"   # "Decrease" / "Increase"
# Roster uses colloquial names that differ from the HGNC symbol the ortholog
# table is keyed on. Aliases are applied only for the control lookup.
POS_CONTROL_ALIASES = {"SCD1": "SCD"}   # human SCD (stearoyl-CoA desaturase) -> mouse Scd1

# Essential genes (Cas13-activity QC; must deplete). PROTOCOL.md §4 names
# RPA1 / RPL9 / POLR2A; we add a small curated pan-essential panel (Hart 2015
# core-essential, ubiquitous translation/transcription/replication machinery)
# to reach the documented 10-15. Human symbols -> mapped to mouse at build time.
ESSENTIAL_GENES_HUMAN = [
    "RPA1", "RPL9", "POLR2A",            # PROTOCOL.md-named
    "RPL3", "RPS3", "RPL5", "RPS5",      # ribosomal core-essential
    "POLR2B", "POLR2L",                  # RNA pol II core
    "EEF2", "PCNA", "RRM1",              # translation elongation / replication
]

# ---------------------------------------------------------------------------
# Guide-design parameters
# ---------------------------------------------------------------------------
N_GUIDES_DEFAULT = 10            # flat guides-per-gene (PI-pending; re-runnable)
GUIDE_LEN = 23                   # RfxCas13d spacer length in this pool
TIGER_MIN = 0.75                 # basic-criteria guard (re-asserted)
CAS13_MIN = 0.75                 # basic-criteria guard (re-asserted: TIGER>=.75 OR cas13>=.75)
HOMOPOLYMERS = ("TTTT", "AAAAA", "CCCCC", "GGGGG")  # disallowed runs
ORTHOLOG_MIN_TIER = "M"          # ortholog_lookup min_tier for control mapping

# Columns kept in the cached index (slim subset of the 25-col genome.csv)
INDEX_COLS = [
    "gene_id", "symbol", "guide_seq", "target_seq", "biotype", "region",
    "tiger_score", "cas13_score", "combined_score", "n_target", "n_appears",
    "tx_id_set", "tx_id_pos", "any_indel", "mismatch", "offtarget_pseudogene_only",
]
