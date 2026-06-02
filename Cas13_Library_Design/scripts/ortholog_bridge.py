"""Mouse↔Human ortholog bridge for the Cas13 MASH in vivo library design.

Problem the previous bridge had:
  * `Analysis/Cas13_Library_Design/reviews/inspect_perdiet_de.py` relied solely on
    the atlas column `mouse_ortholog`.  That column is almost entirely
    protein-coding (1 of 12,671 lncRNAs has a value), so lncRNA cross-species
    analyses came out at 1 / 0 instead of any biologically meaningful number.
  * Mouse DE files carry **versioned** ENSMUSG (`ENSMUSG00000028059.17`).  Any
    merge against the atlas without `str.split('.').str[0]` silently misses
    every gene.
  * Atlas `ensembl_id` is also versioned (`ENSG00000…N`) which fails against
    the unversioned human side of every ortholog table.

This module wraps four independent resources and treats them as fallbacks
in priority order:
    1.  data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz   (canonical
        per project — Ensembl biomaRt symbols + ENSMUSG/ENSG, with
        `ortholog_type` annotation).  ALL 25,439 pairs (one2one + one2many +
        many2many) are kept — paralog fan-in/fan-out is propagated, not
        filtered.
    2.  archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz
        (older Ensembl-ID-only version of the same biomaRt extract).
    3.  data/external/orthologs/layers/L1.5_phasej_synteny.tsv  (Phase J
        ncRNA pipeline output, RNA-seq/55_ncrna_conservation.R — flanking-PCG
        synteny across all 12,671 human lncRNAs.  Recovers 4,674 human↔mouse
        lncRNA syntenic pairs = 40,191 exploded (mouse, human) edges, of which
        14,080 unique mouse lncRNAs.  This is what fills the biomaRt lncRNA
        gap of 1/12,671 → ~4,674).
    4.  Atlas `mouse_ortholog` column (Ensembl-ID, last-resort).

The union of (1)+(2)+(3)+(4) is what the user gets.  Phase J synteny dominates
the ncRNA tier; biomaRt dominates the protein-coding tier; atlas adds a few
back-fills.

Phase 0 quick-win audit (2026-05-24): 243/545 mouse-only-UP lncRNAs in the
current Cas13 candidate set (Analysis/Cas13_Library_Design/data/
candidates_lncrna_independent.csv) recover a human counterpart via the Phase
J layer (44.6%) — matches the audit prediction exactly.

Public API:
    load_mouse_de(diet)           → DataFrame indexed by `gene_id_base`
    load_human_atlas(cols=None)   → atlas subset, ensembl_id_base stripped
    build_ortholog_table()        → mouse↔human bridge (long form;
                                    one row per mouse↔human pair)
    bridge_mouse_to_human(de, atlas)   → joined DE×atlas table
    cross_validate_known_examples()    → 11-pair sanity test (see CLI block)
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PERDIET = ROOT / "RNA-seq/Mouse/Unified_Integration/results/per_diet"
ATLAS = ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
ORTHO_SYM = ROOT / "data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz"
ORTHO_ID = ROOT / "archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz"
ORTHO_PHASEJ = ROOT / "data/external/orthologs/layers/L1.5_phasej_synteny.tsv"
ORTHO_BLAST = ROOT / "data/external/orthologs/layers/L4_blast.tsv"
MASTER_TABLE = ROOT / "data/external/orthologs/master_ortholog_table.tsv.gz"

# Auto-discover diet models from M02 output files in per_diet/
DIETS = sorted(
    f.stem.replace("_de_results", "")
    for f in PERDIET.glob("*_de_results.csv")
) if PERDIET.is_dir() else []

# Atlas columns the library design pipeline needs
ATLAS_COLS_DEFAULT = [
    "human_symbol",
    "ensembl_id",
    "gene_biotype",
    "mouse_ortholog",
    "mouse_ensembl",
    "dream_logFC",
    "dream_padj",
    "coloc_susie_best_pp4",
    "coloc_best_susie_pp4_polyfun",
    "coloc_susie_best_gwas",
    "coloc_best_susie_gwas_polyfun",
]


# ----------------------------- helpers ------------------------------------- #

def strip_version(s: pd.Series) -> pd.Series:
    """Strip `.N` Ensembl version suffix.  Idempotent on already-base IDs."""
    return s.astype(str).str.split(".").str[0]


def _coerce_str(s):
    return s.astype("string").where(s.notna(), pd.NA)


# ----------------------------- loaders ------------------------------------- #

def load_mouse_de(diet: str) -> pd.DataFrame:
    """Load `per_diet/{diet}_de_results.csv`.

    Returns columns:
        gene_id_base (unversioned ENSMUSG; index)
        logFC, AveExpr, t, P.Value, adj.P.Val, B
    """
    de_file = PERDIET / f"{diet}_de_results.csv"
    if not de_file.exists():
        available = DIETS if DIETS else ["(no per_diet files found)"]
        raise ValueError(
            f"No DE results for diet '{diet}' at {de_file}. "
            f"Available diets: {available}"
        )
    df = pd.read_csv(de_file)
    df["gene_id_base"] = strip_version(df["gene"])
    df = df.drop_duplicates("gene_id_base", keep="first")
    return df.set_index("gene_id_base")


def load_human_atlas(cols=None) -> pd.DataFrame:
    """Load multi-evidence atlas with versioned + base Ensembl IDs."""
    use = cols or ATLAS_COLS_DEFAULT
    use = [c for c in use if c]  # tolerate None
    atlas = pd.read_csv(ATLAS, usecols=use, low_memory=False)
    atlas["ensembl_id_base"] = strip_version(atlas["ensembl_id"])
    return atlas


def _load_symbol_ortho() -> pd.DataFrame:
    """biomaRt symbols + IDs.  Columns: human_gene_symbol, human_gene_ensembl,
    mouse_gene_symbol, mouse_gene_ensembl, ortholog_type, conservation_score."""
    df = pd.read_csv(ORTHO_SYM, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    return df


def _load_id_ortho() -> pd.DataFrame:
    """Ensembl-ID-only biomaRt dump.  Columns:
    mouse_ensembl_gene_id, human_ensembl_gene_id, orthology_type."""
    df = pd.read_csv(ORTHO_ID, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    return df


def _load_phasej_ortho() -> pd.DataFrame:
    """Phase J ncRNA synteny output (L1.5).  Common-schema TSV emitted by
    `phasej_to_orthotable.py` — 40,191 mouse↔human lncRNA pairs (14,080 unique
    mouse, 4,674 unique human).  Source is RNA-seq/55_ncrna_conservation.R."""
    df = pd.read_csv(ORTHO_PHASEJ, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    return df


def _load_blast_ortho() -> pd.DataFrame:
    """BLAST sequence-level orthology (L4).  Emitted by blast_to_orthotable.py
    after blastn dc-megablast filtered at pident>=50, length>=50, evalue<=1e-5.
    14,161 gene-level mouse↔human lncRNA pairs (6,094 unique mouse, 5,437
    unique human).  Phase 3a-2 audit verified BLAST recovers Phase J recall
    gaps for canonical lncRNAs (MEG3, HOTAIR, XIST, H19, TERC, PVT1, NORAD,
    HOTTIP, KCNQ1OT1) — 88% spot-check recovery."""
    df = pd.read_csv(ORTHO_BLAST, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    # Numeric promotion for filter columns
    for c in ["blast_pident", "blast_length", "blast_evalue", "blast_bitscore"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _build_from_master(min_tier: str, include_atlas_fallback: bool) -> pd.DataFrame:
    """Build ortholog table from the master table, filtering by confidence tier.

    This is the preferred path when the master table exists and min_tier filtering
    is requested.  The master table already has a curated confidence_tier column
    (H / M / L) that supersedes the layer-by-layer heuristics.

    Returns the same schema as the layer-by-layer fallback: mouse_ensembl,
    human_ensembl, mouse_symbol, human_symbol, ortholog_type, source.
    """
    tier_rank = {"H": 0, "M": 1, "L": 2}
    threshold = tier_rank.get(min_tier, 2)

    master = pd.read_csv(MASTER_TABLE, sep="\t", low_memory=False, dtype=str)
    master["_tier_rank"] = master["confidence_tier"].map(tier_rank)
    master = master[master["_tier_rank"] <= threshold].drop(columns="_tier_rank")

    # Map master-table provenance → ortholog_type + source for downstream compat
    def _classify_row(row):
        prov = str(row.get("provenance_sources", ""))
        if row.get("tier_H_biomart") == "1":
            otype = row.get("biomart_ortholog_type", "unknown")
            # Normalise biomaRt ortholog type names
            if "one2one" in str(otype):
                otype = "ortholog_one2one"
            elif "one2many" in str(otype) or "many2one" in str(otype):
                otype = "ortholog_one2many"
            elif "many2many" in str(otype):
                otype = "ortholog_many2many"
            return otype, "biomaRt_symbols"
        if row.get("tier_H_mirbase") == "1" or row.get("tier_H_mirgenedb") == "1":
            return "mirna_ortholog", "mirbase_mirgenedb"
        if row.get("tier_H_toga") == "1":
            toga_class = str(row.get("toga_orthology_class", ""))
            return f"toga_{toga_class}" if toga_class else "toga", "toga"
        if row.get("tier_M_blast_rbh") == "1":
            return "blastn_rbh", "blast_rbh"
        if row.get("tier_M_blast") == "1":
            return "blastn_sequence", "blast_seq"
        if row.get("tier_M_phasej_synteny") == "1":
            return "synteny_ncrna", "phasej_synteny"
        # Catch-all for any remaining M-tier evidence
        if "lncbook" in prov or "noncode" in prov or "hezroni" in prov:
            return "database_ncrna", prov.split(",")[0] if prov else "unknown"
        return "unknown", prov.split(",")[0] if prov else "unknown"

    classified = master.apply(_classify_row, axis=1, result_type="expand")
    classified.columns = ["ortholog_type", "source"]

    out = pd.DataFrame({
        "mouse_ensembl": strip_version(master["mouse_ensembl"]),
        "human_ensembl": strip_version(master["human_ensembl"]),
        "mouse_symbol": master["mouse_symbol"].where(master["mouse_symbol"].notna(), pd.NA),
        "human_symbol": master["human_symbol"].where(master["human_symbol"].notna(), pd.NA),
        "ortholog_type": classified["ortholog_type"],
        "source": classified["source"],
    })
    out = out.dropna(subset=["mouse_ensembl", "human_ensembl"])

    # Atlas backfill (same as layer-by-layer path)
    if include_atlas_fallback:
        atlas = load_human_atlas(["human_symbol", "ensembl_id", "mouse_ortholog"])
        a = atlas.rename(columns={"mouse_ortholog": "mouse_ensembl"})
        a = a.dropna(subset=["mouse_ensembl"])
        a["mouse_ensembl"] = strip_version(a["mouse_ensembl"])
        a["human_ensembl"] = strip_version(a["ensembl_id"])
        a["mouse_symbol"] = pd.NA
        a["ortholog_type"] = "atlas"
        a["source"] = "atlas_mouse_ortholog"
        a = a[["mouse_ensembl", "human_ensembl", "mouse_symbol",
               "human_symbol", "ortholog_type", "source"]]
        out = pd.concat([out, a], ignore_index=True)

    # Dedup: prefer higher-confidence source (biomaRt > mirbase > toga > blast >
    # phasej > atlas); within tied sources keep first row.
    source_rank = {
        "biomaRt_symbols": 0, "biomaRt_id": 0,
        "mirbase_mirgenedb": 1,
        "toga": 2,
        "blast_rbh": 3, "blast_seq": 4,
        "phasej_synteny": 5,
        "atlas_mouse_ortholog": 9,
    }
    out["_rank"] = out["source"].map(source_rank).fillna(8)
    out = out.sort_values("_rank").drop_duplicates(
        ["mouse_ensembl", "human_ensembl"], keep="first"
    ).drop(columns="_rank")
    return out.reset_index(drop=True)


def build_ortholog_table(
    include_atlas_fallback: bool = True,
    min_tier: str = "M",
) -> pd.DataFrame:
    """Long-form mouse↔human bridge.

    Parameters
    ----------
    include_atlas_fallback : bool
        Whether to include the atlas `mouse_ortholog` column as last-resort.
    min_tier : str, one of "H", "M", "L", or None
        Minimum confidence tier to include.  "H" = biomaRt / TOGA-one2one /
        miRNA databases only.  "M" (default) = also includes BLAST-RBH,
        Phase J + BLAST, and other M-tier evidence.  "L" = includes all pairs
        (Phase J synteny-only positional candidates).  None = legacy behaviour
        (layer-by-layer union without tier filtering).

    Returns columns:
        mouse_ensembl  (unversioned ENSMUSG)
        human_ensembl  (unversioned ENSG)
        mouse_symbol   (optional)
        human_symbol   (optional)
        ortholog_type  (one2one / one2many / many2many / unknown)
        source         ("biomaRt_symbols", "biomaRt_id", "atlas_mouse_ortholog", ...)

    One row per (mouse, human) pair.  Both sides may appear in multiple rows
    when paralog fans-in/fans-out.  Downstream code should call
    `.drop_duplicates()` on the columns it cares about.
    """
    # ---- New path: use curated master table with tier filtering ----
    if min_tier in ("H", "M", "L") and MASTER_TABLE.exists():
        return _build_from_master(min_tier, include_atlas_fallback)

    # ---- Fallback: layer-by-layer assembly (original logic, no tier filter) ----
    pieces = []

    # 1. canonical: symbols file
    s = _load_symbol_ortho()
    s = s.rename(
        columns={
            "mouse_gene_ensembl": "mouse_ensembl",
            "human_gene_ensembl": "human_ensembl",
            "mouse_gene_symbol": "mouse_symbol",
            "human_gene_symbol": "human_symbol",
            "ortholog_type": "ortholog_type",
        }
    )
    s["mouse_ensembl"] = strip_version(s["mouse_ensembl"])
    s["human_ensembl"] = strip_version(s["human_ensembl"])
    s["source"] = "biomaRt_symbols"
    pieces.append(s[["mouse_ensembl", "human_ensembl", "mouse_symbol",
                     "human_symbol", "ortholog_type", "source"]])

    # 2. ID-only file
    e = _load_id_ortho()
    e = e.rename(
        columns={
            "mouse_ensembl_gene_id": "mouse_ensembl",
            "human_ensembl_gene_id": "human_ensembl",
            "orthology_type": "ortholog_type",
        }
    )
    e["mouse_ensembl"] = strip_version(e["mouse_ensembl"])
    e["human_ensembl"] = strip_version(e["human_ensembl"])
    e["mouse_symbol"] = pd.NA
    e["human_symbol"] = pd.NA
    e["source"] = "biomaRt_id"
    pieces.append(e[["mouse_ensembl", "human_ensembl", "mouse_symbol",
                     "human_symbol", "ortholog_type", "source"]])

    # 3. Phase J synteny (L1.5) — fills the ncRNA gap that biomaRt misses
    if ORTHO_PHASEJ.exists():
        pj = _load_phasej_ortho()
        pj["mouse_ensembl"] = strip_version(pj["mouse_ensembl"])
        pj["human_ensembl"] = strip_version(pj["human_ensembl"])
        pj["ortholog_type"] = "synteny_ncrna"
        pj["source"] = "phasej_synteny"
        pieces.append(pj[["mouse_ensembl", "human_ensembl", "mouse_symbol",
                          "human_symbol", "ortholog_type", "source"]])

    # 4. BLAST L4 — sequence-validated lncRNA pairs (Phase 3 of 2026-05-24)
    if ORTHO_BLAST.exists():
        bl = _load_blast_ortho()
        bl["mouse_ensembl"] = strip_version(bl["mouse_ensembl"])
        bl["human_ensembl"] = strip_version(bl["human_ensembl"])
        bl["ortholog_type"] = "blastn_sequence"
        bl["source"] = "blast_seq"
        pieces.append(bl[["mouse_ensembl", "human_ensembl", "mouse_symbol",
                          "human_symbol", "ortholog_type", "source"]])

    # 5. atlas column (last-resort)
    if include_atlas_fallback:
        atlas = load_human_atlas(["human_symbol", "ensembl_id", "mouse_ortholog"])
        a = atlas.rename(
            columns={
                "ensembl_id_base": "human_ensembl",
                "mouse_ortholog": "mouse_ensembl",
            }
        )
        a = a.dropna(subset=["mouse_ensembl"])
        a["mouse_ensembl"] = strip_version(a["mouse_ensembl"])
        a["human_ensembl"] = strip_version(a["ensembl_id"])
        a["mouse_symbol"] = pd.NA
        a["ortholog_type"] = "atlas"
        a["source"] = "atlas_mouse_ortholog"
        pieces.append(a[["mouse_ensembl", "human_ensembl", "mouse_symbol",
                         "human_symbol", "ortholog_type", "source"]])

    out = pd.concat(pieces, ignore_index=True)
    out = out.dropna(subset=["mouse_ensembl", "human_ensembl"])
    # Dedup on the (mouse, human) pair, prefer biomaRt_symbols > biomaRt_id >
    # blast_seq > phasej_synteny > atlas_mouse_ortholog
    # (BLAST ranked above Phase J because it provides direct sequence-level
    # confirmation; Phase J alone is positional and noisy)
    rank = {
        "biomaRt_symbols": 0,
        "biomaRt_id": 1,
        "blast_seq": 2,
        "phasej_synteny": 3,
        "atlas_mouse_ortholog": 4,
    }
    out["_rank"] = out["source"].map(rank).fillna(99)
    out = out.sort_values("_rank").drop_duplicates(
        ["mouse_ensembl", "human_ensembl"], keep="first"
    ).drop(columns="_rank")
    return out.reset_index(drop=True)


# ----------------------------- bridge -------------------------------------- #

def bridge_mouse_to_human(
    mouse_de: pd.DataFrame,
    atlas: pd.DataFrame,
    ortho: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Merge mouse DE with atlas via the ortholog bridge.

    Inputs:
        mouse_de : DataFrame indexed by gene_id_base (output of load_mouse_de).
        atlas    : DataFrame from load_human_atlas (must contain ensembl_id_base).
        ortho    : long ortholog table (default: build_ortholog_table()).

    Output:
        One row per (mouse_ensembl, human_ensembl) successfully matched on
        both sides.  All atlas columns + mouse DE columns are carried.
    """
    if ortho is None:
        ortho = build_ortholog_table()
    de = mouse_de.reset_index().rename(columns={"gene_id_base": "mouse_ensembl"})
    de = de.add_prefix("mouse_").rename(columns={"mouse_mouse_ensembl": "mouse_ensembl"})

    a = atlas.copy()
    a = a.rename(columns={"ensembl_id_base": "human_ensembl"})

    j = (
        ortho.merge(de, on="mouse_ensembl", how="inner")
        .merge(a, on="human_ensembl", how="inner", suffixes=("", "_atlas"))
    )
    return j


# ----------------------------- sanity check -------------------------------- #

KNOWN_PAIRS = [
    ("Col1a1", "COL1A1"),
    ("Mmp12", "MMP12"),
    ("Lcn2", "LCN2"),
    ("Timp1", "TIMP1"),
    ("Tgfb1", "TGFB1"),
    ("Acta2", "ACTA2"),
    ("Cyp7a1", "CYP7A1"),
    ("Pcsk9", "PCSK9"),
    ("Hsd17b13", "HSD17B13"),
    ("Scd1", "SCD"),   # known paralog quirk: biomaRt assigns SCD↔Scd2/Scd4
    ("Fasn", "FASN"),
]


def cross_validate_known_examples() -> dict:
    """Round-trip 11 known fibrotic / metabolic genes through the bridge.

    Returns dict keyed by mouse symbol → dict with keys
        {expected_human, found_human, mouse_ensembl, human_ensembl, ortholog_type, status}
    `status` is one of:
        - "MATCH"           : human symbol matches expectation
        - "PARALOG"         : ortholog table maps to a different human symbol (Ensembl one2many)
        - "MISSING"         : no row in the bridge for this mouse symbol
    """
    sym = _load_symbol_ortho()
    out = {}
    for m_sym, h_sym_expected in KNOWN_PAIRS:
        rows = sym[sym["mouse_gene_symbol"] == m_sym]
        if rows.empty:
            out[m_sym] = dict(
                expected_human=h_sym_expected,
                found_human=None, mouse_ensembl=None,
                human_ensembl=None, ortholog_type=None,
                status="MISSING",
            )
            continue
        # take first row (one2one canonical)
        r = rows.iloc[0]
        found = r["human_gene_symbol"]
        status = "MATCH" if found == h_sym_expected else "PARALOG"
        out[m_sym] = dict(
            expected_human=h_sym_expected,
            found_human=found,
            mouse_ensembl=r["mouse_gene_ensembl"],
            human_ensembl=r["human_gene_ensembl"],
            ortholog_type=r["ortholog_type"],
            status=status,
        )
    return out


# ----------------------------- CLI smoke ----------------------------------- #

if __name__ == "__main__":
    print("=== ortholog_bridge.py smoke test ===")
    res = cross_validate_known_examples()
    n_match = sum(1 for v in res.values() if v["status"] == "MATCH")
    n_paralog = sum(1 for v in res.values() if v["status"] == "PARALOG")
    n_missing = sum(1 for v in res.values() if v["status"] == "MISSING")
    print(f"  matches: {n_match}/{len(res)};  paralogs: {n_paralog};  missing: {n_missing}")
    for m, v in res.items():
        print(f"    {m} → expected {v['expected_human']!s:>10s} | found {v['found_human']!s:>10s} | {v['status']} | {v['ortholog_type']}")
    print()
    for tier in ("H", "M", "L", None):
        kwargs = {"min_tier": tier} if tier else {"min_tier": None}
        ortho = build_ortholog_table(**kwargs)
        label = f"min_tier={tier!r}"
        print(f"Bridge table ({label}): {len(ortho):,} (mouse,human) pairs.")
        print(f"  unique mouse_ensembl: {ortho['mouse_ensembl'].nunique():,}")
        print(f"  unique human_ensembl: {ortho['human_ensembl'].nunique():,}")
        print(f"  source counts:\n{ortho['source'].value_counts()}")
        print()
