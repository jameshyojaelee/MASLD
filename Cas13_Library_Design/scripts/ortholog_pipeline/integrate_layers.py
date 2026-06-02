#!/usr/bin/env python
"""
integrate_layers.py — Build the master ortholog table by merging:
  L0   TOGA Zoonomia mouse-human (layers/L0_toga.tsv) — Kirilenko 2023 Science
  L1   biomaRt Ensembl (already in ortholog_bridge.py)
  L1.5 Phase J synteny (data/external/orthologs/layers/L1.5_phasej_synteny.tsv)
  L4   BLAST sequence (data/external/orthologs/layers/L4_blast.tsv)
  L4r  BLAST reciprocal best hit (layers/L4_blast_rbh.tsv) — when present
  L2   miRBase family-ID matching (layers/L2_mirbase.tsv)
  L5   MMseqs2 reciprocal best hit (layers/L5_mmseqs2_rbh.tsv) — Steinegger & Soeding 2017
  L6   Seekr k-mer functional orthology (layers/L6_seekr_kmer.tsv) — Kirk 2018
  L7   NONCODE v6 curated lncRNA conservation (layers/L7_noncode.tsv) — Zhao 2021
  L8   OrthoFinder gene-tree reconciliation (layers/L8_orthofinder.tsv) — Emms & Kelly 2019
  L9   ortho2align syntenic disambiguation (layers/L9_ortho2align.tsv) — Mylarshchikov 2022
  +    LncBook 2.0 conservation flags (curated_dbs/lncbook_conservation_v2.0.csv.gz)

  L5/L6/L7 added 2026-05-28 (were generated but never integrated). MMseqs2
  joins the sequence evidence class (with BLAST); NONCODE joins curated
  conservation (with LncBook); Seekr is a new independent functional-k-mer
  class. Promote to Tier M (MMseqs2 RBH, NONCODE, Seekr Pearson>0.8).

Tier semantics (v4, 2026-05-24 — adds TOGA L0 layer on top of v3):

  Tier H  Highest confidence — TOGA whole-genome Cactus + ML ortholog,
           OR Ensembl Compara strict 1:1, OR miRBase/MirGeneDB family,
           OR OrthoFinder gene-tree reconciliation (Emms & Kelly 2019).
           Sources allowed:
             - tier_H_toga == 1 (any class: one2one, one2many, many2one,
               many2many; one2many+ flagged in notes/provenance)
             - tier_H_biomart == 1 AND biomart_ortholog_type == "ortholog_one2one"
             - tier_H_orthofinder == 1 (MCL orthogroup + gene-tree reconciliation)
             - tier_H_mirbase == 1
             - tier_H_mirgenedb == 1

  Tier M  Medium confidence — sequence-validated by RBH (Wolf & Koonin 2012)
           OR Ensembl paralog-class orthology OR ortho2align statistical
           syntenic orthology (Mylarshchikov 2022).
           Sources allowed:
             - RBH-confirmed BLAST: tier_M_blast_rbh == 1
             - ortho2align: tier_M_ortho2align == 1 (liftOver + BLASTN within
               syntenic interval + shuffled-genome significance test)
             - Phase J synteny AND tier_M_blast_rbh == 1
             - Phase J synteny AND LncBook Q90 AND any BLAST hit (P1-7 fix:
               Q90 alone flags conservation "somewhere in mouse genome", not
               at the Phase J locus; requires BLAST to confirm co-location)
             - biomaRt paralog classes: one2many / many2one / many2many

  Tier L  Low confidence — positional only OR sequence-only without reciprocity
           Sources allowed:
             - Phase J synteny alone (no sequence confirmation)
             - non-RBH BLAST (sequence similarity without reciprocal best hit;
               paralogs and partial-region matches land here)
             - LncBook conservation alone (no synteny anchor)
             - Phase J + LncBook Q90 without BLAST (P1-7 demoted from Tier M)

TOGA design notes (2026-05-24):
  TOGA classifies all confirmed orthologs (one2one through many2many) as
  Tier H because the underlying signal — whole-genome Cactus chain + CESAR
  exon alignment + intactness ML classifier — is significantly more rigorous
  than biomaRt's gene-tree-based call. The orthology_class merely describes
  the multiplicity of orthologous loci, not the confidence of the calls.
  Paralog-class TOGA pairs carry a "paralog_branch" / "paralog_many2many"
  flag in the notes column for downstream consumers that prefer 1:1 mappings.

Phase 0 pipeline uses ortholog_bridge.py directly. This integration produces
the master table for downstream consumers (atlas, evidence cards, paper).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
LAYERS = PROJECT / "data/external/orthologs/layers"
CURATED = PROJECT / "data/external/orthologs/curated_dbs"
OUT_MASTER = PROJECT / "data/external/orthologs/master_ortholog_table.tsv.gz"


def strip_version(s):
    return s.astype(str).str.split(".").str[0]


def load_toga() -> pd.DataFrame:
    """TOGA Zoonomia mouse-human gold-standard orthology (L0).

    Source: http://genome.senckenberg.de/download/TOGA/human_hg38_reference/
            Rodentia/Mus_musculus__house_mouse__mm39/orthologsClassification.tsv.gz

    Cactus whole-genome alignment + CESAR exon projection + ML intactness
    classifier (Kirilenko et al. 2023 Science 380:eabn3107). Protein-coding
    only (no lncRNA / miRNA).

    Coverage: ~24K (mouse_ensembl, human_ensembl) pairs from ~19.5K human
    Ensembl reference genes.
    """
    p = LAYERS / "L0_toga.tsv"
    if not p.exists():
        print(f"[integrate] L0 TOGA not present at {p}; skipping",
              file=sys.stderr)
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    df["tier_H_toga"] = 1
    # Strip Ensembl version defensively
    for c in ("mouse_ensembl", "human_ensembl"):
        if c in df.columns:
            df[c] = df[c].astype(str).str.split(".").str[0]
            df[c] = df[c].replace({"nan": pd.NA})
    return df


def load_phasej() -> pd.DataFrame:
    p = LAYERS / "L1.5_phasej_synteny.tsv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    df["tier_M_phasej_synteny"] = 1
    return df


def load_blast() -> pd.DataFrame:
    p = LAYERS / "L4_blast.tsv"
    if not p.exists():
        print(f"[integrate] L4 BLAST not yet present at {p}; skipping",
              file=sys.stderr)
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    df["tier_M_blast"] = 1
    for c in ["blast_pident", "blast_length", "blast_evalue", "blast_bitscore"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_blast_rbh() -> pd.DataFrame:
    """L4r — reciprocal best hit BLAST pairs (subset of L4_blast.tsv).

    Returns an empty DataFrame if the file does not exist yet; downstream
    code degrades to v2 (no-RBH) semantics in that case.
    """
    p = LAYERS / "L4_blast_rbh.tsv"
    if not p.exists():
        print(f"[integrate] L4r BLAST RBH not present at {p}; "
              f"RBH-based tiering disabled (all BLAST hits stay at Tier L "
              f"unless they have other evidence).", file=sys.stderr)
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    df["tier_M_blast_rbh"] = 1
    for c in ("blast_pident_h2m", "blast_length_h2m", "blast_bitscore_h2m"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    # Keep only the columns we want to add to master (key cols + h2m fields +
    # the flag); the m->h pident/length/evalue/bitscore are already provided
    # by L4_blast.tsv via load_blast().
    keep = ["mouse_ensembl", "human_ensembl",
            "blast_pident_h2m", "blast_length_h2m", "blast_bitscore_h2m",
            "tier_M_blast_rbh"]
    keep = [c for c in keep if c in df.columns]
    return df[keep]


def load_mirbase() -> pd.DataFrame:
    """miRBase 22 mature-miRNA family-ID matches.

    The source TSV has one row per MATURE miRNA arm (5p / 3p / unknown),
    so a single (mouse_gene, human_gene) pair can appear multiple times.
    We collapse at the GENE level (Bug 2 fix), concatenating arm/MIMAT/
    mature-name fields with comma separators to preserve provenance.
    """
    p = LAYERS / "L2_mirbase.tsv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    df["tier_H_mirbase"] = 1
    # Strip Ensembl version
    for c in ["mouse_ensembl", "human_ensembl"]:
        if c in df.columns:
            df[c] = df[c].astype(str).str.split(".").str[0]
            df[c] = df[c].replace({"nan": pd.NA})
    # Collapse 5p/3p/unknown arm rows for the same (mouse, human) gene pair.
    # We only collapse rows that have BOTH ensembl IDs filled (the others
    # would not merge into master anyway, and may double-count gene-symbol
    # rows that lack ensembl IDs).
    keys = ["mouse_ensembl", "human_ensembl"]
    both = df.dropna(subset=keys, how="any").copy()
    missing = df[df["mouse_ensembl"].isna() | df["human_ensembl"].isna()].copy()
    # Concat arm-level columns; first non-null for symbol/biotype/family.
    agg_concat_cols = [
        c for c in ["mirbase_arm", "mirbase_human_mimat",
                    "mirbase_mouse_mimat",
                    "mirbase_human_mature_name",
                    "mirbase_mouse_mature_name"]
        if c in both.columns]
    agg_first_cols = [
        c for c in ["mouse_symbol", "mouse_biotype",
                    "human_symbol", "human_biotype",
                    "mirbase_family", "notes",
                    "tier_H_mirbase", "tier_H_biomart",
                    "tier_H_mirgenedb"]
        if c in both.columns]
    agg_dict = {c: (lambda s: ",".join(sorted(set(
        x for x in s.dropna().astype(str) if x != "nan"))))
        for c in agg_concat_cols}
    agg_dict.update({c: "first" for c in agg_first_cols})
    if not both.empty:
        collapsed = both.groupby(keys, as_index=False).agg(agg_dict)
    else:
        collapsed = both
    return pd.concat([collapsed, missing], ignore_index=True, sort=False)


def load_mirgenedb() -> pd.DataFrame:
    """MirGeneDB 2.1 strict miRNA family+paralog+arm ortholog pairs (L2b).

    MirGeneDB applies stricter criteria than miRBase: each entry must satisfy
    phylogeny anchoring, sequence support, and read-count thresholds. The
    matched pairs across species are therefore considered gold-standard
    miRNA orthologs (Fromm et al. 2022 NAR).

    Like the miRBase loader, we collapse 5p/3p arm rows at the gene level
    so a single (mouse_ensembl, human_ensembl) pair is not duplicated.
    """
    p = LAYERS / "L2b_mirgenedb.tsv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    df["tier_H_mirgenedb"] = 1
    for c in ["mouse_ensembl", "human_ensembl"]:
        if c in df.columns:
            df[c] = df[c].astype(str).str.split(".").str[0]
            df[c] = df[c].replace({"nan": pd.NA})
    keys = ["mouse_ensembl", "human_ensembl"]
    both = df.dropna(subset=keys, how="any").copy()
    missing = df[df["mouse_ensembl"].isna() | df["human_ensembl"].isna()].copy()
    agg_concat_cols = [
        c for c in ["mirgenedb_arm", "mirgenedb_human_hairpin",
                    "mirgenedb_mouse_hairpin",
                    "mirgenedb_human_mature",
                    "mirgenedb_mouse_mature",
                    "mirgenedb_human_mirbase_hairpin",
                    "mirgenedb_mouse_mirbase_hairpin"]
        if c in both.columns]
    agg_first_cols = [
        c for c in ["mouse_symbol", "mouse_biotype",
                    "human_symbol", "human_biotype",
                    "mirgenedb_family", "mirgenedb_paralog",
                    "mirgenedb_variant", "notes",
                    "tier_H_mirgenedb"]
        if c in both.columns]
    agg_dict = {c: (lambda s: ",".join(sorted(set(
        x for x in s.dropna().astype(str) if x != "nan"))))
        for c in agg_concat_cols}
    agg_dict.update({c: "first" for c in agg_first_cols})
    if not both.empty:
        collapsed = both.groupby(keys, as_index=False).agg(agg_dict)
    else:
        collapsed = both
    return pd.concat([collapsed, missing], ignore_index=True, sort=False)


def load_ortho2align() -> pd.DataFrame:
    """L9 — ortho2align syntenic lncRNA ortholog disambiguation.

    Source: ortho2align v1.0.5 (Mylarshchikov 2022 BMC Bioinformatics).
    LiftOver + BLASTN within syntenic intervals + shuffled-genome statistical
    significance test. Picks the single best mouse ortholog per human lncRNA
    anchor. Upgrades Phase J "all candidates in syntenic interval" to
    "best candidate per anchor with statistical significance."

    Output: one row per (human_ensembl, mouse_ensembl) pair where mouse_ensembl
    is the best-annotated ortho2align hit. Unannotated hits (ortholog region
    doesn't overlap a known mouse gene) are included with mouse_ensembl = NA.
    """
    p = LAYERS / "L9_ortho2align.tsv"
    if not p.exists():
        print(f"[integrate] L9 ortho2align not present at {p}; skipping",
              file=sys.stderr)
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    df["tier_M_ortho2align"] = 1
    for c in ("mouse_ensembl", "human_ensembl"):
        if c in df.columns:
            df[c] = df[c].astype(str).str.split(".").str[0]
            df[c] = df[c].replace({"nan": pd.NA, "<NA>": pd.NA})
    for c in ("ortho2align_score", "ortho2align_min_qvalue",
              "ortho2align_aligned_bp", "ortho2align_coverage",
              "ortho2align_jaccard", "ortho2align_overlap_coeff"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_mmseqs2() -> pd.DataFrame:
    """L5 — MMseqs2 reciprocal-best-hit sequence orthology.

    Source: MMseqs2 easy-rbh (Steinegger & Soeding 2017 Nat Biotechnol).
    Ultra-fast profile-based nucleotide RBH. Independent sequence aligner to
    BLAST (different seeding/scoring), but the same evidence CLASS (sequence
    similarity) — counted with BLAST in independent_evidence_count, not added.
    """
    p = LAYERS / "L5_mmseqs2_rbh.tsv"
    if not p.exists():
        print(f"[integrate] L5 MMseqs2 not present at {p}; skipping",
              file=sys.stderr)
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    df["tier_M_mmseqs2_rbh"] = 1
    for c in ("mouse_ensembl", "human_ensembl"):
        if c in df.columns:
            df[c] = df[c].astype(str).str.split(".").str[0]
            df[c] = df[c].replace({"nan": pd.NA, "<NA>": pd.NA})
    for c in ("mmseqs2_pident", "mmseqs2_length", "mmseqs2_evalue",
              "mmseqs2_bitscore"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    # Rename repeat flag to avoid collision with BLAST's repeat_risk_flag.
    if "repeat_risk_flag" in df.columns:
        df = df.rename(columns={"repeat_risk_flag": "mmseqs2_repeat_risk_flag"})
    return df


def load_seekr() -> pd.DataFrame:
    """L6 — Seekr k-mer functional/compositional orthology.

    Source: Seekr (Kirk 2018 Mol Cell). Pearson correlation of k-mer (k=5,6)
    frequency profiles between mouse/human lncRNAs. Captures compositional
    orthology (shared RBP-binding motif content) for alignment-dark lncRNAs.
    A genuinely independent evidence class (not sequence alignment, not
    synteny). Layer pre-assigns tier_M_seekr (Pearson > 0.8) vs tier_L_seekr.
    """
    p = LAYERS / "L6_seekr_kmer.tsv"
    if not p.exists():
        print(f"[integrate] L6 Seekr not present at {p}; skipping",
              file=sys.stderr)
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    for c in ("mouse_ensembl", "human_ensembl"):
        if c in df.columns:
            df[c] = df[c].astype(str).str.split(".").str[0]
            df[c] = df[c].replace({"nan": pd.NA, "<NA>": pd.NA})
    for c in ("tier_M_seekr", "tier_L_seekr"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    for c in ("seekr_pearson_k5", "seekr_pearson_k6", "seekr_pearson_max"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_noncode() -> pd.DataFrame:
    """L7 — NONCODE v6 curated cross-species lncRNA orthology.

    Source: NONCODE v6 (Zhao 2021 NAR) via RNAcentral ID mapping. Curated
    multi-species lncRNA conservation. Same evidence CLASS as LncBook
    (curated conservation DB) — counted with LncBook, not added.
    """
    p = LAYERS / "L7_noncode.tsv"
    if not p.exists():
        print(f"[integrate] L7 NONCODE not present at {p}; skipping",
              file=sys.stderr)
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    df["tier_M_noncode"] = 1
    for c in ("mouse_ensembl", "human_ensembl"):
        if c in df.columns:
            df[c] = df[c].astype(str).str.split(".").str[0]
            df[c] = df[c].replace({"nan": pd.NA, "<NA>": pd.NA})
    return df


def load_orthofinder() -> pd.DataFrame:
    """OrthoFinder gene-tree reconciliation pairs (L8).

    Source: OrthoFinder 3.1.2 (Emms & Kelly 2019 Genome Biology 20:238).
    DIAMOND all-vs-all + MCL orthogroup inference + MAFFT MSA + FastTree
    gene trees + STRIDE species-tree rooting + gene-tree/species-tree
    reconciliation.

    Protein-coding only (input = GENCODE longest isoform per gene).
    Pairs classified by orthogroup multiplicity: one2one, one2many,
    many2one, many2many.
    """
    p = LAYERS / "L8_orthofinder.tsv"
    if not p.exists():
        print(f"[integrate] L8 OrthoFinder not present at {p}; skipping",
              file=sys.stderr)
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    df["tier_H_orthofinder"] = 1
    for c in ("mouse_ensembl", "human_ensembl"):
        if c in df.columns:
            df[c] = df[c].astype(str).str.split(".").str[0]
            df[c] = df[c].replace({"nan": pd.NA})
    return df


def _load_human_biotype() -> pd.DataFrame:
    """GENCODE v49 human gene metadata → ensembl_base, gene_biotype."""
    p = PROJECT / "data/gencode_v49_gene_metadata.tsv.gz"
    if not p.exists():
        return pd.DataFrame(columns=["human_ensembl", "human_biotype"])
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    df = df.rename(columns={"ensembl_base": "human_ensembl",
                            "gene_biotype": "human_biotype"})
    return df[["human_ensembl", "human_biotype"]].dropna(
        subset=["human_ensembl"]).drop_duplicates(subset=["human_ensembl"])


def _load_mouse_biotype() -> pd.DataFrame:
    """GENCODE vM38 mouse gene metadata → mouse_ensembl_base, mouse_biotype."""
    p = (PROJECT / "Cas13_Library_Design/data"
         / "mouse_gencode_vM38_gene_metadata.csv")
    if not p.exists():
        return pd.DataFrame(columns=["mouse_ensembl", "mouse_biotype"])
    df = pd.read_csv(p, dtype=str, low_memory=False)
    df = df.rename(columns={"mouse_ensembl_base": "mouse_ensembl"})
    return df[["mouse_ensembl", "mouse_biotype"]].dropna(
        subset=["mouse_ensembl"]).drop_duplicates(subset=["mouse_ensembl"])


def load_biomart() -> pd.DataFrame:
    p = PROJECT / "data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz"
    df = pd.read_csv(p, sep="\t", dtype=str)
    df = df.rename(columns={
        "human_gene_symbol": "human_symbol",
        "human_gene_ensembl": "human_ensembl",
        "mouse_gene_symbol": "mouse_symbol",
        "mouse_gene_ensembl": "mouse_ensembl",
    })
    df["mouse_ensembl"] = strip_version(df["mouse_ensembl"])
    df["human_ensembl"] = strip_version(df["human_ensembl"])
    df["tier_H_biomart"] = 1
    df["biomart_ortholog_type"] = df["ortholog_type"]
    # Join biotypes from GENCODE metadata (Bug 3 fix — source TSV does not
    # carry biotype, so without this join all biomaRt rows have NaN biotype
    # and any downstream biotype filter silently drops them).
    hb = _load_human_biotype()
    mb = _load_mouse_biotype()
    if not hb.empty:
        df = df.merge(hb, on="human_ensembl", how="left")
    else:
        df["human_biotype"] = pd.NA
    if not mb.empty:
        df = df.merge(mb, on="mouse_ensembl", how="left")
    else:
        df["mouse_biotype"] = pd.NA
    return df[["mouse_ensembl", "mouse_symbol", "mouse_biotype",
               "human_ensembl", "human_symbol", "human_biotype",
               "tier_H_biomart", "biomart_ortholog_type"]]


def load_lncbook_conservation() -> pd.DataFrame:
    """LncBook conservation table (mouse rows only). Returns one row per
    LncBook human lncRNA with its conservation Q-level vs mouse."""
    p = CURATED / "lncbook_conservation_v2.0.csv.gz"
    if not p.exists():
        return pd.DataFrame()
    # Read only mouse rows
    df = pd.read_csv(p, low_memory=False)
    df = df[df["Species"].astype(str).str.contains("Mus musculus", na=False)]
    df = df.rename(columns={
        "Gene ID": "lncbook_gene_id",
        "Symbol": "lncbook_symbol",
        "Conservation": "lncbook_mouse_conservation",
        "Identity": "lncbook_mouse_identity",
        "Alignment Length": "lncbook_mouse_alignment_length",
    })
    keep = ["lncbook_gene_id", "lncbook_symbol", "lncbook_mouse_conservation",
            "lncbook_mouse_identity", "lncbook_mouse_alignment_length"]
    return df[keep]


def assign_tier(row) -> str:
    """Assign confidence tier per row.

    Tier H: TOGA Zoonomia (any class), biomaRt strict 1:1 (one2one),
            OrthoFinder gene-tree reconciliation (Emms & Kelly 2019),
            miRBase family, or MirGeneDB family.
    Tier M: RBH-confirmed BLAST (when L4r layer available)
            OR strong BLAST (pident ≥ 70 AND length ≥ 200) as the no-RBH
            fallback (Bug 1 spec from user)
            OR ortho2align statistically significant syntenic orthology
            (Mylarshchikov 2022; liftOver + BLASTN + shuffled-genome test)
            OR Phase J + LncBook Q90 + any BLAST hit (P1-7 fix: Q90 alone
              flags conservation "somewhere in mouse genome", not at the
              specific Phase J locus; 90.5% of Q90-only pairs are
              Gm-prefixed predicted genes)
            OR biomaRt paralog classes (one2many / many2one / many2many).
    Tier L: Phase J alone, weak BLAST alone, LncBook conservation alone,
            Phase J + LncBook Q90 without BLAST confirmation.

    Notes:
      - RBH (Wolf & Koonin 2012) takes precedence over strong-BLAST when the
        L4r layer is present (it removes paralog inflation that pident≥70
        alone cannot detect). When L4r is absent, strong-BLAST is the M-tier
        anchor.
      - TOGA (Kirilenko 2023) is the gold-standard mammalian ortholog caller
        — Cactus whole-genome alignment + CESAR exon projection + ML
        intactness classifier — trusted at Tier H for all four classes;
        paralog status is preserved in notes / provenance.
    """
    toga = row.get("tier_H_toga", 0) == 1
    biomart = row.get("tier_H_biomart", 0) == 1
    biomart_type = str(row.get("biomart_ortholog_type", "") or "")
    orthofinder = row.get("tier_H_orthofinder", 0) == 1
    mirbase = row.get("tier_H_mirbase", 0) == 1
    mirgenedb = row.get("tier_H_mirgenedb", 0) == 1
    blast = row.get("tier_M_blast", 0) == 1
    blast_rbh = row.get("tier_M_blast_rbh", 0) == 1
    phasej = row.get("tier_M_phasej_synteny", 0) == 1
    ortho2align = row.get("tier_M_ortho2align", 0) == 1
    mmseqs2_rbh = row.get("tier_M_mmseqs2_rbh", 0) == 1
    mmseqs2_repeat_suspect = (
        str(row.get("mmseqs2_repeat_risk_flag", "") or "") == "high_risk")
    seekr_m = row.get("tier_M_seekr", 0) == 1
    seekr_l = row.get("tier_L_seekr", 0) == 1
    noncode = row.get("tier_M_noncode", 0) == 1
    lncbook_q = str(row.get("lncbook_mouse_conservation", ""))
    lncbook_q90 = lncbook_q == "Q90"
    # P0-6 fix: repeat-risk flag. At 90 Myr divergence, short high-identity
    # fragments are likely LINE/SINE repeats, not true orthology. Pairs
    # flagged high_risk should NOT promote to Tier M via BLAST-only evidence
    # (RBH or strong-BLAST). Non-BLAST evidence (TOGA, biomaRt, miRBase) is
    # unaffected since it does not rely on sequence alignment length.
    repeat_risk = str(row.get("repeat_risk_flag", "") or "")
    blast_is_repeat_suspect = (repeat_risk == "high_risk")
    strong_blast = blast and not blast_is_repeat_suspect and (
        float(row.get("blast_pident", 0) or 0) >= 70
        and float(row.get("blast_length", 0) or 0) >= 200
    )
    # ---- Tier H ---------------------------------------------------------
    if toga:
        return "H"
    if biomart and biomart_type == "ortholog_one2one":
        return "H"
    if orthofinder:
        return "H"
    if mirbase:
        return "H"
    if mirgenedb:
        return "H"
    # ---- Tier M ---------------------------------------------------------
    # biomaRt paralog classes: real homology but ambiguous mapping
    if biomart and biomart_type in (
            "ortholog_one2many", "ortholog_many2one",
            "ortholog_many2many"):
        return "M"
    if blast_rbh and not blast_is_repeat_suspect:
        return "M"  # RBH-confirmed sequence orthology (alone is enough)
    if mmseqs2_rbh and not mmseqs2_repeat_suspect:
        return "M"  # MMseqs2 RBH: independent fast-aligner sequence orthology
    if strong_blast:
        return "M"  # strong-BLAST fallback when L4r unavailable
    if ortho2align:
        return "M"  # ortho2align: statistically significant syntenic orthology
    if noncode:
        return "M"  # NONCODE v6 curated cross-species lncRNA conservation
    if seekr_m:
        return "M"  # Seekr k-mer functional orthology (Pearson > 0.8)
    # P1-7 fix: Phase J + LncBook Q90 requires BLAST confirmation.
    # Without BLAST, Q90 flags conservation somewhere in the mouse genome
    # but not necessarily at the syntenic locus. Demote to Tier L.
    if phasej and lncbook_q90 and blast and not blast_is_repeat_suspect:
        return "M"  # syntenic + LncBook Q90 + BLAST-confirmed at locus
    # ---- Tier L ---------------------------------------------------------
    # Phase J only, weak BLAST only, LncBook only, Phase J+Q90 no BLAST,
    # or BLAST-only with high_risk repeat flag.
    return "L"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=OUT_MASTER)
    p.add_argument("--require-blast", action="store_true",
                   help="abort if L4_blast.tsv not present (default: continue without)")
    args = p.parse_args()

    print("[integrate] loading biomaRt baseline...", file=sys.stderr)
    bm = load_biomart()
    print(f"[integrate]  biomaRt pairs: {len(bm)}", file=sys.stderr)

    print("[integrate] loading TOGA Zoonomia (L0)...", file=sys.stderr)
    tg = load_toga()
    print(f"[integrate]  TOGA pairs: {len(tg)}", file=sys.stderr)

    print("[integrate] loading Phase J synteny...", file=sys.stderr)
    pj = load_phasej()
    print(f"[integrate]  Phase J pairs: {len(pj)}", file=sys.stderr)

    print("[integrate] loading BLAST L4...", file=sys.stderr)
    bl = load_blast()
    if bl.empty and args.require_blast:
        print("[integrate] ERROR: --require-blast set but L4_blast.tsv missing",
              file=sys.stderr)
        return 1
    print(f"[integrate]  BLAST pairs: {len(bl)}", file=sys.stderr)

    print("[integrate] loading LncBook conservation...", file=sys.stderr)
    lb = load_lncbook_conservation()
    print(f"[integrate]  LncBook (mouse-conserved) rows: {len(lb)}", file=sys.stderr)

    # Outer-merge all layers on (mouse_ensembl, human_ensembl)
    keys = ["mouse_ensembl", "human_ensembl"]
    # Bug 3 fix: carry biotype from biomaRt seed (load_biomart joined GENCODE
    # metadata) so biomaRt rows aren't NaN-biotype.
    master = bm[keys + ["tier_H_biomart", "biomart_ortholog_type",
                        "mouse_symbol", "human_symbol",
                        "mouse_biotype", "human_biotype"]].copy()

    if not tg.empty:
        tg_keep = ["mouse_ensembl", "human_ensembl",
                   "mouse_symbol", "human_symbol",
                   "mouse_biotype",
                   "tier_H_toga", "toga_orthology_class"]
        tg_keep = [c for c in tg_keep if c in tg.columns]
        # Drop any TOGA rows missing either Ensembl ID (can't merge)
        tg_use = tg.dropna(subset=keys, how="any")[tg_keep].copy()
        master = master.merge(tg_use, on=keys, how="outer",
                              suffixes=("", "_tg"))
        for col in ("mouse_symbol", "human_symbol", "mouse_biotype"):
            if f"{col}_tg" in master.columns:
                master[col] = master[col].fillna(master[f"{col}_tg"])
                master = master.drop(columns=f"{col}_tg")

    if not pj.empty:
        pj_use = pj[keys + ["mouse_symbol", "human_symbol",
                            "mouse_biotype", "human_biotype",
                            "tier_M_phasej_synteny"]].copy()
        master = master.merge(pj_use, on=keys, how="outer",
                              suffixes=("", "_pj"))
        # consolidate symbols / biotypes (fill biomaRt-side NaNs from phasej)
        for col in ["mouse_symbol", "human_symbol",
                    "mouse_biotype", "human_biotype"]:
            if f"{col}_pj" in master.columns:
                master[col] = master[col].fillna(master[f"{col}_pj"])
                master = master.drop(columns=f"{col}_pj")

    if not bl.empty:
        bl_keep = ["mouse_ensembl", "human_ensembl",
                   "mouse_symbol", "human_symbol",
                   "mouse_biotype", "human_biotype",
                   "tier_M_blast",
                   "blast_pident", "blast_length", "blast_evalue",
                   "blast_bitscore", "repeat_risk_flag"]
        # Only keep columns that actually exist (repeat_risk_flag added in
        # P0-6 fix; older L4_blast.tsv files may lack it).
        bl_keep = [c for c in bl_keep if c in bl.columns]
        master = master.merge(bl[bl_keep], on=keys, how="outer",
                              suffixes=("", "_bl"))
        # Consolidate symbols / biotypes
        for col in ["mouse_symbol", "human_symbol", "mouse_biotype",
                    "human_biotype"]:
            if f"{col}_bl" in master.columns:
                master[col] = master[col].fillna(master[f"{col}_bl"])
                master = master.drop(columns=f"{col}_bl")

    print("[integrate] loading BLAST L4r (RBH)...", file=sys.stderr)
    rbh = load_blast_rbh()
    if not rbh.empty:
        print(f"[integrate]  BLAST RBH pairs: {len(rbh)}", file=sys.stderr)
        master = master.merge(rbh, on=keys, how="left")
    else:
        # L4_blast_rbh.tsv missing — keep schema columns so downstream code
        # can rely on them existing even on a no-RBH run.
        for c in ("blast_pident_h2m", "blast_length_h2m",
                  "blast_bitscore_h2m", "tier_M_blast_rbh"):
            if c not in master.columns:
                master[c] = pd.NA

    print("[integrate] loading miRBase 22 family-ID pairs...", file=sys.stderr)
    mb = load_mirbase()
    if not mb.empty:
        print(f"[integrate]  miRBase pairs: {len(mb)}", file=sys.stderr)
        mb_keep = ["mouse_ensembl", "human_ensembl",
                   "mouse_symbol", "human_symbol",
                   "mouse_biotype", "human_biotype",
                   "tier_H_mirbase",
                   "mirbase_family", "mirbase_arm"]
        # Drop rows missing both Ensembl IDs (can't merge)
        mb = mb.dropna(subset=keys, how="any")
        master = master.merge(mb[mb_keep], on=keys, how="outer",
                              suffixes=("", "_mb"))
        for col in ["mouse_symbol", "human_symbol", "mouse_biotype",
                    "human_biotype"]:
            if f"{col}_mb" in master.columns:
                master[col] = master[col].fillna(master[f"{col}_mb"])
                master = master.drop(columns=f"{col}_mb")

    print("[integrate] loading MirGeneDB 2.1 strict pairs (L2b)...",
          file=sys.stderr)
    mg = load_mirgenedb()
    if not mg.empty:
        print(f"[integrate]  MirGeneDB pairs: {len(mg)}", file=sys.stderr)
        mg_keep = ["mouse_ensembl", "human_ensembl",
                   "mouse_symbol", "human_symbol",
                   "mouse_biotype", "human_biotype",
                   "tier_H_mirgenedb",
                   "mirgenedb_family", "mirgenedb_paralog",
                   "mirgenedb_variant", "mirgenedb_arm"]
        mg_keep = [c for c in mg_keep if c in mg.columns]
        # Drop rows missing both Ensembl IDs (can't merge into master keyed on
        # ensembl); we still record symbol-only pairs through the standalone
        # L2b_mirgenedb.tsv if downstream needs them.
        mg_use = mg.dropna(subset=keys, how="any")[mg_keep]
        master = master.merge(mg_use, on=keys, how="outer",
                              suffixes=("", "_mg"))
        for col in ["mouse_symbol", "human_symbol", "mouse_biotype",
                    "human_biotype"]:
            if f"{col}_mg" in master.columns:
                master[col] = master[col].fillna(master[f"{col}_mg"])
                master = master.drop(columns=f"{col}_mg")

    print("[integrate] loading OrthoFinder gene-tree reconciliation (L8)...",
          file=sys.stderr)
    of = load_orthofinder()
    if not of.empty:
        print(f"[integrate]  OrthoFinder pairs: {len(of)}", file=sys.stderr)
        of_keep = ["mouse_ensembl", "human_ensembl",
                   "mouse_symbol", "human_symbol",
                   "mouse_biotype", "human_biotype",
                   "tier_H_orthofinder",
                   "orthofinder_orthogroup",
                   "orthofinder_orthology_type"]
        of_keep = [c for c in of_keep if c in of.columns]
        of_use = of.dropna(subset=keys, how="any")[of_keep]
        master = master.merge(of_use, on=keys, how="outer",
                              suffixes=("", "_of"))
        for col in ["mouse_symbol", "human_symbol", "mouse_biotype",
                    "human_biotype"]:
            if f"{col}_of" in master.columns:
                master[col] = master[col].fillna(master[f"{col}_of"])
                master = master.drop(columns=f"{col}_of")
    else:
        # Ensure schema columns exist even when layer is absent
        for c in ("tier_H_orthofinder", "orthofinder_orthogroup",
                  "orthofinder_orthology_type"):
            if c not in master.columns:
                master[c] = pd.NA

    print("[integrate] loading ortho2align syntenic disambiguation (L9)...",
          file=sys.stderr)
    o2a = load_ortho2align()
    if not o2a.empty:
        print(f"[integrate]  ortho2align pairs: {len(o2a)}", file=sys.stderr)
        o2a_keep = ["mouse_ensembl", "human_ensembl",
                    "mouse_symbol", "human_symbol",
                    "mouse_biotype", "human_biotype",
                    "tier_M_ortho2align",
                    "ortho2align_score", "ortho2align_min_qvalue",
                    "ortho2align_aligned_bp", "ortho2align_coverage",
                    "ortho2align_jaccard", "ortho2align_overlap_coeff"]
        o2a_keep = [c for c in o2a_keep if c in o2a.columns]
        # Only merge rows that have at least human_ensembl
        o2a_use = o2a.dropna(subset=["human_ensembl"])[o2a_keep].copy()
        # ortho2align unannotated hits (mouse_ensembl=NA) can't merge on both
        # keys; merge annotated rows on (mouse, human), unannotated on human only
        o2a_ann = o2a_use.dropna(subset=keys)
        o2a_noann = o2a_use[o2a_use["mouse_ensembl"].isna()]
        if not o2a_ann.empty:
            master = master.merge(o2a_ann, on=keys, how="outer",
                                  suffixes=("", "_o2a"))
            for col in ["mouse_symbol", "human_symbol", "mouse_biotype",
                        "human_biotype"]:
                if f"{col}_o2a" in master.columns:
                    master[col] = master[col].fillna(master[f"{col}_o2a"])
                    master = master.drop(columns=f"{col}_o2a")
        # Store count of unannotated hits for stats (don't merge — no mouse_ensembl)
        if not o2a_noann.empty:
            print(f"[integrate]  ortho2align unannotated hits (novel syntenic, "
                  f"no known mouse gene): {len(o2a_noann)}",
                  file=sys.stderr)
    else:
        for c in ("tier_M_ortho2align", "ortho2align_score",
                  "ortho2align_min_qvalue", "ortho2align_aligned_bp",
                  "ortho2align_coverage", "ortho2align_jaccard",
                  "ortho2align_overlap_coeff"):
            if c not in master.columns:
                master[c] = pd.NA

    print("[integrate] loading MMseqs2 RBH (L5)...", file=sys.stderr)
    mm = load_mmseqs2()
    if not mm.empty:
        print(f"[integrate]  MMseqs2 RBH pairs: {len(mm)}", file=sys.stderr)
        mm_keep = ["mouse_ensembl", "human_ensembl",
                   "mouse_symbol", "human_symbol",
                   "mouse_biotype", "human_biotype",
                   "tier_M_mmseqs2_rbh",
                   "mmseqs2_pident", "mmseqs2_length", "mmseqs2_evalue",
                   "mmseqs2_bitscore", "mmseqs2_repeat_risk_flag"]
        mm_keep = [c for c in mm_keep if c in mm.columns]
        mm_use = mm.dropna(subset=keys, how="any")[mm_keep]
        master = master.merge(mm_use, on=keys, how="outer", suffixes=("", "_mm"))
        for col in ["mouse_symbol", "human_symbol", "mouse_biotype",
                    "human_biotype"]:
            if f"{col}_mm" in master.columns:
                master[col] = master[col].fillna(master[f"{col}_mm"])
                master = master.drop(columns=f"{col}_mm")
    else:
        for c in ("tier_M_mmseqs2_rbh", "mmseqs2_pident", "mmseqs2_length",
                  "mmseqs2_evalue", "mmseqs2_bitscore",
                  "mmseqs2_repeat_risk_flag"):
            if c not in master.columns:
                master[c] = pd.NA

    print("[integrate] loading Seekr k-mer (L6)...", file=sys.stderr)
    sk = load_seekr()
    if not sk.empty:
        print(f"[integrate]  Seekr pairs: {len(sk)}", file=sys.stderr)
        sk_keep = ["mouse_ensembl", "human_ensembl",
                   "mouse_symbol", "human_symbol",
                   "mouse_biotype", "human_biotype",
                   "tier_M_seekr", "tier_L_seekr",
                   "seekr_pearson_k5", "seekr_pearson_k6", "seekr_pearson_max",
                   "seekr_is_rbh"]
        sk_keep = [c for c in sk_keep if c in sk.columns]
        sk_use = sk.dropna(subset=keys, how="any")[sk_keep]
        master = master.merge(sk_use, on=keys, how="outer", suffixes=("", "_sk"))
        for col in ["mouse_symbol", "human_symbol", "mouse_biotype",
                    "human_biotype"]:
            if f"{col}_sk" in master.columns:
                master[col] = master[col].fillna(master[f"{col}_sk"])
                master = master.drop(columns=f"{col}_sk")
    else:
        for c in ("tier_M_seekr", "tier_L_seekr", "seekr_pearson_k5",
                  "seekr_pearson_k6", "seekr_pearson_max", "seekr_is_rbh"):
            if c not in master.columns:
                master[c] = pd.NA

    print("[integrate] loading NONCODE v6 (L7)...", file=sys.stderr)
    nc = load_noncode()
    if not nc.empty:
        print(f"[integrate]  NONCODE pairs: {len(nc)}", file=sys.stderr)
        nc_keep = ["mouse_ensembl", "human_ensembl",
                   "mouse_symbol", "human_symbol",
                   "mouse_biotype", "human_biotype",
                   "tier_M_noncode", "noncode_human_id", "noncode_mouse_id"]
        nc_keep = [c for c in nc_keep if c in nc.columns]
        nc_use = nc.dropna(subset=keys, how="any")[nc_keep]
        master = master.merge(nc_use, on=keys, how="outer", suffixes=("", "_nc"))
        for col in ["mouse_symbol", "human_symbol", "mouse_biotype",
                    "human_biotype"]:
            if f"{col}_nc" in master.columns:
                master[col] = master[col].fillna(master[f"{col}_nc"])
                master = master.drop(columns=f"{col}_nc")
    else:
        for c in ("tier_M_noncode", "noncode_human_id", "noncode_mouse_id"):
            if c not in master.columns:
                master[c] = pd.NA

    # Annotate LncBook conservation by joining on human gene symbol
    # (LncBook has its own IDs but Symbol column matches HGNC for ~4k pairs)
    if not lb.empty and "lncbook_symbol" in lb.columns:
        # LncBook entries with a real symbol
        lb_sym = lb[lb["lncbook_symbol"].astype(str) != "-"].copy()
        # Take strongest conservation per symbol
        lb_sym["_q_rank"] = lb_sym["lncbook_mouse_conservation"].map(
            {"Q90": 0, "Q75": 1, "Q50": 2, "Not conserved": 3}).fillna(99)
        lb_sym = lb_sym.sort_values("_q_rank").drop_duplicates(
            subset=["lncbook_symbol"], keep="first").drop(columns="_q_rank")
        lb_sym = lb_sym.rename(columns={"lncbook_symbol": "human_symbol"})
        master = master.merge(
            lb_sym[["human_symbol", "lncbook_mouse_conservation",
                    "lncbook_mouse_identity", "lncbook_mouse_alignment_length"]],
            on="human_symbol", how="left"
        )

    # Fill missing tier flags with 0
    for c in ["tier_H_toga", "tier_H_biomart", "tier_H_orthofinder",
              "tier_H_mirbase", "tier_H_mirgenedb",
              "tier_M_phasej_synteny", "tier_M_blast", "tier_M_blast_rbh",
              "tier_M_ortho2align", "tier_M_mmseqs2_rbh",
              "tier_M_seekr", "tier_L_seekr", "tier_M_noncode"]:
        if c in master.columns:
            master[c] = master[c].fillna(0).astype(int)
        else:
            master[c] = 0

    # Drop rows where both ensembl IDs are missing
    master = master.dropna(subset=keys, how="any").reset_index(drop=True)

    # Safety dedup on (mouse_ensembl, human_ensembl). The miRBase loader
    # already collapses 5p/3p arms at the gene level (Bug 2 fix), but this
    # catches any residual duplication from layer overlap.
    before_dedup = len(master)
    master = master.drop_duplicates(
        subset=keys, keep="first").reset_index(drop=True)
    if before_dedup != len(master):
        print(f"[integrate] safety dedup removed "
              f"{before_dedup - len(master)} duplicate (mouse, human) pairs",
              file=sys.stderr)

    # Tier assignment
    print("[integrate] assigning confidence tiers...", file=sys.stderr)
    master["confidence_tier"] = master.apply(assign_tier, axis=1)

    # ---- raw_evidence_count (backward compat, renamed from evidence_count) --
    master["raw_evidence_count"] = (
        master[["tier_H_toga", "tier_H_biomart", "tier_H_orthofinder",
                "tier_H_mirbase", "tier_H_mirgenedb",
                "tier_M_phasej_synteny", "tier_M_blast", "tier_M_blast_rbh",
                "tier_M_ortho2align", "tier_M_mmseqs2_rbh",
                "tier_M_seekr", "tier_M_noncode"]]
        .sum(axis=1).astype(int)
    )

    # ---- independent_evidence_count (P0-5 fix, updated v5 for L8) ----------
    # Counts truly independent evidence CLASSES, not individual sources:
    #   Class 1: Ensembl-derived (biomaRt OR TOGA — shared genome assembly lineage)
    #   Class 2: miRNA family-ID (miRBase OR MirGeneDB — shared precursor DB)
    #   Class 3: Synteny (Phase J OR ortho2align — ortho2align refines
    #            Phase J with alignment + statistical test, same evidence class)
    #   Class 4: Sequence (BLAST, BLAST RBH, or MMseqs2 RBH — all measure
    #            sequence similarity; different aligners, same evidence class)
    #   Class 5: Curated conservation (LncBook Q75/Q90 or NONCODE v6)
    #   Class 6: Gene-tree reconciliation (OrthoFinder — independent DIAMOND
    #            all-vs-all + MCL + gene-tree/species-tree reconciliation;
    #            orthogonal to assembly-based TOGA/biomaRt and raw BLAST)
    #   Class 7: Functional k-mer composition (Seekr — RBP-motif content
    #            similarity; orthogonal to alignment, synteny, and curation)
    _cls1 = ((master["tier_H_biomart"] == 1) | (master["tier_H_toga"] == 1)).astype(int)
    _cls2 = ((master["tier_H_mirbase"] == 1) | (master["tier_H_mirgenedb"] == 1)).astype(int)
    _cls3 = ((master["tier_M_phasej_synteny"] == 1) | (master["tier_M_ortho2align"] == 1)).astype(int)
    _cls4 = ((master["tier_M_blast"] == 1) | (master["tier_M_blast_rbh"] == 1)
             | (master["tier_M_mmseqs2_rbh"] == 1)).astype(int)
    _cls5 = (master["lncbook_mouse_conservation"].isin(["Q75", "Q90"])
             | (master["tier_M_noncode"] == 1)).astype(int)
    _cls6 = (master["tier_H_orthofinder"] == 1).astype(int)
    _cls7 = (master["tier_M_seekr"] == 1).astype(int)
    master["independent_evidence_count"] = (
        _cls1 + _cls2 + _cls3 + _cls4 + _cls5 + _cls6 + _cls7)

    # ---- is_one2one (P1-9 fix, updated v5 for L8) ---------------------------
    # Flag pairs where the mapping source indicates strict 1:1 orthology.
    # biomaRt: ortholog_type == "ortholog_one2one"
    # TOGA: toga_orthology_class == "one2one"
    # OrthoFinder: orthofinder_orthology_type == "one2one"
    # miRNA / lncRNA / BLAST-only pairs: False by convention (no 1:1 class).
    _biomart_121 = (master["biomart_ortholog_type"] == "ortholog_one2one")
    _toga_121 = (master.get("toga_orthology_class", pd.Series(dtype=str)) == "one2one") \
        if "toga_orthology_class" in master.columns \
        else pd.Series(False, index=master.index)
    _of_121 = (master.get("orthofinder_orthology_type",
                           pd.Series(dtype=str)) == "one2one") \
        if "orthofinder_orthology_type" in master.columns \
        else pd.Series(False, index=master.index)
    master["is_one2one"] = (_biomart_121 | _toga_121 | _of_121)

    # ---- v49_absent flag (P1-13 fix) -----------------------------------------
    # 27 human Ensembl gene IDs are absent from GENCODE v49 metadata
    # (retired/merged IDs including 8 MTRNR2L pseudogenes). Flag them
    # for downstream exclusion without deleting rows (preserve provenance).
    hb_v49 = _load_human_biotype()
    if not hb_v49.empty:
        v49_ids = set(hb_v49["human_ensembl"].dropna())
        master["v49_absent"] = (~master["human_ensembl"].isin(v49_ids)).astype(int)
        n_absent = master["v49_absent"].sum()
        n_absent_genes = master.loc[master["v49_absent"] == 1,
                                     "human_ensembl"].nunique()
        print(f"[integrate] v49_absent flag: {n_absent} rows "
              f"({n_absent_genes} unique human genes) absent from GENCODE v49",
              file=sys.stderr)
    else:
        master["v49_absent"] = 0

    # ---- provenance_sources -------------------------------------------------
    prov = []
    for _, r in master.iterrows():
        s = []
        if r.get("tier_H_toga"):
            cls = str(r.get("toga_orthology_class", "") or "")
            if cls and cls != "one2one":
                s.append(f"toga({cls})")
            else:
                s.append("toga")
        if r.get("tier_H_biomart"): s.append("biomart")
        if r.get("tier_H_orthofinder"):
            of_cls = str(r.get("orthofinder_orthology_type", "") or "")
            if of_cls and of_cls != "one2one":
                s.append(f"orthofinder({of_cls})")
            else:
                s.append("orthofinder")
        if r.get("tier_H_mirbase"): s.append("mirbase")
        if r.get("tier_H_mirgenedb"): s.append("mirgenedb")
        if r.get("tier_M_phasej_synteny"): s.append("phasej")
        if r.get("tier_M_ortho2align"): s.append("ortho2align")
        # Distinguish RBH-confirmed BLAST from plain BLAST hits in provenance.
        # P0-6: append repeat risk suffix when flagged high_risk.
        _rr = str(r.get("repeat_risk_flag", "") or "")
        _rr_suffix = ":repeat_suspect" if _rr == "high_risk" else ""
        if r.get("tier_M_blast_rbh"):
            s.append(f"blast_rbh{_rr_suffix}")
        elif r.get("tier_M_blast"):
            s.append(f"blast{_rr_suffix}")
        if r.get("tier_M_mmseqs2_rbh"):
            _mm_rr = str(r.get("mmseqs2_repeat_risk_flag", "") or "")
            s.append("mmseqs2_rbh" + (":repeat_suspect" if _mm_rr == "high_risk" else ""))
        if r.get("tier_M_seekr"):
            s.append("seekr")
        elif r.get("tier_L_seekr"):
            s.append("seekr_weak")
        if r.get("tier_M_noncode"):
            s.append("noncode")
        if str(r.get("lncbook_mouse_conservation", "")) in ("Q75", "Q90"):
            s.append("lncbook_" + str(r["lncbook_mouse_conservation"]).lower())
        prov.append(",".join(s) if s else "")
    master["provenance_sources"] = prov

    print("[integrate] tier distribution:", file=sys.stderr)
    print(master["confidence_tier"].value_counts().to_string(), file=sys.stderr)
    if "repeat_risk_flag" in master.columns:
        blast_rows = master[master["tier_M_blast"] == 1]
        if not blast_rows.empty:
            print("[integrate] repeat_risk_flag (BLAST pairs only):",
                  file=sys.stderr)
            print(blast_rows["repeat_risk_flag"].value_counts().to_string(),
                  file=sys.stderr)
            demoted = ((blast_rows["repeat_risk_flag"] == "high_risk")
                       & (blast_rows["confidence_tier"] == "L")).sum()
            print(f"[integrate]   high_risk pairs demoted to Tier L: {demoted}",
                  file=sys.stderr)
    print(f"[integrate] total pairs: {len(master)}", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    master.to_csv(args.out, sep="\t", index=False, compression="gzip"
                  if str(args.out).endswith(".gz") else None)
    print(f"[integrate] wrote {args.out}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
