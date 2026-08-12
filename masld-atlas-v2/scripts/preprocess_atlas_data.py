#!/usr/bin/env python3
"""
preprocess_atlas_data.py
========================
LEGACY PRE-RESOURCE PRODUCER. This rank-first data model is blocked by default;
set ALLOW_LEGACY_PORTAL_REBUILD=true only for provenance-only regeneration.
Publication data must come from the MASLD Gene Catalog contract and final
promotion.

Reads source CSVs from the MASLD project and produces the web-app data files
for the MASLD Atlas v2 web application:

  1. atlas.parquet       -- Full 27,187-gene atlas (all cols) for /downloads
  2. atlas_core.parquet  -- Compact ~48-col client table (contract §2): atlas
                            subset + joined treat_* (canonical DEG gate) +
                            joined convergence (rank/score/tier/concordance)
  3. gene_index.json     -- Legacy compact search index (retired from gene page)
  4. gene_symbols.json   -- Slim [{symbol, biotype}] for fuse.js autocomplete
  5. atlas_summary.json  -- Landing page statistics
  6. featured_genes.json -- Featured gene cards

Usage:
  python preprocess_atlas_data.py --output-dir ../public/data
"""

import argparse
import json
import math
import os
import sys
from pathlib import Path

if (
    __name__ == "__main__"
    and os.environ.get("ALLOW_LEGACY_PORTAL_REBUILD", "false").lower() != "true"
):
    raise SystemExit(
        "REFUSED: legacy atlas preprocessing is outside the standalone "
        "Resource contract. Use Plans 50/60."
    )

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


# ---------------------------------------------------------------------------
# Paths relative to PROJECT_ROOT
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)

ATLAS_CSV = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
# C2: canonical human bulk DEGs = limma-voom-qw C2 (was dream_results_ashr.csv, retired)
CANONICAL_DEG_CSV = (
    PROJECT_ROOT
    / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"
)
COLOC_GENE_CSV = (
    PROJECT_ROOT / "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"
)
CLINICAL_DRUGS_CSV = (
    PROJECT_ROOT
    / "RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv"
)
CONVERGENCE_CSV = (
    PROJECT_ROOT / "RNA-seq/results/multi_evidence/convergence_evidence.csv"
)
# Tier-1/2 liver-specific COLOC (canonical headline effector set: 473 SuSiE /
# 1,031 union). Keyed by ensembl (no symbol column).
TIER12_COLOC_CSV = (
    PROJECT_ROOT / "GWAS/finemapping/results/susie_coloc/gene_level_coloc_tier12.csv"
)
# Canonical drug-target classification (source of drug_dev_status headline
# 2 approved / 208 clinical / 1,504 preclinical over the full drug universe).
# Keyed by symbol; carries the authoritative drug_dev_status per gene.
DRUG_CLASS_TSV = (
    PROJECT_ROOT / "data/external/drug_targets/drug_target_classification.tsv"
)


# ---------------------------------------------------------------------------
# Featured gene definitions
# ---------------------------------------------------------------------------

FEATURED_GENES = {
    "THRB": {
        "tagline": "FDA-approved target (Resmetirom). GGT COLOC PP.H4 = 0.999.",
        "category": "Drug Target",
    },
    "PNPLA3": {
        "tagline": "Strongest MASLD genetic risk factor. GWAS-validated.",
        "category": "Genetic Risk",
    },
    "RORA": {
        "tagline": "Circadian-metabolic TF. 11 GWAS-ATAC motif disruptions.",
        "category": "Regulator",
    },
    "HSD17B13": {
        "tagline": "Protective variant. Multi-ancestry replicated.",
        "category": "Protective",
    },
    "PPARA": {
        "tagline": "Master metabolic regulator. OxPhos program hub.",
        "category": "Metabolic Regulator",
    },
    "DGAT2": {
        "tagline": "Onset-stage drug target. COLOC PP.H4 = 0.901.",
        "category": "Drug Target",
    },
    "ACSL4": {
        "tagline": "Ferroptosis driver. Upregulated at F2 switch.",
        "category": "Ferroptosis",
    },
    "HKDC1": {
        "tagline": "Top F2→F3 progression driver. UKBB-ALT PP.H4 = 0.993.",
        "category": "Progression",
    },
    "CFLAR": {
        "tagline": "Hepatocyte apoptosis modulator. Death-receptor signalling.",
        "category": "Cell Death",
    },
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def safe_json_value(v):
    """Convert NaN/Inf/-Inf to None for JSON serialization."""
    if v is None:
        return None
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return None
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
        if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
            return None
    if isinstance(v, np.bool_):
        return bool(v)
    return v


def normalize_to_percentile(series: pd.Series, pct: float = 99.0) -> pd.Series:
    """Normalize a series to 0-1 using the pct-th percentile as ceiling."""
    abs_vals = series.abs()
    ceiling = np.nanpercentile(abs_vals.dropna().values, pct)
    if ceiling == 0:
        return pd.Series(0.0, index=series.index)
    normed = abs_vals / ceiling
    return normed.clip(upper=1.0).fillna(0.0)


def coerce_bool(series: pd.Series) -> pd.Series:
    """Convert R-style TRUE/FALSE strings and mixed types to Python bool."""
    return series.map(
        lambda x: (
            True
            if str(x).strip().upper() in ("TRUE", "1", "1.0")
            else (
                False
                if str(x).strip().upper() in ("FALSE", "0", "0.0", "NAN", "NONE", "")
                else False
            )
        )
    )


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def load_data():
    """Load all source CSVs."""
    print(f"Loading atlas from {ATLAS_CSV} ...")
    atlas = pd.read_csv(ATLAS_CSV, low_memory=False)
    print(f"  Atlas shape: {atlas.shape}")
    assert {"bulk_padj", "bulk_logFC"} <= set(atlas.columns), (
        "C2: atlas missing bulk_* — rebuild 27a"
    )

    print(f"Loading canonical DEGs from {CANONICAL_DEG_CSV} ...")
    deg = pd.read_csv(CANONICAL_DEG_CSV)
    print(f"  Canonical DEG shape: {deg.shape}")

    print(f"Loading COLOC gene-level from {COLOC_GENE_CSV} ...")
    coloc = pd.read_csv(COLOC_GENE_CSV)
    print(f"  Coloc shape: {coloc.shape}")

    print(f"Loading clinical drugs from {CLINICAL_DRUGS_CSV} ...")
    drugs = pd.read_csv(CLINICAL_DRUGS_CSV, low_memory=False)
    print(f"  Drugs shape: {drugs.shape}")

    print(f"Loading convergence evidence from {CONVERGENCE_CSV} ...")
    convergence = pd.read_csv(CONVERGENCE_CSV, low_memory=False)
    print(f"  Convergence shape: {convergence.shape}")

    return atlas, deg, coloc, drugs, convergence


def compute_deg_count(deg: pd.DataFrame) -> int:
    """Canonical DEG count = effect-size-aware interval-null FDR gate at
    treat_fdr < 0.05 (McCarthy & Smyth 2009; lfc=0.25). This counts canonical
    DEG *rows* (transcript/id level) = 1,918. NB: collapsed to unique symbols it
    is 1,915 (3 symbols carry 2 significant ids each), which is what the symbol-
    keyed atlas / atlas_core is_deg flag reports.
    """
    if "treat_fdr" in deg.columns:
        return int((deg["treat_fdr"] < 0.05).sum())
    # Legacy fallback (should not trigger for canonical_deg_results.csv)
    return int(((deg["padj"] < 0.05) & (deg["logFC"].abs() > 0.3)).sum())


def compute_deg_symbol_count(deg: pd.DataFrame) -> int:
    """Unique-symbol DEG count = the symbol-level collapse of the canonical gate
    (treat_fdr<0.05). = 1,915 (the 1,918 canonical DEG *rows* map to 1,915 unique
    symbols; 3 symbols carry 2 significant ids). This is what the symbol-keyed
    atlas / atlas_core is_deg flag reports, distinct from the 1,918 headline.
    """
    if "treat_fdr" not in deg.columns or "symbol" not in deg.columns:
        return 0
    return int(deg.loc[deg["treat_fdr"] < 0.05, "symbol"].dropna().nunique())


def merge_gate_columns(atlas: pd.DataFrame, deg: pd.DataFrame) -> pd.DataFrame:
    """Merge the canonical effect-size-aware interval-null FDR gate columns
    (treat_lfc / treat_p / treat_fdr) from canonical_deg_results.csv onto the
    atlas, keyed by gene symbol, so the emitted parquet exposes the numeric
    DEG-gate fields the web app queries. Numeric only; never surfaced by name
    in the UI. The atlas itself only carries the pre-summarized bulk_treat_fdr,
    so these three per-gene columns must come from the canonical DEG table.
    """
    gate_cols = ["treat_lfc", "treat_p", "treat_fdr"]
    present = [c for c in gate_cols if c in deg.columns]
    if (
        not present
        or "human_symbol" not in atlas.columns
        or "symbol" not in deg.columns
    ):
        print(f"  [warn] gate columns not merged (present={present})")
        return atlas
    # The DEG table is keyed on gene id; multiple ids collapse to one symbol
    # (451 duplicated symbols). Collapse to one row per symbol keeping the most
    # significant transcript (lowest treat_fdr), matching the "a symbol is a DEG
    # if any of its transcripts qualifies" set-membership semantics used for the
    # gene index. keep="first" would arbitrarily drop ~7 DEGs.
    sort_key = "treat_fdr" if "treat_fdr" in present else present[0]
    gate = (
        deg[["symbol"] + present]
        .dropna(subset=["symbol"])
        .sort_values(sort_key, kind="stable", na_position="last")
        .drop_duplicates(subset="symbol", keep="first")
        .rename(columns={"symbol": "human_symbol"})
    )
    merged = atlas.merge(gate, on="human_symbol", how="left")
    n_hit = int(merged["treat_fdr"].notna().sum()) if "treat_fdr" in merged.columns else 0
    print(f"  Merged gate columns {present}: {n_hit} atlas rows matched a canonical DEG")
    return merged


def compute_evidence_strengths(atlas: pd.DataFrame) -> pd.DataFrame:
    """Compute 7 evidence source strengths (0-1) for each gene."""

    evidence = pd.DataFrame(index=atlas.index)

    # S1 Human Bulk: |bulk_tstat| normalized to 99th percentile
    if "bulk_tstat" in atlas.columns:
        evidence["s1_human"] = normalize_to_percentile(
            atlas["bulk_tstat"].astype(float), 99.0
        )
    else:
        evidence["s1_human"] = 0.0

    # S2 Genetic Causal: max of all *_coloc_pp4 columns (already 0-1)
    pp4_cols = [c for c in atlas.columns if c.endswith("_coloc_pp4")]
    if pp4_cols:
        pp4_df = atlas[pp4_cols].apply(pd.to_numeric, errors="coerce")
        evidence["s2_genetic"] = pp4_df.max(axis=1).fillna(0.0).clip(0.0, 1.0)
    else:
        evidence["s2_genetic"] = 0.0

    # S3 Essentiality: |essentiality_chronos| normalized to 99th percentile
    if "essentiality_chronos" in atlas.columns:
        evidence["s3_essential"] = normalize_to_percentile(
            atlas["essentiality_chronos"].astype(float), 99.0
        )
    else:
        evidence["s3_essential"] = 0.0

    # S4 Epigenomic: binary from is_disease_regulon_target (if exists, else 0)
    if "is_disease_regulon_target" in atlas.columns:
        evidence["s4_epigenomic"] = coerce_bool(atlas["is_disease_regulon_target"]).astype(float)
    else:
        evidence["s4_epigenomic"] = 0.0

    # S5 Spatial is retained for one compatibility release only. Spatial
    # evidence is categorical and dataset-qualified in the new release tables;
    # it must not contribute a numeric strength or an active-layer vote.
    evidence["s5_spatial"] = 0.0

    # S6 Single-Cell: binary from sceqtl_n_cell_types > 0
    if "sceqtl_n_cell_types" in atlas.columns:
        sceqtl = pd.to_numeric(atlas["sceqtl_n_cell_types"], errors="coerce").fillna(0)
        evidence["s6_singlecell"] = (sceqtl > 0).astype(float)
    else:
        evidence["s6_singlecell"] = 0.0

    # S7 Mouse: translatability_score (already 0-1)
    if "translatability_score" in atlas.columns:
        evidence["s7_mouse"] = (
            pd.to_numeric(atlas["translatability_score"], errors="coerce")
            .fillna(0.0)
            .clip(0.0, 1.0)
        )
    else:
        evidence["s7_mouse"] = 0.0

    return evidence


def build_atlas_parquet(atlas: pd.DataFrame, output_dir: Path):
    """Write the full atlas as a parquet file."""
    out_path = output_dir / "atlas.parquet"
    print(f"Writing {out_path} ...")

    # Replace inf values with NaN before writing
    numeric_cols = atlas.select_dtypes(include=[np.number]).columns
    atlas[numeric_cols] = atlas[numeric_cols].replace([np.inf, -np.inf], np.nan)

    atlas.to_parquet(out_path, index=False, engine="pyarrow")
    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(f"  atlas.parquet: {size_mb:.1f} MB, {atlas.shape[0]} rows x {atlas.shape[1]} cols")


# ---------------------------------------------------------------------------
# atlas_core.parquet (contract §2) — the client's primary table
# ---------------------------------------------------------------------------

# Output column -> atlas source column. Where the spec name differs from the
# actual atlas column, the mapping is applied on output (verified 2026-07-08:
# spec `coloc_best_susie_pp4`/`coloc_best_susie_gwas` are transpositions of the
# real atlas columns `coloc_susie_best_pp4`/`coloc_susie_best_gwas` — the
# inclusive-max SuSiE columns, most populated, THRB=0.9999 / HKDC1=0.9915).
ATLAS_CORE_FROM_ATLAS = {
    "human_symbol": "human_symbol",
    "ensembl_id": "ensembl_id",
    "gene_biotype": "gene_biotype",
    "mouse_ortholog": "mouse_ortholog",
    "bulk_logFC": "bulk_logFC",
    "bulk_padj": "bulk_padj",
    "bulk_tstat": "bulk_tstat",
    "bulk_shrunk_logFC": "bulk_shrunk_logFC",
    "bulk_lfsr": "bulk_lfsr",
    "bulk_sig": "bulk_sig",
    "treat_lfc": "treat_lfc",   # already joined onto atlas by merge_gate_columns()
    "treat_p": "treat_p",
    "treat_fdr": "treat_fdr",
    "bulk_logFC_M": "bulk_logFC_M",
    "bulk_logFC_F": "bulk_logFC_F",
    "sex_class": "sex_class",
    "sex_interaction_padj": "sex_interaction_padj",
    "coloc_best_susie_pp4": "coloc_susie_best_pp4",   # mapped
    "coloc_best_susie_gwas": "coloc_susie_best_gwas",  # mapped
    "coloc_abf_best_pp4": "coloc_abf_best_pp4",
    "coloc_abf_best_gwas": "coloc_abf_best_gwas",
    "n_coloc_sources": "n_coloc_sources",
    "n_ancestry_gwas": "n_ancestry_gwas",
    "coloc_cross_ancestry_replicated": "coloc_cross_ancestry_replicated",
    "coloc_susie_conf_tier": "coloc_susie_conf_tier",
    "twas_z": "twas_z",
    "twas_pval": "twas_pval",
    "essentiality_chronos": "essentiality_chronos",
    "is_essential": "is_essential",
    "spatial_is_svg": "spatial_is_svg",
    "spatial_morans_i": "spatial_morans_i",
    "spatial_consensus_direction": "spatial_consensus_direction",
    "zonation_class": "zonation_class",
    "ferroptosis_class": "ferroptosis_class",
    "dgidb_druggable": "dgidb_druggable",
    "opentargets_drug": "opentargets_drug",
    "max_phase_masld": "max_phase_masld",
    "drug_dev_status": "drug_dev_status",
    "pharos_tdl": "pharos_tdl",
    "dominant_program_for_gene": "dominant_program_for_gene",
    "dominant_program_logFC": "dominant_program_logFC",
    "layers_active": "layers_active",
    "is_conserved": "is_conserved",
}

# Columns to coerce to a clean boolean dtype in atlas_core.
ATLAS_CORE_BOOL_COLS = [
    "bulk_sig", "is_deg", "coloc_cross_ancestry_replicated", "is_essential",
    "spatial_is_svg", "dgidb_druggable", "opentargets_drug", "is_conserved",
]


def build_atlas_core(atlas: pd.DataFrame, convergence: pd.DataFrame, output_dir: Path):
    """Build atlas_core.parquet (contract §2): a compact ~48-col client table =
    subset of atlas + joined treat_* (already on atlas) + derived is_deg +
    joined convergence fields (rank/score/tier/concordance_state, excluding
    excluded_from_ranking rows). Keyed by human_symbol (unique in the atlas).
    """
    out_path = output_dir / "atlas_core.parquet"
    print(f"Building {out_path.name} ...")

    # Select + rename atlas-sourced columns (missing sources -> NaN column).
    core = pd.DataFrame(index=atlas.index)
    for out_col, src_col in ATLAS_CORE_FROM_ATLAS.items():
        if src_col in atlas.columns:
            core[out_col] = atlas[src_col].values
        else:
            print(f"  [warn] atlas missing '{src_col}' for core col '{out_col}' -> NaN")
            core[out_col] = np.nan

    # Derived DEG flag = canonical interval-null gate on the joined treat_fdr.
    treat_fdr = pd.to_numeric(core.get("treat_fdr"), errors="coerce")
    core["is_deg"] = (treat_fdr < 0.05).fillna(False)
    n_deg = int(core["is_deg"].sum())
    print(f"  is_deg (treat_fdr<0.05) = {n_deg} unique symbols "
          f"(canonical DEG *rows* = 1,918; 3 dup-symbol collisions)")

    # Join convergence fields (exclude excluded_from_ranking rows before joining).
    conv_cols_out = {
        "convergence_rank": "convergence_rank",
        "convergence_score": "convergence_score",
        "tier": "convergence_tier",
        "concordance_state": "concordance_state",
    }
    have = [c for c in conv_cols_out if c in convergence.columns]
    if "human_symbol" in convergence.columns and have:
        conv = convergence.copy()
        if "excluded_from_ranking" in conv.columns:
            conv = conv[~coerce_bool(conv["excluded_from_ranking"])]
        conv = (
            conv[["human_symbol"] + have]
            .dropna(subset=["human_symbol"])
            .drop_duplicates(subset="human_symbol", keep="first")
            .rename(columns=conv_cols_out)
        )
        core = core.merge(conv, on="human_symbol", how="left")
        if "convergence_tier" in core.columns:
            t1 = int((core["convergence_tier"] == "1_Genetic_validated").sum())
            print(f"  convergence_tier == 1_Genetic_validated = {t1} (target 677)")
    else:
        print("  [warn] convergence join skipped (missing key/columns)")
        for c in conv_cols_out.values():
            core[c] = np.nan

    # Tier-1/2 liver-specific COLOC flags (canonical headline effector set).
    # The tier12 file is keyed by ensembl only; to maximize recovery onto the
    # symbol-keyed atlas we match a gene by atlas ensembl_id OR by mapped symbol
    # (ensembl->symbol via the main gene_level_coloc.csv). ~24/56 tier12 effector
    # genes are coloc-only (not in the 27,187-gene atlas universe), so the atlas
    # ceiling is ~449 SuSiE / ~975 union, not the file totals 473 / 1,031.
    #   coloc_tier12_pass  = tier12 SuSiE PP.H4 > 0.5             (canonical SuSiE set)
    #   coloc_tier12_union = tier12 any_main (susie|abf PP.H4>0.5) (union set)
    # Two flags because one boolean can reproduce only one of the two counts.
    core["coloc_tier12_pass"] = False
    core["coloc_tier12_union"] = False
    if TIER12_COLOC_CSV.exists() and "ensembl_id" in core.columns:
        t12 = pd.read_csv(TIER12_COLOC_CSV)
        if "ensembl" in t12.columns:
            t12 = t12.assign(_ens=t12["ensembl"].astype(str).str.split(".").str[0])
            ens2sym = {}
            if COLOC_GENE_CSV.exists():
                m = pd.read_csv(COLOC_GENE_CSV, usecols=["gene", "ensembl"])
                ens2sym = dict(
                    zip(m["ensembl"].astype(str).str.split(".").str[0], m["gene"])
                )

            def _ens_sym_sets(mask):
                sub = t12.loc[mask, "_ens"].dropna()
                ens = set(sub)
                syms = set(pd.Series(list(ens)).map(ens2sym).dropna())
                return ens, syms

            su = pd.to_numeric(t12.get("coloc_best_susie_pp4"), errors="coerce")
            su_ens, su_sym = _ens_sym_sets(su > 0.5)
            if "any_main" in t12.columns:
                un_ens, un_sym = _ens_sym_sets(t12["any_main"] == True)
            else:
                un_ens, un_sym = su_ens, su_sym
            core_ens = core["ensembl_id"].astype(str).str.split(".").str[0]
            core_sym = core["human_symbol"].astype(str)
            core["coloc_tier12_pass"] = core_ens.isin(su_ens) | core_sym.isin(su_sym)
            core["coloc_tier12_union"] = core_ens.isin(un_ens) | core_sym.isin(un_sym)
            print(f"  coloc_tier12_pass  (tier12 SuSiE>0.5, atlas-matched) = {int(core['coloc_tier12_pass'].sum())} (file total 473)")
            print(f"  coloc_tier12_union (tier12 any_main, atlas-matched)  = {int(core['coloc_tier12_union'].sum())} (file total 1031)")
        else:
            print(f"  [warn] tier12 file has no 'ensembl' column; tier12 flags left False")
    else:
        print(f"  [warn] tier12 file not found ({TIER12_COLOC_CSV}); tier12 flags left False")

    # Re-derive drug_dev_status from the CANONICAL drug-target classification TSV
    # (joined by symbol), overriding the atlas-sourced value which used a stale
    # classification (e.g. the atlas marked SLC5A2 masld_approved; the canonical
    # TSV does not — the 2 canonical approved targets are THRB + GLP1R).
    # NOTE: the atlas-restricted GROUP BY = 1 approved / 186 clinical / 1,345
    # preclinical, NOT the full-universe headline 2 / 208 / 1,504 — GLP1R (approved)
    # and ~22 clinical / ~159 preclinical drug targets are not in the 27,187-gene
    # atlas. The full-universe headline must come from the TSV, not this column.
    if DRUG_CLASS_TSV.exists() and "human_symbol" in core.columns:
        dt = pd.read_csv(
            DRUG_CLASS_TSV, sep="\t",
            usecols=["symbol", "drug_dev_status", "max_phase_masld"],
            low_memory=False,
        ).dropna(subset=["symbol"])
        _rank = {
            "masld_approved": 0, "masld_clinical": 1, "masld_preclinical": 2,
            "masld_discontinued": 3, "drugged_other_indication": 4,
            "discovery": 5, "undetermined": 6,
        }
        dt = dt.assign(_r=dt["drug_dev_status"].map(_rank).fillna(9))
        dt = dt.sort_values("_r").drop_duplicates("symbol", keep="first")
        # Verbatim copy of the TSV's canonical drug_dev_status + max_phase_masld.
        # Every atlas symbol is in the TSV (a 27,943-gene superset; 0 unmatched),
        # so the drug_dev_status fallback is defensive only — 'undetermined'
        # rather than the non-canonical atlas value. max_phase_masld is NaN where
        # the gene has no MASLD-indication phase (canonical).
        status_map = dict(zip(dt["symbol"], dt["drug_dev_status"]))
        phase_map = dict(zip(dt["symbol"], dt["max_phase_masld"]))
        core["drug_dev_status"] = (
            core["human_symbol"].map(status_map).fillna("undetermined")
        )
        if "max_phase_masld" in core.columns:
            core["max_phase_masld"] = core["human_symbol"].map(phase_map)
        vc = core["drug_dev_status"].value_counts()
        print(
            f"  drug_dev_status (canonical TSV, atlas-restricted): "
            f"approved={int(vc.get('masld_approved',0))} "
            f"clinical={int(vc.get('masld_clinical',0))} "
            f"preclinical={int(vc.get('masld_preclinical',0))} "
            f"(full-universe headline = 2 / 208 / 1504)"
        )
    else:
        print(f"  [warn] drug class TSV not found ({DRUG_CLASS_TSV}); drug_dev_status unchanged")

    # Clean dtypes: bools, and ±Inf -> NaN on numerics.
    for c in ATLAS_CORE_BOOL_COLS:
        if c in core.columns and c != "is_deg":
            core[c] = coerce_bool(core[c])
    num_cols = core.select_dtypes(include=[np.number]).columns
    core[num_cols] = core[num_cols].replace([np.inf, -np.inf], np.nan)

    # Order columns per contract §2 (atlas-sourced order, with is_deg after
    # treat_fdr and convergence block after coloc_susie_conf_tier region).
    ordered = [
        "human_symbol", "ensembl_id", "gene_biotype", "mouse_ortholog",
        "bulk_logFC", "bulk_padj", "bulk_tstat", "bulk_shrunk_logFC", "bulk_lfsr",
        "bulk_sig", "treat_lfc", "treat_p", "treat_fdr", "is_deg",
        "bulk_logFC_M", "bulk_logFC_F", "sex_class", "sex_interaction_padj",
        "coloc_best_susie_pp4", "coloc_best_susie_gwas", "coloc_abf_best_pp4",
        "coloc_abf_best_gwas", "n_coloc_sources", "n_ancestry_gwas",
        "coloc_cross_ancestry_replicated", "coloc_susie_conf_tier",
        "coloc_tier12_pass", "coloc_tier12_union",
        "twas_z", "twas_pval",
        "convergence_rank", "convergence_score", "convergence_tier",
        "concordance_state",
        "essentiality_chronos", "is_essential",
        "spatial_is_svg", "spatial_morans_i", "spatial_consensus_direction",
        "zonation_class", "ferroptosis_class",
        "dgidb_druggable", "opentargets_drug", "max_phase_masld",
        "drug_dev_status", "pharos_tdl",
        "dominant_program_for_gene", "dominant_program_logFC",
        "layers_active", "is_conserved",
    ]
    ordered = [c for c in ordered if c in core.columns]
    core = core[ordered]

    core.to_parquet(out_path, index=False, engine="pyarrow")
    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(f"  atlas_core.parquet: {size_mb:.2f} MB, {core.shape[0]} rows x {core.shape[1]} cols")
    return core


def build_gene_symbols_json(atlas: pd.DataFrame, output_dir: Path):
    """Emit gene_symbols.json: slim [{symbol, biotype}] for fuse.js autocomplete
    (replaces the 10.6 MB gene_index.json dependency on the gene page)."""
    out_path = output_dir / "gene_symbols.json"
    print(f"Building {out_path.name} ...")
    df = atlas[["human_symbol", "gene_biotype"]].copy()
    df = df[df["human_symbol"].notna() & (df["human_symbol"].astype(str).str.strip() != "")]
    records = [
        {"symbol": str(s).strip(), "biotype": safe_json_value(b)}
        for s, b in zip(df["human_symbol"], df["gene_biotype"])
    ]
    records.sort(key=lambda r: r["symbol"])
    with open(out_path, "w") as f:
        json.dump(records, f, separators=(",", ":"), default=safe_json_value)
    size_kb = out_path.stat().st_size / 1024
    print(f"  gene_symbols.json: {size_kb:.1f} KB, {len(records)} symbols")
    return records


def build_gene_index(
    atlas: pd.DataFrame, evidence: pd.DataFrame, deg: pd.DataFrame, output_dir: Path
):
    """Build gene_index.json: compact search index with evidence strengths.

    Optimizations for file size:
    - Null optional fields (sex_class, zonation_class, ferroptosis_class) are omitted
    - Boolean false values for is_conserved and dgidb_druggable are omitted
    - Evidence values of 0.0 are omitted from the evidence dict
    - bulk_logfc/bulk_padj rounded to 3/4 decimal places
    - Evidence rounded to 2 decimal places
    """
    out_path = output_dir / "gene_index.json"
    print(f"Building gene index ...")

    # Canonical DEG status = effect-size-aware interval-null FDR gate
    # (treat_fdr < 0.05). Symbol-level set membership → 1,915 unique symbols.
    if "treat_fdr" in deg.columns:
        deg_mask = deg["treat_fdr"] < 0.05
    else:
        deg_mask = (deg["padj"] < 0.05) & (deg["logFC"].abs() > 0.3)
    canonical_degs = set(deg.loc[deg_mask, "symbol"].dropna())

    records = []
    for idx, row in atlas.iterrows():
        symbol = row.get("human_symbol")
        if pd.isna(symbol) or str(symbol).strip() == "":
            continue

        symbol = str(symbol).strip()
        ensembl_id = safe_json_value(row.get("ensembl_id"))
        biotype = safe_json_value(row.get("gene_biotype"))
        bulk_logfc = safe_json_value(row.get("bulk_logFC"))
        bulk_padj = safe_json_value(row.get("bulk_padj"))
        is_deg = symbol in canonical_degs

        # Boolean / categorical fields
        is_conserved_raw = row.get("is_conserved")
        is_conserved = bool(
            str(is_conserved_raw).strip().upper() in ("TRUE", "1", "1.0")
        )

        sex_class_raw = row.get("sex_class")
        sex_class = (
            safe_json_value(sex_class_raw)
            if pd.notna(sex_class_raw)
            and str(sex_class_raw).strip().lower()
            not in ("", "nan", "not_significant")
            else None
        )

        zonation_class_raw = row.get("zonation_class")
        zonation_class = (
            safe_json_value(zonation_class_raw)
            if pd.notna(zonation_class_raw)
            and str(zonation_class_raw).strip().lower() not in ("", "nan")
            else None
        )

        ferroptosis_class_raw = row.get("ferroptosis_class")
        ferroptosis_class = (
            safe_json_value(ferroptosis_class_raw)
            if pd.notna(ferroptosis_class_raw)
            and str(ferroptosis_class_raw).strip().lower() not in ("", "nan")
            else None
        )

        dgidb_raw = row.get("dgidb_druggable")
        dgidb_druggable = bool(
            str(dgidb_raw).strip().upper() in ("TRUE", "1", "1.0")
        )

        layers_active_raw = row.get("layers_active")
        layers_active = int(
            pd.to_numeric(layers_active_raw, errors="coerce")
            if pd.notna(layers_active_raw)
            else 0
        )

        # Evidence strengths: omit zero-valued entries to save space
        ev = evidence.loc[idx]
        evidence_dict = {}
        for key in ("s1_human", "s2_genetic", "s3_essential", "s4_epigenomic",
                     "s5_spatial", "s6_singlecell", "s7_mouse"):
            val = round(float(ev[key]), 2)
            if val > 0:
                evidence_dict[key] = val

        # Build record, omitting null/false optional fields for compactness
        record = {
            "symbol": symbol,
            "ensembl_id": ensembl_id,
            "biotype": biotype,
            "bulk_logfc": round(float(bulk_logfc), 3) if bulk_logfc is not None else None,
            "bulk_padj": round(float(bulk_padj), 4) if bulk_padj is not None else None,
            "is_deg": is_deg,
            "is_conserved": is_conserved,
            "sex_class": sex_class,
            "zonation_class": zonation_class,
            "ferroptosis_class": ferroptosis_class,
            "dgidb_druggable": dgidb_druggable,
            "layers_active": layers_active,
            "evidence": evidence_dict,
        }

        # Strip null optional fields and false booleans to reduce size
        if sex_class is None:
            del record["sex_class"]
        if zonation_class is None:
            del record["zonation_class"]
        if ferroptosis_class is None:
            del record["ferroptosis_class"]
        if not is_conserved:
            del record["is_conserved"]
        if not dgidb_druggable:
            del record["dgidb_druggable"]
        if not is_deg:
            del record["is_deg"]
        if layers_active == 0:
            del record["layers_active"]

        records.append(record)

    # Sort by symbol
    records.sort(key=lambda r: r["symbol"])

    print(f"  Writing {len(records)} gene entries to {out_path} ...")
    with open(out_path, "w") as f:
        json.dump(records, f, separators=(",", ":"), default=safe_json_value)

    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(f"  gene_index.json: {size_mb:.2f} MB")
    return records


def build_atlas_summary(
    atlas: pd.DataFrame,
    deg: pd.DataFrame,
    coloc: pd.DataFrame,
    drugs: pd.DataFrame,
    output_dir: Path,
):
    """Build atlas_summary.json with landing page statistics."""
    out_path = output_dir / "atlas_summary.json"
    print(f"Building atlas summary ...")

    total_degs = compute_deg_count(deg)
    n_deg_symbols = compute_deg_symbol_count(deg)
    coloc_genes = int((coloc["coloc_best_pp4"] > 0.5).sum())
    drug_targets = int(drugs.shape[0])

    # Conserved core count
    if "is_conserved" in atlas.columns:
        cc_count = int(coerce_bool(atlas["is_conserved"]).sum())
    else:
        cc_count = 0

    # total_degs   = canonical interval-null FDR gate, transcript/gene-model
    #                level = 1,918 (the headline).
    # n_deg_symbols = same gate collapsed to unique symbols = 1,915 (what the
    #                symbol-keyed atlas / atlas_core is_deg flag reports). Both
    #                are exposed so the two levels are explicit and honest.
    # total_cohorts/total_samples = the 5 control-bearing cohorts / 846 samples
    # in the canonical Disease-vs-Control pooled analysis (CLAUDE.md), i.e. the
    # basis for the headline DEG number (not the full 9-cohort / 1,277-sample atlas).
    summary = {
        "total_genes": int(atlas.shape[0]),
        "total_degs": total_degs,
        "n_deg_symbols": n_deg_symbols,
        "total_cohorts": 5,
        "total_samples": 846,
        "conserved_count": cc_count,
        "coloc_genes": coloc_genes,
        "drug_targets": drug_targets,
        "mouse_datasets": 5,
        "evidence_sources": 7,
    }

    with open(out_path, "w") as f:
        json.dump(summary, f, separators=(",", ":"))

    size_kb = out_path.stat().st_size / 1024
    print(f"  atlas_summary.json: {size_kb:.1f} KB")
    print(f"  Stats: {summary}")
    return summary


def build_featured_genes(
    atlas: pd.DataFrame,
    evidence: pd.DataFrame,
    drugs: pd.DataFrame,
    gene_index: list,
    output_dir: Path,
):
    """Build featured_genes.json for 6 highlighted gene cards."""
    out_path = output_dir / "featured_genes.json"
    print(f"Building featured genes ...")

    # Index gene_index by symbol for quick lookup
    gi_lookup = {r["symbol"]: r for r in gene_index}

    # Index drug data by target_gene
    drug_lookup = {}
    for _, drow in drugs.iterrows():
        tg = str(drow.get("target_gene", "")).strip()
        if tg and tg != "nan":
            if tg not in drug_lookup:
                drug_lookup[tg] = []
            drug_lookup[tg].append(
                {
                    "drug": safe_json_value(drow.get("drug")),
                    "stage": safe_json_value(drow.get("stage")),
                    "moa": safe_json_value(drow.get("moa")),
                    "atlas_support": safe_json_value(drow.get("atlas_support")),
                }
            )

    featured = []
    for symbol, meta in FEATURED_GENES.items():
        gi = gi_lookup.get(symbol)
        if gi is None:
            print(f"  WARNING: Featured gene {symbol} not found in atlas, skipping")
            continue

        # Featured gene cards always include all fields (with defaults)
        # for the web UI, even if omitted in compact gene_index
        ev_full = {}
        ev_compact = gi.get("evidence", {})
        for k in ("s1_human", "s2_genetic", "s3_essential", "s4_epigenomic",
                   "s5_spatial", "s6_singlecell", "s7_mouse"):
            ev_full[k] = ev_compact.get(k, 0.0)

        card = {
            "symbol": symbol,
            "ensembl_id": gi.get("ensembl_id"),
            "tagline": meta["tagline"],
            "category": meta["category"],
            "bulk_logfc": gi.get("bulk_logfc"),
            "bulk_padj": gi.get("bulk_padj"),
            "is_deg": gi.get("is_deg", False),
            "is_conserved": gi.get("is_conserved", False),
            "layers_active": gi.get("layers_active", 0),
            "evidence": ev_full,
            "drugs": drug_lookup.get(symbol, []),
        }
        featured.append(card)

    with open(out_path, "w") as f:
        json.dump(featured, f, separators=(",", ":"), indent=None, default=safe_json_value)

    size_kb = out_path.stat().st_size / 1024
    print(f"  featured_genes.json: {size_kb:.1f} KB, {len(featured)} genes")
    return featured


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main():
    if os.environ.get("ALLOW_LEGACY_PORTAL_REBUILD", "false").lower() != "true":
        raise SystemExit(
            "REFUSED: legacy atlas preprocessing is outside the standalone "
            "Resource contract. Use Plans 50/60."
        )
    parser = argparse.ArgumentParser(
        description="Preprocess MASLD atlas data for the web application."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=os.path.join(os.path.dirname(__file__), "..", "public", "data"),
        help="Output directory for generated files (default: ../public/data)",
    )
    parser.add_argument(
        "--only-core",
        action="store_true",
        help="Rebuild ONLY atlas_core.parquet (skip atlas.parquet / gene_index / "
        "gene_symbols / summary / featured). Use for a minimal atlas_core refresh "
        "that must NOT clobber gene_index.json (post-processed by "
        "generate_convergence_data.py).",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Load data
    atlas, deg, coloc, drugs, convergence = load_data()

    # Merge the numeric per-gene DEG-gate columns onto the atlas.
    atlas = merge_gate_columns(atlas, deg)

    if args.only_core:
        print("--only-core: rebuilding atlas_core.parquet only ...")
        build_atlas_core(atlas, convergence, output_dir)
        print("\nDone (only-core). atlas_core.parquet written to:", output_dir)
        return

    # NOTE: dream_robustness_flag is intentionally kept as-is on atlas.parquet.
    # It is a legitimate dream-sensitivity-arm provenance flag (27a-owned; it
    # names the retired method it benchmarks against), NOT a stale canonical-DEG
    # column, so it is not renamed. It is excluded from atlas_core.parquet and is
    # never surfaced as a UI label.

    # Compute evidence strengths
    print("Computing evidence strengths ...")
    evidence = compute_evidence_strengths(atlas)
    atlas["layers_active"] = (
        evidence[[
            "s1_human", "s2_genetic", "s3_essential", "s4_epigenomic",
            "s6_singlecell", "s7_mouse",
        ]] > 0
    ).sum(axis=1)

    # 1. atlas.parquet (full, /downloads)
    build_atlas_parquet(atlas, output_dir)

    # 2. atlas_core.parquet (contract §2, client primary table)
    build_atlas_core(atlas, convergence, output_dir)

    # 3. gene_index.json (legacy; retired from gene page, kept for search index)
    gene_index = build_gene_index(atlas, evidence, deg, output_dir)

    # 4. gene_symbols.json (slim fuse.js autocomplete)
    build_gene_symbols_json(atlas, output_dir)

    # 5. atlas_summary.json
    build_atlas_summary(atlas, deg, coloc, drugs, output_dir)

    # 6. featured_genes.json
    build_featured_genes(atlas, evidence, drugs, gene_index, output_dir)

    print("\nDone. All output files written to:", output_dir)


if __name__ == "__main__":
    main()
