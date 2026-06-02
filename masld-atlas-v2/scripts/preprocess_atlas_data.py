#!/usr/bin/env python3
"""
preprocess_atlas_data.py
========================
Reads source CSVs from the MASLD project and produces 4 output files
for the MASLD Atlas v2 web application:

  1. atlas.parquet      -- Full 33,943-gene atlas for DuckDB-WASM queries
  2. gene_index.json    -- Compact search index with 7-axis evidence per gene
  3. atlas_summary.json -- Landing page statistics
  4. featured_genes.json -- 6 featured gene cards

Usage:
  python preprocess_atlas_data.py --output-dir ../public/data
"""

import argparse
import json
import math
import os
import sys
from pathlib import Path

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
DREAM_ASHR_CSV = (
    PROJECT_ROOT
    / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv"
)
COLOC_GENE_CSV = (
    PROJECT_ROOT / "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"
)
CLINICAL_DRUGS_CSV = (
    PROJECT_ROOT
    / "RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv"
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

    print(f"Loading dream_ashr from {DREAM_ASHR_CSV} ...")
    dream = pd.read_csv(DREAM_ASHR_CSV)
    print(f"  Dream shape: {dream.shape}")

    print(f"Loading COLOC gene-level from {COLOC_GENE_CSV} ...")
    coloc = pd.read_csv(COLOC_GENE_CSV)
    print(f"  Coloc shape: {coloc.shape}")

    print(f"Loading clinical drugs from {CLINICAL_DRUGS_CSV} ...")
    drugs = pd.read_csv(CLINICAL_DRUGS_CSV, low_memory=False)
    print(f"  Drugs shape: {drugs.shape}")

    return atlas, dream, coloc, drugs


def compute_deg_count(dream: pd.DataFrame) -> int:
    """Compute DEG count from dream: padj < 0.05 AND |logFC| > 0.3."""
    return int(
        ((dream["padj"] < 0.05) & (dream["logFC"].abs() > 0.3)).sum()
    )


def compute_evidence_strengths(atlas: pd.DataFrame) -> pd.DataFrame:
    """Compute 7 evidence source strengths (0-1) for each gene."""

    evidence = pd.DataFrame(index=atlas.index)

    # S1 Human Bulk: |dream_tstat| normalized to 99th percentile
    if "dream_tstat" in atlas.columns:
        evidence["s1_human"] = normalize_to_percentile(
            atlas["dream_tstat"].astype(float), 99.0
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

    # S5 Spatial: binary from zonation_class being non-null
    if "zonation_class" in atlas.columns:
        evidence["s5_spatial"] = (
            atlas["zonation_class"].notna()
            & (atlas["zonation_class"].astype(str).str.strip() != "")
            & (atlas["zonation_class"].astype(str).str.lower() != "nan")
        ).astype(float)
    else:
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


def build_gene_index(
    atlas: pd.DataFrame, evidence: pd.DataFrame, dream: pd.DataFrame, output_dir: Path
):
    """Build gene_index.json: compact search index with evidence strengths.

    Optimizations for file size:
    - Null optional fields (sex_class, zonation_class, ferroptosis_class) are omitted
    - Boolean false values for is_conserved and dgidb_druggable are omitted
    - Evidence values of 0.0 are omitted from the evidence dict
    - dream_logfc/dream_padj rounded to 3/4 decimal places
    - Evidence rounded to 2 decimal places
    """
    out_path = output_dir / "gene_index.json"
    print(f"Building gene index ...")

    # Merge dream DEG status into atlas
    dream_degs = set(
        dream.loc[
            (dream["padj"] < 0.05) & (dream["logFC"].abs() > 0.3), "symbol"
        ].dropna()
    )

    records = []
    for idx, row in atlas.iterrows():
        symbol = row.get("human_symbol")
        if pd.isna(symbol) or str(symbol).strip() == "":
            continue

        symbol = str(symbol).strip()
        ensembl_id = safe_json_value(row.get("ensembl_id"))
        biotype = safe_json_value(row.get("gene_biotype"))
        dream_logfc = safe_json_value(row.get("dream_logFC"))
        dream_padj = safe_json_value(row.get("dream_padj"))
        is_deg = symbol in dream_degs

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
            "dream_logfc": round(float(dream_logfc), 3) if dream_logfc is not None else None,
            "dream_padj": round(float(dream_padj), 4) if dream_padj is not None else None,
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
    dream: pd.DataFrame,
    coloc: pd.DataFrame,
    drugs: pd.DataFrame,
    output_dir: Path,
):
    """Build atlas_summary.json with landing page statistics."""
    out_path = output_dir / "atlas_summary.json"
    print(f"Building atlas summary ...")

    total_degs = compute_deg_count(dream)
    coloc_genes = int((coloc["coloc_best_pp4"] > 0.5).sum())
    drug_targets = int(drugs.shape[0])

    # Conserved core count
    if "is_conserved" in atlas.columns:
        cc_count = int(coerce_bool(atlas["is_conserved"]).sum())
    else:
        cc_count = 0

    summary = {
        "total_genes": int(atlas.shape[0]),
        "total_degs": total_degs,
        "total_cohorts": 10,
        "total_samples": 1444,
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
            "dream_logfc": gi.get("dream_logfc"),
            "dream_padj": gi.get("dream_padj"),
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
    parser = argparse.ArgumentParser(
        description="Preprocess MASLD atlas data for the web application."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=os.path.join(os.path.dirname(__file__), "..", "public", "data"),
        help="Output directory for generated files (default: ../public/data)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Load data
    atlas, dream, coloc, drugs = load_data()

    # Compute evidence strengths
    print("Computing evidence strengths ...")
    evidence = compute_evidence_strengths(atlas)

    # 1. atlas.parquet
    build_atlas_parquet(atlas, output_dir)

    # 2. gene_index.json
    gene_index = build_gene_index(atlas, evidence, dream, output_dir)

    # 3. atlas_summary.json
    build_atlas_summary(atlas, dream, coloc, drugs, output_dir)

    # 4. featured_genes.json
    build_featured_genes(atlas, evidence, drugs, gene_index, output_dir)

    print("\nDone. All 4 output files written to:", output_dir)


if __name__ == "__main__":
    main()
