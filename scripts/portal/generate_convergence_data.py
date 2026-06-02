#!/usr/bin/env python3
"""Generate portal data files for the Convergence tab (Phase 7B items 7-8).

Produces:
  1. masld-atlas-v2/public/data/convergence_matrix.json
     - per-gene 7-modality vector (0-1 per modality) + flags
  2. masld-atlas-v2/public/data/bayesian_ranking.json
     - legacy 46b archetype-similarity ranking (from bayesian_posterior.csv, or atlas fallback)
     - NOTE: canonical 46d convergence ranking is now convergence_ranking.json
       (written by generate_convergence_evidence_json.py)
  3. Updates masld-atlas-v2/public/data/gene_index.json IN PLACE
     - renames s3_essential -> s4_essential
             s4_epigenomic -> s5_epigenomic
             s5_spatial    -> s6_spatial
             s6_singlecell -> s7_singlecell
             s2_genetic    -> s3_genetic
             s7_mouse      -> s2_mouse
     - adds missing s5_epigenomic field (atlas-derived)

Modality numbering (applied consistently across this script and gene_index):
  M1 s1_human       = Human bulk DE (|dream_logFC| / max, gated by is_deg)
  M2 s2_mouse       = Mouse conserved (mouse_meta_padj < 0.05 -> 1.0)
  M3 s3_genetic     = Genetic causal (coloc_susie_best_pp4, 0-1)
  M4 s4_essential   = DepMap essentiality (clamp(-essentiality_chronos/2, 0, 1))
  M5 s5_epigenomic  = ATAC/regulatory flag (gwas_variant_in_peak OR any peak/motif)
  M6 s6_spatial     = Spatial SVG (clamp(spatial_morans_i, 0, 1); SVG -> 1.0)
  M7 s7_singlecell  = Cross-celltype signal (min(sc_n_celltypes_sig/5, 1) or LIANA)

Run:
  micromamba run -n spatial python scripts/portal/generate_convergence_data.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
ATLAS_CSV = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
BAYES_CSV = PROJECT_ROOT / "RNA-seq/results/multi_evidence/bayesian_posterior.csv"
PORTAL_DATA = PROJECT_ROOT / "masld-atlas-v2/public/data"
CONVERGENCE_OUT = PORTAL_DATA / "convergence_matrix.json"
BAYES_OUT = PORTAL_DATA / "bayesian_ranking.json"
GENE_INDEX = PORTAL_DATA / "gene_index.json"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _round(val, n=4):
    if val is None:
        return 0.0
    if isinstance(val, float) and (math.isnan(val) or math.isinf(val)):
        return 0.0
    return round(float(val), n)


def clean_for_json(obj):
    if isinstance(obj, dict):
        return {k: clean_for_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [clean_for_json(v) for v in obj]
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
    return obj


def coerce_bool(series: pd.Series) -> pd.Series:
    return series.map(
        lambda x: True
        if str(x).strip().upper() in ("TRUE", "1", "1.0")
        else False
    )


def col_or(df: pd.DataFrame, name: str, default=np.nan) -> pd.Series:
    """Return column if present, else constant-default Series with correct index."""
    if name in df.columns:
        return df[name]
    print(f"  [warn] column missing: {name} (using default)")
    return pd.Series(default, index=df.index, dtype=float)


# ---------------------------------------------------------------------------
# Modality computation
# ---------------------------------------------------------------------------

def compute_modalities(atlas: pd.DataFrame) -> pd.DataFrame:
    """Compute 7 modality strengths (0-1) per gene."""
    out = pd.DataFrame(index=atlas.index)

    # M1 Human bulk DE: |dream_logFC| / max_abs_lfc, gated by is_deg
    lfc = pd.to_numeric(col_or(atlas, "dream_logFC"), errors="coerce")
    padj = pd.to_numeric(col_or(atlas, "dream_padj"), errors="coerce")
    is_deg = (padj < 0.05) & (lfc.abs() > 0.3)
    max_abs_lfc = float(lfc.abs().max())
    if not math.isfinite(max_abs_lfc) or max_abs_lfc == 0:
        max_abs_lfc = 1.0
    m1 = (lfc.abs() / max_abs_lfc).clip(0.0, 1.0).fillna(0.0)
    m1 = m1.where(is_deg, 0.0)
    out["s1_human"] = m1

    # M2 Mouse conserved: 1.0 if mouse_meta_padj < 0.05, else 0
    m_padj = pd.to_numeric(col_or(atlas, "mouse_meta_padj"), errors="coerce")
    out["s2_mouse"] = (m_padj < 0.05).astype(float)

    # M3 Genetic causal: coloc_susie_best_pp4
    pp4 = pd.to_numeric(col_or(atlas, "coloc_susie_best_pp4"), errors="coerce")
    out["s3_genetic"] = pp4.fillna(0.0).clip(0.0, 1.0)

    # M4 Essentiality: clamp(-essentiality_chronos / 2, 0, 1)
    chr_ = pd.to_numeric(col_or(atlas, "essentiality_chronos"), errors="coerce")
    out["s4_essential"] = (-chr_ / 2.0).clip(0.0, 1.0).fillna(0.0)

    # M5 Epigenomic: union of ATAC/GWAS-ATAC signals
    #   gwas_variant_in_peak (bool) OR gwas_variant_da_peak OR
    #   gwas_max_pip_in_peak > 0 OR gwas_motif_disrupted non-null
    ep_sources = []
    for c in ("gwas_variant_in_peak", "gwas_variant_da_peak"):
        if c in atlas.columns:
            ep_sources.append(coerce_bool(atlas[c].astype(str)))
    if "gwas_max_pip_in_peak" in atlas.columns:
        ep_sources.append(
            pd.to_numeric(atlas["gwas_max_pip_in_peak"], errors="coerce")
            .fillna(0.0)
            .gt(0.0)
        )
    if "gwas_motif_disrupted" in atlas.columns:
        gmd = atlas["gwas_motif_disrupted"].astype(str).fillna("")
        ep_sources.append(gmd.str.strip().ne("") & gmd.str.lower().ne("nan"))
    if "gwas_atac_regulatory_score" in atlas.columns:
        ep_sources.append(
            pd.to_numeric(atlas["gwas_atac_regulatory_score"], errors="coerce")
            .fillna(0.0)
            .gt(0.0)
        )
    if ep_sources:
        m5 = ep_sources[0].astype(bool)
        for s in ep_sources[1:]:
            m5 = m5 | s.astype(bool)
        out["s5_epigenomic"] = m5.astype(float)
    else:
        print("  [warn] no epigenomic columns found; s5=0")
        out["s5_epigenomic"] = 0.0

    # M6 Spatial: clamp(spatial_morans_i, 0, 1) OR 1.0 if spatial_is_svg
    mi = pd.to_numeric(col_or(atlas, "spatial_morans_i"), errors="coerce")
    svg_flag = pd.Series(False, index=atlas.index)
    if "spatial_is_svg" in atlas.columns:
        svg = atlas["spatial_is_svg"]
        if svg.dtype == bool:
            svg_flag = svg.fillna(False).astype(bool)
        else:
            svg_flag = coerce_bool(svg.astype(str))
    m6_base = mi.clip(0.0, 1.0).fillna(0.0)
    out["s6_spatial"] = np.where(svg_flag, 1.0, m6_base)

    # M7 Single-cell: sc_n_celltypes_sig/5 OR liana_n_diff_interactions > 0 -> 1.0
    sc = pd.to_numeric(col_or(atlas, "sc_n_celltypes_sig"), errors="coerce").fillna(0.0)
    liana = pd.to_numeric(col_or(atlas, "liana_n_diff_interactions"), errors="coerce").fillna(0.0)
    m7_sc = (sc / 5.0).clip(0.0, 1.0)
    m7_liana = (liana > 0).astype(float)
    out["s7_singlecell"] = np.maximum(m7_sc, m7_liana)

    # M8 Proteomics: best_protein_logFC + best_protein_padj
    prot_lfc = pd.to_numeric(col_or(atlas, "best_protein_logFC"), errors="coerce")
    prot_padj = pd.to_numeric(col_or(atlas, "best_protein_padj"), errors="coerce")
    prot_sig = (prot_padj < 0.05)
    m8_sig = (prot_lfc.abs() / 2.0).clip(0.0, 1.0)
    m8_nonsig = (prot_lfc.abs() * 0.3).clip(0.0, 0.15)
    m8 = np.where(prot_sig.fillna(False), m8_sig.fillna(0.0), m8_nonsig.fillna(0.0))
    out["s8_proteomics"] = pd.Series(m8, index=atlas.index).fillna(0.0)

    return out


# ---------------------------------------------------------------------------
# Convergence matrix
# ---------------------------------------------------------------------------

def build_convergence_matrix(atlas: pd.DataFrame, modalities: pd.DataFrame) -> list:
    # Flags
    lfc = pd.to_numeric(col_or(atlas, "dream_logFC"), errors="coerce")
    padj = pd.to_numeric(col_or(atlas, "dream_padj"), errors="coerce")
    is_deg = ((padj < 0.05) & (lfc.abs() > 0.3)).fillna(False)
    pp4 = pd.to_numeric(col_or(atlas, "coloc_susie_best_pp4"), errors="coerce").fillna(0.0)
    is_coloc = (pp4 >= 0.5)
    is_druggable = pd.Series(False, index=atlas.index)
    if "dgidb_druggable" in atlas.columns:
        d = atlas["dgidb_druggable"]
        if d.dtype == bool:
            is_druggable = d.fillna(False)
        else:
            is_druggable = coerce_bool(d.astype(str))

    col_order = ["s1_human", "s2_mouse", "s3_genetic", "s4_essential",
                 "s5_epigenomic", "s6_spatial", "s7_singlecell", "s8_proteomics"]
    mod_vals = modalities[col_order].fillna(0.0).clip(0.0, 1.0).values

    # count = number of modalities with value > 0.1
    count = (mod_vals > 0.1).sum(axis=1)

    rows = []
    symbols = atlas["human_symbol"].values
    for i, sym in enumerate(symbols):
        if not isinstance(sym, str) or sym == "" or sym.lower() == "nan":
            continue
        rows.append({
            "gene": sym,
            "modalities": [round(float(x), 4) for x in mod_vals[i]],
            "count": int(count[i]),
            "is_deg": bool(is_deg.iloc[i]),
            "is_coloc": bool(is_coloc.iloc[i]),
            "is_druggable": bool(is_druggable.iloc[i]),
        })

    # sort by count desc, then by s1_human desc
    rows.sort(key=lambda r: (-r["count"], -r["modalities"][0]))
    return rows


# ---------------------------------------------------------------------------
# Bayesian ranking
# ---------------------------------------------------------------------------

def build_legacy_archetype_ranking(
    atlas: pd.DataFrame,
    modalities: pd.DataFrame,
) -> dict:
    # Try to load bayesian_posterior.csv (canonical Script 46b output)
    score_source = None
    bayes = None
    if BAYES_CSV.exists():
        try:
            bayes = pd.read_csv(BAYES_CSV)
            if {"human_symbol", "posterior_odds"}.issubset(bayes.columns):
                score_source = (
                    f"atlas file: {BAYES_CSV.relative_to(PROJECT_ROOT)} "
                    f"(column: posterior_odds)"
                )
        except Exception as e:
            print(f"  [warn] failed to load bayesian_posterior.csv: {e}")
            bayes = None

    # Flags
    lfc = pd.to_numeric(col_or(atlas, "dream_logFC"), errors="coerce")
    padj = pd.to_numeric(col_or(atlas, "dream_padj"), errors="coerce")
    is_deg = ((padj < 0.05) & (lfc.abs() > 0.3)).fillna(False)
    pp4 = pd.to_numeric(col_or(atlas, "coloc_susie_best_pp4"), errors="coerce").fillna(0.0)
    is_coloc = (pp4 >= 0.5)
    is_druggable = pd.Series(False, index=atlas.index)
    if "dgidb_druggable" in atlas.columns:
        d = atlas["dgidb_druggable"]
        if d.dtype == bool:
            is_druggable = d.fillna(False)
        else:
            is_druggable = coerce_bool(d.astype(str))

    # Build a per-atlas-row score
    if bayes is not None and score_source is not None:
        score_map = dict(zip(bayes["human_symbol"], bayes["posterior_odds"]))
        score_series = atlas["human_symbol"].map(score_map).astype(float)
        # Rows with no posterior: fill with min observed so they sort last
        min_obs = pd.to_numeric(score_series, errors="coerce").min()
        score_series = score_series.fillna(min_obs - 1.0 if pd.notnull(min_obs) else 0.0)
    else:
        # Fallback formula (documented in metadata)
        score_source = (
            "fallback formula: 0.35*is_deg + 0.25*coloc_susie_best_pp4 + "
            "0.15*mouse_conserved + 0.1*spatial_morans_i + 0.15*proteomics"
        )
        mouse_flag = modalities["s2_mouse"]
        spatial = modalities["s6_spatial"]
        proteomics = modalities["s8_proteomics"]
        score_series = (
            0.35 * is_deg.astype(float)
            + 0.25 * pp4.clip(0.0, 1.0)
            + 0.15 * mouse_flag
            + 0.1 * spatial
            + 0.15 * proteomics
        )

    n = len(atlas)
    # Rank (1-based, descending score)
    order = np.argsort(-score_series.values, kind="stable")
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    # percentile: 100 * (1 - (rank-1)/n)
    pcts = 100.0 * (1.0 - (ranks - 1) / n)

    genes = []
    symbols = atlas["human_symbol"].values
    for i in range(n):
        sym = symbols[i]
        if not isinstance(sym, str) or sym == "" or sym.lower() == "nan":
            continue
        score_val = float(score_series.iloc[i]) if pd.notnull(score_series.iloc[i]) else 0.0
        genes.append({
            "rank": int(ranks[i]),
            "symbol": sym,
            "score": round(score_val, 5),
            "percentile": round(float(pcts[i]), 2),
            "is_deg": bool(is_deg.iloc[i]),
            "is_coloc": bool(is_coloc.iloc[i]),
            "is_druggable": bool(is_druggable.iloc[i]),
        })
    genes.sort(key=lambda r: r["rank"])

    return {
        "genes": genes,
        "metadata": {
            "n_genes": len(genes),
            "score_source": score_source,
            "notes": (
                "rank is 1-based, descending by score. percentile = "
                "100 * (n - rank + 1) / n (higher = more evidence). "
                "Score is posterior_odds from the Bayesian integration "
                "(Script 46b) when available; otherwise a fallback "
                "formula based on DEG / COLOC / mouse conservation / spatial."
            ),
        },
    }


# ---------------------------------------------------------------------------
# gene_index.json rename + s5_epigenomic inject
# ---------------------------------------------------------------------------

KEY_RENAME = {
    "s2_genetic": "s3_genetic",
    "s3_essential": "s4_essential",
    "s4_epigenomic": "s5_epigenomic",
    "s5_spatial": "s6_spatial",
    "s6_singlecell": "s7_singlecell",
    "s7_mouse": "s2_mouse",
}


def update_gene_index(modalities: pd.DataFrame, atlas: pd.DataFrame):
    if not GENE_INDEX.exists():
        print(f"  [warn] gene_index.json not found at {GENE_INDEX}; skipping")
        return
    print(f"Updating {GENE_INDEX} ...")
    with open(GENE_INDEX) as fh:
        index = json.load(fh)

    # Build symbol -> s5_epigenomic lookup (float 0.0 or 1.0)
    s5_map = dict(zip(atlas["human_symbol"].values,
                      modalities["s5_epigenomic"].astype(float).values))
    s8_map = dict(zip(atlas["human_symbol"].values,
                      modalities["s8_proteomics"].astype(float).values))

    n_renamed = 0
    n_s5_added = 0
    n_s8_added = 0
    for entry in index:
        if not isinstance(entry, dict):
            continue
        ev = entry.get("evidence")
        if not isinstance(ev, dict):
            continue
        # Apply renames (preserve values, drop old key)
        new_ev = {}
        for k, v in ev.items():
            nk = KEY_RENAME.get(k, k)
            if nk != k:
                n_renamed += 1
            new_ev[nk] = v
        # Inject s5_epigenomic if atlas flag is set and field missing
        sym = entry.get("symbol")
        s5_val = s5_map.get(sym)
        if s5_val is not None and s5_val > 0 and "s5_epigenomic" not in new_ev:
            new_ev["s5_epigenomic"] = round(float(s5_val), 4)
            n_s5_added += 1
        # Inject s8_proteomics for every gene (always write; 0.0 if inactive)
        s8_val = s8_map.get(sym)
        if s8_val is not None:
            new_ev["s8_proteomics"] = round(float(s8_val), 4)
            n_s8_added += 1
        entry["evidence"] = new_ev

    with open(GENE_INDEX, "w") as fh:
        json.dump(index, fh, indent=2)
    print(f"  renamed {n_renamed} evidence keys across {len(index)} entries")
    print(f"  added s5_epigenomic to {n_s5_added} entries")
    print(f"  added s8_proteomics to {n_s8_added} entries")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if not ATLAS_CSV.exists():
        print(f"ERROR: atlas not found at {ATLAS_CSV}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading atlas from {ATLAS_CSV} ...")
    atlas = pd.read_csv(ATLAS_CSV, low_memory=False)
    print(f"  shape: {atlas.shape}")

    print("Computing 8-modality matrix ...")
    modalities = compute_modalities(atlas)
    n_active = (modalities > 0.1).sum(axis=1)
    print(f"  count distribution: {dict(n_active.value_counts().sort_index())}")
    print(f"  max count: {n_active.max()}")

    # --- Convergence matrix ---
    print("Building convergence_matrix.json ...")
    conv = build_convergence_matrix(atlas, modalities)
    conv = clean_for_json(conv)
    with open(CONVERGENCE_OUT, "w") as fh:
        json.dump(conv, fh, separators=(",", ":"))
    size_kb = CONVERGENCE_OUT.stat().st_size / 1024
    print(f"  -> {CONVERGENCE_OUT.name}: {size_kb:.1f} KB, {len(conv)} genes")

    # --- Legacy 46b archetype ranking (bayesian_ranking.json) ---
    # NOTE: canonical 46d convergence ranking is convergence_ranking.json
    print("Building bayesian_ranking.json (legacy 46b archetype similarity) ...")
    bayes = build_legacy_archetype_ranking(atlas, modalities)
    bayes = clean_for_json(bayes)
    with open(BAYES_OUT, "w") as fh:
        json.dump(bayes, fh, separators=(",", ":"))
    size_kb = BAYES_OUT.stat().st_size / 1024
    print(f"  -> {BAYES_OUT.name}: {size_kb:.1f} KB, "
          f"{len(bayes['genes'])} genes, source: {bayes['metadata']['score_source']}")

    # Verification: THRB / NR1H4 / PPARA ranks + counts
    sym_to_row = {r["gene"]: r for r in conv}
    sym_to_rank = {r["symbol"]: r for r in bayes["genes"]}
    for g in ("THRB", "NR1H4", "PPARA"):
        if g in sym_to_rank:
            r = sym_to_rank[g]
            print(f"  {g}: rank={r['rank']} pct={r['percentile']} "
                  f"score={r['score']}  count={sym_to_row.get(g,{}).get('count','?')}/8")

    # --- Update gene_index.json ---
    update_gene_index(modalities, atlas)
    size_kb = GENE_INDEX.stat().st_size / 1024
    print(f"  -> {GENE_INDEX.name}: {size_kb:.1f} KB")


if __name__ == "__main__":
    main()
