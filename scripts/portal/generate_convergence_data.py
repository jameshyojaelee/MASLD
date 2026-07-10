#!/usr/bin/env python3
"""Generate portal data files for the Convergence tab (Fig 5 six-channel matrix).

Produces:
  1. masld-atlas-v2/public/data/convergence_matrix.json
     - per-gene SIX-channel strength vector (0-1 per channel) + canonical tier/score.
       The six channels are the MAIN atlas modalities (Fig 5A legend), HUMAN-ONLY:
         bulk RNA-seq, GWAS-eQTL COLOC, regulatory/ATAC, spatial, single-cell, proteomics.
       Essentiality (DepMap) and mouse/cross-species are SUPPLEMENTARY, NOT convergence
       channels (2026-07-09 IA reconciliation) and are intentionally excluded here.
  2. masld-atlas-v2/public/data/bayesian_ranking.json
     - canonical 46d convergence ranking (convergence_score / convergence_rank /
       log_convergence_odds from convergence_evidence.csv). The stale March
       bayesian_posterior.csv dependency has been REMOVED.
  3. Updates masld-atlas-v2/public/data/gene_index.json IN PLACE
     - evidence-key normalisation for the gene-card fingerprint (unchanged behaviour).

Channel provenance (46d_convergence_evidence.R):
  S1  = human bulk DE            -> "bulk"       (|bulk_logFC|/max, gated by is_deg)
  S2  = genetic-causal (COLOC)   -> "coloc"      (coloc_susie_best_pp4)
  S4  = epigenomic / ATAC        -> "atac"       (GWAS-ATAC regulatory union flag)
  S5  = spatial                  -> "spatial"    (spatial_morans_i / SVG)
  S6  = single-cell              -> "singlecell" (sc_n_celltypes_sig/5 or LIANA)
  S7  = proteomics               -> "proteomics" (best_protein_logFC)
  (S3 essentiality + S8 mouse are EXCLUDED from convergence — supplementary only.)

NOTE on counting: 46d's canonical `n_modalities_active` (and the paper's Fig 5B
"387 at >=3") counts over {S1,S2,S3,S4,S5,S7} — it INCLUDES essentiality (S3) and
drops single-cell (S6, whose Wakefield-ABF BF is mathematically unreachable). This
portal deliberately counts over the SIX MAIN modalities actually displayed (SET B:
S1,S2,S4,S5,S6,S7), so the per-gene "count" is internally consistent with the dots
shown. The authoritative Tier-1 classification (677) is read verbatim from the
`tier` column and is independent of the channel set.

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
# Canonical 46d convergence ranking (score/rank/tier). Replaces the March
# bayesian_posterior.csv, which used the WRONG 6-channel set (had essentiality,
# no proteomics) and was stale.
CONVERGENCE_CSV = PROJECT_ROOT / "RNA-seq/results/multi_evidence/convergence_evidence.csv"
PORTAL_DATA = PROJECT_ROOT / "masld-atlas-v2/public/data"
CONVERGENCE_OUT = PORTAL_DATA / "convergence_matrix.json"
BAYES_OUT = PORTAL_DATA / "bayesian_ranking.json"
GENE_INDEX = PORTAL_DATA / "gene_index.json"

# Display-strength threshold for the per-gene "N channels active" count. This is a
# display heuristic over the six shown strengths, NOT the canonical BF>log(3) gate.
ACTIVE_STRENGTH = 0.1

# The six MAIN convergence channels, in Fig 5A legend order. (key, human title)
CHANNELS = [
    ("bulk", "Bulk RNA-seq (disease logFC)"),
    ("coloc", "GWAS-eQTL colocalization (SuSiE PP.H4)"),
    ("atac", "Regulatory / ATAC"),
    ("spatial", "Spatial (Moran's I / SVG)"),
    ("singlecell", "Single-cell (cross-cell-type DE / LIANA)"),
    ("proteomics", "Proteomics (best protein logFC)"),
]
CHANNEL_KEYS = [k for k, _ in CHANNELS]

# Canonical tier -> (tier_number, clean label)
TIER_MAP = {
    "1_Genetic_validated": (1, "Tier 1 · genetically validated"),
    "2_Multimodal": (2, "Tier 2 · multimodal"),
    "3_Suggestive": (3, "Tier 3 · suggestive"),
    "4_Weak": (4, "Tier 4 · weak"),
    "Excluded": (0, "Excluded"),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def clean_for_json(obj):
    if isinstance(obj, dict):
        return {k: clean_for_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [clean_for_json(v) for v in obj]
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return None if (math.isnan(v) or math.isinf(v)) else v
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
        lambda x: str(x).strip().upper() in ("TRUE", "1", "1.0")
    )


def col_or(df: pd.DataFrame, name: str, default=np.nan) -> pd.Series:
    if name in df.columns:
        return df[name]
    print(f"  [warn] column missing: {name} (using default)")
    return pd.Series(default, index=df.index, dtype=float)


# ---------------------------------------------------------------------------
# Six-channel strength computation (SET B — main modalities only)
# ---------------------------------------------------------------------------

def compute_channels(atlas: pd.DataFrame) -> pd.DataFrame:
    """Compute the six MAIN-modality strengths (0-1) per gene. Essentiality and
    mouse are intentionally NOT computed — they are supplementary, off-convergence."""
    out = pd.DataFrame(index=atlas.index)

    # bulk (S1): |bulk_logFC| / max_abs, gated by is_deg
    lfc = pd.to_numeric(col_or(atlas, "bulk_logFC"), errors="coerce")
    padj = pd.to_numeric(col_or(atlas, "bulk_padj"), errors="coerce")
    is_deg = (padj < 0.05) & (lfc.abs() > 0.3)
    max_abs = float(lfc.abs().max())
    if not math.isfinite(max_abs) or max_abs == 0:
        max_abs = 1.0
    bulk = (lfc.abs() / max_abs).clip(0.0, 1.0).fillna(0.0).where(is_deg, 0.0)
    out["bulk"] = bulk

    # coloc (S2): SuSiE PP.H4
    pp4 = pd.to_numeric(col_or(atlas, "coloc_susie_best_pp4"), errors="coerce")
    out["coloc"] = pp4.fillna(0.0).clip(0.0, 1.0)

    # atac (S4): union of GWAS-ATAC regulatory signals
    ep = []
    for c in ("gwas_variant_in_peak", "gwas_variant_da_peak"):
        if c in atlas.columns:
            ep.append(coerce_bool(atlas[c].astype(str)))
    if "gwas_max_pip_in_peak" in atlas.columns:
        ep.append(pd.to_numeric(atlas["gwas_max_pip_in_peak"], errors="coerce").fillna(0.0).gt(0.0))
    if "gwas_motif_disrupted" in atlas.columns:
        gmd = atlas["gwas_motif_disrupted"].astype(str).fillna("")
        ep.append(gmd.str.strip().ne("") & gmd.str.lower().ne("nan"))
    if "gwas_atac_regulatory_score" in atlas.columns:
        ep.append(pd.to_numeric(atlas["gwas_atac_regulatory_score"], errors="coerce").fillna(0.0).gt(0.0))
    if ep:
        m = ep[0].astype(bool)
        for s in ep[1:]:
            m = m | s.astype(bool)
        out["atac"] = m.astype(float)
    else:
        print("  [warn] no epigenomic columns found; atac=0")
        out["atac"] = 0.0

    # spatial (S5): clamp(spatial_morans_i,0,1) OR 1.0 if SVG
    mi = pd.to_numeric(col_or(atlas, "spatial_morans_i"), errors="coerce")
    svg_flag = pd.Series(False, index=atlas.index)
    if "spatial_is_svg" in atlas.columns:
        svg = atlas["spatial_is_svg"]
        svg_flag = svg.fillna(False).astype(bool) if svg.dtype == bool else coerce_bool(svg.astype(str))
    out["spatial"] = np.where(svg_flag, 1.0, mi.clip(0.0, 1.0).fillna(0.0))

    # singlecell (S6): sc_n_celltypes_sig/5 OR liana>0 -> 1.0
    sc = pd.to_numeric(col_or(atlas, "sc_n_celltypes_sig"), errors="coerce").fillna(0.0)
    liana = pd.to_numeric(col_or(atlas, "liana_n_diff_interactions"), errors="coerce").fillna(0.0)
    out["singlecell"] = np.maximum((sc / 5.0).clip(0.0, 1.0), (liana > 0).astype(float))

    # proteomics (S7): best_protein_logFC (+ significance)
    plfc = pd.to_numeric(col_or(atlas, "best_protein_logFC"), errors="coerce")
    ppadj = pd.to_numeric(col_or(atlas, "best_protein_padj"), errors="coerce")
    psig = (ppadj < 0.05)
    m_sig = (plfc.abs() / 2.0).clip(0.0, 1.0)
    m_nonsig = (plfc.abs() * 0.3).clip(0.0, 0.15)
    out["proteomics"] = pd.Series(
        np.where(psig.fillna(False), m_sig.fillna(0.0), m_nonsig.fillna(0.0)),
        index=atlas.index,
    ).fillna(0.0)

    return out[CHANNEL_KEYS]


# ---------------------------------------------------------------------------
# Canonical tier / score merge (46d convergence_evidence.csv)
# ---------------------------------------------------------------------------

def load_convergence_ranking() -> pd.DataFrame:
    if not CONVERGENCE_CSV.exists():
        print(f"ERROR: convergence_evidence.csv not found at {CONVERGENCE_CSV}", file=sys.stderr)
        sys.exit(1)
    cev = pd.read_csv(CONVERGENCE_CSV, low_memory=False)
    keep = ["human_symbol", "convergence_score", "convergence_rank",
            "log_convergence_odds", "tier", "n_modalities_active",
            "excluded_from_ranking"]
    keep = [c for c in keep if c in cev.columns]
    cev = cev[keep].copy()
    # Collapse to one row per symbol (best convergence_rank wins)
    if "convergence_rank" in cev.columns:
        cev = cev.sort_values("convergence_rank").drop_duplicates("human_symbol", keep="first")
    else:
        cev = cev.drop_duplicates("human_symbol", keep="first")
    return cev


def build_convergence_matrix(atlas: pd.DataFrame, channels: pd.DataFrame,
                             cev: pd.DataFrame) -> list:
    lfc = pd.to_numeric(col_or(atlas, "bulk_logFC"), errors="coerce")
    padj = pd.to_numeric(col_or(atlas, "bulk_padj"), errors="coerce")
    is_deg = ((padj < 0.05) & (lfc.abs() > 0.3)).fillna(False)
    pp4 = pd.to_numeric(col_or(atlas, "coloc_susie_best_pp4"), errors="coerce").fillna(0.0)
    is_coloc = (pp4 >= 0.5)
    is_druggable = pd.Series(False, index=atlas.index)
    if "dgidb_druggable" in atlas.columns:
        d = atlas["dgidb_druggable"]
        is_druggable = d.fillna(False) if d.dtype == bool else coerce_bool(d.astype(str))

    mod_vals = channels[CHANNEL_KEYS].fillna(0.0).clip(0.0, 1.0).values
    count = (mod_vals > ACTIVE_STRENGTH).sum(axis=1)

    # Canonical tier/score by symbol
    cev_by_sym = cev.set_index("human_symbol")
    symbols = atlas["human_symbol"].values

    rows = []
    for i, sym in enumerate(symbols):
        if not isinstance(sym, str) or sym == "" or sym.lower() == "nan":
            continue
        tier_num, tier_label, score, rank = 4, "Tier 4 · weak", None, None
        if sym in cev_by_sym.index:
            r = cev_by_sym.loc[sym]
            tstr = str(r.get("tier", "")) if pd.notnull(r.get("tier", None)) else ""
            tier_num, tier_label = TIER_MAP.get(tstr, (4, "Tier 4 · weak"))
            sc = r.get("convergence_score", None)
            rk = r.get("convergence_rank", None)
            score = float(sc) if pd.notnull(sc) else None
            rank = int(rk) if pd.notnull(rk) else None
        rows.append({
            "gene": sym,
            "modalities": [round(float(x), 4) for x in mod_vals[i]],
            "count": int(count[i]),
            "tier": tier_num,
            "tier_label": tier_label,
            "score": round(score, 5) if score is not None else None,
            "rank": rank,
            "is_deg": bool(is_deg.iloc[i]),
            "is_coloc": bool(is_coloc.iloc[i]),
            "is_druggable": bool(is_druggable.iloc[i]),
        })

    # Sort by canonical rank (None last), then by count desc
    rows.sort(key=lambda r: (r["rank"] if r["rank"] is not None else 10**9, -r["count"]))
    return rows


# ---------------------------------------------------------------------------
# Canonical convergence ranking (bayesian_ranking.json)
# ---------------------------------------------------------------------------

def build_ranking(atlas: pd.DataFrame, cev: pd.DataFrame, conv_rows: list) -> dict:
    flags = {r["gene"]: r for r in conv_rows}
    cev_valid = cev.dropna(subset=["convergence_rank"]) if "convergence_rank" in cev.columns else cev
    n = len(cev_valid)
    genes = []
    for _, r in cev_valid.iterrows():
        sym = r["human_symbol"]
        if not isinstance(sym, str) or sym == "" or sym.lower() == "nan":
            continue
        rank = int(r["convergence_rank"]) if pd.notnull(r.get("convergence_rank")) else None
        score = float(r["convergence_score"]) if pd.notnull(r.get("convergence_score")) else 0.0
        pct = round(100.0 * (1.0 - (rank - 1) / n), 2) if rank is not None and n else 0.0
        f = flags.get(sym, {})
        genes.append({
            "rank": rank,
            "symbol": sym,
            "score": round(score, 5),
            "percentile": pct,
            "tier": f.get("tier", 4),
            "is_deg": bool(f.get("is_deg", False)),
            "is_coloc": bool(f.get("is_coloc", False)),
            "is_druggable": bool(f.get("is_druggable", False)),
        })
    genes.sort(key=lambda r: (r["rank"] if r["rank"] is not None else 10**9))
    return {
        "genes": genes,
        "metadata": {
            "n_genes": len(genes),
            "score_source": (
                "canonical 46d convergence ranking "
                "(convergence_evidence.csv: convergence_score / convergence_rank). "
                "Tier-1 = genetically-validated convergence gate."
            ),
            "notes": (
                "rank is 1-based ascending (1 = strongest). percentile = "
                "100*(n-rank+1)/n. Score/tier are the canonical Script-46d "
                "multi-evidence convergence values (six MAIN human channels; "
                "essentiality and mouse are supplementary, off-convergence)."
            ),
        },
    }


# ---------------------------------------------------------------------------
# gene_index.json evidence-key normalisation (unchanged behaviour)
# ---------------------------------------------------------------------------

KEY_RENAME = {
    "s2_genetic": "s3_genetic",
    "s3_essential": "s4_essential",
    "s4_epigenomic": "s5_epigenomic",
    "s5_spatial": "s6_spatial",
    "s6_singlecell": "s7_singlecell",
    "s7_mouse": "s2_mouse",
}


def update_gene_index(channels: pd.DataFrame, atlas: pd.DataFrame):
    if not GENE_INDEX.exists():
        print(f"  [warn] gene_index.json not found at {GENE_INDEX}; skipping")
        return
    print(f"Updating {GENE_INDEX} ...")
    with open(GENE_INDEX) as fh:
        index = json.load(fh)
    s5_map = dict(zip(atlas["human_symbol"].values, channels["atac"].astype(float).values))
    s8_map = dict(zip(atlas["human_symbol"].values, channels["proteomics"].astype(float).values))
    n_renamed = n_s5 = n_s8 = 0
    for entry in index:
        if not isinstance(entry, dict):
            continue
        ev = entry.get("evidence")
        if not isinstance(ev, dict):
            continue
        new_ev = {}
        for k, v in ev.items():
            nk = KEY_RENAME.get(k, k)
            if nk != k:
                n_renamed += 1
            new_ev[nk] = v
        sym = entry.get("symbol")
        s5v = s5_map.get(sym)
        if s5v is not None and s5v > 0 and "s5_epigenomic" not in new_ev:
            new_ev["s5_epigenomic"] = round(float(s5v), 4)
            n_s5 += 1
        s8v = s8_map.get(sym)
        if s8v is not None:
            new_ev["s8_proteomics"] = round(float(s8v), 4)
            n_s8 += 1
        entry["evidence"] = new_ev
    with open(GENE_INDEX, "w") as fh:
        json.dump(index, fh, separators=(",", ":"))
    print(f"  renamed {n_renamed} keys; +s5 {n_s5}; +s8 {n_s8} across {len(index)} entries")


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

    print(f"Loading canonical convergence ranking from {CONVERGENCE_CSV.name} ...")
    cev = load_convergence_ranking()
    print(f"  {len(cev)} ranked genes; tier counts: {dict(cev['tier'].value_counts()) if 'tier' in cev else 'n/a'}")

    print("Computing SIX main-modality channels (bulk, coloc, atac, spatial, singlecell, proteomics) ...")
    channels = compute_channels(atlas)
    n_active = (channels > ACTIVE_STRENGTH).sum(axis=1)
    print(f"  display count distribution (over 6 shown channels): {dict(n_active.value_counts().sort_index())}")

    print("Building convergence_matrix.json ...")
    conv = clean_for_json(build_convergence_matrix(atlas, channels, cev))
    with open(CONVERGENCE_OUT, "w") as fh:
        json.dump(conv, fh, separators=(",", ":"))
    print(f"  -> {CONVERGENCE_OUT.name}: {CONVERGENCE_OUT.stat().st_size/1024:.1f} KB, {len(conv)} genes")

    print("Building bayesian_ranking.json (canonical 46d convergence ranking) ...")
    ranking = clean_for_json(build_ranking(atlas, cev, conv))
    with open(BAYES_OUT, "w") as fh:
        json.dump(ranking, fh, separators=(",", ":"))
    print(f"  -> {BAYES_OUT.name}: {BAYES_OUT.stat().st_size/1024:.1f} KB, {len(ranking['genes'])} genes")

    # --- CORRECTNESS GATE ---
    tier1 = sum(1 for r in conv if r["tier"] == 1)
    print(f"\n  [GATE] Tier-1 genes in matrix: {tier1}  (canonical target = 677)")
    for g in ("THRB", "RORA", "HKDC1", "NR1H4", "PPARA"):
        r = next((x for x in conv if x["gene"] == g), None)
        if r:
            print(f"    {g}: tier={r['tier']} rank={r['rank']} score={r['score']} count={r['count']}/6")

    update_gene_index(channels, atlas)
    print(f"  -> {GENE_INDEX.name}: {GENE_INDEX.stat().st_size/1024:.1f} KB")


if __name__ == "__main__":
    main()
