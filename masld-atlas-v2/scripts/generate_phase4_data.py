#!/usr/bin/env python3
"""Generate Phase 4 JSON data files for the MASLD Atlas v2 web app.

Produces:
  1. gwas_atac_browser.json  -- Regulatory variant data for GWAS-ATAC browser
  2. knowledge_graph.json    -- Node-edge graph for interactive network
  3. pathway_genesets.json   -- MSigDB Hallmark gene sets for client-side enrichment

Usage:
  micromamba run -n spatial python scripts/generate_phase4_data.py --output-dir public/data
"""

import argparse
import json
import math
import os
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")


def clean_nan(obj):
    """Recursively replace NaN/inf with None for JSON serialization."""
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    if isinstance(obj, dict):
        return {k: clean_nan(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [clean_nan(v) for v in obj]
    return obj


def safe_float(val, default=None):
    """Convert a value to float, returning default on failure."""
    if val is None or val == "" or (isinstance(val, float) and math.isnan(val)):
        return default
    try:
        f = float(val)
        return default if (math.isnan(f) or math.isinf(f)) else f
    except (ValueError, TypeError):
        return default


# ---------------------------------------------------------------------------
# File 1: GWAS-ATAC Browser
# ---------------------------------------------------------------------------

def generate_gwas_atac_browser(output_dir: Path):
    print("\n=== Generating gwas_atac_browser.json ===")

    # Load variant annotations (variants in scATAC peaks)
    variant_ann = pd.read_csv(
        PROJECT_ROOT / "GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv"
    )
    print(f"  Variant annotations: {len(variant_ann)} rows")

    # Load motif disruption scores
    motif_df = pd.read_csv(
        PROJECT_ROOT / "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"
    )
    print(f"  Motif disruptions: {len(motif_df)} rows")

    # Load variant overlap summary (all credible-set variants)
    overlap_summary = pd.read_csv(
        PROJECT_ROOT / "GWAS/finemapping/results/gwas_atac/variant_overlap_summary.csv"
    )
    print(f"  Overlap summary: {len(overlap_summary)} rows")

    # Load enrichment statistics
    enrichment_df = pd.read_csv(
        PROJECT_ROOT / "GWAS/finemapping/results/gwas_atac/enrichment_statistics.csv"
    )

    # --- Build motif lookup: variant_id -> list of motif disruption records ---
    motif_lookup = {}
    for _, row in motif_df.iterrows():
        vid = row["SNP_id"]
        entry = {
            "tf": row["tf_name"],
            "effect": row["effect"],
            "in_regulon": bool(row.get("motif_in_disease_regulon", False)),
        }
        ad = safe_float(row.get("alleleDiff"))
        if ad is not None:
            entry["alleleDiff"] = round(ad, 4)
        motif_lookup.setdefault(vid, []).append(entry)

    # --- Build per-variant records, keeping the best (highest PIP) entry per variant ---
    # variant_ann may have duplicates (one per cell_type overlap).
    # We want one record per variant with all overlapping cell types.
    variant_groups = variant_ann.groupby("variant_id")

    variant_records = []
    for vid, grp in variant_groups:
        best_row = grp.loc[grp["max_pip"].idxmax()]
        pip = safe_float(best_row["max_pip"], 0)
        rec_pip = safe_float(best_row["max_rec_pip"], 0)

        cell_types = sorted(grp["cell_type"].dropna().unique().tolist())

        rec = {
            "id": vid,
            "chr": best_row.get("chr_hg38", ""),
            "pos": int(best_row["pos_hg38"]) if pd.notna(best_row.get("pos_hg38")) else int(best_row["position"]),
            "pip": round(pip, 6),
            "rec_pip": round(rec_pip, 6),
            "cell_types": cell_types,
            "nearest_gene": best_row.get("nearest_gene", ""),
            "distance": int(best_row["distance_to_tss"]) if pd.notna(best_row.get("distance_to_tss")) else None,
            "motifs_disrupted": motif_lookup.get(vid, []),
        }

        # Add COLOC info if present
        coloc_pp4 = safe_float(best_row.get("coloc_best_pp4"))
        if coloc_pp4 is not None:
            rec["coloc_pp4"] = round(coloc_pp4, 4)
            rec["coloc_gwas"] = best_row.get("coloc_best_gwas", "")

        variant_records.append(rec)

    # Sort by PIP descending, take top 200
    variant_records.sort(key=lambda x: -x["pip"])
    top_variants = variant_records[:200]

    # --- Summary stats ---
    total_variants = len(overlap_summary)
    in_peaks = int(overlap_summary["overlaps_any_peak"].sum()) if "overlaps_any_peak" in overlap_summary.columns else len(variant_ann["variant_id"].unique())
    with_motif = len(set(motif_df["SNP_id"].unique()))
    # Disease regulon TFs
    regulon_tfs = motif_df[motif_df["motif_in_disease_regulon"] == True]["tf_name"].nunique() if "motif_in_disease_regulon" in motif_df.columns else 0

    # --- Cell type enrichment from precomputed enrichment_statistics.csv ---
    cell_type_enrichment = []
    for _, row in enrichment_df.iterrows():
        ct = row["cell_type"]
        if ct == "Circulating_NK_NKT":
            continue  # Skip cell type with 0 overlaps
        cell_type_enrichment.append({
            "cell_type": ct,
            "count": int(row["n_overlapping"]),
            "enrichment": round(safe_float(row["fold_enrichment"], 0), 4),
            "padj": safe_float(row.get("fisher_padj")),
        })
    cell_type_enrichment.sort(key=lambda x: -x["enrichment"])

    # --- Assemble output ---
    output = {
        "variants": top_variants,
        "summary": {
            "total_variants": int(total_variants),
            "in_peaks": int(in_peaks),
            "with_motif_disruption": int(with_motif),
            "disease_regulon_tfs": int(regulon_tfs),
        },
        "cell_type_enrichment": cell_type_enrichment,
    }

    output = clean_nan(output)
    out_path = output_dir / "gwas_atac_browser.json"
    with open(out_path, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    size_kb = out_path.stat().st_size / 1024
    print(f"  Wrote {out_path} ({size_kb:.1f} KB)")
    print(f"  Top 200 variants (by PIP), {len(cell_type_enrichment)} cell types")
    print(f"  Summary: {output['summary']}")
    return output


# ---------------------------------------------------------------------------
# File 2: Knowledge Graph
# ---------------------------------------------------------------------------

def generate_knowledge_graph(output_dir: Path):
    print("\n=== Generating knowledge_graph.json ===")

    # --- Load atlas for gene nodes ---
    atlas = pd.read_csv(
        PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv",
        low_memory=False,
    )
    print(f"  Atlas: {len(atlas)} genes x {len(atlas.columns)} columns")
    assert {"bulk_padj", "bulk_logFC"} <= set(atlas.columns), (
        "C2: atlas missing bulk_* — rebuild 27a"
    )

    # --- Load drug validation ---
    drug_val = pd.read_csv(
        PROJECT_ROOT / "RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv"
    )
    print(f"  Drug validation: {len(drug_val)} rows")

    # --- Load LINCS ranked ---
    lincs = pd.read_csv(
        PROJECT_ROOT / "RNA-seq/results/drug_repurposing/lincs_final_ranked.csv"
    )
    print(f"  LINCS ranked: {len(lincs)} rows")

    # --- Load SCENIC regulons ---
    regulons = pd.read_csv(
        PROJECT_ROOT / "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv"
    )
    print(f"  Disease regulons: {len(regulons)} rows")

    nodes = []
    edges = []
    node_ids = set()

    # -----------------------------------------------------------------------
    # 1. Gene nodes: top 200 by layers_active
    # -----------------------------------------------------------------------
    atlas["layers_active"] = pd.to_numeric(atlas["layers_active"], errors="coerce").fillna(0).astype(int)
    atlas_sorted = atlas.sort_values("layers_active", ascending=False).head(200)

    # Determine is_coloc: any coloc PP4 > 0.5 across available columns
    coloc_cols = [
        "coloc_pp4", "broadaway_coloc_pp4", "ast_coloc_pp4", "ggt_coloc_pp4",
        "pdff_coloc_pp4", "ukbb_alt_coloc_pp4",
        "bbj_alt_coloc_pp4", "bbj_ast_coloc_pp4", "bbj_ggt_coloc_pp4",
    ]
    existing_coloc_cols = [c for c in coloc_cols if c in atlas.columns]

    # Gene set for edge building
    top_gene_symbols = set()

    # Helper to build a gene node from an atlas row
    def make_gene_node(row, symbol_override=None):
        symbol = symbol_override if symbol_override else row["human_symbol"]
        padj = safe_float(row.get("bulk_padj"))
        is_deg = padj is not None and padj < 0.05
        is_coloc = False
        for cc in existing_coloc_cols:
            v = safe_float(row.get(cc))
            if v is not None and v > 0.5:
                is_coloc = True
                break
        node = {
            "id": f"gene_{symbol}",
            "type": "gene",
            "label": symbol,
            "evidence_score": int(row["layers_active"]),
            "is_deg": is_deg,
            "is_coloc": is_coloc,
        }
        lfc = safe_float(row.get("bulk_logFC"))
        if lfc is not None:
            node["logFC"] = round(lfc, 3)
        return node

    for _, row in atlas_sorted.iterrows():
        symbol = row["human_symbol"]
        if not symbol or pd.isna(symbol):
            continue
        nid = f"gene_{symbol}"
        if nid in node_ids:
            continue
        nodes.append(make_gene_node(row))
        node_ids.add(nid)
        top_gene_symbols.add(symbol)

    # Build atlas lookup for adding drug-target / regulon-target genes later
    atlas_by_symbol = atlas.set_index("human_symbol")

    # Collect gene symbols needed for edges that may not be in top 200
    extra_gene_symbols = set()

    # Drug targets from clinical validation
    for _, row in drug_val.iterrows():
        t = row.get("target_gene", "")
        if t and pd.notna(t) and t not in top_gene_symbols:
            extra_gene_symbols.add(t)

    # LINCS targets (top 20)
    lincs_sorted_pre = lincs.sort_values("composite_score", ascending=False).head(20)
    for _, row in lincs_sorted_pre.iterrows():
        target_raw = row.get("target.x", "")
        if target_raw and pd.notna(target_raw) and target_raw not in ('', '""""', '""""""""""'):
            target_clean = target_raw.strip('"')
            if target_clean and target_clean not in top_gene_symbols:
                extra_gene_symbols.add(target_clean)

    # Regulon target genes
    for _, row in regulons.iterrows():
        targets_raw = row.get("target_genes", "")
        for t in str(targets_raw).split(";"):
            t = t.strip()
            if t and not t.startswith("ENSG") and t not in top_gene_symbols:
                extra_gene_symbols.add(t)

    # Add extra gene nodes from atlas (if they exist there)
    n_extra = 0
    for sym in extra_gene_symbols:
        nid = f"gene_{sym}"
        if nid in node_ids:
            continue
        if sym in atlas_by_symbol.index:
            arow = atlas_by_symbol.loc[sym]
            if isinstance(arow, pd.DataFrame):
                arow = arow.iloc[0]
            nodes.append(make_gene_node(arow, symbol_override=sym))
            node_ids.add(nid)
            top_gene_symbols.add(sym)
            n_extra += 1

    print(f"  Gene nodes: {len(top_gene_symbols)} (200 top + {n_extra} edge-connected)")

    # -----------------------------------------------------------------------
    # 2. Drug nodes: clinical validation + top 20 LINCS
    # -----------------------------------------------------------------------
    drug_target_map = {}  # drug_id -> set of target gene symbols

    # Clinical drug validation drugs (deduplicate by drug name)
    seen_drugs = set()
    for _, row in drug_val.iterrows():
        drug_name = row["drug"]
        if drug_name in seen_drugs:
            # Same drug may appear multiple times (different targets)
            drug_id = f"drug_{drug_name.replace(' ', '_').replace('(', '').replace(')', '')}"
            target = row.get("target_gene", "")
            if target and pd.notna(target):
                drug_target_map.setdefault(drug_id, set()).add(target)
            continue

        seen_drugs.add(drug_name)
        drug_id = f"drug_{drug_name.replace(' ', '_').replace('(', '').replace(')', '')}"

        node = {
            "id": drug_id,
            "type": "drug",
            "label": drug_name,
            "stage": row.get("stage", ""),
            "moa": row.get("moa", ""),
            "atlas_support": row.get("atlas_support", ""),
        }
        nodes.append(node)
        node_ids.add(drug_id)

        target = row.get("target_gene", "")
        if target and pd.notna(target):
            drug_target_map.setdefault(drug_id, set()).add(target)

    # Top 20 LINCS drugs (by composite_score, not already in clinical)
    lincs_sorted = lincs.sort_values("composite_score", ascending=False)
    lincs_added = 0
    for _, row in lincs_sorted.iterrows():
        if lincs_added >= 20:
            break
        drug_name = row.get("display_name", row.get("pert_iname", ""))
        if not drug_name or pd.isna(drug_name) or drug_name in seen_drugs:
            continue

        seen_drugs.add(drug_name)
        drug_id = f"drug_{drug_name.replace(' ', '_').replace('(', '').replace(')', '')}"

        node = {
            "id": drug_id,
            "type": "drug",
            "label": drug_name,
            "stage": "LINCS",
            "composite_score": round(safe_float(row.get("composite_score"), 0), 3),
        }
        nodes.append(node)
        node_ids.add(drug_id)
        lincs_added += 1

        # Parse target
        target_raw = row.get("target.x", "")
        if target_raw and pd.notna(target_raw) and target_raw not in ('', '""""', '""""""""""'):
            # Clean R-escaped quotes
            target_clean = target_raw.strip('"')
            if target_clean:
                drug_target_map.setdefault(drug_id, set()).add(target_clean)

    print(f"  Drug nodes: {len(seen_drugs)} ({len(drug_val['drug'].unique())} clinical + {lincs_added} LINCS)")

    # -----------------------------------------------------------------------
    # 3. TF nodes: all disease regulons
    # -----------------------------------------------------------------------
    tf_target_map = {}  # tf_name -> list of target gene symbols

    for _, row in regulons.iterrows():
        tf = row["tf_name"]
        tf_id = f"tf_{tf}"

        targets_raw = row.get("target_genes", "")
        targets = [t.strip() for t in str(targets_raw).split(";") if t.strip() and not t.strip().startswith("ENSG")]

        node = {
            "id": tf_id,
            "type": "tf",
            "label": tf,
            "n_targets": int(row.get("n_target_genes", len(targets))),
            "activity_diff": round(safe_float(row.get("regulon_activity_diff"), 0), 4),
        }
        nodes.append(node)
        node_ids.add(tf_id)
        tf_target_map[tf] = targets

    print(f"  TF nodes: {len(regulons)}")

    # -----------------------------------------------------------------------
    # 4. Pathway nodes: top 20 most frequent pathways among all gene nodes
    # -----------------------------------------------------------------------
    pathway_counts = {}
    gene_pathway_map = {}  # symbol -> list of pathway names

    # Iterate over all gene nodes (top 200 + extra edge-connected)
    for sym in top_gene_symbols:
        if sym not in atlas_by_symbol.index:
            continue
        arow = atlas_by_symbol.loc[sym]
        if isinstance(arow, pd.DataFrame):
            arow = arow.iloc[0]
        pw_raw = arow.get("top_pathways", "")
        if not pw_raw or (isinstance(pw_raw, float) and pd.isna(pw_raw)) or pw_raw == "":
            continue
        pathways = [p.strip() for p in str(pw_raw).split(";") if p.strip()]
        gene_pathway_map[sym] = pathways
        for pw in pathways:
            pathway_counts[pw] = pathway_counts.get(pw, 0) + 1

    # Top 20 pathways
    top_pathways = sorted(pathway_counts.items(), key=lambda x: -x[1])[:20]
    pathway_set = set()

    for pw_name, count in top_pathways:
        pw_id = f"pathway_{pw_name}"
        # Clean display name
        display = pw_name
        for prefix in ("HALLMARK_", "REACTOME_", "KEGG_", "GO_"):
            if display.startswith(prefix):
                display = display[len(prefix):]
                break
        display = display.replace("_", " ").title()

        node = {
            "id": pw_id,
            "type": "pathway",
            "label": display,
            "count": count,
        }
        nodes.append(node)
        node_ids.add(pw_id)
        pathway_set.add(pw_name)

    print(f"  Pathway nodes: {len(top_pathways)}")

    # -----------------------------------------------------------------------
    # 5. Build edges
    # -----------------------------------------------------------------------

    # Drug -> gene edges
    for drug_id, targets in drug_target_map.items():
        for target in targets:
            gene_id = f"gene_{target}"
            if gene_id in node_ids:
                edges.append({
                    "source": drug_id,
                    "target": gene_id,
                    "type": "drug_target",
                })

    # TF -> gene edges (regulon membership)
    for tf, targets in tf_target_map.items():
        tf_id = f"tf_{tf}"
        for target in targets:
            gene_id = f"gene_{target}"
            if gene_id in node_ids:
                edges.append({
                    "source": tf_id,
                    "target": gene_id,
                    "type": "regulon",
                })

    # Gene -> pathway edges
    for symbol, pathways in gene_pathway_map.items():
        gene_id = f"gene_{symbol}"
        if gene_id not in node_ids:
            continue
        for pw in pathways:
            if pw in pathway_set:
                edges.append({
                    "source": gene_id,
                    "target": f"pathway_{pw}",
                    "type": "pathway_member",
                })

    print(f"  Edges: {len(edges)} (drug_target + regulon + pathway_member)")

    # -----------------------------------------------------------------------
    # Assemble and write
    # -----------------------------------------------------------------------
    output = {
        "nodes": nodes,
        "edges": edges,
        "stats": {
            "n_nodes": len(nodes),
            "n_edges": len(edges),
            "n_genes": sum(1 for n in nodes if n["type"] == "gene"),
            "n_drugs": sum(1 for n in nodes if n["type"] == "drug"),
            "n_tfs": sum(1 for n in nodes if n["type"] == "tf"),
            "n_pathways": sum(1 for n in nodes if n["type"] == "pathway"),
        },
    }

    output = clean_nan(output)
    out_path = output_dir / "knowledge_graph.json"
    with open(out_path, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    size_kb = out_path.stat().st_size / 1024
    print(f"  Wrote {out_path} ({size_kb:.1f} KB)")
    print(f"  Stats: {output['stats']}")
    return output


# ---------------------------------------------------------------------------
# File 3: Pathway Gene Sets
# ---------------------------------------------------------------------------

def generate_pathway_genesets(output_dir: Path):
    print("\n=== Generating pathway_genesets.json ===")

    gmt_path = PROJECT_ROOT / "Analysis/downstream_analysis/pathway_analysis/data/genesets/hallmark.gmt"

    gene_sets = []
    all_genes = set()

    with open(gmt_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            set_id = parts[0]
            # parts[1] is description (often "NA")
            genes = [g for g in parts[2:] if g]

            # Clean display name
            display = set_id
            if display.startswith("HALLMARK_"):
                display = display[len("HALLMARK_"):]
            display = display.replace("_", " ").title()

            gene_sets.append({
                "id": set_id,
                "name": display,
                "genes": genes,
                "size": len(genes),
            })
            all_genes.update(genes)

    # Universe size from atlas
    atlas_genes = 33943  # from atlas dimensions documented in CLAUDE.md

    output = {
        "collections": [
            {
                "name": "Hallmark",
                "gene_sets": gene_sets,
            }
        ],
        "total_sets": len(gene_sets),
        "total_genes_in_sets": len(all_genes),
        "universe_size": atlas_genes,
    }

    output = clean_nan(output)
    out_path = output_dir / "pathway_genesets.json"
    with open(out_path, "w") as f:
        json.dump(output, f, separators=(",", ":"))

    size_kb = out_path.stat().st_size / 1024
    print(f"  Wrote {out_path} ({size_kb:.1f} KB)")
    print(f"  {len(gene_sets)} Hallmark gene sets, {len(all_genes)} unique genes")
    return output


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate Phase 4 data files for MASLD Atlas v2")
    parser.add_argument("--output-dir", required=True, help="Output directory for JSON files")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Output directory: {output_dir}")
    print(f"Project root: {PROJECT_ROOT}")

    generate_gwas_atac_browser(output_dir)
    generate_knowledge_graph(output_dir)
    generate_pathway_genesets(output_dir)

    print("\n=== Phase 4 data generation complete ===")


if __name__ == "__main__":
    main()
