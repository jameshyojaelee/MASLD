#!/usr/bin/env python
"""
312a_hepatocyte_literature_benchmark.py
Benchmark hepatocyte subclustering results against published competitor signatures.

Scores published hepatocyte programs (Tzouanas, Li, Govaere, zonation) on our
100K hepatocyte subset, computes Jaccard overlap with our meta-subtype markers,
and validates HMGCS2 trajectory + SOX4/RELB regulon status.

Outputs:
  benchmark/literature_program_scores.csv
  benchmark/literature_overlap_matrix.csv
  benchmark/hmgcs2_trajectory.csv
  benchmark/sox4_relb_check.txt
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

warnings.filterwarnings("ignore")

# ── Paths ────────────────────────────────────────────────────────────────────
PROJECT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
HEPATOCYTE_H5AD = os.path.join(
    PROJECT,
    "Analysis/SingleCell/results_gpu_v2/pseudotime/Hepatocytes_subset.h5ad",
)
META_SUBTYPE_MAPPING = os.path.join(
    PROJECT,
    "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/meta_subtype_mapping.csv",
)
HEPATOCYTE_METADATA = os.path.join(
    PROJECT,
    "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv",
)
SUBTYPE_MARKERS = os.path.join(
    PROJECT,
    "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/subtype_markers.csv",
)
DISEASE_REGULONS = os.path.join(
    PROJECT,
    "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv",
)
OUTDIR = os.path.join(
    PROJECT,
    "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/benchmark",
)
os.makedirs(OUTDIR, exist_ok=True)

# ── Published gene sets ──────────────────────────────────────────────────────
# Tzouanas et al. (Cell 2026) — 4 hepatocyte stress programs
PUBLISHED_PROGRAMS = {
    "Tzouanas_LI": [
        "TNF", "IL6", "CXCL8", "CCL2", "COL1A1", "TGFB1", "ACTA2", "TIMP1",
        "IL1B", "CXCL1", "CCL20", "SERPINE1", "MMP9", "ICAM1", "VCAM1",
    ],
    "Tzouanas_LD": [
        "HMGCS2", "CYP2E1", "CYP3A4", "ALB", "APOB", "PCK1", "G6PC",
        "ALDOB", "HAL", "ASS1", "TAT", "HPD", "SDS", "AGXT", "OTC",
        "ARG1", "CPS1", "ASGR1", "HGD", "ABCB11",
    ],
    "Tzouanas_SU": [
        "LGALS3", "SPP1", "CD9", "TREM2", "GPNMB", "FABP5", "CTSB",
        "CTSD", "LAMP1", "LIPA", "NPC2", "GRN", "APOE", "LPL",
    ],
    "Tzouanas_SD": [
        "GLUL", "AXIN2", "WNT2", "RSPO3", "LGR5", "TBX3", "RNF43",
        "ODAM", "CYP1A2", "CYP2A6",
    ],
    # Li et al. (Nat Genet 2025)
    "Li_LAM": [
        "MITF", "TREM2", "SPP1", "CD9", "GPNMB", "LGALS3", "FABP5",
        "LIPA", "CTSD", "CTSB", "LPL", "APOE",
    ],
    "Li_profibrotic": [
        "RSPO3", "LGR6", "WNT2", "WNT9B", "COL1A1", "COL3A1", "COL1A2",
        "LUM", "DCN", "PDGFRA", "ACTA2",
    ],
    # Govaere et al. (Nat Metab 2023) — NMF subtypes
    "Govaere_metabolic": [
        "FASN", "SCD", "ACACA", "PLIN2", "DGAT2", "FABP1", "APOB",
        "MTTP", "APOA1", "APOC3", "HMGCR", "SREBF1", "PPARA", "ACOX1",
    ],
    "Govaere_inflammatory": [
        "ADAMTSL2", "AKR1B10", "CFHR4", "TREM2", "COL1A1", "LUM",
        "VCAN", "THBS2", "CCL2", "CCL20", "CXCL8", "IL1B", "SPP1",
    ],
    # Zonation markers (Halpern / Aizarani)
    "Zonation_periportal": [
        "ASS1", "HAL", "CPS1", "ALDOB", "PCK1", "SLC1A2", "HAMP",
        "SDS", "TAT", "ARG1", "OTC", "AGXT",
    ],
    "Zonation_pericentral": [
        "GLUL", "CYP2E1", "CYP1A2", "CYP3A4", "AXIN2", "ADH4",
        "CYP2A6", "CYP2C9", "CYP2D6",
    ],
}

print("=" * 70)
print("312a — Hepatocyte Literature Benchmark")
print("=" * 70)

# ── 1. Load data ─────────────────────────────────────────────────────────────
print("\n[1] Loading hepatocyte h5ad (100K cells, all genes)...")
adata = sc.read_h5ad(HEPATOCYTE_H5AD)
print(f"    Shape: {adata.shape}")

# Ensure we have gene names in var index
print(f"    Gene index examples: {list(adata.var_names[:5])}")

# ── 2. Load meta-subtype mapping and transfer labels ─────────────────────────
print("\n[2] Loading meta-subtype mapping and transferring labels...")
mapping = pd.read_csv(META_SUBTYPE_MAPPING)
subtype_to_meta = dict(zip(mapping["subtype"].astype(str), mapping["meta_subtype"]))

# Load full metadata to get subtype assignments for the 100K subset
full_meta = pd.read_csv(HEPATOCYTE_METADATA, index_col=0)
print(f"    Full metadata: {full_meta.shape[0]} cells")

# Match barcodes between pseudotime subset and full metadata
common_barcodes = adata.obs.index.intersection(full_meta.index)
print(f"    Common barcodes: {len(common_barcodes)} / {adata.n_obs}")

# Transfer hepatocyte_subtype and disease_stage_coarse
adata.obs["hepatocyte_subtype"] = np.nan
adata.obs.loc[common_barcodes, "hepatocyte_subtype"] = (
    full_meta.loc[common_barcodes, "hepatocyte_subtype"].astype(str)
)

adata.obs["disease_stage_coarse"] = np.nan
adata.obs.loc[common_barcodes, "disease_stage_coarse"] = (
    full_meta.loc[common_barcodes, "disease_stage_coarse"]
)

# Map subtype -> meta_subtype
adata.obs["meta_subtype"] = adata.obs["hepatocyte_subtype"].map(subtype_to_meta)

n_mapped = adata.obs["meta_subtype"].notna().sum()
print(f"    Cells with meta-subtype label: {n_mapped}")
print(f"    Meta-subtype distribution:")
print(adata.obs["meta_subtype"].value_counts().to_string())

# ── 3. Score published gene sets ─────────────────────────────────────────────
print("\n[3] Scoring published gene sets with sc.tl.score_genes()...")
available_genes = set(adata.var_names)
coverage_report = {}

for program_name, gene_list in PUBLISHED_PROGRAMS.items():
    present = [g for g in gene_list if g in available_genes]
    missing = [g for g in gene_list if g not in available_genes]
    coverage_report[program_name] = {
        "total": len(gene_list),
        "present": len(present),
        "missing_genes": ", ".join(missing) if missing else "none",
    }
    print(f"    {program_name}: {len(present)}/{len(gene_list)} genes present")
    if missing:
        print(f"      Missing: {', '.join(missing)}")

    if len(present) >= 2:
        sc.tl.score_genes(adata, gene_list=present, score_name=program_name)
        # Ensure float dtype for anndata compatibility
        adata.obs[program_name] = np.array(adata.obs[program_name].values, dtype=float)
    else:
        print(f"      SKIPPED (fewer than 2 genes present)")
        adata.obs[program_name] = np.nan

# ── 4. Compute per-meta-subtype mean scores ─────────────────────────────────
print("\n[4] Computing per-meta-subtype mean scores...")
score_cols = list(PUBLISHED_PROGRAMS.keys())
cells_with_meta = adata.obs.dropna(subset=["meta_subtype"])

program_scores = cells_with_meta.groupby("meta_subtype")[score_cols].mean()
program_scores = program_scores.round(4)

# Add cell count
program_scores.insert(0, "n_cells", cells_with_meta.groupby("meta_subtype").size())

# Sort meta-subtypes in biological order
meta_order = ["Healthy", "Neutral", "Disease-Neutral", "Disease-Associated", "Disease-Progressor"]
program_scores = program_scores.reindex([m for m in meta_order if m in program_scores.index])

out_scores = os.path.join(OUTDIR, "literature_program_scores.csv")
program_scores.to_csv(out_scores)
print(f"    Saved: {out_scores}")
print(program_scores.to_string())

# ── 5. Jaccard overlap with our meta-subtype markers ─────────────────────────
print("\n[5] Computing Jaccard overlap with our top-50 markers per subtype...")
markers_df = pd.read_csv(SUBTYPE_MARKERS)

# Get top 50 markers per subtype by score
our_markers = {}
for st, grp in markers_df.groupby("subtype"):
    top50 = grp.nlargest(50, "scores")["names"].tolist()
    our_markers[str(st)] = set(top50)

# Build meta-subtype marker sets (union of constituent subtypes)
meta_markers = {}
for st_str, meta in subtype_to_meta.items():
    if meta not in meta_markers:
        meta_markers[meta] = set()
    if st_str in our_markers:
        meta_markers[meta].update(our_markers[st_str])

# Compute Jaccard: published program vs meta-subtype markers
jaccard_matrix = pd.DataFrame(
    index=list(PUBLISHED_PROGRAMS.keys()),
    columns=[m for m in meta_order if m in meta_markers],
    dtype=float,
)

for prog_name, gene_list in PUBLISHED_PROGRAMS.items():
    prog_set = set(gene_list)
    for meta_name in jaccard_matrix.columns:
        marker_set = meta_markers[meta_name]
        intersection = len(prog_set & marker_set)
        union = len(prog_set | marker_set)
        jaccard_matrix.loc[prog_name, meta_name] = (
            round(intersection / union, 4) if union > 0 else 0.0
        )

# Add overlap gene columns
valid_meta_names = ["Healthy", "Neutral", "Disease-Neutral", "Disease-Associated", "Disease-Progressor"]
valid_meta_names = [m for m in valid_meta_names if m in meta_markers]
overlap_data = {}
for prog_name, gene_list in PUBLISHED_PROGRAMS.items():
    prog_set = set(gene_list)
    for meta_name in valid_meta_names:
        marker_set = meta_markers[meta_name]
        overlap_genes = sorted(prog_set & marker_set)
        col_name = f"{meta_name}_overlap_genes"
        if col_name not in overlap_data:
            overlap_data[col_name] = {}
        overlap_data[col_name][prog_name] = ";".join(overlap_genes) if overlap_genes else ""

for col_name, vals in overlap_data.items():
    jaccard_matrix[col_name] = pd.Series(vals)

out_jaccard = os.path.join(OUTDIR, "literature_overlap_matrix.csv")
jaccard_matrix.to_csv(out_jaccard)
print(f"    Saved: {out_jaccard}")
# Print just the Jaccard values
print(jaccard_matrix[[c for c in jaccard_matrix.columns if "_overlap_genes" not in c]].to_string())

# ── 6. HMGCS2 trajectory ────────────────────────────────────────────────────
print("\n[6] Extracting HMGCS2 trajectory...")
hmgcs2_results = []

if "HMGCS2" in available_genes:
    # Get HMGCS2 expression
    hmgcs2_idx = list(adata.var_names).index("HMGCS2")
    if sparse.issparse(adata.X):
        hmgcs2_expr = np.asarray(adata.X[:, hmgcs2_idx].todense()).flatten()
    else:
        hmgcs2_expr = np.asarray(adata.X[:, hmgcs2_idx]).flatten()

    adata.obs["HMGCS2_expr"] = np.array(hmgcs2_expr, dtype=float)

    # Per meta-subtype
    for meta in meta_order:
        mask = adata.obs["meta_subtype"] == meta
        if mask.sum() > 0:
            vals = adata.obs.loc[mask, "HMGCS2_expr"]
            hmgcs2_results.append({
                "group_type": "meta_subtype",
                "group": meta,
                "n_cells": int(mask.sum()),
                "mean_expr": round(float(vals.mean()), 4),
                "median_expr": round(float(vals.median()), 4),
                "pct_expressing": round(float((vals > 0).mean() * 100), 2),
            })

    # Per disease_stage_coarse
    stage_order = ["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"]
    for stage in stage_order:
        mask = adata.obs["disease_stage_coarse"] == stage
        if mask.sum() > 0:
            vals = adata.obs.loc[mask, "HMGCS2_expr"]
            hmgcs2_results.append({
                "group_type": "disease_stage",
                "group": stage,
                "n_cells": int(mask.sum()),
                "mean_expr": round(float(vals.mean()), 4),
                "median_expr": round(float(vals.median()), 4),
                "pct_expressing": round(float((vals > 0).mean() * 100), 2),
            })
else:
    print("    WARNING: HMGCS2 not found in gene set")

hmgcs2_df = pd.DataFrame(hmgcs2_results)
out_hmgcs2 = os.path.join(OUTDIR, "hmgcs2_trajectory.csv")
hmgcs2_df.to_csv(out_hmgcs2, index=False)
print(f"    Saved: {out_hmgcs2}")
print(hmgcs2_df.to_string(index=False))

# ── 7. SOX4/RELB check in SCENIC+ disease regulons ──────────────────────────
print("\n[7] Checking SOX4/RELB in SCENIC+ disease regulons...")
sox4_relb_lines = []

if os.path.exists(DISEASE_REGULONS):
    regulons = pd.read_csv(DISEASE_REGULONS)
    sox4_relb_lines.append(f"Disease regulons file: {DISEASE_REGULONS}")
    sox4_relb_lines.append(f"Total disease regulons: {len(regulons)}")
    sox4_relb_lines.append(f"TF names: {', '.join(sorted(regulons['tf_name'].tolist()))}")
    sox4_relb_lines.append("")

    for tf in ["SOX4", "RELB"]:
        match = regulons[regulons["tf_name"] == tf]
        if len(match) > 0:
            row = match.iloc[0]
            sox4_relb_lines.append(f"{tf}: FOUND in disease regulons")
            sox4_relb_lines.append(f"  Regulon ID: {row['regulon_id']}")
            sox4_relb_lines.append(f"  N target genes: {row['n_target_genes']}")
            sox4_relb_lines.append(f"  Activity diff (MASLD - Normal): {row['activity_pval']:.2e}")
            sox4_relb_lines.append(f"  Adjusted p-value: {row['activity_padj']:.2e}")
        else:
            sox4_relb_lines.append(f"{tf}: NOT FOUND in disease regulons")

        # Also check if TF appears as a target gene in any regulon
        target_hits = []
        for _, r in regulons.iterrows():
            targets = str(r.get("target_genes", "")).split(";")
            if tf in targets:
                target_hits.append(r["tf_name"])
        if target_hits:
            sox4_relb_lines.append(f"  {tf} is a target gene of: {', '.join(target_hits)}")
        sox4_relb_lines.append("")
else:
    sox4_relb_lines.append(f"Disease regulons file NOT FOUND: {DISEASE_REGULONS}")

out_sox4 = os.path.join(OUTDIR, "sox4_relb_check.txt")
with open(out_sox4, "w") as f:
    f.write("\n".join(sox4_relb_lines))
print(f"    Saved: {out_sox4}")
for line in sox4_relb_lines:
    print(f"    {line}")

# ── 8. Coverage summary ─────────────────────────────────────────────────────
print("\n[8] Gene coverage summary:")
coverage_df = pd.DataFrame(coverage_report).T
coverage_df.index.name = "program"
out_coverage = os.path.join(OUTDIR, "gene_coverage_report.csv")
coverage_df.to_csv(out_coverage)
print(f"    Saved: {out_coverage}")
print(coverage_df[["total", "present"]].to_string())

print("\n" + "=" * 70)
print("DONE. All outputs in:", OUTDIR)
print("=" * 70)
