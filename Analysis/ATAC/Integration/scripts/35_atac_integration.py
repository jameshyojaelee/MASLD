#!/usr/bin/env python3
"""
35_atac_integration.py — L8 ATAC-seq integration into the MASLD Gene Catalog

Integrates chromatin accessibility evidence from three upstream ATAC-seq modules
into the MASLD Gene Catalog as Layer 8 (L8_chromatin):

  Module 1 (Mouse bulk ATAC): Promoter accessibility, differential accessibility
  Module 2 (Human scATAC):    Per-cell-type DA, chromVAR TF motif enrichment
  Module 3 (SCENIC+ GRN):     Regulon targets, enhancer-gene links

Adds 13 new columns (see L8_COLUMNS below) and updates layers_active count.
Handles partial data gracefully — any subset of modules can be missing.

Inputs:
  - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
  - Analysis/ATAC/Mouse_Bulk/results/promoter_accessibility.csv
  - Analysis/ATAC/Human_Multiome/results/scatac_da_results.csv
  - Analysis/ATAC/Human_Multiome/results/chromvar/chromvar_tf_activity.csv
  - Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv
  - Analysis/ATAC/Human_Multiome/scenic_plus/enhancer_gene_links.csv
  - Ortholog map (for mouse-to-human gene symbol mapping)

Outputs:
  - Analysis/ATAC/Integration/results/l8_atac_columns.csv
  - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (updated, with backup)
  - Analysis/ATAC/Integration/results/l8_summary_stats.txt
"""

import argparse
import hashlib
import logging
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BASE = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)

# Default input paths
DEFAULT_ATLAS = BASE / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
DEFAULT_MOUSE_DA = BASE / "Analysis/ATAC/Mouse_Bulk/results/promoter_accessibility.csv"
# Prefer corrected (hepatocyte-only, depth-adjusted logistic regression) over original
_SCATAC_DA_CORRECTED = BASE / "Analysis/ATAC/Human_Multiome/results/l8_annotated_corrected/scatac_da_gene_annotated.csv"
_SCATAC_DA_ORIGINAL = BASE / "Analysis/ATAC/Human_Multiome/results/l8_annotated/scatac_da_gene_annotated.csv"
DEFAULT_SCATAC_DA = _SCATAC_DA_CORRECTED if _SCATAC_DA_CORRECTED.exists() else _SCATAC_DA_ORIGINAL
_CHROMVAR_CORRECTED = BASE / "Analysis/ATAC/Human_Multiome/results/l8_annotated_corrected/chromvar_gene_annotated.csv"
_CHROMVAR_ORIGINAL = BASE / "Analysis/ATAC/Human_Multiome/results/l8_annotated/chromvar_gene_annotated.csv"
DEFAULT_CHROMVAR = _CHROMVAR_CORRECTED if _CHROMVAR_CORRECTED.exists() else _CHROMVAR_ORIGINAL
DEFAULT_REGULONS = BASE / "Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv"
DEFAULT_ENHANCER_LINKS = BASE / "Analysis/ATAC/Human_Multiome/scenic_plus/enhancer_gene_links.csv"
DEFAULT_ORTHOLOG_MAP = (
    BASE
    / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration"
    / "human_mouse_ortholog_comparison.csv"
)

# Output paths
OUTDIR = BASE / "Analysis/ATAC/Integration/results"
L8_COLUMNS_OUT = OUTDIR / "l8_atac_columns.csv"
SUMMARY_OUT = OUTDIR / "l8_summary_stats.txt"

# L8 column definitions: name -> (dtype, default)
L8_COLUMNS = {
    # Module 1: Mouse bulk ATAC (mapped via ortholog)
    "mouse_promoter_accessible": (bool, False),
    "mouse_da_logFC": (float, np.nan),
    "mouse_da_padj": (float, np.nan),
    # Module 2: Human scATAC
    "human_promoter_accessible": (bool, False),
    "hepatocyte_da_logFC": (float, np.nan),
    "hepatocyte_da_padj": (float, np.nan),
    "cell_type_specific_access": (str, ""),
    "chromvar_top_tf": (str, ""),
    # Module 3: SCENIC+ GRN
    "scenic_grn_target": (bool, False),
    "scenic_regulon_tf": (str, ""),
    "scenic_enhancer_link": (bool, False),
    "scenic_regulon_activity_diff": (float, np.nan),
    # Cross-species
    "cross_species_promoter_conserved": (bool, False),
}

ATAC_V3_RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
ATAC_V3_COLUMNS = (
    "atac_v3_release_id",
    "atac_v3_promoter_da_states",
    "atac_v3_promoter_da_source_dependent",
    "atac_v3_program_uids",
    "atac_v3_program_promoter_coverage_states",
    "atac_v3_program_score_states",
    "atac_v3_program_contrast_states",
    "atac_v3_program_cross_cohort_states",
    "atac_v3_testability_reasons",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def setup_logging():
    """Configure logging with timestamps."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


def load_if_exists(path, description):
    """Load a CSV if the file exists; return None with a warning otherwise."""
    path = Path(path)
    if path.exists():
        df = pd.read_csv(path)
        logging.info("Loaded %s: %d rows x %d cols", description, len(df), len(df.columns))
        return df
    logging.warning("%s not found at %s — skipping", description, path)
    return None


def sha256_file(path):
    """Return a streaming SHA256 for a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_atac_v3_root(path):
    """Validate the fixed candidate root, readiness seal, and adapter hashes."""
    supplied = Path(path)
    root = supplied.resolve()
    expected = (
        BASE
        / "Analysis/Multimodal_Program_Projection/candidates"
        / ATAC_V3_RELEASE_ID
    ).resolve()
    if root != expected or supplied.is_symlink():
        raise RuntimeError(f"unsafe ATAC v3 candidate root: {root}")

    gate_path = root / "NON_GENETIC_READY"
    gate = pd.read_csv(gate_path, sep="\t", dtype=str, keep_default_na=False)
    if gate.empty or set(gate["release_id"]) != {ATAC_V3_RELEASE_ID}:
        raise RuntimeError("invalid ATAC v3 NON_GENETIC_READY release")
    if set(gate["gate"]) != {"NON_GENETIC_READY"} or set(gate["status"]) != {"READY"}:
        raise RuntimeError("ATAC v3 NON_GENETIC_READY is not fully ready")
    for row in gate.to_dict("records"):
        artifact = (root / row["artifact"]).resolve()
        if root not in artifact.parents or not artifact.is_file():
            raise RuntimeError(f"unsafe or missing sealed ATAC v3 artifact: {artifact}")
        if sha256_file(artifact) != row["sha256"]:
            raise RuntimeError(f"sealed ATAC v3 artifact hash drift: {row['artifact']}")

    integration = pd.read_csv(
        root / "integration/integration_manifest.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    if integration.empty or set(integration["release_id"]) != {ATAC_V3_RELEASE_ID}:
        raise RuntimeError("invalid ATAC v3 integration manifest")
    for row in integration.to_dict("records"):
        if row["role"] == "integration_artifact":
            artifact = (root / "integration" / row["artifact"]).resolve()
        elif row["role"] == "sealed_non_genetic_input_manifest":
            artifact = (root / row["artifact"]).resolve()
        elif row["role"] == "post_gate_integration_producer":
            artifact = (BASE / row["artifact"]).resolve()
        else:
            raise RuntimeError(f"unknown ATAC v3 integration role: {row['role']}")
        if not artifact.is_file() or sha256_file(artifact) != row["sha256"]:
            raise RuntimeError(f"ATAC v3 integration artifact hash drift: {row['artifact']}")
    return root


def add_atac_v3_columns(l8, root):
    """Add candidate-only observability fields without changing evidence-layer votes."""
    integration = root / "integration"
    promoter = pd.read_csv(
        integration / "gene_catalog_promoter_da_adapter.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    programs = pd.read_csv(
        integration / "gene_catalog_program_atac_adapter.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    required_promoter = {
        "release_id", "gene_symbol", "lineage", "promoter_da_state",
        "n_overlapping_tested_peaks", "source_dependence",
    }
    required_program = {
        "release_id", "gene_symbol", "program_uid", "cohort", "lineage",
        "program_atac_promoter_coverage_state", "program_atac_score_state",
        "program_atac_contrast_state", "program_atac_cross_cohort_state",
        "testability_reason",
    }
    if not required_promoter.issubset(promoter.columns):
        raise RuntimeError("ATAC v3 promoter adapter schema mismatch")
    if not required_program.issubset(programs.columns):
        raise RuntimeError("ATAC v3 program adapter schema mismatch")
    if set(promoter["release_id"]) != {ATAC_V3_RELEASE_ID}:
        raise RuntimeError("ATAC v3 promoter adapter release mismatch")
    if set(programs["release_id"]) != {ATAC_V3_RELEASE_ID}:
        raise RuntimeError("ATAC v3 program adapter release mismatch")

    for column in ATAC_V3_COLUMNS:
        l8[column] = ""
    atlas_genes = set(l8["human_symbol"])
    promoter = promoter[promoter["gene_symbol"].isin(atlas_genes)].copy()
    programs = programs[programs["gene_symbol"].isin(atlas_genes)].copy()

    promoter["state_token"] = (
        promoter["lineage"] + ":" + promoter["promoter_da_state"]
        + ":n=" + promoter["n_overlapping_tested_peaks"]
    )
    promoter_states = promoter.groupby("gene_symbol")["state_token"].apply(
        lambda values: ";".join(sorted(set(values)))
    ).to_dict()
    source_dependent = promoter.groupby("gene_symbol")["source_dependence"].apply(
        lambda values: "TRUE" if "TRUE" in set(values) else "FALSE"
    ).to_dict()

    def program_tokens(frame, value_column):
        tokens = (
            frame["cohort"] + ":" + frame["lineage"] + ":" + frame[value_column]
        )
        return tokens.groupby(frame["gene_symbol"]).apply(
            lambda values: ";".join(sorted(set(values)))
        ).to_dict()

    program_uids = programs.groupby("gene_symbol")["program_uid"].apply(
        lambda values: ";".join(sorted(set(values)))
    ).to_dict()
    reasons = programs[programs["testability_reason"].ne("")].groupby("gene_symbol")[
        "testability_reason"
    ].apply(lambda values: ";".join(sorted(set(values)))).to_dict()

    mappings = {
        "atac_v3_promoter_da_states": promoter_states,
        "atac_v3_promoter_da_source_dependent": source_dependent,
        "atac_v3_program_uids": program_uids,
        "atac_v3_program_promoter_coverage_states": program_tokens(
            programs, "program_atac_promoter_coverage_state"
        ),
        "atac_v3_program_score_states": program_tokens(
            programs, "program_atac_score_state"
        ),
        "atac_v3_program_contrast_states": program_tokens(
            programs, "program_atac_contrast_state"
        ),
        "atac_v3_program_cross_cohort_states": program_tokens(
            programs, "program_atac_cross_cohort_state"
        ),
        "atac_v3_testability_reasons": reasons,
    }
    touched = set(promoter["gene_symbol"]) | set(programs["gene_symbol"])
    l8.loc[l8["human_symbol"].isin(touched), "atac_v3_release_id"] = ATAC_V3_RELEASE_ID
    for column, mapping in mappings.items():
        mask = l8["human_symbol"].isin(mapping)
        l8.loc[mask, column] = l8.loc[mask, "human_symbol"].map(mapping)
    logging.info(
        "ATAC v3 candidate fields: %d catalog genes; no layers_active vote added",
        len(touched),
    )
    return l8


def load_ortholog_map(path):
    """
    Build mouse_symbol -> human_symbol mapping from the ortholog comparison
    file produced by the integration pipeline.

    Returns a dict mapping lowercase mouse_symbol to human_symbol.
    """
    path = Path(path)
    if not path.exists():
        logging.warning("Ortholog map not found at %s", path)
        return {}

    ortho = pd.read_csv(path, usecols=["human_symbol", "mouse_symbol"])
    ortho = ortho.dropna(subset=["human_symbol", "mouse_symbol"])
    ortho = ortho[ortho["human_symbol"].str.strip().ne("")]
    ortho = ortho[ortho["mouse_symbol"].str.strip().ne("")]
    ortho = ortho.drop_duplicates(subset="mouse_symbol", keep="first")

    m2h = dict(
        zip(ortho["mouse_symbol"].str.strip().str.lower(), ortho["human_symbol"].str.strip())
    )
    logging.info("Loaded ortholog map: %d mouse-to-human mappings", len(m2h))
    return m2h


# ---------------------------------------------------------------------------
# Module processors
# ---------------------------------------------------------------------------

def process_module1(mouse_da_path, ortholog_map, atlas_genes):
    """
    Module 1: Mouse bulk ATAC promoter accessibility.

    Maps mouse gene symbols to human orthologs, then extracts:
      - mouse_promoter_accessible
      - mouse_da_logFC
      - mouse_da_padj

    Parameters
    ----------
    mouse_da_path : str or Path
        Path to promoter_accessibility.csv.
    ortholog_map : dict
        mouse_symbol (lowercase) -> human_symbol.
    atlas_genes : set
        Set of human gene symbols in the atlas.

    Returns
    -------
    pd.DataFrame or None
        DataFrame indexed on human_symbol with Module 1 columns.
    """
    df = load_if_exists(mouse_da_path, "Module 1 — Mouse bulk ATAC DA")
    if df is None:
        return None

    # Determine the gene symbol column
    gene_col = None
    for candidate in ("gene_symbol", "gene_name", "symbol"):
        if candidate in df.columns:
            gene_col = candidate
            break
    if gene_col is None:
        logging.error("Module 1: no gene symbol column found (tried gene_symbol, gene_name, symbol)")
        return None

    # Map mouse symbols to human
    df["human_symbol"] = df[gene_col].str.strip().str.lower().map(ortholog_map)
    df = df.dropna(subset=["human_symbol"])
    df = df[df["human_symbol"].isin(atlas_genes)]

    if df.empty:
        logging.warning("Module 1: no genes mapped to atlas after ortholog conversion")
        return None

    # Build output
    out = pd.DataFrame({"human_symbol": df["human_symbol"].values})

    if "mouse_promoter_accessible" in df.columns:
        out["mouse_promoter_accessible"] = df["mouse_promoter_accessible"].values
    else:
        # If column not present, infer: any row in this file = accessible
        out["mouse_promoter_accessible"] = True

    if "mouse_da_logFC" in df.columns:
        out["mouse_da_logFC"] = df["mouse_da_logFC"].astype(float).values
    elif "logFC" in df.columns:
        out["mouse_da_logFC"] = df["logFC"].astype(float).values

    if "mouse_da_padj" in df.columns:
        out["mouse_da_padj"] = df["mouse_da_padj"].astype(float).values
    elif "padj" in df.columns:
        out["mouse_da_padj"] = df["padj"].astype(float).values
    elif "FDR" in df.columns:
        out["mouse_da_padj"] = df["FDR"].astype(float).values

    out = out.drop_duplicates(subset="human_symbol", keep="first")
    logging.info("Module 1: %d genes with mouse ATAC data mapped to atlas", len(out))
    return out


def process_module2(scatac_path, atlas_genes):
    """
    Module 2: Human scATAC differential accessibility.

    Extracts per-cell-type DA and determines overall promoter accessibility.
      - human_promoter_accessible (any cell type)
      - hepatocyte_da_logFC, hepatocyte_da_padj
      - cell_type_specific_access (comma-separated list)

    Parameters
    ----------
    scatac_path : str or Path
        Path to scatac_da_results.csv.
    atlas_genes : set
        Set of human gene symbols in the atlas.

    Returns
    -------
    pd.DataFrame or None
        DataFrame indexed on human_symbol with Module 2 DA columns.
    """
    df = load_if_exists(scatac_path, "Module 2 — Human scATAC DA")
    if df is None:
        return None

    # Determine gene and cell_type columns
    gene_col = None
    for candidate in ("gene_symbol", "gene_name", "gene", "symbol", "human_symbol"):
        if candidate in df.columns:
            gene_col = candidate
            break
    if gene_col is None:
        logging.error("Module 2: no gene symbol column found")
        return None

    cell_type_col = None
    for candidate in ("cell_type", "celltype", "cluster", "annotation"):
        if candidate in df.columns:
            cell_type_col = candidate
            break
    if cell_type_col is None:
        logging.warning("Module 2: no cell_type column; treating all rows as bulk")
        df["cell_type"] = "bulk"
        cell_type_col = "cell_type"

    df = df[df[gene_col].str.strip().isin(atlas_genes)].copy()
    if df.empty:
        logging.warning("Module 2: no genes overlap with atlas")
        return None

    df[gene_col] = df[gene_col].str.strip()

    # Determine logFC / padj column names
    lfc_col = next((c for c in ("logFC", "log2FC", "avg_log2FC", "da_logFC") if c in df.columns), None)
    padj_col = next((c for c in ("padj", "p_val_adj", "FDR", "da_padj") if c in df.columns), None)

    # --- Hepatocyte-specific DA ---
    hep_rows = df[df[cell_type_col].str.lower().str.contains("hepatocyte", na=False)]
    hep_data = {}
    if not hep_rows.empty and lfc_col and padj_col:
        grp = hep_rows.groupby(gene_col)
        hep_agg = grp.agg({lfc_col: "mean", padj_col: "min"}).rename(
            columns={lfc_col: "hepatocyte_da_logFC", padj_col: "hepatocyte_da_padj"}
        )
        # Per-gene multiple testing (review 2026-05-30): taking min-padj over a
        # gene's peaks inflates the chance of a small value for genes near many
        # peaks. Keep the raw min (back-compat) but add the peak count and a
        # Sidak-corrected padj = 1-(1-min_padj)^n_peaks as the defensible column.
        hep_agg["hepatocyte_da_n_peaks"] = grp.size()
        hep_agg["hepatocyte_da_padj_sidak"] = (
            1.0 - (1.0 - hep_agg["hepatocyte_da_padj"]).pow(hep_agg["hepatocyte_da_n_peaks"])
        ).clip(upper=1.0)
        hep_data = hep_agg.to_dict("index")

    # --- Per-gene: accessible cell types ---
    accessible = df.groupby(gene_col)[cell_type_col].apply(
        lambda x: ",".join(sorted(x.unique()))
    )

    # --- Build output ---
    genes_in_data = sorted(set(df[gene_col]))
    records = []
    for gene in genes_in_data:
        rec = {
            "human_symbol": gene,
            "human_promoter_accessible": True,
            "cell_type_specific_access": accessible.get(gene, ""),
        }
        if gene in hep_data:
            rec["hepatocyte_da_logFC"] = hep_data[gene]["hepatocyte_da_logFC"]
            rec["hepatocyte_da_padj"] = hep_data[gene]["hepatocyte_da_padj"]
            rec["hepatocyte_da_n_peaks"] = hep_data[gene]["hepatocyte_da_n_peaks"]
            rec["hepatocyte_da_padj_sidak"] = hep_data[gene]["hepatocyte_da_padj_sidak"]
        records.append(rec)

    out = pd.DataFrame(records)
    out = out.drop_duplicates(subset="human_symbol", keep="first")
    logging.info("Module 2 DA: %d genes with human scATAC data", len(out))
    return out


def process_chromvar(chromvar_path, atlas_genes):
    """
    Module 2c: chromVAR TF motif enrichment.

    For each gene, identifies the most enriched TF motif at its promoter peak.
      - chromvar_top_tf

    Parameters
    ----------
    chromvar_path : str or Path
        Path to chromvar_tf_activity.csv.
    atlas_genes : set
        Set of human gene symbols in the atlas.

    Returns
    -------
    pd.DataFrame or None
        DataFrame indexed on human_symbol with chromvar_top_tf.
    """
    df = load_if_exists(chromvar_path, "Module 2c — chromVAR TF activity")
    if df is None:
        return None

    # Expect columns: tf_name / TF, cell_type, deviation / deviation_score, target_gene / gene
    tf_col = next((c for c in ("tf_name", "TF", "tf", "motif") if c in df.columns), None)
    gene_col = next(
        (c for c in ("target_gene", "gene", "gene_symbol", "nearest_gene") if c in df.columns),
        None,
    )
    dev_col = next(
        (c for c in ("logFC_deviation", "deviation", "deviation_score",
                     "mean_deviation", "score") if c in df.columns),
        None,
    )

    if tf_col is None or gene_col is None:
        logging.warning(
            "Module 2c: required columns not found (tf=%s, gene=%s). Columns: %s",
            tf_col, gene_col, list(df.columns),
        )
        return None

    df = df[df[gene_col].str.strip().isin(atlas_genes)].copy()
    if df.empty:
        logging.warning("Module 2c: no genes overlap with atlas")
        return None

    df[gene_col] = df[gene_col].str.strip()

    # Pick the TF with the highest absolute deviation per gene
    if dev_col:
        df["abs_dev"] = df[dev_col].abs()
        idx = df.groupby(gene_col)["abs_dev"].idxmax()
        top = df.loc[idx, [gene_col, tf_col]].rename(
            columns={gene_col: "human_symbol", tf_col: "chromvar_top_tf"}
        )
    else:
        # No deviation score — take the first TF per gene
        top = (
            df.drop_duplicates(subset=gene_col, keep="first")[[gene_col, tf_col]]
            .rename(columns={gene_col: "human_symbol", tf_col: "chromvar_top_tf"})
        )

    top = top.drop_duplicates(subset="human_symbol", keep="first")
    logging.info("Module 2c: %d genes with chromVAR TF annotation", len(top))
    return top


def process_scenic(regulon_path, enhancer_path, atlas_genes):
    """
    Module 3: SCENIC+ gene regulatory network.

    Extracts regulon target status, regulating TFs, enhancer-gene links,
    and regulon activity differences.
      - scenic_grn_target
      - scenic_regulon_tf
      - scenic_enhancer_link
      - scenic_regulon_activity_diff

    Parameters
    ----------
    regulon_path : str or Path
        Path to hepatocyte_regulons.csv.
    enhancer_path : str or Path
        Path to enhancer_gene_links.csv.
    atlas_genes : set
        Set of human gene symbols in the atlas.

    Returns
    -------
    pd.DataFrame or None
        DataFrame indexed on human_symbol with Module 3 columns.
    """
    regulons = load_if_exists(regulon_path, "Module 3 — SCENIC+ hepatocyte regulons")
    enhancers = load_if_exists(enhancer_path, "Module 3 — SCENIC+ enhancer-gene links")

    if regulons is None and enhancers is None:
        return None

    records = {}  # human_symbol -> dict

    # --- Regulons ---
    if regulons is not None:
        # Expected columns: TF, target_gene, [activity_diff / regulon_activity_diff]
        tf_col = next((c for c in ("TF", "tf", "tf_name", "regulator") if c in regulons.columns), None)
        target_col = next(
            (c for c in ("target_gene", "target", "gene", "gene_symbol") if c in regulons.columns),
            None,
        )
        act_col = next(
            (
                c
                for c in (
                    "regulon_activity_diff",
                    "activity_diff",
                    "mean_activity_diff",
                    "diff",
                )
                if c in regulons.columns
            ),
            None,
        )

        if tf_col is None or target_col is None:
            logging.warning(
                "Module 3 regulons: required columns not found (tf=%s, target=%s). Columns: %s",
                tf_col, target_col, list(regulons.columns),
            )
        else:
            reg = regulons[regulons[target_col].str.strip().isin(atlas_genes)].copy()
            reg[target_col] = reg[target_col].str.strip()

            # Group TFs per target gene
            tf_agg = reg.groupby(target_col)[tf_col].apply(
                lambda x: ",".join(sorted(x.dropna().unique()))
            )

            # Activity diff: take the mean across TFs for each target
            act_agg = {}
            if act_col:
                act_agg = reg.groupby(target_col)[act_col].mean().to_dict()

            for gene in tf_agg.index:
                records.setdefault(gene, {})
                records[gene]["scenic_grn_target"] = True
                records[gene]["scenic_regulon_tf"] = tf_agg[gene]
                if gene in act_agg:
                    records[gene]["scenic_regulon_activity_diff"] = act_agg[gene]

            logging.info("Module 3 regulons: %d target genes", len(tf_agg))

    # --- Enhancer-gene links ---
    if enhancers is not None:
        gene_col = next(
            (c for c in ("gene", "gene_symbol", "target_gene", "gene_name") if c in enhancers.columns),
            None,
        )
        if gene_col is None:
            logging.warning(
                "Module 3 enhancers: no gene column found. Columns: %s",
                list(enhancers.columns),
            )
        else:
            enh = enhancers[enhancers[gene_col].str.strip().isin(atlas_genes)].copy()
            enh_genes = set(enh[gene_col].str.strip().unique())
            for gene in enh_genes:
                records.setdefault(gene, {})
                records[gene]["scenic_enhancer_link"] = True
            logging.info("Module 3 enhancers: %d genes with enhancer links", len(enh_genes))

    if not records:
        return None

    rows = []
    for gene, vals in records.items():
        row = {"human_symbol": gene}
        row.update(vals)
        rows.append(row)

    out = pd.DataFrame(rows)
    out = out.drop_duplicates(subset="human_symbol", keep="first")
    logging.info("Module 3 total: %d genes with SCENIC+ evidence", len(out))
    return out


# ---------------------------------------------------------------------------
# Main integration
# ---------------------------------------------------------------------------

def build_l8(
    atlas_path,
    mouse_da_path,
    scatac_path,
    chromvar_path,
    regulon_path,
    enhancer_path,
    ortholog_path,
):
    """
    Build the L8 ATAC-seq columns and return a DataFrame keyed on human_symbol.

    Returns
    -------
    tuple of (pd.DataFrame, pd.DataFrame)
        (l8_columns DataFrame, atlas DataFrame)
    """
    # Load atlas
    atlas = pd.read_csv(atlas_path)
    logging.info("Loaded atlas: %d genes x %d cols", len(atlas), len(atlas.columns))
    atlas_genes = set(atlas["human_symbol"].dropna().unique())

    # Load ortholog map
    ortholog_map = load_ortholog_map(ortholog_path)

    # --- Process each module ---
    mod1 = process_module1(mouse_da_path, ortholog_map, atlas_genes)
    mod2 = process_module2(scatac_path, atlas_genes)
    chromvar = process_chromvar(chromvar_path, atlas_genes)
    scenic = process_scenic(regulon_path, enhancer_path, atlas_genes)

    # --- Initialize L8 frame with defaults ---
    l8 = pd.DataFrame({"human_symbol": sorted(atlas_genes)})
    for col, (dtype, default) in L8_COLUMNS.items():
        if dtype == bool:
            l8[col] = default
        elif dtype == float:
            l8[col] = np.nan
        else:
            l8[col] = ""

    # --- Merge Module 1 ---
    if mod1 is not None:
        for col in ["mouse_promoter_accessible", "mouse_da_logFC", "mouse_da_padj"]:
            if col in mod1.columns:
                merge_col(l8, mod1, "human_symbol", col)

    # --- Merge Module 2 ---
    if mod2 is not None:
        for col in ["human_promoter_accessible", "hepatocyte_da_logFC", "hepatocyte_da_padj",
                     "cell_type_specific_access"]:
            if col in mod2.columns:
                merge_col(l8, mod2, "human_symbol", col)

    # --- Merge chromVAR ---
    if chromvar is not None:
        if "chromvar_top_tf" in chromvar.columns:
            merge_col(l8, chromvar, "human_symbol", "chromvar_top_tf")

    # --- Merge SCENIC+ ---
    if scenic is not None:
        for col in ["scenic_grn_target", "scenic_regulon_tf", "scenic_enhancer_link",
                     "scenic_regulon_activity_diff"]:
            if col in scenic.columns:
                merge_col(l8, scenic, "human_symbol", col)

    # --- Cross-species conservation ---
    l8["cross_species_promoter_conserved"] = (
        l8["mouse_promoter_accessible"].fillna(False).astype(bool)
        & l8["human_promoter_accessible"].fillna(False).astype(bool)
    )

    # Clean up boolean columns: ensure proper types
    for col in ["mouse_promoter_accessible", "human_promoter_accessible",
                "scenic_grn_target", "scenic_enhancer_link",
                "cross_species_promoter_conserved"]:
        l8[col] = l8[col].fillna(False).astype(bool)

    # Clean up string columns: replace NaN with empty string
    for col in ["cell_type_specific_access", "chromvar_top_tf", "scenic_regulon_tf"]:
        l8[col] = l8[col].fillna("")

    return l8, atlas


def merge_col(target, source, key, col):
    """
    Merge a single column from source into target on key, updating in place.
    Overwrites target[col] where source has non-null values.
    """
    mapping = source.set_index(key)[col].to_dict()
    mask = target[key].isin(mapping)
    target.loc[mask, col] = target.loc[mask, key].map(mapping)


def count_l8_active(row):
    """
    Determine whether L8 provides evidence for a gene.
    Returns 1 if any L8 evidence is present, 0 otherwise.
    """
    has_mouse_da = pd.notna(row.get("mouse_da_padj")) and row.get("mouse_da_padj") < 0.1
    has_human_da = pd.notna(row.get("hepatocyte_da_padj")) and row.get("hepatocyte_da_padj") < 0.1
    has_scenic = row.get("scenic_grn_target", False)
    has_conserved = row.get("cross_species_promoter_conserved", False)
    return int(has_mouse_da or has_human_da or has_scenic or has_conserved)


def generate_summary(l8, atlas, outpath):
    """
    Write coverage and enrichment summary statistics.
    """
    lines = []
    lines.append("=" * 70)
    lines.append("L8 ATAC-seq Integration — Summary Statistics")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("=" * 70)

    total = len(l8)
    lines.append(f"\nTotal atlas genes: {total}")

    # --- Per-column coverage ---
    lines.append("\n--- Coverage per L8 column ---")
    lines.append(f"{'Column':<40} {'Non-default':>12} {'Pct':>8}")
    lines.append("-" * 62)

    for col, (dtype, default) in L8_COLUMNS.items():
        if dtype == bool:
            n = l8[col].sum()
        elif dtype == float:
            n = l8[col].notna().sum()
        else:
            n = (l8[col].str.len() > 0).sum()
        pct = 100.0 * n / total if total > 0 else 0.0
        lines.append(f"{col:<40} {n:>12,} {pct:>7.1f}%")

    # --- Module summary ---
    n_mouse = l8["mouse_promoter_accessible"].sum()
    n_human = l8["human_promoter_accessible"].sum()
    n_conserved = l8["cross_species_promoter_conserved"].sum()
    n_scenic = l8["scenic_grn_target"].sum()
    n_enhancer = l8["scenic_enhancer_link"].sum()

    lines.append("\n--- Module summary ---")
    lines.append(f"Module 1 (Mouse bulk ATAC):  {n_mouse:,} genes with promoter peaks")
    lines.append(f"Module 2 (Human scATAC):     {n_human:,} genes with promoter peaks")
    lines.append(f"Module 3 (SCENIC+ GRN):      {n_scenic:,} regulon targets, {n_enhancer:,} enhancer links")
    lines.append(f"Cross-species conserved:     {n_conserved:,} genes accessible in both species")

    # --- L8 active count ---
    l8_active = l8.apply(count_l8_active, axis=1)
    n_l8_active = l8_active.sum()
    lines.append(f"\nGenes with any L8 evidence (active): {n_l8_active:,} ({100.0 * n_l8_active / total:.1f}%)")

    # --- Enrichment: Conserved vs promoter accessibility ---
    if "is_conserved" in atlas.columns:
        merged = l8.merge(
            atlas[["human_symbol", "is_conserved"]],
            on="human_symbol",
            how="left",
        )
        merged["is_conserved"] = merged["is_conserved"].fillna(False).infer_objects(copy=False).astype(bool)

        for access_col, label in [
            ("human_promoter_accessible", "human promoter accessible"),
            ("mouse_promoter_accessible", "mouse promoter accessible"),
            ("cross_species_promoter_conserved", "cross-species conserved"),
            ("scenic_grn_target", "SCENIC+ regulon target"),
        ]:
            a = merged["is_conserved"].astype(int)
            b = merged[access_col].astype(int)
            table = pd.crosstab(a, b)

            if table.shape == (2, 2):
                odds, pval = stats.fisher_exact(table, alternative="greater")
                lines.append(
                    f"\nConserved enrichment for {label}:"
                    f" OR={odds:.2f}, p={pval:.2e}"
                    f" ({table.iloc[1, 1]:,}/{merged['is_conserved'].sum():,}"
                    f" Core genes are {label})"
                )
            else:
                lines.append(f"\nConserved enrichment for {label}: insufficient data")

    # --- Correlation with other layers ---
    lines.append("\n--- Spearman correlation: L8 columns vs bulk_logFC ---")
    if "bulk_logFC" in atlas.columns:
        merged_corr = l8.merge(atlas[["human_symbol", "bulk_logFC"]], on="human_symbol", how="left")
        for col in ["mouse_da_logFC", "hepatocyte_da_logFC", "scenic_regulon_activity_diff"]:
            valid = merged_corr[[col, "bulk_logFC"]].dropna()
            if len(valid) >= 10:
                rho, pval = stats.spearmanr(valid[col], valid["bulk_logFC"])
                lines.append(f"  {col} vs bulk_logFC: rho={rho:.3f}, p={pval:.2e} (n={len(valid):,})")
            else:
                lines.append(f"  {col} vs bulk_logFC: insufficient data (n={len(valid)})")

    summary_text = "\n".join(lines)
    with open(outpath, "w") as f:
        f.write(summary_text)
    logging.info("Summary written to %s", outpath)

    # Also print to stdout
    print("\n" + summary_text)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="L8 ATAC-seq integration into the MASLD Gene Catalog",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--atlas", default=str(DEFAULT_ATLAS),
                   help="Path to multi_evidence_atlas.csv")
    p.add_argument("--mouse-da", default=str(DEFAULT_MOUSE_DA),
                   help="Path to mouse bulk ATAC promoter_accessibility.csv")
    p.add_argument("--scatac-da", default=str(DEFAULT_SCATAC_DA),
                   help="Path to human scATAC DA results")
    p.add_argument("--chromvar", default=str(DEFAULT_CHROMVAR),
                   help="Path to chromVAR TF activity results")
    p.add_argument("--regulons", default=str(DEFAULT_REGULONS),
                   help="Path to SCENIC+ hepatocyte regulons")
    p.add_argument("--enhancer-links", default=str(DEFAULT_ENHANCER_LINKS),
                   help="Path to SCENIC+ enhancer-gene links")
    p.add_argument("--ortholog-map", default=str(DEFAULT_ORTHOLOG_MAP),
                   help="Path to ortholog comparison CSV")
    p.add_argument("--outdir", default=str(OUTDIR),
                   help="Output directory for L8 columns and summary")
    p.add_argument("--atac-v3-root", default=None,
                   help="Sealed candidate-only ATAC Context v3 root")
    p.add_argument("--dry-run", action="store_true",
                   help="Show what would be merged without modifying the atlas")
    p.add_argument("--force", action="store_true",
                   help="Overwrite existing L8 columns if they exist in the atlas")
    return p.parse_args()


def main():
    setup_logging()
    args = parse_args()

    logging.info("=== L8 ATAC-seq Integration (Script 35) ===")
    logging.info("Atlas:         %s", args.atlas)
    logging.info("Mouse DA:      %s", args.mouse_da)
    logging.info("scATAC DA:     %s", args.scatac_da)
    logging.info("chromVAR:      %s", args.chromvar)
    logging.info("Regulons:      %s", args.regulons)
    logging.info("Enhancer links:%s", args.enhancer_links)
    logging.info("Ortholog map:  %s", args.ortholog_map)
    logging.info("Dry run:       %s", args.dry_run)
    logging.info("Force:         %s", args.force)

    atac_v3_root = None
    if args.atac_v3_root:
        if not args.dry_run:
            logging.error("ATAC v3 candidate mode requires --dry-run; canonical writes are disabled")
            sys.exit(1)
        atac_v3_root = validate_atac_v3_root(args.atac_v3_root)
        outdir_resolved = Path(args.outdir).resolve()
        if atac_v3_root not in outdir_resolved.parents:
            logging.error("ATAC v3 candidate outputs must remain beneath %s", atac_v3_root)
            sys.exit(1)
        logging.info("ATAC v3 root:  %s", atac_v3_root)

    # Validate atlas exists
    if not Path(args.atlas).exists():
        logging.error("Atlas file not found: %s", args.atlas)
        sys.exit(1)

    # Create output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Build L8 columns
    l8, atlas = build_l8(
        atlas_path=args.atlas,
        mouse_da_path=args.mouse_da,
        scatac_path=args.scatac_da,
        chromvar_path=args.chromvar,
        regulon_path=args.regulons,
        enhancer_path=args.enhancer_links,
        ortholog_path=args.ortholog_map,
    )
    if atac_v3_root is not None:
        l8 = add_atac_v3_columns(l8, atac_v3_root)

    # Save L8 columns as standalone file
    l8_out = outdir / "l8_atac_columns.csv"
    l8.to_csv(l8_out, index=False)
    logging.info("L8 columns saved to %s (%d genes)", l8_out, len(l8))

    # Generate summary
    summary_path = outdir / "l8_summary_stats.txt"
    generate_summary(l8, atlas, summary_path)

    # --- Merge into atlas ---
    if args.dry_run:
        n_candidate_columns = len(ATAC_V3_COLUMNS) if atac_v3_root is not None else 0
        logging.info(
            "DRY RUN: would merge %d L8 columns%s into atlas (%d genes)",
            len(L8_COLUMNS),
            f" plus {n_candidate_columns} candidate-only ATAC v3 fields"
            if n_candidate_columns else "",
            len(atlas),
        )
        logging.info("DRY RUN: atlas would go from %d to %d columns",
                     len(atlas.columns), len(atlas.columns) + len(L8_COLUMNS) + n_candidate_columns)
        # Check how many genes would get L8 evidence
        l8_active = l8.apply(count_l8_active, axis=1)
        logging.info("DRY RUN: %d genes would have active L8 evidence", l8_active.sum())
        return

    # Check for existing L8 columns
    existing_l8 = [c for c in L8_COLUMNS if c in atlas.columns]
    if existing_l8 and not args.force:
        logging.error(
            "Atlas already contains L8 columns: %s. Use --force to overwrite.",
            existing_l8,
        )
        sys.exit(1)
    elif existing_l8 and args.force:
        logging.info("Removing existing L8 columns: %s", existing_l8)
        atlas = atlas.drop(columns=existing_l8)

    # Backup atlas
    atlas_path = Path(args.atlas)
    backup_path = atlas_path.with_suffix(".csv.bak")
    shutil.copy2(atlas_path, backup_path)
    logging.info("Atlas backed up to %s", backup_path)

    # Merge L8 columns into atlas
    l8_cols_only = l8.drop(columns=["human_symbol"]).columns.tolist()
    atlas = atlas.merge(
        l8[["human_symbol"] + l8_cols_only],
        on="human_symbol",
        how="left",
    )

    # Fill NaN defaults for genes not in L8
    for col, (dtype, default) in L8_COLUMNS.items():
        if col in atlas.columns:
            if dtype == bool:
                atlas[col] = atlas[col].fillna(False).astype(bool)
            elif dtype == str:
                atlas[col] = atlas[col].fillna("")

    # Update layers_active count
    if "layers_active" in atlas.columns:
        l8_active_map = l8.set_index("human_symbol").apply(count_l8_active, axis=1).to_dict()
        atlas["layers_active"] = atlas.apply(
            lambda row: (
                row["layers_active"]
                + l8_active_map.get(row["human_symbol"], 0)
            ),
            axis=1,
        )
        logging.info("Updated layers_active counts")

    # Save updated atlas
    atlas.to_csv(atlas_path, index=False)
    logging.info(
        "Updated atlas saved: %d genes x %d cols -> %s",
        len(atlas),
        len(atlas.columns),
        atlas_path,
    )
    logging.info("=== L8 integration complete ===")


if __name__ == "__main__":
    main()
