#!/usr/bin/env python3
"""226b_tissue_plasma_bridge.py

Map top-20 plasma proteins (from 226a SHAP importance) to the tissue
transcriptomic atlas, F2 switch signatures, and drug targets.
Pure data joining — no ML.

Inputs
------
1. 226a_protein_importance.csv   – SHAP importance per protein × target
2. multi_evidence_atlas.csv      – 33,943-gene tissue atlas
3. plasma_switch_signatures.csv  – onset / switch / late gene classes
4. clinical_drug_validation_table.csv – MASLD clinical drug targets
5. dream_results_ashr.csv        – tissue DEGs (ashr shrinkage)

Outputs (all to results/multiprogram/)
------
1. 226b_tissue_bridge.csv        – per-protein × target annotation
2. 226b_enrichment.csv           – Fisher's exact DEG enrichment per target
3. 226b_drug_overlaps.csv        – proteins that are MASLD drug targets
"""

import logging
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INTEGRATION = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration",
)
RESULTS_MP = os.path.join(INTEGRATION, "results/multiprogram")
RESULTS_ME = os.path.join(BASE, "RNA-seq/results/multi_evidence")
RESULTS_DR = os.path.join(BASE, "RNA-seq/results/drug_repurposing")
RESULTS_INT = os.path.join(INTEGRATION, "results/integration")

PATH_IMPORTANCE = os.path.join(RESULTS_MP, "226a_protein_importance.csv")
PATH_ATLAS = os.path.join(RESULTS_ME, "multi_evidence_atlas.csv")
PATH_SWITCH = os.path.join(RESULTS_MP, "plasma_switch_signatures.csv")
PATH_DRUGS = os.path.join(RESULTS_DR, "clinical_drug_validation_table.csv")
PATH_DREAM = os.path.join(RESULTS_INT, "dream_results_ashr.csv")

OUT_BRIDGE = os.path.join(RESULTS_MP, "226b_tissue_bridge.csv")
OUT_ENRICHMENT = os.path.join(RESULTS_MP, "226b_enrichment.csv")
OUT_DRUG_OVERLAPS = os.path.join(RESULTS_MP, "226b_drug_overlaps.csv")


def load_importance(path: str) -> pd.DataFrame:
    """Load 226a protein importance and keep top-20 per target."""
    log.info("Loading protein importance from %s", path)
    df = pd.read_csv(path)
    log.info("  %d rows, %d targets, columns: %s",
             len(df), df["target"].nunique(), list(df.columns))
    top20 = df[df["rank"] <= 20].copy()
    log.info("  Top-20 per target: %d rows", len(top20))
    return top20


def load_atlas(path: str) -> pd.DataFrame:
    """Load multi-evidence atlas, keep relevant columns."""
    log.info("Loading atlas from %s", path)
    keep_cols = [
        "human_symbol",
        "dream_logFC", "dream_padj", "dream_tstat",
        "is_conserved", "attribution_class",
        "dgidb_druggable", "opentargets_drug",
        "layers_active", "sources_active",
        "f2_inflection_logFC", "f2_inflection_padj",
    ]
    df = pd.read_csv(path, usecols=keep_cols)
    log.info("  %d genes loaded", len(df))
    return df


def load_switch_signatures(path: str) -> pd.DataFrame:
    """Load F2 switch signatures (onset / switch / late)."""
    log.info("Loading switch signatures from %s", path)
    df = pd.read_csv(path)
    log.info("  %d entries, signatures: %s",
             len(df), df["signature"].value_counts().to_dict())
    # Keep unique gene_symbol → signature mapping (take first if duplicated)
    df = df.drop_duplicates(subset="gene_symbol", keep="first")
    return df[["gene_symbol", "signature"]].rename(
        columns={"signature": "f2_switch_class"}
    )


def load_drug_targets(path: str) -> pd.DataFrame:
    """Load clinical drug validation table and build gene → drug list."""
    log.info("Loading drug targets from %s", path)
    df = pd.read_csv(path)
    log.info("  %d drug entries", len(df))
    # Build gene → comma-separated drug list
    drug_map = (
        df.groupby("target_gene")["drug"]
        .apply(lambda x: "; ".join(sorted(x.unique())))
        .reset_index()
        .rename(columns={"drug": "drug_target_of"})
    )
    log.info("  %d unique drug-target genes", len(drug_map))
    return drug_map


def load_dream_degs(path: str) -> pd.DataFrame:
    """Load dream DEG results (ashr). Mark DEGs: padj<0.05, |shrunk_logFC|>0.2."""
    log.info("Loading dream DEGs from %s", path)
    df = pd.read_csv(path, usecols=["symbol", "padj", "shrunk_logFC", "lfsr"])
    df = df.dropna(subset=["symbol"]).drop_duplicates(subset="symbol", keep="first")
    df["is_tissue_deg"] = (df["padj"] < 0.05) & (df["shrunk_logFC"].abs() > 0.2)
    log.info("  %d genes, %d DEGs (padj<0.05, |logFC|>0.2)",
             len(df), df["is_tissue_deg"].sum())
    return df.rename(columns={
        "shrunk_logFC": "dream_shrunk_logFC_ashr",
        "padj": "dream_padj_ashr",
    })


def run_enrichment(bridge: pd.DataFrame, dream: pd.DataFrame,
                   all_proteins: pd.Series) -> pd.DataFrame:
    """Fisher's exact test: are top-20 enriched for tissue DEGs vs all Olink proteins?"""
    log.info("Running Fisher's exact enrichment tests")

    # Build DEG status for all Olink proteins
    all_prot_set = set(all_proteins.unique())
    deg_set = set(
        dream.loc[dream["is_tissue_deg"], "symbol"].values
    )

    n_all_total = len(all_prot_set)
    n_all_deg = len(all_prot_set & deg_set)

    rows = []
    for target, grp in bridge.groupby("target"):
        top20_proteins = set(grp["protein"].values)
        n_top20_total = len(top20_proteins)
        n_top20_deg = len(top20_proteins & deg_set)

        # 2×2: [[top20_deg, top20_nondeg], [other_deg, other_nondeg]]
        n_other_deg = n_all_deg - n_top20_deg
        n_other_nondeg = (n_all_total - n_top20_total) - n_other_deg

        table = [
            [n_top20_deg, n_top20_total - n_top20_deg],
            [n_other_deg, n_other_nondeg],
        ]
        odds_ratio, pvalue = fisher_exact(table, alternative="greater")

        rows.append({
            "target": target,
            "n_top20_deg": n_top20_deg,
            "n_top20_total": n_top20_total,
            "n_all_deg": n_all_deg,
            "n_all_total": n_all_total,
            "odds_ratio": round(odds_ratio, 3),
            "fisher_pvalue": pvalue,
        })
        log.info("  %s: %d/%d top-20 are DEGs (OR=%.2f, p=%.3g)",
                 target, n_top20_deg, n_top20_total, odds_ratio, pvalue)

    return pd.DataFrame(rows)


def main():
    log.info("=" * 60)
    log.info("226b — Tissue–Plasma Bridge")
    log.info("=" * 60)

    # ------------------------------------------------------------------
    # 1. Load all inputs
    # ------------------------------------------------------------------
    importance = load_importance(PATH_IMPORTANCE)
    atlas = load_atlas(PATH_ATLAS)
    switch = load_switch_signatures(PATH_SWITCH)
    drug_map = load_drug_targets(PATH_DRUGS)
    dream = load_dream_degs(PATH_DREAM)

    # All unique proteins across all targets (for enrichment denominator)
    all_proteins = pd.read_csv(PATH_IMPORTANCE, usecols=["protein"])["protein"]
    log.info("Total unique Olink proteins across all targets: %d",
             all_proteins.nunique())

    # ------------------------------------------------------------------
    # 2. Build bridge table via successive left joins
    # ------------------------------------------------------------------
    bridge = importance.copy()

    # Rename for clarity
    bridge = bridge.rename(columns={"rank": "shap_rank"})

    # Join atlas
    bridge = bridge.merge(
        atlas, left_on="protein", right_on="human_symbol", how="left"
    ).drop(columns=["human_symbol"], errors="ignore")

    # Join switch signatures
    bridge = bridge.merge(
        switch, left_on="protein", right_on="gene_symbol", how="left"
    ).drop(columns=["gene_symbol"], errors="ignore")
    bridge["f2_switch_class"] = bridge["f2_switch_class"].fillna("none")

    # Join drug targets
    bridge = bridge.merge(
        drug_map, left_on="protein", right_on="target_gene", how="left"
    ).drop(columns=["target_gene"], errors="ignore")
    bridge["drug_target_of"] = bridge["drug_target_of"].fillna("")

    # Join dream DEG status (ashr)
    bridge = bridge.merge(
        dream[["symbol", "is_tissue_deg"]],
        left_on="protein", right_on="symbol", how="left"
    ).drop(columns=["symbol"], errors="ignore")
    bridge["is_tissue_deg"] = bridge["is_tissue_deg"].fillna(False)

    # Rename boolean for output clarity
    bridge = bridge.rename(columns={"is_tissue_deg": "is_deg"})

    # ------------------------------------------------------------------
    # 3. Summary statistics
    # ------------------------------------------------------------------
    for target, grp in bridge.groupby("target"):
        n_deg = grp["is_deg"].sum()
        n_switch = (grp["f2_switch_class"] != "none").sum()
        n_drug = (grp["drug_target_of"] != "").sum()
        n_cc = grp["is_conserved"].sum() if "is_conserved" in grp else 0
        log.info(
            "  %s top-20: %d DEGs, %d in F2 switch, %d drug targets, %d conserved core",
            target, n_deg, n_switch, n_drug, n_cc,
        )

    # ------------------------------------------------------------------
    # 4. Fisher's exact enrichment
    # ------------------------------------------------------------------
    enrichment = run_enrichment(bridge, dream, all_proteins)

    # ------------------------------------------------------------------
    # 5. Drug overlaps (subset of bridge where drug_target_of is non-empty)
    # ------------------------------------------------------------------
    drug_overlaps = bridge[bridge["drug_target_of"] != ""][
        ["protein", "target", "shap_rank", "mean_abs_shap", "drug_target_of",
         "dream_logFC", "dream_padj", "is_deg", "f2_switch_class"]
    ].copy()
    log.info("Drug-target overlaps: %d rows", len(drug_overlaps))

    # ------------------------------------------------------------------
    # 6. Order columns for the bridge output
    # ------------------------------------------------------------------
    col_order = [
        "protein", "target", "shap_rank", "fold_stability", "mean_abs_shap",
        "dream_logFC", "dream_padj", "is_deg", "is_conserved",
        "f2_switch_class", "f2_inflection_logFC", "f2_inflection_padj",
        "dgidb_druggable", "opentargets_drug",
        "attribution_class", "layers_active", "sources_active",
        "drug_target_of",
    ]
    # Only keep columns that exist (graceful for missing atlas columns)
    col_order = [c for c in col_order if c in bridge.columns]
    bridge_out = bridge[col_order].sort_values(["target", "shap_rank"])

    # ------------------------------------------------------------------
    # 7. Save outputs
    # ------------------------------------------------------------------
    os.makedirs(RESULTS_MP, exist_ok=True)

    bridge_out.to_csv(OUT_BRIDGE, index=False)
    log.info("Saved bridge table: %s (%d rows)", OUT_BRIDGE, len(bridge_out))

    enrichment.to_csv(OUT_ENRICHMENT, index=False)
    log.info("Saved enrichment: %s (%d rows)", OUT_ENRICHMENT, len(enrichment))

    drug_overlaps.to_csv(OUT_DRUG_OVERLAPS, index=False)
    log.info("Saved drug overlaps: %s (%d rows)", OUT_DRUG_OVERLAPS, len(drug_overlaps))

    log.info("=" * 60)
    log.info("226b complete")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
