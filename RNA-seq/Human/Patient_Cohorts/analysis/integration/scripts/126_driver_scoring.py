#!/usr/bin/env python3
"""
126_driver_scoring.py
Multi-Evidence Driver Scoring: Unify all evidence into per-gene, per-transition scores.

For each gene x transition, computes a composite driver score from:
  1. Transition specificity (tau from Script 115)
  2. Genetic causality (MR/TWAS/COLOC from existing atlas)
  3. Cell-type resolution (BayesPrism attribution from Script 118-119)
  4. Bifurcation involvement (divergence from Script 117)
  5. Pseudotime transport cost (from Script 116)
  6. Cross-species concordance (from existing atlas)
  7. Communication rewiring (CCC L-R DE + pseudotime dynamics from Scripts 124-125)

Integration approach: Bayesian-inspired scoring with positive control calibration.
Known drug targets (THRB, NR1H4, PPARs, GLP1R) serve as validation benchmarks.

Gene name harmonization:
  - Scripts 115/116 output Ensembl IDs (versioned, e.g. ENSG00000000003.17)
  - Script 117 has both Ensembl IDs and gene_symbol column
  - Script 119 uses gene symbols (from BayesPrism)
  - Atlas uses human_symbol + unversioned ensembl_id
  - Strategy: build Ensembl-to-symbol map from atlas, convert everything to symbols

Input:
  - results/progression/transition_programs.csv (Script 115; long format, Ensembl IDs)
  - results/progression/transport_gene_costs.csv (Script 116; NAS transitions, Ensembl IDs)
  - results/progression/divergence_genes.csv (Script 117; Ensembl + gene_symbol)
  - results/progression/celltype_transition_programs.csv (Script 119; gene symbols)
  - results/progression/ccc_transition_de.csv (Script 124; L-R pair DE per transition)
  - results/progression/ccc_pseudotime_dynamics.csv (Script 125; L-R pair pseudotime dynamics)
  - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (human_symbol + ensembl_id)

Output (to results/progression/):
  - driver_scores.csv          (gene x transition x evidence components)
  - driver_summary.csv         (top 50 drivers per transition)
  - therapeutic_roadmap.csv    (druggable targets by transition)

SLURM: io, 4 CPUs, 16G RAM, 48h
Env:   micromamba activate spatial
"""
import os
import sys
import logging
import numpy as np
import pandas as pd
from scipy.stats import rankdata

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# -- Paths -----------------------------------------------------------------
BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INTEG = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PROG = os.path.join(INTEG, "results/progression")
RESULTS = os.path.join(BASE, "RNA-seq/results")
OUTDIR = PROG
os.makedirs(OUTDIR, exist_ok=True)

# -- Known positive controls -----------------------------------------------
POSITIVE_CONTROLS = {
    "THRB":   {"drug": "resmetirom", "stage": "approved", "transition": "F2_to_F3"},
    "NR1H4":  {"drug": "obeticholic acid", "stage": "phase3", "transition": "F2_to_F3"},
    "PPARA":  {"drug": "elafibranor", "stage": "phase3", "transition": "F1_to_F2"},
    "PPARD":  {"drug": "seladelpar", "stage": "phase3", "transition": "F1_to_F2"},
    "GLP1R":  {"drug": "semaglutide", "stage": "phase3", "transition": "F0_to_F1"},
    "LIPE":   {"drug": "lanifibranor (pan-PPAR)", "stage": "phase3", "transition": "F1_to_F2"},
    "FGF21":  {"drug": "pegozafermin", "stage": "phase3", "transition": "F1_to_F2"},
    "ACC1":   {"drug": "firsocostat", "stage": "phase2", "transition": "F0_to_F1"},
    "DGAT2":  {"drug": "ervogastat", "stage": "phase2", "transition": "F0_to_F1"},
    "ASK1":   {"drug": "selonsertib", "stage": "failed", "transition": "F2_to_F3"},
}

# Fibrosis transitions (from Script 115)
FIB_TRANSITIONS = ["F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4"]
# NAS transitions (from Script 116 transport)
NAS_TRANSITIONS = ["NAS_0_1_to_NAS_2_4", "NAS_2_4_to_NAS_5", "NAS_5_to_NAS_6_8"]


def build_ensembl_to_symbol(atlas):
    """Build Ensembl ID -> gene symbol mapping from atlas.

    Atlas has unversioned ensembl_id. Transition data has versioned IDs.
    Returns dict mapping both versioned and unversioned IDs to symbols.
    """
    mapping = {}
    for _, row in atlas[["human_symbol", "ensembl_id"]].dropna().iterrows():
        symbol = str(row["human_symbol"]).strip()
        ens = str(row["ensembl_id"]).strip()
        if symbol and ens and symbol != "nan" and ens != "nan":
            mapping[ens] = symbol
    log.info(f"  Ensembl-to-symbol mapping: {len(mapping)} genes")
    return mapping


def ensembl_to_symbol(gene_id, mapping):
    """Convert a (possibly versioned) Ensembl ID to symbol."""
    gene_id = str(gene_id).strip()
    # Try exact match first (for unversioned IDs)
    if gene_id in mapping:
        return mapping[gene_id]
    # Strip version and retry
    base = gene_id.split(".")[0]
    if base in mapping:
        return mapping[base]
    return None


def harmonize_gene_column(df, gene_col, mapping, source_name):
    """Add a 'gene_symbol' column to df by mapping Ensembl IDs to symbols.

    Returns df with added 'gene_symbol' column, dropping rows that fail to map.
    """
    original_len = len(df)
    # Strip version from Ensembl IDs for lookup
    df = df.copy()
    df["_ensembl_base"] = df[gene_col].astype(str).str.split(".").str[0]
    df["gene_symbol"] = df["_ensembl_base"].map(mapping)
    mapped = df["gene_symbol"].notna().sum()
    log.info(f"  {source_name}: {mapped}/{original_len} genes mapped to symbols "
             f"({mapped/max(original_len,1)*100:.1f}%)")
    df = df[df["gene_symbol"].notna()].drop(columns=["_ensembl_base"])
    return df


def load_evidence(ens_to_sym):
    """Load all evidence sources, harmonizing gene names to symbols."""
    evidence = {}

    # 1. Transition programs (Script 115) -- long format, Ensembl IDs
    #    Columns: gene, logFC, AveExpr, t, P.Value, adj.P.Val, z.std,
    #             transition, n_lower, n_higher, n_datasets, tau,
    #             peak_transition, max_abs_t, peak_logFC, peak_padj
    tp_file = os.path.join(PROG, "transition_programs.csv")
    if os.path.exists(tp_file):
        tp = pd.read_csv(tp_file)
        tp = harmonize_gene_column(tp, "gene", ens_to_sym, "Transition programs")
        evidence["transitions"] = tp
        log.info(f"  Transition programs: {len(tp)} rows, "
                 f"transitions: {sorted(tp['transition'].unique())}")
    else:
        log.warning("  Transition programs not found (Script 115 still running?)")

    # 2. Transport gene costs (Script 116) -- NAS transitions, Ensembl IDs
    #    Columns: gene, transition, transport_cost, transport_direction, cost_rank
    tc_file = os.path.join(PROG, "transport_gene_costs.csv")
    if os.path.exists(tc_file):
        tc = pd.read_csv(tc_file)
        tc = harmonize_gene_column(tc, "gene", ens_to_sym, "Transport costs")
        evidence["transport"] = tc
        log.info(f"  Transport costs: {len(tc)} rows, "
                 f"transitions: {sorted(tc['transition'].unique())}")
    else:
        log.warning("  Transport costs not found")

    # 3. Divergence genes (Script 117) -- has both gene (Ensembl) and gene_symbol
    #    Columns: gene, pval, cohens_d, mean_S1, mean_S2, mean_diff, padj, gene_symbol
    div_file = os.path.join(PROG, "divergence_genes.csv")
    if os.path.exists(div_file):
        div = pd.read_csv(div_file)
        # Script 117 already has gene_symbol column
        if "gene_symbol" in div.columns:
            div = div.rename(columns={"gene_symbol": "gene_symbol"})  # no-op, just explicit
            mapped = div["gene_symbol"].notna().sum()
            log.info(f"  Divergence genes: {mapped}/{len(div)} have gene_symbol")
            div = div[div["gene_symbol"].notna()]
        else:
            div = harmonize_gene_column(div, "gene", ens_to_sym, "Divergence genes")
        evidence["divergence"] = div
        log.info(f"  Divergence genes: {len(div)} rows")
    else:
        log.warning("  Divergence genes not found")

    # 4. Cell-type transitions (Script 119) -- already gene symbols
    #    Columns: gene, cell_type, transition, transition_type, lfc, pval, padj, n_from, n_to
    ct_file = os.path.join(PROG, "celltype_transition_programs.csv")
    if os.path.exists(ct_file):
        ct = pd.read_csv(ct_file)
        # Already uses gene symbols -- rename to gene_symbol for consistency
        ct = ct.rename(columns={"gene": "gene_symbol"})
        evidence["celltype"] = ct
        log.info(f"  Cell-type transitions: {len(ct)} rows, "
                 f"{ct['cell_type'].nunique()} cell types, "
                 f"transitions: {sorted(ct['transition'].unique())}")
    else:
        log.warning("  Cell-type transition programs not found")

    # 5. Multi-evidence atlas -- human_symbol is the gene symbol column
    #    Key columns: human_symbol, ensembl_id, dream_logFC, dream_padj,
    #    mr_*, twas_*, coloc_*, dgidb_druggable, is_conserved, etc.
    atlas_file = os.path.join(RESULTS, "multi_evidence/multi_evidence_atlas.csv")
    if os.path.exists(atlas_file):
        atlas = pd.read_csv(atlas_file)
        # Rename for consistency
        atlas = atlas.rename(columns={"human_symbol": "gene_symbol"})
        evidence["atlas"] = atlas
        log.info(f"  Multi-evidence atlas: {len(atlas)} rows, "
                 f"{len(atlas.columns)} cols")
    else:
        log.warning("  Multi-evidence atlas not found")

    # 6. CCC transition DE (Script 124) -- L-R pair DE per fibrosis transition
    #    Columns: ligand, receptor, axis, lr_pair, transition, lfc, padj, pathway
    #    Gene names are already symbols (ligand/receptor columns).
    ccc_de_file = os.path.join(PROG, "ccc_transition_de.csv")
    if os.path.exists(ccc_de_file):
        ccc_de = pd.read_csv(ccc_de_file)
        evidence["ccc_de"] = ccc_de
        log.info(f"  CCC transition DE: {len(ccc_de)} rows, "
                 f"transitions: {sorted(ccc_de['transition'].unique())}, "
                 f"{ccc_de['ligand'].nunique()} ligands, "
                 f"{ccc_de['receptor'].nunique()} receptors")
    else:
        log.warning("  CCC transition DE not found (Scripts 124-125 not run?)")

    # 7. CCC pseudotime dynamics (Script 125) -- L-R pair pseudotime correlations
    #    Columns: axis, lr_pair, ligand, receptor, spearman_rho, spearman_padj,
    #             trend_adj, cluster, canonical_pathway
    #    Gene names are already symbols (ligand/receptor columns).
    ccc_dyn_file = os.path.join(PROG, "ccc_pseudotime_dynamics.csv")
    if os.path.exists(ccc_dyn_file):
        ccc_dyn = pd.read_csv(ccc_dyn_file)
        evidence["ccc_dynamics"] = ccc_dyn
        log.info(f"  CCC pseudotime dynamics: {len(ccc_dyn)} rows, "
                 f"monotonic (padj<0.05): "
                 f"{len(ccc_dyn[(ccc_dyn['spearman_padj'] < 0.05) & (ccc_dyn['trend_adj'] != 'stable')])}")
    else:
        log.warning("  CCC pseudotime dynamics not found")

    return evidence


def compute_ccc_score(gene_symbol_series, transition, evidence):
    """Compute CCC communication rewiring score for a set of genes at a transition.

    For each gene: count significant L-R interactions where the gene is a ligand
    or receptor (padj < 0.1, |LFC| > 0.25). Normalize to [0, 1] by dividing by
    the max across all genes. Apply a 20% boost for genes whose L-R pairs show
    monotonic pseudotime trends (from ccc_pseudotime_dynamics).

    Returns Series aligned to gene_symbol_series with CCC score in [0, 1].
    """
    genes = gene_symbol_series.values
    n = len(genes)

    # Build set of monotonically changing L-R pairs for the pseudotime boost
    monotonic_lr_pairs = set()
    if "ccc_dynamics" in evidence:
        dyn = evidence["ccc_dynamics"]
        mono = dyn[(dyn["spearman_padj"] < 0.05) & (dyn["trend_adj"] != "stable")]
        monotonic_lr_pairs = set(mono["lr_pair"].unique())

    if "ccc_de" not in evidence:
        return pd.Series(np.zeros(n), index=gene_symbol_series.index)

    ccc_de = evidence["ccc_de"]

    # Filter to this transition and significance thresholds
    sig = ccc_de[(ccc_de["transition"] == transition) &
                 (ccc_de["padj"] < 0.1) &
                 (ccc_de["lfc"].abs() > 0.25)].copy()

    if len(sig) == 0:
        return pd.Series(np.zeros(n), index=gene_symbol_series.index)

    # For each significant L-R pair, check if the pair is also monotonic
    sig["is_monotonic"] = sig["lr_pair"].isin(monotonic_lr_pairs)

    # Count per gene: number of sig LR interactions (with monotonic boost)
    # Each interaction counts as 1.0 base + 0.2 if monotonic = 1.2 max per pair
    gene_scores = {}
    for _, row in sig.iterrows():
        base = 1.0
        boost = 0.2 if row["is_monotonic"] else 0.0
        pair_score = base + boost
        for gene in [row["ligand"], row["receptor"]]:
            if pd.notna(gene):
                gene_scores[gene] = gene_scores.get(gene, 0.0) + pair_score

    # Map to input gene series
    score_values = np.array([gene_scores.get(g, 0.0) for g in genes])

    # Normalize to [0, 1]
    max_score = score_values.max()
    if max_score > 0:
        score_values = score_values / max_score

    return pd.Series(score_values, index=gene_symbol_series.index)


def compute_driver_scores(evidence):
    """Compute composite driver scores per gene x transition.

    Two scoring axes:
      - Fibrosis transitions (F0_to_F1 ... F3_to_F4): from Script 115 + 117 + 119 + atlas
      - NAS transitions: from Script 116 transport + 119 celltype (if available)
    All merged on gene_symbol.
    """
    all_scores = []

    # Build gene universe from all symbol-harmonized sources
    gene_sets = []
    if "transitions" in evidence:
        gene_sets.append(set(evidence["transitions"]["gene_symbol"].dropna()))
    if "atlas" in evidence:
        gene_sets.append(set(evidence["atlas"]["gene_symbol"].dropna()))
    if "transport" in evidence:
        gene_sets.append(set(evidence["transport"]["gene_symbol"].dropna()))
    if "divergence" in evidence:
        gene_sets.append(set(evidence["divergence"]["gene_symbol"].dropna()))
    if "celltype" in evidence:
        gene_sets.append(set(evidence["celltype"]["gene_symbol"].dropna()))
    if "ccc_de" in evidence:
        ccc_genes = (set(evidence["ccc_de"]["ligand"].dropna()) |
                     set(evidence["ccc_de"]["receptor"].dropna()))
        gene_sets.append(ccc_genes)

    if not gene_sets:
        log.error("No evidence sources loaded!")
        return pd.DataFrame()

    all_genes = sorted(set.union(*gene_sets))
    log.info(f"\n  Gene universe (symbols): {len(all_genes)}")

    # ---------- Fibrosis transitions ----------
    for trans in FIB_TRANSITIONS:
        log.info(f"\n  Scoring {trans}...")
        scores = pd.DataFrame({"gene_symbol": all_genes})

        # -- S1: Transition specificity (Script 115) -----------------------
        if "transitions" in evidence:
            tp = evidence["transitions"]
            tp_trans = tp[tp["transition"] == trans].copy()
            if len(tp_trans) > 0:
                # Use actual column names from dream output
                tp_sub = tp_trans[["gene_symbol", "logFC", "adj.P.Val", "tau"]].copy()
                tp_sub = tp_sub.rename(columns={"logFC": "trans_logFC",
                                                "adj.P.Val": "trans_padj"})
                # De-duplicate: keep row with smallest padj per gene
                tp_sub = tp_sub.sort_values("trans_padj").drop_duplicates("gene_symbol", keep="first")
                scores = scores.merge(tp_sub, on="gene_symbol", how="left")
                scores["s1_transition_specificity"] = scores["tau"].fillna(0)
                scores["s1_is_transition_deg"] = (scores["trans_padj"] < 0.1).astype(int).fillna(0)
            else:
                log.warning(f"    No transition data for {trans}")
                scores["s1_transition_specificity"] = 0.0
                scores["s1_is_transition_deg"] = 0
                scores["tau"] = np.nan
                scores["trans_logFC"] = np.nan
                scores["trans_padj"] = np.nan
        else:
            scores["s1_transition_specificity"] = 0.0
            scores["s1_is_transition_deg"] = 0
            scores["tau"] = np.nan
            scores["trans_logFC"] = np.nan
            scores["trans_padj"] = np.nan

        # -- S2: Genetic causality (from atlas) ----------------------------
        if "atlas" in evidence:
            atlas = evidence["atlas"]
            # Select causal evidence columns
            # mr_sig removed 2026-04-22 — MR ditched from paper.
            causal_cols = ["twas_pval", "coloc_pp4", "broadaway_coloc_pp4",
                           "n_coloc_sources", "causal_methods_sig", "causal_robustness",
                           "ctwas_pip", "hyprcoloc_posterior"]
            avail_causal = [c for c in causal_cols if c in atlas.columns]

            if avail_causal:
                atlas_sub = atlas[["gene_symbol"] + avail_causal].copy()
                atlas_sub = atlas_sub.drop_duplicates("gene_symbol", keep="first")
                scores = scores.merge(atlas_sub, on="gene_symbol", how="left")

                # Composite causal score: normalize each indicator to [0,1] and average
                causal_indicators = []

                # MR indicator removed 2026-04-22 — MR ditched from paper.
                # TWAS significant (binary from pval)
                if "twas_pval" in scores.columns:
                    causal_indicators.append(
                        (pd.to_numeric(scores["twas_pval"], errors="coerce") < 0.05)
                        .astype(float).fillna(0)
                    )

                # COLOC (continuous, already 0-1)
                for col in ["coloc_pp4", "broadaway_coloc_pp4"]:
                    if col in scores.columns:
                        causal_indicators.append(
                            pd.to_numeric(scores[col], errors="coerce").fillna(0).clip(0, 1)
                        )

                # n_coloc_sources (normalize by max)
                if "n_coloc_sources" in scores.columns:
                    vals = pd.to_numeric(scores["n_coloc_sources"], errors="coerce").fillna(0)
                    max_val = vals.max()
                    if max_val > 0:
                        causal_indicators.append(vals / max_val)

                # causal_methods_sig (normalize by max)
                if "causal_methods_sig" in scores.columns:
                    vals = pd.to_numeric(scores["causal_methods_sig"], errors="coerce").fillna(0)
                    max_val = vals.max()
                    if max_val > 0:
                        causal_indicators.append(vals / max_val)

                # cTWAS PIP (continuous, already 0-1)
                if "ctwas_pip" in scores.columns:
                    causal_indicators.append(
                        pd.to_numeric(scores["ctwas_pip"], errors="coerce").fillna(0).clip(0, 1)
                    )

                # HyPrColoc posterior (continuous, already 0-1)
                if "hyprcoloc_posterior" in scores.columns:
                    causal_indicators.append(
                        pd.to_numeric(scores["hyprcoloc_posterior"], errors="coerce").fillna(0).clip(0, 1)
                    )

                if causal_indicators:
                    scores["s2_genetic_causality"] = pd.concat(causal_indicators, axis=1).mean(axis=1)
                else:
                    scores["s2_genetic_causality"] = 0.0
            else:
                scores["s2_genetic_causality"] = 0.0

            # Also pull dream LFC and cross-species info
            dream_cols = ["dream_logFC", "dream_padj", "is_conserved"]
            avail_dream = [c for c in dream_cols if c in atlas.columns]
            if avail_dream:
                dream_sub = atlas[["gene_symbol"] + avail_dream].copy()
                dream_sub = dream_sub.drop_duplicates("gene_symbol", keep="first")
                # Only merge columns not already present
                existing = set(scores.columns)
                new_cols = [c for c in avail_dream if c not in existing]
                if new_cols:
                    scores = scores.merge(dream_sub[["gene_symbol"] + new_cols],
                                          on="gene_symbol", how="left")
        else:
            scores["s2_genetic_causality"] = 0.0

        # -- S3: Transport cost (Script 116) -- NAS transitions, not fibrosis
        # Map fibrosis transitions to NAS transitions for a rough correspondence:
        #   F0_to_F1 ~ NAS_0_1_to_NAS_2_4 (early disease)
        #   F1_to_F2 ~ NAS_2_4_to_NAS_5 (mid disease)
        #   F2_to_F3 ~ NAS_5_to_NAS_6_8 (late disease)
        #   F3_to_F4 ~ NAS_5_to_NAS_6_8 (late disease, no better match)
        fib_to_nas_map = {
            "F0_to_F1": "NAS_0_1_to_NAS_2_4",
            "F1_to_F2": "NAS_2_4_to_NAS_5",
            "F2_to_F3": "NAS_5_to_NAS_6_8",
            "F3_to_F4": "NAS_5_to_NAS_6_8",
        }
        if "transport" in evidence:
            tc = evidence["transport"]
            nas_trans = fib_to_nas_map.get(trans)
            if nas_trans:
                tc_trans = tc[tc["transition"] == nas_trans].copy()
            else:
                tc_trans = pd.DataFrame()

            if len(tc_trans) > 0:
                tc_sub = tc_trans[["gene_symbol", "transport_cost"]].copy()
                tc_sub = tc_sub.drop_duplicates("gene_symbol", keep="first")
                scores = scores.merge(tc_sub, on="gene_symbol", how="left")
                valid = scores["transport_cost"].notna()
                if valid.sum() > 0:
                    scores.loc[valid, "s3_transport_cost"] = (
                        rankdata(scores.loc[valid, "transport_cost"].abs().values)
                        / valid.sum()
                    )
                    scores["s3_transport_cost"] = scores["s3_transport_cost"].fillna(0)
                else:
                    scores["s3_transport_cost"] = 0.0
            else:
                scores["s3_transport_cost"] = 0.0
                scores["transport_cost"] = np.nan
        else:
            scores["s3_transport_cost"] = 0.0
            scores["transport_cost"] = np.nan

        # -- S4: Bifurcation divergence (Script 117) -----------------------
        # Divergence is subtype-level (S1 vs S2), not transition-specific.
        # Apply same divergence score to all transitions.
        if "divergence" in evidence:
            div = evidence["divergence"]
            div_sub = div[["gene_symbol", "cohens_d", "padj"]].copy()
            div_sub = div_sub.rename(columns={"cohens_d": "div_cohens_d",
                                               "padj": "div_padj"})
            div_sub = div_sub.drop_duplicates("gene_symbol", keep="first")
            scores = scores.merge(div_sub, on="gene_symbol", how="left")
            scores["s4_bifurcation_divergence"] = (
                (scores["div_padj"].fillna(1) < 0.05).astype(float)
                * scores["div_cohens_d"].abs().fillna(0)
            )
        else:
            scores["s4_bifurcation_divergence"] = 0.0
            scores["div_cohens_d"] = np.nan
            scores["div_padj"] = np.nan

        # -- S5: Cell-type resolution (Script 119) -------------------------
        if "celltype" in evidence:
            ct = evidence["celltype"]
            # Filter for this fibrosis transition + significant
            ct_trans = ct[(ct["transition"] == trans) & (ct["padj"] < 0.1)].copy()
            if len(ct_trans) > 0:
                n_celltypes_total = ct[ct["transition"] == trans]["cell_type"].nunique()
                ct_counts = (ct_trans.groupby("gene_symbol")["cell_type"]
                             .nunique()
                             .reset_index(name="n_celltypes_sig"))
                scores = scores.merge(ct_counts, on="gene_symbol", how="left")
                scores["n_celltypes_sig"] = scores["n_celltypes_sig"].fillna(0)
                scores["s5_celltype_resolution"] = (
                    scores["n_celltypes_sig"] / max(n_celltypes_total, 1)
                )
            else:
                scores["n_celltypes_sig"] = 0
                scores["s5_celltype_resolution"] = 0.0
        else:
            scores["n_celltypes_sig"] = 0
            scores["s5_celltype_resolution"] = 0.0

        # -- S6: Cross-species concordance (from atlas) --------------------
        if "atlas" in evidence and "is_conserved" in scores.columns:
            # is_conserved may be R-style "TRUE"/"FALSE" strings
            cc = scores["is_conserved"].copy()
            if cc.dtype == object:
                cc = cc.map({"TRUE": 1.0, "FALSE": 0.0}).fillna(0.0)
            else:
                cc = cc.fillna(0).astype(float)
            scores["s6_cross_species"] = cc
        else:
            scores["s6_cross_species"] = 0.0

        # -- S7: Communication rewiring (CCC, Scripts 124-125) -------------
        scores["s7_communication_rewiring"] = compute_ccc_score(
            scores["gene_symbol"], trans, evidence
        )
        n_ccc = (scores["s7_communication_rewiring"] > 0).sum()
        log.info(f"    S7 CCC: {n_ccc} genes with non-zero communication rewiring")

        # -- Composite score -----------------------------------------------
        # Weights reduced proportionally from original to accommodate s7 (~0.10)
        weights = {
            "s1_transition_specificity": 0.18,  # Transition-specific tau
            "s1_is_transition_deg":      0.13,  # DEG at this transition
            "s2_genetic_causality":      0.18,  # Genetic causal evidence
            "s3_transport_cost":         0.09,  # OT transport cost (NAS proxy)
            "s4_bifurcation_divergence": 0.13,  # Subtype divergence
            "s5_celltype_resolution":    0.09,  # Cell-type specificity
            "s6_cross_species":          0.10,  # Cross-species conservation
            "s7_communication_rewiring": 0.10,  # CCC L-R rewiring
        }

        scores["driver_score"] = sum(
            scores[col].fillna(0).astype(float) * w for col, w in weights.items()
        )

        # Rank (1 = best)
        scores["driver_rank"] = rankdata(-scores["driver_score"].values)
        scores["transition"] = trans
        scores["transition_type"] = "fibrosis"

        all_scores.append(scores)
        n_nonzero = (scores["driver_score"] > 0).sum()
        log.info(f"    {trans}: {n_nonzero} genes with non-zero score")

    # ---------- NAS transitions (transport-centric) ----------
    for trans in NAS_TRANSITIONS:
        log.info(f"\n  Scoring {trans} (NAS)...")
        scores = pd.DataFrame({"gene_symbol": all_genes})

        # S3: Transport cost (direct match for NAS transitions)
        if "transport" in evidence:
            tc = evidence["transport"]
            tc_trans = tc[tc["transition"] == trans].copy()
            if len(tc_trans) > 0:
                tc_sub = tc_trans[["gene_symbol", "transport_cost"]].copy()
                tc_sub = tc_sub.drop_duplicates("gene_symbol", keep="first")
                scores = scores.merge(tc_sub, on="gene_symbol", how="left")
                valid = scores["transport_cost"].notna()
                if valid.sum() > 0:
                    scores.loc[valid, "s3_transport_cost"] = (
                        rankdata(scores.loc[valid, "transport_cost"].abs().values)
                        / valid.sum()
                    )
                    scores["s3_transport_cost"] = scores["s3_transport_cost"].fillna(0)
                else:
                    scores["s3_transport_cost"] = 0.0
            else:
                scores["s3_transport_cost"] = 0.0
                scores["transport_cost"] = np.nan
        else:
            scores["s3_transport_cost"] = 0.0
            scores["transport_cost"] = np.nan

        # S2: Genetic causality (same for all transitions)
        if "atlas" in evidence:
            atlas = evidence["atlas"]
            # mr_sig removed 2026-04-22 — MR ditched from paper.
            causal_cols = ["twas_pval", "coloc_pp4", "n_coloc_sources",
                           "causal_methods_sig", "ctwas_pip", "hyprcoloc_posterior"]
            avail_causal = [c for c in causal_cols if c in atlas.columns]
            if avail_causal:
                atlas_sub = atlas[["gene_symbol"] + avail_causal].copy()
                atlas_sub = atlas_sub.drop_duplicates("gene_symbol", keep="first")
                scores = scores.merge(atlas_sub, on="gene_symbol", how="left")

                causal_indicators = []
                # MR indicator removed 2026-04-22 — MR ditched from paper.
                if "twas_pval" in scores.columns:
                    causal_indicators.append(
                        (pd.to_numeric(scores["twas_pval"], errors="coerce") < 0.05)
                        .astype(float).fillna(0)
                    )
                for col in ["coloc_pp4"]:
                    if col in scores.columns:
                        causal_indicators.append(
                            pd.to_numeric(scores[col], errors="coerce").fillna(0).clip(0, 1)
                        )
                if "n_coloc_sources" in scores.columns:
                    vals = pd.to_numeric(scores["n_coloc_sources"], errors="coerce").fillna(0)
                    max_val = vals.max()
                    if max_val > 0:
                        causal_indicators.append(vals / max_val)
                if "causal_methods_sig" in scores.columns:
                    vals = pd.to_numeric(scores["causal_methods_sig"], errors="coerce").fillna(0)
                    max_val = vals.max()
                    if max_val > 0:
                        causal_indicators.append(vals / max_val)
                if "ctwas_pip" in scores.columns:
                    causal_indicators.append(
                        pd.to_numeric(scores["ctwas_pip"], errors="coerce").fillna(0).clip(0, 1)
                    )
                if "hyprcoloc_posterior" in scores.columns:
                    causal_indicators.append(
                        pd.to_numeric(scores["hyprcoloc_posterior"], errors="coerce").fillna(0).clip(0, 1)
                    )

                if causal_indicators:
                    scores["s2_genetic_causality"] = pd.concat(causal_indicators, axis=1).mean(axis=1)
                else:
                    scores["s2_genetic_causality"] = 0.0
            else:
                scores["s2_genetic_causality"] = 0.0
        else:
            scores["s2_genetic_causality"] = 0.0

        # S4: Bifurcation divergence (same for all)
        if "divergence" in evidence:
            div = evidence["divergence"]
            div_sub = div[["gene_symbol", "cohens_d", "padj"]].copy()
            div_sub = div_sub.rename(columns={"cohens_d": "div_cohens_d",
                                               "padj": "div_padj"})
            div_sub = div_sub.drop_duplicates("gene_symbol", keep="first")
            scores = scores.merge(div_sub, on="gene_symbol", how="left")
            scores["s4_bifurcation_divergence"] = (
                (scores["div_padj"].fillna(1) < 0.05).astype(float)
                * scores["div_cohens_d"].abs().fillna(0)
            )
        else:
            scores["s4_bifurcation_divergence"] = 0.0

        # S5: Cell-type resolution (NAS transitions exist in celltype data)
        if "celltype" in evidence:
            ct = evidence["celltype"]
            ct_trans = ct[(ct["transition"] == trans) & (ct["padj"] < 0.1)].copy()
            if len(ct_trans) > 0:
                n_celltypes_total = ct[ct["transition"] == trans]["cell_type"].nunique()
                ct_counts = (ct_trans.groupby("gene_symbol")["cell_type"]
                             .nunique()
                             .reset_index(name="n_celltypes_sig"))
                scores = scores.merge(ct_counts, on="gene_symbol", how="left")
                scores["n_celltypes_sig"] = scores["n_celltypes_sig"].fillna(0)
                scores["s5_celltype_resolution"] = (
                    scores["n_celltypes_sig"] / max(n_celltypes_total, 1)
                )
            else:
                scores["n_celltypes_sig"] = 0
                scores["s5_celltype_resolution"] = 0.0
        else:
            scores["n_celltypes_sig"] = 0
            scores["s5_celltype_resolution"] = 0.0

        # NAS transitions have no transition_programs tau or cross-species
        scores["s1_transition_specificity"] = 0.0
        scores["s1_is_transition_deg"] = 0
        scores["s6_cross_species"] = 0.0
        scores["tau"] = np.nan
        scores["trans_logFC"] = np.nan
        scores["trans_padj"] = np.nan

        # S7: Communication rewiring (map NAS transitions to closest fibrosis)
        # CCC transition DE has fibrosis transitions only, so map:
        #   NAS_0_1_to_NAS_2_4 ~ F0_to_F1 (early disease)
        #   NAS_2_4_to_NAS_5   ~ F1_to_F2 (mid disease)
        #   NAS_5_to_NAS_6_8   ~ F2_to_F3 (late disease)
        nas_to_fib_map = {
            "NAS_0_1_to_NAS_2_4": "F0_to_F1",
            "NAS_2_4_to_NAS_5":   "F1_to_F2",
            "NAS_5_to_NAS_6_8":   "F2_to_F3",
        }
        fib_proxy = nas_to_fib_map.get(trans, "F2_to_F3")
        scores["s7_communication_rewiring"] = compute_ccc_score(
            scores["gene_symbol"], fib_proxy, evidence
        )
        n_ccc = (scores["s7_communication_rewiring"] > 0).sum()
        log.info(f"    S7 CCC (via {fib_proxy}): {n_ccc} genes with non-zero rewiring")

        # Composite (re-weighted: transport is primary for NAS, CCC added)
        nas_weights = {
            "s3_transport_cost":         0.30,  # Primary signal for NAS
            "s2_genetic_causality":      0.22,
            "s4_bifurcation_divergence": 0.18,
            "s5_celltype_resolution":    0.18,
            "s7_communication_rewiring": 0.12,  # CCC L-R rewiring
        }

        scores["driver_score"] = sum(
            scores[col].fillna(0).astype(float) * w for col, w in nas_weights.items()
        )

        scores["driver_rank"] = rankdata(-scores["driver_score"].values)
        scores["transition"] = trans
        scores["transition_type"] = "NAS"

        all_scores.append(scores)
        n_nonzero = (scores["driver_score"] > 0).sum()
        log.info(f"    {trans}: {n_nonzero} genes with non-zero score")

    result = pd.concat(all_scores, ignore_index=True)

    # Ensure consistent column ordering in output
    # Identify score columns vs metadata
    core_cols = ["gene_symbol", "transition", "transition_type",
                 "driver_score", "driver_rank"]
    component_cols = [c for c in result.columns if c.startswith("s1_") or
                      c.startswith("s2_") or c.startswith("s3_") or
                      c.startswith("s4_") or c.startswith("s5_") or
                      c.startswith("s6_") or c.startswith("s7_")]
    detail_cols = ["tau", "trans_logFC", "trans_padj", "transport_cost",
                   "div_cohens_d", "div_padj", "n_celltypes_sig"]
    detail_cols = [c for c in detail_cols if c in result.columns]
    other_cols = [c for c in result.columns
                  if c not in core_cols + component_cols + detail_cols]
    ordered = core_cols + sorted(component_cols) + detail_cols + sorted(other_cols)
    ordered = [c for c in ordered if c in result.columns]
    result = result[ordered]

    return result


def validate_against_controls(driver_df):
    """Check if known drug targets rank well."""
    log.info("\n=== Positive Control Validation ===")

    results = []
    for gene, info in POSITIVE_CONTROLS.items():
        trans = info["transition"]
        sub = driver_df[(driver_df["gene_symbol"] == gene) &
                        (driver_df["transition"] == trans)]
        if len(sub) > 0:
            rank = sub["driver_rank"].values[0]
            score = sub["driver_score"].values[0]
            total = len(driver_df[driver_df["transition"] == trans])
            pct = rank / total * 100
            log.info(f"  {gene:8s} ({info['drug']:20s}, {info['stage']:8s}): "
                     f"rank {rank:.0f}/{total} (top {pct:.1f}%), score={score:.3f}")
            results.append({"gene": gene, "rank": rank, "total": total,
                            "pct": pct, "score": score, "found": True})
        else:
            log.info(f"  {gene:8s} ({info['drug']:20s}): NOT FOUND in {trans}")
            results.append({"gene": gene, "rank": np.nan, "total": np.nan,
                            "pct": np.nan, "score": np.nan, "found": False})

    found = sum(1 for r in results if r["found"])
    top10 = sum(1 for r in results if r["found"] and r["pct"] <= 10)
    log.info(f"\n  Summary: {found}/{len(POSITIVE_CONTROLS)} controls found, "
             f"{top10} in top 10%")
    return results


def generate_therapeutic_roadmap(driver_df, evidence):
    """Identify druggable targets per transition."""
    log.info("\n=== Therapeutic Roadmap ===")

    # Get druggable gene set from atlas
    druggable = set()
    if "atlas" in evidence:
        atlas = evidence["atlas"]
        if "dgidb_druggable" in atlas.columns:
            druggable = set(atlas[atlas["dgidb_druggable"].notna()]["gene_symbol"])
            log.info(f"  Druggable genes in atlas: {len(druggable)}")
        if "opentargets_drug" in atlas.columns:
            ot_druggable = set(atlas[atlas["opentargets_drug"].notna()]["gene_symbol"])
            druggable = druggable | ot_druggable
            log.info(f"  Total druggable (DGIdb + OpenTargets): {len(druggable)}")

    roadmap = []
    for trans in driver_df["transition"].unique():
        sub = driver_df[driver_df["transition"] == trans].sort_values("driver_rank")
        top50 = sub.head(50)

        n_druggable = sum(1 for g in top50["gene_symbol"] if g in druggable)
        log.info(f"\n  {trans}: top 50 drivers, {n_druggable} druggable")

        for _, row in top50.iterrows():
            is_drug = row["gene_symbol"] in druggable
            is_control = row["gene_symbol"] in POSITIVE_CONTROLS
            entry = {
                "gene_symbol": row["gene_symbol"],
                "transition": row["transition"],
                "transition_type": row.get("transition_type", "fibrosis"),
                "driver_rank": row["driver_rank"],
                "driver_score": row["driver_score"],
                "is_druggable": is_drug,
                "is_positive_control": is_control,
            }
            # Add component scores
            for col in row.index:
                if col.startswith("s") and "_" in col and col[1].isdigit():
                    entry[col] = row[col]
            roadmap.append(entry)

            if is_drug or is_control:
                log.info(f"    rank {row['driver_rank']:.0f}: {row['gene_symbol']} "
                         f"(score={row['driver_score']:.3f})"
                         f"{' [DRUG TARGET]' if is_drug else ''}"
                         f"{' [+CTRL]' if is_control else ''}")

    return pd.DataFrame(roadmap)


def main():
    log.info("=== 126: Multi-Evidence Driver Scoring ===\n")

    # -- Build Ensembl-to-symbol mapping from atlas ------------------------
    log.info("Building Ensembl-to-symbol mapping from atlas...")
    atlas_file = os.path.join(RESULTS, "multi_evidence/multi_evidence_atlas.csv")
    if not os.path.exists(atlas_file):
        log.error(f"Atlas not found: {atlas_file}")
        sys.exit(1)
    atlas_raw = pd.read_csv(atlas_file, usecols=["human_symbol", "ensembl_id"])
    ens_to_sym = build_ensembl_to_symbol(atlas_raw)
    del atlas_raw

    # -- Load all evidence -------------------------------------------------
    log.info("\nLoading evidence sources...")
    evidence = load_evidence(ens_to_sym)

    if not evidence:
        log.error("No evidence sources found!")
        sys.exit(1)

    log.info(f"\n  Evidence sources loaded: {sorted(evidence.keys())}")

    # -- Compute driver scores ---------------------------------------------
    log.info("\nComputing driver scores...")
    driver_df = compute_driver_scores(evidence)
    log.info(f"  Total scored: {len(driver_df)} (gene x transition entries)")

    # -- Validate ----------------------------------------------------------
    validate_against_controls(driver_df)

    # -- Therapeutic roadmap -----------------------------------------------
    roadmap = generate_therapeutic_roadmap(driver_df, evidence)

    # -- Save --------------------------------------------------------------
    log.info("\n=== Saving Results ===")

    # Full scores
    out_scores = os.path.join(OUTDIR, "driver_scores.csv")
    driver_df.to_csv(out_scores, index=False)
    log.info(f"  Saved driver_scores.csv ({len(driver_df)} rows)")

    # Summary: top 50 per transition
    summary_rows = []
    for trans in driver_df["transition"].unique():
        top = driver_df[driver_df["transition"] == trans].nsmallest(50, "driver_rank")
        summary_rows.append(top)
    if summary_rows:
        summary = pd.concat(summary_rows)
        out_summary = os.path.join(OUTDIR, "driver_summary.csv")
        summary.to_csv(out_summary, index=False)
        log.info(f"  Saved driver_summary.csv ({len(summary)} rows)")
    else:
        log.warning("  No summary rows to save")

    # Therapeutic roadmap
    out_roadmap = os.path.join(OUTDIR, "therapeutic_roadmap.csv")
    roadmap.to_csv(out_roadmap, index=False)
    log.info(f"  Saved therapeutic_roadmap.csv ({len(roadmap)} rows)")

    log.info("\n=== 126: COMPLETE ===")


if __name__ == "__main__":
    main()
