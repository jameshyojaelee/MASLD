#!/usr/bin/env python3
"""
127_spatial_validation.py
Spatial Validation of Progression Signatures using GSE192741 (Guilliams et al., 5 Visium slides).

Tests whether progression-derived signatures show spatial organization:
  1. Transition-specific genes enriched among spatially variable genes (SVGs)
  2. S2 (progressor) divergence genes show pericentral vs periportal spatial bias
  3. BayesPrism cell-type composition changes concordant with cell2location spatial patterns
  4. Disease-emergent SVGs overlap with transition and bifurcation gene sets
  5. Spatial domain distribution of progression signatures

Input:
  Spatial (Analysis/Spatial/results/):
    - svg/svgs_Healthy.csv, svgs_Steatotic.csv, differential_svgs.csv
    - zonation/deg_zonation_classification.csv
    - cell2location/spatial_cell_type_proportions.csv, cross_modality_validation.csv
    - domains/domain_markers.csv, domain_composition.csv
  Progression (results/progression/):
    - transition_programs.csv (gene_symbol column; adj.P.Val for significance)
    - divergence_genes.csv (gene_symbol column; padj for significance)
    - celltype_composition_changes.csv (cell_type, mean_prop_from/to, comp_padj)
    - celltype_attribution.csv (transition attribution summary)
    - bifurcation_summary.csv

Output (to results/progression/):
  - spatial_validation.csv          Per-gene spatial validation results
  - spatial_validation_summary.csv  Summary statistics across all tests

SLURM: io, 4 CPUs, 8G RAM, 48h
Env:   micromamba activate spatial
"""

import os
import sys
import logging
import numpy as np
import pandas as pd
from scipy import stats
from collections import OrderedDict

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SPATIAL_DIR = os.path.join(PROJECT_ROOT, "Analysis/Spatial/results")
PROG_DIR = os.path.join(
    PROJECT_ROOT,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression",
)
OUT_DIR = PROG_DIR  # outputs go alongside other progression results

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(OUT_DIR, "127_spatial_validation.log")),
    ],
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helper: Fisher's exact test (one-sided, greater)
# ---------------------------------------------------------------------------
def fisher_enrichment(query_genes, target_genes, universe_genes):
    """
    Fisher's exact test for enrichment of query_genes in target_genes.
    Returns (odds_ratio, pvalue, n_overlap, n_query, n_target, n_universe).
    """
    query = set(query_genes) & universe_genes
    target = set(target_genes) & universe_genes
    a = len(query & target)
    b = len(query - target)
    c = len(target - query)
    d = len(universe_genes - query - target)
    table = np.array([[a, b], [c, d]])
    odds_ratio, pval = stats.fisher_exact(table, alternative="greater")
    return odds_ratio, pval, a, len(query), len(target), len(universe_genes)


# ---------------------------------------------------------------------------
# 1. Load spatial data
# ---------------------------------------------------------------------------
log.info("Loading spatial data from %s", SPATIAL_DIR)

# SVGs per condition
svg_healthy = pd.read_csv(
    os.path.join(SPATIAL_DIR, "svg/svgs_Healthy.csv"), index_col=0
)
svg_masld = pd.read_csv(
    os.path.join(SPATIAL_DIR, "svg/svgs_Steatotic.csv"), index_col=0
)
diff_svg = pd.read_csv(
    os.path.join(SPATIAL_DIR, "svg/differential_svgs.csv"), index_col=0
)
svg_healthy_genes = set(svg_healthy[svg_healthy["svg"]].index)
svg_masld_genes = set(svg_masld[svg_masld["svg"]].index)
svg_union_genes = svg_healthy_genes | svg_masld_genes
svg_emergent = set(diff_svg[diff_svg["category"] == "disease_emergent_SVG"].index)
svg_universe = set(svg_healthy.index)  # all 3000 tested genes

log.info(
    "SVGs: %d healthy, %d MASLD, %d union, %d emergent, %d universe",
    len(svg_healthy_genes),
    len(svg_masld_genes),
    len(svg_union_genes),
    len(svg_emergent),
    len(svg_universe),
)

# Zonation classification (gene-level spatial data)
zonation = pd.read_csv(
    os.path.join(SPATIAL_DIR, "zonation/deg_zonation_classification.csv"), index_col=0
)
periportal_genes = set(
    zonation[zonation["zonation_class"] == "Periportal-enriched"]["gene"]
)
pericentral_genes = set(
    zonation[zonation["zonation_class"] == "Pericentral-enriched"]["gene"]
)
zonation_universe = set(zonation["gene"])
log.info(
    "Zonation: %d PP, %d PC, %d universe",
    len(periportal_genes),
    len(pericentral_genes),
    len(zonation_universe),
)

# Cell2location spatial proportions (per sample)
c2l_props = pd.read_csv(
    os.path.join(SPATIAL_DIR, "cell2location/spatial_cell_type_proportions.csv")
)
# Clean column names: extract cell type name from long column names
c2l_rename = {}
for col in c2l_props.columns:
    if "mu_fg_" in col:
        ct = col.split("mu_fg_")[1]
        c2l_rename[col] = ct
c2l_props = c2l_props.rename(columns=c2l_rename)

# Domain composition and markers
domain_comp = pd.read_csv(
    os.path.join(SPATIAL_DIR, "domains/domain_composition.csv"), index_col=0
)
domain_markers = pd.read_csv(
    os.path.join(SPATIAL_DIR, "domains/domain_markers.csv"), index_col=0
)
domain_markers_sig = domain_markers[domain_markers["pvals_adj"] < 0.05]

# Moran's I from MASLD condition for spatial autocorrelation strength
morans_masld = svg_masld["I"].to_dict()  # gene -> Moran's I

log.info("Spatial data loaded successfully")


# ---------------------------------------------------------------------------
# 2. Load progression data
# ---------------------------------------------------------------------------
log.info("Loading progression data from %s", PROG_DIR)

# Transition programs (long format: gene x transition)
trans_prog = pd.read_csv(os.path.join(PROG_DIR, "transition_programs.csv"))
log.info("Transition programs: %d rows, %d transitions", len(trans_prog), trans_prog["transition"].nunique())

# Extract significant transition genes per transition
PADJ_THRESH = 0.05
trans_sig = trans_prog[trans_prog["adj.P.Val"] < PADJ_THRESH]
transitions = sorted(trans_sig["transition"].unique())
transition_gene_sets = {}
for tr in transitions:
    genes = set(trans_sig[trans_sig["transition"] == tr]["gene_symbol"].dropna())
    # Remove Ensembl-only names (those that start with ENSG and have no symbol)
    genes = {g for g in genes if not g.startswith("ENSG")}
    transition_gene_sets[tr] = genes
    log.info("  %s: %d significant genes (padj<%g)", tr, len(genes), PADJ_THRESH)

# All transition genes (union across transitions) for summary stats
all_transition_genes = set()
for g in transition_gene_sets.values():
    all_transition_genes |= g

# Divergence genes (S1 vs S2 subtypes)
div_genes = pd.read_csv(os.path.join(PROG_DIR, "divergence_genes.csv"))
div_sig = div_genes[div_genes["padj"] < PADJ_THRESH].copy()
div_sig_symbols = set(div_sig["gene_symbol"].dropna())
div_sig_symbols = {g for g in div_sig_symbols if not g.startswith("ENSG")}
# S2-upregulated (progressor) vs S1-upregulated
div_s2_up = set(
    div_sig[div_sig["mean_diff"] < 0]["gene_symbol"].dropna()
)  # negative mean_diff = S2 > S1
div_s2_up = {g for g in div_s2_up if not g.startswith("ENSG")}
div_s1_up = set(
    div_sig[div_sig["mean_diff"] > 0]["gene_symbol"].dropna()
)
div_s1_up = {g for g in div_s1_up if not g.startswith("ENSG")}
log.info(
    "Divergence: %d sig genes (%d S2-up, %d S1-up)",
    len(div_sig_symbols),
    len(div_s2_up),
    len(div_s1_up),
)

# Cell-type composition changes across transitions
comp_changes = pd.read_csv(
    os.path.join(PROG_DIR, "celltype_composition_changes.csv")
)
log.info("Cell-type composition changes: %d rows", len(comp_changes))

# Peak transition genes (high tau = transition-specific)
trans_peaks = trans_prog.drop_duplicates(subset=["gene"]).copy()
high_tau = trans_peaks[trans_peaks["tau"] > 0.8]
high_tau_genes = set(high_tau["gene_symbol"].dropna())
high_tau_genes = {g for g in high_tau_genes if not g.startswith("ENSG")}
log.info("High tau (>0.8) genes: %d", len(high_tau_genes))

log.info("Progression data loaded successfully")


# ---------------------------------------------------------------------------
# 3. Test 1: Enrichment of transition genes among SVGs
# ---------------------------------------------------------------------------
log.info("=" * 70)
log.info("TEST 1: Transition gene enrichment among SVGs")
log.info("=" * 70)

test1_results = []

# For each transition, test enrichment in SVGs (healthy, MASLD, union, emergent)
for tr in transitions:
    query = transition_gene_sets[tr]
    for svg_name, svg_set in [
        ("SVG_healthy", svg_healthy_genes),
        ("SVG_MASLD", svg_masld_genes),
        ("SVG_union", svg_union_genes),
        ("SVG_emergent", svg_emergent),
    ]:
        OR, pval, n_ovl, n_q, n_t, n_u = fisher_enrichment(
            query, svg_set, svg_universe
        )
        test1_results.append(
            {
                "test": "transition_svg_enrichment",
                "transition": tr,
                "svg_type": svg_name,
                "odds_ratio": OR,
                "pvalue": pval,
                "n_overlap": n_ovl,
                "n_query_in_universe": n_q,
                "n_target": n_t,
                "n_universe": n_u,
            }
        )
        if pval < 0.05:
            log.info(
                "  %s in %s: OR=%.2f, p=%.2e, overlap=%d/%d",
                tr, svg_name, OR, pval, n_ovl, n_q,
            )

# Also test all transition genes (union) and high-tau genes
for name, gene_set in [
    ("all_transition_sig", all_transition_genes),
    ("high_tau_genes", high_tau_genes),
    ("divergence_sig", div_sig_symbols),
    ("divergence_S2_up", div_s2_up),
]:
    for svg_name, svg_set in [
        ("SVG_healthy", svg_healthy_genes),
        ("SVG_MASLD", svg_masld_genes),
        ("SVG_union", svg_union_genes),
        ("SVG_emergent", svg_emergent),
    ]:
        OR, pval, n_ovl, n_q, n_t, n_u = fisher_enrichment(
            gene_set, svg_set, svg_universe
        )
        test1_results.append(
            {
                "test": "gene_set_svg_enrichment",
                "transition": name,
                "svg_type": svg_name,
                "odds_ratio": OR,
                "pvalue": pval,
                "n_overlap": n_ovl,
                "n_query_in_universe": n_q,
                "n_target": n_t,
                "n_universe": n_u,
            }
        )
        if pval < 0.05:
            log.info(
                "  %s in %s: OR=%.2f, p=%.2e, overlap=%d/%d",
                name, svg_name, OR, pval, n_ovl, n_q,
            )

test1_df = pd.DataFrame(test1_results)

# Apply BH correction across all Test 1 p-values
from statsmodels.stats.multitest import multipletests as mt
_, test1_df["padj"], _, _ = mt(test1_df["pvalue"], method="fdr_bh")

log.info("Test 1 complete: %d enrichment tests", len(test1_df))


# ---------------------------------------------------------------------------
# 4. Test 2: Zonation bias of divergence genes
# ---------------------------------------------------------------------------
log.info("=" * 70)
log.info("TEST 2: Zonation bias of S2 divergence genes")
log.info("=" * 70)

test2_results = []

# Test: are divergence genes enriched in pericentral or periportal zones?
for name, gene_set in [
    ("divergence_sig", div_sig_symbols),
    ("divergence_S2_up", div_s2_up),
    ("divergence_S1_up", div_s1_up),
    ("all_transition_sig", all_transition_genes),
    ("high_tau_genes", high_tau_genes),
]:
    for zone_name, zone_genes in [
        ("Periportal", periportal_genes),
        ("Pericentral", pericentral_genes),
    ]:
        OR, pval, n_ovl, n_q, n_t, n_u = fisher_enrichment(
            gene_set, zone_genes, zonation_universe
        )
        test2_results.append(
            {
                "test": "zonation_enrichment",
                "gene_set": name,
                "zone": zone_name,
                "odds_ratio": OR,
                "pvalue": pval,
                "n_overlap": n_ovl,
                "n_query_in_universe": n_q,
                "n_zonated_genes": n_t,
                "n_universe": n_u,
            }
        )
        log.info(
            "  %s in %s: OR=%.2f, p=%.2e, overlap=%d/%d",
            name, zone_name, OR, pval, n_ovl, n_q,
        )

# Also test per-transition zonation bias
for tr in transitions:
    query = transition_gene_sets[tr]
    for zone_name, zone_genes in [
        ("Periportal", periportal_genes),
        ("Pericentral", pericentral_genes),
    ]:
        OR, pval, n_ovl, n_q, n_t, n_u = fisher_enrichment(
            query, zone_genes, zonation_universe
        )
        test2_results.append(
            {
                "test": "transition_zonation_enrichment",
                "gene_set": tr,
                "zone": zone_name,
                "odds_ratio": OR,
                "pvalue": pval,
                "n_overlap": n_ovl,
                "n_query_in_universe": n_q,
                "n_zonated_genes": n_t,
                "n_universe": n_u,
            }
        )
        if pval < 0.05:
            log.info(
                "  %s in %s: OR=%.2f, p=%.2e, overlap=%d",
                tr, zone_name, OR, pval, n_ovl,
            )

# Additional: Moran's I distribution comparison for divergence vs non-divergence genes
div_in_svg = div_sig_symbols & svg_universe
nondiv_in_svg = svg_universe - div_sig_symbols
morans_div = [morans_masld.get(g, np.nan) for g in div_in_svg]
morans_nondiv = [morans_masld.get(g, np.nan) for g in nondiv_in_svg]
morans_div = [x for x in morans_div if not np.isnan(x)]
morans_nondiv = [x for x in morans_nondiv if not np.isnan(x)]

if morans_div and morans_nondiv:
    mw_stat, mw_pval = stats.mannwhitneyu(
        morans_div, morans_nondiv, alternative="greater"
    )
    log.info(
        "  Moran's I: divergence genes median=%.3f vs non-div median=%.3f (MWU p=%.2e)",
        np.median(morans_div),
        np.median(morans_nondiv),
        mw_pval,
    )
    test2_results.append(
        {
            "test": "morans_i_comparison",
            "gene_set": "divergence_sig",
            "zone": "spatial_autocorrelation",
            "odds_ratio": np.median(morans_div) / max(np.median(morans_nondiv), 1e-10),
            "pvalue": mw_pval,
            "n_overlap": len(morans_div),
            "n_query_in_universe": len(div_in_svg),
            "n_zonated_genes": len(nondiv_in_svg),
            "n_universe": len(svg_universe),
        }
    )

test2_df = pd.DataFrame(test2_results)
if len(test2_df) > 0:
    _, test2_df["padj"], _, _ = mt(test2_df["pvalue"], method="fdr_bh")
log.info("Test 2 complete: %d tests", len(test2_df))


# ---------------------------------------------------------------------------
# 5. Test 3: Cell-type composition concordance
# ---------------------------------------------------------------------------
log.info("=" * 70)
log.info("TEST 3: BayesPrism vs cell2location cell-type concordance")
log.info("=" * 70)

test3_results = []

# Build cell-type name mapping between BayesPrism (progression) and cell2location
# BayesPrism names: Hepatocyte, Stellate, Endothelial, Macrophage, Cholangiocyte, etc.
# Cell2location names: Hepatocytes, Fibroblasts, Endothelial cells, Macrophages, Cholangiocytes, etc.
ct_map = {
    "Hepatocyte": "Hepatocytes",
    "Stellate": "Fibroblasts",  # stellate cells = fibroblasts in scRNA reference
    "Endothelial": "Endothelial cells",
    "Macrophage": "Macrophages",
    "Cholangiocyte": "Cholangiocytes",
    "Monocyte": "Mono+mono derived cells",
    "NK_cell": "Circulating NK/NKT",
    "B_cell": "B cells",
    "T_cell": "T cells",
    "DC": "cDC1s",  # approximate mapping
    "Neutrophil": "Neutrophils",
    "Plasma_cell": "Plasma cells",
}

# Get BayesPrism overall composition (mean across F0 baseline = earliest stage)
bp_f0 = comp_changes[comp_changes["transition"] == "F0_to_F1"][
    ["cell_type", "mean_prop_from"]
].copy()
bp_f0 = bp_f0.rename(columns={"mean_prop_from": "bayesprism_proportion"})

# Get cell2location overall composition (mean across healthy slides only)
# JBO018 and JBO022 are Healthy
c2l_healthy = c2l_props[c2l_props["sample_id"].isin(["JBO018", "JBO022"])]
c2l_mean = c2l_healthy.drop(columns=["sample_id"]).mean()

# Build comparison table
concordance_rows = []
bp_vals = []
c2l_vals = []
for bp_ct, c2l_ct in ct_map.items():
    bp_row = bp_f0[bp_f0["cell_type"] == bp_ct]
    if len(bp_row) == 0:
        continue
    bp_prop = bp_row["bayesprism_proportion"].values[0]
    c2l_prop = c2l_mean.get(c2l_ct, np.nan)
    if np.isnan(c2l_prop):
        continue
    concordance_rows.append(
        {
            "bayesprism_celltype": bp_ct,
            "cell2location_celltype": c2l_ct,
            "bayesprism_proportion": bp_prop,
            "cell2location_proportion": c2l_prop,
        }
    )
    bp_vals.append(bp_prop)
    c2l_vals.append(c2l_prop)

concordance_df = pd.DataFrame(concordance_rows)

# Spearman correlation between BayesPrism and cell2location proportions
if len(bp_vals) >= 3:
    rho, pval = stats.spearmanr(bp_vals, c2l_vals)
    log.info(
        "  BayesPrism vs cell2location proportion Spearman rho=%.3f, p=%.2e (n=%d cell types)",
        rho, pval, len(bp_vals),
    )
    test3_results.append(
        {
            "test": "celltype_proportion_correlation",
            "metric": "spearman_rho",
            "value": rho,
            "pvalue": pval,
            "n_celltypes": len(bp_vals),
            "description": "BayesPrism baseline (F0) vs cell2location healthy slides",
        }
    )

    # Also Pearson on log-transformed values (for proportions spanning orders of magnitude)
    log_bp = np.log10(np.array(bp_vals) + 1e-6)
    log_c2l = np.log10(np.array(c2l_vals) + 1e-6)
    r_pearson, p_pearson = stats.pearsonr(log_bp, log_c2l)
    log.info(
        "  Log10 Pearson r=%.3f, p=%.2e",
        r_pearson, p_pearson,
    )
    test3_results.append(
        {
            "test": "celltype_proportion_correlation",
            "metric": "pearson_r_log10",
            "value": r_pearson,
            "pvalue": p_pearson,
            "n_celltypes": len(bp_vals),
            "description": "BayesPrism baseline (F0) vs cell2location healthy slides (log10)",
        }
    )

# Test: direction of significant composition changes matches healthy→steatotic spatial shift
# Compare BayesPrism F0→F1 direction with healthy→steatotic spatial differences
c2l_steatotic = c2l_props[c2l_props["sample_id"].isin(["JBO014", "JBO015", "JBO019"])]
c2l_steatotic_mean = c2l_steatotic.drop(columns=["sample_id"]).mean()
c2l_healthy_mean = c2l_mean

direction_rows = []
bp_directions = []
c2l_directions = []
for bp_ct, c2l_ct in ct_map.items():
    bp_row = comp_changes[
        (comp_changes["transition"] == "F0_to_F1")
        & (comp_changes["cell_type"] == bp_ct)
    ]
    if len(bp_row) == 0:
        continue
    bp_change = bp_row["pct_change"].values[0]
    c2l_h = c2l_healthy_mean.get(c2l_ct, np.nan)
    c2l_s = c2l_steatotic_mean.get(c2l_ct, np.nan)
    if np.isnan(c2l_h) or np.isnan(c2l_s) or c2l_h == 0:
        continue
    c2l_change = (c2l_s - c2l_h) / c2l_h * 100

    direction_rows.append(
        {
            "bayesprism_celltype": bp_ct,
            "cell2location_celltype": c2l_ct,
            "bayesprism_pct_change": bp_change,
            "cell2location_pct_change": c2l_change,
            "direction_concordant": (bp_change > 0) == (c2l_change > 0),
        }
    )
    bp_directions.append(bp_change)
    c2l_directions.append(c2l_change)

direction_df = pd.DataFrame(direction_rows)
if len(direction_df) > 0:
    n_concordant = direction_df["direction_concordant"].sum()
    n_total = len(direction_df)
    # Binomial test: is concordance > 50% (chance)?
    binom_pval = stats.binomtest(n_concordant, n_total, 0.5, alternative="greater").pvalue
    log.info(
        "  Direction concordance: %d/%d (%.1f%%), binomial p=%.3f",
        n_concordant, n_total, 100 * n_concordant / n_total, binom_pval,
    )
    test3_results.append(
        {
            "test": "composition_direction_concordance",
            "metric": "fraction_concordant",
            "value": n_concordant / n_total,
            "pvalue": binom_pval,
            "n_celltypes": n_total,
            "description": f"{n_concordant}/{n_total} cell types change in same direction (BayesPrism F0->F1 vs spatial healthy->steatotic)",
        }
    )

    # Spearman of magnitudes
    if len(bp_directions) >= 3:
        rho_dir, pval_dir = stats.spearmanr(bp_directions, c2l_directions)
        log.info(
            "  Direction magnitude Spearman rho=%.3f, p=%.2e",
            rho_dir, pval_dir,
        )
        test3_results.append(
            {
                "test": "composition_magnitude_correlation",
                "metric": "spearman_rho",
                "value": rho_dir,
                "pvalue": pval_dir,
                "n_celltypes": len(bp_directions),
                "description": "Spearman of pct_change (BayesPrism F0->F1 vs spatial healthy->steatotic)",
            }
        )

test3_df = pd.DataFrame(test3_results)
log.info("Test 3 complete: %d tests", len(test3_df))


# ---------------------------------------------------------------------------
# 6. Test 4: Disease-emergent SVGs overlap with progression signatures
# ---------------------------------------------------------------------------
log.info("=" * 70)
log.info("TEST 4: Disease-emergent SVG overlap with progression signatures")
log.info("=" * 70)

test4_results = []

# Disease-emergent SVGs: genes that become spatially variable in MASLD but not healthy
for name, gene_set in [
    ("all_transition_sig", all_transition_genes),
    ("high_tau_genes", high_tau_genes),
    ("divergence_sig", div_sig_symbols),
    ("divergence_S2_up", div_s2_up),
]:
    OR, pval, n_ovl, n_q, n_t, n_u = fisher_enrichment(
        gene_set, svg_emergent, svg_universe
    )
    test4_results.append(
        {
            "test": "disease_emergent_svg_enrichment",
            "gene_set": name,
            "odds_ratio": OR,
            "pvalue": pval,
            "n_overlap": n_ovl,
            "n_query_in_universe": n_q,
            "n_emergent_svgs": n_t,
            "n_universe": n_u,
        }
    )
    overlap_genes = gene_set & svg_emergent
    log.info(
        "  %s in emergent SVGs: OR=%.2f, p=%.2e, overlap=%d (%s)",
        name, OR, pval, n_ovl,
        ", ".join(sorted(overlap_genes)[:10]) if overlap_genes else "none",
    )

# Per-transition test for emergent SVGs
for tr in transitions:
    query = transition_gene_sets[tr]
    OR, pval, n_ovl, n_q, n_t, n_u = fisher_enrichment(
        query, svg_emergent, svg_universe
    )
    test4_results.append(
        {
            "test": "transition_emergent_svg_enrichment",
            "gene_set": tr,
            "odds_ratio": OR,
            "pvalue": pval,
            "n_overlap": n_ovl,
            "n_query_in_universe": n_q,
            "n_emergent_svgs": n_t,
            "n_universe": n_u,
        }
    )

test4_df = pd.DataFrame(test4_results)
if len(test4_df) > 0:
    _, test4_df["padj"], _, _ = mt(test4_df["pvalue"], method="fdr_bh")
log.info("Test 4 complete: %d tests", len(test4_df))


# ---------------------------------------------------------------------------
# 7. Test 5: Spatial domain distribution of progression signatures
# ---------------------------------------------------------------------------
log.info("=" * 70)
log.info("TEST 5: Spatial domain enrichment of progression signatures")
log.info("=" * 70)

test5_results = []

# Get marker gene sets per spatial domain
domain_gene_sets = {}
for dom in sorted(domain_markers_sig["group"].unique()):
    genes = set(domain_markers_sig[domain_markers_sig["group"] == dom]["names"])
    domain_gene_sets[f"domain_{dom}"] = genes
    log.info("  Domain %s: %d marker genes", dom, len(genes))

# Universe = all domain marker genes tested
domain_universe = set(domain_markers["names"])
log.info("  Domain marker universe: %d genes", len(domain_universe))

for name, gene_set in [
    ("all_transition_sig", all_transition_genes),
    ("high_tau_genes", high_tau_genes),
    ("divergence_sig", div_sig_symbols),
    ("divergence_S2_up", div_s2_up),
    ("divergence_S1_up", div_s1_up),
]:
    for dom_name, dom_genes in domain_gene_sets.items():
        OR, pval, n_ovl, n_q, n_t, n_u = fisher_enrichment(
            gene_set, dom_genes, domain_universe
        )
        test5_results.append(
            {
                "test": "domain_enrichment",
                "gene_set": name,
                "spatial_domain": dom_name,
                "odds_ratio": OR,
                "pvalue": pval,
                "n_overlap": n_ovl,
                "n_query_in_universe": n_q,
                "n_domain_markers": n_t,
                "n_universe": n_u,
            }
        )
        if pval < 0.05:
            log.info(
                "  %s in %s: OR=%.2f, p=%.2e, overlap=%d",
                name, dom_name, OR, pval, n_ovl,
            )

test5_df = pd.DataFrame(test5_results)
if len(test5_df) > 0:
    _, test5_df["padj"], _, _ = mt(test5_df["pvalue"], method="fdr_bh")
log.info("Test 5 complete: %d tests", len(test5_df))


# ---------------------------------------------------------------------------
# 8. Build per-gene spatial validation table
# ---------------------------------------------------------------------------
log.info("=" * 70)
log.info("Building per-gene spatial validation table")
log.info("=" * 70)

# Collect all progression genes with spatial annotations
all_prog_genes = all_transition_genes | div_sig_symbols

per_gene_rows = []
for gene in sorted(all_prog_genes):
    row = {"gene_symbol": gene}

    # SVG status
    row["is_svg_healthy"] = gene in svg_healthy_genes
    row["is_svg_masld"] = gene in svg_masld_genes
    row["is_svg_union"] = gene in svg_union_genes
    row["is_disease_emergent_svg"] = gene in svg_emergent

    # Moran's I in MASLD
    row["morans_i_masld"] = morans_masld.get(gene, np.nan)
    row["morans_i_healthy"] = (
        svg_healthy.loc[gene, "I"] if gene in svg_healthy.index else np.nan
    )

    # Zonation
    zon_row = zonation[zonation["gene"] == gene]
    if len(zon_row) > 0:
        row["zonation_class"] = zon_row["zonation_class"].values[0]
        row["zonation_spearman_rho"] = zon_row["spearman_rho"].values[0]
    else:
        row["zonation_class"] = np.nan
        row["zonation_spearman_rho"] = np.nan

    # Progression membership
    row["is_transition_sig"] = gene in all_transition_genes
    row["is_divergence_sig"] = gene in div_sig_symbols
    row["is_S2_up"] = gene in div_s2_up
    row["is_high_tau"] = gene in high_tau_genes

    # Which transitions significant in
    gene_transitions = []
    for tr in transitions:
        if gene in transition_gene_sets[tr]:
            gene_transitions.append(tr)
    row["n_transitions_sig"] = len(gene_transitions)
    row["transitions_sig"] = ";".join(gene_transitions) if gene_transitions else ""

    # Domain membership
    for dom_name, dom_genes in domain_gene_sets.items():
        row[f"is_{dom_name}_marker"] = gene in dom_genes

    per_gene_rows.append(row)

per_gene_df = pd.DataFrame(per_gene_rows)
log.info("Per-gene table: %d genes, %d columns", len(per_gene_df), len(per_gene_df.columns))

# Summary stats
n_svg = per_gene_df["is_svg_union"].sum()
n_emergent = per_gene_df["is_disease_emergent_svg"].sum()
n_zonated = per_gene_df["zonation_class"].notna().sum()
log.info(
    "  %d/%d progression genes are SVGs (%.1f%%)",
    n_svg, len(per_gene_df), 100 * n_svg / len(per_gene_df),
)
log.info(
    "  %d/%d are disease-emergent SVGs (%.1f%%)",
    n_emergent, len(per_gene_df), 100 * n_emergent / len(per_gene_df),
)


# ---------------------------------------------------------------------------
# 9. Compile summary
# ---------------------------------------------------------------------------
log.info("=" * 70)
log.info("Compiling summary")
log.info("=" * 70)

summary_rows = []

# Test 1 summary: best SVG enrichment per transition
for tr in transitions:
    tr_test = test1_df[
        (test1_df["transition"] == tr) & (test1_df["svg_type"] == "SVG_MASLD")
    ]
    if len(tr_test) > 0:
        row = tr_test.iloc[0]
        summary_rows.append(
            {
                "test": "T1_transition_SVG_enrichment",
                "category": tr,
                "odds_ratio": row["odds_ratio"],
                "pvalue": row["pvalue"],
                "padj": row["padj"],
                "n_overlap": row["n_overlap"],
                "n_tested": row["n_query_in_universe"],
                "description": f"{tr} transition genes in MASLD SVGs",
            }
        )

# Test 1 summary: aggregate gene sets
for name in ["all_transition_sig", "high_tau_genes", "divergence_sig", "divergence_S2_up"]:
    sub = test1_df[
        (test1_df["transition"] == name) & (test1_df["svg_type"] == "SVG_MASLD")
    ]
    if len(sub) > 0:
        row = sub.iloc[0]
        summary_rows.append(
            {
                "test": "T1_geneset_SVG_enrichment",
                "category": name,
                "odds_ratio": row["odds_ratio"],
                "pvalue": row["pvalue"],
                "padj": row["padj"],
                "n_overlap": row["n_overlap"],
                "n_tested": row["n_query_in_universe"],
                "description": f"{name} in MASLD SVGs",
            }
        )

# Test 2 summary: zonation
for _, row in test2_df.iterrows():
    if row["test"] in ("zonation_enrichment", "morans_i_comparison"):
        summary_rows.append(
            {
                "test": f"T2_{row['test']}",
                "category": f"{row['gene_set']}_{row['zone']}",
                "odds_ratio": row["odds_ratio"],
                "pvalue": row["pvalue"],
                "padj": row.get("padj", np.nan),
                "n_overlap": row["n_overlap"],
                "n_tested": row["n_query_in_universe"],
                "description": f"{row['gene_set']} in {row['zone']} zone",
            }
        )

# Test 3 summary
for _, row in test3_df.iterrows():
    summary_rows.append(
        {
            "test": f"T3_{row['test']}",
            "category": row["metric"],
            "odds_ratio": row["value"],
            "pvalue": row["pvalue"],
            "padj": np.nan,
            "n_overlap": np.nan,
            "n_tested": row["n_celltypes"],
            "description": row["description"],
        }
    )

# Test 4 summary
for _, row in test4_df.iterrows():
    summary_rows.append(
        {
            "test": "T4_emergent_SVG",
            "category": row["gene_set"],
            "odds_ratio": row["odds_ratio"],
            "pvalue": row["pvalue"],
            "padj": row.get("padj", np.nan),
            "n_overlap": row["n_overlap"],
            "n_tested": row["n_query_in_universe"],
            "description": f"{row['gene_set']} in disease-emergent SVGs",
        }
    )

# Test 5 summary: only significant domain enrichments
for _, row in test5_df.iterrows():
    if row["pvalue"] < 0.05:
        summary_rows.append(
            {
                "test": "T5_domain_enrichment",
                "category": f"{row['gene_set']}_{row['spatial_domain']}",
                "odds_ratio": row["odds_ratio"],
                "pvalue": row["pvalue"],
                "padj": row.get("padj", np.nan),
                "n_overlap": row["n_overlap"],
                "n_tested": row["n_query_in_universe"],
                "description": f"{row['gene_set']} markers in {row['spatial_domain']}",
            }
        )

summary_df = pd.DataFrame(summary_rows)

# Overall validation score: fraction of tests that are significant
n_sig_tests = summary_df[summary_df["pvalue"] < 0.05].shape[0]
n_total_tests = len(summary_df)
log.info(
    "Overall: %d/%d summary tests significant at p<0.05 (%.1f%%)",
    n_sig_tests, n_total_tests, 100 * n_sig_tests / max(n_total_tests, 1),
)


# ---------------------------------------------------------------------------
# 10. Save outputs
# ---------------------------------------------------------------------------
log.info("Saving outputs to %s", OUT_DIR)

per_gene_df.to_csv(os.path.join(OUT_DIR, "spatial_validation.csv"), index=False)
log.info("  spatial_validation.csv: %d rows x %d cols", *per_gene_df.shape)

summary_df.to_csv(os.path.join(OUT_DIR, "spatial_validation_summary.csv"), index=False)
log.info("  spatial_validation_summary.csv: %d rows", len(summary_df))

# Also save detailed test results for transparency
test1_df.to_csv(os.path.join(OUT_DIR, "spatial_validation_test1_svg_enrichment.csv"), index=False)
test2_df.to_csv(os.path.join(OUT_DIR, "spatial_validation_test2_zonation.csv"), index=False)
test3_df.to_csv(os.path.join(OUT_DIR, "spatial_validation_test3_celltype.csv"), index=False)
test4_df.to_csv(os.path.join(OUT_DIR, "spatial_validation_test4_emergent.csv"), index=False)
test5_df.to_csv(os.path.join(OUT_DIR, "spatial_validation_test5_domains.csv"), index=False)

# Save concordance and direction tables as well
concordance_df.to_csv(
    os.path.join(OUT_DIR, "spatial_validation_celltype_concordance.csv"), index=False
)
direction_df.to_csv(
    os.path.join(OUT_DIR, "spatial_validation_direction_concordance.csv"), index=False
)

log.info("=" * 70)
log.info("Script 127 complete")
log.info("=" * 70)
