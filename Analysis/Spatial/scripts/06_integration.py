#!/usr/bin/env python3
"""
06_integration.py — Integrate spatial results with bulk multi-evidence atlas.

Adds spatial evidence columns (zonation, deconvolution, SVGs, domains) to
the existing multi-evidence atlas. Performs enrichment tests of spatial
findings against key gene sets.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config, load_multi_evidence_atlas,
    load_conserved, load_deconv_scores, build_ensembl_to_symbol_map,
    resolve_path, save_csv, print_header,
)


def load_spatial_results():
    """Load all spatial result CSVs."""
    results = {}

    # Zonation classification
    zon_path = RESULTS_DIR / "zonation" / "deg_zonation_classification.csv"
    if zon_path.exists():
        results["zonation"] = pd.read_csv(zon_path, index_col=0)

    # SVGs (use MASLD if available, else first condition)
    for pattern in RESULTS_DIR.glob("svg/svgs_*.csv"):
        results["svg"] = pd.read_csv(pattern, index_col=0)
        break  # Use first found

    # Differential SVGs
    diff_svg = RESULTS_DIR / "svg" / "differential_svgs.csv"
    if diff_svg.exists():
        results["diff_svg"] = pd.read_csv(diff_svg, index_col=0)

    # Cell type deconvolution
    hep_val = RESULTS_DIR / "cell2location" / "hep_intrinsic_validation.csv"
    if hep_val.exists():
        results["hep_validation"] = pd.read_csv(hep_val, index_col=0)

    # Spatial domains
    dom_path = RESULTS_DIR / "domains" / "spatial_domains.csv"
    if dom_path.exists():
        results["domains"] = pd.read_csv(dom_path, index_col=0)

    # Co-expression modules
    for pattern in RESULTS_DIR.glob("coexpression/modules_*.csv"):
        results["modules"] = pd.read_csv(pattern, index_col=0)
        break

    # Phase 4: Zonation-aware DE (09a)
    zon_de = RESULTS_DIR / "zonation" / "zonation_aware_de.csv"
    if zon_de.exists():
        results["zonation_de"] = pd.read_csv(zon_de, index_col=0)

    # Phase 4: GWAS overlay (09b)
    gwas_zon = RESULTS_DIR / "gwas_spatial" / "coloc_zonation_mapping.csv"
    if gwas_zon.exists():
        results["gwas_spatial"] = pd.read_csv(gwas_zon, index_col=0)

    # Phase 4: Niche targets (09c)
    niche = RESULTS_DIR / "drug_repurposing" / "spatial_niche_targets.csv"
    if niche.exists():
        results["niche_targets"] = pd.read_csv(niche, index_col=0)

    # Phase 4: Cross-modality (09d)
    cm = RESULTS_DIR / "cell2location" / "cross_modality_validation.csv"
    if cm.exists():
        results["cross_modality"] = pd.read_csv(cm, index_col=0)

    # --- Govaere2026 spatial outputs (GeoMx regional + CosMx cell-type DE) ---
    govaere_dir = RESULTS_DIR / "govaere2026"
    gov_specs = {
        "govaere_geomx_sh_vs_pt":  ("geomx_de_sh_vs_pt.csv", "gene_symbol"),
        "govaere_geomx_sh_vs_ls":  ("geomx_de_sh_vs_ls.csv", "gene_symbol"),
        "govaere_cosmx_hep_mash":  ("cosmx_de_Hepatocyte_MASH_vs_noMASH.csv", "gene"),
        # KC pre-subcluster (loader 41 leiden res=0.5 merged MetMac+TransMac into KC)
        "govaere_cosmx_kc_mash":   ("cosmx_de_KC_MASH_vs_noMASH.csv",         "gene"),
        # Sub-clustered macrophage lineage (from 41b at leiden res=1.5,
        # following Govaere 2026 Methods page 14 "resolution = 1-2").
        "govaere_cosmx_metmac_mash":   ("cosmx_de_MetMac_MASH_vs_noMASH.csv",         "gene"),
        "govaere_cosmx_transmac_mash": ("cosmx_de_TransMac_MASH_vs_noMASH.csv",       "gene"),
        "govaere_cosmx_kcpost_mash":   ("cosmx_de_KC_post_subcluster_MASH_vs_noMASH.csv", "gene"),
        "govaere_cosmx_premac_mash":   ("cosmx_de_preMac_MASH_vs_noMASH.csv",         "gene"),
        "govaere_cosmx_mono_mash":     ("cosmx_de_Monocyte_MASH_vs_noMASH.csv",       "gene"),
    }
    for key, (fname, gene_col) in gov_specs.items():
        path = govaere_dir / fname
        if path.exists():
            df = pd.read_csv(path)
            # Stash gene column under a uniform name for downstream merge
            df = df.rename(columns={gene_col: "_gene_key"})
            results[key] = df

    # Govaere2026 signature panels (already gene-keyed, wide)
    sig_path = (PROJECT_ROOT / "RNA-seq" / "results" / "govaere2026"
                / "govaere2026_signatures_wide.tsv")
    if sig_path.exists():
        results["govaere_signatures"] = pd.read_csv(sig_path, sep="\t")

    # Multi-cohort spatial consensus (produced by 44_spatial_consensus.py)
    consensus_path = RESULTS_DIR / "integration" / "spatial_consensus.csv"
    if consensus_path.exists():
        results["spatial_consensus"] = pd.read_csv(consensus_path)

    return results


def build_spatial_evidence(atlas, spatial_results):
    """Add spatial evidence columns to multi-evidence atlas."""
    symbol_col = "human_symbol" if "human_symbol" in atlas.columns else ("symbol" if "symbol" in atlas.columns else atlas.columns[0])

    # Zonation
    if "zonation" in spatial_results:
        zon = spatial_results["zonation"]
        gene_col = "gene" if "gene" in zon.columns else zon.index.name
        if gene_col and gene_col in zon.columns:
            zon_map = zon.set_index(gene_col) if gene_col in zon.columns else zon
        else:
            zon_map = zon
        for col in ["zonation_class", "spearman_rho"]:
            if col in zon_map.columns:
                atlas[f"spatial_{col}"] = atlas[symbol_col].map(
                    zon_map[col].to_dict()
                )

    # SVGs
    if "svg" in spatial_results:
        svg = spatial_results["svg"]
        atlas["spatial_morans_i"] = atlas[symbol_col].map(
            svg["I"].to_dict() if "I" in svg.columns else {}
        )
        atlas["spatial_is_svg"] = atlas[symbol_col].map(
            svg["svg"].to_dict() if "svg" in svg.columns else {}
        )

    # Differential SVGs
    if "diff_svg" in spatial_results:
        diff = spatial_results["diff_svg"]
        atlas["spatial_svg_category"] = atlas[symbol_col].map(
            diff["category"].to_dict() if "category" in diff.columns else {}
        )

    # Hepatocyte validation
    if "hep_validation" in spatial_results:
        hep = spatial_results["hep_validation"]
        gene_col = "gene" if "gene" in hep.columns else hep.index.name
        if gene_col and gene_col in hep.columns:
            hep_map = hep.set_index(gene_col)
        else:
            hep_map = hep
        for col in ["spatial_fc", "validated", "wilcoxon_padj_bh"]:
            if col in hep_map.columns:
                atlas[f"spatial_hep_{col}"] = atlas[symbol_col].map(
                    hep_map[col].to_dict()
                )

    # Co-expression modules
    if "modules" in spatial_results:
        mod = spatial_results["modules"]
        if "Module" in mod.columns:
            atlas["spatial_coexpr_module"] = atlas[symbol_col].map(
                mod["Module"].to_dict()
            )

    # --- Phase 4 columns ---

    # Zonation-aware DE (09a)
    if "zonation_de" in spatial_results:
        zde = spatial_results["zonation_de"]
        gene_col = "gene" if "gene" in zde.columns else zde.index.name
        if gene_col and gene_col in zde.columns:
            zde_map = zde.set_index(gene_col)
        else:
            zde_map = zde
        if "classification" in zde_map.columns:
            atlas["spatial_zone_de_class"] = atlas[symbol_col].map(
                zde_map["classification"].to_dict()
            )
        if "n_sig_bins" in zde_map.columns:
            atlas["spatial_zone_de_n_sig_bins"] = atlas[symbol_col].map(
                zde_map["n_sig_bins"].to_dict()
            )

    # GWAS overlay (09b)
    if "gwas_spatial" in spatial_results:
        gwas = spatial_results["gwas_spatial"]
        gene_col = "gene" if "gene" in gwas.columns else gwas.index.name
        if gene_col and gene_col in gwas.columns:
            gwas_map = gwas.set_index(gene_col)
        else:
            gwas_map = gwas
        if "zonation_class" in gwas_map.columns:
            atlas["spatial_coloc_zone"] = atlas[symbol_col].map(
                gwas_map["zonation_class"].to_dict()
            )

    # Niche targets (09c)
    if "niche_targets" in spatial_results:
        niche = spatial_results["niche_targets"]
        gene_col = "gene" if "gene" in niche.columns else niche.index.name
        if gene_col and gene_col in niche.columns:
            niche_map = niche.set_index(gene_col)
        else:
            niche_map = niche
        if "disease_niche_domain" in niche_map.columns:
            atlas["spatial_niche_domain"] = atlas[symbol_col].map(
                niche_map["disease_niche_domain"].to_dict()
            )
        if "niche_fc" in niche_map.columns:
            atlas["spatial_niche_fc"] = atlas[symbol_col].map(
                niche_map["niche_fc"].to_dict()
            )

    # Cross-modality validation (09d)
    if "cross_modality" in spatial_results:
        cm = spatial_results["cross_modality"]
        # This is aggregated per-cell-type, map by cell_type proportions
        # Add as metadata note, not per-gene
        pass

    # --- Govaere2026 spatial DE (GeoMx + CosMx) ---
    # GeoMx column schema: gene_symbol, logFC, AveExpr, t_or_z_stat, pval,
    #                      padj_bh, n_segs_group1, n_segs_group2
    # CosMx column schema: gene, logfoldchange, pval, pval_adj, score,
    #                      pct_nz_MASH, pct_nz_no_MASH, n_cells_MASH,
    #                      n_cells_no_MASH, cell_type
    govaere_map = [
        # (results-key,                  output_prefix,                            logfc_col,        padj_col)
        ("govaere_geomx_sh_vs_pt",       "spatial_govaere2026_geomx_sh_vs_pt",     "logFC",         "padj_bh"),
        ("govaere_geomx_sh_vs_ls",       "spatial_govaere2026_geomx_sh_vs_ls",     "logFC",         "padj_bh"),
        ("govaere_cosmx_hep_mash",       "spatial_govaere2026_cosmx_hep_mash",     "logfoldchange", "pval_adj"),
        ("govaere_cosmx_kc_mash",        "spatial_govaere2026_cosmx_kc_mash",      "logfoldchange", "pval_adj"),
        # Macrophage sub-cluster CosMx outputs (41b) — adds the GPNMB+/MetMac
        # axis directly. Paper's primary macrophage finding.
        ("govaere_cosmx_metmac_mash",    "spatial_govaere2026_cosmx_metmac_mash",  "logfoldchange", "pval_adj"),
        ("govaere_cosmx_transmac_mash",  "spatial_govaere2026_cosmx_transmac_mash","logfoldchange", "pval_adj"),
        ("govaere_cosmx_kcpost_mash",    "spatial_govaere2026_cosmx_kcpost_mash",  "logfoldchange", "pval_adj"),
        ("govaere_cosmx_premac_mash",    "spatial_govaere2026_cosmx_premac_mash",  "logfoldchange", "pval_adj"),
        ("govaere_cosmx_mono_mash",      "spatial_govaere2026_cosmx_mono_mash",    "logfoldchange", "pval_adj"),
    ]
    for key, prefix, lfc_col, padj_col in govaere_map:
        if key not in spatial_results:
            continue
        df = spatial_results[key]
        if "_gene_key" not in df.columns or lfc_col not in df.columns or padj_col not in df.columns:
            print(f"  WARNING: {key} missing required columns; skipping")
            continue
        # Aggregate duplicate symbols (some panels have isoform rows). Keep
        # the row with the smallest padj per gene (most significant evidence).
        df_agg = (df.sort_values(padj_col)
                    .drop_duplicates(subset="_gene_key", keep="first"))
        lfc_map  = dict(zip(df_agg["_gene_key"], df_agg[lfc_col]))
        padj_map = dict(zip(df_agg["_gene_key"], df_agg[padj_col]))
        atlas[f"{prefix}_logfc"] = atlas[symbol_col].map(lfc_map)
        atlas[f"{prefix}_padj"]  = atlas[symbol_col].map(padj_map)

    # --- Govaere2026 signature panels (wide, already gene-keyed) ---
    if "govaere_signatures" in spatial_results:
        sig = spatial_results["govaere_signatures"]
        if "human_symbol" in sig.columns:
            # Aggregate any duplicate symbols (defensive); use first occurrence
            sig = sig.drop_duplicates(subset="human_symbol", keep="first")
            sig_cols = [c for c in sig.columns if c.startswith("signature_govaere2026_")]
            sig_indexed = sig.set_index("human_symbol")[sig_cols]
            # Idempotency: drop any pre-existing signature_govaere2026_* cols from
            # the atlas before the merge (otherwise pd.merge appends _x/_y suffixes
            # when 06_integration is re-run on an atlas that already has them).
            preexisting = [c for c in sig_cols if c in atlas.columns]
            if preexisting:
                atlas = atlas.drop(columns=preexisting)
            # Merge via left-join on the symbol column without inflating rows.
            atlas = atlas.merge(
                sig_indexed,
                how="left",
                left_on=symbol_col,
                right_index=True,
            )

    # --- Multi-cohort spatial consensus (produced by 44_spatial_consensus.py) ---
    if "spatial_consensus" in spatial_results:
        cons = spatial_results["spatial_consensus"]
        if "gene" in cons.columns:
            cons_cols = [c for c in cons.columns if c.startswith("spatial_consensus_")]
            cons = cons.drop_duplicates(subset="gene", keep="first")
            cons_indexed = cons.set_index("gene")[cons_cols]
            # Idempotency: drop any pre-existing spatial_consensus_* cols from atlas.
            preexisting = [c for c in cons_cols if c in atlas.columns]
            if preexisting:
                atlas = atlas.drop(columns=preexisting)
            atlas = atlas.merge(
                cons_indexed,
                how="left",
                left_on=symbol_col,
                right_index=True,
            )

    return atlas


def enrichment_tests(atlas):
    """Test enrichment of spatial features in key gene sets (6+ tests)."""
    config = load_config()
    symbol_col = "human_symbol" if "human_symbol" in atlas.columns else ("symbol" if "symbol" in atlas.columns else atlas.columns[0])
    tests = []

    # Define SVG set
    svg_genes = set(atlas[atlas.get("spatial_is_svg", False) == True][symbol_col])
    if len(svg_genes) == 0:
        print("  No SVGs found, skipping enrichment tests")
        return pd.DataFrame()

    # Fisher background = genes actually TESTED for spatial autocorrelation, NOT
    # the full atlas (F071/F184). SVGs can only be called among the ~3,000 HVGs
    # squidpy tested per condition; `spatial_is_svg` is non-null exactly for
    # atlas genes present in that tested panel (True/False), and NaN for the
    # ~24k genes never eligible to be an SVG. Using the full atlas as the
    # background inflates the Fisher 'd' cell and every odds ratio. Restricting
    # to the tested universe gives the correctly-conditioned enrichment.
    if "spatial_is_svg" in atlas.columns:
        tested_mask = atlas["spatial_is_svg"].notna()
        all_genes = set(atlas.loc[tested_mask, symbol_col].dropna())
    else:
        all_genes = set(atlas[symbol_col].dropna())
    # SVGs must lie within the tested universe (defensive intersection).
    svg_genes = svg_genes & all_genes
    print(f"  SVG-tested background universe: {len(all_genes)} genes")

    # Build gene sets to test against SVGs
    gene_sets = {}

    # 1. Dream DEGs (intersect with the SVG-tested universe so the Fisher 2x2
    #    is conditioned on the same background — F071/F184)
    if "dream_padj" in atlas.columns:
        gene_sets["Dream_DEGs"] = set(atlas[atlas["dream_padj"] < 0.1][symbol_col]) & all_genes

    # 2. Conserved
    conserved = load_conserved()
    if conserved:
        gene_sets["Conserved"] = set(conserved) & all_genes

    # 3. Hepatocyte_intrinsic (from deconv attribution scores)
    try:
        deconv_scores = load_deconv_scores()
        hep_genes = deconv_scores[deconv_scores["category"] == "Hepatocyte_intrinsic"]
        ensembl_to_symbol = build_ensembl_to_symbol_map()
        hep_symbols = set(hep_genes["gene"].map(ensembl_to_symbol).dropna()) & all_genes
        if hep_symbols:
            gene_sets["Hepatocyte_intrinsic"] = hep_symbols
    except Exception as e:
        print(f"  WARNING: Could not load hep-intrinsic gene set: {e}")

    # 4. Sex-dimorphic DEGs (Male_biased + Female_biased + Divergent)
    try:
        sex_path = resolve_path(config["paths"]["sex_degs"])
        sex_df = pd.read_csv(sex_path)
        sex_classes = ["Male_biased", "Female_biased", "Divergent",
                       "Male_specific", "Female_specific", "Sex_divergent"]
        sex_genes_df = sex_df[sex_df["sex_class"].isin(sex_classes)]
        sex_ensembl = sex_genes_df["gene"].tolist()
        if not hasattr(enrichment_tests, "_e2s"):
            enrichment_tests._e2s = build_ensembl_to_symbol_map()
        sex_symbols = set(pd.Series(sex_ensembl).map(enrichment_tests._e2s).dropna()) & all_genes
        if sex_symbols:
            gene_sets["Sex_dimorphic_DEGs"] = sex_symbols
    except Exception as e:
        print(f"  WARNING: Could not load sex DEG gene set: {e}")

    # 5. Druggable genes (from multi-evidence atlas)
    if "dgidb_druggable" in atlas.columns:
        druggable = set(atlas[atlas["dgidb_druggable"] == True][symbol_col]) & all_genes
        if druggable:
            gene_sets["Druggable_genes"] = druggable

    # 6. Top multi-evidence (top 5% by active evidence layers)
    if "layers_active" in atlas.columns:
        threshold = atlas["layers_active"].quantile(0.95)
        top_me = set(atlas[atlas["layers_active"] >= threshold][symbol_col]) & all_genes
        if top_me:
            gene_sets["Top5pct_multi_evidence"] = top_me

    print(f"  Testing SVGs ({len(svg_genes)}) against {len(gene_sets)} gene sets")

    # Run Fisher's exact tests
    for name, gene_set in gene_sets.items():
        a = len(svg_genes & gene_set)
        b = len(svg_genes - gene_set)
        c = len(gene_set - svg_genes)
        d = len(all_genes - svg_genes - gene_set)
        if min(a + b, c + d, a + c, b + d) > 0:
            odds, pval = fisher_exact([[a, b], [c, d]])
            tests.append({
                "gene_set": name, "n_overlap": a,
                "n_svg": len(svg_genes), "n_set": len(gene_set),
                "odds_ratio": odds, "fisher_pval": pval,
            })

    tests_df = pd.DataFrame(tests)
    if len(tests_df) == 0:
        return tests_df

    # Also test zonation-classified genes if available
    zon_tests = _zonation_enrichment_tests(atlas, gene_sets, symbol_col, all_genes)
    if len(zon_tests) > 0:
        tests_df = pd.concat([tests_df, zon_tests], ignore_index=True)

    # Apply BH correction across all tests
    _, padj, _, _ = multipletests(tests_df["fisher_pval"].values, method="fdr_bh")
    tests_df["fisher_padj_bh"] = padj

    return tests_df


def _zonation_enrichment_tests(atlas, gene_sets, symbol_col, all_genes):
    """Test enrichment of zonation classes in gene sets."""
    if "spatial_zonation_class" not in atlas.columns:
        return pd.DataFrame()

    tests = []
    for zon_class in ["Periportal-enriched", "Pericentral-enriched"]:
        zon_genes = set(atlas[atlas["spatial_zonation_class"] == zon_class][symbol_col]) & all_genes
        if len(zon_genes) < 5:
            continue
        for name, gene_set in gene_sets.items():
            a = len(zon_genes & gene_set)
            b = len(zon_genes - gene_set)
            c = len(gene_set - zon_genes)
            d = len(all_genes - zon_genes - gene_set)
            if min(a + b, c + d, a + c, b + d) > 0:
                odds, pval = fisher_exact([[a, b], [c, d]])
                tests.append({
                    "gene_set": f"{zon_class}_in_{name}",
                    "n_overlap": a, "n_svg": len(zon_genes),
                    "n_set": len(gene_set),
                    "odds_ratio": odds, "fisher_pval": pval,
                })

    return pd.DataFrame(tests)


def main():
    print_header("06: Integration with Multi-Evidence Atlas")

    config = load_config()

    # Load existing atlas
    atlas = load_multi_evidence_atlas()
    n_cols_before = len(atlas.columns)
    n_rows_before = len(atlas)
    print(f"  Atlas: {n_rows_before} genes × {n_cols_before} columns")

    # Load spatial results
    spatial_results = load_spatial_results()
    print(f"  Spatial result sets loaded: {list(spatial_results.keys())}")

    # Add spatial evidence
    atlas = build_spatial_evidence(atlas, spatial_results)
    n_new = len(atlas.columns) - n_cols_before
    print(f"  Added {n_new} spatial columns")

    # Count genes with spatial annotations (existing + new govaere/consensus)
    new_family_prefixes = (
        "spatial_govaere2026_",
        "signature_govaere2026_",
        "spatial_consensus_",
    )
    for prefix in new_family_prefixes:
        cols = [c for c in atlas.columns if c.startswith(prefix)]
        for col in cols:
            n_annotated = atlas[col].notna().sum()
            print(f"    {col}: {n_annotated} genes annotated")

    # Save updated atlas
    output_path = RESULTS_DIR / "integration" / "multi_evidence_atlas_with_spatial.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    atlas.to_csv(output_path, index=False)
    print(f"\n  Saved: {output_path} ({len(atlas)} genes × {len(atlas.columns)} columns)")

    # ---------- Validation assertions ----------
    print("\n  Running validation assertions...")
    failures = []

    new_govaere_spatial = [c for c in atlas.columns if c.startswith("spatial_govaere2026_")]
    new_signature      = [c for c in atlas.columns if c.startswith("signature_govaere2026_")]
    new_consensus      = [c for c in atlas.columns if c.startswith("spatial_consensus_")]

    # Minimum 8 cols are the pre-subcluster contrasts (4 contrasts × {logfc, padj}).
    # MetMac/TransMac/KCpost from 41b add 6 more → 14 expected when 41b ran;
    # preMac + Monocyte are conditional on per-group n>=50 in 41b.
    if len(new_govaere_spatial) >= 8:
        print(f"    PASS: {len(new_govaere_spatial)} spatial_govaere2026_* columns (>=8 required)")
    else:
        msg = f"FAIL: only {len(new_govaere_spatial)} spatial_govaere2026_* columns (>=8 required)"
        failures.append(msg); print(f"    {msg}")

    # Hard requirement: the MetMac axis must be present once 41b has run.
    metmac_cols = [c for c in atlas.columns
                   if c.startswith("spatial_govaere2026_cosmx_metmac_mash_")]
    metmac_csv = RESULTS_DIR / "govaere2026" / "cosmx_de_MetMac_MASH_vs_noMASH.csv"
    if metmac_csv.exists():
        if len(metmac_cols) >= 2:
            print(f"    PASS: MetMac cols present ({metmac_cols})")
        else:
            msg = f"FAIL: MetMac DE CSV present but {len(metmac_cols)} MetMac cols (expected 2)"
            failures.append(msg); print(f"    {msg}")
    else:
        print(f"    NOTE: MetMac DE CSV not yet produced (41b not run) — skipping MetMac col check")

    if len(new_signature) >= 12:
        print(f"    PASS: {len(new_signature)} signature_govaere2026_* columns (>=12 required)")
    else:
        msg = f"FAIL: only {len(new_signature)} signature_govaere2026_* columns (>=12 required)"
        failures.append(msg); print(f"    {msg}")

    if len(new_consensus) >= 3:
        print(f"    PASS: {len(new_consensus)} spatial_consensus_* columns (>=3 required)")
    else:
        msg = f"FAIL: only {len(new_consensus)} spatial_consensus_* columns (>=3 required)"
        failures.append(msg); print(f"    {msg}")

    for col in ("spatial_morans_i", "spatial_is_svg"):
        if col in atlas.columns:
            print(f"    PASS: existing column {col} preserved")
        else:
            msg = f"FAIL: existing column {col} missing from output"
            failures.append(msg); print(f"    {msg}")

    if len(atlas) == n_rows_before:
        print(f"    PASS: row count unchanged ({n_rows_before})")
    else:
        msg = f"FAIL: row count changed {n_rows_before} -> {len(atlas)}"
        failures.append(msg); print(f"    {msg}")

    # GPNMB sanity: positive logFC in cosmx KC MASH
    if "spatial_govaere2026_cosmx_kc_mash_logfc" in atlas.columns:
        sym_col = "human_symbol" if "human_symbol" in atlas.columns else atlas.columns[0]
        gpnmb_row = atlas[atlas[sym_col] == "GPNMB"]
        if len(gpnmb_row) == 1:
            val = gpnmb_row["spatial_govaere2026_cosmx_kc_mash_logfc"].iloc[0]
            if pd.notna(val) and val > 0:
                print(f"    PASS: GPNMB cosmx_kc_mash_logfc = {val:.4f} (>0)")
            else:
                msg = f"FAIL: GPNMB cosmx_kc_mash_logfc = {val} (expected >0)"
                failures.append(msg); print(f"    {msg}")
        else:
            msg = f"FAIL: GPNMB row count = {len(gpnmb_row)} (expected 1)"
            failures.append(msg); print(f"    {msg}")

    # Spot checks: print values
    print("\n  Spot checks:")
    sym_col = "human_symbol" if "human_symbol" in atlas.columns else atlas.columns[0]
    spot_cols = [c for c in atlas.columns if (c.startswith("spatial_govaere2026_")
                                              or c.startswith("signature_govaere2026_")
                                              or c.startswith("spatial_consensus_"))]
    for g in ("GPNMB", "LPL", "FABP5"):
        row = atlas[atlas[sym_col] == g]
        if row.empty:
            print(f"    {g}: NOT in atlas")
            continue
        r = row.iloc[0]
        non_null = {c: r[c] for c in spot_cols if pd.notna(r[c])}
        print(f"    {g}: {len(non_null)} non-null new columns")
        for c, v in list(non_null.items())[:8]:
            if isinstance(v, float):
                print(f"      {c} = {v:.4g}")
            else:
                print(f"      {c} = {v}")

    if failures:
        print(f"\n  {len(failures)} validation FAILURES:")
        for f in failures:
            print(f"    - {f}")
        sys.exit(1)
    else:
        print("\n  All validation assertions PASSED")

    # Enrichment tests
    print("\n  Running enrichment tests...")
    enrichment = enrichment_tests(atlas)
    if len(enrichment) > 0:
        save_csv(enrichment, "spatial_enrichment_tests.csv", subdir="integration")
        for _, row in enrichment.iterrows():
            sig = "***" if row["fisher_pval"] < 0.001 else "**" if row["fisher_pval"] < 0.01 else "*" if row["fisher_pval"] < 0.05 else "ns"
            padj_sig = f" (BH: {row['fisher_padj_bh']:.2e})" if "fisher_padj_bh" in row.index else ""
            print(f"    {row['gene_set']}: OR={row['odds_ratio']:.2f}, "
                  f"p={row['fisher_pval']:.2e}{padj_sig} {sig} "
                  f"({row['n_overlap']}/{row['n_svg']})")

    print_header("06: Complete")


if __name__ == "__main__":
    main()
