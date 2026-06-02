"""Build machine-readable column dictionary for multi_evidence_atlas.csv.

Output: RNA-seq/results/multi_evidence/atlas_columns.tsv

Curated annotation (column -> source/method/LD panel/deprecation/etc.) is
combined with auto-computed schema (dtype, missingness, value range, n_unique).

RP7 Top-1 deliverable from 2026-05-12 deliberation; resolves the five+
`coloc_*_pp4` variants by tagging each with its method (ABF/SuSiE) x LD panel
(PolyFun / sghatan-v1 / Broadaway / sc-eQTL / per-GWAS / legacy-ambiguous),
flags the eight retired MR columns, and anchors the join key on GENCODE v49.

Run under rnaseq env:
    micromamba activate rnaseq
    python scripts/build_atlas_column_dictionary.py
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
ATLAS_PATH = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
OUT_PATH = PROJECT_ROOT / "RNA-seq/results/multi_evidence/atlas_columns.tsv"


# Hand-curated annotation. Schema:
#   source         : S1 bulk human dream / S2 mouse meta / S3 genetic causal /
#                    S4 essentiality / S5 epigenomic / S6 spatial /
#                    S7 single-cell / proteomics / drug / progression /
#                    sex / pathway / annotation / ortholog / identifier
#   method         : free-text method tag (e.g., dream limma-voom, coloc.abf(),
#                    coloc.susie(), SuSiEx, MESuSiE, MR-IVW, OTTERS S-PrediXcan,
#                    cTWAS, INTACT, HyPrColoc, MuSiC deconvolution, etc.).
#                    Empty string = no specific method tag.
#   ld_panel       : PolyFun-EUR / sghatan-v1 / Broadaway-EUR /
#                    1KG-EAS / TopLD / PanUKBB-AFR / PanUKBB-CSA / NA
#   gwas_scope     : best-of-portfolio / per-GWAS:<gwas> / NA
#   units          : free-text physical/statistical units
#   deprecated     : True for retired columns (e.g., MR retired 2026-04-22)
#   provenance     : script that writes the column (best effort).
#   notes          : ambiguity flags, sensitivity caveats, paper references.
#
# Every numeric column is assumed to live in real-valued domain unless the
# `units` field implies otherwise (probability [0,1], count >=0, etc.).
COL_ANNOT: dict[str, dict[str, str | bool]] = {
    # -------- identifiers (1-3) --------
    "human_symbol": dict(
        source="identifier", method="HGNC", ld_panel="NA", gwas_scope="NA",
        units="HGNC gene symbol", deprecated=False,
        provenance="27a_assemble_evidence_atlas.R",
        notes="Primary human gene symbol; not version-pinned. Use ensembl_id for joins.",
    ),
    "ensembl_id": dict(
        source="identifier", method="GENCODE v49 / GRCh38.p14",
        ld_panel="NA", gwas_scope="NA",
        units="Ensembl gene ID (versioned)", deprecated=False,
        provenance="27a_assemble_evidence_atlas.R",
        notes="JOIN KEY. GENCODE v49 anchor. Pre-Phase-10 atlas snapshots may carry unversioned IDs.",
    ),
    "gene_biotype": dict(
        source="annotation", method="GENCODE v49 biotype",
        ld_panel="NA", gwas_scope="NA",
        units="categorical (protein_coding, lncRNA, pseudogene, ...)",
        deprecated=False, provenance="build_gencode_metadata.R",
        notes="Filter atlas to specific biotype via this column.",
    ),
    # -------- ortholog (4, 136) --------
    "mouse_ortholog": dict(
        source="ortholog", method="Ensembl biomaRt 1:1 ortholog",
        ld_panel="NA", gwas_scope="NA",
        units="MGI gene symbol", deprecated=False,
        provenance="streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz",
        notes="Empty when no 1:1 mouse-human ortholog. 1:many collapsed silently.",
    ),
    "mouse_ensembl": dict(
        source="ortholog", method="Ensembl biomaRt 1:1 ortholog",
        ld_panel="NA", gwas_scope="NA",
        units="Ensembl mouse gene ID (GENCODE vM38)",
        deprecated=False,
        provenance="streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz",
        notes="Empty when no 1:1 ortholog.",
    ),
    # -------- S1 bulk human dream DEG (5-7, 86-87) --------
    "dream_logFC": dict(
        source="S1", method="dream limma-voom mega-analysis",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change (Disease vs Control)", deprecated=False,
        provenance="05_dream_mega_analysis.R / 05b output dream_results_ashr.csv",
        notes="5-cohort canonical mega (n=847). Tier 1 = padj<0.05 & |logFC|>0.5.",
    ),
    "dream_padj": dict(
        source="S1", method="dream limma-voom (BH-adjusted)",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="05_dream_mega_analysis.R",
        notes="Two-tier system: padj<0.05 (Tier 1); padj<0.1 (Tier 2/exploratory).",
    ),
    "dream_tstat": dict(
        source="S1", method="dream limma-voom",
        ld_panel="NA", gwas_scope="NA",
        units="t-statistic", deprecated=False,
        provenance="05_dream_mega_analysis.R",
        notes="Signed test statistic; direction matches dream_logFC.",
    ),
    "human_consensus_tier": dict(
        source="S1", method="multi-threshold tier assignment",
        ld_panel="NA", gwas_scope="NA",
        units="categorical (Tier1/Tier2/Tier3/NotDE)", deprecated=False,
        provenance="27a_assemble_evidence_atlas.R",
        notes="Tier 1: padj<0.05 & |logFC|>0.5; Tier 2: padj<0.05 no LFC filter; Tier 3: padj<0.1.",
    ),
    "dream_logFC_M": dict(
        source="sex", method="dream sex-stratified (male only)",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change (M only)", deprecated=False,
        provenance="26_sex_stratified_analysis.R",
        notes="Stratified-DE artifact-prone; cross-check sex_interaction_padj.",
    ),
    "dream_logFC_F": dict(
        source="sex", method="dream sex-stratified (female only)",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change (F only)", deprecated=False,
        provenance="26_sex_stratified_analysis.R",
        notes="Stratified-DE artifact-prone; cross-check sex_interaction_padj.",
    ),
    # -------- S2 mouse meta-analysis (8-17) --------
    "mouse_meta_logFC": dict(
        source="S2", method="dream pooled / metafor across 5 diet models",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change (mouse Disease vs Control)",
        deprecated=False,
        provenance="M03-M04 mouse integration",
        notes="5 mouse diet models (CDA, CDHFD, FPC, NASH, WD-CCl4). 463 samples.",
    ),
    "mouse_meta_padj": dict(
        source="S2", method="dream pooled BH-adjusted",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="M03-M04",
        notes="",
    ),
    "n_diets_sig": dict(
        source="S2", method="diet-wise DE count",
        ld_panel="NA", gwas_scope="NA",
        units="integer count (0-5)", deprecated=False,
        provenance="M03-M04",
        notes="Number of mouse diet models where gene is DE (padj<0.05).",
    ),
    "mouse_consensus_tier": dict(
        source="S2", method="diet-wise consensus tier",
        ld_panel="NA", gwas_scope="NA",
        units="categorical", deprecated=False,
        provenance="M03-M04",
        notes="",
    ),
    "primary_category": dict(
        source="annotation", method="cross-species concordance classification",
        ld_panel="NA", gwas_scope="NA",
        units="categorical (Conserved / Human_only / Mouse_only / Discordant / ...)",
        deprecated=False,
        provenance="Analysis/Cross_Species_Concordance/",
        notes="Conserved_Core = primary_category==Conserved & n_concordant_diets>=3.",
    ),
    "translatability_score": dict(
        source="annotation", method="cross-species composite score",
        ld_panel="NA", gwas_scope="NA",
        units="continuous (higher = more translatable)",
        deprecated=False,
        provenance="Analysis/Cross_Species_Concordance/",
        notes="",
    ),
    "n_concordant_diets": dict(
        source="annotation", method="diet-wise direction concordance",
        ld_panel="NA", gwas_scope="NA",
        units="integer count (0-5)", deprecated=False,
        provenance="Analysis/Cross_Species_Concordance/",
        notes="",
    ),
    "is_conserved": dict(
        source="annotation", method="binary conservation flag",
        ld_panel="NA", gwas_scope="NA",
        units="boolean", deprecated=False,
        provenance="Analysis/Cross_Species_Concordance/",
        notes="TRUE for Conserved_Core membership.",
    ),
    "best_mouse_model": dict(
        source="S2", method="strongest-LFC diet model per gene",
        ld_panel="NA", gwas_scope="NA",
        units="categorical (CDA/CDHFD/FPC/NASH/WD-CCl4)",
        deprecated=False, provenance="M03-M04",
        notes="",
    ),
    # -------- S3 TWAS (18-19, 180-184) --------
    "twas_z": dict(
        source="S3", method="S-PrediXcan TWAS Z-score (GTEx Liver)",
        ld_panel="GTEx-Liver", gwas_scope="best-of-portfolio",
        units="Z-statistic", deprecated=False,
        provenance="19 (archived) / OTTERS replacement see otters_*",
        notes="GTEx v8 Liver N=208 (underpowered for formal replication).",
    ),
    "twas_pval": dict(
        source="S3", method="S-PrediXcan TWAS p-value",
        ld_panel="GTEx-Liver", gwas_scope="best-of-portfolio",
        units="p-value", deprecated=False,
        provenance="19 (archived)",
        notes="",
    ),
    "otters_broadaway_pval": dict(
        source="S3", method="OTTERS TWAS (Broadaway eQTL)",
        ld_panel="Broadaway-EUR", gwas_scope="best-of-portfolio",
        units="p-value", deprecated=False,
        provenance="63_otters_twas.sh",
        notes="Broadaway liver eQTL N=1,183 (EUR-only, hg19).",
    ),
    "otters_broadaway_z": dict(
        source="S3", method="OTTERS TWAS Z-statistic",
        ld_panel="Broadaway-EUR", gwas_scope="best-of-portfolio",
        units="Z-statistic", deprecated=False,
        provenance="63_otters_twas.sh",
        notes="",
    ),
    "otters_n_gwas_sig": dict(
        source="S3", method="OTTERS TWAS portfolio count",
        ld_panel="Broadaway-EUR", gwas_scope="all-tested",
        units="integer count", deprecated=False,
        provenance="63_otters_twas.sh",
        notes="Number of GWAS in portfolio with significant OTTERS TWAS Z.",
    ),
    "ctwas_pip": dict(
        source="S3", method="cTWAS posterior inclusion probability",
        ld_panel="Broadaway-EUR", gwas_scope="best-of-portfolio",
        units="PIP [0,1]", deprecated=False,
        provenance="Phase L causal overhaul",
        notes="cTWAS = causal TWAS (Zhao et al. 2024).",
    ),
    "ctwas_pval": dict(
        source="S3", method="cTWAS p-value",
        ld_panel="Broadaway-EUR", gwas_scope="best-of-portfolio",
        units="p-value", deprecated=False,
        provenance="Phase L",
        notes="",
    ),
    # -------- S3 COLOC: LEGACY AMBIGUOUS (20) --------
    "coloc_pp4": dict(
        source="S3", method="LEGACY COLOC (ambiguous)",
        ld_panel="legacy", gwas_scope="best-of-portfolio (legacy)",
        units="PP.H4 posterior probability [0,1]",
        deprecated=False,
        provenance="pre-Phase-10 atlas snapshot",
        notes="LEGACY COLUMN -- predates Phase 10 PolyFun production swap (2026-05-06). Method (ABF vs SuSiE) and LD panel ambiguous. Canonical replacements: coloc_best_pp4_polyfun (ABF) and coloc_best_susie_pp4_polyfun (SuSiE). Retain for reproducibility of pre-2026-05-06 figures only.",
    ),
    # -------- S3 COLOC ABF method (21-27, 112-113) --------
    "coloc_abf_best_pp4": dict(
        source="S3", method="coloc.abf() Wakefield single-causal",
        ld_panel="multi-panel (best across portfolio)",
        gwas_scope="best-of-portfolio",
        units="PP.H4 posterior [0,1]",
        deprecated=False, provenance="35_coloc_pipeline.R",
        notes="ABF assumes single causal variant per region.",
    ),
    "coloc_abf_best_gwas": dict(
        source="S3", method="ABF best GWAS identifier",
        ld_panel="multi-panel", gwas_scope="best-of-portfolio",
        units="GWAS name", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="Which GWAS achieves the highest ABF PP4 per gene.",
    ),
    "coloc_abf_n_gwas_h4_05": dict(
        source="S3", method="ABF portfolio count at PP4>0.5",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "coloc_abf_n_gwas_h4_08": dict(
        source="S3", method="ABF portfolio count at PP4>0.8",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "coloc_n_gwas_tested": dict(
        source="S3", method="COLOC portfolio test count",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="Number of GWAS tested against this eGene.",
    ),
    "coloc_n_groups_h4_05": dict(
        source="S3", method="ABF group-collapsed PP4>0.5 count",
        ld_panel="multi-panel", gwas_scope="grouped",
        units="integer count", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="Collapsing related GWAS (e.g., overlapping cohorts) into groups.",
    ),
    "coloc_n_groups_h4_08": dict(
        source="S3", method="ABF group-collapsed PP4>0.8 count",
        ld_panel="multi-panel", gwas_scope="grouped",
        units="integer count", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "coloc_n_gwas_h4_08": dict(
        source="S3", method="ABF portfolio count at PP4>0.8 (duplicate)",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="27a_assemble_evidence_atlas.R",
        notes="DUPLICATE of coloc_abf_n_gwas_h4_08; kept for backward compatibility.",
    ),
    # -------- S3 COLOC SuSiE method (28-34, 185-187) --------
    "coloc_susie_n_gwas_h4_05": dict(
        source="S3", method="coloc.susie() portfolio count at PP4>0.5",
        ld_panel="legacy/sghatan-v1 era",
        gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="35s_susie_coloc_broadaway.R",
        notes="Pre-Phase-10 LD panel; for PolyFun-canonical use coloc_n_gwas_susie_h4_05_polyfun.",
    ),
    "coloc_susie_n_gwas_h4_08": dict(
        source="S3", method="coloc.susie() portfolio count at PP4>0.8",
        ld_panel="legacy/sghatan-v1 era",
        gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="35s_susie_coloc_broadaway.R",
        notes="",
    ),
    "coloc_susie_n_gwas_h4_09": dict(
        source="S3", method="coloc.susie() portfolio count at PP4>0.9",
        ld_panel="legacy/sghatan-v1 era",
        gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="35s_susie_coloc_broadaway.R",
        notes="",
    ),
    "coloc_susie_n_pairs_total": dict(
        source="S3", method="coloc.susie() total eGene-GWAS pair count",
        ld_panel="legacy/sghatan-v1 era",
        gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="35s_susie_coloc_broadaway.R",
        notes="",
    ),
    "coloc_susie_success_rate": dict(
        source="S3", method="SuSiE convergence success rate",
        ld_panel="legacy/sghatan-v1 era",
        gwas_scope="all-portfolio",
        units="fraction [0,1]", deprecated=False,
        provenance="35s_susie_coloc_broadaway.R",
        notes="Atlas-wide SuSiE convergence ceiling = 18.3% across 18,975 eGenes.",
    ),
    "coloc_n_groups_susie_h4_05": dict(
        source="S3", method="SuSiE group-collapsed PP4>0.5 count",
        ld_panel="legacy/sghatan-v1 era",
        gwas_scope="grouped",
        units="integer count", deprecated=False,
        provenance="35s_susie_coloc_broadaway.R",
        notes="",
    ),
    "coloc_n_groups_susie_h4_08": dict(
        source="S3", method="SuSiE group-collapsed PP4>0.8 count",
        ld_panel="legacy/sghatan-v1 era",
        gwas_scope="grouped",
        units="integer count", deprecated=False,
        provenance="35s_susie_coloc_broadaway.R",
        notes="",
    ),
    "coloc_susie_best_pp4": dict(
        source="S3", method="coloc.susie() best PP4",
        ld_panel="legacy/sghatan-v1 era",
        gwas_scope="best-of-portfolio",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35s_susie_coloc_broadaway.R",
        notes="LEGACY for pre-Phase-10. Canonical = coloc_best_susie_pp4_polyfun.",
    ),
    "coloc_susie_n_signals": dict(
        source="S3", method="SuSiE credible-set count per eGene",
        ld_panel="legacy/sghatan-v1 era",
        gwas_scope="best-of-portfolio",
        units="integer count", deprecated=False,
        provenance="35s_susie_coloc_broadaway.R",
        notes="Multi-causal indicator.",
    ),
    "coloc_susie_best_gwas": dict(
        source="S3", method="SuSiE best GWAS identifier",
        ld_panel="legacy/sghatan-v1 era",
        gwas_scope="best-of-portfolio",
        units="GWAS name", deprecated=False,
        provenance="35s_susie_coloc_broadaway.R",
        notes="",
    ),
    # -------- S3 COLOC PolyFun-canonical (112-115, 118-124) --------
    "coloc_best_gwas_polyfun": dict(
        source="S3", method="ABF best GWAS (PolyFun-EUR LD)",
        ld_panel="PolyFun-EUR",
        gwas_scope="best-of-portfolio (17 EUR GWAS)",
        units="GWAS name", deprecated=False,
        provenance="77b_add_polyfun_atlas_columns.R",
        notes="PolyFun EUR LD reference (Weissbrod 2020, ~337K UKBB EUR).",
    ),
    "coloc_best_pp4_polyfun": dict(
        source="S3", method="coloc.abf() with PolyFun-EUR LD",
        ld_panel="PolyFun-EUR",
        gwas_scope="best-of-portfolio (17 EUR GWAS)",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="77b_add_polyfun_atlas_columns.R",
        notes="CANONICAL ABF column post Phase 10 swap (2026-05-06).",
    ),
    "coloc_best_susie_gwas_polyfun": dict(
        source="S3", method="SuSiE best GWAS (PolyFun-EUR LD)",
        ld_panel="PolyFun-EUR",
        gwas_scope="best-of-portfolio (17 EUR GWAS)",
        units="GWAS name", deprecated=False,
        provenance="77b_add_polyfun_atlas_columns.R",
        notes="",
    ),
    "coloc_best_susie_pp4_polyfun": dict(
        source="S3", method="coloc.susie() with PolyFun-EUR LD",
        ld_panel="PolyFun-EUR",
        gwas_scope="best-of-portfolio (17 EUR GWAS)",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="77b_add_polyfun_atlas_columns.R",
        notes="CANONICAL SuSiE column post Phase 10 swap (2026-05-06). 538 genes >0.5; 322 >0.8; 226 >0.9.",
    ),
    "coloc_is_mhc": dict(
        source="S3", method="MHC region flag",
        ld_panel="NA", gwas_scope="NA",
        units="boolean", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="True for genes in MHC region (chr6:25-35Mb); fine-mapping unreliable here.",
    ),
    "coloc_ld_cluster_flag": dict(
        source="S3", method="LD cluster annotation",
        ld_panel="NA", gwas_scope="NA",
        units="categorical/string", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="LD-cluster membership for multi-gene credible sets.",
    ),
    "coloc_n_gwas_h4_05_polyfun": dict(
        source="S3", method="ABF portfolio count PP4>0.5 (PolyFun)",
        ld_panel="PolyFun-EUR", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="77b_add_polyfun_atlas_columns.R",
        notes="",
    ),
    "coloc_n_gwas_h4_08_polyfun": dict(
        source="S3", method="ABF portfolio count PP4>0.8 (PolyFun)",
        ld_panel="PolyFun-EUR", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="77b_add_polyfun_atlas_columns.R",
        notes="",
    ),
    "coloc_n_gwas_susie_h4_05_polyfun": dict(
        source="S3", method="SuSiE portfolio count PP4>0.5 (PolyFun)",
        ld_panel="PolyFun-EUR", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="77b_add_polyfun_atlas_columns.R",
        notes="",
    ),
    "coloc_n_gwas_susie_h4_08_polyfun": dict(
        source="S3", method="SuSiE portfolio count PP4>0.8 (PolyFun)",
        ld_panel="PolyFun-EUR", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="77b_add_polyfun_atlas_columns.R",
        notes="",
    ),
    "coloc_n_gwas_susie_h4_09_polyfun": dict(
        source="S3", method="SuSiE portfolio count PP4>0.9 (PolyFun)",
        ld_panel="PolyFun-EUR", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="77b_add_polyfun_atlas_columns.R",
        notes="",
    ),
    "coloc_susie_success_rate_polyfun": dict(
        source="S3", method="SuSiE convergence success rate (PolyFun)",
        ld_panel="PolyFun-EUR", gwas_scope="all-portfolio",
        units="fraction [0,1]", deprecated=False,
        provenance="77b_add_polyfun_atlas_columns.R",
        notes="",
    ),
    # -------- S3 COLOC sghatan-v1 archived (215-216) --------
    "coloc_susie_best_pp4_v1_sghatan": dict(
        source="S3", method="coloc.susie() with sghatan-v1 UKBB EUR LD",
        ld_panel="sghatan-v1 (archived)",
        gwas_scope="best-of-portfolio (17 EUR GWAS)",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="archive/ld_panel_v1_sghatan_2026-05-06/",
        notes="ARCHIVED frozen sensitivity reference (Phase 10 swap 2026-05-06). PolyFun-vs-v1 paired r=0.986, rho=0.970, Jaccard 0.456 at PP4>0.5. v1 retention 83.3%. 324 genes >0.5 in v1.",
    ),
    "coloc_susie_n_signals_v1_sghatan": dict(
        source="S3", method="SuSiE credible-set count (sghatan-v1)",
        ld_panel="sghatan-v1 (archived)",
        gwas_scope="best-of-portfolio",
        units="integer count", deprecated=False,
        provenance="archive/ld_panel_v1_sghatan_2026-05-06/",
        notes="ARCHIVED sensitivity reference.",
    ),
    # -------- S3 sc-eQTL COLOC (35-38) --------
    "sceqtl_coloc_pp4_hep": dict(
        source="S3", method="sc-eQTL COLOC (hepatocyte)",
        ld_panel="see-source", gwas_scope="best-of-portfolio",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="sc-TWAS pipeline",
        notes="Hepatocyte-specific sc-eQTL COLOC.",
    ),
    "sceqtl_coloc_best_pp4": dict(
        source="S3", method="sc-eQTL COLOC (best cell type)",
        ld_panel="see-source", gwas_scope="best-of-portfolio",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="sc-TWAS pipeline",
        notes="Best PP4 across cell types.",
    ),
    "sceqtl_coloc_cell_type": dict(
        source="S3", method="sc-eQTL COLOC (best cell type identifier)",
        ld_panel="see-source", gwas_scope="best-of-portfolio",
        units="cell type name", deprecated=False,
        provenance="sc-TWAS pipeline",
        notes="",
    ),
    "sceqtl_n_cell_types": dict(
        source="S3", method="sc-eQTL COLOC cell-type count",
        ld_panel="see-source", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="sc-TWAS pipeline",
        notes="",
    ),
    # -------- S3 per-GWAS COLOC PP4 (39-57) --------
    "broadaway_coloc_pp4": dict(
        source="S3", method="COLOC with Broadaway eQTL",
        ld_panel="Broadaway-EUR", gwas_scope="best-of-portfolio",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35s_susie_coloc_broadaway.R",
        notes="Broadaway liver eQTL N=1,183 EUR (hg19).",
    ),
    "ast_coloc_pp4": dict(
        source="S3", method="COLOC vs AST GWAS",
        ld_panel="see-source", gwas_scope="per-GWAS:AST",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "ggt_coloc_pp4": dict(
        source="S3", method="COLOC vs GGT GWAS",
        ld_panel="see-source", gwas_scope="per-GWAS:GGT",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "pdff_coloc_pp4": dict(
        source="S3", method="COLOC vs PDFF GWAS",
        ld_panel="see-source", gwas_scope="per-GWAS:PDFF",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="Sample overlap: 3 PDFF GWAS = same UKBB MRI cohort.",
    ),
    "ukbb_alt_coloc_pp4": dict(
        source="S3", method="COLOC vs UKBB ALT GWAS",
        ld_panel="see-source", gwas_scope="per-GWAS:UKBB_ALT",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "ukbb_alt_coloc_cell_type": dict(
        source="S3", method="UKBB ALT cell-type assignment",
        ld_panel="see-source", gwas_scope="per-GWAS:UKBB_ALT",
        units="cell type name", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "n_liver_enzyme_coloc": dict(
        source="S3", method="liver-enzyme GWAS COLOC count",
        ld_panel="multi-panel",
        gwas_scope="liver-enzyme-subset (ALT/AST/GGT)",
        units="integer count", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "best_liver_enzyme_pp4": dict(
        source="S3", method="best liver-enzyme COLOC PP4",
        ld_panel="multi-panel",
        gwas_scope="liver-enzyme-subset",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "finngen_nafld_coloc_pp4": dict(
        source="S3", method="COLOC vs FinnGen R12 NAFLD",
        ld_panel="see-source", gwas_scope="per-GWAS:FinnGen_NAFLD",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="run_finngen_coloc.sh",
        notes="",
    ),
    "finngen_nash_coloc_pp4": dict(
        source="S3", method="COLOC vs FinnGen R12 NASH",
        ld_panel="see-source", gwas_scope="per-GWAS:FinnGen_NASH",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="run_finngen_coloc.sh",
        notes="",
    ),
    "finngen_hcc_coloc_pp4": dict(
        source="S3", method="COLOC vs FinnGen R12 HCC",
        ld_panel="see-source", gwas_scope="per-GWAS:FinnGen_HCC",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="run_finngen_coloc.sh",
        notes="",
    ),
    "bbj_alt_coloc_pp4": dict(
        source="S3", method="COLOC vs BBJ ALT",
        ld_panel="1KG-EAS (N=504)", gwas_scope="per-GWAS:BBJ_ALT",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="run_bbj_coloc.sh",
        notes="BBJ EAS GWAS x EUR eQTL (known LD mismatch). ABF fallback retained.",
    ),
    "bbj_ast_coloc_pp4": dict(
        source="S3", method="COLOC vs BBJ AST",
        ld_panel="1KG-EAS (N=504)", gwas_scope="per-GWAS:BBJ_AST",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="run_bbj_coloc.sh",
        notes="",
    ),
    "bbj_ggt_coloc_pp4": dict(
        source="S3", method="COLOC vs BBJ GGT",
        ld_panel="1KG-EAS (N=504)", gwas_scope="per-GWAS:BBJ_GGT",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="run_bbj_coloc.sh",
        notes="THRB anchor: PP4=1.000 (resmetirom-validated target).",
    ),
    "ghouse_cirrhosis_coloc_pp4": dict(
        source="S3", method="COLOC vs Ghouse cirrhosis",
        ld_panel="see-source",
        gwas_scope="per-GWAS:Ghouse_cirrhosis",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "ghouse_hcc_coloc_pp4": dict(
        source="S3", method="COLOC vs Ghouse HCC",
        ld_panel="see-source", gwas_scope="per-GWAS:Ghouse_HCC",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "decode_nafl_coloc_pp4": dict(
        source="S3", method="COLOC vs deCODE NAFL",
        ld_panel="see-source", gwas_scope="per-GWAS:deCODE_NAFL",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "decode_cirrhosis_coloc_pp4": dict(
        source="S3", method="COLOC vs deCODE cirrhosis",
        ld_panel="see-source",
        gwas_scope="per-GWAS:deCODE_cirrhosis",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "decode_hcc_coloc_pp4": dict(
        source="S3", method="COLOC vs deCODE HCC",
        ld_panel="see-source", gwas_scope="per-GWAS:deCODE_HCC",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    # -------- S3 Pan-UKBB ancestry-specific (159-164) --------
    "panukbb_afr_alt_coloc_pp4": dict(
        source="S3", method="ABF fallback (no AFR liver eQTL)",
        ld_panel="PanUKBB-AFR", gwas_scope="per-GWAS:PanUKBB_AFR_ALT",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="AFR Pan-UKBB GWAS x EUR Broadaway eQTL. ABF-only. Methodologically NOT cross-ancestry validation; reframe as directional sensitivity.",
    ),
    "panukbb_afr_ast_coloc_pp4": dict(
        source="S3", method="ABF fallback (no AFR liver eQTL)",
        ld_panel="PanUKBB-AFR", gwas_scope="per-GWAS:PanUKBB_AFR_AST",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "panukbb_afr_ggt_coloc_pp4": dict(
        source="S3", method="ABF fallback (no AFR liver eQTL)",
        ld_panel="PanUKBB-AFR", gwas_scope="per-GWAS:PanUKBB_AFR_GGT",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "panukbb_csa_alt_coloc_pp4": dict(
        source="S3", method="ABF fallback (no SAS liver eQTL)",
        ld_panel="PanUKBB-CSA", gwas_scope="per-GWAS:PanUKBB_CSA_ALT",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="South Asian Pan-UKBB x EUR Broadaway. ABF-only.",
    ),
    "panukbb_csa_ast_coloc_pp4": dict(
        source="S3", method="ABF fallback (no SAS liver eQTL)",
        ld_panel="PanUKBB-CSA", gwas_scope="per-GWAS:PanUKBB_CSA_AST",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    "panukbb_csa_ggt_coloc_pp4": dict(
        source="S3", method="ABF fallback (no SAS liver eQTL)",
        ld_panel="PanUKBB-CSA", gwas_scope="per-GWAS:PanUKBB_CSA_GGT",
        units="PP.H4 posterior [0,1]", deprecated=False,
        provenance="35_coloc_pipeline.R",
        notes="",
    ),
    # -------- S3 SuSiEx + MESuSiE multi-ancestry (58-69) --------
    "susiex_max_pip": dict(
        source="S3", method="SuSiEx joint cross-ancestry fine-mapping",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="PIP [0,1]", deprecated=False,
        provenance="10_run_susiex.py",
        notes="SuSiEx = ancestry-aware joint SuSiE (Yuan 2024).",
    ),
    "susiex_n_loci": dict(
        source="S3", method="SuSiEx locus count",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="integer count", deprecated=False,
        provenance="10_run_susiex.py",
        notes="",
    ),
    "susiex_n_cs": dict(
        source="S3", method="SuSiEx credible-set count",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="integer count", deprecated=False,
        provenance="10_run_susiex.py",
        notes="",
    ),
    "susiex_cs_size_joint": dict(
        source="S3", method="SuSiEx joint credible-set size",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="integer count", deprecated=False,
        provenance="10_run_susiex.py",
        notes="",
    ),
    "susiex_trait_pairs": dict(
        source="S3", method="SuSiEx trait-pair identifiers",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="trait-pair string", deprecated=False,
        provenance="10_run_susiex.py",
        notes="",
    ),
    "mesusie_max_pip": dict(
        source="S3", method="MESuSiE max PIP across credible sets",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="PIP [0,1]", deprecated=False,
        provenance="Phase L mesusie",
        notes="MESuSiE = Multi-Ethnic SuSiE (Cai 2024).",
    ),
    "mesusie_max_pip_shared": dict(
        source="S3", method="MESuSiE shared-causal max PIP",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="PIP [0,1]", deprecated=False,
        provenance="Phase L mesusie",
        notes="Shared-causal variant probability across ancestries.",
    ),
    "mesusie_in_shared_cs": dict(
        source="S3", method="MESuSiE shared-credible-set membership",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="boolean", deprecated=False,
        provenance="Phase L mesusie",
        notes="",
    ),
    "mesusie_in_eur_cs": dict(
        source="S3", method="MESuSiE EUR-credible-set membership",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="boolean", deprecated=False,
        provenance="Phase L mesusie",
        notes="",
    ),
    "mesusie_in_eas_cs": dict(
        source="S3", method="MESuSiE EAS-credible-set membership",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="boolean", deprecated=False,
        provenance="Phase L mesusie",
        notes="",
    ),
    "mesusie_n_loci": dict(
        source="S3", method="MESuSiE locus count",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="integer count", deprecated=False,
        provenance="Phase L mesusie",
        notes="",
    ),
    "mesusie_trait_pairs": dict(
        source="S3", method="MESuSiE trait-pair identifiers",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="trait-pair string", deprecated=False,
        provenance="Phase L mesusie",
        notes="",
    ),
    # -------- S3 ancestry/multi-coloc summaries (70-72, 145-147) --------
    "n_coloc_sources": dict(
        source="S3", method="COLOC source count",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="27a_assemble_evidence_atlas.R",
        notes="Number of distinct COLOC data sources supporting this gene.",
    ),
    "n_ancestry_gwas": dict(
        source="S3", method="ancestry-stratified GWAS count",
        ld_panel="NA", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="48_cross_ancestry_replication.R",
        notes="",
    ),
    "cross_ancestry_coloc_replication": dict(
        source="S3", method="cross-ancestry COLOC replication flag",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="boolean/categorical", deprecated=False,
        provenance="48_cross_ancestry_replication.R",
        notes="67 EUR+EAS cross-ancestry COLOC genes.",
    ),
    "n_ancestry_gwas_dedup": dict(
        source="S3", method="ancestry GWAS count (deduplicated)",
        ld_panel="NA", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="48_cross_ancestry_replication.R",
        notes="Deduplicates overlapping cohorts (e.g., 3 UKBB-PDFF same cohort).",
    ),
    "n_validated_ancestry_gwas": dict(
        source="S3", method="ancestry-validated GWAS count",
        ld_panel="NA", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="48_cross_ancestry_replication.R",
        notes="",
    ),
    "n_validated_ancestry_gwas_dedup": dict(
        source="S3", method="ancestry-validated GWAS (deduplicated)",
        ld_panel="NA", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="48_cross_ancestry_replication.R",
        notes="",
    ),
    "cross_ancestry_coloc_replication_dedup": dict(
        source="S3", method="cross-ancestry replication (dedup)",
        ld_panel="multi-ancestry", gwas_scope="cross-ancestry pairs",
        units="boolean/categorical", deprecated=False,
        provenance="48_cross_ancestry_replication.R",
        notes="",
    ),
    # -------- S3 progression COLOC / Zenodo (73-76) --------
    "has_progression_coloc": dict(
        source="S3", method="progression-stratified COLOC flag",
        ld_panel="multi-panel", gwas_scope="progression-subset",
        units="boolean", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="See progression_coloc_class for granularity.",
    ),
    "zenodo_nafld_coloc": dict(
        source="S3", method="Zenodo external NAFLD COLOC overlay",
        ld_panel="external", gwas_scope="external",
        units="boolean/PP4", deprecated=False,
        provenance="external Zenodo import",
        notes="",
    ),
    "zenodo_coloc_cell_types": dict(
        source="S3", method="Zenodo cell-type COLOC overlay",
        ld_panel="external", gwas_scope="external",
        units="cell-type string", deprecated=False,
        provenance="external Zenodo import",
        notes="",
    ),
    "zenodo_n_traits_coloc": dict(
        source="S3", method="Zenodo trait COLOC count",
        ld_panel="external", gwas_scope="external",
        units="integer count", deprecated=False,
        provenance="external Zenodo import",
        notes="",
    ),
    # -------- S3 ieQTL (77-79) --------
    "ieqtl_disease_interaction": dict(
        source="S3", method="interaction eQTL disease-interaction term",
        ld_panel="see-source", gwas_scope="NA",
        units="boolean/significant flag", deprecated=False,
        provenance="ieQTL pipeline",
        notes="Disease-state-dependent eQTL.",
    ),
    "ieqtl_interaction_pval": dict(
        source="S3", method="ieQTL interaction p-value",
        ld_panel="see-source", gwas_scope="NA",
        units="p-value", deprecated=False,
        provenance="ieQTL pipeline",
        notes="",
    ),
    "ieqtl_cell_type": dict(
        source="S3", method="ieQTL cell-type assignment",
        ld_panel="see-source", gwas_scope="NA",
        units="cell-type string", deprecated=False,
        provenance="ieQTL pipeline",
        notes="",
    ),
    # -------- S3 pleiotropy (80-83) --------
    "pleiotropy_class": dict(
        source="S3", method="pleiotropy classification",
        ld_panel="multi-panel", gwas_scope="cross-disease",
        units="categorical", deprecated=False,
        provenance="pleiotropy_pipeline",
        notes="",
    ),
    "pleiotropy_n_traits": dict(
        source="S3", method="pleiotropy trait count",
        ld_panel="multi-panel", gwas_scope="cross-disease",
        units="integer count", deprecated=False,
        provenance="pleiotropy_pipeline",
        notes="",
    ),
    "pleiotropy_n_domains": dict(
        source="S3", method="pleiotropy disease-domain count",
        ld_panel="multi-panel", gwas_scope="cross-disease",
        units="integer count", deprecated=False,
        provenance="pleiotropy_pipeline",
        notes="",
    ),
    "pleiotropy_domains": dict(
        source="S3", method="pleiotropy disease-domain list",
        ld_panel="multi-panel", gwas_scope="cross-disease",
        units="domain-list string", deprecated=False,
        provenance="pleiotropy_pipeline",
        notes="",
    ),
    # -------- S3 HyPrColoc (188-190) --------
    "hyprcoloc_posterior": dict(
        source="S3", method="HyPrColoc joint posterior",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="posterior [0,1]", deprecated=False,
        provenance="HyPrColoc pipeline",
        notes="HyPrColoc = hypothesis-prioritized COLOC (Foley 2021).",
    ),
    "hyprcoloc_n_traits": dict(
        source="S3", method="HyPrColoc cluster size",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="integer count", deprecated=False,
        provenance="HyPrColoc pipeline",
        notes="",
    ),
    "hyprcoloc_traits": dict(
        source="S3", method="HyPrColoc cluster trait list",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="trait-list string", deprecated=False,
        provenance="HyPrColoc pipeline",
        notes="",
    ),
    # -------- DEPRECATED: MR retired 2026-04-22 (137-144) --------
    "mr_beta": dict(
        source="S3", method="MR causal beta estimate (RETIRED)",
        ld_panel="see-source", gwas_scope="see-source",
        units="MR causal effect estimate", deprecated=True,
        provenance="ARCHIVED archive/mr_ditched_2026-04-22/",
        notes="DEPRECATED. MR permanently ditched from paper 2026-04-22. Use TWAS + COLOC + INTACT instead.",
    ),
    "mr_bidirectional_pval": dict(
        source="S3", method="bidirectional MR p-value (RETIRED)",
        ld_panel="see-source", gwas_scope="see-source",
        units="p-value", deprecated=True,
        provenance="ARCHIVED",
        notes="DEPRECATED. See note on mr_beta.",
    ),
    "mr_ivw_beta": dict(
        source="S3", method="MR-IVW beta (RETIRED)",
        ld_panel="see-source", gwas_scope="see-source",
        units="MR causal effect estimate", deprecated=True,
        provenance="ARCHIVED",
        notes="DEPRECATED.",
    ),
    "mr_ivw_pval": dict(
        source="S3", method="MR-IVW p-value (RETIRED)",
        ld_panel="see-source", gwas_scope="see-source",
        units="p-value", deprecated=True,
        provenance="ARCHIVED",
        notes="DEPRECATED.",
    ),
    "mr_n_gwas_sig": dict(
        source="S3", method="MR significant GWAS count (RETIRED)",
        ld_panel="see-source", gwas_scope="all-portfolio",
        units="integer count", deprecated=True,
        provenance="ARCHIVED",
        notes="DEPRECATED.",
    ),
    "mr_n_instruments": dict(
        source="S3", method="MR genetic-instrument count (RETIRED)",
        ld_panel="see-source", gwas_scope="see-source",
        units="integer count", deprecated=True,
        provenance="ARCHIVED",
        notes="DEPRECATED.",
    ),
    "mr_pval": dict(
        source="S3", method="MR p-value (RETIRED)",
        ld_panel="see-source", gwas_scope="see-source",
        units="p-value", deprecated=True,
        provenance="ARCHIVED",
        notes="DEPRECATED.",
    ),
    "mr_sig": dict(
        source="S3", method="MR significance flag (RETIRED)",
        ld_panel="see-source", gwas_scope="see-source",
        units="boolean", deprecated=True,
        provenance="ARCHIVED",
        notes="DEPRECATED.",
    ),
    # -------- SuSiE backed flag (177) --------
    "susie_backed": dict(
        source="S3", method="SuSiE-backed credible-set indicator",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="boolean", deprecated=False,
        provenance="27a_assemble_evidence_atlas.R",
        notes="True if at least one SuSiE PP4 column populated (any LD panel).",
    ),
    # -------- NAFLD specificity / progression flags (84, 134-135, 165) --------
    "is_nafld_specific": dict(
        source="annotation", method="trait specificity",
        ld_panel="NA", gwas_scope="cross-disease",
        units="boolean", deprecated=False,
        provenance="27a_assemble_evidence_atlas.R",
        notes="True if COLOC signal is NAFLD/MASLD-specific vs cross-disease pleiotropic.",
    ),
    "is_onset_specific": dict(
        source="progression",
        method="onset-specific stratified classification",
        ld_panel="NA", gwas_scope="NA",
        units="boolean", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "is_progression_specific": dict(
        source="progression",
        method="progression-specific stratified classification",
        ld_panel="NA", gwas_scope="NA",
        units="boolean", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "prog_is_onset_specific": dict(
        source="progression", method="onset-specific (post-rebuild)",
        ld_panel="NA", gwas_scope="NA",
        units="boolean", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="Refresh after 5-cohort canonical rebuild 2026-05-01.",
    ),
    # -------- sex stratification (85-88, 174-175, 191) --------
    "sex_class": dict(
        source="sex", method="sex-dimorphic interaction classification",
        ld_panel="NA", gwas_scope="NA",
        units="categorical (Female_biased / Male_biased / Divergent / Concordant)",
        deprecated=False,
        provenance="26_sex_stratified_analysis.R",
        notes="Interaction-model output (2026-04-03). 3,080 sex-dimorphic; 1,978 F-biased + 675 M-biased + 427 Divergent + 31,373 Concordant. Old significance-based 7,886/300/1,021 scheme retired.",
    ),
    "sex_interaction_padj": dict(
        source="sex", method="sex-interaction term BH-adjusted",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="26_sex_stratified_analysis.R",
        notes="",
    ),
    "sex_interaction_padj_from_cls": dict(
        source="sex", method="sex-interaction padj (classifier-derived)",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="26_sex_stratified_analysis.R",
        notes="",
    ),
    "sex_progression_class": dict(
        source="sex", method="sex x progression interaction class",
        ld_panel="NA", gwas_scope="NA",
        units="categorical", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "sex_causal_score": dict(
        source="sex", method="sex-stratified causal evidence score",
        ld_panel="NA", gwas_scope="NA",
        units="composite score", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    # -------- pathway (89-90) --------
    "n_leading_edge_pathways": dict(
        source="pathway", method="fgsea leading-edge pathway count",
        ld_panel="NA", gwas_scope="NA",
        units="integer count", deprecated=False,
        provenance="downstream_analysis/pathway_analysis/",
        notes="MSigDB v2025.1. fgsea against Hallmark + C2 / C5.",
    ),
    "top_pathways": dict(
        source="pathway", method="fgsea top pathway list",
        ld_panel="NA", gwas_scope="NA",
        units="pathway-name string", deprecated=False,
        provenance="downstream_analysis/pathway_analysis/",
        notes="",
    ),
    # -------- S4 essentiality (91-93) --------
    "essentiality_chronos": dict(
        source="S4", method="DepMap CHRONOS essentiality score",
        ld_panel="NA", gwas_scope="NA",
        units="CHRONOS score (continuous)", deprecated=False,
        provenance="DepMap external download",
        notes="More negative = more essential. Liver-line-restricted version available via n_liver_lines.",
    ),
    "is_essential": dict(
        source="S4", method="DepMap pan-essentiality flag",
        ld_panel="NA", gwas_scope="NA",
        units="boolean", deprecated=False,
        provenance="DepMap external download",
        notes="True for pan-essential genes (safety liability if drug target).",
    ),
    "n_liver_lines": dict(
        source="S4", method="DepMap liver cell-line count",
        ld_panel="NA", gwas_scope="NA",
        units="integer count", deprecated=False,
        provenance="DepMap external download",
        notes="",
    ),
    # -------- drug repurposing (94-97, 167, 213-214) --------
    "attribution_class": dict(
        source="annotation", method="DEG cell-type attribution",
        ld_panel="NA", gwas_scope="NA",
        units="categorical (Hepatocyte_intrinsic / Composition_driven / Unmasked / ...)",
        deprecated=False,
        provenance="25_celltype_attribution.R (MuSiC + BayesPrism)",
        notes="MuSiC: 1,188 hepatocyte-intrinsic, 717 Composition_driven, 87 Unmasked.",
    ),
    "dgidb_druggable": dict(
        source="drug", method="DGIdb druggable-gene category",
        ld_panel="NA", gwas_scope="NA",
        units="categorical/string", deprecated=False,
        provenance="DGIdb external download",
        notes="",
    ),
    "opentargets_drug": dict(
        source="drug", method="OpenTargets drug-target evidence",
        ld_panel="NA", gwas_scope="NA",
        units="drug-name string", deprecated=False,
        provenance="OpenTargets external download",
        notes="",
    ),
    "lincs_reversal": dict(
        source="drug", method="LINCS L1000 reversal score",
        ld_panel="NA", gwas_scope="NA",
        units="reversal score (continuous)", deprecated=False,
        provenance="run_strategy11_v2.sh",
        notes="1,107 reversal compounds; 56-75% stage-specific.",
    ),
    "progression_drug_reversal_n": dict(
        source="drug", method="progression-stratified LINCS reversal count",
        ld_panel="NA", gwas_scope="NA",
        units="integer count", deprecated=False,
        provenance="run_strategy11_v2.sh",
        notes="",
    ),
    "drug_target_stratification": dict(
        source="drug", method="progression-stratified drug-target class",
        ld_panel="NA", gwas_scope="NA",
        units="categorical", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "lincs_reversal_stratum": dict(
        source="drug", method="LINCS stratum assignment",
        ld_panel="NA", gwas_scope="NA",
        units="categorical", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    # -------- proteomics (98-100) --------
    "n_prot_datasets": dict(
        source="proteomics", method="protein dataset coverage",
        ld_panel="NA", gwas_scope="NA",
        units="integer count", deprecated=False,
        provenance="proteomics_validation/",
        notes="",
    ),
    "best_protein_logFC": dict(
        source="proteomics", method="protein DE logFC (best dataset)",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="proteomics_validation/",
        notes="mRNA-protein concordance rho=0.317 overall / 0.563 Conserved.",
    ),
    "best_protein_padj": dict(
        source="proteomics", method="protein DE BH-adjusted",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="proteomics_validation/",
        notes="",
    ),
    # -------- S6 spatial (101-104, 212) --------
    "spatial_max_I": dict(
        source="S6", method="Visium Moran's I (best dataset)",
        ld_panel="NA", gwas_scope="NA",
        units="Moran's I statistic", deprecated=False,
        provenance="Analysis/Spatial/",
        notes="GSE192741 + Vu 2025 (validation).",
    ),
    "spatial_sig": dict(
        source="S6", method="Visium spatial significance flag",
        ld_panel="NA", gwas_scope="NA",
        units="boolean (BH-adjusted)", deprecated=False,
        provenance="Analysis/Spatial/",
        notes="",
    ),
    "spatial_n_datasets": dict(
        source="S6", method="Visium dataset coverage",
        ld_panel="NA", gwas_scope="NA",
        units="integer count", deprecated=False,
        provenance="Analysis/Spatial/",
        notes="",
    ),
    "layers_active": dict(
        source="annotation", method="multi-source evidence layer count",
        ld_panel="NA", gwas_scope="NA",
        units="integer count", deprecated=False,
        provenance="27a_assemble_evidence_atlas.R",
        notes="Number of evidence sources (S1-S7) flagging this gene.",
    ),
    "spatial_coloc_overlap": dict(
        source="S6", method="spatial x COLOC overlap flag",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="boolean", deprecated=False,
        provenance="gwas_spatial_convergence/",
        notes="",
    ),
    # -------- progression / contrasts (105-111, 126-133, 148-158, 168-173, 176, 217-219) --------
    "adv_fib_logFC": dict(
        source="progression", method="advanced-fibrosis contrast logFC",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="F3-F4 vs F0-F1 binary.",
    ),
    "adv_fib_padj": dict(
        source="progression", method="advanced-fibrosis padj",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="",
    ),
    "ballooning_ordinal_coef": dict(
        source="progression", method="ballooning ordinal regression coef",
        ld_panel="NA", gwas_scope="NA",
        units="ordinal coefficient", deprecated=False,
        provenance="15_progression_ordinal.R",
        notes="",
    ),
    "causal_methods_sig": dict(
        source="S3",
        method="multi-method causal-inference significance count",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="integer count (0-N methods)", deprecated=False,
        provenance="27a_assemble_evidence_atlas.R",
        notes="Counts {TWAS, COLOC-ABF, COLOC-SuSiE, INTACT, HyPrColoc, cTWAS}.",
    ),
    "causal_robustness": dict(
        source="S3", method="causal-inference robustness composite",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="composite score", deprecated=False,
        provenance="27a_assemble_evidence_atlas.R",
        notes="",
    ),
    "cirrhosis_logFC": dict(
        source="progression", method="cirrhosis contrast logFC",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="F4 vs F0 binary.",
    ),
    "cirrhosis_padj": dict(
        source="progression", method="cirrhosis contrast padj",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="",
    ),
    "early_late_nash_logFC": dict(
        source="progression", method="early-vs-late NASH contrast logFC",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="",
    ),
    "early_late_nash_padj": dict(
        source="progression", method="early-vs-late NASH padj",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="",
    ),
    "extreme_logFC": dict(
        source="progression", method="extreme-stage contrast logFC",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="F4 vs F0 extreme contrast.",
    ),
    "extreme_padj": dict(
        source="progression", method="extreme-stage contrast padj",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="",
    ),
    "f2_inflection_logFC": dict(
        source="progression", method="F2-inflection contrast logFC",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="ANCHOR for legacy F2-switch claim (RETIRED 2026-05-09). Multi-step cascade reframing: F0-F2 vs F3-F4 binary contrast; treat as one of four transitions, not a switch anchor.",
    ),
    "f2_inflection_padj": dict(
        source="progression", method="F2-inflection contrast padj",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="See f2_inflection_logFC notes.",
    ),
    "fibrosis_ordinal_coef": dict(
        source="progression", method="fibrosis ordinal coef (F0-F4)",
        ld_panel="NA", gwas_scope="NA",
        units="ordinal coefficient", deprecated=False,
        provenance="15_progression_ordinal.R",
        notes="",
    ),
    "inflammation_ordinal_coef": dict(
        source="progression", method="inflammation ordinal coef",
        ld_panel="NA", gwas_scope="NA",
        units="ordinal coefficient", deprecated=False,
        provenance="15_progression_ordinal.R",
        notes="",
    ),
    "nafl_vs_ctrl_logFC": dict(
        source="progression", method="NAFL-vs-Ctrl logFC",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="13_nafl_vs_nash_de.R",
        notes="",
    ),
    "nafl_vs_ctrl_padj": dict(
        source="progression", method="NAFL-vs-Ctrl padj",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="13_nafl_vs_nash_de.R",
        notes="",
    ),
    "nafl_vs_nash_logFC": dict(
        source="progression",
        method="NAFL-vs-NASH harmonized logFC",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="13_nafl_vs_nash_de.R",
        notes="Diagnosis-harmonized; NAS-based reclassification.",
    ),
    "nafl_vs_nash_padj": dict(
        source="progression", method="NAFL-vs-NASH padj",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="13_nafl_vs_nash_de.R",
        notes="",
    ),
    "nafl_vs_nash_tstat": dict(
        source="progression", method="NAFL-vs-NASH t-statistic",
        ld_panel="NA", gwas_scope="NA",
        units="t-statistic", deprecated=False,
        provenance="13_nafl_vs_nash_de.R",
        notes="",
    ),
    "nas_ge5_logFC": dict(
        source="progression", method="NAS>=5 contrast logFC",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="NAS Kleiner-2005 >=5 vs <5.",
    ),
    "nas_ge5_padj": dict(
        source="progression", method="NAS>=5 padj",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="",
    ),
    "nash_vs_ctrl_logFC": dict(
        source="progression", method="NASH-vs-Ctrl logFC",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="13_nafl_vs_nash_de.R",
        notes="",
    ),
    "nash_vs_ctrl_padj": dict(
        source="progression", method="NASH-vs-Ctrl padj",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="13_nafl_vs_nash_de.R",
        notes="",
    ),
    "nash_vs_nafl_fibadj_logFC": dict(
        source="progression",
        method="NASH-vs-NAFL fibrosis-adjusted logFC",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="13_nafl_vs_nash_de.R",
        notes="Fibrosis-adjusted on matched 615-sample subset. 90.2% (87.3-96.1%) DEG reduction vs unadjusted.",
    ),
    "nash_vs_nafl_fibadj_padj": dict(
        source="progression",
        method="NASH-vs-NAFL fibrosis-adjusted padj",
        ld_panel="NA", gwas_scope="NA",
        units="BH-adjusted p-value", deprecated=False,
        provenance="13_nafl_vs_nash_de.R",
        notes="",
    ),
    "progression_deconv_class": dict(
        source="progression",
        method="cell-type deconvolution progression class",
        ld_panel="NA", gwas_scope="NA",
        units="categorical", deprecated=False,
        provenance="118_bayesprism_deconv.R",
        notes="BayesPrism cell-type composition shift class.",
    ),
    "progression_gene_class": dict(
        source="progression",
        method="progression-trajectory gene classification",
        ld_panel="NA", gwas_scope="NA",
        units="categorical", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "progression_n_contrasts_sig": dict(
        source="progression",
        method="progression-contrast significance count",
        ld_panel="NA", gwas_scope="NA",
        units="integer count", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="",
    ),
    "progression_peak_contrast": dict(
        source="progression",
        method="progression peak-contrast identifier",
        ld_panel="NA", gwas_scope="NA",
        units="contrast-name string", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="",
    ),
    "progression_tau": dict(
        source="progression", method="progression Kendall's tau",
        ld_panel="NA", gwas_scope="NA",
        units="Kendall's tau", deprecated=False,
        provenance="14_progression_contrasts.R",
        notes="",
    ),
    "progression_twas_direction": dict(
        source="progression", method="progression-TWAS direction",
        ld_panel="multi-panel", gwas_scope="progression-subset",
        units="sign / direction", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "progression_twas_n_contrasts": dict(
        source="progression", method="progression-TWAS contrast count",
        ld_panel="multi-panel", gwas_scope="progression-subset",
        units="integer count", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "steatosis_ordinal_coef": dict(
        source="progression", method="steatosis ordinal coef",
        ld_panel="NA", gwas_scope="NA",
        units="ordinal coefficient", deprecated=False,
        provenance="15_progression_ordinal.R",
        notes="",
    ),
    "progression_coloc_class": dict(
        source="S3", method="progression-stratified COLOC class",
        ld_panel="multi-panel", gwas_scope="progression-subset",
        units="categorical", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "progression_causal_score": dict(
        source="S3", method="progression causal evidence score",
        ld_panel="multi-panel", gwas_scope="progression-subset",
        units="composite score", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "stage_specific_coloc_gwas": dict(
        source="S3", method="stage-specific COLOC GWAS list",
        ld_panel="multi-panel", gwas_scope="progression-subset",
        units="GWAS-list string", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    # -------- annotation classes (178-179) --------
    "ferroptosis_class": dict(
        source="annotation", method="ferroptosis-gene classification",
        ld_panel="NA", gwas_scope="NA",
        units="categorical", deprecated=False,
        provenance="ferroptosis_pipeline (Phase F)",
        notes="88 ferroptosis DEGs.",
    ),
    "zonation_class": dict(
        source="annotation", method="hepatocyte-zonation classification",
        ld_panel="NA", gwas_scope="NA",
        units="categorical (Periportal / Pericentral / Pan-lobular / Other)",
        deprecated=False,
        provenance="zonation_pipeline (Phase F)",
        notes="Periportal loss in disease; ADH4 lone pericentral causal locus.",
    ),
    # -------- NMF k=6 program scores (192-205) --------
    "prog_Innate_immune_inflammation_causal_score": dict(
        source="progression", method="NMF-program causal score (P2)",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="composite score", deprecated=False,
        provenance="217_stratified_causal_layer.R / 95_nmf_program_atlas.R",
        notes="NMF k=6 program P2 = Innate-immune.",
    ),
    "prog_Progression_Inflammatory_causal_score": dict(
        source="progression", method="NMF-program causal score (P1)",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="composite score", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="NMF k=6 program P1 = Pro-inflammatory.",
    ),
    "prog_Progression_Other_causal_score": dict(
        source="progression", method="NMF-program causal score (Other)",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="composite score", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "prog_Quiescent_Parenchyma_causal_score": dict(
        source="progression", method="NMF-program causal score (P3)",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="composite score", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="NMF k=6 program P3 = Parenchymal.",
    ),
    "prog_Stable_causal_score": dict(
        source="progression", method="NMF-program causal score (Stable)",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="composite score", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "prog_Vascular_smooth_muscle_pericyte_causal_score": dict(
        source="progression",
        method="NMF-program causal score (P6 stellate-myofibroblast)",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="composite score", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="NMF k=6 program P6 = Stellate-myofibroblast / vascular-smooth-muscle.",
    ),
    "prog_Innate_immune_inflammation_logFC": dict(
        source="progression", method="NMF-program logFC (P2)",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="95_nmf_program_atlas.R",
        notes="",
    ),
    "prog_Progression_Inflammatory_logFC": dict(
        source="progression", method="NMF-program logFC (P1)",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="95_nmf_program_atlas.R",
        notes="",
    ),
    "prog_Progression_Other_logFC": dict(
        source="progression", method="NMF-program logFC (Other)",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="95_nmf_program_atlas.R",
        notes="",
    ),
    "prog_Quiescent_Parenchyma_logFC": dict(
        source="progression", method="NMF-program logFC (P3)",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="95_nmf_program_atlas.R",
        notes="",
    ),
    "prog_Stable_logFC": dict(
        source="progression", method="NMF-program logFC (Stable)",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="95_nmf_program_atlas.R",
        notes="",
    ),
    "prog_Vascular_smooth_muscle_pericyte_logFC": dict(
        source="progression", method="NMF-program logFC (P6)",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="95_nmf_program_atlas.R",
        notes="",
    ),
    "dominant_program_for_gene": dict(
        source="progression", method="dominant NMF-program assignment",
        ld_panel="NA", gwas_scope="NA",
        units="program-name string", deprecated=False,
        provenance="95_nmf_program_atlas.R",
        notes="",
    ),
    "dominant_program_logFC": dict(
        source="progression", method="dominant-program logFC",
        ld_panel="NA", gwas_scope="NA",
        units="log2 fold change", deprecated=False,
        provenance="95_nmf_program_atlas.R",
        notes="",
    ),
    # -------- S2 fate / subtypes (206-209) --------
    "s2_causal_score": dict(
        source="progression",
        method="S2-fate (Progressor-female) causal evidence score",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="composite score", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="S2 fate = female-enriched Progressor signature (legacy k=2 framing, retired but kept as comparator).",
    ),
    "subtype_coloc_enriched": dict(
        source="S3", method="hepatocyte-subtype COLOC enrichment flag",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="boolean", deprecated=False,
        provenance="312_hepatocyte_subtype_coloc.R",
        notes="Disease-Associated subtype OR=4.50 at PP4>0.5.",
    ),
    "s2_genetic_risk_prob": dict(
        source="S3", method="S2-fate genetic-risk probability",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="probability [0,1]", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "progression_driver_validated": dict(
        source="S3", method="progression-driver validation flag",
        ld_panel="multi-panel", gwas_scope="progression-subset",
        units="boolean", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    # -------- zonation causal (210-211) --------
    "zone_causal_score": dict(
        source="S3", method="zonation causal evidence score",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="composite score", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="",
    ),
    "zone_coloc_enriched": dict(
        source="S3", method="zonation COLOC enrichment flag",
        ld_panel="multi-panel", gwas_scope="all-portfolio",
        units="boolean", deprecated=False,
        provenance="217_stratified_causal_layer.R",
        notes="Pericentral Hep-ATAC x Spatial_Pericentral OR = 3.86, p = 0.048 nominal.",
    ),
}


def auto_schema(df: pd.DataFrame) -> pd.DataFrame:
    """Auto-compute dtype, missingness, value range, n_unique per column."""
    rows = []
    for col in df.columns:
        s = df[col]
        n = len(s)
        n_missing = int(s.isna().sum())
        miss_frac = n_missing / n if n else 0.0
        dtype = str(s.dtype)
        n_unique = int(s.nunique(dropna=True))
        if pd.api.types.is_numeric_dtype(s):
            non_na = s.dropna()
            if len(non_na):
                vmin = float(non_na.min())
                vmax = float(non_na.max())
            else:
                vmin = vmax = float("nan")
            value_range = f"[{vmin:.4g}, {vmax:.4g}]"
            example = ""
        else:
            value_range = ""
            sample = s.dropna()
            example = str(sample.iloc[0])[:80] if len(sample) else ""
        rows.append(dict(
            column=col, dtype=dtype, n_missing=n_missing,
            missingness=round(miss_frac, 4), n_unique=n_unique,
            value_range=value_range, example=example,
        ))
    return pd.DataFrame(rows)


def merge_annotation(schema_df: pd.DataFrame) -> pd.DataFrame:
    """Merge auto-schema with hand-curated annotation."""
    out_rows = []
    for _, row in schema_df.iterrows():
        col = row["column"]
        annot = COL_ANNOT.get(col, {})
        out_rows.append(dict(
            column=col,
            dtype=row["dtype"],
            n_missing=row["n_missing"],
            missingness=row["missingness"],
            n_unique=row["n_unique"],
            value_range=row["value_range"],
            example=row["example"],
            source=annot.get("source", ""),
            method=annot.get("method", ""),
            ld_panel=annot.get("ld_panel", ""),
            gwas_scope=annot.get("gwas_scope", ""),
            units=annot.get("units", ""),
            deprecated=annot.get("deprecated", ""),
            provenance=annot.get("provenance", ""),
            ensembl_version_anchor="GENCODE v49 / GRCh38.p14",
            notes=annot.get("notes", ""),
        ))
    return pd.DataFrame(out_rows)


def main() -> None:
    print(f"Reading {ATLAS_PATH}")
    df = pd.read_csv(ATLAS_PATH, low_memory=False)
    print(f"Atlas dims: {df.shape[0]:,} rows x {df.shape[1]} cols")
    schema_df = auto_schema(df)
    full_df = merge_annotation(schema_df)

    n_total = len(full_df)
    n_annotated = int((full_df["source"] != "").sum())
    n_deprecated = int((full_df["deprecated"] == True).sum())  # noqa: E712
    print(f"Annotated: {n_annotated}/{n_total} ({n_annotated/n_total:.1%})")
    print(f"Deprecated (e.g., MR retired 2026-04-22): {n_deprecated}")
    unannotated = full_df.loc[full_df["source"] == "", "column"].tolist()
    if unannotated:
        print(f"UNANNOTATED columns ({len(unannotated)}): {unannotated[:10]} "
              f"{'...' if len(unannotated) > 10 else ''}")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    full_df.to_csv(OUT_PATH, sep="\t", index=False)
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
