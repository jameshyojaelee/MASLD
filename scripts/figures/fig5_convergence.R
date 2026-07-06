##############################################################################
# fig5_convergence.R (was fig6_evidence_matrix_v2.R; renamed 2026-04-15)
# Fig 5: Multi-Modal Evidence Matrix (ComplexHeatmap v2)
# Top genes ranked by causal evidence, showing continuous scores across
# all modalities: bulk RNA, multiple COLOC phenotypes, INTACT, TWAS, cTWAS,
# HyPrColoc, cross-species, GWAS-ATAC, progression, druggability
# Uses atlas v12 (156 cols) + convergence summary (INTACT)
# Output: figures/main/fig5_convergence/panel_evidence_matrix.pdf (FIG5_DIR)
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUTDIR <- FIG5_DIR
PANDIR <- file.path(OUTDIR, "panels")
dir.create(PANDIR, showWarnings = FALSE, recursive = TRUE)

# Layout selector. The default is the canonical evidence-ranked panel; the
# therapeutic_axes variant preserves the same selected genes/evidence but
# groups rows by one curated primary therapeutic interpretation.
FIG5A_LAYOUT <- Sys.getenv("FIG5A_LAYOUT", "standard")
VALID_LAYOUTS <- c("standard", "therapeutic_axes")
if (!FIG5A_LAYOUT %in% VALID_LAYOUTS) {
  stop(sprintf("Invalid FIG5A_LAYOUT='%s'; expected one of: %s",
               FIG5A_LAYOUT, paste(VALID_LAYOUTS, collapse = ", ")))
}
OUTPDF <- file.path(
  PANDIR,
  if (FIG5A_LAYOUT == "therapeutic_axes") "fig5a_therapeutic_axes.pdf" else "fig5a.pdf"
)

# Optional: render a truncated top-N variant with taller cells.
# `FIG5A_TOPN=20 Rscript fig5_convergence.R` -> panels/fig5a_top20.pdf
TOPN_ENV <- Sys.getenv("FIG5A_TOPN", "")
TOPN     <- if (nzchar(TOPN_ENV)) as.integer(TOPN_ENV) else NA_integer_
SQ_CELL_MM <- 8L  # uniform square cell size (mm) for the top-N variant
if (FIG5A_LAYOUT == "therapeutic_axes" && !is.na(TOPN)) {
  stop("FIG5A_LAYOUT=therapeutic_axes cannot be combined with FIG5A_TOPN")
}
if (!is.na(TOPN)) {
  OUTPDF <- file.path(PANDIR, sprintf("fig5a_top%d.pdf", TOPN))
}

# ── Per-phenotype-category max PP.H4 (ABF or SuSiE) ─────────────────────────
# Reads per-GWAS susie_coloc/<GWAS>/susie_coloc_chr*.csv (each row has both
# PP.H4.abf and PP.H4.susie for that gene × GWAS), takes the per-row max,
# classifies GWAS to category (enzyme / disease / pdff), and aggregates the
# max across GWAS within each category. Result is cached to a CSV; rebuild
# is triggered if the cache is older than any source file.
build_category_pp4 <- function() {
  cache_f    <- file.path(BASE, "RNA-seq/results/multi_evidence/susie_abf_max_pp4_by_category.csv")
  susie_root <- file.path(BASE, "GWAS/finemapping/results/susie_coloc")
  gwas_dirs  <- list.dirs(susie_root, recursive = FALSE, full.names = TRUE)
  gwas_dirs  <- gwas_dirs[!grepl("/archived", gwas_dirs)]

  src_files <- list.files(susie_root, recursive = TRUE, full.names = TRUE,
                          pattern = "^susie_coloc_chr.*\\.csv$")
  if (file.exists(cache_f) && length(src_files) > 0) {
    if (file.info(cache_f)$mtime > max(file.info(src_files)$mtime, na.rm = TRUE)) {
      cat("  Using cached per-category PP4:", cache_f, "\n")
      return(fread(cache_f))
    }
  }

  classify_gwas <- function(name) {
    if (grepl("ALT|AST|GGT", name, ignore.case = TRUE)) return("enzyme")
    if (grepl("PDFF",       name, ignore.case = TRUE)) return("pdff")
    if (grepl("NAFLD|NASH|HCC|Cirrhosis", name, ignore.case = TRUE)) return("disease")
    NA_character_
  }

  cat("  Building per-category PP4 from", length(gwas_dirs), "GWAS dirs ...\n")
  per_gwas <- list()
  for (d in gwas_dirs) {
    gwas_name <- basename(d)
    cat_v     <- classify_gwas(gwas_name)
    if (is.na(cat_v)) next
    chr_files <- list.files(d, pattern = "^susie_coloc_chr.*\\.csv$", full.names = TRUE)
    if (length(chr_files) == 0) next
    dt <- tryCatch(rbindlist(lapply(chr_files, fread), fill = TRUE),
                   error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) next
    if (!"PP.H4.abf" %in% names(dt)) next
    abf <- suppressWarnings(as.numeric(dt$PP.H4.abf))
    sus <- if ("PP.H4.susie" %in% names(dt)) suppressWarnings(as.numeric(dt$PP.H4.susie)) else rep(NA_real_, nrow(dt))
    dt[, max_pp4 := pmax(abf, sus, na.rm = TRUE)]
    dt[is.infinite(max_pp4), max_pp4 := NA_real_]
    dt[, category := cat_v]
    per_gwas[[gwas_name]] <- dt[!is.na(gene), .(gene, category, max_pp4)]
  }

  agg <- rbindlist(per_gwas, fill = TRUE)
  out <- agg[, .(max_pp4 = suppressWarnings(max(max_pp4, na.rm = TRUE))),
             by = .(gene, category)]
  out[is.infinite(max_pp4), max_pp4 := NA_real_]
  out_wide <- dcast(out, gene ~ category, value.var = "max_pp4")
  setnames(out_wide,
           intersect(c("enzyme", "disease", "pdff"), names(out_wide)),
           paste0("category_pp4_", intersect(c("enzyme", "disease", "pdff"), names(out_wide))))
  fwrite(out_wide, cache_f)
  cat("  Cached per-category PP4 to", cache_f, "(", nrow(out_wide), "genes)\n")
  out_wide
}

# ── Load data ────────────────────────────────────────────────────────────────
cat("Loading atlas ...\n")
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
cat(sprintf("  Atlas: %d genes x %d cols\n", nrow(atlas), ncol(atlas)))

# Merge per-category max PP4 (ABF or SuSiE) into atlas
category_pp4 <- build_category_pp4()
atlas <- merge(atlas, category_pp4, by.x = "human_symbol", by.y = "gene", all.x = TRUE)
cat(sprintf("  Per-category PP4 merged (%d genes)\n",
            sum(!is.na(atlas$category_pp4_enzyme) |
                !is.na(atlas$category_pp4_disease) |
                !is.na(atlas$category_pp4_pdff))))

# INTACT load removed 2026-06-19 (INTACT dropped; genetic signal = COLOC via rs_coloc/rs_susie).

# Load progression driver data
prog_file <- file.path(BASE, "RNA-seq/results/stratified_causal/progression_driver_genetics.csv")
if (file.exists(prog_file)) {
  prog <- fread(prog_file)
  prog_dedup <- prog[!duplicated(gene_symbol), .(gene_symbol, progression_causal_score_driver = progression_causal_score,
                                                    genetically_validated)]
  atlas <- merge(atlas, prog_dedup, by.x = "human_symbol", by.y = "gene_symbol", all.x = TRUE)
  cat("  Progression driver scores merged:", sum(!is.na(atlas$progression_causal_score_driver)), "genes\n")
}

# Load broader progression classification (5,121 genes with onset/progression/pan_stage)
prog_class_file <- file.path(BASE, "RNA-seq/results/stratified_causal/progression_coloc_by_transition.csv")
if (file.exists(prog_class_file)) {
  prog_class <- fread(prog_class_file)
  prog_class_dedup <- prog_class[!duplicated(symbol), .(symbol, progression_class,
                                                          n_stages_sig, max_abs_logfc)]
  atlas <- merge(atlas, prog_class_dedup, by.x = "human_symbol", by.y = "symbol", all.x = TRUE)
  cat("  Progression classification merged:", sum(!is.na(atlas$progression_class)), "genes\n")
  # Create composite progression score: use causal score if available, else derive from class
  # progression_class genes get a scaled score based on class + n_stages_sig
  atlas[, progression_composite := fifelse(
    !is.na(progression_causal_score_driver), progression_causal_score_driver,
    fifelse(progression_class == "progression", 0.3 * pmin(n_stages_sig / 4, 1),
    fifelse(progression_class == "pan_stage", 0.15 * pmin(n_stages_sig / 4, 1),
    fifelse(progression_class == "onset", 0.05, NA_real_))))]
  cat("  Composite progression score:", sum(!is.na(atlas$progression_composite)), "genes\n")
} else {
  # Fallback: use driver score directly as composite
  atlas[, progression_composite := progression_causal_score_driver]
}

# Load SuSiE-COLOC gene-level (may not be in atlas after rebuild)
susie_file <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
if (file.exists(susie_file) && !"coloc_susie_best_pp4" %in% names(atlas)) {
  susie <- fread(susie_file)
  susie <- susie[gene != "" & !is.na(gene), .(gene, coloc_susie_best_pp4 = coloc_best_pp4)]
  susie <- susie[!duplicated(gene)]
  atlas <- merge(atlas, susie, by.x = "human_symbol", by.y = "gene", all.x = TRUE)
  cat("  SuSiE-COLOC PP.H4 merged\n")
}

# Load GWAS-ATAC gene-level (may not be in atlas)
gwas_atac_file <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/gene_level_gwas_atac.csv")
if (file.exists(gwas_atac_file) && !"gwas_max_pip_in_peak" %in% names(atlas)) {
  ga <- fread(gwas_atac_file)
  if ("assigned_gene" %in% names(ga)) {
    ga_dedup <- ga[assigned_gene != "" & !is.na(assigned_gene)]
    ga_dedup <- ga_dedup[!duplicated(assigned_gene)]
    atlas <- merge(atlas, ga_dedup, by.x = "human_symbol", by.y = "assigned_gene", all.x = TRUE)
    cat("  GWAS-ATAC gene-level merged\n")
  }
}

# Load Proteomics validation
# Direction concordance with bulk RNA is NOT required — discordance is valid biology
# (post-transcriptional regulation, secretion, protein stability, etc.)
if (!"best_protein_logFC" %in% names(atlas)) {
  prot_file <- file.path(BASE, "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv")
  if (file.exists(prot_file)) {
    prot <- fread(prot_file)
    gene_prot <- prot[protein_padj < 0.05,
      .(n_prot_datasets = uniqueN(dataset),
        best_protein_logFC = protein_logFC[which.min(protein_padj)],
        best_protein_padj = min(protein_padj)),
      by = gene]
    gene_prot <- gene_prot[!duplicated(gene)]
    atlas <- merge(atlas, gene_prot, by.x = "human_symbol", by.y = "gene", all.x = TRUE)
  }
}
cat("  Proteomics validation merged:", sum(!is.na(atlas$best_protein_logFC)), "genes\n")

# Load Spatial validation (Moran's I from both Visium datasets)
sp_guilliams <- file.path(BASE, "Analysis/Spatial/results/validation_bulk/spatial_bulk_merged_Guilliams_et_al.csv")
sp_vu <- file.path(BASE, "Analysis/Spatial/results/validation_bulk/spatial_bulk_merged_Vu_et_al.csv")
if (file.exists(sp_guilliams)) {
  sg <- fread(sp_guilliams)[, .(symbol, morans_I_guilliams = I, spatial_padj_guilliams = padj_bh)]
  sg <- sg[!duplicated(symbol)]
  atlas <- merge(atlas, sg, by.x = "human_symbol", by.y = "symbol", all.x = TRUE)
  cat("  Spatial (Guilliams) merged:", sum(!is.na(atlas$morans_I_guilliams)), "genes\n")
}
if (file.exists(sp_vu)) {
  sv <- fread(sp_vu)[, .(symbol, morans_I_vu = I, spatial_padj_vu = padj_bh)]
  sv <- sv[!duplicated(symbol)]
  atlas <- merge(atlas, sv, by.x = "human_symbol", by.y = "symbol", all.x = TRUE)
  cat("  Spatial (Vu) merged:", sum(!is.na(atlas$morans_I_vu)), "genes\n")
}
# Composite spatial metric
if ("morans_I_guilliams" %in% names(atlas)) {
  atlas[, spatial_max_I := pmax(morans_I_guilliams, morans_I_vu, na.rm = TRUE)]
  atlas[, spatial_sig := (spatial_padj_guilliams < 0.05 | spatial_padj_vu < 0.05)]
  atlas[is.na(spatial_sig), spatial_sig := FALSE]
}

# Load scATAC DA results (hepatocyte DA peaks annotated to nearest gene)
da_file <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/l8_annotated_corrected/scatac_da_gene_annotated.csv")
if (file.exists(da_file)) {
  da <- fread(da_file)
  da_gene <- da[padj < 0.05, .(n_da_peaks = .N,
    best_da_logFC = logFC[which.min(padj)],
    best_da_padj = min(padj)), by = gene_symbol]
  da_gene <- da_gene[!duplicated(gene_symbol)]
  atlas <- merge(atlas, da_gene, by.x = "human_symbol", by.y = "gene_symbol", all.x = TRUE)
  cat("  DA peaks merged:", sum(!is.na(atlas$n_da_peaks)), "genes\n")
}

# Load SCENIC+ disease regulon targets (is gene regulated by a disease-relevant TF?)
dis_reg_file <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv")
if (file.exists(dis_reg_file)) {
  dis_reg <- fread(dis_reg_file)
  # Guard (2026-06-01): disease_regulons.csv may be header-only (empty) -> target_genes
  # reads as logical and strsplit() errors. Generate "with what we got": skip the
  # regulon annotation if empty/non-character, leaving is_disease_reg_target = FALSE.
  if (nrow(dis_reg) > 0 && is.character(dis_reg$target_genes)) {
    reg_targets <- dis_reg[, .(target = unlist(strsplit(target_genes, ";"))), by = tf_name]
    reg_targets <- reg_targets[target != ""]
    reg_gene <- reg_targets[, .(is_disease_reg_target = TRUE,
      n_regulating_tfs = uniqueN(tf_name),
      regulating_tfs = paste(unique(tf_name), collapse = ";")), by = target]
    reg_gene <- reg_gene[!duplicated(target)]
    atlas <- merge(atlas, reg_gene, by.x = "human_symbol", by.y = "target", all.x = TRUE)
    atlas[is.na(is_disease_reg_target), is_disease_reg_target := FALSE]
    cat("  Disease regulon targets merged:", sum(atlas$is_disease_reg_target), "genes\n")
  } else {
    atlas[, is_disease_reg_target := FALSE]
    cat("  Disease regulon file empty/non-character; skipping regulon annotation (generating with available data)\n")
  }
}

# Disrupted master-regulator targets (the dense ATAC/regulatory column, 2026-06-22):
# # of disease master-regulators that are motif-disrupted by fine-mapped MASLD risk
# variants (fig4 disease_master_regulators) AND regulate this gene (CollecTRI).
# Precomputed by scripts/figures/compute_disrupted_mr_targets.R.
mrt_file <- file.path(BASE, "RNA-seq/results/multi_evidence/disrupted_mr_targets.csv")
if (file.exists(mrt_file)) {
  mrt <- fread(mrt_file)[, .(human_symbol, n_disrupted_mr_targets)]
  atlas <- merge(atlas, mrt, by = "human_symbol", all.x = TRUE)
  atlas[is.na(n_disrupted_mr_targets), n_disrupted_mr_targets := 0L]
  cat("  Disrupted master-regulator targets merged:", sum(atlas$n_disrupted_mr_targets > 0), "genes\n")
}

# Clinical drug-development status (replaces binary Druggable, 2026-06-22): the
# atlas already carries drug_dev_status (same source as fig4c / fig5c) — use it.
cat("  Clinical drug-dev status (atlas):",
    sum(atlas$drug_dev_status %in% c("masld_approved","masld_clinical","masld_preclinical","masld_discontinued")),
    "MASLD-pipeline genes\n")

# Load clinical drug validation for target genes — prefer v2 (2026-04-24+).
drug_val_v2 <- file.path(BASE, "RNA-seq/results/drug_repurposing/clinical_drug_validation_table_v2.csv")
drug_val_v1 <- file.path(BASE, "RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv")
drug_val_file <- if (file.exists(drug_val_v2)) drug_val_v2 else drug_val_v1
drug_val_best <- NULL
if (file.exists(drug_val_file)) {
  drug_val <- fread(drug_val_file)
  # Map atlas_support to ordinal: Strong=3, Moderate=2, Weak=1, Absent=0
  drug_val[, clin_level := fcase(
    atlas_support == "Strong", 3L,
    atlas_support == "Moderate", 2L,
    atlas_support == "Weak", 1L,
    default = 0L)]
  # Keep best evidence per target gene
  drug_val_best <- drug_val[, .(clin_level = max(clin_level, na.rm = TRUE),
                                 drug_name = drug[which.max(clin_level)],
                                 drug_stage = stage[which.max(clin_level)]),
                             by = target_gene]
  cat("  Clinical drug validation loaded:", nrow(drug_val_best), "target genes\n")
}

# Load pathway functional annotation
pathway_file <- file.path(BASE, "RNA-seq/results/multi_evidence/functional_activity/gene_functional_annotation.csv")
if (file.exists(pathway_file)) {
  pw <- fread(pathway_file, select = c("human_symbol", "n_pathway_memberships", "functional_score"))
  pw <- pw[!duplicated(human_symbol)]
  atlas <- merge(atlas, pw, by = "human_symbol", all.x = TRUE)
  cat("  Pathway memberships merged:", sum(!is.na(atlas$n_pathway_memberships)), "genes\n")
}

# ── Gene selection ───────────────────────────────────────────────────────────
# Must be DEG + have at least some causal evidence
# Rank by composite: max COLOC PP.H4 + INTACT + bulk DEG significance
# ARCHIVED 2026-04-08: FinnGen COLOC columns removed
coloc_pp4_cols <- intersect(c("broadaway_coloc_pp4", "ukbb_alt_coloc_pp4",
  "ast_coloc_pp4", "ggt_coloc_pp4", "pdff_coloc_pp4",
  "bbj_alt_coloc_pp4"), names(atlas))

atlas[, max_coloc_pp4 := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = coloc_pp4_cols]

# Composite ranking score
# Rank using available columns (some may not exist after atlas rebuild)
# Build score component by component
atlas[, rs_coloc := fifelse(is.na(max_coloc_pp4), 0, max_coloc_pp4) * 3]
atlas[, rs_deg := fifelse(is.na(bulk_treat_fdr) | bulk_treat_fdr >= 0.05, 0, -log10(pmax(bulk_padj, 1e-300)) / 10)]  # TREAT canonical; effect floor (lfc=0.25) IS in the test
if ("coloc_susie_best_pp4" %in% names(atlas)) {
  atlas[, rs_susie := fifelse(is.na(coloc_susie_best_pp4), 0, coloc_susie_best_pp4)]
} else {
  atlas[, rs_susie := 0]
}
# Progression component: use composite score (causal score or class-derived)
if ("progression_composite" %in% names(atlas)) {
  atlas[, rs_prog := fifelse(is.na(progression_composite), 0, progression_composite)]
} else {
  atlas[, rs_prog := 0]
}
atlas[, rank_score := rs_coloc + rs_deg + rs_susie + rs_prog]
atlas[, c("rs_coloc", "rs_deg", "rs_susie", "rs_prog") := NULL]

# ── Merge 46d convergence evidence ──────────────────────────────────────────
# Coexists with the heuristic rank_score above. 46d's empirically-calibrated
# posterior is used as the primary tie-breaker within each n_convergence level
# (replacing the legacy 46b archetype score), and tier + concordance state are
# rendered as left-annotation stripes on the heatmap.
bayev_f <- file.path(BASE, "RNA-seq/results/multi_evidence/convergence_evidence.csv")
if (file.exists(bayev_f)) {
  bayev <- fread(bayev_f, select = c("human_symbol", "convergence_score",
                                       "convergence_rank", "tier",
                                       "concordance_state", "dominant_stage_S1"))
  setnames(bayev, c("convergence_score","convergence_rank","tier",
                     "concordance_state","dominant_stage_S1"),
                  c("bayev_score","bayev_rank","bayev_tier",
                     "bayev_state","bayev_stage"))
  atlas <- merge(atlas, bayev, by = "human_symbol", all.x = TRUE)
  cat(sprintf("  46d convergence evidence merged (%d genes have score)\n",
              sum(!is.na(atlas$bayev_score))))
} else {
  cat("  WARN: 46d output missing; falling back to legacy rank_score\n")
  atlas[, bayev_score := NA_real_]
  atlas[, bayev_rank  := NA_integer_]
  atlas[, bayev_tier  := NA_character_]
  atlas[, bayev_state := NA_character_]
  atlas[, bayev_stage := NA_character_]
}

# ── Build GWAS-ATAC gene-level motif disruption from variant annotations ─────
# This maps variant-level motif disruptions to nearest genes (184 genes total)
ann_file <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv")
motif_file <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv")
atac_gene_dt <- data.table()

if (file.exists(ann_file) && file.exists(motif_file)) {
  ann_atac <- fread(ann_file)
  motif_atac <- fread(motif_file)
  # Map variants to nearest genes via hepatocyte peak overlaps
  hep_ann <- ann_atac[cell_type == "Hepatocytes" & !is.na(nearest_gene) & nearest_gene != ""]
  var_gene <- unique(hep_ann[, .(variant_id, nearest_gene, max_pip)])
  # Merge motif hits
  mg <- merge(motif_atac[, .(SNP_id, tf_name, motif_in_disease_regulon, alleleDiff)],
              var_gene, by.x = "SNP_id", by.y = "variant_id", allow.cartesian = TRUE)
  atac_gene_dt <- mg[, .(
    n_motifs_disrupted = uniqueN(tf_name),
    n_regulon_hits     = sum(motif_in_disease_regulon, na.rm = TRUE),
    regulon_tfs        = paste(unique(tf_name[motif_in_disease_regulon == TRUE]), collapse = ";"),
    max_pip_in_peak    = max(max_pip, na.rm = TRUE)
  ), by = nearest_gene]
  cat("  GWAS-ATAC: ", nrow(atac_gene_dt), "genes with motif disruption data\n")
}

# Notable genes: include GWAS-ATAC hits with regulon disruptions
notable <- c("THRB", "HKDC1", "RORA", "GCKR", "TM6SF2", "SLC39A8", "DGAT2",
             "CFLAR", "SYNJ2", "MTTP", "EPHA2", "FABP1", "HNF1A", "CDK6",
             "ADH4", "SGCB", "PNPLA3", "HSD17B13", "NR1H4", "PPARA",
             "SORT1", "CELSR2", "ACTG1", "GNMT",
             # GWAS-ATAC regulon disruption genes
             "HNF1B", "TM4SF4", "CYP26A1", "ALDH2", "TRIB1AL", "CYFIP2",
             # Druggable DEGs
             "AKR1B10", "LPL", "MMP9", "OLR1", "SPP1", "FAP", "COL1A1",
             "FMO1", "CXCL10", "ACHE",
             # Proteomics-validated + spatial
             "PLIN2", "ANXA2", "CYP2E1", "NNMT", "CCL19", "ADH4",
             # Top progression drivers (causal score > 0.1)
             "ADCY8", "ID4", "FCN2", "LRRC1", "SH3YL1", "SEMA3G",
             "HSBP1L1", "CDH16", "GLS", "MGP",
             # Clinical drug target genes (from clinical_drug_validation_table)
             "FGFR1", "GCGR", "GLP1R", "KLB", "PPARD", "PPARG", "SCD", "LGALS3", "CCR2",
             "HAL", "HGD", "ETNPPL", "C1QC", "CORO1A", "COL14A1", "EMILIN1",
             "FASN", "KHK", "ACACA", "ACACB", "GHR", "HMGCR", "MAP3K5", "CNR1", "ADGRE5", "LOXL2", "ADORA3", "FGF21")

# Determine DEG filter column (canonical: bulk_treat_fdr<0.05; TREAT lfc=0.25, effect floor in test -> NO separate |shrunk_logFC| filter)
deg_filter <- !is.na(atlas$bulk_treat_fdr) & atlas$bulk_treat_fdr < 0.05

# Gene selection: ALL notable genes guaranteed, remaining slots filled by rank
# 1. Get all notable genes first (bypass DEG filter)
notable_rows <- atlas[human_symbol %in% notable & human_symbol != "" & !is.na(human_symbol)]
notable_rows <- notable_rows[!duplicated(human_symbol)]
# 2. Fill remaining slots with top-ranked DEGs not already in notable
max_genes <- 160
n_fill <- max(0, max_genes - nrow(notable_rows))
top_ranked <- head(atlas[human_symbol != "" & !is.na(human_symbol) &
                           deg_filter & !(human_symbol %in% notable_rows$human_symbol)
                          ][order(-bayev_score, -rank_score, na.last = TRUE)], n_fill)
top_genes <- rbind(notable_rows, top_ranked)
top_genes <- top_genes[!duplicated(human_symbol)]

# Merge GWAS-ATAC data (before computing convergence)
if (nrow(atac_gene_dt) > 0) {
  top_genes <- merge(top_genes, atac_gene_dt, by.x = "human_symbol", by.y = "nearest_gene", all.x = TRUE)
}

# Define DEG contrast list here so it is available for both convergence scoring
# and matrix construction below.
deg_contrast_defs <- list(
  list(lfc = "bulk_shrunk_logFC", padj = "bulk_treat_fdr",   lfc_thresh = 0.0, label = "MASLD\nvs Ctrl"),  # canonical: bulk_treat_fdr<0.05 (TREAT lfc=0.25; effect floor IN the test). lfc kept for z-score display only.
  list(lfc = "nafl_vs_ctrl_logFC", padj = "nafl_vs_ctrl_padj", lfc_thresh = 0.0, label = "MAFL\nvs Ctrl"),
  list(lfc = "nash_vs_ctrl_logFC", padj = "nash_vs_ctrl_padj", lfc_thresh = 0.0, label = "MASH\nvs Ctrl"),
  list(lfc = "nafl_vs_nash_logFC", padj = "nafl_vs_nash_padj", lfc_thresh = 0.0, label = "MASH\nvs MAFL"),  # coeff=nafl_nashNASH: positive = up in MASH
  list(lfc = "adv_fib_logFC",      padj = "adv_fib_padj",      lfc_thresh = 0.0, label = "Adv Fib\n(F3-F4)"),
  list(lfc = "cirrhosis_logFC",    padj = "cirrhosis_padj",    lfc_thresh = 0.0, label = "Cirrhosis\n(F4)")
)

# ── Compute convergence score (N independent modalities with evidence) ──────
# Max = 5. Each modality counted once (0 or 1).
# Progression is excluded — it is derived from RNA, not an independent source.
# Mouse is excluded from the convergence count per user request.
top_genes[, n_convergence := 0L]

# 1. Human DEG: significant in ≥1 of the 6 RNA contrasts
deg_any_sig <- Reduce(`|`, lapply(deg_contrast_defs, function(d) {
  if (!d$lfc %in% names(top_genes)) return(rep(FALSE, nrow(top_genes)))
  lfc  <- as.numeric(top_genes[[d$lfc]])
  padj <- if (d$padj %in% names(top_genes)) as.numeric(top_genes[[d$padj]]) else rep(1, top_n)
  !is.na(padj) & padj < 0.05 & abs(lfc) > d$lfc_thresh
}))
top_genes[, n_convergence := n_convergence + as.integer(deg_any_sig)]

# 2. Genetic causal (any COLOC PP.H4 > 0.8)
top_genes[, n_convergence := n_convergence + fifelse(!is.na(max_coloc_pp4) & max_coloc_pp4 > 0.8, 1L, 0L)]
# 3. Regulatory (target of >=1 motif-disrupted master-regulator — matches the
#    displayed "Targets of disrupted TFs" column).
if ("n_disrupted_mr_targets" %in% names(top_genes)) {
  top_genes[, n_convergence := n_convergence + fifelse(!is.na(n_disrupted_mr_targets) & n_disrupted_mr_targets > 0, 1L, 0L)]
}
# 4. Proteomics validation
if ("best_protein_logFC" %in% names(top_genes)) {
  top_genes[, n_convergence := n_convergence + fifelse(!is.na(best_protein_logFC), 1L, 0L)]
}
# 5. Spatial
if ("spatial_max_I" %in% names(top_genes)) {
  top_genes[, n_convergence := n_convergence + fifelse(!is.na(spatial_max_I) & spatial_max_I > 0.05, 1L, 0L)]
}

# Sort by convergence (descending), then by 46d convergence_score, falling back
# to the heuristic rank_score for genes without a 46d posterior.
top_genes <- top_genes[order(-n_convergence, -bayev_score, -rank_score,
                              na.last = TRUE)]

label_genes <- c("THRB", "HKDC1", "RORA", "DGAT2", "PPARA", "NR1H4", "PNPLA3",
                 "HSD17B13", "GCKR",
                 "MMP9", "FGFR1", "GLP1R", "PPARG", "KLB",
                 "SCD", "LPL", "EPHA2", "AKR1B10", "FAP",
                 "GAS6", "TM4SF4", "CDK6", "LGALS3", "CCR2",
                 "HAL", "HGD", "ETNPPL", "C1QC", "CORO1A", "COL14A1", "EMILIN1",
                 "FASN", "KHK", "ACACA", "ACACB", "GHR", "HMGCR", "MAP3K5", "CNR1", "ADGRE5", "LOXL2", "ADORA3", "FGF21")

if (FIG5A_LAYOUT == "therapeutic_axes") {
  axis_file <- file.path(BASE, "scripts/figures/fig5a_therapeutic_axes.tsv")
  if (file.exists(axis_file)) {
    label_genes <- unique(c(label_genes, fread(axis_file)$gene))
  }
}

# Readability thinning (2026-06-22): the regulatory modality became dense, so the
# conv==4 tier ballooned (~37 genes) and crowded the matrix. Keep only the genuinely
# high-convergence tier (conv>=5) plus every labeled paper anchor at any tier.
top_genes <- top_genes[n_convergence >= 3 | human_symbol %in% label_genes]
top_genes <- top_genes[n_convergence >= 4 | human_symbol %in% label_genes]

# Optional top-N truncation (genes already sorted by convergence/46d/rank above)
if (!is.na(TOPN)) {
  top_genes <- head(top_genes, TOPN)
  cat(sprintf("FIG5A_TOPN set: truncating to top %d genes\n", TOPN))
}

# Optional row grouping by a single curated primary therapeutic axis. These
# categories are interpretive therapeutic bins, not claims of exclusive
# molecular-pathway membership. Preserve the existing convergence/rank order
# within each axis.
axis_levels <- c(
  "Metabolic dysfunction & lipid handling",
  "Hepatocellular injury & inflammation",
  "Fibrogenesis & matrix remodeling"
)
axis_display <- c(
  "Metabolic dysfunction\n& lipid handling",
  "Hepatocellular injury\n& inflammation",
  "Fibrogenesis\n& matrix remodeling"
)
axis_colors <- structure(
  c("#E69F00", "#C9265E", "#7B1FA2"),
  names = axis_levels
)

if (FIG5A_LAYOUT == "therapeutic_axes") {
  axis_file <- file.path(BASE, "scripts/figures/fig5a_therapeutic_axes.tsv")
  if (!file.exists(axis_file)) stop("Missing therapeutic-axis annotation: ", axis_file)

  axis_annot <- fread(axis_file)
  required_axis_cols <- c("gene", "primary_therapeutic_axis", "rationale")
  missing_axis_cols <- setdiff(required_axis_cols, names(axis_annot))
  if (length(missing_axis_cols)) {
    stop("Therapeutic-axis annotation is missing columns: ",
         paste(missing_axis_cols, collapse = ", "))
  }
  duplicate_axis_genes <- unique(axis_annot$gene[duplicated(axis_annot$gene)])
  if (length(duplicate_axis_genes)) {
    stop("Genes assigned more than once in therapeutic-axis annotation: ",
         paste(duplicate_axis_genes, collapse = ", "))
  }
  invalid_axes <- setdiff(unique(axis_annot$primary_therapeutic_axis), axis_levels)
  if (length(invalid_axes)) {
    stop("Invalid primary therapeutic axes: ", paste(invalid_axes, collapse = ", "))
  }
  top_genes <- top_genes[human_symbol %in% axis_annot$gene]
  unassigned_genes <- setdiff(top_genes$human_symbol, axis_annot$gene)
  if (length(unassigned_genes)) {
    stop("Displayed genes lack a primary therapeutic-axis assignment: ",
         paste(unassigned_genes, collapse = ", "))
  }

  top_genes[, primary_therapeutic_axis :=
    axis_annot$primary_therapeutic_axis[match(human_symbol, axis_annot$gene)]]

  # Filter out genes missing MASLD Therapeutics info (not undetermined/other)
  valid_statuses <- c("masld_approved", "masld_clinical", "masld_discontinued", "masld_preclinical", "discovery")
  top_genes <- top_genes[drug_dev_status %in% valid_statuses]

  # Prioritize approved, clinical, discontinued (failed), and discovery (novel) genes
  top_genes[, is_priority_status := drug_dev_status %in% c("masld_approved", "masld_clinical", "masld_discontinued", "discovery")]

  # Sort by priority status first, then by convergence and other scores
  setorder(top_genes, -is_priority_status, -n_convergence, -bayev_score, -rank_score, na.last = TRUE)

  # Equalize counts so that the difference is within +/- 3 genes
  counts <- table(top_genes$primary_therapeutic_axis)
  min_cnt <- min(counts)
  max_allowed <- min_cnt + 3

  # Take top max_allowed within each category
  top_genes_list <- split(top_genes, top_genes$primary_therapeutic_axis)
  top_genes_list <- lapply(top_genes_list, function(sub_dt) {
    head(sub_dt, max_allowed)
  })
  top_genes <- rbindlist(top_genes_list)

  top_genes[, is_priority_status := NULL]

  top_genes[, primary_therapeutic_axis :=
    factor(primary_therapeutic_axis, levels = axis_levels)]
  setorder(top_genes, primary_therapeutic_axis, -n_convergence, -bayev_score, -rank_score, na.last = TRUE)

  cat("Therapeutic-axis groups:\n")
  print(table(top_genes$primary_therapeutic_axis, useNA = "ifany"))
}

top_n <- nrow(top_genes)
genes <- top_genes$human_symbol
convergence_scores <- top_genes$n_convergence
therapeutic_axis_split <- if (FIG5A_LAYOUT == "therapeutic_axes") {
  factor(as.character(top_genes$primary_therapeutic_axis),
         levels = axis_levels, labels = axis_display)
} else {
  NULL
}
cat(sprintf("Selected %d genes (convergence >= 3; %d notable forced in)\n", top_n, sum(notable_rows$human_symbol %in% genes)))
cat(sprintf("  Convergence range: %d - %d\n", min(convergence_scores), max(convergence_scores)))
cat("  Displayed conv tiers:\n"); print(table(convergence_scores))
cat(sprintf("  Labeled genes shown: %d / %d\n", sum(genes %in% label_genes), length(label_genes)))

# ── Build matrices ───────────────────────────────────────────────────────────

zscore_col <- function(x) {
  x_num <- as.numeric(x)
  mu <- mean(x_num, na.rm = TRUE)
  s  <- sd(x_num, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x_num)
  (x_num - mu) / s
}
# Robust z-score using median + MAD (resistant to skewed distributions)
zscore_robust <- function(x) {
  x_num <- as.numeric(x)
  med <- median(x_num, na.rm = TRUE)
  mad_val <- mad(x_num, na.rm = TRUE)
  if (is.na(mad_val) || mad_val == 0) return(zscore_col(x_num))
  (x_num - med) / mad_val
}
clamp <- function(x, lo, hi) {
  out <- pmin(pmax(x, lo), hi)
  out[is.na(x)] <- NA
  out
}

# Column 1: Human logFC (z-score normalized) — kept for backward compat
mat_human_z <- matrix(zscore_col(as.numeric(top_genes$bulk_logFC)), ncol = 1,
                      dimnames = list(genes, "bulk_logFC"))

# Multi-contrast DEG matrices: z-scored over selected genes, NA for non-significant
# (deg_contrast_defs defined above, before convergence scoring)
mat_deg_z_sig <- do.call(cbind, lapply(deg_contrast_defs, function(d) {
  if (!d$lfc %in% names(top_genes)) {
    return(matrix(NA_real_, nrow = top_n, ncol = 1))
  }
  lfc  <- as.numeric(top_genes[[d$lfc]])
  padj <- if (d$padj %in% names(top_genes)) as.numeric(top_genes[[d$padj]]) else rep(1, top_n)
  sig  <- !is.na(padj) & padj < 0.05 & abs(lfc) > d$lfc_thresh
  z    <- zscore_col(lfc)
  z[!sig] <- NA_real_
  matrix(z, ncol = 1)
}))
colnames(mat_deg_z_sig) <- sapply(deg_contrast_defs, `[[`, "label")
rownames(mat_deg_z_sig) <- genes
n_deg_cols <- ncol(mat_deg_z_sig)  # = 6
cat(sprintf("DEG contrast columns built: %d\n", n_deg_cols))

# COLOC grouped by phenotype category — unified PP.H4 = max(ABF, SuSiE) per
# GWAS, then max across GWAS within the category. Built from per-GWAS
# susie_coloc CSVs via build_category_pp4() at the top of this script.
get_cat <- function(col) {
  if (col %in% names(top_genes)) as.numeric(top_genes[[col]]) else rep(NA_real_, top_n)
}
mat_enzyme_pp4  <- matrix(get_cat("category_pp4_enzyme"),  ncol = 1,
                           dimnames = list(genes, "Enzyme\nCOLOC"))
mat_disease_pp4 <- matrix(get_cat("category_pp4_disease"), ncol = 1,
                           dimnames = list(genes, "Disease\nCOLOC"))
mat_pdff        <- matrix(get_cat("category_pp4_pdff"),    ncol = 1,
                           dimnames = list(genes, "PDFF\nCOLOC"))

# Ancestry COUNT: # of ancestries (EUR / EAS / AFR / CSA) with any COLOC PP.H4 > 0.5.
# ARCHIVED 2026-04-08: FinnGen COLOC columns removed from EUR classification.
eur_coloc_cols <- intersect(c("broadaway_coloc_pp4", "ukbb_alt_coloc_pp4", "ast_coloc_pp4",
  "ggt_coloc_pp4", "pdff_coloc_pp4", "decode_nafl_coloc_pp4"), names(top_genes))
eas_coloc_cols <- intersect(c("bbj_alt_coloc_pp4", "bbj_ast_coloc_pp4", "bbj_ggt_coloc_pp4"), names(top_genes))
afr_coloc_cols <- intersect(c("panukbb_afr_alt_coloc_pp4", "panukbb_afr_ast_coloc_pp4", "panukbb_afr_ggt_coloc_pp4"), names(top_genes))
csa_coloc_cols <- intersect(c("panukbb_csa_alt_coloc_pp4", "panukbb_csa_ast_coloc_pp4", "panukbb_csa_ggt_coloc_pp4"), names(top_genes))
.has_anc <- function(cols) { if (!length(cols)) return(rep(FALSE, top_n))
  apply(as.matrix(top_genes[, ..cols]), 1, function(x) any(!is.na(x) & x > 0.5)) }
n_ancestry <- as.integer(.has_anc(eur_coloc_cols)) + as.integer(.has_anc(eas_coloc_cols)) +
              as.integer(.has_anc(afr_coloc_cols)) + as.integer(.has_anc(csa_coloc_cols))

# Column 7: Multi-ancestry gene-level finemapping max PIP
# pmax over MESuSiE (multi-ancestry shared signal) and SuSiEx (multi-ancestry joint).
finemap_pip <- rep(NA_real_, top_n)
if ("mesusie_max_pip" %in% names(top_genes)) finemap_pip <- as.numeric(top_genes$mesusie_max_pip)
if ("susiex_max_pip"  %in% names(top_genes)) {
  finemap_pip <- pmax(finemap_pip, as.numeric(top_genes$susiex_max_pip), na.rm = TRUE)
}
finemap_pip[is.infinite(finemap_pip)] <- NA
mat_finemap <- matrix(finemap_pip, ncol = 1, dimnames = list(genes, "Finemap_PIP"))

# Column 9: scEQTL COLOC (best PP.H4 across cell types)
mat_sceqtl <- matrix(
  if ("sceqtl_coloc_best_pp4" %in% names(top_genes)) as.numeric(top_genes$sceqtl_coloc_best_pp4) else NA,
  ncol = 1, dimnames = list(genes, "scEQTL_COLOC"))

# Column 10: ieQTL (disease-interaction eQTL, binary)
mat_ieqtl <- matrix(0L, nrow = top_n, ncol = 1, dimnames = list(genes, "ieQTL"))
if ("ieqtl_disease_interaction" %in% names(top_genes)) {
  mat_ieqtl[, 1] <- as.integer(top_genes$ieqtl_disease_interaction == TRUE & !is.na(top_genes$ieqtl_disease_interaction))
}

# Column 11: TWAS z
mat_twas <- matrix(
  if ("twas_z" %in% names(top_genes)) as.numeric(top_genes$twas_z) else NA,
  ncol = 1, dimnames = list(genes, "TWAS_z"))

# Column 10: Mouse logFC (z-score normalized, same scale as human)
mat_mouse_z <- matrix(zscore_col(
  if ("mouse_meta_logFC" %in% names(top_genes)) as.numeric(top_genes$mouse_meta_logFC) else rep(NA, top_n)),
  ncol = 1, dimnames = list(genes, "mouse_logFC"))

# Column 11: Progression score (composite: causal score where available, class-derived otherwise)
mat_prog <- matrix(
  if ("progression_composite" %in% names(top_genes)) as.numeric(top_genes$progression_composite) else
  if ("progression_causal_score_driver" %in% names(top_genes)) as.numeric(top_genes$progression_causal_score_driver) else NA,
  ncol = 1, dimnames = list(genes, "prog_score"))

# ATAC: DA logFC (diverging) + N motifs disrupted (count) + Disease regulon target (binary)
# DA logFC: chromatin accessibility change (sig only, non-sig = NA)
mat_da_lfc <- matrix(NA, nrow = top_n, ncol = 1, dimnames = list(genes, "DA\nlogFC"))
if ("best_da_logFC" %in% names(top_genes)) {
  sig_da <- !is.na(top_genes$best_da_padj) & top_genes$best_da_padj < 0.05
  mat_da_lfc[sig_da, 1] <- as.numeric(top_genes$best_da_logFC[sig_da])
}

# Disrupted master-regulator targets (count) — ONE dense regulatory column that
# replaces the old sparse per-gene "N motif disrupted" + empty "Disease regulon".
# = # of motif-disrupted disease master-regulators (CollecTRI) regulating the gene.
mat_mrtarget <- matrix(NA, nrow = top_n, ncol = 1, dimnames = list(genes, "Targets of\ndisrupted TFs"))
if ("n_disrupted_mr_targets" %in% names(top_genes)) {
  v <- as.numeric(top_genes$n_disrupted_mr_targets)
  mat_mrtarget[!is.na(v) & v > 0, 1] <- v[!is.na(v) & v > 0]   # 0 -> NA (white = no evidence)
}

# Proteomics: protein logFC (z-score normalized) + n_datasets validated
mat_prot <- matrix(NA, nrow = top_n, ncol = 2,
  dimnames = list(genes, c("Protein\nlogFC (z)", "N datasets\nvalidated")))
if ("best_protein_logFC" %in% names(top_genes)) {
  mat_prot[, 1] <- zscore_robust(as.numeric(top_genes$best_protein_logFC))
}
if ("n_prot_datasets" %in% names(top_genes)) {
  mat_prot[, 2] <- as.numeric(top_genes$n_prot_datasets)
}

# Spatial: Moran's I (max across datasets)
mat_spatial <- matrix(NA, nrow = top_n, ncol = 1,
  dimnames = list(genes, "Moran's I"))
if ("spatial_max_I" %in% names(top_genes)) {
  mat_spatial[, 1] <- as.numeric(top_genes$spatial_max_I)
}

# Clinical drug-development status (replaces binary Druggable; categories from fig4h)
# 0: None, 1: Approved, 2: In trials, 3: Failed / discontinued, 4: Preclinical, 5: Discovery
mat_drug <- matrix(0L, nrow = top_n, ncol = 1, dimnames = list(genes, "Clinical\nstatus"))
if ("drug_dev_status" %in% names(top_genes)) {
  s <- top_genes$drug_dev_status
  mat_drug[!is.na(s) & s == "masld_approved",     1] <- 1L
  mat_drug[!is.na(s) & s == "masld_clinical",      1] <- 2L
  mat_drug[!is.na(s) & s == "masld_discontinued",  1] <- 3L
  mat_drug[!is.na(s) & s == "masld_preclinical",   1] <- 4L
  mat_drug[!is.na(s) & s == "discovery",           1] <- 5L
}

# Pathway memberships count
mat_pathway <- matrix(NA, nrow = top_n, ncol = 1,
  dimnames = list(genes, "N_pathways"))
pw_col <- intersect(c("n_pathway_memberships", "n_leading_edge_pathways"), names(top_genes))
if (length(pw_col) > 0) {
  mat_pathway[, 1] <- as.numeric(top_genes[[pw_col[1]]])
}

# Sex dimorphism: based on sex INTERACTION test from Script 26 v2
# sex_class now uses interaction-based names (Female_biased, Male_biased, Concordant, Divergent)
# Only genes with significant sex interaction (padj < 0.05) are truly sex-dimorphic
# Sex x disease interaction: interaction-based classification (v2)
mat_sex <- matrix(NA_character_, nrow = top_n, ncol = 1, dimnames = list(genes, "Sex\nclass"))
if ("sex_class" %in% names(top_genes)) {
  sc <- as.character(top_genes$sex_class)
  mat_sex[sc == "Female_biased", 1] <- "F"
  mat_sex[sc == "Male_biased", 1]   <- "M"
  mat_sex[sc == "Divergent", 1]     <- "D"
  # Concordant = white (leave as NA)
}

# Clinical drug validation level (ordinal)
mat_clin <- matrix(NA_real_, nrow = top_n, ncol = 1, dimnames = list(genes, "Clinical\nevidence"))
if (!is.null(drug_val_best)) {
  idx <- match(genes, drug_val_best$target_gene)
  matched <- !is.na(idx)
  mat_clin[matched, 1] <- drug_val_best$clin_level[idx[matched]]
}

# ── Significance vectors for DEG stars ───────────────────────────────────────
# Human DEG significance: bulk_treat_fdr < 0.05 (canonical TREAT; lfc=0.25 floor in test)
human_sig <- !is.na(top_genes$bulk_treat_fdr) & top_genes$bulk_treat_fdr < 0.05
# Mouse DEG significance
mouse_sig <- if ("mouse_meta_padj" %in% names(top_genes)) {
  !is.na(top_genes$mouse_meta_padj) & top_genes$mouse_meta_padj < 0.05
} else { rep(FALSE, top_n) }

# ── Color scales (publication_color_themes.R gradients) ──────────────────────
source(file.path(BASE, "scripts/publication_color_themes.R"))
na_bg <- background1_gradient[1]  # light lavender-gray for NA cells

# Sunset: Genetic causal (PP.H4), Clinical evidence — core causal (KEEP AS-IS)
col_pp4     <- colorRamp2(c(0, 0.5, 1), c("white", sunset_gradient[4], sunset_gradient[11]))
col_clin    <- colorRamp2(c(0, 1, 2, 3), c("white", sunset_gradient[2], sunset_gradient[6], sunset_gradient[11]))

# Sea/Sunset diverging: z-scores — brighter endpoints
col_zscore  <- colorRamp2(c(-2.5, 0, 2.5), c(sea_gradient[9], "white", sunset_gradient[7]))
col_twas    <- col_zscore  # share the same scale (TWAS z clamped to same range)

# Progression score: green gradient (distinct from blue/purple/sunset)
col_prog    <- colorRamp2(c(0, 0.3, 1), c("white", green_gradient[4], green_gradient[10]))
# Unified count scale: purple gradient for N pathways, N motifs, N datasets
col_count <- colorRamp2(c(0, 5, 20), c("white", purple_gradient[3], purple_gradient[9]))
col_morans  <- colorRamp2(c(0, 0.05, 0.15), c("white", forest_gradient[3], forest_gradient[9]))

# ATAC counts in purple; Drug/EAS stay gray
col_binary  <- structure(c("white", gray_gradient[7]), names = c("0", "1"))
# Clinical drug-development status colors: categories and colors from fig4h calibration plot
clin_status_cols <- c("white", "#2E7D32", "#F9A825", "#C62828", "#5B9BD5", "#9E9E9E")
clin_status_labs <- c("Approved", "In trials", "Failed / discontinued", "Preclinical", "Discovery")
col_ancestry <- structure(c(blue_gradient[7], sunset_gradient[6], purple_gradient[8]), names = c("EUR", "EAS", "Both"))
col_sex     <- structure(c(sunset_gradient[7], sea_gradient[9], purple_gradient[8]), names = c("F", "M", "D"))

ht_opt(
  heatmap_row_names_gp    = gpar(fontsize = 7, fontface = "italic", fontfamily = "Helvetica"),
  heatmap_column_names_gp = gpar(fontsize = 7, fontfamily = "Helvetica"),
  heatmap_column_title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "Helvetica")
)

# mat_deg_z_sig already has NA for non-significant entries (built above)
# mat_human_z_sig: first column of multi-contrast matrix for backward compat
mat_human_z_sig <- mat_deg_z_sig[, 1, drop = FALSE]
mat_mouse_z_sig <- mat_mouse_z
mat_mouse_z_sig[!mouse_sig, 1] <- NA

# ── Row annotations: convergence bar (left) + optional Primary axis (left). ──────────
if (FIG5A_LAYOUT == "therapeutic_axes") {
  ha_left <- rowAnnotation(
    " " = therapeutic_axis_split,
    "Conv." = anno_barplot(convergence_scores,
      width = unit(13, "mm"), border = FALSE,
      gp = gpar(fill = gray_gradient[9], col = NA),
      axis_param = list(gp = gpar(fontsize = 6), direction = "reverse"),
      bar_width = 0.7, direction = "reverse"),
    col = list(" " = structure(axis_colors, names = axis_display)),
    simple_anno_size = unit(2.5, "mm"),
    show_legend = FALSE,
    annotation_name_gp = gpar(fontsize = 6.5, fontface = "bold"),
    annotation_name_rot = 0,
    annotation_name_side = "top",
    annotation_name_offset = unit(2, "mm"))
} else {
  ha_left <- rowAnnotation(
    "Conv." = anno_barplot(convergence_scores,
      width = unit(13, "mm"), border = FALSE,
      gp = gpar(fill = gray_gradient[9], col = NA),
      axis_param = list(gp = gpar(fontsize = 6), direction = "reverse"),
      bar_width = 0.7, direction = "reverse"),
    annotation_name_gp = gpar(fontsize = 6.5, fontface = "bold"),
    annotation_name_rot = 0,
    annotation_name_side = "top",
    annotation_name_offset = unit(2, "mm"))
}

# ── Build heatmaps — ALL COLUMNS FLAT (no group headers) ────────────────────

# Human logFC (z-score, * for significant DEGs)
# Consolidated 2026-06-22: show only 3 key contrasts (disease / steatohepatitis
# progression / advanced fibrosis) instead of all 6 — the rest were redundant.
mat_rna_combined <- clamp(mat_deg_z_sig, -2.5, 2.5)
rna_show <- c("MASLD\nvs Ctrl", "MASH\nvs MAFL", "Adv Fib\n(F3-F4)")
mat_rna_combined <- mat_rna_combined[, intersect(rna_show, colnames(mat_rna_combined)), drop = FALSE]

h_human <- Heatmap(mat_rna_combined, name = "z-score", col = col_zscore,
  column_labels = colnames(mat_rna_combined), column_title = "Transcriptomics",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  row_split = therapeutic_axis_split,
  cluster_row_slices = FALSE,
  row_gap = if (FIG5A_LAYOUT == "therapeutic_axes") unit(2.2, "mm") else unit(0, "mm"),
  row_title_gp = gpar(
    fontsize = 7, fontface = "bold", fontfamily = "Helvetica",
    col = "black"),
  row_title_side = "left",
  row_title_rot = 0,
  left_annotation = ha_left,
  na_col = "white",
  width = unit((if (!is.na(TOPN)) SQ_CELL_MM else 9) * ncol(mat_rna_combined), "mm"),
  height = if (!is.na(TOPN)) unit(top_n * SQ_CELL_MM, "mm") else NULL,
  border = TRUE,
  heatmap_legend_param = list(title = "z-score", title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5), legend_height = unit(2, "cm"), at = c(-2, 0, 2)))

# Consolidated 2026-06-22: 3 COLOC categories (enzyme/disease/PDFF) -> one best
# PP.H4 column. Genetics block = COLOC | Finemap PIP.
mat_coloc_best <- matrix(
  suppressWarnings(pmax(mat_enzyme_pp4[, 1], mat_disease_pp4[, 1], mat_pdff[, 1], na.rm = TRUE)),
  ncol = 1, dimnames = list(genes, "COLOC\n(PP.H4)"))
mat_coloc_best[is.infinite(mat_coloc_best)] <- NA
mat_genetic <- cbind(
  clamp(mat_coloc_best, 0, 1),
  clamp(mat_finemap, 0, 1)
)
genetic_labels <- c("COLOC\n(PP.H4)", "Finemap\nPIP")

h_genetic <- Heatmap(mat_genetic, name = "PP.H4", col = col_pp4,
  column_labels = genetic_labels, column_title = "Genetics",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  na_col = "white",
  width = unit((if (!is.na(TOPN)) SQ_CELL_MM else 10) * ncol(mat_genetic), "mm"),
  border = TRUE,
  heatmap_legend_param = list(title = "PP.H4", title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5), legend_height = unit(2, "cm"), at = c(0, 0.5, 1)))

# TWAS: -log10(p) * sign(z) — significance-weighted direction
mat_twas_signed <- matrix(NA, nrow = top_n, ncol = 1, dimnames = list(genes, "TWAS"))
if ("twas_z" %in% names(top_genes)) {
  tz <- as.numeric(top_genes$twas_z)
  twas_p <- 2 * pnorm(-abs(tz))  # two-sided p from z
  mat_twas_signed[, 1] <- -log10(pmax(twas_p, 1e-20)) * sign(tz)
}
col_twas_sig <- colorRamp2(c(-10, 0, 10), c(sea_gradient[9], "white", sunset_gradient[7]))
h_twas <- Heatmap(clamp(mat_twas_signed, -10, 10), name = "TWAS", col = col_twas_sig,
  column_labels = "TWAS\n-log10(p)", column_title = "TWAS",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  na_col = "white", width = unit(10, "mm"), border = TRUE,
  heatmap_legend_param = list(title = "TWAS", title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5), legend_height = unit(1.5, "cm"), at = c(-10, 0, 10)))

# Mouse logFC (z-score, * for significant)
h_mouse <- Heatmap(clamp(mat_mouse_z_sig, -2.5, 2.5), name = "Mouse z", col = col_zscore,
  column_labels = "Mouse\nlogFC (z)", column_title = "Mouse",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  na_col = "white", width = unit(if (!is.na(TOPN)) SQ_CELL_MM else 12, "mm"), border = TRUE,
  show_heatmap_legend = FALSE)

# ATAC/regulatory: ONE dense column — disrupted master-regulator targets (count).
mat_atac_combined <- mat_mrtarget

h_atac <- Heatmap(mat_atac_combined, name = "ATAC", col = col_count,
  column_labels = colnames(mat_atac_combined), column_title = "Regulatory",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  na_col = "white",
  width = unit(if (!is.na(TOPN)) SQ_CELL_MM * ncol(mat_atac_combined) else 14, "mm"),
  border = TRUE,
  show_heatmap_legend = FALSE)  # uses unified Count legend

# Proteomics validation (protein logFC z-score)
mat_prot_combined <- clamp(mat_prot[, 1, drop = FALSE], -2.5, 2.5)

h_prot <- Heatmap(mat_prot_combined,
  name = "Protein", col = col_zscore,
  column_labels = "Protein\nlogFC (z)",
  column_title = "Proteomics",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  na_col = "white",
  width = unit(if (!is.na(TOPN)) SQ_CELL_MM * ncol(mat_prot_combined) else 12, "mm"),
  border = TRUE,
  show_heatmap_legend = FALSE)  # shares z-score legend

# Spatial validation (Moran's I)
h_spatial <- Heatmap(mat_spatial, name = "Spatial", col = col_morans,
  column_labels = "Moran's I", column_title = "Spatial",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  na_col = "white", width = unit(if (!is.na(TOPN)) SQ_CELL_MM else 12, "mm"), border = TRUE,
  heatmap_legend_param = list(title = "Moran's I", title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5), legend_height = unit(1.5, "cm"), at = c(0, 0.05, 0.1, 0.15)))

# Clinical drug validation level
h_clin <- Heatmap(mat_clin, name = "ClinVal", col = col_clin,
  column_labels = "Clinical\nevidence", column_title = "Clinical",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  na_col = "white", width = unit(10, "mm"), border = TRUE,
  heatmap_legend_param = list(title = "Clinical", title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5), at = c(0, 2, 3),
    labels = c("Absent", "Moderate", "Strong"),
    legend_height = unit(1.5, "cm")))

# Clinical status (rightmost). With the thinned, taller layout every gene fits,
# so label ALL genes directly as right-side row names (no anno_mark leader lines).
# Categorical fill via cell_fun (status colors); custom legend in annotation_legend_list.
h_drug <- Heatmap(mat_drug, name = "Clinical", col = col_binary,
  column_labels = "Clinical\nstatus", column_title = "Therapeutics",
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, row_names_side = "right", row_labels = genes,
  row_names_gp = gpar(fontsize = 7, fontface = "italic", fontfamily = "Helvetica"),
  na_col = "white", width = unit(if (!is.na(TOPN)) SQ_CELL_MM else 12, "mm"), border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- mat_drug[i, j]; if (is.na(v)) v <- 0
    grid.rect(x, y, width, height, gp = gpar(fill = clin_status_cols[v + 1], col = NA))
  },
  show_heatmap_legend = FALSE)

# ── Assemble ─────────────────────────────────────────────────────────────────
cat(sprintf("Saving to %s ...\n", OUTPDF))
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf

if (!is.na(TOPN)) {
  # ── Transposed layout (top-N variant) ──────────────────────────────────────
  # Genes -> columns, modalities -> rows. Group names become horizontal row
  # titles on the left (room for "GWAS + eQTL", "Proteomics"), and the gene
  # names label the columns along the bottom.
  cw    <- unit(SQ_CELL_MM * top_n, "mm")             # shared gene-axis width
  rh    <- function(n) unit(SQ_CELL_MM * n, "mm")     # modality-block height
  rt_gp <- gpar(fontsize = 8, fontface = "bold", fontfamily = "Helvetica")

  # Convergence barplot now sits on top (per-gene = per-column)
  ha_top <- HeatmapAnnotation(
    "Conv." = anno_barplot(convergence_scores, height = unit(10, "mm"),
      border = FALSE, gp = gpar(fill = gray_gradient[9], col = NA),
      axis_param = list(gp = gpar(fontsize = 6)), bar_width = 0.7),
    annotation_name_gp = gpar(fontsize = 6.5, fontface = "bold"),
    annotation_name_side = "left")

  # ── Per-modality row cap (transposed top-N variant): max 3 rows/modality ─────
  # RNA and Genetics already consolidated to key cols; match those labels.
  rna_keep_labels <- c("MASLD\nvs Ctrl", "MASH\nvs MAFL", "Adv Fib\n(F3-F4)")
  gen_keep_labels <- c("COLOC\n(PP.H4)", "Finemap\nPIP")
  rna_keep_idx <- match(rna_keep_labels, colnames(mat_rna_combined)); rna_keep_idx <- rna_keep_idx[!is.na(rna_keep_idx)]
  gen_keep_idx <- match(gen_keep_labels, genetic_labels);             gen_keep_idx <- gen_keep_idx[!is.na(gen_keep_idx)]
  mat_rna_sel    <- mat_rna_combined[, rna_keep_idx, drop = FALSE]
  mat_gen_sel    <- mat_genetic[,      gen_keep_idx, drop = FALSE]
  gen_labels_sel <- genetic_labels[gen_keep_idx]

  t_rna <- t(mat_rna_sel)
  th_rna <- Heatmap(t_rna, name = "z-score", col = col_zscore,
    row_title = "Transcriptomics", top_annotation = ha_top, height = rh(nrow(t_rna)),
    cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6.5, fontfamily = "Helvetica"),
    show_column_names = FALSE, row_title_rot = 0, row_title_gp = rt_gp,
    na_col = "white", border = TRUE, width = cw,
    heatmap_legend_param = list(title = "z-score", title_gp = gpar(fontsize = 7, fontface = "bold"),
      labels_gp = gpar(fontsize = 6.5), legend_height = unit(2, "cm"), at = c(-2, 0, 2)))

  t_gen <- t(mat_gen_sel)
  th_gen <- Heatmap(t_gen, name = "PP.H4", col = col_pp4,
    row_title = "GWAS + eQTL", row_labels = gen_labels_sel, height = rh(nrow(t_gen)),
    cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6.5, fontfamily = "Helvetica"),
    show_column_names = FALSE, row_title_rot = 0, row_title_gp = rt_gp,
    na_col = "white", border = TRUE, width = cw,
    heatmap_legend_param = list(title = "PP.H4", title_gp = gpar(fontsize = 7, fontface = "bold"),
      labels_gp = gpar(fontsize = 6.5), legend_height = unit(2, "cm"), at = c(0, 0.5, 1)))

  t_prot <- t(mat_prot_combined)
  th_prot <- Heatmap(t_prot, name = "Protein", col = col_zscore,
    row_title = "Proteomics", row_labels = "Protein\nlogFC (z)", height = rh(nrow(t_prot)),
    cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6.5, fontfamily = "Helvetica"),
    show_column_names = FALSE, row_title_rot = 0, row_title_gp = rt_gp,
    na_col = "white", border = TRUE, width = cw, show_heatmap_legend = FALSE)

  t_atac <- t(mat_atac_combined)
  th_atac <- Heatmap(t_atac, name = "ATAC", col = col_count,
    row_title = "Regulatory", height = rh(nrow(t_atac)),
    cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6.5, fontfamily = "Helvetica"),
    show_column_names = FALSE, row_title_rot = 0, row_title_gp = rt_gp,
    na_col = "white", border = TRUE, width = cw,
    cell_fun = function(j, i, x, y, width, height, fill) {
      v <- t_atac[i, j]
      if (is.na(v)) return()
      if (i == 2) {  # disease regulon row: binary gray
        col_val <- if (v == 1) gray_gradient[7] else "white"
        grid.rect(x, y, width, height, gp = gpar(fill = col_val, col = NA))
      }
    },
    show_heatmap_legend = FALSE)

  t_spa <- t(mat_spatial)
  th_spa <- Heatmap(t_spa, name = "Moran's I", col = col_morans,
    row_title = "Spatial", row_labels = "Moran's I", height = rh(nrow(t_spa)),
    cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6.5, fontfamily = "Helvetica"),
    show_column_names = FALSE, row_title_rot = 0, row_title_gp = rt_gp,
    na_col = "white", border = TRUE, width = cw,
    heatmap_legend_param = list(title = "Moran's I", title_gp = gpar(fontsize = 7, fontface = "bold"),
      labels_gp = gpar(fontsize = 6.5), legend_height = unit(1.5, "cm"), at = c(0, 0.05, 0.1, 0.15)))

  t_mou <- t(mat_mouse_z_sig)
  th_mou <- Heatmap(t_mou, name = "Mouse z", col = col_zscore,
    row_title = "Mouse", row_labels = "Mouse\nlogFC (z)", height = rh(nrow(t_mou)),
    cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6.5, fontfamily = "Helvetica"),
    show_column_names = FALSE, row_title_rot = 0, row_title_gp = rt_gp,
    na_col = "white", border = TRUE, width = cw, show_heatmap_legend = FALSE)

  t_dru <- t(mat_drug)
  th_dru <- Heatmap(t_dru, name = "Clinical", col = col_binary,
    row_title = "Therapeutics", row_labels = "Clinical\nstatus", height = rh(nrow(t_dru)),
    cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left",
    row_names_gp = gpar(fontsize = 6.5, fontfamily = "Helvetica"),
    show_column_names = TRUE, column_labels = genes, column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 7, fontface = "italic", fontfamily = "Helvetica"),
    column_names_rot = 90, row_title_rot = 0, row_title_gp = rt_gp,
    na_col = "white", border = TRUE, width = cw,
    cell_fun = function(j, i, x, y, width, height, fill) {
      v <- t_dru[i, j]; if (is.na(v)) v <- 0
      grid.rect(x, y, width, height, gp = gpar(fill = clin_status_cols[v + 1], col = NA))
    },
    show_heatmap_legend = FALSE)

  ht_list <- th_rna %v% th_gen %v% th_prot %v% th_atac %v% th_spa %v% th_dru
  total_rows <- nrow(t_rna) + nrow(t_gen) + nrow(t_prot) + nrow(t_atac) +
                nrow(t_spa) + nrow(t_dru)
  dev_w <- SQ_CELL_MM * top_n / 25.4 + 5.5
  dev_h <- total_rows * SQ_CELL_MM / 25.4 + 3
} else {
  # Order: RNA | Genetic Causal | Proteomics | ATAC | Spatial | Drug
  # (ATAC moved after Proteomics 2026-06-17)
  ht_list <- h_human + h_genetic + h_prot + h_atac + h_spatial + h_drug
  dev_w <- if (FIG5A_LAYOUT == "therapeutic_axes") 16.9 else 15.0
  # Taller gene boxes (2026-06-22): fixed per-gene row height + fixed margin so
  # each box is the same height regardless of how many genes survive the filter.
  dev_h <- top_n * 0.165 + if (FIG5A_LAYOUT == "therapeutic_axes") 1.7 else 1.4
}

pdf_device(OUTPDF, width = dev_w, height = dev_h)

# Progression score legend (green)
score_legend <- Legend(
  col_fun = col_prog,
  title = "Score",
  title_gp = gpar(fontsize = 7, fontface = "bold"),
  labels_gp = gpar(fontsize = 6.5),
  legend_height = unit(1.5, "cm"),
  at = c(0, 0.5, 1))

# Unified Count legend (N pathways, N motifs disrupted, N datasets — all use purple scale)
count_legend <- Legend(
  col_fun = col_count,
  title = "Count",
  title_gp = gpar(fontsize = 7, fontface = "bold"),
  labels_gp = gpar(fontsize = 6.5),
  legend_height = unit(1.5, "cm"),
  at = c(0, 5, 10))

ancestry_legend <- Legend(
  labels = c("1", "2", "3", "4"),
  legend_gp = gpar(fill = purple_gradient[c(3, 5, 7, 9)], col = NA),
  title = "COLOC ancestry (n)",
  title_gp = gpar(fontsize = 7, fontface = "bold"),
  labels_gp = gpar(fontsize = 6.5),
  grid_height = unit(3, "mm"),
  grid_width = unit(3, "mm"))

clinical_legend <- Legend(
  labels = clin_status_labs,
  legend_gp = gpar(fill = clin_status_cols[2:6], col = NA),
  title = "Clinical status",
  title_gp = gpar(fontsize = 7, fontface = "bold"),
  labels_gp = gpar(fontsize = 6.5),
  grid_height = unit(3, "mm"),
  grid_width = unit(3, "mm"))

draw(ht_list,
  heatmap_legend_side = "right",
  merge_legends = TRUE,
  annotation_legend_list = list(count_legend, clinical_legend),
  padding = unit(c(3, 3, 5, 3), "mm"))

dev.off()
ht_opt(RESET = TRUE)

n_total_cols <- ncol(mat_genetic) + 1 + 1 + 1 + 1 + 1 + 2 + 2 + 1 + 1 + 1 + 1 + 1 + 1  # genetic + EAS + prog + pathway + atac(2) + prot(2) + spatial + mouse + CC + sex + clin + drug
cat("Done.\n")
cat(sprintf("  %d genes, %d evidence columns\n", top_n, n_total_cols))
cat(sprintf("  Notable genes included: %s\n",
  paste(intersect(notable, genes), collapse = ", ")))
