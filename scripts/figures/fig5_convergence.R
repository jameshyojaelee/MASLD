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
OUTPDF <- file.path(PANDIR, "fig5a.pdf")

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

# Load INTACT scores — prefer the newer intact_scores.csv (2026-05-06+),
# fall back to convergence_summary.csv for older runs.
intact_file <- file.path(BASE, "RNA-seq/results/gwas_rna_integration/intact_scores.csv")
conv_file   <- file.path(BASE, "RNA-seq/results/gwas_rna_integration/convergence_summary.csv")
if (file.exists(intact_file)) {
  intact_dt <- fread(intact_file)
  if ("intact_score_bulk" %in% names(intact_dt)) {
    atlas <- merge(atlas, intact_dt[, .(gene, intact_score_bulk)],
                   by.x = "human_symbol", by.y = "gene", all.x = TRUE)
    cat("  INTACT scores merged from intact_scores.csv\n")
  }
} else if (file.exists(conv_file)) {
  conv <- fread(conv_file)
  if ("intact_score_bulk" %in% names(conv)) {
    atlas <- merge(atlas, conv[, .(gene, intact_score_bulk)],
                   by.x = "human_symbol", by.y = "gene", all.x = TRUE)
    cat("  INTACT scores merged from convergence_summary.csv (legacy)\n")
  }
}

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
# Rank by composite: max COLOC PP.H4 + INTACT + dream significance
# ARCHIVED 2026-04-08: FinnGen COLOC columns removed
coloc_pp4_cols <- intersect(c("broadaway_coloc_pp4", "ukbb_alt_coloc_pp4",
  "ast_coloc_pp4", "ggt_coloc_pp4", "pdff_coloc_pp4",
  "bbj_alt_coloc_pp4"), names(atlas))

atlas[, max_coloc_pp4 := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = coloc_pp4_cols]

# Composite ranking score
# Rank using available columns (some may not exist after atlas rebuild)
# Build score component by component
atlas[, rs_coloc := fifelse(is.na(max_coloc_pp4), 0, max_coloc_pp4) * 3]
atlas[, rs_intact := fifelse(is.na(intact_score_bulk), 0, intact_score_bulk) * 2]
atlas[, rs_deg := fifelse(is.na(dream_padj) | dream_padj >= 0.05 | abs(dream_logFC) <= 0.5, 0, -log10(pmax(dream_padj, 1e-300)) / 10)]
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
atlas[, rank_score := rs_coloc + rs_intact + rs_deg + rs_susie + rs_prog]
atlas[, c("rs_coloc", "rs_intact", "rs_deg", "rs_susie", "rs_prog") := NULL]

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
             "FGFR1", "GCGR", "GLP1R", "KLB", "PPARD", "PPARG", "SCD")

# Determine DEG filter column
deg_filter <- !is.na(atlas$dream_padj) & atlas$dream_padj < 0.05 & abs(atlas$dream_logFC) > 0.5

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
  list(lfc = "dream_logFC",        padj = "dream_padj",        lfc_thresh = 0.5, label = "MASLD\nvs Ctrl"),
  list(lfc = "nafl_vs_ctrl_logFC", padj = "nafl_vs_ctrl_padj", lfc_thresh = 0.0, label = "MAFL\nvs Ctrl"),
  list(lfc = "nash_vs_ctrl_logFC", padj = "nash_vs_ctrl_padj", lfc_thresh = 0.0, label = "MASH\nvs Ctrl"),
  list(lfc = "nafl_vs_nash_logFC", padj = "nafl_vs_nash_padj", lfc_thresh = 0.0, label = "MASH\nvs MAFL"),  # coeff=nafl_nashNASH: positive = up in MASH
  list(lfc = "adv_fib_logFC",      padj = "adv_fib_padj",      lfc_thresh = 0.0, label = "Adv Fib\n(F3-F4)"),
  list(lfc = "cirrhosis_logFC",    padj = "cirrhosis_padj",    lfc_thresh = 0.0, label = "Cirrhosis\n(F4)")
)

# ── Compute convergence score (N independent modalities with evidence) ──────
# Max = 6. Each modality counted once (0 or 1).
# Progression is excluded — it is derived from RNA, not an independent source.
top_genes[, n_convergence := 0L]

# 1. Human DEG: significant in ≥1 of the 6 RNA contrasts
deg_any_sig <- Reduce(`|`, lapply(deg_contrast_defs, function(d) {
  if (!d$lfc %in% names(top_genes)) return(rep(FALSE, nrow(top_genes)))
  lfc  <- as.numeric(top_genes[[d$lfc]])
  padj <- if (d$padj %in% names(top_genes)) as.numeric(top_genes[[d$padj]]) else rep(1, top_n)
  !is.na(padj) & padj < 0.05 & abs(lfc) > d$lfc_thresh
}))
top_genes[, n_convergence := n_convergence + as.integer(deg_any_sig)]

# 2. Mouse DEG
if ("mouse_meta_padj" %in% names(top_genes)) {
  top_genes[, n_convergence := n_convergence + fifelse(!is.na(mouse_meta_padj) & mouse_meta_padj < 0.05, 1L, 0L)]
}
# 3. Genetic causal (any COLOC PP.H4 > 0.8)
top_genes[, n_convergence := n_convergence + fifelse(!is.na(max_coloc_pp4) & max_coloc_pp4 > 0.8, 1L, 0L)]
# 4. Regulatory (ATAC motif disruption)
if ("n_motifs_disrupted" %in% names(top_genes)) {
  top_genes[, n_convergence := n_convergence + fifelse(!is.na(n_motifs_disrupted) & n_motifs_disrupted > 0, 1L, 0L)]
}
# 5. Proteomics validation
if ("best_protein_logFC" %in% names(top_genes)) {
  top_genes[, n_convergence := n_convergence + fifelse(!is.na(best_protein_logFC), 1L, 0L)]
}
# 6. Spatial
if ("spatial_max_I" %in% names(top_genes)) {
  top_genes[, n_convergence := n_convergence + fifelse(!is.na(spatial_max_I) & spatial_max_I > 0.05, 1L, 0L)]
}

# Sort by convergence (descending), then by 46d convergence_score, falling back
# to the heuristic rank_score for genes without a 46d posterior.
top_genes <- top_genes[order(-n_convergence, -bayev_score, -rank_score,
                              na.last = TRUE)]

# Keep only genes with ≥3 independent modalities converging
top_genes <- top_genes[n_convergence >= 3]

top_n <- nrow(top_genes)
genes <- top_genes$human_symbol
convergence_scores <- top_genes$n_convergence
cat(sprintf("Selected %d genes (convergence >= 3; %d notable forced in)\n", top_n, sum(notable_rows$human_symbol %in% genes)))
cat(sprintf("  Convergence range: %d - %d\n", min(convergence_scores), max(convergence_scores)))

# Key genes to label (discussed in paper)
label_genes <- c("THRB", "HKDC1", "RORA", "DGAT2", "PPARA", "NR1H4", "PNPLA3",
                  "HSD17B13", "GCKR",
                  "MMP9", "FGFR1", "GLP1R", "PPARG", "KLB",
                  "SCD", "LPL", "EPHA2", "AKR1B10", "FAP",
                  "GAS6", "TM4SF4", "CDK6")

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
mat_human_z <- matrix(zscore_col(as.numeric(top_genes$dream_logFC)), ncol = 1,
                      dimnames = list(genes, "dream_logFC"))

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

# Ancestry classification: EUR only / EAS only / Both / None
# ARCHIVED 2026-04-08: FinnGen COLOC columns removed from EUR classification
eur_coloc_cols <- intersect(c("broadaway_coloc_pp4", "ukbb_alt_coloc_pp4", "ast_coloc_pp4",
  "ggt_coloc_pp4", "pdff_coloc_pp4",
  "ghouse_cirrhosis_coloc_pp4", "ghouse_hcc_coloc_pp4",
  "decode_nafl_coloc_pp4", "decode_cirrhosis_coloc_pp4", "decode_hcc_coloc_pp4"), names(top_genes))
eas_coloc_cols <- intersect(c("bbj_alt_coloc_pp4", "bbj_ast_coloc_pp4", "bbj_ggt_coloc_pp4"), names(top_genes))
mat_ancestry <- matrix(NA_character_, nrow = top_n, ncol = 1, dimnames = list(genes, "Ancestry"))
if (length(eur_coloc_cols) > 0 && length(eas_coloc_cols) > 0) {
  eur_mat <- as.matrix(top_genes[, ..eur_coloc_cols])
  eas_mat <- as.matrix(top_genes[, ..eas_coloc_cols])
  has_eur <- apply(eur_mat, 1, function(x) any(!is.na(x) & x > 0.5))
  has_eas <- apply(eas_mat, 1, function(x) any(!is.na(x) & x > 0.5))
  mat_ancestry[has_eur & has_eas, 1]  <- "Both"
  mat_ancestry[has_eur & !has_eas, 1] <- "EUR"
  mat_ancestry[!has_eur & has_eas, 1] <- "EAS"
}

# Column 7: Multi-ancestry gene-level finemapping max PIP
# pmax over MESuSiE (multi-ancestry shared signal) and SuSiEx (multi-ancestry joint).
finemap_pip <- rep(NA_real_, top_n)
if ("mesusie_max_pip" %in% names(top_genes)) finemap_pip <- as.numeric(top_genes$mesusie_max_pip)
if ("susiex_max_pip"  %in% names(top_genes)) {
  finemap_pip <- pmax(finemap_pip, as.numeric(top_genes$susiex_max_pip), na.rm = TRUE)
}
finemap_pip[is.infinite(finemap_pip)] <- NA
mat_finemap <- matrix(finemap_pip, ncol = 1, dimnames = list(genes, "Finemap_PIP"))

# Column 8: Multi-INTACT
mat_intact <- matrix(
  if ("intact_score_bulk" %in% names(top_genes)) as.numeric(top_genes$intact_score_bulk) else NA,
  ncol = 1, dimnames = list(genes, "Multi_INTACT"))

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

# N motifs disrupted (GWAS-ATAC, count)
mat_motif <- matrix(NA, nrow = top_n, ncol = 1, dimnames = list(genes, "N motif\ndisrupted"))
if ("n_motifs_disrupted" %in% names(top_genes)) {
  mat_motif[, 1] <- as.numeric(top_genes$n_motifs_disrupted)
}

# Disease regulon target (binary)
mat_reg <- matrix(0L, nrow = top_n, ncol = 1, dimnames = list(genes, "Disease\nregulon"))
if ("is_disease_reg_target" %in% names(top_genes)) {
  mat_reg[, 1] <- as.integer(top_genes$is_disease_reg_target == TRUE)
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

# Druggable (binary)
# Druggable: combine DGIdb + OpenTargets + clinical drug validation targets
mat_drug <- matrix(0L, nrow = top_n, ncol = 1, dimnames = list(genes, "Druggable"))
if ("dgidb_druggable" %in% names(top_genes)) {
  mat_drug[, 1] <- fifelse(top_genes$dgidb_druggable == TRUE & !is.na(top_genes$dgidb_druggable), 1L, 0L)
}
if ("opentargets_drug" %in% names(top_genes)) {
  mat_drug[, 1] <- pmax(mat_drug[, 1], fifelse(top_genes$opentargets_drug == TRUE & !is.na(top_genes$opentargets_drug), 1L, 0L))
}
# Also mark clinical drug validation targets as druggable
if (!is.null(drug_val_best)) {
  clin_idx <- genes %in% drug_val_best$target_gene
  mat_drug[clin_idx, 1] <- 1L
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
# Human DEG significance: padj < 0.05 + |logFC| > 0.5
human_sig <- !is.na(top_genes$dream_padj) & top_genes$dream_padj < 0.05 &
             abs(top_genes$dream_logFC) > 0.5
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

# ── Row annotations: convergence bar (left) + gene labels (right). ──────────
ha_left <- rowAnnotation(
  "Conv." = anno_barplot(convergence_scores,
    width = unit(13, "mm"), border = FALSE,
    gp = gpar(fill = gray_gradient[9], col = NA),
    axis_param = list(gp = gpar(fontsize = 6), direction = "reverse"),
    bar_width = 0.7, direction = "reverse"),
  annotation_name_gp = gpar(fontsize = 6.5, fontface = "bold"),
  annotation_name_rot = 90,
  annotation_name_side = "top",
  annotation_name_offset = unit(2, "mm"))

# Right: selective gene labels using anno_mark
label_idx <- which(genes %in% label_genes)
ha_right <- rowAnnotation(
  gene = anno_mark(at = label_idx, labels = genes[label_idx],
    labels_gp = gpar(fontsize = 7.5, fontface = "italic", fontfamily = "Helvetica"),
    link_width = unit(8, "mm"),
    padding = unit(1, "mm"),
    link_gp = gpar(lwd = 0.5, col = gray_gradient[7])))

# ── Build heatmaps — ALL COLUMNS FLAT (no group headers) ────────────────────

# Human logFC (z-score, * for significant DEGs)
# Combined RNA matrix: 6 DEG contrast cols
mat_rna_combined <- clamp(mat_deg_z_sig, -2.5, 2.5)
colnames(mat_rna_combined) <- sapply(deg_contrast_defs, `[[`, "label")

h_human <- Heatmap(mat_rna_combined, name = "z-score", col = col_zscore,
  column_labels = colnames(mat_rna_combined), column_title = "RNA",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  left_annotation = ha_left,
  na_col = "white", width = unit(7 * ncol(mat_rna_combined), "mm"), border = TRUE,
  heatmap_legend_param = list(title = "z-score", title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5), legend_height = unit(2, "cm"), at = c(-2, 0, 2)))

# Encode ancestry as integer so it can join the numeric mat_genetic matrix.
# 0 = none, 1 = EUR, 2 = EAS, 3 = Both — rendered via cell_fun override.
mat_ancestry_num <- matrix(0L, nrow = top_n, ncol = 1,
                           dimnames = list(genes, "COLOC\nancestry"))
if (!all(is.na(mat_ancestry[, 1]))) {
  mat_ancestry_num[mat_ancestry[, 1] == "EUR",  1] <- 1L
  mat_ancestry_num[mat_ancestry[, 1] == "EAS",  1] <- 2L
  mat_ancestry_num[mat_ancestry[, 1] == "Both", 1] <- 3L
}

# Combine genetic columns: Enzyme COLOC | Disease COLOC | PDFF | SuSiE | INTACT | ieQTL | Ancestry
mat_genetic <- cbind(
  clamp(mat_enzyme_pp4, 0, 1),
  clamp(mat_disease_pp4, 0, 1),
  clamp(mat_pdff, 0, 1),
  clamp(mat_finemap, 0, 1),
  mat_ancestry_num
)
genetic_labels <- c("ALT/AST/GGT\nCOLOC", "NAFLD/NASH\nCOLOC", "PDFF\nCOLOC",
                     "Finemap\nPIP", "COLOC\nancestry")

h_genetic <- Heatmap(mat_genetic, name = "PP.H4", col = col_pp4,
  column_labels = genetic_labels, column_title = "GWAS + eQTL",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  na_col = "white", width = unit(8 * ncol(mat_genetic), "mm"), border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- mat_genetic[i, j]
    if (is.na(v)) return()
    if (j == ncol(mat_genetic)) {
      # Ancestry (last): categorical — overrides pp4 color scale
      col_val <- if (v == 1) col_ancestry["EUR"] else if (v == 2) col_ancestry["EAS"] else if (v == 3) col_ancestry["Both"] else "white"
      grid.rect(x, y, width, height, gp = gpar(fill = col_val, col = NA))
    }
  },
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
  na_col = "white", width = unit(10, "mm"), border = TRUE,
  show_heatmap_legend = FALSE)

# GWAS-ATAC: motif disruption counts + regulon hits
# Normalize motif count to max for color scale
# Combined ATAC matrix: N motifs disrupted | Disease regulon — single box
mat_atac_combined <- cbind(mat_motif, mat_reg)
colnames(mat_atac_combined) <- c("N motif\ndisrupted", "Disease\nregulon")

h_atac <- Heatmap(mat_atac_combined, name = "ATAC", col = col_count,
  column_labels = colnames(mat_atac_combined), column_title = "ATAC",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  na_col = "white", width = unit(18, "mm"), border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- mat_atac_combined[i, j]
    if (is.na(v)) return()
    if (j == 1) {
      # N motifs: purple count scale (default)
    } else {
      # Disease regulon: binary gray
      col_val <- if (v == 1) gray_gradient[7] else "white"
      grid.rect(x, y, width, height, gp = gpar(fill = col_val, col = NA))
    }
  },
  show_heatmap_legend = FALSE)  # uses unified Count legend

# Proteomics validation (protein logFC z-score)
mat_prot_combined <- clamp(mat_prot[, 1, drop = FALSE], -2.5, 2.5)

h_prot <- Heatmap(mat_prot_combined,
  name = "Protein", col = col_zscore,
  column_labels = "Protein\nlogFC (z)",
  column_title = "Proteomics",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  na_col = "white", width = unit(10, "mm"), border = TRUE,
  show_heatmap_legend = FALSE)  # shares z-score legend

# Spatial validation (Moran's I)
h_spatial <- Heatmap(mat_spatial, name = "Spatial", col = col_morans,
  column_labels = "Moran's I", column_title = "Spatial",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  na_col = "white", width = unit(10, "mm"), border = TRUE,
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

# Druggable (binary, rightmost — gene labels via anno_mark instead)
h_drug <- Heatmap(mat_drug, name = "Binary", col = col_binary,
  column_labels = "Druggable", column_title = "Drug",
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE,
  na_col = "white", width = unit(10, "mm"), border = TRUE,
  right_annotation = ha_right,
  heatmap_legend_param = list(title = "Binary", title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5), at = c(0, 1), labels = c("No", "Yes"),
    legend_height = unit(1, "cm")))

# ── Assemble ─────────────────────────────────────────────────────────────────
# Order: RNA | Progression | Pathways | ATAC | Genetic Causal | Ancestry | TWAS | Proteomics | Spatial | Mouse RNA | Sex | Clinical | Drug
ht_list <- h_human + h_atac + h_genetic + h_prot + h_spatial + h_mouse + h_drug

cat(sprintf("Saving to %s ...\n", OUTPDF))
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
pdf_device(OUTPDF, width = 14, height = max(5, top_n * 0.038))

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
  labels = c("EUR only", "EAS only", "Both"),
  legend_gp = gpar(fill = c(col_ancestry["EUR"], col_ancestry["EAS"], col_ancestry["Both"]), col = NA),
  title = "COLOC ancestry",
  title_gp = gpar(fontsize = 7, fontface = "bold"),
  labels_gp = gpar(fontsize = 6.5),
  grid_height = unit(3, "mm"),
  grid_width = unit(3, "mm"))

draw(ht_list,
  heatmap_legend_side = "right",
  annotation_legend_list = list(score_legend, count_legend, ancestry_legend),
  column_title = "Multi-Modal Convergence Matrix",
  column_title_gp = gpar(fontsize = 11, fontface = "bold", fontfamily = "Helvetica"),
  padding = unit(c(3, 3, 5, 3), "mm"))

dev.off()
ht_opt(RESET = TRUE)

n_total_cols <- ncol(mat_genetic) + 1 + 1 + 1 + 1 + 1 + 2 + 2 + 1 + 1 + 1 + 1 + 1 + 1  # genetic + EAS + prog + pathway + atac(2) + prot(2) + spatial + mouse + CC + sex + clin + drug
cat("Done.\n")
cat(sprintf("  %d genes, %d evidence columns\n", top_n, n_total_cols))
cat(sprintf("  Notable genes included: %s\n",
  paste(intersect(notable, genes), collapse = ", ")))
