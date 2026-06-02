#!/usr/bin/env Rscript
# 30_gwas_overlay_v2.R
# ---------------------------------------------------------------------------
# Strategy 12 (v2): GWAS Gene Overlay with Human-Only Mapping + Continuous
# Enrichment.  Fixes three problems in the original 21_gwas_spatial_overlay.R:
#   1. ~50% mapping failure (atlas requires both species) → added human-only tier
#   2. Pipe-delimited GWAS entries never matched       → split on "|"
#   3. Fisher exact test was under-powered              → Wilcoxon + permutation
# Follows Sanjana lab publication theme: PDF, Helvetica, lab palette.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggnewscale)
  library(scales)
  library(biomaRt)
  library(patchwork)
})

pdf.options(useDingbats = FALSE)
set.seed(42)

cat("=== Strategy 12 v2: GWAS Overlay (Human-Only + Continuous Enrichment) ===\n\n")

# ============================================================
#  Paths
# ============================================================
PROJECT  <- Sys.getenv("MASLD_PROJECT_ROOT",
              "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RNA_DIR  <- file.path(PROJECT, "RNA-seq")
OUT_DIR  <- file.path(RNA_DIR, "results", "gwas_spatial_convergence")
FIG_DIR  <- file.path(PROJECT, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_FILE  <- file.path(PROJECT, "GWAS", "Closest_genes.csv")
# Use ashr-shrunk dream results if available (lfsr-based significance)
DREAM_FILE <- file.path(RNA_DIR, "Human", "Patient_Cohorts", "analysis",
                        "integration", "results", "integration",
                        "dream_results_ashr.csv")
if (!file.exists(DREAM_FILE)) {
  DREAM_FILE <- file.path(RNA_DIR, "Human", "Patient_Cohorts", "analysis",
                          "integration", "results", "integration",
                          "dream_results.csv")
  cat("WARNING: dream_results_ashr.csv not found; falling back to dream_results.csv\n")
}
ATLAS_FILE <- file.path(RNA_DIR, "Analysis", "Cross_Species_Concordance",
                        "results", "concordance_atlas_unified.csv")
# MR_FILE removed 2026-04-22 — MR ditched from paper.

# ============================================================
#  Theme & palette (Sanjana lab standard)
# ============================================================
theme_pub <- theme_minimal(base_family = "Helvetica", base_size = 7) +
  theme(
    plot.title       = element_text(size = 8, face = "bold"),
    plot.subtitle    = element_text(size = 7, color = "grey40"),
    axis.title       = element_text(size = 8),
    axis.text        = element_text(size = 6),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 7),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(size = 7, face = "bold")
  )

pal <- list(
  magenta = "#e14b9d", pink = "#e35070", purple = "#d358c7",
  blue = "#4baeef", orange = "#e1b172", green = "#30d796",
  teal = "#2bbfbd", grey = "#808080"
)

concordance_colors <- c(
  Conserved       = pal$magenta,
  Moderate_Concordance = pal$pink,
  Human_Enriched       = pal$purple,
  Mouse_Specific       = pal$blue,
  Diet_Selective       = pal$orange,
  Species_Discordant   = pal$teal,
  Not_Significant      = "#cccccc",
  Unclassified         = "#888888"
)

# ============================================================
#  1. Clean GWAS input: split pipe-delimited entries, deduplicate
# ============================================================
cat("[1] Loading and cleaning GWAS genes...\n")
gwas_raw <- fread(GWAS_FILE, header = FALSE, col.names = "raw")
gwas_raw <- gwas_raw[raw != "" & !is.na(raw)]
gwas_raw[, raw := trimws(raw)]

# Split pipe-delimited entries (e.g. "PNPLA3|SAMM50")
gwas_split <- unlist(strsplit(gwas_raw$raw, "\\|"))
gwas_split <- trimws(gwas_split)

# Remove Ensembl IDs that crept in (e.g. ENSG00000258674) - keep only HGNC symbols
# HGNC symbols: alphabetical start, alphanumeric characters, possible hyphens
is_ensembl <- grepl("^ENSG[0-9]", gwas_split)
cat(sprintf("    Dropping %d Ensembl IDs from pipe-split entries\n", sum(is_ensembl)))
gwas_symbols <- unique(gwas_split[!is_ensembl])
gwas_symbols <- gwas_symbols[gwas_symbols != ""]

cat(sprintf("    Raw entries: %d | After pipe-split & dedup: %d unique HGNC symbols\n",
            nrow(gwas_raw), length(gwas_symbols)))

# ============================================================
#  2. Human-only mapping (Tier 1): biomaRt Ensembl → HGNC, then
#     join to dream_results.csv
# ============================================================
cat("\n[2] Human-only mapping via dream results...\n")
dream <- fread(DREAM_FILE)
cat(sprintf("    Dream results: %d genes\n", nrow(dream)))

# If ashr file: use shrunk_logFC as primary logFC, lfsr as significance metric
if ("shrunk_logFC" %in% names(dream)) {
  dream[, logFC_raw := logFC]
  dream[, logFC := shrunk_logFC]
  cat("    Using ashr shrunk_logFC and lfsr for significance\n")
}

# Dream has versioned Ensembl IDs (e.g. ENSG00000154978.14)
dream[, ensembl_base := sub("\\..*$", "", gene)]

# --- Ensembl → HGNC mapping: try local cache first, then biomaRt ---
ann_cache_path <- file.path(PROJECT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")
ensembl_to_hgnc <- NULL

if (file.exists(ann_cache_path)) {
  cat("    Using local annotation cache for Ensembl → HGNC mapping...\n")
  ann_cache <- fread(ann_cache_path)
  bm <- ann_cache[symbol != "" & !is.na(symbol), .(ensembl_base = gene_base, hgnc_symbol = symbol)]
  bm <- bm[!duplicated(ensembl_base)]
  cat(sprintf("    Local cache: %d Ensembl → HGNC mappings\n", nrow(bm)))
  ensembl_to_hgnc <- bm
} else {
  cat("    Local annotation cache not found. Querying biomaRt...\n")
  ensembl_to_hgnc <- tryCatch({
    mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
                       mirror = "useast")
    ids <- unique(dream$ensembl_base)
    # Query in batches to avoid timeout on large sets
    batch_size <- 5000
    n_batches  <- ceiling(length(ids) / batch_size)
    bm_list <- lapply(seq_len(n_batches), function(i) {
      start <- (i - 1) * batch_size + 1
      end   <- min(i * batch_size, length(ids))
      as.data.table(getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters    = "ensembl_gene_id",
        values     = ids[start:end],
        mart       = mart
      ))
    })
    bm <- rbindlist(bm_list)
    setnames(bm, c("ensembl_base", "hgnc_symbol"))
    bm <- bm[hgnc_symbol != ""]
    bm <- bm[!duplicated(ensembl_base)]
    cat(sprintf("    biomaRt returned %d Ensembl → HGNC mappings\n", nrow(bm)))
    bm
  }, error = function(e) {
    cat(sprintf("    WARNING: biomaRt failed: %s\n", conditionMessage(e)))
    NULL
  })
}

if (!is.null(ensembl_to_hgnc)) {
  dream <- merge(dream, ensembl_to_hgnc, by = "ensembl_base", all.x = TRUE)
} else {
  # If biomaRt totally failed, we still proceed but with NAs
  dream[, hgnc_symbol := NA_character_]
  cat("    WARNING: proceeding without Ensembl→HGNC mapping. Human-only tier will be sparse.\n")
}

# Map GWAS symbols against dream
gwas_in_dream <- dream[hgnc_symbol %in% gwas_symbols]
gwas_in_dream <- gwas_in_dream[order(padj)]
gwas_in_dream <- gwas_in_dream[!duplicated(hgnc_symbol)]  # best padj per symbol

human_mapping <- gwas_in_dream[, .(
  gwas_symbol    = hgnc_symbol,
  ensembl_id     = gene,
  dream_logFC    = logFC,
  dream_padj     = padj,
  dream_AveExpr  = AveExpr,
  dream_t        = t
)]

n_mapped_dream <- nrow(human_mapping)
n_unmapped_dream <- length(gwas_symbols) - n_mapped_dream
cat(sprintf("    GWAS genes mapped to dream: %d / %d  (%d unmapped)\n",
            n_mapped_dream, length(gwas_symbols), n_unmapped_dream))

# Report unmapped GWAS genes
unmapped_dream <- setdiff(gwas_symbols, human_mapping$gwas_symbol)
if (length(unmapped_dream) > 0) {
  cat(sprintf("    Unmapped examples: %s\n",
              paste(head(unmapped_dream, 15), collapse = ", ")))
  if (length(unmapped_dream) > 15)
    cat(sprintf("    ... and %d more\n", length(unmapped_dream) - 15))
}

# Save human-only mapping
fwrite(human_mapping, file.path(OUT_DIR, "gwas_human_mapping.csv"))
cat(sprintf("    Saved: %s\n", file.path(OUT_DIR, "gwas_human_mapping.csv")))

# ============================================================
#  3. Cross-species overlay (Tier 2): concordance atlas
# ============================================================
cat("\n[3] Cross-species overlay via concordance atlas...\n")
atlas <- fread(ATLAS_FILE)
cat(sprintf("    Atlas: %d genes, %d categories\n",
            nrow(atlas), uniqueN(atlas$primary_category)))

atlas[, is_gwas := human_symbol %in% gwas_symbols]
gwas_in_atlas <- atlas[is_gwas == TRUE]
n_atlas_mapped <- uniqueN(gwas_in_atlas$human_symbol)

cat(sprintf("    GWAS genes in atlas: %d / %d\n",
            n_atlas_mapped, length(gwas_symbols)))

# Build Tier-2 table
atlas_mapping <- gwas_in_atlas[, .(
  gwas_symbol            = human_symbol,
  primary_category       = primary_category,
  best_n_concordant      = best_n_concordant,
  n_diets_sig            = n_diets_sig,
  mean_h_lfc             = mean_h_lfc,
  translatability_score  = translatability_score,
  translatability_tier   = translatability_tier,
  h_significant          = h_significant,
  n_concordant           = n_concordant,
  n_discordant           = n_discordant
)]

# ============================================================
#  4. Continuous enrichment tests (Wilcoxon + Permutation)
# ============================================================
cat("\n[4] Running continuous enrichment tests...\n")

# --- 4a. Wilcoxon rank-sum: translatability_score GWAS vs non-GWAS ---
cat("    [4a] Wilcoxon rank-sum on translatability_score...\n")
scores_gwas    <- atlas[is_gwas == TRUE  & !is.na(translatability_score), translatability_score]
scores_nongwas <- atlas[is_gwas == FALSE & !is.na(translatability_score), translatability_score]

if (length(scores_gwas) >= 5 && length(scores_nongwas) >= 5) {
  wt_global <- wilcox.test(scores_gwas, scores_nongwas, alternative = "greater")
  # Effect size: rank-biserial correlation r = 1 - 2U/(n1*n2)
  U_global  <- wt_global$statistic
  n1 <- length(scores_gwas)
  n2 <- length(scores_nongwas)
  r_global  <- 1 - (2 * U_global) / (n1 * n2)
  cat(sprintf("    Global Wilcoxon: W=%d, p=%.4g, rank-biserial r=%.4f\n",
              as.integer(U_global), wt_global$p.value, r_global))
  cat(sprintf("    Median trans-score: GWAS=%.3f, non-GWAS=%.3f\n",
              median(scores_gwas), median(scores_nongwas)))
} else {
  cat("    WARNING: Insufficient observations for global Wilcoxon test.\n")
  wt_global <- list(p.value = NA_real_, statistic = NA_real_)
  r_global  <- NA_real_
}

# Per-category Wilcoxon: for each category, compare translatability_score of
# GWAS-in-category vs GWAS-not-in-category (among atlas genes only)
categories <- sort(unique(atlas$primary_category))

wilcox_per_cat <- rbindlist(lapply(categories, function(cat_name) {
  in_cat    <- atlas[primary_category == cat_name & !is.na(translatability_score)]
  not_in_cat <- atlas[primary_category != cat_name & !is.na(translatability_score)]

  gwas_in   <- in_cat[is_gwas == TRUE, translatability_score]
  bg_in     <- in_cat[is_gwas == FALSE, translatability_score]

  n_gwas_cat <- length(gwas_in)
  n_bg_cat   <- length(bg_in)

  if (n_gwas_cat >= 3 && n_bg_cat >= 3) {
    wt <- wilcox.test(gwas_in, bg_in, alternative = "two.sided")
    U  <- wt$statistic
    r  <- 1 - (2 * U) / (n_gwas_cat * n_bg_cat)
  } else {
    wt <- list(p.value = NA_real_, statistic = NA_real_)
    r  <- NA_real_
  }

  data.table(
    category         = cat_name,
    n_gwas_in_cat    = n_gwas_cat,
    n_bg_in_cat      = n_bg_cat,
    median_gwas      = ifelse(n_gwas_cat > 0, median(gwas_in), NA_real_),
    median_bg        = ifelse(n_bg_cat > 0, median(bg_in), NA_real_),
    wilcox_W         = as.numeric(wt$statistic),
    wilcox_p         = wt$p.value,
    rank_biserial_r  = r
  )
}))

wilcox_per_cat[, wilcox_padj := p.adjust(wilcox_p, method = "fdr")]
wilcox_per_cat <- wilcox_per_cat[order(wilcox_p)]

cat("    Per-category Wilcoxon results:\n")
print(wilcox_per_cat[, .(category, n_gwas_in_cat, median_gwas, median_bg,
                         rank_biserial_r, wilcox_padj)])

# --- 4b. Permutation test: Conserved enrichment ---
cat("\n    [4b] Permutation test: Conserved enrichment...\n")
n_gwas_in_atlas  <- sum(atlas$is_gwas)
n_obs_core       <- sum(atlas$is_gwas & atlas$primary_category == "Conserved")
n_atlas_total    <- nrow(atlas)
N_PERM           <- 10000

perm_counts <- replicate(N_PERM, {
  idx <- sample.int(n_atlas_total, n_gwas_in_atlas, replace = FALSE)
  sum(atlas$primary_category[idx] == "Conserved")
})

perm_p <- (sum(perm_counts >= n_obs_core) + 1) / (N_PERM + 1)
cat(sprintf("    Observed GWAS in Conserved: %d / %d atlas-mapped\n",
            n_obs_core, n_gwas_in_atlas))
cat(sprintf("    Permutation p-value (N=%d): %.4g\n", N_PERM, perm_p))
cat(sprintf("    Permutation mean: %.1f, sd: %.2f\n",
            mean(perm_counts), sd(perm_counts)))

# Compile continuous enrichment results
enrichment_results <- data.table(
  test = c("wilcoxon_global", "permutation_conserved",
           paste0("wilcoxon_", wilcox_per_cat$category)),
  statistic = c(as.numeric(wt_global$statistic), n_obs_core,
                wilcox_per_cat$wilcox_W),
  effect_size = c(r_global, n_obs_core / n_gwas_in_atlas,
                  wilcox_per_cat$rank_biserial_r),
  p_value = c(wt_global$p.value, perm_p,
              wilcox_per_cat$wilcox_p),
  n_gwas = c(length(scores_gwas), n_gwas_in_atlas,
             wilcox_per_cat$n_gwas_in_cat),
  description = c(
    sprintf("GWAS vs non-GWAS translatability score (median %.3f vs %.3f)",
            median(scores_gwas, na.rm = TRUE), median(scores_nongwas, na.rm = TRUE)),
    sprintf("Observed %d in Conserved vs permutation mean %.1f",
            n_obs_core, mean(perm_counts)),
    sprintf("Within %s: GWAS median=%.3f vs bg median=%.3f",
            wilcox_per_cat$category,
            wilcox_per_cat$median_gwas,
            wilcox_per_cat$median_bg)
  )
)

fwrite(enrichment_results,
       file.path(OUT_DIR, "gwas_concordance_continuous_enrichment.csv"))
cat(sprintf("    Saved: %s\n",
            file.path(OUT_DIR, "gwas_concordance_continuous_enrichment.csv")))

# ============================================================
#  5. Evidence cards for 14 priority loci
# ============================================================
cat("\n[5] Building evidence cards for priority loci...\n")

priority_loci <- c("PNPLA3", "TM6SF2", "GCKR", "MBOAT7", "HSD17B13",
                    "FTO", "SERPINA1", "GPAM", "MTTP", "ALDH2",
                    "APOE", "MTARC1", "GATAD2A", "MAU2")

# Start with all priority loci as the scaffold
cards <- data.table(gwas_symbol = priority_loci)

# Merge human-only mapping (dream)
cards <- merge(cards, human_mapping, by = "gwas_symbol", all.x = TRUE)

# Merge atlas mapping
cards <- merge(cards, atlas_mapping, by = "gwas_symbol", all.x = TRUE)

# MR/TWAS merge REMOVED 2026-04-22 — MR ditched from paper.
# Initialize placeholder columns (kept as all-NA so downstream figure/tally
# code that references these names continues to function without errors).
# TWAS + COLOC + INTACT causal layers are surfaced via the main atlas; this
# script now reports transcriptomic + concordance evidence only.
cards[, c("twas_z", "twas_p", "twas_fdr", "mr_beta", "mr_p", "mr_fdr",
          "coloc_PP4", "mr_score") :=
         .(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
           NA_real_, NA_real_)]

# Add evidence summary flags
cards[, in_dream         := !is.na(dream_logFC)]
cards[, in_atlas         := !is.na(primary_category)]
# Standard padj-based significance (padj < 0.05, |logFC| > 0.5)
# dream_padj already in cards from line 184 column selection
cards[, dream_sig := !is.na(dream_padj) & dream_padj < 0.05 & !is.na(dream_logFC) & abs(dream_logFC) > 0.5]
cards[, in_conserved := primary_category == "Conserved" & !is.na(primary_category)]
# has_mr is now always FALSE (MR ditched 2026-04-22); column retained for compat
cards[, has_mr := FALSE]

# Count evidence lines
cards[, n_evidence_types := as.integer(in_dream) + as.integer(in_atlas)]

# Order by priority list
cards[, priority_order := match(gwas_symbol, priority_loci)]
cards <- cards[order(priority_order)]
cards[, priority_order := NULL]

fwrite(cards, file.path(OUT_DIR, "gwas_convergence_evidence_cards.csv"))
cat(sprintf("    Saved: %s\n",
            file.path(OUT_DIR, "gwas_convergence_evidence_cards.csv")))

# Print summary
cat("\n    Priority loci evidence summary:\n")
for (i in seq_len(nrow(cards))) {
  row <- cards[i]
  cat(sprintf("    - %-12s | dream_logFC=%6s | dream_padj=%8s | atlas=%s | n_concord=%s | trans=%.3s\n",
              row$gwas_symbol,
              ifelse(is.na(row$dream_logFC), "  --  ",
                     sprintf("%+.3f", row$dream_logFC)),
              ifelse(is.na(row$dream_padj), "   --   ",
                     sprintf("%.2e", row$dream_padj)),
              ifelse(is.na(row$primary_category), "NOT_IN_ATLAS",
                     row$primary_category),
              ifelse(is.na(row$best_n_concordant), "-",
                     as.character(row$best_n_concordant)),
              ifelse(is.na(row$translatability_score), "-- ",
                     sprintf("%.3f", row$translatability_score))))
}

# ============================================================
#  6. Figures
# ============================================================
cat("\n[6] Generating figures...\n")

# --- 6a. Enrichment bar plot: Wilcoxon effect sizes per category ---
cat("    [6a] Wilcoxon effect-size bar plot...\n")
plot_wilcox <- copy(wilcox_per_cat)
plot_wilcox[, sig_label := fifelse(is.na(wilcox_padj), "",
                           fifelse(wilcox_padj < 0.001, "***",
                           fifelse(wilcox_padj < 0.01,  "**",
                           fifelse(wilcox_padj < 0.05,  "*", "ns"))))]
# Order by effect size
plot_wilcox[, category := factor(category,
              levels = category[order(rank_biserial_r, na.last = TRUE)])]

# Direction color: positive r = GWAS scores higher
plot_wilcox[, direction := fifelse(rank_biserial_r > 0, "higher", "lower")]

p_wilcox <- ggplot(plot_wilcox[!is.na(rank_biserial_r)],
                   aes(x = category, y = rank_biserial_r, fill = category)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40",
             linewidth = 0.3) +
  geom_text(aes(label = sig_label,
                y = rank_biserial_r + sign(rank_biserial_r) * 0.03),
            size = 2.5) +
  scale_fill_manual(values = concordance_colors, drop = FALSE) +
  coord_flip() +
  labs(
    title    = "GWAS Gene Translatability Within Concordance Categories",
    subtitle = sprintf("Wilcoxon rank-biserial r; %d GWAS genes in atlas (FDR-corrected)",
                       n_gwas_in_atlas),
    x = NULL,
    y = "Rank-Biserial r (GWAS vs background)"
  ) +
  theme_pub

# --- 6b. Category distribution bar chart for GWAS genes ---
cat("    [6b] Category distribution bar chart...\n")
cat_dist <- gwas_in_atlas[, .N, by = primary_category]
setnames(cat_dist, "N", "count")
all_cats <- data.table(primary_category = names(concordance_colors))
cat_dist <- merge(all_cats, cat_dist, by = "primary_category", all.x = TRUE)
cat_dist[is.na(count), count := 0]
cat_dist[, pct := round(count / sum(count) * 100, 1)]
cat_dist[, primary_category := factor(primary_category,
           levels = cat_dist$primary_category[order(-count)])]

p_dist <- ggplot(cat_dist[count > 0],
                 aes(x = primary_category, y = count, fill = primary_category)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", count, pct)),
            vjust = -0.3, size = 2, lineheight = 0.8) +
  scale_fill_manual(values = concordance_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title    = "Distribution of GWAS Genes Across Concordance Categories",
    subtitle = sprintf("%d / %d GWAS genes mapped",
                       nrow(gwas_in_atlas), length(gwas_symbols)),
    x = NULL, y = "Number of GWAS genes"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.5))

# Combine Wilcoxon enrichment and Category Distribution into a single figure
p_combined_enrichment <- p_dist + p_wilcox + plot_layout(ncol = 2, widths = c(1, 1.2))
ggsave(file.path(FIG_DIR, "gwas_category_enrichment_combined_v3.pdf"),
       p_combined_enrichment, width = 7.5, height = 3.5, device = cairo_pdf)
cat(sprintf("    Saved: %s\n", file.path(FIG_DIR, "gwas_category_enrichment_combined_v3.pdf")))

# --- 6c. Convergence dot plot for 14 priority loci ---
cat("    [6c] Priority loci convergence dot plot...\n")

# Build a long-form evidence matrix for the dot plot
evidence_long <- rbindlist(list(
  cards[, .(gwas_symbol,
            evidence = "Human DE (dream)",
            value    = fifelse(in_dream, 1, 0),
            lfc      = dream_logFC,
            nlog10p  = -log10(dream_padj),
            detail   = fifelse(!is.na(dream_padj), sprintf("p=%.1e", dream_padj), ""))],
  cards[, .(gwas_symbol,
            evidence = "Concordance Atlas",
            value    = fifelse(in_atlas, 1, 0),
            lfc      = NA_real_,
            nlog10p  = NA_real_,
            detail   = fifelse(!is.na(primary_category), primary_category, ""))],
  cards[, .(gwas_symbol,
            evidence = "Conserved",
            value    = fifelse(in_conserved, 1, 0),
            lfc      = NA_real_,
            nlog10p  = NA_real_,
            detail   = "")],
  cards[, .(gwas_symbol,
            evidence = "Dream padj < 0.05",
            value    = fifelse(dream_sig, 1, 0),
            lfc      = NA_real_,
            nlog10p  = NA_real_,
            detail   = "")]
  # "MR/TWAS" evidence layer dropped 2026-04-22 — MR ditched from paper.
))

# Factor ordering: priority list for genes, logical order for evidence
evidence_long[, gwas_symbol := factor(gwas_symbol, levels = rev(priority_loci))]
evidence_long[, evidence := factor(evidence,
                levels = c("Human DE (dream)", "Dream padj < 0.05",
                           "Concordance Atlas", "Conserved"))]

# Plot mapping:
# For Human DE (dream), map color to LFC and size to significance
# For binary evidence layers, keep simple filled/unfilled circles
p_dot_base <- ggplot(evidence_long, aes(x = evidence, y = gwas_symbol)) +
  geom_point(data = evidence_long[evidence != "Human DE (dream)"],
             aes(fill = factor(value)), shape = 21, size = 3.5, stroke = 0.3, color = "grey30") +
  scale_fill_manual(
    values = c("0" = "white", "1" = pal$magenta),
    labels = c("Absent", "Present"),
    name   = "Binary Evidence"
  )

# Overlay continuous layer for Human DE (ggnewscale loaded at top)
p_dot <- ggplot(evidence_long, aes(x = evidence, y = gwas_symbol)) +
  # Background tiles for grid structure
  geom_tile(color = "white", fill = "transparent") +
  # 1. Binary points (all except Human DE)
  geom_point(data = evidence_long[evidence != "Human DE (dream)"],
             aes(fill = factor(value)), shape = 21, size = 3.5, stroke = 0.3, color = "grey30") +
  scale_fill_manual(values = c("0" = "white", "1" = pal$magenta), name = "Evidence") +
  # 2. Continuous points (Human DE) -> color by LFC, size by significance.
  # Cap nlog10p for point scaling to avoid huge blobs
  geom_point(data = evidence_long[evidence == "Human DE (dream)" & !is.na(lfc)],
             aes(color = lfc, size = pmin(nlog10p, 10)), shape = 16) +
  scale_color_gradient2(low = pal$teal, mid = "grey90", high = pal$magenta, midpoint = 0, name = "LogFC", limits=c(-1,1), oob=squish) +
  scale_size_continuous(range = c(1, 5), name = "-log10(p) [cap 10]") +
  labs(
    title    = "Multi-Evidence Convergence for Priority GWAS Loci",
    subtitle = "Integrated Clinical, Cross-Species, and Causal Triangulation",
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x  = element_text(angle = 40, hjust = 1, size = 6),
    axis.text.y  = element_text(face = "italic", size = 6),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    legend.position  = "right",
    legend.box = "vertical",
    legend.margin = margin(t=0, b=0, l=0, r=0)
  )

# Marginal bar plot for translatability score
trans_dt <- cards[, .(gwas_symbol, translatability_score)]
trans_dt[, gwas_symbol := factor(gwas_symbol, levels = rev(priority_loci))]

p_trans_bar <- ggplot(trans_dt, aes(x = translatability_score, y = gwas_symbol)) +
  geom_col(fill = pal$purple, width = 0.7) +
  scale_x_continuous(expand=c(0,0)) +
  labs(x = "Translatability Score", y = NULL, title = "") +
  theme_pub +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

p_combined_dot <- p_dot + p_trans_bar + plot_layout(ncol = 2, widths = c(4, 1))

ggsave(file.path(FIG_DIR, "gwas_priority_convergence_dotplot_v3.pdf"),
       p_combined_dot, width = 6.5, height = 4.5, device = cairo_pdf)
cat(sprintf("    Saved: %s\n",
            file.path(FIG_DIR, "gwas_priority_convergence_dotplot_v3.pdf")))

# ============================================================
#  Session summary
# ============================================================
cat("\n=== Summary ===\n")
cat(sprintf("  GWAS raw entries:              %d\n", nrow(gwas_raw)))
cat(sprintf("  After pipe-split + dedup:      %d unique HGNC symbols\n",
            length(gwas_symbols)))
cat(sprintf("  Mapped to dream (human-only):  %d / %d\n",
            n_mapped_dream, length(gwas_symbols)))
cat(sprintf("  Mapped to concordance atlas:   %d / %d\n",
            n_atlas_mapped, length(gwas_symbols)))
cat(sprintf("  Wilcoxon global p-value:       %.4g\n", wt_global$p.value))
cat(sprintf("  Permutation Conserved p:  %.4g (obs=%d, mean=%.1f)\n",
            perm_p, n_obs_core, mean(perm_counts)))
cat(sprintf("  Priority loci with dream data: %d / %d\n",
            sum(cards$in_dream), length(priority_loci)))
cat(sprintf("  Priority loci in atlas:        %d / %d\n",
            sum(cards$in_atlas), length(priority_loci)))
# MR/TWAS tally line removed 2026-04-22 — MR ditched from paper.
cat(sprintf("\nOutputs in: %s\n", OUT_DIR))
cat(sprintf("Figures in: %s\n", FIG_DIR))
cat(sprintf("\nCompleted: %s\n", Sys.time()))
