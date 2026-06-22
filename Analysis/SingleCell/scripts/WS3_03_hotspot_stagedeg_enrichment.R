#!/usr/bin/env Rscript
# =============================================================================
# WS3_03 — HOTSPOT module axis of the stage-DEG routing analysis.
#
# QUESTION (DATA-DECIDES framing): do the stage-specific bulk DEGs cluster into
# stage-associated per-cell-type Hotspot autocorrelation modules, and which cell
# types' modules carry them?
#
# DESIGN
#   - REUSE the existing per-cell-type Hotspot modules (501 GPU discovery is NOT
#     re-run). Modules are stage-AGNOSTIC autocorrelation structure.
#   - RE-FIT the stage-association statistic (505 logic) on TWO axes:
#       * disease_stage_coarse (ordinal 0-3)  -> PRIMARY axis
#       * F_stage_augmented                   -> SECONDARY axis
#     via lmer(score ~ axis + (1|dataset)) on donor module scores (CPU, cheap;
#     reuses donor_scores_all.tsv, NO GPU). Same protocol-contamination filter
#     as canonical 505.
#   - A module is "stage-associated / progression-module" on an axis if its
#     stage-beta q < 0.05 AND beta > 0 (score strengthens with disease stage).
#   - HYPERGEOMETRIC enrichment of each stage-DEG set in each module; background
#     = that cell type's Hotspot autocorrelation gene universe (autocorr.tsv,
#     ~9-13k genes per CT). HEADLINE restricted to progression-modules.
#   - Direction check via pct_concordant_up / mean_bulk_logFC (all_modules.tsv).
#
# OUTPUTS (under stagedeg_routing/)
#   - phenotype_correlations_stagedeg.tsv     : re-fit stage-beta, both axes
#   - hotspot_module_stagedeg_enrichment.csv  : CT x module x stage_set enrichment
#   - hotspot_celltype_stagedeg_capture.csv   : CT x stage_set capture fractions
#
# ENV: rnaseq    SLURM: CPU only
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(purrr); library(stringr)
  library(lme4); library(lmerTest)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HS   <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
OUT  <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/disease_signatures/stagedeg_routing")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

STAGE_DEG <- file.path(OUT, "stage_deg_sets.csv")

# Cell types with canonical (post-protocol-fix) donor scores + module_genes +
# autocorr universe. donor_scores_all.tsv carries exactly these 5 after the
# 2026-05-22 protocol-contamination remediation (global / endothelial were not
# re-scored into the canonical file). We report which CTs have modules.
CELL_TYPES <- c("hepatocytes", "macrophages", "fibroblasts",
                "cholangiocytes", "tcells")

# 9-lineage denominator (from CLAUDE.md taxonomy) for honest reporting of which
# lineages have Hotspot modules available.
NINE_LINEAGES <- c("Hepatocytes", "Macrophages", "Fibroblasts", "Cholangiocytes",
                   "Endothelial cells", "T cells", "B cells", "Plasma cells",
                   "NK cells")
CT_TO_LINEAGE <- c(hepatocytes = "Hepatocytes", macrophages = "Macrophages",
                   fibroblasts = "Fibroblasts", cholangiocytes = "Cholangiocytes",
                   tcells = "T cells")

# -----------------------------------------------------------------------------
# 1. RE-FIT stage-association (505 logic) on disease_stage_coarse + F_stage_augmented
# -----------------------------------------------------------------------------
message("== Step 1: re-fit module stage-beta on coarse + augmented axes ==")

donor_meta <- read_tsv(
  file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"),
  show_col_types = FALSE) |>
  mutate(disease_stage_ordinal = case_when(
    disease_stage_coarse == "Healthy"         ~ 0,
    disease_stage_coarse == "Steatosis"       ~ 1,
    disease_stage_coarse == "Steatohepatitis" ~ 2,
    disease_stage_coarse == "Cirrhosis"       ~ 3,
    TRUE ~ NA_real_
  ))

# Backwards-compat (mirror 505): treat absent protocol flag as not-excluded.
if (!"exclude_stage_analysis" %in% names(donor_meta)) {
  donor_meta$exclude_stage_analysis <- FALSE
}

scores <- read_tsv(file.path(HS, "donor_scores_all.tsv"), show_col_types = FALSE)

joined <- scores |>
  inner_join(donor_meta, by = "sample") |>
  filter(!is.na(score))

# Protocol-contamination exclusion (GSE136103 / Liver_Atlas CD45+ enrichment).
n_pre <- dplyr::n_distinct(joined$sample)
joined <- joined |>
  filter(is.na(exclude_stage_analysis) | exclude_stage_analysis == FALSE)
n_post <- dplyr::n_distinct(joined$sample)
message(sprintf("  protocol filter: %d donors retained (was %d, dropped %d)",
                n_post, n_pre, n_pre - n_post))

# Two axes requested by the task.
axes <- list(
  disease_stage_coarse = "disease_stage_ordinal",  # PRIMARY
  F_stage_augmented    = "F_stage_augmented"        # SECONDARY
)

# Ensure both axis columns are numeric in `joined` (disease_stage_ordinal is
# already numeric from case_when; F_stage_augmented may read as character).
for (axis_col in unname(unlist(axes))) {
  if (axis_col %in% names(joined)) {
    joined[[axis_col]] <- suppressWarnings(as.numeric(joined[[axis_col]]))
  }
}

fit_one <- function(df, axis_col) {
  if (sum(!is.na(df[[axis_col]])) < 20 || dplyr::n_distinct(df$dataset) < 2) {
    return(tibble(beta = NA_real_, SE = NA_real_, t = NA_real_, p = NA_real_,
                  n = sum(!is.na(df[[axis_col]]))))
  }
  fm <- as.formula(sprintf("score ~ %s + (1|dataset)", axis_col))
  m <- tryCatch(lmerTest::lmer(fm, data = df, REML = FALSE), error = function(e) NULL)
  if (is.null(m)) return(tibble(beta = NA_real_, SE = NA_real_, t = NA_real_,
                                p = NA_real_, n = nrow(df)))
  cf <- tryCatch(summary(m)$coefficients, error = function(e) NULL)
  if (is.null(cf) || !(axis_col %in% rownames(cf))) {
    return(tibble(beta = NA_real_, SE = NA_real_, t = NA_real_, p = NA_real_,
                  n = nrow(df)))
  }
  row <- cf[axis_col, ]
  p_val <- if ("Pr(>|t|)" %in% colnames(cf)) row["Pr(>|t|)"] else NA_real_
  tibble(beta = unname(row["Estimate"]), SE = unname(row["Std. Error"]),
         t = unname(row["t value"]), p = unname(p_val), n = nrow(df))
}

pheno <- tibble()
for (ax in names(axes)) {
  axis_col <- axes[[ax]]
  per_module <- joined |>
    group_by(cell_type, module) |>
    group_modify(~ fit_one(.x, axis_col)) |>
    ungroup() |>
    mutate(axis = ax, q = p.adjust(p, method = "BH"))
  pheno <- bind_rows(pheno, per_module)
}

# Strip the "<cell_type>__" namespace so module is a plain integer (matches
# module_genes.tsv). 508 does the same.
pheno <- pheno |>
  mutate(module_int = as.integer(sub("^.*__", "", module)))

write_tsv(pheno, file.path(OUT, "phenotype_correlations_stagedeg.tsv"))
message(sprintf("  wrote phenotype_correlations_stagedeg.tsv: %d rows", nrow(pheno)))

# Per-(axis, ct, module) stage-association flag.
# progression_module := stage-beta q < 0.05 AND beta > 0 (strengthens with stage).
pheno_flag <- pheno |>
  transmute(axis, cell_type, module = module_int,
            stage_beta = beta, stage_q = q, stage_n = n,
            is_progression_module = !is.na(q) & q < 0.05 & !is.na(beta) & beta > 0)

# -----------------------------------------------------------------------------
# 2. Load module gene sets, autocorr universes, all_modules direction columns
# -----------------------------------------------------------------------------
message("== Step 2: load module gene sets + autocorr universes ==")

# A hotspot gene matches a stage_deg_set row if it equals the row's `symbol`
# (real symbol when one exists) OR its unversioned ENSG `gene`. Both files use
# the identical "symbol-when-available, else ENSG" namespace.
module_genes <- map_dfr(CELL_TYPES, function(ct) {
  p <- file.path(HS, ct, "module_genes.tsv")
  if (!file.exists(p)) { message(sprintf("  [WARN] missing %s", p)); return(tibble()) }
  read_tsv(p, show_col_types = FALSE) |>
    transmute(cell_type = ct, module = as.integer(module), gene = gene)
})

universe <- map(CELL_TYPES, function(ct) {
  p <- file.path(HS, ct, "autocorr.tsv")
  read_tsv(p, show_col_types = FALSE)$gene
})
names(universe) <- CELL_TYPES
for (ct in CELL_TYPES) message(sprintf("  %s: universe=%d genes, modules=%d, module_genes=%d",
  ct, length(universe[[ct]]), dplyr::n_distinct(module_genes$module[module_genes$cell_type==ct]),
  sum(module_genes$cell_type==ct)))

# all_modules.tsv direction columns (mean_bulk_logFC = C2/LVQW canonical bulk
# logFC; pct_concordant_up = fraction of module genes up in bulk).
all_modules <- read_tsv(file.path(HS, "all_modules.tsv"), show_col_types = FALSE) |>
  transmute(cell_type, module = as.integer(module),
            mean_bulk_logFC, pct_concordant_up, bulk_replicated, is_novel)

# -----------------------------------------------------------------------------
# 3. Stage-DEG sets
# -----------------------------------------------------------------------------
message("== Step 3: stage-DEG sets ==")
sdeg <- read_csv(STAGE_DEG, show_col_types = FALSE)

# Build the two matchable keys per row.
sdeg <- sdeg |>
  mutate(key_symbol = symbol, key_ensg = gene)

stage_sets <- sort(unique(sdeg$stage_set))
message(sprintf("  %d stage sets: %s", length(stage_sets), paste(stage_sets, collapse=", ")))

# For a given cell-type universe (a character vector of hotspot genes), determine
# which universe genes are DEGs in a given stage_set. A universe gene `g` is a
# DEG if any stage_set row with is_deg==TRUE has key_symbol==g OR key_ensg==g.
deg_lookup <- function(stage_set_df) {
  degs <- stage_set_df |> filter(is_deg)
  unique(c(degs$key_symbol, degs$key_ensg))
}
# Also a "tested" lookup (the stage-set gene universe = all rows, for reference).
tested_lookup <- function(stage_set_df) {
  unique(c(stage_set_df$key_symbol, stage_set_df$key_ensg))
}

# Pre-split stage sets.
sdeg_split <- split(sdeg, sdeg$stage_set)
deg_keys   <- map(sdeg_split, deg_lookup)
deg_dir    <- map(sdeg_split, function(d) {
  dd <- d |> filter(is_deg)
  setNames(dd$direction, dd$key_symbol)  # direction by symbol key (sufficient for summary)
})

# -----------------------------------------------------------------------------
# 4. Hypergeometric enrichment: CT x module x stage_set
# -----------------------------------------------------------------------------
message("== Step 4: hypergeometric enrichment (CT x module x stage_set) ==")

# Restrict every test to the cell type's autocorr universe (the genes Hotspot
# actually tested for autocorrelation). DEGs and module genes are intersected
# with this universe before the hypergeometric test, so the background is
# strictly the ~9-13k autocorr-tested genes per CT.
enrich_rows <- list()
ri <- 0L
for (ct in CELL_TYPES) {
  uni <- universe[[ct]]
  N   <- length(uni)
  ct_mods <- module_genes |> filter(cell_type == ct)
  mod_ids <- sort(unique(ct_mods$module))
  for (ss in stage_sets) {
    deg_in_uni <- intersect(deg_keys[[ss]], uni)  # white balls in urn
    K <- length(deg_in_uni)
    if (K == 0) next
    for (mid in mod_ids) {
      mod_genes <- ct_mods$gene[ct_mods$module == mid]
      mod_in_uni <- intersect(mod_genes, uni)
      n_draw <- length(mod_in_uni)            # balls drawn
      if (n_draw == 0) next
      hits <- length(intersect(mod_in_uni, deg_in_uni))  # white balls drawn
      # hypergeometric upper-tail (enrichment): P(X >= hits)
      p_enrich <- phyper(hits - 1L, K, N - K, n_draw, lower.tail = FALSE)
      # odds ratio (2x2): module-vs-not x DEG-vs-not within universe
      a <- hits; b <- n_draw - hits; c <- K - hits; d <- N - K - (n_draw - hits)
      OR <- (a * d) / (b * c)
      ri <- ri + 1L
      enrich_rows[[ri]] <- tibble(
        cell_type = ct, module = mid, stage_set = ss,
        universe_N = N, deg_in_universe_K = K, module_size = n_draw,
        overlap = hits, expected = n_draw * K / N, OR = OR, p = p_enrich
      )
    }
  }
}
enrich <- bind_rows(enrich_rows)

# BH within each axis x stage_set across all CT-modules tested (multiplicity is
# over modules for a given stage set). Axis matters because progression-module
# status differs per axis; attach both axes via a cross-join on the flags.
enrich <- enrich |>
  left_join(all_modules, by = c("cell_type", "module"))

# Attach the per-axis progression flag (each enrichment row x each axis).
enrich_axis <- pheno_flag |>
  select(axis, cell_type, module, stage_beta, stage_q, is_progression_module) |>
  inner_join(enrich, by = c("cell_type", "module"),
             relationship = "many-to-many")

# q-value: BH across all (CT, module) tested for a given (axis, stage_set).
enrich_axis <- enrich_axis |>
  group_by(axis, stage_set) |>
  mutate(q = p.adjust(p, method = "BH")) |>
  ungroup() |>
  arrange(axis, stage_set, q)

# Direction concordance: stage-set direction (majority of its DEGs) vs module
# bulk direction (pct_concordant_up). Report module-level direction signal.
enrich_axis <- enrich_axis |>
  mutate(module_dir = case_when(
            is.na(pct_concordant_up) ~ NA_character_,
            pct_concordant_up >= 0.5 ~ "up",
            TRUE ~ "down"),
         module_mean_logFC = mean_bulk_logFC)

write_csv(enrich_axis |>
            select(axis, cell_type, module, stage_set, universe_N,
                   deg_in_universe_K, module_size, overlap, expected, OR, p, q,
                   stage_beta, stage_q, is_progression_module,
                   module_dir, module_mean_logFC, pct_concordant_up,
                   bulk_replicated, is_novel),
          file.path(OUT, "hotspot_module_stagedeg_enrichment.csv"))
message(sprintf("  wrote hotspot_module_stagedeg_enrichment.csv: %d rows", nrow(enrich_axis)))

# -----------------------------------------------------------------------------
# 5. Per-cell-type "module carries stage-DEGs" capture statistics
# -----------------------------------------------------------------------------
message("== Step 5: per-cell-type stage-DEG capture ==")

SIG_Q <- 0.05

# For each (axis, cell_type, stage_set):
#   - n_prog_modules            : # progression-modules of this CT (axis-specific)
#   - n_prog_modules_enriched   : # progression-modules sig enriched (q<0.05) for set
#   - frac_prog_modules_enriched: enriched / progression-modules
#   - n_deg_in_universe         : DEGs of this set within CT universe (K)
#   - n_deg_captured_prog       : # of those DEGs in >=1 enriched progression-module
#   - frac_deg_captured_prog    : captured / K   <- HEADLINE capture fraction
#   - n_any_modules_enriched    : # of ALL modules (not just prog) enriched (context)
#   - frac_deg_captured_any     : DEGs in >=1 enriched module of any kind

# Pre-compute, per (axis, ct, stage_set), the set of enriched progression-modules
# and the genes they capture.
ct_universe_chr <- universe

capture_rows <- list()
ci <- 0L
for (ax in names(axes)) {
  flags_ax <- pheno_flag |> filter(axis == ax)
  for (ct in CELL_TYPES) {
    uni <- ct_universe_chr[[ct]]
    ct_mods <- module_genes |> filter(cell_type == ct)
    prog_mods <- flags_ax$module[flags_ax$cell_type == ct & flags_ax$is_progression_module]
    n_prog <- length(prog_mods)
    n_all  <- dplyr::n_distinct(ct_mods$module)
    for (ss in stage_sets) {
      deg_in_uni <- intersect(deg_keys[[ss]], uni)
      K <- length(deg_in_uni)
      if (K == 0) next
      sub <- enrich_axis |>
        filter(axis == ax, cell_type == ct, stage_set == ss)
      enr_prog <- sub |> filter(is_progression_module, q < SIG_Q)
      enr_any  <- sub |> filter(q < SIG_Q)
      # genes captured by enriched progression modules
      cap_prog_genes <- ct_mods$gene[ct_mods$module %in% enr_prog$module]
      cap_any_genes  <- ct_mods$gene[ct_mods$module %in% enr_any$module]
      n_cap_prog <- length(intersect(cap_prog_genes, deg_in_uni))
      n_cap_any  <- length(intersect(cap_any_genes,  deg_in_uni))
      ci <- ci + 1L
      capture_rows[[ci]] <- tibble(
        axis = ax, cell_type = ct, lineage = unname(CT_TO_LINEAGE[ct]),
        stage_set = ss,
        n_prog_modules = n_prog, n_all_modules = n_all,
        n_prog_modules_enriched = nrow(enr_prog),
        frac_prog_modules_enriched = if (n_prog > 0) nrow(enr_prog) / n_prog else NA_real_,
        n_any_modules_enriched = nrow(enr_any),
        deg_in_universe_K = K,
        n_deg_captured_prog = n_cap_prog,
        frac_deg_captured_prog = n_cap_prog / K,
        n_deg_captured_any = n_cap_any,
        frac_deg_captured_any = n_cap_any / K
      )
    }
  }
}
capture <- bind_rows(capture_rows) |>
  arrange(axis, stage_set, desc(frac_deg_captured_prog))

write_csv(capture, file.path(OUT, "hotspot_celltype_stagedeg_capture.csv"))
message(sprintf("  wrote hotspot_celltype_stagedeg_capture.csv: %d rows", nrow(capture)))

# -----------------------------------------------------------------------------
# 6. Console verification summary
# -----------------------------------------------------------------------------
message("\n================= VERIFICATION SUMMARY =================")

message("\n[Cell types with Hotspot modules] (vs 9-lineage denominator)")
message(sprintf("  available: %d / 9 -> %s",
  length(CELL_TYPES), paste(unname(CT_TO_LINEAGE[CELL_TYPES]), collapse=", ")))
message(sprintf("  missing lineages: %s",
  paste(setdiff(NINE_LINEAGES, unname(CT_TO_LINEAGE[CELL_TYPES])), collapse=", ")))

for (ax in names(axes)) {
  message(sprintf("\n[Axis = %s] progression-modules per cell type (q<0.05 & beta>0)", ax))
  pm <- pheno_flag |> filter(axis == ax, is_progression_module) |>
    count(cell_type, name = "n_prog")
  allm <- pheno_flag |> filter(axis == ax) |> count(cell_type, name = "n_mod")
  pm_full <- allm |> left_join(pm, by = "cell_type") |>
    mutate(n_prog = coalesce(n_prog, 0L))
  for (i in seq_len(nrow(pm_full))) {
    message(sprintf("  %-15s %d/%d modules stage-associated",
      pm_full$cell_type[i], pm_full$n_prog[i], pm_full$n_mod[i]))
  }
}

message("\n[HEADLINE: which CTs' progression-modules are most enriched for stage-DEGs]")
message("  (primary axis = disease_stage_coarse; mean frac of progression-modules enriched, across stage sets)")
hl <- capture |>
  filter(axis == "disease_stage_coarse", n_prog_modules > 0) |>
  group_by(cell_type) |>
  summarise(mean_frac_prog_enriched = mean(frac_prog_modules_enriched, na.rm = TRUE),
            mean_frac_deg_captured = mean(frac_deg_captured_prog, na.rm = TRUE),
            n_prog_modules = first(n_prog_modules), .groups = "drop") |>
  arrange(desc(mean_frac_deg_captured))
print(as.data.frame(hl))

message("\n[Per-transition capture: frac of stage-DEGs in >=1 progression-module, primary axis]")
pt <- capture |>
  filter(axis == "disease_stage_coarse") |>
  select(stage_set, cell_type, n_prog_modules, frac_deg_captured_prog) |>
  pivot_wider(names_from = cell_type, values_from = frac_deg_captured_prog)
print(as.data.frame(pt))

message("\n[Hepatocyte vs non-hepatocyte module-carrying comparison, primary axis]")
hep_vs <- capture |>
  filter(axis == "disease_stage_coarse", n_prog_modules > 0) |>
  mutate(grp = ifelse(cell_type == "hepatocytes", "hepatocyte", "non-hepatocyte")) |>
  group_by(grp, stage_set) |>
  summarise(frac_deg_captured_prog = mean(frac_deg_captured_prog, na.rm = TRUE),
            .groups = "drop") |>
  group_by(grp) |>
  summarise(mean_frac_deg_captured = mean(frac_deg_captured_prog, na.rm = TRUE),
            .groups = "drop")
print(as.data.frame(hep_vs))

message("\n[Known-biology recovery: fibroblast ECM-type module enrichment]")
# Top fibroblast progression-modules enriched for advanced-fibrosis stage sets.
fib_chk <- enrich_axis |>
  filter(axis == "disease_stage_coarse", cell_type == "fibroblasts",
         is_progression_module,
         stage_set %in% c("coarse_Cirrhosis", "fine_F3_F4", "prog_c3_adv_fib")) |>
  arrange(q) |>
  select(stage_set, module, overlap, module_size, OR, q, module_dir, module_mean_logFC) |>
  head(10)
print(as.data.frame(fib_chk))

message("\n================= DONE =================")
