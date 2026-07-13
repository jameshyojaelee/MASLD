#!/usr/bin/env Rscript
# ============================================================================
# T0.15 / H6-inv audit: Conserved vs positive-control / drug-target /
# multi-evidence top-5% enrichment, under two null models and a selection-
# effect check on mRNA-protein rho.
#
# Concern (Team 2 adversarial): the reported ORs (16.47 DGIdb, 27.85 top-5%,
# ~18-ish vs positive-control panel) may be inflated by protein-coding-vs-
# lncRNA curation bias. Conserved is 100% protein-coding by construction
# (cross-species orthology requirement); if the *universe* in the Fisher test
# is also 100% protein-coding, biotype alone cannot drive the OR — but if the
# universe were the full 33,943-gene GENCODE atlas (12,671 lncRNAs etc.),
# biotype could be a large or dominant confounder.
#
# This audit:
#   (1) Inventories biotype distribution of each gene set.
#   (2) Recomputes OR under Null A (full 33,943-gene atlas universe)
#       and Null B (biotype-matched random draws, same universe).
#   (3) Repeats on the original "concordance atlas" universe (11,611 genes,
#       99.9% protein-coding) for direct comparison with the published numbers.
#   (4) Checks the selection-effect on Conserved mRNA-protein rho
#       (628 of 1,108 CC have proteomics; is 0.563 inflated?).
#
# Output:
#   RNA-seq/results/audit_sensitivity/H6_inv_cc_positive_control.csv
#   RNA-seq/results/audit_sensitivity/H6_inv_cc_positive_control.rds
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

N_DRAWS <- 1000L

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
atlas <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select = c("human_symbol", "gene_biotype", "is_conserved"))
atlas <- atlas[!is.na(human_symbol) & human_symbol != ""]
atlas <- atlas[!duplicated(human_symbol)]

conc <- fread(file.path(BASE,
  "Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv"),
  select = c("human_symbol", "primary_category"))
conc <- conc[!is.na(human_symbol) & human_symbol != ""]
conc <- conc[!duplicated(human_symbol)]

me_scored <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/multi_evidence_scored_genes.csv"),
  select = c("human_symbol", "composite_score"))
me_scored <- me_scored[!is.na(human_symbol) & human_symbol != ""]

dgidb <- fread(file.path(BASE,
  "RNA-seq/results/drug_repurposing/dgidb_drug_gene_interactions.csv"),
  select = "gene")
drug_targets <- unique(dgidb$gene[!is.na(dgidb$gene) & dgidb$gene != ""])

pos_ctrl <- fread(file.path(BASE, "results/library/positive_control.csv"))
pos_ctrl_genes <- unique(pos_ctrl[[1]])
pos_ctrl_genes <- pos_ctrl_genes[!is.na(pos_ctrl_genes) & pos_ctrl_genes != ""]

# ---------------------------------------------------------------------------
# Derive gene sets
# ---------------------------------------------------------------------------
cc_genes <- conc[primary_category == "Conserved"]$human_symbol
top5_thresh <- quantile(me_scored$composite_score, 0.95, na.rm = TRUE)
top5_genes <- unique(me_scored[composite_score >= top5_thresh]$human_symbol)

# ---------------------------------------------------------------------------
# (1) Biotype distribution audit
# ---------------------------------------------------------------------------
collapse_biotype <- function(bt) {
  bt2 <- rep("other", length(bt))
  bt2[is.na(bt)] <- "unknown"
  bt2[bt == "protein_coding"] <- "protein_coding"
  bt2[bt == "lncRNA"] <- "lncRNA"
  bt2[grepl("pseudogene", bt, ignore.case = TRUE)] <- "pseudogene"
  bt2
}
atlas[, biotype_class := collapse_biotype(gene_biotype)]

biotype_of <- function(gene_set) {
  dt <- atlas[human_symbol %in% gene_set]
  tab <- table(dt$biotype_class)
  # also count genes not in atlas
  miss <- length(setdiff(gene_set, atlas$human_symbol))
  c(as.list(tab), list(not_in_atlas = miss, total = length(gene_set)))
}

biotype_report <- rbindlist(list(
  data.table(set = "full_atlas_universe",
             protein_coding = sum(atlas$biotype_class == "protein_coding"),
             lncRNA         = sum(atlas$biotype_class == "lncRNA"),
             pseudogene     = sum(atlas$biotype_class == "pseudogene"),
             other          = sum(atlas$biotype_class == "other"),
             unknown        = sum(atlas$biotype_class == "unknown"),
             total          = nrow(atlas)),
  data.table(set = "concordance_atlas_universe",
             protein_coding = sum(atlas[human_symbol %in% conc$human_symbol]$biotype_class == "protein_coding"),
             lncRNA         = sum(atlas[human_symbol %in% conc$human_symbol]$biotype_class == "lncRNA"),
             pseudogene     = sum(atlas[human_symbol %in% conc$human_symbol]$biotype_class == "pseudogene"),
             other          = sum(atlas[human_symbol %in% conc$human_symbol]$biotype_class == "other"),
             unknown        = sum(atlas[human_symbol %in% conc$human_symbol]$biotype_class == "unknown"),
             total          = uniqueN(conc$human_symbol)),
  data.table(set = "conserved",
             protein_coding = sum(atlas[human_symbol %in% cc_genes]$biotype_class == "protein_coding"),
             lncRNA         = sum(atlas[human_symbol %in% cc_genes]$biotype_class == "lncRNA"),
             pseudogene     = sum(atlas[human_symbol %in% cc_genes]$biotype_class == "pseudogene"),
             other          = sum(atlas[human_symbol %in% cc_genes]$biotype_class == "other"),
             unknown        = sum(atlas[human_symbol %in% cc_genes]$biotype_class == "unknown"),
             total          = length(cc_genes)),
  data.table(set = "positive_control_65",
             protein_coding = sum(atlas[human_symbol %in% pos_ctrl_genes]$biotype_class == "protein_coding"),
             lncRNA         = sum(atlas[human_symbol %in% pos_ctrl_genes]$biotype_class == "lncRNA"),
             pseudogene     = sum(atlas[human_symbol %in% pos_ctrl_genes]$biotype_class == "pseudogene"),
             other          = sum(atlas[human_symbol %in% pos_ctrl_genes]$biotype_class == "other"),
             unknown        = sum(atlas[human_symbol %in% pos_ctrl_genes]$biotype_class == "unknown"),
             total          = length(pos_ctrl_genes)),
  data.table(set = "dgidb_drug_targets",
             protein_coding = sum(atlas[human_symbol %in% drug_targets]$biotype_class == "protein_coding"),
             lncRNA         = sum(atlas[human_symbol %in% drug_targets]$biotype_class == "lncRNA"),
             pseudogene     = sum(atlas[human_symbol %in% drug_targets]$biotype_class == "pseudogene"),
             other          = sum(atlas[human_symbol %in% drug_targets]$biotype_class == "other"),
             unknown        = sum(atlas[human_symbol %in% drug_targets]$biotype_class == "unknown"),
             total          = length(drug_targets)),
  data.table(set = "multi_evidence_top5",
             protein_coding = sum(atlas[human_symbol %in% top5_genes]$biotype_class == "protein_coding"),
             lncRNA         = sum(atlas[human_symbol %in% top5_genes]$biotype_class == "lncRNA"),
             pseudogene     = sum(atlas[human_symbol %in% top5_genes]$biotype_class == "pseudogene"),
             other          = sum(atlas[human_symbol %in% top5_genes]$biotype_class == "other"),
             unknown        = sum(atlas[human_symbol %in% top5_genes]$biotype_class == "unknown"),
             total          = length(top5_genes))
))
fwrite(biotype_report, file.path(OUTDIR, "H6_inv_biotype_distribution.csv"))
cat("\n=== Biotype distributions ===\n"); print(biotype_report)

# ---------------------------------------------------------------------------
# (2) Fisher OR helper
# ---------------------------------------------------------------------------
fisher_or <- function(drawn, target, universe) {
  drawn    <- intersect(drawn, universe)
  target_u <- intersect(target, universe)
  a <- length(intersect(drawn, target_u))
  b <- length(drawn) - a
  c <- length(target_u) - a
  d <- length(universe) - a - b - c
  m <- matrix(c(a, b, c, d), nrow = 2)
  # Haldane-Anscombe correction when cells == 0
  if (any(c(a, b, c, d) == 0)) m <- m + 0.5
  ft <- tryCatch(fisher.test(m), error = function(e) NULL)
  if (is.null(ft)) return(list(OR = NA_real_, p = NA_real_, a = a, b = b, c = c, d = d))
  list(OR = unname(ft$estimate), p = ft$p.value, a = a, b = b, c = c, d = d)
}

# ---------------------------------------------------------------------------
# (3) Null A (full-universe random) + Null B (biotype-matched random)
#     Computed per universe (concordance atlas = published setting;
#     full atlas = what an lncRNA-inclusive universe would imply).
# ---------------------------------------------------------------------------
cc_biotype_vec <- atlas[human_symbol %in% cc_genes]$biotype_class
cc_in_atlas <- intersect(cc_genes, atlas$human_symbol)

run_null <- function(universe, comparator_name, target_set, n_draws = N_DRAWS) {
  k <- length(cc_in_atlas)
  obs <- fisher_or(cc_in_atlas, target_set, universe)

  # Null A: uniform random
  null_A <- numeric(n_draws)
  for (i in seq_len(n_draws)) {
    drawn <- sample(universe, size = k, replace = FALSE)
    null_A[i] <- fisher_or(drawn, target_set, universe)$OR
  }

  # Null B: biotype-matched
  u_bt <- atlas[human_symbol %in% universe, .(human_symbol, biotype_class)]
  bt_tbl <- table(cc_biotype_vec)
  null_B <- numeric(n_draws)
  for (i in seq_len(n_draws)) {
    drawn <- character(0)
    for (bt in names(bt_tbl)) {
      pool <- u_bt[biotype_class == bt]$human_symbol
      k_bt <- as.integer(bt_tbl[[bt]])
      if (length(pool) < k_bt) {
        # fall back: sample with replacement if pool is smaller than needed
        drawn <- c(drawn, sample(pool, size = k_bt, replace = TRUE))
      } else {
        drawn <- c(drawn, sample(pool, size = k_bt, replace = FALSE))
      }
    }
    null_B[i] <- fisher_or(drawn, target_set, universe)$OR
  }

  emp_p_A <- mean(null_A >= obs$OR, na.rm = TRUE)
  emp_p_B <- mean(null_B >= obs$OR, na.rm = TRUE)

  list(
    comparator = comparator_name,
    universe   = attr(universe, "name") %||% "universe",
    observed_OR = obs$OR,
    observed_p  = obs$p,
    a = obs$a, b = obs$b, c = obs$c, d = obs$d,
    null_A_median = median(null_A, na.rm = TRUE),
    null_A_95     = as.numeric(quantile(null_A, 0.95, na.rm = TRUE)),
    null_A_2.5    = as.numeric(quantile(null_A, 0.025, na.rm = TRUE)),
    null_A_97.5   = as.numeric(quantile(null_A, 0.975, na.rm = TRUE)),
    null_A_p      = emp_p_A,
    null_B_median = median(null_B, na.rm = TRUE),
    null_B_95     = as.numeric(quantile(null_B, 0.95, na.rm = TRUE)),
    null_B_2.5    = as.numeric(quantile(null_B, 0.025, na.rm = TRUE)),
    null_B_97.5   = as.numeric(quantile(null_B, 0.975, na.rm = TRUE)),
    null_B_p      = emp_p_B,
    n_draws       = n_draws,
    k_draw_size   = k,
    n_cc_in_universe = length(intersect(cc_in_atlas, universe))
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

universe_full  <- atlas$human_symbol
attr(universe_full, "name") <- "full_atlas"
universe_conc  <- intersect(atlas$human_symbol, conc$human_symbol)
attr(universe_conc, "name") <- "concordance_atlas"

comparators <- list(
  pos_ctrl     = pos_ctrl_genes,
  drug_target  = drug_targets,
  top5         = top5_genes
)

cat("\n=== Running null simulations (", N_DRAWS, "draws) ===\n")
all_rows <- list()
for (uni_name in c("concordance_atlas", "full_atlas")) {
  U <- if (uni_name == "concordance_atlas") universe_conc else universe_full
  attr(U, "name") <- uni_name
  for (cmp in names(comparators)) {
    cat(sprintf("  [%s | %s]\n", uni_name, cmp))
    res <- run_null(U, cmp, comparators[[cmp]], n_draws = N_DRAWS)
    res$universe <- uni_name
    all_rows[[paste(uni_name, cmp, sep = "_")]] <- res
  }
}

# ---------------------------------------------------------------------------
# (4) Selection-effect check on mRNA-protein rho
# ---------------------------------------------------------------------------
cat("\n=== Selection-effect check: mRNA-protein rho ===\n")
prot_sum <- fread(file.path(BASE,
  "Analysis/Proteomics/results/mrna_protein_concordance_summary.csv"))
# Use the main contrast (advanced vs early fibrosis from PXD052937 if present,
# otherwise fall back to the first dataset so the scalar is comparable to the
# published number). Use the union of all rows to estimate rho.
# The published rho=0.563 is from mrna_protein_concordance_stratified.csv for
# stratum = "Conserved"; the underlying data is mrna_protein_concordance_summary.csv
# NOTE: the mRNA effect column from mrna_protein_concordance_summary.csv (proteomics
# area producer, not migrated to bulk_*) is the mRNA side of the mRNA-protein
# concordance check, not the multi-evidence atlas. Reads below are kept against the
# on-disk column name and flagged for the CI gate. # C2-OK-sensitivity
prot <- prot_sum[!is.na(dream_logFC) & !is.na(protein_logFC)]  # C2-OK-sensitivity
# Collapse duplicates (same gene, multiple datasets) -> take mean per gene
prot_gene <- prot[, .(dream_logFC = mean(dream_logFC, na.rm = TRUE),  # C2-OK-sensitivity
                      protein_logFC = mean(protein_logFC, na.rm = TRUE),
                      is_conserved = any(isTRUE(is_conserved) |
                                              is_conserved %in% c("TRUE", TRUE, 1L))),
                  by = gene]
prot_gene <- prot_gene[!is.na(dream_logFC) & !is.na(protein_logFC)]  # C2-OK-sensitivity
cat("  Proteomics genes available:", nrow(prot_gene), "\n")

# Attach biotype
prot_gene <- merge(prot_gene, atlas[, .(human_symbol, biotype_class)],
                    by.x = "gene", by.y = "human_symbol", all.x = TRUE)

# (a) Published: CC subset rho
cc_prot <- prot_gene[gene %in% cc_genes]
rho_cc <- suppressWarnings(cor(cc_prot$dream_logFC, cc_prot$protein_logFC, method = "spearman"))  # C2-OK-sensitivity
cat(sprintf("  (a) CC subset rho (n=%d): %.3f\n", nrow(cc_prot), rho_cc))

# (b) Non-CC protein-coding — size-matched random draw (100 reps)
nonCC_pc <- prot_gene[!(gene %in% cc_genes) & biotype_class == "protein_coding"]
set.seed(42)
rho_b_vec <- numeric(N_DRAWS)
for (i in seq_len(N_DRAWS)) {
  smp <- nonCC_pc[sample(.N, min(nrow(cc_prot), .N))]
  rho_b_vec[i] <- suppressWarnings(cor(smp$dream_logFC, smp$protein_logFC,  # C2-OK-sensitivity
                                        method = "spearman"))
}
cat(sprintf("  (b) Non-CC protein-coding size-matched rho median=%.3f, 2.5-97.5%% [%.3f, %.3f]\n",
            median(rho_b_vec), quantile(rho_b_vec, 0.025), quantile(rho_b_vec, 0.975)))

# (c) Biotype-matched random draw (100 reps)
#     Resample prot_gene to match the biotype distribution of cc_prot.
cc_bt_tbl <- table(cc_prot$biotype_class)
rho_c_vec <- numeric(N_DRAWS)
for (i in seq_len(N_DRAWS)) {
  drawn <- character(0)
  for (bt in names(cc_bt_tbl)) {
    pool <- prot_gene[biotype_class == bt & !(gene %in% cc_genes)]$gene
    k_bt <- as.integer(cc_bt_tbl[[bt]])
    if (length(pool) == 0) next
    if (length(pool) < k_bt) {
      drawn <- c(drawn, sample(pool, k_bt, replace = TRUE))
    } else {
      drawn <- c(drawn, sample(pool, k_bt, replace = FALSE))
    }
  }
  smp <- prot_gene[gene %in% drawn]
  rho_c_vec[i] <- suppressWarnings(cor(smp$dream_logFC, smp$protein_logFC,  # C2-OK-sensitivity
                                        method = "spearman"))
}
cat(sprintf("  (c) Biotype-matched rho median=%.3f, 2.5-97.5%% [%.3f, %.3f]\n",
            median(rho_c_vec), quantile(rho_c_vec, 0.025), quantile(rho_c_vec, 0.975)))

# empirical p: how often did a random draw exceed CC's rho?
emp_p_b <- mean(rho_b_vec >= rho_cc, na.rm = TRUE)
emp_p_c <- mean(rho_c_vec >= rho_cc, na.rm = TRUE)

selection_row <- data.table(
  test = "C_selection",
  comparator = "mRNA_protein_rho",
  observed_OR = rho_cc,  # reusing OR column for rho (scalar)
  null_median_OR = median(rho_c_vec),
  null_2.5_OR = quantile(rho_c_vec, 0.025),
  null_97.5_OR = quantile(rho_c_vec, 0.975),
  emp_p = emp_p_c,
  n_draws = N_DRAWS,
  n_cc_genes = nrow(cc_prot),
  k_draw_size = nrow(cc_prot),
  universe = "proteomics_detected",
  notes = sprintf("(b) non-CC pc rho med=%.3f; (c) biotype-matched rho med=%.3f; published CC rho=0.563 (n=628 in Proteomics/mrna_protein_concordance_stratified.csv)",
                  median(rho_b_vec), median(rho_c_vec))
)

# ---------------------------------------------------------------------------
# Assemble final CSV
# ---------------------------------------------------------------------------
rows <- rbindlist(lapply(all_rows, function(r) {
  data.table(
    test = c("A_fullrand", "B_biotypematched"),
    comparator = r$comparator,
    observed_OR = r$observed_OR,
    null_median_OR = c(r$null_A_median, r$null_B_median),
    null_2.5_OR = c(r$null_A_2.5, r$null_B_2.5),
    null_97.5_OR = c(r$null_A_97.5, r$null_B_97.5),
    emp_p = c(r$null_A_p, r$null_B_p),
    n_draws = r$n_draws,
    n_cc_genes = r$n_cc_in_universe,
    k_draw_size = r$k_draw_size,
    universe = r$universe,
    notes = sprintf("observed_p=%.3g; a=%d,b=%d,c=%d,d=%d",
                    r$observed_p, r$a, r$b, r$c, r$d)
  )
}), use.names = TRUE, fill = TRUE)

rows <- rbind(rows, selection_row, use.names = TRUE, fill = TRUE)
fwrite(rows, file.path(OUTDIR, "H6_inv_cc_positive_control.csv"))
saveRDS(list(rows = rows, biotype_report = biotype_report,
             rho_b_vec = rho_b_vec, rho_c_vec = rho_c_vec,
             rho_cc = rho_cc, n_cc_prot = nrow(cc_prot)),
        file.path(OUTDIR, "H6_inv_cc_positive_control.rds"))

cat("\n=== Output CSV ===\n"); print(rows)
cat("\nSaved:\n  ", file.path(OUTDIR, "H6_inv_cc_positive_control.csv"),
    "\n  ", file.path(OUTDIR, "H6_inv_cc_positive_control.rds"), "\n")
cat("End time:", format(Sys.time()), "\n")
