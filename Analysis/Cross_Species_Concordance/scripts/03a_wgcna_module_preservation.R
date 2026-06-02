#!/usr/bin/env Rscript
# 03a_wgcna_module_preservation.R
# ---------------------------------------------------------------------------
# WGCNA module preservation: build co-expression network on human, test
# whether modules preserve in each mouse diet subgroup
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
  library(ggplot2)
})

pdf.options(useDingbats = FALSE)
allowWGCNAThreads(nThreads = 16)

cat("=== Phase 3a: WGCNA Module Preservation ===\n\n")

BASE    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
H_INT   <- file.path(BASE, "Human/Patient_Cohorts/analysis/integration")
ANNOT   <- file.path(H_INT, "results/gene_annotation")
INT_DIR <- file.path(H_INT, "results/integration")
MOUSE_UI <- file.path(BASE, "Mouse/Unified_Integration")
WD      <- file.path(BASE, "Analysis/Cross_Species_Concordance")
RES     <- file.path(WD, "results")

DIETS <- c("MCD", "HFD", "CDAHFD", "FPC")  # 2026-05-29: LIDPAD archived/dropped

# ============================================================
#  Load ortholog mapping + expression data
# ============================================================
ortho <- fread(file.path(ANNOT, "ortholog_mapping.tsv"))
cat("1:1 orthologs:", nrow(ortho), "\n")

# Human expression (logCPM from merged_dge.rds, 1,444 samples, 10-cohort)
h_dge    <- readRDS(file.path(INT_DIR, "merged_dge.rds"))
h_logcpm <- edgeR::cpm(h_dge, log = TRUE, prior.count = 1)
h_meta   <- as.data.table(readRDS(file.path(INT_DIR, "meta_matched.rds")))

# Mouse expression
m_logcpm <- readRDS(file.path(MOUSE_UI, "results/corrected_logcpm.rds"))
m_meta   <- as.data.table(readRDS(file.path(MOUSE_UI, "results/meta_matched.rds")))

# ============================================================
#  Restrict to 1:1 orthologs (common gene space)
# ============================================================
h_base <- gsub("\\..*", "", rownames(h_logcpm))
m_base <- gsub("\\..*", "", rownames(m_logcpm))

# Map: human ENSG base → ortholog row
h_in_ortho <- h_base %in% ortho$human_gene_id
m_in_ortho <- m_base %in% ortho$mouse_gene_id

h_sub <- h_logcpm[h_in_ortho, ]
rownames(h_sub) <- h_base[h_in_ortho]

m_sub <- m_logcpm[m_in_ortho, ]
rownames(m_sub) <- m_base[m_in_ortho]

# Create common mapping
ortho_map <- ortho[human_gene_id %in% rownames(h_sub) & mouse_gene_id %in% rownames(m_sub)]
# Use human symbol as common ID
ortho_map[, common_id := human_symbol]
ortho_map <- ortho_map[common_id != "" & !is.na(common_id)]
ortho_map <- ortho_map[!duplicated(common_id)]  # one symbol per gene

h_common <- h_sub[ortho_map$human_gene_id, ]
rownames(h_common) <- ortho_map$common_id

m_common <- m_sub[ortho_map$mouse_gene_id, ]
rownames(m_common) <- ortho_map$common_id

cat(sprintf("Common gene space: %d genes\n", nrow(h_common)))

# ============================================================
#  Filter to most variable genes (top 5000 by variance)
# ============================================================
n_top <- 5000
h_var <- apply(h_common, 1, var)
top_genes <- names(sort(h_var, decreasing = TRUE))[1:min(n_top, length(h_var))]

h_wgcna <- t(h_common[top_genes, ])  # WGCNA needs samples × genes
m_wgcna <- t(m_common[top_genes, ])

cat(sprintf("WGCNA matrix: %d genes × %d human samples × %d mouse samples\n",
  length(top_genes), nrow(h_wgcna), nrow(m_wgcna)))

# ============================================================
#  Build human WGCNA network
# ============================================================
cat("\n=== Building Human Network ===\n")

# Pick soft-threshold
powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
sft <- pickSoftThreshold(h_wgcna, powerVector = powers, verbose = 0)
# Use first power where SFT R² > 0.80
best_power <- sft$fitIndices$Power[which(sft$fitIndices$SFT.R.sq > 0.80)[1]]
if (is.na(best_power)) best_power <- 6  # fallback
cat(sprintf("  Soft threshold power: %d\n", best_power))

# Build network
h_net <- blockwiseModules(
  h_wgcna, power = best_power,
  TOMType = "signed", networkType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  maxBlockSize = n_top, verbose = 0
)

h_modules <- h_net$colors
names(h_modules) <- colnames(h_wgcna)
n_mods <- length(unique(h_modules)) - 1  # exclude module 0 (unassigned)
cat(sprintf("  Human modules detected: %d (plus grey/unassigned)\n", n_mods))
print(table(h_modules))

# Save module assignments
mod_dt <- data.table(gene = names(h_modules), module = h_modules)
fwrite(mod_dt, file.path(RES, "wgcna_module_assignments.csv"))

# ============================================================
#  Module preservation in mouse (per diet subgroup)
# ============================================================
cat("\n=== Testing Module Preservation ===\n")

# Prepare multiExpr and multiColor
multiExpr <- list(human = list(data = h_wgcna))
multiColor <- list(human = h_modules)

preservation_results <- data.table()

for (diet in DIETS) {
  # Get mouse samples for this diet (disease + matching controls)
  diet_datasets <- m_meta[diet_model == diet, unique(dataset)]
  diet_samples <- m_meta[dataset %in% diet_datasets, sample_id]
  diet_samples <- intersect(diet_samples, rownames(m_wgcna))

  if (length(diet_samples) < 10) {
    cat(sprintf("  %s: only %d samples, skipping\n", diet, length(diet_samples)))
    next
  }

  m_diet <- m_wgcna[diet_samples, ]
  cat(sprintf("  %s: testing with %d samples...\n", diet, nrow(m_diet)))

  multiExpr_test <- list(
    human = list(data = h_wgcna),
    mouse = list(data = m_diet)
  )
  multiColor_test <- list(human = h_modules)

  mp <- tryCatch({
    modulePreservation(
      multiExpr_test, multiColor_test,
      referenceNetworks = 1, testNetworks = 2,
      nPermutations = 500,  # 500+ recommended for publishable Zsummary statistics
      randomSeed = 42, verbose = 0
    )
  }, error = function(e) {
    cat(sprintf("  %s: modulePreservation failed: %s\n", diet, e$message))
    NULL
  })

  if (!is.null(mp)) {
    # Extract Zsummary
    ref <- 1; test <- 2
    stats <- mp$preservation$Z[[ref]][[test]]
    p_stats <- mp$preservation$observed[[ref]][[test]]

    for (mod in rownames(stats)) {
      if (mod == "gold") next
      row <- data.table(
        diet = diet,
        module = as.integer(mod),
        module_size = p_stats[mod, "moduleSize"],
        Zsummary = stats[mod, "Zsummary.pres"],
        Zdensity = stats[mod, "Zdensity.pres"],
        Zconnectivity = stats[mod, "Zconnectivity.pres"],
        medianRank = mp$preservation$log.p[[ref]][[test]][mod, "log.p.medianRank.pres"]
      )
      preservation_results <- rbindlist(list(preservation_results, row), fill = TRUE)
    }
  }
}

if (nrow(preservation_results) > 0) {
  # Classify preservation
  preservation_results[, preservation := fifelse(
    Zsummary > 10, "Highly_Preserved",
    fifelse(Zsummary > 2, "Moderately_Preserved", "Not_Preserved")
  )]

  fwrite(preservation_results, file.path(RES, "wgcna_preservation_stats.csv"))
  cat("\nPreservation summary:\n")
  print(preservation_results[, .N, by = .(diet, preservation)])

  # Hub genes per preserved module
  hub_genes <- data.table()
  for (mod in unique(preservation_results[preservation == "Highly_Preserved", module])) {
    mod_genes <- names(h_modules[h_modules == mod])
    # Get connectivity (kME)
    kme <- cor(h_wgcna[, mod_genes], h_net$MEs[, paste0("ME", mod)], use = "pairwise.complete.obs")
    top_hub <- mod_genes[order(-abs(kme[, 1]))][1:min(20, length(mod_genes))]
    hub_genes <- rbindlist(list(hub_genes,
      data.table(module = mod, gene = top_hub, kME = kme[top_hub, 1])))
  }
  if (nrow(hub_genes) > 0) {
    fwrite(hub_genes, file.path(RES, "wgcna_hub_genes.csv"))
    cat("\nTop hub genes from preserved modules:\n")
    print(head(hub_genes, 20))
  }
}

cat("\nSaved: wgcna_preservation_stats.csv, wgcna_module_assignments.csv\n")
cat("=== Phase 3a complete ===\n")
