# WGCNA Pipeline - Step 05: Cross-Dataset Comparison
# Systematic comparison of modules across Patient, Other_MCD, and In-House_MCD datasets.

library(WGCNA)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# ==============================================================================
# CONFIGURATION
# ==============================================================================
BASE_DIR <- "../results"
OUTPUT_DIR <- "../results/05_comparison"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Datasets to compare
datasets <- c("Patient", "Other_MCD", "InHouse_MCD")

# ==============================================================================
# LOAD DATA
# ==============================================================================
# We need Expression Data (datExpr) and Module Colors (moduleColors) for each
multiExpr <- list()
multiColor <- list()

for (ds in datasets) {
  cat("Loading", ds, "...\n")
  # Load Expression
  load(file.path(BASE_DIR, "01_preprocessing", paste0(ds, ".RData")))
  # Renaming to generic
  if (ds == "Patient") {
      d <- datExpr_patient
  } else if (ds == "Other_MCD") {
      d <- datExpr_other
  } else if (ds == "InHouse_MCD") {
      d <- datExpr_inhouse
  }
  
  # Normalize gene names for cross-species (ToUpper)
  colnames(d) <- toupper(colnames(d))
  
  multiExpr[[ds]] <- list(data = d)
  
  # Load Network Results
  # Note: You must have run Step 03 for each dataset! 
  # For now, we assume the output file contains 'moduleColors'
  # We might need to adjust paths if you run them in subfolders.
  # Let's assume standard names or we need to standardize Step 03 output names.
  net_file <- file.path(BASE_DIR, "03_network_construction", paste0("WGCNA_network_", ds, ".RData"))
  if (file.exists(net_file)) {
      load(net_file)
      multiColor[[ds]] <- moduleColors
  } else {
      warning(paste("Network file not found for", ds, "- skipping module loading."))
      multiColor[[ds]] <- NULL
  }
}

# Check if we have data
if (length(multiExpr) < 2) stop("Need at least 2 datasets to compare.")

# ==============================================================================
# MODULE PRESERVATION (Z-STATISTICS)
# ==============================================================================
# Calculate preservation of Human (Patient) modules in Mouse datasets
cat("Calculating module preservation (Patient as Reference)...\n")

# Only run if Patient modules exist
if (!is.null(multiColor[["Patient"]])) {
    mp <- modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1, # 1 is Patient
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 3)
    
    save(mp, file = file.path(OUTPUT_DIR, "module_preservation_stats.RData"))
    
    # Plotting Z-summary
    ref <- 1
    for (test in 2:length(datasets)) {
        statsObs <- mp$quality$observed[[ref]][[test]]
        statsZ <- mp$quality$Z[[ref]][[test]]
        
        pdf(file.path(OUTPUT_DIR, paste0("Preservation_Patient_vs_", datasets[test], ".pdf")), 
            width = 7, height = 7)
        par(mfrow = c(1,1))
        # Z-summary
        modColors <- rownames(statsZ)
        moduleSizes <- statsZ$moduleSize
        plot(moduleSizes, statsZ$Zsummary, pch = 21, bg = modColors,
             main = paste("Preservation Z-summary: Patient vs", datasets[test]),
             xlab = "Module Size", ylab = "Zsummary")
        abline(h = 0, col = "grey")
        abline(h = 2, col = "blue", lty = 2) # Weak evidence
        abline(h = 10, col = "darkgreen", lty = 2) # Strong evidence
        dev.off()
    }
}

# ==============================================================================
# MODULE OVERLAP (HYPERGEOMETRIC TEST)
# ==============================================================================
cat("Calculating module overlap...\n")

# Define overlap function
calculate_overlap <- function(colors1, colors2, genes1, genes2) {
  common_genes <- intersect(genes1, genes2)
  if (length(common_genes) == 0) return(NULL)
  
  # Subset to common
  c1 <- colors1[match(common_genes, genes1)]
  c2 <- colors2[match(common_genes, genes2)]
  
  tbl <- table(c1, c2)
  pvals <- matrix(NA, nrow=nrow(tbl), ncol=ncol(tbl))
  rownames(pvals) <- rownames(tbl)
  colnames(pvals) <- colnames(tbl)
  
  universe <- length(common_genes)
  
  for(i in rownames(tbl)) {
    for(j in colnames(tbl)) {
      # Phyper(q, m, n, k)
      # q: overlap size - 1
      # m: size of module 1
      # n: universe - m
      # k: size of module 2
      q <- tbl[i,j]
      m <- sum(c1 == i)
      n <- universe - m
      k <- sum(c2 == j)
      pvals[i,j] <- phyper(q - 1, m, n, k, lower.tail = FALSE)
    }
  }
  return(list(pvals = pvals, counts = tbl))
}

# Pairwise comparisons
for (i in 1:(length(datasets)-1)) {
    for (j in (i+1):length(datasets)) {
        ds1 <- datasets[i]
        ds2 <- datasets[j]
        
        if (is.null(multiColor[[ds1]]) | is.null(multiColor[[ds2]])) next
        
        # Get genes
        g1 <- colnames(multiExpr[[ds1]]$data)
        g2 <- colnames(multiExpr[[ds2]]$data)
        
        # Upper Case for orthology matching (Simple approach)
        # Assuming Patient is Human and others Mouse
        if (ds1 == "Patient") g1 <- toupper(g1) else g1 <- toupper(g1)
        if (ds2 == "Patient") g2 <- toupper(g2) else g2 <- toupper(g2)
        
        res <- calculate_overlap(multiColor[[ds1]], multiColor[[ds2]], g1, g2)
        
        if (!is.null(res)) {
            # Plot Heatmap of -log10 P-values
            pdf(file.path(OUTPUT_DIR, paste0("Overlap_", ds1, "_vs_", ds2, ".pdf")), width = 10, height = 10)
            
            p_mat <- -log10(res$pvals + 1e-50) # Avoid inf
            Heatmap(p_mat, 
                    name = "-log10(P)", 
                    column_title = paste(ds2, "Modules"), 
                    row_title = paste(ds1, "Modules"),
                    col = colorRamp2(c(0, 10, 50), c("white", "yellow", "red")))
            dev.off()
        }
    }
}

cat("Comparison complete.\n")
