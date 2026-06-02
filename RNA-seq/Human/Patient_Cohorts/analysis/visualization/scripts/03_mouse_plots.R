
# 03_mouse_plots.R
# Generates visualizations for Mouse Models (Trajectory PCA, Cross-Model comparisons)
# Updated: Feb 2026 (Sanjana Lab Theme)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
  library(grid)
  library(RColorBrewer)
})

# Paths
BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse"
MCD_DIR <- file.path(BASE_DIR, "InHouse_MCD")
INT_DIR <- file.path(BASE_DIR, "Unified_Integration")
RESULTS_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/visualization/plots/mouse"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Source plotting functions & theme
FUNC_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/visualization/functions"
source(file.path(FUNC_DIR, "theme_publication.R"))

# --- 1. In-House MCD Trajectory PCA ---
message("Generating In-House MCD Trajectory PCA...")

# Load Counts
counts_file <- file.path(MCD_DIR, "results/normalized_counts_all_samples.csv")
meta_file <- file.path(MCD_DIR, "metadata/samples.tsv")

if (file.exists(counts_file) && file.exists(meta_file)) {
    counts_df <- fread(counts_file)
    meta_df <- fread(meta_file)
    setDF(counts_df)
    setDF(meta_df)
    
    # Process counts (Symbol as rownames)
    # Check for duplicates
    counts_df <- counts_df[!duplicated(counts_df$gene_name), ]
    rownames(counts_df) <- counts_df$gene_name
    
    # Filter numeric columns (samples)
    sample_cols <- meta_df$sample_id
    sample_cols <- intersect(sample_cols, colnames(counts_df))
    
    counts_mat <- as.matrix(counts_df[, sample_cols])
    meta_df <- meta_df[meta_df$sample_id %in% sample_cols, ]
    rownames(meta_df) <- meta_df$sample_id
    
    # PCA
    # Log2 transform (counts are likely generic normalized, maybe TPM or similar? 
    # normalized_counts_all_samples.csv suggests normalized. 
    # Check max value. If > 1000, probably length corrected. Log2(x+1) for safety.
    pca_mat <- log2(counts_mat + 1)
    
    pca <- prcomp(t(pca_mat), center=TRUE, scale.=TRUE)
    percentVar <- round(100 * summary(pca)$importance[2, 1:2])
    
    d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], meta_df)
    d$week <- factor(d$week)
    d$diet <- factor(d$diet)
    
    # Trajectory Plot: Color by Diet, Shape by Week? Or Color by Week?
    # Guidelines: "Visualize disease progression trajectory relative to control."
    # Let's Color by Week (Gradient?) and Shape by Diet.
    # Or Color by Diet (Control=Blue, MCD=Orange/Magenta) and connect centroids by week.
    
    # Let's calculate centroids for centroids path
    centroids <- aggregate(cbind(PC1, PC2) ~ diet + week, data=d, FUN=mean)
    
    p_traj <- ggplot(d, aes(x=PC1, y=PC2, color=week, shape=diet)) +
        geom_point(size=3, alpha=0.8) +
        geom_path(data=centroids, aes(group=diet, color=week), linewidth=1, arrow=arrow(length=unit(0.2, "cm"))) + # Path doesn't make sense if color varies by week along path.
        # Better: Color by Diet, label Week.
        theme_publication()
        
    # Revised: Color by Diet, Text Label for Week (Centroids)
    p_traj <- ggplot(d, aes(x=PC1, y=PC2, color=diet, shape=as.factor(sex))) +
        geom_point(size=2, alpha=0.6) +
        scale_color_manual(values=c("Control"=sanjana_colors[["Blue"]], "MCD"=sanjana_colors[["Magenta"]])) +
        theme_publication() +
        labs(title="In-House MCD Timecourse PCA", shape="Sex", color="Diet") +
        xlab(paste0("PC1: ", percentVar[1], "%")) +
        ylab(paste0("PC2: ", percentVar[2], "%"))

    # Add centroids and trajectory
    centroids <- aggregate(cbind(PC1, PC2) ~ diet + week, data=d, FUN=mean)
    p_traj <- p_traj + 
        geom_point(data=centroids, aes(shape=NULL), size=5, alpha=1, color="black") +
        geom_point(data=centroids, aes(shape=NULL, color=diet), size=4, alpha=1) +
        geom_path(data=centroids, aes(group=diet, shape=NULL, color=diet), 
                  arrow=arrow(type="closed", length=unit(0.15, "inches")), linewidth=1.2) +
        geom_text_repel(data=centroids, aes(label=paste0("W", week), shape=NULL), 
                        color="black", size=4, fontface="bold", box.padding=0.6, point.padding=0.5)
        
    save_pdf(p_traj, file.path(RESULTS_DIR, "trajectory_pca.pdf"), width=6, height=5)
}

# --- 2. Public Model Comparison ---
message("Generating Public Model Comparison Plots...")

# Load Per-Diet Results
diet_dir <- file.path(INT_DIR, "results/per_diet")
files <- list.files(diet_dir, pattern="_de_results.csv", full.names=TRUE)

diets <- list()
diet_names <- c()

for (f in files) {
    name <- gsub("_de_results.csv", "", basename(f))
    if (name == "de_summary") next
    
    dt <- fread(f)
    setDF(dt)
    
    # Ensure columns
    # Likely 'gene', 'logFC', 'adj.P.Val' (limma) or 'log2FoldChange', 'padj' (DESeq2)
    # Map to standard
    if ("logFC" %in% names(dt)) dt $log2FoldChange <- dt$logFC
    if ("adj.P.Val" %in% names(dt)) dt$padj <- dt$adj.P.Val
    if (!"symbol" %in% names(dt) && "gene" %in% names(dt)) dt$symbol <- dt$gene
    
    diets[[name]] <- dt
    diet_names <- c(diet_names, name)
}

# Load Mapping for Mouse Symbols
map_file <- file.path(INT_DIR, "results/integration/human_mouse_ortholog_comparison.csv")
mouse_map <- NULL
if (file.exists(map_file)) {
    map_df <- fread(map_file)
    # Use mouse_gene_id -> mouse_symbol
    # mouse_gene_id in map seems to be without version (ENSMUSG...) based on standard Ensembl.
    # Let's check map content from previous head: ENSMUSG00000078652 (no version).
    mouse_map <- map_df[, c("mouse_gene_id", "mouse_symbol")]
    mouse_map <- mouse_map[!duplicated(mouse_map$mouse_gene_id), ]
    mouse_map <- mouse_map[mouse_map$mouse_gene_id != "", ]
}

# Apply mapping to all loaded diet DFs
for (name in diet_names) {
    dt <- diets[[name]]
    if ("gene" %in% names(dt)) {
        dt$gene_base <- gsub("\\..*", "", dt$gene)
        if (!is.null(mouse_map)) {
            dt <- merge(dt, mouse_map, by.x="gene_base", by.y="mouse_gene_id", all.x=TRUE)
            dt$symbol <- ifelse(is.na(dt$mouse_symbol) | dt$mouse_symbol == "", dt$gene, dt$mouse_symbol)
        } else {
            dt$symbol <- dt$gene
        }
    }
    diets[[name]] <- dt
}

# Jaccard Heatmap
message("  Generating Jaccard Heatmap...")
jaccard_mat <- matrix(NA, nrow=length(diet_names), ncol=length(diet_names), dimnames=list(diet_names, diet_names))

get_degs <- function(df) {
    # loose cutoff for overlap
    return(unique(df$symbol[df$padj < 0.05 & abs(df$log2FoldChange) > 1]))
}

deg_lists <- lapply(diets, get_degs)

for (i in 1:length(diet_names)) {
    for (j in 1:length(diet_names)) {
        s1 <- deg_lists[[diet_names[i]]]
        s2 <- deg_lists[[diet_names[j]]]
        
        uni <- length(union(s1, s2))
        inte <- length(intersect(s1, s2))
        
        jaccard_mat[i, j] <- ifelse(uni > 0, inte/uni, 0)
    }
}

pdf(file.path(RESULTS_DIR, "jaccard_heatmap.pdf"), width=5, height=5, useDingbats = FALSE)
Heatmap(jaccard_mat, 
        name="Jaccard Index", 
        col=colorRamp2(c(0, 1), c("white", sanjana_colors[["Magenta"]])),
        cluster_rows=TRUE, 
        cluster_columns=TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", jaccard_mat[i, j]), x, y, gp = gpar(fontsize = 8, fontfamily="Helvetica"))
        },
        column_title = "DEG Overlap (Jaccard)",
        column_title_gp = gpar(fontsize = 10, fontfamily="Helvetica", fontface="bold"))
dev.off()

# LFC Correlation (MCD vs HFD example)
# Find MCD and HFD
if ("MCD" %in% diet_names && "HFD" %in% diet_names) {
    message("  Generating LFC Scatter (MCD vs HFD)...")
    mcd <- diets[["MCD"]]
    hfd <- diets[["HFD"]]
    
    merged <- merge(mcd, hfd, by="symbol", suffixes=c(".MCD", ".HFD"))
    
    # Filter for significant in at least one
    sig_either <- merged[merged$padj.MCD < 0.05 | merged$padj.HFD < 0.05, ]
    
    cor_val <- cor(sig_either$log2FoldChange.MCD, sig_either$log2FoldChange.HFD, method="spearman")
    
    p_scatter <- ggplot(sig_either, aes(x=log2FoldChange.MCD, y=log2FoldChange.HFD)) +
        geom_point(alpha=0.3, size=1, color=sanjana_colors[["DarkGrey"]]) +
        geom_smooth(method="lm", color=sanjana_colors[["Blue"]], linewidth=0.5) +
        labs(title=paste0("LFC Correlation: MCD vs HFD (rho=", round(cor_val, 2), ")"),
             subtitle="Genes significant in at least one model") +
        theme_publication() +
        geom_hline(yintercept=0, linetype="dashed", color="grey") +
        geom_vline(xintercept=0, linetype="dashed", color="grey")
        
    save_pdf(p_scatter, file.path(RESULTS_DIR, "lfc_scatter_mcd_hfd.pdf"), width=5, height=5)
}

message("Done.")
