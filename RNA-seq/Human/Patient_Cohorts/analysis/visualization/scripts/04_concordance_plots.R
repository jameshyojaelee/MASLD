
# 04_concordance_plots.R
# Generates Cross-Species Concordance visualizations
# Updated: Feb 2026 (Sanjana Lab Theme)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggrepel)
})

# Paths
BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
CONC_DIR <- file.path(BASE_DIR, "Analysis/Cross_Species_Concordance/results")
HUMAN_DIR <- file.path(BASE_DIR, "Human/Patient_Cohorts/analysis/integration/results/integration")
MOUSE_DIR <- file.path(BASE_DIR, "Mouse/Unified_Integration/results/per_diet")
RESULTS_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/visualization/plots/concordance"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Source plotting functions & theme
FUNC_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/visualization/functions"
source(file.path(FUNC_DIR, "theme_publication.R"))

# --- 1. Concordance Matrix (Correlation Heatmap) ---
message("Generating Concordance Matrix Heatmap...")

conc_file <- file.path(CONC_DIR, "gene_concordance_matrix_20x.csv")
if (file.exists(conc_file)) {
    conc_df <- fread(conc_file)
    setDF(conc_df)
    
    # Reshape 'rho_sig' (Spearman rho of significant genes) to matrix
    # Columns: human_signature, diet, rho_sig
    
    mat_df <- reshape(conc_df[, c("human_signature", "diet", "rho_sig")], 
                      idvar = "diet", timevar = "human_signature", direction = "wide")
    rownames(mat_df) <- mat_df$diet
    mat_df$diet <- NULL
    colnames(mat_df) <- gsub("rho_sig.", "", colnames(mat_df))
    
    mat <- as.matrix(mat_df)
    mat[is.na(mat)] <- 0
    
    # Plot
    pdf(file.path(RESULTS_DIR, "concordance_matrix_rho.pdf"), width=5, height=5, useDingbats=FALSE)
    draw(Heatmap(mat,
            name = "Spearman Rho\n(Sig Genes)",
            col = colorRamp2(c(-0.5, 0, 0.5), c(sanjana_colors[["Blue"]], "white", sanjana_colors[["Magenta"]])),
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 8, fontfamily="Helvetica"))
            },
            column_title = "Cross-Species Concordance",
            column_title_gp = gpar(fontsize = 10, fontfamily="Helvetica", fontface="bold")
    ))
    dev.off()
}

# --- 2. Gene-Level Quadrant Scatter (Human vs MCD) ---
message("Generating Quadrant Scatterplot (Human vs MCD)...")

# --- 2. Gene-Level Quadrant Scatter (Human vs MCD) ---
message("Generating Quadrant Scatterplot (Human vs MCD)...")

conc_gene_file <- file.path(CONC_DIR, "gene_concordance_per_gene.csv")
mouse_file <- file.path(MOUSE_DIR, "MCD_de_results.csv") # Representative model

if (file.exists(conc_gene_file) && file.exists(mouse_file)) {
    c_res <- fread(conc_gene_file)
    m_res <- fread(mouse_file)
    setDF(c_res)
    setDF(m_res)
    
    # c_res has: mouse_gene_id, human_symbol, mean_h_lfc (Human LFC from meta-analysis is likely this or close to it)
    # m_res has: gene (ENSMUSG), logFC, adj.P.Val
    
    # Clean mouse data
    if ("logFC" %in% names(m_res)) m_res$logFC_Mouse <- m_res$logFC
    if ("adj.P.Val" %in% names(m_res)) m_res$padj_Mouse <- m_res$adj.P.Val
    # Use 'gene' as key - Strip version
    m_res$gene <- gsub("\\..*", "", m_res$gene)
    
    # Merge on mouse_gene_id
    merged <- merge(c_res, m_res, by.x="mouse_gene_id", by.y="gene")
    
    # We need significant genes in Human. 
    # c_res has 'h_significant' boolean? 
    # Checked head earlier: "h_significant,TRUE/FALSE".
    # Perfect.
    
    sig_human <- merged[merged$h_significant == TRUE, ]
    
    if (nrow(sig_human) > 0) {
        # Define Quadrants
        sig_human$Quadrant <- "Discordant"
        sig_human$Quadrant[sig_human$mean_h_lfc > 0 & sig_human$logFC_Mouse > 0] <- "Concordant Up"
        sig_human$Quadrant[sig_human$mean_h_lfc < 0 & sig_human$logFC_Mouse < 0] <- "Concordant Down"
        
        # Colors
        quad_colors <- c(
            "Concordant Up" = sanjana_colors[["Magenta"]], 
            "Concordant Down" = sanjana_colors[["Blue"]],
            "Discordant" = sanjana_colors[["Grey"]]
        )
        
        p_quad <- ggplot(sig_human, aes(x=mean_h_lfc, y=logFC_Mouse, color=Quadrant)) +
            geom_point(alpha=0.6, size=1.5) +
            geom_hline(yintercept=0, linetype="dashed", color="grey50") +
            geom_vline(xintercept=0, linetype="dashed", color="grey50") +
            theme_publication() +
            scale_color_manual(values=quad_colors) +
            labs(title="Human vs Mouse (MCD) Concordance",
                 subtitle=paste0("Genes Significant in Human Meta-Analysis (n=", nrow(sig_human), ")"),
                 x="Human log2 Fold Change",
                 y="Mouse MCD log2 Fold Change") +
            theme(legend.position="bottom")
            
        # Label top genes in each quadrant (highest combined LFC magnitude)
        sig_human$dist <- sqrt(sig_human$mean_h_lfc^2 + sig_human$logFC_Mouse^2)
        top_genes <- sig_human %>% group_by(Quadrant) %>% top_n(5, dist)
        
        # Need dplyr for pivot/group_by or base R
        # Use base R for safety if dplyr not loaded (it isn't)
        # Split by quadrant, sort by distance, take top 5
        top_genes_indices <- c()
        for (q in unique(sig_human$Quadrant)) {
             sub <- sig_human[sig_human$Quadrant == q, ]
             sub <- sub[order(sub$dist, decreasing=TRUE), ]
             top_genes_indices <- c(top_genes_indices, rownames(sub)[1:min(5, nrow(sub))])
        }
        label_data <- sig_human[top_genes_indices, ]
        
        p_quad <- p_quad + geom_text_repel(data=label_data, aes(label=human_symbol), 
                                           size=3, fontface="bold", color="black", max.overlaps=20)
            
        save_pdf(p_quad, file.path(RESULTS_DIR, "quadrant_scatter_mcd.pdf"), width=5, height=6)
    } else {
        message("No significant human genes found in overlap.")
    }
}

# --- 3. Pathway NES Scatter ---
message("Generating Pathway Concordance Scatter...")

h_fgsea_file <- file.path(CONC_DIR, "fgsea_human_results.csv")
m_fgsea_file <- file.path(CONC_DIR, "fgsea_mouse_results.csv")

if (file.exists(h_fgsea_file) && file.exists(m_fgsea_file)) {
    h_fgsea <- fread(h_fgsea_file)
    m_fgsea <- fread(m_fgsea_file)
    setDF(h_fgsea)
    setDF(m_fgsea)
    
    # Filter Human for 'disease_vs_ctrl'
    h_sub <- h_fgsea[h_fgsea$source == "disease_vs_ctrl", ]
    
    # Filter Mouse for 'MCD' (Representative)
    # Could also loop through diets
    m_sub <- m_fgsea[m_fgsea$source == "MCD", ]
    
    # Merge
    path_df <- merge(h_sub[, c("pathway", "NES", "padj")], 
                     m_sub[, c("pathway", "NES", "padj")], 
                     by="pathway", suffixes=c("_human", "_mouse"))
    
    # Clean names
    path_df$pathway_clean <- gsub("HALLMARK_", "", path_df$pathway)
    path_df$pathway_clean <- gsub("_", " ", path_df$pathway_clean)
    
    # Calculate Correlation
    cor_val <- cor(path_df$NES_human, path_df$NES_mouse, method="pearson", use="complete.obs")
    
    p_path <- ggplot(path_df, aes(x=NES_human, y=NES_mouse, label=pathway_clean)) +
        geom_point(color=sanjana_colors[["Purple"]], size=2, alpha=0.8) +
        geom_text_repel(size=2, max.overlaps=15, box.padding=0.3) +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="grey") +
        geom_hline(yintercept=0, linetype="dotted", color="grey") +
        geom_vline(xintercept=0, linetype="dotted", color="grey") +
        theme_publication() +
        labs(title=paste0("Pathway Concordance (MCD vs Human)"), 
             subtitle=paste0("Pearson r = ", round(cor_val, 2)),
             x="Human NES", y="Mouse MCD NES")
            
    save_pdf(p_path, file.path(RESULTS_DIR, "pathway_scatter_mcd.pdf"), width=6, height=6)
}

message("Done.")
