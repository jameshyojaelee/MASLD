#!/usr/bin/env Rscript

# Script to analyze gene biotypes in the GWAS closest genes file
# This will categorize genes as protein-coding, lncRNA, pseudogene, etc.

# Load required libraries
library(tidyverse)
library(biomaRt)
library(ggplot2)

# Define file paths
input_file <- "/gpfs/commons/home/jameslee/Cas13/GWAS/Closest_genes.csv"
output_dir <- "/gpfs/commons/home/jameslee/Cas13/GWAS"

# Read the input file
cat("Reading input file:", input_file, "\n")
gene_list <- read.csv(input_file, header = FALSE, stringsAsFactors = FALSE)
colnames(gene_list) <- c("gene_symbol")
cat("Total genes in list:", nrow(gene_list), "\n")

# Clean gene symbols by removing everything after the pipe symbol (|)
gene_list$gene_symbol_clean <- sapply(strsplit(as.character(gene_list$gene_symbol), "\\|"), `[`, 1)
cat("After cleaning gene names:", length(unique(gene_list$gene_symbol_clean)), "unique gene symbols\n")

# Connect to Ensembl BioMart database
cat("Connecting to Ensembl BioMart...\n")
ensembl <- tryCatch({
  useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
}, error = function(e) {
  cat("Error connecting to Ensembl:", e$message, "\n")
  cat("Trying alternative Ensembl mirror...\n")
  # Try an alternative mirror
  tryCatch({
    useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
  }, error = function(e2) {
    cat("Error connecting to alternative Ensembl mirror:", e2$message, "\n")
    return(NULL)
  })
})

if (!is.null(ensembl)) {
  cat("Retrieving gene biotype information...\n")
  
  # Define which attributes to retrieve
  attributes <- c("ensembl_gene_id", "external_gene_name", "gene_biotype", "chromosome_name", "description")
  
  # Query biomaRt for the gene symbols
  gene_info <- tryCatch({
    getBM(
      attributes = attributes,
      filters = "external_gene_name",
      values = unique(gene_list$gene_symbol_clean),
      mart = ensembl
    )
  }, error = function(e) {
    cat("Error querying Ensembl:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(gene_info) && nrow(gene_info) > 0) {
    cat("Successfully retrieved biotype information for", nrow(gene_info), "genes\n")
    
    # Check for missing genes
    found_genes <- gene_info$external_gene_name
    missing_genes <- setdiff(unique(gene_list$gene_symbol_clean), found_genes)
    
    if (length(missing_genes) > 0) {
      cat("Warning:", length(missing_genes), "genes not found in Ensembl. These might be:\n")
      cat("- Deprecated gene symbols\n")
      cat("- Ensembl IDs rather than gene symbols\n")
      cat("- Non-standard identifiers\n")
      
      # If there are genes starting with ENSG, they might be Ensembl IDs
      ensg_ids <- missing_genes[grepl("^ENSG", missing_genes)]
      if (length(ensg_ids) > 0) {
        cat("Found", length(ensg_ids), "potential Ensembl IDs, querying them directly...\n")
        
        # Query biomaRt for the Ensembl IDs
        ensg_info <- tryCatch({
          getBM(
            attributes = attributes,
            filters = "ensembl_gene_id",
            values = ensg_ids,
            mart = ensembl
          )
        }, error = function(e) {
          cat("Error querying Ensembl with IDs:", e$message, "\n")
          return(NULL)
        })
        
        if (!is.null(ensg_info) && nrow(ensg_info) > 0) {
          cat("Retrieved information for", nrow(ensg_info), "Ensembl IDs\n")
          # Add this information to the main gene_info dataframe
          gene_info <- rbind(gene_info, ensg_info)
        }
      }
      
      # Write missing genes to file
      write.csv(
        data.frame(gene_symbol = missing_genes), 
        file.path(output_dir, "missing_genes.csv"), 
        row.names = FALSE
      )
    }
    
    # Count the number of genes in each biotype category
    biotype_counts <- gene_info %>%
      group_by(gene_biotype) %>%
      summarize(
        count = n(),
        percent = round(n() / nrow(gene_info) * 100, 2),
        .groups = "drop"
      ) %>%
      arrange(desc(count))
    
    # Print summary
    cat("\nBiotype summary:\n")
    print(biotype_counts)
    
    # Create a lookup table for gene symbols to join with original list
    gene_lookup <- gene_info %>%
      dplyr::select(gene_symbol = external_gene_name, ensembl_id = ensembl_gene_id, biotype = gene_biotype, chromosome = chromosome_name)
    
    # Join with original gene list
    gene_list_with_biotype <- gene_list %>%
      left_join(gene_lookup, by = c("gene_symbol_clean" = "gene_symbol"))
    
    # For genes with missing biotypes, try to find them by Ensembl ID if they match the pattern
    missing_info_genes <- gene_list_with_biotype %>% filter(is.na(biotype))
    for (i in 1:nrow(missing_info_genes)) {
      if (grepl("^ENSG", missing_info_genes$gene_symbol[i])) {
        # Extract the ENSG ID before any pipe symbol
        ensg_id <- gsub("\\|.*$", "", missing_info_genes$gene_symbol[i])
        # Look up in gene_info
        matching_row <- gene_info %>% filter(ensembl_gene_id == ensg_id)
        if (nrow(matching_row) > 0) {
          # Update the main dataframe
          idx <- which(gene_list_with_biotype$gene_symbol == missing_info_genes$gene_symbol[i])
          gene_list_with_biotype$ensembl_id[idx] <- matching_row$ensembl_gene_id[1]
          gene_list_with_biotype$biotype[idx] <- matching_row$gene_biotype[1]
          gene_list_with_biotype$chromosome[idx] <- matching_row$chromosome_name[1]
        }
      }
    }
    
    # Count genes with and without biotype information
    cat("Genes with biotype information:", sum(!is.na(gene_list_with_biotype$biotype)), "/", nrow(gene_list_with_biotype), "\n")
    
    # Write detailed results to file
    detailed_output_file <- file.path(output_dir, "gwas_gene_biotype_analysis.csv")
    cat("Writing detailed results to:", detailed_output_file, "\n")
    write.csv(gene_list_with_biotype, detailed_output_file, row.names = FALSE)
    
    # Write biotype summary to file
    summary_file <- file.path(output_dir, "gwas_gene_biotype_summary.csv")
    cat("Writing biotype summary to:", summary_file, "\n")
    write.csv(biotype_counts, summary_file, row.names = FALSE)
    
    # Create pie chart of biotypes
    pdf(file.path(output_dir, "gwas_gene_biotype_piechart.pdf"), width = 10, height = 8)
    # Combine small categories into "Other" for better visualization
    biotype_for_plot <- biotype_counts %>%
      mutate(
        biotype_group = case_when(
          count < 5 ~ "Other",
          TRUE ~ gene_biotype
        )
      ) %>%
      group_by(biotype_group) %>%
      summarize(
        count = sum(count),
        percent = sum(percent),
        .groups = "drop"
      )
    
    # Create custom color palette for primary biotypes
    biotype_colors <- c(
      "protein_coding" = "#4285F4",         # Blue
      "lncRNA" = "#EA4335",                 # Red
      "processed_pseudogene" = "#FBBC05",   # Yellow
      "antisense" = "#34A853",              # Green
      "TEC" = "#8F9BBA",                    # Lavender
      "sense_intronic" = "#F6BDC0",         # Light pink
      "transcribed_unprocessed_pseudogene" = "#D5A6BD", # Light purple
      "Other" = "#CCCCCC"                  # Grey
    )
    
    # Assign colors to other biotypes not in our predefined set
    for (bt in unique(biotype_for_plot$biotype_group)) {
      if (!(bt %in% names(biotype_colors))) {
        biotype_colors[bt] <- "#AAAAAA"  # Default gray for other categories
      }
    }
    
    # Plot
    pie_data <- biotype_for_plot %>%
      arrange(desc(count)) %>%
      mutate(
        ypos = cumsum(percent) - 0.5*percent,
        label = paste0(biotype_group, "\n(", percent, "%)")
      )
    
    ggplot(pie_data, aes(x = "", y = percent, fill = biotype_group)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      theme_void() +
      scale_fill_manual(values = biotype_colors) +
      labs(
        title = "Distribution of Gene Biotypes in GWAS Hit Genes",
        fill = "Biotype"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")
      )
    dev.off()
    
    # Create barplot of top biotypes
    pdf(file.path(output_dir, "gwas_gene_biotype_barplot.pdf"), width = 12, height = 8)
    
    # Take top 10 biotypes for barplot
    top_biotypes <- biotype_counts %>%
      arrange(desc(count)) %>%
      head(10)
    
    ggplot(top_biotypes, aes(x = reorder(gene_biotype, count), y = count, fill = gene_biotype)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = biotype_colors) +
      labs(
        title = "Top Gene Biotypes in GWAS Hit Genes",
        x = "Biotype",
        y = "Count"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none"
      )
    dev.off()
    
    # Create chromosome distribution plot
    pdf(file.path(output_dir, "gwas_gene_chromosome_distribution.pdf"), width = 12, height = 8)
    
    # Count genes per chromosome
    chrom_count <- gene_list_with_biotype %>%
      filter(!is.na(chromosome)) %>%
      group_by(chromosome) %>%
      summarize(
        count = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(count))
    
    # Plot chromosome distribution
    ggplot(chrom_count, aes(x = reorder(chromosome, -count), y = count, fill = chromosome)) +
      geom_bar(stat = "identity") +
      labs(
        title = "Chromosome Distribution of GWAS Hit Genes",
        x = "Chromosome",
        y = "Number of Genes"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    dev.off()
    
    # Compare with general gene distribution
    # Get a count of all human genes by biotype for comparison
    cat("Retrieving general gene biotype distribution for comparison...\n")
    
    all_human_biotypes <- tryCatch({
      # Sample a subset of all human genes (full dataset would be too large)
      sample_genes <- getBM(
        attributes = c("ensembl_gene_id", "gene_biotype"),
        filters = "chromosome_name",
        values = c(1:22, "X", "Y"), # Standard chromosomes
        mart = ensembl,
        uniqueRows = TRUE
      )
      
      # Summarize biotype distribution
      sample_genes %>%
        group_by(gene_biotype) %>%
        summarize(
          count = n(),
          percent = round(n() / nrow(sample_genes) * 100, 2),
          .groups = "drop"
        ) %>%
        arrange(desc(count))
    }, error = function(e) {
      cat("Error retrieving general gene distribution:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(all_human_biotypes) && nrow(all_human_biotypes) > 0) {
      cat("Successfully retrieved general gene distribution\n")
      
      # Write to file
      write.csv(all_human_biotypes, file.path(output_dir, "human_gene_biotype_distribution.csv"), row.names = FALSE)
      
      # Create comparison barplot for top biotypes
      pdf(file.path(output_dir, "gwas_vs_general_biotype_comparison.pdf"), width = 12, height = 8)
      
      # Prepare data for comparison
      top_biotypes_general <- all_human_biotypes %>%
        arrange(desc(count)) %>%
        head(5) %>%
        mutate(source = "Human Genome")
      
      top_biotypes_gwas <- biotype_counts %>%
        arrange(desc(count)) %>%
        head(5) %>%
        rename(gene_biotype = gene_biotype) %>%
        mutate(source = "GWAS Hits")
      
      # Combine data
      comparison_data <- rbind(
        top_biotypes_general %>% dplyr::select(gene_biotype, percent, source),
        top_biotypes_gwas %>% dplyr::select(gene_biotype, percent, source)
      )
      
      # Plot
      ggplot(comparison_data, aes(x = gene_biotype, y = percent, fill = source)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(
          title = "GWAS Genes vs. Human Genome Biotype Distribution",
          x = "Biotype",
          y = "Percentage (%)",
          fill = "Source"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
      dev.off()
    }
    
  } else {
    cat("Failed to retrieve biotype information. Using manual classification.\n")
    
    # Attempt a manual classification based on gene name patterns
    # (This is a fallback and not very accurate)
    gene_list$biotype <- "unknown"
    
    # Try to guess based on gene symbol patterns (very rough)
    gene_list$biotype[grepl("^LINC", gene_list$gene_symbol_clean)] <- "lncRNA"
    gene_list$biotype[grepl("-AS", gene_list$gene_symbol_clean)] <- "antisense"
    gene_list$biotype[grepl("^MT-", gene_list$gene_symbol_clean)] <- "Mt_rRNA"
    
    # Count the number of genes in each biotype category
    biotype_counts <- gene_list %>%
      group_by(biotype) %>%
      summarize(
        count = n(),
        percent = round(n() / nrow(gene_list) * 100, 2),
        .groups = "drop"
      ) %>%
      arrange(desc(count))
    
    # Print summary
    cat("\nEstimated biotype summary (limited accuracy):\n")
    print(biotype_counts)
    
    # Write results to file
    output_file <- file.path(output_dir, "gwas_gene_biotype_analysis_estimated.csv")
    cat("Writing results with estimated biotypes to:", output_file, "\n")
    write.csv(gene_list, output_file, row.names = FALSE)
  }
} else {
  cat("Failed to connect to Ensembl. Cannot determine gene biotypes.\n")
}

cat("Process complete.\n") 