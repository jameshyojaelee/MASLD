#!/usr/bin/env Rscript
# meta_05_library_comparison.R
# Compare combined MCD DEGs with final_core_degs.csv library

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

root <- normalizePath(".")
out_dir <- file.path(root, "meta_analysis")
plot_dir <- file.path(root, "plots")

message("=== Library Comparison: Combined MCD vs final_core_degs.csv ===\n")

# Load gene sets
gene_sets <- readRDS(file.path(out_dir, "upregulated_gene_sets.rds"))
all_degs <- readRDS(file.path(out_dir, "all_degs_list.rds"))

# Load library
library_path <- file.path(root, "..", "..", "final_core_degs.csv")
if (!file.exists(library_path)) {
  library_path <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/results/library/final_core_degs.csv"
}
library_df <- read_csv(library_path, col_types = cols())

message("Library genes: ", nrow(library_df))

# Parse source analyses to identify MCD vs Patient evidence
library_df <- library_df %>%
  mutate(
    has_mcd = str_detect(source_analyses, "MCD"),
    has_inhouse_mcd = str_detect(source_analyses, "MCD \\(in-house\\)"),
    has_external_mcd = str_detect(source_analyses, "MCD \\(external\\)"),
    has_patient = str_detect(source_analyses, "Patient"),
    has_gwas = str_detect(source_analyses, "GWAS")
  )

message("Library with MCD evidence: ", sum(library_df$has_mcd, na.rm = TRUE))
message("Library with Patient evidence: ", sum(library_df$has_patient, na.rm = TRUE))

# Get library gene IDs
library_genes <- library_df$mouse_gene_id

# Compare with combined upregulated DEGs
combined_lenient <- gene_sets$Combined_lenient
combined_stringent <- gene_sets$Combined_stringent

# Coverage analysis
in_combined_lenient <- intersect(library_genes, combined_lenient)
in_combined_stringent <- intersect(library_genes, combined_stringent)

message("\n=== Library Coverage by Combined Analysis ===")
message("Library genes in Combined (lenient): ", length(in_combined_lenient), 
        " (", round(100 * length(in_combined_lenient) / length(library_genes), 1), "%)")
message("Library genes in Combined (stringent): ", length(in_combined_stringent),
        " (", round(100 * length(in_combined_stringent) / length(library_genes), 1), "%)")

# Also check individual datasets
inhouse_lenient <- gene_sets$InHouse_lenient
gse156918_lenient <- gene_sets$GSE156918_lenient
gse205974_lenient <- gene_sets$GSE205974_lenient

coverage_df <- tibble(
  dataset = c("Combined", "InHouse", "GSE156918", "GSE205974"),
  lenient_coverage = c(
    length(intersect(library_genes, combined_lenient)),
    length(intersect(library_genes, inhouse_lenient)),
    length(intersect(library_genes, gse156918_lenient)),
    length(intersect(library_genes, gse205974_lenient))
  ),
  stringent_coverage = c(
    length(intersect(library_genes, combined_stringent)),
    length(intersect(library_genes, gene_sets$InHouse_stringent)),
    length(intersect(library_genes, gene_sets$GSE156918_stringent)),
    length(intersect(library_genes, gene_sets$GSE205974_stringent))
  )
) %>%
  mutate(
    lenient_pct = round(100 * lenient_coverage / nrow(library_df), 1),
    stringent_pct = round(100 * stringent_coverage / nrow(library_df), 1)
  )

write_csv(coverage_df, file.path(out_dir, "library_coverage.csv"))
message("\n=== Coverage Summary ===")
print(coverage_df)

# Which library genes are NOT in combined but ARE in individuals?
individual_union <- unique(c(inhouse_lenient, gse156918_lenient, gse205974_lenient))
library_in_individual_not_combined <- setdiff(
  intersect(library_genes, individual_union),
  combined_lenient
)
message("\nLibrary genes in individuals but NOT in combined: ", length(library_in_individual_not_combined))

# Which combined genes are NOT in library? (potential additions)
combined_not_in_library <- setdiff(combined_stringent, library_genes)

# Get gene info for potential additions
potential_additions <- all_degs$Combined %>%
  filter(gene_id %in% combined_not_in_library) %>%
  select(gene_id, gene_name, log2FoldChange, padj, baseMean) %>%
  arrange(padj) %>%
  head(500)

write_csv(potential_additions, file.path(out_dir, "potential_library_additions.csv"))
message("Potential library additions (stringent, not in library): ", 
        length(combined_not_in_library), " → saved top 500")

# Create bar chart comparing coverage
p <- ggplot(coverage_df, aes(x = reorder(dataset, -lenient_coverage), y = lenient_coverage)) +
  geom_bar(stat = "identity", fill = "#3498DB", alpha = 0.8) +
  geom_text(aes(label = paste0(lenient_pct, "%")), vjust = -0.5, size = 4) +
  labs(
    title = "Library Gene Coverage by MCD Dataset",
    subtitle = paste0("Out of ", nrow(library_df), " library genes (lenient threshold: padj<0.1, LFC>0)"),
    x = "Dataset",
    y = "Library Genes Detected as Upregulated"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "library_coverage_bar.pdf"), p, width = 7, height = 6)
message("\nSaved: ", file.path(plot_dir, "library_coverage_bar.pdf"))

# Create Venn-style summary
message("\n=== Library Gene Summary ===")
message("Total library genes: ", nrow(library_df))
message("Detected in Combined (lenient): ", length(in_combined_lenient))
message("Detected in Combined (stringent): ", length(in_combined_stringent))
message("NOT detected in Combined (lenient): ", length(setdiff(library_genes, combined_lenient)))
