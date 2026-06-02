#!/usr/bin/env Rscript
# 01_prepare_combined_metadata.R
# Merge metadata from in-house and external MCD datasets

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

root <- normalizePath(".")
out_dir <- file.path(root, "metadata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("=== Preparing Combined MCD Metadata ===")

# --- In-house MCD ---
inhouse_meta <- read_tsv(

  file.path(root, "..", "InHouse_MCD", "metadata", "samples.tsv"),
  col_types = cols()
) %>%
  mutate(
    batch = "InHouse",
    dataset_label = "In-House MCD",
    genotype = NA_character_,
    layout = "PAIRED",
    source_counts = file.path(root, "..", "InHouse_MCD", "counts", "featurecounts", "gene_counts.txt")
  ) %>%
  select(sample_id, diet, batch, dataset_label, genotype, layout, source_counts)

message("In-house samples: ", nrow(inhouse_meta))

# --- GSE156918 (from Public_MCD) ---
# Note: Only Cre (wild-type) samples were processed; FLCN-KO samples excluded
gse156918_meta <- read_tsv(
  file.path(root, "..", "Public_MCD", "metadata", "samples.tsv"),
  col_types = cols()
) %>%
  filter(grepl("GSM47481", sample_id)) %>%
  filter(genotype == "Cre") %>%
  mutate(
    batch = "GSE156918",
    dataset_label = "External MCD",
    layout = "SINGLE",
    source_counts = file.path(root, "..", "Public_MCD", "counts", "featurecounts", "gene_counts.txt")
  ) %>%
  select(sample_id, diet, batch, dataset_label, genotype, layout, source_counts)

message("GSE156918 samples (Cre-only): ", nrow(gse156918_meta))

# --- GSE205974 (from Public_MCD) ---
gse205974_meta <- read_tsv(
  file.path(root, "..", "Public_MCD", "metadata", "samples.tsv"),
  col_types = cols()
) %>% 
  filter(grepl("GSM623625", sample_id)) %>%
  mutate(
    batch = "GSE205974",
    dataset_label = "External MCD",
    layout = ifelse(grepl("GSM623625", sample_id), "PAIRED", "SINGLE"),
    source_counts = file.path(root, "..", "Public_MCD", "counts", "featurecounts", "gene_counts.txt")
  ) %>%
  select(sample_id, diet, batch, dataset_label, genotype, layout, source_counts)

message("GSE205974 samples: ", nrow(gse205974_meta))

# --- Combine all ---
combined <- bind_rows(inhouse_meta, gse156918_meta, gse205974_meta) %>%
  mutate(
    diet = factor(diet, levels = c("Control", "MCD")),
    batch = factor(batch),
    dataset_label = factor(dataset_label)
  )

message("\n=== Combined Dataset Summary ===")
message("Total samples: ", nrow(combined))
message("By diet: ", paste(table(combined$diet), collapse = " / "))
message("By batch: ")
print(table(combined$batch, combined$diet))
message("By label:")
print(table(combined$dataset_label, combined$diet))

# Write combined metadata
write_tsv(combined, file.path(out_dir, "samples.tsv"))
message("\nWrote combined metadata to: ", file.path(out_dir, "samples.tsv"))
