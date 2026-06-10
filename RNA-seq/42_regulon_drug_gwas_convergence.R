#!/usr/bin/env Rscript
# =============================================================================
# 42_regulon_drug_gwas_convergence.R
# Module A3: Regulon → Drug → GWAS Convergence Analysis
#
# Identifies TFs whose disease regulons converge on GWAS-supported, druggable
# targets — the "triple convergence" analysis.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# --- Configuration ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

atlas_path     <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
hep_reg_path   <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv")
dis_reg_path   <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv")
enh_path       <- file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/enhancer_gene_links.csv")

# COLOC results (gene-level)
coloc_alt_path <- file.path(BASE, "RNA-seq/results/causal_inference/broadaway_ukbb/coloc_results.csv")
coloc_ast_path <- file.path(BASE, "RNA-seq/results/causal_inference/broadaway_ukbb_ast/coloc_results.csv")
coloc_ggt_path <- file.path(BASE, "RNA-seq/results/causal_inference/broadaway_ukbb_ggt/coloc_results.csv")
coloc_pdff_path <- file.path(BASE, "RNA-seq/results/causal_inference/broadaway_pdff/coloc_results.csv")

# Drug data
dgidb_path     <- file.path(BASE, "RNA-seq/results/drug_repurposing/dgidb_drug_gene_interactions.csv")
lincs_path     <- file.path(BASE, "RNA-seq/results/drug_repurposing/lincs_annotated_compounds.csv")
conv_drug_path <- file.path(BASE, "RNA-seq/results/drug_repurposing/convergent_drug_targets.csv")

# ieQTL
ieqtl_path     <- file.path(BASE, "RNA-seq/results/causal_inference/sceqtl/ieqtl_disease_genes.csv")

out_dir <- file.path(BASE, "RNA-seq/results/convergence")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Module A3: Regulon-Drug-GWAS Convergence ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. Load All Data
# =============================================================================
cat("--- Loading data ---\n")

atlas <- fread(atlas_path)
cat("  Atlas: ", nrow(atlas), " genes\n")

# Disease regulons
dis_regs <- fread(dis_reg_path)
cat("  Disease regulons: ", nrow(dis_regs), " TFs\n")

# Hepatocyte regulon edges (TF → target)
hep_regs <- fread(hep_reg_path)
cat("  Hepatocyte regulon edges: ", nrow(hep_regs), "\n")

# Enhancer-gene links
enhancers <- fread(enh_path)
cat("  Enhancer-gene links: ", nrow(enhancers), "\n")

# COLOC results - load all 4 GWAS
load_coloc <- function(path, gwas_name) {
  if (!file.exists(path)) {
    cat(sprintf("  WARNING: %s not found\n", path))
    return(data.frame(gene = character(), PP.H4 = numeric(), gwas = character()))
  }
  df <- fread(path)
  df$gwas <- gwas_name
  df %>% select(gene, ensembl, PP.H4, gwas)
}

coloc_alt <- load_coloc(coloc_alt_path, "UKBB_ALT")
coloc_ast <- load_coloc(coloc_ast_path, "UKBB_AST")
coloc_ggt <- load_coloc(coloc_ggt_path, "UKBB_GGT")
coloc_pdff <- load_coloc(coloc_pdff_path, "PDFF")

coloc_all <- bind_rows(coloc_alt, coloc_ast, coloc_ggt, coloc_pdff)
coloc_sig <- coloc_all %>% filter(PP.H4 > 0.5)
cat(sprintf("  COLOC genes (PP4>0.5): ALT=%d, AST=%d, GGT=%d, PDFF=%d\n",
            sum(coloc_alt$PP.H4 > 0.5), sum(coloc_ast$PP.H4 > 0.5),
            sum(coloc_ggt$PP.H4 > 0.5), sum(coloc_pdff$PP.H4 > 0.5)))

# DGIdb drug-gene interactions
dgidb <- if (file.exists(dgidb_path)) {
  df <- fread(dgidb_path)
  cat("  DGIdb interactions: ", nrow(df), "\n")
  df
} else {
  cat("  DGIdb: not found\n")
  data.frame(gene = character())
}

# LINCS annotated compounds
lincs <- if (file.exists(lincs_path)) {
  df <- fread(lincs_path)
  cat("  LINCS compounds: ", nrow(df), "\n")
  df
} else {
  cat("  LINCS: not found\n")
  data.frame()
}

# Convergent drug targets (pre-computed)
conv_drugs <- if (file.exists(conv_drug_path)) {
  df <- fread(conv_drug_path)
  cat("  Convergent drug targets: ", nrow(df), "\n")
  df
} else {
  cat("  Convergent drugs: not found\n")
  data.frame(symbol = character())
}

# ieQTL disease interactions
ieqtl <- if (file.exists(ieqtl_path)) {
  df <- fread(ieqtl_path)
  cat("  ieQTL genes: ", nrow(df), "\n")
  df
} else {
  cat("  ieQTL: not found\n")
  data.frame(gene = character())
}

# =============================================================================
# 2. Build Disease Regulon → Target Gene Map
# =============================================================================
cat("\n--- Building regulon target maps ---\n")

disease_tfs <- dis_regs$tf_name

# Disease regulon edges
disease_edges <- hep_regs %>%
  filter(tf_name %in% disease_tfs)

# Collect all unique target genes from disease regulons
all_regulon_targets <- unique(disease_edges$target_gene)
cat(sprintf("  Disease regulon target genes: %d (from %d TFs)\n",
            length(all_regulon_targets), length(disease_tfs)))

# Also include TFs themselves as they self-regulate
all_regulon_genes <- unique(c(all_regulon_targets, disease_tfs))
cat(sprintf("  All regulon-associated genes (targets + TFs): %d\n",
            length(all_regulon_genes)))

# =============================================================================
# 3. Layer 1 — GWAS Support for Regulon Genes
# =============================================================================
cat("\n--- Layer 1: GWAS support ---\n")

# Check which regulon genes have COLOC support
gwas_support <- coloc_sig %>%
  filter(gene %in% all_regulon_genes) %>%
  group_by(gene) %>%
  summarise(
    gwas_coloc = paste(gwas, collapse = ";"),
    n_gwas_coloc = n_distinct(gwas),
    max_pp4 = max(PP.H4),
    .groups = "drop"
  )

cat(sprintf("  Regulon genes with COLOC (PP4>0.5): %d / %d\n",
            nrow(gwas_support), length(all_regulon_genes)))

# Also check ieQTL disease interaction
ieqtl_support <- if (nrow(ieqtl) > 0) {
  ieqtl %>%
    filter(gene %in% all_regulon_genes) %>%
    group_by(gene) %>%
    summarise(
      ieqtl_cell_types = paste(unique(cell_type), collapse = ";"),
      min_ieqtl_pval = min(interaction_pval, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  data.frame(gene = character())
}

cat(sprintf("  Regulon genes with ieQTL: %d\n", nrow(ieqtl_support)))

# Check TWAS from atlas (MR column removed 2026-04-22 — MR ditched from paper)
mr_twas_support <- atlas %>%
  filter(human_symbol %in% all_regulon_genes) %>%
  filter(!is.na(twas_pval) & twas_pval < 0.05) %>%
  select(human_symbol, twas_pval) %>%
  rename(gene = human_symbol)

cat(sprintf("  Regulon genes with TWAS support: %d\n", nrow(mr_twas_support)))

# =============================================================================
# 4. Layer 2 — Druggability for Regulon Genes
# =============================================================================
cat("\n--- Layer 2: Druggability ---\n")

# DGIdb druggable genes
dgidb_druggable <- if (nrow(dgidb) > 0) {
  dgidb %>%
    filter(gene %in% all_regulon_genes) %>%
    group_by(gene) %>%
    summarise(
      n_drugs = n_distinct(drug_name),
      n_approved = sum(approved == TRUE, na.rm = TRUE),
      drugs = paste(unique(drug_name), collapse = ";"),
      .groups = "drop"
    )
} else {
  data.frame(gene = character())
}

cat(sprintf("  Regulon genes in DGIdb: %d\n", nrow(dgidb_druggable)))

# LINCS reversal from atlas
lincs_reversal <- atlas %>%
  filter(human_symbol %in% all_regulon_genes & !is.na(lincs_reversal) & lincs_reversal == TRUE) %>%
  select(human_symbol) %>%
  rename(gene = human_symbol) %>%
  mutate(lincs_hit = TRUE)

cat(sprintf("  Regulon genes with LINCS reversal: %d\n", nrow(lincs_reversal)))

# OpenTargets drugs from atlas
ot_drugs <- atlas %>%
  filter(human_symbol %in% all_regulon_genes & !is.na(opentargets_drug) & opentargets_drug == TRUE) %>%
  select(human_symbol) %>%
  rename(gene = human_symbol) %>%
  mutate(opentargets_hit = TRUE)

cat(sprintf("  Regulon genes with OpenTargets drugs: %d\n", nrow(ot_drugs)))

# =============================================================================
# 5. Layer 3 — Clinical Pipeline Cross-Reference
# =============================================================================
cat("\n--- Layer 3: Clinical pipeline ---\n")

clinical_drugs <- data.frame(
  drug = c("Resmetirom", "Semaglutide", "Efruxifermin", "Survodutide",
           "ION224", "Rapirosiran", "Obeticholic acid", "Elafibranor",
           "Lanifibranor", "Aramchol"),
  target = c("THRB", "GLP1R", "FGFR1;KLB", "GLP1R;GCGR", "DGAT2",
             "HSD17B13", "NR1H4", "PPARA;PPARD", "PPARA;PPARD;PPARG", "SCD"),
  stringsAsFactors = FALSE
) %>%
  mutate(target_list = strsplit(target, ";")) %>%
  unnest(target_list) %>%
  rename(gene = target_list)

clinical_match <- clinical_drugs %>%
  filter(gene %in% all_regulon_genes)

cat(sprintf("  Regulon genes matching clinical drug targets: %d\n", nrow(clinical_match)))
if (nrow(clinical_match) > 0) {
  print(clinical_match %>% select(gene, drug), row.names = FALSE)
}

# =============================================================================
# 6. Build Convergence Table
# =============================================================================
cat("\n--- Building convergence table ---\n")

# Start with all disease regulon edges
convergence <- disease_edges %>%
  select(tf_name, target_gene, regulon_activity_diff, activity_padj) %>%
  rename(gene = target_gene)

# Also add TFs themselves
tf_self <- data.frame(
  tf_name = disease_tfs,
  gene = disease_tfs,
  stringsAsFactors = FALSE
) %>%
  left_join(dis_regs %>% select(tf_name, regulon_activity_diff, activity_padj), by = "tf_name")

convergence <- bind_rows(convergence, tf_self) %>% distinct()

# Join atlas evidence
atlas_slim <- atlas %>%
  select(human_symbol, bulk_logFC, bulk_padj, bulk_tstat,
         is_conserved, attribution_class) %>%
  rename(gene = human_symbol)

convergence <- convergence %>%
  left_join(atlas_slim, by = "gene")

# Join GWAS support
convergence <- convergence %>%
  left_join(gwas_support, by = "gene") %>%
  mutate(has_gwas = !is.na(gwas_coloc))

# Join ieQTL
convergence <- convergence %>%
  left_join(ieqtl_support, by = "gene") %>%
  mutate(has_ieqtl = !is.na(ieqtl_cell_types))

# Join TWAS (MR ditched 2026-04-22; has_mr_twas retained as TWAS-only flag)
convergence <- convergence %>%
  left_join(mr_twas_support, by = "gene", suffix = c("", "_mr")) %>%
  mutate(has_mr_twas = !is.na(twas_pval))

# Join druggability
convergence <- convergence %>%
  left_join(dgidb_druggable, by = "gene") %>%
  left_join(lincs_reversal, by = "gene") %>%
  left_join(ot_drugs, by = "gene") %>%
  mutate(
    has_dgidb = !is.na(n_drugs),
    has_lincs = !is.na(lincs_hit) & lincs_hit == TRUE,
    has_opentargets = !is.na(opentargets_hit) & opentargets_hit == TRUE,
    is_druggable = has_dgidb | has_lincs | has_opentargets
  )

# Join clinical drugs
clinical_lookup <- clinical_match %>%
  group_by(gene) %>%
  summarise(clinical_drug = paste(drug, collapse = ";"), .groups = "drop")

convergence <- convergence %>%
  left_join(clinical_lookup, by = "gene") %>%
  mutate(has_clinical = !is.na(clinical_drug))

# Compute convergence flags
convergence <- convergence %>%
  mutate(
    is_deg = !is.na(bulk_padj) & bulk_padj < 0.1,  # Exploratory annotation threshold; primary DEGs: padj<0.05 + |logFC|>0.5 (Script 05b)
    has_genetic = has_gwas | has_ieqtl | has_mr_twas,
    # Triple convergence: TF-regulated + GWAS + druggable
    is_triple = has_genetic & is_druggable,
    # Double convergence
    is_double_genetic_drug = has_genetic & is_druggable & !is_deg,
    is_double_genetic_deg = has_genetic & is_deg,
    is_double_drug_deg = is_druggable & is_deg,
    # Convergence score (count of evidence layers)
    convergence_layers = rowSums(cbind(
      is_deg, has_gwas, has_ieqtl, has_mr_twas,
      has_dgidb, has_lincs, has_opentargets, has_clinical
    ), na.rm = TRUE)
  )

cat("\nConvergence summary:\n")
cat(sprintf("  Total regulon gene entries: %d (unique genes: %d)\n",
            nrow(convergence), n_distinct(convergence$gene)))
cat(sprintf("  With GWAS support: %d\n", sum(convergence$has_gwas)))
cat(sprintf("  With ieQTL: %d\n", sum(convergence$has_ieqtl)))
cat(sprintf("  Druggable: %d\n", sum(convergence$is_druggable)))
cat(sprintf("  Triple convergence (TF × Genetic × Drug): %d\n", sum(convergence$is_triple)))
cat(sprintf("  With clinical drug: %d\n", sum(convergence$has_clinical)))

# =============================================================================
# 7. Triple Convergence Targets
# =============================================================================
cat("\n--- Triple convergence targets ---\n")

triple_targets <- convergence %>%
  filter(is_triple) %>%
  select(tf_name, gene, regulon_activity_diff, activity_padj,
         bulk_logFC, bulk_padj, is_deg,
         gwas_coloc, n_gwas_coloc, max_pp4,
         ieqtl_cell_types,
         n_drugs, drugs, clinical_drug,
         lincs_hit, opentargets_hit,
         convergence_layers) %>%
  arrange(desc(convergence_layers), activity_padj)

if (nrow(triple_targets) > 0) {
  cat("Triple convergence targets:\n")
  print(triple_targets %>%
    select(tf_name, gene, gwas_coloc, n_drugs, clinical_drug, convergence_layers),
    row.names = FALSE)
}

# =============================================================================
# 8. Convergence Network (Edge List)
# =============================================================================
cat("\n--- Building convergence network ---\n")

edges <- list()

# TF → gene regulatory edges
reg_edges <- convergence %>%
  filter(tf_name != gene) %>%
  select(from = tf_name, to = gene) %>%
  mutate(edge_type = "regulatory", weight = 1)
edges[["regulatory"]] <- reg_edges

# Drug → gene pharmacological edges
if (nrow(dgidb_druggable) > 0) {
  drug_genes <- dgidb %>%
    filter(gene %in% all_regulon_genes) %>%
    select(from = drug_name, to = gene) %>%
    mutate(edge_type = "pharmacological", weight = 1)
  edges[["pharmacological"]] <- drug_genes
}

# GWAS → gene genetic edges
if (nrow(gwas_support) > 0) {
  gwas_edges <- coloc_sig %>%
    filter(gene %in% all_regulon_genes) %>%
    select(from = gwas, to = gene) %>%
    mutate(edge_type = "genetic", weight = 1)
  edges[["genetic"]] <- gwas_edges
}

network_edges <- bind_rows(edges)
cat(sprintf("  Network: %d edges (%d regulatory, %d pharmacological, %d genetic)\n",
            nrow(network_edges),
            sum(network_edges$edge_type == "regulatory"),
            sum(network_edges$edge_type == "pharmacological"),
            sum(network_edges$edge_type == "genetic")))

# =============================================================================
# 9. Permutation Test
# =============================================================================
cat("\n--- Permutation testing ---\n")

# Test: are disease regulon targets more likely to be GWAS-supported than random?
n_perm <- 10000
set.seed(42)

observed_gwas <- sum(all_regulon_genes %in% gwas_support$gene)
observed_drug <- sum(all_regulon_genes %in% dgidb_druggable$gene)
observed_triple <- sum(convergence$is_triple)

# Universe: all genes in atlas
all_genes <- atlas$human_symbol
all_coloc_genes <- unique(coloc_sig$gene)
all_dgidb_genes <- if (nrow(dgidb) > 0) unique(dgidb$gene) else character()

perm_gwas <- numeric(n_perm)
perm_drug <- numeric(n_perm)
perm_triple <- numeric(n_perm)

n_regulon <- length(all_regulon_genes)

for (i in seq_len(n_perm)) {
  random_genes <- sample(all_genes, n_regulon)
  perm_gwas[i] <- sum(random_genes %in% all_coloc_genes)
  perm_drug[i] <- sum(random_genes %in% all_dgidb_genes)
  perm_triple[i] <- sum(random_genes %in% all_coloc_genes &
                         random_genes %in% all_dgidb_genes)
}

p_gwas <- (sum(perm_gwas >= observed_gwas) + 1) / (n_perm + 1)
p_drug <- (sum(perm_drug >= observed_drug) + 1) / (n_perm + 1)
p_triple <- (sum(perm_triple >= observed_triple) + 1) / (n_perm + 1)

perm_stats <- data.frame(
  test = c("GWAS_enrichment", "Drug_enrichment", "Triple_convergence"),
  observed = c(observed_gwas, observed_drug, observed_triple),
  expected_mean = c(mean(perm_gwas), mean(perm_drug), mean(perm_triple)),
  expected_sd = c(sd(perm_gwas), sd(perm_drug), sd(perm_triple)),
  fold_enrichment = c(
    observed_gwas / max(mean(perm_gwas), 0.001),
    observed_drug / max(mean(perm_drug), 0.001),
    observed_triple / max(mean(perm_triple), 0.001)
  ),
  perm_pvalue = c(p_gwas, p_drug, p_triple),
  n_permutations = n_perm,
  stringsAsFactors = FALSE
)

cat("\nPermutation results:\n")
print(perm_stats, row.names = FALSE)

# =============================================================================
# 10. Per-TF Convergence Summary
# =============================================================================
cat("\n--- Per-TF convergence summary ---\n")

tf_summary <- convergence %>%
  group_by(tf_name) %>%
  summarise(
    n_targets = n_distinct(gene),
    n_deg_targets = sum(is_deg, na.rm = TRUE),
    n_gwas_targets = sum(has_gwas, na.rm = TRUE),
    n_druggable_targets = sum(is_druggable, na.rm = TRUE),
    n_triple = sum(is_triple, na.rm = TRUE),
    n_clinical = sum(has_clinical, na.rm = TRUE),
    max_convergence = max(convergence_layers, na.rm = TRUE),
    gwas_genes = paste(gene[has_gwas], collapse = ";"),
    drug_genes = paste(gene[is_druggable], collapse = ";"),
    .groups = "drop"
  ) %>%
  left_join(dis_regs %>% select(tf_name, regulon_activity_diff, activity_padj),
            by = "tf_name") %>%
  arrange(desc(n_triple), desc(n_gwas_targets))

cat("\nTop converging TFs:\n")
print(tf_summary %>%
  select(tf_name, n_targets, n_gwas_targets, n_druggable_targets, n_triple,
         regulon_activity_diff) %>%
  head(15), row.names = FALSE)

# =============================================================================
# 11. Write Outputs
# =============================================================================
cat("\n--- Writing outputs ---\n")

fwrite(convergence,
       file.path(out_dir, "regulon_drug_gwas_convergence.csv"))
cat("  Written: regulon_drug_gwas_convergence.csv\n")

fwrite(triple_targets,
       file.path(out_dir, "triple_convergence_targets.csv"))
cat("  Written: triple_convergence_targets.csv\n")

fwrite(network_edges,
       file.path(out_dir, "convergence_network.csv"))
cat("  Written: convergence_network.csv\n")

fwrite(perm_stats,
       file.path(out_dir, "convergence_permutation_stats.csv"))
cat("  Written: convergence_permutation_stats.csv\n")

fwrite(tf_summary,
       file.path(out_dir, "tf_convergence_summary.csv"))
cat("  Written: tf_convergence_summary.csv\n")

cat("\n=== Module A3 Complete ===\n")
cat("End:", format(Sys.time()), "\n")
