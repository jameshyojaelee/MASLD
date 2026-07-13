#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_250_nodes
#SBATCH --output=logs/net_250_nodes_%j.out
#SBATCH --error=logs/net_250_nodes_%j.err
# 250_define_network_nodes.R — Define gene node set V for the Bayesian multiplex gene network.
# Filters the multi-evidence atlas to genes with meaningful evidence (not just tested).

suppressPackageStartupMessages(library(data.table))

t0 <- Sys.time()

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS_PATH <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
OUTDIR <- file.path(BASE, "RNA-seq/results/network")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 250: Define Network Node Set V ===\n")
cat("Start:", format(t0, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Atlas:", ATLAS_PATH, "\n\n")

atlas <- fread(ATLAS_PATH)
cat("Atlas loaded:", nrow(atlas), "genes x", ncol(atlas), "columns\n\n")
stopifnot(all(c("bulk_padj","bulk_logFC") %in% names(atlas)))

has_col <- function(x) x %in% names(atlas)

crit_deg        <- !is.na(atlas$bulk_padj) & atlas$bulk_padj < 0.05
crit_coloc      <- !is.na(atlas$coloc_susie_best_pp4) & atlas$coloc_susie_best_pp4 > 0.5
crit_cc         <- !is.na(atlas$is_conserved) & atlas$is_conserved == TRUE
crit_essential  <- !is.na(atlas$essentiality_chronos) & atlas$essentiality_chronos < -0.5
crit_liana      <- !is.na(atlas$liana_n_diff_interactions) & atlas$liana_n_diff_interactions > 0
crit_spatial    <- if (has_col("spatial_sig")) !is.na(atlas$spatial_sig) & atlas$spatial_sig == TRUE else rep(FALSE, nrow(atlas))
crit_gwas_atac  <- if (has_col("gwas_variant_in_peak")) !is.na(atlas$gwas_variant_in_peak) & atlas$gwas_variant_in_peak == TRUE else rep(FALSE, nrow(atlas))
# crit_mr dropped 2026-04-22 — MR ditched from paper; TWAS+COLOC+INTACT only.
crit_mouse      <- !is.na(atlas$mouse_meta_padj) & atlas$mouse_meta_padj < 0.05
crit_sceqtl     <- !is.na(atlas$sceqtl_coloc_best_pp4) & atlas$sceqtl_coloc_best_pp4 > 0.5

any_criterion <- crit_deg | crit_coloc | crit_cc | crit_essential | crit_liana |
                 crit_spatial | crit_gwas_atac | crit_mouse | crit_sceqtl

nodes <- atlas[any_criterion]

nodes[, is_deg := !is.na(bulk_padj) & bulk_padj < 0.05 & !is.na(bulk_logFC) & abs(bulk_logFC) > 0.3]

keep_cols <- intersect(
  c("human_symbol", "ensembl_id", "gene_biotype", "bulk_logFC", "bulk_padj",
    "bulk_tstat", "is_deg", "coloc_susie_best_pp4", "is_conserved",
    "sex_class", "zonation_class", "ferroptosis_class", "dgidb_druggable",
    "n_coloc_sources", "attribution_class", "mouse_ortholog"),
  c(names(nodes), "is_deg")
)
nodes <- nodes[, .SD, .SDcols = keep_cols]

fwrite(nodes, file.path(OUTDIR, "network_nodes.csv"))

cat("--- Inclusion Criteria Breakdown ---\n")
criteria <- data.table(
  criterion = c("bulk_padj < 0.05", "coloc_susie_best_pp4 > 0.5",
                "is_conserved", "essentiality_chronos < -0.5",
                "liana_n_diff_interactions > 0", "spatial_sig",
                "gwas_variant_in_peak",
                "mouse_meta_padj < 0.05", "sceqtl_coloc_best_pp4 > 0.5"),
  # mr_sig row dropped 2026-04-22 — MR ditched from paper.
  n = c(sum(crit_deg), sum(crit_coloc), sum(crit_cc), sum(crit_essential),
        sum(crit_liana), sum(crit_spatial), sum(crit_gwas_atac),
        sum(crit_mouse), sum(crit_sceqtl))
)
criteria[, pct := sprintf("%.1f%%", 100 * n / nrow(atlas))]
print(criteria, row.names = FALSE)

cat("\n--- Summary ---\n")
cat("Atlas genes:", nrow(atlas), "\n")
cat("Node set V:", nrow(nodes), sprintf("(%.1f%% of atlas)\n", 100 * nrow(nodes) / nrow(atlas)))
cat("Output:", file.path(OUTDIR, "network_nodes.csv"), "\n")

t1 <- Sys.time()
cat("\nElapsed:", format(round(difftime(t1, t0, units = "secs"), 1)), "\n")
cat("Done.\n")
