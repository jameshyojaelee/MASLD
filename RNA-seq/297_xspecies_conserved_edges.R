#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_297_dxs
#SBATCH --output=logs/net_297_dxs_%j.out
#SBATCH --error=logs/net_297_dxs_%j.err
# ===========================================================================
# Script 297: D-XS edges -- cross-species conserved pairs (Phase 3 rehab).
# ---------------------------------------------------------------------------
# Currently there is NO source for D-XS edges: both edges_d_xs.csv and the
# portal's "xspecies" layer are empty. This script builds two channels
# from the Cross_Species_Concordance pipeline outputs:
#
#   Channel 1: WGCNA module preservation.
#     For each WGCNA module with max_Zsummary (across diet models) >= 5
#     (well-preserved) OR >= 2 (moderately, as a fallback floor), emit all
#     gene pairs within that module. Annotate with module_id and
#     preservation score.
#
#   Channel 2: Cell-type conserved gene cliques.
#     From celltype_conserved_genes_all.csv, build a per-cell-type clique
#     of conserved genes (concordant LFC, both sig). Cross-cell-type pairs
#     excluded; within-cell-type pairs only. This prevents an N^2 explosion
#     across the full 1,108-Conserved set.
#
# Output: edges_d_xs_conserved.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
})

t0 <- Sys.time()

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

XS_DIR   <- file.path(PROJ, "Analysis/Cross_Species_Concordance/results")
PRES     <- file.path(XS_DIR, "wgcna_preservation_stats.csv")
MOD_ASSG <- file.path(XS_DIR, "wgcna_module_assignments.csv")
CT_CONS  <- file.path(XS_DIR, "celltype_conserved_genes_all.csv")
NODES    <- file.path(PROJ, "RNA-seq/results/network/network_nodes.csv")
OUT      <- file.path(PROJ, "RNA-seq/results/network/edges_d_xs_conserved.csv")

Z_STRICT   <- 5
Z_FALLBACK <- 2
MAX_MODULE_SIZE <- 300  # guard against a giant module producing N^2 edges

message("[297] Cross-species D-XS edge construction starting at ", format(t0))

stopifnot(file.exists(PRES), file.exists(MOD_ASSG), file.exists(CT_CONS),
          file.exists(NODES))

nodes <- fread(NODES); V <- unique(nodes$human_symbol)
message("  Node set V: ", length(V), " genes")

# ==========================================================================
#  Channel 1: WGCNA module pairs
# ==========================================================================
pres <- fread(PRES)
message("  preservation rows (diet x module): ", nrow(pres))
module_max <- pres[, .(max_Z = max(Zsummary, na.rm = TRUE),
                       mean_Z = mean(Zsummary, na.rm = TRUE),
                       module_size = module_size[1]),
                   by = module]
message("  unique modules: ", nrow(module_max))
message("  max Zsummary distribution:")
print(round(quantile(module_max$max_Z, c(0.5, 0.75, 0.9, 1.0),
                     na.rm = TRUE), 2))

well_preserved <- module_max[max_Z >= Z_STRICT]
threshold_used <- Z_STRICT
if (nrow(well_preserved) == 0) {
  message("  no modules at Zsummary >= ", Z_STRICT,
          "; falling back to Z >= ", Z_FALLBACK)
  well_preserved <- module_max[max_Z >= Z_FALLBACK]
  threshold_used <- Z_FALLBACK
}
message("  modules retained (Zsummary >= ", threshold_used, "): ",
        nrow(well_preserved))

mod_assg <- fread(MOD_ASSG)
stopifnot(all(c("gene", "module") %in% colnames(mod_assg)))
mod_assg <- mod_assg[gene %in% V]
mod_assg <- mod_assg[module %in% well_preserved$module]
message("  gene-module rows after V filter + well-preserved modules: ",
        nrow(mod_assg))

# Emit within-module pairs, guarding against oversized modules
channel1 <- mod_assg[, {
  if (.N >= 2 && .N <= MAX_MODULE_SIZE) {
    idx <- combn(.N, 2)
    ga <- gene[idx[1, ]]; gb <- gene[idx[2, ]]
    swap <- ga > gb
    tmp <- ga[swap]; ga[swap] <- gb[swap]; gb[swap] <- tmp
    list(gene_a = ga, gene_b = gb, module_id = as.integer(module[1]))
  } else {
    if (.N > MAX_MODULE_SIZE) {
      message("    skipping oversized module ", module[1], " (", .N, " genes)")
    }
    NULL
  }
}, by = module]
channel1[, channel := "wgcna_module"]
channel1 <- merge(channel1, well_preserved[, .(module, max_Z, mean_Z)],
                  by.x = "module_id", by.y = "module", all.x = TRUE)
setnames(channel1, c("max_Z", "mean_Z"),
         c("preservation_Z", "preservation_Z_mean"))
channel1 <- channel1[, .(gene_a, gene_b, channel, module_id,
                         preservation_Z, preservation_Z_mean,
                         celltype = NA_character_)]
message("  channel 1 WGCNA module edges: ", nrow(channel1))

# ==========================================================================
#  Channel 2: Cell-type conserved gene cliques
# ==========================================================================
ct <- fread(CT_CONS)
message("  cell-type conserved gene rows: ", nrow(ct))

# Require both sig flag and V membership
if (!"bulk_sig" %in% colnames(ct)) ct[, bulk_sig := TRUE]
ct <- ct[bulk_sig == TRUE & human_symbol %in% V]

# Build per-celltype cliques (concordant direction)
if ("human_logFC" %in% colnames(ct) && "mouse_logFC" %in% colnames(ct)) {
  ct[, concordant_direction := sign(human_logFC) == sign(mouse_logFC)]
  ct <- ct[concordant_direction == TRUE]
}

message("  cell-type conserved genes post-filter: ", nrow(ct))
message("  unique cell types: ", uniqueN(ct$cell_type))

channel2 <- ct[, {
  if (.N >= 2) {
    idx <- combn(.N, 2)
    ga <- human_symbol[idx[1, ]]; gb <- human_symbol[idx[2, ]]
    swap <- ga > gb
    tmp <- ga[swap]; ga[swap] <- gb[swap]; gb[swap] <- tmp
    list(gene_a = ga, gene_b = gb, celltype = cell_type[1])
  } else NULL
}, by = cell_type]
channel2[, cell_type := NULL]
channel2[, channel := "celltype_conserved"]
channel2[, module_id := NA_integer_]
channel2[, preservation_Z := NA_real_]
channel2[, preservation_Z_mean := NA_real_]
channel2 <- channel2[, .(gene_a, gene_b, channel, module_id,
                         preservation_Z, preservation_Z_mean, celltype)]
message("  channel 2 celltype-conserved edges: ", nrow(channel2))

# ==========================================================================
#  Combine + emit
# ==========================================================================
all_edges <- rbind(channel1, channel2)
all_edges <- all_edges[gene_a != gene_b]

# Dedup: same pair in multiple channels -> aggregate channel tags + keep best Z
setorder(all_edges, gene_a, gene_b, -preservation_Z)
dedup <- all_edges[, .(
  channel              = paste(sort(unique(channel)), collapse = ";"),
  module_id            = module_id[1],
  preservation_Z       = max(preservation_Z, na.rm = TRUE),
  preservation_Z_mean  = max(preservation_Z_mean, na.rm = TRUE),
  celltype             = paste(sort(unique(na.omit(celltype))), collapse = ";"),
  raw_score            = max(preservation_Z / 10, 0.2, na.rm = TRUE)
), by = .(gene_a, gene_b)]
# When preservation_Z is all NA (celltype-only), raw_score falls to 0.2 floor
dedup[!is.finite(preservation_Z), preservation_Z := NA_real_]
dedup[!is.finite(preservation_Z_mean), preservation_Z_mean := NA_real_]
dedup[celltype == "", celltype := NA_character_]

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(dedup, OUT)
message("[", Sys.time(), "] Wrote ", nrow(dedup),
        " D-XS conserved edges -> ", OUT)
message("  Genes touched: ", length(unique(c(dedup$gene_a, dedup$gene_b))))
message("  Channel breakdown:")
print(dedup[, .N, by = channel])
message("[297] Elapsed: ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
        " min")
