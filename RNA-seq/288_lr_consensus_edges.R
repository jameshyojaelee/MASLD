#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --job-name=net_288_dlr_consensus
#SBATCH --output=logs/net_288_dlr_consensus_%j.out
#SBATCH --error=logs/net_288_dlr_consensus_%j.err
# ===========================================================================
# Script 288: Expanded D-LR edges (Phase 3 sparse-layer rehab).
# ---------------------------------------------------------------------------
# 284 (strict differential) yielded <10 LR edges — too few to justify a
# portal tab. 253 already built a broader `edges_lr.csv` (~1,500 LIANA
# consensus edges across cell-type pairs) from
# Analysis/SingleCell/results_gpu_v2/fig2_data/liana_differential_interactions.csv
# but that file is keyed on symbols+metadata-strings rather than the
# structured schema the portal atlas merge (290) expects.
#
# 288 reads edges_lr.csv, parses the pipe-delimited metadata column into
# structured (cell_type_pairs, direction, n_contexts), and adds
# `differential_at_F2` if the same gene-pair appears in 284's edges_d_lr.csv
# as a differential hit -- preserving the F2-dynamic subset as a filter.
#
# Optional annotation sources (present if available):
#   Analysis/SingleCell/results_gpu_v2/sex_ccc/sex_specific_ccc_axes.csv
#   Analysis/SingleCell/results_gpu_v2/metabolic_ccc/hep_to_mac_metabolic_ligand_axes.csv
#
# Output: edges_d_lr_consensus.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
})

t0 <- Sys.time()

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

LR_BASE   <- file.path(PROJ, "RNA-seq/results/network/edges_lr.csv")
D_LR_DIFF <- file.path(PROJ, "RNA-seq/results/network/edges_d_lr.csv")
SEX_CCC   <- file.path(PROJ, "Analysis/SingleCell/results_gpu_v2/sex_ccc/sex_specific_ccc_axes.csv")
MET_CCC   <- file.path(PROJ, "Analysis/SingleCell/results_gpu_v2/metabolic_ccc/hep_to_mac_metabolic_ligand_axes.csv")
NODES     <- file.path(PROJ, "RNA-seq/results/network/network_nodes.csv")
OUT       <- file.path(PROJ, "RNA-seq/results/network/edges_d_lr_consensus.csv")

message("[288] LR consensus edge construction starting at ", format(t0))

stopifnot(file.exists(LR_BASE), file.exists(NODES))

nodes <- fread(NODES); V <- unique(nodes$human_symbol)
message("  Node set V: ", length(V), " genes")

# ==========================================================================
#  Base table: 253's edges_lr.csv -- gene_a, gene_b, raw_score, metadata
# ==========================================================================
base <- fread(LR_BASE)
stopifnot(all(c("gene_a", "gene_b", "raw_score", "metadata") %in% colnames(base)))
message("  base LR edges (from 253): ", nrow(base))

# Parse metadata string: "lr|celltype_pairs=A->B;C->D|direction=MASLD_up|n_contexts=42"
extract_field <- function(metastr, key) {
  pat <- paste0("(?:^|[|])", key, "=([^|]*)")
  m <- regmatches(metastr, regexpr(pat, metastr, perl = TRUE))
  out <- rep(NA_character_, length(metastr))
  hit <- regexpr(pat, metastr, perl = TRUE)
  start <- as.integer(hit)
  matches <- regmatches(metastr, hit)
  has <- start > 0 & nchar(matches) > 0
  val <- sub(paste0("^(?:^|[|])", key, "="), "", matches, perl = TRUE)
  out[has] <- val[has]
  out
}

base[, cell_type_pairs := extract_field(metadata, "celltype_pairs")]
base[, direction       := extract_field(metadata, "direction")]
base[, n_contexts      := suppressWarnings(as.integer(
                                 extract_field(metadata, "n_contexts")))]

# Restrict to V on both endpoints
base <- base[gene_a %in% V & gene_b %in% V]
message("  after V filter: ", nrow(base))

# ==========================================================================
#  Differential flag: gene-pair present in 284's edges_d_lr.csv
# ==========================================================================
differential_pairs <- data.table()
if (file.exists(D_LR_DIFF)) {
  d_lr <- fread(D_LR_DIFF)
  if ("gene_a" %in% colnames(d_lr) && "gene_b" %in% colnames(d_lr)) {
    differential_pairs <- d_lr[, .(gene_a, gene_b)]
    differential_pairs[, key := paste(pmin(gene_a, gene_b),
                                      pmax(gene_a, gene_b), sep = "|")]
    message("  differential F2-dynamic LR pairs from 284: ", nrow(differential_pairs))
  }
}
base[, key := paste(pmin(gene_a, gene_b), pmax(gene_a, gene_b), sep = "|")]
base[, differential_at_F2 := key %in% differential_pairs$key]

# ==========================================================================
#  Optional annotations
# ==========================================================================
sex_pairs <- data.table()
if (file.exists(SEX_CCC)) {
  sex_tab <- tryCatch(fread(SEX_CCC), error = function(e) NULL)
  if (!is.null(sex_tab) && all(c("ligand", "receptor") %in% colnames(sex_tab))) {
    sex_pairs <- sex_tab[, .(ligand, receptor)]
    sex_pairs[, key := paste(pmin(ligand, receptor),
                             pmax(ligand, receptor), sep = "|")]
    message("  sex-specific LR axes: ", nrow(sex_pairs))
  }
}
base[, sex_specific := key %in% sex_pairs$key]

metabolic_pairs <- data.table()
if (file.exists(MET_CCC)) {
  met_tab <- tryCatch(fread(MET_CCC), error = function(e) NULL)
  if (!is.null(met_tab)) {
    # Heuristic column picks
    lg_col <- colnames(met_tab)[grepl("ligand", tolower(colnames(met_tab)))][1]
    rg_col <- colnames(met_tab)[grepl("receptor", tolower(colnames(met_tab)))][1]
    if (!is.na(lg_col) && !is.na(rg_col)) {
      metabolic_pairs <- met_tab[, c(lg_col, rg_col), with = FALSE]
      setnames(metabolic_pairs, c("ligand", "receptor"))
      metabolic_pairs[, key := paste(pmin(ligand, receptor),
                                     pmax(ligand, receptor), sep = "|")]
      message("  hep->macrophage metabolic ligand axes: ", nrow(metabolic_pairs))
    }
  }
}
base[, metabolic_hep_mac := key %in% metabolic_pairs$key]

# ==========================================================================
#  Canonicalize + emit
# ==========================================================================
swap <- base$gene_a > base$gene_b
a2 <- base$gene_a; b2 <- base$gene_b
base[, gene_a := ifelse(swap, b2, a2)]
base[, gene_b := ifelse(swap, a2, b2)]
base <- base[gene_a != gene_b]
base <- unique(base, by = c("gene_a", "gene_b"))
base[, key := NULL]

out <- base[, .(gene_a, gene_b,
                score_diff_lr   = raw_score,
                cell_type_pairs,
                direction,
                n_contexts,
                differential_at_F2,
                sex_specific,
                metabolic_hep_mac)]

setorder(out, -score_diff_lr, gene_a, gene_b)

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(out, OUT)
message("[", Sys.time(), "] Wrote ", nrow(out),
        " D-LR consensus edges -> ", OUT)

message("  Genes touched: ", length(unique(c(out$gene_a, out$gene_b))),
        " / ", length(V))
message("  Direction breakdown:")
print(out[, .N, by = direction])
message("  Flags:")
message("    differential_at_F2: ", sum(out$differential_at_F2, na.rm = TRUE))
message("    sex_specific:       ", sum(out$sex_specific, na.rm = TRUE))
message("    metabolic_hep_mac:  ", sum(out$metabolic_hep_mac, na.rm = TRUE))
message("[288] Elapsed: ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
        " min")
