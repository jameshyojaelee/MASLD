#!/usr/bin/env Rscript
#SBATCH --job-name=net_293_dreg
#SBATCH --partition=cpu
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --output=logs/net_293_dreg_%j.out
#SBATCH --error=logs/net_293_dreg_%j.err
# ===========================================================================
# Script 293: D-REG edges -- directed TF -> target edges from SCENIC+ regulons.
# ---------------------------------------------------------------------------
# Reads hepatocyte_regulons.csv (pre-exploded TF-target rows) plus
# disease_regulons.csv (target_genes as semicolon-joined strings) and emits
# edges_d_regulon.csv in the same schema family as edges_d_coloc /
# edges_d_lr (gene_a, gene_b, + layer-specific columns). Edge type = D-REG,
# directed (gene_a = TF, gene_b = target; no canonical swap).
#
# Inputs:
#   Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv
#   Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv
#   RNA-seq/results/network/network_nodes.csv     (node set V + symbol map)
#
# Output:
#   RNA-seq/results/network/edges_d_regulon.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
})

t0 <- Sys.time()

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

NETDIR    <- file.path(PROJ, "RNA-seq/results/network")
NODE_PATH <- file.path(NETDIR, "network_nodes.csv")
HEP_PATH  <- file.path(PROJ,
                       "Analysis/ATAC/Human_Multiome/scenic_plus/hepatocyte_regulons.csv")
DIS_PATH  <- file.path(PROJ,
                       "Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv")
OUT       <- file.path(NETDIR, "edges_d_regulon.csv")

cat("[293] D-REG edge construction\n")
cat("  Start:", format(t0, "%Y-%m-%d %H:%M:%S"), "\n")

stopifnot(file.exists(NODE_PATH), file.exists(HEP_PATH), file.exists(DIS_PATH))

nodes <- fread(NODE_PATH)
V <- unique(nodes$human_symbol)
cat("[293] Node set V:", length(V), "genes\n")

# ---- Hepatocyte regulons: already exploded (tf_name, target_gene rows) ----
hep <- fread(HEP_PATH)
hep_edges <- hep[, .(
  tf                      = tf_name,
  target                  = target_gene,
  regulon_id              = regulon_id,
  regulon_source          = "hepatocyte",
  regulon_activity_diff   = regulon_activity_diff,
  activity_padj           = activity_padj
)]
cat("[293] Hepatocyte regulon rows:", nrow(hep_edges), "\n")

# ---- Disease regulons: target_genes is a semicolon-joined string ----
dis <- fread(DIS_PATH)
dis_exploded <- dis[, {
  tgts <- strsplit(target_genes, ";", fixed = TRUE)[[1]]
  .(
    tf                      = rep(tf_name, length(tgts)),
    target                  = tgts,
    regulon_id              = rep(regulon_id, length(tgts)),
    regulon_source          = rep("disease", length(tgts)),
    regulon_activity_diff   = rep(regulon_activity_diff, length(tgts)),
    activity_padj           = rep(activity_padj, length(tgts))
  )
}, by = seq_len(nrow(dis))]
dis_exploded[, seq_len := NULL]
cat("[293] Disease regulon edges (post-explode):", nrow(dis_exploded), "\n")

# ---- Union ----
all_edges <- rbind(hep_edges, dis_exploded)
cat("[293] Union rows:", nrow(all_edges), "\n")

# Filter to V on both endpoints. Some targets are raw ensembl_ids (e.g.,
# ENSG00000237773) rather than symbols -- drop those since the portal is
# symbol-keyed.
before <- nrow(all_edges)
all_edges <- all_edges[tf %in% V & target %in% V & tf != target]
cat("[293] After V filter (both endpoints symbols in V, self-edges removed):",
    nrow(all_edges), "(dropped", before - nrow(all_edges), ")\n")

# Deduplicate (TF, target): a TF-target pair may appear in both hepatocyte
# and disease regulons. Keep the row with strongest |activity_diff|.
all_edges[, abs_diff := abs(regulon_activity_diff)]
setorder(all_edges, tf, target, -abs_diff)
dedup <- all_edges[!duplicated(all_edges, by = c("tf", "target"))]
# Merge source tag across duplicates
src_tag <- all_edges[, .(regulon_source = paste(sort(unique(regulon_source)),
                                                collapse = ";")),
                     by = .(tf, target)]
dedup[, regulon_source := NULL]
dedup <- merge(dedup, src_tag, by = c("tf", "target"))

# Score: 1 - activity_padj (higher = more confident), floored at 0
dedup[, raw_score := pmax(0, pmin(1, 1 - activity_padj))]

# Emit in edges_d_* schema:
#   gene_a = TF, gene_b = target (directed; no canonical swap)
out <- dedup[, .(
  gene_a                = tf,
  gene_b                = target,
  tf                    = tf,
  target                = target,
  regulon_id            = regulon_id,
  regulon_source        = regulon_source,
  regulon_activity_diff = regulon_activity_diff,
  activity_padj         = activity_padj,
  raw_score             = raw_score,
  directed              = TRUE
)]

setorder(out, -raw_score, gene_a, gene_b)

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(out, OUT)
cat("[293] Wrote", nrow(out), "D-REG edges to", OUT, "\n")

# ---- Diagnostics ----
cat("[293] Unique TFs:", uniqueN(out$tf), "\n")
cat("[293] Unique targets:", uniqueN(out$gene_b), "\n")
cat("[293] Genes touched (tf U target):",
    length(unique(c(out$gene_a, out$gene_b))),
    "/", length(V),
    sprintf("(%.1f%%)", 100 * length(unique(c(out$gene_a, out$gene_b))) / length(V)),
    "\n")
cat("[293] Source breakdown:\n")
print(out[, .N, by = regulon_source])
cat("[293] Top 10 TFs by out-degree:\n")
print(out[, .N, by = tf][order(-N)][1:10])

cat("[293] Elapsed:",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
    "min\n")
