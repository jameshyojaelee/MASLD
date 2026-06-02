#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --job-name=net_289_dcerna_triplet
#SBATCH --output=logs/net_289_dcerna_triplet_%j.out
#SBATCH --error=logs/net_289_dcerna_triplet_%j.err
# ===========================================================================
# Script 289: ceRNA triplet -> pairwise edge expansion (Phase 3 rehab).
# ---------------------------------------------------------------------------
# The ncRNA pipeline (53-57) emitted ~11 curated ceRNA triplets
# (cerna_network_full.csv) and ~6 validated axes
# (cerna_validated_axes.csv). Both are lncRNA-mRNA pairs sharing a set of
# miRNAs.
#
# This script emits edges_d_cerna_triplet.csv in the edges_d_* schema
# (gene_a, gene_b, n_shared_mirnas, shared_mirnas + provenance columns).
# Due to the tiny size (<< 500-edge floor in the plan), the atlas merge
# step will collapse this layer to a node-level badge in the portal
# rather than a full layer tab. This script just produces the substrate;
# the layer-vs-badge decision lives in 290_edge_annotation_atlas.py.
#
# Output: edges_d_cerna_triplet.csv
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
})

t0 <- Sys.time()

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

NETDIR    <- file.path(PROJ, "RNA-seq/results/network")
FULL_PATH <- file.path(PROJ, "RNA-seq/results/ncrna/cerna_network_full.csv")
MASLD_PATH <- file.path(PROJ, "RNA-seq/results/ncrna/cerna_network_masld.csv")
VALID_PATH <- file.path(PROJ, "RNA-seq/results/ncrna/cerna_validated_axes.csv")
NODES      <- file.path(NETDIR, "network_nodes.csv")
OUT        <- file.path(NETDIR, "edges_d_cerna_triplet.csv")

message("[289] ceRNA triplet -> edge expansion starting at ", format(t0))

stopifnot(file.exists(NODES))
nodes <- fread(NODES); V <- unique(nodes$human_symbol)
message("  Node set V: ", length(V), " genes")

# Union of inputs, each annotated with provenance
combine_one <- function(path, src) {
  if (!file.exists(path)) {
    message("  missing: ", path, " -- skipping")
    return(data.table())
  }
  dt <- fread(path)
  required <- c("lncrna", "mrna", "shared_mirnas", "cerna_score")
  if (!all(required %in% colnames(dt))) {
    # cerna_validated_axes.csv may use different columns; be permissive
    if (!all(c("lncrna", "mrna") %in% colnames(dt))) return(data.table())
    for (col in setdiff(required, colnames(dt))) dt[, (col) := NA]
  }
  dt[, provenance := src]
  dt[, .(lncrna, mrna, shared_mirnas, cerna_score, provenance)]
}

all_triplets <- rbindlist(list(
  combine_one(FULL_PATH,  "ncrna_full"),
  combine_one(MASLD_PATH, "ncrna_masld_expressed"),
  combine_one(VALID_PATH, "ncrna_validated")
), fill = TRUE)
message("  raw triplet rows (all sources): ", nrow(all_triplets))

# Restrict to both endpoints in V
all_triplets <- all_triplets[lncrna %in% V & mrna %in% V & lncrna != mrna]
message("  after V filter: ", nrow(all_triplets))

# Canonicalize
triplets_c <- all_triplets[, {
  ga <- pmin(lncrna, mrna); gb <- pmax(lncrna, mrna)
  list(gene_a = ga, gene_b = gb,
       shared_mirnas = shared_mirnas,
       n_shared_mirnas = suppressWarnings(as.integer(cerna_score)),
       provenance = provenance)
}]

# For duplicate pairs across sources, keep the highest shared_mirna count
# and aggregate provenance.
setorder(triplets_c, gene_a, gene_b, -n_shared_mirnas)
dedup <- triplets_c[, .(
  n_shared_mirnas = n_shared_mirnas[1],
  shared_mirnas   = shared_mirnas[1],
  provenance      = paste(sort(unique(provenance)), collapse = ";"),
  raw_score       = n_shared_mirnas[1] /
                    max(1L, max(triplets_c$n_shared_mirnas, na.rm = TRUE))
), by = .(gene_a, gene_b)]

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(dedup, OUT)
message("[", Sys.time(), "] Wrote ", nrow(dedup),
        " ceRNA triplet edges -> ", OUT)

message("  Genes touched: ", length(unique(c(dedup$gene_a, dedup$gene_b))))
message("  Provenance breakdown:")
print(dedup[, .N, by = provenance])
message("[289] Elapsed: ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
        " min")
