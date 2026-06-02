#!/usr/bin/env Rscript
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_285_dcerna
#SBATCH --output=logs/net_285_dcerna_%j.out
#SBATCH --error=logs/net_285_dcerna_%j.err
#
# 285_cerna_masld_expressed.R
# Phase 2 of MASLD network v2 (Architecture C): Build D-ceRNA edges
#
# Refactor of the ceRNA portion of 253_regulon_lr_cerna_edges.R to:
#   1. Use ALL ceRNA database sources (ENCORI, TargetScan, miRTarBase, LncBase)
#   2. Require MASLD expression of both lncRNA and mRNA endpoints
#   3. Preserve lncRNA-mRNA provenance (biotype tags on each endpoint)
#
# Output: RNA-seq/results/network/edges_d_cerna.csv

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

t0 <- Sys.time()
cat("=== 285: D-ceRNA Edges (MASLD-expressed) ===\n")
cat("Start:", format(t0, "%Y-%m-%d %H:%M:%S"), "\n\n")

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

OUTDIR    <- file.path(BASE, "RNA-seq/results/network")
NODE_PATH <- file.path(OUTDIR, "network_nodes.csv")
DB_DIR    <- file.path(BASE, "data/ncrna_databases")
DGE_PATH  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")

ENCORI_LNC_F  <- file.path(DB_DIR, "encori_mirna_lncrna.csv")
ENCORI_MRNA_F <- file.path(DB_DIR, "encori_mirna_mrna.csv")
TARGETSCAN_F  <- file.path(DB_DIR, "targetscan_conserved.csv")
MIRTARBASE_F  <- file.path(DB_DIR, "mirtarbase_validated.csv")
LNCBASE_F     <- file.path(DB_DIR, "lncbase_interactions.csv")

OUT_PATH <- file.path(OUTDIR, "edges_d_cerna.csv")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Load node set V (with biotype) and build MASLD-expression filter
# ----------------------------------------------------------------------------
stopifnot(file.exists(NODE_PATH))
nodes <- fread(NODE_PATH)
V <- unique(nodes$human_symbol)
biotype_lookup <- setNames(nodes$gene_biotype, nodes$human_symbol)
cat("Node set V:", length(V), "genes\n")

cat("Loading merged DGE for expression filter ...\n")
stopifnot(file.exists(DGE_PATH))
dge <- readRDS(DGE_PATH)

# Compute logCPM and require >1 in >=20% of samples
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
prop_expr <- rowMeans(logcpm > 1)
expr_rownames <- rownames(logcpm)

# Try to map from rownames (may be Ensembl or symbol) to symbols in V
expressed_mask <- prop_expr >= 0.20
expressed_ids  <- expr_rownames[expressed_mask]

# Attempt to translate via nodes$ensembl_id
ens_to_sym <- setNames(nodes$human_symbol, nodes$ensembl_id)
expressed_symbols <- unique(c(
  expressed_ids[expressed_ids %in% V],                    # if rownames are symbols
  ens_to_sym[expressed_ids[expressed_ids %in% names(ens_to_sym)]]  # if rownames are ensembl
))
expressed_symbols <- expressed_symbols[!is.na(expressed_symbols) & expressed_symbols != ""]

cat("  logCPM matrix:", nrow(logcpm), "x", ncol(logcpm), "\n")
cat("  Genes passing expression filter (>1 logCPM in >=20% samples):",
    sum(expressed_mask), "\n")
cat("  Expressed symbols in V:", length(intersect(expressed_symbols, V)), "\n\n")

EXPRESSED_V <- intersect(expressed_symbols, V)

# ----------------------------------------------------------------------------
# Load ceRNA database sources
# ----------------------------------------------------------------------------
read_if <- function(path) {
  if (!file.exists(path)) {
    cat("  MISSING:", path, "\n"); return(NULL)
  }
  fread(path)
}

cat("--- Loading ceRNA databases ---\n")
encori_lnc  <- read_if(ENCORI_LNC_F)
encori_mrna <- read_if(ENCORI_MRNA_F)
targetscan  <- read_if(TARGETSCAN_F)
mirtarbase  <- read_if(MIRTARBASE_F)
lncbase     <- read_if(LNCBASE_F)

cat("  ENCORI miRNA-lncRNA:", if (is.null(encori_lnc)) 0 else nrow(encori_lnc), "\n")
cat("  ENCORI miRNA-mRNA:  ", if (is.null(encori_mrna)) 0 else nrow(encori_mrna), "\n")
cat("  TargetScan:         ", if (is.null(targetscan)) 0 else nrow(targetscan), "\n")
cat("  miRTarBase:         ", if (is.null(mirtarbase)) 0 else nrow(mirtarbase), "\n")
cat("  LncBase:            ", if (is.null(lncbase)) 0 else nrow(lncbase), "\n\n")

# ----------------------------------------------------------------------------
# Union of miRNA-lncRNA and miRNA-mRNA interactions
# ----------------------------------------------------------------------------
lnc_tab <- rbindlist(list(
  if (!is.null(encori_lnc)) encori_lnc[, .(mirna, lncrna, source)] else NULL,
  if (!is.null(lncbase))    lncbase[,    .(mirna, lncrna, source)] else NULL
), use.names = TRUE, fill = TRUE)

mrna_tab <- rbindlist(list(
  if (!is.null(encori_mrna)) encori_mrna[, .(mirna, mrna, source)] else NULL,
  if (!is.null(targetscan))  targetscan[,  .(mirna, mrna, source)] else NULL,
  if (!is.null(mirtarbase))  mirtarbase[,  .(mirna, mrna, source)] else NULL
), use.names = TRUE, fill = TRUE)

# Normalise: drop missing, de-duplicate on (mirna, gene)
lnc_tab  <- unique(lnc_tab[!is.na(mirna) & !is.na(lncrna) & mirna != "" & lncrna != ""])
mrna_tab <- unique(mrna_tab[!is.na(mirna) & !is.na(mrna)  & mirna != "" & mrna  != ""])

cat("Union miRNA-lncRNA interactions: ", nrow(lnc_tab),
    " (", uniqueN(lnc_tab$mirna), " miRNAs,", uniqueN(lnc_tab$lncrna), " lncRNAs)\n", sep = "")
cat("Union miRNA-mRNA interactions:   ", nrow(mrna_tab),
    " (", uniqueN(mrna_tab$mirna), " miRNAs,", uniqueN(mrna_tab$mrna),  " mRNAs)\n\n", sep = "")

# ----------------------------------------------------------------------------
# Filter endpoints to MASLD-expressed genes in V
# ----------------------------------------------------------------------------
lnc_tab  <- lnc_tab[lncrna %in% EXPRESSED_V]
mrna_tab <- mrna_tab[mrna  %in% EXPRESSED_V]

cat("After V + MASLD-expression filter:\n")
cat("  lncRNA side:", nrow(lnc_tab),  "(", uniqueN(lnc_tab$lncrna), "lncRNAs)\n")
cat("  mRNA side:  ", nrow(mrna_tab), "(", uniqueN(mrna_tab$mrna),  "mRNAs)\n\n")

# ----------------------------------------------------------------------------
# Build lncRNA -> miRNA set map and count total targets per lncRNA
# ----------------------------------------------------------------------------
lnc_to_mirs <- lnc_tab[, .(mirs = list(unique(mirna)),
                           n_mirs = uniqueN(mirna)), by = lncrna]

# Total miRNA targets per lncRNA (size of miRNA target set) for normalisation
max_lnc_mirs <- if (nrow(lnc_to_mirs) > 0) max(lnc_to_mirs$n_mirs) else 1L

# ----------------------------------------------------------------------------
# For each miRNA, enumerate mRNA targets for the join
# ----------------------------------------------------------------------------
setkey(mrna_tab, mirna)

# Expand lnc_to_mirs to (lncrna, mirna) long form, then join to mRNA targets
lnc_long <- lnc_tab[, .(mirna = unique(mirna)), by = lncrna]

# lncRNA x mRNA triplets sharing each miRNA
triplets <- merge(lnc_long, mrna_tab[, .(mirna, mrna)],
                  by = "mirna", allow.cartesian = TRUE)
triplets <- triplets[lncrna != mrna]
cat("Candidate (lncRNA, mRNA, miRNA) triplets:", nrow(triplets), "\n")

if (nrow(triplets) == 0) {
  cat("\nNo triplets found — writing empty output.\n")
  empty <- data.table(gene_a = character(), gene_b = character(),
                      n_shared_mirnas = integer(), shared_mirnas = character(),
                      cerna_score = numeric(),
                      biotype_a = character(), biotype_b = character(),
                      score_norm = numeric())
  fwrite(empty, OUT_PATH)
  cat("Wrote empty file to", OUT_PATH, "\n")
  quit(status = 0)
}

# ----------------------------------------------------------------------------
# Aggregate per (lncRNA, mRNA) pair
# ----------------------------------------------------------------------------
pair_tab <- triplets[, .(
  n_shared_mirnas = uniqueN(mirna),
  shared_mirnas   = paste(sort(unique(mirna)), collapse = ";")
), by = .(lncrna, mrna)]

# Normalise: n_shared / total miRNAs targeting the lncRNA (i.e. fraction of
# the lncRNA's miRNA "sponge capacity" shared with the mRNA). Also provide a
# global score_norm = n_shared / max(total_miRNA_targets_of_lncRNA).
lnc_total <- setNames(lnc_to_mirs$n_mirs, lnc_to_mirs$lncrna)
pair_tab[, cerna_score := n_shared_mirnas / lnc_total[lncrna]]
pair_tab[, score_norm  := n_shared_mirnas / max_lnc_mirs]

# Filter: n_shared >= 1 (by construction) and score > 0.01
pair_tab <- pair_tab[n_shared_mirnas >= 1 & cerna_score > 0.01]

# ----------------------------------------------------------------------------
# Canonical ordering (alphabetical) and biotype tags
# ----------------------------------------------------------------------------
pair_tab[, c("gene_a", "gene_b") := .(pmin(lncrna, mrna), pmax(lncrna, mrna))]

# Biotypes from node table
pair_tab[, biotype_a := biotype_lookup[gene_a]]
pair_tab[, biotype_b := biotype_lookup[gene_b]]

out <- pair_tab[, .(gene_a, gene_b,
                    n_shared_mirnas, shared_mirnas,
                    cerna_score, biotype_a, biotype_b, score_norm)]

# Deduplicate — a pair could appear once (lnc-mrna); collapse if any duplication
setorder(out, gene_a, gene_b, -n_shared_mirnas)
out <- out[!duplicated(out, by = c("gene_a", "gene_b"))]

# ----------------------------------------------------------------------------
# Write output + report
# ----------------------------------------------------------------------------
fwrite(out, OUT_PATH)

cat("\n--- Per-source summary ---\n")
cat("  n unique lncRNAs in edges:", uniqueN(pair_tab$lncrna), "\n")
cat("  n unique mRNAs in edges:  ", uniqueN(pair_tab$mrna),   "\n")
cat("  Total D-ceRNA edges:      ", nrow(out),                "\n")
if (nrow(out) > 0) {
  cat("  Mean n_shared_mirnas:   ", round(mean(out$n_shared_mirnas), 3), "\n")
  cat("  Mean cerna_score:       ", round(mean(out$cerna_score),     4), "\n")
  cat("  Max n_shared_mirnas:    ", max(out$n_shared_mirnas),          "\n")
  cat("\n  Top 5 edges by n_shared_mirnas:\n")
  top <- out[order(-n_shared_mirnas, -cerna_score)][1:min(5, .N)]
  for (r in seq_len(nrow(top))) {
    cat("    ", top$gene_a[r], " -- ", top$gene_b[r],
        " (n_shared=", top$n_shared_mirnas[r],
        ", score=", round(top$cerna_score[r], 3),
        ", ", top$biotype_a[r], "/", top$biotype_b[r], ")\n", sep = "")
  }
}

cat("\nWrote:", OUT_PATH, "\n")
cat("Done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
    "minutes.\n")
