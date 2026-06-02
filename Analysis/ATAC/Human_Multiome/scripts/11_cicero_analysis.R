#!/usr/bin/env Rscript
# ==========================================================================
# Cicero co-accessibility analysis on hepatocyte scATAC-seq data
#
# Input:  Analysis/ATAC/Human_Multiome/results/cicero/ (from 10_cicero_prepare.py)
# Output: Analysis/ATAC/Human_Multiome/results/cicero/ (CSVs + RDS)
#         figures/supplementary/epigenomic/ (panels 17-21)
# ==========================================================================

suppressPackageStartupMessages({
  library(cicero)
  library(monocle)
  library(VGAM)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# --- Paths ----------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts", "figures", "publication_theme.R"))

CICERO_DIR <- file.path(BASE, "Analysis", "ATAC", "Human_Multiome", "results", "cicero")
FIG_OUT    <- file.path(BASE, "figures", "supplementary", "figS05_epigenomic_spatial")
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(p, name, w = fig_half_width, h = 3) {
  path <- file.path(FIG_OUT, name)
  save_fig(p, path, width = w, height = h)
  message("  Saved ", path)
}

cat("==============================================================\n")
cat("Cicero Co-accessibility Analysis\n")
cat("==============================================================\n")

# ==========================================================================
# 1. Load data and create CDS via make_atac_cds
# ==========================================================================
message("\n1. Loading Cicero input...")

# Read triplet format directly (peak \t cell \t count)
cicero_input <- fread(file.path(CICERO_DIR, "cicero_input.tsv"),
                      header = FALSE, col.names = c("Peak", "Cell", "Count"))
peaks    <- readLines(file.path(CICERO_DIR, "peaks.tsv"))
barcodes <- readLines(file.path(CICERO_DIR, "barcodes.tsv"))
meta     <- fread(file.path(CICERO_DIR, "cell_metadata.csv"))

message(sprintf("  Triplets: %s entries, %d peaks, %d cells",
                formatC(nrow(cicero_input), big.mark = ","),
                length(peaks), length(barcodes)))

# Create CDS manually (make_atac_cds has R 4.x class() bug)
# Build sparse matrix from triplets
peak_ids <- factor(cicero_input$Peak, levels = peaks)
cell_ids <- factor(cicero_input$Cell, levels = barcodes)
mat <- sparseMatrix(
  i = as.integer(peak_ids),
  j = as.integer(cell_ids),
  x = rep(1L, nrow(cicero_input)),  # binarized
  dims = c(length(peaks), length(barcodes)),
  dimnames = list(peaks, barcodes)
)
# Remove any zero-count peaks
keep <- rowSums(mat) > 0
mat <- mat[keep, ]
message(sprintf("  Sparse matrix: %d peaks x %d cells", nrow(mat), ncol(mat)))

# Build CDS (monocle v2 CellDataSet)
fd <- data.frame(site_name = rownames(mat), row.names = rownames(mat))
pd <- data.frame(cells = colnames(mat), row.names = colnames(mat))
fd <- new("AnnotatedDataFrame", data = fd)
pd <- new("AnnotatedDataFrame", data = pd)
input_cds <- newCellDataSet(mat, phenoData = pd, featureData = fd,
                             expressionFamily = VGAM::negbinomial.size())
message("  CDS created successfully")

# Add cell metadata
pData(input_cds)$condition <- meta$condition[match(colnames(input_cds), barcodes)]
pData(input_cds)$donor_id  <- meta$donor_id[match(colnames(input_cds), barcodes)]

# ==========================================================================
# 2. Dimensionality reduction
# ==========================================================================
message("\n2. Dimensionality reduction...")

# Use pre-computed UMAP from SnapATAC2
umap_coords <- data.frame(
  UMAP_1 = meta$UMAP1,
  UMAP_2 = meta$UMAP2,
  row.names = barcodes
)
# Subset to cells in CDS
umap_coords <- umap_coords[colnames(input_cds), , drop = FALSE]

# Cicero needs a reduced dimension embedding for KNN aggregation
# We need to run detect_genes and dimensionality reduction on the CDS
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)

# Reduce dimensions with LSI (needed for make_cicero_cds internals)
set.seed(42)
input_cds <- reduceDimension(input_cds, max_components = 2,
                              num_dim = 6,
                              reduction_method = "tSNE",
                              norm_method = "none")

message("  Dimensionality reduction complete")

# ==========================================================================
# 3. Run Cicero
# ==========================================================================
message("\n3. Running Cicero co-accessibility...")

# Make Cicero CDS (aggregates cells into neighborhoods)
cicero_cds <- make_cicero_cds(input_cds,
                               reduced_coordinates = as.matrix(umap_coords),
                               k = 50)
message("  Cicero CDS created (k=50 neighborhoods)")

# GRCh38 chromosome sizes
chr_sizes <- data.frame(
  V1 = paste0("chr", c(1:22, "X", "Y")),
  V2 = c(248956422, 242193529, 198295559, 190214555, 181538259,
         170805979, 159345973, 145138636, 138394717, 133797422,
         135086622, 133275309, 114364328, 107043718, 101991189,
         90338345, 83257441, 80373285, 58617616, 64444167,
         46709983, 50818468, 156040895, 57227415)
)

# Run Cicero
message("  Computing co-accessibility (this may take 30-60 minutes)...")
conns <- run_cicero(cicero_cds, chr_sizes, window = 500000,
                    sample_num = 100, silent = FALSE)

message(sprintf("  Connections computed: %d peak pairs", nrow(conns)))
fwrite(conns, file.path(CICERO_DIR, "cicero_connections.csv"))

# Deduplicate reciprocal pairs (Cicero returns both A->B and B->A)
conns_dedup <- conns[Peak1 < Peak2]
message(sprintf("  Unique connections (deduplicated): %d", nrow(conns_dedup)))

# ==========================================================================
# 4. Identify CCANs
# ==========================================================================
message("\n4. Identifying CCANs...")

CCANS <- generate_ccans(conns, coaccess_cutoff_override = 0.1)
n_ccans <- if (nrow(CCANS) > 0 && "CCAN" %in% names(CCANS))
  length(unique(CCANS$CCAN[!is.na(CCANS$CCAN)])) else 0
message(sprintf("  CCANs identified: %d", n_ccans))
fwrite(CCANS, file.path(CICERO_DIR, "ccans.csv"))

# ==========================================================================
# 5. Link peaks to genes (Cicero gene activity)
# ==========================================================================
message("\n5. Computing Cicero gene activity scores...")

# Use Cicero's built-in gene annotation sample or build from GENCODE
gene_anno <- NULL  # will be built from GENCODE v49

# Better: try to use our own annotation
gencode_gtf <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
if (file.exists(gencode_gtf)) {
  message("  Building gene annotation from GENCODE v49...")
  # Read a subset of the GTF for gene-level annotations
  gtf_cmd <- paste0("zcat ", gencode_gtf,
                     " | awk -F'\\t' '$3==\"gene\"'")
  gtf_lines <- system(gtf_cmd, intern = TRUE)
  if (length(gtf_lines) > 0) {
    gene_anno_dt <- rbindlist(lapply(gtf_lines, function(line) {
      fields <- strsplit(line, "\t")[[1]]
      if (length(fields) < 9) return(NULL)
      attrs <- fields[9]
      # Extract gene_name
      gname <- sub('.*gene_name "([^"]+)".*', '\\1', attrs)
      if (gname == attrs) return(NULL)  # no match
      data.table(
        chromosome = fields[1],
        start = as.integer(fields[4]),
        end = as.integer(fields[5]),
        strand = fields[7],
        gene = gname
      )
    }))
    if (nrow(gene_anno_dt) > 0) {
      gene_anno <- as.data.frame(gene_anno_dt)
      message(sprintf("  Gene annotation: %d genes", nrow(gene_anno)))
    }
  }
}

# Add 'gene' column to featureData (required by build_gene_activity_matrix)
# For tiles, use tile coordinate as placeholder
fData(input_cds)$gene <- fData(input_cds)$site_name

# Build gene activity matrix using Cicero connections + gene annotation
gene_act <- tryCatch(
  build_gene_activity_matrix(input_cds, conns, coaccess_cutoff = 0.25),
  error = function(e) {
    message("  WARNING: build_gene_activity_matrix failed: ", e$message)
    message("  Skipping gene activity scores (connections + CCANs still valid)")
    NULL
  }
)
if (!is.null(gene_act)) {
  message(sprintf("  Gene activity matrix: %d genes x %d cells",
                  nrow(gene_act), ncol(gene_act)))
  saveRDS(gene_act, file.path(CICERO_DIR, "cicero_gene_activity.rds"))
} else {
  message("  Gene activity matrix skipped")
}

# ==========================================================================
# 6. Extract distal links
# ==========================================================================
message("\n6. Extracting distal links...")

conns_dt <- as.data.table(conns)
conns_dt <- conns_dt[!is.na(coaccess) & coaccess > 0.1]

# Parse coordinates
conns_dt[, c("chr1", "start1", "end1") := tstrsplit(Peak1, "_", type.convert = TRUE)]
conns_dt[, c("chr2", "start2", "end2") := tstrsplit(Peak2, "_", type.convert = TRUE)]
conns_dt[, distance := abs((start1 + end1)/2 - (start2 + end2)/2)]

distal_links <- conns_dt[distance > 10000 & coaccess > 0.2]
distal_links <- distal_links[Peak1 < Peak2]  # deduplicate reciprocal pairs
message(sprintf("  Distal links (>10kb, coaccess>0.2): %d", nrow(distal_links)))
fwrite(distal_links, file.path(CICERO_DIR, "distal_links.csv"))

# ==========================================================================
# 7. Disease-differential co-accessibility
# ==========================================================================
message("\n7. Differential co-accessibility (MASLD vs Normal)...")

conditions <- unique(meta$condition)
message(sprintf("  Conditions: %s", paste(conditions, collapse = ", ")))

if (length(conditions) >= 2) {
  cond_conns_list <- list()

  for (cond in conditions) {
    cond_safe <- gsub("[^A-Za-z0-9]", "_", cond)
    triplet_file <- file.path(CICERO_DIR, paste0("cicero_input_", cond_safe, ".tsv"))

    if (!file.exists(triplet_file)) {
      message(sprintf("  Skipping %s (file not found: %s)", cond, triplet_file))
      next
    }

    input_c <- fread(triplet_file, header = FALSE,
                     col.names = c("Peak", "Cell", "Count"))
    n_cells_c <- length(unique(input_c$Cell))

    if (n_cells_c < 500) {
      message(sprintf("  Skipping %s (%d cells < 500)", cond, n_cells_c))
      next
    }

    # Build CDS manually (R 4.x compat)
    peaks_c <- unique(input_c$Peak)
    cells_c <- unique(input_c$Cell)
    mat_c <- sparseMatrix(
      i = match(input_c$Peak, peaks_c),
      j = match(input_c$Cell, cells_c),
      x = rep(1L, nrow(input_c)),
      dims = c(length(peaks_c), length(cells_c)),
      dimnames = list(peaks_c, cells_c)
    )
    fd_c <- new("AnnotatedDataFrame",
                data = data.frame(site_name = peaks_c, row.names = peaks_c))
    pd_c <- new("AnnotatedDataFrame",
                data = data.frame(cells = cells_c, row.names = cells_c))
    cds_c <- newCellDataSet(mat_c, phenoData = pd_c, featureData = fd_c,
                             expressionFamily = VGAM::negbinomial.size())
    cds_c <- detectGenes(cds_c)
    cds_c <- estimateSizeFactors(cds_c)

    set.seed(42)
    cds_c <- reduceDimension(cds_c, max_components = 2, num_dim = 6,
                              reduction_method = "tSNE", norm_method = "none")

    # Use tSNE coords for make_cicero_cds
    tsne_c <- t(reducedDimA(cds_c))
    colnames(tsne_c) <- c("tSNE_1", "tSNE_2")

    cicero_c <- make_cicero_cds(cds_c, reduced_coordinates = tsne_c, k = 30)

    message(sprintf("  Running Cicero for %s (%d cells)...", cond, ncol(mat_c)))
    conns_c <- run_cicero(cicero_c, chr_sizes, window = 500000,
                          sample_num = 100, silent = TRUE)
    conns_c$condition <- cond
    cond_conns_list[[cond]] <- as.data.table(conns_c)
  }

  if (length(cond_conns_list) >= 2) {
    all_cond <- rbindlist(cond_conns_list, use.names = TRUE, fill = TRUE)
    fwrite(all_cond, file.path(CICERO_DIR, "condition_connections.csv"))

    # Compare by EXPLICIT disease-vs-control contrast (review 2026-05-30, F043).
    # The previous positional cn[1]/cn[2] used first-appearance order of
    # unique(meta$condition) = NORMAL, MASL, MASH, so it silently computed
    # NORMAL - MASL and dropped MASH, contradicting the "MASLD - Normal" label.
    # Select the disease (MASH) and control (NORMAL) arms by name and assert
    # both are present; delta = disease - control to match the stated direction.
    disease_cond <- "MASH"
    control_cond <- "NORMAL"
    missing_conds <- setdiff(c(disease_cond, control_cond), names(cond_conns_list))
    if (length(missing_conds) > 0) {
      stop(sprintf("Cicero differential: required condition(s) [%s] not present; have: %s",
                   paste(missing_conds, collapse = ", "),
                   paste(names(cond_conns_list), collapse = ", ")))
    }
    dis_col <- paste0("coaccess_", disease_cond)
    ctrl_col <- paste0("coaccess_", control_cond)
    diff <- merge(
      cond_conns_list[[disease_cond]][, .(Peak1, Peak2, coaccess)],
      cond_conns_list[[control_cond]][, .(Peak1, Peak2, coaccess)],
      by = c("Peak1", "Peak2"),
      suffixes = c(paste0("_", disease_cond), paste0("_", control_cond)),
      all = TRUE
    )
    for (col in c(dis_col, ctrl_col)) set(diff, which(is.na(diff[[col]])), col, 0)
    diff[, delta_coaccess := get(dis_col) - get(ctrl_col)]  # disease - control
    diff <- diff[order(-abs(delta_coaccess))]
    fwrite(diff, file.path(CICERO_DIR, "differential_coaccessibility.csv"))
    message(sprintf("  Differential connections (%s - %s): %d pairs",
                    disease_cond, control_cond, nrow(diff)))
  }
}

# ==========================================================================
# 8. Generate supplementary figure panels
# ==========================================================================
message("\n8. Generating supplementary figure panels...")

# --- Panel 17: Co-accessibility score distribution ---
message("  Panel 17: Co-accessibility distribution...")
conns_all <- as.data.table(conns)[!is.na(coaccess)]

p17 <- ggplot(conns_all, aes(x = coaccess)) +
  geom_histogram(bins = 100, fill = masld_colors$down, color = NA, alpha = 0.8) +
  geom_vline(xintercept = c(0.1, 0.25), linetype = "dashed",
             linewidth = 0.3, color = c("gray40", masld_colors$up)) +
  labs(x = "Co-accessibility score",
       y = "Peak pairs",
       title = paste0("Cicero co-accessibility (",
                      formatC(nrow(conns_all), big.mark = ","), " pairs)")) +
  theme_masld()

save_panel(p17, "panel_17_coaccess_distribution.pdf", w = fig_half_width, h = 2.5)

# --- Panel 18: Co-accessibility vs distance ---
message("  Panel 18: Co-accessibility vs distance...")
conns_pos <- conns_all[coaccess > 0]
conns_pos[, c("chr1", "s1", "e1") := tstrsplit(Peak1, "_", type.convert = TRUE)]
conns_pos[, c("chr2", "s2", "e2") := tstrsplit(Peak2, "_", type.convert = TRUE)]
conns_pos[, dist_kb := abs((s1+e1)/2 - (s2+e2)/2) / 1000]
conns_pos <- conns_pos[dist_kb > 0 & dist_kb < 500]

p18 <- ggplot(conns_pos[sample(.N, min(.N, 50000))],
              aes(x = dist_kb, y = coaccess)) +
  geom_point(size = 0.1, alpha = 0.1, color = masld_colors$down) +
  geom_smooth(method = "loess", se = TRUE, color = masld_colors$up,
              linewidth = 0.8, span = 0.3) +
  labs(x = "Distance (kb)", y = "Co-accessibility",
       title = "Co-accessibility decays with genomic distance") +
  theme_masld()

save_panel(p18, "panel_18_coaccess_vs_distance.pdf", w = fig_half_width, h = 2.5)

# --- Panel 19: CCAN size distribution ---
message("  Panel 19: CCAN sizes...")
if (nrow(CCANS) > 0 && "CCAN" %in% names(CCANS)) {
  ccan_sizes <- as.data.table(CCANS)[!is.na(CCAN), .N, by = CCAN]
  p19 <- ggplot(ccan_sizes, aes(x = N)) +
    geom_histogram(bins = 30, fill = masld_colors$down, color = NA, alpha = 0.8) +
    labs(x = "Peaks per CCAN", y = "Number of CCANs",
         title = paste0(nrow(ccan_sizes), " cis-co-accessibility networks")) +
    theme_masld()
  save_panel(p19, "panel_19_ccan_sizes.pdf", w = fig_half_width, h = 2.5)
}

# --- Panel 20: Top distal links ---
message("  Panel 20: Top distal links...")
if (nrow(distal_links) > 0) {
  top_links <- distal_links[order(-coaccess)][1:min(25, .N)]
  top_links[, label := paste0(chr1, ":", round(start1/1e6, 2), "-",
                               round(start2/1e6, 2), "Mb")]

  p20 <- ggplot(top_links, aes(x = reorder(label, coaccess), y = coaccess)) +
    geom_col(aes(fill = distance/1000), width = 0.7) +
    scale_fill_gradient(low = masld_colors$down, high = masld_colors$up,
                        name = "Dist (kb)") +
    labs(x = NULL, y = "Co-accessibility",
         title = "Top distal peak pairs (>10kb)") +
    coord_flip() +
    theme_masld() +
    theme(axis.text.y = element_text(size = 4))

  save_panel(p20, "panel_20_top_distal_links.pdf", w = fig_half_width, h = 4)
}

# --- Panel 21: Differential co-accessibility ---
message("  Panel 21: Differential co-accessibility...")
diff_file <- file.path(CICERO_DIR, "differential_coaccessibility.csv")
if (file.exists(diff_file)) {
  diff <- fread(diff_file)
  diff <- diff[!is.na(delta_coaccess)]
  p21 <- ggplot(diff[sample(.N, min(.N, 50000))], aes(x = delta_coaccess)) +
    geom_histogram(bins = 100, fill = masld_colors$down, color = NA, alpha = 0.8) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "gray40") +
    labs(x = expression(Delta~"co-accessibility (MASLD - Normal)"),
         y = "Peak pairs",
         title = "Differential co-accessibility") +
    annotate("text", x = max(diff$delta_coaccess, na.rm=TRUE)*0.5, y = Inf, vjust = 2,
             label = paste0("Gained: ", sum(diff$delta_coaccess > 0.1, na.rm=TRUE)),
             size = 2.2, color = masld_colors$up) +
    annotate("text", x = min(diff$delta_coaccess, na.rm=TRUE)*0.5, y = Inf, vjust = 2,
             label = paste0("Lost: ", sum(diff$delta_coaccess < -0.1, na.rm=TRUE)),
             size = 2.2, color = masld_colors$down) +
    theme_masld()
  save_panel(p21, "panel_21_differential_coaccess.pdf", w = fig_half_width, h = 2.5)
}

cat("\n==============================================================\n")
cat("Cicero analysis complete.\n")
cat(sprintf("Results: %s\n", CICERO_DIR))
cat(sprintf("Figures: %s\n", FIG_OUT))
cat("==============================================================\n")
