#!/usr/bin/env Rscript
# =============================================================================
# 55_ncrna_conservation.R
# Module 4: Synteny-Based lncRNA Conservation Analysis
#
# Assesses evolutionary conservation of lncRNA DEGs via:
# - Syntenic block conservation (flanking PCG ortholog mapping)
# - liftOver coordinate mapping (hg38 -> mm39)
# - PhastCons 100-way multi-species alignment scores
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(biomaRt)
})

select <- dplyr::select
filter <- dplyr::filter

# --- Configuration ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

ncrna_deg_path  <- file.path(BASE, "RNA-seq/results/ncrna/ncrna_deg_annotated.csv")
atlas_path      <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
human_gtf_path  <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
cons_dir        <- file.path(BASE, "data/ncrna_conservation")
out_dir         <- file.path(BASE, "RNA-seq/results/ncrna")

# Data files (must be pre-downloaded)
mouse_gtf_path  <- file.path(cons_dir, "gencode.vM38.annotation.gtf.gz")
chain_path      <- file.path(cons_dir, "hg38ToMm39.over.chain")
phastcons_path  <- file.path(cons_dir, "hg38.phastCons100way.bw")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Module 4: Synteny-Based lncRNA Conservation ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. Check Required Data Files
# =============================================================================
cat("--- 1. Checking required data files ---\n")

required_files <- c(
  "Module 1 output" = ncrna_deg_path,
  "Human GTF" = human_gtf_path,
  "Mouse GTF" = mouse_gtf_path,
  "liftOver chain" = chain_path,
  "PhastCons BigWig" = phastcons_path
)

for (nm in names(required_files)) {
  exists <- file.exists(required_files[nm])
  cat(sprintf("  %-20s %s [%s]\n", nm, basename(required_files[nm]),
              ifelse(exists, "OK", "MISSING")))
}

missing <- !file.exists(required_files)
if (any(missing)) {
  cat("\n  WARNING: Missing files:\n")
  for (nm in names(required_files)[missing]) {
    cat(sprintf("    %s: %s\n", nm, required_files[nm]))
  }
  cat("  Will skip analyses requiring missing files.\n")
}

# =============================================================================
# 2. Load lncRNA DEG Data
# =============================================================================
cat("\n--- 2. Loading ncRNA DEG annotations ---\n")

ncrna <- fread(ncrna_deg_path)
lncrna <- ncrna[gene_biotype == "lncRNA"]
cat(sprintf("  lncRNAs loaded: %d (DEGs: %d)\n", nrow(lncrna),
            sum(lncrna$is_deg, na.rm = TRUE)))

# =============================================================================
# 3. Extract lncRNA Coordinates from Human GTF
# =============================================================================
cat("\n--- 3. Extracting lncRNA coordinates from GENCODE v49 ---\n")

human_gtf <- import(human_gtf_path)
human_genes <- human_gtf[human_gtf$type == "gene"]
human_lnc <- human_genes[human_genes$gene_type == "lncRNA"]
human_pc <- human_genes[human_genes$gene_type == "protein_coding"]

cat(sprintf("  Human lncRNA genes: %d\n", length(human_lnc)))
cat(sprintf("  Human protein-coding genes: %d\n", length(human_pc)))

# Build lncRNA GRanges with gene names
lnc_gr <- human_lnc
names(lnc_gr) <- lnc_gr$gene_name

# Filter to lncRNAs in our dataset
lnc_in_data <- lnc_gr[lnc_gr$gene_name %in% lncrna$human_symbol]
cat(sprintf("  lncRNAs matched to dataset: %d\n", length(lnc_in_data)))

# =============================================================================
# 4. Find Flanking Protein-Coding Genes
# =============================================================================
cat("\n--- 4. Finding flanking protein-coding genes ---\n")

# For each lncRNA, find the nearest upstream and downstream PCG
pc_gr <- human_pc
names(pc_gr) <- pc_gr$gene_name

# Separate by chromosome and strand for proper flanking
flanking_list <- lapply(seq_along(lnc_in_data), function(i) {
  lnc <- lnc_in_data[i]
  lnc_chr <- as.character(seqnames(lnc))
  lnc_start <- start(lnc)
  lnc_end <- end(lnc)

  # PCGs on same chromosome
  pc_same_chr <- pc_gr[seqnames(pc_gr) == lnc_chr]
  if (length(pc_same_chr) == 0) {
    return(data.table(lncrna = lnc$gene_name,
                      upstream_pcg = NA_character_, upstream_dist = NA_integer_,
                      downstream_pcg = NA_character_, downstream_dist = NA_integer_))
  }

  # Upstream: PCGs ending before lncRNA start
  upstream_idx <- which(end(pc_same_chr) < lnc_start)
  if (length(upstream_idx) > 0) {
    dists <- lnc_start - end(pc_same_chr[upstream_idx])
    best_up <- upstream_idx[which.min(dists)]
    up_gene <- pc_same_chr$gene_name[best_up]
    up_dist <- min(dists)
  } else {
    up_gene <- NA_character_
    up_dist <- NA_integer_
  }

  # Downstream: PCGs starting after lncRNA end
  downstream_idx <- which(start(pc_same_chr) > lnc_end)
  if (length(downstream_idx) > 0) {
    dists <- start(pc_same_chr[downstream_idx]) - lnc_end
    best_down <- downstream_idx[which.min(dists)]
    down_gene <- pc_same_chr$gene_name[best_down]
    down_dist <- min(dists)
  } else {
    down_gene <- NA_character_
    down_dist <- NA_integer_
  }

  data.table(lncrna = lnc$gene_name,
             upstream_pcg = up_gene, upstream_dist = as.integer(up_dist),
             downstream_pcg = down_gene, downstream_dist = as.integer(down_dist))
})

flanking_dt <- rbindlist(flanking_list)
cat(sprintf("  lncRNAs with both flanking PCGs: %d / %d\n",
            sum(!is.na(flanking_dt$upstream_pcg) & !is.na(flanking_dt$downstream_pcg)),
            nrow(flanking_dt)))

# =============================================================================
# 5. Map Flanking PCGs to Mouse Orthologs
# =============================================================================
cat("\n--- 5. Mapping flanking PCGs to mouse orthologs ---\n")

all_pcgs <- unique(c(flanking_dt$upstream_pcg, flanking_dt$downstream_pcg))
all_pcgs <- all_pcgs[!is.na(all_pcgs)]

# Use biomaRt with retry logic
get_orthologs <- function(genes, max_retries = 3) {
  for (attempt in 1:max_retries) {
    result <- tryCatch({
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                         host = "https://useast.ensembl.org")
      getLDS(attributes = "hgnc_symbol",
             filters = "hgnc_symbol",
             values = genes,
             mart = ensembl,
             attributesL = c("mgi_symbol"),
             martL = useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                             host = "https://useast.ensembl.org"))
    }, error = function(e) {
      cat(sprintf("  biomaRt attempt %d failed: %s\n", attempt, e$message))
      if (attempt < max_retries) Sys.sleep(5)
      NULL
    })
    if (!is.null(result)) return(result)
  }

  # Fallback: use existing ortholog mapping from atlas
  # W2 fix: atlas mouse_ortholog stores Ensembl IDs (ENSMUSG*), not MGI symbols.
  # Convert to MGI symbols using the mouse GTF gene_id -> gene_name mapping.
  cat("  biomaRt failed, using atlas ortholog mapping as fallback\n")
  atlas_orth <- fread(atlas_path, select = c("human_symbol", "mouse_ortholog"))
  atlas_orth <- atlas_orth[human_symbol %in% genes & !is.na(mouse_ortholog) & mouse_ortholog != ""
                           & mouse_ortholog != '""']

  # Build Ensembl ID -> MGI symbol lookup from mouse GTF
  if (exists("mouse_gtf_path") && file.exists(mouse_gtf_path)) {
    cat("  Building Ensembl->MGI lookup from mouse GTF...\n")
    m_gtf <- import(mouse_gtf_path)
    m_genes <- m_gtf[m_gtf$type == "gene"]
    ens_to_mgi <- data.table(
      ensembl_id = sub("\\..*", "", m_genes$gene_id),
      mgi_symbol = m_genes$gene_name
    )
    ens_to_mgi <- unique(ens_to_mgi)
    cat(sprintf("  Ensembl->MGI lookup: %d entries\n", nrow(ens_to_mgi)))

    # Convert atlas Ensembl IDs to MGI symbols
    atlas_orth <- merge(atlas_orth, ens_to_mgi,
                        by.x = "mouse_ortholog", by.y = "ensembl_id", all.x = TRUE)
    # Use MGI symbol where available, fall back to Ensembl ID
    atlas_orth[, final_symbol := fifelse(!is.na(mgi_symbol) & mgi_symbol != "",
                                          mgi_symbol, mouse_ortholog)]
    n_converted <- sum(!is.na(atlas_orth$mgi_symbol) & atlas_orth$mgi_symbol != "")
    cat(sprintf("  Converted %d / %d Ensembl IDs to MGI symbols\n",
                n_converted, nrow(atlas_orth)))
    result <- atlas_orth[, .(human_symbol, final_symbol)]
    setnames(result, c("HGNC.symbol", "MGI.symbol"))
  } else {
    cat("  WARNING: Mouse GTF not available for Ensembl->MGI conversion\n")
    result <- atlas_orth[, .(human_symbol, mouse_ortholog)]
    setnames(result, c("HGNC.symbol", "MGI.symbol"))
  }
  return(result)
}

orthologs <- get_orthologs(all_pcgs)
if (!is.null(orthologs) && nrow(orthologs) > 0) {
  setDT(orthologs)
  setnames(orthologs, c("human_gene", "mouse_gene"))
  orthologs <- unique(orthologs)
  # Deduplicate: keep one mouse ortholog per human gene (first match)
  orthologs <- orthologs[!duplicated(human_gene)]
  cat(sprintf("  PCGs with mouse orthologs: %d / %d (1:1 after dedup)\n",
              nrow(orthologs), length(all_pcgs)))
} else {
  cat("  WARNING: No orthologs retrieved. Synteny analysis will be limited.\n")
  orthologs <- data.table(human_gene = character(), mouse_gene = character())
}

# Annotate flanking with mouse orthologs (1:1 mapping ensures no row expansion)
flanking_dt <- merge(flanking_dt,
                      orthologs[, .(human_gene, mouse_upstream = mouse_gene)],
                      by.x = "upstream_pcg", by.y = "human_gene", all.x = TRUE)
flanking_dt <- merge(flanking_dt,
                      orthologs[, .(human_gene, mouse_downstream = mouse_gene)],
                      by.x = "downstream_pcg", by.y = "human_gene", all.x = TRUE)

flanking_dt[, has_upstream_ortholog := !is.na(mouse_upstream)]
flanking_dt[, has_downstream_ortholog := !is.na(mouse_downstream)]
cat(sprintf("  flanking_dt after ortholog merge: %d rows (expect %d)\n",
            nrow(flanking_dt), length(lnc_in_data)))

# =============================================================================
# 6. Check Syntenic Conservation with Mouse Genome
# =============================================================================
cat("\n--- 6. Checking syntenic conservation ---\n")

if (file.exists(mouse_gtf_path)) {
  mouse_gtf <- import(mouse_gtf_path)
  mouse_genes <- mouse_gtf[mouse_gtf$type == "gene"]
  mouse_lnc <- mouse_genes[mouse_genes$gene_type == "lncRNA"]
  mouse_pc <- mouse_genes[mouse_genes$gene_type == "protein_coding"]

  cat(sprintf("  Mouse lncRNA genes: %d\n", length(mouse_lnc)))
  cat(sprintf("  Mouse protein-coding genes: %d\n", length(mouse_pc)))

  # For each lncRNA with both flanking orthologs, check if a mouse lncRNA
  # exists between the orthologous flanking genes
  has_both <- flanking_dt[has_upstream_ortholog == TRUE & has_downstream_ortholog == TRUE]
  cat(sprintf("  lncRNAs with both flanking orthologs: %d\n", nrow(has_both)))

  synteny_results <- lapply(seq_len(nrow(has_both)), function(i) {
    row <- has_both[i]
    up_mouse <- row$mouse_upstream
    down_mouse <- row$mouse_downstream

    # Find mouse PCG coordinates
    up_match <- mouse_pc[mouse_pc$gene_name == up_mouse]
    down_match <- mouse_pc[mouse_pc$gene_name == down_mouse]

    if (length(up_match) == 0 | length(down_match) == 0) {
      return(data.table(lncrna = row$lncrna, synteny_status = "No_flanking_match"))
    }

    # Check if on same chromosome
    up_chr <- as.character(seqnames(up_match[1]))
    down_chr <- as.character(seqnames(down_match[1]))

    if (up_chr != down_chr) {
      return(data.table(lncrna = row$lncrna, synteny_status = "Synteny_absent",
                        mouse_lnc_between = NA_character_))
    }

    # Define region between flanking orthologs
    region_start <- min(end(up_match[1]), end(down_match[1]))
    region_end <- max(start(up_match[1]), start(down_match[1]))

    if (region_start >= region_end) {
      # Genes overlap or are inverted
      return(data.table(lncrna = row$lncrna, synteny_status = "Synteny_partial",
                        mouse_lnc_between = NA_character_))
    }

    region_gr <- GRanges(up_chr, IRanges(region_start, region_end))

    # Find mouse lncRNAs in this region
    overlaps <- findOverlaps(mouse_lnc, region_gr)
    if (length(overlaps) > 0) {
      mouse_lncs <- mouse_lnc$gene_name[queryHits(overlaps)]
      return(data.table(lncrna = row$lncrna, synteny_status = "Synteny_conserved",
                        mouse_lnc_between = paste(unique(mouse_lncs), collapse = ";")))
    } else {
      return(data.table(lncrna = row$lncrna, synteny_status = "Synteny_absent",
                        mouse_lnc_between = NA_character_))
    }
  })

  synteny_dt <- rbindlist(synteny_results, fill = TRUE)

  # Merge back to flanking table
  flanking_dt <- merge(flanking_dt, synteny_dt, by = "lncrna", all.x = TRUE)
  flanking_dt[is.na(synteny_status), synteny_status := "No_flanking_ortholog"]

  cat("  Synteny classification:\n")
  print(flanking_dt[, .N, by = synteny_status][order(-N)])
} else {
  cat("  SKIPPED: Mouse GTF not found\n")
  flanking_dt[, synteny_status := "Not_assessed"]
  flanking_dt[, mouse_lnc_between := NA_character_]
}

# =============================================================================
# 7. liftOver Analysis
# =============================================================================
cat("\n--- 7. liftOver coordinate mapping ---\n")

if (file.exists(chain_path)) {
  chain <- import.chain(chain_path)

  # liftOver lncRNA coordinates
  lnc_for_lift <- lnc_in_data
  lifted <- liftOver(lnc_for_lift, chain)

  # Compute liftover success rate
  lift_results <- data.table(
    lncrna = lnc_for_lift$gene_name,
    n_human_ranges = 1,
    n_mouse_ranges = lengths(lifted),
    liftover_success = lengths(lifted) > 0
  )

  # For successful lifts, compute fraction of bases mapped
  for (i in seq_along(lifted)) {
    if (length(lifted[[i]]) > 0) {
      human_width <- width(lnc_for_lift[i])
      mouse_width <- sum(width(lifted[[i]]))
      lift_results$frac_bases_mapped[i] <- min(mouse_width / human_width, 1.0)
    } else {
      lift_results$frac_bases_mapped[i] <- 0
    }
  }

  cat(sprintf("  liftOver success: %d / %d (%.1f%%)\n",
              sum(lift_results$liftover_success),
              nrow(lift_results),
              100 * mean(lift_results$liftover_success)))
  cat(sprintf("  Median fraction bases mapped (successful): %.3f\n",
              median(lift_results$frac_bases_mapped[lift_results$liftover_success], na.rm = TRUE)))

  # Merge liftover results
  flanking_dt <- merge(flanking_dt, lift_results[, .(lncrna, liftover_success, frac_bases_mapped)],
                        by = "lncrna", all.x = TRUE)
} else {
  cat("  SKIPPED: liftOver chain not found\n")
  flanking_dt[, liftover_success := NA]
  flanking_dt[, frac_bases_mapped := NA_real_]
}

# =============================================================================
# 8. PhastCons Conservation Scores
# =============================================================================
cat("\n--- 8. PhastCons 100-way conservation scores ---\n")

if (file.exists(phastcons_path)) {
  bw <- import(phastcons_path, format = "BigWig",
               which = GRanges(seqlevels(lnc_in_data)[seqlevels(lnc_in_data) %in%
                                paste0("chr", c(1:22, "X", "Y"))],
                               IRanges(1, 1)))

  # Score at TSS +/- 500bp
  tss_gr <- resize(lnc_in_data, width = 1, fix = "start")
  tss_window <- resize(tss_gr, width = 1001, fix = "center")
  # Constrain to valid chromosomes
  valid_chr <- paste0("chr", c(1:22, "X", "Y"))
  tss_window <- tss_window[as.character(seqnames(tss_window)) %in% valid_chr]

  # Import PhastCons scores in the TSS windows
  cat("  Computing PhastCons scores at TSS +/- 500bp...\n")
  phastcons_scores <- tryCatch({
    scores <- import(phastcons_path, format = "BigWig", which = tss_window)
    # Average score per lncRNA
    olaps <- findOverlaps(tss_window, scores)
    score_by_gene <- data.table(
      idx = queryHits(olaps),
      score = scores$score[subjectHits(olaps)]
    )
    mean_scores <- score_by_gene[, .(phastcons_tss = mean(score, na.rm = TRUE)), by = idx]
    result <- data.table(
      lncrna = tss_window$gene_name[mean_scores$idx],
      phastcons_tss = mean_scores$phastcons_tss
    )
    result
  }, error = function(e) {
    cat(sprintf("  WARNING: PhastCons scoring failed: %s\n", e$message))
    data.table(lncrna = character(), phastcons_tss = numeric())
  })

  if (nrow(phastcons_scores) > 0) {
    cat(sprintf("  PhastCons scored: %d lncRNAs\n", nrow(phastcons_scores)))
    cat(sprintf("  Median PhastCons at TSS: %.4f\n",
                median(phastcons_scores$phastcons_tss, na.rm = TRUE)))

    # Merge
    flanking_dt <- merge(flanking_dt, phastcons_scores, by = "lncrna", all.x = TRUE)
  } else {
    flanking_dt[, phastcons_tss := NA_real_]
  }

  # Also score gene body
  cat("  Computing PhastCons scores at gene body...\n")
  gene_body_scores <- tryCatch({
    lnc_body <- lnc_in_data[as.character(seqnames(lnc_in_data)) %in% valid_chr]
    scores <- import(phastcons_path, format = "BigWig", which = lnc_body)
    olaps <- findOverlaps(lnc_body, scores)
    score_by_gene <- data.table(
      idx = queryHits(olaps),
      score = scores$score[subjectHits(olaps)]
    )
    mean_scores <- score_by_gene[, .(phastcons_body = mean(score, na.rm = TRUE)), by = idx]
    data.table(
      lncrna = lnc_body$gene_name[mean_scores$idx],
      phastcons_body = mean_scores$phastcons_body
    )
  }, error = function(e) {
    cat(sprintf("  WARNING: Gene body scoring failed: %s\n", e$message))
    data.table(lncrna = character(), phastcons_body = numeric())
  })

  if (nrow(gene_body_scores) > 0) {
    flanking_dt <- merge(flanking_dt, gene_body_scores, by = "lncrna", all.x = TRUE)
    cat(sprintf("  Median PhastCons gene body: %.4f\n",
                median(flanking_dt$phastcons_body, na.rm = TRUE)))
  } else {
    flanking_dt[, phastcons_body := NA_real_]
  }
} else {
  cat("  SKIPPED: PhastCons BigWig not found\n")
  flanking_dt[, phastcons_tss := NA_real_]
  flanking_dt[, phastcons_body := NA_real_]
}

# =============================================================================
# 9. Merge DEG status, Deduplicate, and Save
# =============================================================================
cat("\n--- 9. Merging DEG status and saving ---\n")

# Merge DEG status — deduplicate lncrna data first to prevent cartesian product
lncrna_dedup <- lncrna[!duplicated(human_symbol), .(human_symbol, is_deg, bulk_logFC, bulk_padj)]
cat(sprintf("  DEG data: %d unique lncRNAs (from %d rows)\n", nrow(lncrna_dedup), nrow(lncrna)))
flanking_dt <- merge(flanking_dt, lncrna_dedup,
                      by.x = "lncrna", by.y = "human_symbol", all.x = TRUE)

# Deduplicate by lncrna before saving — merges in steps 5-8 may create
# multiple rows per lncRNA from 1:many ortholog/liftOver/PhastCons mappings.
# Keep the row with the best synteny status per lncRNA.
if (anyDuplicated(flanking_dt$lncrna) > 0) {
  n_before <- nrow(flanking_dt)
  # Priority: Synteny_conserved > Synteny_partial > Synteny_absent > No_flanking_match > No_flanking_ortholog
  status_order <- c("Synteny_conserved" = 1, "Synteny_partial" = 2,
                     "Synteny_absent" = 3, "No_flanking_match" = 4,
                     "No_flanking_ortholog" = 5, "Not_assessed" = 6)
  flanking_dt[, status_rank := status_order[synteny_status]]
  flanking_dt[is.na(status_rank), status_rank := 99]
  setorder(flanking_dt, lncrna, status_rank)
  flanking_dt <- flanking_dt[!duplicated(lncrna)]
  flanking_dt[, status_rank := NULL]
  cat(sprintf("  Deduplicated flanking_dt: %d -> %d rows\n", n_before, nrow(flanking_dt)))
}

# Full conservation table
fwrite(flanking_dt, file.path(out_dir, "lncrna_synteny_conservation.csv"))

# Summary table
cons_summary <- flanking_dt[, .(
  n_total = .N,
  n_with_both_flanking = sum(!is.na(upstream_pcg) & !is.na(downstream_pcg)),
  n_synteny_conserved = sum(synteny_status == "Synteny_conserved", na.rm = TRUE),
  n_synteny_partial = sum(synteny_status == "Synteny_partial", na.rm = TRUE),
  n_synteny_absent = sum(synteny_status == "Synteny_absent", na.rm = TRUE),
  n_liftover_success = sum(liftover_success == TRUE, na.rm = TRUE),
  median_phastcons_tss = median(phastcons_tss, na.rm = TRUE),
  median_phastcons_body = median(phastcons_body, na.rm = TRUE),
  n_deg = sum(is_deg, na.rm = TRUE),
  n_deg_conserved = sum(is_deg == TRUE & synteny_status == "Synteny_conserved", na.rm = TRUE)
)]

fwrite(cons_summary, file.path(out_dir, "lncrna_conservation_summary.csv"))

cat(sprintf("  Saved lncrna_synteny_conservation.csv (%d rows)\n", nrow(flanking_dt)))
cat(sprintf("  Saved lncrna_conservation_summary.csv\n"))

# =============================================================================
# 11. Summary
# =============================================================================
cat("\n--- 11. Summary ---\n")
cat(sprintf("  lncRNAs analyzed: %d\n", nrow(flanking_dt)))
cat(sprintf("  Synteny-conserved: %d (%.1f%%)\n",
            cons_summary$n_synteny_conserved,
            100 * cons_summary$n_synteny_conserved / max(cons_summary$n_total, 1)))
cat(sprintf("  liftOver success: %d (%.1f%%)\n",
            cons_summary$n_liftover_success,
            100 * cons_summary$n_liftover_success / max(cons_summary$n_total, 1)))
cat(sprintf("  Median PhastCons (TSS): %.4f\n", cons_summary$median_phastcons_tss))

# =============================================================================
# 11. Enrichment Tests (after save — non-critical)
# =============================================================================
cat("\n--- 11. Conservation enrichment among DEGs ---\n")

if ("synteny_status" %in% names(flanking_dt) && any(flanking_dt$synteny_status == "Synteny_conserved")) {
  deg_syn <- flanking_dt[!is.na(is_deg)]
  mat <- matrix(c(
    sum(deg_syn$is_deg == TRUE & deg_syn$synteny_status == "Synteny_conserved"),
    sum(deg_syn$is_deg == FALSE & deg_syn$synteny_status == "Synteny_conserved"),
    sum(deg_syn$is_deg == TRUE & deg_syn$synteny_status != "Synteny_conserved"),
    sum(deg_syn$is_deg == FALSE & deg_syn$synteny_status != "Synteny_conserved")
  ), nrow = 2)
  fisher_syn <- fisher.test(mat)
  cat(sprintf("  Synteny-conserved enrichment in DEGs: OR=%.2f, p=%.4f\n",
              fisher_syn$estimate, fisher_syn$p.value))
}

if ("phastcons_tss" %in% names(flanking_dt)) {
  deg_scores <- flanking_dt[is_deg == TRUE & !is.na(phastcons_tss)]$phastcons_tss
  nondeg_scores <- flanking_dt[is_deg == FALSE & !is.na(phastcons_tss)]$phastcons_tss
  if (length(deg_scores) > 5 & length(nondeg_scores) > 5) {
    wt <- wilcox.test(deg_scores, nondeg_scores)
    cat(sprintf("  PhastCons TSS: DEG median=%.4f, non-DEG median=%.4f, p=%.4f\n",
                median(deg_scores), median(nondeg_scores), wt$p.value))
  }
}

if ("liftover_success" %in% names(flanking_dt)) {
  lo_dt <- flanking_dt[!is.na(liftover_success) & !is.na(is_deg)]
  if (nrow(lo_dt) > 0) {
    mat_lo <- matrix(c(
      sum(lo_dt$is_deg == TRUE & lo_dt$liftover_success == TRUE),
      sum(lo_dt$is_deg == FALSE & lo_dt$liftover_success == TRUE),
      sum(lo_dt$is_deg == TRUE & lo_dt$liftover_success == FALSE),
      sum(lo_dt$is_deg == FALSE & lo_dt$liftover_success == FALSE)
    ), nrow = 2)
    fisher_lo <- fisher.test(mat_lo)
    cat(sprintf("  liftOver success enrichment in DEGs: OR=%.2f, p=%.4f\n",
                fisher_lo$estimate, fisher_lo$p.value))
  }
}

cat("\n=== Module 4 Complete ===\n")
cat("End:", format(Sys.time()), "\n")
