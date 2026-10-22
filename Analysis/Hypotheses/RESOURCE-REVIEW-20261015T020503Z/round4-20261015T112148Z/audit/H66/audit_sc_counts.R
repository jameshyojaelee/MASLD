# H66 audit: recount single-cell runs, donors and eligible donors per lineage.
# Reads pseudobulk column sums (library size) only; no per-gene detection is computed.
suppressPackageStartupMessages(library(data.table))
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
m <- build_srr_to_donor_map(BASE)
cc <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/lineage_cell_counts.tsv"))
meta <- fread(file.path(BASE, "Analysis/SingleCell/metadata/snrna_donor_sex.csv"))[, .(sample, dataset)]
cc <- merge(cc, meta, by = "sample", all.x = TRUE)
cat("runs", nrow(cc), "unique runs", uniqueN(cc$sample), "runs without dataset", sum(is.na(cc$dataset)), "\n")
cat("GSE289173 in cell-count datasets:", any(grepl("GSE289173", cc$dataset)), "\n")
cc[, donor := ifelse(sample %in% names(m), m[sample], sample)]
cat("donors", uniqueN(cc$donor), "\n")
print(cc[, .(donors = uniqueN(donor), runs = .N), by = dataset][order(-donors)])
lins <- c("B_cells", "Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Hepatocytes",
          "Macrophages", "Plasma_cells", "Resident_NK", "T_cells")
out <- list()
for (L in lins) {
  f <- file.path(BASE, sprintf("Analysis/SingleCell/results_gpu_v2/pseudobulk/%s_pseudobulk.csv", L))
  pb <- fread(f)
  genes <- pb[[1]]
  v <- as.matrix(pb[, -1])
  runs <- colnames(v)
  lib <- colSums(v)
  in_cc <- runs %in% cc$sample
  ds <- meta$dataset[match(runs, meta$sample)]
  donor_of_run <- ifelse(runs %in% names(m), m[runs], runs)
  dlib <- tapply(lib, donor_of_run, sum)
  col <- paste0("n_", L)
  # cells summed over every run in the cell-count table (prespec text)
  cells_all <- cc[, .(cells = sum(get(col))), by = donor]
  # cells summed only over runs that are columns of this lineage's pseudobulk
  cells_pb <- cc[sample %in% runs, .(cells = sum(get(col))), by = donor]
  d <- data.table(donor = names(dlib), lib = as.numeric(dlib))
  d <- merge(d, cells_all, by = "donor", all.x = TRUE)
  setnames(d, "cells", "cells_all")
  d <- merge(d, cells_pb, by = "donor", all.x = TRUE)
  setnames(d, "cells", "cells_pbruns")
  d[is.na(cells_all), cells_all := 0]; d[is.na(cells_pbruns), cells_pbruns := 0]
  runs_cc_pos <- cc[get(col) > 0, sample]
  out[[L]] <- data.table(
    lineage = L, pb_genes = length(genes), genes_unique = uniqueN(genes),
    pb_cols = length(runs), cols_in_cc = sum(in_cc), cols_gse289173 = sum(grepl("GSE289173", ds), na.rm = TRUE),
    cols_no_dataset = sum(is.na(ds)), integer_valued = all(abs(v - round(v)) < 1e-6),
    pb_donors = length(dlib),
    cc_runs_cells_gt0 = length(runs_cc_pos), cc_runs_gt0_not_in_pb = sum(!runs_cc_pos %in% runs),
    pb_cols_with_0_cells = sum(cc[match(runs, sample)][[col]] == 0, na.rm = TRUE),
    donors_ge20_allruns = sum(d$cells_all >= 20),
    donors_ge20_pbruns = sum(d$cells_pbruns >= 20),
    donors_lib_ge1e4 = sum(d$lib >= 1e4),
    elig_text_rule = sum(d$cells_all >= 20 & d$lib >= 1e4),
    elig_pbruns_rule = sum(d$cells_pbruns >= 20 & d$lib >= 1e4),
    elig_pbruns_50cells = sum(d$cells_pbruns >= 50 & d$lib >= 1e4),
    min_lib_eligible = min(d$lib[d$cells_pbruns >= 20 & d$lib >= 1e4])
  )
  # donors per dataset among eligible donors (pbruns rule)
  dd <- unique(cc[, .(donor, dataset)])
  e <- merge(d[cells_pbruns >= 20 & lib >= 1e4], dd, by = "donor", all.x = TRUE)
  cat(L, "eligible donors by dataset:", paste(names(table(e$dataset)), table(e$dataset), sep = "=", collapse = " "), "\n")
  rm(pb, v); invisible(gc())
}
res <- rbindlist(out)
print(res, width = 250)
fwrite(res, "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H66/audit_sc_counts.tsv", sep = "\t")
# gene name format check
pb <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk/Hepatocytes_pseudobulk.csv"), select = 1)
g <- pb[[1]]
cat("pseudobulk rows", length(g), "ENSG-like", sum(grepl("^ENSG", g)), "examples", head(g, 5), "\n")
writeLines(capture.output(sessionInfo()), "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H66/sessionInfo_audit_sc.txt")
