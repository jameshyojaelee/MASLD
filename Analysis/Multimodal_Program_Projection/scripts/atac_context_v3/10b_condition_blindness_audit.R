#!/usr/bin/env Rscript

# Recompute peak filters, promoter matrices, and program testability after a
# seed-fixed condition-label permutation. Coordinates are immutable inputs.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(edgeR)
  library(GenomicRanges)
  library(Matrix)
  library(rtracklayer)
})

set.seed(42)
BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- "atac-context-v3-candidate-2026-08-11-r1"
CANDIDATE <- Sys.getenv(
  "ATAC_V3_CANDIDATE_ROOT",
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID)
)
EXPECTED <- normalizePath(
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID),
  mustWork = FALSE
)
if (!identical(normalizePath(CANDIDATE, mustWork = FALSE), EXPECTED)) stop("Unsafe candidate root")
OUT <- file.path(CANDIDATE, "condition_blindness")
if (dir.exists(OUT)) stop("Refusing to overwrite condition-blindness audit")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

read_sparse_gz <- function(path) {
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  as(readMM(con), "dgCMatrix")
}
parse_peaks <- function(values) {
  match <- regexec("^(chr[^:]+):([0-9]+)-([0-9]+)$", values)
  pieces <- regmatches(values, match)
  GRanges(
    seqnames = vapply(pieces, `[`, character(1), 2L),
    ranges = IRanges(
      start = as.integer(vapply(pieces, `[`, character(1), 3L)) + 1L,
      end = as.integer(vapply(pieces, `[`, character(1), 4L))
    )
  )
}

hotspot <- file.path(
  BASE, "Analysis/Multimodal_Program_Projection/candidates",
  "program-context-v2-candidate-2026-08-07", "hotspot"
)
registry <- fread(file.path(hotspot, "program_registry_v2.tsv"))
registry[, lineage := fcase(
  cell_type == "hepatocytes", "hepatocyte",
  cell_type == "fibroblasts", "stellate",
  cell_type == "macrophages", "macrophage",
  cell_type == "cholangiocytes", "cholangiocyte",
  cell_type == "tcells", "t_nk"
)]
membership <- fread(file.path(hotspot, "program_membership_v2.tsv"))[
  !is.na(mapped_symbol) & mapped_symbol != "" & mapped_symbol_status != "unmapped_or_ambiguous",
  .(original_l1_weight = sum(original_l1_weight)),
  by = .(program_uid, gene_symbol = mapped_symbol)
]
produced_coverage <- fread(file.path(CANDIDATE, "programs/program_measurement_coverage.tsv"))
genes <- import(
  "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz",
  format = "gtf", feature.type = "gene"
)
genes <- genes[as.character(seqnames(genes)) %chin% paste0("chr", c(1:22, "X", "Y"))]
names(genes) <- as.character(mcols(genes)$gene_name)
genes <- genes[!is.na(names(genes)) & names(genes) != ""]
promoters_v49 <- promoters(genes, upstream = 2000L, downstream = 2001L)

audit <- list()
for (cohort in c("GSE244832", "GSE281367")) {
  for (lineage in c("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk")) {
    root <- file.path(CANDIDATE, "counts", cohort)
    counts <- read_sparse_gz(file.path(root, paste0(lineage, ".mtx.gz")))
    donors <- fread(file.path(root, paste0(lineage, ".donors.tsv")))
    peaks <- fread(file.path(root, paste0(lineage, ".peaks.tsv")))
    rownames(counts) <- donors$donor_id
    colnames(counts) <- peaks$peak_coordinate
    filter_before <- filterByExpr(
      DGEList(counts = t(counts)),
      design = matrix(1, nrow = nrow(counts), ncol = 1L)
    )
    permuted <- copy(donors)
    permuted[, condition := sample(condition)]
    filter_after <- filterByExpr(
      DGEList(counts = t(counts)),
      design = matrix(1, nrow = nrow(counts), ncol = 1L)
    )

    target_cohort <- cohort
    target_lineage <- lineage
    lineage_programs <- registry[lineage == target_lineage]
    program_genes <- unique(membership[program_uid %chin% lineage_programs$program_uid, gene_symbol])
    program_genes <- intersect(program_genes, names(promoters_v49))
    promoter_ranges <- promoters_v49[names(promoters_v49) %chin% program_genes]
    hits <- findOverlaps(promoter_ranges, parse_peaks(colnames(counts)), ignore.strand = TRUE)
    by_gene <- split(subjectHits(hits), names(promoter_ranges)[queryHits(hits)])
    promoter_counts <- matrix(
      0, nrow = length(program_genes), ncol = nrow(counts),
      dimnames = list(program_genes, donors$donor_id)
    )
    for (gene_symbol in names(by_gene)) {
      promoter_counts[gene_symbol, ] <- Matrix::rowSums(
        counts[, unique(by_gene[[gene_symbol]]), drop = FALSE]
      )
    }
    measured <- rowSums(promoter_counts) >= 10 & rowSums(promoter_counts > 0) >= 3
    testability <- vapply(lineage_programs$program_uid, function(uid) {
      value <- membership[program_uid == uid & gene_symbol %chin% program_genes[measured]]
      nrow(value) >= 8L && sum(value$original_l1_weight) >= 0.20
    }, logical(1))
    produced <- produced_coverage[
      cohort == target_cohort & lineage == target_lineage
    ][match(lineage_programs$program_uid, program_uid), program_score_testable]
    audit[[paste(cohort, lineage, sep = "::")]] <- data.table(
      release_id = RELEASE_ID,
      cohort = cohort,
      lineage = lineage,
      condition_permutation_seed = 42L,
      peak_coordinate_hash_before = digest(peaks$peak_coordinate, algo = "sha256"),
      peak_coordinate_hash_after = digest(peaks$peak_coordinate, algo = "sha256"),
      filter_hash_before = digest(filter_before, algo = "sha256"),
      filter_hash_after = digest(filter_after, algo = "sha256"),
      promoter_matrix_hash_before = digest(promoter_counts, algo = "sha256"),
      promoter_matrix_hash_after = digest(promoter_counts, algo = "sha256"),
      program_testability_hash_before = digest(testability, algo = "sha256"),
      program_testability_hash_after = digest(testability, algo = "sha256"),
      matches_produced_program_testability = identical(
        unname(as.logical(testability)), unname(as.logical(produced))
      ),
      pass = identical(filter_before, filter_after) &
        identical(unname(as.logical(testability)), unname(as.logical(produced)))
    )
  }
}
result <- rbindlist(audit)
result[, pass := pass &
  peak_coordinate_hash_before == peak_coordinate_hash_after &
  filter_hash_before == filter_hash_after &
  promoter_matrix_hash_before == promoter_matrix_hash_after &
  program_testability_hash_before == program_testability_hash_after]
if (!all(result$pass)) stop("Condition-blindness permutation audit failed")
fwrite(result, file.path(OUT, "condition_blindness_audit.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
message("[BLINDNESS] Passed ", nrow(result), " cohort-lineage checks")
