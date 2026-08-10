#!/usr/bin/env Rscript

# Two-layer ATAC projection for Figure 4.
#
# Layer 1 (static, descriptive): available credible-set rows from the
# prespecified 35-study Tier-1/2 registry that lie in lineage-matched scATAC
# peaks and within +/-2 kb of a frozen program-gene TSS. Direct MASLD/PDFF and
# liver-enzyme loci are retained as separate scopes; no stale mixed-ancestry
# enrichment null is reused.
#
# Layer 2 (dynamic): donor-pseudobulk promoter accessibility scores for the 22
# frozen programs in external GSE281367 (primary) and GSE244832 multiome
# (internal cross-modal consistency; the donor contrast itself is unpaired).
# Disease labels are never used for feature
# filtering or standardization. Exact donor-label permutations provide p values.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(GenomicRanges)
  library(rtracklayer)
})
set.seed(42)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
ROOT <- file.path(BASE, "Analysis/Multimodal_Program_Projection")
SMOKE <- toupper(Sys.getenv("FIG4_ATAC_SMOKE", "FALSE")) %in% c("TRUE", "1", "YES")
OUT <- if (SMOKE) file.path(ROOT, "results/smoke/atac") else file.path(ROOT, "results/atac")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

REGISTRY_FILE <- file.path(ROOT, "results/frozen_programs.tsv")
MEMBERSHIP_FILE <- file.path(ROOT, "results/frozen_program_membership.tsv")
GTF_FILE <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
G244 <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/snapatac2")
G244_META <- file.path(BASE, "Analysis/ATAC/Human_Multiome/metadata/donor_metadata_curated.tsv")
G281 <- file.path(BASE, "Analysis/ATAC/Human_External/pseudobulk")
PEAK_DIR <- file.path(
  BASE, "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2"
)
CREDIBLE_FILE <- file.path(BASE, "GWAS/finemapping/results/credible_sets.csv")
TIER_FILE <- file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv")
CHAIN_FILE <- file.path(BASE, "data/broadaway_eqtl/hg19ToHg38.over.chain")

required <- c(
  REGISTRY_FILE, MEMBERSHIP_FILE, GTF_FILE, G244_META, CREDIBLE_FILE,
  TIER_FILE, CHAIN_FILE
)
if (any(!file.exists(required))) {
  stop("Missing ATAC input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

registry <- fread(REGISTRY_FILE)
membership <- fread(MEMBERSHIP_FILE)
membership <- membership[
  mapped_symbol == TRUE & !is.na(gene_symbol) & gene_symbol != "",
  .(original_l1_weight = sum(original_l1_weight)),
  by = .(program_id, cell_type, module, gene_symbol)
]

ct_map <- data.table(
  cell_type = c("hepatocytes", "fibroblasts", "macrophages", "cholangiocytes"),
  atac_prefix = c("hep", "stellate", "macrophage", "cholangiocyte"),
  peak_label = c("Hepatocytes", "Fibroblasts", "Macrophages", "Cholangiocytes")
)

message("[ATAC] importing GENCODE v49 gene coordinates")
genes <- import(GTF_FILE, format = "gtf", feature.type = "gene")
standard_chr <- paste0("chr", c(1:22, "X", "Y"))
genes <- genes[as.character(seqnames(genes)) %in% standard_chr]
gene_symbol <- as.character(mcols(genes)$gene_name)
names(genes) <- gene_symbol
genes <- genes[!is.na(names(genes)) & names(genes) != ""]
genes <- genes[!duplicated(names(genes))]
promoter_gr <- promoters(genes, upstream = 2000L, downstream = 2001L)

fast_read_counts <- function(path) {
  con <- gzfile(path, "r")
  lines <- readLines(con)
  close(con)
  peaks <- strsplit(lines[1], "\t", fixed = TRUE)[[1]][-1]
  body <- lines[-1]
  mat <- matrix(0, nrow = length(body), ncol = length(peaks))
  donors <- character(length(body))
  for (i in seq_along(body)) {
    fields <- strsplit(body[i], "\t", fixed = TRUE)[[1]]
    donors[i] <- fields[1]
    mat[i, ] <- as.numeric(fields[-1])
  }
  rownames(mat) <- donors
  colnames(mat) <- peaks
  mat
}

parse_peaks <- function(x) {
  mt <- regexec("^(chr[^:]+):([0-9]+)-([0-9]+)$", x)
  bits <- regmatches(x, mt)
  good <- lengths(bits) == 4L
  out <- data.table(
    peak_index = which(good),
    peak_name = x[good],
    chr = vapply(bits[good], `[`, character(1), 2),
    start0 = as.integer(vapply(bits[good], `[`, character(1), 3)),
    end = as.integer(vapply(bits[good], `[`, character(1), 4))
  )
  out[, start := start0 + 1L]
  out
}

exact_two_group_p <- function(score, case) {
  score <- as.numeric(score)
  case <- as.logical(case)
  ok <- is.finite(score) & !is.na(case)
  score <- score[ok]
  case <- case[ok]
  n <- length(score)
  n_case <- sum(case)
  if (n < 6L || n_case == 0L || n_case == n) return(c(effect = NA_real_, pvalue = NA_real_))
  observed <- mean(score[case]) - mean(score[!case])
  cmb <- combn(n, n_case)
  total <- sum(score)
  case_means <- colSums(matrix(score[cmb], nrow = n_case)) / n_case
  perm <- case_means - (total - case_means * n_case) / (n - n_case)
  p <- mean(abs(perm) >= abs(observed) - 1e-12)
  c(effect = observed, pvalue = p)
}

score_module <- function(z, m, samples, case) {
  measured <- m[gene_symbol %in% rownames(z)]
  n_measured <- nrow(measured)
  retained <- sum(measured$original_l1_weight)
  testable <- n_measured >= 8L && retained >= 0.20
  if (!testable) {
    return(list(
      summary = data.table(
        n_measured = n_measured, retained_l1_weight = retained,
        testable = FALSE, effect = NA_real_, pvalue = NA_real_,
        equal_effect = NA_real_, leave_top_effect = NA_real_
      ),
      scores = NULL
    ))
  }
  measured[, weight := original_l1_weight / sum(original_l1_weight)]
  zz <- z[measured$gene_symbol, samples, drop = FALSE]
  observed_weight <- colSums(is.finite(zz) * measured$weight)
  weighted <- colSums(zz * measured$weight, na.rm = TRUE) / observed_weight
  weighted[!is.finite(weighted)] <- NA_real_
  equal <- colMeans(zz, na.rm = TRUE)
  top_gene <- measured$gene_symbol[which.max(measured$weight)]
  leave <- measured[gene_symbol != top_gene]
  leave[, weight := original_l1_weight / sum(original_l1_weight)]
  leave_zz <- z[leave$gene_symbol, samples, drop = FALSE]
  leave_observed_weight <- colSums(is.finite(leave_zz) * leave$weight)
  leave_score <- colSums(leave_zz * leave$weight, na.rm = TRUE) / leave_observed_weight
  leave_score[!is.finite(leave_score)] <- NA_real_
  primary <- exact_two_group_p(weighted, case)
  eq <- exact_two_group_p(equal, case)
  lo <- exact_two_group_p(leave_score, case)
  list(
    summary = data.table(
      n_measured = n_measured,
      retained_l1_weight = retained,
      testable = TRUE,
      effect = unname(primary["effect"]),
      pvalue = unname(primary["pvalue"]),
      equal_effect = unname(eq["effect"]),
      leave_top_effect = unname(lo["effect"])
    ),
    scores = data.table(sample_id = samples, score = weighted)
  )
}

dynamic_results <- list()
dynamic_scores <- list()
coverage_rows <- list()

cohorts_to_run <- if (SMOKE) "GSE281367" else c("GSE281367", "GSE244832")
ct_to_run <- if (SMOKE) ct_map[1] else ct_map
for (cohort in cohorts_to_run) {
  message("[ATAC] dynamic promoter accessibility: ", cohort)
  for (j in seq_len(nrow(ct_to_run))) {
    ct <- ct_to_run$cell_type[j]
    prefix <- ct_to_run$atac_prefix[j]
    if (cohort == "GSE281367") {
      counts_file <- file.path(G281, sprintf("%s_pseudobulk_counts_GSE281367.tsv.gz", prefix))
      coldata_file <- file.path(G281, sprintf("%s_pseudobulk_coldata_GSE281367.tsv", prefix))
      if (!file.exists(counts_file) || !file.exists(coldata_file)) next
      coldata <- fread(coldata_file)
      coldata[, condition_curated := toupper(condition)]
      coldata <- coldata[condition_curated %in% c("NORMAL", "MASH")]
    } else {
      counts_file <- file.path(G244, sprintf("%s_pseudobulk_counts.tsv.gz", prefix))
      if (!file.exists(counts_file)) next
      coldata <- fread(G244_META, select = c("donor_id", "condition"))
      coldata[, condition_curated := toupper(condition)]
      coldata <- coldata[condition_curated %in% c("NORMAL", "MASH")]
    }

    counts <- fast_read_counts(counts_file)
    counts <- counts[, !duplicated(colnames(counts)), drop = FALSE]
    samples <- intersect(coldata$donor_id, rownames(counts))
    coldata <- coldata[match(samples, donor_id)]
    counts <- counts[samples, , drop = FALSE]
    if (cohort == "GSE281367" && !identical(c(sum(coldata$condition_curated == "NORMAL"), sum(coldata$condition_curated == "MASH")), c(6L, 6L))) {
      stop("GSE281367 must contain 6 NORMAL and 6 MASH donors for ", ct)
    }
    if (cohort == "GSE244832" && !identical(c(sum(coldata$condition_curated == "NORMAL"), sum(coldata$condition_curated == "MASH")), c(5L, 9L))) {
      stop("GSE244832 must contain 5 NORMAL and 9 MASH donors for ", ct)
    }

    peaks <- parse_peaks(colnames(counts))
    counts <- counts[, peaks$peak_index, drop = FALSE]
    peak_gr <- GRanges(
      seqnames = peaks$chr,
      ranges = IRanges(start = peaks$start, end = peaks$end)
    )
    names(peak_gr) <- peaks$peak_name

    genes_ct <- unique(membership[cell_type == ct, gene_symbol])
    genes_ct <- intersect(genes_ct, names(promoter_gr))
    prom <- promoter_gr[genes_ct]
    hits <- findOverlaps(prom, peak_gr, ignore.strand = TRUE)
    by_gene <- split(subjectHits(hits), queryHits(hits))
    gene_counts <- matrix(
      0,
      nrow = length(genes_ct),
      ncol = length(samples),
      dimnames = list(genes_ct, samples)
    )
    for (idx in names(by_gene)) {
      peak_idx <- unique(by_gene[[idx]])
      gene_counts[as.integer(idx), ] <- rowSums(counts[, peak_idx, drop = FALSE])
    }

    # Group-blind measured-gene rule.
    measured <- rowSums(gene_counts) >= 10 & rowSums(gene_counts > 0) >= 3
    gene_counts <- gene_counts[measured, , drop = FALSE]
    dge <- DGEList(counts = t(counts))
    dge <- calcNormFactors(dge, method = "TMM")
    eff_lib <- dge$samples$lib.size * dge$samples$norm.factors
    logcpm <- cpm(gene_counts, lib.size = eff_lib, log = TRUE, prior.count = 0.5)
    z <- t(scale(t(logcpm)))
    z[!is.finite(z)] <- NA_real_

    for (pid in registry[cell_type == ct, program_id]) {
      ans <- score_module(
        z,
        membership[program_id == pid],
        samples,
        coldata$condition_curated == "MASH"
      )
      ans$summary[, `:=`(
        program_id = pid,
        cohort = cohort,
        contrast = "MASH_vs_NORMAL",
        cell_type = ct,
        n_normal = sum(coldata$condition_curated == "NORMAL"),
        n_mash = sum(coldata$condition_curated == "MASH")
      )]
      dynamic_results[[paste(cohort, pid, sep = "::")]] <- ans$summary
      if (!is.null(ans$scores)) {
        ans$scores[, `:=`(
          program_id = pid,
          cohort = cohort,
          condition = coldata$condition_curated
        )]
        dynamic_scores[[paste(cohort, pid, sep = "::")]] <- ans$scores
      }
    }
    coverage_rows[[paste(cohort, ct, sep = "::")]] <- data.table(
      cohort = cohort,
      cell_type = ct,
      n_donors = length(samples),
      n_unique_peaks = ncol(counts),
      n_program_genes_with_promoter_activity = nrow(gene_counts)
    )
  }
}

dynamic <- rbindlist(dynamic_results, fill = TRUE)
dynamic[testable == TRUE, padj := p.adjust(pvalue, method = "BH"), by = cohort]
dynamic[, sensitivity_sign_agree := testable & is.finite(effect) &
  sign(effect) == sign(equal_effect) & sign(effect) == sign(leave_top_effect)]
dynamic[, robust := testable & is.finite(padj) & padj < 0.05 & sensitivity_sign_agree]
dynamic <- merge(
  registry[, .(program_id, display_order, module, program_name)],
  dynamic,
  by = "program_id",
  all.y = TRUE,
  sort = FALSE
)
setorder(dynamic, cohort, display_order)
fwrite(dynamic, file.path(OUT, "dynamic_program_accessibility.tsv"), sep = "\t", quote = FALSE)
if (length(dynamic_scores)) {
  fwrite(
    rbindlist(dynamic_scores),
    file.path(OUT, "dynamic_program_scores.tsv"),
    sep = "\t",
    quote = FALSE
  )
}
fwrite(
  rbindlist(coverage_rows),
  file.path(OUT, "dynamic_coverage.tsv"),
  sep = "\t",
  quote = FALSE
)

# -------------------------------------------------------------------------
# Static fine-mapped-variant/promoter/open-chromatin layer.
# -------------------------------------------------------------------------
message("[ATAC] static 35-study-registry promoter/open-chromatin layer")
tiers <- fread(TIER_FILE)
tiers <- tiers[placement == "main" & tier %in% c(1L, 2L)]
if (nrow(tiers) != 35L || uniqueN(tiers$study_name) != 35L) {
  stop("Expected exactly 35 prespecified Tier-1/2 studies")
}
cs <- fread(CREDIBLE_FILE)
study_coverage <- copy(tiers)
study_coverage[, credible_set_rows_present := study_name %in% unique(cs$study)]
fwrite(
  study_coverage,
  file.path(OUT, "static_study_coverage.tsv"),
  sep = "\t",
  quote = FALSE
)
cs <- merge(cs, tiers, by.x = "study", by.y = "study_name", all = FALSE)
if ("trait.y" %in% names(cs)) cs[, trait := trait.y]
if (uniqueN(cs$study) != 35L) {
  warning(
    "Credible-set artifact represents ", uniqueN(cs$study),
    " of 35 primary studies; missing studies are retained in static_study_coverage.tsv"
  )
}
cs <- cs[
  (!is.na(either_in_cs) & either_in_cs == TRUE) |
    (!is.na(recommended_pip) & recommended_pip > 0.1) |
    (is.na(recommended_pip) & !is.na(max_pip) & max_pip > 0.1)
]
if (SMOKE && nrow(cs) > 5000L) {
  cs <- cs[unique(round(seq(1, .N, length.out = 5000L)))]
}
cs[, row_id := .I]
gr19 <- GRanges(
  seqnames = paste0("chr", cs$chromosome),
  ranges = IRanges(start = as.integer(cs$position), width = 1L)
)
lifted <- liftOver(gr19, import.chain(CHAIN_FILE))
one <- which(lengths(lifted) == 1L)
gr38 <- unlist(lifted[one])
cs38 <- cs[one]
cs38[, `:=`(
  chr_hg38 = as.character(seqnames(gr38)),
  pos_hg38 = start(gr38),
  source_locus = as.character(locus)
)]
variant_gr <- GRanges(
  seqnames = cs38$chr_hg38,
  ranges = IRanges(start = cs38$pos_hg38, width = 1L)
)

static_hits <- list()
static_ct <- if (SMOKE) ct_map[1] else ct_map
for (j in seq_len(nrow(static_ct))) {
  ct <- static_ct$cell_type[j]
  peak_file <- file.path(PEAK_DIR, paste0(static_ct$peak_label[j], "_peaks.bed"))
  if (!file.exists(peak_file)) stop("Missing lineage peak file: ", peak_file)
  bed <- fread(peak_file, header = FALSE)
  if (ncol(bed) == 4L) bed <- bed[, 2:4]
  setnames(bed, c("chr", "start0", "end"))
  bed <- unique(bed)
  peaks <- GRanges(
    seqnames = bed$chr,
    ranges = IRanges(start = as.integer(bed$start0) + 1L, end = as.integer(bed$end))
  )
  accessible_variant <- unique(queryHits(findOverlaps(variant_gr, peaks, ignore.strand = TRUE)))
  if (!length(accessible_variant)) next

  genes_ct <- unique(membership[cell_type == ct, gene_symbol])
  genes_ct <- intersect(genes_ct, names(promoter_gr))
  prom <- promoter_gr[genes_ct]
  hp <- findOverlaps(variant_gr[accessible_variant], prom, ignore.strand = TRUE)
  if (!length(hp)) next
  v_idx <- accessible_variant[queryHits(hp)]
  h <- copy(cs38[v_idx])
  h[, gene_symbol := names(prom)[subjectHits(hp)]]
  h[, cell_type := ct]
  h <- merge(
    h,
    unique(membership[cell_type == ct, .(program_id, gene_symbol)]),
    by = "gene_symbol",
    allow.cartesian = TRUE,
    all = FALSE
  )
  static_hits[[ct]] <- h[, .(
    program_id, cell_type, gene_symbol, study, ancestry, trait,
    trait_scope = fifelse(tier_label == "direct_MASLD", "direct_disease_PDFF", "liver_enzyme"),
    source_locus, chromosome, locus, variant_id, recommended_pip,
    chr_hg38, pos_hg38
  )]
}
hits <- rbindlist(static_hits, fill = TRUE)
if (nrow(hits)) {
  hits <- unique(hits)
  fwrite(hits, file.path(OUT, "static_gwas_atac_hits.tsv"), sep = "\t", quote = FALSE)
  static_global_summary <- hits[, .(
    n_source_rows = .N,
    n_programs = uniqueN(program_id),
    n_program_gene_incidences = uniqueN(paste(program_id, gene_symbol, sep = "::")),
    n_unique_genes = uniqueN(gene_symbol),
    n_program_variant_incidences = uniqueN(paste(program_id, variant_id, sep = "::")),
    n_unique_variants = uniqueN(variant_id),
    n_unique_studies = uniqueN(study)
  ), by = trait_scope]
  fwrite(
    static_global_summary,
    file.path(OUT, "static_gwas_atac_global_summary.tsv"),
    sep = "\t",
    quote = FALSE
  )
  static_summary <- hits[, .(
    n_source_loci = uniqueN(source_locus),
    n_study_loci = uniqueN(paste(study, source_locus, sep = "::")),
    n_variants = uniqueN(variant_id),
    n_genes = uniqueN(gene_symbol),
    n_studies = uniqueN(study),
    n_ancestries = uniqueN(ancestry)
  ), by = .(program_id, trait_scope)]
} else {
  fwrite(data.table(), file.path(OUT, "static_gwas_atac_hits.tsv"), sep = "\t")
  fwrite(data.table(), file.path(OUT, "static_gwas_atac_global_summary.tsv"), sep = "\t")
  static_summary <- data.table()
}
template <- CJ(
  program_id = registry$program_id,
  trait_scope = c("direct_disease_PDFF", "liver_enzyme"),
  unique = TRUE
)
static_summary <- merge(template, static_summary, by = c("program_id", "trait_scope"), all.x = TRUE)
for (v in c("n_source_loci", "n_study_loci", "n_variants", "n_genes", "n_studies", "n_ancestries")) {
  static_summary[is.na(get(v)), (v) := 0L]
}
static_summary <- merge(
  registry[, .(program_id, display_order, cell_type, module, program_name)],
  static_summary,
  by = "program_id",
  all.y = TRUE,
  sort = FALSE
)
scope_coverage <- study_coverage[, .(
  n_registry_studies = .N,
  n_studies_with_credible_sets = sum(credible_set_rows_present)
), by = .(
  trait_scope = fifelse(tier_label == "direct_MASLD", "direct_disease_PDFF", "liver_enzyme")
)]
static_summary <- merge(static_summary, scope_coverage, by = "trait_scope", all.x = TRUE, sort = FALSE)
setorder(static_summary, display_order, trait_scope)
fwrite(static_summary, file.path(OUT, "static_gwas_atac_summary.tsv"), sep = "\t", quote = FALSE)

cat("[ATAC] dynamic testable/robust by cohort:\n")
print(dynamic[, .(testable = sum(testable), robust = sum(robust)), by = cohort])
cat("[ATAC] static promoter/open-peak hits:", nrow(hits), "rows across",
    uniqueN(hits$program_id), "programs\n")
if (SMOKE) cat("[ATAC] smoke test complete; static layer used a 5,000-row credible-set subset\n")
