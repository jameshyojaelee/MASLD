#!/usr/bin/env Rscript
# 56_motif_disruption_v2.R
# Variant of 56_motif_disruption.R with selectable nucleotide background:
#   MB_BACKGROUND=uniform  (legacy 0.25 each)
#   MB_BACKGROUND=genome   (autosomal A/C/G/T from BSgenome)
#   MB_BACKGROUND=peak     (merged scATAC peak A/C/G/T)
#
# All other logic identical to Script 56. Output filenames suffixed with
# `_bg{uniform|genome|peak}` so the 3 runs land side-by-side.
#
# Prereq: src/56e_gc_background.R has been run (produces *_bg_{genome,peaks}.rds)
# Usage:  MB_BACKGROUND=peak Rscript 56_motif_disruption_v2.R
# Env:    motifbreakr (+ module load meme/5.5.9 for FIMO)

library(data.table)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifbreakR)
library(TFBSTools)
library(MotifDb)

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE_DIR, "GWAS/finemapping")
ATAC_DIR <- file.path(BASE_DIR, "Analysis/ATAC/Human_Multiome")
OUT_DIR  <- file.path(FM_DIR, "results/gwas_atac")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Resolve background ───────────────────────────────────────────────────────
BG_MODE <- tolower(Sys.getenv("MB_BACKGROUND", unset = "uniform"))
stopifnot(BG_MODE %in% c("uniform", "genome", "peak"))

if (BG_MODE == "uniform") {
  bg_vec <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
} else if (BG_MODE == "genome") {
  bg_file <- file.path(OUT_DIR, "motifbreakr_bg_genome.rds")
  if (!file.exists(bg_file)) stop("Missing ", bg_file, " — run 56e_gc_background.R first.")
  bg_vec <- readRDS(bg_file)
} else {
  bg_file <- file.path(OUT_DIR, "motifbreakr_bg_peaks.rds")
  if (!file.exists(bg_file)) stop("Missing ", bg_file, " — run 56e_gc_background.R first.")
  bg_vec <- readRDS(bg_file)
}
# Ensure A/C/G/T order, finite, sums to 1
bg_vec <- bg_vec[c("A","C","G","T")]
stopifnot(all(is.finite(bg_vec)), abs(sum(bg_vec) - 1) < 1e-6)

cat("============================================================\n")
cat("56_motif_disruption_v2.R  (background =", BG_MODE, ")\n")
cat("============================================================\n")
cat(sprintf("Background: A=%.4f C=%.4f G=%.4f T=%.4f  (GC=%.4f)\n\n",
            bg_vec["A"], bg_vec["C"], bg_vec["G"], bg_vec["T"],
            bg_vec["C"] + bg_vec["G"]))

# Filename suffix for this run
SUFFIX <- paste0("_bg", BG_MODE)

# ── 1. Load peak-overlapping variants ────────────────────────────────────────
ann_file <- file.path(OUT_DIR, "gwas_atac_variant_annotation.csv")
if (!file.exists(ann_file)) {
  stop("Variant annotation not found. Run 55_gwas_atac_variant_overlap.R first.")
}
ann <- fread(ann_file)
cat("Loaded", nrow(ann), "variant-peak annotations\n")

variants <- ann %>%
  distinct(chromosome, position, allele1, allele2, chr_hg38, pos_hg38, max_pip) %>%
  arrange(desc(max_pip))
cat("Unique variants in peaks:", nrow(variants), "\n")

# ── 2. Load disease regulons for cross-referencing ───────────────────────────
regulon_file <- file.path(ATAC_DIR, "scenic_plus/disease_regulons.csv")
if (file.exists(regulon_file)) {
  regulons <- fread(regulon_file)
  disease_tfs <- toupper(regulons$tf_name)
  cat("Loaded", length(disease_tfs), "SCENIC+ disease regulon TFs\n")
} else {
  cat("WARNING: disease_regulons.csv not found — skipping regulon cross-reference\n")
  disease_tfs <- character(0)
  regulons <- data.table()
}

# ── 3. Prepare variant input for motifbreakR ─────────────────────────────────
cat("\n--- Preparing variants for motifbreakR ---\n")
genome <- BSgenome.Hsapiens.UCSC.hg38

gr_snps <- GRanges(
  seqnames = variants$chr_hg38,
  ranges   = IRanges(start = variants$pos_hg38, width = 1),
  strand   = "*"
)
ref_alleles <- as.character(getSeq(genome, gr_snps))
cat("Retrieved reference alleles for", length(ref_alleles), "variants\n")

variants$ref_genome <- ref_alleles
variants$REF <- ifelse(variants$allele1 == ref_alleles, variants$allele1,
                 ifelse(variants$allele2 == ref_alleles, variants$allele2, NA_character_))
variants$ALT <- ifelse(variants$allele1 == ref_alleles, variants$allele2,
                 ifelse(variants$allele2 == ref_alleles, variants$allele1, NA_character_))

n_mismatch <- sum(is.na(variants$REF))
if (n_mismatch > 0) {
  cat("WARNING:", n_mismatch, "variants where neither allele matches hg38 reference — skipping\n")
  variants <- variants %>% filter(!is.na(REF))
}

variants <- variants %>% filter(nchar(REF) == 1 & nchar(ALT) == 1)
cat("SNV variants for motif analysis:", nrow(variants), "\n")

if (nrow(variants) == 0) {
  cat("WARNING: No SNVs remain. Writing empty results.\n")
  fwrite(data.table(),
         file.path(OUT_DIR, paste0("motif_disruption_scores", SUFFIX, ".csv")))
  quit(save = "no", status = 0)
}

gr_mb <- GRanges(
  seqnames = variants$chr_hg38,
  ranges   = IRanges(start = variants$pos_hg38, width = 1),
  strand   = "*"
)
mcols(gr_mb)$SNP_id <- paste(variants$chromosome, variants$position,
                              variants$allele1, variants$allele2, sep = ":")
mcols(gr_mb)$REF <- DNAStringSet(variants$REF)
mcols(gr_mb)$ALT <- DNAStringSet(variants$ALT)
names(gr_mb) <- mcols(gr_mb)$SNP_id
attr(gr_mb, "genome.package") <- "BSgenome.Hsapiens.UCSC.hg38"
cat("Prepared", length(gr_mb), "SNVs for motifbreakR\n")

# ── 4. Run motifbreakR ──────────────────────────────────────────────────────
cat("\n--- Running motifbreakR (bg =", BG_MODE, ") ---\n")

jaspar_motifs <- query(MotifDb, andStrings = c("hsapiens"),
                       orStrings = c("jaspar2022", "jaspar2024"))
if (length(jaspar_motifs) == 0) {
  jaspar_motifs <- query(MotifDb, "jaspar")
  jaspar_motifs <- jaspar_motifs[grep("Hsapiens", names(jaspar_motifs))]
}
if (length(jaspar_motifs) == 0) {
  jaspar_motifs <- query(MotifDb, "HOCOMOCOv11-core-A")
  if (length(jaspar_motifs) == 0) jaspar_motifs <- query(MotifDb, "HOCOMOCO")
}
cat("Using", length(jaspar_motifs), "TF motifs from MotifDb\n")

tryCatch({
  mb_results <- motifbreakR(
    snpList   = gr_mb,
    filterp   = TRUE,
    pwmList   = jaspar_motifs,
    threshold = 1e-4,
    method    = "ic",
    bkg       = bg_vec,
    BPPARAM   = BiocParallel::MulticoreParam(
      workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
  )
  cat("motifbreakR returned", length(mb_results), "variant-motif results\n")
}, error = function(e) {
  cat("ERROR in motifbreakR:", conditionMessage(e), "\nAttempting serial ...\n")
  mb_results <<- motifbreakR(
    snpList   = gr_mb,
    filterp   = TRUE,
    pwmList   = jaspar_motifs,
    threshold = 1e-4,
    method    = "ic",
    bkg       = bg_vec,
    BPPARAM   = BiocParallel::SerialParam()
  )
  cat("motifbreakR (serial) returned", length(mb_results), "variant-motif results\n")
})

if (length(mb_results) == 0) {
  cat("WARNING: motifbreakR returned no results. Writing empty output.\n")
  fwrite(data.table(),
         file.path(OUT_DIR, paste0("motif_disruption_scores", SUFFIX, ".csv")))
  quit(save = "no", status = 0)
}

mb_dt <- as.data.table(mb_results)
cat("Raw motifbreakR hits:", nrow(mb_dt), "\n")

mb_dt <- mb_dt %>% filter(effect == "strong" | abs(alleleDiff) > 0.3)
cat("Strong effect / |alleleDiff| > 0.3:", nrow(mb_dt), "variant-motif pairs\n")

if ("geneSymbol" %in% colnames(mb_dt)) {
  mb_dt$tf_name <- mb_dt$geneSymbol
} else if ("providerName" %in% colnames(mb_dt)) {
  mb_dt$tf_name <- mb_dt$providerName
} else {
  mb_dt$tf_name <- mb_dt$providerId
}

mb_dt$motif_in_disease_regulon <- toupper(mb_dt$tf_name) %in% disease_tfs

if (nrow(regulons) > 0) {
  regulon_lookup <- regulons %>%
    mutate(tf_upper = toupper(tf_name)) %>%
    select(tf_upper, regulon_activity_diff, activity_padj) %>%
    distinct(tf_upper, .keep_all = TRUE)

  mb_dt <- mb_dt %>%
    mutate(tf_upper = toupper(tf_name)) %>%
    left_join(regulon_lookup, by = "tf_upper") %>%
    select(-tf_upper)
}

n_regulon <- sum(mb_dt$motif_in_disease_regulon)
cat("  Disrupting SCENIC+ disease regulon motifs:", n_regulon, "\n")

# ── 5. (Skip FIMO validation in the 3-background sweep) ─────────────────────
# FIMO validation is identical regardless of motifbreakR background, so we
# don't burn 3× the wall time recomputing it. The canonical Script 56 run
# already produced fimo_comparison.csv. We keep an `fimo_concordant = NA`
# column for schema compatibility with downstream consumers.
mb_dt$fimo_concordant <- NA
cat("Skipping FIMO validation (background-independent; see Script 56 baseline)\n")

# ── 6. Write final output ───────────────────────────────────────────────────
cat("\n=== Writing outputs ===\n")

out_cols <- c("SNP_id", "seqnames", "start", "end",
              "tf_name", "effect", "alleleDiff", "Refpvalue", "Altpvalue",
              "motif_in_disease_regulon", "fimo_concordant")
if ("regulon_activity_diff" %in% colnames(mb_dt)) {
  out_cols <- c(out_cols, "regulon_activity_diff", "activity_padj")
}
out_cols <- intersect(out_cols, colnames(mb_dt))
motif_out <- mb_dt %>% select(all_of(out_cols))

pip_lookup <- variants %>%
  mutate(SNP_id = paste(chromosome, position, allele1, allele2, sep = ":")) %>%
  # FIX (review F100): carry the BSgenome-resolved REF/ALT (re-resolved at v2:97-101) so 56c can
  # align the Broadaway/Currin beta to the genome ALT (the motifbreakR alleleDiff axis), not the
  # unordered GWAS allele2. SNP_id is allele1:allele2 order, which is NOT REF/ALT-ordered.
  select(SNP_id, max_pip, ref_genome = REF, alt_genome = ALT)
motif_out <- motif_out %>% left_join(pip_lookup, by = "SNP_id")

motif_out <- motif_out %>%
  mutate(priority_score = max_pip * abs(alleleDiff)) %>%
  arrange(desc(priority_score))

motif_out$bg_mode <- BG_MODE

out_file <- file.path(OUT_DIR, paste0("motif_disruption_scores", SUFFIX, ".csv"))
fwrite(motif_out, out_file)
cat("Motif disruption scores:", nrow(motif_out), "variant-motif pairs ->\n  ", out_file, "\n")

# FIX (review Script-56 drift): when run as the canonical producer (chain via run_56.sbatch), also
# write the unsuffixed motif_disruption_scores.csv that 56c/57/56i consume — so v2 (BSgenome-resolved
# REF/ALT, the F100 dependency) is the single source of truth instead of the original 56_motif_disruption.R.
if (tolower(Sys.getenv("MB_WRITE_CANONICAL", unset = "false")) %in% c("1", "true", "yes")) {
  canon_file <- file.path(OUT_DIR, "motif_disruption_scores.csv")
  fwrite(motif_out, canon_file)
  cat("  Also wrote canonical (review Script-56 drift):", canon_file, "\n")
}

cat("\n=== Summary (bg =", BG_MODE, ") ===\n")
cat("Total variant-motif disruptions:", nrow(motif_out), "\n")
cat("Unique variants with disruptions:", n_distinct(motif_out$SNP_id), "\n")
cat("Unique TF motifs disrupted:", n_distinct(motif_out$tf_name), "\n")
cat("  Strong effect:", sum(motif_out$effect == "strong", na.rm = TRUE), "\n")
cat("  In disease regulon:", sum(motif_out$motif_in_disease_regulon, na.rm = TRUE), "\n")

top_tfs <- motif_out %>% count(tf_name, sort = TRUE) %>% head(10)
cat("\nTop 10 disrupted TF motifs:\n"); print(as.data.frame(top_tfs))

if (any(motif_out$motif_in_disease_regulon)) {
  reg_hits <- motif_out %>%
    filter(motif_in_disease_regulon) %>%
    select(SNP_id, tf_name, alleleDiff, max_pip) %>%
    arrange(desc(max_pip))
  cat("\nDisease regulon disruptions (top 20):\n")
  print(as.data.frame(head(reg_hits, 20)))
}

cat("\n============================================================\n")
cat("Script 56_v2 (bg =", BG_MODE, ") complete.\n")
cat("============================================================\n")
