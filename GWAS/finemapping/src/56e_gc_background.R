#!/usr/bin/env Rscript
# 56e_gc_background.R
# Compute genome-wide AND peak-based nucleotide backgrounds for motifbreakR.
# Liver chromatin is GC-heterogeneous (~45-55% GC in open peaks vs ~41% genome-wide).
# Uniform A/C/G/T = 0.25 background mis-prices motif PWM scoring.
#
# Outputs (in results/gwas_atac/):
#   motifbreakr_bg_genome.rds   — genome-wide autosomal A/C/G/T frequencies
#   motifbreakr_bg_peaks.rds    — merged liver-atlas peak A/C/G/T frequencies
#   motifbreakr_bg_summary.csv  — comparison table for paper supplement
#
# Usage: Rscript 56e_gc_background.R
# Env:   motifbreakr

suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(data.table)
  library(dplyr)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE_DIR, "GWAS/finemapping")
ATAC_DIR <- file.path(BASE_DIR, "Analysis/ATAC/Human_Multiome")
PEAK_DIR <- file.path(ATAC_DIR, "results/label_transfer/cell_type_peak_sets_v2")
OUT_DIR  <- file.path(FM_DIR, "results/gwas_atac")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("56e_gc_background.R\n")
cat("Compute nucleotide backgrounds for motifbreakR\n")
cat("============================================================\n\n")

genome <- BSgenome.Hsapiens.UCSC.hg38
autosomes <- paste0("chr", 1:22)

# ── 1. Genome-wide background (autosomes, N-stripped) ───────────────────────
cat("--- 1. Genome-wide background (autosomes) ---\n")

genome_counts <- numeric(4)
names(genome_counts) <- c("A", "C", "G", "T")
total_n_count <- 0
total_length  <- 0

for (chrom in autosomes) {
  seq <- genome[[chrom]]
  freq <- alphabetFrequency(seq, baseOnly = TRUE)
  # baseOnly=TRUE returns A,C,G,T,other
  genome_counts <- genome_counts + as.numeric(freq[c("A","C","G","T")])
  total_n_count <- total_n_count + as.numeric(freq["other"])
  total_length  <- total_length  + as.numeric(length(seq))
  cat(sprintf("  %s: length=%.0f, N=%.0f\n",
              chrom, as.numeric(length(seq)), as.numeric(freq["other"])))
}

bg_genome <- genome_counts / sum(genome_counts)
gc_genome <- bg_genome["C"] + bg_genome["G"]

cat(sprintf("\nGenome-wide (autosomes, N-stripped): total length=%.0f, N=%.0f (%.2f%%)\n",
            total_length, total_n_count, 100 * total_n_count / total_length))
cat(sprintf("  A=%.4f C=%.4f G=%.4f T=%.4f  | GC=%.4f\n",
            bg_genome["A"], bg_genome["C"], bg_genome["G"], bg_genome["T"], gc_genome))

saveRDS(bg_genome, file.path(OUT_DIR, "motifbreakr_bg_genome.rds"))

# ── 2. Peak-based background (9 cell-type BEDs merged) ──────────────────────
cat("\n--- 2. Peak-based background ---\n")

bed_files <- list.files(PEAK_DIR, pattern = "_peaks\\.bed$", full.names = TRUE)
cat("Found", length(bed_files), "cell-type peak BEDs\n")

all_gr <- GRanges()
for (bf in bed_files) {
  ct <- sub("_peaks\\.bed$", "", basename(bf))
  bed <- fread(bf, header = FALSE, col.names = c("chr", "start", "end"))
  # BED is 0-based half-open; GRanges is 1-based closed
  gr <- GRanges(
    seqnames = bed$chr,
    ranges   = IRanges(start = bed$start + 1L, end = bed$end)
  )
  # Restrict to autosomes
  gr <- gr[seqnames(gr) %in% autosomes]
  cat(sprintf("  %-24s n_peaks=%d  total_bp=%d\n",
              ct, length(gr), sum(width(gr))))
  all_gr <- c(all_gr, gr)
}

# Set seqinfo for safe getSeq
seqlevels(all_gr, pruning.mode = "coarse") <- autosomes
seqinfo(all_gr) <- seqinfo(genome)[autosomes]

# Merge overlapping peaks across cell types
merged_gr <- reduce(all_gr)
cat(sprintf("\nMerged peaks across 9 cell types: n=%d total_bp=%d\n",
            length(merged_gr), sum(width(merged_gr))))

# Extract sequences
cat("Extracting peak sequences from BSgenome ...\n")
peak_seqs <- getSeq(genome, merged_gr)

# Alphabet frequencies
cat("Computing alphabet frequencies ...\n")
peak_freq_mat <- alphabetFrequency(peak_seqs, baseOnly = TRUE)
peak_counts <- as.numeric(colSums(peak_freq_mat))
names(peak_counts) <- colnames(peak_freq_mat)
# keep A,C,G,T only
peak_counts_acgt <- peak_counts[c("A","C","G","T")]
n_count_peaks    <- peak_counts["other"]

bg_peaks <- peak_counts_acgt / sum(peak_counts_acgt)
gc_peaks <- bg_peaks["C"] + bg_peaks["G"]

cat(sprintf("\nPeak-based (merged 9-CT): total length=%.0f, N=%.0f (%.2f%%)\n",
            sum(peak_counts), n_count_peaks,
            100 * n_count_peaks / sum(peak_counts)))
cat(sprintf("  A=%.4f C=%.4f G=%.4f T=%.4f  | GC=%.4f\n",
            bg_peaks["A"], bg_peaks["C"], bg_peaks["G"], bg_peaks["T"], gc_peaks))

saveRDS(bg_peaks, file.path(OUT_DIR, "motifbreakr_bg_peaks.rds"))

# ── 3. Summary table ────────────────────────────────────────────────────────
summary_dt <- data.table(
  background = c("uniform", "genome", "peak"),
  A  = c(0.25, bg_genome["A"], bg_peaks["A"]),
  C  = c(0.25, bg_genome["C"], bg_peaks["C"]),
  G  = c(0.25, bg_genome["G"], bg_peaks["G"]),
  T  = c(0.25, bg_genome["T"], bg_peaks["T"]),
  GC = c(0.50, gc_genome,      gc_peaks),
  source = c("flat 0.25 prior",
             "BSgenome.Hsapiens.UCSC.hg38 autosomes (N-stripped)",
             sprintf("Merged scATAC peaks across %d cell types (n=%d, %.1f Mb)",
                     length(bed_files), length(merged_gr),
                     sum(width(merged_gr)) / 1e6))
)

fwrite(summary_dt, file.path(OUT_DIR, "motifbreakr_bg_summary.csv"))

cat("\n=== Background summary table ===\n")
print(summary_dt)

cat("\nWrote:\n")
cat("  ", file.path(OUT_DIR, "motifbreakr_bg_genome.rds"), "\n")
cat("  ", file.path(OUT_DIR, "motifbreakr_bg_peaks.rds"), "\n")
cat("  ", file.path(OUT_DIR, "motifbreakr_bg_summary.csv"), "\n")

cat("\n============================================================\n")
cat("Script 56e_gc_background complete.\n")
cat("============================================================\n")
