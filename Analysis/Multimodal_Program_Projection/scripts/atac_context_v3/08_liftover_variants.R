#!/usr/bin/env Rscript

# Unique hg19-to-hg38 one-base mapping for replayed SuSiE posterior variants.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

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
OUT <- file.path(CANDIDATE, "genetics", "liftover")
if (dir.exists(OUT)) stop("Refusing to overwrite liftover output: ", OUT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

files <- Sys.glob(file.path(
  CANDIDATE, "genetics", "replay_exports", "*", "chr*", "*", "variant_posteriors.tsv.gz"
))
if (length(files) == 0L) stop("No promoted-COLOC replay posterior files")
tables <- lapply(files, function(path) {
  value <- fread(path)
  value[, replay_source := substring(path, nchar(CANDIDATE) + 2L)]
  value
})
posterior <- rbindlist(tables, fill = TRUE)
required <- c(
  "gwas_name", "gene", "ensembl", "chr", "signal_pair_index", "snp",
  "SNP.PP.H4", "hg19_position", "allele1", "allele2"
)
if (!all(required %in% names(posterior))) stop("Replay posterior schema incomplete")
posterior[, hg19_chrom := paste0("chr", chr)]
posterior[, variant_key := paste(hg19_chrom, hg19_position, allele1, allele2, sep = ":")]
variants <- unique(
  posterior[, .(variant_key, hg19_chrom, hg19_position, allele1, allele2)],
  by = "variant_key"
)
ranges <- GRanges(
  seqnames = variants$hg19_chrom,
  ranges = IRanges(start = variants$hg19_position, width = 1L)
)
names(ranges) <- variants$variant_key
chain <- import.chain(file.path(BASE, "data/broadaway_eqtl/hg19ToHg38.over.chain"))
lifted <- liftOver(ranges, chain)
mapping <- rbindlist(lapply(seq_along(lifted), function(index) {
  value <- lifted[[index]]
  if (length(value) == 0L) {
    return(data.table(
      variant_key = variants$variant_key[index], n_liftover_mappings = 0L,
      hg38_chrom = NA_character_, hg38_position_1based = NA_integer_
    ))
  }
  data.table(
    variant_key = variants$variant_key[index],
    n_liftover_mappings = length(value),
    hg38_chrom = as.character(seqnames(value))[1L],
    hg38_position_1based = start(value)[1L]
  )
}))
posterior <- merge(posterior, mapping, by = "variant_key", all.x = TRUE, sort = FALSE)
posterior[, unique_liftover := n_liftover_mappings == 1L]
fwrite(posterior, file.path(OUT, "variant_liftover_raw.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
message("[LIFTOVER] Wrote ", nrow(posterior), " posterior rows")
