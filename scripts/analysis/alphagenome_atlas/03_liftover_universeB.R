#!/usr/bin/env Rscript
# Step 03: hg19 -> hg38 one-base liftover for Universe B credible-set members.
# Copied rule from Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/08_liftover_variants.R:
# per variant, the number of liftover mappings and the FIRST mapping's chrom/position.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT <- Sys.getenv("AGA_OUT_ROOT", "")
if (!nzchar(OUT)) stop("AGA_OUT_ROOT is unset")
chain_path <- file.path(BASE, "data/broadaway_eqtl/hg19ToHg38.over.chain")
in_path <- file.path(OUT, "tables/variant_posteriors_B_hg19.tsv")
out_path <- file.path(OUT, "tables/universeB_liftover_raw.tsv")
if (file.exists(out_path)) stop("Refusing to overwrite: ", out_path)

posterior <- fread(in_path)
posterior[, hg19_chrom := paste0("chr", chromosome)]
posterior[, variant_key := paste(hg19_chrom, position, effect_allele, other_allele, sep = ":")]
variants <- unique(posterior[, .(variant_key, hg19_chrom, position)], by = "variant_key")
ranges <- GRanges(seqnames = variants$hg19_chrom, ranges = IRanges(start = variants$position, width = 1L))
names(ranges) <- variants$variant_key
lifted <- liftOver(ranges, import.chain(chain_path))
n_liftover_mappings <- as.integer(elementNROWS(lifted))
flat <- unlist(lifted, use.names = FALSE)
first_index <- head(cumsum(c(1L, n_liftover_mappings)), -1L)
mapped <- n_liftover_mappings > 0L
hg38_chrom <- rep(NA_character_, length(lifted))
hg38_position_1based <- rep(NA_integer_, length(lifted))
if (any(mapped)) {
  hg38_chrom[mapped] <- as.character(seqnames(flat))[first_index[mapped]]
  hg38_position_1based[mapped] <- start(flat)[first_index[mapped]]
}
mapping <- data.table(variant_key = variants$variant_key, n_liftover_mappings, hg38_chrom, hg38_position_1based)
fwrite(mapping, out_path, sep = "\t")
message(sprintf("liftover: %d unique variants, %d mapped once, %d unmapped, %d multimapped",
  nrow(mapping), sum(n_liftover_mappings == 1L), sum(n_liftover_mappings == 0L), sum(n_liftover_mappings > 1L)))
