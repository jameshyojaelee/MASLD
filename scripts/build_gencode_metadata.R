#!/usr/bin/env Rscript
# build_gencode_metadata.R — Parse GENCODE v49 GTF into a gene-level metadata
# TSV (ensembl_id, gene_name, chromosome, gene_biotype) for confounder-gene
# classification in NMF + other downstream scripts.
#
# Writes: data/gencode_v49_gene_metadata.tsv.gz (~60K rows)

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
gtf_path <- Sys.getenv("GENCODE_GTF",
  "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz")
out_path <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

cat("Reading GTF (gene entries only):", gtf_path, "\n")
# Read gzipped GTF, keep only 'gene' feature rows (avoid exons/transcripts)
con <- gzfile(gtf_path, "rt")
on.exit(close(con), add = TRUE)
lines <- character(0)
chunk <- 1e6
repeat {
  x <- readLines(con, n = chunk)
  if (length(x) == 0) break
  x <- x[!startsWith(x, "#")]
  # Only lines where 3rd tab-field == "gene"
  cols3 <- vapply(strsplit(x, "\t", fixed = TRUE), `[`, character(1), 3)
  lines <- c(lines, x[!is.na(cols3) & cols3 == "gene"])
}
cat(sprintf("  %d gene-level entries parsed\n", length(lines)))

split_tab <- strsplit(lines, "\t", fixed = TRUE)
chrom <- vapply(split_tab, `[`, character(1), 1)
attr_str <- vapply(split_tab, `[`, character(1), 9)

# Extract fields from GTF attributes
extract_attr <- function(field) {
  pat <- sprintf('%s "([^"]+)"', field)
  m <- regmatches(attr_str, regexpr(pat, attr_str))
  out <- rep(NA_character_, length(attr_str))
  hit <- nzchar(m)
  out[hit] <- sub(pat, "\\1", m[hit])
  out
}

gene_id    <- extract_attr("gene_id")
gene_name  <- extract_attr("gene_name")
gene_type  <- extract_attr("gene_type")    # GENCODE uses gene_type (Ensembl: gene_biotype)

dt <- data.table(
  gene_id    = gene_id,
  gene_name  = gene_name,
  chromosome = chrom,
  gene_biotype = gene_type
)
# Strip Ensembl version suffix on gene_id for matching
dt[, ensembl_base := sub("\\.[0-9]+(_PAR_Y)?$", "", gene_id)]
# Deduplicate (PAR_Y entries have same ensembl_base; keep first = autosomal)
dt <- unique(dt, by = c("ensembl_base", "chromosome"))

cat(sprintf("  Unique gene rows: %d\n", nrow(dt)))
cat("  Chromosome distribution:\n")
print(dt[, .N, by = chromosome][order(-N)][1:25])
cat("  Biotype distribution (top 10):\n")
print(dt[, .N, by = gene_biotype][order(-N)][1:10])

fwrite(dt, out_path, sep = "\t", compress = "gzip")
cat(sprintf("Wrote %s (%d rows)\n", out_path, nrow(dt)))
