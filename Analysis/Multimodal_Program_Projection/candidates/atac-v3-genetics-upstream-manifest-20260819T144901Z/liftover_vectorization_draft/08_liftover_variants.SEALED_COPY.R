#!/usr/bin/env Rscript

# Unique hg19-to-hg38 one-base mapping for replayed SuSiE posterior variants.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
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

files <- Sys.glob(file.path(
  CANDIDATE, "genetics", "replay_execution", "batch_*", "exports",
  "*", "chr*", "*", "variant_posteriors.tsv.gz"
))
if (length(files) == 0L) stop("No promoted-COLOC replay posterior files")
batch_manifests <- Sys.glob(file.path(
  CANDIDATE, "genetics", "replay_execution", "batch_*", "batch_manifest.tsv"
))
expected_manifests <- file.path(
  CANDIDATE, "genetics", "replay_execution", paste0("batch_", 1:5),
  "batch_manifest.tsv"
)
if (!setequal(normalizePath(batch_manifests), normalizePath(expected_manifests))) {
  stop("Exactly five named replay batch manifests are required")
}
manifested_posteriors <- character()
manifest_input_rows <- list()
for (manifest_path in batch_manifests) {
  manifest <- fread(manifest_path)
  batch_root <- dirname(manifest_path)
  expected_batch <- as.integer(sub("batch_", "", basename(batch_root)))
  required_manifest <- c(
    "release_id", "batch_id", "gwas_name", "ensembl", "artifact", "bytes", "sha256"
  )
  if (!all(required_manifest %in% names(manifest)) ||
      any(manifest$release_id != RELEASE_ID) ||
      any(manifest$batch_id != expected_batch)) {
    stop("Replay batch manifest schema or identity is invalid: ", manifest_path)
  }
  for (index in seq_len(nrow(manifest))) {
    artifact <- file.path(batch_root, manifest$artifact[[index]])
    artifact_norm <- normalizePath(artifact, mustWork = FALSE)
    batch_norm <- normalizePath(batch_root)
    if (!startsWith(artifact_norm, paste0(batch_norm, .Platform$file.sep))) {
      stop("Replay manifest artifact escapes its batch root: ", artifact)
    }
    if (!file.exists(artifact) ||
        file.info(artifact)$size != manifest$bytes[[index]] ||
        digest(artifact, algo = "sha256", file = TRUE, serialize = FALSE) !=
          manifest$sha256[[index]]) {
      stop("Validated replay artifact changed: ", artifact)
    }
    if (basename(artifact) == "variant_posteriors.tsv.gz") {
      manifested_posteriors <- c(manifested_posteriors, artifact_norm)
    }
  }
  manifest_input_rows[[length(manifest_input_rows) + 1L]] <- data.table(
    release_id = RELEASE_ID,
    role = paste0("replay_batch_manifest_", expected_batch),
    path = manifest_path,
    sha256 = digest(manifest_path, algo = "sha256", file = TRUE, serialize = FALSE)
  )
}
if (!setequal(normalizePath(files), manifested_posteriors)) {
  stop("Replay posterior files do not exactly match validated batch manifests")
}
frozen <- fread(file.path(CANDIDATE, "input_manifest.tsv"))
chain_row <- frozen[role == "hg19_to_hg38_chain"]
if (nrow(chain_row) != 1L) stop("Frozen hg19-to-hg38 chain is not unique")
chain_path <- chain_row$relative_or_absolute_path[[1L]]
if (!grepl("^/", chain_path)) chain_path <- file.path(BASE, chain_path)
if (!file.exists(chain_path) ||
    digest(chain_path, algo = "sha256", file = TRUE, serialize = FALSE) !=
      chain_row$sha256[[1L]]) {
  stop("Frozen hg19-to-hg38 chain hash mismatch")
}
PENDING <- file.path(
  CANDIDATE, "genetics", paste0(".liftover.pending.", Sys.getpid())
)
stale_pending <- Sys.glob(file.path(CANDIDATE, "genetics", ".liftover.pending.*"))
if (length(stale_pending) > 0L) stop("Unresolved prior liftover stage: ", stale_pending[[1L]])
if (dir.exists(PENDING) || file.exists(PENDING)) stop("Pending liftover path exists")
dir.create(PENDING, recursive = TRUE, showWarnings = FALSE)
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
chain <- import.chain(chain_path)
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
fwrite(posterior, file.path(PENDING, "variant_liftover_raw.tsv.gz"), sep = "\t")
liftover_inputs <- rbindlist(c(list(data.table(
  release_id = RELEASE_ID,
  role = "hg19_to_hg38_chain",
  path = chain_path,
  sha256 = chain_row$sha256[[1L]]
)), manifest_input_rows), use.names = TRUE)
fwrite(liftover_inputs, file.path(PENDING, "liftover_input_manifest.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(PENDING, "sessionInfo.txt"))
if (!file.rename(PENDING, OUT)) stop("Atomic liftover stage rename failed")
message("[LIFTOVER] Wrote ", nrow(posterior), " posterior rows")
