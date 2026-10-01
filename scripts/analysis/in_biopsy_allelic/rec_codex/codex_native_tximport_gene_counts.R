#!/usr/bin/env Rscript
# Native Salmon expected-count aggregation only; no missing-row imputation.

fail <- function(message) stop(message, call. = FALSE)
require_ok <- function(ok, message) if (!isTRUE(ok)) fail(message)
argv <- commandArgs(trailingOnly = TRUE)
required <- c("--quant", "--tx2gene", "--sample", "--out", "--report")
require_ok(length(argv) == 2L * length(required), "Expected --quant --tx2gene --sample --out --report pairs")
keys <- argv[seq.int(1L, length(argv), by = 2L)]
require_ok(!anyDuplicated(keys) && setequal(keys, required), "Unknown, missing or duplicate CLI option")
opt <- setNames(argv[seq.int(2L, length(argv), by = 2L)], keys)
require_ok(nzchar(Sys.getenv("SLURM_JOB_ID")), "Native quantification import must run under SLURM")
require_ok(grepl("^[A-Za-z0-9_.-]+$", opt[["--sample"]]), "Invalid sample identifier")
require_ok(requireNamespace("tximport", quietly = TRUE) &&
             as.character(utils::packageVersion("tximport")) == "1.34.0", "Requires installed tximport1.34.0")
require_ok(requireNamespace("jsonlite", quietly = TRUE) &&
             as.character(utils::packageVersion("jsonlite")) == "2.0.0", "Requires installed jsonlite2.0.0")

sha256 <- function(path) {
  result <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  require_ok(is.null(attr(result, "status")) && length(result) == 1L, "SHA256 command failed")
  sub(" .*", "", result)
}
mapping_sha <- "0ac3dd2c4cc9b6923afea3041f5d4e94b67b9ded0f33acc4c028e6f3845eb8d2"
require_ok(file.exists(opt[["--quant"]]) && file.exists(opt[["--tx2gene"]]), "Missing quantification or mapping input")
require_ok(sha256(opt[["--tx2gene"]]) == mapping_sha, "Complete frozen native transcript mapping changed")
outputs <- c(counts = opt[["--out"]], report = opt[["--report"]],
             native_zero_transcripts = paste0(opt[["--out"]], ".native_zero_transcripts.tsv"))
require_ok(!anyDuplicated(normalizePath(outputs, mustWork = FALSE)), "Output paths must differ")
require_ok(!any(file.exists(outputs)), "Refusing output overwrite")
require_ok(all(dir.exists(dirname(outputs))), "Output parent directories must already exist")

# Select only identity columns, avoiding the large exon-geometry metadata.
header <- strsplit(readLines(opt[["--tx2gene"]], n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
require_ok(!anyDuplicated(header) && all(c("transcript_id", "raw_gene_id", "length") %in% header), "Invalid native mapping columns")
classes <- rep("NULL", length(header))
classes[header %in% c("transcript_id", "raw_gene_id")] <- "character"
classes[header == "length"] <- "numeric"
mapping <- utils::read.delim(opt[["--tx2gene"]], colClasses = classes,
                            check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
require_ok(nrow(mapping) == 227368L && !anyDuplicated(mapping$transcript_id), "Requires all227368 unique native versioned transcripts")
require_ok(all(!is.na(mapping$transcript_id) & nzchar(mapping$transcript_id) &
                 !is.na(mapping$raw_gene_id) & nzchar(mapping$raw_gene_id)), "Empty native transcript/gene identity")
require_ok(all(grepl("\\.[0-9]+$", mapping$transcript_id)), "Versioned transcript identities must be preserved")
require_ok(is.numeric(mapping$length) && all(is.finite(mapping$length) & mapping$length > 0 &
                                             mapping$length == floor(mapping$length)),
           "Native mapping lengths must be finite positive integers")
tx2gene <- data.frame(TXNAME = mapping$transcript_id, GENEID = mapping$raw_gene_id,
                      stringsAsFactors = FALSE)

quant <- utils::read.delim(opt[["--quant"]], check.names = FALSE, stringsAsFactors = FALSE)
quant_columns <- c("Name", "Length", "EffectiveLength", "TPM", "NumReads")
require_ok(identical(names(quant), quant_columns), "Unexpected quant.sf schema")
require_ok(nrow(quant) == 227368L && !anyDuplicated(quant$Name) &&
             !anyNA(quant$Name) && setequal(quant$Name, tx2gene$TXNAME),
           "Native quantifier must return the complete exact227368 versioned transcript axis")
for (column in quant_columns[-1L]) {
  require_ok(is.numeric(quant[[column]]) && all(is.finite(quant[[column]]) & quant[[column]] >= 0),
             paste("Nonfinite, negative or nonnumeric quant.sf", column))
}
require_ok(all(quant$Length > 0 & quant$Length == floor(quant$Length)),
           "Native quantifier lengths must be positive integers")
# The admitted index uses no clipping; quantifier lengths must retain the
# frozen reconstructed transcript lengths, matched by exact versioned ID.
require_ok(all(quant$Length == mapping$length[match(quant$Name, mapping$transcript_id)]),
           "Native quantifier Length differs from the frozen no-clip transcript length")
require_ok(sum(quant$NumReads) > 0, "Quantification has zero aggregate estimated count")
quant_sha <- sha256(opt[["--quant"]])

tx <- tximport::tximport(setNames(opt[["--quant"]], opt[["--sample"]]), type = "salmon",
                        txOut = TRUE, countsFromAbundance = "no", ignoreTxVersion = FALSE,
                        ignoreAfterBar = FALSE, dropInfReps = TRUE)
require_ok(nrow(tx$counts) == 227368L && setequal(rownames(tx$counts), quant$Name),
           "tximport altered or omitted native transcript identities")
require_ok(identical(tx$countsFromAbundance, "no"), "tximport changed expected-count units")
native <- quant$NumReads[match(rownames(tx$counts), quant$Name)]
require_ok(all(abs(tx$counts[, 1L] - native) <= 1e-10 + 1e-12 * abs(native)),
           "Imported transcript estimates differ from native NumReads")
gene <- tximport::summarizeToGene(tx, tx2gene = tx2gene, countsFromAbundance = "no",
                                 ignoreTxVersion = FALSE, ignoreAfterBar = FALSE)
require_ok(setequal(rownames(gene$counts), unique(tx2gene$GENEID)) &&
             !anyDuplicated(rownames(gene$counts)), "Gene aggregation omitted, padded or duplicated native gene rows")
require_ok(identical(gene$countsFromAbundance, "no") &&
             all(is.finite(gene$counts) & gene$counts >= 0), "Invalid aggregated expected counts")
expected <- rowsum(matrix(quant$NumReads, ncol = 1L),
                   tx2gene$GENEID[match(quant$Name, tx2gene$TXNAME)])
expected <- expected[match(rownames(gene$counts), rownames(expected)), 1L]
require_ok(all(abs(gene$counts[, 1L] - expected) <= 1e-8 + 1e-12 * abs(expected)),
           "Gene counts differ from direct native transcript sums")
require_ok(abs(sum(gene$counts) - sum(quant$NumReads)) <= 1e-8 + 1e-12 * sum(quant$NumReads),
           "Gene aggregation does not conserve native estimated-count sum")
require_ok(sha256(opt[["--quant"]]) == quant_sha && sha256(opt[["--tx2gene"]]) == mapping_sha,
           "Input changed during import")

zero <- quant[quant$NumReads == 0, , drop = FALSE]
zero$raw_gene_id <- tx2gene$GENEID[match(zero$Name, tx2gene$TXNAME)]
zero$length_le_k31 <- zero$Length <= 31
zero$effective_length_zero <- zero$EffectiveLength == 0
zero$structural_unobservability_candidate <- zero$length_le_k31 | zero$effective_length_zero
# These flags are algorithmic diagnostics, not proof of biological absence;
# unclipped Length alone does not reconstruct Salmon reference preprocessing.
result <- data.frame(gene_id = rownames(gene$counts), count = gene$counts[, 1L],
                     check.names = FALSE, stringsAsFactors = FALSE)
names(result)[2L] <- opt[["--sample"]]
require_ok(!any(file.exists(outputs)), "Output appeared during import; refusing overwrite")
utils::write.table(result, outputs[["counts"]], sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(zero, outputs[["native_zero_transcripts"]], sep = "\t", quote = FALSE, row.names = FALSE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
require_ok(length(script_arg) == 1L, "Cannot record script path")
script <- sub("^--file=", "", script_arg)
report <- list(sample = opt[["--sample"]], native_transcripts = nrow(quant), native_genes = nrow(result),
               transcript_count_sum = sum(quant$NumReads), gene_count_sum = sum(gene$counts),
               native_zero_transcripts = nrow(zero), native_zero_genes = sum(gene$counts == 0),
               structural_candidate_zero_transcripts = sum(zero$structural_unobservability_candidate),
               countsFromAbundance = "no", ignoreTxVersion = FALSE, ignoreAfterBar = FALSE,
               inferential_replicates_imported = FALSE, rows_padded_or_removed = FALSE,
               full_model_gene_extraction_delegated = TRUE, private_PISCES_equivalence_asserted = FALSE,
               biological_absence_asserted = FALSE,
               structural_flag_limit = "Zero estimates retained natively; length<=31 or zero effective length is an algorithmic candidate, not biological absence or a complete post-preprocessing eligibility test",
               inputs = list(quant = normalizePath(opt[["--quant"]]), quant_sha256 = quant_sha,
                             tx2gene = normalizePath(opt[["--tx2gene"]]), tx2gene_sha256 = mapping_sha),
               script_sha256 = sha256(script), slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
               output_sha256 = list(gene_counts = sha256(outputs[["counts"]]),
                                    native_zero_transcripts = sha256(outputs[["native_zero_transcripts"]])),
               versions = list(R = R.version.string, tximport = as.character(utils::packageVersion("tximport")),
                               jsonlite = as.character(utils::packageVersion("jsonlite"))),
               sessionInfo = capture.output(sessionInfo()))
jsonlite::write_json(report, outputs[["report"]], auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(sprintf("sample=%s native_transcripts=%d native_genes=%d count_sum=%.12g\n",
            opt[["--sample"]], nrow(quant), nrow(result), sum(gene$counts)))
