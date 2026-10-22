# H66 step 0: run-to-donor table for the single-cell lineage field.
# Sources Analysis/SingleCell/scripts/lib_donor_collapse.R (read-only, not edited) after
# verifying its sha256 and the sha256 of every file it reads against PRESPEC_H66.md.
# Draws no random numbers.
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
out_i <- which(args == "--out")
if (length(out_i) != 1) stop("usage: 00_donor_map.R --out <run dir>")
OUT <- args[out_i + 1]
if (!dir.exists(OUT)) stop("output directory does not exist: ", OUT)

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
H66 <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
PRESPEC <- file.path(H66, "PRESPEC_H66.md")
PRESPEC_SHA <- "66be821e9195362a44c837835ed2b2e2c730da96e83e5ce6c3040b5594fd6b6c"

sha256 <- function(p) {
  r <- system2("sha256sum", shQuote(p), stdout = TRUE)
  strsplit(r, " +")[[1]][1]
}

log_lines <- character(0)
say <- function(...) {
  s <- paste0(...)
  message(s)
  log_lines <<- c(log_lines, s)
}

# Prespec hash
got <- sha256(PRESPEC)
if (got != PRESPEC_SHA) {
  writeLines(c(sprintf("S1 STOP: PRESPEC_H66.md sha256 %s != %s", got, PRESPEC_SHA)), file.path(OUT, "00_donor_map_STOP.txt"))
  stop("S1: prespec sha256 mismatch")
}
say("prespec sha256 OK ", got)
prespec_txt <- paste(readLines(PRESPEC, warn = FALSE), collapse = "\n")

# Inputs read by this step (sha256 copied from the prespec input table; each is also
# checked to appear verbatim in the prespec text)
inp <- data.table(
  path = c(
    "Analysis/SingleCell/scripts/lib_donor_collapse.R",
    "data/GSE244832/metadata/donor_pairing.csv",
    "data/GSE185477/metadata/donor_pairing.csv",
    "data/GSE202379/metadata/donor_pairing.csv",
    "data/GSE136103/metadata/donor_pairing.csv",
    "Liver_Atlas/metadata/SraRunTable.csv",
    "Liver_Atlas/metadata/GSE192740_sampleInfo_scRNAseq.tsv",
    "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/lineage_cell_counts.tsv",
    "Analysis/SingleCell/metadata/snrna_donor_sex.csv"),
  sha = c(
    "e85072df76335e4d15580b0aeac7c669465c3cb9c39294f560a5cb505fe81d0f",
    "822764cb0740dadb55617ba850e592d9884e254edd56cca846f5083f75e7c426",
    "240b8427db9b5c080e8875f3d9a4f1a91f9daf13f092b595e8b14a329e71027b",
    "3ce1f50fe863fab210ed097c46c99ee73abf4acded86f2586c02efcf402918fb",
    "1c17ef26c7d64cab253ca77b1c45a33f6f0718af723ed9c8f548181acdb93d8d",
    "327ade4c42b4759b3d17544f391484736958984ca641f4a008e10d2dba159452",
    "3f36dbc7414bd7a408fb163ef12466e8192502498c544f47a1bff9fe9f06e35f",
    "f1e98eceb95c050455e73bffe835c05f21fe2df001df8cf8842ee757733dec8b",
    "1e0c2a8862eeb8fbc69da173d907fcac22384860265f7cc0b99fb0eca32f0f1d"))
inp[, in_prespec := vapply(sha, function(s) grepl(s, prespec_txt, fixed = TRUE), logical(1))]
if (!all(inp$in_prespec)) stop("internal: an expected sha256 is not in the prespec text")
inp[, got := vapply(file.path(BASE, path), sha256, character(1))]
inp[, ok := got == sha]
fwrite(inp, file.path(OUT, "00_donor_map_input_sha256.tsv"), sep = "\t")
if (!all(inp$ok)) {
  bad <- inp[ok == FALSE, path]
  writeLines(sprintf("S1 STOP: sha256 mismatch for %s", bad), file.path(OUT, "00_donor_map_STOP.txt"))
  stop("S1: input sha256 mismatch: ", paste(bad, collapse = ", "))
}
say("00_donor_map: ", nrow(inp), " input sha256 values match the prespec")

source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
m <- build_srr_to_donor_map(BASE)
say("srr_to_donor map entries: ", length(m), "; donors in map: ", uniqueN(unname(m)))

cc <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/lineage_cell_counts.tsv"))
meta <- fread(file.path(BASE, "Analysis/SingleCell/metadata/snrna_donor_sex.csv"))[, .(sample, dataset)]
runs <- data.table(run = cc$sample)
runs[, in_map := run %in% names(m)]
runs[, donor := ifelse(in_map, unname(m[run]), run)]
runs <- merge(runs, meta, by.x = "run", by.y = "sample", all.x = TRUE, sort = FALSE)
say("runs ", nrow(runs), "; donors ", uniqueN(runs$donor), "; runs without dataset ", sum(is.na(runs$dataset)))
# a donor must belong to a single dataset
dd <- runs[, .(n_ds = uniqueN(dataset)), by = donor]
if (any(dd$n_ds != 1)) stop("a donor spans more than one dataset")
fwrite(runs[order(dataset, donor, run)], file.path(OUT, "run_to_donor.tsv"), sep = "\t")
fwrite(data.table(run = names(m), donor = unname(m)), file.path(OUT, "srr_to_donor_map_full.tsv"), sep = "\t")
writeLines(log_lines, file.path(OUT, "00_donor_map.log"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_00_donor_map.txt"))
