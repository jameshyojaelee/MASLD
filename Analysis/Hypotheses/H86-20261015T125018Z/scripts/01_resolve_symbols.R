#!/usr/bin/env Rscript
# H86 step 1 (prespec PRESPEC_H86.md, sha256 7e96f97f...): resolve PXD051911 protein symbols with
# HGNChelper 0.8.15 and org.Hs.eg.db 3.20.0. Reads only the Catalog `symbol` and `ensembl_id`
# columns (never `primary_evidence_class`). Writes:
#   inputs_sha256.tsv          D1 hashes (S1 check; on mismatch writes STEP01_STATUS.txt=STOP_S1 and exits 0
#                              so that step 2 writes RESULT_H86.md with the S1 verdict)
#   versions.tsv               D1 package versions and HGNChelper map date
#   hgnc_all_symbols.tsv       checkGeneSymbols on every DE-table and deposit symbol (D5, D7)
#   unmatched_119.tsv          checkGeneSymbols + alias2SymbolTable on DE symbols with no Catalog symbol (D4)
#   nc2_alias_ensembl.tsv      org.Hs.eg.db ALIAS -> ENSEMBL rows for the remapped old symbols (NC-2)
#   deposit_symbols.tsv        deposit Genes field split on ';' (D7)
# No step is random; the seed is set as the prespec requires.

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k) {
  i <- match(k, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", k)
  args[i + 1]
}
OUT <- getarg("--out")
CAT <- getarg("--catalog")
SEED <- as.integer(getarg("--seed"))
set.seed(SEED)

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
H86 <- normalizePath(file.path(dirname(script_path), ".."))
stopifnot(dir.exists(OUT))
setwd(ROOT)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

PRESPEC_SHA <- "7e96f97f0185657daa34eb561898e44358d1c03a8a2a9e6bb399e9678b2ebbc1"
R4 <- "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z"
expected <- data.table(
  role = c("prespec_copy", "prespec_original", "catalog_gene_index", "catalog_quarantine",
           "catalog_domain_coverage", "protein_de", "deposit_liver", "deposit_meta",
           "projection_script", "gate_cov_py", "writer_remap_R"),
  path = c(file.path(H86, "PRESPEC_H86.md"),
           "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H86.md",
           file.path(CAT, "passport_gene_index.parquet"),
           file.path(CAT, "passport_identity_quarantine.tsv"),
           file.path(CAT, "passport_domain_coverage.tsv"),
           "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv",
           "data/PXD051911/liver_protein_quant.txt",
           "data/PXD051911/meta_data.txt",
           "Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R",
           file.path(R4, "gate/R11-09/cov.py"),
           file.path(R4, "write/H86/remap_feasibility.R")),
  size_expected = c(NA, NA, 468456, 243387, 3868, 1361989, 5599042, 35867, 29975, 1059, 2092),
  sha256_expected = c(PRESPEC_SHA, PRESPEC_SHA,
    "ffd3c4650c77a899e9ebed7327697398c4b6462cfffcd7f0e635cd0d023db06b",
    "f696105dee84a5837f080feadf8dc8f44164b2febaafce3f676883c758e50e6d",
    "d5c541b71700f8ec2365c0ccb9b9875ffbf934cf7f8c3e24fa2416e75e3f194e",
    "4cab14f405049964a441d45a298e6209f5067fb48ec233873b3c839800cb6df3",
    "c4fafdfb2339d54d17a91779098012e08448aa5ebcc8a865988dcc49f92c969e",
    "340874d4f720f909fcfca428fed6da16ba3f7cfdf4d9fd53dd6e0378198ba2bc",
    "20a3c06cfa4e1a726339cedc4b701ab0d9d276585791d1a3c7202045dc19dd6b",
    "ce6af28bf2431a8fe13a388d1337ec566cb8df8f875278847f684b5d6924e0f0",
    "d1332e26541e4c779eae34c881ae79f3da7596dd056e62b01abcd6fff7c4e528")
)
expected[, exists := file.exists(path)]
expected[, size_observed := ifelse(exists, file.size(path), NA_real_)]
expected[, sha256_observed := vapply(seq_len(.N), function(i)
  if (exists[i]) digest(file = path[i], algo = "sha256") else NA_character_, character(1))]
expected[, match := exists & sha256_observed == sha256_expected &
           (is.na(size_expected) | size_observed == size_expected)]
fwrite(expected, file.path(OUT, "inputs_sha256.tsv"), sep = "\t")
cat("[01] D1 input hashes:\n"); print(expected[, .(role, match, sha256_observed)])
if (!all(expected$match)) {
  writeLines("STOP_S1", file.path(OUT, "STEP01_STATUS.txt"))
  cat("[01] S1: hash mismatch for:", expected[match == FALSE]$role, "\n")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(HGNChelper)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(limma)
})

# D1 versions and HGNChelper map date
map_date <- readLines(system.file(file.path("extdata", "date_of_last_update.txt"), package = "HGNChelper"))
versions <- data.table(
  item = c("HGNChelper", "HGNChelper_map_date", "org.Hs.eg.db", "AnnotationDbi", "limma", "data.table", "R"),
  value = c(as.character(packageVersion("HGNChelper")), paste(map_date, collapse = " "),
            as.character(packageVersion("org.Hs.eg.db")), as.character(packageVersion("AnnotationDbi")),
            as.character(packageVersion("limma")), as.character(packageVersion("data.table")),
            R.version.string))
fwrite(versions, file.path(OUT, "versions.tsv"), sep = "\t")
print(versions)

# Catalog symbols and Ensembl IDs only (R library has no parquet reader; pyarrow reads two columns).
py <- normalizePath(file.path(R.home(), "..", "..", "bin", "python"))
cat_tsv <- file.path(OUT, "tmp", "catalog_symbols_noclass.tsv")
code <- sprintf(paste0(
  "import pyarrow.parquet as pq; ",
  "t=pq.read_table('%s', columns=['ensembl_id','symbol']).to_pandas(); ",
  "t.to_csv('%s', sep='\\t', index=False)"),
  file.path(ROOT, CAT, "passport_gene_index.parquet"), cat_tsv)
rc <- system2(py, c("-c", shQuote(code)))
if (rc != 0) stop("catalog symbol export failed")
g <- fread(cat_tsv, na.strings = "")
stopifnot(nrow(g) == 28233L, uniqueN(g$symbol) == 28233L)

p <- fread("Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv", na.strings = "")
dep <- fread("data/PXD051911/liver_protein_quant.txt", select = c("ProteinAccessions", "Genes"), na.strings = "")
dep_long <- dep[!is.na(Genes) & Genes != "", .(symbol = trimws(unlist(strsplit(Genes, ";", fixed = TRUE)))),
                by = .(ProteinAccessions, Genes)]
dep_long <- dep_long[symbol != ""]
dep_long[, multi_gene_group := grepl(";", Genes, fixed = TRUE)]
fwrite(dep_long, file.path(OUT, "deposit_symbols.tsv"), sep = "\t")
cat("[01] deposit groups", nrow(dep), "empty Genes", dep[is.na(Genes) | Genes == "", .N],
    "multi-gene groups", dep[grepl(";", Genes, fixed = TRUE), .N],
    "unique single-gene symbols", dep[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE), uniqueN(Genes)],
    "unique symbols after split", uniqueN(dep_long$symbol), "\n")

# checkGeneSymbols on every DE-table and deposit symbol (elementwise lookup)
all_syms <- sort(unique(c(p$gene, dep_long$symbol)))
hmsg <- character(0)
h_all <- withCallingHandlers(
  suppressWarnings(checkGeneSymbols(all_syms, unmapped.as.na = TRUE)),
  message = function(m) { hmsg <<- c(hmsg, conditionMessage(m)); invokeRestart("muffleMessage") })
h_all <- as.data.table(h_all)
h_all[, n_sugg := ifelse(is.na(Suggested.Symbol), 0L,
                         lengths(strsplit(Suggested.Symbol, " /// ", fixed = TRUE)))]
h_all[, `:=`(in_de = x %in% p$gene, in_deposit = x %in% dep_long$symbol)]
fwrite(h_all, file.path(OUT, "hgnc_all_symbols.tsv"), sep = "\t")
cat("[01] checkGeneSymbols message:", hmsg, "\n")

# D4: the DE symbols with no Catalog symbol, exactly as the writer and audit computed them
um <- setdiff(p$gene, g$symbol)
cat("[01] unmatched DE symbols", length(um), "\n")
h <- as.data.table(suppressMessages(suppressWarnings(checkGeneSymbols(um, unmapped.as.na = TRUE))))
h[, n_sugg := ifelse(is.na(Suggested.Symbol), 0L, lengths(strsplit(Suggested.Symbol, " /// ", fixed = TRUE)))]
a <- data.table(x = um, org_symbol = suppressWarnings(alias2SymbolTable(um, species = "Hs")))
u <- merge(h, a, by = "x", sort = TRUE)
fwrite(u, file.path(OUT, "unmatched_119.tsv"), sep = "\t")

# NC-2 identifier rows: org.Hs.eg.db ALIAS -> ENSEMBL for the old symbols whose single suggestion is a
# Catalog symbol (the remap candidates). Step 2 decides which of these are used.
cand <- h[n_sugg == 1L & Suggested.Symbol %in% g$symbol]$x
cat("[01] remap candidates for NC-2 query", length(cand), "\n")
ens <- as.data.table(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = cand,
                                                            keytype = "ALIAS", columns = "ENSEMBL")))
fwrite(ens, file.path(OUT, "nc2_alias_ensembl.tsv"), sep = "\t")
cat("[01] NC-2 old symbols with >=1 Ensembl ID", ens[!is.na(ENSEMBL), uniqueN(ALIAS)], "/", length(cand), "\n")

writeLines("OK", file.path(OUT, "STEP01_STATUS.txt"))
cat("[01] done\n")
