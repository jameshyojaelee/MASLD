#!/usr/bin/env Rscript
# 61b_combine_metabolite_coloc.R
# ---------------------------------------------------------------------------
# Combine per-shard metabolite/lipid COLOC outputs (60b) into a gene-level
# atlas-column table for the multi-evidence atlas.
#
# COMPLEMENTARY (supporting) layer — NOT a primary genetic-causal tier.
#
# Locus -> gene assignment (many-to-many; metabolites are enzyme/transporter
# products, so a locus can map to several genes and a gene to several
# metabolites — we do NOT collapse to a single 1:1 per-gene PP4):
#   (1) PRIMARY: the Broadaway eGene(s) that themselves colocalize at the same
#       MASLD locus (existing eGene COLOC, susie_coloc/<gwas>/). A shared causal
#       variant between metabolite-QTL and eQTL implicates that gene.
#   (2) FALLBACK: nearest Broadaway eGene by hg19 cis-window midpoint to the
#       metabolite-QTL top SNP (when no eGene colocalizes at the locus).
#
# Known-MASLD-biomarker flag: BCAAs / ceramides / bile acids / acylcarnitines /
# PC-PE (regex on metabolite trait name).
#
# Output: RNA-seq/results/multi_evidence/mqtl_atlas_columns.tsv  (key human_symbol)
#   mqtl_n_colocalizing_metabolites, mqtl_best_pp4, mqtl_metabolite_list,
#   mqtl_known_masld_biomarker  (+ supporting detail cols)
# Also writes a long per-(metabolite x locus x gene) evidence table.
# ---------------------------------------------------------------------------
suppressMessages({library(data.table)})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR  <- file.path(BASE_DIR, "GWAS/finemapping")
PP4_TH  <- as.numeric(Sys.getenv("MQTL_PP4_THRESH", "0.5"))
EGENE_LOCUS_PAD <- 5e5   # eGene COLOC top_snp within +/-500kb of metabolite locus

OUT_TSV <- file.path(BASE_DIR, "RNA-seq/results/multi_evidence/mqtl_atlas_columns.tsv")
LONG_TSV <- file.path(BASE_DIR, "RNA-seq/results/multi_evidence/mqtl_coloc_long.tsv")
dir.create(dirname(OUT_TSV), recursive = TRUE, showWarnings = FALSE)

# --- gather shard outputs from both arms ---
shard_files <- c(
  list.files(file.path(BASE_DIR, "data/external/chen2023_mqtl/coloc_out"),
             pattern = "\\.tsv$", full.names = TRUE),
  list.files(file.path(BASE_DIR, "data/external/lipidqtl/coloc_out"),
             pattern = "\\.tsv$", full.names = TRUE)
)
shard_files <- shard_files[file.exists(shard_files)]
cat("Shard files:", length(shard_files), "\n")
if (length(shard_files) == 0) stop("No shard outputs found — run 60b first.")

res <- rbindlist(lapply(shard_files, function(f)
  tryCatch(fread(f), error = function(e) NULL)), fill = TRUE)
res <- res[!is.na(PP.H4)]
cat("Total locus-metabolite tests (PP.H4>0.1 kept by 60b):", nrow(res), "\n")

# --- colocalizing metabolite-locus pairs ---
hits <- res[PP.H4 > PP4_TH]
cat("Colocalizing metabolite-locus pairs (PP.H4 >", PP4_TH, "):", nrow(hits), "\n")
if (nrow(hits) == 0) {
  cat("No colocalizing pairs above threshold; writing empty atlas columns.\n")
  fwrite(data.table(human_symbol=character(), mqtl_n_colocalizing_metabolites=integer(),
                    mqtl_best_pp4=double(), mqtl_metabolite_list=character(),
                    mqtl_known_masld_biomarker=logical()), OUT_TSV, sep="\t")
  quit(save = "no", status = 0)
}

# --- known-MASLD-biomarker flag ---
biomarker_pat <- paste(c(
  "leucine","isoleucine","\\bvaline\\b","ceramide",
  "bile acid","cholate","cholic","deoxycholat","chenodeoxychol","glycochol",
  "taurochol","glycodeoxychol","ursodeoxychol","lithochol",
  "carnitine","acylcarnitine","acyl carnitine",
  "phosphatidylcholine","phosphatidylethanolamine","lysophosphatidyl",
  "glycerophosphocholine","glycerophosphoethanolamine"
), collapse = "|")
hits[, known_biomarker := grepl(biomarker_pat, metabolite, ignore.case = TRUE)]

# --- locus -> gene assignment ---
# (1) eGene COLOC at the same locus: load per-GWAS susie_coloc top_snp positions.
egene <- list()
for (gw in unique(hits$gwas)) {
  d <- file.path(FM_DIR, "results/susie_coloc", gw)
  if (!dir.exists(d)) next
  fs <- list.files(d, pattern = "^susie_coloc_chr.*\\.csv$", full.names = TRUE)
  if (length(fs) == 0) next
  e <- rbindlist(lapply(fs, function(f) tryCatch(fread(f), error=function(z) NULL)), fill = TRUE)
  if (nrow(e) == 0) next
  # parse top_snp "chr:pos"
  e <- e[!is.na(top_snp) & grepl(":", top_snp)]
  e[, e_pos := as.integer(sub("^[0-9]+:", "", top_snp))]
  e[, gwas := gw]
  egene[[gw]] <- e[, .(gene, ensembl, chr, gwas, e_pos, egene_pp4 = PP.H4.abf)]
}
egene <- rbindlist(egene, fill = TRUE)
cat("eGene COLOC top-SNP records loaded:", nrow(egene), "\n")

# (2) Broadaway hg19 gene cis-windows for nearest-gene fallback
gw_win <- fread(file.path(BASE_DIR, "data/external/chen2023_mqtl/broadaway_gene_hg19_window.tsv"))

assign_genes <- function(hrow) {
  # primary: eGenes that colocalize at the same locus (same gwas, chr, top_snp near metabolite locus)
  cands <- egene[gwas == hrow$gwas & chr == hrow$chr &
                 abs(e_pos - hrow$cs_lead_pos) <= EGENE_LOCUS_PAD &
                 egene_pp4 > 0.3]
  if (nrow(cands) > 0) {
    return(data.table(gene = cands$gene, ensembl = cands$ensembl,
                      assign_method = "egene_coloc", egene_pp4 = cands$egene_pp4))
  }
  # fallback: nearest Broadaway gene cis-window midpoint to metabolite top SNP (hg19)
  ts_pos <- suppressWarnings(as.integer(sub("^[0-9]+:", "", hrow$top_snp)))
  if (is.na(ts_pos)) ts_pos <- hrow$cs_lead_pos
  near <- gw_win[chr == hrow$chr]
  if (nrow(near) == 0) return(NULL)
  near[, dist := pmin(abs(cis_mid - ts_pos), abs(cis_start - ts_pos), abs(cis_end - ts_pos))]
  best <- near[order(dist)][1]
  data.table(gene = best$gene, ensembl = best$ENSG,
             assign_method = "nearest_egene", egene_pp4 = NA_real_)
}

long_list <- list()
for (i in seq_len(nrow(hits))) {
  hrow <- hits[i]
  ga <- assign_genes(hrow)
  if (is.null(ga) || nrow(ga) == 0) next
  long_list[[i]] <- cbind(
    ga,
    data.table(human_symbol = ga$gene, metabolite = hrow$metabolite,
               source = hrow$source, accession = hrow$accession,
               gwas = hrow$gwas, locus_id = hrow$locus_id,
               chr = hrow$chr, cs_lead_pos = hrow$cs_lead_pos,
               mqtl_pp4 = hrow$PP.H4, known_biomarker = hrow$known_biomarker))
}
long <- rbindlist(long_list, fill = TRUE)
long <- long[!is.na(human_symbol) & human_symbol != ""]
fwrite(long, LONG_TSV, sep = "\t")
cat("Wrote long per-(metabolite x locus x gene) table:", nrow(long), "rows ->", LONG_TSV, "\n")

# --- collapse to gene-level atlas columns (many-to-many preserved as a list) ---
atlas <- long[, .(
  mqtl_n_colocalizing_metabolites = uniqueN(metabolite),
  mqtl_best_pp4 = max(mqtl_pp4, na.rm = TRUE),
  mqtl_metabolite_list = paste(sort(unique(metabolite)), collapse = "; "),
  mqtl_known_masld_biomarker = any(known_biomarker, na.rm = TRUE),
  mqtl_n_biomarker_metabolites = uniqueN(metabolite[known_biomarker]),
  mqtl_best_gwas = gwas[which.max(mqtl_pp4)],
  mqtl_assign_method = paste(sort(unique(assign_method)), collapse = ";"),
  mqtl_sources = paste(sort(unique(source)), collapse = ";")
), by = human_symbol]
setorder(atlas, -mqtl_best_pp4)
fwrite(atlas, OUT_TSV, sep = "\t")

cat("\n============================================================\n")
cat("Genes with metabolite/lipid-QTL COLOC (PP.H4>", PP4_TH, "):", nrow(atlas), "\n")
cat("  via eGene-COLOC:", long[assign_method=="egene_coloc", uniqueN(human_symbol)],
    "| via nearest-eGene:", long[assign_method=="nearest_egene", uniqueN(human_symbol)], "\n")
cat("  genes with a KNOWN-MASLD-biomarker hit:", atlas[mqtl_known_masld_biomarker==TRUE, .N], "\n")
cat("  top genes:\n")
print(head(atlas[, .(human_symbol, mqtl_n_colocalizing_metabolites, mqtl_best_pp4,
                     mqtl_known_masld_biomarker, mqtl_best_gwas)], 20))
cat("Wrote", OUT_TSV, "\n")
cat("============================================================\n")
