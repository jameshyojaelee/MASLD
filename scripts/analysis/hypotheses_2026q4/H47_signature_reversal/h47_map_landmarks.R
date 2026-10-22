#!/usr/bin/env Rscript
# H47: one-to-one map of the 978 L1000 landmark genes to refit Ensembl IDs.
# Rule (PRESPEC_H47.md): Entrez -> Ensembl via org.Hs.eg.db, kept if in the
# refit universe; otherwise (or if ambiguous) exact symbol match to GENCODE v49
# gene_name within the refit universe; Ensembl IDs hit by two landmarks dropped.
# Usage: Rscript h47_map_landmarks.R <gene_info.txt.gz> <deg_results.csv> <gencode_meta.tsv.gz> <out_dir>
suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
gene_info_path <- args[1]; deg_path <- args[2]; meta_path <- args[3]; out_dir <- args[4]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gi <- read.delim(gzfile(gene_info_path), colClasses = "character")
lm <- gi[gi$pr_is_lm == "1", c("pr_gene_id", "pr_gene_symbol")]
stopifnot(nrow(lm) == 978, !anyDuplicated(lm$pr_gene_id))

deg <- read.csv(deg_path, colClasses = c(gene = "character"))
deg$ensembl <- sub("\\.[0-9]+(_PAR_Y)?$", "", deg$gene)
stopifnot(!anyDuplicated(deg$ensembl))
universe <- deg$ensembl

meta <- read.delim(gzfile(meta_path), colClasses = "character")
meta <- meta[meta$ensembl_base %in% universe, c("ensembl_base", "gene_name")]

eg2ens <- AnnotationDbi::select(org.Hs.eg.db, keys = lm$pr_gene_id,
                                keytype = "ENTREZID", columns = "ENSEMBL")
eg2ens <- eg2ens[!is.na(eg2ens$ENSEMBL) & eg2ens$ENSEMBL %in% universe, ]

res <- lapply(seq_len(nrow(lm)), function(i) {
  eid <- lm$pr_gene_id[i]; sym <- lm$pr_gene_symbol[i]
  cand <- unique(eg2ens$ENSEMBL[eg2ens$ENTREZID == eid])
  sym_hits <- unique(meta$ensembl_base[meta$gene_name == sym])
  if (length(cand) == 1) {
    return(data.frame(pr_gene_id = eid, pr_gene_symbol = sym, ensembl = cand,
                      map_step = "entrez", n_entrez_candidates = 1))
  }
  pool <- if (length(cand) > 1) intersect(cand, sym_hits) else sym_hits
  if (length(pool) == 1) {
    return(data.frame(pr_gene_id = eid, pr_gene_symbol = sym, ensembl = pool,
                      map_step = if (length(cand) > 1) "entrez_ambiguous_symbol" else "symbol",
                      n_entrez_candidates = length(cand)))
  }
  data.frame(pr_gene_id = eid, pr_gene_symbol = sym, ensembl = NA_character_,
             map_step = "unmapped", n_entrez_candidates = length(cand))
})
res <- do.call(rbind, res)
dup <- res$ensembl[!is.na(res$ensembl) & duplicated(res$ensembl)]
res$map_step[!is.na(res$ensembl) & res$ensembl %in% dup] <- "dropped_shared_ensembl"
res$ensembl[res$map_step == "dropped_shared_ensembl"] <- NA_character_
res$gencode_v49_name <- meta$gene_name[match(res$ensembl, meta$ensembl_base)]
write.table(res, file.path(out_dir, "landmark_map.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "")
cat("landmark map steps:\n"); print(table(res$map_step))
cat("mapped:", sum(!is.na(res$ensembl)), "of", nrow(res), "\n")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_map.txt"))
