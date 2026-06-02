#!/usr/bin/env Rscript
# Find DEGs in all 6 analyses of fig1d (5 per-study + dream) and check library coverage.
suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

dream <- load_dream_results()
per_study <- load_per_study_de()

upset_studies <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
padj_thr <- 0.05; lfc_thr <- 0.3

deg_lists <- list()
for (ds in upset_studies) {
  ds_data <- per_study[dataset == ds]
  pcol <- intersect(c("padj", "adj.P.Val"), names(ds_data))[1]
  lcol <- intersect(c("logFC", "study_logFC"), names(ds_data))[1]
  degs <- ds_data[get(pcol) < padj_thr & abs(get(lcol)) > lfc_thr, gene]
  deg_lists[[ds]] <- unique(sub("\\..*", "", degs))
}
deg_lists[["Integrated"]] <- unique(sub("\\..*", "", dream[is_dream_deg(dream), gene]))

# Intersection across all 6 sets
all6 <- Reduce(intersect, deg_lists)
cat("6-way intersection size:", length(all6), "\n")

# Load library + map
lib <- fread(file.path(BASE, "results/library/final_core_degs.csv"))
lib_human_ensg <- unique(unlist(strsplit(paste(lib$human_ortholog_ids, collapse = ";"), ";")))
lib_human_sym  <- unique(unlist(strsplit(paste(lib$human_ortholog_symbols, collapse = ";"), ";")))
lib_human_ensg <- lib_human_ensg[lib_human_ensg != "" & !is.na(lib_human_ensg)]
lib_human_sym  <- lib_human_sym[lib_human_sym != "" & !is.na(lib_human_sym)]

# fig1d genes are ENSG-stripped (version removed); check against library human_ortholog_ids
in_lib_ensg <- all6 %in% lib_human_ensg
cat("In library (by ENSG):", sum(in_lib_ensg), "/", length(all6), "\n")

# Also try symbol mapping via dream metadata
sym_col <- intersect(c("gene_name", "symbol", "hgnc_symbol"), names(dream))
if (length(sym_col) > 0) {
  dream_map <- unique(dream[, .(gene = sub("\\..*", "", gene), sym = get(sym_col[1]))])
  all6_sym <- dream_map[gene %in% all6, unique(sym)]
  in_lib_sym <- all6_sym %in% lib_human_sym
  cat("Symbol-resolved hits:", length(all6_sym), "  In library (by symbol):", sum(in_lib_sym), "/", length(all6_sym), "\n")
  out <- data.table(ensg = all6,
                    symbol = dream_map[match(all6, dream_map$gene), sym],
                    in_library = in_lib_ensg)
} else {
  out <- data.table(ensg = all6, in_library = in_lib_ensg)
}

out_path <- file.path(BASE, "results/library/fig1d_6way_intersection_vs_library.csv")
fwrite(out[order(-in_library)], out_path)
cat("Wrote:", out_path, "\n")
cat("\n--- Missing from library ---\n")
miss <- out[in_library == FALSE]
if (nrow(miss) > 0) print(miss) else cat("(none — all present)\n")
