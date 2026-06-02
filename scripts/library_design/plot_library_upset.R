#!/usr/bin/env Rscript
# UpSet plot replacing Panel A — aesthetic multi-way Venn alternative
suppressPackageStartupMessages({
  library(data.table); library(stringr); library(ComplexHeatmap); library(grid); library(circlize)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/publication_theme.R"))

lib_dir    <- file.path(BASE, "results/library")
panels_dir <- file.path(lib_dir, "overview_panels")
dir.create(panels_dir, recursive = TRUE, showWarnings = FALSE)
out_pdf    <- file.path(panels_dir, "panel_A_upset.pdf")

core <- fread(file.path(lib_dir, "final_core_degs.csv"))

classify <- function(tok) {
  tok <- trimws(tok)
  if (grepl("^GWAS\\s*\\|", tok))       return("GWAS (human)")
  if (grepl("^Perturb\\s*\\|", tok))    return("Perturb (mouse)")
  if (grepl("^MCD.*\\|\\s*mouse", tok)) return("Mouse RNA-seq")
  if (grepl("^Patient\\s*\\|", tok))    return("Human RNA-seq")
  "Other"
}
dataset_rename <- c("GSE135251" = "Govaere", "GSE130970" = "Hoang",
                    "GSE156918" = "Paquette MCD", "GSE205974" = "Yue MCD",
                    "MCD Week pooled" = "Cas13 MCD mouse")
extract_dataset <- function(tok) {
  tok <- trimws(tok); parts <- str_split(tok, "\\s*\\|\\s*")[[1]]
  if (grepl("^GWAS",    tok)) return("GWAS closest")
  if (grepl("^Perturb", tok)) return("Perturb prior")
  if (grepl("^MCD",     tok)) {
    raw <- sub("\\s*\\(.*", "", parts[3])
    return(if (raw %in% names(dataset_rename)) unname(dataset_rename[raw]) else raw)
  }
  if (grepl("^Patient", tok)) {
    author <- if (parts[2] %in% names(dataset_rename)) unname(dataset_rename[parts[2]]) else parts[2]
    contrast <- if (grepl("NAS", parts[3], ignore.case = TRUE)) "NAS1+"
                else if (grepl("Fibrosis", parts[3], ignore.case = TRUE)) "Fibrosis"
                else sub("\\s*\\(.*", "", parts[3])
    return(paste(author, contrast))
  }
  parts[1]
}

rows <- core[, .(mouse_gene_id, toks = str_split(source_analyses, ";"))]
long <- rows[, .(tok = unlist(toks)), by = mouse_gene_id]
long[, tok := trimws(tok)]
long[, source_type := sapply(tok, classify)]
long[, dataset     := sapply(tok, extract_dataset)]
gd <- unique(long[, .(mouse_gene_id, dataset, source_type)])

# Build gene-membership list per dataset
set_order <- c(
  "Govaere NAS1+", "Govaere Fibrosis", "Hoang NAS1+", "Hoang Fibrosis",
  "Cas13 MCD mouse", "Paquette MCD", "Yue MCD",
  "GWAS closest", "Perturb prior"
)
set_order <- intersect(set_order, unique(gd$dataset))
sets <- lapply(set_order, function(d) unique(gd[dataset == d, mouse_gene_id]))
names(sets) <- set_order

# Dataset -> source_type map (for row coloring)
ds_type <- unique(gd[, .(dataset, source_type)])
ds_type_v <- setNames(ds_type$source_type, ds_type$dataset)[set_order]

src_palette <- c(
  "Human RNA-seq"   = "#C2185B",
  "Mouse RNA-seq"   = "#42A5F5",
  "GWAS (human)"    = "#00695C",
  "Perturb (mouse)" = "#7B1FA2"
)
row_cols <- unname(src_palette[ds_type_v])

m <- make_comb_mat(sets, mode = "distinct")
# Keep singletons (degree=1) + higher-degree intersections >= threshold for clarity
cs <- comb_size(m); degs <- comb_degree(m)
keep <- (degs == 1) | (cs >= 50)
m <- m[, keep]
cs <- comb_size(m); degs <- comb_degree(m)
# Singletons first (sorted by set_order); then multi-set by degree asc, size desc
singleton_idx <- which(degs == 1)
singleton_order <- singleton_idx[order(match(
  sapply(singleton_idx, function(i) set_order[which(as.integer(strsplit(comb_name(m)[i], "")[[1]]) == 1)]),
  set_order))]
multi_idx <- which(degs > 1)
multi_order <- multi_idx[order(degs[multi_idx], -cs[multi_idx])]
m <- m[, c(singleton_order, multi_order)]

set_sizes  <- set_size(m)
comb_sizes <- comb_size(m)

top_ann <- HeatmapAnnotation(
  "Genes" = anno_barplot(comb_sizes, border = FALSE, gp = gpar(fill = "#880E4F", col = NA),
                         height = unit(3.2, "cm"),
                         add_numbers = TRUE,
                         numbers_gp = gpar(fontsize = 6, fontfamily = "Helvetica")),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 7, fontfamily = "Helvetica", fontface = "bold"),
  annotation_name_rot = 90
)
right_ann <- rowAnnotation(
  "Set size" = anno_barplot(set_sizes, border = FALSE,
                            gp = gpar(fill = row_cols, col = NA),
                            width = unit(2.5, "cm"),
                            add_numbers = TRUE,
                            numbers_gp = gpar(fontsize = 6, fontfamily = "Helvetica")),
  annotation_name_gp = gpar(fontsize = 7, fontfamily = "Helvetica", fontface = "bold")
)

ht <- UpSet(m,
            set_order  = set_order,
            comb_order = seq_len(ncol(m)),
            top_annotation   = top_ann,
            right_annotation = right_ann,
            pt_size = unit(3.2, "mm"),
            lwd = 1.4,
            comb_col = "#880E4F",
            bg_col = "#F5F5F5",
            bg_pt_col = "#E0E0E0",
            row_names_gp = gpar(fontsize = 7, fontfamily = "Helvetica", col = row_cols))

src_lgd <- Legend(labels = names(src_palette), title = "Source",
                  legend_gp = gpar(fill = unname(src_palette)),
                  labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica"),
                  title_gp  = gpar(fontsize = 7, fontfamily = "Helvetica", fontface = "bold"),
                  grid_height = unit(3, "mm"), grid_width = unit(3, "mm"))

pdf(out_pdf, width = 9.2, height = 4.8, useDingbats = FALSE)
draw(ht,
     padding = unit(c(3, 3, 3, 3), "mm"),
     annotation_legend_list = list(src_lgd),
     annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()

cat("Wrote:", out_pdf, "\n")
cat("Intersections kept (>=20 genes):", ncol(m), "\n")
