#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(stringr)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/publication_theme.R"))

lib_dir <- file.path(BASE, "results/library")
out_pdf <- file.path(lib_dir, "library_overview.pdf")

core <- fread(file.path(lib_dir, "final_core_degs.csv"))

# ---- Palette: magenta (human/disease) + cyan (mouse) + teal (GWAS) + violet (Perturb) ----
src_palette <- c(
  "Human RNA-seq"   = "#C2185B",
  "Mouse RNA-seq"   = "#42A5F5",
  "GWAS (human)"    = "#00695C",
  "Perturb (mouse)" = "#7B1FA2"
)
cross_species_palette <- c(
  "Both species"  = "#880E4F",
  "Human only"    = "#C2185B",
  "Mouse only"    = "#42A5F5",
  "Non-RNA-seq"   = "#00695C"
)
share_palette <- c(
  "Unique to dataset" = "#880E4F",
  "Shared"            = "#42A5F5"
)
biotype_palette <- c(
  protein_coding = "#C2185B", lncRNA = "#E91E63", miRNA = "#880E4F",
  IG_V_gene = "#42A5F5", IG_C_gene = "#1565C0", TR_C_gene = "#0D47A1",
  snoRNA = "#00695C", misc_RNA = "#7B1FA2"
)

# ---- Source classification + dataset name remap ----
classify <- function(tok) {
  tok <- trimws(tok)
  if (grepl("^GWAS\\s*\\|", tok))       return("GWAS (human)")
  if (grepl("^Perturb\\s*\\|", tok))    return("Perturb (mouse)")
  if (grepl("^MCD.*\\|\\s*mouse", tok)) return("Mouse RNA-seq")
  if (grepl("^Patient\\s*\\|", tok))    return("Human RNA-seq")
  "Other"
}
dataset_rename <- c(
  "GSE135251"       = "Govaere",
  "GSE130970"       = "Hoang",
  "GSE156918"       = "Paquette MCD",
  "GSE205974"       = "Yue MCD",
  "MCD Week pooled" = "Cas13 MCD mouse"
)
extract_dataset <- function(tok) {
  tok <- trimws(tok)
  parts <- str_split(tok, "\\s*\\|\\s*")[[1]]
  if (grepl("^GWAS",    tok)) return("GWAS closest genes")
  if (grepl("^Perturb", tok)) return("Perturb-Multimodal")
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

# ---- Parse source_analyses to gene x (dataset, source_type) ----
rows <- core[, .(mouse_gene_id, mouse_biotype,
                 toks = str_split(source_analyses, ";"))]
long <- rows[, .(tok = unlist(toks)), by = .(mouse_gene_id, mouse_biotype)]
long[, tok := trimws(tok)]
long[, source_type := sapply(tok, classify)]
long[, dataset     := sapply(tok, extract_dataset)]
gd <- unique(long[, .(mouse_gene_id, dataset, source_type)])

# How many *distinct* datasets back each gene (shared vs unique across whole library)
n_ds_per_gene <- gd[, .(n_ds = uniqueN(dataset)), by = mouse_gene_id]
gd <- merge(gd, n_ds_per_gene, by = "mouse_gene_id")
gd[, shared_class := ifelse(n_ds == 1, "Unique to dataset", "Shared")]
gd[, shared_class := factor(shared_class, levels = c("Unique to dataset", "Shared"))]

# ================================================================
# Panel A — per-dataset stacked bars, faceted by source type
# ================================================================
pA_dt <- gd[, .N, by = .(source_type, dataset, shared_class)]
dataset_order <- gd[, .N, by = dataset][order(N), dataset]
pA_dt[, dataset := factor(dataset, levels = dataset_order)]
pA_dt[, source_type := factor(source_type,
                              levels = c("Human RNA-seq","Mouse RNA-seq","GWAS (human)","Perturb (mouse)"))]
totals <- pA_dt[, .(N = sum(N)), by = .(source_type, dataset)]

pA <- ggplot(pA_dt, aes(N, dataset, fill = shared_class)) +
  geom_col(width = 0.75) +
  geom_text(data = totals, aes(x = N, y = dataset, label = N),
            inherit.aes = FALSE, hjust = -0.15, size = 2.3, family = "Helvetica") +
  facet_grid(source_type ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = as_labeller(c(
               "Human RNA-seq" = "Human RNA-seq",
               "Mouse RNA-seq" = "Mouse RNA-seq",
               "GWAS (human)"    = "",
               "Perturb (mouse)" = ""))) +
  scale_fill_manual(values = share_palette, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.17))) +
  labs(title = sprintf("Final library composition (n = %s DEGs)",
                       format(nrow(core), big.mark = ",")),
       x = "Genes", y = NULL) +
  theme_masld() +
  theme(legend.position = "top",
        strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1),
        strip.placement = "outside",
        panel.spacing.y = unit(0.25, "lines"))

# ================================================================
# Panel B — convergence histogram stacked by cross-species category
# ================================================================
gene_src <- dcast(unique(long[, .(mouse_gene_id, source_type)]),
                  mouse_gene_id ~ source_type, value.var = "source_type",
                  fun.aggregate = function(x) length(x) > 0, fill = FALSE)
has_h <- gene_src[["Human RNA-seq"]]; has_m <- gene_src[["Mouse RNA-seq"]]
gene_src[, cross := ifelse(has_h & has_m, "Both species",
                    ifelse(has_h, "Human only",
                    ifelse(has_m, "Mouse only", "Non-RNA-seq")))]
core_x <- merge(core[, .(mouse_gene_id, n_analyses)],
                gene_src[, .(mouse_gene_id, cross)], by = "mouse_gene_id")
core_x[, n_bin := ifelse(n_analyses >= 7, "7+", as.character(n_analyses))]
core_x[, n_bin := factor(n_bin, levels = c("1","2","3","4","5","6","7+"))]
core_x[, cross := factor(cross, levels = names(cross_species_palette))]

pB_totals <- core_x[, .N, by = n_bin]
pB <- ggplot(core_x, aes(n_bin, fill = cross)) +
  geom_bar(width = 0.7) +
  geom_text(data = pB_totals, aes(x = n_bin, y = N, label = N),
           inherit.aes = FALSE, vjust = -0.4, size = 2.3, family = "Helvetica") +
  scale_fill_manual(values = cross_species_palette, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "B. Evidence convergence × cross-species support",
       x = "Supporting analyses (n)", y = "Genes") +
  theme_masld() + theme(legend.position = "top")

# ================================================================
# Panel C — biotype composition
# ================================================================
pC_dt <- core[, .N, by = mouse_biotype][order(-N)]
pC_dt[, mouse_biotype := factor(mouse_biotype, levels = mouse_biotype)]
pC <- ggplot(pC_dt, aes(mouse_biotype, N, fill = mouse_biotype)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), vjust = -0.4, size = 2.2, family = "Helvetica") +
  scale_y_log10(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = biotype_palette, guide = "none") +
  labs(title = "C. Gene biotypes (log10)", x = NULL, y = "Genes") +
  theme_masld() + theme(axis.text.x = element_text(angle = 25, hjust = 1))

# ================================================================
# Save each panel individually so the user can scale them
# ================================================================
panels_dir <- file.path(lib_dir, "overview_panels")
dir.create(panels_dir, recursive = TRUE, showWarnings = FALSE)

save_fig(pA, file.path(panels_dir, "panel_A_datasets_unique_shared.pdf"),
         width = fig_full_width * 0.75, height = 3.0)
save_fig(pB, file.path(panels_dir, "panel_B_convergence_cross_species.pdf"),
         width = fig_half_width, height = 2.8)
save_fig(pC, file.path(panels_dir, "panel_C_biotypes.pdf"),
         width = fig_half_width, height = 2.5)

cat("Wrote individual panels to: ", panels_dir, "\n", sep = "")
