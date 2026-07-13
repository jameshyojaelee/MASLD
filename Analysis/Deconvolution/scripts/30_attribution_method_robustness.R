#!/usr/bin/env Rscript
# 30_attribution_method_robustness.R
# ---------------------------------------------------------------------------
# Is the Fig3 "hepatocyte-intrinsic vs composition-driven" split robust to the
# deconvolution method used for the composition adjustment?
#
# Loads the three composition-adjusted attribution tables and compares the gene
# membership of the Hepatocyte_intrinsic (cell-autonomous / controlled-direct-
# effect) and Composition_driven (composition-mediated) classes across methods
# via Jaccard overlap.
#
# Inputs (guarded: warn + skip if not yet produced):
#   Rectangle : RNA-seq/results/causal_inference/rectangle/deconv_attribution_scores.csv   (from 25b; class col = category)   <- produced later
#   MuSiC     : RNA-seq/results/causal_inference/c2_recount/deconv_attribution_scores.csv   (from 25 ; class col = category)
#   BayesPrism: RNA-seq/Human/Patient_Cohorts/analysis/integration/results/deconvolution/bayesprism/bayesprism_attribution_scores.csv (class col = bp_class; "NS" == Not_significant)
#
# Outputs (NEW dir Analysis/Deconvolution/rectangle_comparison/):
#   attribution_class_jaccard.csv     (per set: method-pair Jaccard + intersect/union/sizes)
#   attribution_class_sizes.csv       (class-set size per method)
#   attribution_jaccard_heatmap.pdf   (method x method, one facet per class set)
#
# Figure rules: PDF only (useDingbats=FALSE); base 6 pt Helvetica, no bold; all
# text black; captions via message(); control/reference = #9E9E9E; no lollipops; no 3D.
# Env: rnaseq (ggplot2 / data.table).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})
pdf.options(useDingbats = FALSE)

project_root <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTPUT_DIR <- Sys.getenv("RECT_CMP_OUTDIR",
  unset = file.path(project_root, "Analysis/Deconvolution/rectangle_comparison"))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

CONTROL_GREY <- "#9E9E9E"

# Method -> attribution-scores path. Ordered Rectangle, MuSiC, BayesPrism.
ATTRIB_PATHS <- list(
  Rectangle  = file.path(project_root, "RNA-seq/results/causal_inference/rectangle/deconv_attribution_scores.csv"),
  MuSiC      = file.path(project_root, "RNA-seq/results/causal_inference/c2_recount/deconv_attribution_scores.csv"),
  BayesPrism = file.path(project_root, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/deconvolution/bayesprism/bayesprism_attribution_scores.csv")
)

# Class sets whose robustness we test.
CLASS_SETS <- c("Hepatocyte_intrinsic", "Composition_driven")

# ---------------------------------------------------------------------------
# Loader: return data.table(gene, class) with a normalised class vocabulary.
#   - class column is "category" (MuSiC/Rectangle, from script 25/25b) or
#     "bp_class" (BayesPrism); "NS" is normalised to "Not_significant".
# ---------------------------------------------------------------------------
load_attrib_classes <- function(path, method) {
  if (!file.exists(path)) {
    warning(sprintf("Missing %s attribution scores: %s", method, path))
    return(NULL)
  }
  dt <- data.table::fread(path)
  cls_col <- if ("category" %in% names(dt)) "category" else if ("bp_class" %in% names(dt)) "bp_class" else NA_character_
  if (is.na(cls_col) || !"gene" %in% names(dt)) {
    warning(sprintf("%s: no recognised class column (category/bp_class) or gene col -> skip", method))
    return(NULL)
  }
  cl <- as.character(dt[[cls_col]])
  cl[cl == "NS"] <- "Not_significant"
  # Harmonise gene IDs: MuSiC/Rectangle (script 25/25b) carry a GENCODE version
  # suffix (ENSG...x.NN); BayesPrism is unversioned. Strip the suffix so the
  # three universes are comparable (else every BayesPrism Jaccard is a false 0).
  gene <- sub("\\..*$", "", as.character(dt$gene))
  data.table(gene = gene, class = cl, method = method)
}

jaccard <- function(a, b) {
  if (length(a) == 0L && length(b) == 0L) return(NA_real_)
  u <- length(union(a, b))
  if (u == 0L) return(NA_real_)
  length(intersect(a, b)) / u
}

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
attrib_list <- Filter(Negate(is.null),
                      Map(load_attrib_classes, ATTRIB_PATHS, names(ATTRIB_PATHS)))
methods_present <- names(attrib_list)
if (length(methods_present) < 2L) {
  message("Only ", length(methods_present), " attribution table(s) present (",
          paste(methods_present, collapse = ", "),
          "); need >=2 for a cross-method comparison. ",
          "Rectangle is produced later by 25b_deconv_attribution_rectangle.R.")
  message("Nothing to compare yet -> exiting cleanly.")
  quit(save = "no", status = 0)
}
message("Attribution tables present: ", paste(methods_present, collapse = ", "))

# Gene sets per method per class.
gene_sets <- list()
size_rows <- list()
for (m in methods_present) {
  dt <- attrib_list[[m]]
  for (cs in CLASS_SETS) {
    g <- unique(dt[class == cs, gene])
    gene_sets[[paste(m, cs, sep = "||")]] <- g
    size_rows[[paste(m, cs, sep = "||")]] <- data.table(
      method = m, class_set = cs, n_genes = length(g),
      n_tested = nrow(dt))
  }
}
sizes <- data.table::rbindlist(size_rows, use.names = TRUE, fill = TRUE)
data.table::fwrite(sizes, file.path(OUTPUT_DIR, "attribution_class_sizes.csv"))
message("Wrote attribution_class_sizes.csv.")
print(sizes)

# ---------------------------------------------------------------------------
# Pairwise Jaccard per class set (full method x method matrix incl. diagonal=1)
# ---------------------------------------------------------------------------
jac_rows <- list()
for (cs in CLASS_SETS) {
  for (ma in methods_present) {
    for (mb in methods_present) {
      ga <- gene_sets[[paste(ma, cs, sep = "||")]]
      gb <- gene_sets[[paste(mb, cs, sep = "||")]]
      jac_rows[[paste(cs, ma, mb, sep = "||")]] <- data.table(
        class_set = cs,
        method_a  = ma,
        method_b  = mb,
        jaccard   = jaccard(ga, gb),
        n_intersect = length(intersect(ga, gb)),
        n_union     = length(union(ga, gb)),
        n_a = length(ga),
        n_b = length(gb)
      )
    }
  }
}
jac <- data.table::rbindlist(jac_rows, use.names = TRUE, fill = TRUE)
data.table::fwrite(jac, file.path(OUTPUT_DIR, "attribution_class_jaccard.csv"))
message("Wrote attribution_class_jaccard.csv.")

# Report the off-diagonal (distinct-method) overlaps.
offdiag <- jac[method_a != method_b]
if (nrow(offdiag) > 0) {
  for (i in seq_len(nrow(offdiag))) {
    r <- offdiag[i]
    message(sprintf("  %-22s %s vs %s: Jaccard=%.3f (|A|=%d |B|=%d, shared=%d)",
                    r$class_set, r$method_a, r$method_b,
                    ifelse(is.na(r$jaccard), NA_real_, r$jaccard),
                    r$n_a, r$n_b, r$n_intersect))
  }
}

# ---------------------------------------------------------------------------
# Heatmap: method x method Jaccard, one facet per class set
# ---------------------------------------------------------------------------
jac[, method_a := factor(method_a, levels = methods_present)]
jac[, method_b := factor(method_b, levels = rev(methods_present))]
jac[, class_lab := factor(class_set,
      levels = CLASS_SETS,
      labels = c(Hepatocyte_intrinsic = "Hepatocyte-intrinsic\n(cell-autonomous)",
                 Composition_driven   = "Composition-driven\n(composition-mediated)")[CLASS_SETS])]

p <- ggplot(jac, aes(x = method_a, y = method_b, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(jaccard), "NA", sprintf("%.2f", jaccard))),
            size = 2.1, color = "black") +
  scale_fill_gradient(low = "white", high = "#2E7D32", limits = c(0, 1),
                      na.value = CONTROL_GREY, name = "Jaccard") +
  facet_wrap(~ class_lab) +
  coord_fixed() +
  labs(title = "Attribution-class robustness across deconvolution methods",
       x = NULL, y = NULL) +
  theme_minimal(base_family = "Helvetica", base_size = 6) +
  theme(
    plot.title  = element_text(size = 7, face = "plain", color = "black"),
    strip.text  = element_text(size = 6, color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    legend.title = element_text(size = 6, color = "black"),
    legend.text  = element_text(size = 6, color = "black"),
    panel.grid   = element_blank()
  )
pdf(file.path(OUTPUT_DIR, "attribution_jaccard_heatmap.pdf"), width = 6.0, height = 3.4)
print(p)
invisible(dev.off())
message("Wrote attribution_jaccard_heatmap.pdf.")

message("\n30_attribution_method_robustness.R complete. Outputs in: ", OUTPUT_DIR)
