#!/usr/bin/env Rscript
# KEY MESSAGE: The overlap is close to the chance expectation for these set sizes.
# Retains the reviewed Figure 1B two-bar layout on the native 774-unit interface.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(scales)})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
stopifnot(!file.exists(out))
dir.create(out, recursive = TRUE)
source("scripts/figures/publication_theme.R")
input <- "RNA-seq/results/genetics_state_interface/20261009T142027Z_coloc_r2_corrected_controls_r8/overlap_by_arm.tsv"
audit <- fread(input)
p <- audit[genetic_definition == "susie_pp4_0.5" & bulk_arm == "fragment_23370"]
stopifnot(nrow(p) == 1L, p$n_joint_testable == 14704L,
          p$n_genetic_joint == 426L, p$n_deg_joint == 1155L,
          p$n_overlap_observed == 26L, p$n_genes_bulk == 24115L)
view <- rbind(
  data.table(map = "Colocalized genes", segment = c("convergent", "genetic_only"),
             n = c(26L, 400L)),
  data.table(map = "Canonical DEGs\n(BH q<0.05, |log2FC|>0.50)",
             segment = c("convergent", "disease_state_only"), n = c(26L, 1129L)))
view[, fraction := n / sum(n), by = map]
view[, segment := factor(segment, levels = c("convergent", "genetic_only", "disease_state_only"))]
setorder(view, map, segment)
view[, label := sprintf("%s\n%s", comma(n), percent(fraction, accuracy = .1))]
view[segment == "convergent", label := NA_character_]
view[, map := factor(map, levels = c("Canonical DEGs\n(BH q<0.05, |log2FC|>0.50)", "Colocalized genes"))]
expected <- 426 * 1155 / 14704
segs <- data.table(map = levels(view$map), x = c(expected / 1155, expected / 426))
segs[, map := factor(map, levels = levels(view$map))]
shared <- merge(view[segment == "convergent", .(map, fraction, n)], segs[, .(map, xexp = x)], by = "map")
shared[, xlab := pmax(fraction, xexp)]
shared[, txt := sprintf("%s observed (%s)", comma(n), percent(fraction, accuracy = .1))]
colors <- c(convergent = "#00695C", genetic_only = "#1565C0", disease_state_only = "#C9265E")
plot <- ggplot(view, aes(fraction, map, fill = segment)) +
  geom_col(width = .56, color = "white", linewidth = .35, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = label), position = position_stack(vjust = .5, reverse = TRUE),
            size = 6 / .pt, color = "white", lineheight = .9, na.rm = TRUE) +
  geom_segment(data = segs, inherit.aes = FALSE,
               aes(x = x, xend = x, y = as.numeric(map) - .40, yend = as.numeric(map) + .40),
               linetype = "22", linewidth = .4, color = "black") +
  geom_text(data = segs, inherit.aes = FALSE,
            aes(x = x, y = as.numeric(map) + .46, label = sprintf("expected %s", percent(x, accuracy = .1))),
            size = 6 / .pt, color = "black", hjust = 0, nudge_x = .008) +
  geom_text(data = shared, inherit.aes = FALSE, aes(x = xlab, y = map, label = txt),
            size = 6 / .pt, color = "black", hjust = 0, nudge_x = .014) +
  scale_fill_manual(values = colors, name = NULL,
                    breaks = c("convergent", "disease_state_only", "genetic_only"),
                    labels = c("Both", "Disease-state only", "Colocalized only")) +
  scale_x_continuous(labels = percent, expand = expansion(mult = c(0, .02))) +
  scale_y_discrete(expand = expansion(add = .75)) +
  labs(x = "Fraction of native evidence interface", y = NULL) + theme_masld(base_size = 6) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
        legend.key.width = unit(.38, "cm"))
ggsave(file.path(out, "fig1b_complementary_maps.pdf"), plot, width = 4.8,
       height = 2.25, device = cairo_pdf, bg = "white")
fwrite(p, file.path(out, "fig1b_interface_summary.tsv"), sep = "\t")
fwrite(view, file.path(out, "fig1b_display_source.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(sprintf("Figure 1B: 426 colocalized and 1,155 canonical disease-state genes share 26 in the native 14,704-gene stable-Ensembl interface; %.6f expected. Canonical calls use BH q<0.05 and absolute log2FC>0.5. Seven SuSiE-untestable genes are retained as not positive; they are not eligible SuSiE negatives. The categorical 14,720-gene and Catalog 14,682-gene universes are distinct.", expected), file.path(out, "fig1b_caption.txt"))
writeLines(system2("sha256sum", shQuote(input), stdout = TRUE), file.path(out, "input_sha256.txt"))
cat("FIG1B_CURRENT_INTERFACE_PASS\n")
