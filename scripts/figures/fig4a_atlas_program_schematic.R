#!/usr/bin/env Rscript
# KEY MESSAGE: donor-aware liver data define fixed programs for cross-sectional comparisons.
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
project_root <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Descriptive atlas census from the current paper ledger, not inferential n.
ledger <- readLines(file.path(project_root, "docs/RESULTS.md"), warn = FALSE)
required <- c(
  "| Prespecified Hotspot programs | 117 |",
  "| Integrated libraries / biological donors / cells | 275 / 104 / 1,232,318 |",
  "| Analyzed libraries / biological donors / cells | 273 / 102 / 1,232,285 |"
)
stopifnot(all(vapply(required, function(value) any(grepl(value, ledger, fixed = TRUE)), logical(1))))
facts <- data.frame(
  metric = c("single_cell_datasets", "integrated_donors", "integrated_cells",
             "analyzed_donors", "prespecified_programs", "tested_lineages"),
  value = c(7, 104, 1232318, 102, 117, 5),
  unit = c("dataset", "biological donor", "cell", "biological donor", "program", "lineage"),
  state = c("descriptive atlas coverage", "descriptive atlas census", "descriptive atlas census",
            "analysis census; inferential n varies by contrast", "fixed definitions and weights",
            "hepatocyte;macrophage;fibroblast;cholangiocyte;T cell"),
  source = "docs/RESULTS.md#single-cell-programs-and-composition"
)
get_value <- function(metric_name) facts$value[match(metric_name, facts$metric)]
first <- sprintf("%d datasets  ·  %d donors  ·  %.2fM cells",
                 get_value("single_cell_datasets"), get_value("integrated_donors"),
                 get_value("integrated_cells") / 1e6)
second_left <- sprintf("%d analyzed donors", get_value("analyzed_donors"))
second_right <- sprintf("%d programs / %d lineages", get_value("prespecified_programs"),
                        get_value("tested_lineages"))
third <- "Disease  ·  lineage  ·  cross-sectional transport"

grDevices::pdf.options(useDingbats = FALSE)
stem <- "fig4a_atlas_program_schematic"
cairo_pdf(file.path(out, paste0(stem, ".pdf")), width = 180 / 72,
          height = 29 / 72, family = "Helvetica", onefile = FALSE)
grid.newpage()
pushViewport(viewport(xscale = c(0, 180), yscale = c(0, 29), default.units = "points"))
gp <- gpar(fontfamily = "Helvetica", fontsize = 6, fontface = "plain", col = "#4D4D4D")
text_width <- function(value) convertWidth(grobWidth(textGrob(value, gp = gp)), "points", valueOnly = TRUE)
stopifnot(text_width(first) < 174, text_width(third) < 174)
left_x <- 2
arrow_start <- left_x + text_width(second_left) + 5
arrow_end <- arrow_start + 12
right_x <- arrow_end + 5
stopifnot(right_x + text_width(second_right) < 178)
label <- function(value, x, y, color = "#4D4D4D") {
  grid.text(value, x = unit(x, "native"), y = unit(y, "native"), just = "left",
            gp = gpar(fontfamily = "Helvetica", fontsize = 6, fontface = "plain", col = color))
}
label(first, 2, 24)
label(second_left, left_x, 14.5)
grid.segments(unit(arrow_start, "native"), unit(14.5, "native"),
              unit(arrow_end, "native"), unit(14.5, "native"),
              arrow = arrow(length = unit(2.5, "points"), type = "closed"),
              gp = gpar(col = "#9E9E9E", lwd = 0.5))
label(second_right, right_x, 14.5, "#C9265E")
label(third, 2, 5)
popViewport()
dev.off()
write.table(facts, file.path(out, paste0(stem, "_source.tsv")), sep = "\t",
            row.names = FALSE, quote = FALSE)
writeLines(paste(
  "Descriptive overview of seven human liver single-cell datasets: 104 biological donors and",
  "1,232,318 cells in the integrated atlas, with 102 donors and 1,232,285 cells analyzed.",
  "The 117 prespecified programs across five tested lineages retain fixed definitions and weights.",
  "Program disease association, lineage localization and cross-sectional bulk transport are separate",
  "comparisons with their own eligible biological units and testing families. Cells and libraries",
  "are not independent donors. The arrow denotes analysis, not a causal or temporal relationship."
), file.path(out, paste0(stem, "_caption.txt")))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
