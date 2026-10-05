#!/usr/bin/env Rscript
# KEY MESSAGE: Shared-signal genetic candidates retain distinct phenotype support.
# Genetics-only candidate; does not read or adopt a bulk disease-effect call.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: fig1c_corrected_genetics.R <evidence_class_table.tsv> <r2_summary.json> <output.pdf>")
}
project_root <- Sys.getenv("MASLD_PROJECT_ROOT", getwd())
source(file.path(project_root, "scripts/figures/publication_theme.R"))
classes <- fread(args[[1L]])
summary <- fromJSON(args[[2L]])
required <- c("symbol", "genetic_confidence", "genetic_trait_scope", "analysis_release_id")
stopifnot(all(required %in% names(classes)))
genes <- classes[genetic_confidence == "susie"]
stopifnot(!anyNA(genes$symbol), all(nzchar(genes$symbol)), !anyDuplicated(genes$symbol),
          all(genes$genetic_trait_scope %in% c("enzyme", "direct_disease", "both")),
          nrow(genes) == summary$tier12$tier12_susie_gt_0.5,
          uniqueN(genes$analysis_release_id) == 1L)
trait <- genes[, .(n_genes = .N), by = genetic_trait_scope]
trait <- merge(data.table(genetic_trait_scope = c("enzyme", "both", "direct_disease")),
               trait, by = "genetic_trait_scope", all.x = TRUE, sort = FALSE)
trait[is.na(n_genes), n_genes := 0L]
trait[, label := fcase(genetic_trait_scope == "enzyme", "Liver enzyme only",
                       genetic_trait_scope == "both", "Both",
                       default = "Direct disease/PDFF only")]
trait[, label := factor(label, levels = c("Liver enzyme only", "Both", "Direct disease/PDFF only"))]
trait[, `:=`(analysis_release_id = unique(genes$analysis_release_id),
             support_definition = "any_eligible_trait_support", release_state = "candidate_not_adopted")]
stopifnot(sum(trait$n_genes) == nrow(genes))

out <- args[[3L]]
source_path <- file.path(dirname(out), "fig1c_trait_scope_source.tsv")
manifest_path <- file.path(dirname(out), "fig1c_input_manifest.tsv")
caption_path <- file.path(dirname(out), "fig1c_caption.txt")
if (any(file.exists(c(out, source_path, manifest_path, caption_path)))) {
  stop("Refusing to overwrite an existing Figure 1C candidate")
}
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
plot <- ggplot(trait, aes(x = n_genes, y = label, fill = label)) +
  geom_col(width = 0.58, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(n_genes)), hjust = -0.15, size = 6 / .pt) +
  scale_fill_manual(values = c("Liver enzyme only" = "#F4A674", "Both" = "#00695C",
                               "Direct disease/PDFF only" = "#1565C0")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Primary SuSiE genes", y = NULL) +
  theme_masld(base_size = 6) + theme(panel.grid.major.y = element_blank())
ggsave(out, plot, width = 4, height = 2.25, device = grDevices::pdf,
       useDingbats = FALSE, bg = "white")
fwrite(trait, source_path, sep = "\t")
sha256 <- function(path) strsplit(system2("sha256sum", shQuote(normalizePath(path)),
                                        stdout = TRUE)[[1L]], "[[:space:]]+")[[1L]][[1L]]
inputs <- data.table(role = c("evidence_classes", "coloc_summary", "producer"),
                     path = normalizePath(c(args[[1L]], args[[2L]],
                                            file.path(project_root, "scripts/figures/fig1c_corrected_genetics.R"))))
inputs[, sha256 := vapply(path, sha256, character(1))]
fwrite(inputs, manifest_path, sep = "\t")
writeLines(sprintf(
  "Figure 1C candidate: %d primary tier-1/2 SuSiE genes partitioned by any eligible direct-disease/liver-fat and liver-enzyme support. This is not the strongest-trait-driver partition. Corrected genetics candidate; no synchronized release adoption.",
  nrow(genes)), caption_path)
writeLines(capture.output(sessionInfo()), file.path(dirname(out), "sessionInfo.txt"))
cat("FIG1C_CANDIDATE\tPASS\n")
print(trait)
