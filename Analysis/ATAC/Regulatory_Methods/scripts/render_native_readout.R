#!/usr/bin/env Rscript
# Explicitly bounded native-method diagnostics; no analysis is fitted here.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, !file.exists(args[2L]))
input <- args[1L]
state <- jsonlite::fromJSON(file.path(input, "completion.json"))
stopifnot(identical(state$status, "completed"),
  !length(list.files(input, pattern = "^INVALIDATED.*\\.json$")))
d <- fread(file.path(input, "exploratory_comparisons_context.tsv"))
stopifnot(all(d$gene_id == "ENSG00000000419"), all(d$cohort == "GSE296875"),
  all(d$lineage == "hepatocyte"), !any(d$predictive_improvement))
dir.create(args[2L], recursive = TRUE)
palette <- c(LINGER = "#C2185B", TRIPOD = "#1565C0")
manifest <- list()
for (target in c("promoter_accessibility", "rna")) {
  baseline_name <- if (target == "rna") "promoter_only" else "intercept"
  choices <- list(
    LINGER_native_bounded = c(baseline_name, "LINGER_native_bounded_pairing_shuffle",
      "LINGER_reference_RNA_only", "LINGER_reference_prior_only"),
    TRIPOD_Xt_L1 = c(baseline_name, "TRIPOD_Xt_L1_pairing_shuffle",
      "TRIPOD_reference_RNA_only", "TRIPOD_reference_prior_only"))
  selected <- rbindlist(lapply(names(choices), function(candidate)
    d[endpoint == target & method == candidate & baseline %in% choices[[candidate]]]))
  stopifnot(nrow(selected) == 8L, all(selected$n_donors == 39L),
    all(selected$status == "evaluated"), all(is.finite(selected$ci_low)),
    all(is.finite(selected$ci_high)), all(selected$ci_low <= selected$improvement),
    all(selected$ci_high >= selected$improvement))
  selected[, algorithm := fifelse(method == "LINGER_native_bounded", "LINGER", "TRIPOD")]
  selected[, comparator := fcase(baseline == baseline_name,
    if (target == "rna") "Promoter only" else "Intercept",
    grepl("shuffle", baseline), "Pairing shuffle",
    grepl("RNA_only", baseline), "RNA-only proxy",
    default = "Prior-only proxy")]
  selected[, label := paste(algorithm, "vs", comparator)]
  selected[, label := factor(label, levels = rev(label))]
  p <- ggplot(selected, aes(improvement, label, color = algorithm)) +
    geom_vline(xintercept = 0, color = "#9E9E9E", linewidth = .3) +
    geom_segment(aes(x = ci_low, xend = ci_high, yend = label), linewidth = .4) +
    geom_point(size = 1.3) + scale_color_manual(values = palette) +
    labs(x = "Comparator - candidate standardized squared error", y = NULL,
      title = paste("DPM1 held-out donor prediction:", target),
      subtitle = "GSE296875 hepatocytes; one gene; 39 donors; five outer folds\nBounded technical comparison; no method-selection claim",
      caption = "2,000 donor-bootstrap 95% intervals; positive values favor candidate.\nRNA/prior references are proxies, not native architectural ablations.\nNative source/API/attribution adaptations are recorded separately.") +
    theme_classic(base_size = 6, base_family = "Helvetica") +
    theme(text = element_text(size = 6, face = "plain"),
      plot.title = element_text(size = 6, face = "plain"),
      plot.subtitle = element_text(size = 6, face = "plain"),
      plot.caption = element_text(size = 6, face = "plain", hjust = 0),
      axis.text = element_text(size = 6, color = "black"),
      legend.position = "none", plot.caption.position = "plot")
  file <- file.path(args[2L], paste0("native_DPM1_", target, ".pdf"))
  ggsave(file, p, width = 4.8, height = 3.4, device = grDevices::pdf, useDingbats = FALSE)
  fwrite(selected[, .(cohort, lineage, lineage_state, gene_id, endpoint, method,
    baseline, n_donors, improvement, ci_low, ci_high)],
    file.path(args[2L], paste0("panel_values_", target, ".tsv")), sep = "\t")
  manifest[[target]] <- data.table(panel = basename(file), comparisons = nrow(selected),
    donors = 39L, genes = 1L, inference = "exploratory_only")
}
fwrite(rbindlist(manifest), file.path(args[2L], "render_manifest.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(args[2L], "sessionInfo.txt"))
