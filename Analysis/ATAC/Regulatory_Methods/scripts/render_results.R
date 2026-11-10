#!/usr/bin/env Rscript
# Candidate diagnostics only; no scientific analysis is inferred from missing tables.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
grDevices::pdf.options(useDingbats = FALSE)
source(file.path(Sys.getenv("HOME"), "publication_color_themes.R"))

require_true <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
FONT_SIZE <- 6
BASELINE <- "#9E9E9E"
METHOD_COLORS <- c(voom_qw = BASELINE, edger_ql = "#1565C0")

panel_theme <- function() {
  theme_classic(base_size = FONT_SIZE, base_family = "Helvetica") +
    theme(text = element_text(size = FONT_SIZE, face = "plain"),
      plot.title = element_text(size = FONT_SIZE, face = "plain"),
      plot.subtitle = element_text(size = FONT_SIZE, face = "plain"),
      plot.caption = element_text(size = FONT_SIZE, face = "plain", hjust = 0),
      plot.caption.position = "plot",
      axis.title = element_text(size = FONT_SIZE, face = "plain"),
      axis.text = element_text(size = FONT_SIZE, face = "plain", color = "black"),
      legend.title = element_text(size = FONT_SIZE, face = "plain"),
      legend.text = element_text(size = FONT_SIZE, face = "plain"),
      strip.text = element_text(size = FONT_SIZE, face = "plain"),
      strip.background = element_blank(), plot.margin = margin(4, 4, 4, 4))
}

read_existing <- function(root, name, columns = NULL, optional_columns = character()) {
  path <- file.path(root, name)
  if (!file.exists(path) || file.info(path)$size < 2) return(NULL)
  header <- names(fread(path, nrows = 0L))
  if (!is.null(columns) && !all(columns %in% header)) return(NULL)
  selected <- if (is.null(columns)) NULL else unique(c(columns, intersect(optional_columns, header)))
  result <- fread(path, select = selected)
  attr(result, "source_path") <- normalizePath(path)
  result
}

prepare_effect_scatter <- function(d, x, y, maximum_display = 5000L) {
  require_true(!anyDuplicated(d$peak_coordinate), "Duplicated peak in effect scatter")
  valid <- is.finite(d[[x]]) & is.finite(d[[y]])
  paired <- d[valid]
  n <- nrow(paired)
  if (!n) return(NULL)
  rho <- if (n > 2L && sd(paired[[x]]) > 0 && sd(paired[[y]]) > 0)
    cor(paired[[x]], paired[[y]], method = "spearman") else NA_real_
  same_sign <- mean(sign(paired[[x]]) == sign(paired[[y]]))
  setorder(paired, peak_coordinate)
  shown <- paired[unique(round(seq(1L, n, length.out = min(n, maximum_display))))]
  list(shown = shown, n = n, rho = rho, same_sign = same_sign)
}

parse_render_args <- function(args) {
  values <- list()
  require_true(length(args) >= 4L && length(args) %% 2L == 0L,
    "Usage: render_results.R [--da ROOT] [--network ROOT] --output NEWDIR")
  for (i in seq.int(1L, length(args), by = 2L)) {
    require_true(args[i] %in% c("--da", "--network", "--output"), "Unknown argument")
    values[[substring(args[i], 3L)]] <- args[i + 1L]
  }
  require_true((!is.null(values$da) || !is.null(values$network)) && !is.null(values$output),
    "--output and at least one of --da/--network required")
  values
}

render_results <- function(args) {
  if (!is.null(args$da)) {
    require_true(dir.exists(args$da), "Missing DA directory")
    require_true(file.exists(file.path(args$da, "run_summary.json")), "DA run has no completed summary")
  }
  require_true(!file.exists(args$output), "Refusing to overwrite output")
  dir.create(args$output, recursive = TRUE)
  manifest <- list()
  record <- function(panel, source, status, reason = "") {
    manifest[[length(manifest) + 1L]] <<- data.table(panel = panel, source = source,
      status = status, reason = reason)
  }
  save_panel <- function(p, name, source, height = 2.8, width = 3.5) {
    path <- file.path(args$output, paste0(name, ".pdf"))
    grDevices::pdf(path, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
    tryCatch(print(p + panel_theme()), finally = grDevices::dev.off())
    record(basename(path), source, "rendered")
  }
  groups <- function(d, columns) split(d, interaction(d[, ..columns], drop = TRUE))
  suffix <- function(d, columns) paste(vapply(columns, function(x) as.character(d[[x]][1L]), ""), collapse = "_")

  if (is.null(args$da)) {
    render_network(args$network, save_panel, record)
    fwrite(rbindlist(manifest), file.path(args$output, "render_manifest.tsv"), sep = "\t")
    writeLines(capture.output(sessionInfo()), file.path(args$output, "sessionInfo.txt"))
    return(invisible(NULL))
  }

  summary <- read_existing(args$da, "cohort_summary.tsv")
  if (!is.null(summary) && all(c("cohort", "n_tested_peaks", "method", "lineage", "representation", "filter") %in% names(summary))) {
    source <- attr(summary, "source_path")
    for (d in groups(summary[cohort != "cross_cohort" & is.finite(n_tested_peaks)], c("representation", "filter"))) {
      # KEY MESSAGE: Report the number of peaks actually tested for each cohort and lineage.
      p <- ggplot(d, aes(n_tested_peaks, lineage, color = method, shape = method)) +
        geom_point(position = position_dodge(width = 0.35), size = 1.5) + facet_wrap(~cohort, ncol = 1) +
        scale_color_manual(values = METHOD_COLORS) + labs(x = "Tested peaks (count)", y = NULL,
          title = paste("Accessibility test coverage:", d$representation[1L], d$filter[1L]),
          caption = "Coverage describes assay observability; it is not accuracy.")
      save_panel(p, paste0("da_coverage_", suffix(d, c("representation", "filter"))), source)
    }
  } else record("da_coverage", file.path(args$da, "cohort_summary.tsv"), "skipped", "missing_table_or_schema")

  null <- read_existing(args$da, "null_diagnostics.tsv")
  needed <- c("n_tested_peaks", "n_nominal_p05", "n_bh_q05", "method", "lineage", "cohort", "representation", "filter")
  if (!is.null(null) && all(needed %in% names(null))) {
    source <- attr(null, "source_path")
    for (d in groups(null[is.finite(n_tested_peaks) & n_tested_peaks > 0], c("representation", "filter"))) {
      label <- if ("null_model" %in% names(d) && all(grepl("NB", d$null_model)))
        "Fitted negative-binomial simulation; model diagnostic, not biological ground truth." else
        "Donor relabeling diagnostic; exchangeability assumptions apply, not biological ground truth."
      for (metric in c("n_nominal_p05", "n_bh_q05")) {
        d[, fraction := get(metric) / n_tested_peaks]
        # KEY MESSAGE: Display the observed rejection fractions under the recorded diagnostic null.
        p <- ggplot(d, aes(method, fraction, color = method, shape = method)) +
          geom_point(position = position_dodge(width = 0.2), size = 1, alpha = 0.55) +
          facet_grid(lineage ~ cohort) + scale_color_manual(values = METHOD_COLORS) +
          labs(x = NULL, y = if (metric == "n_nominal_p05") "Fraction with nominal p < 0.05" else "Fraction with BH q < 0.05",
            title = "Model calibration diagnostic", caption = label) + theme(legend.position = "none")
        if (metric == "n_nominal_p05") p <- p + geom_hline(yintercept = 0.05, linetype = 2, color = BASELINE, linewidth = 0.3)
        save_panel(p, paste0("da_null_", metric, "_", suffix(d, c("representation", "filter"))), source, width = 4.4, height = 3.4)
      }
    }
  } else record("da_null", file.path(args$da, "null_diagnostics.tsv"), "skipped", "missing_table_or_no_completed_null")

  needed <- c("peak_coordinate", "lineage", "representation", "filter", "method", "logFC_gse244832", "logFC_gse281367")
  replication <- read_existing(args$da, "replication.tsv.gz", needed,
    optional_columns = c("qvalue_gse244832", "qvalue_gse281367"))
  if (!is.null(replication)) {
    source <- attr(replication, "source_path")
    for (d in groups(replication, c("lineage", "representation", "filter"))) {
      for (cohort in c("gse244832", "gse281367")) {
        effect <- paste0("logFC_", cohort)
        x <- d[method == "voom_qw", .(peak_coordinate, voom = get(effect))]
        y <- d[method == "edger_ql", .(peak_coordinate, edger = get(effect))]
        require_true(!anyDuplicated(x$peak_coordinate) && !anyDuplicated(y$peak_coordinate), "Duplicated method/feature rows")
        paired <- merge(x, y, by = "peak_coordinate")[is.finite(voom) & is.finite(edger)]
        if (!nrow(paired)) next
        n <- nrow(paired)
        rho <- if (n > 2 && sd(paired$voom) > 0 && sd(paired$edger) > 0) cor(paired$voom, paired$edger, method = "spearman") else NA_real_
        # Coordinate order gives a reproducible display sample independent of effect or significance.
        setorder(paired, peak_coordinate)
        shown <- paired[unique(round(seq(1, n, length.out = min(n, 5000L))))]
        # KEY MESSAGE: Compare raw method effects on the same peaks without significance selection.
        p <- ggplot(shown, aes(voom, edger)) + geom_abline(slope = 1, intercept = 0, color = BASELINE, linewidth = 0.3) +
          geom_point(size = 0.3, alpha = 0.35, color = "#1565C0") + coord_equal() +
          labs(x = "voom quality weights: MASH - normal log2 effect", y = "edgeR QL: MASH - normal log2 effect",
            title = paste0("Raw accessibility effect agreement: ", cohort, "\n", d$lineage[1L], " | ", d$representation[1L], " | ", d$filter[1L]),
            caption = sprintf("All shared peaks: %s; Spearman rho: %.3f.\nDisplayed: %s; %s counts.\nMethod agreement is not ground-truth validation.",
              format(n, big.mark = ","), rho, nrow(shown), d$representation[1L]))
        save_panel(p, paste0("da_method_effects_", cohort, "_", suffix(d, c("lineage", "representation", "filter"))), source)
      }
    }
    for (d in groups(replication, c("lineage", "representation", "filter", "method"))) {
      stats <- prepare_effect_scatter(d, "logFC_gse244832", "logFC_gse281367")
      if (is.null(stats)) next
      # KEY MESSAGE: Show whether independent donor cohorts agree on disease effects over all shared peaks.
      p <- ggplot(stats$shown, aes(logFC_gse244832, logFC_gse281367)) +
        geom_abline(slope = 1, intercept = 0, color = BASELINE, linewidth = 0.3) +
        geom_point(size = 0.3, alpha = 0.4, color = METHOD_COLORS[[d$method[1L]]]) + coord_equal() +
        labs(x = "GSE244832: MASH - normal log2 effect", y = "GSE281367: MASH - normal log2 effect",
          title = paste0("Candidate cross-cohort accessibility effects\n", d$lineage[1L], " | ", d$representation[1L], " | ", d$filter[1L], " | ", d$method[1L]),
          caption = sprintf("All shared peaks: %s; Spearman rho: %.3f.\nSame sign: %.1f%%; displayed: %s (coordinate sample).\nNo significance selection; candidate comparison only.\nBoth cohorts contributed to outcome-blind peak ascertainment.",
            format(stats$n, big.mark = ","), stats$rho, 100 * stats$same_sign, format(nrow(stats$shown), big.mark = ",")))
      save_panel(p, paste0("da_cohort_effects_", suffix(d, c("lineage", "representation", "filter", "method"))), source)
    }
    adjusted <- read_existing(args$da, "covariate_sensitivity.tsv.gz",
      c("peak_coordinate", "adjusted_logFC", "adjusted_qvalue", "covariate", "status", "cohort", "lineage", "representation", "filter", "method"))
    if (!is.null(adjusted)) {
      adjusted_source <- paste(source, attr(adjusted, "source_path"), sep = ";")
      batch <- adjusted[status == "testable" & grepl("^batch", covariate, ignore.case = TRUE)]
      if (!nrow(batch)) record("da_batch_effects", attr(adjusted, "source_path"), "skipped", "no_testable_recorded_batch")
      for (d in groups(batch, c("cohort", "lineage", "representation", "filter", "method", "covariate"))) {
        keys <- c("lineage", "representation", "filter", "method")
        lookup <- unique(d[, ..keys])
        effect <- paste0("logFC_", tolower(d$cohort[1L]))
        qvalue <- paste0("qvalue_", tolower(d$cohort[1L]))
        if (!effect %in% names(replication)) next
        raw <- replication[lookup, on = keys, nomatch = 0L]
        raw <- raw[, .(peak_coordinate, condition_only_logFC = get(effect),
          condition_only_qvalue = if (qvalue %in% names(raw)) get(qvalue) else NA_real_)]
        paired <- merge(raw, d[, .(peak_coordinate, adjusted_logFC, adjusted_qvalue)], by = "peak_coordinate")
        stats <- prepare_effect_scatter(paired, "condition_only_logFC", "adjusted_logFC")
        if (is.null(stats)) next
        finite <- paired[is.finite(condition_only_logFC) & is.finite(adjusted_logFC)]
        count_label <- if (qvalue %in% names(replication)) sprintf("BH q < 0.05: condition only %s; adjusted %s.",
          format(sum(finite$condition_only_qvalue < .05, na.rm = TRUE), big.mark = ","),
          format(sum(finite$adjusted_qvalue < .05, na.rm = TRUE), big.mark = ",")) else "BH comparison unavailable in the supplied source."
        # KEY MESSAGE: Show the change in disease effects after adjustment for recorded technical batch.
        p <- ggplot(stats$shown, aes(condition_only_logFC, adjusted_logFC)) +
          geom_abline(slope = 1, intercept = 0, color = BASELINE, linewidth = 0.3) +
          geom_point(size = 0.3, alpha = 0.4, color = METHOD_COLORS[[d$method[1L]]]) + coord_equal() +
          labs(x = "Condition only: MASH - normal log2 effect", y = "Batch adjusted: MASH - normal log2 effect",
            title = paste0("Candidate batch sensitivity: ", d$cohort[1L], "\n", d$lineage[1L], " | ", d$representation[1L], " | ", d$filter[1L], " | ", d$method[1L]),
            caption = sprintf("Shared peaks: %s; Spearman rho: %.3f.\n%s\nAdjustment: %s; displayed: %s (coordinate sample).\nOriginal BH families retained; no significance selection.",
              format(stats$n, big.mark = ","), stats$rho, count_label, d$covariate[1L], format(nrow(stats$shown), big.mark = ",")))
        save_panel(p, paste0("da_batch_effects_", suffix(d, c("cohort", "lineage", "representation", "filter", "method", "covariate"))), adjusted_source)
      }
    } else record("da_batch_effects", file.path(args$da, "covariate_sensitivity.tsv.gz"), "skipped", "missing_table_or_no_testable_batch")
    rm(replication)
  } else record("da_method_effects", file.path(args$da, "replication.tsv.gz"), "skipped", "missing_table")

  lodo <- read_existing(args$da, "lodo_summary.tsv")
  needed <- c("status", "effect_spearman", "method", "cohort", "lineage", "representation", "filter")
  if (!is.null(lodo) && all(needed %in% names(lodo))) {
    source <- attr(lodo, "source_path")
    for (d in groups(lodo[status == "testable" & is.finite(effect_spearman)], c("representation", "filter"))) {
      # KEY MESSAGE: Show how much accessibility effects change when one biological donor is omitted.
      p <- ggplot(d, aes(method, effect_spearman, color = method, shape = method)) + geom_point(size = 1, alpha = 0.6) +
        facet_grid(lineage ~ cohort) + scale_color_manual(values = METHOD_COLORS) +
        labs(x = NULL, y = "Spearman rho with full fit", title = "Donor-omission sensitivity",
          caption = "Each point omits one donor; overlapping subsets are not independent replicates.") + theme(legend.position = "none")
      save_panel(p, paste0("da_donor_influence_", suffix(d, c("representation", "filter"))), source, width = 4.4, height = 3.4)
    }
  } else record("da_donor_influence", file.path(args$da, "lodo_summary.tsv"), "skipped", "missing_table_or_not_run")

  if (is.null(args$network)) record("network", "", "skipped", "no_network_directory_supplied")
  else render_network(args$network, save_panel, record)
  fwrite(rbindlist(manifest), file.path(args$output, "render_manifest.tsv"), sep = "\t")
  writeLines(capture.output(sessionInfo()), file.path(args$output, "sessionInfo.txt"))
}

render_network <- function(root, save_panel, record) {
  # Network table readers are deliberately separate from DA effect units.
  require_true(dir.exists(root), "Missing network directory")
  require_true(file.exists(file.path(root, "completion.json")), "Network run has no completed summary")
  completion <- jsonlite::fromJSON(file.path(root, "completion.json"))
  require_true(identical(completion$status, "completed") &&
    !length(list.files(root, pattern = "^INVALIDATED.*\\.json$")), "Network run is incomplete or invalidated")
  settings <- if (file.exists(file.path(root, "run_settings.json"))) jsonlite::fromJSON(file.path(root, "run_settings.json")) else list()
  maximum <- settings$max_genes
  bounded <- length(maximum) == 1L && is.finite(maximum) && maximum > 0
  scope <- if (isTRUE(settings$pilot) || isTRUE(completion$pilot) || identical(settings$analysis_scope, "technical_pilot")) {
    "Technical pilot"
  } else if (bounded || identical(settings$analysis_scope, "bounded_gene_development") || identical(completion$interpretation, "bounded_gene_development")) {
    "Bounded development"
  } else if (identical(settings$analysis_scope, "genome_wide")) "Candidate genome-wide analysis" else "Candidate analysis; scope unrecorded"
  if (isTRUE(settings$technical_comparison_only)) {
    scope <- paste0("Bounded technical comparison; no method-selection claim\n",
      settings$cicero_distance_samples, " Cicero calibration windows")
  }
  if (bounded) scope <- paste0(scope, "; maximum ", maximum, " genes")
  comparisons <- read_existing(root, "comparisons.tsv")
  needed <- c("cohort", "lineage", "endpoint", "method", "baseline", "improvement", "ci_low", "ci_high", "n_donors", "status")
  if (!is.null(comparisons) && all(needed %in% names(comparisons))) {
    source <- attr(comparisons, "source_path")
    d <- comparisons[status == "evaluated" & is.finite(improvement) & is.finite(ci_low) & is.finite(ci_high)]
    if (nrow(d)) {
      method_label <- function(x) {
        labels <- c(cicero = "Donor-weighted Cicero", directnet_adapted = "Adapted DIRECT-NET",
          distance = "Distance", intercept = "Intercept", promoter_only = "Promoter only")
        result <- unname(labels[x])
        result[is.na(result)] <- x[is.na(result)]
        result
      }
      d[, comparison := paste(method_label(method), "vs", method_label(baseline))]
      d[, lineage_label := fifelse(cohort == "GSE296875" & lineage == "stellate",
        "Mesenchymal (source label)", lineage)]
      for (g in split(d, interaction(d$cohort, d$endpoint, drop = TRUE))) {
      # KEY MESSAGE: Report held-out donor prediction improvements with their recorded donor bootstrap intervals.
      p <- ggplot(g, aes(improvement, comparison)) + geom_vline(xintercept = 0, color = BASELINE, linewidth = 0.3) +
        geom_segment(aes(x = ci_low, xend = ci_high, yend = comparison), linewidth = 0.4, color = "#1565C0") +
        geom_point(size = 1.3, color = "#1565C0") + facet_wrap(~lineage_label, ncol = 1, scales = "free_y") +
        labs(x = "Baseline - candidate standardized squared error", y = NULL,
          title = paste("Held-out donor prediction:", g$cohort[1L], g$endpoint[1L]),
          subtitle = scope,
          caption = "Recorded donor bootstrap 95% intervals; positive values favor candidate.\nPrediction does not establish causal regulation.")
      name <- paste(g$cohort[1L], g$endpoint[1L], sep = "_")
      save_panel(p, paste0("network_prediction_gain_", name), source, width = 4.4, height = 3.2)
      # KEY MESSAGE: Show the biological-donor coverage supporting each prediction comparison.
      p <- ggplot(g, aes(n_donors, comparison)) + geom_point(size = 1.3, color = BASELINE) +
        scale_x_continuous(limits = c(0, max(g$n_donors) * 1.1),
          breaks = function(x) { b <- pretty(x); b[b == floor(b)] }) +
        facet_wrap(~lineage_label, ncol = 1, scales = "free_y") + labs(x = "Held-out biological donors (count)", y = NULL,
          title = paste("Prediction coverage:", g$cohort[1L], g$endpoint[1L]), subtitle = scope)
      save_panel(p, paste0("network_donor_coverage_", name), source, width = 4.4, height = 3.2)
      }
    } else record("network_prediction_gain", source, "skipped", "no_evaluated_comparisons_with_intervals")
  } else record("network_prediction_gain", file.path(root, "comparisons.tsv"), "skipped", "missing_table_or_schema")
  for (directory in list.dirs(root, recursive = FALSE, full.names = TRUE)) {
    coverage <- read_existing(directory, "coverage.tsv")
    if (is.null(coverage) || !all(c("outer_fold", "gene_id", "status", "n_candidates") %in% names(coverage))) next
    source <- attr(coverage, "source_path")
    d <- coverage[, .(gene_fold_count = .N), by = status]
    context_label <- if (basename(directory) == "GSE296875_stellate")
      "GSE296875 Mesenchymal (source label)" else basename(directory)
    # KEY MESSAGE: Distinguish tested from untestable gene-fold combinations.
    p <- ggplot(d, aes(gene_fold_count, reorder(status, gene_fold_count))) + geom_point(size = 1.4, color = BASELINE) +
      scale_x_continuous(limits = c(0, max(d$gene_fold_count) * 1.1),
        breaks = function(x) { b <- pretty(x); b[b == floor(b)] }) +
      labs(x = "Gene-fold combinations (count)", y = NULL, title = paste("Cis-context coverage:", context_label),
        subtitle = scope,
        caption = "Folds overlap in training donors; these counts are not independent biological replicates.")
    save_panel(p, paste0("network_gene_coverage_", basename(directory)), source)
  }
}

if (sys.nframe() == 0L) render_results(parse_render_args(commandArgs(trailingOnly = TRUE)))
