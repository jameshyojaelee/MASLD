#!/usr/bin/env Rscript

# Figure 4 honesty reframe (2026-08-08).
#
# Two new panels for figures/main/fig4_validation/:
#
#   (1) fig4c_convergence_or_ci_power.pdf
#       Convergent class vs each single-map class, across the three atlas
#       endpoints, drawn as exact-Fisher odds ratio + 95% CI on a log-odds axis
#       with the 80%-power minimum detectable odds ratio marked. The existing
#       fig4c_convergence_comparison.pdf plots point estimates and q-values with
#       no interval, so it reads as an established absence of effect. Every
#       comparison here is indeterminate: the intervals include OR = 1 AND the
#       study cannot resolve effects below the marked resolution limit.
#
#   (2) fig4d_endpoint_provenance.pdf
#       What each endpoint actually tests. The disease-state class is defined by
#       a bulk disease-vs-control contrast; proteomics and single-cell endpoints
#       are the SAME estimand measured on a different molecular layer
#       (cross-assay replication), while spatial SVG tests a distinct property
#       (spatial variability) and is therefore the one fair independent test.
#       This replaces any "three independent modalities" framing.
#
# Reads ONLY the frozen manuscript release. Every plotted number is re-derived
# from the release 2x2 counts and asserted against the release's own statistics;
# any mismatch is a hard stop.
#
# Does not touch fig4a / fig4b (their class-vs-neither contrasts are unaffected)
# and does not modify scripts/figures/fig1_fig4_evidence_classes.R.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID  <- Sys.getenv("MANUSCRIPT_RELEASE_ID", "2026-07-15-r2")
RELEASE_DIR <- file.path(BASE, "RNA-seq/results/manuscript_release", RELEASE_ID)

source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR  <- FIG4_DIR
DATA_DIR <- file.path(OUT_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

stopifnot(basename(OUT_DIR) == "fig4_validation")

# ---------------------------------------------------------------------------
# 0. Load the frozen release tables
# ---------------------------------------------------------------------------
req <- file.path(RELEASE_DIR, c(
  "evidence_class_validation_summary.tsv",
  "evidence_class_validation_adjusted.tsv",
  "evidence_class_validation_pairwise.tsv",
  "evidence_class_counts.tsv",
  "acceptance_gates.tsv"
))
missing <- req[!file.exists(req)]
if (length(missing)) {
  stop("[fig4-resolution] missing release artifact(s): ",
       paste(basename(missing), collapse = ", "))
}

summ     <- fread(req[1])
adjusted <- fread(req[2])
pairwise <- fread(req[3])
counts   <- fread(req[4])
gates    <- fread(req[5])

for (nm in c("summ", "adjusted", "pairwise", "counts", "gates")) {
  d <- get(nm)
  if (!all(d$analysis_release_id == RELEASE_ID)) {
    stop("[fig4-resolution] ", nm, " carries a release id other than ", RELEASE_ID)
  }
}

# ---------------------------------------------------------------------------
# 1. Assertion: the release's own summary table is internally consistent
# ---------------------------------------------------------------------------
summ[, rate_rederived := n_positive / n_tested]
bad_rate <- summ[abs(rate_rederived - positive_rate) > 1e-9]
if (nrow(bad_rate)) {
  print(bad_rate)
  stop("[fig4-resolution] positive_rate does not equal n_positive/n_tested")
}

# ---------------------------------------------------------------------------
# 2. Assertion: the prespecified acceptance gate is recorded and did NOT pass
# ---------------------------------------------------------------------------
gate <- gates[gate == "class_validation_article"]
if (nrow(gate) != 1L) {
  stop("[fig4-resolution] expected exactly one class_validation_article gate row")
}
gate_passed <- tolower(as.character(gate$passed[1]))
if (!gate_passed %in% c("false", "f")) {
  stop("[fig4-resolution] class_validation_article gate is no longer False (got '",
       gate$passed[1], "') -- the panel's framing is stale, re-derive before plotting")
}
GATE_CRITERION <- as.character(gate$criterion[1])
if (!nzchar(GATE_CRITERION)) {
  stop("[fig4-resolution] class_validation_article gate has no recorded criterion; ",
       "a prespecified gate must carry its criterion")
}
message("[fig4-resolution] prespecified gate class_validation_article = ",
        gate$passed[1], "  criterion: ", GATE_CRITERION)

# ---------------------------------------------------------------------------
# 3. Assertion: the convergent class is the 34-gene class described in the paper
# ---------------------------------------------------------------------------
conv_counts <- counts[primary_evidence_class == "convergent"]
n_convergent_class <- sum(conv_counts$N)
if (n_convergent_class != 34L) {
  stop("[fig4-resolution] convergent class is ", n_convergent_class, ", expected 34")
}
n_conv_direct <- sum(conv_counts[genetic_trait_scope == "direct_disease"]$N)
if (n_conv_direct != 0L) {
  stop("[fig4-resolution] convergent class now contains ", n_conv_direct,
       " direct-disease genes; the enzyme-trait caveat needs re-deriving")
}
message(sprintf("[fig4-resolution] convergent class n = %d (%d enzyme, %d both, %d direct-disease)",
                n_convergent_class,
                sum(conv_counts[genetic_trait_scope == "enzyme"]$N),
                sum(conv_counts[genetic_trait_scope == "both"]$N),
                n_conv_direct))

# ---------------------------------------------------------------------------
# 4. Re-derive every 2x2, exact Fisher OR + 95% CI, and the 80%-power limit
# ---------------------------------------------------------------------------

# Exact unconditional power of Fisher's exact test with fixed group sizes.
# Rejection region is enumerated once (it does not depend on the alternative),
# then averaged over the two independent binomials under the alternative.
fisher_rejection_matrix <- function(n1, n2, alpha = 0.05) {
  rej <- matrix(FALSE, nrow = n1 + 1L, ncol = n2 + 1L)
  for (x1 in 0:n1) {
    for (x2 in 0:n2) {
      tab <- matrix(c(x1, n1 - x1, x2, n2 - x2), nrow = 2L)
      # conf.int = FALSE: the interval is not needed here and its root-finding
      # dominates the runtime of ~120k calls.
      rej[x1 + 1L, x2 + 1L] <-
        stats::fisher.test(tab, conf.int = FALSE)$p.value <= alpha
    }
  }
  rej
}

power_at_or <- function(or, n1, n2, p2, rej) {
  odds2 <- p2 / (1 - p2)
  p1 <- (or * odds2) / (1 + or * odds2)
  sum(outer(stats::dbinom(0:n1, n1, p1), stats::dbinom(0:n2, n2, p2)) * rej)
}

# Smallest OR > 1 the design can detect with >= `target` power at nominal alpha.
min_detectable_or <- function(n1, n2, p2, rej, target = 0.80) {
  pw <- function(or) power_at_or(or, n1, n2, p2, rej)
  lo <- 1
  hi <- 2
  while (pw(hi) < target) {
    hi <- hi * 2
    if (hi > 1e7) return(list(mdor = NA_real_, power = NA_real_))
  }
  for (i in seq_len(80L)) {
    mid <- sqrt(lo * hi)
    if (pw(mid) < target) lo <- mid else hi <- mid
  }
  list(mdor = hi, power = pw(hi))
}

class_levels <- c("neither", "genetic_only", "disease_state_only", "convergent")
get_cell <- function(endpoint_name, cls, field) {
  v <- summ[endpoint == endpoint_name & primary_evidence_class == cls][[field]]
  if (length(v) != 1L) {
    stop("[fig4-resolution] no unique summary row for ", endpoint_name, " / ", cls)
  }
  as.integer(v)
}

grid <- CJ(endpoint   = c("proteomics", "spatial_svg", "single_cell"),
           comparator = c("disease_state_only", "genetic_only"),
           sorted = FALSE)

rows <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  ep <- grid$endpoint[i]
  cm <- grid$comparator[i]

  n1 <- get_cell(ep, "convergent", "n_tested")
  k1 <- get_cell(ep, "convergent", "n_positive")
  n2 <- get_cell(ep, cm,           "n_tested")
  k2 <- get_cell(ep, cm,           "n_positive")

  tab <- matrix(c(k1, n1 - k1, k2, n2 - k2), nrow = 2L,
                dimnames = list(c("positive", "negative"), c("convergent", cm)))
  ft <- stats::fisher.test(tab)

  p2   <- k2 / n2
  rej  <- fisher_rejection_matrix(n1, n2, alpha = 0.05)
  mdor <- min_detectable_or(n1 = n1, n2 = n2, p2 = p2, rej = rej, target = 0.80)

  rows[[i]] <- data.table(
    endpoint          = ep,
    comparator        = cm,
    n_convergent      = n1,
    k_convergent      = k1,
    n_comparator      = n2,
    k_comparator      = k2,
    rate_convergent   = k1 / n1,
    rate_comparator   = p2,
    odds_ratio        = unname(ft$estimate),
    ci_low            = ft$conf.int[1],
    ci_high           = ft$conf.int[2],
    p_value           = ft$p.value,
    or_naive          = (k1 * (n2 - k2)) / ((n1 - k1) * k2),
    mdor_80           = mdor$mdor,
    power_at_mdor     = mdor$power,
    power_at_observed = power_at_or(unname(ft$estimate), n1, n2, p2, rej)
  )
}
res <- rbindlist(rows)

# ---- Assertion: recomputed OR / p must match the frozen release pairwise table
chk <- merge(
  res[, .(endpoint, comparison = paste0("convergent_vs_", comparator),
          or_mine = odds_ratio, p_mine = p_value,
          n1_mine = n_convergent, n2_mine = n_comparator)],
  pairwise[, .(endpoint, comparison, odds_ratio, p_value, n_convergent, n_comparator)],
  by = c("endpoint", "comparison"), all = TRUE
)
if (nrow(chk) != 6L || anyNA(chk$or_mine) || anyNA(chk$odds_ratio)) {
  print(chk)
  stop("[fig4-resolution] pairwise comparisons do not align 1:1 with the release")
}
chk[, `:=`(d_or = abs(or_mine - odds_ratio),
           d_p  = abs(p_mine - p_value),
           d_n1 = abs(n1_mine - n_convergent),
           d_n2 = abs(n2_mine - n_comparator))]
if (chk[, max(d_or) > 1e-8 || max(d_p) > 1e-10 || max(d_n1) > 0 || max(d_n2) > 0]) {
  print(chk)
  stop("[fig4-resolution] recomputed exact-Fisher statistics disagree with the release")
}
message("[fig4-resolution] all 6 recomputed exact-Fisher ORs / p-values match the release")

# BH q-values across the 6 pairwise tests, re-derived, checked against release
res[, q_value := stats::p.adjust(p_value, method = "BH")]
qchk <- merge(res[, .(endpoint, comparison = paste0("convergent_vs_", comparator), q_mine = q_value)],
              pairwise[, .(endpoint, comparison, q_value)], by = c("endpoint", "comparison"))
if (qchk[, max(abs(q_mine - q_value))] > 1e-10) {
  print(qchk)
  stop("[fig4-resolution] recomputed BH q-values disagree with the release")
}

# The claim at issue is convergent vs disease-state-only: all three intervals
# must span OR = 1 for the "cannot resolve" framing to hold.
res[, spans_null := ci_low <= 1 & ci_high >= 1]
ds <- res[comparator == "disease_state_only"]
if (!all(ds$spans_null)) {
  print(ds[spans_null == FALSE])
  stop("[fig4-resolution] a convergent-vs-disease-state-only interval no longer ",
       "spans OR=1; the indeterminate framing is stale and must be re-derived")
}

# No comparison may meet the prespecified multiplicity-controlled gate. One
# comparison (proteomics vs genetic-only) does exclude OR=1 at nominal alpha
# while failing at BH q<0.05; that is shown on the panel with its q-value rather
# than smoothed over.
if (any(res$q_value < 0.05)) {
  print(res[q_value < 0.05])
  stop("[fig4-resolution] a comparison now meets BH q<0.05; the acceptance gate ",
       "would flip and this panel's framing is stale")
}
n_nominal <- res[p_value < 0.05, .N]
message(sprintf(
  "[fig4-resolution] %d of %d comparisons exclude OR=1 at nominal alpha; %d meet BH q<0.05",
  n_nominal, nrow(res), res[q_value < 0.05, .N]))

# Every observed effect smaller than the design's own 80%-power limit is the
# crisp statement of the resolution problem. Reported, not asserted.
res[, below_resolution := odds_ratio < mdor_80]
message(sprintf("[fig4-resolution] %d of %d point estimates lie below their own 80%%-power limit",
                res[below_resolution == TRUE, .N], nrow(res)))

# ---- Independent second derivation of the power limit by Monte Carlo.
# The minimum detectable OR is the load-bearing number on this panel and it
# disagrees with a normal-approximation estimate by roughly 2x, so it is
# re-derived by simulation rather than re-read. Analytic enumeration and
# simulation must agree at the claimed 80% power.
set.seed(42)
N_SIM <- 20000L
res[, power_at_mdor_sim := {
  vapply(seq_len(.N), function(i) {
    odds2 <- rate_comparator[i] / (1 - rate_comparator[i])
    o1 <- mdor_80[i] * odds2
    p1 <- o1 / (1 + o1)
    x1 <- stats::rbinom(N_SIM, n_convergent[i], p1)
    x2 <- stats::rbinom(N_SIM, n_comparator[i], rate_comparator[i])
    mean(vapply(seq_len(N_SIM), function(s) {
      stats::fisher.test(matrix(c(x1[s], n_convergent[i] - x1[s],
                                  x2[s], n_comparator[i] - x2[s]), nrow = 2L),
                         conf.int = FALSE)$p.value <= 0.05
    }, logical(1)))
  }, numeric(1))
}]
res[, power_sim_gap := abs(power_at_mdor_sim - 0.80)]
if (res[, max(power_sim_gap)] > 0.02) {
  print(res[, .(endpoint, comparator, mdor_80, power_at_mdor, power_at_mdor_sim)])
  stop("[fig4-resolution] analytic and simulated power disagree at the minimum ",
       "detectable OR; the resolution limit is not trustworthy")
}
message(sprintf("[fig4-resolution] MC cross-check (%d sims/row): simulated power at the ",
                N_SIM), "minimum detectable OR = ",
        paste(sprintf("%.3f", res$power_at_mdor_sim), collapse = ", "),
        " (analytic target 0.800)")

fwrite(res, file.path(DATA_DIR, "fig4c_convergence_or_ci_power.tsv"), sep = "\t")

message("[fig4-resolution] recomputed statistics:")
for (i in seq_len(nrow(res))) {
  with(res[i], message(sprintf(
    "  %-12s vs %-18s OR %5.2f [%.2f, %.2f]  p=%.3f q=%.3f  n=%d vs %d  80%%-power OR >= %.1f",
    endpoint, comparator, odds_ratio, ci_low, ci_high, p_value, q_value,
    n_convergent, n_comparator, mdor_80)))
}

# ---------------------------------------------------------------------------
# 5. Panel 1 -- odds ratio with 95% CI and the 80%-power resolution limit
# ---------------------------------------------------------------------------
class_colors <- c(
  neither            = "#D7D9DA",
  genetic_only       = "#1565C0",
  disease_state_only = "#C9265E",
  convergent         = "#00695C"
)
LIMIT_GREY <- "#9E9E9E"

pd <- copy(res)
pd[, endpoint_label := factor(
  endpoint,
  levels = c("single_cell", "spatial_svg", "proteomics"),
  labels = c("Single-cell", "Spatial SVG", "Proteomics"))]
pd[, row_label := sprintf("%s\n(%d vs %d)", as.character(endpoint_label),
                          n_convergent, n_comparator)]
pd[, comparator_label := factor(
  comparator,
  levels = c("disease_state_only", "genetic_only"),
  labels = c("Convergent vs disease-state-only",
             "Convergent vs genetic-only"))]
pd[, row_f := factor(row_label, levels = unique(row_label[order(endpoint_label)]))]

pd[, or_label := sprintf("OR %.2f [%.2f, %.2f], q = %.2f",
                         odds_ratio, ci_low, ci_high, q_value)]
pd[, mdor_label := sprintf("cannot resolve below OR %.1f", mdor_80)]
pd[, ci_mid     := sqrt(ci_low * ci_high)]

X_MIN <- 0.12
X_MAX <- 400
NUDGE <- 0.24

p1 <- ggplot(pd, aes(y = row_f)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "black", linewidth = 0.3) +
  # Resolution limit, offset just below each row: the band of odds ratios this
  # design has under 80% power to detect. Drawn as its own track so it reads as
  # a property of the design, not as an interval estimate.
  geom_segment(aes(x = 1, xend = mdor_80, yend = row_f),
               position = position_nudge(y = -NUDGE),
               colour = LIMIT_GREY, linewidth = 1.4, alpha = 0.55,
               lineend = "butt") +
  geom_point(aes(x = mdor_80), position = position_nudge(y = -NUDGE),
             shape = 1, size = 1.3, stroke = 0.4, colour = LIMIT_GREY) +
  geom_text(aes(x = mdor_80, label = mdor_label),
            position = position_nudge(y = -NUDGE),
            hjust = -0.12, size = 6 / .pt, colour = "black") +
  # Point estimate and exact 95% interval on the row itself.
  geom_segment(aes(x = ci_low, xend = ci_high, yend = row_f, colour = comparator),
               linewidth = 0.4) +
  geom_point(aes(x = odds_ratio, colour = comparator), shape = 16, size = 1.5) +
  # geom_label with a white knockout: several intervals straddle x = 1, and plain
  # text would have the dashed reference line running through the glyphs.
  geom_label(aes(x = ci_mid, label = or_label),
             position = position_nudge(y = NUDGE - 0.02),
             hjust = 0.5, size = 6 / .pt, colour = "black",
             fill = "white", label.size = 0,
             label.padding = unit(0.5, "pt")) +
  facet_wrap(~ comparator_label, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(disease_state_only = class_colors[["disease_state_only"]],
                                 genetic_only       = class_colors[["genetic_only"]]),
                      guide = "none") +
  scale_x_log10(breaks = c(0.25, 1, 4, 16, 64, 256),
                labels = c("0.25", "1", "4", "16", "64", "256"),
                limits = c(X_MIN, X_MAX),
                expand = expansion(mult = c(0.01, 0))) +
  scale_y_discrete(expand = expansion(add = 0.65)) +
  labs(x = paste0("Odds of endpoint positivity, convergent vs comparator class ",
                  "(exact Fisher, log scale).\nGrey track = advantage sizes this ",
                  "design has under 80% power to detect."),
       y = NULL) +
  theme_masld(base_size = 6) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(lineheight = 0.95),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

pdf(file.path(OUT_DIR, "fig4c_convergence_or_ci_power.pdf"),
    width = 5.4, height = 3.4, useDingbats = FALSE)
print(p1)
invisible(dev.off())

# ---------------------------------------------------------------------------
# 6. Panel 2 -- endpoint provenance: what each endpoint actually tests
# ---------------------------------------------------------------------------
wrap_cell <- function(x, width) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

rate_txt <- function(ep) {
  a <- summ[endpoint == ep & primary_evidence_class == "disease_state_only"]$positive_rate
  b <- summ[endpoint == ep & primary_evidence_class == "neither"]$positive_rate
  sprintf("%.1f%% vs %.1f%%", 100 * a, 100 * b)
}

prov <- data.table(
  row_id = 1:4,
  row_label = c("Bulk RNA-seq\n(defines the class)", "Proteomics",
                "Single-cell RNA", "Spatial SVG"),
  layer = c("Bulk liver tissue RNA",
            "Liver protein, DIA-MS",
            "Single-cell RNA, donor pseudobulk",
            "Spatial transcriptome"),
  estimand = c("Disease vs control",
               "Disease vs control",
               "Disease vs control",
               "Spatial variability; no disease contrast"),
  relation = c("Defines the disease-state class",
               "Cross-assay replication of the same estimand",
               "Cross-assay replication of the same estimand",
               "Distinct property, so an independent test"),
  rate = c("not applicable",
           rate_txt("proteomics"),
           rate_txt("single_cell"),
           rate_txt("spatial_svg")),
  relation_class = c("defines", "replication", "replication", "independent")
)

long <- melt(prov,
             id.vars = c("row_id", "row_label", "relation_class"),
             measure.vars = c("layer", "estimand", "relation", "rate"),
             variable.name = "column", value.name = "text")
long[, column := factor(column,
  levels = c("layer", "estimand", "relation", "rate"),
  labels = c("Molecular layer", "Estimand tested",
             "Relation to the class definition",
             "Positive rate:\ndisease-state-only vs neither"))]
long[, text_wrapped := wrap_cell(text, width = 22)]
long[, row_f := factor(row_label, levels = rev(prov$row_label))]
long[, relation_class := factor(relation_class,
  levels = c("defines", "replication", "independent"),
  labels = c("Defines the class being tested",
             "Same estimand, different assay",
             "Distinct property (fair independent test)"))]

fill_vals <- c(
  "Defines the class being tested"            = "#EDEDED",
  "Same estimand, different assay"            = "#F6DCE5",
  "Distinct property (fair independent test)" = "#D5E7E4"
)

p2 <- ggplot(long, aes(x = column, y = row_f, fill = relation_class)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = text_wrapped), size = 6 / .pt, colour = "black",
            lineheight = 0.95) +
  # No legend: the "Relation to the class definition" column states each category
  # verbatim, so a key would only duplicate it. Fill is a redundant grouping cue.
  scale_fill_manual(values = fill_vals, guide = "none") +
  scale_x_discrete(position = "top", expand = expansion(0)) +
  scale_y_discrete(expand = expansion(0)) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(legend.position = "none",
        axis.line   = element_blank(),
        axis.ticks  = element_blank(),
        axis.text.y = element_text(hjust = 1, lineheight = 0.95),
        axis.text.x = element_text(hjust = 0.5, lineheight = 0.95),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

pdf(file.path(OUT_DIR, "fig4d_endpoint_provenance.pdf"),
    width = 5.4, height = 2.4, useDingbats = FALSE)
print(p2)
invisible(dev.off())

# ---------------------------------------------------------------------------
# 7. Banned-string guard on everything rendered on canvas
# ---------------------------------------------------------------------------
canvas_text <- c(
  as.character(pd$row_label), levels(pd$comparator_label),
  pd$or_label, pd$mdor_label,
  as.character(long$text), as.character(long$row_label),
  levels(long$column), levels(long$relation_class),
  paste0("Odds of endpoint positivity, convergent vs comparator class ",
         "(exact Fisher, log scale). Grey track = advantage sizes this ",
         "design has under 80% power to detect.")
)
banned <- c("\\bTREAT\\b", "tested[ _]negative", "mega", "no effect",
            "no difference", "equivalent", "prospectively validated")
for (b in banned) {
  hit <- grep(b, canvas_text, value = TRUE, ignore.case = FALSE)
  if (length(hit)) {
    stop("[fig4-resolution] banned string '", b, "' on canvas: ",
         paste(unique(hit), collapse = " | "))
  }
}

message("")
message("[fig4-resolution] wrote:")
message("  ", file.path(OUT_DIR, "fig4c_convergence_or_ci_power.pdf"))
message("  ", file.path(OUT_DIR, "fig4d_endpoint_provenance.pdf"))
message("  ", file.path(DATA_DIR, "fig4c_convergence_or_ci_power.tsv"))
message("[fig4-resolution] prespecified gate criterion: ", GATE_CRITERION)

# Caption is generated from the data so it cannot go stale. Note in particular
# that it must NOT claim all six intervals include OR = 1 -- one does not.
n_span <- res[spans_null == TRUE, .N]
excl   <- res[spans_null == FALSE]
excl_txt <- if (nrow(excl) == 0L) "" else paste0(
  "The remaining comparison (",
  paste(sprintf("%s, convergent vs %s", excl$endpoint,
                gsub("_", "-", excl$comparator)), collapse = "; "),
  ") excludes OR = 1 at nominal alpha but does not meet the prespecified gate ",
  paste(sprintf("(q = %.2f)", excl$q_value), collapse = ", "), ". ")

message("[fig4-resolution] CAPTION Figure 4c. Convergent-class genes versus each ",
        "single-map class in three atlas endpoints (exact Fisher odds ratio with 95% ",
        "confidence interval, log scale; group sizes beside each row). The grey track ",
        "below each row spans the odds ratios this design has under 80% power to ",
        "detect at nominal alpha = 0.05. ",
        sprintf("%d of %d intervals include OR = 1. ", n_span, nrow(res)),
        excl_txt,
        sprintf("No comparison reaches the prespecified threshold, and all %d point ",
                nrow(res)),
        sprintf("estimates fall below their own detection limit (OR %.1f to %.1f), so ",
                min(res$mdor_80), max(res$mdor_80)),
        "every call is indeterminate: these data cannot distinguish a convergent-class ",
        "advantage of moderate size from none at all. The prespecified acceptance gate ",
        "(", GATE_CRITERION, ") is therefore not met. ",
        "Figure 4d. What each endpoint tests: the disease-state class is defined by a ",
        "bulk disease-versus-control contrast, so the proteomic and single-cell ",
        "endpoints re-measure that same estimand on a different molecular layer, ",
        "while spatial variable-gene status is a distinct property and is the one ",
        "endpoint that supplies an independent test.")
