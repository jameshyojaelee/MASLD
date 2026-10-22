# KEY MESSAGE: a question receives one of four words, and which word it receives is decided
# by where the observed statistic falls relative to a threshold fixed before the answer.

source(file.path(Sys.getenv("FIG7_SCRIPTS"), "_fig7_common.R"))
source(file.path(Sys.getenv("FIG7_SCRIPTS"), "_fig7e_guards.R"))

v <- tab("panel_7e_verdicts.tsv")
v$observed <- num(v$observed); v$threshold <- num(v$threshold)
stopifnot(!any(is.na(v$observed)), !any(is.na(v$threshold)), all(v$threshold > 0))

# The class plotted for each row must be the REALISED verdict, never a conditional
# "if it had been negative" field. Both guards run before anything is drawn.
VOCAB <- verdict_vocabulary(tab("panel_7a_framework.tsv"))
check_verdict_rows(v, VOCAB)
VERD <- setNames(gsub("_", " ", VOCAB, fixed = TRUE), VOCAB)

QUES <- c(
  "RNA->chromatin map differs between normal and diseased liver"  = "normal vs diseased map",
  "liver zone recoverable from bulk expression (GSE105127)"       = "zone from expression",
  "graded severity WITHIN diseased liver (normal livers dropped)" = "severity within disease",
  "graded severity within diseased liver, female-only"            = "severity, female only",
  "severity dependence of the RNA->ATAC map (fibrosis any)"       = "ATAC map vs fibrosis",
  "severity dependence of the RNA->ATAC map (steatosis >= 5%)"    = "ATAC map vs steatosis",
  "common-variant heritability of the 117 programs: FinnGen R12 NAFLD" = "NAFLD endpoint",
  "common-variant heritability of the 117 programs: Ghodsian 2021 NAFLD meta-analysis" = "NAFLD meta-analysis",
  "common-variant heritability of the 117 programs: FinnGen R12 cirrhosis" = "cirrhosis endpoint",
  "signed cis-eQTL direction from DNA-LM embeddings"              = "signed eQTL direction",
  "117 programs carry liver-enzyme heritability beyond a lineage-matched background (ukbb_alt)" = "program h2, ALT",
  "117 programs carry liver-enzyme heritability beyond a lineage-matched background (ukbb_ast)" = "program h2, AST",
  "severity scores move correctly when the same person is treated" = "treatment response",
  "a routing rule over model components adds skill"               = "routing adds skill")

v$verdf <- droplevels(factor(short(v$verdict, VERD), levels = unname(VERD)))
v$qshort <- short(v$question, QUES)
v$ratio <- v$observed / v$threshold

v <- v[order(v$verdf, v$ratio), ]
v$key <- paste0(as.integer(v$verdf), "|", v$qshort)
stopifnot(!any(duplicated(v$key)))
v$key <- factor(v$key, levels = v$key)
LABS <- setNames(v$qshort, as.character(v$key))

COLS <- c("supported" = COL_UP, "indeterminate" = COL_OPEN,
          "untestable" = COL_NULL, "tested negative" = COL_DOWN)

p <- ggplot(v, aes(x = ratio, y = key, colour = verdf)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = COL_NULL) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = COLS, guide = "none") +
  scale_y_discrete(labels = LABS) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5), limits = c(-0.55, 1.95), expand = c(0, 0)) +
  facet_grid(rows = vars(verdf), scales = "free_y", space = "free_y") +
  labs(x = "observed / its own decision threshold", y = NULL) +
  theme_fig7() +
  theme(strip.text.y = element_text(size = FS, colour = "black", face = "plain",
                                    angle = 0, hjust = 0, margin = margin(1, 1, 1, 5)),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

save_panel(p, "fig7E_verdict_vocabulary_and_refusal_thresholds.pdf", width = 3.3, height = 2.5)

caption("Fig. 7E. Every question in Figure 7 plotted as its observed statistic divided by ",
        "the threshold fixed for it before the answer was seen; the dashed line is 1. A ",
        "question to the right of the line cleared its own bar, one to the left did not. Each ",
        "statistic is oriented so that positive is the prespecified direction; the gray line is 0. ",
        "The four words are the only verdicts a question can receive. Rows:")
for (vv in levels(v$verdf)) {
  d <- v[v$verdf == vv, ]
  for (i in seq_len(nrow(d))) {
    caption("  [", vv, "] ", d$qshort[i], " - ", d$question[i], "; substrate ", d$unit_and_n[i],
            "; statistic ", d$statistic[i], "; observed ", signif(d$observed[i], 5),
            " vs threshold ", signif(d$threshold[i], 5), " (", d$threshold_name[i],
            "), ratio ", sprintf("%.3f", d$ratio[i]), ". ", d$decision_rule[i],
            " Source ", d$source[i], ".")
  }
}
caption("The vocabulary has four words; ", paste(setdiff(unname(VERD), levels(v$verdf)),
        collapse = ", "), " was available but no question in this figure received it, so it ",
        "is not drawn. The zonation question is classed supported from the REALISED verdict ",
        "field of its producing file (zone_metrics.json::verdict = 'positive'). That file ",
        "also records a conditional negative_result_label = 'indeterminate', the label that ",
        "would have applied had the result come out negative; it never applied, because the ",
        "accuracy cleared both the pre-registered 0.05 bar and the MDE. An earlier build of ",
        "this panel plotted that conditional field as though it were the outcome; guards G1 ",
        "and G2 in _fig7e_guards.R now reject that table (see ",
        "scripts/99_demo_verdict_guard_fires.R).")
