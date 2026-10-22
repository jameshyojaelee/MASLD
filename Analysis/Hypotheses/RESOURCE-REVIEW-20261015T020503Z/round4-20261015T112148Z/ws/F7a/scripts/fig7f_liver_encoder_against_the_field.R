# KEY MESSAGE: the two encoders with the highest pooled score have ZERO clean held-out
# studies, so their lead is not evidence of skill; and the advantage every encoder shows
# over a 50-dimensional PCA of raw counts collapses once the head has capacity.

source(file.path(Sys.getenv("FIG7_SCRIPTS"), "_fig7_common.R"))
suppressPackageStartupMessages(library(patchwork))

ROOT <- Sys.getenv("MASLD_ROOT"); stopifnot(nzchar(ROOT))
MM_REL <- "Analysis/MASLD_Model_Benchmark/executions/encoder-field-score-21314196/evaluation/model_metrics.tsv"
mm <- read.delim(file.path(ROOT, MM_REL), sep = "\t", check.names = FALSE,
                 quote = "", stringsAsFactors = FALSE)
mm$donor_class_balanced_macro_f1 <- num(mm$donor_class_balanced_macro_f1)

g <- tab("panel_7f_gain.tsv")
e <- tab("panel_7f_exposure_grid.tsv")
g$gain <- num(g$gain); g$lo <- num(g$lo); g$hi <- num(g$hi)
e$donors <- num(e$donors)

SHORTM <- c("Geneformer v2 316M"               = "Geneformer",
            "TranscriptFormer TF-Sapiens"      = "TranscriptFormer",
            "SCimilarity v1.1"                 = "SCimilarity",
            "UCE 4-layer"                      = "UCE-4L",
            "UCE 33-layer"                     = "UCE-33L",
            "scVI liver encoder, trained here" = "scVI liver",
            "HVG + PCA of raw counts"          = "HVG+PCA")

# Row order is the pooled donor-class-balanced macro F1 of each block carrying the
# two-layer MLP head, read from the producing file: best at the top.
pool <- do.call(rbind, lapply(unique(c(g$block_id, e$block_id)), function(bid) {
  r <- mm[mm$model_id == paste0(bid, "::two_layer_mlp"), ]
  stopifnot(nrow(r) == 1L)
  data.frame(block_id = bid, pooled = r$donor_class_balanced_macro_f1,
             clean = if (any(e$block_id == bid)) sum(e$block_id == bid & e$head_id == "linear" & e$disposition == "clean") else as.integer(r$clean_folds), stringsAsFactors = FALSE)
}))
pool <- pool[order(-pool$pooled), ]
MORDER <- pool$block_id

# ------------------------------------------------------ block a: the exposure grid
el <- e[e$head_label == "linear head", ]
em <- e[e$head_label == "two-layer MLP head", ]
chk <- merge(el[, c("block_id", "held_out_study", "disposition")],
             em[, c("block_id", "held_out_study", "disposition")],
             by = c("block_id", "held_out_study"))
stopifnot(nrow(chk) == nrow(el), all(chk$disposition.x == chk$disposition.y))

DISP <- c("clean" = "clean", "confounded_by_exposure" = "confounded",
          "exposure_unresolved" = "unresolved")
el$disp <- factor(short(el$disposition, DISP), levels = unname(DISP))
# Order model rows by the pooled ranking, keeping only the models present in this block.
lev_all <- vapply(MORDER, function(bid) {
  lab <- unique(c(e$model_label[e$block_id == bid], g$model_label[g$block_id == bid]))
  stopifnot(length(lab) == 1L)
  short(lab, SHORTM)
}, character(1))
el$mshort <- factor(short(el$model_label, SHORTM),
                    levels = rev(lev_all[lev_all %in% short(el$model_label, SHORTM)]))
stopifnot(!any(is.na(el$mshort)))

sord <- unique(el[, c("held_out_study", "donors")])
sord <- sord[order(-sord$donors), ]
el$study <- factor(el$held_out_study, levels = sord$held_out_study)

FILLD <- c("clean" = "#1565C0", "confounded" = COL_UP, "unresolved" = COL_NULL)

pa <- ggplot(el, aes(x = study, y = mshort, fill = disp)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = FILLD, drop = FALSE) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_fig7() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(size = FS, colour = "black", face = "plain",
                                       angle = 45, hjust = 0, vjust = 0))

# ---------------------------------------- blocks b and c: the head-capacity collapse
BB <- unique(g$block[startsWith(g$block, "b.")]); stopifnot(length(BB) == 1L)
gb <- g[g$block == BB, ]
gb$mshort <- factor(short(gb$model_label, SHORTM),
                    levels = rev(lev_all[lev_all %in% short(gb$model_label, SHORTM)]))
stopifnot(!any(is.na(gb$mshort)))
gb$headf <- factor(gb$head_label, levels = c("linear head", "two-layer MLP head"))
stopifnot(!any(is.na(gb$headf)))
# Review R08-05 (2026-10-15): an encoder with no clean held-out study is drawn as
# uninterpretable (hollow gray), whatever its interval, as Legend A states.
NOCLEAN <- startsWith(gb$eligibility, "uninterpretable")
stopifnot(all(as.integer(gb$clean_studies[NOCLEAN]) == 0L), all(as.integer(gb$clean_studies[!NOCLEAN]) > 0L))
gb$cls <- ifelse(NOCLEAN, "no clean study",
                 ifelse(gb$lo > 0, "above zero", ifelse(gb$hi < 0, "below zero", "crosses zero")))
gb$cls <- factor(gb$cls, levels = c("above zero", "below zero", "crosses zero", "no clean study"))
COLS <- c("above zero" = COL_UP, "below zero" = COL_DOWN, "crosses zero" = COL_NULL,
          "no clean study" = COL_NULL)
SHPF <- c("above zero" = 16, "below zero" = 16, "crosses zero" = 16, "no clean study" = 1)

pb <- ggplot(gb, aes(x = gain, y = mshort, colour = cls, shape = cls)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "black") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4) +
  geom_point(size = 1.1, stroke = 0.4) +
  scale_colour_manual(values = COLS, drop = TRUE) +
  scale_shape_manual(values = SHPF, drop = TRUE) +
  scale_x_continuous(breaks = c(0, 0.02, 0.04, 0.06)) +
  facet_grid(cols = vars(headf)) +
  labs(x = "gain over HVG + PCA", y = NULL) +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme_fig7() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

p <- pa / pb + plot_layout(heights = c(1, 1.15))

save_panel(p, "fig7F_liver_encoder_against_the_field.pdf", width = 4.5, height = 3.1)

nclean <- setNames(pool$clean, pool$block_id)
top2 <- pool$block_id[1:2]
lin <- gb[gb$headf == "linear head" & gb$block_id != "scvi_liver_latent" & !NOCLEAN, ]
mlp <- gb[gb$headf == "two-layer MLP head" & gb$block_id != "scvi_liver_latent" & !NOCLEAN, ]

caption("Fig. 7F. The liver-specific encoder against the field under identical exposure ",
        "rules. Substrate: 102 donors (the unit of inference), 50,000 identical cells, 5 ",
        "classes, 7 studies held out by study; 10,000-resample donor cluster bootstrap with ",
        "shared multiplicities across all models, so deltas are paired on the same draws. ",
        "Source: Analysis/MASLD_Model_Benchmark/executions/encoder-field-score-21314196/evaluation/.")
caption("Top: the training exposure of every held-out study for every encoder. Rows are ",
        "ordered by pooled donor-class-balanced macro F1 (two-layer MLP head), best at the ",
        "top, read from ", MM_REL, ": ",
        paste(sprintf("%s %.4f (%d clean folds)", lev_all, pool$pooled, pool$clean),
              collapse = "; "),
        ". The two best-pooled encoders, ", lev_all[1], " and ", lev_all[2],
        ", have respectively ", pool$clean[1], " and ", pool$clean[2], " source-qualified clean studies. All-study pooled scores and gains retain encoder-seen studies and are not clean-only estimates. ",
        "Column order is by donors held out: ",
        paste(sprintf("%s %d donors", sord$held_out_study, sord$donors), collapse = ", "), ".")
caption("Bottom: paired gain in donor-class-balanced macro F1 over HVG + 50-d PCA of raw ",
        "counts carrying the identical head; point = estimate, whisker = 95% donor cluster ",
        "bootstrap interval; magenta the interval excludes zero, gray it crosses zero; hollow ",
        "gray = no clean held-out study, so the gain is uninterpretable whatever its interval ",
        "(Geneformer: 0 clean, 3 confounded, 4 unresolved). With a linear head the pretrained ",
        "encoders with at least one clean study gain ",
        sprintf("%+.4f to %+.4f", min(lin$gain), max(lin$gain)),
        " and the liver encoder ",
        sprintf("%+.4f", gb$gain[gb$block_id == "scvi_liver_latent" & gb$headf == "linear head"]),
        "; with a two-layer MLP head the same encoders gain only ",
        sprintf("%+.4f to %+.4f", min(mlp$gain), max(mlp$gain)),
        " and the liver encoder ",
        sprintf("%+.4f", gb$gain[gb$block_id == "scvi_liver_latent" & gb$headf == "two-layer MLP head"]),
        ", whose interval then crosses zero. Row values: ",
        paste(sprintf("%s, %s: %+.4f [%+.4f, %+.4f] (clean %s / confounded %s / unresolved %s of 7)",
                      short(gb$model_label, SHORTM), gb$head_label, gb$gain, gb$lo, gb$hi,
                      gb$clean_studies, gb$confounded_studies, gb$unresolved_studies),
              collapse = "; "), ".")
