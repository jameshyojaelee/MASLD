# Shared setup for the rebuilt Figure 7 panels.
#
# The single rule this rebuild enforces: a panel shows data. Every word that explains,
# qualifies or attributes is emitted with message() as caption text and is never drawn
# on the canvas. Y-axis labels carry the short name only.
#
# Standing figure rules enforced here, not per script:
#   font size 6 everywhere, weight normal; all text black; no lollipops; no 3D;
#   control colour #9E9E9E; PDF only with useDingbats = FALSE; one panel per file.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

source(file.path(Sys.getenv("HOME"), "publication_color_themes.R"))

FS <- 6                       # the only font size in the figure
GS <- FS / .pt                # the same 6 pt expressed for ggplot geom text
# Candidate 2026-10-17 (panel-deep-fixes): one colour code across 7B-7F.
#   interval excludes the baseline above: magenta filled circle
#   interval excludes the baseline below: blue filled down-triangle
#   interval includes the baseline:       gray open circle
#   a bound measured from the data:       black open diamond
#   chance / a pre-declared reference:    gray filled square
# The 7E verdicts reuse it (supported / tested negative / untestable). The 7F exposure grid
# uses teal / violet / light gray, which no other mark in the figure uses.
COL_UP    <- "#C9265E"
COL_DOWN  <- "#1565C0"
COL_NULL  <- "#9E9E9E"

OUT_DIR  <- Sys.getenv("FIG7_OUT_DIR")
TAB_DIR  <- Sys.getenv("FIG7_TABLES")
stopifnot(nzchar(OUT_DIR), nzchar(TAB_DIR))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(20261017)  # nothing here is random; set for the record

# Every drawn value is read from the panel's source table; nothing is computed in R.
src_tab <- function(panel) {
  d <- read.delim(file.path(TAB_DIR, paste0("fig", panel, ".tsv")), sep = "\t", quote = "",
                  check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  d$colour <- toupper(d$colour)
  d
}

SHAPE <- c(circle_filled = 21, circle_open = 21, triangle_down_filled = 25,
           diamond_open = 23, square_filled = 22)
glyph_fill <- function(glyph, colour) ifelse(glyph %in% c("circle_open", "diamond_open"), "white", colour)

# Draw one block of marks on its own x axis. `key` names the legend entry of each mark.
forest_block <- function(d, xlab, vlines = NULL, legend = TRUE, nudge = NULL, xlim = NULL,
                         breaks = waiver(), keys = NULL) {
  d$y <- factor(d$row_label, levels = rev(unique(d$row_label[order(d$row_order)])))
  d$yn <- as.numeric(d$y) + if (is.null(nudge)) 0 else nudge
  # draw reference squares and bound diamonds first, so an achieved mark is never hidden
  d <- d[order(match(d$glyph, c("square_filled", "diamond_open", "circle_open",
                                "circle_filled", "triangle_down_filled"))), ]
  if (is.null(keys)) keys <- unique(d[, c("key", "glyph", "colour")])
  stopifnot(!anyDuplicated(keys$key), all(d$key %in% keys$key))
  p <- ggplot(d, aes(x = x, y = yn))
  if (!is.null(vlines)) for (v in vlines)
    p <- p + geom_vline(xintercept = v$x, linetype = v$lty, colour = v$col, linewidth = 0.3)
  di <- d[!is.na(d$lo), ]
  if (nrow(di)) p <- p + geom_segment(data = di, aes(x = lo, xend = hi, y = yn, yend = yn, colour = key),
                                      linewidth = 0.4, show.legend = FALSE)
  p <- p + geom_point(aes(shape = key, fill = key, colour = key), size = 1.3, stroke = 0.45,
                      show.legend = TRUE) +
    scale_shape_manual(values = setNames(SHAPE[keys$glyph], keys$key), breaks = keys$key, limits = keys$key, drop = FALSE, name = NULL) +
    scale_fill_manual(values = setNames(glyph_fill(keys$glyph, keys$colour), keys$key), breaks = keys$key, limits = keys$key, drop = FALSE, name = NULL) +
    scale_colour_manual(values = setNames(keys$colour, keys$key), breaks = keys$key, limits = keys$key, drop = FALSE, name = NULL) +
    scale_y_continuous(breaks = seq_along(levels(d$y)), labels = levels(d$y),
                       expand = expansion(add = 0.6)) +
    scale_x_continuous(breaks = breaks, limits = xlim) +
    labs(x = xlab, y = NULL) +
    guides(shape = guide_legend(nrow = if (nrow(keys) > 3) 2 else 1, byrow = TRUE)) +
    theme_fig7() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = FS, colour = "black", face = "plain", hjust = 1,
                                      margin = margin(t = 2)),
          legend.position = if (legend) "bottom" else "none",
          legend.justification = "right",
          plot.margin = margin(2, 4, 2, 2))
  # A block title is drawn as a facet strip (data label naming the block), never as a caption.
  # Labels are right-aligned to the panel edge and may extend left over the row-label column.
  if ("strip_lab" %in% names(d)) p <- p + facet_grid(cols = vars(strip_lab)) +
      theme(strip.text.x = element_text(size = FS, colour = "black", face = "plain", hjust = 1),
            strip.clip = "off", panel.spacing.x = unit(10, "pt"))
  p
}

num <- function(x) suppressWarnings(as.numeric(x))

# Short-name lookup. A name that is not in the map raises, so a long method sentence
# can never leak back onto an axis by falling through unmatched.
short <- function(x, map) {
  hit <- match(x, names(map))
  if (anyNA(hit)) stop("unmapped label(s): ", paste(unique(x[is.na(hit)]), collapse = " | "))
  unname(map[hit])
}

# Caption text goes to the log, never to the canvas.
caption <- function(...) {
  message("\n[CAPTION] ", ...)
}

# Base theme: 6 pt, no bold, black text, no grid.
theme_fig7 <- function() {
  theme_classic(base_size = FS, base_family = "Helvetica") %+replace%
    theme(
      text              = element_text(size = FS, colour = "black", face = "plain"),
      plot.title        = element_blank(),
      plot.subtitle     = element_blank(),
      plot.caption      = element_blank(),
      axis.title        = element_text(size = FS, colour = "black", face = "plain"),
      axis.title.y      = element_text(size = FS, colour = "black", face = "plain",
                                       angle = 90, margin = margin(r = 2)),
      axis.title.x      = element_text(size = FS, colour = "black", face = "plain",
                                       margin = margin(t = 2)),
      axis.text         = element_text(size = FS, colour = "black", face = "plain"),
      axis.text.y       = element_text(size = FS, colour = "black", face = "plain", hjust = 1),
      legend.title      = element_blank(),
      legend.text       = element_text(size = FS, colour = "black", face = "plain"),
      legend.key.size   = unit(0.20, "cm"),
      legend.position   = "bottom",
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.margin = margin(-4, 0, 0, 0),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      strip.text        = element_text(size = FS, colour = "black", face = "plain", hjust = 0,
                                       margin = margin(1, 1, 1, 1)),
      strip.background  = element_blank(),
      strip.placement   = "outside",
      panel.grid        = element_blank(),
      panel.spacing.y   = unit(3, "pt"),
      panel.spacing.x   = unit(4, "pt"),
      plot.margin       = margin(3, 3, 3, 3),
      axis.line         = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks        = element_line(linewidth = 0.3, colour = "black")
    )
}

save_panel <- function(plot, filename, width, height) {
  p <- file.path(OUT_DIR, filename)
  grDevices::pdf.options(useDingbats = FALSE)
  ggplot2::ggsave(p, plot, width = width, height = height, units = "in",
                  device = grDevices::pdf, useDingbats = FALSE)
  message("wrote ", p)
  invisible(p)
}
