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
COL_UP    <- "#C9265E"        # Liang deep magenta: resolved above zero / achieved
COL_DOWN  <- "#1565C0"        # Deep blue: resolved below zero, the model loses
COL_NULL  <- "#9E9E9E"        # Liang neutral gray: control, unresolved, not resolvable
COL_OPEN  <- "#F4A674"        # Liang warm peach: the indeterminate label

EXEC_DIR <- Sys.getenv("FIG7_EXEC_DIR")
OUT_DIR  <- Sys.getenv("FIG7_OUT_DIR")
stopifnot(nzchar(EXEC_DIR), nzchar(OUT_DIR))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

tab <- function(name) {
  p <- file.path(EXEC_DIR, "tables", name)
  read.delim(p, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE,
             colClasses = "character")
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
