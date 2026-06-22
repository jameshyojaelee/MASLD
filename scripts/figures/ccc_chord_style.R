#!/usr/bin/env Rscript
# ============================================================================
# ccc_chord_style.R
#
# Shared polished-chord renderer for the CCC chord diagrams (the stage-gated
# headline chord in ccc_v3_panels.R and the full-network chord in
# ccc_full_network_chord.R). Centralises the "fancy" aesthetic so both panels
# stay visually consistent.
#
# pretty_chord(adj, nodes, grid_col, ...)
#   adj      : data.table/data.frame with columns from, to, value, col
#   nodes    : ordered character vector of sector ids (short labels)
#   grid_col : named vector sector-id -> sector block colour
# ============================================================================

suppressPackageStartupMessages({
  library(circlize)
  library(data.table)
})

# soft tints of the sector colours for the outer halo ring
.lighten <- function(hex, amt = 0.55) {
  rgb <- grDevices::col2rgb(hex) / 255
  rgb <- rgb + (1 - rgb) * amt
  grDevices::rgb(rgb[1], rgb[2], rgb[3])
}

pretty_chord <- function(adj, nodes, grid_col,
                         start.degree = 90,
                         gap.degree   = 4,
                         label_cex    = 1.05,
                         show_axis    = FALSE,
                         show_total   = TRUE,
                         axis_step    = NULL,
                         link_border  = "white",
                         link_lwd     = 0.5) {

  adj <- as.data.table(adj)
  nodes <- unname(nodes[nodes %in% c(adj$from, adj$to)])

  # per-sector total = sum of all ribbon widths touching the sector
  # (condensed alternative to a full tick ruler). unname() above keeps the
  # sapply names equal to the short sector ids (= sector.index in panel.fun).
  sector_total <- sapply(nodes, function(s)
    sum(adj[from == s | to == s]$value))

  circos.clear()
  circos.par(start.degree            = start.degree,
             gap.degree              = gap.degree,
             track.margin            = c(0.005, 0.005),
             cell.padding            = c(0, 0, 0, 0),
             points.overflow.warning = FALSE,
             canvas.xlim             = c(-1.45, 1.45),
             canvas.ylim             = c(-1.45, 1.45))

  chordDiagram(
    x                  = adj[, .(from, to, value)],
    order              = nodes,
    grid.col           = grid_col,
    col                = adj$col,
    transparency       = NA,            # alpha already baked into adj$col
    directional        = 1,
    direction.type     = c("diffHeight", "arrows"),
    diffHeight         = mm_h(1.6),
    link.arr.type      = "big.arrow",
    link.lwd           = link_lwd,
    link.lty           = 1,
    link.border        = link_border,
    grid.border        = "white",
    link.sort          = TRUE,
    link.largest.ontop = TRUE,
    annotationTrack    = "grid",
    annotationTrackHeight = mm_h(3.5),
    preAllocateTracks  = list(list(track.height = mm_h(9.5)))
  )

  # outer halo: a faint wide ring echoing each sector colour (depth cue)
  circos.track(track.index = 1, ylim = c(0, 1), bg.border = NA,
               track.height = mm_h(1.4), panel.fun = function(x, y) {
    s  <- get.cell.meta.data("sector.index")
    xl <- get.cell.meta.data("xlim")
    circos.rect(xl[1], 0, xl[2], 1, col = .lighten(grid_col[s], 0.6),
                border = NA)
  })

  # axis ticks + sector labels on the grid track (track index 2)
  circos.track(track.index = 2, bg.border = NA, panel.fun = function(x, y) {
    s  <- get.cell.meta.data("sector.index")
    xl <- get.cell.meta.data("xlim")
    yl <- get.cell.meta.data("ylim")
    if (show_axis) {
      span <- xl[2] - xl[1]
      step <- if (!is.null(axis_step)) axis_step else {
        raw <- span / 4
        mag <- 10^floor(log10(max(raw, 1)))
        c(1, 2, 5, 10)[which.min(abs(c(1, 2, 5, 10) * mag - raw))] * mag
      }
      brks <- seq(0, floor(xl[2] / step) * step, by = step)
      circos.axis(h = "top", major.at = brks, labels = brks,
                  labels.cex = 0.5, labels.col = "grey35",
                  col = "grey55", lwd = 0.7,
                  major.tick.length = mm_y(0.9), minor.ticks = 1)
    }
    # condensed count: one total per sector, set below the name with a clear gap
    if (show_total) {
      circos.text(mean(xl), yl[2] + mm_y(1.8),
                  format(sector_total[s], big.mark = ","),
                  facing = "bending.inside", niceFacing = TRUE,
                  cex = label_cex * 0.66, font = 1, col = "grey40")
    }
    circos.text(mean(xl), yl[2] + mm_y(7.2), s,
                facing = "bending.inside", niceFacing = TRUE,
                cex = label_cex, font = 2, col = "grey12")
  })

  circos.clear()
}
