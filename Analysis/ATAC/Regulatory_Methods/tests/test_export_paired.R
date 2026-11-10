#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Pass path to export_paired.R")
source(args[[1L]])
fixture <- data.frame(donor_id = c(rep("a", 100), rep("b", 99), rep("c", 1200)), cell_id = seq_len(1399))
selected <- select_cells(fixture, "pilot", seed = 42L)
stopifnot(length(selected) == 200L,
          identical(sort(unique(fixture$donor_id[selected])), c("a", "c")),
          all(table(fixture$donor_id[selected]) == 100L),
          identical(selected, select_cells(fixture, "pilot", seed = 42L)))
full <- select_cells(fixture, "full", seed = 42L)
stopifnot(length(full) == 1100L, !anyDuplicated(full), all(table(fixture$donor_id[full]) >= 100L))
message("paired donor eligibility, cap, and deterministic selection: PASS")
