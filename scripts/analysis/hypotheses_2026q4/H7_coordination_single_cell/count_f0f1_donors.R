#!/usr/bin/env Rscript
# H7 stop-rule check: count biological donors with documented fibrosis F0 or F1.
# Usage: Rscript count_f0f1_donors.R <OUT_DIR>
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1]]
if (!file.exists(file.path(out_dir, "PRESPEC_H7.md"))) stop("PRESPEC_H7.md missing")

base <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
fmap <- file.path(base, "Analysis/Multimodal_Program_Projection/candidates",
                  "program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot",
                  "documented_fstage_donor_map.tsv")
m <- fread(fmap)
if (anyDuplicated(m$donor)) stop("donor map is not one row per donor")

by_stage <- m[, .(n_donors = .N), by = .(dataset, F_stage_documented)][order(F_stage_documented)]
n_f0f1 <- m[F_stage_documented %in% c(0L, 1L), .N]
threshold <- 15L
decision <- if (n_f0f1 < threshold) "STOP_fewer_than_15_F0_F1_donors" else "PROCEED"

fwrite(by_stage, file.path(out_dir, "documented_fstage_donor_counts.tsv"), sep = "\t")
fwrite(m[F_stage_documented %in% c(0L, 1L), .(donor, dataset, F_stage_documented)],
       file.path(out_dir, "f0f1_donor_roster.tsv"), sep = "\t")
fwrite(data.table(n_documented_donors = nrow(m), n_f0 = m[F_stage_documented == 0L, .N],
                  n_f1 = m[F_stage_documented == 1L, .N], n_f0f1 = n_f0f1,
                  threshold = threshold, decision = decision, source = fmap),
       file.path(out_dir, "stop_rule.tsv"), sep = "\t")
message(decision, ": ", n_f0f1, " documented F0-F1 donors")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
