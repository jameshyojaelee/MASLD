# Independent recompute of SSX-C stage coefficients from per-library axis scores.
suppressPackageStartupMessages(library(data.table))
wd <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z"
out <- file.path(wd, "verify_ssx-20261102T153423Z")
sc <- fread(file.path(wd, "out_c-20261102T134043Z/steatosite_scores_shift.tsv"))
units <- fread(file.path(wd, "../SS-20261025T135237Z/out_s2_qc/units.tsv"))[include == TRUE]
st <- fread(cmd = "tr -d '\\r' < /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv")
st[, stage := as.integer(nash_crn_kleiner_fibrosis_stage)]
d <- merge(sc, units[, .(sample_alias, inferred_sex, batch_u = batch)], by = "sample_alias")
cat("scores rows:", nrow(sc), "; with included units:", nrow(d), "; batch mismatch:", sum(d$batch != d$batch_u), "\n")
d <- merge(d, st[, .(sample_alias, stage)], by = "sample_alias")
d <- d[stage %in% 1:4]
cat("staged units:", nrow(d), "; stage table:", paste(names(table(d$stage)), table(d$stage), collapse = " "), "\n")
d[, bin := factor(c("F0F1", "F2", "F3F4", "F3F4")[stage], levels = c("F0F1", "F2", "F3F4"))]
d[, stage4 := factor(c("F0F1", "F2", "F3", "F4")[stage], levels = c("F0F1", "F2", "F3", "F4"))]
d[, batch := factor(batch)]; d[, inferred_sex := factor(inferred_sex)]
res <- rbindlist(lapply(paste0("P", 1:5), function(a) {
  d[, z := (get(a) - mean(get(a))) / sd(get(a))]
  cb <- coef(lm(z ~ batch + inferred_sex + bin, data = d))
  c4 <- coef(lm(z ~ batch + inferred_sex + stage4, data = d))
  data.table(axis = a, binF3F4 = cb[["binF3F4"]], stage4F2 = c4[["stage4F2"]],
             stage4F3 = c4[["stage4F3"]], stage4F4 = c4[["stage4F4"]])
}))
fwrite(res, file.path(out, "c_coefficients.tsv"), sep = "\t")
print(res, digits = 6)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_c.txt"))
