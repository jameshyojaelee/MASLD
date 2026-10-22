# Diagnostic after K3 stop: in which labelled AOI does each marker set peak, per ROI. Not an outcome.
suppressPackageStartupMessages(library(data.table))
RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V6-20261021T234432Z"
d <- readRDS(file.path(RUN, "derived/geomx_rebuilt.rds")); a <- as.data.table(d$aoi)
Z <- t(scale(t(d$L)))
sets <- list(Albumin = c("ALB", "APOA1"), Vinculin = c("COL1A1", "COL1A2", "COL3A1"), CD45 = "PTPRC")
for (s in names(sets)) a[[paste0("z_", s)]] <- colMeans(Z[sets[[s]], a$dcc, drop = FALSE])
a[, peak := c("Albumin", "Vinculin", "CD45")[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = c("z_Albumin", "z_Vinculin", "z_CD45")]
print(a[, .N, by = .(labelled = stain, marker_peak = peak)][order(labelled, -N)])
r <- rbindlist(lapply(names(sets), function(s) a[, .(set = s, peak_in = stain[which.max(get(paste0("z_", s)))]), by = roi]))
print(dcast(r[, .N, by = .(set, peak_in)], set ~ peak_in, value.var = "N", fill = 0))
print(a[, .(slide, roi, stain, well = sub(".*-", "", dcc), q3 = round(q3), aligned, z_Albumin = round(z_Albumin, 2), z_Vinculin = round(z_Vinculin, 2), z_CD45 = round(z_CD45, 2))][order(roi, stain)][1:24])
fwrite(a, file.path(RUN, "outputs/diag_K3_marker_peaks.tsv"), sep = "\t")
