#!/usr/bin/env Rscript
# DIR6 post hoc diagnostic (not prespecified): the PNPLA3 and TM6SF2 positive controls got lfsr 1 because
# their effects lie outside the range of the halfuniform prior fitted on the H26 instruments (fixed g).
# Here ashr is refitted on instruments plus the six control rows together; control calls and any change in
# instrument calls are reported. Gene classes in dir6_gene_classes.tsv stay the prespecified ones.
# Usage: 03_posthoc_control_refit.R <step2 out dir> <new out dir>
suppressMessages({library(data.table); library(ashr)})
a <- commandArgs(trailingOnly = TRUE); IN <- a[1]; OUT <- a[2]
set.seed(20261111)
w <- fread(file.path(IN, "instrument_effects.tsv"))
res <- list()
for (tr in c("pdff", "cirr")) {
  b <- w[[paste0("beta_al_", tr)]]; s <- w[[paste0("se_", tr)]]
  use <- is.finite(b) & is.finite(s)
  fit <- ash(b[use], s[use], mixcompdist = "halfuniform")
  x <- w[use, .(id, outcome = tr, z = b[use] / s[use], pm_refit = get_pm(fit), lfsr_refit = get_lfsr(fit),
                call_prespec = get(paste0("call_", tr)))]
  x[, call_refit := fifelse(lfsr_refit >= 0.05, "no outcome support", fifelse(pm_refit > 0, "concordant", "reversed"))]
  res[[tr]] <- x
}
r <- rbindlist(res)
fwrite(r, file.path(OUT, "posthoc_control_refit.tsv"), sep = "\t")
print(r[grepl("^ctrl_", id)])
cat("instrument rows whose call changes after refit:", r[!grepl("^ctrl_", id) & call_refit != call_prespec, .N], "\n")
print(r[!grepl("^ctrl_", id) & call_refit != call_prespec])
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
