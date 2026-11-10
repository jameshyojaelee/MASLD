suppressPackageStartupMessages(library(data.table))
task <- '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/NUCLEUS-FIX-20261107'
r <- fread(file.path(task,'fixed_programs/results.tsv'))
saved <- readRDS(file.path(task,'fixed_programs/nulls.rds'))
stopifnot(nrow(r)==117L,length(saved$null)==117L,saved$draws==19999L,saved$seed==20260923L)
for (i in seq_len(117L)) {
  if (!r$testable[i]) { stopifnot(r$p[i]==1); next }
  null <- saved$null[[i]]
  stopifnot(length(null)==19999L,all(is.finite(null)))
  p <- (1+sum(abs(null)>=abs(r$t[i])*(1-1e-12)))/20000
  stopifnot(abs(p-r$p[i])<1e-12)
}
stopifnot(max(abs(p.adjust(r$p,'BH',n=117)-r$q))<1e-12)
writeLines('PASS: saved 19,999-draw nulls reproduce all 117 permutation p values and BH q values.',file.path(task,'checks/permutation_check.txt'))
