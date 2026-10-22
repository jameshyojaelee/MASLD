# Dump the 774-unit sample metadata used by the adopted fit to TSV (read-only on inputs)
m <- readRDS("RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds")
write.table(as.data.frame(m), "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/meta_774.tsv", sep="\t", quote=FALSE, row.names=FALSE)
for (d in unique(m$dataset)) { cat(d, "\n"); print(table(fib=m$fibrosis_stage[m$dataset==d], useNA="always")); print(table(nas=m$nas_score[m$dataset==d], useNA="always")); print(table(sex=m$sex[m$dataset==d], inf=m$inferred_sex[m$dataset==d], useNA="always")); print(table(m$condition[m$dataset==d])) }
