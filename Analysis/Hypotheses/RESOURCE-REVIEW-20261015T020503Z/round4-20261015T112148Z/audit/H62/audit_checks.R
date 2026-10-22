# H62 audit: n of the adult refit units, package versions, power arithmetic. No adolescent expression value read.
m <- readRDS("RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds")
print(table(m$group_binary)); print(table(m$dataset, m$group_binary))
cat("R", R.version.string, "limma", as.character(packageVersion("limma")), "edgeR", as.character(packageVersion("edgeR")), "\n")
f <- function(d, df, a) { ncp <- d/sqrt(1/34+1/33); q <- qt(1-a/2, df); 1-pt(q,df,ncp)+pt(-q,df,ncp) }
cat("MDE80 alpha 0.025 df62", uniroot(function(d) f(d,62,0.025)-.8, c(.1,2))$root, "\n")
cat("normal approx MDE80", (qnorm(.975)+qnorm(.8))*sqrt(1/34+1/33), "\n")
cat("Spearman/Pearson n=52 r for 80% power (Fisher z)", tanh((qnorm(.975)+qnorm(.8))/sqrt(52-3)), "\n")
cat("MDE80 df62 for 33 vs 34 after S4 removal of 3 (e.g. 32 vs 32)", uniroot(function(d){ncp<-d/sqrt(1/32+1/32);q<-qt(.975,59);1-pt(q,59,ncp)+pt(-q,59,ncp)-.8},c(.1,2))$root, "\n")
