suppressPackageStartupMessages({library(data.table); library(metafor)})
D <- "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis"
duc <- "hotspot_hepatocytes_48f39dd4d817a10e"; sec <- "hotspot_hepatocytes_f05c535ae5bbc0b9"
for (arm in c("new774_three_cohorts","new774_plus_GSE174478")) {
  ce <- fread(file.path(D, arm, "cohort_effects.tsv"))
  for (fid in c(duc, sec)) {
    d <- ce[feature_id==fid & axis=="fibrosis" & model_kind=="primary" & estimable==TRUE]
    f <- rma(yi=d$beta, sei=d$se_hc3, method="REML", test="knha")
    cat(arm, ifelse(fid==duc,"Ductular","StromalECM"), "k=", f$k, "df=", f$dfs, "beta=", round(f$b,4), "p=", signif(f$pval,4), "\n")
    print(d[, .(cohort, beta=round(beta,3), se_hc3=round(se_hc3,3), n)], row.names=FALSE)
  }
}
cat("t crit df2", qt(0.975,2), "df3", qt(0.975,3), "\n")
