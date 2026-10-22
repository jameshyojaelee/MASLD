# MDE and power for GSE289173 estimands; inputs are discovery donor scores only.
set.seed(20260923)
R <- "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/"
G <- "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gse289173/gate/"
s <- read.delim(paste0(R,"donor_program_scores_primary.tsv"))[,c("cell_type","module","donor","score")]
m <- read.delim("Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv")
m <- m[toupper(as.character(m$exclude_stage_analysis))!="TRUE" & m$disease_stage_coarse %in% c("Healthy","Steatosis","Steatohepatitis") & m$dataset!="GSE189600",]
d <- merge(s, m[,c("sample","dataset","disease_stage_coarse")], by.x="donor", by.y="sample")
d$stage <- match(d$disease_stage_coarse, c("Healthy","Steatosis","Steatohepatitis"))-1
d <- d[d$dataset %in% c("GSE202379","GSE244832"),]
# standardized group means per dataset (residual SD of stage-factor model), n-weighted across datasets
stdmeans <- function(h){ f <- lm(score ~ factor(stage), h); sdr <- sigma(f); mu <- tapply(h$score,h$stage,mean); (mu-mu["0"])/sdr }
out <- do.call(rbind, lapply(split(d, paste(d$cell_type,d$module)), function(g){
  r <- lapply(split(g,g$dataset), function(h) c(stdmeans(h), n=nrow(h)))
  w <- sapply(r,`[`,"n"); M <- sapply(r, function(x) x[c("1","2")])
  data.frame(cell_type=g$cell_type[1], module=g$module[1], m1=sum(M[1,]*w)/sum(w), m2=sum(M[2,]*w)/sum(w))}))
write.table(out, paste0(G,"discovery_std_stage_means.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(out[out$cell_type=="hepatocytes" & out$module %in% c(8,20),])
# analytic two-sample MDE (equal SD) via noncentral t
mde2 <- function(n1,n2,alpha,power=.8){ df<-n1+n2-2; tc<-qt(1-alpha/2,df); f<-function(dd){ ncp<-dd/sqrt(1/n1+1/n2); 1-pt(tc,df,ncp)+pt(-tc,df,ncp)-power}; uniroot(f,c(0.01,10))$root }
pow2 <- function(n1,n2,dd,alpha){ df<-n1+n2-2; tc<-qt(1-alpha/2,df); ncp<-dd/sqrt(1/n1+1/n2); 1-pt(tc,df,ncp)+pt(-tc,df,ncp)}
# ordinal trend: power for slope test given group sizes and standardized means
powtrend <- function(ns,mu,alpha){ x<-rep(0:2,ns); mux<-rep(mu,ns); b<-cov(x,mux)/var(x); res<-var(mux-b*x)*(length(x)-1)/(length(x)-2)
  s2 <- 1+res; se<-sqrt(s2/sum((x-mean(x))^2)); df<-length(x)-2; tc<-qt(1-alpha/2,df); ncp<-b/se; 1-pt(tc,df,ncp)+pt(-tc,df,ncp)}
mdetrend <- function(ns,alpha){ x<-rep(0:2,ns); df<-length(x)-2; tc<-qt(1-alpha/2,df); sx<-sqrt(sum((x-mean(x))^2)); uniroot(function(b) 1-pt(tc,df,b*sx)+pt(-tc,df,b*sx)-.8,c(.01,10))$root }
a1 <- 0.05; aF <- 0.05/117
cat("\nMDE (SD units, 80% power, two-sided)\n")
tab <- data.frame(estimand=c("E1 MASLD 25 vs ctrl 23","E1b MASLD 25 vs FL-series ctrl 3","E2 MASH 19 vs MASL 6","E3 ordinal ctrl23/MASL6/MASH19 slope per step","E3b ordinal FL ctrl3/MASL6/MASH19 slope per step"),
  alpha_single=c(mde2(25,23,a1),mde2(25,3,a1),mde2(19,6,a1),mdetrend(c(23,6,19),a1),mdetrend(c(3,6,19),a1)),
  alpha_117=c(mde2(25,23,aF),mde2(25,3,aF),mde2(19,6,aF),mdetrend(c(23,6,19),aF),mdetrend(c(3,6,19),aF)))
print(tab, digits=3)
cat("\nExpected power, Ductular injury (hep 20) and Stromal ECM (hep 8), from discovery standardized means\n")
for (mod in c(20,8)) { r <- out[out$cell_type=="hepatocytes"&out$module==mod,]
  dc <- (6*r$m1+19*r$m2)/25  # disease mean vs ctrl in this cohort's mix
  cat(sprintf("hep %d: m1=%.2f m2=%.2f | E1 d=%.2f power a.05=%.2f a117=%.3f | E2 d=%.2f power a.05=%.2f | E3 power a.05=%.2f a117=%.3f\n",
    mod, r$m1, r$m2, dc, pow2(25,23,dc,a1), pow2(25,23,dc,aF), r$m2-r$m1, pow2(19,6,r$m2-r$m1,a1), powtrend(c(23,6,19),c(0,r$m1,r$m2),a1), powtrend(c(23,6,19),c(0,r$m1,r$m2),aF)))}
out$dc <- (6*out$m1+19*out$m2)/25
out$powE1_117 <- sapply(out$dc, function(x) pow2(25,23,abs(x),aF))
out$powE3_117 <- mapply(function(a,b) powtrend(c(23,6,19),c(0,a,b),aF), out$m1, out$m2)
cat("\n117-family: programs with >=80% power at alpha 0.05/117: E1", sum(out$powE1_117>=.8), " E3", sum(out$powE3_117>=.8), "\n")
cat("expected number of 117 programs detected (sum of powers, alpha 0.05/117): E1", round(sum(out$powE1_117),1), " E3", round(sum(out$powE3_117),1), "\n")
cat("|dc| quantiles:"); print(round(quantile(abs(out$dc),c(.5,.75,.9,1)),2))
write.table(tab, paste0(G,"mde_table.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(out, paste0(G,"expected_power_117.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
writeLines(capture.output(sessionInfo()), paste0(G,"sessionInfo.txt"))
