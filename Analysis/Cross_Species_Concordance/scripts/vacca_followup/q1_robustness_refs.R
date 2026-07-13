#!/usr/bin/env Rscript
# Robustness of the Q1 two-arm result to the choice of metabolic reference.
# The default metab ref (Mild-vs-Ctrl + Moderate + EPoS-Moderate) is r=0.89 to the
# fibrotic ref. Test a PURER metabolic ref (Mild-vs-Control ONLY) to confirm the
# MCD-below-HFD metabolic ranking is not an artifact of ref correlation.
suppressPackageStartupMessages({ library(data.table); library(readxl) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE,"data/external/vacca_2024"); CS <- file.path(BASE,"Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(CS,"vacca_benchmark")
pnrm <- function(x) toupper(trimws(gsub("_"," ",sub("^KEGG_","",x))))
nrm  <- function(x) toupper(gsub("[^A-Za-z0-9]","",x))
sp   <- function(a,b) suppressWarnings(cor(a,b,method="spearman",use="complete.obs"))
pr   <- function(a,b) suppressWarnings(cor(a,b,method="pearson", use="complete.obs"))
p6 <- as.data.table(read_excel(file.path(VACCA,"42255_2024_1043_MOESM6_ESM.xlsx"),sheet="NES"))
setnames(p6,1,"pathway"); p6[,pname:=pnrm(pathway)]
nf <- function(c) suppressWarnings(as.numeric(p6[[c]]))
refs <- list(
  metab_pure   = nf("UCAM/VCU: Mild vs Control"),
  metab_default= rowMeans(cbind(nf("UCAM/VCU: Mild vs Control"),nf("UCAM/VCU: Moderate vs Mild"),nf("EPoS: Moderate vs Mild")),na.rm=TRUE),
  fibro        = rowMeans(cbind(nf("UCAM/VCU: Severe vs Mild"),nf("EPoS: Severe vs Mild")),na.rm=TRUE))
hr <- data.table(pname=p6$pname, metab_pure=refs$metab_pure, metab_default=refs$metab_default, fibro=refs$fibro)
mr <- fread(file.path(CS,"fgsea_mouse_results.csv"))[grepl("^KEGG_",pathway)]; mr[,pname:=pnrm(pathway)]
shared <- intersect(unique(mr$pname), hr[is.finite(metab_pure)&is.finite(fibro)]$pname)
hh <- hr[match(shared,pname)]
cat(sprintf("shared=%d ; r(metab_pure,fibro)=%.2f ; r(metab_default,fibro)=%.2f\n",
            length(shared), pr(hh$metab_pure,hh$fibro), pr(hh$metab_default,hh$fibro)))
diets <- c("MCD","HFD","CDAHFD","FPC")
out <- rbindlist(lapply(diets, function(d){
  v <- mr[source==d][match(shared,pname)]$NES
  data.table(diet=d, metab_pure=pr(v,hh$metab_pure), metab_default=pr(v,hh$metab_default), fibro=pr(v,hh$fibro))
}))
cat("\n=== our diets under PURE vs DEFAULT metabolic ref ===\n")
out[, rank_pure := frank(-metab_pure)][, rank_default := frank(-metab_default)]
print(out[order(rank_pure), .(diet, metab_pure=round(metab_pure,3), rank_pure,
                              metab_default=round(metab_default,3), rank_default, fibro=round(fibro,3))])
mcd <- out[diet=="MCD"]; hfd <- out[diet=="HFD"]; cda <- out[diet=="CDAHFD"]
cat(sprintf("\nMCD metab rank pure=%d default=%d ; MCD<HFD pure=%s default=%s ; MCD<CDAHFD pure=%s default=%s\n",
            mcd$rank_pure, mcd$rank_default,
            mcd$metab_pure<hfd$metab_pure, mcd$metab_default<hfd$metab_default,
            mcd$metab_pure<cda$metab_pure, mcd$metab_default<cda$metab_default))
fwrite(out, file.path(OUT,"q1_robustness_refs.csv"))
