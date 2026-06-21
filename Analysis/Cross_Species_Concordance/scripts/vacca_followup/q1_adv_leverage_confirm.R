#!/usr/bin/env Rscript
# Confirm the single-pathway leverage finding & FPC's RELATIVE robustness.
# (1) FPC r with the top contributor (Gly/Ser/Thr metab) removed, vs other diets.
# (2) Rank stability: across all single-pathway drops, does FPC stay rank-1 among OUR 4?
# (3) Bootstrap pathways (resample 28 pathways w/ replacement, 2000x): CI on each diet's r
#     and P(FPC = max of the 4). Robust if FPC stays top even resampling pathways.
suppressPackageStartupMessages({ library(data.table); library(readxl) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE,"data/external/vacca_2024"); CS <- file.path(BASE,"Analysis/Cross_Species_Concordance/results")
RES   <- file.path(BASE,"results/vacca_benchmark"); set.seed(42)
pnrm <- function(x) toupper(trimws(gsub("_"," ",sub("^KEGG_","",x))))
pear <- function(a,b) suppressWarnings(cor(a,b,method="pearson",use="complete.obs"))

p6 <- as.data.table(read_excel(file.path(VACCA,"42255_2024_1043_MOESM6_ESM.xlsx"),sheet="NES")); setnames(p6,1,"pathway")
p6[, pname := pnrm(pathway)]
fib_mat <- sapply(c("UCAM/VCU: Severe vs Mild","EPoS: Severe vs Mild"), function(cc) suppressWarnings(as.numeric(p6[[cc]])))
href <- data.table(pname=p6$pname, h_nes=rowMeans(fib_mat,na.rm=TRUE))[is.finite(h_nes)]
mr <- fread(file.path(CS,"fgsea_mouse_results.csv"))[grepl("^KEGG_",pathway)]; mr[, pname := pnrm(pathway)]
diets <- c("FPC","MCD","HFD","CDAHFD"); shared <- intersect(unique(mr$pname),href$pname)
href_v <- href[match(shared,pname)]$h_nes

# build aligned matrix of NES per diet over shared pathways (complete cases only)
M <- sapply(diets, function(d) mr[source==d][match(shared,pname)]$NES)
keep <- complete.cases(M) & is.finite(href_v)
M <- M[keep,,drop=FALSE]; h <- href_v[keep]; pn <- shared[keep]
cat(sprintf("Aligned pathways (complete across all 4 diets + human): %d\n", nrow(M)))

r_full <- apply(M,2,function(v) pear(v,h)); names(r_full) <- diets
cat("\nFull r per diet:\n"); print(round(r_full,3))

# (1) drop the FPC top contributor pathway from ALL diets, recompute
zc <- scale(M[,"FPC"])[,1]; zh <- scale(h)[,1]; contrib <- zc*zh/(nrow(M)-1)
top_path <- pn[which.max(contrib)]
cat(sprintf("\nFPC top-contributor pathway = '%s' (contrib %.3f of r=%.3f)\n",
            top_path, max(contrib), r_full["FPC"]))
idx <- which(pn != top_path)
r_noTop <- apply(M[idx,,drop=FALSE],2,function(v) pear(v,h[idx]))
cat("r after removing FPC's #1 pathway (all diets):\n"); print(round(r_noTop,3))
cat(sprintf(">> FPC drops %.3f -> %.3f ; still rank-1 of 4? %s\n",
            r_full["FPC"], r_noTop["FPC"], r_noTop["FPC"]==max(r_noTop)))

# (2) single-pathway-drop: does FPC stay rank-1 among the 4 across ALL drops?
stay1 <- sapply(seq_len(nrow(M)), function(i){
  rr <- apply(M[-i,,drop=FALSE],2,function(v) pear(v,h[-i])); which.max(rr)==which(diets=="FPC")
})
cat(sprintf("\n(2) FPC stays rank-1 of 4 in %d/%d single-pathway-drop jackknives (%.0f%%)\n",
            sum(stay1), length(stay1), 100*mean(stay1)))

# (3) bootstrap pathways (resample rows w/ replacement)
B <- 2000
boot <- replicate(B, {
  s <- sample(nrow(M), replace=TRUE)
  apply(M[s,,drop=FALSE],2,function(v) pear(v,h[s]))
})  # 4 x B
fpc_is_max <- mean(boot["FPC",] == apply(boot,2,max))
ci <- t(apply(boot,1,quantile,c(.025,.5,.975),na.rm=TRUE))
cat("\n(3) Pathway-bootstrap (B=2000) 95% CI on r per diet:\n")
print(round(ci,3))
cat(sprintf(">> P(FPC is the MAX of 4 diets across pathway bootstraps) = %.3f\n", fpc_is_max))

out <- data.table(
  metric=c("n_aligned_pathways","fpc_r_full","fpc_top_pathway","fpc_top_contrib",
           "fpc_r_drop_top1","fpc_stays_rank1_after_drop_top1",
           "fpc_rank1_frac_single_drop","fpc_boot_ci_lo","fpc_boot_ci_hi","P_fpc_is_max_boot"),
  value=c(nrow(M), round(r_full["FPC"],3), top_path, round(max(contrib),3),
          round(r_noTop["FPC"],3), r_noTop["FPC"]==max(r_noTop),
          round(mean(stay1),3), round(ci["FPC",1],3), round(ci["FPC",3],3), round(fpc_is_max,3)))
fwrite(out, file.path(RES,"q1_adv_leverage_confirm.csv"))
cat("\nWrote:", file.path(RES,"q1_adv_leverage_confirm.csv"),"\n")
