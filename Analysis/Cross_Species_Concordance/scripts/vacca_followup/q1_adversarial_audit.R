#!/usr/bin/env Rscript
# ============================================================================
# q1_adversarial_audit.R  — ADVERSARIAL CRITIC of the Q1 two-arm resolution
#
# Independently re-derives the load-bearing Q1 numbers and stress-tests the
# three claims a reviewer will attack:
#   C1. "Two-arm PEP reproduces BOTH Vacca arms" — but the METABOLIC arm is weak.
#   C2. "MCD>CDAHFD/HFD is explained by the fibrotic-vs-metabolic arm split."
#   C3. "CDAHFD discrepancy is DATA (GSE162876 outlier), not METHOD."
#
# Tests:
#   T1  Validation strength of each arm on Vacca's OWN 33 mouse models
#       (Spearman PEP-proximity vs published per-arm DHPS) + bootstrap CI.
#   T2  Is the metabolic arm distinguishable from the fibrotic arm? Arm refs
#       inter-correlation; partial Spearman of metab-PEP vs metab-DHPS
#       CONTROLLING for fibro-DHPS (does metab arm carry independent signal?).
#   T3  n=4 fragility: the headline 2x2 ranks 4 diets. Bootstrap the SHARED
#       pathway set (resample pathways) and report how often MCD's metabolic
#       rank stays >= CDAHFD & HFD, and how often the fibrotic ranking is stable.
#   T4  CDAHFD data-not-method: our-CDAHFD vs Vacca-CDAHFD (0.81) vs
#       within-Vacca-CDAHFD ceiling (0.94) vs our-CDAHFD vs random-Vacca (0.52).
#       Also: is our-CDAHFD's BEST Vacca match actually a CDAHFD model, or a
#       different-diet model (=method/normalization artifact)?
# Env: rnaseq
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(readxl) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(BASE, "results/vacca_benchmark")
set.seed(42)
pnrm <- function(x) toupper(trimws(gsub("_", " ", sub("^KEGG_", "", x))))
nrm  <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
sp   <- function(a,b) suppressWarnings(cor(a,b,method="spearman",use="complete.obs"))
pr   <- function(a,b) suppressWarnings(cor(a,b,method="pearson", use="complete.obs"))

## ---- Load Vacca pathway NES + per-arm DHPS -------------------------------
p6 <- as.data.table(read_excel(file.path(VACCA,"42255_2024_1043_MOESM6_ESM.xlsx"),sheet="NES"))
setnames(p6,1,"pathway"); p6[,pname:=pnrm(pathway)]
num <- function(col) suppressWarnings(as.numeric(p6[[col]]))
metab_ref <- rowMeans(cbind(num("UCAM/VCU: Mild vs Control"),
                            num("UCAM/VCU: Moderate vs Mild"),
                            num("EPoS: Moderate vs Mild")), na.rm=TRUE)
fibro_ref <- rowMeans(cbind(num("UCAM/VCU: Severe vs Mild"),
                            num("EPoS: Severe vs Mild")), na.rm=TRUE)
human_cols <- grep("UCAM|EPoS",names(p6),value=TRUE)
model_cols <- setdiff(names(p6)[-1], c(human_cols,"pname"))
for (cc in model_cols) suppressWarnings(p6[[cc]] <- as.numeric(p6[[cc]]))

m9 <- as.data.table(read_excel(file.path(VACCA,"42255_2024_1043_MOESM9_ESM.xlsx"),
                               col_names=FALSE, skip=2))
setnames(m9,1:10,c("model","diet_group","mP","mH","m_DHPS","mM","fP","fH","f_DHPS","fM"))
m9 <- m9[!is.na(model)]; for(cc in c("m_DHPS","f_DHPS")) m9[[cc]]<-as.numeric(m9[[cc]])
m9[,key:=nrm(model)]

mr <- fread(file.path(CS,"fgsea_mouse_results.csv"))[grepl("^KEGG_",pathway)]
mr[,pname:=pnrm(pathway)]

href <- data.table(pname=p6$pname, metab=metab_ref, fibro=fibro_ref)
shared_human <- href[is.finite(metab)&is.finite(fibro)]$pname
shared <- intersect(unique(mr$pname), shared_human)
hr_m <- href[match(shared,pname)]$metab
hr_f <- href[match(shared,pname)]$fibro

cat("================== Q1 ADVERSARIAL AUDIT ==================\n")
cat(sprintf("Shared KEGG pathways (ours x Vacca, both arms finite): %d\n", length(shared)))
cat(sprintf("[T2] arm-reference INDEPENDENCE: Pearson(metab_ref,fibro_ref)=%.3f  Spearman=%.3f\n",
            pr(hr_m,hr_f), sp(hr_m,hr_f)))

## ---- Per-model two-arm proximity on Vacca models -------------------------
vmod <- rbindlist(lapply(model_cols, function(mc){
  v <- p6[match(shared,pname)][[mc]]
  dt <- data.table(model=mc, metab_prox=pr(v,hr_m), fibro_prox=pr(v,hr_f))
  dt[, key := nrm(mc)][]
}))
vmod <- merge(vmod, m9[,.(key,diet_group,m_DHPS,f_DHPS)], by="key")[!grepl("^R[-.]",model)]
N <- nrow(vmod)
cat(sprintf("\n[T1] Vacca mouse models n=%d\n", N))

## T1: arm validation + bootstrap CI over models
boot_sp <- function(x,y,B=2000){
  v <- replicate(B,{i<-sample.int(length(x),replace=TRUE); sp(x[i],y[i])})
  quantile(v, c(.025,.975), na.rm=TRUE)
}
m_val <- sp(vmod$metab_prox, vmod$m_DHPS); m_ci <- boot_sp(vmod$metab_prox, vmod$m_DHPS)
f_val <- sp(vmod$fibro_prox, vmod$f_DHPS); f_ci <- boot_sp(vmod$fibro_prox, vmod$f_DHPS)
mx    <- sp(vmod$metab_prox, vmod$f_DHPS)  # cross
fx    <- sp(vmod$fibro_prox, vmod$m_DHPS)
cat(sprintf("  METAB arm  : metab_prox vs m_DHPS  Spearman=%.3f  95%%CI[%.3f,%.3f]\n", m_val,m_ci[1],m_ci[2]))
cat(sprintf("  FIBRO arm  : fibro_prox vs f_DHPS  Spearman=%.3f  95%%CI[%.3f,%.3f]\n", f_val,f_ci[1],f_ci[2]))
cat(sprintf("  CROSS      : metab_prox vs f_DHPS=%.3f ; fibro_prox vs m_DHPS=%.3f\n", mx, fx))
cat(sprintf("  >> metab CI crosses 0: %s ; metab arm > its OWN cross (fibro_prox vs m_DHPS): %s\n",
            (m_ci[1] < 0), (m_val > fx)))

## T2: does metab-arm carry signal INDEPENDENT of fibrotic? partial corr (rank)
rank_resid <- function(y, given){
  r <- residuals(lm(rank(y) ~ rank(given))); r
}
pm <- cor(rank_resid(vmod$metab_prox, vmod$f_DHPS),
          rank_resid(vmod$m_DHPS,     vmod$f_DHPS))
cat(sprintf("\n[T2] PARTIAL Spearman metab_prox~m_DHPS | f_DHPS = %.3f  (independent metabolic signal?)\n", pm))

## T3: n=4 fragility — bootstrap the shared pathway set, re-rank our 4 diets
diet_nes <- sapply(c("MCD","HFD","CDAHFD","FPC"), function(d)
  mr[source==d][match(shared,pname)]$NES)
boot_rank <- function(B=2000){
  keep_mcd_below <- 0L; fpc_top_fibro <- 0L; mcd_top2_fibro <- 0L
  cd_bottom_fibro <- 0L
  for(b in 1:B){
    i <- sample.int(length(shared), replace=TRUE)
    mp <- sapply(c("MCD","HFD","CDAHFD","FPC"), function(d) pr(diet_nes[i,d], hr_m[i]))
    fp <- sapply(c("MCD","HFD","CDAHFD","FPC"), function(d) pr(diet_nes[i,d], hr_f[i]))
    mr_ <- rank(-mp); fr_ <- rank(-fp)
    if(mr_["MCD"]>mr_["CDAHFD"] & mr_["MCD"]>mr_["HFD"]) keep_mcd_below <- keep_mcd_below+1L
    if(which.max(fp)==which(c("MCD","HFD","CDAHFD","FPC")=="FPC")) fpc_top_fibro <- fpc_top_fibro+1L
    if(fr_["MCD"]<=2) mcd_top2_fibro <- mcd_top2_fibro+1L
    if(which.min(fp)==which(c("MCD","HFD","CDAHFD","FPC")=="CDAHFD")) cd_bottom_fibro <- cd_bottom_fibro+1L
  }
  c(mcd_metab_below_cdahfd_and_hfd = keep_mcd_below/B,
    fpc_is_fibro_top   = fpc_top_fibro/B,
    mcd_fibro_top2     = mcd_top2_fibro/B,
    cdahfd_fibro_bottom= cd_bottom_fibro/B)
}
br <- boot_rank()
cat("\n[T3] n=4 ranking FRAGILITY (bootstrap shared pathways, B=2000): fraction of resamples where...\n")
for(nm in names(br)) cat(sprintf("     %-30s = %.3f\n", nm, br[nm]))

## T4: CDAHFD data-not-method
vacca_cdahfd <- c("6J-CDAHFD-F45-6W","6J-CDAHFD-F45-8W","6J-CDAHFD-F45-12W")
ourNES <- function(d){v<-mr[source==d]; setNames(v$NES,v$pname)}
href_idx <- match(shared, p6$pname)
corr_one <- function(d, mc){
  ov <- ourNES(d)[shared]; vv <- p6[[mc]][href_idx]
  data.table(vacca_model=mc, pearson=pr(ov,vv), spearman=sp(ov,vv))
}
cd_all <- rbindlist(lapply(model_cols, function(mc) corr_one("CDAHFD",mc)))
cd_all[, is_cdahfd := vacca_model %in% vacca_cdahfd]
setorder(cd_all,-pearson); cd_all[,rank:=.I]
best_match <- cd_all[1]
best_cd    <- cd_all[is_cdahfd==TRUE][which.max(pearson)]
cat("\n[T4] CDAHFD: our-CDAHFD best overall match + best CDAHFD-specific match\n")
cat(sprintf("     BEST overall Vacca match for our CDAHFD: %s (r=%.3f, is a Vacca-CDAHFD model: %s)\n",
            best_match$vacca_model, best_match$pearson, best_match$is_cdahfd))
cat(sprintf("     BEST Vacca-CDAHFD match: %s (r=%.3f, overall rank %d/%d)\n",
            best_cd$vacca_model, best_cd$pearson, best_cd$rank, nrow(cd_all)))
cat(sprintf("     within-Vacca-CDAHFD ceiling (mean pairwise Pearson): %.3f\n",
            mean(combn(vacca_cdahfd,2,function(p) pr(p6[[p[1]]][href_idx],p6[[p[2]]][href_idx])))))
cat(sprintf("     our-CDAHFD vs Vacca-CDAHFD mean=%.3f ; vs all-other-Vacca mean=%.3f\n",
            mean(cd_all[is_cdahfd==TRUE]$pearson), mean(cd_all[is_cdahfd==FALSE]$pearson)))

## Also: how does our CDAHFD's gap compare to FPC's gap (is CDAHFD uniquely bad,
## or are ALL our diets ~equally distant = method/normalization, not data)?
ourFPC_vs_WD <- {
  wd <- grep("-WD-", model_cols, value=TRUE)
  mean(rbindlist(lapply(wd, function(mc) corr_one("FPC",mc)))$pearson)
}
ourMCD_vs_MCD <- {
  vm <- grep("MCD", model_cols, value=TRUE)
  if(length(vm)) mean(rbindlist(lapply(vm, function(mc) corr_one("MCD",mc)))$pearson) else NA
}
cat(sprintf("\n[T4b] cross-diet self-match Pearson: CDAHFD->VaccaCDAHFD=%.3f ; FPC->VaccaWD=%.3f ; MCD->VaccaMCD=%.3f\n",
            mean(cd_all[is_cdahfd==TRUE]$pearson), ourFPC_vs_WD, ourMCD_vs_MCD))

## ---- write audit summary -------------------------------------------------
audit <- data.table(
  metric = c("n_shared_pathways","n_vacca_models","arm_ref_independence_pearson",
             "metab_arm_validation_spearman","metab_arm_ci_lo","metab_arm_ci_hi",
             "fibro_arm_validation_spearman","fibro_arm_ci_lo","fibro_arm_ci_hi",
             "cross_metabProx_vs_fDHPS","cross_fibroProx_vs_mDHPS",
             "metab_partial_indep_spearman",
             "boot_MCD_metab_below_both","boot_FPC_fibro_top","boot_MCD_fibro_top2","boot_CDAHFD_fibro_bottom",
             "ourCDAHFD_vs_VaccaCDAHFD","within_VaccaCDAHFD_ceiling","ourCDAHFD_vs_otherVacca",
             "ourCDAHFD_best_match_is_cdahfd","ourCDAHFD_best_cdahfd_overall_rank",
             "ourFPC_vs_VaccaWD"),
  value = c(length(shared), N, round(pr(hr_m,hr_f),3),
            round(m_val,3), round(m_ci[1],3), round(m_ci[2],3),
            round(f_val,3), round(f_ci[1],3), round(f_ci[2],3),
            round(mx,3), round(fx,3),
            round(pm,3),
            round(br["mcd_metab_below_cdahfd_and_hfd"],3), round(br["fpc_is_fibro_top"],3),
            round(br["mcd_fibro_top2"],3), round(br["cdahfd_fibro_bottom"],3),
            round(mean(cd_all[is_cdahfd==TRUE]$pearson),3),
            round(mean(combn(vacca_cdahfd,2,function(p) pr(p6[[p[1]]][href_idx],p6[[p[2]]][href_idx]))),3),
            round(mean(cd_all[is_cdahfd==FALSE]$pearson),3),
            as.integer(best_match$is_cdahfd), best_cd$rank,
            round(ourFPC_vs_WD,3)))
fwrite(audit, file.path(OUT,"q1_adversarial_audit.csv"))
cat(sprintf("\nWrote: %s\n", file.path(OUT,"q1_adversarial_audit.csv")))
cat("=========================================================\n")
