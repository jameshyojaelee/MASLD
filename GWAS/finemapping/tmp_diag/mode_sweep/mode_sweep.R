#!/usr/bin/env Rscript
# TASK 16 + liftover hypothesis.
# (A) palindromic-vs-nonpalindromic "opposite" enrichment, ALL 50 studies x ALL 22 chr,
#     driving the CURRENT shipped helpers (panel_freq_of now pair-matches, 5 args).
# (B) stratify enrichment by AF-build mode: existing / join (hg19 native) / join_lift / none.
#     If enrichment is confined to join_lift -> position-error artifact.
# (C) DIRECT test of the position-error mechanism: a wrong-position join silently
#     REMOVES non-palindromic variants (letter mismatch) while silently ACCEPTING
#     palindromes. So palindromic "opposite" calls should sit in neighbourhoods with
#     an elevated non-palindromic JOIN-FAILURE rate. Measured in 100kb windows.
suppressPackageStartupMessages(library(data.table))
FM_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
OUT <- file.path(FM_DIR, "tmp_diag/mode_sweep")
src <- readLines(file.path(FM_DIR,"src/06_susie_coloc.R"))
i0 <- grep("^PANEL_AF_ROOT <- ", src); cl <- grep("^\\}$", src)
i1 <- min(cl[cl > grep("^strand_verdict_vs_panel <- ", src)])
eval(parse(text=paste(src[i0:i1], collapse="\n")))
stopifnot(length(formals(panel_freq_of)) == 5L)
cat("helpers sourced; panel_freq_of args:", paste(names(formals(panel_freq_of)),collapse=","), "\n")

reg   <- fread(file.path(FM_DIR,"config/gwas_registry.tsv"))
afsrc <- fread(file.path(FM_DIR,"config/gwas_af_sources.tsv"))
reg <- merge(reg[, .(study_name, ancestry, sumstats_path)],
             afsrc[, .(study_name, mode, af_sumstats_path, af_coverage)], by="study_name", all.x=TRUE)
reg[is.na(mode), mode := "absent"]
setorder(reg, ancestry, study_name)      # group by ancestry so the sidecar cache is reused
CHRS <- 1:22

pac <- new.env(); pac_anc <- ""
get_pa <- function(anc, CH) {
  if (!identical(anc, pac_anc)) { rm(list=ls(pac), envir=pac); pac_anc <<- anc; gc(verbose=FALSE) }
  k <- paste0(anc, "_", CH)
  if (is.null(pac[[k]])) pac[[k]] <- load_panel_af(anc, CH)
  pac[[k]]
}
rows <- list(); nb <- list()
for (si in seq_len(nrow(reg))) {
  S <- reg$study_name[si]; anc <- reg$ancestry[si]; md <- reg$mode[si]
  gp <- if (!is.na(reg$af_sumstats_path[si]) && file.exists(reg$af_sumstats_path[si]))
          reg$af_sumstats_path[si] else reg$sumstats_path[si]
  if (!file.exists(gp)) { cat("SKIP", S, "\n"); next }
  hdr <- tryCatch(names(fread(gp, nrows=0L)), error=function(e) NULL)
  if (is.null(hdr) || !all(c("chromosome","position","allele1","allele2") %in% hdr)) { cat("SKIP schema", S, "\n"); next }
  has_af <- "af" %in% hdr
  g <- fread(gp, select=c("chromosome","position","allele1","allele2", if (has_af) "af"), showProgress=FALSE)
  if (!has_af) g[, af := NA_real_]
  for (CH in CHRS) {
    x <- g[chromosome == CH]; if (!nrow(x)) next
    pa <- get_pa(anc, CH); if (is.null(pa)) next
    x <- unique(x[, .(position, a1=toupper(allele1), a2=toupper(allele2), gaf=af)], by="position")
    x[pa, on=.(position), `:=`(p1=i.bim_a1, p2=i.bim_a2, paf=i.af_a1)]
    x[, ref := panel_freq_of(a1, a2, p1, p2, paf)]          # SHIPPED, pair-matching
    x[, v   := strand_verdict_vs_panel(gaf, ref)]
    x[, is_pal := (a1 %in% c("A","T") & a2 %in% c("A","T")) | (a1 %in% c("C","G") & a2 %in% c("C","G"))]
    x[, absent    := is.na(p1)]                      # position not in the panel at all (benign)
    x[, pairmiss  := !is.na(p1) & is.na(ref)]        # position present, alleles do NOT pair-match
    x[, join_fail := is.na(ref)]
    # (C) 100kb-window non-palindromic join-failure rate
    x[, win := position %/% 100000L]
    wf <- x[is_pal == FALSE, .(nonpal_n=.N, nonpal_fail=mean(join_fail),
                               nonpal_pairmiss=mean(pairmiss), nonpal_absent=mean(absent)), by=win]
    x <- merge(x, wf, by="win", all.x=TRUE)
    np <- x[is_pal==FALSE]; pl <- x[is_pal==TRUE]
    rows[[length(rows)+1L]] <- data.table(study=S, anc=anc, mode=md, chr=CH, has_af=has_af,
      n_nonpal=nrow(np), nonpal_verdict=np[!is.na(v),.N], nonpal_opp=np[v %in% "opposite",.N],
      nonpal_joinfail=np[join_fail==TRUE,.N], nonpal_pairmiss=np[pairmiss==TRUE,.N],
      nonpal_absent=np[absent==TRUE,.N],
      n_pal=nrow(pl), pal_verdict=pl[!is.na(v),.N], pal_opp=pl[v %in% "opposite",.N],
      pal_joinfail=pl[join_fail==TRUE,.N], pal_pairmiss=pl[pairmiss==TRUE,.N],
      pal_absent=pl[absent==TRUE,.N])
    if (pl[v %in% "opposite", .N] > 0)
      nb[[length(nb)+1L]] <- pl[!is.na(v), .(study=S, mode=md, chr=CH, position,
            grp=fifelse(v=="opposite","pal_opposite","pal_same"), nonpal_fail,
            nonpal_pairmiss, nonpal_absent, nonpal_n)]
  }
  rm(g); gc(verbose=FALSE)
  cat(sprintf("  done %-40s %-4s %-10s\n", S, anc, md)); flush.console()
}
res <- rbindlist(rows); fwrite(res, file.path(OUT,"mode_sweep_raw.tsv"), sep="\t")
nbd <- rbindlist(nb); if (nrow(nbd)) fwrite(nbd, file.path(OUT,"neighbourhood.tsv"), sep="\t")

a <- res[, .(mode=mode[1], anc=anc[1], has_af=all(has_af),
   nonpal_v=sum(nonpal_verdict), nonpal_opp=sum(nonpal_opp),
   pal_v=sum(pal_verdict), pal_opp=sum(pal_opp)), by=study]
a[, `:=`(ctrl_pct=100*nonpal_opp/nonpal_v, pal_pct=100*pal_opp/pal_v)]
a[, ratio := pal_pct/ctrl_pct]
a[, p := mapply(function(po,pv,no,nv) if (pv>0 && nv>0)
      fisher.test(matrix(c(po,pv-po,no,nv-no),nrow=2))$p.value else NA_real_,
      pal_opp, pal_v, nonpal_opp, nonpal_v)]
a[, p_bonf := p.adjust(p, method="bonferroni")]
setorder(a, p)
fwrite(a, file.path(OUT,"mode_sweep_perstudy.tsv"), sep="\t")
cat("\n=== PER-STUDY ENRICHMENT (all 22 chr) ===\n")
print(a[, .(study, mode, anc, ctrl_pct=round(ctrl_pct,4), pal_pct=round(pal_pct,4),
            ratio=round(ratio,2), p=signif(p,3), sig=p_bonf<0.05)], nrows=60)
cat("\n=== BY AF-BUILD MODE (pooled) ===\n")
m <- a[, .(studies=.N, sig=sum(p_bonf<0.05, na.rm=TRUE),
           nonpal_opp=sum(nonpal_opp), nonpal_v=sum(nonpal_v),
           pal_opp=sum(pal_opp), pal_v=sum(pal_v)), by=mode]
m[, `:=`(ctrl_pct=round(100*nonpal_opp/nonpal_v,4), pal_pct=round(100*pal_opp/pal_v,4))]
m[, ratio := round(pal_pct/ctrl_pct,2)]
print(m[order(-studies)])
cat("\n=== (C) POSITION-ERROR TEST: non-palindromic JOIN-FAILURE rate in the 100kb window ===\n")
if (nrow(nbd)) {
  print(nbd[, .(n=.N, mean_joinfail=round(mean(nonpal_fail,na.rm=TRUE),5),
                mean_PAIRMISS=round(mean(nonpal_pairmiss,na.rm=TRUE),5),
                mean_absent=round(mean(nonpal_absent,na.rm=TRUE),5)), by=.(mode,grp)][order(mode,grp)])
  cat("\n  Wilcoxon pal_opposite vs pal_same (pooled):\n")
  o <- nbd[grp=="pal_opposite", nonpal_pairmiss]; s <- nbd[grp=="pal_same", nonpal_pairmiss]
  if (length(o)>2 && length(s)>2) {
    w <- wilcox.test(o, s)
    cat(sprintf("    [nonpal PAIR-MISMATCH rate] opposite n=%d mean=%.5f | same n=%d mean=%.5f | W p=%.3g\n",
        length(o), mean(o,na.rm=TRUE), length(s), mean(s,na.rm=TRUE), w$p.value))
  }
}
cat("\nwrote:", OUT, "\n")
