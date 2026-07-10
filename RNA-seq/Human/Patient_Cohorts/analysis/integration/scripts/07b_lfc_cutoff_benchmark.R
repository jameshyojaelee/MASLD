#!/usr/bin/env Rscript
# 07b_lfc_cutoff_benchmark.R
# ---------------------------------------------------------------------------
# Data-driven selection of the bulk Tier-1 DEG cutoff: which SCALE
# (raw |logFC| vs ashr-shrunk |logFC| vs lfsr+shrunk-magnitude) and which
# THRESHOLD best recovers known MASLD biology while staying on the CV-stability
# plateau. Read-only on pipeline inputs; writes only its own results dir.
#
# Inputs : canonical_deg_results.csv (logFC, shrunk_logFC, lfsr, padj, AveExpr,
#            symbol) — the limma-voom-qw C2 canonical (promoted 2026-06-08,
#            superseding dream_results_ashr.csv, which is retained on disk only
#            as a retired-method sensitivity arm)
#          positive_control_validation.csv (62 Expression_driven + 3 GWAS_variant)
#          published panels (Govaere25 / Feng / SteatoSITE) — vectors below
#          drug targets (THRB, DGAT2, SCD, HSD17B13, ...) from Script 40
# Output : results/audit_sensitivity/lfc_cutoff_benchmark/
#            cutoff_sweep.csv, auroc_by_scale.csv, matched_N.csv, cutoff_benchmark.pdf
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT  <- file.path(BASE, "RNA-seq/results/audit_sensitivity/lfc_cutoff_benchmark")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- load canonical (limma-voom-qw C2) + DE-DUPLICATED ensembl->symbol map ---
a <- fread(file.path(INT, "canonical_deg_results.csv"))
a[, eb := sub("\\.[0-9]+$", "", gene)]
stopifnot(!anyDuplicated(a$eb))                       # canonical rows must be unique per gene
map <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
             select = c("ensembl_id", "human_symbol"))
map[, eb := sub("\\.[0-9]+$", "", ensembl_id)]
map <- unique(map[!is.na(human_symbol) & human_symbol != ""], by = "eb")   # DEDUP
a[, human_symbol := map[.(a$eb), on = "eb", human_symbol]]
# canonical carries its own `symbol`; fall back to it where the atlas map missed
if ("symbol" %in% names(a))
  a[is.na(human_symbol) & !is.na(symbol) & symbol != "", human_symbol := symbol]
cat(sprintf("canonical genes: %d (%d with symbol); rows unique per gene: %s\n",
            nrow(a), sum(!is.na(a$human_symbol)), !anyDuplicated(a$eb)))

# --- positive sets ---
pc <- fread(file.path(BASE, "RNA-seq/results/validation/positive_control_validation.csv"))
expr_ctrl <- unique(pc[control_type == "Expression_driven", gene])
gw_ctrl   <- unique(pc[control_type == "GWAS_variant", gene])
govaere <- c("AKR1B10","DUSP6","GDF15","THBS2","A2M","CDH2","COL1A1","COL3A1","COL4A1","COL4A2",
             "COL6A3","DCN","FBN1","FSTL1","IGFBP7","LUM","MFAP4","MMP2","POSTN","SPARC","SPP1","TAGLN","THY1","TIMP1","VCAN")
feng <- c("EFHD1","MLIP","TREM2","SPP1","GPNMB","CCL2","CCL20","CXCL1","CXCL6","IL1B","IL32",
          "COL1A1","COL1A2","COL3A1","FN1","LOX","LOXL2","ACTA2","PDGFRB","TGFB1","SERPINE1","MMP9","AKR1B10","GDF15","THY1","THBS2","LUM")
steatosite <- c("MT1F","KCNH7","COL25A1","RASD2","CTGF","STC1","GDNF","PRRX1","FGF7","LCNL1","DPEP1","CHRDL2","LHX6","POU4F1","CDH16")
drug_expr <- c("THRB","DGAT2","SCD","HSD17B13")        # expression-relevant approved/pipeline drug targets
known_union <- unique(c(expr_ctrl, govaere, feng, steatosite, drug_expr))
panels <- list(expr_ctrl=expr_ctrl, govaere25=govaere, feng=feng, steatosite=steatosite,
               drug_expr=drug_expr, known_union=known_union, gwas_neg=gw_ctrl)
for(nm in names(panels)) panels[[nm]] <- intersect(panels[[nm]], a$human_symbol[!is.na(a$human_symbol)])
cat("panel sizes (mapped to atlas):\n"); for(nm in names(panels)) cat(sprintf("  %-12s %d\n", nm, length(panels[[nm]])))

sym_in <- function(degsyms, set) round(100*sum(set %in% degsyms)/max(length(set),1),1)

# --- (1) THRESHOLD-FREE AUROC per scale: does the scale rank known genes high? ---
auroc <- function(score, is_pos){ # Mann-Whitney
  score <- score[!is.na(score)]; is_pos <- is_pos[!is.na(score)]
  r <- rank(score); (sum(r[is_pos]) - sum(is_pos)*(sum(is_pos)+1)/2) / (sum(is_pos)*sum(!is_pos))
}
a[, pos_known := human_symbol %in% panels$known_union]
auroc_dt <- rbindlist(lapply(c("raw_absLFC","shrunk_absLFC","neglog10_lfsr"), function(sc){
  s <- switch(sc, raw_absLFC=abs(a$logFC), shrunk_absLFC=abs(a$shrunk_logFC), neglog10_lfsr=-log10(pmax(a$lfsr,1e-300)))
  data.table(scale=sc, auroc_known=round(auroc(s, a$pos_known),4))
}))
fwrite(auroc_dt, file.path(OUT,"auroc_by_scale.csv"))
cat("\n=== Threshold-free AUROC (rank known-MASLD genes above background) ===\n"); print(auroc_dt)

# --- (2) threshold sweep per scale ---
sweep_rows <- list()
add <- function(scale, thr, mask){
  d <- a[mask & !is.na(human_symbol), human_symbol]
  n <- a[mask, .N]
  sweep_rows[[length(sweep_rows)+1]] <<- data.table(
    scale=scale, threshold=thr, nDEG=n,
    recall_expr=sym_in(d,panels$expr_ctrl), recall_govaere=sym_in(d,panels$govaere25),
    recall_feng=sym_in(d,panels$feng), recall_steato=sym_in(d,panels$steatosite),
    recall_known=sym_in(d,panels$known_union), gwas_capture=sum(panels$gwas_neg %in% d),
    median_AveExpr=round(median(a[mask, AveExpr],na.rm=TRUE),2))
}
for(t in seq(0.10,1.00,0.05)) add("raw_absLFC",   t, a$padj<0.05 & abs(a$logFC)>t)
for(t in seq(0.05,0.70,0.05)) add("shrunk_absLFC",t, a$padj<0.05 & abs(a$shrunk_logFC)>t)
for(t in c(0.005,0.01,0.05)) for(mg in c(0,0.2,0.3,0.5)) add(sprintf("lfsr%.3f",t), mg, a$lfsr<t & abs(a$shrunk_logFC)>mg)
sweep <- rbindlist(sweep_rows)
fwrite(sweep, file.path(OUT,"cutoff_sweep.csv"))

# --- (3) CV-stability plateau per scale (CV of nDEG across padj thresholds) ---
padjs <- c(0.001,0.005,0.01,0.025,0.05,0.1)
cv_rows <- list()
for(sc in c("raw","shrunk")){ col <- if(sc=="raw") a$logFC else a$shrunk_logFC
  for(t in seq(0.1,1.0,0.05)){ ns <- sapply(padjs, function(p) sum(a$padj<p & abs(col)>t))
    cv_rows[[length(cv_rows)+1]] <- data.table(scale=sc, threshold=t, cv=round(sd(ns)/mean(ns)*100,1)) } }
cv <- rbindlist(cv_rows)
plateau <- cv[cv<10, .(first_stable_thr=min(threshold)), by=scale]
cat("\n=== CV<10% stability plateau onset per scale ===\n"); print(plateau)

# Persist the per-scale CV<10% crossing point (was ephemeral stdout only):
#   scale, threshold (first |LFC| with CV<10%), cv (its CV value).
plateau_cv <- merge(plateau, cv,
                    by.x = c("scale","first_stable_thr"),
                    by.y = c("scale","threshold"))[, .(scale, threshold = first_stable_thr, cv)]
fwrite(plateau_cv, file.path(OUT,"cutoff_cv_plateau.csv"))
cat("\n=== Persisted CV plateau crossing points (cutoff_cv_plateau.csv) ===\n"); print(plateau_cv)

# --- (4) matched-N recall (fair cross-scale at equal DEG count) ---
matchN <- function(targetN){
  rbindlist(lapply(c("raw","shrunk"), function(sc){ col <- if(sc=="raw") a$logFC else a$shrunk_logFC
    bt<-NA;bn<-Inf; for(t in seq(0.1,1.2,0.01)){n<-sum(a$padj<0.05 & abs(col)>t); if(abs(n-targetN)<abs(bn-targetN)){bn<-n;bt<-t}}
    d<-a[a$padj<0.05 & abs(col)>bt & !is.na(human_symbol),human_symbol]
    data.table(targetN=targetN, scale=sc, thr=bt, nDEG=bn, recall_known=sym_in(d,panels$known_union),
               recall_govaere=sym_in(d,panels$govaere25), median_AveExpr=round(median(a[a$padj<0.05&abs(col)>bt,AveExpr],na.rm=TRUE),2)) }))}
matched <- rbindlist(lapply(c(1000,1500,1885,2500), matchN))
fwrite(matched, file.path(OUT,"matched_N.csv"))
cat("\n=== Recall at matched DEG count (raw vs shrunk) ===\n"); print(matched)

# --- figure ---
p1 <- ggplot(sweep[grepl("raw|shrunk",scale)], aes(nDEG, recall_known, color=scale)) + geom_line() + geom_point(size=1) +
  geom_vline(xintercept=1885, linetype="dashed", color="grey50") +
  labs(title="Known-MASLD recall vs DEG count", x="# DEGs", y="recall known union (%)") + theme_bw(base_size=10)
p2 <- ggplot(cv, aes(threshold, cv, color=scale)) + geom_line() + geom_hline(yintercept=10, linetype="dashed") +
  labs(title="CV-stability plateau per scale", x="|LFC| threshold", y="CV of DEG count (%)") + theme_bw(base_size=10)
p3 <- ggplot(sweep[grepl("raw|shrunk",scale)], aes(nDEG, median_AveExpr, color=scale)) + geom_line() +
  labs(title="Low-expression composition vs DEG count", x="# DEGs", y="median AveExpr of set") + theme_bw(base_size=10)
ggsave(file.path(OUT,"cutoff_benchmark.pdf"), (p1|p2)/(p3|patchwork::plot_spacer()),
       width=11, height=8, useDingbats=FALSE)

# --- per-scale STANDALONE panels (unambiguous on-disk raw vs shrunk set) ---
SCALE_TITLE <- c(raw="raw |log2FC|", shrunk="ashr-shrunk |log2FC|")
for(sc in c("raw","shrunk")){
  lab <- SCALE_TITLE[[sc]]; scl <- paste0(sc,"_absLFC")
  sw <- sweep[scale==scl]; cvs <- cv[scale==sc]
  q1 <- ggplot(sw, aes(nDEG, recall_known)) + geom_line() + geom_point(size=1) +
    geom_vline(xintercept=1853, linetype="dashed", color="grey50") +
    labs(title=sprintf("Known-MASLD recall vs DEG count (%s)", lab),
         x="# DEGs", y="recall known union (%)") + theme_bw(base_size=10)
  q2 <- ggplot(cvs, aes(threshold, cv)) + geom_line() + geom_point(size=1) +
    geom_hline(yintercept=10, linetype="dashed") +
    labs(title=sprintf("CV-stability plateau (%s)", lab),
         x="|LFC| threshold", y="CV of DEG count (%)") + theme_bw(base_size=10)
  q3 <- ggplot(sw, aes(nDEG, median_AveExpr)) + geom_line() +
    labs(title=sprintf("Low-expression composition (%s)", lab),
         x="# DEGs", y="median AveExpr of set") + theme_bw(base_size=10)
  ggsave(file.path(OUT, sprintf("cutoff_benchmark_%s.pdf", sc)),
         (q1|q2)/(q3|patchwork::plot_spacer()), width=11, height=8, useDingbats=FALSE)
  cat(sprintf("[done] wrote cutoff_benchmark_%s.pdf\n", sc))
}
cat("\n[done] outputs in", OUT, "\n")
