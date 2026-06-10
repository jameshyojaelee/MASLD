#!/usr/bin/env Rscript
# ===========================================================================
# Does ashr shrinkage improve recovery of TRUE MASLD DEGs for limma_voom_qw__C2?
# Applies ashr to the C2 (logFC, SE), then benchmarks RAW vs ashr DEG sets
# against four truth proxies: curated positive controls (best proxy for
# "true MASLD gene"), OpenTargets MASLD enrichment, directional mouse
# replication, and COLOC genetic enrichment. Reports at each method's native
# threshold AND at matched DEG-set size (isolates ranking quality from the
# threshold/size effect). Read-only; writes ashr_vs_raw_truth.csv.
# ===========================================================================
suppressMessages({library(data.table); library(ashr)})
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
DEGX <- file.path(BASE, "RNA-seq/results/degx_factorial")
GENEMETA <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
OT_FILE  <- file.path(BASE, "data/published_gene_panels/opentargets_masld_2025.tsv")
ATLAS    <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
PC_FILE  <- file.path(BASE, "RNA-seq/results/validation/positive_control_validation.csv")
strip_v <- function(x) sub("[.][0-9]+$", "", x)

# ---- C2 cell + ashr ----
d <- fread(file.path(DEGX, "cells/cell_limma_voom_qw__C2__kna.csv"))
d[, gb := strip_v(gene)]
d <- d[!is.na(logFC) & !is.na(SE) & SE > 0]
a <- ashr::ash(d$logFC, d$SE, mixcompdist = "normal")
d[, shrunk_logFC := a$result$PosteriorMean]
d[, lfsr := a$result$lfsr]
cat(sprintf("[ashr] %d genes; mean |logFC| %.3f -> mean |shrunk| %.3f (shrinkage %.1f%%)\n",
            nrow(d), mean(abs(d$logFC)), mean(abs(d$shrunk_logFC)),
            100*(1 - mean(abs(d$shrunk_logFC))/mean(abs(d$logFC)))))

# ---- gene map + truth sources ----
gm <- fread(GENEMETA); gm[, eb := strip_v(gene_id)]
sym2eb <- gm[gene_name != "", .(gene_name, eb)]; eb2sym <- gm[, .(eb, gene_name)]
map_sym <- function(s) sym2eb[match(s, gene_name), eb]
pc <- fread(PC_FILE); pc_eb <- unique(na.omit(map_sym(unique(pc$gene))))
ot <- fread(OT_FILE, skip = "gene_symbol"); ot_eb <- unique(na.omit(map_sym(unique(ot$gene_symbol))))
atlas <- fread(ATLAS, select=c("ensembl_id","mouse_meta_padj","mouse_meta_logFC","coloc_best_susie_pp4_polyfun"))
atlas[, gb := strip_v(ensembl_id)]; atlas <- atlas[!is.na(gb)][!duplicated(gb)]
universe   <- intersect(d$gb, atlas$gb)
pc_u    <- intersect(pc_eb, universe)
ot_u    <- intersect(ot_eb, universe)
coloc_u <- atlas[!is.na(coloc_best_susie_pp4_polyfun) & coloc_best_susie_pp4_polyfun > 0.5, gb]
mouse_sig <- atlas[!is.na(mouse_meta_padj) & mouse_meta_padj < 0.05]
cat(sprintf("[truth] universe=%d  PC=%d  OT=%d  COLOC=%d  mouse-tested=%d\n",
            length(universe), length(pc_u), length(ot_u), length(coloc_u), atlas[!is.na(mouse_meta_padj), .N]))

fisher_or <- function(hit, truth, univ){ hit<-intersect(hit,univ); truth<-intersect(truth,univ)
  A<-length(intersect(hit,truth)); B<-length(setdiff(hit,truth)); C<-length(setdiff(truth,hit)); D<-length(univ)-A-B-C
  ft<-tryCatch(fisher.test(matrix(c(A,B,C,max(D,0)),2)),error=function(e)NULL); if(is.null(ft)) NA_real_ else unname(ft$estimate) }
mouse_dir_prec <- function(degset){ m<-merge(data.table(gb=degset), d[,.(gb,hlfc=logFC)],by="gb")
  m<-merge(m, atlas[,.(gb,mp=mouse_meta_padj,ml=mouse_meta_logFC)],by="gb")[!is.na(mp)]
  if(!nrow(m)) return(NA_real_); mean(m$mp<0.05 & sign(m$hlfc)==sign(m$ml)) }

metrics <- function(degset, up){
  degset<-intersect(degset,universe); up<-intersect(up,universe)
  pc_hit<-intersect(degset,pc_u)
  data.table(n_deg=length(degset),
    pc_recall=length(pc_hit)/length(pc_u), pc_precision=length(pc_hit)/max(length(degset),1),
    pc_n=length(pc_hit),
    ot_or_up=fisher_or(up, ot_u, universe),
    mouse_dir=mouse_dir_prec(degset),
    coloc_or=fisher_or(degset, coloc_u, universe)) }

# ---- native-threshold DEG sets ----
raw_deg <- d[padj<0.05 & abs(logFC)>0.5, gb];  raw_up <- d[padj<0.05 & logFC>0.5, gb]
ash_deg <- d[lfsr<0.05 & abs(shrunk_logFC)>0.5, gb]; ash_up <- d[lfsr<0.05 & shrunk_logFC>0.5, gb]
res <- rbind(
  cbind(method="RAW  (padj<.05,|lfc|>.5)",  metrics(raw_deg, raw_up)),
  cbind(method="ashr (lfsr<.05,|slfc|>.5)", metrics(ash_deg, ash_up)))

# ---- matched-N: rank by raw p vs by lfsr, top-N, same truth ----
matchedN <- function(N){
  setorder(d, pval);  r_top<-d[seq_len(min(N,.N))]; r_up<-r_top[logFC>0,gb]
  setorder(d, lfsr);  a_top<-d[seq_len(min(N,.N))]; a_up<-a_top[shrunk_logFC>0,gb]
  rbind(cbind(method=sprintf("RAW  top-%d (rank by p)",N),    metrics(r_top$gb, r_up)),
        cbind(method=sprintf("ashr top-%d (rank by lfsr)",N), metrics(a_top$gb, a_up))) }
resN <- rbindlist(lapply(c(500,1000,1774), matchedN))

out <- rbind(res, resN, fill=TRUE)
fwrite(out, file.path(DEGX, "ashr_vs_raw_truth.csv"))
cat("\n=== NATIVE THRESHOLD ===\n")
print(res[,.(method, n_deg, pc_n, pc_recall=round(pc_recall,3), pc_prec=round(pc_precision,4),
  ot_or_up=round(ot_or_up,3), mouse_dir=round(mouse_dir,3), coloc_or=round(coloc_or,3))], row.names=FALSE)
cat("\n=== MATCHED-N (isolates ranking quality from threshold) ===\n")
print(resN[,.(method, n_deg, pc_n, pc_recall=round(pc_recall,3),
  ot_or_up=round(ot_or_up,3), mouse_dir=round(mouse_dir,3), coloc_or=round(coloc_or,3))], row.names=FALSE)
cat("\nwrote ashr_vs_raw_truth.csv\n")
