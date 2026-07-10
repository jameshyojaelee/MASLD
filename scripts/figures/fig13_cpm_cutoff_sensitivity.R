#!/usr/bin/env Rscript
# =============================================================================
# Figure 13 -- Library sensitivity plots (final gene count after ALL gates)
# Each panel sweeps one parameter while holding all others at canonical values.
#
# Canonical: cohort>=2, TREAT FDR<0.05, LFC>=0.2 (PCG), padj<0.05 (lncRNA),
#            mouse-hep CPM>=1 (PCG), CPM>=0.3 (lncRNA), guideable in vM38.
#
# Output:
#   13a_lncrna_count_vs_cpm_cutoff.pdf   -- sweep lncRNA CPM (canonical marked)
#   13b_pcg_count_vs_cpm_cutoff.pdf      -- sweep PCG CPM (canonical marked)
#   13c_pcg_count_vs_lfc_cutoff.pdf      -- sweep PCG LFC / TREAT threshold
#   13d_count_vs_padj_cutoff.pdf         -- sweep padj (PCG TREAT FDR + lncRNA raw FDR)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT_DIR <- FIGS_CAS13LIB_DIR

# ---- Canonical parameter values (marked on every plot) ----------------------
CANON_CPM_PC   <- 1.0
CANON_CPM_LNC  <- 0.1
CANON_LFC_PC   <- 0.2
CANON_PADJ     <- 0.05
CANON_COHORT   <- 2L

# ---- Source paths -----------------------------------------------------------
PERDIET   <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet_cas13")
PERSTUDY  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study")
ORTHO     <- file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz")
META      <- file.path(BASE, "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv")
MOUSEHEP  <- file.path(BASE, "Cas13_Library_Design/data/mouse_hep_specificity_vm38.csv")
CANONICAL <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
COLOCFILE <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
POSCTRL   <- file.path(BASE, "results/library/positive_control.csv")
UNGUID    <- file.path(BASE, "Cas13_Library_Design/data/guides/unguideable_vM38.csv")

COHORTS      <- c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621")
DIETS        <- c("MCD","CDAHFD","Western","HFD")
MIN_DIETS    <- 3L
LFSR_THR     <- 0.05;  SHRUNK_LFC_THR_LNC <- 0.0
# (PC effect-size floor = CANON_LFC_PC = 0.2, applied via TREAT; no standalone 0.5 cutoff.)
COLOC_PP4_THR <- 0.5;  COLOC_GATE_EXEMPT <- 0.9
POSCTRL_ALIASES       <- c("SCD1" = "SCD")
POSCTRL_EXCLUDE       <- c("GIPR","GLP1R","FAP")
POSCTRL_MOUSE_OVERRIDE <- c("SCD" = "Scd1")
KEEP_BIOTYPES <- c("protein_coding","lncRNA")
strip_v <- function(x) sub("[.][0-9]+$","",x)

# ---- Analytical TREAT FDR (mirror rebuild_cas13_library.R add_treat_fdr*) ----
# Two-sided shifted-t tail at logFC offset `lfc`, then BH. BH is applied SEPARATELY
# within each biotype (PC vs lncRNA) -- exactly as the canonical rebuild does -- so the
# advertised PC LFC floor binds on protein-coding genes and the BH universe matches.
.treat_fdr <- function(logfc, se, df_use, lfc) {
  se2 <- se; se2[!is.finite(se2) | se2 <= 0] <- NA_real_
  p <- pt((abs(logfc) - lfc) / se2, df = df_use, lower.tail = FALSE) +
       pt((abs(logfc) + lfc) / se2, df = df_use, lower.tail = FALSE)
  p.adjust(p, method = "BH")
}
# Biotype-split TREAT FDR on a data.table carrying logFC, bt, SE/df (or t for SE).
# Returns the input dt with an fdr_treat column (PC at lfc_pc, lncRNA at lfc_lnc=0).
add_treat_fdr_bt <- function(dt, se, df_use, lfc_pc, lfc_lnc = 0.0) {
  out <- copy(dt)
  out[, fdr_treat := NA_real_]
  is_pc  <- out$bt != "lncRNA"
  is_lnc <- out$bt == "lncRNA"
  if (any(is_pc))
    out[is_pc,  fdr_treat := .treat_fdr(logFC, se[is_pc],  df_use[is_pc],  lfc_pc)]
  if (any(is_lnc))
    out[is_lnc, fdr_treat := .treat_fdr(logFC, se[is_lnc], df_use[is_lnc], lfc_lnc)]
  out
}

# =============================================================================
# LOAD ALL STATIC DATA ONCE
# =============================================================================
cat("Loading metadata and orthologs...\n")
meta_full <- fread(META)
setnames(meta_full,
  c("mouse_ensembl_base","mouse_biotype","mouse_symbol_gtf"),
  c("gene_id_mouse","biotype","gene_symbol_mouse"), skip_absent=TRUE)
meta_full[grepl("protein_coding",biotype), biotype := "protein_coding"]
meta_full[grepl("lncRNA|lincRNA",biotype), biotype := "lncRNA"]
mouse_biotype_of <- setNames(meta_full$biotype, meta_full$gene_id_mouse)
pc_ids  <- meta_full[biotype=="protein_coding", gene_id_mouse]
lnc_ids <- meta_full[biotype=="lncRNA",         gene_id_mouse]
pclnc   <- c(pc_ids, lnc_ids)

ortho_raw <- fread(cmd=paste0("zcat ",ORTHO),
  select=c("mouse_ensembl","human_ensembl","human_symbol","confidence_tier","is_one2one"))
ortho_raw[, gene_id_mouse:=strip_v(mouse_ensembl)]; ortho_raw[, human_ensembl:=strip_v(human_ensembl)]
ortho_raw <- ortho_raw[confidence_tier %in% c("H","M")]
ortho_raw[, trank:=match(confidence_tier,c("H","M"))]
ortho_raw[, o2o:=is_one2one %in% c(TRUE,"True","TRUE","true")]
ortho_raw[, one2one_rank:=ifelse(o2o,0L,1L)]   # one2one-preferred (rebuild parity); o2o asc would pick NON-one2one
oh <- copy(ortho_raw); setorder(oh,human_ensembl,trank,one2one_rank,gene_id_mouse); oh<-unique(oh,by="human_ensembl")
om <- copy(ortho_raw); setorder(om,gene_id_mouse,trank,one2one_rank,human_ensembl); om<-unique(om,by="gene_id_mouse")
osym <- copy(ortho_raw); osym[,hsym:=toupper(human_symbol)]; osym<-osym[hsym!=""]
setorder(osym,hsym,trank,one2one_rank,gene_id_mouse); osym<-unique(osym,by="hsym")
h2m  <- function(hv) intersect(unique(oh[human_ensembl %in% hv, gene_id_mouse]), pclnc)
sym2m<- function(sv) intersect(unique(osym[hsym %in% toupper(sv), gene_id_mouse]), pclnc)
human_biotype_of <- setNames(mouse_biotype_of[oh$gene_id_mouse], oh$human_ensembl)

# Canonical DEG file + df_total for TREAT
cat("Loading canonical DEGs...\n")
can <- fread(CANONICAL, select=c("gene","logFC","SE","t","P.Value","padj","symbol"))
can[, hb:=strip_v(gene)]
can_m <- merge(can, oh[,.(human_ensembl,gene_id_mouse)], by.x="hb", by.y="human_ensembl", all.x=TRUE)
can_m[, bt:=mouse_biotype_of[gene_id_mouse]]; can_m[is.na(bt), bt:="protein_coding"]
pr <- can[is.finite(t)&is.finite(P.Value)&P.Value>0&P.Value<1&abs(t)>1e-6]
pi2 <- unique(round(seq(1,nrow(pr),length.out=min(nrow(pr),12))))
df_total_canon <- median(vapply(pi2, function(i)
  uniroot(function(df) 2*pt(-abs(pr$t[i]),df=df)-pr$P.Value[i],c(0.1,1e6))$root, numeric(1)))
cat(sprintf("  canonical df_total = %.1f\n", df_total_canon))

# Canonical human CORE (C1) human-ENSG set at lfc=CANON_LFC_PC -- used to gate the
# mouse-confirmed tier's human directional concordance (rebuild restricts human_logfc to
# hb %in% core_h, NOT any positive human logFC). Recomputed inside build_human_tiers for
# the sweeps; this fixed copy only feeds the mouse-confirmed concordance map below.
.canon_core_dt <- add_treat_fdr_bt(can_m, se=can_m$SE,
                                   df_use=rep(df_total_canon, nrow(can_m)),
                                   lfc_pc=CANON_LFC_PC, lfc_lnc=SHRUNK_LFC_THR_LNC)
core_h_canon <- unique(.canon_core_dt[!is.na(logFC) & !is.na(fdr_treat) &
                                      fdr_treat < CANON_PADJ & logFC > 0, hb])

# Per-study DEG files + per-cohort df_total
cat("Loading per-study cohort DEGs...\n")
perstudy <- lapply(COHORTS, function(coh) {
  d <- fread(file.path(PERSTUDY,paste0(coh,"_de_results.csv")))
  d[, hb:=strip_v(gene)]
  dm <- merge(d, oh[,.(human_ensembl,gene_id_mouse)], by.x="hb",by.y="human_ensembl",all.x=TRUE)
  dm[, bt:=mouse_biotype_of[gene_id_mouse]]; dm[is.na(bt), bt:="protein_coding"]
  dft <- dm$df.total[1]
  if (is.null(dft)||is.na(dft)) {
    pr2 <- dm[is.finite(t)&is.finite(P.Value)&P.Value>0&P.Value<1&abs(t)>1e-6]
    pi3 <- unique(round(seq(1,nrow(pr2),length.out=min(nrow(pr2),12))))
    dft <- median(vapply(pi3,function(i)
      uniroot(function(df)2*pt(-abs(pr2$t[i]),df=df)-pr2$P.Value[i],c(0.1,1e6))$root,numeric(1)))
  }
  dm[, df_study:=dft]
  dm
})
names(perstudy) <- COHORTS

# Mouse cross-diet (fixed)
cat("Loading mouse cross-diet DEGs...\n")
# Mouse cross-diet membership = TREAT FDR per diet (PC lfc=CANON_LFC_PC, lncRNA lfc=0),
# UP only -- identical to the canonical rebuild's add_treat_fdr mouse arm (se=|logFC/t|,
# df=df_total column). Replaces the retired ashr lfsr/shrunk_logFC point-estimate gate.
up_by_diet <- lapply(DIETS, function(d) {
  dt <- fread(file.path(PERDIET,paste0(d,"_de_results.csv")))
  dt[, gb:=strip_v(gene)]; dt[, bt:=mouse_biotype_of[gb]]; dt[is.na(bt),bt:="protein_coding"]
  se_d <- abs(dt$logFC / dt$t)
  df_d <- if ("df_total" %in% names(dt)) dt$df_total else rep(df_total_canon, nrow(dt))
  dt <- add_treat_fdr_bt(dt, se=se_d, df_use=df_d, lfc_pc=CANON_LFC_PC, lfc_lnc=SHRUNK_LFC_THR_LNC)
  dt[!is.na(fdr_treat) & fdr_treat < LFSR_THR & logFC > 0, gb]
})
names(up_by_diet)<-DIETS; all_up<-sort(unique(unlist(up_by_diet)))
cross <- data.table(gene_id_mouse=all_up)
for (d in DIETS) cross[,(d):=gene_id_mouse %in% up_by_diet[[d]]]
cross[,n_diets_up:=rowSums(.SD),.SDcols=DIETS]
# MOUSE-CONFIRMED requires the mouse gene's best-human ortholog to be a CORE (C1) DEG --
# directional concordance against the human disease-vs-control CORE set (core_h_canon),
# exactly as the rebuild does (human_logfc restricted to hb %in% core_h). The prior
# "any positive human logFC" filter admitted ~470 net-new PC passengers.
m2h_map <- merge(data.table(gene_id_mouse=cross[n_diets_up>=MIN_DIETS,gene_id_mouse]),
                 om[,.(gene_id_mouse,human_ensembl)],by="gene_id_mouse",all.x=TRUE)
mouse_conf_base <- m2h_map[human_ensembl %in% core_h_canon, gene_id_mouse]
mouse_conf_base <- intersect(mouse_conf_base, pclnc)

# COLOC tier (fixed)
cl <- fread(COLOCFILE); cl[,hb:=strip_v(ensembl)]
coloc_mouse <- h2m(unique(cl[hb!=""&!is.na(coloc_best_susie_pp4)&coloc_best_susie_pp4>COLOC_PP4_THR,hb]))
coloc_pp4_dt <- merge(cl[hb!=""&!is.na(coloc_best_susie_pp4),.(human_ensembl=hb,pp4=coloc_best_susie_pp4)],
                      oh[,.(human_ensembl,gene_id_mouse)],by="human_ensembl")
coloc_pp4_dt <- coloc_pp4_dt[,.(coloc_pp4=max(pp4,na.rm=TRUE)),by=gene_id_mouse]
coloc_pp4_of <- setNames(coloc_pp4_dt$coloc_pp4, coloc_pp4_dt$gene_id_mouse)

# Positive controls (fixed)
msym <- meta_full[,.(gene_id_mouse,gene_symbol_mouse)]; msym[,msu:=toupper(gene_symbol_mouse)]
msym_to_id <- function(s){v<-msym[msu==toupper(s),gene_id_mouse];if(length(v))v[1] else NA_character_}
pc_raw <- fread(POSCTRL); pc_raw[,hsym:=toupper(`Gene symbol`)]
pc_raw[hsym %in% names(POSCTRL_ALIASES), hsym:=POSCTRL_ALIASES[hsym]]
pc_raw <- pc_raw[!(hsym %in% POSCTRL_EXCLUDE)]
pc_map <- merge(pc_raw[,.(hsym,pos_control_direction=`Steatosis_Change_upon_KD`)],
                osym[,.(hsym,gene_id_mouse)],by="hsym",all.x=TRUE)
for (hs in names(POSCTRL_MOUSE_OVERRIDE)){
  mid<-msym_to_id(POSCTRL_MOUSE_OVERRIDE[[hs]]); if(!is.na(mid)) pc_map[hsym==hs,gene_id_mouse:=mid]
}
pc_map<-unique(pc_map[!is.na(gene_id_mouse)],by="gene_id_mouse")
posctrl_mouse <- pc_map$gene_id_mouse

# Mouse hepatocyte CPM
cat("Loading mouse hepatocyte CPM...\n")
mhs <- fread(MOUSEHEP)
mhs_id <- unique(mhs[!is.na(gene_id)&gene_id!="",.(gene_id_mouse=gene_id,mouse_hep_cpm)],by="gene_id_mouse")
mhs_sym <- unique(mhs[,.(gene_symbol_mouse=gene_symbol,cpm_s=mouse_hep_cpm)],by="gene_symbol_mouse")
cpm_table <- merge(meta_full[,.(gene_id_mouse,gene_symbol_mouse,biotype)], mhs_id, by="gene_id_mouse",all.x=TRUE)
cpm_table <- merge(cpm_table, mhs_sym, by="gene_symbol_mouse", all.x=TRUE)
cpm_table[is.na(mouse_hep_cpm), mouse_hep_cpm:=cpm_s]; cpm_table[,cpm_s:=NULL]
cpm_table[is.na(mouse_hep_cpm), mouse_hep_cpm:=0]
cpm_of <- setNames(cpm_table$mouse_hep_cpm, cpm_table$gene_id_mouse)

# Guidability
ung_ids <- if(file.exists(UNGUID)) strip_v(fread(UNGUID)$gene_id_mouse) else character(0)
cat(sprintf("Unguideable genes: %d\n", length(ung_ids)))

# =============================================================================
# CORE TIER BUILDER -- parameterised by lfc_pc and padj_thr
# =============================================================================
build_human_tiers <- function(lfc_pc = CANON_LFC_PC, padj_thr = CANON_PADJ,
                              cohort_min = CANON_COHORT) {
  # CORE = C1 integrated disease-vs-control ONLY (canonical_deg_results.csv), biotype-
  # split TREAT FDR (PC at lfc_pc, lncRNA at 0), UP only -- identical to the rebuild.
  can_m_local <- add_treat_fdr_bt(can_m, se=can_m$SE,
                                  df_use=rep(df_total_canon, nrow(can_m)),
                                  lfc_pc=lfc_pc, lfc_lnc=SHRUNK_LFC_THR_LNC)
  core_h <- unique(can_m_local[!is.na(logFC) & !is.na(fdr_treat) &
                               fdr_treat < padj_thr & logFC > 0, hb])
  core_mouse <- h2m(core_h)

  # COHORT-REPLICATED: human DEGs UP (TREAT FDR<padj_thr, biotype-split) in >=N cohorts.
  cohort_lists <- lapply(perstudy, function(dm) {
    dml <- add_treat_fdr_bt(dm, se=dm$SE, df_use=dm$df_study,
                            lfc_pc=lfc_pc, lfc_lnc=SHRUNK_LFC_THR_LNC)
    unique(dml[!is.na(logFC) & !is.na(fdr_treat) & fdr_treat < padj_thr & logFC > 0, hb])
  })
  cohort_n <- table(unlist(cohort_lists))
  cohort_mouse <- h2m(names(cohort_n)[cohort_n >= cohort_min])

  list(core=core_mouse, cohort=cohort_mouse)
}

# =============================================================================
# FINAL LIBRARY COUNTER -- applies CPM + exemptions + guidability
# =============================================================================
count_final <- function(tiers,
                        cpm_pc = CANON_CPM_PC, cpm_lnc = CANON_CPM_LNC) {
  pool_ungated <- sort(Reduce(union, list(
    tiers$core, tiers$cohort, mouse_conf_base, coloc_mouse, posctrl_mouse)))
  pool <- pool_ungated[pool_ungated %in% pclnc]

  cpm  <- cpm_of[pool]; cpm[is.na(cpm)] <- 0
  bt   <- mouse_biotype_of[pool]; bt[is.na(bt)] <- "protein_coding"
  exempt <- (pool %in% posctrl_mouse) |
            (!is.na(coloc_pp4_of[pool]) & coloc_pp4_of[pool] >= COLOC_GATE_EXEMPT)
  pass_cpm <- (bt=="protein_coding" & cpm >= cpm_pc) |
              (bt=="lncRNA"         & cpm >= cpm_lnc) | exempt
  pass_guide <- !(pool %in% ung_ids)
  final <- pool[pass_cpm & pass_guide]
  list(
    total = length(final),
    pc    = sum(mouse_biotype_of[final] == "protein_coding", na.rm=TRUE),
    lnc   = sum(mouse_biotype_of[final] == "lncRNA",         na.rm=TRUE)
  )
}

# =============================================================================
# SWEEP HELPERS
# =============================================================================
canonical_tiers <- build_human_tiers()   # computed once for CPM sweeps
cat(sprintf("Canonical tiers: core=%d  cohort=%d\n",
            length(canonical_tiers$core), length(canonical_tiers$cohort)))

sweep_cpm_lnc <- function(cuts) {
  rbindlist(lapply(cuts, function(x) {
    r <- count_final(canonical_tiers, cpm_lnc=x)
    data.table(cutoff=x, n_lnc=r$lnc, n_pc=r$pc, n_total=r$total)
  }))
}

sweep_cpm_pc <- function(cuts) {
  rbindlist(lapply(cuts, function(x) {
    r <- count_final(canonical_tiers, cpm_pc=x)
    data.table(cutoff=x, n_lnc=r$lnc, n_pc=r$pc, n_total=r$total)
  }))
}

sweep_lfc_pc <- function(cuts) {
  rbindlist(lapply(cuts, function(x) {
    cat(sprintf("  LFC sweep: %.2f\n", x))
    t <- build_human_tiers(lfc_pc=x)
    r <- count_final(t)
    data.table(cutoff=x, n_lnc=r$lnc, n_pc=r$pc, n_total=r$total)
  }))
}

sweep_padj <- function(cuts) {
  rbindlist(lapply(cuts, function(x) {
    cat(sprintf("  padj sweep: %.4f\n", x))
    t <- build_human_tiers(padj_thr=x)
    r <- count_final(t)
    data.table(cutoff=x, n_lnc=r$lnc, n_pc=r$pc, n_total=r$total)
  }))
}

# =============================================================================
# RUN SWEEPS
# =============================================================================
cat("Sweeping lncRNA CPM cutoff...\n")
cuts_cpm_lnc <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
s_lnc <- sweep_cpm_lnc(cuts_cpm_lnc)
print(s_lnc)

cat("Sweeping PCG CPM cutoff...\n")
cuts_cpm_pc <- c(0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0)
s_pcg <- sweep_cpm_pc(cuts_cpm_pc)
print(s_pcg)

cat("Sweeping PCG LFC cutoff (TREAT FDR < 0.05)...\n")
cuts_lfc <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0)
s_lfc <- sweep_lfc_pc(cuts_lfc)
print(s_lfc)

cat("Sweeping padj / TREAT FDR cutoff...\n")
cuts_padj <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2)
s_padj <- sweep_padj(cuts_padj)
print(s_padj)

# =============================================================================
# SHARED THEME HELPERS
# =============================================================================
mark_canonical <- function(x_val, y_max, offset_x=0, label_side="right") {
  list(
    geom_vline(xintercept=x_val, linetype="dashed", color="grey40", linewidth=0.5),
    annotate("text", x=x_val + offset_x, y=y_max,
             label=paste0("canonical\n(", x_val, ")"), hjust=ifelse(label_side=="right",0,1),
             size=GEOM_TEXT_6PT, color="black", vjust=1)
  )
}

pdf_dev <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
col_lnc   <- "#C0143C"
col_pc    <- "#2271B2"
col_total <- "#333333"

# =============================================================================
# 13a: lncRNA count vs CPM cutoff
# =============================================================================
p13a <- ggplot(s_lnc, aes(cutoff, n_lnc)) +
  mark_canonical(CANON_CPM_LNC, max(s_lnc$n_lnc)*0.98, offset_x=0.02) +
  geom_line(linewidth=1, color=col_lnc) +
  geom_point(size=2.5, color=col_lnc) +
  geom_text(data=s_lnc[cutoff==CANON_CPM_LNC],
            aes(label=n_lnc), vjust=-1, hjust=0.5, size=GEOM_TEXT_6PT, fontface="plain", color=col_lnc) +
  scale_x_continuous(breaks=cuts_cpm_lnc) +
  scale_y_continuous(labels=comma, expand=expansion(mult=c(0.05,0.12))) +
  labs(x="Mouse hepatocyte CPM cutoff", y="lncRNA genes in final library") +
  theme_masld(base_size=11) + theme_pub()

ggsave(file.path(OUT_DIR,"13a_lncrna_count_vs_cpm_cutoff.pdf"), p13a,
       width=6, height=4.5, device=pdf_dev)
message("[caption] lncRNA: library size vs CPM cutoff")
cat("Wrote 13a\n")

# =============================================================================
# 13b: PCG count vs CPM cutoff
# =============================================================================
p13b <- ggplot(s_pcg, aes(cutoff, n_pc)) +
  mark_canonical(CANON_CPM_PC, max(s_pcg$n_pc)*0.98, offset_x=0.05) +
  geom_line(linewidth=1, color=col_pc) +
  geom_point(size=2.5, color=col_pc) +
  geom_text(data=s_pcg[cutoff==CANON_CPM_PC],
            aes(label=n_pc), vjust=-1, hjust=0.5, size=GEOM_TEXT_6PT, fontface="plain", color=col_pc) +
  scale_x_continuous(breaks=cuts_cpm_pc) +
  scale_y_continuous(labels=comma, expand=expansion(mult=c(0.05,0.12))) +
  labs(x="Mouse hepatocyte CPM cutoff", y="PCG genes in final library") +
  theme_masld(base_size=11) + theme_pub()

ggsave(file.path(OUT_DIR,"13b_pcg_count_vs_cpm_cutoff.pdf"), p13b,
       width=6, height=4.5, device=pdf_dev)
message("[caption] PCG: library size vs CPM cutoff")
cat("Wrote 13b\n")

# =============================================================================
# 13c: library size vs LFC cutoff  (PCG TREAT, lncRNA unchanged -- also show total)
# =============================================================================
s_lfc_long <- melt(s_lfc[,.(cutoff,PCG=n_pc,lncRNA=n_lnc,Total=n_total)],
                   id.vars="cutoff", variable.name="biotype", value.name="n")
s_lfc_long[, biotype := factor(biotype, levels=c("Total","PCG","lncRNA"))]
lfc_cols <- c(Total=col_total, PCG=col_pc, lncRNA=col_lnc)

p13c <- ggplot(s_lfc_long, aes(cutoff, n, color=biotype, group=biotype)) +
  mark_canonical(CANON_LFC_PC, max(s_lfc$n_total)*0.98, offset_x=0.02) +
  geom_line(linewidth=1) + geom_point(size=2.5) +
  geom_text(data=s_lfc_long[cutoff==CANON_LFC_PC],
            aes(label=n), vjust=-1, hjust=0.5, size=GEOM_TEXT_6PT, fontface="plain") +
  scale_color_manual(values=lfc_cols, name=NULL) +
  scale_x_continuous(breaks=cuts_lfc) +
  scale_y_continuous(labels=comma, expand=expansion(mult=c(0.05,0.12))) +
  labs(x="PCG LFC cutoff (TREAT, padj < 0.05)",
       y="Genes in final library") +
  theme_masld(base_size=11) + theme_pub() +
  theme(legend.position="top")

ggsave(file.path(OUT_DIR,"count_vs_lfc_cutoff.pdf"), p13c,
       width=6.5, height=4.5, device=pdf_dev)
message("[caption] Library size vs LFC cutoff")
cat("Wrote 13c\n")

# =============================================================================
# 13d: library size vs padj / TREAT FDR cutoff (PCG + lncRNA + total)
# =============================================================================
s_padj_long <- melt(s_padj[,.(cutoff,PCG=n_pc,lncRNA=n_lnc,Total=n_total)],
                    id.vars="cutoff", variable.name="biotype", value.name="n")
s_padj_long[, biotype := factor(biotype, levels=c("Total","PCG","lncRNA"))]

p13d <- ggplot(s_padj_long, aes(cutoff, n, color=biotype, group=biotype)) +
  mark_canonical(CANON_PADJ, max(s_padj$n_total)*0.98, offset_x=0.002) +
  geom_line(linewidth=1) + geom_point(size=2.5) +
  geom_text(data=s_padj_long[cutoff==CANON_PADJ],
            aes(label=n), vjust=-1, hjust=0.5, size=GEOM_TEXT_6PT, fontface="plain") +
  scale_color_manual(values=lfc_cols, name=NULL) +
  scale_x_continuous(breaks=cuts_padj,
                     labels=function(x) ifelse(x<0.01, formatC(x,format="e",digits=0),
                                               as.character(x))) +
  scale_y_continuous(labels=comma, expand=expansion(mult=c(0.05,0.12))) +
  labs(x="padj cutoff (PCG: TREAT FDR; lncRNA: raw FDR)",
       y="Genes in final library") +
  theme_masld(base_size=11) + theme_pub() +
  theme(legend.position="top")

ggsave(file.path(OUT_DIR,"count_vs_padj_cutoff.pdf"), p13d,
       width=6.5, height=4.5, device=pdf_dev)
message("[caption] Library size vs significance threshold")
cat("Wrote 13d\n")

cat("fig13_cpm_cutoff_sensitivity.R complete.\n")
