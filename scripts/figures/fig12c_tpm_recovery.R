#!/usr/bin/env Rscript
# =============================================================================
# Figure 12c -- Known-MASLD-gene recovery vs mouse disease-liver TPM floor
# Companion to fig12 (which sweeps the LFC cutoff). Here the human ashr spine is
# FIXED at shrunk_logFC > 0.2 (lfsr<0.05, up-only) and we sweep the mouse
# disease-liver TPM gate, to show how aggressively the TPM floor prunes recovery
# of curated / OpenTargets MASLD genes.
#
# Two lines per panel:
#   - Human spine (ashr>0.2, up)                 : panel gene is an ashr>0.2 up-DEG
#   - Full library (+ mouse-confirmed + GWAS)    : panel gene's mouse ortholog is
#                                                  in spine ∪ mouse-conf ∪ GWAS tier
# Recovery requires a mouse ortholog whose pooled disease-liver TPM >= floor.
# Dashed line = canonical PC gate (TPM = 1.0).
#
# Output: 12c_pc_tpm_recovery_curve.pdf , 12c_ot_tpm_recovery_curve.pdf
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE,"scripts/figures/publication_theme.R"))
source(file.path(BASE,"scripts/figures/load_figure_data.R"))
OUT_DIR <- FIGS_CAS13LIB_DIR
strip_v <- function(x) sub("[.][0-9]+$","",x)

INT   <- file.path(BASE,"RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
ORTHO <- file.path(BASE,"data/external/orthologs/master_ortholog_table.tsv.gz")
MOUSE_META <- file.path(BASE,"Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv")
PERDIET <- file.path(BASE,"RNA-seq/Mouse/Unified_Integration/results/per_diet_cas13")  # Cas13 library Western pool; decoupled from paper 4-model per_diet (2026-06-16)
FC_DIR <- file.path(BASE,"RNA-seq/Mouse/Unified_Integration/counts/featurecounts")
PUB_FC <- file.path(BASE,"RNA-seq/Mouse/Public_Diet_Models/counts/featurecounts/gene_counts.txt")
WD_DIR <- file.path(BASE,"RNA-seq/Mouse/Western_Diet_Datasets")
MAIN_META <- file.path(BASE,"RNA-seq/Mouse/Unified_Integration/metadata/unified_mouse_metadata.csv")
CAND_PC  <- file.path(BASE,"Cas13_Library_Design/data/candidates_pc_independent.csv")
CAND_LNC <- file.path(BASE,"Cas13_Library_Design/data/candidates_lncrna_independent.csv")
DIETS <- c("MCD","CDAHFD","Western","HFD"); LFSR <- 0.05; SM <- 0.5; H_THR <- 0.1
EXCLUDE_DATASETS <- c("GSE159911","GSE225616","GSE263273")

# ---- human ashr + panels ---------------------------------------------------
d <- fread(file.path(INT,"meta_results_ashr.csv")); d[, hb := strip_v(gene)]  # v6: metafor ashr (was dream)
pc_sym <- fread(file.path(BASE,"results/library/positive_control.csv"))[["Gene symbol"]]
ot_sym <- unique(fread(file.path(BASE,"data/published_gene_panels/opentargets_masld_2025.tsv"), skip="gene_symbol")$gene_symbol)
d[, is_pc := symbol %in% pc_sym]; d[, is_ot := symbol %in% ot_sym]

# ---- ortholog: best mouse per human ----------------------------------------
ortho <- fread(cmd=paste0("zcat ",ORTHO), select=c("mouse_ensembl","human_ensembl","confidence_tier","is_one2one"))
ortho[, mouse_ensembl := strip_v(mouse_ensembl)]; ortho[, human_ensembl := strip_v(human_ensembl)]
ortho <- ortho[confidence_tier %in% c("H","M")]; ortho[, trank := match(confidence_tier,c("H","M"))]
ortho[, o2o := fifelse(is_one2one %in% c(TRUE,"True","TRUE","true"),0L,1L)]
setorder(ortho, human_ensembl, trank, o2o, mouse_ensembl); ortho_byhuman <- unique(ortho, by="human_ensembl")
setorder(ortho, mouse_ensembl, trank, o2o, human_ensembl); ortho_m2h <- unique(ortho, by="mouse_ensembl")
d <- merge(d, ortho_byhuman[, .(hb=human_ensembl, mgb=mouse_ensembl)], by="hb", all.x=TRUE)

mm <- fread(MOUSE_META)
mm[grepl("protein_coding",mouse_biotype), bt2:="protein_coding"]; mm[grepl("lncRNA|lincRNA",mouse_biotype), bt2:="lncRNA"]
pc_lnc_ids <- mm[bt2 %in% c("protein_coding","lncRNA"), mouse_ensembl_base]

# ---- mouse-confirmed + GWAS tiers (mouse space) ----------------------------
human_logfc <- d[, .(hb, hlfc=logFC)]
up_by_diet <- lapply(DIETS, function(dd){ x<-fread(file.path(PERDIET,paste0(dd,"_de_results.csv")),select=c("gene","shrunk_logFC","lfsr")); x[, gb:=strip_v(gene)]; x[lfsr<LFSR & shrunk_logFC>SM, gb] })
mouse_any<-unique(unlist(up_by_diet)); n_up<-sapply(mouse_any,function(g) sum(sapply(up_by_diet,function(u) g%in%u)))
mc_map<-merge(ortho_m2h[mouse_ensembl %in% mouse_any[n_up>=3], .(mouse_ensembl,hb=human_ensembl)], human_logfc, by="hb", all.x=TRUE)
mouse_conf<-intersect(mc_map[!is.na(hlfc)&hlfc>0,mouse_ensembl], pc_lnc_ids)
.gpc<-fread(CAND_PC)[axis_finemap==1]; .glnc<-fread(CAND_LNC)[axis_finemap==1]
gwas_finemap<-intersect(unique(strip_v(c(.gpc$gene_id_mouse,.glnc$gene_id_mouse))), pc_lnc_ids)

# ---- pooled mouse disease-liver TPM (261 samples) --------------------------
load_fc<-function(p){dt<-fread(p,skip="Geneid");sc<-setdiff(names(dt),c("Geneid","Chr","Start","End","Strand","Length"));len<-setNames(dt[["Length"]],strip_v(dt[["Geneid"]]));mat<-as.matrix(dt[,..sc]);rownames(mat)<-strip_v(dt[["Geneid"]]);colnames(mat)<-basename(dirname(sc));list(counts=mat,lengths=len)}
tpm_of<-function(c,l){g<-intersect(rownames(c),names(l));c<-c[g,,drop=FALSE];rpk<-sweep(c,1,l[g]/1e3,"/");sweep(rpk,2,colSums(rpk)/1e6,"/")}
meta<-fread(MAIN_META); dz<-meta[group_binary!="Control"]; dz<-dz[!dataset%in%EXCLUDE_DATASETS]; dz_main<-dz$sample_id
mcd<-lapply(file.path(FC_DIR,c("gene_counts_inhouse.txt","gene_counts_gse156918.txt","gene_counts_gse205974.txt")),load_fc)
ref<-mcd[[1]]$lengths; mcd_t<-tpm_of(do.call(cbind,lapply(mcd,`[[`,"counts")),ref); pub_t<-tpm_of(load_fc(PUB_FC)$counts,ref)
g220<-tpm_of(load_fc(file.path(WD_DIR,"GSE220575/counts/featurecounts/gene_counts.txt"))$counts,ref)
g246<-tpm_of(load_fc(file.path(WD_DIR,"GSE246088/counts/featurecounts/gene_counts.txt"))$counts,ref)
g305<-tpm_of(load_fc(file.path(WD_DIR,"GSE305484/counts/featurecounts/gene_counts.txt"))$counts,ref)
g328<-tpm_of(load_fc(file.path(WD_DIR,"GSE246328/counts/featurecounts/gene_counts.txt"))$counts,ref)
g565<-tpm_of(load_fc(file.path(WD_DIR,"GSE292565/counts/featurecounts/gene_counts.txt"))$counts,ref)
g220d<-merge(fread(file.path(WD_DIR,"GSE220575/metadata/sample_metadata.csv")),fread(file.path(WD_DIR,"GSE220575/metadata/gsm_to_srr.tsv")),by="gsm")[condition%in%c("MASH","HCC"),srr]
g246mm<-merge(fread(file.path(WD_DIR,"GSE246088/metadata/sample_metadata.csv")),fread(file.path(WD_DIR,"GSE246088/metadata/gsm_to_srr.tsv")),by.x="geo_accession",by.y="gsm")
g246d<-g246mm[genotype=="Plvap_Control"&diet%in%c("Western_Diet","High_Fat_Diet"),srr]
g305d<-merge(fread(file.path(WD_DIR,"GSE305484/metadata/sample_metadata.csv")),fread(file.path(WD_DIR,"GSE305484/metadata/gsm_to_srr.tsv")),by.x="gsm_accession",by.y="gsm")[condition!="Chow",srr]
g328d<-fread(file.path(WD_DIR,"GSE246328/metadata/sample_metadata.csv"))[condition=="Disease",srr]
g565d<-fread(file.path(WD_DIR,"GSE292565/metadata/sample_metadata.csv"))[condition=="Disease",srr]
src<-list(mcd_t[,intersect(colnames(mcd_t),dz_main),drop=FALSE],pub_t[,intersect(colnames(pub_t),dz_main),drop=FALSE],
  g220[,intersect(colnames(g220),g220d),drop=FALSE],g246[,intersect(colnames(g246),g246d),drop=FALSE],
  g305[,intersect(colnames(g305),g305d),drop=FALSE],g328[,intersect(colnames(g328),g328d),drop=FALSE],g565[,intersect(colnames(g565),g565d),drop=FALSE])
allg<-Reduce(intersect,lapply(src,rownames)); pooled<-do.call(cbind,lapply(src,function(m) m[allg,,drop=FALSE]))
tpm_vec<-setNames(apply(pooled,1,median),allg)
cat(sprintf("Pooled disease samples: %d\n", ncol(pooled)))

# ---- per-gene flags --------------------------------------------------------
d[, in_spine := !is.na(lfsr) & lfsr < LFSR & shrunk_logFC > H_THR]
d[, in_lib_tier := in_spine | (mgb %in% mouse_conf) | (mgb %in% gwas_finemap)]
d[, mtpm := tpm_vec[mgb]]

# ---- recovery vs TPM floor -------------------------------------------------
TPM_CUTS <- c(0, 0.25, 0.5, 1, 1.5, 2, 3, 5, 10)
recovery <- function(flagcol, panelcol, Npanel) {
  sub <- d[get(panelcol) == TRUE]
  sapply(TPM_CUTS, function(x) sum(sub[[flagcol]] & !is.na(sub$mtpm) & sub$mtpm >= x) / Npanel)
}
N_PC <- length(pc_sym); N_OT <- length(unique(d[is_ot==TRUE, symbol]))
LINE_LABELS <- c(spine = sprintf("Human DEGs (ashr>%.1f, up)", H_THR), lib = "Full library (+ mouse-conf + GWAS)")
LINE_COLORS <- setNames(c(palette2[1], masld_colors$up), LINE_LABELS)  # teal, magenta

mk <- function(panelcol, Npanel, ylab, title_str) {
  long <- rbindlist(list(
    data.table(tpm=TPM_CUTS, frac=recovery("in_spine",   panelcol, Npanel), line=LINE_LABELS[1]),
    data.table(tpm=TPM_CUTS, frac=recovery("in_lib_tier",panelcol, Npanel), line=LINE_LABELS[2])))
  long[, line := factor(line, levels=LINE_LABELS)]
  ggplot(long, aes(tpm, frac, color=line)) +
    geom_vline(xintercept=1.0, linetype="dashed", color="grey45", linewidth=0.4) +
    # draw full-library thick underneath, then Human-DEGs thinner on top so both
    # stay visible where the curves coincide (no data offset)
    geom_line(data=long[line==LINE_LABELS[2]], linewidth=1.5) +
    geom_point(data=long[line==LINE_LABELS[2]], size=3.2) +
    geom_line(data=long[line==LINE_LABELS[1]], linewidth=0.7) +
    geom_point(data=long[line==LINE_LABELS[1]], size=1.7) +
    scale_color_manual(values=LINE_COLORS, name=NULL) +
    scale_x_continuous(name="Mouse disease-liver TPM floor", breaks=TPM_CUTS,
                       guide=guide_axis(n.dodge=2)) +
    scale_y_continuous(name=ylab, labels=scales::percent_format(accuracy=1),
                       limits=c(0,NA), expand=expansion(mult=c(0,0.12))) +
    labs(title=title_str) +
    theme_masld(base_size=11) + theme(legend.position="top", legend.text=element_text(size=9))
}

pdf_dev <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
pB <- mk("is_pc", N_PC, "Fraction of curated panel recovered",
         sprintf("Curated MASLD panel (n=%d) — recovery vs TPM floor (ashr>%.1f fixed)", N_PC, H_THR))
ggsave(file.path(OUT_DIR,"12c_pc_tpm_recovery_curve.pdf"), pB, width=6.8, height=5, device=pdf_dev)
pC <- mk("is_ot", N_OT, "Fraction of OpenTargets MASLD genes recovered",
         sprintf("OpenTargets MASLD (n=%d) — recovery vs TPM floor (ashr>%.1f fixed)", N_OT, H_THR))
ggsave(file.path(OUT_DIR,"12c_ot_tpm_recovery_curve.pdf"), pC, width=6.8, height=5, device=pdf_dev)
cat("Wrote 12c_pc_tpm_recovery_curve.pdf , 12c_ot_tpm_recovery_curve.pdf\n")

# summary at canonical PC gate
sp1 <- recovery("in_spine","is_pc",N_PC); lb1 <- recovery("in_lib_tier","is_pc",N_PC)
cat(sprintf("\nCurated panel @TPM>=1: spine=%.0f%%  full-lib=%.0f%%  (@TPM>=0: spine=%.0f%% full=%.0f%%)\n",
  100*sp1[TPM_CUTS==1],100*lb1[TPM_CUTS==1],100*sp1[TPM_CUTS==0],100*lb1[TPM_CUTS==0]))
