#!/usr/bin/env Rscript
# gwas_overlap_with_library.R
# ---------------------------------------------------------------------------
# How many of the GWAS nearest-genes (finemap tier) does the current library
# already capture via its human-DEG spine and mouse-confirmed tier, at various
# human ashr cutoffs?  And how many are NET-NEW if we add the GWAS tier to
# Option F (ashr>0.20 + split TPM gate, 1,991 genes)?
#
# All in mouse-ortholog space (the library targets mouse). GWAS genes gated by
# the same biotype-aware TPM floor (PC>=1.0 / lncRNA>=0.5) when added, since a
# knockdown target must be expressed.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(data.table) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FC_DIR <- file.path(BASE,"RNA-seq/Mouse/Unified_Integration/counts/featurecounts")
PUB_FC <- file.path(BASE,"RNA-seq/Mouse/Public_Diet_Models/counts/featurecounts/gene_counts.txt")
WD_DIR <- file.path(BASE,"RNA-seq/Mouse/Western_Diet_Datasets")
MAIN_META <- file.path(BASE,"RNA-seq/Mouse/Unified_Integration/metadata/unified_mouse_metadata.csv")
ASHR  <- file.path(BASE,"RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv")
ORTHO <- file.path(BASE,"data/external/orthologs/master_ortholog_table.tsv.gz")
MOUSE_META <- file.path(BASE,"Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv")
PERDIET <- file.path(BASE,"RNA-seq/Mouse/Unified_Integration/results/per_diet")
CAND_PC  <- file.path(BASE,"Cas13_Library_Design/data/candidates_pc_independent.csv")
CAND_LNC <- file.path(BASE,"Cas13_Library_Design/data/candidates_lncrna_independent.csv")
DATD <- file.path(BASE,"Cas13_Library_Design/data")
DIETS <- c("MCD","CDAHFD","Western","HFD"); LFSR <- 0.05; SM <- 0.5
PC_GATE <- 1.0; LNC_GATE <- 0.5
EXCLUDE_DATASETS <- c("GSE159911","GSE225616","GSE263273")
CHOSEN <- 0.20                       # Option F
strip_v <- function(x) sub("[.][0-9]+$","",x)

# ---- selection inputs ------------------------------------------------------
ash <- fread(ASHR, select=c("gene","symbol","shrunk_logFC","lfsr","logFC")); ash[, hb := strip_v(gene)]
human_logfc <- ash[, .(hb, hlfc = logFC)]
ortho <- fread(cmd=paste0("zcat ",ORTHO), select=c("mouse_ensembl","human_ensembl","confidence_tier","is_one2one"))
ortho[, mouse_ensembl := strip_v(mouse_ensembl)]; ortho[, human_ensembl := strip_v(human_ensembl)]
ortho <- ortho[confidence_tier %in% c("H","M")]; ortho[, trank := match(confidence_tier,c("H","M"))]
ortho[, o2o := fifelse(is_one2one %in% c(TRUE,"True","TRUE","true"),0L,1L)]
setorder(ortho, human_ensembl, trank, o2o, mouse_ensembl); ortho_byhuman <- unique(ortho, by="human_ensembl")
setorder(ortho, mouse_ensembl, trank, o2o, human_ensembl); ortho_m2h <- unique(ortho, by="mouse_ensembl")

mm <- fread(MOUSE_META)
mm[grepl("protein_coding", mouse_biotype), bt2 := "protein_coding"]
mm[grepl("lncRNA|lincRNA", mouse_biotype), bt2 := "lncRNA"]
pc_lnc_ids <- mm[bt2 %in% c("protein_coding","lncRNA"), mouse_ensembl_base]
bt_of <- setNames(mm$bt2, mm$mouse_ensembl_base)

up_by_diet <- lapply(DIETS, function(d){ dt<-fread(file.path(PERDIET,paste0(d,"_de_results.csv")),
  select=c("gene","shrunk_logFC","lfsr")); dt[, gb:=strip_v(gene)]; dt[lfsr<LFSR & shrunk_logFC>SM, gb] })
mouse_any <- unique(unlist(up_by_diet)); n_up <- sapply(mouse_any, function(g) sum(sapply(up_by_diet, function(u) g %in% u)))
mouse_3plus <- mouse_any[n_up>=3]
mc_map <- merge(ortho_m2h[mouse_ensembl %in% mouse_3plus, .(mouse_ensembl, hb=human_ensembl)], human_logfc, by="hb", all.x=TRUE)
mouse_conf <- intersect(mc_map[!is.na(hlfc) & hlfc>0, mouse_ensembl], pc_lnc_ids)
get_spine <- function(thr) intersect(ortho_byhuman[human_ensembl %in% ash[!is.na(lfsr)&lfsr<LFSR&shrunk_logFC>thr, hb], mouse_ensembl], pc_lnc_ids)
get_lib   <- function(thr) sort(union(get_spine(thr), mouse_conf))

# ---- pooled disease-liver TPM (261 samples) --------------------------------
load_fc <- function(p){ dt<-fread(p,skip="Geneid"); sc<-setdiff(names(dt),c("Geneid","Chr","Start","End","Strand","Length"))
  len<-setNames(dt[["Length"]],strip_v(dt[["Geneid"]])); mat<-as.matrix(dt[,..sc]); rownames(mat)<-strip_v(dt[["Geneid"]]); colnames(mat)<-basename(dirname(sc)); list(counts=mat,lengths=len) }
tpm_of <- function(c,l){ g<-intersect(rownames(c),names(l)); c<-c[g,,drop=FALSE]; rpk<-sweep(c,1,l[g]/1e3,"/"); sweep(rpk,2,colSums(rpk)/1e6,"/") }
meta<-fread(MAIN_META); dz<-meta[group_binary!="Control"]; dz<-dz[!dataset %in% EXCLUDE_DATASETS]; dz_main<-dz$sample_id
mcd<-lapply(file.path(FC_DIR,c("gene_counts_inhouse.txt","gene_counts_gse156918.txt","gene_counts_gse205974.txt")),load_fc)
ref_len<-mcd[[1]]$lengths; mcd_tpm<-tpm_of(do.call(cbind,lapply(mcd,`[[`,"counts")),ref_len)
pub_tpm<-tpm_of(load_fc(PUB_FC)$counts,ref_len)
g220<-tpm_of(load_fc(file.path(WD_DIR,"GSE220575/counts/featurecounts/gene_counts.txt"))$counts,ref_len)
g246<-tpm_of(load_fc(file.path(WD_DIR,"GSE246088/counts/featurecounts/gene_counts.txt"))$counts,ref_len)
g305<-tpm_of(load_fc(file.path(WD_DIR,"GSE305484/counts/featurecounts/gene_counts.txt"))$counts,ref_len)
g328<-tpm_of(load_fc(file.path(WD_DIR,"GSE246328/counts/featurecounts/gene_counts.txt"))$counts,ref_len)
g565<-tpm_of(load_fc(file.path(WD_DIR,"GSE292565/counts/featurecounts/gene_counts.txt"))$counts,ref_len)
g220d<-merge(fread(file.path(WD_DIR,"GSE220575/metadata/sample_metadata.csv")),fread(file.path(WD_DIR,"GSE220575/metadata/gsm_to_srr.tsv")),by="gsm")[condition %in% c("MASH","HCC"),srr]
g246mm<-merge(fread(file.path(WD_DIR,"GSE246088/metadata/sample_metadata.csv")),fread(file.path(WD_DIR,"GSE246088/metadata/gsm_to_srr.tsv")),by.x="geo_accession",by.y="gsm")
g246d<-g246mm[genotype=="Plvap_Control"&diet %in% c("Western_Diet","High_Fat_Diet"),srr]
g305d<-merge(fread(file.path(WD_DIR,"GSE305484/metadata/sample_metadata.csv")),fread(file.path(WD_DIR,"GSE305484/metadata/gsm_to_srr.tsv")),by.x="gsm_accession",by.y="gsm")[condition!="Chow",srr]
g328d<-fread(file.path(WD_DIR,"GSE246328/metadata/sample_metadata.csv"))[condition=="Disease",srr]
g565d<-fread(file.path(WD_DIR,"GSE292565/metadata/sample_metadata.csv"))[condition=="Disease",srr]
src<-list(mcd_tpm[,intersect(colnames(mcd_tpm),dz_main),drop=FALSE], pub_tpm[,intersect(colnames(pub_tpm),dz_main),drop=FALSE],
  g220[,intersect(colnames(g220),g220d),drop=FALSE], g246[,intersect(colnames(g246),g246d),drop=FALSE],
  g305[,intersect(colnames(g305),g305d),drop=FALSE], g328[,intersect(colnames(g328),g328d),drop=FALSE], g565[,intersect(colnames(g565),g565d),drop=FALSE])
allg<-Reduce(intersect,lapply(src,rownames)); pooled<-do.call(cbind,lapply(src,function(m) m[allg,,drop=FALSE]))
tpm_vec<-setNames(apply(pooled,1,median),allg)
scoreable <- function(genes){ g<-genes[genes %in% names(tpm_vec)]; thr<-ifelse(bt_of[g]=="lncRNA",LNC_GATE,PC_GATE); g[tpm_vec[g]>=thr] }
cat(sprintf("Pooled disease samples: %d\n", ncol(pooled)))

# ---- GWAS nearest-gene set (finemap tier) ----------------------------------
gpc  <- fread(CAND_PC)[axis_finemap==1]; glnc <- fread(CAND_LNC)[axis_finemap==1]
G_all <- unique(strip_v(c(gpc$gene_id_mouse, glnc$gene_id_mouse)))
G <- intersect(G_all, pc_lnc_ids)
sym_by_mouse <- setNames(c(gpc$gene_symbol_human, glnc$gene_symbol_human), strip_v(c(gpc$gene_id_mouse, glnc$gene_id_mouse)))
cat(sprintf("GWAS finemap mouse-ortholog genes (PC+lncRNA): %d ; pass TPM gate: %d\n",
            length(G), length(scoreable(G))))

# ---- overlap by cutoff -----------------------------------------------------
CUTS <- c(0.10,0.15,0.20,0.30,0.40,0.50)
res <- rbindlist(lapply(CUTS, function(thr){
  sp <- get_spine(thr); lib <- get_lib(thr); libg <- scoreable(lib)
  data.table(
    cutoff               = thr,
    library_size         = length(libg),
    GWAS_in_human_spine  = length(intersect(G, sp)),
    GWAS_in_mouse_conf   = length(intersect(G, mouse_conf)),
    GWAS_in_library_gated= length(intersect(G, libg)),
    GWAS_netnew_gated    = length(setdiff(scoreable(G), libg)),  # pass TPM, not already in lib
    GWAS_netnew_nogate   = length(setdiff(G, lib)))
}))
cat("\n=== GWAS nearest-gene capture by current library, by human ashr cutoff ===\n"); print(res)

# ---- Option F focus (ashr>0.20) --------------------------------------------
libF  <- get_lib(CHOSEN); libFg <- scoreable(libF)
G_sc  <- scoreable(G)
in_F  <- intersect(G, libFg)
new_sc<- setdiff(G_sc, libFg)
cat(sprintf("\n=== OPTION F (ashr>0.20 + TPM gate, %d genes) + GWAS tier ===\n", length(libFg)))
cat(sprintf("  GWAS genes total (PC+lnc, ortholog-mapped): %d\n", length(G)))
cat(sprintf("  - already in Option F: %d\n", length(in_F)))
cat(sprintf("  - NOT in F but pass TPM gate (net-new): %d\n", length(new_sc)))
cat(sprintf("  - NOT in F and fail/again no TPM: %d\n", length(setdiff(G, libFg)) - length(new_sc)))
cat(sprintf("  Augmented library = Option F + net-new scoreable GWAS = %d + %d = %d genes\n",
            length(libFg), length(new_sc), length(union(libFg, new_sc))))

# net-new gene list with human symbols
new_dt <- data.table(mouse_ensembl=new_sc, gene_symbol_human=sym_by_mouse[new_sc],
                     biotype=bt_of[new_sc], tpm=round(tpm_vec[new_sc],2))
setorder(new_dt, -tpm)
fwrite(res,    file.path(DATD,"gwas_overlap_by_cutoff.csv"))
fwrite(new_dt, file.path(DATD,"gwas_netnew_genes_optionF.csv"))
cat(sprintf("\nNet-new GWAS genes added to Option F (%d), top by TPM:\n", nrow(new_dt)))
print(head(new_dt[, .(gene_symbol_human, biotype, tpm)], 30))
cat("\nWrote: gwas_overlap_by_cutoff.csv , gwas_netnew_genes_optionF.csv\n")
