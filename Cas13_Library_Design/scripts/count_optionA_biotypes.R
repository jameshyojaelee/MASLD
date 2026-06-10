#!/usr/bin/env Rscript
# count_optionA_biotypes.R -- biotype breakdown of Option A library
#   A = human ashr spine (lfsr<0.05 & shrunk_logFC>0.20) UNION mouse-confirmed
#       tier (>=3/4 diets, human-concordant), biotype PC/lncRNA,
#       then gated to mouse disease-liver pooled median TPM >= 1.
suppressPackageStartupMessages({ library(data.table) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FC_DIR <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/counts/featurecounts")
PUB_FC <- file.path(BASE, "RNA-seq/Mouse/Public_Diet_Models/counts/featurecounts/gene_counts.txt")
WD_DIR <- file.path(BASE, "RNA-seq/Mouse/Western_Diet_Datasets")
MAIN_META <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/metadata/unified_mouse_metadata.csv")
ASHR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv")
ORTHO <- file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz")
MOUSE_META <- file.path(BASE, "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv")
PERDIET <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
DIETS <- c("MCD","CDAHFD","Western","HFD"); LFSR <- 0.05; SM <- 0.5; THR <- 0.20; TPM_GATE <- 1
strip_v <- function(x) sub("[.][0-9]+$","",x)

ash <- fread(ASHR, select=c("gene","shrunk_logFC","lfsr","logFC")); ash[, hb := strip_v(gene)]
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

up_by_diet <- lapply(DIETS, function(d){ dt <- fread(file.path(PERDIET,paste0(d,"_de_results.csv")),
  select=c("gene","shrunk_logFC","lfsr")); dt[, gb := strip_v(gene)]; dt[lfsr<LFSR & shrunk_logFC>SM, gb] })
mouse_any <- unique(unlist(up_by_diet)); n_up <- sapply(mouse_any, function(g) sum(sapply(up_by_diet, function(u) g %in% u)))
mouse_3plus <- mouse_any[n_up>=3]
mc_map <- merge(ortho_m2h[mouse_ensembl %in% mouse_3plus, .(mouse_ensembl, hb=human_ensembl)], human_logfc, by="hb", all.x=TRUE)
mouse_conf <- intersect(mc_map[!is.na(hlfc) & hlfc>0, mouse_ensembl], pc_lnc_ids)
spine <- intersect(ortho_byhuman[human_ensembl %in% ash[!is.na(lfsr)&lfsr<LFSR&shrunk_logFC>THR, hb], mouse_ensembl], pc_lnc_ids)
lib <- sort(union(spine, mouse_conf))

# pooled disease-liver TPM
load_fc <- function(path){ dt <- fread(path, skip="Geneid"); sc <- setdiff(names(dt),c("Geneid","Chr","Start","End","Strand","Length"))
  len <- setNames(dt[["Length"]], strip_v(dt[["Geneid"]])); mat <- as.matrix(dt[,..sc]); rownames(mat) <- strip_v(dt[["Geneid"]])
  colnames(mat) <- basename(dirname(sc)); list(counts=mat,lengths=len) }
tpm_of <- function(counts,len){ g <- intersect(rownames(counts),names(len)); counts <- counts[g,,drop=FALSE]
  rpk <- sweep(counts,1,len[g]/1e3,"/"); sweep(rpk,2,colSums(rpk)/1e6,"/") }
meta <- fread(MAIN_META); dz <- if("group_binary" %in% names(meta)) meta[group_binary!="Control"] else meta[!grepl("Control|LFD|Chow|Ctrl",condition)]
mcd_fcs <- lapply(file.path(FC_DIR,c("gene_counts_inhouse.txt","gene_counts_gse156918.txt","gene_counts_gse205974.txt")), load_fc)
ref_len <- mcd_fcs[[1]]$lengths; mcd_tpm <- tpm_of(do.call(cbind,lapply(mcd_fcs,`[[`,"counts")),ref_len)
pub_tpm <- tpm_of(load_fc(PUB_FC)$counts, ref_len)
g220 <- load_fc(file.path(WD_DIR,"GSE220575/counts/featurecounts/gene_counts.txt")); g220t <- tpm_of(g220$counts,ref_len)
g246 <- load_fc(file.path(WD_DIR,"GSE246088/counts/featurecounts/gene_counts.txt")); g246t <- tpm_of(g246$counts,ref_len)
g305 <- load_fc(file.path(WD_DIR,"GSE305484/counts/featurecounts/gene_counts.txt")); g305t <- tpm_of(g305$counts,ref_len)
g220m <- merge(fread(file.path(WD_DIR,"GSE220575/metadata/sample_metadata.csv")), fread(file.path(WD_DIR,"GSE220575/metadata/gsm_to_srr.tsv")), by="gsm")[condition %in% c("MASH","HCC"), srr]
g246mm <- merge(fread(file.path(WD_DIR,"GSE246088/metadata/sample_metadata.csv")), fread(file.path(WD_DIR,"GSE246088/metadata/gsm_to_srr.tsv")), by.x="geo_accession", by.y="gsm")
g246d <- g246mm[genotype=="Plvap_Control" & diet %in% c("Western_Diet","High_Fat_Diet"), srr]
g305d <- merge(fread(file.path(WD_DIR,"GSE305484/metadata/sample_metadata.csv")), fread(file.path(WD_DIR,"GSE305484/metadata/gsm_to_srr.tsv")), by.x="gsm_accession", by.y="gsm")[condition!="Chow", srr]
dz_main <- dz$sample_id
src <- list(mcd_tpm[,intersect(colnames(mcd_tpm),dz_main),drop=FALSE], pub_tpm[,intersect(colnames(pub_tpm),dz_main),drop=FALSE],
  g220t[,intersect(colnames(g220t),g220m),drop=FALSE], g246t[,intersect(colnames(g246t),g246d),drop=FALSE], g305t[,intersect(colnames(g305t),g305d),drop=FALSE])
allg <- Reduce(intersect, lapply(src, rownames)); pooled <- do.call(cbind, lapply(src, function(m) m[allg,,drop=FALSE]))
tpm_vec <- setNames(apply(pooled,1,median), allg)

libA <- lib[ lib %in% names(tpm_vec) ]; libA <- libA[ tpm_vec[libA] >= TPM_GATE ]

cat(sprintf("\nOption A library (ashr>%.2f UNION mouse-confirmed, then TPM>=1): %d genes\n", THR, length(libA)))
cat("--- biotype breakdown (mouse GENCODE vM38) ---\n")
bt <- table(factor(bt_of[libA], levels=c("protein_coding","lncRNA")), useNA="ifany")
print(bt)
nl <- sum(bt_of[libA]=="lncRNA", na.rm=TRUE)
cat(sprintf("\nlncRNAs in Option A: %d  (%.1f%% of %d)\n", nl, 100*nl/length(libA), length(libA)))

# how many lncRNAs entered via human spine vs mouse-confirmed; and how many lost to TPM gate
lib_lnc  <- lib[ bt_of[lib]=="lncRNA" ]; lib_lnc <- lib_lnc[!is.na(lib_lnc)]
gated_lnc <- intersect(lib_lnc, libA)
cat(sprintf("lncRNAs before TPM gate: %d  ->  survive TPM>=1: %d  (lost %d to low expression)\n",
            length(lib_lnc), length(gated_lnc), length(lib_lnc)-length(gated_lnc)))
sp_lnc <- intersect(spine[ bt_of[spine]=="lncRNA" ], gated_lnc); sp_lnc <- sp_lnc[!is.na(sp_lnc)]
mc_lnc <- intersect(mouse_conf[ bt_of[mouse_conf]=="lncRNA" ], gated_lnc); mc_lnc <- mc_lnc[!is.na(mc_lnc)]
cat(sprintf("  of surviving lncRNAs: %d via human spine, %d via mouse-confirmed (overlap counted once in total)\n",
            length(sp_lnc), length(mc_lnc)))
