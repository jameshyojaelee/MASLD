suppressPackageStartupMessages({library(data.table); library(limma)})
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
raw <- fread(file.path(BASE, "data/PXD051911/liver_protein_quant.txt"))
meta <- fread(file.path(BASE, "data/PXD051911/meta_data.txt"))
sc <- setdiff(names(raw), c("ProteinAccessions","Genes","ProteinDescriptions"))
meta <- unique(meta[!is.na(liver_proteomics_filename) & liver_proteomics_filename %in% sc], by="liver_proteomics_filename")
meta[, sample_order := match(liver_proteomics_filename, sc)]; setorder(meta, sample_order)
meta[, group := ifelse(saf_diagnosis=="No_MASLD","Control","MASLD")]
meta[, batch := ifelse(grepl("^2019", liver_proteomics_filename),"2019","2020")]
meta[, `:=`(age_z=as.numeric(scale(as.numeric(alder))), bmi_z=as.numeric(scale(as.numeric(bmi))), sex=factor(gender),
            Fib=as.numeric(sub("^F","",kleiner_fibrosis_grade)), NAS=as.numeric(nafld_activity_score))]
s <- meta$liver_proteomics_filename
raw <- raw[!is.na(Genes) & Genes!="" & !grepl(";",Genes,fixed=TRUE)]
x <- as.matrix(raw[, ..s]); storage.mode(x) <- "double"; x[!is.finite(x) | x<=0] <- NA; x <- log2(x)
keep <- rowMeans(is.na(x)) <= 0.5; x <- x[keep,]; g <- raw$Genes[keep]
x <- normalizeBetweenArrays(x, method="quantile")
ps <- function(p, d, f){ X <- model.matrix(~ batch + age_z + bmi_z + sex, data=d)
  r <- cor(residuals(lm.fit(X, rank(p))), residuals(lm.fit(X, rank(d[[f]])))); r }
cat("All-58 on raw quantile abundance (should equal panel4c_histology_partial):\n")
for (gene in c("GNMT","MAT1A","CYP2C19")) { p <- x[which(g==gene),]; cat(gene, round(ps(p, meta, "Fib"),6), round(ps(p, meta, "NAS"),6), "\n") }
adj <- fread(file.path(BASE,"Analysis/Multimodal_Program_Projection/results/proteomics/panel4c_adjusted_abundance.tsv"))
m <- meta[group=="MASLD"]
cat("Within-MASLD on display-adjusted abundance (candidate input):\n")
for (gene in c("GNMT","MAT1A","CYP2C19")) { gg <- gene; p <- unlist(adj[adj$gene==gg, m$liver_proteomics_filename, with=FALSE]); cat(gene, round(ps(p, m, "Fib"),6), round(ps(p, m, "NAS"),6), "\n") }
