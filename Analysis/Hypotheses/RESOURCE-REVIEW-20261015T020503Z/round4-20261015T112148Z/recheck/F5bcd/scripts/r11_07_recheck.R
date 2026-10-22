suppressPackageStartupMessages({library(data.table); library(limma)})
set.seed(20260923)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
raw <- fread(file.path(BASE, "data/PXD051911/liver_protein_quant.txt"))
meta <- fread(file.path(BASE, "data/PXD051911/meta_data.txt"))
sc <- setdiff(names(raw), c("ProteinAccessions","Genes","ProteinDescriptions"))
meta <- unique(meta[!is.na(liver_proteomics_filename) & liver_proteomics_filename %in% sc], by="liver_proteomics_filename")
cat("n samples", nrow(meta), "\n")
meta[, group := ifelse(saf_diagnosis=="No_MASLD","Control","MASLD")]
meta[, batch := ifelse(grepl("^2019", liver_proteomics_filename),"2019","2020")]
meta[, `:=`(age=as.numeric(alder), bmi=as.numeric(bmi), sex=factor(gender),
            Fib=as.numeric(sub("^F","",kleiner_fibrosis_grade)), NAS=as.numeric(nafld_activity_score))]
print(table(meta$group)); print(meta[, .(fib_range=paste(range(Fib,na.rm=T),collapse="-"), nas_range=paste(range(NAS,na.rm=T),collapse="-"), nNA_fib=sum(is.na(Fib)), nNA_nas=sum(is.na(NAS))), by=group])
s <- meta$liver_proteomics_filename
raw <- raw[!is.na(Genes) & Genes!="" & !grepl(";",Genes,fixed=TRUE)]
x <- as.matrix(raw[, ..s]); storage.mode(x) <- "double"; x[!is.finite(x) | x<=0] <- NA; x <- log2(x)
keep <- rowMeans(is.na(x)) <= 0.5; x <- x[keep,]; g <- raw$Genes[keep]
x <- normalizeBetweenArrays(x, method="quantile")
res <- list()
m <- meta[group=="MASLD"]
for (gene in c("GNMT","MAT1A","CYP2C19")) {
  idx <- which(g==gene); cat(gene, "rows", length(idx), "\n")
  p <- if (length(idx)==1) x[idx,] else apply(x[idx,,drop=FALSE],2,median,na.rm=TRUE)
  names(p) <- s
  for (f in c("Fib","NAS")) {
    d <- m[is.finite(p[m$liver_proteomics_filename]) & is.finite(m[[f]])]
    X <- model.matrix(~ batch + age + bmi + sex, data=d)
    pr <- residuals(lm.fit(X, rank(p[d$liver_proteomics_filename])))
    fr <- residuals(lm.fit(X, rank(d[[f]])))
    r <- cor(pr, fr); df <- nrow(d) - qr(X)$rank - 1; t <- r*sqrt(df/(1-r^2))
    res[[length(res)+1]] <- data.table(gene, feature=f, n=nrow(d), rho=r, df=df, p=2*pt(-abs(t),df))
  }
}
o <- rbindlist(res); o[, q := p.adjust(p,"BH")]; print(o)
fwrite(o, file.path(commandArgs(TRUE)[1], "r11_07_recheck_raw_abundance.tsv"), sep="\t")
# Figure 5C histology column: who is in it
h <- fread(file.path(BASE,"Analysis/Multimodal_Program_Projection/results/proteomics/panel4c_histology_partial.tsv"))
print(h[gene %in% c("GNMT","MAT1A","CYP2C19") & feature %in% c("Fibrosis","NAS")])
print(table(h$n))
