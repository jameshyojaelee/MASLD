# H86 writer check: can the 119 unmatched protein symbols be remapped with locally installed
# annotation (no download)? Counts resolution only; no join to Catalog evidence classes.
suppressPackageStartupMessages({library(data.table); library(HGNChelper); library(org.Hs.eg.db); library(limma)})
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
g <- fread(file.path(Sys.getenv("W"),"tmp/catalog_symbols.tsv"))
p <- fread(file.path(root,"Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv"), na.strings="")
cat("protein rows", nrow(p), "unique", uniqueN(p$gene), "padj<0.05", sum(p$padj<0.05), "\n")
cat("BH recompute max abs diff", max(abs(p.adjust(p$pvalue,"BH")-p$padj)), "\n")
um <- setdiff(p$gene, g$symbol); cat("unmatched", length(um), "\n")
cat("HGNChelper", as.character(packageVersion("HGNChelper")), " org.Hs.eg.db", as.character(packageVersion("org.Hs.eg.db")), "\n")
h <- suppressWarnings(checkGeneSymbols(um, unmapped.as.na=TRUE))
h <- as.data.table(h)
h[, n_sugg := ifelse(is.na(Suggested.Symbol), 0L, lengths(strsplit(Suggested.Symbol, " /// ", fixed=TRUE)))]
h[, in_cat := !is.na(Suggested.Symbol) & n_sugg==1L & Suggested.Symbol %in% g$symbol]
h[, collides := in_cat & Suggested.Symbol %in% p$gene]
cat("HGNChelper: unique suggestion", sum(h$n_sugg==1), " ambiguous", sum(h$n_sugg>1), " none", sum(h$n_sugg==0), " unique and in Catalog", sum(h$in_cat), " target already a protein row", sum(h$collides), "\n")
a <- as.data.table(alias2SymbolTable(um, species="Hs"))
a <- data.table(old=um, new=alias2SymbolTable(um, species="Hs"))
a[, in_cat := !is.na(new) & new %in% g$symbol]
cat("org.Hs.eg.db alias2SymbolTable: resolved", sum(!is.na(a$new)), " in Catalog", sum(a$in_cat), "\n")
m <- merge(h[, .(old=x, hgnc=Suggested.Symbol, n_sugg)], a[, .(old, org=new)], by="old")
cat("agree when both unique", m[n_sugg==1 & !is.na(org), sum(hgnc==org)], "/", m[n_sugg==1 & !is.na(org), .N], "\n")
print(m[old %in% c("MARC1","MARC2","AARS","CD97","CERT","ATP5MD")])
fwrite(m, file.path(Sys.getenv("W"), "remap_feasibility_119.tsv"), sep="\t")
