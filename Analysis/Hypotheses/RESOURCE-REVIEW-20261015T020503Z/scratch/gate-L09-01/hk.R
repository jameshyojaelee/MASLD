suppressPackageStartupMessages(library(msigdbr))
cat("msigdbr version", as.character(packageVersion("msigdbr")), "\n")
a <- tryCatch(msigdbr(species="Homo sapiens", collection="C2"), error=function(e) msigdbr(species="Homo sapiens", category="C2"))
h <- a[a$gs_name=="HOUNKPE_HOUSEKEEPING_GENES",]
cat("rows", nrow(h), "unique symbols", length(unique(h$gene_symbol)), "unique ensembl", length(unique(h$ensembl_gene)), "\n")
print(unique(h[,intersect(c("gs_name","gs_subcollection","gs_subcat","gs_pmid","gs_description","gs_url","db_version","gs_source_species"), colnames(h))]))
write.table(unique(h[,c("gene_symbol","ensembl_gene")]), "hk_genes.tsv", sep="\t", quote=FALSE, row.names=FALSE)
