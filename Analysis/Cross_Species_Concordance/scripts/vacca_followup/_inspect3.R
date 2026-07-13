suppressPackageStartupMessages({library(data.table);library(msigdbr)})
cat("msigdbr version:", as.character(packageVersion("msigdbr")), "\n")
co <- as.data.table(msigdbr_collections())
cat("collections with KEGG/CP in name:\n")
print(co[grepl("KEGG|CP", gs_collection) | grepl("KEGG|CP", gs_subcollection),
         .(gs_collection, gs_subcollection, num_genesets)])
