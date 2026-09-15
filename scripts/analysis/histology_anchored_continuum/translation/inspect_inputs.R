suppressPackageStartupMessages(library(data.table))
g <- fread("GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv",
  select = c("gwas_name", "gene", "ensembl", "chr", "PP.H4.susie", "PP.H4.abf"))
print(g[duplicated(g[, .(gwas_name, ensembl)]) | duplicated(g[, .(gwas_name, ensembl)], fromLast = TRUE)][1:20])
cat("Duplicate keys", sum(duplicated(g[, .(gwas_name, ensembl)])), "\n")
cat("Exact duplicates", sum(duplicated(g)), "\n")
m <- fread("data/PXD051911/meta_data.txt")
l <- names(fread("data/PXD051911/liver_protein_quant.txt", nrows = 0))
p <- names(fread("data/PXD051911/plasma_protein_quant.txt", nrows = 0))
paired <- m[liver_proteomics_filename %in% l & plasma_proteomics_filename %in% p]
print(paired[, .N, by = .(sample_group, saf_diagnosis)])
print(paired[, .(unique_identifier, patient_name, liver_proteomics_filename, plasma_proteomics_filename, plasma_batch_effects)][1:5])
cat("Liver metadata rows", nrow(m[liver_proteomics_filename %in% l]), "Paired", nrow(paired), "unique donors", uniqueN(paired$patient_name), "\n")
