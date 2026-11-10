suppressMessages(library(data.table)); a<-commandArgs(TRUE)
g<-rownames(readRDS("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/merged_dge.rds")); G<-length(g)
n<-integer(G); for(b in 1:200){ix<-readRDS(sprintf("%s/draws/perm_%03d.rds",a[1],b))$changed_idx; n[ix]<-n[ix]+1L}
tc<-logical(G); tc[readRDS(file.path(a[2],"draws","t1_true.rds"))$changed_idx]<-TRUE
mine<-sub("[.].*","",g[tc & n<=9])
ct<-fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.tsv.gz",select=c("ensembl_id","rna_ancestry_sensitive_h79"))
pr<-ct[rna_ancestry_sensitive_h79==TRUE,ensembl_id]
cat("G",G,"seedB D1",length(mine),"Catalog",length(pr),"shared",length(intersect(mine,pr)),"seedB-only",length(setdiff(mine,pr)),"Catalog-only",length(setdiff(pr,mine)),"\n")
writeLines(mine,file.path(a[1],"d1_genes_seedB.txt"))
