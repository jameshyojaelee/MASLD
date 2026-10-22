# H55 planning: freeze the acute-phase / immediate-early exclusion set and count overlap with the primary family.
suppressMessages({library(msigdbr); library(data.table)})
S <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/plan-H55/"
h <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
go <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))
x <- rbind(h[gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"], go[gs_name == "GOBP_ACUTE_PHASE_RESPONSE"], fill = TRUE)
cat("db_version:", paste(unique(x$db_version), collapse = ","), "\n")
ex <- unique(x[, .(gs_name, ensembl_gene, gene_symbol)])
fwrite(ex, paste0(S, "exclusion_set_tnfa_acutephase.tsv"), sep = "\t")
fam <- fread(paste0(S, "primary_family_genes.tsv"))
u <- unique(ex$ensembl_gene)
cat("exclusion set genes:", length(u), " in primary family:", sum(fam$ens %in% u), "\n")
cat("family genes removed:", paste(sort(fam$symbol[fam$ens %in% u]), collapse = ","), "\n")
