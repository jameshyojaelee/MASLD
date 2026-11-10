#!/usr/bin/env Rscript
# Verifier C step 3c: full intron list vs list without empty elements, all S2K variants.
suppressPackageStartupMessages({ library(data.table); library(GenomicFeatures); library(TxDb.Hsapiens.UCSC.hg19.knownGene) })
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
a <- fread(file.path(R, "figures/candidates/fig2-coloc-r2-20260924T053244Z/variant_class/coloc_variant_annotation.csv"))
main <- fread(file.path(R, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[placement == "main", study_name]
v <- unique(a[study %in% main & gene_symbol != "", .(chr, pos_hg19)])
g <- GRanges(paste0("chr", v$chr), IRanges(v$pos_hg19, width = 1))
ibt <- intronsByTranscript(txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene)
cat("intron lists:", length(ibt), "; empty:", sum(lengths(ibt) == 0), "\n")
key <- function(m, nm) unique(paste(m$xHits, nm[m$transcriptsHits]))
full <- mapToTranscripts(g, ibt); ne <- ibt[lengths(ibt) > 0]; nonempty <- mapToTranscripts(g, ne)
direct <- findOverlaps(g, ibt, type = "within")
kd <- unique(paste(queryHits(direct), names(ibt)[subjectHits(direct)]))
kf <- key(full, names(ibt)); kn <- key(nonempty, names(ne))
cat(sprintf("variant-transcript pairs: direct=%d full_list=%d nonempty_list=%d\n", length(kd), length(kf), length(kn)))
cat(sprintf("full-list pairs not in direct=%d; nonempty-list pairs not in direct=%d; direct pairs missed by full=%d\n",
            length(setdiff(kf, kd)), length(setdiff(kn, kd)), length(setdiff(kd, kf))))
cat(sprintf("variants with any intron hit: direct=%d full=%d nonempty=%d\n",
            length(unique(queryHits(direct))), length(unique(full$xHits)), length(unique(nonempty$xHits))))
