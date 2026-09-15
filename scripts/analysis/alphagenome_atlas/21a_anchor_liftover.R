#!/usr/bin/env Rscript
# Step 21a (P2): rsID -> hg19 (SuSiE aggregate) -> hg38 (repository liftOver rule, as in 03_liftover_universeB.R) for the anchor panel.
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(rtracklayer); library(jsonlite) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT <- Sys.getenv("AGA_OUT_ROOT", ""); if (!nzchar(OUT)) stop("AGA_OUT_ROOT is unset")
out_path <- file.path(OUT, "tables/anchor_hg38.tsv"); if (file.exists(out_path)) stop("Refusing to overwrite: ", out_path)
pre <- fromJSON(file.path(BASE, "scripts/analysis/alphagenome_atlas/21_anchor_panel_prespec.json"))
rs <- unique(pre$anchors$rsid)
su <- fread(file.path(BASE, "GWAS/finemapping/runs/uniform35_v4_2026-08-04/aggregate/susie_variant_results.tsv.gz"),
            select = c("rsid", "chromosome", "position", "effect_allele", "other_allele", "palindromic"))
x <- unique(su[rsid %in% rs])
pos <- unique(x[, .(rsid, chromosome, position)])
if (any(duplicated(pos$rsid))) stop("an rsID maps to more than one hg19 position in the SuSiE aggregate")
gr <- GRanges(paste0("chr", pos$chromosome), IRanges(pos$position, width = 1L)); names(gr) <- pos$rsid
lifted <- liftOver(gr, import.chain(file.path(BASE, "data/broadaway_eqtl/hg19ToHg38.over.chain")))
n_map <- as.integer(elementNROWS(lifted)); flat <- unlist(lifted, use.names = FALSE); first <- head(cumsum(c(1L, n_map)), -1L)
pos[, n_map := n_map]; pos[, hg38_chrom := NA_character_]; pos[, hg38_pos := NA_integer_]
ok <- n_map > 0L; pos$hg38_chrom[ok] <- as.character(seqnames(flat))[first[ok]]; pos$hg38_pos[ok] <- start(flat)[first[ok]]
al <- x[, .(alleles = paste(sort(unique(c(effect_allele, other_allele))), collapse = "/"), palindromic = any(palindromic == TRUE | palindromic == "TRUE")), by = rsid]
out <- merge(pos, al, by = "rsid"); out[, hg19_chrom := chromosome]; out[, hg19_pos := position]
missing <- setdiff(rs, out$rsid)
fwrite(out[, .(rsid, hg19_chrom, hg19_pos, hg38_chrom, hg38_pos, n_map, alleles, palindromic)], out_path, sep = "\t")
message(sprintf("anchors: %d requested, %d found in the SuSiE aggregate, %d mapped once; missing: %s", length(rs), nrow(out), sum(out$n_map == 1L), paste(missing, collapse = ",")))
