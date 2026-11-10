#!/usr/bin/env Rscript
# A4 step 0: inputs only, no association.
#   work/peaks_annot.tsv     peak index, hg38 coordinates, width, centre, GC (BSgenome hg38)
#   work/genes_tss.tsv       GSE296875 RNA features with a TSS in the cellranger 2024-A GTF
#   work/cs_members_nc.tsv   gated noncoding credible-set members, hg19 and lifted hg38
#   work/loci.tsv            physical loci (lead-member single linkage, gap > 1 Mb)
#   work/variant_peak_cs.tsv peak x credible set noncoding PIP sums
# Usage: 00_inputs.R <A4_DIR>
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
  library(BSgenome.Hsapiens.UCSC.hg38); library(Biostrings)
})
A4 <- commandArgs(TRUE)[1]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
H24 <- file.path(ROOT, "Analysis/Hypotheses/H24-20261008T211750Z/gse296875")
GENES <- file.path(ROOT, "Analysis/MASLD_Model_Benchmark/executions/gse296875-donor-lineage-pseudobulk-21107249/pseudobulk/genes.tsv")
GTF <- "/gpfs/commons/home/jameslee/reference_genome/cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
MEMB <- file.path(ROOT, "GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv")
CHAIN <- file.path(ROOT, "data/broadaway_eqtl/hg19ToHg38.over.chain")
GAP <- 1e6
PRIMARY <- paste0("chr", c(1:22, "X"))
W <- file.path(A4, "work"); dir.create(W, showWarnings = FALSE)

# ---- peaks ----
pk <- fread(file.path(H24, "peaks.tsv"))
pk[, peak_index := .I]
pk[, c("chr", "start", "end") := tstrsplit(peak_id, "-", fixed = TRUE)]
pk[, `:=`(start = as.integer(start), end = as.integer(end))]
pk[, `:=`(width = end - start + 1L, centre = (start + end) %/% 2L, primary = chr %in% PRIMARY)]
gr <- GRanges(pk[primary == TRUE, chr], IRanges(pk[primary == TRUE, start], pk[primary == TRUE, end]))
pk[primary == TRUE, gc := as.numeric(letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, gr), "GC", as.prob = TRUE))]
fwrite(pk, file.path(W, "peaks_annot.tsv"), sep = "\t")
cat("peaks:", nrow(pk), " primary:", sum(pk$primary), "\n")

# ---- genes ----
g <- fread(GENES)
gtf <- import(GTF, feature.type = "gene")
gt <- data.table(ensembl = sub("\\..*$", "", gtf$gene_id), chr = as.character(seqnames(gtf)),
                 start = start(gtf), end = end(gtf), strand = as.character(strand(gtf)),
                 gene_type = gtf$gene_type, gtf_name = gtf$gene_name)
gt[, tss := fifelse(strand == "-", end, start)]
gt <- unique(gt, by = "ensembl")
g <- merge(g[, .(gene_index, ensembl = ensembl_stable_id, symbol)], gt, by = "ensembl", all.x = TRUE, sort = FALSE)
setorder(g, gene_index)
g[, has_tss := !is.na(tss) & chr %in% PRIMARY]
fwrite(g, file.path(W, "genes_tss.tsv"), sep = "\t")
cat("genes:", nrow(g), " with primary TSS:", sum(g$has_tss), "\n")

# ---- credible-set members ----
m <- fread(MEMB, colClasses = list(character = c("chromosome", "locus")))
m <- m[tolower(as.character(pip_sum_gate_passed)) == "true"]
cat("gated credible sets:", uniqueN(m$credible_set_uid), "\n")
print(m[, .(n_cs = uniqueN(credible_set_uid)), by = trait_scope])
lead <- m[order(-normalized_susie_pip)][, .SD[1], by = credible_set_uid][
  , .(credible_set_uid, lead_chr = chromosome, lead_pos = position)]
setorder(lead, lead_chr, lead_pos)
lead[, locus_id := {
  new <- c(TRUE, diff(lead_pos) > GAP)
  paste0(lead_chr[1], ":", cumsum(new))
}, by = lead_chr]
m <- merge(m, lead[, .(credible_set_uid, locus_id)], by = "credible_set_uid")
loci <- m[, .(chrom = chromosome[1], span_start_hg19 = min(position), span_end_hg19 = max(position),
              n_cs = uniqueN(credible_set_uid),
              n_direct_cs = uniqueN(credible_set_uid[trait_scope == "tier1_direct_masld_pdff"]),
              n_enzyme_cs = uniqueN(credible_set_uid[trait_scope == "tier2_liver_enzyme"]),
              studies = paste(sort(unique(study)), collapse = ";")), by = locus_id]
loci[, scope := fifelse(n_direct_cs > 0 & n_enzyme_cs > 0, "shared", fifelse(n_direct_cs > 0, "direct", "enzyme"))]
fwrite(loci, file.path(W, "loci.tsv"), sep = "\t")
cat("loci:", nrow(loci), "\n"); print(loci[, .N, by = scope])

nc <- m[consequence_category == "other_noncoding"]
g19 <- GRanges(paste0("chr", nc$chromosome), IRanges(nc$position, width = 1L))
g19$row <- seq_len(nrow(nc))
lo <- liftOver(g19, import.chain(CHAIN))
n_hits <- lengths(lo)
lo1 <- unlist(lo[n_hits == 1L])
nc[, pos38 := NA_integer_]; nc[, chr38 := NA_character_]
nc[lo1$row, `:=`(pos38 = start(lo1), chr38 = as.character(seqnames(lo1)))]
nc[chr38 != paste0("chr", chromosome), pos38 := NA_integer_]
cat("noncoding members:", nrow(nc), " lifted 1:1:", sum(!is.na(nc$pos38)), "\n")
fwrite(nc[, .(credible_set_uid, study, trait_scope, ancestry, locus_id, variant_id, chromosome,
              position_hg19 = position, chr38, pos38, normalized_susie_pip)],
       file.path(W, "cs_members_nc.tsv"), sep = "\t")

v <- nc[!is.na(pos38)]
ov <- findOverlaps(GRanges(v$chr38, IRanges(v$pos38, width = 1L)),
                   GRanges(pk$chr, IRanges(pk$start, pk$end)))
vp <- data.table(v[queryHits(ov), .(credible_set_uid, locus_id, study, trait_scope, variant_id, normalized_susie_pip)],
                 peak_index = pk$peak_index[subjectHits(ov)], peak_id = pk$peak_id[subjectHits(ov)])
vps <- vp[, .(pip_sum = sum(normalized_susie_pip), n_members = .N,
              top_variant = variant_id[which.max(normalized_susie_pip)],
              top_pip = max(normalized_susie_pip)),
          by = .(peak_index, peak_id, credible_set_uid, locus_id, study, trait_scope)]
fwrite(vps, file.path(W, "variant_peak_cs.tsv"), sep = "\t")
cat("peak x credible-set rows:", nrow(vps), " peaks:", uniqueN(vps$peak_index),
    " with pip_sum >= 0.1:", uniqueN(vps[pip_sum >= 0.1, peak_index]), "\n")
writeLines(capture.output(sessionInfo()), file.path(A4, "sessionInfo_00.txt"))
