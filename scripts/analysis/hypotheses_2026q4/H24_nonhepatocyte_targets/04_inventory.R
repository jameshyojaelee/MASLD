#!/usr/bin/env Rscript
# H24 step 4: test universe and family sizes. Inputs only; no peak-gene
# association, colocalization, or QTL statistic is computed or read here.
#
# Revision v2 (written before any association was computed). The v1 run
# (inventory/) applied edgeR::filterByExpr to ATAC peak counts. Its CPM cutoff
# is calibrated to ~20k genes and removed every variant peak in every
# non-hepatocyte lineage (an RNA rule applied to 306,706 peaks). v2 keeps a
# peak in a lineage if its raw pseudobulk count is >= 5 in >= 50% of that
# lineage's donor units. v2 writes to inventory_v2/; inventory/ is kept as is.
#
# 1. Lift credible-set members (H24 and reference clusters) hg19 -> hg38.
# 2. Physical loci: H24 clusters whose span +/- 500 kb windows overlap on the
#    same chromosome are merged (4:2+4:9, 8:2+8:8); same for reference clusters.
# 3. GSE296875 peaks that contain at least one member ("variant peaks").
# 4. Genes whose TSS (GENCODE v44, refdata-gex-GRCh38-2024-A) lies within
#    500 kb of the variant-peak centre.
# 5. Per lineage, donor units with >= 20 nuclei; RNA filterByExpr (intercept
#    design); ATAC count rule above.
# 6. Families: H24 (primary), REF (separate BH): PNPLA3-region clusters
#    22:3+22:9 all lineages; Zhu 2026 CRISPRi DAV peaks (10:15, 19:2+19:10)
#    hepatocyte only.
# 7. N_perm = min(2e6, max(1e5, 40 * m_H24 rounded up to 1e4)).
# Usage: Rscript 04_inventory.R <OUT_DIR>
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer); library(edgeR); library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE); OUT <- args[[1]]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CHAIN <- file.path(ROOT, "data/broadaway_eqtl/hg19ToHg38.over.chain")
GTF <- "/gpfs/commons/home/jameslee/reference_genome/cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
PB <- file.path(ROOT, "Analysis/MASLD_Model_Benchmark/executions/gse296875-donor-lineage-pseudobulk-21107249/pseudobulk")
ATAC244 <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/counts/GSE244832")
WIN <- 500000L; MIN_NUC <- 20L; ATAC_MIN_COUNT <- 5; ATAC_MIN_FRAC <- 0.5
REF_PNPLA3 <- c("22:3", "22:9"); REF_ZHU_CLUSTERS <- c("10:15", "19:2", "19:10")
inv <- file.path(OUT, "inventory_v2"); dir.create(inv, showWarnings = FALSE)
atac_keep <- function(m) rowSums(m >= ATAC_MIN_COUNT) >= ATAC_MIN_FRAC * ncol(m)

loci <- fread(file.path(OUT, "loci_direct_clusters.tsv"))
memb_all <- fread(file.path(OUT, "direct_members_hg19.tsv"), colClasses = list(character = "chromosome"))
h24 <- loci[h24_primary == TRUE]

# physical loci: merge clusters whose +/- WIN windows overlap
merge_physical <- function(cl) {
  cl <- cl[order(chrom, span_start_hg19)]
  cl[, `:=`(ws = span_start_hg19 - WIN, we = span_end_hg19 + WIN)]
  grp <- integer(nrow(cl)); g <- 0L; cur_chr <- ""; cur_end <- -Inf
  for (i in seq_len(nrow(cl))) {
    if (cl$chrom[i] != cur_chr || cl$ws[i] > cur_end) { g <- g + 1L; cur_chr <- cl$chrom[i]; cur_end <- cl$we[i] }
    else cur_end <- max(cur_end, cl$we[i])
    grp[i] <- g
  }
  cl[, grp := grp]
  cl[, physical_locus := paste(sort(cluster_id), collapse = "+"), by = grp]
  cl[, .(cluster_id, physical_locus)]
}
phys_h24 <- merge_physical(copy(h24))
phys_ref <- merge_physical(copy(loci[cluster_id %in% c(REF_PNPLA3, REF_ZHU_CLUSTERS)]))
phys <- rbind(phys_h24[, set := "H24"], phys_ref[, set := "REF"])
strict_ids <- h24[h24_strict_no_any_coloc == TRUE, cluster_id]
phys_tab <- merge(phys, loci[, .(cluster_id, chrom, span_start_hg19, span_end_hg19, direct_studies, ancestries,
                                 top_member, top_member_pip, h24_strict_no_any_coloc)], by = "cluster_id")
fwrite(phys_tab, file.path(inv, "physical_loci.tsv"), sep = "\t")

# 1. liftover (unique 1:1 mappings only)
memb <- memb_all[cluster_id %in% phys$cluster_id]
chain <- import.chain(CHAIN)
lift <- function(chr, pos) {
  lo <- liftOver(GRanges(paste0("chr", chr), IRanges(pos, width = 1)), chain)
  ok <- lengths(lo) == 1L; out <- rep(NA_integer_, length(pos)); out[ok] <- start(unlist(lo[ok])); out
}
memb[, pos_hg38 := lift(chromosome, position)]
memb[, lifted := !is.na(pos_hg38)]
memb <- merge(memb, phys[, .(cluster_id, physical_locus, set)], by = "cluster_id")
fwrite(memb, file.path(inv, "members_hg19_hg38.tsv"), sep = "\t")

# 3. variant peaks
peaks <- fread(file.path(OUT, "gse296875/peaks.tsv"))
pp <- tstrsplit(peaks$peak_id, "-")
pgr <- GRanges(pp[[1]], IRanges(as.integer(pp[[2]]), as.integer(pp[[3]])))
m38 <- memb[lifted == TRUE]
hit <- findOverlaps(GRanges(paste0("chr", m38$chromosome), IRanges(m38$pos_hg38, width = 1)), pgr)
vp <- data.table(m38[queryHits(hit)], peak_index = subjectHits(hit))
vpeak <- vp[, .(peak_id = peaks$peak_id[peak_index],
                cluster_id = paste(sort(unique(cluster_id)), collapse = ";"),
                physical_locus = physical_locus[1],
                n_members = uniqueN(variant_id),
                max_cs_pip_in_peak = max(.SD[, .(s = sum(normalized_susie_pip)), by = credible_set_uid]$s),
                variants = paste(sort(unique(variant_id)), collapse = ";")),
            by = .(peak_index, set)]
fwrite(vpeak, file.path(inv, "variant_peaks.tsv"), sep = "\t")

# Zhu 2026 CRISPRi DAVs (positions: Ensembl GRCh37 REST lookup, inputs/)
zhu <- fread(file.path(OUT, "inputs/zhu2026_crispri_positions_hg19.tsv"), colClasses = list(character = "chrom"))
zhu[, pos_hg38 := NA_integer_]
zhu[!is.na(pos_hg19), pos_hg38 := lift(chrom, pos_hg19)]
zhu[, `:=`(any_direct_cluster = NA_character_, h24_cluster = NA_character_, peak_index = NA_integer_)]
for (i in seq_len(nrow(zhu))) {
  if (is.na(zhu$pos_hg19[i])) next
  w <- loci[chrom == zhu$chrom[i] & span_start_hg19 - WIN <= zhu$pos_hg19[i] & span_end_hg19 + WIN >= zhu$pos_hg19[i]]
  if (nrow(w)) {
    zhu$any_direct_cluster[i] <- paste(w$cluster_id, collapse = ";")
    if (any(w$h24_primary)) zhu$h24_cluster[i] <- paste(w[h24_primary == TRUE]$cluster_id, collapse = ";")
  }
  if (!is.na(zhu$pos_hg38[i])) {
    h <- subjectHits(findOverlaps(GRanges(paste0("chr", zhu$chrom[i]), IRanges(zhu$pos_hg38[i], width = 1)), pgr))
    if (length(h)) zhu$peak_index[i] <- h[1]
  }
}
fwrite(zhu, file.path(inv, "zhu2026_crispri_overlap.tsv"), sep = "\t")
zpk <- zhu[!is.na(peak_index) & !is.na(any_direct_cluster) &
             vapply(strsplit(any_direct_cluster, ";"), function(x) any(x %in% REF_ZHU_CLUSTERS), logical(1))]
zpeak <- if (nrow(zpk)) data.table(peak_index = zpk$peak_index, peak_id = peaks$peak_id[zpk$peak_index],
                                   cluster_id = zpk$any_direct_cluster, physical_locus = paste0("ZHU_", zpk$rsid),
                                   set = "REF_ZHU", n_members = NA_integer_, max_cs_pip_in_peak = NA_real_,
                                   variants = zpk$rsid) else NULL

# 4. genes within window of peak centre
gtf <- import(GTF, feature.type = "gene")
gtf$gene_id <- sub("\\..*$", "", gtf$gene_id)
tss <- resize(gtf, 1L, fix = "start")
genes <- fread(file.path(PB, "genes.tsv"))
tss <- tss[tss$gene_id %in% genes$ensembl_stable_id]
allpk <- rbind(vpeak[set == "H24", .(peak_index, peak_id, cluster_id, physical_locus, family = "H24", max_cs_pip_in_peak)],
               vpeak[set == "REF", .(peak_index, peak_id, cluster_id, physical_locus, family = "REF_PNPLA3", max_cs_pip_in_peak)],
               if (!is.null(zpeak)) zpeak[, .(peak_index, peak_id, cluster_id, physical_locus, family = "REF_ZHU", max_cs_pip_in_peak)])
allpk <- allpk[!(family == "REF_PNPLA3" & !grepl("22:", cluster_id))]  # REF set holds PNPLA3 and Zhu clusters; keep 22:x only here
centre <- GRanges(seqnames(pgr[allpk$peak_index]),
                  IRanges((start(pgr[allpk$peak_index]) + end(pgr[allpk$peak_index])) %/% 2L, width = 1))
h2 <- findOverlaps(resize(centre, 2L * WIN + 1L, fix = "center"), tss, ignore.strand = TRUE)
pg <- data.table(allpk[queryHits(h2)], ensembl = tss$gene_id[subjectHits(h2)],
                 symbol = tss$gene_name[subjectHits(h2)],
                 distance = start(tss)[subjectHits(h2)] - start(centre)[queryHits(h2)])
pg <- unique(pg, by = c("family", "peak_index", "ensembl"))

# 5. filters per lineage
units <- fread(file.path(OUT, "gse296875/units.tsv"))
dims <- fromJSON(file.path(OUT, "gse296875/dims.json"))
rna <- matrix(readBin(file.path(OUT, "gse296875/rna_gene_by_unit.f32"), "numeric",
                      n = prod(dims$rna), size = 4, endian = "little"), nrow = dims$rna[1])
atac <- matrix(readBin(file.path(OUT, "gse296875/atac_peak_by_unit.f32"), "numeric",
                       n = prod(dims$atac), size = 4, endian = "little"), nrow = dims$atac[1])
gidx <- match(pg$ensembl, genes$ensembl_stable_id)
fam <- list(); lin_tab <- list()
for (lin in sort(unique(units$lineage_id))) {
  u <- units[lineage_id == lin & nuclei >= MIN_NUC]
  cols <- u$unit_index + 1L
  keep_g <- filterByExpr(DGEList(rna[, cols]), design = matrix(1, length(cols), 1))
  keep_p <- atac_keep(atac[, cols, drop = FALSE])
  tt <- pg[keep_p[peak_index] & keep_g[gidx]]
  if (lin != "hepatocyte") tt <- tt[family != "REF_ZHU"]
  fam[[lin]] <- tt[, lineage := lin]
  vh <- vpeak[set == "H24"]$peak_index
  lin_tab[[lin]] <- data.table(lineage = lin, donors = length(cols), wells = uniqueN(u$well_id),
                               median_nuclei = median(u$nuclei), expressed_genes = sum(keep_g),
                               accessible_peaks = sum(keep_p), h24_variant_peaks_accessible = sum(keep_p[vh]),
                               h24_tests = nrow(tt[family == "H24"]),
                               h24_physical_loci = uniqueN(tt[family == "H24"]$physical_locus),
                               ref_tests = nrow(tt[family != "H24"]))
}
fam <- rbindlist(fam); lin_tab <- rbindlist(lin_tab)
fwrite(fam, file.path(inv, "l1_test_universe.tsv"), sep = "\t")
fwrite(lin_tab, file.path(inv, "l1_lineage_summary.tsv"), sep = "\t")

# GSE244832 donors per lineage (RNA and ATAC nuclei >= 20)
r244 <- fread(file.path(OUT, "gse244832/rna_units.tsv"))
a244 <- rbindlist(lapply(c("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"), function(l)
  fread(file.path(ATAC244, paste0(l, ".donors.tsv")))[, .(donor_id, lineage, atac_nuclei = n_cells)]))
both <- merge(r244[, .(donor_id, lineage = lineage_id, rna_nuclei = nuclei)], a244, by = c("donor_id", "lineage"))
l2_tab <- both[, .(donors_both_ge20 = sum(rna_nuclei >= MIN_NUC & atac_nuclei >= MIN_NUC)), by = lineage]
fwrite(l2_tab, file.path(inv, "l2_gse244832_donors.tsv"), sep = "\t")

m_h24 <- nrow(fam[family == "H24"])
n_perm <- as.integer(min(2e6, max(1e5, ceiling(40 * m_h24 / 1e4) * 1e4)))
rep <- list(revision = "v2_atac_count_rule", atac_rule = "raw count >= 5 in >= 50% of the lineage's donor units",
            h24_clusters = nrow(h24), h24_physical_loci = uniqueN(phys_h24$physical_locus),
            h24_members = nrow(memb[set == "H24"]), h24_members_lifted = sum(memb[set == "H24"]$lifted),
            h24_variant_peaks = nrow(vpeak[set == "H24"]),
            h24_physical_loci_with_variant_peak = uniqueN(vpeak[set == "H24"]$physical_locus),
            h24_peak_gene_pairs_in_window = nrow(pg[family == "H24"]),
            family_size_m_h24 = m_h24, family_size_ref_pnpla3 = nrow(fam[family == "REF_PNPLA3"]),
            family_size_ref_zhu = nrow(fam[family == "REF_ZHU"]), n_perm = n_perm,
            lineages = lin_tab, l2 = l2_tab,
            zhu_dav_in_h24_cluster = sum(!is.na(zhu$h24_cluster)),
            zhu_dav_in_any_direct_cluster = sum(!is.na(zhu$any_direct_cluster)),
            zhu_dav_peaks_in_ref = if (!is.null(zpeak)) zpeak$variants else character(0))
write_json(rep, file.path(inv, "inventory.json"), pretty = TRUE, auto_unbox = TRUE)
print(rep)
writeLines(capture.output(sessionInfo()), file.path(inv, "sessionInfo.txt"))
