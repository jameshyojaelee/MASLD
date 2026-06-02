#!/usr/bin/env Rscript
# 56h_roadmap_replication.R
# Cross-validate GWAS-ATAC + scATAC findings against Roadmap Epigenomics E066 (Liver).
#
# Track A (A8): variant x chromHMM state overlap
#   - For 8,329 finemapped variants, compute Roadmap chromHMM state distribution
#   - Compare variants overlapping our scATAC peaks vs not
#   - Compute Fisher OR for "concordant" calls (in our peaks AND in Roadmap active states)
#   - Disease-regulon replication
#
# Track B (B4-Roadmap-addon): peak x chromHMM state overlap
#   - For each cell-type scATAC peak BED, compute % overlap with Roadmap active states
#   - Hepatocyte disease-up vs disease-down peaks

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(rtracklayer)
})

PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
ROAD_DIR <- file.path(PROJ, "data/external/roadmap_liver_E066")
GWAS_ATAC_DIR <- file.path(PROJ, "GWAS/finemapping/results/gwas_atac")
PEAKS_DIR <- file.path(PROJ, "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2")
DA_HEP <- file.path(PROJ, "Analysis/ATAC/Human_Multiome/results/snapatac2/scatac_da_corrected_hep.csv")
OUT_DIR <- GWAS_ATAC_DIR
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

# Roadmap 15-state chromHMM model labels (BED file uses E1..E15 short labels).
# Mapping per Roadmap publication (Roadmap Epigenomics Consortium, Nature 2015):
state_label <- c(
    "E1"="TssA",        "E2"="TssAFlnk",  "E3"="TxFlnk",
    "E4"="Tx",          "E5"="TxWk",      "E6"="EnhG",
    "E7"="Enh",         "E8"="ZNF/Rpts",  "E9"="Het",
    "E10"="TssBiv",     "E11"="BivFlnk",  "E12"="EnhBiv",
    "E13"="ReprPC",     "E14"="ReprPCWk", "E15"="Quies"
)
# Per task: "active regulatory" = TssA, TssFlnk, TssFlnkU, EnhG1, EnhG2, Enh, EnhWk
# In the 15-state core model: E1=TssA, E2=TssAFlnk(=TssFlnk), E6=EnhG(=EnhG1), E7=Enh
# (TssFlnkU/EnhG2/EnhWk are 18-state labels — map to same promoter/enhancer family
#  whose 15-state equivalents are already included).
ACTIVE_STATES <- c("E1","E2","E6","E7")

cat("================================================================\n")
cat(" Script 56h: Roadmap E066 replication of GWAS-ATAC + scATAC peaks\n")
cat("================================================================\n")

# -----------------------------------------------------------------------------
# 1. Load Roadmap E066 chromHMM segments
# -----------------------------------------------------------------------------
chmm_f <- file.path(ROAD_DIR, "E066_15_chromHMM_hg38.bed.gz")
stopifnot(file.exists(chmm_f))
cat(sprintf("[%s] Loading chromHMM E066: %s\n", format(Sys.time()), chmm_f))
chmm <- fread(chmm_f, sep="\t", header=FALSE,
              col.names=c("chrom","start","end","state"))
cat(sprintf("  Loaded %d segments. State distribution:\n", nrow(chmm)))
print(chmm[, .N, by=state][order(-N)])

chmm_gr <- GRanges(seqnames=chmm$chrom,
                   ranges=IRanges(start=chmm$start+1L, end=chmm$end),
                   state=chmm$state)

# -----------------------------------------------------------------------------
# 2. Load Roadmap E066 H3K27ac (broadPeak, hg38) and H3K4me3 (narrowPeak, hg38).
#    NOTE: E066 (Liver) was not assayed for DNase by Roadmap. We use H3K27ac
#    (active enhancers/promoters) and H3K4me3 (active TSS) as empirical
#    regulatory annotations. The original task mentioned "DNase" — substituting
#    H3K27ac+H3K4me3 because DNase is unavailable for E066.
# -----------------------------------------------------------------------------
h3k27ac_f <- file.path(ROAD_DIR, "E066_H3K27ac_hg38.broadPeak.gz")
h3k4me3_f <- file.path(ROAD_DIR, "E066_H3K4me3_hg38.narrowPeak.gz")
stopifnot(file.exists(h3k27ac_f), file.exists(h3k4me3_f))

read_peak <- function(f) {
    pk <- fread(f, sep="\t", header=FALSE)
    setnames(pk, c("chrom","start","end","name","score","strand",
                   "signalValue","pValue","qValue")[seq_len(ncol(pk))])
    pk
}
h3k27ac <- read_peak(h3k27ac_f)
h3k4me3 <- read_peak(h3k4me3_f)
cat(sprintf("  H3K27ac peaks (hg38): %d\n", nrow(h3k27ac)))
cat(sprintf("  H3K4me3 peaks (hg38): %d\n", nrow(h3k4me3)))

h3k27ac_gr <- GRanges(seqnames=h3k27ac$chrom,
                      ranges=IRanges(start=h3k27ac$start+1L, end=h3k27ac$end))
h3k4me3_gr <- GRanges(seqnames=h3k4me3$chrom,
                      ranges=IRanges(start=h3k4me3$start+1L, end=h3k4me3$end))
# "Empirical active" = H3K27ac OR H3K4me3
empirical_gr <- reduce(c(h3k27ac_gr, h3k4me3_gr))
cat(sprintf("  Merged H3K27ac+H3K4me3 active regions: %d\n", length(empirical_gr)))

# =============================================================================
# A8 — VARIANT x chromHMM OVERLAP
# =============================================================================
cat("\n================ A8: Variant x chromHMM overlap ================\n")

vov <- fread(file.path(GWAS_ATAC_DIR, "variant_overlap_summary.csv"))
cat(sprintf("Loaded %d finemapped variants\n", nrow(vov)))

# Drop NA hg38 coords
vov <- vov[!is.na(chr_hg38) & !is.na(pos_hg38)]
cat(sprintf("After NA filter: %d\n", nrow(vov)))

v_gr <- GRanges(seqnames=vov$chr_hg38,
                ranges=IRanges(start=vov$pos_hg38, end=vov$pos_hg38))

# Assign chromHMM state per variant
hits <- findOverlaps(v_gr, chmm_gr)
vov[, chromhmm_state := NA_character_]
vov[queryHits(hits), chromhmm_state := chmm_gr$state[subjectHits(hits)]]
# Variants without any state -> assign "0_NoCov" (rare; chromHMM segments
# usually tile the genome but the 15-state model can leave alt contigs blank)
vov[is.na(chromhmm_state), chromhmm_state := "0_NoCov"]

vov[, in_active := chromhmm_state %in% ACTIVE_STATES]

# Per-state count
state_counts <- vov[, .N, by=chromhmm_state][order(-N)]
state_counts[, frac := round(N/sum(N), 4)]
cat("\nVariant chromHMM state distribution:\n")
print(state_counts)
fwrite(state_counts, file.path(OUT_DIR, "roadmap_E066_overlap.csv"))

n_total  <- nrow(vov)
n_active <- sum(vov$in_active)
cat(sprintf("\n[A8] Variants in Roadmap active states: %d / %d (%.2f%%)\n",
            n_active, n_total, 100*n_active/n_total))

# Variants in H3K27ac broadPeak
vov[, in_h3k27ac := overlapsAny(v_gr, h3k27ac_gr)]
vov[, in_h3k4me3 := overlapsAny(v_gr, h3k4me3_gr)]
vov[, in_empirical_active := in_h3k27ac | in_h3k4me3]
cat(sprintf("[A8] Variants in Roadmap H3K27ac peaks: %d / %d (%.2f%%)\n",
            sum(vov$in_h3k27ac), n_total, 100*sum(vov$in_h3k27ac)/n_total))
cat(sprintf("[A8] Variants in Roadmap H3K4me3 peaks: %d / %d (%.2f%%)\n",
            sum(vov$in_h3k4me3), n_total, 100*sum(vov$in_h3k4me3)/n_total))
cat(sprintf("[A8] Variants in H3K27ac OR H3K4me3: %d / %d (%.2f%%)\n",
            sum(vov$in_empirical_active), n_total,
            100*sum(vov$in_empirical_active)/n_total))

# Concordance: variants overlapping our scATAC peaks vs in Roadmap active
# overlaps_any_peak == our scATAC peak overlap (provided in variant_overlap_summary)
tab <- table(scatac=vov$overlaps_any_peak, roadmap_active=vov$in_active)
cat("\nContingency (rows=scATAC, cols=Roadmap active):\n"); print(tab)
# Fisher test - 2x2 only if at least 2 levels each
if (all(dim(tab) >= c(2,2))) {
    ft <- fisher.test(tab)
    cat(sprintf("Fisher OR scATAC vs Roadmap_active: %.3f (95%% CI %.3f-%.3f), p = %.3g\n",
                ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
    or_val <- as.numeric(ft$estimate); ci <- ft$conf.int; pv <- ft$p.value
} else {
    or_val <- NA_real_; ci <- c(NA_real_, NA_real_); pv <- NA_real_
}

# Per-cell-type breakdown: hepatocyte-peak-overlapping variants vs other CT
# cell_types_overlapping is a ";"-separated list (or "")
vov[, hep_peak := grepl("Hepatocyte", cell_types_overlapping, ignore.case=TRUE)]
vov[, any_peripheral := (n_cell_types_overlapping > 0) & !hep_peak]

ct_table <- rbindlist(list(
    data.table(group="hep_peak_overlapping",
               n=sum(vov$hep_peak),
               n_active=sum(vov$hep_peak & vov$in_active),
               pct_active=round(100*sum(vov$hep_peak & vov$in_active)/max(sum(vov$hep_peak),1),3)),
    data.table(group="peripheral_only_overlapping",
               n=sum(vov$any_peripheral),
               n_active=sum(vov$any_peripheral & vov$in_active),
               pct_active=round(100*sum(vov$any_peripheral & vov$in_active)/max(sum(vov$any_peripheral),1),3)),
    data.table(group="no_scatac_peak",
               n=sum(!vov$overlaps_any_peak),
               n_active=sum(!vov$overlaps_any_peak & vov$in_active),
               pct_active=round(100*sum(!vov$overlaps_any_peak & vov$in_active)/max(sum(!vov$overlaps_any_peak),1),3))
))
cat("\n[A8] Per-cell-type breakdown:\n"); print(ct_table)
fwrite(ct_table, file.path(OUT_DIR, "roadmap_concordance.csv"))

# 2x2 Fisher for hepatocyte-overlapping vs no-peak
tab_hep <- matrix(c(
    sum(vov$hep_peak & vov$in_active),
    sum(vov$hep_peak & !vov$in_active),
    sum(!vov$overlaps_any_peak & vov$in_active),
    sum(!vov$overlaps_any_peak & !vov$in_active)
), nrow=2, byrow=TRUE,
    dimnames=list(c("hep_peak","no_peak"), c("active","not_active")))
cat("\n[A8] Hepatocyte-peak vs no-peak in Roadmap active:\n"); print(tab_hep)
ft_hep <- fisher.test(tab_hep)
cat(sprintf("Fisher OR hep-peak vs no-peak (active): %.3f (95%% CI %.3f-%.3f), p = %.3g\n",
            ft_hep$estimate, ft_hep$conf.int[1], ft_hep$conf.int[2], ft_hep$p.value))

# Save summary stats CSV
summary_dt <- data.table(
    metric=c(
        "total_variants_hg38",
        "in_chromhmm_active",
        "frac_in_chromhmm_active_pct",
        "in_h3k27ac",
        "frac_in_h3k27ac_pct",
        "in_h3k4me3",
        "frac_in_h3k4me3_pct",
        "in_empirical_active_h3k27ac_or_h3k4me3",
        "frac_in_empirical_active_pct",
        "in_scatac_peak",
        "fisher_or_scatac_vs_chromhmm_active",
        "fisher_or_scatac_vs_chromhmm_active_ci_lower",
        "fisher_or_scatac_vs_chromhmm_active_ci_upper",
        "fisher_or_scatac_vs_chromhmm_active_pvalue",
        "fisher_or_hep_vs_nopeak_chromhmm_active",
        "fisher_or_hep_vs_nopeak_pvalue"
    ),
    value=c(
        n_total,
        n_active,
        round(100*n_active/n_total, 3),
        sum(vov$in_h3k27ac),
        round(100*sum(vov$in_h3k27ac)/n_total, 3),
        sum(vov$in_h3k4me3),
        round(100*sum(vov$in_h3k4me3)/n_total, 3),
        sum(vov$in_empirical_active),
        round(100*sum(vov$in_empirical_active)/n_total, 3),
        sum(vov$overlaps_any_peak),
        round(or_val, 3),
        round(ci[1], 3),
        round(ci[2], 3),
        signif(pv, 3),
        round(as.numeric(ft_hep$estimate), 3),
        signif(ft_hep$p.value, 3)
    )
)
fwrite(summary_dt, file.path(OUT_DIR, "roadmap_A8_summary.csv"))
cat("\nWrote roadmap_A8_summary.csv\n")

# -----------------------------------------------------------------------------
# A8 sub-module: disease-regulon variant replication
# -----------------------------------------------------------------------------
cat("\n--- A8 disease-regulon replication ---\n")
md <- fread(file.path(GWAS_ATAC_DIR, "motif_disruption_scores.csv"))
# motif_in_disease_regulon TRUE rows
dis_var <- unique(md[motif_in_disease_regulon == TRUE, SNP_id])
cat(sprintf("Disease-regulon disrupting variants (unique): %d\n", length(dis_var)))

# Parse "chr:pos:ref:alt" -> chr+pos (hg38, since motif_disruption uses hg38 coords)
# md$seqnames already has "chr*" form
dis_md <- unique(md[motif_in_disease_regulon == TRUE,
                    .(SNP_id, seqnames, start, tf_name, alleleDiff)])
dis_gr <- GRanges(seqnames=dis_md$seqnames,
                  ranges=IRanges(start=dis_md$start, end=dis_md$start))

# Roadmap chromHMM state per dis-regulon variant
dis_hits <- findOverlaps(dis_gr, chmm_gr)
dis_md[, chromhmm_state := NA_character_]
dis_md[queryHits(dis_hits), chromhmm_state := chmm_gr$state[subjectHits(dis_hits)]]
dis_md[is.na(chromhmm_state), chromhmm_state := "0_NoCov"]
dis_md[, in_active := chromhmm_state %in% ACTIVE_STATES]

# Roadmap H3K27ac + H3K4me3
dis_md[, in_h3k27ac := overlapsAny(dis_gr, h3k27ac_gr)]
dis_md[, in_h3k4me3 := overlapsAny(dis_gr, h3k4me3_gr)]
dis_md[, in_empirical_active := in_h3k27ac | in_h3k4me3]

# Per-variant: collapse to unique SNP_id (motifs can hit one variant multiple times)
dis_uniq <- unique(dis_md[, .(SNP_id, seqnames, start, chromhmm_state, in_active,
                              in_h3k27ac, in_h3k4me3, in_empirical_active)])
cat(sprintf("\nUnique disease-regulon variants: %d\n", nrow(dis_uniq)))
cat(sprintf("  in Roadmap chromHMM active state: %d (%.1f%%)\n",
            sum(dis_uniq$in_active), 100*mean(dis_uniq$in_active)))
cat(sprintf("  in Roadmap H3K27ac peak:           %d (%.1f%%)\n",
            sum(dis_uniq$in_h3k27ac), 100*mean(dis_uniq$in_h3k27ac)))
cat(sprintf("  in Roadmap H3K4me3 peak:           %d (%.1f%%)\n",
            sum(dis_uniq$in_h3k4me3), 100*mean(dis_uniq$in_h3k4me3)))
cat(sprintf("  in chromHMM active OR H3K27ac/H3K4me3 peak: %d (%.1f%%)\n",
            sum(dis_uniq$in_active | dis_uniq$in_empirical_active),
            100*mean(dis_uniq$in_active | dis_uniq$in_empirical_active)))

fwrite(dis_uniq, file.path(OUT_DIR, "roadmap_disease_regulon_replication.csv"))

# =============================================================================
# B4-Roadmap-addon — PEAK x chromHMM OVERLAP
# =============================================================================
cat("\n============ B4-Roadmap-addon: Peak x chromHMM overlap ============\n")

ct_files <- list.files(PEAKS_DIR, pattern="_peaks\\.bed$", full.names=TRUE)
cat(sprintf("Found %d cell-type peak BEDs\n", length(ct_files)))

ct_results <- lapply(ct_files, function(f) {
    ct <- sub("_peaks\\.bed$", "", basename(f))
    pk <- fread(f, header=FALSE, sep="\t",
                col.names=c("chrom","start","end"),
                fill=TRUE, select=1:3)
    pk <- pk[grepl("^chr", chrom)]
    gr <- GRanges(seqnames=pk$chrom,
                  ranges=IRanges(start=pk$start+1L, end=pk$end))
    h <- findOverlaps(gr, chmm_gr)
    pk_state <- rep("NoChromHMM", length(gr))
    if (length(h) > 0) {
        ov <- pintersect(gr[queryHits(h)], chmm_gr[subjectHits(h)])
        dt <- data.table(q=queryHits(h), s=chmm_gr$state[subjectHits(h)],
                         w=width(ov))
        best <- dt[, .SD[which.max(w)], by=q]
        pk_state[best$q] <- best$s
    }
    pk[, state := pk_state]
    pk[, in_active := state %in% ACTIVE_STATES]
    pk[, in_h3k27ac := overlapsAny(gr, h3k27ac_gr)]
    pk[, in_h3k4me3 := overlapsAny(gr, h3k4me3_gr)]
    summary_row <- data.table(cell_type=ct, n_peaks=nrow(pk),
                              n_chromhmm_active=sum(pk$in_active),
                              pct_chromhmm_active=round(100*mean(pk$in_active), 3),
                              n_h3k27ac=sum(pk$in_h3k27ac),
                              pct_h3k27ac=round(100*mean(pk$in_h3k27ac), 3),
                              n_h3k4me3=sum(pk$in_h3k4me3),
                              pct_h3k4me3=round(100*mean(pk$in_h3k4me3), 3))
    state_break <- pk[, .N, by=state][order(-N)]
    state_break[, pct := round(100*N/sum(N), 3)]
    state_break[, cell_type := ct]
    list(summary=summary_row, per_state=state_break)
})

peak_summary <- rbindlist(lapply(ct_results, `[[`, "summary"))
peak_per_state <- rbindlist(lapply(ct_results, `[[`, "per_state"))
cat("\n[B4] Per-CT peak x chromHMM summary:\n"); print(peak_summary)
fwrite(peak_summary, file.path(OUT_DIR, "peak_chromhmm_overlap.csv"))
fwrite(peak_per_state, file.path(OUT_DIR, "peak_chromhmm_per_state.csv"))

# -----------------------------------------------------------------------------
# B4 sub: hepatocyte DA peaks (disease-up vs disease-down) in Roadmap active
# -----------------------------------------------------------------------------
cat("\n--- B4-Roadmap-addon: Hepatocyte DA up vs down in chromHMM ---\n")
da <- fread(DA_HEP)
cat(sprintf("Hepatocyte DA peaks (Cas13): %d\n", nrow(da)))
# Format: chr1:6785000-6785500
parse_peak <- function(x) {
    m <- regmatches(x, regexec("^(chr[0-9XYM]+):([0-9]+)-([0-9]+)$", x))
    data.table(chrom=sapply(m, function(z) z[2]),
               start=as.integer(sapply(m, function(z) z[3])),
               end=as.integer(sapply(m, function(z) z[4])))
}
da_parsed <- parse_peak(da[["feature name"]])
da <- cbind(da, da_parsed)
# Drop unparsed
da <- da[!is.na(chrom)]
cat(sprintf("After coord parsing: %d\n", nrow(da)))

# DA categorization: significant + direction
da[, sig := `adjusted p-value` < 0.05]
da[, direction := fifelse(`log2(fold_change)` > 0, "up",
                          fifelse(`log2(fold_change)` < 0, "down", "ns"))]
cat("DA distribution:\n"); print(da[, .N, by=.(sig, direction)])

da_gr <- GRanges(seqnames=da$chrom,
                 ranges=IRanges(start=da$start+1L, end=da$end))
h <- findOverlaps(da_gr, chmm_gr)
da_state <- rep("NoChromHMM", length(da_gr))
if (length(h) > 0) {
    ov <- pintersect(da_gr[queryHits(h)], chmm_gr[subjectHits(h)])
    dt2 <- data.table(q=queryHits(h), s=chmm_gr$state[subjectHits(h)],
                      w=width(ov))
    best2 <- dt2[, .SD[which.max(w)], by=q]
    da_state[best2$q] <- best2$s
}
da[, state := da_state]
da[, in_active := state %in% ACTIVE_STATES]
da[, in_h3k27ac := overlapsAny(da_gr, h3k27ac_gr)]
da[, in_h3k4me3 := overlapsAny(da_gr, h3k4me3_gr)]

# Distribution per direction
ha_dist <- da[sig == TRUE, .(n_peaks=.N,
                              n_active=sum(in_active),
                              pct_active=round(100*mean(in_active), 3),
                              n_h3k27ac=sum(in_h3k27ac),
                              pct_h3k27ac=round(100*mean(in_h3k27ac), 3),
                              n_h3k4me3=sum(in_h3k4me3),
                              pct_h3k4me3=round(100*mean(in_h3k4me3), 3)),
                by=direction]
cat("\n[B4] Hepatocyte DA in Roadmap by direction:\n"); print(ha_dist)
fwrite(ha_dist, file.path(OUT_DIR, "hep_da_chromhmm_distribution.csv"))

# Per-state breakdown for hep DA
hep_per_state <- da[sig == TRUE, .N, by=.(direction, state)][order(direction, -N)]
hep_per_state[, pct_within_dir := round(100*N/sum(N), 3), by=direction]
fwrite(hep_per_state, file.path(OUT_DIR, "hep_da_chromhmm_per_state.csv"))

# Fisher: up vs down enrichment in active states
if (nrow(da[sig==TRUE & direction=="up"]) > 0 &&
    nrow(da[sig==TRUE & direction=="down"]) > 0) {
    mtx <- matrix(c(
        sum(da$sig & da$direction=="up"   & da$in_active),
        sum(da$sig & da$direction=="up"   & !da$in_active),
        sum(da$sig & da$direction=="down" & da$in_active),
        sum(da$sig & da$direction=="down" & !da$in_active)
    ), nrow=2, byrow=TRUE,
        dimnames=list(c("up","down"), c("active","not_active")))
    cat("\nHep DA up vs down x active state:\n"); print(mtx)
    ft_da <- fisher.test(mtx)
    cat(sprintf("Fisher OR (up vs down in active): %.3f (95%% CI %.3f-%.3f), p = %.3g\n",
                ft_da$estimate, ft_da$conf.int[1], ft_da$conf.int[2], ft_da$p.value))
}

cat("\n================================================================\n")
cat(sprintf("[%s] Script 56h done. Outputs in %s\n",
            format(Sys.time()), OUT_DIR))
cat("================================================================\n")
