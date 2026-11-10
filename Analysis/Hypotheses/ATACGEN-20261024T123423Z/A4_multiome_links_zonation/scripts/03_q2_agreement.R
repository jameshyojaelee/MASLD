#!/usr/bin/env Rscript
# A4 step 3 (Q2): do stable peak-gene links from credible-set variant peaks name the
# colocalized eQTL gene or the nearest gene more often than matched random genes?
#   Locus: work/loci.tsv (lead-member single linkage, gap > 1 Mb).
#   L: stably linked genes (Q2 family, q < 0.05, z > 0, bootstrap reproduction >= 0.8)
#      over the locus's variant peaks. One representative link per (locus, gene): highest z.
#   C: COLOC r2 tier-1/2 genes, SuSiE eligible PP.H4 > 0.5 (hit1) or ABF PP.H4 > 0.5
#      (top_snp) within the locus member span +/- 500 kb (hg19); testable if tested in the
#      lineage with TSS within 500 kb of a variant peak of the locus (hg38).
#   N: nearest tested gene TSS to each variant peak centre (hg38).
#   A1 = #loci with L and C sharing a gene (loci with |L| >= 1 and testable C);
#   A2 = #loci with L and N sharing a gene (loci with |L| >= 1); A3 = C vs N, descriptive.
#   Null: 10,000 draws, each representative link's gene replaced by a random tested gene in
#   the same peak's +/- 500 kb window and TSS-distance bin; p = (1 + #null >= obs) / 10,001.
#   BH over {A1, A2} per lineage and analysis.
# Analyses: primary (PIP sum >= 0.1, SuSiE + ABF coloc), sens_any_member, sens_susie_only.
# Seed: 20261024 + 9000 + 10 * lineage index + analysis index.
# Usage: 03_q2_agreement.R <A4_DIR> <lineages,comma> <N_DRAW>
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })
args <- commandArgs(TRUE)
A4 <- args[1]; LINS <- strsplit(args[2], ",")[[1]]; N_DRAW <- as.integer(args[3])
SEED <- 20261024L; WIN <- 5e5; STABLE <- 0.8; PIP_MIN <- 0.1
DBIN <- c(0, 1e4, 5e4, 1e5, 2.5e5, 5e5 + 1)
LIN_INDEX <- c(hepatocyte = 1L, fibroblast = 2L, macrophage = 3L)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CO <- file.path(ROOT, "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z")
TIER <- file.path(ROOT, "GWAS/finemapping/config/gwas_trait_tier.tsv")
W <- file.path(A4, "work"); RES <- file.path(A4, "results")

loci <- fread(file.path(W, "loci.tsv"), colClasses = list(character = "chrom"))
vpc <- fread(file.path(W, "variant_peak_cs.tsv"))
pk <- fread(file.path(W, "peaks_annot.tsv"))
genes <- fread(file.path(W, "genes_tss.tsv"))
lsum <- fread(file.path(RES, "link_summary.tsv"))

tier <- fread(TIER); t12 <- tier[tier %in% 1:2]
sp <- fread(file.path(CO, "susie_coloc_signal_pairs.tsv.gz"))
sp <- sp[as.character(eligible) %in% c("TRUE", "True", "true") & PP.H4 > 0.5 & gwas %in% t12$study_name]
sp[, c("chrom", "pos") := tstrsplit(hit1, ":", fixed = TRUE)]
ab <- fread(file.path(CO, "susie_coloc_all_gwas.csv"), select = c("gwas_name", "gene", "ensembl", "PP.H4.abf", "top_snp"))
ab <- ab[PP.H4.abf > 0.5 & gwas_name %in% t12$study_name & !is.na(top_snp) & top_snp != ""]
ab[, c("chrom", "pos") := tstrsplit(as.character(top_snp), ":", fixed = TRUE)]
coloc <- rbind(sp[, .(ensembl, gene, study = gwas, src = "susie", chrom, pos = as.numeric(pos))],
               ab[, .(ensembl, gene, study = gwas_name, src = "abf", chrom, pos = as.numeric(pos))])
loc_c <- loci[, {
  w <- coloc[chrom == .BY$chrom & pos >= span_start_hg19 - WIN & pos <= span_end_hg19 + WIN]
  .(ensembl = w$ensembl, src = w$src)
}, by = .(locus_id, chrom)]
loc_c <- unique(loc_c[, .(locus_id, ensembl, susie = src == "susie")])[, .(susie = any(susie)), by = .(locus_id, ensembl)]

agree_rows <- list(); locus_rows <- list(); lnc_rows <- list()
for (lin in LINS) {
  k <- LIN_INDEX[[lin]]
  info <- lsum[lineage == lin]
  if (!nrow(info) || info$status != "tested" || !isTRUE(info$informative)) {
    agree_rows[[lin]] <- data.table(lineage = lin, analysis = "all", statistic = c("A1", "A2", "A3"),
                                    note = "not run: lineage uninformative (stop rule S1) or untested")
    next
  }
  tg <- fread(file.path(W, paste0("tested_genes_", lin, ".tsv")))
  tpk <- fread(file.path(W, paste0("tested_peaks_", lin, ".tsv")))
  q2 <- fread(file.path(RES, paste0("q2_pairs_", lin, ".tsv.gz")))
  tpk <- merge(tpk[, .(peak_index)], pk[, .(peak_index, peak_id, chr, centre)], by = "peak_index")
  for (a in c("primary", "sens_any_member", "sens_susie_only")) {
    ai <- match(a, c("primary", "sens_any_member", "sens_susie_only"))
    anymem <- a == "sens_any_member"
    lp <- unique(vpc[peak_index %in% tpk$peak_index & (anymem | pip_sum >= PIP_MIN), .(locus_id, peak_index)])
    lp <- merge(lp, tpk, by = "peak_index")
    if (anymem) q2[, `:=`(lk = link_q2any %in% TRUE, rp = repro_q2any)] else q2[, `:=`(lk = link_q2 %in% TRUE, rp = repro_q2)]
    st <- q2[lk & rp >= STABLE, .(peak_index, gene_index, ensembl, symbol, dist, z, q = if (anymem) q2any else q2, rp)]
    # L: representative link per (locus, gene)
    Lr <- merge(lp[, .(locus_id, peak_index, chr, centre)], st, by = "peak_index", allow.cartesian = TRUE)
    Lr <- Lr[order(locus_id, gene_index, -z)][, .SD[1], by = .(locus_id, gene_index)]
    # C (testable)
    cc <- loc_c[if (a == "sens_susie_only") susie else TRUE]
    cc <- merge(cc, tg[, .(ensembl, gene_index, symbol, chr_g = chr, tss)], by = "ensembl")
    cc <- merge(cc, lp[, .(locus_id, chr, centre)], by = "locus_id", allow.cartesian = TRUE)
    cc <- unique(cc[chr == chr_g & abs(centre - tss) <= WIN, .(locus_id, gene_index, symbol)])
    call_all <- unique(merge(loc_c[if (a == "sens_susie_only") susie else TRUE], genes[, .(ensembl, symbol)], by = "ensembl", all.x = TRUE)[
      , .(locus_id, csym = fifelse(is.na(symbol), ensembl, symbol))])
    # N: nearest tested gene per variant peak
    gtr <- GRanges(tg$chr, IRanges(tg$tss, width = 1L))
    pgr <- GRanges(lp$chr, IRanges(lp$centre, width = 1L))
    nn <- nearest(pgr, gtr, ignore.strand = TRUE)
    Nt <- unique(data.table(locus_id = lp$locus_id, gene_index = tg$gene_index[nn])[!is.na(gene_index)])
    pcg <- genes[gene_type == "protein_coding" & has_tss == TRUE]
    nn2 <- nearest(pgr, GRanges(pcg$chr, IRanges(pcg$tss, width = 1L)), ignore.strand = TRUE)
    Npc <- unique(data.table(locus_id = lp$locus_id, sym = pcg$symbol[nn2])[!is.na(sym)])

    lociL <- unique(Lr$locus_id); lociC <- unique(cc$locus_id)
    inC <- Lr[, paste(locus_id, gene_index)] %in% cc[, paste(locus_id, gene_index)]
    inN <- Lr[, paste(locus_id, gene_index)] %in% Nt[, paste(locus_id, gene_index)]
    e1 <- intersect(lociL, lociC)
    A1 <- uniqueN(Lr$locus_id[inC & Lr$locus_id %in% e1]); A2 <- uniqueN(Lr$locus_id[inN])
    A3 <- uniqueN(cc[paste(locus_id, gene_index) %in% Nt[, paste(locus_id, gene_index)], locus_id])

    # matched null
    set.seed(SEED + 9000L + 10L * k + ai)
    tg_chr <- split(seq_len(nrow(tg)), tg$chr)
    pc_ <- tpk$centre[match(Lr$peak_index, tpk$peak_index)]; pch <- tpk$chr[match(Lr$peak_index, tpk$peak_index)]
    D <- matrix(NA_integer_, nrow(Lr), N_DRAW)
    for (j in seq_len(nrow(Lr))) {
      ii <- tg_chr[[pch[j]]]; dd <- abs(tg$tss[ii] - pc_[j])
      b <- findInterval(abs(Lr$dist[j]), DBIN)
      cand <- tg$gene_index[ii][dd <= WIN & findInterval(dd, DBIN) == b]
      D[j, ] <- if (length(cand) == 1L) cand else sample(cand, N_DRAW, replace = TRUE)
    }
    keyC <- cc[, paste(locus_id, gene_index)]; keyN <- Nt[, paste(locus_id, gene_index)]
    hitC <- matrix(paste(Lr$locus_id, D) %in% keyC, nrow(Lr)); hitN <- matrix(paste(Lr$locus_id, D) %in% keyN, nrow(Lr))
    lf <- factor(Lr$locus_id)
    null_count <- function(H, eligible) {
      m <- rowsum(H * 1L, lf) > 0
      colSums(m[rownames(m) %in% eligible, , drop = FALSE])
    }
    if (nrow(Lr)) { n1 <- null_count(hitC, e1); n2 <- null_count(hitN, lociL) } else { n1 <- n2 <- integer(N_DRAW) }
    pA1 <- (1 + sum(n1 >= A1)) / (1 + N_DRAW); pA2 <- (1 + sum(n2 >= A2)) / (1 + N_DRAW)
    qq <- p.adjust(c(pA1, pA2), "BH")
    agree_rows[[paste(lin, a)]] <- data.table(lineage = lin, analysis = a, statistic = c("A1", "A2", "A3"),
      eligible_loci = c(length(e1), length(lociL), length(lociC)), observed = c(A1, A2, A3),
      observed_frac = c(A1 / length(e1), A2 / length(lociL), A3 / length(lociC)),
      null_mean = c(mean(n1), mean(n2), NA), null_q025 = c(quantile(n1, 0.025), quantile(n2, 0.025), NA),
      null_q975 = c(quantile(n1, 0.975), quantile(n2, 0.975), NA),
      p_empirical = c(pA1, pA2, NA), q_bh = c(qq, NA), n_draw = N_DRAW,
      loci_with_tested_variant_peak = uniqueN(lp$locus_id), stable_links = nrow(st), representative_links = nrow(Lr),
      note = c("", "", "descriptive, no null"))

    if (a == "primary") {
      sy <- function(x) paste(sort(unique(x)), collapse = ";")
      gsym <- setNames(tg$symbol, tg$gene_index)
      lr_s <- loci[locus_id %in% unique(lp$locus_id), .(locus_id, chrom, span_start_hg19, span_end_hg19, scope, n_cs)]
      lr_s <- merge(lr_s, lp[, .(variant_peaks = uniqueN(peak_index)), by = locus_id], by = "locus_id")
      lr_s <- merge(lr_s, Lr[, .(L_genes = sy(symbol)), by = locus_id], by = "locus_id", all.x = TRUE)
      lr_s <- merge(lr_s, cc[, .(C_testable = sy(symbol)), by = locus_id], by = "locus_id", all.x = TRUE)
      lr_s <- merge(lr_s, call_all[, .(C_all = sy(csym)), by = locus_id], by = "locus_id", all.x = TRUE)
      lr_s <- merge(lr_s, Nt[, .(N_tested = sy(gsym[as.character(gene_index)])), by = locus_id], by = "locus_id", all.x = TRUE)
      lr_s <- merge(lr_s, Npc[, .(N_protein_coding = sy(sym)), by = locus_id], by = "locus_id", all.x = TRUE)
      lr_s[, `:=`(L_in_C = locus_id %in% Lr$locus_id[inC], L_in_N = locus_id %in% Lr$locus_id[inN],
                  C_in_N = locus_id %in% cc[paste(locus_id, gene_index) %in% keyN, locus_id])]
      locus_rows[[lin]] <- cbind(lineage = lin, lr_s)
      lnc <- Lr[!inC]
      lnc <- merge(lnc, lr_s[, .(locus_id, scope, C_all, C_testable, N_tested)], by = "locus_id")
      lnc[, `:=`(is_nearest_tested = paste(locus_id, gene_index) %in% keyN,
                 top_variant = vpc[match(paste(lnc$peak_index, lnc$locus_id), paste(vpc$peak_index, vpc$locus_id)), top_variant])]
      lnc_rows[[lin]] <- cbind(lineage = lin, lnc[, .(locus_id, scope, symbol, ensembl, peak_index,
        peak_id = pk$peak_id[match(peak_index, pk$peak_index)], top_variant, dist, z, q, bootstrap_reproduction = rp,
        is_nearest_tested, C_all, C_testable)])
    }
  }
}
fwrite(rbindlist(agree_rows, fill = TRUE), file.path(RES, "q2_agreement.tsv"), sep = "\t")
fwrite(rbindlist(locus_rows, fill = TRUE), file.path(RES, "q2_loci.tsv"), sep = "\t")
fwrite(rbindlist(lnc_rows, fill = TRUE), file.path(RES, "q2_link_not_coloc.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(A4, "sessionInfo_03.txt"))
message("done")
