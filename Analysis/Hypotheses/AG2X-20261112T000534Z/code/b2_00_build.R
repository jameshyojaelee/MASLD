#!/usr/bin/env Rscript
# AG2X Stage B2 step 0: members of the credible sets Stage B left out, matched background, counts. No hosted calls.
# Read-only on all inputs; writes only stageB2/inputs/ and stageB2/tmp/.
# Target loci: Stage B loci (stageB/inputs/loci.tsv, B0 job 23229228) with scored_as ==
#   "cs_outside_substrate_lead_proxies": 51 tier 2 (DIG1 ALT) loci and 502 tier 3 loci. Their member rows come
#   from stageB/inputs/credible_set_member_rows.tsv.gz (same set and locus assignment as Stage B).
# Members: hg19 -> hg38 liftOver (unique hit, + strand); REF = GRCh38 FASTA bases at the hg38 position with the
#   length of one allele; the allele that matches is REF. SNVs and indels; failures excluded and counted.
# Member ancestry: ancestry of the member's highest-PIP set at the locus (ties: first by set ID).
# Background per locus: one per member (1:1), matched on MAF in the member's own ancestry panel
#   (1000 Genomes AFR, AMR, EAS, SAS, EUR; panel_af/<anc>); pool 1 = combined_finemapping SNVs in the locus
#   windows with max PIP < 0.001 over all studies, not at a credible-set position, present in that ancestry
#   panel; pool 2 (only when pool 1 is exhausted) = that ancestry panel's SNVs in the window not at a
#   PIP >= 0.001 position or credible-set position. Nearest MAF without replacement, seeded tiebreak.
# Variants already scored in Stage B (stageB/inputs/variants_scored_final.tsv) are not re-called; their Stage B
#   scores are reused.
suppressPackageStartupMessages({library(data.table); library(rtracklayer); library(GenomicRanges); library(Rsamtools)})
set.seed(20261111)
RP <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
PK <- file.path(RP, "Analysis/Hypotheses/AG2X-20261112T000534Z")
O <- file.path(PK, "stageB2/inputs")
SB <- file.path(PK, "stageB/inputs")
FM <- file.path(RP, "GWAS/finemapping/results")
FA <- file.path(RP, "Analysis/MASLD_Model_Benchmark/executions/sequence-fasta-21063808/GRCh38.p14.sequence_model.fa")
CHAIN <- file.path(RP, "data/broadaway_eqtl/hg19ToHg38.over.chain")
LDREF <- file.path(RP, "GWAS/finemapping/data/ld_ref")
WIN <- 500000L
ANC <- c(EUR = "eur", EAS = "eas", AFR = "afr", AMR = "amr", SAS = "sas")
cnt <- list()
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n"); flush.console() }
akey <- function(c, p, a, b) paste(c, p, pmin(a, b), pmax(a, b), sep = ":")
chain <- import.chain(CHAIN)
fa <- FaFile(FA)
lift <- function(chr, pos) {
  gr <- GRanges(paste0("chr", chr), IRanges(pos, pos), strand = "+")
  lo <- liftOver(gr, chain)
  n <- elementNROWS(lo)
  out <- data.table(hits = n, c38 = NA_character_, p38 = NA_integer_, s38 = NA_character_)
  u <- unlist(lo[n == 1]); idx <- which(n == 1)
  out[idx, `:=`(c38 = as.character(seqnames(u)), p38 = start(u), s38 = as.character(strand(u)))]
  out
}
seqat <- function(c38, p38, len) {
  r <- rep(NA_character_, length(c38)); ok <- !is.na(p38)
  if (any(ok)) r[ok] <- as.character(getSeq(fa, GRanges(c38[ok], IRanges(p38[ok], p38[ok] + len[ok] - 1L))))
  toupper(r)
}

loci <- fread(file.path(SB, "loci.tsv"))
tl <- loci[scored_as == "cs_outside_substrate_lead_proxies"]
cnt$target_loci_by_tier <- as.list(tl[, .N, by = tier][, setNames(N, tier)])
mem <- fread(file.path(SB, "credible_set_member_rows.tsv.gz"))
mem <- mem[locus_id %in% tl$locus_id]
cnt$member_rows <- nrow(mem); cnt$member_keys <- uniqueN(mem$k)
cnt$member_rows_with_uid_already <- sum(!is.na(mem$uid) & mem$uid != "")
cnt$member_keys_by_ancestry <- as.list(mem[, .(n = uniqueN(k)), by = ancestry][, setNames(n, ancestry)])

# ---------------- lift and REF-check every member key
kk <- unique(mem[, .(k)])
sp <- tstrsplit(kk$k, ":")
kk[, `:=`(chr = as.integer(sp[[1]]), pos19 = as.integer(sp[[2]]), a = toupper(sp[[3]]), b = toupper(sp[[4]]))]
lo <- lift(kk$chr, kk$pos19)
kk[, `:=`(c38 = lo$c38, p38 = lo$p38, s38 = lo$s38, hits = lo$hits)]
kk[, lift_ok := hits == 1 & s38 %in% "+"]
kk[, sa := seqat(c38, ifelse(lift_ok, p38, NA_integer_), nchar(a))]
kk[, sb := seqat(c38, ifelse(lift_ok, p38, NA_integer_), nchar(b))]
kk[, ref := fifelse(lift_ok & sa == a, a, fifelse(lift_ok & sb == b, b, NA_character_))]
kk[, alt := fifelse(!is.na(ref) & ref == a, b, fifelse(!is.na(ref), a, NA_character_))]
kk[, uid := fifelse(is.na(ref), NA_character_, paste0(c38, ":", p38, ":", ref, ":", alt))]
cnt$member_keys_lift_fail <- sum(!kk$lift_ok, na.rm = TRUE) + sum(is.na(kk$lift_ok))
cnt$member_keys_ref_fail <- sum(kk$lift_ok %in% TRUE & is.na(kk$uid))
cnt$member_keys_ok <- sum(!is.na(kk$uid))
cnt$member_keys_ok_indel <- sum(!is.na(kk$uid) & (nchar(kk$a) > 1 | nchar(kk$b) > 1))
mem[, uid := NULL]
mem <- merge(mem, kk[, .(k, uid)], by = "k")
memok <- mem[!is.na(uid)]

# member table: primary locus (first by locus ID), ancestry of the highest-PIP set at that locus, max PIP
setorder(memok, uid, locus_id, -pip, set_id)
mvar <- memok[, .(primary_locus = locus_id[1], loci = paste(sort(unique(locus_id)), collapse = ";"),
                  tier = min(tier), max_pip = max(pip, na.rm = TRUE), k = k[1],
                  studies = paste(sort(unique(study)), collapse = ","),
                  ancestries = paste(sort(unique(ancestry)), collapse = ",")), by = uid]
an <- memok[, .SD[1], by = .(uid, locus_id)][, .(uid, locus_id, member_ancestry = ancestry)]
mvar <- merge(mvar, an, by.x = c("uid", "primary_locus"), by.y = c("uid", "locus_id"))
mvar[, panel := ifelse(member_ancestry %in% names(ANC), ANC[member_ancestry], "eur")]
cnt$members_panel_fallback_eur <- sum(!(mvar$member_ancestry %in% names(ANC)))
cnt$member_ancestry_values <- as.list(table(mvar$member_ancestry))
sbv <- fread(file.path(PK, "stageB/inputs/variants_scored_final.tsv"))
mvar[, scored_in_stageB := uid %in% sbv$uid]
cnt$member_variants <- nrow(mvar)
cnt$member_variants_scored_in_stageB <- sum(mvar$scored_in_stageB)
cnt$member_variants_by_tier <- as.list(mvar[, .N, by = tier][, setNames(N, tier)])
cnt$member_variants_by_ancestry <- as.list(mvar[, .N, by = member_ancestry][, setNames(N, member_ancestry)])
cons <- fread(file.path(FM, "credible_set_variant_consequences.csv"))
cons[, k := akey(chr, pos, ref, alt)]
mvar <- merge(mvar, unique(cons[, .(k, Consequence, class)], by = "k"), by = "k", all.x = TRUE)
mvar[, coding_pip05 := class %in% "coding_protein_altering" & max_pip >= 0.5]
cnt$members_with_consequence <- sum(!is.na(mvar$class))

# ---------------- background candidates
fmr <- fread(file.path(FM, "combined_finemapping.csv"),
             select = c("chromosome", "position", "allele1", "allele2", "study", "locus", "max_pip"))
say("combined_finemapping rows %d", nrow(fmr))
wins <- fmr[, .(wchr = chromosome[1], wmin = min(position), wmax = max(position)), by = .(study, locus)]
pmx <- fmr[, .(pmax = max(max_pip, na.rm = TRUE)), by = .(chromosome, position)]
cand <- unique(fmr[nchar(allele1) == 1 & nchar(allele2) == 1, .(chromosome, position, allele1, allele2)])
rm(fmr); gc()
cand[, npos := .N, by = .(chromosome, position)]
cand <- cand[npos == 1][, npos := NULL]
cand <- merge(cand, pmx[is.finite(pmax) & pmax < 0.001, .(chromosome, position)], by = c("chromosome", "position"))
csall <- fread(file.path(FM, "credible_sets.csv"), select = c("chromosome", "position", "either_in_cs"))
csk <- unique(csall[either_in_cs == TRUE, .(chromosome, position)]); rm(csall)
cand <- cand[!csk, on = .(chromosome, position)]
cand[, k := akey(chromosome, position, allele1, allele2)]
sl <- unique(mem[, .(locus_id, study, locus)])
lw <- merge(sl, wins, by = c("study", "locus"))[, .(wchr = wchr[1], wmin = min(wmin), wmax = max(wmax)), by = locus_id]
lw <- rbind(lw, tl[!(locus_id %in% lw$locus_id), .(locus_id, wchr = chr, wmin = pos19 - WIN, wmax = pos19 + WIN)])
cc1 <- cand[lw, on = .(chromosome = wchr, position >= wmin, position <= wmax), nomatch = 0,
            .(locus_id, chromosome = x.chromosome, position = x.position, allele1, allele2, k)]
pm_bad <- pmx[!is.finite(pmax) | pmax >= 0.001, .(chromosome, position)]
rm(cand, pmx); gc()

refcheck <- function(d) {
  lo <- lift(d$chromosome, d$position)
  d[, `:=`(c38 = lo$c38, p38 = lo$p38, ok = lo$hits == 1 & lo$s38 %in% "+")]
  d <- d[ok == TRUE]
  d[, ref38 := seqat(c38, p38, rep(1L, .N))]
  d <- d[ref38 == allele1 | ref38 == allele2]
  d[, uid := paste0(c38, ":", p38, ":", ref38, ":", ifelse(ref38 == allele1, allele2, allele1))]
  d[!(uid %in% mvar$uid)]
}
bg <- list(); pool_n <- list()
for (A in sort(unique(mvar$panel))) {
  lwA0 <- lw[locus_id %in% mvar[panel == A]$primary_locus]
  af <- rbindlist(lapply(sort(unique(lwA0$wchr)), function(ch) {   # keep only panel rows inside this panel's windows
    x <- fread(file.path(LDREF, sprintf("panel_af/%s/chr%d.af.tsv.gz", A, ch)),
               select = c("chromosome", "position", "bim_a1", "bim_a2", "af_a1"))
    w <- lwA0[wchr == ch]
    x[unique(x[w, on = .(position >= wmin, position <= wmax), nomatch = 0, which = TRUE])]
  }))
  af[, k := akey(chromosome, position, bim_a1, bim_a2)]
  af[, maf := pmin(af_a1, 1 - af_a1)]
  af <- unique(af, by = "k")
  need <- mvar[panel == A]
  need <- merge(need, af[, .(k, maf)], by = "k", all.x = TRUE)
  lwA <- lw[locus_id %in% need$primary_locus]
  p1 <- merge(cc1[locus_id %in% lwA$locus_id], af[, .(k, maf)], by = "k")
  p1 <- refcheck(p1); p1[, src := 1L]
  p2 <- af[nchar(bim_a1) == 1 & nchar(bim_a2) == 1, .(chromosome, position, allele1 = bim_a1, allele2 = bim_a2, maf, k)]
  p2[, npos := .N, by = .(chromosome, position)]
  p2 <- p2[npos == 1][, npos := NULL][!csk, on = .(chromosome, position)][!pm_bad, on = .(chromosome, position)]
  p2 <- p2[lwA, on = .(chromosome = wchr, position >= wmin, position <= wmax), nomatch = 0,
           .(locus_id, chromosome = x.chromosome, position = x.position, allele1, allele2, maf, k)]
  p2 <- p2[!p1, on = .(locus_id, chromosome, position)]
  rm(af); gc()
  p2 <- refcheck(p2); p2[, src := 2L]
  pool <- rbind(p1, p2, fill = TRUE)
  setorder(pool, locus_id, chromosome, position)
  pool[, rk := runif(.N)]
  pool_n[[A]] <- list(pool1 = nrow(p1), pool2 = nrow(p2))
  for (L in sort(unique(need$primary_locus))) {
    m <- need[primary_locus == L][order(uid)]
    pl <- pool[locus_id == L]
    used <- character(0)
    for (j in seq_len(nrow(m))) {
      p <- pl[!(uid %in% used)]
      if (!nrow(p)) break
      if (any(p$src == 1L)) p <- p[src == 1L]
      dm <- if (is.na(m$maf[j])) p$rk else abs(p$maf - m$maf[j])
      pick <- p[order(dm, rk)][1]
      used <- c(used, pick$uid)
      bg[[length(bg) + 1]] <- data.table(locus_id = L, member_uid = m$uid[j], uid = pick$uid, panel = A,
                                         member_maf = m$maf[j], bg_maf = pick$maf,
                                         maf_diff = abs(pick$maf - m$maf[j]), bg_source = pick$src)
    }
  }
  say("panel %s: %d members, background so far %d", A, nrow(need), length(bg))
}
bg <- rbindlist(bg)
bg <- merge(bg, tl[, .(locus_id, tier)], by = "locus_id")
cnt$background_pools <- pool_n
cnt$background_rows <- nrow(bg); cnt$background_unique <- uniqueN(bg$uid)
cnt$background_by_panel <- as.list(bg[, .N, by = panel][, setNames(N, panel)])
cnt$background_by_source <- as.list(table(bg$bg_source))
cnt$background_median_maf_diff <- median(bg$maf_diff, na.rm = TRUE)
cnt$members_without_background <- nrow(mvar) - nrow(bg)
cnt$members_without_af_in_panel <- sum(is.na(bg$member_maf)) + cnt$members_without_background

# ---------------- variants to call (new uids only) and order: tier 2 (ALT) loci first, then tier 3
v <- rbind(mvar[, .(uid, role = "member", tier, locus_id = primary_locus)],
           bg[, .(uid, role = "background", tier, locus_id)])
vu <- v[, .(tier = min(tier), roles = paste(sort(unique(role)), collapse = ";"),
            loci = paste(sort(unique(locus_id)), collapse = ";")), by = uid]
vu[, scored_in_stageB := uid %in% sbv$uid]
cnt$variants_total <- nrow(vu)
cnt$variants_new_by_tier <- as.list(vu[scored_in_stageB == FALSE, .N, by = tier][, setNames(N, tier)])
cnt$variants_reused_from_stageB <- sum(vu$scored_in_stageB)
fwrite(tl, file.path(O, "loci_b2.tsv"), sep = "\t")
fwrite(memok, file.path(O, "credible_set_member_rows_b2.tsv.gz"), sep = "\t")
fwrite(kk, file.path(O, "member_key_liftover.tsv.gz"), sep = "\t")
fwrite(mvar, file.path(O, "member_variants_b2.tsv"), sep = "\t")
fwrite(bg, file.path(O, "background_b2.tsv"), sep = "\t")
fwrite(vu, file.path(O, "variants_b2.tsv"), sep = "\t")
write(jsonlite::toJSON(cnt, auto_unbox = TRUE, pretty = TRUE, na = "null"), file.path(O, "b2_counts.json"))
capture.output(sessionInfo(), file = file.path(O, "sessionInfo_b2_build.txt"))
say("done")
