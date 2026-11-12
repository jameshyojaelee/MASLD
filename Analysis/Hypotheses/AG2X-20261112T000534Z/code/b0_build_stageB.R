#!/usr/bin/env Rscript
# AG2X Stage B step B0: variant sets, tiers, background and budget. No hosted calls.
# Read-only on all inputs; writes only stageB/inputs/ and stageB/tmp/.
# Rules (to be fixed in PRESPEC_AG2X_stageB.md; this script only builds sets and counts):
# - Credible-set rows: GWAS/finemapping/results/credible_sets.csv, rows with either_in_cs == TRUE.
#   A set is study|trait|locus|susie<cs> (susie_cs > 0, PIP = susie_pip) or study|trait|locus|carma<cs>
#   (carma_cs > 0, PIP = carma_pip). Variants map to the REF-checked hg38 substrate
#   (seqfunc/variant_substrate_hg38.tsv) by hg19 chr, pos and unordered allele pair.
# - Locus assignment: a set belongs to the DIG1 locus (169 ALT loci, GRCh37 lead) whose lead lies within 500 kb
#   of the set's top-PIP member; nearest lead if several. Extra named loci: PNPLA3 lead 22:44324727 (DIG1
#   script 73 benchmark constant); named genes not found in DIG1 gene lists are looked up in the substrate gene
#   columns, then by GENCODE gene start (hg38) +-500 kb.
# - Tier 1: DIG1 loci with lfsr < 0.05 for liver fat, Ghouse cirrhosis or FinnGen NAFLD in ENZFAT1 S2, plus loci
#   of TNFSF10, PANX1, GPT (ENZFAT1 locus 68), GPAM, HSD17B13, MTARC1, PNPLA3, TM6SF2 (ENZFAT1 locus 148).
#   Tier 2: other DIG1 loci. Tier 3: all other sets, merged into loci by single linkage of top positions within
#   500 kb per chromosome.
# - Loci without a set: lead plus r2 > 0.8 proxies within 500 kb in the 1000 Genomes EUR panel used for
#   fine-mapping (GWAS/finemapping/data/ld_ref/1kg_eur, PLINK 1.9).
# - Background per locus: as many variants as the locus has scored members; candidates are biallelic SNVs in
#   the locus fine-mapping windows (combined_finemapping.csv rows of the locus's study x locus pairs; lead +-500 kb
#   for loci without a set), maximum PIP over all studies < 0.001, never a credible-set member, present in the
#   1KG EUR AF panel, unique hg19->hg38 liftover on +, hg38 REF equal to one allele. Matched to members on EUR MAF,
#   nearest without replacement, ties by a seeded random key (seed 20261111).
suppressPackageStartupMessages({library(data.table); library(rtracklayer); library(GenomicRanges); library(Rsamtools)})
set.seed(20261111)
RP <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
O <- file.path(RP, "Analysis/Hypotheses/AG2X-20261112T000534Z/stageB/inputs")
TMP <- file.path(RP, "Analysis/Hypotheses/AG2X-20261112T000534Z/stageB/tmp")
FM <- file.path(RP, "GWAS/finemapping/results")
FA <- file.path(RP, "Analysis/MASLD_Model_Benchmark/executions/sequence-fasta-21063808/GRCh38.p14.sequence_model.fa")
CHAIN <- file.path(RP, "data/broadaway_eqtl/hg19ToHg38.over.chain")
LDREF <- file.path(RP, "GWAS/finemapping/data/ld_ref")
PLINK <- "/nfs/sw/plink/plink-1.90b6.24/plink"
GTF <- "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/genes/genes.gtf.gz"
WIN <- 500000L
NAMED <- c("TNFSF10", "PANX1", "GPT", "GPAM", "HSD17B13", "MTARC1", "PNPLA3", "TM6SF2")
# old symbols used in the substrate and consequence tables (team-lead 2026-11-12: MTARC1 is MARC1 on disk)
ALIAS <- list(MTARC1 = c("MTARC1", "MARC1"))
syms <- function(g) if (is.null(ALIAS[[g]])) g else ALIAS[[g]]
cnt <- list()
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n"); flush.console() }
akey <- function(c, p, a, b) paste(c, p, pmin(a, b), pmax(a, b), sep = ":")

# ---------------------------------------------------------------- liftover and REF helpers
chain <- import.chain(CHAIN)
fa <- FaFile(FA)
lift <- function(chr, pos) {
  gr <- GRanges(paste0("chr", chr), IRanges(pos, pos), strand = "+")
  lo <- liftOver(gr, chain)
  n <- elementNROWS(lo)
  out <- data.table(i = seq_along(gr), hits = n, c38 = NA_character_, p38 = NA_integer_, s38 = NA_character_)
  u <- unlist(lo[n == 1]); idx <- which(n == 1)
  out[idx, `:=`(c38 = as.character(seqnames(u)), p38 = start(u), s38 = as.character(strand(u)))]
  out
}
refbase <- function(c38, p38) {
  r <- rep(NA_character_, length(c38)); ok <- !is.na(p38)
  if (any(ok)) r[ok] <- as.character(getSeq(fa, GRanges(c38[ok], IRanges(p38[ok], p38[ok]))))
  toupper(r)
}

# ---------------------------------------------------------------- substrate and credible-set rows
sub <- fread(file.path(FM, "seqfunc/variant_substrate_hg38.tsv"))
sp <- tstrsplit(sub$variant_id_hg19, ":")
sub[, `:=`(chr19 = as.integer(sp[[1]]), pos19 = as.integer(sp[[2]]), a19 = sp[[3]], b19 = sp[[4]])]
sub[, k := akey(chr19, pos19, a19, b19)]
sub[, uid := paste0(chr, ":", pos_hg38, ":", ref_hg38, ":", alt_hg38)]
cnt$substrate_rows <- nrow(sub); cnt$substrate_unique_keys <- uniqueN(sub$k)
cnt$substrate_ref_status <- as.list(table(sub$ref_status, useNA = "ifany"))

cs <- fread(file.path(FM, "credible_sets.csv"),
            select = c("chromosome", "position", "allele1", "allele2", "trait", "locus", "study", "ancestry",
                       "susie_pip", "susie_cs", "carma_pip", "carma_cs", "either_in_cs", "variant_id"))
cnt$cs_rows <- nrow(cs)
cnt$cs_either_in_cs <- as.list(table(cs$either_in_cs, useNA = "ifany"))
cs <- cs[either_in_cs == TRUE]
cs[, k := akey(chromosome, position, allele1, allele2)]
cnt$cs_member_unique_variants <- uniqueN(cs$k)
s1 <- cs[!is.na(susie_cs) & susie_cs > 0, .(k, chromosome, position, trait, locus, study, ancestry,
                                            set_id = paste(study, trait, locus, paste0("susie", susie_cs), sep = "|"),
                                            method = "susie", pip = susie_pip)]
s2 <- cs[!is.na(carma_cs) & carma_cs > 0, .(k, chromosome, position, trait, locus, study, ancestry,
                                            set_id = paste(study, trait, locus, paste0("carma", carma_cs), sep = "|"),
                                            method = "carma", pip = carma_pip)]
mem <- unique(rbind(s1, s2), by = c("set_id", "k"))
cnt$cs_rows_in_no_set <- nrow(cs[!(k %in% mem$k)])
cnt$sets_susie <- uniqueN(s1$set_id); cnt$sets_carma <- uniqueN(s2$set_id)
mem <- merge(mem, unique(sub[, .(k, uid, ref_status, is_indel, gene, primary_gene)], by = "k"), by = "k", all.x = TRUE)
cnt$member_rows_without_substrate_row <- sum(is.na(mem$uid))
cnt$member_variants_without_substrate_row <- uniqueN(mem[is.na(uid)]$k)
tops <- mem[order(set_id, -pip)][, .SD[1], by = set_id][, .(set_id, top_chr = chromosome, top_pos = position, top_pip = pip)]

# ---------------------------------------------------------------- loci
dig <- fread(file.path(RP, "Analysis/Hypotheses/DIG1-20261024T113951Z/out_release_vs_disease_loci/loci.tsv"))
enz <- fread(file.path(RP, "Analysis/Hypotheses/ENZFAT1-20261110T204401Z/out_S2/locus_effects.tsv"))
enz <- enz[!grepl("^ctrl_", locus)][, locus := as.integer(locus)]
enz[, disease_effect := (!is.na(lfsr_pdff) & lfsr_pdff < 0.05) | (!is.na(lfsr_cirr_ghouse) & lfsr_cirr_ghouse < 0.05) |
                        (!is.na(lfsr_nafld_finngen) & lfsr_nafld_finngen < 0.05)]
loci <- dig[, .(locus_id = paste0("DIG1_", locus), dig_locus = locus, chr, pos19 = pos, a1, a2, genes, source = "DIG1")]
loci <- merge(loci, enz[, .(dig_locus = locus, lfsr_pdff, lfsr_cirr_ghouse, lfsr_nafld_finngen, disease_effect)],
              by = "dig_locus", all.x = TRUE)
loci[is.na(disease_effect), disease_effect := FALSE]
cnt$dig1_loci <- nrow(loci); cnt$dig1_disease_effect_loci <- sum(loci$disease_effect)
named_hit <- list()
for (g in NAMED) {
  h <- loci[grepl(paste0("(^|,)", g, "(,|$)"), genes)]$locus_id
  if (g == "GPT") h <- union(h, "DIG1_68")         # ENZFAT1 S2 controls: GPT locus = 68
  if (g == "TM6SF2") h <- union(h, "DIG1_148")     # ENZFAT1: locus 148 is rs58542926
  named_hit[[g]] <- h
}
# PNPLA3 lead (DIG1 script 73 BENCH constant, hg19)
loci <- rbind(loci, data.table(locus_id = "named_PNPLA3", dig_locus = NA_integer_, chr = 22L, pos19 = 44324727L, a1 = "C",
                               a2 = "G", genes = "PNPLA3", source = "named", disease_effect = FALSE), fill = TRUE)
named_hit[["PNPLA3"]] <- union(named_hit[["PNPLA3"]], "named_PNPLA3")

# assign sets to the nearest locus lead within 500 kb
asg <- tops[loci[, .(locus_id, chr, pos19)], on = .(top_chr = chr), allow.cartesian = TRUE, nomatch = 0]
asg[, d := abs(top_pos - pos19)]
asg <- asg[d <= WIN][order(set_id, d)][, .SD[1], by = set_id]
setmap <- asg[, .(set_id, locus_id)]
# named genes still without a locus: substrate gene columns, then GENCODE gene start +-500 kb (hg38)
gtf <- fread(cmd = sprintf("zcat %s | awk -F'\\t' '$3==\"gene\"'", GTF), header = FALSE, sep = "\t", quote = "")
gtf[, gname := sub('.*gene_name "([^"]+)".*', "\\1", V9)]
gtf[, start38 := ifelse(V7 == "+", V4, V5)]
memh38 <- unique(mem[!is.na(uid), .(set_id, uid, gene, primary_gene)])
memh38[, `:=`(c38 = sub(":.*", "", uid), p38 = as.integer(tstrsplit(uid, ":")[[2]]))]
for (g in NAMED) {
  if (length(named_hit[[g]])) next
  ids <- unique(memh38[gene %in% syms(g) | primary_gene %in% syms(g)]$set_id)
  how <- paste("substrate gene column", paste(syms(g), collapse = "/"))
  if (!length(ids)) {
    gg <- gtf[gname == g]
    if (nrow(gg)) ids <- unique(memh38[c38 == gg$V1[1] & abs(p38 - gg$start38[1]) <= WIN]$set_id)
    how <- "GENCODE gene start +-500 kb"
  }
  # sets already assigned to a DIG1 locus: that locus becomes the named locus (as GPT = DIG1_68)
  already <- unique(setmap[set_id %in% ids]$locus_id)
  if (length(already)) {
    named_hit[[g]] <- already
    cnt[[paste0("named_", g, "_already_in")]] <- already
  }
  ids <- setdiff(ids, setmap$set_id)
  cnt[[paste0("named_", g, "_rule")]] <- if (length(ids) || length(already)) how else "not found"
  if (length(ids)) {
    lid <- paste0("named_", g)
    t0 <- tops[set_id %in% ids][order(-top_pip)][1]
    loci <- rbind(loci, data.table(locus_id = lid, chr = t0$top_chr, pos19 = t0$top_pos, genes = g, source = "named",
                                   disease_effect = FALSE), fill = TRUE)
    setmap <- rbind(setmap, data.table(set_id = ids, locus_id = lid))
    named_hit[[g]] <- union(named_hit[[g]], lid)
  }
}
# old-symbol check for every named gene: substrate gene/primary_gene rows under each symbol
cnt$named_symbol_rows_in_substrate <- lapply(setNames(NAMED, NAMED), function(g)
  as.list(sapply(syms(g), function(x) sum(sub$gene == x | sub$primary_gene == x, na.rm = TRUE))))
cnt$named_gene_loci <- lapply(named_hit, function(x) if (length(x)) x else "not found")
tier1 <- union(loci[disease_effect == TRUE]$locus_id, unlist(named_hit))
loci[, tier := ifelse(locus_id %in% tier1, 1L, ifelse(source == "DIG1", 2L, 3L))]
# tier 3 loci: remaining sets, single linkage on top positions within 500 kb
rest <- tops[!(set_id %in% setmap$set_id)][order(top_chr, top_pos)]
rest[, grp := cumsum(c(TRUE, diff(top_pos) > WIN | diff(top_chr) != 0))]
rest[, locus_id := paste0("T3_", top_chr, "_", min(top_pos)), by = grp]
t3 <- rest[, .(chr = top_chr[1], pos19 = top_pos[which.max(top_pip)], source = "tier3_cluster", tier = 3L,
               disease_effect = FALSE), by = locus_id]
loci <- rbind(loci, t3, fill = TRUE)
setmap <- rbind(setmap, rest[, .(set_id, locus_id)])
mem <- merge(mem, setmap, by = "set_id")
mem <- merge(mem, loci[, .(locus_id, tier)], by = "locus_id")
mem[, member_source := fifelse(is.na(uid), NA_character_, "substrate")]
# Tier 1 loci whose credible sets have no substrate member: lift and REF-check those members here (SNVs only)
sub_loci <- unique(mem[!is.na(uid)]$locus_id)
ix <- which(is.na(mem$uid) & mem$tier == 1L & !(mem$locus_id %in% sub_loci))
cnt$tier1_outside_substrate_rows <- length(ix)
if (length(ix)) {
  kk <- tstrsplit(mem$k[ix], ":")
  a <- kk[[3]]; b <- kk[[4]]
  lo <- lift(as.integer(kk[[1]]), as.integer(kk[[2]]))
  r38 <- refbase(lo$c38, lo$p38)
  ok <- lo$hits == 1 & lo$s38 %in% "+" & nchar(a) == 1 & nchar(b) == 1 & (r38 == a | r38 == b)
  ok[is.na(ok)] <- FALSE
  u <- ifelse(ok, paste0(lo$c38, ":", lo$p38, ":", r38, ":", ifelse(r38 == a, b, a)), NA_character_)
  mem[ix, `:=`(uid = u, member_source = ifelse(ok, "lifted_here", NA_character_))]
  cnt$tier1_outside_substrate_lifted_ok_variants <- uniqueN(u[ok])
  cnt$tier1_outside_substrate_failed_variants <- uniqueN(mem$k[ix][!ok])
}
# TRAIL4 ALT credible-set top variants (hg19; members beyond the top were not saved by TRAIL4)
tr <- fread(file.path(RP, "Analysis/Hypotheses/TRAIL4-20261111T025501Z/out_23190150/susie_credible_sets.tsv"))
tr <- tr[which == "trait1" & grepl("^(T1_|T1r_|T1s_|T2_ALT|S_ALT|C3_)", test)]
tr[, gene := fifelse(grepl("PNPLA3", test), "PNPLA3", fifelse(grepl("PANX1", test), "PANX1", "TNFSF10"))]
tp <- tstrsplit(tr$top_snp, ":")
tr[, `:=`(chromosome = as.integer(tp[[1]]), position = as.integer(tp[[2]]), a = tp[[3]], b = tp[[4]])]
lo <- lift(tr$chromosome, tr$position)
tr[, `:=`(c38 = lo$c38, p38 = lo$p38, r38 = refbase(lo$c38, lo$p38), hits = lo$hits, s38 = lo$s38)]
tr[, ok := hits == 1 & s38 %in% "+" & nchar(a) == 1 & nchar(b) == 1 & (r38 == a | r38 == b)]
tr[is.na(ok), ok := FALSE]
tr[, locus_id := fifelse(gene == "PNPLA3", "named_PNPLA3", fifelse(gene == "PANX1", named_hit[["PANX1"]][1], named_hit[["TNFSF10"]][1]))]
cnt$trail4_top_rows <- nrow(tr); cnt$trail4_top_ref_fail <- sum(!tr$ok)
tr <- tr[ok == TRUE]
mem <- rbind(mem, tr[, .(locus_id, set_id = paste("TRAIL4", test, cs, sep = "|"), k = akey(chromosome, position, a, b), chromosome, position,
                         trait = "ALT", locus = NA_real_, study = paste0("TRAIL4_", test), ancestry = "EUR", method = "trail4_top",
                         pip = top_pip, uid = paste0(c38, ":", p38, ":", r38, ":", ifelse(r38 == a, b, a)), tier = 1L,
                         member_source = "trail4_top")], fill = TRUE)
nset <- setmap[, .(n_sets = .N), by = locus_id]
loci <- merge(loci, nset, by = "locus_id", all.x = TRUE)
loci[is.na(n_sets), n_sets := 0L]
loci[, has_cs := n_sets > 0]
nsub <- mem[!is.na(uid), .(n_substrate_members = uniqueN(uid[member_source == "substrate"]), n_scored_members = uniqueN(uid)), by = locus_id]
nout <- mem[is.na(uid), .(n_members_outside_substrate = uniqueN(k)), by = locus_id]
loci <- merge(merge(loci, nsub, by = "locus_id", all.x = TRUE), nout, by = "locus_id", all.x = TRUE)
loci[is.na(n_substrate_members), n_substrate_members := 0L]
loci[is.na(n_scored_members), n_scored_members := 0L]
loci[is.na(n_members_outside_substrate), n_members_outside_substrate := 0L]
loci[, scored_as := fifelse(n_scored_members > 0, "credible_set_members",
                     fifelse(has_cs, "cs_outside_substrate_lead_proxies", "no_cs_lead_proxies"))]
cnt$loci_by_tier_and_scoring <- loci[, .N, by = .(tier, scored_as)][order(tier, scored_as)]
cnt$members_outside_substrate_tier12_by_ancestry <- mem[is.na(uid) & tier <= 2, .(variants = uniqueN(k)), by = ancestry]
cnt$members_outside_substrate_tier12_variants <- uniqueN(mem[is.na(uid) & tier <= 2]$k)
cnt$members_outside_substrate_tier12_variants_EUR_EAS <- uniqueN(mem[is.na(uid) & tier <= 2 & ancestry %in% c("EUR", "EAS")]$k)
cnt$tier1_disease_effect_loci <- sum(loci$disease_effect)
cnt$tier1_disease_effect_loci_with_cs <- sum(loci$disease_effect & loci$has_cs)
cnt$half_rule_met <- cnt$tier1_disease_effect_loci_with_cs >= cnt$tier1_disease_effect_loci / 2

# ---------------------------------------------------------------- 1KG EUR AF panel (hg19)
af <- rbindlist(lapply(1:22, function(ch) fread(file.path(LDREF, sprintf("panel_af/eur/chr%d.af.tsv.gz", ch)),
                                                 select = c("chromosome", "position", "bim_a1", "bim_a2", "af_a1"))))
af[, k := akey(chromosome, position, bim_a1, bim_a2)]
af[, maf := pmin(af_a1, 1 - af_a1)]
af <- unique(af[, .(k, chromosome, position, bim_a1, bim_a2, maf)], by = "k")
afm <- af[, .(k, maf)]

# ---------------------------------------------------------------- loci without a set: lead + r2 > 0.8 proxies
nocs <- loci[scored_as != "credible_set_members" & source %in% c("DIG1", "named")]
prox <- list()
for (i in seq_len(nrow(nocs))) {
  L <- nocs[i]
  bim <- fread(file.path(LDREF, sprintf("1kg_eur/chr%d_eur.bim", L$chr)), header = FALSE)
  hit <- bim[V4 == L$pos19]
  rows <- data.table(locus_id = L$locus_id, chr = L$chr, pos19 = L$pos19, a = L$a1, b = L$a2, r2 = 1, role = "lead")
  if (nrow(hit) == 1) {
    of <- file.path(TMP, paste0("ld_", L$locus_id))
    system2(PLINK, c("--bfile", file.path(LDREF, sprintf("1kg_eur/chr%d_eur", L$chr)), "--r2", "--ld-snp", hit$V2,
                     "--ld-window-kb", "500", "--ld-window", "999999", "--ld-window-r2", "0.8", "--memory", "4000",
                     "--out", of), stdout = FALSE, stderr = FALSE)
    if (file.exists(paste0(of, ".ld"))) {
      ld <- fread(paste0(of, ".ld"))
      ld <- ld[SNP_B != SNP_A]
      if (nrow(ld)) {
        pa <- tstrsplit(ld$SNP_B, ":")
        rows <- rbind(rows, data.table(locus_id = L$locus_id, chr = L$chr, pos19 = ld$BP_B, a = pa[[3]], b = pa[[4]],
                                       r2 = ld$R2, role = "proxy"))
      }
    }
    rows[, lead_in_panel := TRUE]
  } else rows[, lead_in_panel := FALSE]
  prox[[i]] <- rows
}
prox <- rbindlist(prox, fill = TRUE)
cnt$loci_without_cs <- nrow(nocs)
cnt$loci_without_cs_lead_in_panel <- if (nrow(prox)) uniqueN(prox[lead_in_panel == TRUE]$locus_id) else 0
if (nrow(prox)) {
  lo <- lift(prox$chr, prox$pos19)
  prox[, `:=`(c38 = lo$c38, p38 = lo$p38, s38 = lo$s38, hits = lo$hits)]
  prox[, ref38 := refbase(c38, p38)]
  prox[, ok := hits == 1 & s38 == "+" & nchar(a) == 1 & nchar(b) == 1 & (ref38 == a | ref38 == b)]
  prox[, alt38 := ifelse(ref38 == a, b, a)]
  prox[, uid := ifelse(ok, paste0(c38, ":", p38, ":", ref38, ":", alt38), NA_character_)]
  prox[, k := akey(chr, pos19, a, b)]
  cnt$lead_proxy_rows <- nrow(prox); cnt$lead_proxy_failed_lift_or_ref <- sum(!prox$ok)
}

# ---------------------------------------------------------------- members to score
mv <- mem[!is.na(uid)]
cnt$member_rows_by_tier <- as.list(mv[, .N, by = tier][order(tier)][, setNames(N, tier)])
mvar <- mv[, .(tier = min(tier), loci = paste(sort(unique(locus_id)), collapse = ","),
               traits = paste(sort(unique(trait)), collapse = ","), studies = paste(sort(unique(study)), collapse = ","),
               n_sets = uniqueN(set_id), max_pip = max(pip, na.rm = TRUE), k = k[1],
               member_sources = paste(sort(unique(member_source)), collapse = ",")), by = uid]
cons <- fread(file.path(FM, "credible_set_variant_consequences.csv"))
cons[, k := akey(chr, pos, ref, alt)]
mvar <- merge(mvar, unique(cons[, .(k, Consequence, class)], by = "k"), by = "k", all.x = TRUE)
mvar[, coding_pip05 := class == "coding_protein_altering" & max_pip >= 0.5]
mvar[is.na(coding_pip05), coding_pip05 := FALSE]
mvar <- merge(mvar, afm, by = "k", all.x = TRUE)
# primary locus of a variant for background accounting: its lowest-tier locus, first by name
prim <- mv[order(uid, tier, locus_id)][, .(primary_locus = locus_id[1]), by = uid]
mvar <- merge(mvar, prim, by = "uid")
cnt$member_variants_by_tier <- as.list(mvar[, .N, by = tier][order(tier)][, setNames(N, tier)])

# ---------------------------------------------------------------- background
fmr <- fread(file.path(FM, "combined_finemapping.csv"),
             select = c("chromosome", "position", "allele1", "allele2", "study", "locus", "max_pip"))
say("combined_finemapping rows %d", nrow(fmr))
cnt$combined_finemapping_rows <- nrow(fmr)
wins <- fmr[, .(wchr = chromosome[1], wmin = min(position), wmax = max(position)), by = .(study, locus)]
pmx <- fmr[, .(pmax = max(max_pip, na.rm = TRUE)), by = .(chromosome, position)]
csk <- unique(cs[, .(chromosome, position)])
cand <- unique(fmr[nchar(allele1) == 1 & nchar(allele2) == 1, .(chromosome, position, allele1, allele2)])
rm(fmr); gc()
cand[, npos := .N, by = .(chromosome, position)]
cand <- cand[npos == 1][, npos := NULL]           # one allele pair per position
cand <- merge(cand, pmx[is.finite(pmax) & pmax < 0.001], by = c("chromosome", "position"))
cand <- cand[!csk, on = .(chromosome, position)]
cand[, k := akey(chromosome, position, allele1, allele2)]
cand <- merge(cand, afm, by = "k")
say("background candidates after PIP, CS and AF filters: %d", nrow(cand))
# locus windows: union of the windows of the locus's study x locus pairs; lead +-500 kb if none
sl <- unique(mem[, .(locus_id, study, locus)])
lw <- merge(sl, wins, by = c("study", "locus"))[, .(wchr = wchr[1], wmin = min(wmin), wmax = max(wmax)), by = locus_id]
lw <- rbind(lw, loci[!(locus_id %in% lw$locus_id), .(locus_id, wchr = chr, wmin = pos19 - WIN, wmax = pos19 + WIN)])
# members (and leads/proxies) per locus that need a matched background
need <- rbind(mvar[, .(locus_id = primary_locus, uid, maf)],
              if (nrow(prox)) merge(prox[ok == TRUE, .(locus_id, uid, k)], afm, by = "k", all.x = TRUE)[, .(locus_id, uid, maf)])
need <- unique(need, by = "uid")
cc <- cand[lw, on = .(chromosome = wchr, position >= wmin, position <= wmax), nomatch = 0,
           .(locus_id, chromosome = x.chromosome, position = x.position, allele1, allele2, maf, k)]
lo <- lift(cc$chromosome, cc$position)
cc[, `:=`(c38 = lo$c38, p38 = lo$p38, ok = lo$hits == 1 & lo$s38 == "+")]
cc <- cc[ok == TRUE]
cc[, ref38 := refbase(c38, p38)]
cc <- cc[ref38 == allele1 | ref38 == allele2]
cc[, uid := paste0(c38, ":", p38, ":", ref38, ":", ifelse(ref38 == allele1, allele2, allele1))]
cc <- cc[!(uid %in% need$uid)]
cc[, src := 1L]
# second pool, used only when the first is exhausted: 1KG EUR panel SNVs in the window that are not in any
# fine-mapping row with PIP >= 0.001, not at a credible-set position, and not already in the first pool
c2 <- af[nchar(bim_a1) == 1 & nchar(bim_a2) == 1, .(chromosome, position, allele1 = bim_a1, allele2 = bim_a2, maf, k)]
c2[, npos := .N, by = .(chromosome, position)]
c2 <- c2[npos == 1][, npos := NULL]
c2 <- c2[!csk, on = .(chromosome, position)]
c2 <- c2[!pmx[!is.finite(pmax) | pmax >= 0.001], on = .(chromosome, position)]
c2 <- c2[lw, on = .(chromosome = wchr, position >= wmin, position <= wmax), nomatch = 0,
         .(locus_id, chromosome = x.chromosome, position = x.position, allele1, allele2, maf, k)]
c2 <- c2[!cc, on = .(locus_id, chromosome, position)]
lo <- lift(c2$chromosome, c2$position)
c2[, `:=`(c38 = lo$c38, p38 = lo$p38, ok = lo$hits == 1 & lo$s38 == "+")]
c2 <- c2[ok == TRUE]
c2[, ref38 := refbase(c38, p38)]
c2 <- c2[ref38 == allele1 | ref38 == allele2]
c2[, uid := paste0(c38, ":", p38, ":", ref38, ":", ifelse(ref38 == allele1, allele2, allele1))]
c2 <- c2[!(uid %in% need$uid)]
c2[, src := 2L]
cc <- rbind(cc, c2, fill = TRUE)
setorder(cc, locus_id, chromosome, position)
cc[, rk := runif(.N)]
bg <- list()
for (L in unique(need$locus_id)) {
  m <- need[locus_id == L][order(uid)]
  pool <- cc[locus_id == L]
  used <- character(0)
  for (j in seq_len(nrow(m))) {
    p <- pool[!(uid %in% used)]
    if (!nrow(p)) break
    if (any(p$src == 1L)) p <- p[src == 1L]
    dm <- if (is.na(m$maf[j])) p$rk else abs(p$maf - m$maf[j])
    pick <- p[order(dm, rk)][1]
    used <- c(used, pick$uid)
    bg[[length(bg) + 1]] <- data.table(locus_id = L, member_uid = m$uid[j], uid = pick$uid, member_maf = m$maf[j],
                                       bg_maf = pick$maf, maf_diff = abs(pick$maf - m$maf[j]), bg_source = pick$src)
  }
}
bg <- rbindlist(bg)
bg <- merge(bg, loci[, .(locus_id, tier)], by = "locus_id")
cnt$background_rows <- nrow(bg); cnt$background_unique_variants <- uniqueN(bg$uid)
cnt$background_median_maf_diff <- median(bg$maf_diff, na.rm = TRUE)
cnt$background_members_without_match <- nrow(need) - nrow(bg)
cnt$members_without_af <- sum(is.na(need$maf))
cnt$background_by_source <- as.list(table(bg$bg_source))
cnt$unmatched_by_locus <- merge(need[, .(n_need = .N), by = locus_id], bg[, .(n_bg = .N), by = locus_id], all.x = TRUE)[is.na(n_bg) | n_bg < n_need]

# ---------------------------------------------------------------- controls (hg19 -> hg38, REF-checked)
# SORT1 position = DIG1 locus 6 lead; PNPLA3/TM6SF2/GCKR = DIG1 script 73 constants; HSD17B13 rs72613567 hg19
# 4:88231392 (T>TA insertion) is recalled, not read from a repo file, and is checked by its hg38 REF base here.
ctl <- data.table(control = c("SORT1_rs12740374", "PNPLA3_rs738409", "TM6SF2_rs58542926", "GCKR_rs1260326",
                              "HSD17B13_rs72613567"),
                  chr = c(1L, 22L, 19L, 2L, 4L), pos19 = c(109817590L, 44324727L, 19379549L, 27730940L, 88231392L))
lo <- lift(ctl$chr, ctl$pos19)
ctl[, `:=`(c38 = lo$c38, p38 = lo$p38, ref38 = refbase(lo$c38, lo$p38))]
ctl <- merge(ctl, afm[, c("kc", "kp") := .(as.integer(sub(":.*", "", k)), as.integer(tstrsplit(k, ":")[[2]]))],
             by.x = c("chr", "pos19"), by.y = c("kc", "kp"), all.x = TRUE)
ctl[, in_member_set := sapply(paste0(c38, ":", p38, ":"), function(p) any(startsWith(mvar$uid, p)))]

# ---------------------------------------------------------------- variant list and budget
v <- rbind(mvar[, .(uid, role = "member", tier, loci, traits, coding_pip05, max_pip)],
           if (nrow(prox)) prox[ok == TRUE, .(uid, role = ifelse(role == "lead", "lead_no_cs", "proxy_no_cs"), tier = loci[match(prox[ok == TRUE]$locus_id, loci$locus_id)]$tier,
                                         loci = locus_id, traits = NA_character_, coding_pip05 = FALSE, max_pip = NA_real_)],
           bg[, .(uid, role = "background", tier, loci = locus_id, traits = NA_character_, coding_pip05 = FALSE, max_pip = NA_real_)],
           fill = TRUE)
v <- v[order(uid, tier, role)]
vu <- v[, .(tier = min(tier), roles = paste(sort(unique(role)), collapse = ","), loci = paste(sort(unique(unlist(strsplit(loci, ",")))), collapse = ",")), by = uid]
cnt$variants_to_score_by_tier <- as.list(vu[, .N, by = tier][order(tier)][, setNames(N, tier)])
fwrite(loci, file.path(O, "loci.tsv"), sep = "\t")
fwrite(mem, file.path(O, "credible_set_member_rows.tsv.gz"), sep = "\t")
fwrite(mvar, file.path(O, "member_variants.tsv"), sep = "\t")
fwrite(prox, file.path(O, "lead_proxy_no_cs.tsv"), sep = "\t")
fwrite(bg, file.path(O, "background.tsv"), sep = "\t")
fwrite(v, file.path(O, "variant_roles.tsv"), sep = "\t")
fwrite(vu, file.path(O, "variants_to_score.tsv"), sep = "\t")
fwrite(ctl, file.path(O, "controls.tsv"), sep = "\t")
write(jsonlite::toJSON(cnt, auto_unbox = TRUE, pretty = TRUE, na = "null"), file.path(O, "b0_counts.json"))
capture.output(sessionInfo(), file = file.path(O, "sessionInfo_b0.txt"))
print(cnt)
print(loci[tier == 1, .(locus_id, genes, disease_effect, n_sets, has_cs)])
print(ctl)
