#!/usr/bin/env Rscript
# Blind re-check (round 4, workstream F2abce): re-derive Fig 2B, 2C ownership,
# 2E locus grouping and the R16-08 COLOC r2 row count from release files.
suppressPackageStartupMessages({
  library(data.table); library(VariantAnnotation)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(GenomicRanges)
  library(org.Hs.eg.db); library(AnnotationDbi)
})
set.seed(20260923)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT <- file.path(BASE, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F2abce/out")
dir.create(OUT, showWarnings = FALSE)
FM <- file.path(BASE, "GWAS/finemapping/results")
say <- function(...) cat(sprintf(...), "\n", sep = "")

tier <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))
MAIN <- tier[placement == "main", study_name]
say("main strata %d; MVP EUR main: %s", length(MAIN), paste(grep("^MVP_.*_EUR$", MAIN, value = TRUE), collapse = ","))

cr <- fread(file.path(FM, "susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"),
            select = c("gwas_name", "gene", "chr", "PP.H4.susie", "PP.H4.abf", "top_snp", "ancestry"))
# ---- R16-08: rows for UKBB EUR NAFLD chr1 in COLOC r2 (all rows, before any filter)
u <- cr[gwas_name == "2023_36280732_NAFLD_UKBB_EUR" & as.character(chr) == "1"]
s4 <- suppressWarnings(as.numeric(u$PP.H4.susie))
say("R16-08 COLOC r2 rows study chr1: %d; finite PP.H4.susie %d; >0.5 %d; max %.4f",
    nrow(u), sum(is.finite(s4)), sum(s4 > 0.5, na.rm = TRUE), max(s4, na.rm = TRUE))
say("R16-08 all chromosomes for study: %d", cr[gwas_name == "2023_36280732_NAFLD_UKBB_EUR", .N])

co <- cr[gwas_name %in% MAIN & !is.na(gene) & gene != "" & !is.na(top_snp) & top_snp != ""]
co[, s4 := suppressWarnings(as.numeric(PP.H4.susie))]
co[, a4 := suppressWarnings(as.numeric(PP.H4.abf))]
co[, row_id := .I]

# ---------------- Fig 2B ----------------
fin <- co[is.finite(s4)]
say("genes with finite SuSiE PP.H4 in main strata: %d", uniqueN(fin$gene))
mx <- fin[, .(m = max(s4)), by = gene]
fin <- merge(fin, mx, by = "gene")
atmax <- fin[s4 == m]
ties <- atmax[, .(n_keys = uniqueN(paste(gwas_name, top_snp))), by = gene][n_keys > 1]
say("genes whose max PP.H4 is tied across distinct (GWAS, lead): %d", nrow(ties))
y <- atmax[order(gene, row_id)][, .SD[1], by = gene]
y[, c("cchr", "cpos") := tstrsplit(top_snp, ":", fixed = TRUE)]
y[, cpos := as.integer(cpos)]
say("genes whose max comes from MVP EUR main strata: %d", y[grepl("^MVP_.*_EUR$", gwas_name), .N])

keys <- unique(y[, .(study = gwas_name, chromosome = cchr, position = cpos)])
rd <- function(f, src) {
  x <- fread(f, select = c("study", "chromosome", "position", "susie_pip"),
             colClasses = list(character = c("chromosome", "susie_pip")))
  x <- x[keys, on = .(study, chromosome, position), nomatch = 0L]
  x[, pip := suppressWarnings(as.numeric(susie_pip))]
  x[, src := src]
  x
}
fm <- rbind(rd(file.path(FM, "combined_finemapping.csv"), "canonical"),
            rd(file.path(FM, "combined_finemapping_coloc_targeted.csv"), "targeted"))
fwrite(fm, file.path(OUT, "fig2b_lead_rows_recheck.tsv"), sep = "\t")
mvp_eur <- grep("^MVP_.*_EUR$", MAIN, value = TRUE)
st_all <- c(unique(fread(file.path(FM, "combined_finemapping.csv"), select = "study")$study),
            unique(fread(file.path(FM, "combined_finemapping_coloc_targeted.csv"), select = "study")$study))
say("MVP EUR main strata present in either fine-mapping file: %d (%s)", sum(mvp_eur %in% st_all),
    paste(intersect(mvp_eur, st_all), collapse = ","))
agg <- fm[, .(n_rec = .N, n_fin = sum(is.finite(pip)),
              xmax = if (any(is.finite(pip))) max(pip, na.rm = TRUE) else NA_real_,
              xmin = if (any(is.finite(pip))) min(pip, na.rm = TRUE) else NA_real_),
          by = .(study, chromosome, position)]
b <- merge(y[, .(gene, gwas_name, top_snp, cchr, cpos, pp4 = s4)], agg,
           by.x = c("gwas_name", "cchr", "cpos"), by.y = c("study", "chromosome", "position"), all.x = TRUE)
say("genes with any fine-mapping record at y's lead: %d; with finite SuSiE PIP: %d",
    b[!is.na(n_rec), .N], b[is.finite(xmax), .N])
rule <- function(x) {
  pl <- is.finite(x) & (b$pp4 > 0.05 | x >= 0.99)
  p4 <- pmax(b$pp4, 0)
  c(shown = sum(pl), corner = sum(pl & x >= 0.9 & p4 < 0.5), topright = sum(pl & x >= 0.9 & p4 > 0.5))
}
print(rbind(max = rule(b$xmax), min = rule(b$xmin)))
b[, shown := is.finite(xmax) & (pp4 > 0.05 | xmax >= 0.99)]
sh <- b[shown == TRUE]
say("shown with >1 record: %d; xmax != xmin exact: %d; |diff|>1e-6: %d; differ at 4 dp: %d; differ at 3 dp: %d",
    sh[n_rec > 1, .N], sh[xmax != xmin, .N], sh[abs(xmax - xmin) > 1e-6, .N],
    sh[round(xmax, 4) != round(xmin, 4), .N], sh[round(xmax, 3) != round(xmin, 3), .N])
say("shown with >1 FINITE record: %d", sh[n_fin > 1, .N])
print(b[gene %in% c("RORA", "POLL", "HSD17B13", "TM6SF2", "MBOAT7", "GCKR", "PNPLA3"),
        .(gene, gwas_name, top_snp, pp4 = round(pp4, 4), n_rec, xmax = round(xmax, 4), xmin = round(xmin, 4), shown)])
min_shown <- b[is.finite(xmin) & (pp4 > 0.05 | xmin >= 0.99), gene]
say("dropped under min rule: %s", paste(setdiff(sh$gene, min_shown), collapse = ","))

# Class of the plotted lead (hg19 TxDb knownGene), own-gene coding ownership
annot <- function(vk) {
  vk <- unique(vk)
  p <- tstrsplit(vk, ":", fixed = TRUE)
  gr <- GRanges(paste0("chr", p[[1]]), IRanges(as.integer(p[[2]]), width = 1), vk = vk)
  loc <- locateVariants(gr, TxDb.Hsapiens.UCSC.hg19.knownGene,
                        AllVariants(promoter = PromoterVariants(upstream = 2000, downstream = 200)))
  data.table(vk = mcols(gr)$vk[loc$QUERYID], loc = as.character(loc$LOCATION), eg = as.character(loc$GENEID))
}
rank <- c(coding = 1, spliceSite = 2, fiveUTR = 3, threeUTR = 4, promoter = 5, intron = 6, intergenic = 7)
classify <- function(genes, vks) {
  L <- annot(vks)
  L[!loc %in% names(rank), loc := "intergenic"]
  top <- L[, .(cls = loc[which.min(rank[loc])]), by = vk]
  codeg <- L[loc %in% c("coding", "spliceSite") & !is.na(eg), .(eg = list(unique(eg))), by = vk]
  d <- data.table(gene = genes, vk = vks)
  d <- merge(d, top, by = "vk", all.x = TRUE); d[is.na(cls), cls := "intergenic"]
  d <- merge(d, codeg, by = "vk", all.x = TRUE)
  sym2eg <- unique(rbind(
    as.data.table(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, unique(d$gene), "ENTREZID", "SYMBOL"))),
    setnames(as.data.table(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, unique(d$gene), "ENTREZID", "ALIAS"))), "ALIAS", "SYMBOL")))
  alleg <- unique(unlist(d$eg)); alleg <- alleg[!is.na(alleg)]
  eg2sym <- as.data.table(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, alleg, "SYMBOL", "ENTREZID")))
  d[, own := mapply(function(g, e) length(e) > 0 && any(e %in% sym2eg[SYMBOL == g, ENTREZID]), gene, eg)]
  d[, coding_syms := vapply(eg, function(e) if (is.null(e)) "" else paste(sort(eg2sym[ENTREZID %in% e, SYMBOL]), collapse = ";"), "")]
  d[, cls3 := fifelse(cls %in% c("coding", "spliceSite"), fifelse(own, "coding_this", "coding_other"), "noncoding")]
  d[, eg := NULL]
  d
}
cb <- classify(sh$gene, sh$top_snp)
print(cb[, .N, by = cls3])
print(cb[cls3 != "noncoding"][order(cls3, gene)])
print(cb[gene %in% c("HSD17B13", "RORA", "TM6SF2")])
fwrite(merge(sh, cb[, .(gene, cls, cls3, own, coding_syms)], by = "gene"), file.path(OUT, "fig2b_recheck_shown.tsv"), sep = "\t")

# ---------------- Fig 2C ownership ----------------
co[, best := pmax(s4, a4, na.rm = TRUE)]
co[!is.finite(best), best := NA_real_]
setorder(co, gene, -best, gwas_name)
cbest <- co[, .SD[1], by = gene][is.finite(best) & best > 0.5]
say("Fig2C union: %d genes", nrow(cbest))
cc <- classify(cbest$gene, cbest$top_snp)
cc[, c5 := fifelse(cls %in% c("coding", "spliceSite"), "coding",
           fifelse(cls %in% c("fiveUTR", "threeUTR"), "UTR", cls))]
print(cc[, .N, by = c5])
cod <- cc[c5 == "coding"]
say("coding genes %d; own %d; other %d; distinct leads %d", nrow(cod), cod[own == TRUE, .N], cod[own == FALSE, .N], uniqueN(cod$vk))
print(cod[order(own, gene), .(gene, vk, cls, own, coding_syms)])
fwrite(cod, file.path(OUT, "fig2c_recheck_coding.tsv"), sep = "\t")

# ---------------- Fig 2E locus grouping ----------------
pos <- co[is.finite(s4) & s4 > 0.5]
pos[, c("vchr", "vpos") := tstrsplit(top_snp, ":", fixed = TRUE)]
pos[, vpos := as.integer(vpos)]
say("multi-signal positive genes: %d; ancestry values: %s", uniqueN(pos$gene), paste(unique(pos$ancestry), collapse = ","))
eur_genes <- unique(pos[ancestry == "EUR", gene])
say("EUR multi-signal positive genes: %d", length(eur_genes))
noneur_only <- setdiff(unique(pos$gene), eur_genes)
say("non-EUR-only genes: %d", length(noneur_only))
eur_leads <- unique(pos[ancestry == "EUR", .(egene = gene, vchr, vpos)])
eur_maxlead <- pos[ancestry == "EUR"][order(-s4)][, .SD[1], by = gene][, .(egene = gene, vchr, vpos)]
mem <- fread(file.path(BASE, "figures/main/fig2_genetics/panels/Fig2E_ancestry_unique_coloc_GWS_membership.tsv"))
mem <- mem[gene %in% noneur_only]
say("membership non-EUR-only rows: %d; GWS pass: %d", nrow(mem), mem[non_eur_gws_gate_pass == TRUE, .N])
ml <- pos[gene %in% noneur_only][order(-s4)][, .SD[1], by = gene][, .(gene, ml_chr = vchr, ml_pos = vpos)]
mem <- merge(mem, ml, by = "gene")
mem[, c("mp_chr", "mp_pos") := tstrsplit(best_snp, ":", fixed = TRUE)]
mem[, mp_pos := as.integer(mp_pos)]
nb <- function(ch, p, L) L[vchr == ch & abs(vpos - p) <= 1e6, paste(sort(unique(egene)), collapse = ",")]
mem[, nb_maxpp4_anyEURlead := mapply(nb, ml_chr, ml_pos, MoreArgs = list(L = eur_leads))]
mem[, nb_minp_anyEURlead := mapply(nb, mp_chr, mp_pos, MoreArgs = list(L = eur_leads))]
mem[, nb_maxpp4_EURmaxlead := mapply(nb, ml_chr, ml_pos, MoreArgs = list(L = eur_maxlead))]
mem[, nb_minp_EURmaxlead := mapply(nb, mp_chr, mp_pos, MoreArgs = list(L = eur_maxlead))]
for (cl in grep("^nb_", names(mem), value = TRUE)) say("%s: %d of %d have a neighbour", cl, sum(mem[[cl]] != ""), nrow(mem))
g <- mem[non_eur_gws_gate_pass == TRUE][order(ml_chr, ml_pos)]
# single-linkage 1 Mb loci among the 9 GWS genes on max-PP.H4 lead
g[, locus := cumsum(c(TRUE, ml_chr[-1] != ml_chr[-.N] | diff(ml_pos) > 1e6))]
print(g[, .(gene, best_gwas, ml = paste(ml_chr, ml_pos, sep = ":"), best_snp, locus, nb_maxpp4_anyEURlead, nb_minp_anyEURlead)])
say("GWS genes %d in %d loci; genes without EUR neighbour %d in %d loci", nrow(g), uniqueN(g$locus),
    g[nb_maxpp4_anyEURlead == "", .N], g[nb_maxpp4_anyEURlead == "", uniqueN(locus)])
say("ancestry of GWS genes: %s", paste(names(table(g$supported_ancestries)), table(g$supported_ancestries), collapse = "; "))
fwrite(mem, file.path(OUT, "fig2e_recheck_neighbours.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
say("done")
