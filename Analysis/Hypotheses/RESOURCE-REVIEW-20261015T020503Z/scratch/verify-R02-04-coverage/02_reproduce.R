suppressPackageStartupMessages(library(data.table))
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
S <- paste0(B, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R02-04-coverage/")
coloc <- fread(paste0(B, "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"),
               select = c("gwas_name","gene","PP.H4.susie","PP.H4.abf","top_snp","top_snp_PP"))
MAIN <- fread(paste0(B, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[placement == "main", study_name]
coloc <- coloc[gwas_name %in% MAIN]
setnames(coloc, c("PP.H4.susie","PP.H4.abf"), c("pp4_susie","pp4_abf"))
coloc <- coloc[!is.na(gene) & gene != ""]
coloc[, pp4_susie := suppressWarnings(as.numeric(pp4_susie))]
coloc[, pp4_abf := suppressWarnings(as.numeric(pp4_abf))]
coloc <- coloc[!is.na(top_snp) & top_snp != ""]
coloc[, c("snp_chr","snp_pos") := tstrsplit(top_snp, ":", fixed = TRUE, type.convert = TRUE)]

cn <- fread(paste0(S, "fm_canon_filtered.csv"),
            select = c("study","chromosome","position","ancestry","susie_pip","susie_converged","carma_pip","recommended_pip","susiex_pip","mesusie_pip"))
cn[, x_source := "canonical_finemapping"]
tg <- fread(paste0(S, "fm_targ_filtered.csv"),
            select = c("study","chromosome","position","ancestry","susie_pip","susie_converged","carma_pip","recommended_pip"))
tg[, `:=`(susiex_pip = NA_real_, mesusie_pip = NA_real_, x_source = "coloc_targeted_finemapping")]
fm0 <- rbindlist(list(cn, tg), use.names = TRUE)
num <- function(x) suppressWarnings(as.numeric(x))
for (v in c("susie_pip","carma_pip","recommended_pip","susiex_pip","mesusie_pip")) set(fm0, j = v, value = num(fm0[[v]]))
fm0[, source_rank := fifelse(x_source == "canonical_finemapping", 1L, 2L)]

max_or_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
label_at_max <- function(value, label) if (all(is.na(value))) NA_character_ else label[which.max(fifelse(is.na(value), -Inf, value))]
gene_y <- coloc[, .(pp4_susie = max_or_na(pp4_susie),
                    susie_gwas = label_at_max(pp4_susie, gwas_name),
                    y_top_snp = label_at_max(pp4_susie, top_snp)), by = gene]

build <- function(xcol_expr, tag) {
  fm <- copy(fm0)
  fm[, causal_pip := eval(xcol_expr)]
  fm[!is.finite(causal_pip), causal_pip := NA_real_]
  fm <- fm[is.finite(causal_pip)]
  fm <- fm[order(study, chromosome, position, -causal_pip, source_rank)]
  fu <- unique(fm, by = c("study","chromosome","position"))
  cm <- merge(coloc, fu[, .(study, chromosome, position, causal_pip, x_source, ancestry, susie_pip, susie_converged, carma_pip, recommended_pip, susiex_pip, mesusie_pip)],
              by.x = c("gwas_name","snp_chr","snp_pos"), by.y = c("study","chromosome","position"))
  gx <- cm[order(gene, -causal_pip, x_source), .SD[1], by = gene][,
          .(gene, causal_pip, x_source, x_top_snp = top_snp, x_gwas = gwas_name, ancestry, susie_pip, susie_converged, carma_pip, recommended_pip, susiex_pip, mesusie_pip)]
  dt <- merge(gene_y, gx, by = "gene", all.x = TRUE)
  pl <- dt[!is.na(pp4_susie) & is.finite(causal_pip)]
  pl[, pp4 := pmax(pp4_susie, 0)]
  roster <- c("HSD17B13","TM6SF2","MBOAT7","GCKR","RORA")
  pl[, plotted := (pp4 > 0.05) | (causal_pip >= 0.99) | gene %in% roster]
  pl[, variant := tag]
  pl
}
base <- build(quote(pmax(recommended_pip, susie_pip, na.rm = TRUE)), "baseline")
cat("baseline rows", nrow(base), "plotted", sum(base$plotted), "\n")
# compare with installed source CSV
src <- fread(paste0(B, "figures/main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc_source.csv"))
m <- merge(src, base[, .(gene, cp = causal_pip, xs = x_top_snp, pg = susie_gwas, pp = pp4, pl = plotted)], by = "gene", all = TRUE)
cat("genes only in src:", sum(is.na(m$cp)), " only in rederived:", sum(is.na(m$causal_pip)), "\n")
cat("max |x diff|:", max(abs(m$causal_pip - round(m$cp, 4)), na.rm = TRUE),
    " max |y diff|:", max(abs(m$pp4 - round(m$pp, 4)), na.rm = TRUE),
    " x_top_snp mismatches:", sum(m$x_top_snp != m$xs, na.rm = TRUE),
    " pp4_gwas mismatches:", sum(m$pp4_gwas != m$pg, na.rm = TRUE),
    " plotted mismatches:", sum(m$plotted != m$pl, na.rm = TRUE), "\n")
saveRDS(list(base = base, fm0 = fm0, coloc = coloc, gene_y = gene_y), paste0(S, "state.rds"))
