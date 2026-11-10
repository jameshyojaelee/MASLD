#!/usr/bin/env Rscript
# U4 (a): independent recount of the gene-level COLOC calls in release r2 from the
# per-pair master table. Written from Methods/legend text, not from the producer.
#   Tier 1 = direct NAFLD/NASH/PDFF; tier 2 = ALT/AST/GGT (Methods "35 primary strata").
#   SuSiE positive gene = any tier-1/2 row with PP.H4.susie > 0.5.
#   Untestable row = method "susie_untestable_insufficient_shared_posterior".
#   Gene state priority: positive > untestable > tested (finite PP.H4.susie <= 0.5) > none.
#   ABF positive gene = any tier-1/2 row PP.H4.abf > 0.5. Union = SuSiE pos OR ABF pos.
#   Driver = tier-1/2 study with the highest SuSiE PP.H4; ties -> alphabetically first study id.
suppressPackageStartupMessages(library(data.table))
D <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
OUT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/recount"
UNT <- "susie_untestable_insufficient_shared_posterior"
res <- list()
put <- function(k, v) { res[[k]] <<- v; cat(sprintf("%-60s %s\n", k, paste(v, collapse = " "))) }

m <- fread(file.path(D, "susie_coloc_all_gwas.csv"),
           select = c("gwas_name", "gene", "ensembl", "chr", "PP.H4.abf", "PP.H3.susie", "PP.H4.susie",
                      "method", "ancestry", "n_cs_pairs", "n_cs_pairs_eligible", "eligibility_source"),
           colClasses = list(character = c("gene", "ensembl", "gwas_name")))
put("master_rows", nrow(m))
put("master_unique_gene_study_rows(ensembl+symbol+chr)", uniqueN(m, by = c("gwas_name", "ensembl", "gene", "chr")))
put("master_unique_gwas_ensembl_rows", uniqueN(m[ensembl != ""], by = c("gwas_name", "ensembl")))
put("n_studies", uniqueN(m$gwas_name))
put("rows_blank_ensembl", m[ensembl == "", .N])
put("ensembl_with_version_suffix", m[grepl("\\.", ensembl), .N])
print(m[, .N, by = method])
print(m[, .N, by = eligibility_source])

studies <- sort(unique(m$gwas_name))
trait_of <- function(s) {
  fifelse(grepl("NAFLD|NASH|PDFF", s) & !grepl("ChronLiver|Cirrhosis", s), "direct",
  fifelse(grepl("(^|_)(ALT|AST|GGT)(_|$)", s), "enzyme", "other"))
}
st <- data.table(gwas_name = studies, cls = trait_of(studies))
print(st[, .N, by = cls])
put("tier12_studies", st[cls != "other", .N])
put("tier1_direct", st[cls == "direct", .N]); put("tier2_enzyme", st[cls == "enzyme", .N])
m <- merge(m, st, by = "gwas_name")

# rerun corrections (eligibility_source)
put("rerun_rows_status_all", paste(capture.output(print(m[eligibility_source == "rerun_shared_posterior_check",
     .(pos = sum(PP.H4.susie > 0.5, na.rm = TRUE), le = sum(PP.H4.susie <= 0.5, na.rm = TRUE),
       unt = sum(method == UNT)), ])), collapse = " | "))
put("rerun_rows_status_t12", paste(capture.output(print(m[eligibility_source == "rerun_shared_posterior_check" & cls != "other",
     .(n = .N, pos = sum(PP.H4.susie > 0.5, na.rm = TRUE), le = sum(PP.H4.susie <= 0.5, na.rm = TRUE),
       unt = sum(method == UNT))])), collapse = " | "))

# exact-0.5 values (> vs >= rule)
put("rows_PP.H4.susie_exactly_0.5", m[PP.H4.susie == 0.5, .N])
put("rows_PP.H4.abf_exactly_0.5", m[PP.H4.abf == 0.5, .N])

gene_state <- function(x) {
  # x: rows for a set of studies; gene-level states on Ensembl id
  x <- x[ensembl != ""]
  x[, .(susie_max = suppressWarnings(max(PP.H4.susie, na.rm = TRUE)),
        any_finite = any(is.finite(PP.H4.susie)),
        any_unt = any(method == UNT),
        abf_max = max(PP.H4.abf, na.rm = TRUE)), by = ensembl][
    , state := fifelse(any_finite & susie_max > 0.5, "positive",
               fifelse(any_unt, "untestable",
               fifelse(any_finite, "tested_le_0.5", "no_susie_result")))][]
}

t12 <- m[cls != "other"]
g <- gene_state(t12)
print(g[, .N, by = state])
for (s in c("positive", "tested_le_0.5", "untestable", "no_susie_result")) put(paste0("t12_state_", s), g[state == s, .N])
put("t12_n_ensembl", nrow(g))
put("t12_susie_tested_universe(pos+le)", g[state %in% c("positive", "tested_le_0.5"), .N])
g[, abf_pos := abf_max > 0.5]
put("t12_abf_pos", g[abf_pos == TRUE, .N])
put("t12_union", g[abf_pos | state == "positive", .N])
put("t12_abf_only(union minus susie_pos)", g[abf_pos & state != "positive", .N])
print(g[abf_pos & state != "positive", .N, by = state])
put("untestable_genes", paste(sort(unique(t12[ensembl %in% g[state == "untestable", ensembl], gene])), collapse = ","))
put("untestable_abf_pos", g[state == "untestable" & abf_pos, .N])
put("untestable_best_finite_susie_max", round(g[state == "untestable", max(susie_max[is.finite(susie_max)])], 4))

# p12 sensitivity: PP4' = k PP4 / (k PP4 + 1 - PP4); k = 0.2 (p12 1e-6), k = 2 (p12 1e-5)
p12x <- function(p, k) k * p / (k * p + 1 - p)
put("p12_1e-6_positive_falling_to_le_0.5", g[state == "positive" & p12x(susie_max, 0.2) <= 0.5, .N])
put("p12_1e-5_tested_le_rising_above_0.5", g[state == "tested_le_0.5" & p12x(susie_max, 2) > 0.5, .N])
p12up <- g[state == "tested_le_0.5" & p12x(susie_max, 2) > 0.5, .(ensembl, susie_max)][order(susie_max)]
p12up <- merge(p12up, unique(t12[ensembl != "", .(ensembl, gene)]), by = "ensembl")
fwrite(p12up, file.path(OUT, "p12_1e-5_rising_genes.tsv"), sep = "\t")
put("p12_1e-6_positive_max_falling_value",round(g[state == "positive" & p12x(susie_max, 0.2) <= 0.5, max(susie_max)], 4))

# >= instead of > rule
put("t12_susie_ge_0.5", g[is.finite(susie_max) & susie_max >= 0.5, .N])
put("t12_abf_ge_0.5", g[abf_max >= 0.5, .N])

# symbol-level counts (Fig 2A counts symbols)
sym <- t12[gene != "", .(susie_max = suppressWarnings(max(PP.H4.susie, na.rm = TRUE)), abf_max = max(PP.H4.abf, na.rm = TRUE)), by = gene]
put("sym_susie_pos", sym[is.finite(susie_max) & susie_max > 0.5, .N])
put("sym_union", sym[(is.finite(susie_max) & susie_max > 0.5) | abf_max > 0.5, .N])
put("sym_abf_only", sym[!(is.finite(susie_max) & susie_max > 0.5) & abf_max > 0.5, .N])
# symbols mapping to >1 ensembl among union genes
su <- unique(t12[ensembl %in% g[abf_pos | state == "positive", ensembl], .(gene, ensembl)])
put("union_symbols_with_multiple_ensembl", paste(su[, .N, by = gene][N > 1, gene], collapse = ","))
put("union_ensembl_with_multiple_symbols", paste(su[, .N, by = ensembl][N > 1, ensembl], collapse = ","))

# driver split
pos <- t12[ensembl %in% g[state == "positive", ensembl] & is.finite(PP.H4.susie)]
setorder(pos, ensembl, -PP.H4.susie, gwas_name)
drv <- pos[, .SD[1], by = ensembl]
put("driver_enzyme", drv[cls == "enzyme", .N]); put("driver_direct", drv[cls == "direct", .N])
# ties across classes
tie <- pos[, .(mx = max(PP.H4.susie)), by = ensembl][pos, on = "ensembl"][PP.H4.susie == mx, .(ncls = uniqueN(cls), genes = gene[1]), by = ensembl][ncls > 1]
put("tie_across_classes_genes", paste(tie$genes, collapse = ","))
# direct-preferring rule
setorder(pos, ensembl, -PP.H4.susie, cls, gwas_name)   # "direct" < "enzyme"
drv2 <- pos[, .SD[1], by = ensembl]
put("driver_direct_preferred_enzyme_direct", paste(drv2[cls == "enzyme", .N], drv2[cls == "direct", .N]))

# Fig 1C: support by class
cs <- t12[ensembl != "", .(pos = any(PP.H4.susie > 0.5, na.rm = TRUE), fin = any(is.finite(PP.H4.susie)),
                           unt = any(method == UNT)), by = .(ensembl, cls)]
w <- dcast(cs, ensembl ~ cls, value.var = c("pos", "fin", "unt"), fill = FALSE)
put("fig1c_enzyme_only", w[pos_enzyme & !pos_direct, .N])
put("fig1c_direct_only", w[!pos_enzyme & pos_direct, .N])
put("fig1c_both", w[pos_enzyme & pos_direct, .N])
put("fig1c_enzyme_supported", w[pos_enzyme == TRUE, .N]); put("fig1c_direct_supported", w[pos_direct == TRUE, .N])
put("fig1c_enzyme_sup_direct_tested_le", w[pos_enzyme & !pos_direct & fin_direct, .N])
put("fig1c_enzyme_sup_direct_untested", w[pos_enzyme & !pos_direct & !fin_direct, .N])
put("fig1c_direct_sup_enzyme_tested_le", w[pos_direct & !pos_enzyme & fin_enzyme, .N])
put("fig1c_direct_sup_enzyme_untested", w[pos_direct & !pos_enzyme & !fin_enzyme, .N])
put("fig1c_tested_any_enzyme", w[fin_enzyme == TRUE, .N]); put("fig1c_tested_any_direct", w[fin_direct == TRUE, .N])
put("fig1c_tested_both", w[fin_enzyme & fin_direct, .N])
put("fig1c_tested_both_enzyme_sup", w[fin_enzyme & fin_direct & pos_enzyme, .N])
put("fig1c_tested_both_direct_sup", w[fin_enzyme & fin_direct & pos_direct, .N])

# Fig 2E: ancestry of supporting GWAS for the positive genes
anc <- t12[ensembl %in% g[state == "positive", ensembl]]
pa <- anc[, .(eur_pos = any(PP.H4.susie > 0.5 & ancestry == "EUR", na.rm = TRUE),
              noneur_pos = any(PP.H4.susie > 0.5 & ancestry != "EUR", na.rm = TRUE),
              noneur_fin = any(is.finite(PP.H4.susie) & ancestry != "EUR"),
              noneur_abf = any(PP.H4.abf > 0.5 & ancestry != "EUR")), by = ensembl]
put("fig2e_eur_only", pa[eur_pos & !noneur_pos, .N])
put("fig2e_eur_only_noneur_tested", pa[eur_pos & !noneur_pos & noneur_fin, .N])
put("fig2e_eur_only_noneur_untested", pa[eur_pos & !noneur_pos & !noneur_fin, .N])
put("fig2e_eur_only_noneur_abf_pos", pa[eur_pos & !noneur_pos & noneur_abf, .N])
put("fig2e_eur_only_noneur_abf_pos_tested/untested", paste(pa[eur_pos & !noneur_pos & noneur_abf & noneur_fin, .N],
                                                          pa[eur_pos & !noneur_pos & noneur_abf & !noneur_fin, .N]))
put("fig2e_eur_and_noneur", pa[eur_pos & noneur_pos, .N])
put("fig2e_noneur_only", pa[!eur_pos & noneur_pos, .N])
# Fig 2A split of single-signal-only genes by SuSiE state, symbol level (state priority as above)
symst <- t12[gene != "", .(susie_max = suppressWarnings(max(PP.H4.susie, na.rm = TRUE)),
                           any_finite = any(is.finite(PP.H4.susie)), any_unt = any(method == UNT),
                           abf_max = max(PP.H4.abf, na.rm = TRUE)), by = gene][
  , state := fifelse(any_finite & susie_max > 0.5, "positive", fifelse(any_unt, "untestable",
             fifelse(any_finite, "tested_le_0.5", "no_susie_result")))]
put("fig2a_sym_abf_only_by_state", paste(capture.output(print(symst[abf_max > 0.5 & state != "positive", .N, by = state])), collapse = " | "))
put("fig2a_sym_untestable_abf_pos", symst[state == "untestable" & abf_max > 0.5, .N])

# full 50-study portfolio
gf <- gene_state(m)
for (s in c("positive", "tested_le_0.5", "untestable", "no_susie_result")) put(paste0("full_state_", s), gf[state == s, .N])
put("full_abf_pos", gf[abf_max > 0.5, .N]); put("full_union", gf[abf_max > 0.5 | state == "positive", .N])

# producer table comparison (Ensembl level)
p <- fread(file.path(D, "gene_level_coloc_tier12.csv"))
cmp <- merge(g[, .(ensembl, state, susie_max, abf_max)], p[, .(ensembl, susie_state_t12, coloc_best_susie_pp4, coloc_best_abf_pp4, driving_gwas, susie_driving_gwas)],
             by = "ensembl", all = TRUE)
put("cmp_rows_only_mine", cmp[is.na(susie_state_t12), .N]); put("cmp_rows_only_producer", cmp[is.na(state), .N])
put("cmp_state_disagree", cmp[state != susie_state_t12, .N])
put("cmp_abf_maxdiff", max(abs(cmp$abf_max - cmp$coloc_best_abf_pp4), na.rm = TRUE))
put("cmp_susie_maxdiff", max(abs(cmp[is.finite(susie_max)]$susie_max - cmp[is.finite(susie_max)]$coloc_best_susie_pp4), na.rm = TRUE))
dd <- merge(drv[, .(ensembl, my_driver = gwas_name)], p[, .(ensembl, susie_driving_gwas)], by = "ensembl")
put("cmp_driver_disagree", dd[my_driver != susie_driving_gwas, .N])
fwrite(g, file.path(OUT, "gene_state_t12_independent.tsv"), sep = "\t")
fwrite(drv[, .(ensembl, gene, driver = gwas_name, cls, PP.H4.susie)], file.path(OUT, "driver_independent.tsv"), sep = "\t")
fwrite(data.table(metric = names(res), value = sapply(res, function(v) paste(v, collapse = " "))), file.path(OUT, "recount_results.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
