#!/usr/bin/env Rscript
# Prep data tables for the Fig 3 intro/summary panel sketches.
# Outputs CSVs in scripts/figures/sketches_fig3_intro/data/
suppressPackageStartupMessages({ library(data.table) })

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT  <- file.path(BASE, "scripts/figures/sketches_fig3_intro/data")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. GWAS portfolio --------------------------------------------
reg <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_registry.tsv"))
reg[, trait := fcase(
  grepl("NAFLD|NASH",       study_name, ignore.case = TRUE), "NAFLD/NASH",
  grepl("HCC",              study_name, ignore.case = TRUE), "HCC",
  grepl("Cirrhos|Cirrh",    study_name, ignore.case = TRUE), "Cirrhosis",
  grepl("PDFF",             study_name, ignore.case = TRUE), "PDFF (MRI)",
  grepl("ALT|AST|GGT",      study_name),                     "Liver enzymes",
  default                                                   = "Other")]
reg[, trait := factor(trait, levels = c(
  "NAFLD/NASH", "HCC", "Cirrhosis", "PDFF (MRI)", "Liver enzymes", "Other"))]
reg[, ancestry := factor(ancestry, levels = c("EUR","EAS","AFR","SAS"))]
reg[, label := gsub("^20[0-9]{2}_[0-9]+_", "", study_name)]
reg[, label := gsub("_EUR|_EAS|_AFR|_SAS|_CSA", "", label)]
reg[, label := gsub("PanUKBB_", "PanUKBB ", label)]

# ---------- 2. Combined finemapping (per-study counts) ------------------
cf <- fread(file.path(BASE, "GWAS/finemapping/results/combined_finemapping.csv"),
            select = c("study","locus","susie_pip","susie_cs","susie_converged",
                       "susiex_pip","mesusie_pip"))

susie_per_study <- cf[, .(
    susie_n_loci_attempted = uniqueN(locus),
    susie_n_loci_converged = uniqueN(locus[susie_converged == TRUE]),
    susie_n_variants_pip09 = sum(susie_pip >= 0.9, na.rm = TRUE),
    susie_n_variants_pip05 = sum(susie_pip >= 0.5, na.rm = TRUE),
    susiex_n_variants_pip09 = sum(susiex_pip  >= 0.9, na.rm = TRUE),
    mesusie_n_variants_pip09 = sum(mesusie_pip >= 0.9, na.rm = TRUE)),
  by = study]
setnames(susie_per_study, "study", "study_name")

# ---------- 3. SuSiE-COLOC per-study (gene-level) ------------------------
gc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"),
  select = c("gene","coloc_best_susie_pp4","coloc_best_susie_gwas"))

coloc_per_study <- gc[!is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 >= 0.5,
  .(coloc_n_genes_pp4_05 = uniqueN(gene[coloc_best_susie_pp4 >= 0.5]),
    coloc_n_genes_pp4_08 = uniqueN(gene[coloc_best_susie_pp4 >= 0.8]),
    coloc_n_genes_pp4_09 = uniqueN(gene[coloc_best_susie_pp4 >= 0.9])),
  by = coloc_best_susie_gwas]
setnames(coloc_per_study, "coloc_best_susie_gwas", "study_name")

# ---------- 4. Merge per-GWAS counts into registry -----------------------
portfolio <- merge(reg, susie_per_study, by = "study_name", all.x = TRUE)
portfolio <- merge(portfolio, coloc_per_study, by = "study_name", all.x = TRUE)
num_cols <- c("susie_n_loci_attempted","susie_n_loci_converged",
              "susie_n_variants_pip09","susie_n_variants_pip05",
              "susiex_n_variants_pip09","mesusie_n_variants_pip09",
              "coloc_n_genes_pp4_05","coloc_n_genes_pp4_08","coloc_n_genes_pp4_09")
for (c in num_cols) portfolio[is.na(get(c)), (c) := 0]

fwrite(portfolio, file.path(OUT, "gwas_portfolio.csv"))
cat("\nPortfolio (28 GWAS):\n")
print(portfolio[, .(study_name, ancestry, trait, N_tot, N_cases,
                    susie_n_loci_converged, susie_n_variants_pip09,
                    coloc_n_genes_pp4_05)])

# ---------- 5. Method-level overall counts -------------------------------
sx_gene <- fread(file.path(BASE,
  "GWAS/finemapping/results/susiex_1kg/susiex_gene_summary_1kg.csv"))
sx_var  <- fread(file.path(BASE,
  "GWAS/finemapping/results/susiex_1kg/susiex_variant_summary_1kg.csv"))
me_gene <- fread(file.path(BASE,
  "GWAS/finemapping/results/mesusie/mesusie_gene_summary.csv"))
me_var  <- fread(file.path(BASE,
  "GWAS/finemapping/results/mesusie/mesusie_variant_summary.csv"))

total_loci_attempted   <- uniqueN(cf$locus)
total_loci_susie_conv  <- uniqueN(cf[susie_converged == TRUE, locus])
total_loci_susiex      <- uniqueN(sx_var$locus_id)
total_loci_mesusie     <- uniqueN(me_var$locus_id)

n_genes_pp4_05 <- gc[!is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 >= 0.5, uniqueN(gene)]
n_genes_pp4_08 <- gc[!is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 >= 0.8, uniqueN(gene)]
n_genes_pp4_09 <- gc[!is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 >= 0.9, uniqueN(gene)]

cascade <- data.table(
  step  = c("Lead loci screened",
            "SuSiE converged loci",
            "SuSiE-X joint loci (cross-ancestry)",
            "meSuSiE loci (multi-ethnic)",
            "Genes colocalized PP4>=0.5",
            "Genes colocalized PP4>=0.8",
            "Genes colocalized PP4>=0.9"),
  count = c(total_loci_attempted, total_loci_susie_conv,
            total_loci_susiex, total_loci_mesusie,
            n_genes_pp4_05, n_genes_pp4_08, n_genes_pp4_09))
fwrite(cascade, file.path(OUT, "method_cascade.csv"))
cat("\nCascade:\n"); print(cascade)

method_totals <- data.table(
  method = c("SuSiE (single ancestry)",
             "SuSiE-X (cross-ancestry)",
             "meSuSiE (multi-ethnic)"),
  n_loci_converged    = c(total_loci_susie_conv, total_loci_susiex, total_loci_mesusie),
  n_genes_with_pip05  = c(NA_integer_,
                          sx_gene[susiex_max_pip >= 0.5, uniqueN(GeneSymbol)],
                          me_gene[mesusie_max_pip >= 0.5, uniqueN(GeneSymbol)]),
  n_genes_with_pip09  = c(NA_integer_,
                          sx_gene[susiex_max_pip >= 0.9, uniqueN(GeneSymbol)],
                          me_gene[mesusie_max_pip >= 0.9, uniqueN(GeneSymbol)]))
fwrite(method_totals, file.path(OUT, "method_totals.csv"))
cat("\nMethod totals:\n"); print(method_totals)

# ---------- 6. eQTL panels ------------------------------------------------
broadaway_lead <- fread(file.path(BASE,
  "data/broadaway_eqtl/Liver_eQTL_Meta_Leads_ST3_20240530.tsv"))
broadaway_n_egenes <- uniqueN(broadaway_lead$Gene)

eqtl_panels <- data.table(
  panel    = c("Broadaway liver", "GTEx v8 Liver"),
  ancestry = c("EUR", "EUR"),
  n_samples = c(1183L, 208L),
  n_egenes  = c(broadaway_n_egenes, 8657L))   # GTEx liver eGenes (PredictDB)
fwrite(eqtl_panels, file.path(OUT, "eqtl_panels.csv"))
cat("\neQTL panels:\n"); print(eqtl_panels)

# ---------- 7. Ancestry-level rollup -------------------------------------
ancestry_summary <- portfolio[, .(
    n_gwas = .N,
    n_samples_total = sum(N_tot, na.rm = TRUE),
    n_loci_converged = sum(susie_n_loci_converged, na.rm = TRUE),
    n_variants_pip09 = sum(susie_n_variants_pip09, na.rm = TRUE),
    n_coloc_genes_pp4_05 = sum(coloc_n_genes_pp4_05, na.rm = TRUE)),
  by = ancestry]
fwrite(ancestry_summary, file.path(OUT, "ancestry_summary.csv"))
cat("\nAncestry rollup:\n"); print(ancestry_summary)

cat("\n[prep_data] DONE\n")
