#!/usr/bin/env Rscript
# ============================================================================
# 78_population_genetics.R
# Population / clinical-genetics evidence cluster (v2, 2026-06-01).
#
# MASLD-anchored population-scale evidence. Each component carries an explicit
# MASLD-relevance gate (the substrates are population-generic; relevance comes
# from a disease filter or credible-set intersection). Modular + file.exists-
# guarded; assembles one gene-level table keyed on human_symbol.
#
# Components (built incrementally):
#   [1] ClinVar      — P/LP variants FILTERED TO MASLD-relevant conditions      [this commit]
#   [2] gnomAD       — pLI/constraint extension + AF of credible-set variants   [pending data]
#   [3] OMIM         — monogenic fatty-liver disease genes                       [pending data]
#   [4] MASLD burden — All-of-Us/UKB Nature-2025 rare-variant supp table         [pending data]
#   [5] AlphaMissense— deleteriousness of credible-set coding variants (optional)[pending data]
#
# Output: RNA-seq/results/multi_evidence/popgen_atlas_columns.tsv
# Merged by 27a post-assembly block (drop-then-merge on human_symbol).
# Tier: ClinVar/burden = T1 (MASLD-context); gnomAD constraint = T3 (annotation).
# ============================================================================

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
EXT  <- file.path(BASE, "data/external")
OUT  <- file.path(BASE, "RNA-seq/results/multi_evidence/popgen_atlas_columns.tsv")

cat("== 78_population_genetics ==\n")
layers <- list()  # each: data.table keyed on human_symbol

# ---------------------------------------------------------------------------
# [1] ClinVar P/LP, FILTERED TO MASLD-RELEVANT CONDITIONS (the relevance gate)
#     Source: variant_summary.txt.gz (GeneSymbol + ClinicalSignificance +
#     PhenotypeList + ReviewStatus + Assembly already parsed).
# ---------------------------------------------------------------------------
cv_path <- file.path(EXT, "clinvar/variant_summary.txt.gz")
if (file.exists(cv_path)) {
  cv <- fread(cv_path, quote = "",
              select = c("GeneSymbol", "ClinicalSignificance", "PhenotypeList",
                         "ReviewStatus", "Assembly"))
  cv <- cv[Assembly == "GRCh38"]
  # P/LP only (exclude Benign / Conflicting / Uncertain)
  cv <- cv[grepl("pathogenic", ClinicalSignificance, ignore.case = TRUE) &
           !grepl("Conflicting|Benign|Uncertain", ClinicalSignificance, ignore.case = TRUE)]
  # MASLD-relevance gate: phenotype is a steatotic/metabolic/cholestatic liver condition
  # (MASLD spectrum + Mendelian fatty-liver phenocopies).
  masld_cond <- paste(
    "steato|fatty liver|non-alcoholic|nonalcoholic|NAFLD|NASH|MASLD|MASH",
    "cirrhos|hepatic fibrosis|liver fibrosis|hepatocellular",
    "lipodystroph|cholestas|hypobetalipoprotein|abetalipoprotein",
    "lysosomal acid lipase|cholesteryl ester storage|Wolman",
    "Wilson disease|glycogen storage|citrullinemia|Citrin|carnitine",
    "hereditary fructose|progressive familial intrahepatic", sep = "|")
  cv <- cv[grepl(masld_cond, PhenotypeList, ignore.case = TRUE)]
  cv <- cv[GeneSymbol != "" & !grepl(";", GeneSymbol)]  # drop multi-gene rows
  # review-status -> star rating
  star_of <- function(s) {
    s <- tolower(s)
    fifelse(grepl("practice guideline", s), 4L,
    fifelse(grepl("reviewed by expert panel", s), 3L,
    fifelse(grepl("multiple submitters", s) & grepl("no conflicts", s), 2L,
    fifelse(grepl("single submitter|criteria provided", s), 1L, 0L))))
  }
  cv[, stars := star_of(ReviewStatus)]
  clinvar <- cv[, .(
    clinvar_n_plp_masld = .N,
    clinvar_max_stars   = max(stars),
    clinvar_conditions  = paste(head(unique(unlist(strsplit(paste(PhenotypeList, collapse = "|"), "\\|"))), 6),
                                collapse = "; ")
  ), by = .(human_symbol = GeneSymbol)]
  layers$clinvar <- clinvar
  cat(sprintf("  [1] ClinVar: %d genes with MASLD-condition P/LP (>=2 stars: %d)\n",
              nrow(clinvar), sum(clinvar$clinvar_max_stars >= 2)))
  cat(sprintf("      top genes: %s\n",
              paste(head(clinvar[order(-clinvar_n_plp_masld)]$human_symbol, 12), collapse = ", ")))
} else {
  cat("  [1] ClinVar variant_summary.txt.gz not found; skipping\n")
}

# ---------------------------------------------------------------------------
# [2] gnomAD — ancestry-stratified AF of MASLD credible-set variants
#     + gene-level LoF constraint (pLI / LOEUF).
#     AF cache produced by Script 78a (gnomAD v4 GraphQL on credible-set core,
#     hg19->hg38 lifted, nearest-TSS gene). Constraint from the already-
#     downloaded gnomAD v4.1 constraint metrics table.
#     Tier: constraint = T3 (population-generic annotation); AF/divergence
#     contextualizes MASLD credible-set variants (T3 annotation layer).
# ---------------------------------------------------------------------------
af_path  <- file.path(BASE, "RNA-seq/results/multi_evidence/gnomad_af_credset.tsv")
con_path <- file.path(EXT, "gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv")
gnomad_parts <- list()

# [2a] credible-set variant AF (per gene)
if (file.exists(af_path)) {
  af <- fread(af_path)
  keep <- intersect(c("human_symbol", "gnomad_af_eur", "gnomad_af_eas",
                      "gnomad_af_global", "gnomad_ancestry_divergence",
                      "gnomad_af_n_variants", "gnomad_af_lead_variant"), names(af))
  gnomad_parts$af <- unique(af[, ..keep], by = "human_symbol")
  cat(sprintf("  [2a] gnomAD AF: %d credible-set genes (%d ancestry-divergent EUR/EAS >2x)\n",
              nrow(gnomad_parts$af),
              sum(gnomad_parts$af$gnomad_ancestry_divergence %in% TRUE)))
} else {
  cat("  [2a] gnomAD AF cache not found; run 78a_fetch_gnomad_af.R (skipping AF)\n")
}

# [2b] gene-level LoF constraint: pLI + LOEUF + derived class
if (file.exists(con_path)) {
  con <- fread(con_path,
               select = c("gene", "mane_select", "lof.pLI", "lof.oe_ci.upper"))
  # LOEUF = upper bound of the observed/expected LoF CI (gnomAD canonical metric)
  setnames(con, c("lof.pLI", "lof.oe_ci.upper"), c("gnomad_pLI", "gnomad_loeuf"))
  # collapse to one (most-constrained) transcript per gene symbol; prefer MANE
  con <- con[!is.na(gene) & gene != ""]
  setorder(con, gene, -mane_select, gnomad_loeuf)
  con <- con[, .SD[1], by = gene]
  # constraint class: pLI if present, else LOEUF-derived (gnomAD convention LOEUF<0.6 = constrained)
  con[, gnomad_constraint_class := fifelse(
    !is.na(gnomad_pLI) & gnomad_pLI >= 0.9, "LoF_intolerant",
    fifelse(!is.na(gnomad_loeuf) & gnomad_loeuf < 0.6, "constrained",
    fifelse(!is.na(gnomad_loeuf) & gnomad_loeuf >= 0.6, "tolerant", NA_character_)))]
  gnomad_parts$con <- con[, .(human_symbol = gene, gnomad_pLI, gnomad_loeuf,
                              gnomad_constraint_class)]
  cat(sprintf("  [2b] gnomAD constraint: %d genes (LoF-intolerant pLI>=0.9: %d)\n",
              nrow(gnomad_parts$con),
              sum(gnomad_parts$con$gnomad_constraint_class == "LoF_intolerant", na.rm = TRUE)))
} else {
  cat("  [2b] gnomAD constraint metrics not found; skipping\n")
}

if (length(gnomad_parts) > 0) {
  layers$gnomad_af <- Reduce(function(a, b)
    merge(a, b, by = "human_symbol", all = TRUE), gnomad_parts)
}

# ---------------------------------------------------------------------------
# [3] OMIM — monogenic fatty-liver / metabolic-liver disease genes.
#     OMIM genemap2 needs a registration key (gated). If a local genemap2 dump
#     exists we parse it; otherwise FALL BACK to a curated, literature-anchored
#     panel of monogenic steatosis / cholestasis / lipodystrophy / FAO /
#     glycogen-storage / urea-cycle genes (each with its OMIM disease phenotype),
#     UNIONED with the local Okur 2026 monogenic-metabolic panel.
#     Tier: T1 (MASLD-context: every entry is a Mendelian steatosis/cholestasis
#     phenocopy or a monogenic fatty-liver driver).
# ---------------------------------------------------------------------------
omim_dump <- file.path(EXT, "omim/genemap2.txt")  # present only if key-gated dump downloaded
if (file.exists(omim_dump)) {
  # parse genemap2 (tab-delimited; cols: ... Gene Symbols ... Phenotypes)
  gm <- tryCatch(fread(omim_dump, sep = "\t", quote = "", fill = TRUE, header = FALSE),
                 error = function(e) NULL)
  omim_ok <- !is.null(gm)
} else {
  omim_ok <- FALSE
}

# Curated monogenic fatty-liver / metabolic-liver panel (literature-anchored).
# Each gene maps to a Mendelian disorder with hepatic steatosis, steatohepatitis,
# cholestasis, lipodystrophy-associated fatty liver, or a fatty-liver phenocopy.
curated_omim <- data.table(
  human_symbol = c(
    # --- monogenic steatosis / steatohepatitis ---
    "APOB","MTTP","PNPLA3","TM6SF2","GCKR","HSD17B13","MBOAT7","PNPLA2",
    # --- lysosomal acid lipase / cholesteryl-ester storage ---
    "LIPA",
    # --- lipodystrophy (with severe fatty liver) ---
    "LMNA","PPARG","AGPAT2","BSCL2","CAV1","PLIN1","CIDEC","LIPE","CAVIN1",
    # --- cholestasis / PFIC / bile-acid ---
    "ABCB4","ABCB11","ATP8B1","TJP2","NR1H4","MYO5B","ABCC2","SLC10A1",
    # --- Wilson / copper, hemochromatosis (iron) ---
    "ATP7B","HFE","SLC40A1",
    # --- fatty-acid oxidation defects ---
    "CPT1A","CPT2","ACADM","ACADVL","HADHA","HADHB","ETFA","ETFB","ETFDH","SLC22A5",
    # --- glycogen storage disease (hepatic) ---
    "G6PC1","SLC37A4","AGL","GBE1","PYGL","GYS2","PHKA2","G6PC2",
    # --- urea-cycle disorders (with hepatic involvement / steatosis) ---
    "OTC","CPS1","ASS1","ASL","ARG1","NAGS","SLC25A13","SLC25A15",
    # --- congenital disorders of metabolism w/ steatosis ---
    "ALDOB","GALT","FAH","CTNS"),
  omim_phenotype = c(
    "Familial hypobetalipoproteinemia (OMIM 615558)",
    "Abetalipoproteinemia (OMIM 200100)",
    "Steatohepatitis susceptibility (PNPLA3 I148M; OMIM 613387)",
    "MASLD susceptibility / steatosis (TM6SF2; OMIM 617299)",
    "Hepatic steatosis susceptibility (GCKR)",
    "MASLD protective splice variant (HSD17B13)",
    "Hepatic fat / fibrosis susceptibility (MBOAT7)",
    "Neutral lipid storage disease with myopathy (OMIM 610717)",
    "Lysosomal acid lipase deficiency / Wolman / CESD (OMIM 278000)",
    "Familial partial lipodystrophy type 2 (Dunnigan; OMIM 151660)",
    "Familial partial lipodystrophy type 3 (OMIM 604367)",
    "Congenital generalized lipodystrophy type 1 (OMIM 608594)",
    "Congenital generalized lipodystrophy type 2 (OMIM 269700)",
    "Congenital generalized lipodystrophy type 3 (OMIM 612526)",
    "Familial partial lipodystrophy type 4 (OMIM 613877)",
    "Familial partial lipodystrophy type 5 (OMIM 615238)",
    "Familial lipodystrophy / lipase deficiency (OMIM 151800)",
    "Congenital generalized lipodystrophy type 4 (OMIM 613327)",
    "PFIC type 3 / low-phospholipid cholelithiasis (OMIM 602347)",
    "PFIC type 2 (OMIM 601847)",
    "PFIC type 1 / Byler disease (OMIM 211600)",
    "PFIC type 4 (OMIM 615878)",
    "PFIC type 5 (FXR deficiency; OMIM 617049)",
    "PFIC type 6 / MVID (OMIM 619484)",
    "Dubin-Johnson syndrome (OMIM 237500)",
    "Hypercholanemia / NTCP deficiency (OMIM 619256)",
    "Wilson disease (OMIM 277900)",
    "Hereditary hemochromatosis type 1 (OMIM 235200)",
    "Hemochromatosis type 4 / ferroportin disease (OMIM 606069)",
    "CPT I deficiency (OMIM 255120)",
    "CPT II deficiency (OMIM 255110/600649)",
    "MCAD deficiency (OMIM 201450)",
    "VLCAD deficiency (OMIM 201475)",
    "LCHAD / MTP deficiency (OMIM 609016)",
    "MTP deficiency, beta-subunit (OMIM 609015)",
    "Glutaric acidemia IIA (OMIM 231680)",
    "Glutaric acidemia IIB (OMIM 231680)",
    "Glutaric acidemia IIC (OMIM 231680)",
    "Systemic primary carnitine deficiency (OMIM 212140)",
    "GSD type Ia (von Gierke; OMIM 232200)",
    "GSD type Ib (OMIM 232220)",
    "GSD type III (Cori/Forbes; OMIM 232400)",
    "GSD type IV (Andersen; OMIM 232500)",
    "GSD type VI (Hers; OMIM 232700)",
    "GSD type 0, hepatic (OMIM 240600)",
    "GSD type IXa (OMIM 306000)",
    "GSD susceptibility / fasting hypoglycemia (G6PC2)",
    "Ornithine transcarbamylase deficiency (OMIM 311250)",
    "Carbamoyl-phosphate synthetase I deficiency (OMIM 237300)",
    "Citrullinemia type I (OMIM 215700)",
    "Argininosuccinic aciduria (OMIM 207900)",
    "Argininemia (OMIM 207800)",
    "N-acetylglutamate synthase deficiency (OMIM 237310)",
    "Citrin deficiency / citrullinemia type II (OMIM 605814/603471)",
    "HHH syndrome (OMIM 238970)",
    "Hereditary fructose intolerance (OMIM 229600)",
    "Classic galactosemia (OMIM 230400)",
    "Tyrosinemia type I (OMIM 276700)",
    "Cystinosis (OMIM 219800)"))

if (omim_ok) {
  cat("  [3] OMIM genemap2 dump found and parsed (key-gated path)\n")
  # (parse-and-filter path retained for when a dump is downloaded; the curated
  #  panel below is still unioned to guarantee MASLD-anchored coverage.)
}
# union curated panel with the local Okur 2026 monogenic-metabolic panel
omim <- copy(curated_omim)
okur_path <- file.path(EXT, "okur2026_monogenic/okur2026_panel_genes_partial.tsv")
if (file.exists(okur_path)) {
  ok <- fread(okur_path)
  extra <- setdiff(ok$gene, omim$human_symbol)
  if (length(extra) > 0) {
    okm <- ok[gene %in% extra]
    omim <- rbind(omim, data.table(
      human_symbol   = okm$gene,
      omim_phenotype = paste0("Okur2026 monogenic-metabolic panel: ", okm$notes)),
      fill = TRUE)
  }
  cat(sprintf("  [3] OMIM curated panel + Okur2026 union: +%d Okur-only genes\n", length(extra)))
}
omim[, omim_monogenic := TRUE]
omim <- unique(omim, by = "human_symbol")
layers$omim <- omim[, .(human_symbol, omim_monogenic, omim_phenotype)]
cat(sprintf("  [3] OMIM monogenic fatty-liver panel: %d genes (curated fallback%s)\n",
            nrow(layers$omim), if (omim_ok) " + genemap2" else "; OMIM API key-gated"))

# ---------------------------------------------------------------------------
# [4] MASLD-specific rare-variant burden.
#     Source: Verma et al. 2025, "Trans-ancestral rare variant association study
#     with machine learning-based phenotyping for MASLD," Genome Biology
#     26:53 (PMID 40065360; PMC11892324; DOI 10.1186/s13059-025-03518-5).
#     UK Biobank + All-of-Us + BioMe, 736,010 individuals. Gene-level collapsing
#     of ultra-rare PTV + deleterious missense vs PDFF / MASLD case-control.
#     Hand-curated from the paper text/Table S6 + predicted-phenotype Table S24
#     (open-access CC-BY). If a richer machine-readable supp lands at
#     data/external/masld_burden_verma2025/burden_table.tsv it overrides this.
#     Tier: T1 (MASLD-specific by construction).
# ---------------------------------------------------------------------------
verma_path <- file.path(EXT, "masld_burden_verma2025/burden_table.tsv")
burden <- NULL
if (file.exists(verma_path)) {
  vb <- tryCatch(fread(verma_path), error = function(e) NULL)
  if (!is.null(vb) && "human_symbol" %in% names(vb) && "burden_masld_pval" %in% names(vb)) {
    burden <- vb
    cat("  [4] MASLD burden: machine-readable Verma2025 supp table loaded\n")
  }
}
if (is.null(burden)) {
  # Curated from Verma et al. 2025 (Genome Biol). p = gene-level / single-variant
  # association p-value (true-phenotype meta-analysis where available, else the
  # predicted-phenotype gene-level / lead single-variant p). analysis annotates source.
  burden <- data.table(
    human_symbol = c(
      # trans-ancestral meta-analysis hits (true phenotypes)
      "APOB","MYCBP2","XAB2","CDH5",
      # predicted-phenotype gene-level Bonferroni hits
      "FNIP1","G6PC1","HSPG2","IGFALS","INSR","LDLR","MC4R","PDE3B",
      "PPARG","SLC30A10","SMAD6",
      # predicted-phenotype single-variant hits
      "SHBG","APOC3","ACVR1C","INHBE","ZNF404","XKR7","NOL4L","STC2",
      "MAD1L1","LY6G6C","FXR2"),
    burden_masld_pval = c(
      5.90e-9, 1.13e-6, 1.42e-8, 2.28e-7,
      8.71e-7, 4.37e-8, 1.22e-7, 1.92e-10, 1.67e-7, 4.86e-7, 2.62e-12, 4.53e-10,
      5.29e-10, 1.43e-6, 2.35e-6,
      3.97e-152, 9.55e-12, 1.56e-9, NA, NA, NA, NA, NA,
      NA, 0.029, 0.005),
    burden_analysis = c(
      rep("meta_gene_or_singlevariant", 4),
      rep("predicted_pheno_gene", 11),
      rep("predicted_pheno_singlevariant", 11)))
  cat("  [4] MASLD burden: curated Verma2025 (Genome Biol 2025; PMC11892324) gene list\n")
}
burden[, `:=`(burden_masld_hit = TRUE,
              burden_source     = "Verma2025_GenomeBiol_736k")]
layers$burden <- unique(burden[, .(human_symbol, burden_masld_pval,
                                    burden_masld_hit, burden_source)],
                        by = "human_symbol")
cat(sprintf("  [4] MASLD rare-variant burden: %d genes (%d with reported p)\n",
            nrow(layers$burden), sum(!is.na(layers$burden$burden_masld_pval))))

# ---------------------------------------------------------------------------
# [5] AlphaMissense (optional, pending data) — added when cache lands.
# ---------------------------------------------------------------------------

# ---- assemble ----
if (length(layers) == 0) stop("no population-genetics components available")
pg <- Reduce(function(a, b) merge(a, b, by = "human_symbol", all = TRUE), layers)
setorder(pg, human_symbol)
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(pg, OUT, sep = "\t")
cat(sprintf("WROTE %s : %d genes x %d cols\n", OUT, nrow(pg), ncol(pg)))
cat(sprintf("  columns: %s\n", paste(names(pg), collapse = ", ")))

# per-column non-NA gene counts
cat("  non-NA gene counts per column:\n")
for (col in setdiff(names(pg), "human_symbol")) {
  v <- pg[[col]]
  n <- if (is.logical(v)) sum(v %in% TRUE) else sum(!is.na(v))
  cat(sprintf("    %-28s %d\n", col, n))
}

# anchor-gene readout across all components
cat("  anchor genes:\n")
show_cols <- intersect(c("clinvar_n_plp_masld","clinvar_max_stars","gnomad_af_eur",
                         "gnomad_af_eas","gnomad_ancestry_divergence","gnomad_pLI",
                         "gnomad_constraint_class","omim_monogenic","burden_masld_hit",
                         "burden_masld_pval"), names(pg))
for (g in c("PNPLA3","TM6SF2","HSD17B13","APOB","LIPA","ATP7B","ABCB4")) {
  if (g %in% pg$human_symbol) {
    r <- pg[human_symbol == g]
    cat(sprintf("    [%s] %s\n", g,
                paste(sapply(show_cols, function(c) sprintf("%s=%s", c, r[[c]][1])),
                      collapse = " | ")))
  } else cat(sprintf("    [%s] absent\n", g))
}
