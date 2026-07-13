#!/usr/bin/env Rscript
# ============================================================================
# 07_vacca_benchmark_concordance.R
#
# Benchmark OUR cross-species concordance against the consortium-grade murine
# model human-proximity ranking of Vacca et al. 2024 Nat Metab (LITMUS).
#   DOI 10.1038/s42255-024-01043-6 ; data data/external/vacca_2024/ (CC BY).
#
# We DO NOT re-run our cross-species pipeline or re-derive a DHPS. We compare
# our EXISTING per-diet concordance + Conserved_Core to Vacca's PUBLISHED
# DHPS rankings and 951-gene human-MASH signature.
#
#   (A) MODEL-LEVEL  (more independent): our 4 diets' direction-concordance
#       vs Vacca's per-model DHPS, aggregated to matched diet groups.
#   (B) GENE-LEVEL   (CONSISTENCY check, NOT independent — see caveat): our
#       1,355-gene Conserved_Core vs Vacca's 951-gene human-MASH signature,
#       Fisher OR + biotype-matched permutation null (mirrors audit H6).
#
# ⚠ SHARED-COHORT CAVEAT: Vacca's human benchmark = UCAM/VCU + EPoS. EPoS =
#   Govaere GSE135251, which IS one of our 5 mega cohorts. So the human side
#   is partially shared → the gene-level overlap is a CONSISTENCY check, not
#   independent validation. The model-level (mouse) axis is the more
#   independent comparison.
#
# Outputs: Analysis/Cross_Species_Concordance/results/vacca_benchmark/
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA  <- file.path(BASE, "data/external/vacca_2024")
CS     <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUTDIR <- file.path(CS, "vacca_benchmark")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

N_DRAWS <- 1000L
set.seed(42)

cat(strrep("=", 78), "\n")
cat("Vacca 2024 (LITMUS) cross-species benchmark concordance\n")
cat(strrep("=", 78), "\n")

# ============================================================================
# (A) MODEL-LEVEL: our per-diet concordance vs Vacca per-model DHPS
# ⚠ SUPERSEDED 2026-06-21: this gene-level model comparison ANTI-correlates with
#   DHPS because DHPS is a PATHWAY-level (DSEA/PEP) construct, not gene-level
#   (gene-level reproduces DHPS at only Spearman 0.22 vs 0.73 for the pathway
#   metric — scripts 08/09/10). The figures (figS06_vacca_benchmark, fig4g) now
#   use the PATHWAY-level proximity from script 10. This block is RETAINED only as
#   the record of the gene-level attempt; do NOT cite its model ranking.
# ============================================================================
cat("\n--- (A) Model-level (SUPERSEDED — gene-level; see 08/09/10 for pathway) ---\n")

# Vacca MOESM9: per-model PHPS/HHPS/DHPS/MHPS in two arms. 2-row header →
# read raw and assign by position (cols: 1 Model, 2 DietGroup,
# 3-6 metabolic PHPS/HHPS/DHPS/MHPS, 7-10 fibrosis PHPS/HHPS/DHPS/MHPS, ...).
m9 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM9_ESM.xlsx"),
                               col_names = FALSE, skip = 2))
setnames(m9, seq_len(ncol(m9)),
  c("model", "diet_group", "m_PHPS", "m_HHPS", "m_DHPS", "m_MHPS",
    "f_PHPS", "f_HHPS", "f_DHPS", "f_MHPS",
    head(c("BW_l2fc", "BW_p", "fibrosis_avg"), ncol(m9) - 10)))
m9 <- m9[!is.na(model)]
m9 <- m9[!grepl("^R-", model)]   # drop RAT models (R-CDAA…); compare mouse-only
for (cc in c("m_DHPS", "f_DHPS")) m9[[cc]] <- as.numeric(m9[[cc]])

# Map Vacca "Diet Group" → OUR diet. Transparent, auditable:
#   MCD   ← CDD            (choline-deficient: DB-MCD, 6J-MCD)
#   CDAHFD← CDHFD          (choline-deficient high-fat)
#   HFD   ← HFD
#   FPC   ← WD*            (Western incl. GAN/AMLN, all coded WD0.2%/WD2%/…)
# Unmapped (CDHFHCD, CCl4, GMCHOW, HFDAMLD, rat CDAA) → reported, excluded.
m9[, our_diet := fcase(
  diet_group == "CDD",                  "MCD",
  grepl("^CDHFD$", diet_group),         "CDAHFD",
  diet_group == "HFD",                  "HFD",
  grepl("^WD", diet_group),             "FPC",
  default = NA_character_)]
cat("  Vacca diet-group → our-diet mapping (model membership):\n")
print(m9[!is.na(our_diet), .(n_models = .N,
        models = paste(model, collapse = "; ")), by = .(our_diet, diet_group)][order(our_diet)])
cat("  Vacca models NOT mapped to our 4 diets:\n")
print(m9[is.na(our_diet), .(diet_group, model)])

# Aggregate to our-diet groups. Mapping is many-to-one (esp. FPC ← 20 WD models),
# so keep the DHPS RANGE, not just the mean, for honest presentation.
vacca_diet <- m9[!is.na(our_diet), .(
  vacca_n_models       = .N,
  vacca_DHPS_metab     = mean(m_DHPS, na.rm = TRUE),
  vacca_DHPS_metab_min = min(m_DHPS,  na.rm = TRUE),
  vacca_DHPS_metab_max = max(m_DHPS,  na.rm = TRUE),
  vacca_DHPS_fibro     = mean(f_DHPS, na.rm = TRUE)), by = our_diet]

# Our per-diet concordance (disease-vs-ctrl + nafl-vs-nash anchors)
mx <- fread(file.path(CS, "gene_concordance_matrix_20x.csv"))
our_diet <- dcast(mx[human_signature %in% c("disease_vs_ctrl", "nafl_vs_nash"),
                     .(diet, human_signature, concordance_pct, rho_sig)],
                  diet ~ human_signature, value.var = c("concordance_pct", "rho_sig"))
setnames(our_diet, "diet", "our_diet")

# Our DSEA-ANALOG: GSEA NES of the human disease signature in each diet's ranked
# transcriptome — the metric fig4g (crossspecies_conservation.R) plots and the
# DIRECT methodological analog of Vacca's DHPS/DSEA (replicated here so we compare
# like-with-like, NOT the both-significant direction-concordance which differs).
suppressPackageStartupMessages(library(fgsea))
ANNOT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation")
INT_I <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
MPD   <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
ortho <- fread(file.path(ANNOT, "ortholog_mapping.tsv"))
humd  <- fread(file.path(INT_I, "canonical_deg_results.csv"))
humd[, gene_base := gsub("\\..*", "", gene)]
hsym  <- merge(humd[, .(gene_base, h_lfc = logFC, h_padj = padj)],
               ortho[, .(human_gene_id, mouse_gene_id)],
               by.x = "gene_base", by.y = "human_gene_id")
gs <- list(human_up   = unique(hsym[h_padj < 0.05 & h_lfc >  0.5, mouse_gene_id]),
           human_down = unique(hsym[h_padj < 0.05 & h_lfc < -0.5, mouse_gene_id]))
nes_rows <- rbindlist(lapply(c("MCD", "HFD", "CDAHFD", "FPC"), function(d) {
  md <- fread(file.path(MPD, paste0(d, "_de_results.csv")))
  md[, mouse_base := gsub("\\..*", "", gene)]; md <- md[is.finite(t)]
  st <- md[, setNames(t, mouse_base)]; st <- st[!duplicated(names(st))]
  fg <- as.data.table(suppressWarnings(fgsea(gs, st, scoreType = "std", nPermSimple = 10000)))
  data.table(our_diet = d, nes_up = fg[pathway == "human_up", NES],
             nes_down = fg[pathway == "human_down", NES])
}))
# DHPS-analog = mean(NES_up, -NES_down): high when human-UP genes are up AND
# human-DOWN genes are down in the model (same construction as Vacca's DSEA).
nes_rows[, our_dsea_nes := (nes_up - nes_down) / 2]

model_level <- Reduce(function(a, b) merge(a, b, by = "our_diet", all.x = TRUE),
                      list(our_diet, vacca_diet, nes_rows))
setcolorder(model_level, "our_diet")
model_level <- model_level[order(-vacca_DHPS_metab)]
cat("\n  MODEL-LEVEL table (our GSEA DSEA-NES vs Vacca DHPS):\n")
print(model_level[, .(our_diet, our_dsea_nes = round(our_dsea_nes, 2),
                      nes_up = round(nes_up, 2), nes_down = round(nes_down, 2),
                      vacca_DHPS_metab = round(vacca_DHPS_metab, 2),
                      vacca_DHPS_fibro = round(vacca_DHPS_fibro, 2))])

# Agreement: our DSEA-NES vs Vacca DHPS (the matched DSEA-style metric)
sp_metab <- suppressWarnings(cor(model_level$our_dsea_nes, model_level$vacca_DHPS_metab,
                                 method = "spearman", use = "complete.obs"))
sp_fibro <- suppressWarnings(cor(model_level$our_dsea_nes, model_level$vacca_DHPS_fibro,
                                 method = "spearman", use = "complete.obs"))
pe_metab <- suppressWarnings(cor(model_level$our_dsea_nes, model_level$vacca_DHPS_metab,
                                 method = "pearson", use = "complete.obs"))
cat(sprintf("\n  Agreement (4 diets) — our DSEA-NES vs Vacca DHPS:\n    Spearman vs metabolic = %.2f ; vs fibrotic = %.2f ; Pearson(metab) = %.2f\n",
            sp_metab, sp_fibro, pe_metab))

mcd_rank <- which(model_level[order(-our_dsea_nes)]$our_diet == "MCD")
mcd <- model_level[our_diet == "MCD"]
cat(sprintf("\n  MCD: our DSEA-NES = %.2f (rank %d/4); Vacca metabolic DHPS = %.2f (lowest of 4).\n    MCD strongly perturbs the disease transcriptome yet is the least metabolically\n    human-proximal model (lean, non-obese) — translatability caveat for our MCD-heavy pool.\n",
            mcd$our_dsea_nes, mcd_rank, mcd$vacca_DHPS_metab))
cat("  CAVEAT: our CDAHFD + FPC are both arms of GSE162876 (shared controls); not independent.\n")

fwrite(model_level, file.path(OUTDIR, "model_level_concordance.csv"))

# ============================================================================
# (B) GENE-LEVEL: Conserved_Core vs Vacca 951-gene human-MASH signature
# ============================================================================
cat("\n--- (B) Gene-level: Conserved_Core vs Vacca human-MASH signature ---\n")

s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"),
                               sheet = "Table S4"))
setnames(s4,
  c("GeneSymbol", "Gene_used_in_DSEA(0:No/1:Yes)",
    "Early/All/Late disease stage (if used in DSEA otherwise NA)"),
  c("gene", "in_dsea", "vacca_stage"), skip_absent = TRUE)
# Representative human progression logFC = mean(UCAM/VCU + EPoS Severe-vs-Mild).
# read_excel parses these as character (mixed cells) → coerce numeric first.
hcols <- c("L2FC_UCAM/VCU: Severe vs Mild", "L2FC_EPoS: Severe vs Mild")
s4[, (hcols) := lapply(.SD, as.numeric), .SDcols = hcols]
s4[, vacca_human_logfc := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = hcols]
s4 <- s4[!is.na(gene) & gene != ""][!duplicated(gene)]

vacca_sig  <- s4[in_dsea == 1]$gene                                   # 951
vacca_prog <- s4[in_dsea == 1 & vacca_stage == "Disease progression"]$gene  # 526
cat(sprintf("  Vacca human-MASH signature: %d genes (DSEA); progression subset %d.\n",
            length(vacca_sig), length(vacca_prog)))

# Our Conserved_Core (NAFL-vs-NASH anchor, >=3/4 diets) + human logFC
conc <- fread(file.path(CS, "concordance_atlas_unified.csv"),
              select = c("human_symbol", "primary_category", "mean_h_lfc"))
conc <- conc[!is.na(human_symbol) & human_symbol != ""][!duplicated(human_symbol)]
cc_genes <- conc[primary_category == "Conserved"]$human_symbol
cat(sprintf("  Our Conserved_Core: %d genes.\n", length(cc_genes)))

# Biotype for the matched null (from multi-evidence atlas, as audit H6 does)
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("human_symbol", "gene_biotype"))
atlas <- atlas[!is.na(human_symbol) & human_symbol != ""][!duplicated(human_symbol)]
collapse_biotype <- function(bt) {
  bt2 <- rep("other", length(bt)); bt2[is.na(bt)] <- "unknown"
  bt2[bt == "protein_coding"] <- "protein_coding"
  bt2[bt == "lncRNA"] <- "lncRNA"
  bt2[grepl("pseudogene", bt, ignore.case = TRUE)] <- "pseudogene"; bt2 }
atlas[, biotype_class := collapse_biotype(gene_biotype)]

# Fair universe = genes testable in BOTH pipelines (our concordance atlas
# ∩ Vacca's S4 tested universe).
universe <- intersect(conc$human_symbol, s4$gene)
cat(sprintf("  Shared testable universe (our concordance atlas ∩ Vacca S4): %d genes.\n",
            length(universe)))

# ---- Fisher OR + biotype-matched null (mirrors RNA-seq/audit/H6) ----------
fisher_or <- function(drawn, target, universe) {
  drawn <- intersect(drawn, universe); target_u <- intersect(target, universe)
  a <- length(intersect(drawn, target_u)); b <- length(drawn) - a
  cc <- length(target_u) - a; d <- length(universe) - a - b - cc
  m <- matrix(c(a, b, cc, d), nrow = 2); if (any(c(a, b, cc, d) == 0)) m <- m + 0.5
  ft <- tryCatch(fisher.test(m), error = function(e) NULL)
  if (is.null(ft)) return(list(OR = NA_real_, p = NA_real_, a = a, b = b, c = cc, d = d))
  list(OR = unname(ft$estimate), p = ft$p.value, a = a, b = b, c = cc, d = d)
}
cc_in_u   <- intersect(cc_genes, universe)
cc_bt_vec <- atlas[human_symbol %in% cc_in_u]$biotype_class
u_bt      <- atlas[human_symbol %in% universe, .(human_symbol, biotype_class)]
# universe genes missing from biotype atlas → "unknown"
miss_u <- setdiff(universe, u_bt$human_symbol)
if (length(miss_u)) u_bt <- rbind(u_bt, data.table(human_symbol = miss_u, biotype_class = "unknown"))

run_enrich <- function(target_set, label) {
  obs <- fisher_or(cc_in_u, target_set, universe)
  k <- length(cc_in_u)
  null_A <- numeric(N_DRAWS); null_B <- numeric(N_DRAWS)
  bt_tbl <- table(cc_bt_vec)
  for (i in seq_len(N_DRAWS)) {
    null_A[i] <- fisher_or(sample(universe, k), target_set, universe)$OR
    drawn <- character(0)
    for (bt in names(bt_tbl)) {
      pool <- u_bt[biotype_class == bt]$human_symbol; kbt <- as.integer(bt_tbl[[bt]])
      drawn <- c(drawn, if (length(pool) < kbt) sample(pool, kbt, TRUE) else sample(pool, kbt))
    }
    null_B[i] <- fisher_or(drawn, target_set, universe)$OR
  }
  data.table(comparator = label, universe = "concordance_atlas∩vacca_S4",
    observed_OR = obs$OR, fisher_p = obs$p, a = obs$a, b = obs$b, c = obs$c, d = obs$d,
    null_A_median = median(null_A, na.rm = TRUE), null_A_p = mean(null_A >= obs$OR, na.rm = TRUE),
    null_B_median = median(null_B, na.rm = TRUE), null_B_p = mean(null_B >= obs$OR, na.rm = TRUE),
    n_draws = N_DRAWS, k_conserved = k, n_target = length(intersect(target_set, universe)))
}
enr <- rbindlist(list(
  run_enrich(vacca_sig,  "vacca_human_mash_951"),
  run_enrich(vacca_prog, "vacca_progression_526")))
cat("\n  GENE-LEVEL enrichment (Conserved_Core in Vacca signature):\n"); print(enr)
fwrite(enr, file.path(OUTDIR, "gene_level_enrichment.csv"))

# ---- Direction concordance + per-gene overlap table -----------------------
ov <- merge(conc[, .(human_symbol, our_mean_h_lfc = mean_h_lfc,
                     in_conserved_core = primary_category == "Conserved")],
            s4[, .(human_symbol = gene, in_vacca_signature = in_dsea == 1,
                   vacca_stage, vacca_human_logfc)],
            by = "human_symbol", all = TRUE)
ov <- ov[human_symbol %in% universe]
ov[is.na(in_conserved_core),  in_conserved_core  := FALSE]
ov[is.na(in_vacca_signature), in_vacca_signature := FALSE]
ov[, direction_match := sign(our_mean_h_lfc) == sign(vacca_human_logfc)]

shared <- ov[in_conserved_core & in_vacca_signature & !is.na(direction_match)]
dir_pct <- 100 * mean(shared$direction_match, na.rm = TRUE)
cat(sprintf("\n  Overlap: Conserved_Core ∩ Vacca-951 = %d genes (of %d core in universe).\n",
            nrow(ov[in_conserved_core & in_vacca_signature]), nrow(ov[in_conserved_core == TRUE])))
cat(sprintf("  Conserved_Core genes BEYOND Vacca signature (our extension): %d.\n",
            nrow(ov[in_conserved_core & !in_vacca_signature])))
cat(sprintf("  Direction concordance of shared genes (human logFC sign): %.1f%% (n=%d).\n",
            dir_pct, nrow(shared)))
fwrite(ov[order(-in_conserved_core, -in_vacca_signature)],
       file.path(OUTDIR, "conserved_core_vs_vacca_overlap.csv"))

# Spot-check anchors
cat("\n  Spot-check anchors:\n")
print(ov[human_symbol %in% c("THRB", "RORA", "FGF21", "HKDC1", "COL1A1"),
         .(human_symbol, in_conserved_core, in_vacca_signature, vacca_stage,
           our_mean_h_lfc = round(our_mean_h_lfc, 2),
           vacca_human_logfc = round(vacca_human_logfc, 2), direction_match)])

cat("\n", strrep("=", 78), "\n")
cat("DONE. Outputs in ", OUTDIR, "\n")
cat("⚠ Gene-level = CONSISTENCY check (Vacca EPoS = our Govaere cohort, shared human input).\n")
cat("  Model-level (mouse) = the more independent axis.\n")
cat(strrep("=", 78), "\n")
