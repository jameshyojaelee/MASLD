#!/usr/bin/env Rscript
# 03_aggregate_and_deg_overlap.R
# ----------------------------------------------------------------------------
# Join variants_long ⨝ variant_classification, attach DEG context (Tier 1) and
# GWAS-ATAC overlay, emit master CSVs + per-panel summary tables.
suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT  <- file.path(BASE, "RNA-seq/results/coloc_variant_classes")

vlong <- fread(file.path(OUT, "variants_long.csv"))
ann   <- fread(file.path(OUT, "variant_classification.csv"))
cs_members <- fread(file.path(OUT, "cs_members_long.csv"))

cat(sprintf("variants_long: %d rows | classification: %d rows\n",
            nrow(vlong), nrow(ann)))

# Inner-join classification onto every variant_long row (keeps definition rows)
v <- merge(vlong, ann, by = "variant_key", all.x = TRUE)

# ---- DEG (Tier 1) -----------------------------------------------------------
deg <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
  select = c("gene","logFC","padj","symbol","treat_fdr","treat_p","treat_lfc"))
deg_tier1 <- deg[!is.na(treat_fdr) & treat_fdr < 0.05]
cat(sprintf("Tier 1 DEGs (TREAT FDR<0.05 at lfc=0.25): %d (expected 1,918 rows)\n",
            nrow(deg_tier1)))
if (nrow(deg_tier1) != 1918L)
  warning("TREAT Tier-1 DEG count differs from the frozen 1,918-row release")
deg_tier1[, deg_direction := fifelse(logFC > 0, "up", "down")]
setnames(deg_tier1, c("logFC","padj","treat_fdr"),
         c("deg_logFC","deg_padj","deg_treat_fdr"))
# Keep one row per non-empty symbol (most-significant by padj)
deg_tier1 <- deg_tier1[symbol != "" & !is.na(symbol)]
setorder(deg_tier1, symbol, deg_padj)
deg_tier1_short <- deg_tier1[!duplicated(symbol),
                              .(symbol, deg_logFC, deg_padj, deg_direction)]

# Definition C only has gene context — annotate is_DEG, deg_*
v[, is_DEG := definition == "C_coloc_top" &
              gene_symbol %in% deg_tier1_short$symbol]
v <- merge(v, deg_tier1_short, by.x = "gene_symbol", by.y = "symbol",
           all.x = TRUE)

# ---- GWAS-ATAC overlay -----------------------------------------------------
atac_path <- file.path(BASE,
  "GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv")
if (file.exists(atac_path)) {
  atac <- fread(atac_path)
  atac_keys <- intersect(c("variant_id","SNP","variant_key"), names(atac))
  if (length(atac_keys)) {
    atac[, variant_key := get(atac_keys[1])]
    # Some ATAC variant_ids may carry alleles ("chr:pos:a1:a2") — strip
    atac[, variant_key := sub("^([0-9]+):([0-9]+).*$", "\\1:\\2", variant_key)]
    atac_summary <- atac[, .(in_atac_peak = TRUE,
                              atac_cell_types = paste(unique(cell_type),
                                                       collapse = ";"),
                              n_atac_peaks = .N),
                          by = variant_key]
    v <- merge(v, atac_summary, by = "variant_key", all.x = TRUE)
    v[is.na(in_atac_peak), in_atac_peak := FALSE]
  } else {
    cat("ATAC file has unexpected schema; skipping overlay\n")
    v[, in_atac_peak := NA]
  }
} else {
  cat("ATAC file not found — skipping overlay\n")
  v[, in_atac_peak := NA]
}

# ---- Emit three master CSVs ------------------------------------------------
defC_master <- v[definition == "C_coloc_top"]
defAB_master <- v[definition %in% c("A_lead_gwas","B_finemapped")]
fwrite(defC_master,
       file.path(OUT, "coloc_variant_annotation.csv"))
fwrite(defAB_master,
       file.path(OUT, "lead_causal_annotation.csv"))

# CS-member annotation
cs_members[, variant_key := paste0(chromosome, ":", position)]
cs_ann <- merge(cs_members, ann, by = "variant_key", all.x = TRUE)
fwrite(cs_ann, file.path(OUT, "cs_member_annotation.csv"))

cat(sprintf("\nMaster CSVs written:\n"))
cat(sprintf("  Definition C (coloc_top):  %d rows | %d unique genes | %d coding | %d non-coding\n",
            nrow(defC_master), uniqueN(defC_master$gene_symbol),
            sum(defC_master$coarse_class == "coding"),
            sum(defC_master$coarse_class == "non-coding")))
cat(sprintf("  Definition A+B (leads):    %d rows\n", nrow(defAB_master)))
cat(sprintf("  CS members:                %d rows across %d (study, locus, CS)\n",
            nrow(cs_ann),
            if (nrow(cs_ann)) uniqueN(cs_ann[, paste(study, locus, susie_cs)]) else 0L))

# ---- Per-figure summary tables ---------------------------------------------
# F1: All COLOC genes — counts by method × class.
# Per-gene coarse-class assignment with coding-precedence: a gene called
# coding in any (GWAS, study) is coding-overall; otherwise non-coding.
# Per-gene fine-class: take the highest-precedence fine class observed.
defC_master[, method := fcase(
  !is.na(pp4_susie) & pp4_susie >= 0.5, "SuSiE-COLOC",
  !is.na(pp4_abf)   & pp4_abf   >= 0.5, "ABF-COLOC",
  default = "other")]
fine_rank <- c("coding"=1L,"spliceSite"=2L,"fiveUTR"=3L,
               "threeUTR"=4L,"promoter"=5L,"intron"=6L,"intergenic"=7L)
defC_per_gene_method <- defC_master[, {
    cls <- na.omit(fine_class)
    if (!length(cls)) data.table(fine_class = "intergenic",
                                 coarse_class = "non-coding")
    else {
      best <- cls[which.min(fine_rank[cls])]
      data.table(fine_class = best,
                 coarse_class = if (best %in% c("coding","spliceSite"))
                   "coding" else "non-coding")
    }
  }, by = .(method, gene_symbol)]
panel_F1 <- defC_per_gene_method[, .(n_genes = .N),
                                  by = .(method, fine_class, coarse_class)]
fwrite(panel_F1, file.path(OUT, "panel_F1_data.csv"))

# F2: DEG funnel
n_tier1   <- nrow(deg_tier1)
deg_coloc <- defC_master[is_DEG == TRUE, unique(gene_symbol)]
n_overlap <- length(deg_coloc)
# For coding/non-coding split: take per-gene worst-case class (coding wins
# precedence so a gene that's coding in any GWAS gets called coding overall).
deg_gene_class <- defC_master[is_DEG == TRUE, .(
    is_coding_anywhere = any(coarse_class == "coding", na.rm = TRUE)),
    by = gene_symbol]
n_deg_coding    <- deg_gene_class[, sum(is_coding_anywhere)]
n_deg_noncoding <- n_overlap - n_deg_coding

deg_coloc_split <- defC_master[is_DEG == TRUE,
  .(n_genes = uniqueN(gene_symbol)), by = .(coarse_class, fine_class)]
panel_F2 <- list(
  funnel = data.table(
    step = c("Tier-1 DEGs (TREAT FDR<0.05; lfc=0.25)",
             "DEGs ∩ COLOC (PP.H4 >= 0.5)",
             "Coding among DEG-COLOC",
             "Non-coding among DEG-COLOC"),
    count = c(n_tier1, n_overlap, n_deg_coding, n_deg_noncoding)),
  fine_split = deg_coloc_split)
fwrite(panel_F2$funnel,    file.path(OUT, "panel_F2_funnel.csv"))
fwrite(panel_F2$fine_split, file.path(OUT, "panel_F2_fine_split.csv"))

# Hypergeometric: is coding-fraction in DEG-COLOC different from coding-fraction in all COLOC?
n_all_coding    <- defC_master[coarse_class == "coding", uniqueN(gene_symbol)]
n_all_total     <- defC_master[, uniqueN(gene_symbol)]
n_deg_coding    <- defC_master[is_DEG == TRUE & coarse_class == "coding",
                               uniqueN(gene_symbol)]
n_deg_total     <- length(deg_coloc)
ht <- if (n_deg_total > 0) {
  phyper(n_deg_coding - 1, n_all_coding, n_all_total - n_all_coding,
         n_deg_total, lower.tail = FALSE)
} else NA_real_
cat(sprintf("\nF2 hypergeometric: %d/%d DEG-COLOC genes are coding (vs %d/%d overall) p=%.3g\n",
            n_deg_coding, n_deg_total, n_all_coding, n_all_total, ht))
fwrite(data.table(
  n_deg_coloc = n_deg_total, n_deg_coloc_coding = n_deg_coding,
  n_all_coloc = n_all_total, n_all_coloc_coding = n_all_coding,
  hypergeom_p = ht),
  file.path(OUT, "panel_F2_hypergeom.csv"))

# F3: Definition concordance — exact variant match only (study, variant_key).
# For each defC variant, check whether the same variant_key appears as a defA
# (GWAS lead) or defB (top finemapped) within the same study.
defA <- unique(v[definition == "A_lead_gwas",
                  .(study, variant_key, A_class = fine_class,
                    A_coarse = coarse_class)])
defB <- unique(v[definition == "B_finemapped",
                  .(study, variant_key, B_class = fine_class,
                    B_coarse = coarse_class)])
defC_for_pair <- unique(defC_master[, .(study, gene_symbol,
                                         variant_key,
                                         C_class = fine_class,
                                         C_coarse = coarse_class,
                                         pp4_best, is_DEG)])
ac <- merge(defC_for_pair, defA, by = c("study","variant_key"), all.x = TRUE)
ac[, AC_same_variant := !is.na(A_class)]
ac[, A_match_class := fifelse(AC_same_variant, A_coarse, "no_match")]
panel_F3_AC <- ac[, .(n = .N), by = .(C_coarse, A_match_class)]
fwrite(panel_F3_AC, file.path(OUT, "panel_F3_AC.csv"))

bc <- merge(defC_for_pair, defB, by = c("study","variant_key"), all.x = TRUE)
bc[, BC_same_variant := !is.na(B_class)]
bc[, B_match_class := fifelse(BC_same_variant, B_coarse, "no_match")]
panel_F3_BC <- bc[, .(n = .N), by = .(C_coarse, B_match_class)]
fwrite(panel_F3_BC, file.path(OUT, "panel_F3_BC.csv"))

# F4: Per-ancestry coding fraction × definition
panel_F4 <- v[!is.na(coarse_class),
              .(n = .N), by = .(definition, ancestry, coarse_class)]
fwrite(panel_F4, file.path(OUT, "panel_F4_ancestry.csv"))

# F5: Per-trait
panel_F5 <- v[!is.na(coarse_class),
              .(n = .N), by = .(definition, trait_cat, coarse_class)]
fwrite(panel_F5, file.path(OUT, "panel_F5_trait.csv"))

# F6: DEG-COLOC LFC × PP4
panel_F6 <- defC_master[is_DEG == TRUE,
  .(gene_symbol, study, ancestry, trait_cat, pp4_best, deg_logFC, deg_padj,
    deg_direction, coarse_class, fine_class, variant_key,
    distance_to_tss_kb, in_atac_peak)]
panel_F6 <- panel_F6[!duplicated(panel_F6, by = "gene_symbol")]
fwrite(panel_F6, file.path(OUT, "panel_F6_lfc_pp4.csv"))

# F7: Canonical examples
canonical_genes <- c("PNPLA3","TM6SF2","GCKR","HSD17B13","MTARC1","MARC1",
                     "APOE","RORA","HKDC1","SAMM50","TRIB1","SERPINA1",
                     "ADH4","CETP","EFHD1","MLIP")
panel_F7 <- defC_master[gene_symbol %in% canonical_genes,
  .(gene_symbol, study, ancestry, trait_cat, pp4_best, fine_class, coarse_class,
    deg_logFC, deg_padj, definition, variant_key)]
panel_F7 <- panel_F7[!duplicated(panel_F7, by = "gene_symbol")]
missing <- setdiff(canonical_genes, panel_F7$gene_symbol)
if (length(missing)) {
  panel_F7 <- rbind(panel_F7,
    data.table(gene_symbol = missing, study = NA, ancestry = NA,
               trait_cat = NA, pp4_best = NA_real_, fine_class = NA,
               coarse_class = NA, deg_logFC = NA_real_,
               deg_padj = NA_real_, definition = NA, variant_key = NA),
    fill = TRUE)
}
fwrite(panel_F7, file.path(OUT, "panel_F7_canonical.csv"))

# F8: CS-member coding fraction per locus
if (nrow(cs_ann)) {
  cs_ann[, is_coding := coarse_class == "coding"]
  panel_F8 <- cs_ann[, .(
    n_members = .N,
    n_coding  = sum(is_coding, na.rm = TRUE),
    pip_weighted_coding = sum(is_coding * recommended_pip, na.rm = TRUE) /
                           pmax(sum(recommended_pip, na.rm = TRUE), 1e-12),
    coloc_genes = first(coloc_genes),
    coloc_pp4_best = first(coloc_pp4_best)),
    by = .(study, locus, susie_cs)]
  panel_F8[, frac_coding_unweighted := n_coding / n_members]
  fwrite(panel_F8, file.path(OUT, "panel_F8_cs_coding_fraction.csv"))
} else {
  cat("No CS members — skipping F8 panel data\n")
}

cat("\n[03_aggregate_and_deg_overlap] DONE\n")
