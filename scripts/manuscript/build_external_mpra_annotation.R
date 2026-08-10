#!/usr/bin/env Rscript
# build_external_mpra_annotation.R
#
# Join the external Zhu / Hu et al. 2026 (Nat Genet; DOI 10.1038/s41588-026-02617-8)
# MASLD MPRA differential-activity-variant (DAV) + featured target-gene calls to our
# FROZEN expression-QTL + disease-state evidence-class table (2026-07-15-r2) and
# quantify, per gene, how the experimental variant-function layer relates to our
# genetic + disease-state classes. ADDITIVE / external annotation only -- does NOT
# mutate the frozen release and is NOT an input to the convergence/heuristic score.
# See data/external/zhu2026_masld_mpra/README.md and the number-verification firewall.
#
# Analog of build_external_scchromatin_annotation.R (Elison sc-chromatin) and
# build_external_pqtl_annotation.R (Gobeil pQTL), one modality over (experimental MPRA).
# The Zhu-side calls are DERIVED deterministically from the on-disk script-92 parse
# (GWAS/finemapping/results/seqfunc/mpra_benchmark/) -- NOT hand-typed, NOT re-downloaded.
# Our per-gene values are pulled LIVE from the frozen table (never hardcoded).
#
# Out: RNA-seq/results/multi_evidence/external_mpra/mpra_external_annotation.tsv
#      RNA-seq/results/multi_evidence/external_mpra/mpra_external_summary.txt
#      RNA-seq/results/multi_evidence/external_mpra/figS2R_source.tsv          (figure-ready)
#      data/external/zhu2026_masld_mpra/zhu2026_masld_mpra_calls.csv           (build artifact)
# Env: rnaseq

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
sub_f    <- file.path(BASE, "GWAS/finemapping/results/seqfunc/mpra_benchmark/mpra_substrate_truth.tsv")
gencode_f<- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
frozen_f <- file.path(BASE, "RNA-seq/results/manuscript_release/2026-07-15-r2/evidence_class_table.tsv")
calls_dir<- file.path(BASE, "data/external/zhu2026_masld_mpra")
out_dir  <- file.path(BASE, "RNA-seq/results/multi_evidence/external_mpra")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(sub_f), file.exists(frozen_f), file.exists(gencode_f))

num <- function(v) suppressWarnings(as.numeric(v))
DOI <- "Zhu 2026 NatGenet, DOI 10.1038/s41588-026-02617-8"

# --- (1) derive the DAV calls from the script-92 parse (firewall: no hand-typed values) ---
s <- fread(sub_f)
dav <- s[num(dav) == 1]
# gene = our finemapping/coloc locus assignment for the tested variant (made explicit in source_note)
dav[, g := fifelse(!is.na(primary_gene) & primary_gene != "" & primary_gene != "-",
                   primary_gene, gene)]
dav <- dav[!is.na(g) & g != "" & g != "-"]
dav <- dav[!grepl("^(RP11-|RP13-|RP[0-9]|CTD-|CTC-|AC[0-9]{6}|AL[0-9]{6}|LINC|MIR|SNOR)", g)]  # drop clone/ncRNA labels
# strongest allelic effect per gene
dav[, abslfc := abs(num(mpra_log2fc))]
setorder(dav, g, -abslfc)
dav1 <- dav[, .SD[1], by = g]
dav_calls <- dav1[, .(
  gene            = g,
  zhu_locus       = rsid,
  zhu_finding_type= "mpra_dav",
  zhu_layer       = "MPRA",
  zhu_celltype    = cell_model,
  zhu_lead_variant= rsid,
  zhu_direction   = fifelse(num(mpra_log2fc) > 0, "activity_up", "activity_down"),
  zhu_mpra_log2fc = round(num(mpra_log2fc), 3),
  zhu_context     = paste0(cell_model, "_", stimulus),
  zhu_detail      = sprintf("significant MPRA DAV (log2FC=%.2f, FDR=%.1e) in %s_%s at the variant our map assigns to this locus",
                            num(mpra_log2fc), num(mpra_fdr), cell_model, stimulus),
  source_note     = sprintf("GSE281364 MPRA DAV %s (gene=our coloc/finemap locus label); via 92_hu2025_mpra_benchmark.R substrate parse; %s",
                            rsid, DOI))]

# --- (2) featured main-text targets (only hand-entered symbols; evidence pulled live below) ---
feat_calls <- data.table(
  gene            = c("SLC22A3", "APOA5", "ANGPTL3", "LPL"),
  zhu_locus       = c("6q25", "11q23", "1p31", "8p21"),
  zhu_finding_type= "featured_target",
  zhu_layer       = "MPRA+eQTL+looping+CRISPRi",
  zhu_celltype    = c("HepG2", "HepG2", "HepG2", "HepG2/LX2"),
  zhu_lead_variant= NA_character_,
  zhu_direction   = NA_character_,
  zhu_mpra_log2fc = NA_real_,
  zhu_context     = "featured",
  zhu_detail      = c(
    "featured DAV target; modulates SLC22A3 expression",
    "featured triglyceride-metabolism DAV target",
    "featured triglyceride-metabolism DAV target",
    "featured triglyceride-metabolism / HSC-activation DAV target"),
  source_note     = paste0("main-text featured target; ", DOI))

calls <- rbindlist(list(feat_calls, dav_calls[!gene %in% feat_calls$gene]), use.names = TRUE)

# protein-coding restriction (drop residual pseudogene/ncRNA labels) using gencode v49
gc <- fread(gencode_f)
pc <- unique(gc[gene_biotype == "protein_coding", gene_name])
calls <- calls[gene %in% pc | zhu_finding_type == "featured_target"]

fwrite(calls, file.path(calls_dir, "zhu2026_masld_mpra_calls.csv"))

# --- (3) pull our frozen per-gene expression + disease-state evidence (LIVE) ----------------
e <- fread(frozen_f)
keep <- c("symbol", "primary_evidence_class", "sensitivity_evidence_class",
          "genetic_trait_scope", "bulk_logFC", "bulk_treat_fdr", "in_treat_deg",
          "max_susie_pp4_all", "max_abf_pp4_all",
          "max_susie_pp4_direct", "max_abf_pp4_direct",
          "max_susie_pp4_enzyme", "max_abf_pp4_enzyme")
e <- e[, ..keep]
setnames(e, keep,
         c("symbol", "our_primary_class", "our_sensitivity_class",
           "our_genetic_trait_scope", "our_bulk_logFC", "our_bulk_treat_fdr",
           "our_in_treat_deg", "our_susie_all", "our_abf_all",
           "our_susie_direct", "our_abf_direct",
           "our_susie_enzyme", "our_abf_enzyme"))

m <- merge(calls, e, by.x = "gene", by.y = "symbol", all.x = TRUE, sort = FALSE)

GATE <- 0.5
sA <- fcoalesce(num(m$our_susie_all), 0); aA <- fcoalesce(num(m$our_abf_all), 0)
sD <- fcoalesce(num(m$our_susie_direct), 0); aD <- fcoalesce(num(m$our_abf_direct), 0)
m[, our_genetic_any    := (sA > GATE) | (aA > GATE)]
m[, our_genetic_direct := (sD > GATE) | (aD > GATE)]
m[, our_is_deg         := our_in_treat_deg == TRUE]

# --- (4) relationship call (their DAV/featured + our frozen values) -------------------------
classify <- function(r) {
  if (isTRUE(r$our_genetic_direct))                    return("corroborated_genetic_direct")   # direct-MASLD coloc >0.5
  if (isTRUE(r$our_is_deg))                            return("corroborated_disease_state")    # canonical DEG
  if (isTRUE(r$our_genetic_any))                       return("corroborated_genetic_enzyme")   # enzyme/indirect coloc >0.5
  return("expression_miss")                                                                    # Zhu functional, our maps null (gap)
}
m[, relationship := vapply(seq_len(.N), function(i) classify(m[i]), character(1))]

# honest caveats: enzyme-trait anchoring, DEG direction, coding-mechanism blind spots
coding <- c("PNPLA3","TM6SF2","HSD17B13","GCKR")
m[, caveat := fifelse(our_genetic_any == TRUE & our_genetic_direct == FALSE & our_genetic_trait_scope == "enzyme",
                      "coloc is enzyme/liver-injury-trait anchored; direct-MASLD sub-gate",
              fifelse(gene %in% coding & relationship == "expression_miss",
                      "coding/post-translational MASLD gene: blind to expression-COLOC (MPRA still flags the cis element)",
              fifelse(relationship == "expression_miss",
                      "MPRA-functional but our expression-genetics null: experimental layer adds signal we miss",
                      "")))]

ordv <- c("gene", "zhu_locus", "zhu_finding_type", "zhu_layer", "zhu_celltype",
          "zhu_lead_variant", "zhu_direction", "zhu_mpra_log2fc", "zhu_context",
          "our_primary_class", "our_sensitivity_class", "our_genetic_trait_scope",
          "our_susie_all", "our_abf_all", "our_susie_direct", "our_abf_direct",
          "our_susie_enzyme", "our_abf_enzyme", "our_is_deg", "our_bulk_logFC",
          "our_bulk_treat_fdr", "our_genetic_any", "our_genetic_direct",
          "relationship", "caveat", "zhu_detail", "source_note")
m <- m[, intersect(ordv, names(m)), with = FALSE]
out_tsv <- file.path(out_dir, "mpra_external_annotation.tsv")
fwrite(m, out_tsv, sep = "\t")

# figure-ready long table (FigS2R)
fig <- m[, .(gene, zhu_finding_type, zhu_celltype, zhu_direction, zhu_mpra_log2fc,
             relationship, our_primary_class, our_susie_direct, our_abf_direct,
             our_susie_enzyme, our_abf_enzyme, our_susie_all, our_abf_all,
             our_is_deg, our_bulk_treat_fdr, our_bulk_logFC, caveat)]
fwrite(fig, file.path(out_dir, "figS2R_source.tsv"), sep = "\t")

# --- (5) headline counts -------------------------------------------------------------------
sink(file.path(out_dir, "mpra_external_summary.txt"))
cat("External experimental MPRA/CRISPRi (Zhu/Hu et al. 2026 Nat Genet) vs our frozen maps\n")
cat("====================================================================================\n\n")
cat("Rows:", nrow(m), " genes  (", sum(m$zhu_finding_type=="featured_target"),
    "featured +", sum(m$zhu_finding_type=="mpra_dav"), "MPRA-DAV )\n\n")
cat("Per-gene relationship to our expression-genetic + disease-state maps:\n")
print(m[order(zhu_finding_type, relationship),
        .(gene, zhu_finding_type,
          susieDir = round(fcoalesce(num(our_susie_direct),0),3),
          abfAll   = round(fcoalesce(num(our_abf_all),0),3),
          deg = our_is_deg, scope = our_genetic_trait_scope, relationship)])
cat("\nHeadline counts:\n"); print(m[, .N, by = relationship][order(-N)])
cat("\nFeatured-target verdicts:\n")
print(m[zhu_finding_type=="featured_target", .(gene, relationship, caveat)])
cat("\nCorroborated at our GENETIC map (direct-MASLD coloc > 0.5):\n  ",
    paste(m[our_genetic_direct == TRUE]$gene, collapse = ", "), "\n", sep = "")
cat("Corroborated at our DISEASE-STATE map (canonical DEG):\n  ",
    paste(m[our_is_deg == TRUE]$gene, collapse = ", "), "\n", sep = "")
cat("Expression-miss gaps (Zhu MPRA-functional, our expression-genetics null):\n  ",
    paste(m[relationship == "expression_miss"]$gene, collapse = ", "), "\n", sep = "")
sink()

cat("Wrote:", out_tsv, "\n")
cat(readLines(file.path(out_dir, "mpra_external_summary.txt")), sep = "\n")
