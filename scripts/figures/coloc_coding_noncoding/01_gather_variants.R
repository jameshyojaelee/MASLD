#!/usr/bin/env Rscript
# 01_gather_variants.R
# ----------------------------------------------------------------------------
# Build a long-format variant table for the coding/non-coding pipeline.
# Three operational definitions of the "lead variant":
#   A. Lead causal (GWAS lead)        — pre-finemapping; data/lead_snps/*.tsv
#   B. Top finemapped causal          — credible_sets_1kg.csv, max recommended_pip
#   C. Top gene-level COLOC           — promoted aggregate COLOC release
#
# Outputs (RNA-seq/results/coloc_variant_classes/):
#   variants_long.csv      — stacked rows from A, B, C with variant_key
#   cs_members_long.csv    — every CS member at COLOC-positive (study, locus)
suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT  <- Sys.getenv(
  "FIG2_VARIANT_CLASS_DIR",
  file.path(BASE, "RNA-seq/results/coloc_variant_classes")
)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

POLYFUN_DIR <- file.path(BASE, "GWAS/finemapping/results/susie_coloc_polyfun")
KG_DIR      <- file.path(BASE, "GWAS/finemapping/results/susie_coloc_1kg")
OTHER_DIR   <- file.path(BASE, "GWAS/finemapping/results/susie_coloc")
COLOC_FILE  <- Sys.getenv(
  "FIG2_COLOC_INPUT",
  file.path(OTHER_DIR, "susie_coloc_all_gwas.csv")
)
LEAD_DIR    <- file.path(BASE, "GWAS/finemapping/data/lead_snps")
CS_FILE     <- file.path(BASE, "GWAS/finemapping/results/credible_sets.csv")

# Active studies are defined by the promoted aggregate, not by leftover
# per-study directories from older releases.
stopifnot(file.exists(COLOC_FILE))
active_studies <- unique(fread(COLOC_FILE, select = "gwas_name")$gwas_name)

# ---- Canonical study → method/ancestry/trait map ---------------------------
# LD-panel priority (per Phase 10 swap, 2026-05-06): PolyFun > 1kg > other.
# Any study with PolyFun output is routed to PolyFun automatically below.


trait_of <- function(s) {
  # Tokens are NOT end-anchored: MVP strata carry the trait mid-name
  # (MVP_ALT_AFR / MVP_AST_EAS / MVP_Albumin_EUR), so the prior "ALT$/AST$/GGT$"
  # anchors dropped every MVP quant-trait stratum into "Other". Albumin/Platelet/
  # ChronLiver tokens were also absent (MVP-only traits). ALT/AST/GGT do not
  # substring-collide with any other registry study_name.
  fcase(
    grepl("ALT",        s), "Liver enzymes",
    grepl("AST",        s), "Liver enzymes",
    grepl("GGT",        s), "Liver enzymes",
    grepl("Albumin",    s, ignore.case = TRUE), "Albumin",
    grepl("Platelet",   s, ignore.case = TRUE), "Platelet",
    grepl("PDFF",       s), "PDFF",
    grepl("ChronLiver", s, ignore.case = TRUE), "Chronic liver disease",
    grepl("Cirrhosis|CHIRHEP", s, ignore.case = TRUE), "Cirrhosis",
    grepl("HCC",        s), "HCC",
    grepl("NAFLD|NASH", s, ignore.case = TRUE), "NAFLD/NASH",
    default = "Other")
}
# Registry-driven ancestry (canonical source). The prior regex heuristic had no
# AMR branch, so MVP *_AMR strata fell through to the EUR default. The finemapping
# registry carries an explicit per-study ancestry column (EUR/AFR/AMR/EAS/SAS);
# use it as the authority, keeping the regex only as a defensive fallback for any
# study not present in the registry (e.g. lead-SNP files without a registry row).
REGISTRY <- file.path(BASE, "GWAS/finemapping/config/gwas_registry.tsv")
reg_ancestry <- NULL
if (file.exists(REGISTRY)) {
  reg <- fread(REGISTRY)
  if (all(c("study_name", "ancestry") %in% names(reg)))
    reg_ancestry <- setNames(as.character(reg$ancestry), reg$study_name)
}
ancestry_of <- function(s) {
  fallback <- fcase(
    grepl("EAS$|^BBJ_", s), "EAS",
    grepl("AMR$",       s), "AMR",
    grepl("AFR",        s), "AFR",
    grepl("CSA",        s), "SAS",
    grepl("EUR$|^UKBB|^FinnGen|^Ghouse|^2019|^2020|^2021|^2022|^2023", s), "EUR",
    default = "EUR")
  if (is.null(reg_ancestry)) return(fallback)
  reg_hit <- unname(reg_ancestry[s])
  fifelse(is.na(reg_hit), fallback, reg_hit)
}

# ---- 1. Definition C: top gene-level COLOC variant per (gene, GWAS) ---------
gather_definition_C <- function() {
  d <- fread(COLOC_FILE)
  required <- c("gene", "ensembl", "chr", "gwas_name", "PP.H4.abf",
                "PP.H4.susie", "top_snp", "top_snp_PP")
  stopifnot(all(required %in% names(d)))
  d[, `:=`(
    PP.H4.susie = suppressWarnings(as.numeric(PP.H4.susie)),
    PP.H4.abf = suppressWarnings(as.numeric(PP.H4.abf))
  )]
  d[, pp4_best := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
  d[!is.finite(pp4_best), pp4_best := NA_real_]
  d <- d[is.finite(pp4_best) & pp4_best > 0.5 &
         !is.na(gene) & gene != "" & !is.na(top_snp) & top_snp != ""]
  setorder(d, gene, gwas_name, -pp4_best, -top_snp_PP)
  out <- d[!duplicated(d, by = c("gene", "gwas_name"))]
  if (!nrow(out)) return(data.table())
  out[, study := gwas_name]
  out[, ancestry  := ancestry_of(study)]
  out[, trait_cat := trait_of(study)]
  out[, definition := "C_coloc_top"]
  out[, c("chr_p","pos_p") := tstrsplit(top_snp, ":", fixed = TRUE)]
  out[, chr := suppressWarnings(as.integer(chr_p))]
  out[, pos_hg19 := suppressWarnings(as.integer(pos_p))]
  out[, locus := paste0(chr, ".", pos_hg19)]    # filled by GWAS lead later if needed
  out[, variant_key := top_snp]
  out[, c("chr_p","pos_p") := NULL]
  setnames(out, c("PP.H4.susie","PP.H4.abf","gene","ensembl"),
                 c("pp4_susie","pp4_abf","gene_symbol","ensembl"))
  cat(sprintf("Definition C: %d (gene, study) rows; %d unique genes\n",
              nrow(out), uniqueN(out$gene_symbol)))
  out[, .(definition, study, ancestry, trait_cat,
          chr, pos_hg19, locus,
          gene_symbol, ensembl,
          pp4_susie, pp4_abf, pp4_best,
          top_snp_PP,
          recommended_pip = NA_real_, n_cs_members = NA_integer_,
          variant_key)]
}

# ---- 2. Definition B: top finemapped causal per converged credible set ------
gather_definition_B <- function() {
  cs <- fread(CS_FILE, select = c("chromosome","position","allele1","allele2",
                                   "locus","study","ancestry",
                                   "susie_pip","susie_cs","susie_converged",
                                   "susie_reliable","recommended_pip",
                                   "either_in_cs"))
  cs <- cs[!is.na(susie_cs) & susie_cs > 0 &
           susie_converged == TRUE & susie_reliable == TRUE]
  if (!nrow(cs)) return(data.table())
  setorder(cs, study, locus, susie_cs, -recommended_pip)
  top <- cs[!duplicated(cs, by = c("study","locus","susie_cs"))]
  cs_size <- cs[, .(n_cs_members = .N), by = .(study, locus, susie_cs)]
  top <- merge(top, cs_size, by = c("study","locus","susie_cs"), all.x = TRUE)
  top[, definition := "B_finemapped"]
  top[, trait_cat  := trait_of(study)]
  top[, ancestry   := ifelse(is.na(ancestry) | ancestry == "",
                              ancestry_of(study), ancestry)]
  top[, variant_key := paste0(chromosome, ":", position)]
  top[, .(definition, study, ancestry, trait_cat,
          chr = chromosome, pos_hg19 = position, locus,
          gene_symbol = NA_character_, ensembl = NA_character_,
          pp4_susie = NA_real_, pp4_abf = NA_real_, pp4_best = NA_real_,
          top_snp_PP = NA_real_,
          recommended_pip, n_cs_members,
          variant_key)]
}

# ---- 3. Definition A: GWAS lead per study × locus --------------------------
gather_definition_A <- function() {
  files <- list.files(LEAD_DIR, pattern = "_leadSNPs.*\\.tsv$",
                      full.names = TRUE)
  # Prefer "_merged" version when both present (matches finemapping pipeline)
  files_dt <- data.table(path = files, base = basename(files))
  files_dt[, study := sub("(_(infr_)?preprocessed)?_leadSNPs.*", "", base)]
  # Filter for active studies only (prevents dropped studies from being loaded)
  files_dt <- files_dt[study %in% active_studies]
  files_dt[, has_merged := grepl("_merged\\.tsv$", base)]
  setorder(files_dt, study, -has_merged)
  files_dt <- files_dt[!duplicated(study)]
  rows <- list()
  for (i in seq_len(nrow(files_dt))) {
    d <- tryCatch(fread(files_dt$path[i]), error = function(e) NULL)
    if (is.null(d) || !nrow(d)) next
    chr_col <- intersect(c("chr","CHR","chromosome"), names(d))[1]
    pos_col <- intersect(c("pos","BP","position","POS","bp"), names(d))[1]
    if (is.na(chr_col) || is.na(pos_col)) next
    d[, study := files_dt$study[i]]
    d[, chr  := as.integer(sub("^chr","", get(chr_col)))]
    d[, pos_hg19 := as.integer(get(pos_col))]
    d <- d[!is.na(chr) & !is.na(pos_hg19)]
    rows[[files_dt$study[i]]] <- d[, .(study, chr, pos_hg19,
                                       locus = paste0(chr, ".", pos_hg19))]
  }
  out <- rbindlist(rows, fill = TRUE)
  if (!nrow(out)) return(data.table())
  out[, definition := "A_lead_gwas"]
  out[, ancestry   := ancestry_of(study)]
  out[, trait_cat  := trait_of(study)]
  out[, variant_key := paste0(chr, ":", pos_hg19)]
  out[, .(definition, study, ancestry, trait_cat,
          chr, pos_hg19, locus,
          gene_symbol = NA_character_, ensembl = NA_character_,
          pp4_susie = NA_real_, pp4_abf = NA_real_, pp4_best = NA_real_,
          top_snp_PP = NA_real_,
          recommended_pip = NA_real_, n_cs_members = NA_integer_,
          variant_key)]
}

# ---- 4. CS members at COLOC-positive (study, locus) ------------------------
# For each COLOC-positive gene, locate the CS that contains its top SNP, then
# pull all variants in that CS. This gives the actual credible set behind the
# colocalization signal — the uncertainty layer for downstream F8.
gather_cs_members <- function(defC) {
  if (!nrow(defC)) return(data.table())
  cs <- fread(CS_FILE)
  cs[, variant_key := paste0(chromosome, ":", position)]
  # Match COLOC top SNPs to CS rows on (study, variant_key)
  defC_keys <- unique(defC[, .(study, variant_key,
                                gene_symbol, ensembl, pp4_best)])
  hit <- merge(cs, defC_keys, by = c("study","variant_key"))
  if (!nrow(hit)) return(data.table())
  # Each hit reveals the (study, locus, susie_cs) of the COLOC-driving CS
  cs_keys <- unique(hit[, .(study, locus, susie_cs)])
  # Now pull all CS members at those (study, locus, susie_cs)
  members <- merge(cs, cs_keys, by = c("study","locus","susie_cs"))
  # Attach the colocalized gene(s) at this CS
  gene_lookup <- hit[, .(coloc_genes = paste(unique(gene_symbol), collapse = ";"),
                          coloc_pp4_best = max(pp4_best, na.rm = TRUE)),
                      by = .(study, locus, susie_cs)]
  members <- merge(members, gene_lookup,
                   by = c("study","locus","susie_cs"), all.x = TRUE)
  members
}

# ---- Run ---------------------------------------------------------------------
cat("Gathering Definition C (top gene-level COLOC)...\n")
defC <- gather_definition_C()
cat("Gathering Definition B (top finemapped causal)...\n")
defB <- gather_definition_B()
cat("Gathering Definition A (GWAS lead)...\n")
defA <- gather_definition_A()

variants_long <- rbindlist(list(defA, defB, defC), use.names = TRUE, fill = TRUE)
fwrite(variants_long, file.path(OUT, "variants_long.csv"))
cat(sprintf("\nvariants_long: %d rows | A=%d, B=%d, C=%d | %d unique variant_keys\n",
            nrow(variants_long), nrow(defA), nrow(defB), nrow(defC),
            uniqueN(variants_long$variant_key)))

# CS members at the credible sets that contain each COLOC top SNP
cs_members <- gather_cs_members(defC)
fwrite(cs_members, file.path(OUT, "cs_members_long.csv"))
cat(sprintf("cs_members_long: %d rows across %d (study, locus, CS) triples\n",
            nrow(cs_members),
            if (nrow(cs_members)) uniqueN(cs_members[, paste(study, locus, susie_cs)]) else 0L))

cat("\n[01_gather_variants] DONE\n")
