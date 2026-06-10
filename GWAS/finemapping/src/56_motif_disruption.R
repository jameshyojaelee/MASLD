#!/usr/bin/env Rscript
# 56_motif_disruption.R
# Predicts TF motif disruption for GWAS variants in scATAC peaks
# Uses motifbreakR (primary) + FIMO (validation)
# Requires: 55_gwas_atac_variant_overlap.R output
# Usage: Rscript 56_motif_disruption.R

library(data.table)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifbreakR)
library(TFBSTools)
library(MotifDb)

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE_DIR, "GWAS/finemapping")
ATAC_DIR <- file.path(BASE_DIR, "Analysis/ATAC/Human_Multiome")
OUT_DIR  <- file.path(FM_DIR, "results/gwas_atac")

cat("============================================================\n")
cat("56_motif_disruption.R\n")
cat("TF motif disruption analysis for GWAS-ATAC variants\n")
cat("============================================================\n\n")

# ── 1. Load peak-overlapping variants ────────────────────────────────────────
ann_file <- file.path(OUT_DIR, "gwas_atac_variant_annotation.csv")
if (!file.exists(ann_file)) {
  stop("Variant annotation not found. Run 55_gwas_atac_variant_overlap.R first.")
}
ann <- fread(ann_file)
cat("Loaded", nrow(ann), "variant-peak annotations\n")

# Deduplicate to unique variants (same variant can appear in multiple cell types)
variants <- ann %>%
  distinct(chromosome, position, allele1, allele2, chr_hg38, pos_hg38, max_pip) %>%
  arrange(desc(max_pip))
cat("Unique variants in peaks:", nrow(variants), "\n")

# ── 2. Load disease master regulators for cross-referencing ──────────────────
# Cross-modality definition (Analysis/ATAC/Human_Multiome/scripts/04c_*):
# hepatocyte SCENIC+ regulon TFs that are well-powered MASLD disease signals
# (bulk DEG at n=846 OR genetic COLOC PP.H4>0.5). This is INDEPENDENT of the
# motif-disruption test below, so the cross-reference is non-circular.
# NOTE: the legacy FDR-gated `disease_regulons.csv` is EMPTY — the donor-level
# (n=18) regulon-activity DE returns 0 at FDR; the old "12/24 TFs" came from a
# cell-level (pseudoreplicated) test and is deliberately NOT used. Override the
# source with REGULON_FILE if needed.
regulon_file <- Sys.getenv("REGULON_FILE",
  unset = file.path(ATAC_DIR, "scenic_plus/disease_master_regulators.csv"))
if (file.exists(regulon_file)) {
  regulons <- fread(regulon_file)
  disease_tfs <- toupper(regulons$tf_name)
  cat("Loaded", length(disease_tfs), "SCENIC+ disease master-regulator TFs from",
      basename(regulon_file), "\n")
} else {
  cat("WARNING:", basename(regulon_file), "not found — skipping regulon cross-reference\n")
  disease_tfs <- character(0)
  regulons <- data.table()
}

# ── 3. Prepare variant input for motifbreakR ─────────────────────────────────
cat("\n--- Preparing variants for motifbreakR ---\n")

# motifbreakR accepts GRanges with REF and ALT columns
# Our variants have allele1/allele2 but we need to determine which is REF
# Look up the reference allele from BSgenome
genome <- BSgenome.Hsapiens.UCSC.hg38

# Create GRanges
gr_snps <- GRanges(
  seqnames = variants$chr_hg38,
  ranges   = IRanges(start = variants$pos_hg38, width = 1),
  strand   = "*"
)

# Get reference allele from genome
ref_alleles <- as.character(getSeq(genome, gr_snps))
cat("Retrieved reference alleles for", length(ref_alleles), "variants\n")

# Assign REF/ALT based on genome reference
variants$ref_genome <- ref_alleles
variants$REF <- ifelse(variants$allele1 == ref_alleles, variants$allele1,
                 ifelse(variants$allele2 == ref_alleles, variants$allele2, NA_character_))
variants$ALT <- ifelse(variants$allele1 == ref_alleles, variants$allele2,
                 ifelse(variants$allele2 == ref_alleles, variants$allele1, NA_character_))

# Filter out variants where neither allele matches reference (could be strand issues)
n_mismatch <- sum(is.na(variants$REF))
if (n_mismatch > 0) {
  cat("WARNING:", n_mismatch, "variants where neither allele matches hg38 reference — skipping these\n")
  variants <- variants %>% filter(!is.na(REF))
}

# Also filter to SNVs only (motifbreakR handles SNVs, not indels)
variants <- variants %>%
  filter(nchar(REF) == 1 & nchar(ALT) == 1)
cat("SNV variants for motif analysis:", nrow(variants), "\n")

if (nrow(variants) == 0) {
  cat("WARNING: No SNVs remain for motif analysis. Writing empty results.\n")
  fwrite(data.table(), file.path(OUT_DIR, "motif_disruption_scores.csv"))
  quit(save = "no", status = 0)
}

# Create motifbreakR-compatible input
# motifbreakR needs a GRanges with SNP_id, REF, and ALT in mcols
gr_mb <- GRanges(
  seqnames = variants$chr_hg38,
  ranges   = IRanges(start = variants$pos_hg38, width = 1),
  strand   = "*"
)
mcols(gr_mb)$SNP_id <- paste(variants$chromosome, variants$position,
                              variants$allele1, variants$allele2, sep = ":")
mcols(gr_mb)$REF <- DNAStringSet(variants$REF)
mcols(gr_mb)$ALT <- DNAStringSet(variants$ALT)
names(gr_mb) <- mcols(gr_mb)$SNP_id

# Critical: motifbreakR reads genome from this attribute (not from genome())
attr(gr_mb, "genome.package") <- "BSgenome.Hsapiens.UCSC.hg38"

cat("Prepared", length(gr_mb), "SNVs for motifbreakR\n")

# ── 4. Run motifbreakR ──────────────────────────────────────────────────────
cat("\n--- Running motifbreakR ---\n")

# Use JASPAR vertebrate core motifs via MotifDb
# Filter to human TF motifs
jaspar_motifs <- query(MotifDb, andStrings = c("hsapiens"), orStrings = c("jaspar2022", "jaspar2024"))
if (length(jaspar_motifs) == 0) {
  # Fallback: use all JASPAR motifs
  jaspar_motifs <- query(MotifDb, "jaspar")
  jaspar_motifs <- jaspar_motifs[grep("Hsapiens", names(jaspar_motifs))]
}

if (length(jaspar_motifs) == 0) {
  # Broader fallback: use HOCOMOCO
  jaspar_motifs <- query(MotifDb, "HOCOMOCOv11-core-A")
  if (length(jaspar_motifs) == 0) {
    jaspar_motifs <- query(MotifDb, "HOCOMOCO")
  }
}

cat("Using", length(jaspar_motifs), "TF motifs from MotifDb\n")

# Run motifbreakR — this is the main computation
tryCatch({
  mb_results <- motifbreakR(
    snpList    = gr_mb,
    filterp    = TRUE,
    pwmList    = jaspar_motifs,
    threshold  = 1e-4,
    method     = "ic",
    bkg        = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    BPPARAM    = BiocParallel::MulticoreParam(workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
  )
  cat("motifbreakR returned", length(mb_results), "variant-motif results\n")
}, error = function(e) {
  cat("ERROR in motifbreakR:", conditionMessage(e), "\n")
  cat("Attempting without parallel...\n")
  mb_results <<- motifbreakR(
    snpList    = gr_mb,
    filterp    = TRUE,
    pwmList    = jaspar_motifs,
    threshold  = 1e-4,
    method     = "ic",
    bkg        = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    BPPARAM    = BiocParallel::SerialParam()
  )
  cat("motifbreakR (serial) returned", length(mb_results), "variant-motif results\n")
})

if (length(mb_results) == 0) {
  cat("WARNING: motifbreakR returned no results. Writing empty output.\n")
  fwrite(data.table(), file.path(OUT_DIR, "motif_disruption_scores.csv"))
  quit(save = "no", status = 0)
}

# Convert to data.table
mb_dt <- as.data.table(mb_results)
cat("Raw motifbreakR hits:", nrow(mb_dt), "\n")

# Filter to strong effects
mb_dt <- mb_dt %>%
  filter(effect == "strong" | abs(alleleDiff) > 0.3)
cat("Strong effect / alleleDiff > 0.3:", nrow(mb_dt), "variant-motif pairs\n")

# Extract TF name from motif providerName/geneSymbol
if ("geneSymbol" %in% colnames(mb_dt)) {
  mb_dt$tf_name <- mb_dt$geneSymbol
} else if ("providerName" %in% colnames(mb_dt)) {
  mb_dt$tf_name <- mb_dt$providerName
} else {
  mb_dt$tf_name <- mb_dt$providerId
}

# Cross-reference with SCENIC+ disease regulons
mb_dt$motif_in_disease_regulon <- toupper(mb_dt$tf_name) %in% disease_tfs

# Add regulon activity diff for matching TFs
if (nrow(regulons) > 0) {
  regulon_lookup <- regulons %>%
    mutate(tf_upper = toupper(tf_name)) %>%
    select(tf_upper, regulon_activity_diff, activity_padj) %>%
    distinct(tf_upper, .keep_all = TRUE)

  mb_dt <- mb_dt %>%
    mutate(tf_upper = toupper(tf_name)) %>%
    left_join(regulon_lookup, by = "tf_upper") %>%
    select(-tf_upper)
}

n_regulon <- sum(mb_dt$motif_in_disease_regulon)
cat("  Disrupting SCENIC+ disease regulon motifs:", n_regulon, "\n")

# ── 5. FIMO validation ──────────────────────────────────────────────────────
cat("\n--- FIMO validation ---\n")

# Only validate strong motifbreakR hits
strong_hits <- mb_dt %>% filter(effect == "strong")
cat("Running FIMO validation on", nrow(strong_hits), "strong disruptions\n")

# Check if FIMO is available
fimo_bin <- Sys.which("fimo")
if (fimo_bin == "") {
  # Hardcoded fallback path (module load in sbatch sets PATH, but just in case)
  fimo_bin <- "/nfs/sw/easybuild/software/meme/5.5.9/bin/fimo"
}

fimo_available <- file.exists(fimo_bin)
if (!fimo_available) {
  cat("WARNING: FIMO not found. Ensure 'module load meme/5.5.9' was run before launching R.\n")
}

if (fimo_available && nrow(strong_hits) > 0) {
  # Create temp directory for FIMO
  fimo_dir <- file.path(OUT_DIR, "fimo_tmp")
  dir.create(fimo_dir, showWarnings = FALSE)

  # Get unique variant positions from strong hits
  fimo_variants <- strong_hits %>%
    distinct(SNP_id, seqnames, start, end) %>%
    head(500)  # Limit to top 500 to keep compute reasonable

  # Extract ±25bp flanking sequences
  flank_size <- 25L
  gr_fimo <- GRanges(
    seqnames = fimo_variants$seqnames,
    ranges   = IRanges(
      start = as.integer(fimo_variants$start) - flank_size,
      end   = as.integer(fimo_variants$end) + flank_size
    )
  )

  # Get REF sequences
  ref_seqs <- getSeq(genome, gr_fimo)
  names(ref_seqs) <- fimo_variants$SNP_id

  # Create ALT sequences (substitute at position flank_size + 1)
  alt_seqs <- ref_seqs
  # Need to get ALT alleles
  var_lookup <- variants %>%
    mutate(SNP_id = paste(chromosome, position, allele1, allele2, sep = ":")) %>%
    select(SNP_id, ALT)
  fimo_variants <- fimo_variants %>% left_join(var_lookup, by = "SNP_id")

  for (i in seq_len(nrow(fimo_variants))) {
    alt_allele <- fimo_variants$ALT[i]
    if (!is.na(alt_allele) && nchar(alt_allele) == 1) {
      pos_in_seq <- flank_size + 1L
      alt_seq <- as.character(alt_seqs[[i]])
      substr(alt_seq, pos_in_seq, pos_in_seq) <- alt_allele
      alt_seqs[[i]] <- DNAString(alt_seq)
    }
  }
  alt_seqs <- DNAStringSet(alt_seqs)
  names(alt_seqs) <- paste0(fimo_variants$SNP_id, "_alt")

  # Write FASTA files
  ref_fa <- file.path(fimo_dir, "ref_sequences.fa")
  alt_fa <- file.path(fimo_dir, "alt_sequences.fa")
  writeXStringSet(ref_seqs, ref_fa)
  writeXStringSet(alt_seqs, alt_fa)

  # Find JASPAR motif database file
  # motifbreakR uses MotifDb PWMs, for FIMO we need MEME format
  # Export PWMs to MEME format
  meme_file <- file.path(fimo_dir, "motifs.meme")

  # Get the PWMs used by motifbreakR for the strong hits
  hit_motif_ids <- unique(strong_hits$providerId)
  hit_motifs <- jaspar_motifs[names(jaspar_motifs) %in% hit_motif_ids]
  if (length(hit_motifs) == 0) {
    # Try matching on partial names
    hit_motifs <- jaspar_motifs[grep(paste(head(hit_motif_ids, 10), collapse = "|"), names(jaspar_motifs))]
  }

  if (length(hit_motifs) > 0) {
    # Write MEME format manually
    sink(meme_file)
    cat("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n")
    cat("Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")

    for (m_name in names(hit_motifs)) {
      pwm <- as.matrix(hit_motifs[[m_name]])
      cat("MOTIF", m_name, "\n")
      cat("letter-probability matrix: alength= 4 w=", nrow(pwm), "\n")
      write.table(pwm, row.names = FALSE, col.names = FALSE, sep = "\t")
      cat("\n")
    }
    sink()
    cat("Exported", length(hit_motifs), "motifs to MEME format\n")

    # Run FIMO on ref and alt
    ref_out <- file.path(fimo_dir, "fimo_ref.tsv")
    alt_out <- file.path(fimo_dir, "fimo_alt.tsv")

    ref_cmd <- sprintf("%s --thresh 1e-4 --text --no-qvalue %s %s > %s 2>/dev/null",
                       fimo_bin, meme_file, ref_fa, ref_out)
    alt_cmd <- sprintf("%s --thresh 1e-4 --text --no-qvalue %s %s > %s 2>/dev/null",
                       fimo_bin, meme_file, alt_fa, alt_out)

    system(ref_cmd)
    system(alt_cmd)

    # Parse FIMO results
    if (file.exists(ref_out) && file.info(ref_out)$size > 0) {
      fimo_ref <- tryCatch(fread(ref_out, skip = 1), error = function(e) data.table())
      fimo_alt <- tryCatch(fread(alt_out, skip = 1), error = function(e) data.table())

      if (nrow(fimo_ref) > 0 || nrow(fimo_alt) > 0) {
        # Compare: which motifs are found in ref but not alt (disrupted) or alt but not ref (created)
        if (nrow(fimo_ref) > 0) {
          colnames(fimo_ref)[1:3] <- c("motif_id", "motif_alt_id", "sequence_name")
          fimo_ref$source <- "ref"
        }
        if (nrow(fimo_alt) > 0) {
          colnames(fimo_alt)[1:3] <- c("motif_id", "motif_alt_id", "sequence_name")
          fimo_alt$source <- "alt"
          # Strip _alt suffix from sequence names
          fimo_alt$sequence_name <- gsub("_alt$", "", fimo_alt$sequence_name)
        }

        # Count motif hits per variant for ref vs alt
        ref_counts <- fimo_ref %>%
          group_by(sequence_name, motif_id) %>%
          summarise(ref_n_hits = n(), ref_min_pval = min(as.numeric(V8), na.rm = TRUE), .groups = "drop")
        alt_counts <- fimo_alt %>%
          group_by(sequence_name, motif_id) %>%
          summarise(alt_n_hits = n(), alt_min_pval = min(as.numeric(V8), na.rm = TRUE), .groups = "drop")

        fimo_compare <- full_join(ref_counts, alt_counts,
                                   by = c("sequence_name", "motif_id")) %>%
          mutate(
            ref_n_hits = replace_na(ref_n_hits, 0L),
            alt_n_hits = replace_na(alt_n_hits, 0L),
            fimo_disrupted = ref_n_hits > 0 & alt_n_hits == 0,
            fimo_created   = ref_n_hits == 0 & alt_n_hits > 0,
            fimo_changed   = fimo_disrupted | fimo_created
          )

        fwrite(fimo_compare, file.path(OUT_DIR, "fimo_comparison.csv"))
        cat("FIMO comparison:", nrow(fimo_compare), "variant-motif pairs\n")
        cat("  Disrupted (ref only):", sum(fimo_compare$fimo_disrupted), "\n")
        cat("  Created (alt only):", sum(fimo_compare$fimo_created), "\n")

        # Add FIMO concordance to motifbreakR results
        mb_dt <- mb_dt %>%
          mutate(fimo_concordant = NA)

        # Match on SNP_id + motif name
        for (j in seq_len(nrow(mb_dt))) {
          snp_id <- mb_dt$SNP_id[j]
          # Try to match motif
          fimo_match <- fimo_compare %>%
            filter(sequence_name == snp_id & grepl(mb_dt$tf_name[j], motif_id, ignore.case = TRUE))
          if (nrow(fimo_match) > 0) {
            mb_dt$fimo_concordant[j] <- any(fimo_match$fimo_changed)
          }
        }

        n_concordant <- sum(mb_dt$fimo_concordant == TRUE, na.rm = TRUE)
        n_tested <- sum(!is.na(mb_dt$fimo_concordant))
        if (n_tested > 0) {
          cat("FIMO concordance:", n_concordant, "/", n_tested,
              "(", round(100 * n_concordant / n_tested, 1), "%)\n")
        }
      } else {
        cat("FIMO returned no hits\n")
      }
    } else {
      cat("FIMO output empty or missing\n")
    }
  } else {
    cat("WARNING: Could not match motifs for FIMO export\n")
  }

  # Clean up temp files
  unlink(fimo_dir, recursive = TRUE)
} else {
  if (!fimo_available) cat("FIMO not available — skipping validation\n")
  mb_dt$fimo_concordant <- NA
}

# ── 6. Write final output ───────────────────────────────────────────────────
cat("\n=== Writing outputs ===\n")

# Select relevant columns for output
out_cols <- c("SNP_id", "seqnames", "start", "end",
              "tf_name", "effect", "alleleDiff", "Refpvalue", "Altpvalue",
              "motif_in_disease_regulon", "fimo_concordant")

# Add regulon columns if available
if ("regulon_activity_diff" %in% colnames(mb_dt)) {
  out_cols <- c(out_cols, "regulon_activity_diff", "activity_padj")
}

# Keep only existing columns
out_cols <- intersect(out_cols, colnames(mb_dt))
motif_out <- mb_dt %>% select(all_of(out_cols))

# Add variant PIP
pip_lookup <- variants %>%
  mutate(SNP_id = paste(chromosome, position, allele1, allele2, sep = ":")) %>%
  # FIX (review A11#1): carry genome REF/ALT. SNP_id is in GWAS allele1:allele2 order,
  # but motifbreakR's alleleDiff sign is REF->ALT (genome). Exporting ref_genome/alt_genome
  # lets 56c align the eQTL Beta to the genome ALT instead of the GWAS allele2.
  select(SNP_id, max_pip, ref_genome = REF, alt_genome = ALT)
motif_out <- motif_out %>% left_join(pip_lookup, by = "SNP_id")

# Sort by PIP × |alleleDiff|
motif_out <- motif_out %>%
  mutate(priority_score = max_pip * abs(alleleDiff)) %>%
  arrange(desc(priority_score))

fwrite(motif_out, file.path(OUT_DIR, "motif_disruption_scores.csv"))
cat("Motif disruption scores:", nrow(motif_out), "variant-motif pairs written\n")

# Summary stats
cat("\n=== Summary ===\n")
cat("Total variant-motif disruptions:", nrow(motif_out), "\n")
cat("Unique variants with disruptions:", n_distinct(motif_out$SNP_id), "\n")
cat("Unique TF motifs disrupted:", n_distinct(motif_out$tf_name), "\n")
cat("  Strong effect:", sum(motif_out$effect == "strong", na.rm = TRUE), "\n")
cat("  In disease regulon:", sum(motif_out$motif_in_disease_regulon, na.rm = TRUE), "\n")

# Top disrupted TFs
top_tfs <- motif_out %>%
  count(tf_name, sort = TRUE) %>%
  head(10)
cat("\nTop 10 disrupted TF motifs:\n")
print(as.data.frame(top_tfs))

# Disease regulon disruptions
if (any(motif_out$motif_in_disease_regulon)) {
  reg_hits <- motif_out %>%
    filter(motif_in_disease_regulon) %>%
    select(SNP_id, tf_name, alleleDiff, max_pip) %>%
    arrange(desc(max_pip))
  cat("\nDisease regulon disruptions:\n")
  print(as.data.frame(head(reg_hits, 20)))
}

cat("\n============================================================\n")
cat("Script 56 complete. Results in:", OUT_DIR, "\n")
cat("============================================================\n")
