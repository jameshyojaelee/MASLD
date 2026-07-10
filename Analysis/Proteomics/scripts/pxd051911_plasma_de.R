#!/usr/bin/env Rscript
# pxd051911_plasma_de.R — MASLD-vs-control differential proteomics in the PXD051911 *plasma* arm.
# WHY: enables a within-COHORT liver->plasma translation analysis (the plasma twin of the liver DE that
# produced protein_transcript_concordance_v3.csv). PXD051911 deposited a plasma protein matrix from the
# SAME Spectronaut pipeline as the liver matrix (PRIDE 2025/03; provenance checksum-verified). Of the 143
# V1 plasma samples, 41 are the SAME patients (same visit) as liver samples -> genuine paired cohort,
# far more rigorous than the cross-cohort PXD052937 plasma (7 controls) used before.
#
# METHOD (modality-appropriate, mirrors differential_proteomics.R liver DE): limma on log2 DIA-MS
# intensities. Raw plasma matrix (linear) log2'd; design ~0 + group + age + bmi + sex + batch (plasma has
# 5 acquisition batches, and disease is NOT confounded with batch -> model batch as covariate rather than
# using the deposited batch-corrected matrix, keeping handling symmetric with liver). Contrast MASLD-Control.
# OUT: Analysis/Proteomics/results/pxd051911_plasma_differential.csv
suppressPackageStartupMessages({ library(data.table); library(limma) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DATA <- file.path(BASE, "data/PXD051911"); RES <- file.path(BASE, "Analysis/Proteomics/results")

# ── plasma matrix (raw linear intensities; cols = .raw filenames) ─────────────
ab <- fread(file.path(DATA, "plasma_protein_quant.txt"))
ab <- ab[Genes != "" & !is.na(Genes)]
ab[, gene := tstrsplit(Genes, ";", fixed = TRUE)[[1]]]
ab <- ab[!duplicated(gene)]
scols <- setdiff(names(ab), c("ProteinAccessions", "Genes", "ProteinDescriptions", "gene"))
mat <- as.matrix(ab[, ..scols]); rownames(mat) <- ab$gene; mode(mat) <- "numeric"
mat <- log2(mat)                                                  # linear -> log2 (matches liver handling)

# ── plasma sample metadata (V1 plasma arm) ────────────────────────────────────
md <- fread(file.path(DATA, "meta_data.txt"))
md <- md[plasma_proteomics_filename != "" & !is.na(plasma_proteomics_filename) &
         saf_diagnosis %in% c("MASH", "MASL", "No_MASLD")]
md[, `:=`(raw = plasma_proteomics_filename,
          group = factor(fifelse(saf_diagnosis == "No_MASLD", "Control", "MASLD"),
                         levels = c("Control", "MASLD")),
          age = as.numeric(alder), bmi = as.numeric(bmi),
          sex = factor(gender), batch = factor(plasma_batch_effects))]
md <- md[raw %in% colnames(mat)]
md <- md[stats::complete.cases(md[, .(group, age, bmi, sex, batch)])]  # limma needs complete design
mat <- mat[, md$raw, drop = FALSE]

# ── 50% per-protein NA filter (as liver) ──────────────────────────────────────
keep <- rowMeans(is.na(mat)) < 0.5
mat <- mat[keep, , drop = FALSE]

# paired-with-liver flag (same patient AND liver sample present)
n_paired <- md[liver_proteomics_filename != "" & !is.na(liver_proteomics_filename), .N]
cat(sprintf("[plasma-de] %d proteins x %d plasma samples (MASLD %d / Control %d); %d also liver-paired; %d batches\n",
            nrow(mat), ncol(mat), sum(md$group == "MASLD"), sum(md$group == "Control"),
            n_paired, nlevels(droplevels(md$batch))))

# ── limma DE: MASLD vs Control, + age/bmi/sex/batch ───────────────────────────
design <- model.matrix(~ 0 + group + age + bmi + sex + batch, data = md)
colnames(design) <- gsub("^group", "", colnames(design))
fit <- lmFit(mat, design)
cont <- makeContrasts(MASLD_vs_Control = MASLD - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont), robust = TRUE)
tt <- as.data.table(topTable(fit2, coef = "MASLD_vs_Control", number = Inf), keep.rownames = "gene")
setnames(tt, c("logFC", "adj.P.Val", "t", "P.Value"),
             c("plasma_logFC", "plasma_padj", "plasma_t", "plasma_p"), skip_absent = TRUE)
out <- tt[, .(gene, plasma_logFC, plasma_p, plasma_padj, plasma_t)]
fwrite(out, file.path(RES, "pxd051911_plasma_differential.csv"))
cat(sprintf("[plasma-de] wrote %d genes; sig (padj<0.05) = %d\n",
            nrow(out), out[plasma_padj < 0.05, .N]))
cat("[plasma-de] done\n")
