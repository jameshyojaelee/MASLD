#!/usr/bin/env Rscript
# 95_nmf_clean_ksweep.R — NMF k-sweep with comprehensive confounder strip
#
# Scope: remove standard NMF-confounding gene classes BEFORE top-IQR gene
# selection, so the decomposition allocates rank to MASLD disease biology
# rather than to sex, mitochondrial, ribosomal, or immunoglobulin axes.
#
# Exclusions (always on):
#   - chrY (all genes + Y-linked pseudogenes — full chromosome drop)
#   - chrM / MT-* (mitochondrial)
#   - RPS*/RPL*/MRPS*/MRPL* + rRNA biotypes (cytoplasmic + mito ribosomal)
#   - IG[HKL][VJCD] (immunoglobulin variable + constant regions)
#   - Hemoglobin (HBA*/HBB/HBD/HBE1/HBG*/HBM/HBQ1/HBZ)
#   - X-inactivation escapees: Oliva 2020 (PMID 32913072) + Tukiainen 2017 union
#     {XIST, TSIX, KDM6A, DDX3X, EIF1AX, UBA1, RLIM, ZFX, RPS4X, EIF2S3, TXLNG,
#      KDM5C, STS, SMC1A, NAA10, MAP7D2, PUDP, OFD1} — UTX removed (alias of KDM6A)
#   - MT-pseudogene symbols: MT[A-Z0-9]*P[0-9]+ (chr-1-resident mito pseudogenes
#     like MTCO1P20, MTATP6P1, MTND4P11 that escape the chrM chromosome drop)
#   - Pseudogene biotypes (processed/unprocessed/transcribed/unitary/IG/TR/rRNA_pseudogene)
#
# Opt-in exclusions (env var):
#   NMF_STRIP_HLA=1           → remove HLA-* (high-polymorphism donor variance)
#   NMF_PROTEIN_CODING_ONLY=1 → restrict input to gene_biotype=="protein_coding"
#
# Remediation run (bravo-coordinator, 2026-04-23):
#   Output cache path overridable via NMF_CACHE_PATH so PROT_ONLY=TRUE / FALSE
#   variants write side-by-side without clobbering.  K_SWEEP fit inline rather
#   than via shard split.
#
# Emits (same schema as Script 94):
#   RNA-seq/results/subtypes/nmf_ksweep_biology.md      (overwrites)
#   RNA-seq/results/subtypes/nmf_ksweep_program_atlas.csv
#   RNA-seq/results/subtypes/nmf_ksweep_rubric_scores.csv
#   RNA-seq/results/subtypes/nmf_results_cache_clean.rds
#   figures/misc/nmf_ksweep_diagnostic.pdf

suppressPackageStartupMessages({
  library(data.table); library(dplyr, warn.conflicts = FALSE)
  library(edgeR); library(limma); library(NMF)
  library(fgsea); library(msigdbr)
  library(ggplot2); library(patchwork); library(RColorBrewer)
  library(doParallel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

dge_path   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/merged_dge.rds")
meta_path  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
qc_path    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv")
atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
meta_gtf   <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
out_dir    <- file.path(BASE, "RNA-seq/results/subtypes")
# NMF_CACHE_PATH override lets remediation variants (prot-only vs not) write
# side-by-side without clobbering the canonical clean cache.
cache_path <- Sys.getenv("NMF_CACHE_PATH",
                         file.path(out_dir, "nmf_results_cache_clean.rds"))

# Output-file suffix derived from cache filename so PROT_ONLY variants don't
# clobber each other's biology.md / atlas.csv / rubric.csv. Default (canonical
# nmf_results_cache_clean.rds) keeps the original unsuffixed names.
.cache_base <- sub("\\.rds$", "", basename(cache_path))
.suffix <- sub("^nmf_results_cache_clean", "", .cache_base)  # "" or e.g. "_protonly"
md_out     <- file.path(out_dir, sprintf("nmf_ksweep_biology%s.md", .suffix))
atlas_out  <- file.path(out_dir, sprintf("nmf_ksweep_program_atlas%s.csv", .suffix))
rubric_out <- file.path(out_dir, sprintf("nmf_ksweep_rubric_scores%s.csv", .suffix))
pdf_out    <- file.path(FIG_MISC, sprintf("nmf_ksweep_diagnostic%s.pdf", .suffix))

K_SWEEP  <- 3:8
NMF_RUNS <- as.integer(Sys.getenv("NMF_RUNS", "20"))
N_TOP_GENES <- as.integer(Sys.getenv("NMF_N_TOP_GENES", "5000"))
NMF_CPUS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
STRIP_HLA   <- toupper(Sys.getenv("NMF_STRIP_HLA", "0")) %in% c("1","TRUE","T","YES")
PROT_ONLY   <- toupper(Sys.getenv("NMF_PROTEIN_CODING_ONLY", "0")) %in% c("1","TRUE","T","YES")
# NMF_INLINE_ALL_K=1 fits all K_SWEEP inline (used by remediation runs on bigmem).
# Default behaviour keeps legacy mode: inline only k=3,4 and expect shards for k>=5.
INLINE_ALL_K <- toupper(Sys.getenv("NMF_INLINE_ALL_K", "0")) %in% c("1","TRUE","T","YES")
cat(sprintf("=== 95 clean k-sweep (CPU=%d, nrun=%d, STRIP_HLA=%s, PROT_ONLY=%s) ===\n",
            NMF_CPUS, NMF_RUNS, STRIP_HLA, PROT_ONLY))

cl <- makeCluster(NMF_CPUS); registerDoParallel(cl)
on.exit(try(stopCluster(cl), silent = TRUE))

# ============================================================
# Canonical gene sets (for labeling only, not for stripping)
# ============================================================
FIBROTIC_CANON <- c(
  "COL1A1","COL1A2","COL3A1","COL5A1","COL6A3","COL15A1",
  "LUM","DCN","FBN1","LOX","SPARC","ELN","FN1",
  "ACTA2","DES","MYOCD","CNN1","TAGLN","MYL9","MYH11","FLNC",
  "PDGFRB","THY1","VIM","CDH11","LTBP2","CTSK",
  "TIMP1","PDPN","SMOC2","EDIL3","POSTN","RCN3","BGN",
  "LRRC15","SFRP4","PTPRQ")

# Lineage-marker signatures for tag_program() disambiguation (bravo-label-fixer, Task 9).
# Historical bug: pure stage_rho thresholding mis-labeled neutrophil programs
# (S100A8/A9-dominant) as "Quiescent-Parenchyma" and vSMC/pericyte programs
# (DES/CNN1/MYOCD with NO COL1A1/COL3A1 co-expression) as "Fibrogenic".
NEUTROPHIL_MARKERS <- c("S100A8","S100A9","CXCR1","CXCR2","FCGR3B",
                        "MPO","ELANE","LCN2","LTF","CSF3R")
VSMC_PERICYTE_MARKERS <- c("DES","CNN1","MYOCD","ACTG2","ACTA2","MYH11",
                           "TAGLN","MYL9","RGS5")
HSC_COLLAGEN_MARKERS <- c("COL1A1","COL3A1","CTHRC1","LUM","BGN")

# Biology-anchored labeling helper: inspects lineage markers BEFORE stage_rho rule.
tag_program <- function(top50_syms, fgsea_tbl, stage_rho, switch_ratio,
                        pct_dom) {
  sig <- as.data.table(fgsea_tbl)[padj < 0.05 & NES > 0][order(-NES)]
  emt_NES <- sig[grepl("EPITHELIAL_MESENCHYMAL|TGF_BETA|MATRISOME|COAGULATION",
                       pathway), max(NES, na.rm = TRUE)]
  inf_NES <- sig[grepl("TNFA|INFLAMMATORY|INTERFERON|IL6|IL2|COMPLEMENT",
                       pathway), max(NES, na.rm = TRUE)]
  if (!is.finite(emt_NES)) emt_NES <- 0
  if (!is.finite(inf_NES)) inf_NES <- 0

  top15 <- head(top50_syms, 15)
  n_neut <- sum(top15 %in% NEUTROPHIL_MARKERS)
  n_vsmc_top50 <- sum(top50_syms %in% VSMC_PERICYTE_MARKERS)
  n_hsc_col <- sum(top50_syms %in% HSC_COLLAGEN_MARKERS)
  n_fib <- sum(top50_syms %in% FIBROTIC_CANON)

  # --- Lineage overrides (evaluated BEFORE stage_rho rule) ---

  # (1) Neutrophil signature — >=3 markers in top-15 → Innate-immune-inflammation
  if (n_neut >= 3) return("Innate-immune-inflammation")

  # (2) vSMC/pericyte WITHOUT HSC collagen → Vascular-smooth-muscle / pericyte
  if (n_vsmc_top50 >= 3 && n_hsc_col == 0) {
    return("Vascular-smooth-muscle / pericyte")
  }

  # (3) vSMC + HSC collagen co-present → Activated-HSC / myofibroblast
  if (n_vsmc_top50 >= 2 && n_hsc_col >= 1) {
    return("Activated-HSC / myofibroblast")
  }

  # (4) HSC activation with EMT + rising along fibrosis (legacy Fibrogenic rule)
  if (n_fib >= 3 && emt_NES >= 2 && !is.na(stage_rho) && stage_rho > 0.5) {
    return("Fibrogenic")
  }

  # --- Fall-through: stage_rho directional rule ---
  if (is.na(stage_rho)) return("Unlabeled")

  if (stage_rho > 0.6) {
    if (inf_NES >= 2) return("Progression-Inflammatory")
    return("Progression-Other")
  }
  if (stage_rho < -0.6) return("Quiescent-Parenchyma")
  # Weakly monotonic or non-monotonic
  if (!is.na(switch_ratio) && switch_ratio > 1.3) return("Transition-Late")
  if (!is.na(switch_ratio) && switch_ratio < 0.77) return("Transition-Early")
  "Stable"
}

# ============================================================
# 1. Build clean matrix
# ============================================================
cat("\n=== Step 1: Build clean input matrix ===\n")
if (file.exists(cache_path)) {
  cat("  Loading existing clean cache\n")
  cached <- readRDS(cache_path)
  mat_nn <- cached$mat_nn; mat <- cached$mat
  sample_ids <- colnames(mat_nn); gene_ids <- rownames(mat_nn)
} else {
  dge <- readRDS(dge_path); meta <- fread(meta_path); qc <- fread(qc_path)
  meta <- meta %>% left_join(qc %>% select(sample_id, pass_technical),
                              by = "sample_id")
  disease_samples <- meta %>%
    filter(group_binary == "Disease" & pass_technical == TRUE) %>%
    pull(sample_id)
  available <- intersect(disease_samples, colnames(dge))
  dge_sub <- dge[, available]
  meta_sub <- meta %>% filter(sample_id %in% available)

  dge_sub <- calcNormFactors(dge_sub, method = "TMM")
  v <- voom(dge_sub, design = NULL, plot = FALSE)
  logcpm <- removeBatchEffect(v$E, batch = meta_sub$dataset)
  cat(sprintf("  Full log-CPM: %d genes x %d samples\n",
              nrow(logcpm), ncol(logcpm)))

  # Load GENCODE metadata + atlas symbol map
  gm <- fread(meta_gtf)
  atlas <- fread(atlas_path, select = c("ensembl_id","human_symbol"))
  atlas[, ens_base := sub("\\..*", "", ensembl_id)]
  id2sym <- setNames(atlas$human_symbol, atlas$ens_base)

  row_bases <- sub("\\..*", "", rownames(logcpm))
  row_syms  <- id2sym[row_bases]

  # Join GENCODE metadata by ensembl_base
  gm_u <- unique(gm, by = "ensembl_base")
  row_chr <- gm_u$chromosome[match(row_bases, gm_u$ensembl_base)]
  row_bio <- gm_u$gene_biotype[match(row_bases, gm_u$ensembl_base)]

  # -------- Confounder masks --------
  drop_chrY <- !is.na(row_chr) & row_chr %in% c("chrY")
  drop_chrM <- !is.na(row_chr) & row_chr %in% c("chrM","chrMT")
  drop_rRNA_biotype <- !is.na(row_bio) & row_bio %in% c("rRNA","Mt_rRNA","rRNA_pseudogene")

  # Symbol-based drops (use symbol if available, else skip)
  s <- row_syms
  s[is.na(s)] <- ""
  drop_rps_rpl  <- grepl("^RPS[0-9]+[A-Z]*$|^RPL[0-9]+[A-Z]*$", s) |
                   s %in% c("RPSA","RPLP0","RPLP1","RPLP2") |
                   grepl("^RPL[0-9]+P[0-9]+$|^RPS[0-9]+P[0-9]+$", s)
  drop_mrps_mrpl <- grepl("^MRPS[0-9]+[A-Z]*$|^MRPL[0-9]+[A-Z]*$", s)
  drop_ig <- grepl("^IG[HKL][VJCD][0-9]", s)
  drop_hb <- s %in% c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")

  # X-escape list: Oliva 2020 (PMID 32913072) + Tukiainen 2017 union.
  # UTX removed (alias of KDM6A — was duplicate). Added: KDM5C, STS, SMC1A,
  # NAA10, MAP7D2, PUDP, OFD1 (established escapees with >0 allelic expression
  # in reference cohorts).
  X_ESCAPE_GENES <- c(
    # Original 11 members (UTX dropped as KDM6A alias)
    "XIST","TSIX","KDM6A","DDX3X","EIF1AX","UBA1","RLIM","ZFX","RPS4X",
    "EIF2S3","TXLNG",
    # Oliva 2020 + Tukiainen 2017 additions
    "KDM5C","STS","SMC1A","NAA10","MAP7D2","PUDP","OFD1"
  )
  drop_xesc <- s %in% X_ESCAPE_GENES

  # MT-pseudogene strip: chr-resident mito pseudogenes (MTCO1P20, MTATP6P1,
  # MTND4P11, MTCYBP45, ...) that are named after mito genes but sit on
  # autosomes and so escape the chrM chromosome drop.  Regex matches
  # MT[letters][optional digits]P[digits]$ (all MT-* pseudogene families).
  drop_mt_pseudo <- grepl("^MT[A-Z]+[0-9]*L?P[0-9]+$", s)

  # Generic pseudogene biotype strip — applied to ALL pseudogene biotypes
  # per GENCODE v49.  When PROT_ONLY=TRUE this is redundant with drop_nonpc
  # but harmless; when PROT_ONLY=FALSE it removes the N>21k pseudogenes that
  # would otherwise dominate top-IQR via low-expression stochasticity.
  PSEUDOGENE_BIOTYPES <- c("processed_pseudogene","unprocessed_pseudogene",
                            "transcribed_processed_pseudogene",
                            "transcribed_unprocessed_pseudogene",
                            "translated_processed_pseudogene",
                            "translated_unprocessed_pseudogene",
                            "transcribed_unitary_pseudogene",
                            "unitary_pseudogene","pseudogene",
                            "IG_V_pseudogene","IG_C_pseudogene","IG_J_pseudogene",
                            "IG_pseudogene","TR_V_pseudogene","TR_J_pseudogene",
                            "polymorphic_pseudogene","rRNA_pseudogene")
  drop_pseudogene_biotype <- !is.na(row_bio) & row_bio %in% PSEUDOGENE_BIOTYPES

  # Belt-and-suspenders: P[0-9]+$-terminated symbols whose biotype is missing
  # OR non-protein-coding.  Gate by biotype so we don't accidentally strip
  # protein-coding genes like NLRP3/NLRP6/AQP1 that happen to match the regex.
  drop_P_suffix_symbol <- grepl("^.*P[0-9]+$", s) &
                          (is.na(row_bio) | row_bio != "protein_coding")

  # Optional masks
  drop_hla <- if (STRIP_HLA) grepl("^HLA-", s) else rep(FALSE, length(s))
  drop_nonpc <- if (PROT_ONLY) (!is.na(row_bio) & row_bio != "protein_coding") else rep(FALSE, length(s))

  drop_any <- drop_chrY | drop_chrM | drop_rRNA_biotype |
              drop_rps_rpl | drop_mrps_mrpl | drop_ig | drop_hb | drop_xesc |
              drop_mt_pseudo | drop_pseudogene_biotype | drop_P_suffix_symbol |
              drop_hla | drop_nonpc

  cat("  Dropping by class:\n")
  cat(sprintf("    chrY                %d\n", sum(drop_chrY)))
  cat(sprintf("    chrM                %d\n", sum(drop_chrM)))
  cat(sprintf("    rRNA biotype        %d\n", sum(drop_rRNA_biotype)))
  cat(sprintf("    RPS*/RPL*           %d\n", sum(drop_rps_rpl)))
  cat(sprintf("    MRPS*/MRPL*         %d\n", sum(drop_mrps_mrpl)))
  cat(sprintf("    IG[HKL][VJCD]       %d\n", sum(drop_ig)))
  cat(sprintf("    Hemoglobin          %d\n", sum(drop_hb)))
  cat(sprintf("    X-escapees (18)     %d\n", sum(drop_xesc)))
  cat(sprintf("    MT-pseudogenes      %d\n", sum(drop_mt_pseudo)))
  cat(sprintf("    pseudogene biotype  %d\n", sum(drop_pseudogene_biotype)))
  cat(sprintf("    P\\d+ symbol suffix  %d  (biotype-gated)\n", sum(drop_P_suffix_symbol)))
  cat(sprintf("    HLA (opt)           %d  [STRIP_HLA=%s]\n", sum(drop_hla), STRIP_HLA))
  cat(sprintf("    non-PC (opt)        %d  [PROT_ONLY=%s]\n", sum(drop_nonpc), PROT_ONLY))
  cat(sprintf("    TOTAL DROP          %d\n", sum(drop_any)))

  logcpm <- logcpm[!drop_any, ]
  cat(sprintf("  Clean log-CPM: %d genes x %d samples\n",
              nrow(logcpm), ncol(logcpm)))

  # Top IQR
  gene_iqr <- apply(logcpm, 1, IQR)
  top_genes <- names(sort(gene_iqr, decreasing = TRUE))[1:min(N_TOP_GENES, length(gene_iqr))]
  mat <- logcpm[top_genes, ]
  mat_nn <- mat - apply(mat, 1, min); mat_nn[mat_nn < 0] <- 0
  sample_ids <- colnames(mat_nn); gene_ids <- rownames(mat_nn)

  # Sanity assertions — use same regexes as the strip, so we verify strip
  # actually removed the targeted classes (kinases like RPS6KA1/RPS6KB* and
  # IG superfamily members like IGF1 are NOT ribosomal/immunoglobulin).
  kept_syms <- id2sym[sub("\\..*", "", gene_ids)]
  kept_syms[is.na(kept_syms)] <- ""
  stopifnot("Y pseudogene leaked" =
    sum(kept_syms %in% c("BCORP1","CDY4P","ANOS2P","AGKP1")) == 0)
  stopifnot("XIST/TSIX leaked" = sum(kept_syms %in% c("XIST","TSIX")) == 0)
  stopifnot("MT- leaked"  = sum(startsWith(kept_syms, "MT-")) == 0)
  stopifnot("MT-pseudogene leaked" =
    sum(grepl("^MT[A-Z]+[0-9]*L?P[0-9]+$", kept_syms)) == 0)
  stopifnot("X-escape leaked" = sum(kept_syms %in% X_ESCAPE_GENES) == 0)
  stopifnot("RPS/RPL subunit leaked" =
    sum(grepl("^RPS[0-9]+[A-Z]*$|^RPL[0-9]+[A-Z]*$", kept_syms)) == 0)
  stopifnot("Ig variable leaked" =
    sum(grepl("^IG[HKL][VJCD][0-9]", kept_syms)) == 0)
  # Positive-control assertions: ensure we retained key protein-coding
  # genes that end in Pn but are NOT pseudogenes (NLRP3, NLRP6, AQP1, AQP4)
  for (g in c("NLRP3","NLRP6","AQP1","AQP4")) {
    if (g %in% row_syms && !(g %in% kept_syms)) {
      warning(sprintf("POSITIVE CONTROL %s incorrectly stripped — check regex gating", g))
    }
  }
  cat("  Sanity asserts passed.\n")
  # Report retained kinase/superfamily cousins for transparency
  n_rps_kin <- sum(grepl("^RPS[0-9]+K", kept_syms))   # RPS6KA*, RPS6KB* kinases
  if (n_rps_kin > 0)
    cat(sprintf("  (retained %d RPS-kinase protein-coding genes e.g. RPS6KA1)\n", n_rps_kin))

  saveRDS(list(nmf_results = list(), metrics = data.frame(),
               mat_nn = mat_nn, mat = mat), cache_path)
  cached <- readRDS(cache_path)
}

meta <- fread(meta_path)[sample_id %in% sample_ids]; setkey(meta, sample_id)
atlas <- fread(atlas_path, select = c("ensembl_id","human_symbol"))
id2sym <- setNames(atlas$human_symbol, atlas$ensembl_id)
strip_ver <- function(x) sub("\\.[0-9]+$", "", x)
sym_W <- id2sym[strip_ver(gene_ids)]; names(sym_W) <- gene_ids

# ============================================================
# 2. NMF fits (only k=3,4 here; k=5..8 via parallel shards)
# ============================================================
cat("\n=== Step 2: NMF fits ===\n")
nmf_list <- cached$nmf_results
metrics <- if (!is.null(cached$metrics) && nrow(cached$metrics)) cached$metrics else data.frame()
inline_ks <- if (INLINE_ALL_K) K_SWEEP else c(3, 4)
for (k in inline_ks) {
  if (as.character(k) %in% names(nmf_list)) {
    cat(sprintf("  k=%d: reusing cached\n", k)); next
  }
  cat(sprintf("  k=%d: fitting (nrun=%d)\n", k, NMF_RUNS))
  t0 <- Sys.time()
  set.seed(42 + k)
  fit <- nmf(mat_nn, rank = k, nrun = NMF_RUNS,
             method = "brunet", seed = "random",
             .pbackend = "par", .options = list(verbose = FALSE))
  nmf_list[[as.character(k)]] <- fit
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("    done in %.1f min\n", el))
  metrics <- rbind(metrics, data.frame(
    k = k, cophenetic = cophcor(fit), dispersion = dispersion(fit),
    silhouette = mean(silhouette(fit, what = "consensus")[, "sil_width"])
  ))
  saveRDS(list(nmf_results = nmf_list, metrics = metrics,
               mat_nn = mat_nn, mat = mat), cache_path)
}

# If shards for k=5..8 haven't been run yet, bail out of rubric step with a note
needed_k <- setdiff(as.character(K_SWEEP), names(nmf_list))
if (length(needed_k)) {
  cat(sprintf("\nMissing k values: %s\n", paste(needed_k, collapse = ",")))
  cat("Submit parallel shards with NMF_K=<k> NMF_CACHE_PATH=", cache_path, " ",
      "via run_94a_single_k.sbatch, then run 94b merge, then rerun 95 for rubric.\n", sep = "")
  cat("Exiting cleanly (Step 3+ rubric skipped).\n")
  quit(save = "no", status = 0)
}

# ============================================================
# 3. Per-k characterization + rubric (runs only once all k cached)
# ============================================================
cat("\n=== Step 3: Per-k characterization ===\n")
hallmark  <- msigdbr(species = "Homo sapiens", collection = "H")
hall_list <- split(hallmark$gene_symbol, hallmark$gs_name)

atlas_rows <- list(); rubric_rows <- list(); per_k_summary <- list()

for (k in K_SWEEP) {
  cat(sprintf("  -- k=%d --\n", k))
  res_k <- nmf_list[[as.character(k)]]
  W <- basis(res_k); H <- coef(res_k)
  rownames(W) <- gene_ids; colnames(H) <- sample_ids
  W_centered <- sapply(seq_len(k), function(j) W[, j] - rowMeans(W[, -j, drop = FALSE]))
  rownames(W_centered) <- rownames(W)
  dom <- apply(H, 2, which.max)

  # Stage-level dominance rate per program
  fib_int <- as.integer(meta$fibrosis_stage[match(sample_ids, meta$sample_id)])
  stage_rho <- rep(NA_real_, k)
  for (j in seq_len(k)) {
    prop <- vapply(0:4, function(s) {
      n_s <- sum(fib_int == s, na.rm = TRUE)
      if (n_s == 0) NA_real_ else mean(dom[fib_int == s] == j, na.rm = TRUE)
    }, numeric(1))
    if (sum(!is.na(prop)) >= 3) {
      stage_rho[j] <- suppressWarnings(cor(prop, 0:4, method = "spearman", use = "pairwise"))
    }
  }

  prog_rows <- list()
  for (j in seq_len(k)) {
    w_j <- W_centered[, j]
    ord <- order(w_j, decreasing = TRUE)
    top50_ids <- gene_ids[ord][1:50]
    top50_syms <- unname(sym_W[top50_ids])
    top50_syms <- top50_syms[!is.na(top50_syms) & top50_syms != ""]

    has_sym <- !is.na(sym_W) & sym_W != ""
    ranked <- w_j[has_sym]; names(ranked) <- sym_W[has_sym]
    ranked <- ranked[!duplicated(names(ranked))]
    ranked <- sort(ranked, decreasing = TRUE)
    set.seed(42)
    fg <- suppressWarnings(fgsea(pathways = hall_list, stats = ranked,
                                 minSize = 10, maxSize = 500,
                                 nPermSimple = 10000, scoreType = "std"))

    h_j <- H[j, ]
    rho_fib <- suppressWarnings(cor(h_j, fib_int, method = "spearman", use = "pairwise"))
    ok_f <- !is.na(fib_int)
    m_lo <- mean(h_j[ok_f & fib_int <= 2], na.rm = TRUE)
    m_hi <- mean(h_j[ok_f & fib_int >= 3], na.rm = TRUE)
    sw_ratio <- m_hi / (m_lo + 1e-10)
    pct_dom <- mean(dom == j) * 100
    ds_vec <- meta$dataset[match(sample_ids, meta$sample_id)]
    kw <- suppressWarnings(kruskal.test(h_j ~ factor(ds_vec)))
    n_ds <- length(unique(ds_vec))
    eta2 <- max(0, (as.numeric(kw$statistic) - n_ds + 1) / (length(h_j) - n_ds))
    sex_vec <- meta$sex[match(sample_ids, meta$sample_id)]
    ok_sx <- sex_vec %in% c("M","F")
    rho_sex <- if (sum(ok_sx) > 10)
      suppressWarnings(cor(h_j[ok_sx], as.integer(sex_vec[ok_sx] == "F"),
                           method = "spearman")) else NA_real_

    tag <- tag_program(top50_syms, fg, stage_rho[j], sw_ratio, pct_dom)

    sigp <- as.data.table(fg)[padj < 0.05 & NES > 0][order(-NES)][1:3, .(pathway, NES)]
    top_path <- paste(sprintf("%s(%.2f)", sigp$pathway, sigp$NES), collapse = ";")

    prog_rows[[j]] <- data.table(
      k = k, program = paste0("P", j), label = tag,
      top10_genes = paste(head(top50_syms, 10), collapse = ","),
      top_pathways = top_path,
      rho_fibrosis = rho_fib, stage_rho = stage_rho[j],
      rho_sex_F = rho_sex, eta2_dataset = eta2,
      switch_ratio = sw_ratio, pct_dominant = pct_dom,
      n_fibrotic_canon = sum(top50_syms %in% FIBROTIC_CANON)
    )
  }
  prog_dt <- rbindlist(prog_rows); atlas_rows[[as.character(k)]] <- prog_dt

  # Rubric
  Hcor <- cor(t(H)); max_off <- max(abs(Hcor[lower.tri(Hcor)]))
  labels_k <- prog_dt$label
  fib_rows <- prog_dt[label == "Fibrogenic"]
  n_fib <- nrow(fib_rows)
  fib_clean <- if (n_fib == 1) {
    abs(fib_rows$rho_sex_F) < 0.2 & fib_rows$eta2_dataset < 0.2 &
    !is.na(fib_rows$stage_rho) & fib_rows$stage_rho > 0.5
  } else FALSE
  fib_canon_ok <- if (n_fib == 1) fib_rows$n_fibrotic_canon >= 3 else FALSE

  rubric <- data.table(
    k = k,
    crit1_all_labeled = (sum(labels_k == "Unlabeled") == 0),
    crit2_one_fibrogenic = (n_fib == 1),
    crit3_fib_clean = (n_fib == 1 && fib_clean && fib_canon_ok),
    crit4_any_progression = any(grepl("Progression", labels_k)),
    crit5_any_quiescent = any(labels_k == "Quiescent-Parenchyma"),
    crit6_orthogonal = (max_off < 0.6),
    max_H_cor = max_off,
    min_coverage_pct = min(prog_dt$pct_dominant),
    n_programs = k, score = NA_integer_
  )
  rubric[, score := sum(unlist(.SD)),
         .SDcols = grep("^crit", names(rubric), value = TRUE)]
  rubric_rows[[as.character(k)]] <- rubric
  per_k_summary[[as.character(k)]] <- list(prog_dt = prog_dt, Hcor = Hcor)
  cat(sprintf("    rubric=%d/6, max|H|=%.2f, min_cov=%.1f%%\n",
              rubric$score, max_off, rubric$min_coverage_pct))
}

atlas_dt  <- rbindlist(atlas_rows)
rubric_dt <- rbindlist(rubric_rows)
fwrite(atlas_dt,  atlas_out)
fwrite(rubric_dt, rubric_out)

# Recommendation
full_pass <- rubric_dt[score == 6][order(k)]
if (nrow(full_pass) > 0) {
  chosen_k <- full_pass$k[1]
  rationale <- sprintf("Full pass: smallest k with 6/6 = k=%d", chosen_k)
} else {
  best <- rubric_dt[order(-score, k)]; chosen_k <- best$k[1]
  rationale <- sprintf("No k fully passes. Best k=%d with %d/6; next: %s",
                       chosen_k, best$score[1],
                       paste(sprintf("k=%d(%d/6)", best$k[2:min(nrow(best),4)],
                                     best$score[2:min(nrow(best),4)]),
                             collapse=", "))
}
cat(sprintf("\nRecommended k = %d\n  %s\n", chosen_k, rationale))

# Markdown + PDF (same shape as Script 94)
md <- c(
  "# NMF Clean k-Sweep (comprehensive confounder strip)",
  sprintf("**Generated:** %s  ", format(Sys.time())),
  sprintf("**Samples:** %d disease-only, QC-passing", length(sample_ids)),
  sprintf("**Gene pool:** %d top-IQR after strip  (STRIP_HLA=%s, PROT_ONLY=%s)",
          length(gene_ids), STRIP_HLA, PROT_ONLY),
  sprintf("**k range:** %d..%d, NMF nrun=%d", min(K_SWEEP), max(K_SWEEP), NMF_RUNS),
  "", sprintf("## RECOMMENDED k: **%d**", chosen_k), "", rationale, "",
  "## Rubric scores by k", "",
  paste(capture.output(print(rubric_dt[, .(k, score, crit1_all_labeled,
                                           crit2_one_fibrogenic,
                                           crit3_fib_clean,
                                           crit4_any_progression,
                                           crit5_any_quiescent,
                                           crit6_orthogonal,
                                           min_coverage_pct = round(min_coverage_pct,1),
                                           max_H_cor = round(max_H_cor,2))],
                            row.names = FALSE)), collapse = "\n"),
  "", "---", "", "## Per-k program atlas"
)
for (k in K_SWEEP) {
  prog_dt <- per_k_summary[[as.character(k)]]$prog_dt
  md <- c(md, sprintf("### k=%d", k), "")
  tbl <- prog_dt[, .(program, label, pct_dominant = round(pct_dominant,1),
                     stage_rho = round(stage_rho,3),
                     rho_sex_F = round(rho_sex_F,3),
                     eta2_dataset = round(eta2_dataset,3),
                     switch_ratio = round(switch_ratio,2),
                     fib_canon = n_fibrotic_canon,
                     top_genes = substr(top10_genes, 1, 70))]
  md <- c(md, paste(capture.output(print(tbl, row.names = FALSE)), collapse = "\n"), "")
}
writeLines(md, md_out); cat(sprintf("Wrote %s\n", md_out))

# PDF
rh <- melt(rubric_dt[, .(k, crit1_all_labeled, crit2_one_fibrogenic,
                         crit3_fib_clean, crit4_any_progression,
                         crit5_any_quiescent, crit6_orthogonal)],
           id.vars = "k", variable.name = "criterion", value.name = "pass")
rh[, criterion := sub("^crit[0-9]+_", "", criterion)]
rh[, pass := as.integer(pass)]
pA <- ggplot(rh, aes(factor(k), criterion, fill = factor(pass))) +
  geom_tile(color = "white") + geom_text(aes(label = ifelse(pass==1,"o","x")), size = 3) +
  scale_fill_manual(values = c("0"="#EF9A9A","1"="#A5D6A7"), guide = "none") +
  labs(x = "k", y = NULL, title = sprintf("Rubric (clean strip); recommended k=%d", chosen_k)) +
  theme_masld()

label_long <- rbindlist(lapply(K_SWEEP, function(k) {
  per_k_summary[[as.character(k)]]$prog_dt[, .(k = k, label)]
}))
pB <- ggplot(label_long[, .N, by = .(k, label)], aes(factor(k), N, fill = label)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3) +
  labs(x = "k", y = "programs", title = "Label composition per k") +
  theme_masld()

fib_track <- atlas_dt[label == "Fibrogenic"]
pC <- if (nrow(fib_track)) {
  ggplot(fib_track, aes(factor(k), stage_rho, fill = switch_ratio)) +
    geom_col(color = "white") +
    geom_text(aes(label = sprintf("ρ=%.2f\nsw=%.2f", stage_rho, switch_ratio)),
              vjust = -0.2, size = 2.3) +
    labs(x = "k", y = "stage ρ (dominance vs F stage)",
         title = "Fibrogenic program vs k") + theme_masld()
} else ggplot() + labs(title = "No Fibrogenic at any k") + theme_masld()

cov_track <- atlas_dt[, .(min_cov = min(pct_dominant),
                          max_cov = max(pct_dominant)), by = k]
pD <- ggplot(cov_track, aes(factor(k))) +
  geom_linerange(aes(ymin = min_cov, ymax = max_cov), color = "#1565C0", linewidth = 1.5) +
  geom_point(aes(y = min_cov), color = "#C2185B", size = 2) +
  geom_point(aes(y = max_cov), color = "#1565C0", size = 2) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "grey40") +
  labs(x = "k", y = "% coverage", title = "Coverage floor") + theme_masld()

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(title = sprintf("NMF clean k-sweep — recommended k=%d", chosen_k))
ggsave(pdf_out, fig, width = 14, height = 10, device = cairo_pdf)
cat("Wrote", pdf_out, "\n")
cat("Done.\n")
