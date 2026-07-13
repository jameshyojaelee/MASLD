#!/usr/bin/env Rscript
# rebuild_cas13_library.R
# ---------------------------------------------------------------------------
# Build cas13_library.csv under the v9 definition (2026-06-27). v9 changes
# the target-gene strategy from v6: the human arm broadens from the integrated
# spine alone to integrated UNION cohort-replicated DEGs, COLOC genetic targets
# and positive controls are folded in, and the MASH + miRNA tiers are removed.
#
#   Target selection moved to TREAT FDR (2026-06-27). Every DE arm is now gated by an
#   ANALYTICAL TREAT reconstruction (add_treat_fdr): the moderated-t statistic is
#   re-tested against an explicit logFC offset `lfc` (a two-sided shifted-t tail using
#   the topTable moderated SE and a median-inferred residual+prior df), then BH-adjusted;
#   the gate is fdr_treat < 0.05 & logFC > 0 (UP). This REPLACES the prior raw-padj /
#   ashr-lfsr selection (the old "padj<0.05 & logFC>0.5" / "lfsr<0.05 & shrunk_logFC>0.5"
#   point-estimate cutoffs). The `lfc` offset, NOT a point-estimate floor, now encodes
#   effect-size stringency, applied BIOTYPE-WISE (add_treat_fdr_bt): protein-coding genes
#   use lfc = SHRUNK_LFC_THR (env CAS13_PC_LFC) on the human C1 + per-cohort axes and the
#   mouse PC arm; lncRNAs use lfc = 0.0 everywhere (TREAT@0 == ordinary moderated-t
#   FDR<0.05, no effect-size floor). The floor is keyed on BIOTYPE, not on which contrast a
#   gene enters through, so the advertised PC LFC floor binds on protein-coding genes
#   regardless of axis (this also closed an earlier C2/C3 lfc=0-for-all leak; C2/C3 have
#   since been dropped as entry axes entirely, 2026-06-29).
#
#   LIBRARY = CORE              (mouse orthologs of UP human DEGs from the integrated
#                                disease-vs-control axis ONLY -- canonical limma-voom QW
#                                C2 -- gated TREAT FDR<0.05 & logFC>0, biotype-split lfc
#                                (PC=SHRUNK_LFC_THR, lncRNA=0). The MASH-vs-MASL (C2) and
#                                advanced-vs-early-fibrosis (C3) progression axes were
#                                DROPPED 2026-06-29 (PI); they are no longer entry axes.)
#             UNION
#             COHORT-REPLICATED (mouse orthologs of human DEGs significant [TREAT FDR<0.05,
#                                lfc=0.5, UP] in >=N of the 5 control-bearing cohorts;
#                                N = COHORT_MIN, default 2)
#             UNION
#             MOUSE-CONFIRMED   (mouse cross-diet UP, TREAT FDR<0.05 per diet [PC lfc=0.5,
#                                lncRNA lfc=0.0] in >=3 of 4 diets, AND human ortholog
#                                logFC>0 -- human directional concordance, never mouse-alone)
#             UNION
#             COLOC             (canonical SuSiE COLOC PP.H4 > 0.5, gene_level_coloc.csv;
#                                genetic causal support, NOT the permissive ABF fallback)
#             UNION
#             POSITIVE CONTROLS (65 curated MASLD-biology genes, folded in as targets)
#   restricted to protein_coding + lncRNA. NO miRNA, NO MASH tier. CORE = the single C1
#   integrated disease-vs-control axis (PC at lfc=SHRUNK_LFC_THR, lncRNA no-floor); the
#   human spine is CORE UNION COHORT-REPLICATED only (C2/C3 dropped 2026-06-29).
#   Screen-expression gate (v9): the screen is hepatocyte-autonomous, so gate on the
#   MOUSE-HEPATOCYTE scRNA pseudobulk substrate (mouse_hep_cpm) with biotype-aware
#   floors: protein_coding >=1.0 CPM, lncRNA >=0.1 CPM. The metric is CPM (10x UMI
#   counts / library-size x 1e6, hepatocyte pseudobulk) -- the correct unit for droplet
#   scRNA, NOT length-normalized TPM. A whole-liver bulk gate (the prior MCD-TPM gate)
#   admits non-parenchymal passengers (e.g. Col1a1, stellate) that cannot move a cell-
#   autonomous hepatocyte-lipid readout, so the hepatocyte substrate is the correct gate.
#   Positive controls AND high-COLOC genes (SuSiE PP.H4 >= CAS13_COLOC_GATE_EXEMPT_PP4,
#   default 0.9) remain exempt and are kept even if hep-absent, but flagged
#   mouse_untestable. In-house MCD bulk TPM (mcd_mean_tpm) is retained as an ANNOTATION
#   column only; human hep scRNA stays a separate SOFT tag (hep_substrate). lncRNA is
#   still demoted to an EXPLORATORY arm (library_arm); human lncRNA remains raw-FDR/no-LFC.
#
# tier priority on overlap: positive_control > core > cohort_replicated >
#   mouse_confirmed > coloc. has_human_de = core OR cohort. is_mash_deg kept as a soft
#   annotation only (adds no genes).
#
# v7 SOURCE CHANGE: the human arm is sourced from the paper-canonical limma-voom
# quality-weighted C2 analysis (canonical_deg_results.csv) + the per-study cohort
# DEGs, replacing the v6 metafor spine -- this matches the integrated-vs-per-study
# venn and the project-wide canonical DEG method.
#
# Consumed by figS_cas13_library_*.R. Output filename kept as cas13_library.csv
# for consumer compatibility; library_version column = v9.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

BASE     <- Sys.getenv("MASLD_PROJECT_ROOT",
                       "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PERDIET  <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet_cas13")  # Cas13 library Western pool; decoupled from paper 4-model per_diet (2026-06-16)
PERSTUDY <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study")
OUT      <- file.path(BASE, "Cas13_Library_Design/data/cas13_library.csv")
BACKUP   <- file.path(BASE, "Cas13_Library_Design/data/cas13_library_v3.0_pre_v9.csv")
V7REF    <- file.path(BASE, "Cas13_Library_Design/data/cas13_library_v8_canonical.csv")  # frozen v8, baseline for v8->v9 diff
DIFFOUT  <- file.path(BASE, "Cas13_Library_Design/data/cas13_library_v8_to_v9_diff.csv")
MANIFEST <- file.path(BASE, "Cas13_Library_Design/data/BUILD_MANIFEST.txt")
ATLAS    <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
ORTHO    <- file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz")
META     <- file.path(BASE, "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv")
HEPSPEC  <- file.path(BASE, "Cas13_Library_Design/data/hep_specificity.csv")
# Mouse MCD mean TPM -> hard screen-expression gate + annotation.
MCDTPM   <- file.path(BASE, "RNA-seq/Mouse/InHouse_MCD/results/mean_tpm_mcd.csv")
# Mouse-hepatocyte scRNA substrate -> annotation only. Now the vM38 re-quantification
# (2026-06-25): full GENCODE vM38 atlas (114 samples, custom permissive-lncRNA cellranger
# ref), native vM38 gene_ids matching the library, all loci measured. vM37<->vM38 hep_cpm
# Spearman rho = 0.972 on shared genes. Old vM37 substrate: mouse_hep_specificity.csv.
MOUSEHEP <- file.path(BASE, "Cas13_Library_Design/data/mouse_hep_specificity_vm38.csv")
# v7: human disease-vs-control arm = canonical limma-voom QW C2 (integrated).
CANONICAL<- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
MASHMASL <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/mash_vs_masl_dream.csv")
ADVFIB   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/adv_vs_early_fibrosis_dream.csv")
# COLOC tier = canonical SuSiE PP.H4 > 0.5 (gene-level, EUR/EAS/AFR/SAS portfolio).
COLOCFILE<- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
# Positive controls (human symbols + steatosis-direction label).
POSCTRL  <- file.path(BASE, "results/library/positive_control.csv")
# MASH-vs-MASL kept for the is_mash_deg ANNOTATION only (no longer an entry tier).
MASH     <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/nafl_vs_nash_meta_ashr.csv")

# 5 control-bearing cohorts (disease-vs-control); per-study cohort replication.
COHORTS    <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
COHORT_MIN <- as.integer(Sys.getenv("CAS13_COHORT_MIN", "2"))  # >=N cohorts for the cohort tier

DIETS          <- c("MCD", "CDAHFD", "Western", "HFD")
MIN_DIETS      <- 3L       # mouse cross-diet replication
LFSR_THR       <- 0.05     # FDR cutoff applied to fdr_treat (name is historical; not an lfsr now)
# NOTE: despite the name, SHRUNK_LFC_THR is now the TREAT lfc OFFSET (the |logFC| the
# moderated-t is tested against), NOT a point-estimate cutoff on a shrunk logFC. It is
# reused as the lfc offset for every PROTEIN-CODING arm (all human axes C1/C2/C3 + per-
# cohort + mouse cross-diet PC); lncRNAs use SHRUNK_LFC_THR_LNC=0. Kept as-is for ref
# stability; rename deferred (touches many call sites).
SHRUNK_LFC_THR <- as.numeric(Sys.getenv("CAS13_PC_LFC", "0.2"))  # TREAT lfc offset, PC arms, UP only (env-tunable; 0.2 -> ~2,177 lib, 2026-06-29 PI sizing)
SHRUNK_LFC_THR_LNC <- 0.0  # TREAT lfc offset, lncRNA (no effect-size floor)
HUMAN_PADJ_THR <- 0.05     # human arm FDR cutoff (applied to fdr_treat)
COLOC_PP4_THR  <- 0.5      # canonical SuSiE PP.H4 threshold (tier + annotation)
KEEP_BIOTYPES  <- c("protein_coding", "lncRNA")
HEP_CPM_PC_GATE  <- as.numeric(Sys.getenv("CAS13_PC_HEP_CPM_GATE", "1.0"))   # mouse-hepatocyte scRNA pseudobulk CPM floor (PC)
HEP_CPM_LNC_GATE <- as.numeric(Sys.getenv("CAS13_LNC_HEP_CPM_GATE", "0.1"))  # mouse-hepatocyte scRNA pseudobulk CPM floor (lncRNA)
# high-COLOC genes (SuSiE PP.H4 >= this) are exempt (kept, flagged mouse_untestable).
COLOC_GATE_EXEMPT_PP4 <- as.numeric(Sys.getenv("CAS13_COLOC_GATE_EXEMPT_PP4", "0.9"))
LIB_VERSION    <- "v9"

# Positive-control roster aliases (roster colloquial name -> HGNC symbol).
POSCTRL_ALIASES <- c("SCD1" = "SCD")
# Drop systemic / non-parenchymal "controls": they act through non-hepatocyte or
# whole-body mechanisms (incretin receptors, fibroblast FAP) and cannot move a
# cell-autonomous hepatocyte-lipid readout regardless of expression -- keeping them
# would only deflate the positive-control AUROC / dynamic range.
POSCTRL_EXCLUDE <- c("GIPR", "GLP1R", "FAP")
# Override mis-mapped control orthologs: human SCD -> mouse Scd1 (the liver paralog;
# the ortholog table's deterministic pick is Scd3, a skin paralog ~0 TPM in liver).
POSCTRL_MOUSE_OVERRIDE <- c("SCD" = "Scd1")   # values are MOUSE symbols

strip_v <- function(x) sub("[.][0-9]+$", "", x)
infer_df_total <- function(dt) {
  pr <- dt[is.finite(t) & is.finite(P.Value) & P.Value > 0 & P.Value < 1 & abs(t) > 1e-6]
  idx <- unique(round(seq(1, nrow(pr), length.out = min(nrow(pr), 12))))
  median(vapply(idx, function(i)
    uniroot(function(df) 2 * pt(-abs(pr$t[i]), df = df) - pr$P.Value[i], c(0.1, 1e6))$root,
    numeric(1)))
}
add_treat_fdr <- function(dt, lfc, se_col = NULL, df_col = NULL) {
  out <- copy(dt)
  se <- if (!is.null(se_col) && se_col %in% names(out)) out[[se_col]] else abs(out$logFC / out$t)
  se[!is.finite(se) | se <= 0] <- NA_real_
  df_use <- if (!is.null(df_col) && df_col %in% names(out)) out[[df_col]] else infer_df_total(out)
  out[, p_treat := pt((abs(logFC) - lfc) / se, df = df_use, lower.tail = FALSE) +
                   pt((abs(logFC) + lfc) / se, df = df_use, lower.tail = FALSE)]
  out[, fdr_treat := p.adjust(p_treat, method = "BH")]
  out
}
# Biotype-aware TREAT: protein-coding genes are tested against the PC effect-size
# offset (lfc_pc), lncRNAs against no floor (lfc_lnc=0) -- the same biotype split the
# mouse cross-diet arm already uses, applied uniformly to EVERY human axis (C1/C2/C3 +
# per-cohort) so the advertised "PC LFC>=0.5" floor actually binds on PROTEIN-CODING
# genes regardless of which contrast they enter through. `dt` must already carry `bt`.
add_treat_fdr_bt <- function(dt, se_col = NULL, df_col = NULL,
                             lfc_pc = SHRUNK_LFC_THR, lfc_lnc = SHRUNK_LFC_THR_LNC) {
  pc  <- add_treat_fdr(dt[bt != "lncRNA"], lfc = lfc_pc,  se_col = se_col, df_col = df_col)
  lnc <- add_treat_fdr(dt[bt == "lncRNA"], lfc = lfc_lnc, se_col = se_col, df_col = df_col)
  rbindlist(list(pc, lnc), fill = TRUE)
}

# Biotype maps retained for explicit PCG/lncRNA threshold provenance.
meta_bt0 <- fread(META, select = c("mouse_ensembl_base", "mouse_biotype"))
setnames(meta_bt0, c("mouse_ensembl_base", "mouse_biotype"), c("gene_id_mouse", "biotype"))
meta_bt0[grepl("protein_coding", biotype), biotype := "protein_coding"]
meta_bt0[grepl("lncRNA|lincRNA", biotype), biotype := "lncRNA"]
mouse_biotype_of <- setNames(meta_bt0$biotype, meta_bt0$gene_id_mouse)

ortho_bt0 <- fread(cmd = paste0("zcat ", ORTHO),
                   select = c("mouse_ensembl", "human_ensembl",
                              "confidence_tier", "is_one2one"))
ortho_bt0[, gene_id_mouse := strip_v(mouse_ensembl)]
ortho_bt0[, human_ensembl := strip_v(human_ensembl)]
ortho_bt0 <- ortho_bt0[confidence_tier %in% c("H", "M")]
ortho_bt0[, trank := match(confidence_tier, c("H", "M"))]
ortho_bt0[, o2o := is_one2one %in% c(TRUE, "True", "TRUE", "true")]
ortho_bt0[, one2one_rank := ifelse(o2o, 0L, 1L)]
setorder(ortho_bt0, human_ensembl, trank, one2one_rank, mouse_ensembl)
ortho_bt0 <- unique(ortho_bt0, by = "human_ensembl")
human_biotype_of <- setNames(mouse_biotype_of[ortho_bt0$gene_id_mouse],
                             ortho_bt0$human_ensembl)
# --- 0. Freeze the current (v7) library, then back it up ---------------------
if (file.exists(OUT)) {
  if (!file.exists(V7REF)) { file.copy(OUT, V7REF); cat("froze v7 reference ->", V7REF, "\n") }
  file.copy(OUT, BACKUP, overwrite = TRUE)
  cat("backed up current library ->", BACKUP, "\n")
}

# --- 1. Mouse cross-diet membership (TREAT FDR, consistent with human) -------
up_by_diet <- lapply(DIETS, function(d) {
  dt <- fread(file.path(PERDIET, paste0(d, "_de_results.csv")))
  dt[, gene_base := strip_v(gene)]
  dt[, bt := unname(mouse_biotype_of[gene_base])]
  dt[is.na(bt), bt := "protein_coding"]
  dt_pc <- add_treat_fdr(dt[bt != "lncRNA"], lfc = SHRUNK_LFC_THR, df_col = "df_total")
  dt_lnc <- add_treat_fdr(dt[bt == "lncRNA"], lfc = 0.0, df_col = "df_total")
  rbindlist(list(dt_pc, dt_lnc), fill = TRUE)[fdr_treat < LFSR_THR & logFC > 0, gene_base]
})
names(up_by_diet) <- DIETS
all_up <- sort(unique(unlist(up_by_diet)))
cross <- data.table(gene_id_mouse = all_up)
for (d in DIETS) cross[, (d) := gene_id_mouse %in% up_by_diet[[d]]]
cross[, n_diets_up := rowSums(.SD), .SDcols = DIETS]
cross[, diet_list := apply(.SD, 1, function(r) paste(DIETS[as.logical(r)], collapse = ",")),
      .SDcols = DIETS]

# --- 2. Human CORE: union of therapeutic-axis human DEGs ---------------------
can <- fread(CANONICAL, select = c("gene", "logFC", "padj", "symbol", "SE", "t", "P.Value"))
can[, hb := strip_v(gene)]
universe_n <- nrow(can)
col1a1_lfc <- can[symbol == "COL1A1", logFC][1]
cat(sprintf("PROVENANCE canonical-C2: universe=%d  COL1A1 logFC=%.4f\n", universe_n, col1a1_lfc))
if (universe_n > 30000)
  stop(sprintf("canonical universe=%d > 30000 (mis-pointed file). Expected C2 (~27,638).", universe_n))
if (is.na(col1a1_lfc) || col1a1_lfc < 0.5)
  stop(sprintf("COL1A1 canonical logFC=%.3f unexpected (expected ~1.3, up in disease). Check CANONICAL file.", col1a1_lfc))

can[, bt := human_biotype_of[hb]]
can[is.na(bt), bt := "protein_coding"]
can <- add_treat_fdr_bt(can, se_col = "SE")  # CORE = C1 integrated disease-vs-control ONLY: PC lfc=0.5, lncRNA 0
# CORE is the integrated disease-vs-control DEG set ALONE. The MASH-vs-MASL (C2) and
# advanced-fibrosis (C3) progression axes were dropped (2026-06-29, PI): the human arm
# is CORE integrated DEGs UNION COHORT-REPLICATED DEGs only (see spine_h below).
core_h <- unique(can[!is.na(logFC) & fdr_treat < HUMAN_PADJ_THR & logFC > 0, hb])
cat(sprintf("Human CORE (C1 MASLD-vs-Control integrated only, TREAT FDR<%.2f, PC LFC>=%.1f / lnc no-floor): %d genes\n",
            HUMAN_PADJ_THR, SHRUNK_LFC_THR, length(core_h)))
human_support <- can[hb %in% core_h, .(hb, human_symbol = symbol, human_logFC = logFC, contrast = "C1_disease_vs_ctrl")]
human_logfc <- human_support[, .(human_logFC = max(human_logFC, na.rm = TRUE)), by = hb]

# --- 2b. COHORT-REPLICATED: human DEGs significant in >=COHORT_MIN cohorts ----
cohort_lists <- lapply(COHORTS, function(c) {
  d <- fread(file.path(PERSTUDY, paste0(c, "_de_results.csv")))
  d[, hb := strip_v(gene)]
  d[, bt := human_biotype_of[hb]]
  d[is.na(bt), bt := "protein_coding"]
  
  d <- add_treat_fdr_bt(d, se_col = "SE", df_col = "df.total")   # cohort: PC at lfc=0.5, lncRNA at 0
  unique(d[!is.na(logFC) & fdr_treat < HUMAN_PADJ_THR & logFC > 0, hb])
})
cohort_n_tab <- table(unlist(cohort_lists))                  # per human gene: #cohorts significant
cohort_h <- names(cohort_n_tab)[cohort_n_tab >= COHORT_MIN]
cat(sprintf("Cohort-replicated (UP in >=%d of %d cohorts): %d human genes\n",
            COHORT_MIN, length(COHORTS), length(cohort_h)))
human_cohort_n <- data.table(hb = names(cohort_n_tab), n_cohorts_sig = as.integer(cohort_n_tab))

spine_h <- union(core_h, cohort_h)                           # full human spine (human ENSG)

# --- 3. Ortholog table (deterministic, one2one-preferred) --------------------
ortho <- fread(cmd = paste0("zcat ", ORTHO),
               select = c("mouse_ensembl", "human_ensembl", "human_symbol",
                          "confidence_tier", "is_one2one"))
ortho[, gene_id_mouse := strip_v(mouse_ensembl)]
ortho[, human_ensembl := strip_v(human_ensembl)]
ortho <- ortho[confidence_tier %in% c("H", "M")]
ortho[, trank := match(confidence_tier, c("H", "M"))]
ortho[, o2o := is_one2one %in% c(TRUE, "True", "TRUE", "true")]
ortho[, one2one_rank := ifelse(o2o, 0L, 1L)]

nho <- ortho[, .(n_human_orthologs = uniqueN(human_ensembl)), by = gene_id_mouse]

# (a) best HUMAN per mouse -- annotation map
setorder(ortho, gene_id_mouse, trank, one2one_rank, human_ensembl)
m2h <- ortho[, .(m2h_human_ensembl = human_ensembl[1],
                 m2h_human_symbol  = human_symbol[1]), by = gene_id_mouse]
# (b) best MOUSE per human ENSG -- backbone map
ortho_byhuman <- copy(ortho); setorder(ortho_byhuman, human_ensembl, trank, one2one_rank, mouse_ensembl)
ortho_byhuman <- unique(ortho_byhuman, by = "human_ensembl")
# (c) best MOUSE per human SYMBOL -- for positive-control mapping
ortho_bysym <- copy(ortho); ortho_bysym[, hsym := toupper(human_symbol)]
ortho_bysym <- ortho_bysym[hsym != "" & !is.na(hsym)]
setorder(ortho_bysym, hsym, trank, one2one_rank, mouse_ensembl)
ortho_bysym <- unique(ortho_bysym, by = "hsym")

map_h2m <- function(hset) unique(ortho_byhuman[human_ensembl %in% hset, gene_id_mouse])

# --- 4a. HUMAN SPINE: map (core U cohort) human ENSG -> mouse, 1 entry/mouse ---
core_mouse   <- map_h2m(core_h)
cohort_mouse <- map_h2m(cohort_h)
spine_mouse  <- union(core_mouse, cohort_mouse)
cat(sprintf("Human spine (mouse orthologs): core=%d cohort=%d union=%d\n",
            length(core_mouse), length(cohort_mouse), length(spine_mouse)))
# entry human (strongest logFC among the contributing therapeutic-axis contrasts)
spine_hu <- merge(human_support[hb %in% spine_h, .(human_ensembl = hb, human_symbol, human_logFC)],
                  ortho_byhuman[, .(human_ensembl, gene_id_mouse)], by = "human_ensembl")
setorder(spine_hu, gene_id_mouse, -human_logFC, human_ensembl)
spine <- spine_hu[, .(entry_human_ensembl = human_ensembl[1],
                      gene_symbol_human   = human_symbol[1]), by = gene_id_mouse]

# --- 4b. MOUSE-CONFIRMED: cross-diet UP + human directional concordance -------
mouse_confirmed_at <- function(min_diets) {
  mc_genes <- cross[n_diets_up >= min_diets, gene_id_mouse]
  mc <- merge(data.table(gene_id_mouse = mc_genes), m2h, by = "gene_id_mouse", all.x = TRUE)
  mc <- merge(mc, human_logfc, by.x = "m2h_human_ensembl", by.y = "hb", all.x = TRUE)
  mc[!is.na(human_logFC) & human_logFC > 0, gene_id_mouse]
}
mouse_confirmed <- mouse_confirmed_at(MIN_DIETS)
cat(sprintf("Mouse-confirmed tier (>=%d diets + human concordance): %d genes\n",
            MIN_DIETS, length(mouse_confirmed)))

# --- 4c. COLOC tier: canonical SuSiE PP.H4 > 0.5 -> mouse --------------------
cl <- fread(COLOCFILE)
cl[, hb := strip_v(ensembl)]
coloc_h <- unique(cl[hb != "" & !is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 > COLOC_PP4_THR, hb])
coloc_mouse <- map_h2m(coloc_h)
# best SuSiE PP4 per mouse gene (annotation)
coloc_pp4_dt <- merge(cl[hb %in% coloc_h, .(human_ensembl = hb, pp4 = coloc_best_susie_pp4)],
                      ortho_byhuman[, .(human_ensembl, gene_id_mouse)], by = "human_ensembl")
coloc_pp4_dt <- coloc_pp4_dt[, .(coloc_best_susie_pp4 = max(pp4, na.rm = TRUE)), by = gene_id_mouse]
cat(sprintf("COLOC tier (SuSiE PP.H4>%.1f): %d human -> %d mouse genes\n",
            COLOC_PP4_THR, length(coloc_h), length(coloc_mouse)))

# --- 4d. POSITIVE CONTROLS: human symbol -> mouse, folded in -----------------
# Mechanism-clean the roster (drop systemic targets) + fix mis-mapped orthologs.
msym <- fread(META, select = c("mouse_ensembl_base", "mouse_symbol_gtf"))
msym[, msu := toupper(mouse_symbol_gtf)]
msym_to_id <- function(s) { v <- msym[msu == toupper(s), mouse_ensembl_base]; if (length(v)) v[1] else NA_character_ }
pc_raw <- fread(POSCTRL)
pc_raw[, hsym := toupper(`Gene symbol`)]
pc_raw[hsym %in% names(POSCTRL_ALIASES), hsym := POSCTRL_ALIASES[hsym]]
n_excl <- sum(pc_raw$hsym %in% POSCTRL_EXCLUDE)
pc_raw <- pc_raw[!(hsym %in% POSCTRL_EXCLUDE)]
pc_map <- merge(pc_raw[, .(hsym, pos_control_direction = `Steatosis_Change_upon_KD`)],
                ortho_bysym[, .(hsym, gene_id_mouse)], by = "hsym", all.x = TRUE)
for (hs in names(POSCTRL_MOUSE_OVERRIDE)) {       # remap e.g. SCD -> Scd1
  mid <- msym_to_id(POSCTRL_MOUSE_OVERRIDE[[hs]])
  if (!is.na(mid)) pc_map[hsym == hs, gene_id_mouse := mid]
}
pc_map <- unique(pc_map[!is.na(gene_id_mouse)], by = "gene_id_mouse")
posctrl_mouse <- pc_map$gene_id_mouse
cat(sprintf("Positive controls: roster - %d excluded (%s) -> %d mapped (SCD->%s)\n",
            n_excl, paste(POSCTRL_EXCLUDE, collapse = ","), length(posctrl_mouse),
            POSCTRL_MOUSE_OVERRIDE[["SCD"]]))

# --- 4e. MASH annotation (no longer a tier; flag only) -----------------------
mash <- fread(MASH, select = c("gene", "shrunk_logFC", "lfsr"))
mash[, hb := strip_v(gene)]
mash_h  <- unique(mash[!is.na(lfsr) & lfsr < 0.05 & shrunk_logFC > 0.2, hb])
mash_mouse <- map_h2m(mash_h)

# --- 5. Assemble library: union of 5 tiers, PC+lncRNA only -------------------
lib_genes <- sort(Reduce(union, list(spine_mouse, mouse_confirmed, coloc_mouse, posctrl_mouse)))
mem <- data.table(gene_id_mouse = lib_genes)
# tier priority: positive_control > core > cohort_replicated > mouse_confirmed > coloc
mem[, tier := fifelse(gene_id_mouse %in% posctrl_mouse,    "positive_control",
              fifelse(gene_id_mouse %in% core_mouse,       "core",
              fifelse(gene_id_mouse %in% cohort_mouse,     "cohort_replicated",
              fifelse(gene_id_mouse %in% mouse_confirmed,  "mouse_confirmed",
                                                           "coloc"))))]
mem[, is_core              := gene_id_mouse %in% core_mouse]
mem[, is_cohort_replicated := gene_id_mouse %in% cohort_mouse]
mem[, has_human_de         := gene_id_mouse %in% spine_mouse]   # core OR cohort
mem[, is_mash_deg          := gene_id_mouse %in% mash_mouse]    # annotation (overlap-aware)
mem[, is_positive_control  := gene_id_mouse %in% posctrl_mouse]
cat("Library union (pre biotype filter):", nrow(mem), "genes\n")

# diet membership
mem <- merge(mem, cross[, .(gene_id_mouse, n_diets_up, diet_list)], by = "gene_id_mouse", all.x = TRUE)
mem[is.na(n_diets_up), n_diets_up := 0L]; mem[is.na(diet_list), diet_list := ""]

# entry-human (spine) + m2h (annotation)
mem <- merge(mem, spine[, .(gene_id_mouse, gene_symbol_human, entry_human_ensembl)],
             by = "gene_id_mouse", all.x = TRUE)
mem <- merge(mem, m2h, by = "gene_id_mouse", all.x = TRUE)
mem[is.na(gene_symbol_human), gene_symbol_human := m2h_human_symbol]
mem[, human_ensembl_annot := ifelse(has_human_de, entry_human_ensembl, m2h_human_ensembl)]

# COLOC pp4 + positive-control direction
mem <- merge(mem, coloc_pp4_dt, by = "gene_id_mouse", all.x = TRUE)
mem[, has_coloc := !is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 > COLOC_PP4_THR]
mem <- merge(mem, pc_map[, .(gene_id_mouse, pos_control_direction)], by = "gene_id_mouse", all.x = TRUE)

# cohort replication count (annotation, via representative human ortholog)
mem <- merge(mem, human_cohort_n, by.x = "human_ensembl_annot", by.y = "hb", all.x = TRUE)
mem[is.na(n_cohorts_sig), n_cohorts_sig := 0L]

# ortholog ambiguity + direction-conflict flags
mem <- merge(mem, nho, by = "gene_id_mouse", all.x = TRUE)
mem[is.na(n_human_orthologs), n_human_orthologs := 0L]
mem[, ortholog_ambiguous := n_human_orthologs > 1L]
mem <- merge(mem, can[, .(hb, m2h_logFC = logFC)], by.x = "m2h_human_ensembl", by.y = "hb", all.x = TRUE)
mem[, direction_conflict := has_human_de & !is.na(m2h_human_symbol) &
      m2h_human_symbol != gene_symbol_human & !is.na(m2h_logFC) & m2h_logFC < 0]

# biotype + mouse symbol
meta <- fread(META)
setnames(meta, c("mouse_ensembl_base", "mouse_symbol_gtf", "mouse_biotype"),
         c("gene_id_mouse", "gene_symbol_mouse", "biotype"), skip_absent = TRUE)
mem <- merge(mem, meta[, .(gene_id_mouse, gene_symbol_mouse, biotype)], by = "gene_id_mouse", all.x = TRUE)
mem[grepl("protein_coding", biotype), biotype := "protein_coding"]
mem[grepl("lncRNA|lincRNA", biotype), biotype := "lncRNA"]
n_before <- nrow(mem)
mem <- mem[biotype %in% KEEP_BIOTYPES]
cat("Biotype filter PC+lncRNA:", n_before, "->", nrow(mem), "genes\n")

# --- 5b. Screen-expression gate: mouse-HEPATOCYTE scRNA pseudobulk CPM (v9) ---
#         The screen is hepatocyte-autonomous, so a target can only score if its
#         mouse transcript is present in mouse HEPATOCYTES -- not a whole-liver-bulk
#         passenger from another lineage (e.g. Col1a1, stellate). Gate on mouse_hep_cpm
#         with biotype-aware floors: PC >= 1.0 CPM, lncRNA >= 0.1 CPM. The metric is
#         CPM (10x UMI counts / library-size x 1e6, hepatocyte pseudobulk) -- the
#         correct unit for droplet scRNA, NOT length-normalized TPM. Positive controls
#         and high-COLOC genes remain exempt (kept, flagged mouse_untestable). In-house
#         MCD bulk TPM (mcd_mean_tpm) is retained as an ANNOTATION column only.
mhs <- fread(MOUSEHEP)                                   # carries gene_id (GRCm39-2024-A) + gene_symbol
# Join the substrate on STABLE Ensembl gene_id (2026-06-24): survives symbol drift
# between the library roster and the scRNA atlas (e.g. Akr1b1<->Akr1b3, Atp5f1d<->
# Atp5d, Lars1<->Lars, Fyb1<->Fyb) that previously zeroed ~19 real hepatocyte genes.
# Fall back to gene_symbol for the few substrate rows lacking a ref gene_id.
mhs_id <- unique(mhs[!is.na(gene_id) & gene_id != "",
                     .(gene_id_mouse = gene_id, mouse_hep_cpm,
                       mouse_hep_ratio, mouse_hep_substrate)], by = "gene_id_mouse")
mem <- merge(mem, mhs_id, by = "gene_id_mouse", all.x = TRUE)
mhs_sym <- unique(mhs[, .(gene_symbol_mouse = gene_symbol, cpm_s = mouse_hep_cpm,
                          ratio_s = mouse_hep_ratio, sub_s = mouse_hep_substrate)],
                  by = "gene_symbol_mouse")
mem <- merge(mem, mhs_sym, by = "gene_symbol_mouse", all.x = TRUE)
mem[is.na(mouse_hep_substrate), `:=`(mouse_hep_cpm = cpm_s, mouse_hep_ratio = ratio_s,
                                     mouse_hep_substrate = sub_s)]
mem[, c("cpm_s", "ratio_s", "sub_s") := NULL]
mem[is.na(mouse_hep_cpm), mouse_hep_cpm := 0]

# bulk MCD mean TPM (v9: ANNOTATION only; the hard gate is mouse_hep_cpm above)
mcd <- fread(MCDTPM, select = c("gene_id", "mean_tpm"))
mcd[, gene_id_mouse := strip_v(gene_id)]
mcd <- mcd[, .(mcd_mean_tpm = max(mean_tpm)), by = gene_id_mouse]
mem <- merge(mem, mcd, by = "gene_id_mouse", all.x = TRUE)
mem[is.na(mcd_mean_tpm), mcd_mean_tpm := 0]

gate_fail <- (mem$biotype == "protein_coding" & mem$mouse_hep_cpm < HEP_CPM_PC_GATE) |
             (mem$biotype == "lncRNA" & mem$mouse_hep_cpm < HEP_CPM_LNC_GATE)

gate_exempt <- mem$is_positive_control |
               (mem$has_coloc & !is.na(mem$coloc_best_susie_pp4) &
                mem$coloc_best_susie_pp4 >= COLOC_GATE_EXEMPT_PP4)

drop_gene  <- gate_fail & !gate_exempt
n_exempt <- sum(gate_fail & gate_exempt)
n_pre <- nrow(mem)
mem <- mem[!drop_gene]

cat(sprintf("mouse-hep CPM gate (drop PC <%.1f CPM, lncRNA <%.1f CPM; %d exempt [ctrl|COLOC>=%.2f]): %d -> %d genes\n",
            HEP_CPM_PC_GATE, HEP_CPM_LNC_GATE, n_exempt, COLOC_GATE_EXEMPT_PP4, n_pre, nrow(mem)))

mem[, mouse_untestable := (biotype == "protein_coding" & mouse_hep_cpm < HEP_CPM_PC_GATE) |
                          (biotype == "lncRNA" & mouse_hep_cpm < HEP_CPM_LNC_GATE)]
mem[, hep_substrate_confident := !mouse_untestable]
mem[, control_benchmark_eligible := is_positive_control & !mouse_untestable]

# --- 5c. Pool-guideability gate (2026-06-24): drop genes with NO qualifying Cas13
#         guide in the upstream vM38 pool. A target with zero designable guides cannot
#         be in the screen at all -- a HARD physical constraint, so it applies to EVERY
#         tier INCLUDING positive controls and high-COLOC exemptions (the hep exemption
#         only concerns the lipid readout, not guide existence). The list is maintained
#         as a UNION across guide-build runs by build_library_guides.py
#         (coverage status == 'no_guides_in_source'); reasons: mitochondrial, wrong
#         biotype, no basic-CCDS PC transcript, or every spacer failing the off-target/
#         uniqueness filter (single-exon or large paralog families -- Sox4, Rpl36,
#         Tuba1b, Ugt1a1, ...). genome_n_guides_passed == 0 confirmed for all (not a
#         symbol/gene_id join artifact).
UNGUIDEABLE <- file.path(BASE, "Cas13_Library_Design/data/guides/unguideable_vM38.csv")
if (file.exists(UNGUIDEABLE)) {
  ung_ids <- strip_v(fread(UNGUIDEABLE)$gene_id_mouse)
  in_ung    <- mem$gene_id_mouse %in% ung_ids
  n_ctrl_g  <- sum(in_ung & mem$is_positive_control)
  n_coloc_g <- sum(in_ung & mem$has_coloc & !is.na(mem$coloc_best_susie_pp4) &
                   mem$coloc_best_susie_pp4 >= COLOC_GATE_EXEMPT_PP4)
  n_pre_g <- nrow(mem); mem <- mem[!in_ung]
  cat(sprintf("pool-guideability gate (drop no-guide genes in vM38 pool): %d -> %d  [-%d; incl %d pos-ctrl, %d COLOC>=%.2f]\n",
              n_pre_g, nrow(mem), sum(in_ung), n_ctrl_g, n_coloc_g, COLOC_GATE_EXEMPT_PP4))
} else {
  cat("pool-guideability gate: unguideable_vM38.csv not found (run build_library_guides.py first) -- no genes dropped\n")
}

# --- 6. Readout-aware soft tags (hepatocyte substrate + disease cell type) ----
if (file.exists(HEPSPEC)) {
  hep <- fread(HEPSPEC)
  mem <- merge(mem, hep[, .(gene_symbol_human = gene_symbol, hep_mean_cpm,
                            hep_ratio, hep_substrate, hep_top_other = top_other_celltype)],
               by = "gene_symbol_human", all.x = TRUE)
  mem[is.na(hep_substrate), hep_substrate := "unknown"]
} else {
  warning("hep_specificity.csv not found; run compute_hep_specificity.R first. Tagging 'unknown'.")
  mem[, `:=`(hep_mean_cpm = NA_real_, hep_ratio = NA_real_,
             hep_substrate = "unknown", hep_top_other = NA_character_)]
}

# disease-DE cell type from the multi-evidence atlas (annotation)
atlas_sc <- fread(ATLAS, select = c("ensembl_id", "sc_best_celltype"))
atlas_sc[, hb := strip_v(ensembl_id)]
mem <- merge(mem, atlas_sc[, .(hb, sc_disease_celltype = sc_best_celltype)],
             by.x = "human_ensembl_annot", by.y = "hb", all.x = TRUE)
mem[, hep_expressed_causal := has_coloc & hep_substrate == "high"]
mem[, has_human_evidence := has_human_de | has_coloc]

# evidence-axis count (for lncRNA prioritisation): human DE, mouse cross-diet,
# COLOC, positive control -- each an independent line of support.
mem[, n_evidence_axes := as.integer(has_human_de) + as.integer(n_diets_up >= MIN_DIETS) +
      as.integer(has_coloc) + as.integer(is_positive_control)]

# library arm: PC tiers = primary discovery; lncRNA = exploratory secondary (decided
# v9 -- kept in full but de-prioritised; weakest orthology + low mouse expression).
mem[, library_arm := ifelse(biotype == "lncRNA", "exploratory_lncrna", "primary")]

# lncRNA flags
mem[, lnc_human_ortholog := biotype == "lncRNA" & !is.na(gene_symbol_human) &
      gene_symbol_human != "" & !grepl("^ENSG", gene_symbol_human)]
mem[, lnc_hep_expressed := biotype == "lncRNA" & hep_substrate %in% c("high", "ambient_suspect")]
# higher-confidence exploratory lncRNA: present in mouse hepatocytes AND multi-axis.
mem[, lnc_priority := biotype == "lncRNA" &
      mouse_hep_substrate %in% c("high", "ambient_suspect") & n_evidence_axes >= 2]

mem[, library_version := LIB_VERSION]

# --- 7. Write canonical schema (v9) -----------------------------------------
out <- mem[, .(gene_id_mouse, gene_symbol_mouse, gene_symbol_human, biotype, tier,
               library_arm,
               is_core, is_cohort_replicated, n_cohorts_sig, n_diets_up, diet_list,
               has_human_de, has_coloc, coloc_best_susie_pp4, is_positive_control,
               pos_control_direction, control_benchmark_eligible, mouse_untestable,
               is_mash_deg, has_human_evidence, n_evidence_axes,
               mouse_hep_cpm, mouse_hep_ratio, mouse_hep_substrate, hep_substrate_confident,
               hep_substrate, hep_mean_cpm, hep_ratio, mcd_mean_tpm, sc_disease_celltype,
               n_human_orthologs, ortholog_ambiguous, direction_conflict,
               hep_expressed_causal, lnc_human_ortholog, lnc_hep_expressed, lnc_priority,
               library_version)]
setorder(out, tier, -n_cohorts_sig, -n_diets_up, gene_symbol_mouse)
fwrite(out, OUT)

# --- 8. Verification + reported outcomes ------------------------------------
n_total <- nrow(out)
n_pc    <- out[biotype == "protein_coding", .N]
n_lnc   <- out[biotype == "lncRNA", .N]
N_SGRNA <- 4L; N_NT <- 500L; N_EXTRA_CTRL <- 30L               # 4 gRNA/target (build_library_guides.py --n 4, 2026-06-29); NT + safe-harbor/essential
n_sgrna <- (n_pc + n_lnc) * N_SGRNA + N_NT + N_EXTRA_CTRL * N_SGRNA
HEP <- 1e7; TRANSD <- 0.30; CRE_EFF <- 0.80; FACS_EFF <- 0.60
IN_FRAC <- 0.13; GATE <- 0.15; MORTALITY <- 1/8
input_cells <- HEP * TRANSD * CRE_EFF * IN_FRAC
arm_cells   <- HEP * TRANSD * CRE_EFF * (1 - IN_FRAC) * GATE * FACS_EFF
eff_cells   <- min(input_cells, arm_cells)
mice_cov <- ceiling(n_sgrna * 500 / eff_cells)
mice     <- ceiling(mice_cov / (1 - MORTALITY))

cat("\n=== LIBRARY", LIB_VERSION, "(", OUT, ") ===\n")
cat(sprintf("TOTAL: %d  (PC=%d, lncRNA=%d)\n", n_total, n_pc, n_lnc))
cat("tier distribution (positive_control>core>cohort_replicated>mouse_confirmed>coloc):\n")
print(table(out$tier))
cat(sprintf("  is_core=%d  is_cohort_replicated=%d  has_human_de(core|cohort)=%d\n",
            sum(out$is_core), sum(out$is_cohort_replicated), sum(out$has_human_de)))
cat(sprintf("  has_coloc=%d  is_positive_control=%d  is_mash_deg(annotation)=%d\n",
            sum(out$has_coloc), sum(out$is_positive_control), sum(out$is_mash_deg)))
cat("library_arm distribution:\n"); print(table(out$library_arm))
stopifnot("miRNA leaked into the library" = out[biotype == "miRNA", .N] == 0)
stopifnot("MASH tier should not exist" = out[tier == "mash_progression", .N] == 0)
# v9 gate diagnostics: mouse-hepatocyte substrate + exemption book-keeping
cat("mouse_hep_substrate distribution (PC only; gate substrate):\n")
print(table(out[biotype == "protein_coding", mouse_hep_substrate]))
cat(sprintf("  mouse_untestable (fails mouse-hep CPM gate, kept via exemption)=%d  control_benchmark_eligible=%d\n",
            sum(out$mouse_untestable), sum(out$control_benchmark_eligible)))
cat(sprintf("  scoreable PC (mouse-hep CPM >= %.1f): %d (%.1f%% of PC)\n",
            HEP_CPM_PC_GATE,
            out[biotype == "protein_coding" & mouse_hep_cpm >= HEP_CPM_PC_GATE, .N],
            100 * out[biotype == "protein_coding" & mouse_hep_cpm >= HEP_CPM_PC_GATE, .N] / max(n_pc, 1)))
# flagship verification: HKDC1 rescued-but-flagged; Gck (mouse hepatic hexokinase) clean
hk <- out[gene_symbol_mouse == "Hkdc1"]
gk <- out[gene_symbol_mouse == "Gck"]
cat(sprintf("  CHECK Hkdc1: present=%s tier=%s mouse_hep_substrate=%s mouse_hep_cpm=%.3f mouse_untestable=%s\n",
            nrow(hk) > 0, if (nrow(hk)) hk$tier else "NA",
            if (nrow(hk)) hk$mouse_hep_substrate else "NA",
            if (nrow(hk)) hk$mouse_hep_cpm else NA_real_,
            if (nrow(hk)) hk$mouse_untestable else "NA"))
cat(sprintf("  CHECK Gck:   present=%s mouse_hep_substrate=%s\n",
            nrow(gk) > 0, if (nrow(gk)) gk$mouse_hep_substrate else "NA"))
if (nrow(hk) == 0)
  warning("HKDC1 absent from v9 library -- expected present via COLOC exemption (PP4 0.992 >= ",
          COLOC_GATE_EXEMPT_PP4, "). Check coloc tier / exemption threshold.")
cat("human hep_substrate distribution (soft tag):\n"); print(table(out$hep_substrate))
cat("disease-DE cell type (top 6):\n"); print(head(sort(table(out$sc_disease_celltype), decreasing = TRUE), 6))
cat(sprintf("lncRNA: %d total | human-orthologous=%d | hep-expressed=%d\n",
            n_lnc, out[lnc_human_ortholog == TRUE, .N], out[lnc_hep_expressed == TRUE, .N]))
cat(sprintf("sgRNAs: %d (%d/target, +%d NT, +%d safe-harbor/essential x%d)\n",
            n_sgrna, N_SGRNA, N_NT, N_EXTRA_CTRL, N_SGRNA))
cat(sprintf("mice: %d surviving @500x; INJECT %d (1/8 mortality)\n", mice_cov, mice))

# COHORT_MIN sensitivity (FULL library AFTER the v9 mouse-hep CPM gate)
cat("\nCOHORT_MIN sensitivity (FULL library after mouse-hep CPM gate):\n")
pc_set  <- meta[grepl("protein_coding", biotype), gene_id_mouse]
lnc_set <- meta[grepl("lncRNA|lincRNA", biotype), gene_id_mouse]
id2sym    <- setNames(meta$gene_symbol_mouse, meta$gene_id_mouse)
# gene_id-keyed mouse-hep CPM lookup (matches the actual gate)
hep_cpm_of <- setNames(mhs_id$mouse_hep_cpm, mhs_id$gene_id_mouse)
fail_of   <- function(ids) {
  cpm <- unname(hep_cpm_of[ids])
  cpm[is.na(cpm)] <- 0.0
  bt <- unname(mouse_biotype_of[ids])
  bt[is.na(bt)] <- "protein_coding"
  (bt == "protein_coding" & cpm < HEP_CPM_PC_GATE) | (bt == "lncRNA" & cpm < HEP_CPM_LNC_GATE)
}
coloc_pp4_of <- setNames(coloc_pp4_dt$coloc_best_susie_pp4, coloc_pp4_dt$gene_id_mouse)
exempt_pc <- function(ids) (ids %in% posctrl_mouse) |
             (!is.na(coloc_pp4_of[ids]) & coloc_pp4_of[ids] >= COLOC_GATE_EXEMPT_PP4)
gated_n <- function(g) {
  pc <- intersect(g, pc_set); lnc <- intersect(g, lnc_set)
  pc_ok <- pc[!fail_of(pc) | exempt_pc(pc)]
  lnc_ok <- lnc[!fail_of(lnc) | exempt_pc(lnc)]
  length(union(pc_ok, lnc_ok))
}
for (cm in c(2L, 3L, 4L)) {
  sm <- union(core_mouse, map_h2m(names(cohort_n_tab)[cohort_n_tab >= cm]))
  full <- Reduce(union, list(sm, mouse_confirmed, coloc_mouse, posctrl_mouse))
  cat(sprintf("  >=%d cohorts: FULL (mouse-hep CPM-gated) = %d\n", cm, gated_n(full)))
}

# --- 9. v8 -> v9 diff --------------------------------------------------------
diffref <- if (file.exists(V7REF)) V7REF else BACKUP
if (file.exists(diffref)) {
  v7 <- fread(diffref); v7g <- v7$gene_id_mouse; v8g <- out$gene_id_mouse
  added   <- setdiff(v8g, v7g); dropped <- setdiff(v7g, v8g)
  diff <- rbind(data.table(gene_id_mouse = added,   change = "added"),
                data.table(gene_id_mouse = dropped, change = "dropped"))
  diff <- merge(diff, out[, .(gene_id_mouse, gene_symbol_mouse, gene_symbol_human, tier, biotype)],
                by = "gene_id_mouse", all.x = TRUE)
  fwrite(diff, DIFFOUT)
  cat(sprintf("\nv8->v9 diff: +%d added / -%d dropped (v8=%d, v9=%d) -> %s\n",
              length(added), length(dropped), length(v7g), length(v8g), DIFFOUT))
}

# --- 10. BUILD_MANIFEST ------------------------------------------------------
md5 <- function(p) if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_
git_sha <- tryCatch(system2("git", c("-C", BASE, "rev-parse", "--short", "HEAD"),
                            stdout = TRUE, stderr = FALSE), error = function(e) "NA")
manifest <- c(
  sprintf("library_version: %s", LIB_VERSION),
  sprintf("built_utc: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("git_sha: %s", paste(git_sha, collapse = "")),
  "human_deg_source: CORE = C1 MASLD-vs-Control integrated DEGs only, UNION per-study cohort-replicated DEGs (C2 MASH-vs-MASL + C3 fibrosis axes dropped 2026-06-29)",
  sprintf("cohort_min: %d of %d  |  human TREAT FDR<%.2f (PC LFC>=%.1f biotype-split on C1 core + per-cohort axes; lncRNA no-LFC, UP)",
          COHORT_MIN, length(COHORTS), HUMAN_PADJ_THR, SHRUNK_LFC_THR),
  sprintf("mouse_confirmed: >=%d diets  |  TREAT FDR<%.2f per diet (PC LFC>=%.1f; lncRNA lfc=0.0), UP",
          MIN_DIETS, LFSR_THR, SHRUNK_LFC_THR),
  sprintf("pc_gate: mouse-HEPATOCYTE scRNA pseudobulk CPM (drop PC <%.1f CPM, lncRNA <%.1f CPM; positive controls + COLOC SuSiE>=%.2f exempt, flagged mouse_untestable). MCD bulk TPM = annotation only.", HEP_CPM_PC_GATE, HEP_CPM_LNC_GATE, COLOC_GATE_EXEMPT_PP4),
  sprintf("mouse_untestable(fails mouse-hep CPM gate, kept via exemption): %d | control_benchmark_eligible: %d",
          sum(out$mouse_untestable), sum(out$control_benchmark_eligible)),
  sprintf("library_arm: primary=%d exploratory_lncrna=%d",
          out[library_arm == "primary", .N], out[library_arm == "exploratory_lncrna", .N]),
  sprintf("n_genes: %d (PC=%d lncRNA=%d)", n_total, n_pc, n_lnc),
  sprintf("tier: core=%d cohort_replicated=%d mouse_confirmed=%d coloc=%d positive_control=%d",
          out[tier == "core", .N], out[tier == "cohort_replicated", .N],
          out[tier == "mouse_confirmed", .N], out[tier == "coloc", .N],
          out[tier == "positive_control", .N]),
  sprintf("is_positive_control: %d | has_coloc: %d | is_mash_deg(annotation): %d",
          sum(out$is_positive_control), sum(out$has_coloc), sum(out$is_mash_deg)),
  sprintf("n_sgrna: %d  mice_surviving_500x: %d  mice_to_inject: %d", n_sgrna, mice_cov, mice),
  "inputs (path | mtime | md5):",
  sprintf("  human_canonical_C2: %s | %s | %s", CANONICAL, format(file.mtime(CANONICAL)), md5(CANONICAL)),
  sprintf("  per_study_dir: %s | cohorts=%s", PERSTUDY, paste(COHORTS, collapse = ",")),
  sprintf("  coloc: %s | %s | %s", COLOCFILE, format(file.mtime(COLOCFILE)), md5(COLOCFILE)),
  sprintf("  pos_control: %s | %s | %s", POSCTRL, format(file.mtime(POSCTRL)), md5(POSCTRL)),
  sprintf("  ortho: %s | %s | %s", ORTHO, format(file.mtime(ORTHO)), md5(ORTHO)),
  sprintf("  atlas: %s | %s | %s", ATLAS, format(file.mtime(ATLAS)), md5(ATLAS)),
  sprintf("  hepspec_human: %s | %s | %s", HEPSPEC, format(file.mtime(HEPSPEC)), md5(HEPSPEC)),
  sprintf("  mouse_hep_gate: %s | %s | %s", MOUSEHEP, format(file.mtime(MOUSEHEP)), md5(MOUSEHEP)),
  sprintf("  mcd_tpm_annot: %s | %s | %s", MCDTPM, format(file.mtime(MCDTPM)), md5(MCDTPM)),
  sprintf("  canonical_universe: %d  COL1A1_logFC: %.4f", universe_n, col1a1_lfc),
  paste0("  per_diet: ", paste(DIETS, collapse = ","))
)
writeLines(manifest, MANIFEST)
cat("wrote", MANIFEST, "\n")
