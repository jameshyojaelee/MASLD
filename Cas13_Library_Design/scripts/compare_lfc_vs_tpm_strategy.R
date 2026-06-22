#!/usr/bin/env Rscript
# compare_lfc_vs_tpm_strategy.R
# ---------------------------------------------------------------------------
# For a ~2,000-gene Cas13 KD screen, compare two ways to reach that size:
#   STRATEGY A : lower human ashr LFC cutoff + biotype-aware mouse disease-liver
#                TPM gate (protein_coding >=1.0, lncRNA >=0.5 -- canonical floors)
#   STRATEGY B : higher human ashr LFC cutoff, NO TPM gate
#
# Cas13 knockdown requires the transcript to be expressed in the target tissue,
# so "scoreable fraction" (passes the split TPM gate) is a hard prerequisite.
# This script quantifies what each strategy trades:
#   - scoreable fraction (TPM>=1)             -> live vs dead library slots
#   - mouse disease-liver TPM distribution
#   - human disease association (mean shrunk LFC)
#   - OpenTargets MASLD recovery + enrichment -> disease relevance
#   - corr(human LFC, mouse TPM)              -> does high-LFC pick low-TPM genes?
#
# Mouse disease-liver TPM = row-median across ALL disease samples pooled from
# every featurecounts source (per-source TPM normalisation, then pooled median).
# Library = human spine (ashr>thr) UNION mouse-confirmed tier, biotype PC/lncRNA.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

FC_DIR     <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/counts/featurecounts")
PUB_FC     <- file.path(BASE, "RNA-seq/Mouse/Public_Diet_Models/counts/featurecounts/gene_counts.txt")
WD_DIR     <- file.path(BASE, "RNA-seq/Mouse/Western_Diet_Datasets")
MAIN_META  <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/metadata/unified_mouse_metadata.csv")
ASHR       <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                        "results/integration/meta_results_ashr.csv")  # v6: metafor ashr (was dream)
ORTHO      <- file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz")
MOUSE_META <- file.path(BASE, "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv")
PERDIET    <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet_cas13")  # Cas13 library Western pool; decoupled from paper 4-model per_diet (2026-06-16)
OT_FILE    <- file.path(BASE, "data/published_gene_panels/opentargets_masld_2025.tsv")
DATDIR     <- file.path(BASE, "Cas13_Library_Design/data")

DIETS         <- c("MCD", "CDAHFD", "Western", "HFD")
KEEP_BIOTYPES <- c("protein_coding", "lncRNA")
LFSR_THR      <- 0.05
SHRUNK_MOUSE  <- 0.5
PC_GATE       <- 1.0   # protein-coding TPM floor (canonical)
LNC_GATE      <- 0.5   # lncRNA TPM floor (canonical; lncRNAs run lower)
TARGET        <- 2000
# Datasets NOT in the mouse DEG analysis (excluded in config/mouse_datasets.yaml:
# no matched controls / problematic) -> must NOT feed the disease-liver TPM reference.
EXCLUDE_DATASETS <- c("GSE159911",  # LIDPAD (excluded 2026-05-28)
                      "GSE225616",  # no matched controls (excluded 2026-05-27)
                      "GSE263273")  # AMLN_ob, ob/ob, no matched controls (excluded 2026-05-27)
# TPM reference scope: "all" = pooled across the 4 DEG diet groups (default);
# "HFD" = HFD diet alone (GSE224069 + GSE274914 + GSE246088 HFD arm).
TPM_REF <- Sys.getenv("TPM_REF", "all")
SUF     <- if (TPM_REF == "HFD") "_hfd" else ""
strip_v <- function(x) sub("[.][0-9]+$", "", x)

# ============================== 1. selection inputs =========================
ash <- fread(ASHR, select = c("gene", "symbol", "shrunk_logFC", "lfsr", "logFC"))
ash[, hb := strip_v(gene)]
human_logfc <- ash[, .(hb, hlfc = logFC)]

ortho <- fread(cmd = paste0("zcat ", ORTHO),
               select = c("mouse_ensembl", "human_ensembl", "confidence_tier", "is_one2one"))
ortho[, mouse_ensembl := strip_v(mouse_ensembl)]; ortho[, human_ensembl := strip_v(human_ensembl)]
ortho <- ortho[confidence_tier %in% c("H", "M")]
ortho[, trank := match(confidence_tier, c("H", "M"))]
ortho[, o2o  := fifelse(is_one2one %in% c(TRUE, "True", "TRUE", "true"), 0L, 1L)]
setorder(ortho, human_ensembl, trank, o2o, mouse_ensembl)
ortho_byhuman <- unique(ortho, by = "human_ensembl")
setorder(ortho, mouse_ensembl, trank, o2o, human_ensembl)
ortho_m2h <- unique(ortho, by = "mouse_ensembl")

mm <- fread(MOUSE_META)
mm[grepl("protein_coding", mouse_biotype), bt2 := "protein_coding"]
mm[grepl("lncRNA|lincRNA",  mouse_biotype), bt2 := "lncRNA"]
pc_lnc_ids <- mm[bt2 %in% KEEP_BIOTYPES, mouse_ensembl_base]
bt_of <- setNames(mm$bt2, mm$mouse_ensembl_base)   # mouse ensembl -> {protein_coding, lncRNA}

# human symbol per mouse-ensembl (for OpenTargets membership)
sym_by_mouse <- merge(ortho_m2h[, .(mgb = mouse_ensembl, hb = human_ensembl)],
                      ash[, .(hb, symbol, hlfc = logFC, h_shrunk = shrunk_logFC)],
                      by = "hb", all.x = TRUE)
setkey(sym_by_mouse, mgb)

# mouse-confirmed tier (constant)
up_by_diet <- lapply(DIETS, function(d) {
  dt <- fread(file.path(PERDIET, paste0(d, "_de_results.csv")),
              select = c("gene", "shrunk_logFC", "lfsr"))
  dt[, gb := strip_v(gene)]
  dt[lfsr < LFSR_THR & shrunk_logFC > SHRUNK_MOUSE, gb]
})
mouse_any   <- unique(unlist(up_by_diet))
n_up        <- sapply(mouse_any, function(g) sum(sapply(up_by_diet, function(u) g %in% u)))
mouse_3plus <- mouse_any[n_up >= 3]
mc_map <- merge(ortho_m2h[mouse_ensembl %in% mouse_3plus, .(mouse_ensembl, hb = human_ensembl)],
                human_logfc, by = "hb", all.x = TRUE)
mouse_conf <- intersect(mc_map[!is.na(hlfc) & hlfc > 0, mouse_ensembl], pc_lnc_ids)

# GWAS finemap tier: nearest-gene of recommended_pip>0.5 finemapped variants (constant)
.gpc  <- fread(file.path(BASE, "Cas13_Library_Design/data/candidates_pc_independent.csv"))[axis_finemap == 1]
.glnc <- fread(file.path(BASE, "Cas13_Library_Design/data/candidates_lncrna_independent.csv"))[axis_finemap == 1]
gwas_finemap <- intersect(unique(strip_v(c(.gpc$gene_id_mouse, .glnc$gene_id_mouse))), pc_lnc_ids)

get_spine <- function(thr)
  intersect(ortho_byhuman[human_ensembl %in% ash[!is.na(lfsr) & lfsr < LFSR_THR & shrunk_logFC > thr, hb],
                          mouse_ensembl], pc_lnc_ids)
get_lib <- function(thr) sort(union(union(get_spine(thr), mouse_conf), gwas_finemap))  # + GWAS tier

# ============================== 2. mouse disease-liver TPM (pooled) =========
load_fc <- function(path) {
  dt <- fread(path, skip = "Geneid")
  sc <- setdiff(names(dt), c("Geneid","Chr","Start","End","Strand","Length"))
  len <- setNames(dt[["Length"]], strip_v(dt[["Geneid"]]))
  mat <- as.matrix(dt[, ..sc]); rownames(mat) <- strip_v(dt[["Geneid"]])
  colnames(mat) <- basename(dirname(sc)); list(counts = mat, lengths = len)
}
tpm_of <- function(counts, len) {
  g <- intersect(rownames(counts), names(len)); counts <- counts[g,,drop=FALSE]
  rpk <- sweep(counts, 1, len[g]/1e3, "/"); sweep(rpk, 2, colSums(rpk)/1e6, "/")
}
meta <- fread(MAIN_META)
dz   <- if ("group_binary" %in% names(meta)) meta[group_binary != "Control"] else
          meta[!grepl("Control|LFD|Chow|Ctrl", condition)]
dz <- dz[!dataset %in% EXCLUDE_DATASETS]   # drop non-DEG datasets (LIDPAD etc.)
cat("Disease datasets feeding TPM (unified):\n"); print(dz[, .N, by = dataset][order(dataset)])

message("Pooling mouse disease-liver TPM ...")
mcd_fcs <- lapply(file.path(FC_DIR, c("gene_counts_inhouse.txt","gene_counts_gse156918.txt",
                                       "gene_counts_gse205974.txt")), load_fc)
ref_len <- mcd_fcs[[1]]$lengths
mcd_tpm <- tpm_of(do.call(cbind, lapply(mcd_fcs, `[[`, "counts")), ref_len)
pub_tpm <- tpm_of(load_fc(PUB_FC)$counts, ref_len)
g220 <- load_fc(file.path(WD_DIR,"GSE220575/counts/featurecounts/gene_counts.txt")); g220t <- tpm_of(g220$counts, ref_len)
g246 <- load_fc(file.path(WD_DIR,"GSE246088/counts/featurecounts/gene_counts.txt")); g246t <- tpm_of(g246$counts, ref_len)
g305 <- load_fc(file.path(WD_DIR,"GSE305484/counts/featurecounts/gene_counts.txt")); g305t <- tpm_of(g305$counts, ref_len)
g328 <- load_fc(file.path(WD_DIR,"GSE246328/counts/featurecounts/gene_counts.txt")); g328t <- tpm_of(g328$counts, ref_len)  # GAN (NASH)
g565 <- load_fc(file.path(WD_DIR,"GSE292565/counts/featurecounts/gene_counts.txt")); g565t <- tpm_of(g565$counts, ref_len)  # FFC (NASH)

g220m <- merge(fread(file.path(WD_DIR,"GSE220575/metadata/sample_metadata.csv")),
               fread(file.path(WD_DIR,"GSE220575/metadata/gsm_to_srr.tsv")), by="gsm")[condition %in% c("MASH","HCC"), srr]
g246mm <- merge(fread(file.path(WD_DIR,"GSE246088/metadata/sample_metadata.csv")),
                fread(file.path(WD_DIR,"GSE246088/metadata/gsm_to_srr.tsv")), by.x="geo_accession", by.y="gsm")
g246d    <- g246mm[genotype=="Plvap_Control" & diet %in% c("Western_Diet","High_Fat_Diet"), srr]
g246_hfd <- g246mm[genotype=="Plvap_Control" & diet == "High_Fat_Diet", srr]   # HFD arm only
g305d  <- merge(fread(file.path(WD_DIR,"GSE305484/metadata/sample_metadata.csv")),
                fread(file.path(WD_DIR,"GSE305484/metadata/gsm_to_srr.tsv")), by.x="gsm_accession", by.y="gsm")[condition!="Chow", srr]
g328d  <- fread(file.path(WD_DIR,"GSE246328/metadata/sample_metadata.csv"))[condition=="Disease", srr]  # GAN (excl. Reversal)
g565d  <- fread(file.path(WD_DIR,"GSE292565/metadata/sample_metadata.csv"))[condition=="Disease", srr]  # FFC

dz_main <- dz$sample_id
src <- list(
  mcd_tpm[, intersect(colnames(mcd_tpm), dz_main),  drop=FALSE],
  pub_tpm[, intersect(colnames(pub_tpm), dz_main),  drop=FALSE],
  g220t[,   intersect(colnames(g220t),   g220m),    drop=FALSE],
  g246t[,   intersect(colnames(g246t),   g246d),    drop=FALSE],
  g305t[,   intersect(colnames(g305t),   g305d),    drop=FALSE],
  g328t[,   intersect(colnames(g328t),   g328d),    drop=FALSE],
  g565t[,   intersect(colnames(g565t),   g565d),    drop=FALSE])
if (TPM_REF == "HFD") {                       # HFD diet alone as TPM reference
  hfd_ids <- dz[diet_model == "HFD", sample_id]
  src <- list(
    pub_tpm[, intersect(colnames(pub_tpm), hfd_ids),  drop=FALSE],   # GSE224069 + GSE274914
    g246t[,   intersect(colnames(g246t),   g246_hfd), drop=FALSE])   # GSE246088 HFD arm
  cat("TPM reference = HFD diet ALONE (GSE224069 + GSE274914 + GSE246088 HFD arm)\n")
}
allg <- Reduce(intersect, lapply(src, rownames))
pooled <- do.call(cbind, lapply(src, function(m) m[allg,,drop=FALSE]))
tpm_vec <- setNames(apply(pooled, 1, median), allg)
cat(sprintf("Pooled disease samples: %d ; genes with TPM: %d\n", ncol(pooled), length(tpm_vec)))

# biotype-aware scoreability gate: PC TPM>=1.0, lncRNA TPM>=0.5 (canonical floors)
scoreable <- function(genes) {
  g   <- genes[genes %in% names(tpm_vec)]
  thr <- ifelse(bt_of[g] == "lncRNA", LNC_GATE, PC_GATE)
  g[ tpm_vec[g] >= thr ]                              # returns the genes that pass
}

# ============================== 3. OpenTargets ==============================
ot <- fread(OT_FILE, skip = "gene_symbol"); ot_genes <- unique(ot$gene_symbol)
N_total <- length(ot_genes)

# ============================== 4. size sweep ===============================
GRID <- c(0.10,0.15,0.20,0.30,0.40,0.50,0.60,0.75,1.00)
sweep_dt <- rbindlist(lapply(GRID, function(thr) {
  lib  <- get_lib(thr)
  pres <- intersect(lib, names(tpm_vec))
  ns   <- length(scoreable(pres))                    # split-gate scoreable count
  data.table(thr = thr,
             n_lib = length(lib),
             n_lib_tpm1 = ns,
             pct_tpm1 = round(100 * ns / length(pres), 1))
}))
cat("\n=== library size by human ashr LFC cutoff ===\n"); print(sweep_dt)

# pick thr_B (no TPM) closest to TARGET, thr_A (with TPM) closest to TARGET
thr_B <- GRID[which.min(abs(sweep_dt$n_lib      - TARGET))]
thr_A <- GRID[which.min(abs(sweep_dt$n_lib_tpm1 - TARGET))]
cat(sprintf("\nChosen: STRATEGY A thr=%.2f (+ split TPM gate PC>=%.1f/lncRNA>=%.1f)  |  STRATEGY B thr=%.2f (no TPM)\n",
            thr_A, PC_GATE, LNC_GATE, thr_B))

# ============================== 5. head-to-head =============================
characterise <- function(genes, label) {
  g    <- intersect(genes, names(tpm_vec))
  tg   <- tpm_vec[g]
  info <- sym_by_mouse[.(genes)]
  syms <- unique(na.omit(info$symbol))
  n_ot <- sum(syms %in% ot_genes)
  base_ot <- N_total / length(unique(na.omit(sym_by_mouse$symbol)))   # prevalence in ortholog universe
  data.table(
    strategy        = label,
    n_genes         = length(genes),
    pct_scoreable   = round(100 * length(scoreable(g)) / length(g), 1),  # split-gate scoreable
    median_TPM       = round(median(tg, na.rm = TRUE), 2),
    mean_human_LFC   = round(mean(info$h_shrunk, na.rm = TRUE), 3),
    median_human_LFC = round(median(info$h_shrunk, na.rm = TRUE), 3),
    n_OpenTargets   = n_ot,
    OT_recovery_pct = round(100*n_ot/N_total, 1),
    OT_enrichment   = round((n_ot/length(syms)) / base_ot, 2))
}
libA <- scoreable(get_lib(thr_A))   # split-gate (PC>=1.0 / lncRNA>=0.5)
libB <- get_lib(thr_B)

cmp <- rbindlist(list(
  characterise(libA, sprintf("A: ashr>%.2f + TPM gate", thr_A)),
  characterise(libB, sprintf("B: ashr>%.2f, no TPM", thr_B))))
cat("\n=== HEAD-TO-HEAD (both ~", TARGET, "genes) ===\n", sep=""); print(cmp)

ov <- length(intersect(libA, libB))
cat(sprintf("\nOverlap A vs B: %d genes (%.0f%% of A, %.0f%% of B)\n",
            ov, 100*ov/length(libA), 100*ov/length(libB)))

# does a higher LFC cutoff pull in lower-TPM genes?
spine_all <- get_spine(0.10)
sp <- data.table(mgb = spine_all)
sp <- merge(sp, sym_by_mouse[, .(mgb, h_shrunk)], by = "mgb", all.x = TRUE)
sp[, tpm := tpm_vec[mgb]]
sp[, gate := ifelse(bt_of[mgb] == "lncRNA", LNC_GATE, PC_GATE)]
sp[, scor := !is.na(tpm) & tpm >= gate]
rho_lfc_tpm <- cor(sp$h_shrunk, log2(sp$tpm + 1), method = "spearman", use = "complete.obs")
cat(sprintf("\nSpearman corr(human shrunk LFC, log2 mouse disease TPM) over spine: %.3f\n", rho_lfc_tpm))
hi <- sp[h_shrunk > thr_B]; lo <- sp[h_shrunk > thr_A & h_shrunk <= thr_B]
cat(sprintf("  scoreable(split gate): LFC>%.2f tier = %.1f%%  vs  %.2f<LFC<=%.2f tier = %.1f%%\n",
            thr_B, 100*mean(hi$scor, na.rm=TRUE), thr_A, thr_B, 100*mean(lo$scor, na.rm=TRUE)))

fwrite(sweep_dt, file.path(DATDIR, paste0("lfc_vs_tpm_size_sweep",  SUF, ".csv")))
fwrite(cmp,      file.path(DATDIR, paste0("lfc_vs_tpm_headtohead", SUF, ".csv")))
cat(sprintf("\nWrote: lfc_vs_tpm_size_sweep%s.csv , lfc_vs_tpm_headtohead%s.csv\n", SUF, SUF))
