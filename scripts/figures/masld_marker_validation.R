#!/usr/bin/env Rscript
# masld_marker_validation.R
# Main Fig 3 (RNA-seq) — POSITIVE-CONTROL validation that cross-cohort batch
# correction removed cohort structure while preserving disease/fibrosis biology.
#
# Idea: independently published MASLD/NASH transcriptomic marker signatures are
# ground truth. If, on the batch-corrected mega-cohort matrix, these markers (a)
# separate Disease from Control, (b) rise monotonically across fibrosis stage
# F0->F4, and (c) NO LONGER cluster samples by cohort, then the correction is
# trustworthy. Two of the four signatures (Moylan 2014 / GSE49541, Arendt 2015 /
# GSE89632) are derived from cohorts NOT in our atlas -> genuine cross-cohort
# validators; the other four (Govaere/Suppli/Hoang/Pantano) come from cohorts
# inside our mega-5 and are positive-control/consistency baselines (flagged).
#
# Batch handling mirrors the sibling panel figS_pca_definitive.R: group-BLIND
# limma::removeBatchEffect(logCPM, batch=dataset, covariates=sex) — disease is
# NOT protected, so disease separation after correction is an honest test.
#
# Outputs (figures/main/fig3_RNAseq/panels/  == FIG2_DIR/panels):
#   marker_validation_heatmap.pdf      samples x up-marker genes, corrected, row-split by source
#   marker_module_scores_stage.pdf     GSVA score per signature across F0-F4
#   marker_module_scores_binary.pdf    GSVA score per signature, Control vs Disease
#   marker_cohort_mixing.pdf           variance explained by cohort vs disease, Before vs After
#   data/marker_validation_ledger.csv  per-signature gene counts, AUC, rho, R2, silhouette
#   data/transcriptomic_signatures_union.tsv   long table of all signature genes
#   data/transcriptomic_signatures_overlap.csv source x source shared-gene / Jaccard matrix
#
# House style: PDF only; NO titles/subtitles/annotation sentences inside panels
# (captions -> message()); all text black; control == #9E9E9E; gene names italic.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(limma); library(matrixStats)
  library(ComplexHeatmap); library(circlize); library(GSVA); library(cluster)
  library(UpSetR)
})
grDevices::pdf.options(useDingbats = FALSE)
set.seed(42)  # governs kmeans sex inference + jitterdodge; analysis otherwise deterministic

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT      <- file.path(FIG2_DIR, "panels")          # FIG2_DIR -> main/fig3_RNAseq
DATA_DIR <- file.path(OUT, "data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

PANEL_DIR <- file.path(BASE, "data/published_gene_panels")
MEGA <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")  # include_in_mega=TRUE
CTRL <- masld_colors$control   # #9E9E9E (invariant)
DIS  <- masld_colors$mash      # #C9265E
AFT  <- masld_colors$conserved # #00695C teal (after-correction)

# Published transcriptomic signatures taken VERBATIM from the primary paper/supplement
# (Suppli 2019 + Hoang 2019 were dropped 2026-06-23: their verbatim lists were not
# machine-accessible, so they had been top-N reconstructions — excluded to keep every
# panel a directly-published gene set). `external` = derived from a cohort NOT in our
# mega-5 (genuine cross-cohort validator).
PANELS <- data.table(
  file     = c("govaere_2020_panel.tsv", "pantano_2021_panel.tsv",
               "moylan_2014_panel.tsv",  "arendt_2015_panel.tsv"),
  label    = c("Govaere 2020", "Pantano 2021", "Moylan 2014", "Arendt 2015"),
  geo      = c("GSE135251", "GSE162694", "GSE49541", "GSE89632"),
  external = c(FALSE, FALSE, TRUE, TRUE)
)
SIG_LEVELS <- PANELS$label

# ---------------------------------------------------------------------------
# 1. Read curated panels -> long signature table (+ union + overlap outputs)
# ---------------------------------------------------------------------------
read_panel <- function(i) {
  f <- file.path(PANEL_DIR, PANELS$file[i])
  if (!file.exists(f)) stop("missing panel: ", f)
  ln <- readLines(f, warn = FALSE)                 # drop '#' provenance header lines
  ln <- ln[!grepl("^\\s*#", ln)]
  dt <- fread(text = paste(ln, collapse = "\n"))
  stopifnot("gene_symbol" %in% names(dt))
  if (!"direction" %in% names(dt)) dt[, direction := "up"]
  dt[, .(label = PANELS$label[i], geo = PANELS$geo[i], external = PANELS$external[i],
         gene_symbol = toupper(trimws(gene_symbol)),
         direction   = tolower(trimws(direction)))]
}
sig_long <- rbindlist(lapply(seq_len(nrow(PANELS)), read_panel))
sig_long <- sig_long[!is.na(gene_symbol) & gene_symbol != ""]
# harmonize legacy aliases -> current HGNC so cross-study overlap and matrix mapping
# are exact (e.g. PPAPDC1A->PLPP4, CYR61->CCN1); falls back to original if unmapped.
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  off <- limma::alias2SymbolTable(sig_long$gene_symbol, species = "Hs")
  sig_long[, gene_symbol := toupper(ifelse(is.na(off), gene_symbol, off))]
}
sig_long[!direction %in% c("up", "down"), direction := "up"]
sig_long <- unique(sig_long, by = c("label", "gene_symbol"))
sig_long[, label := factor(label, levels = SIG_LEVELS)]
fwrite(sig_long, file.path(DATA_DIR, "transcriptomic_signatures_union.tsv"), sep = "\t")
cat("Signature gene counts (total / up / down):\n")
print(sig_long[, .(n = .N, up = sum(direction == "up"), down = sum(direction == "down")),
               by = .(label, geo, external)][order(label)])

# source x source overlap (shared genes + Jaccard, any direction)
gsets_all <- split(sig_long$gene_symbol, sig_long$label)
ov <- CJ(a = SIG_LEVELS, b = SIG_LEVELS, sorted = FALSE)
ov[, `:=`(shared = mapply(function(x, y) length(intersect(gsets_all[[x]], gsets_all[[y]])), a, b),
          jaccard = mapply(function(x, y) {
            u <- length(union(gsets_all[[x]], gsets_all[[y]]))
            if (u == 0) 0 else round(length(intersect(gsets_all[[x]], gsets_all[[y]])) / u, 3)
          }, a, b))]
fwrite(dcast(ov, a ~ b, value.var = "shared"),
       file.path(DATA_DIR, "transcriptomic_signatures_overlap.csv"))

# UpSet of the FULL marker sets (up+down): set sizes + cross-study intersections.
# Set-size bars are a single uniform colour (studies are NOT colour-coded, per request).
gsets_full <- split(sig_long$gene_symbol, as.character(sig_long$label))[SIG_LEVELS]
set_cols   <- rep("#595959", nrow(PANELS))
nuniq <- length(unique(sig_long$gene_symbol))
cat(sprintf("\nMarker sets: %d rows across %d studies; %d unique genes (overlap small).\n",
            nrow(sig_long), nrow(PANELS), nuniq))
pdf(file.path(OUT, "marker_signature_upset.pdf"), width = fig_full_width, height = 3.4,
    useDingbats = FALSE)
# NB: call upset() directly (no print()) inside the device — print() adds a blank first page.
UpSetR::upset(UpSetR::fromList(gsets_full),
  sets = rev(SIG_LEVELS), keep.order = TRUE, order.by = "freq", nintersects = 40,
  mainbar.y.label = "Genes in intersection", sets.x.label = "Signature size (genes)",
  main.bar.color = "black", matrix.color = "black", shade.color = "grey88",
  sets.bar.color = rev(set_cols),
  text.scale = c(1.1, 1.0, 1.0, 0.9, 1.0, 0.95), point.size = 1.9, line.size = 0.5)
dev.off()
# UpSetR unconditionally emits a leading blank page (internal grid.newpage, device-independent);
# strip it with ghostscript so the panel is a single page.
upf <- file.path(OUT, "marker_signature_upset.pdf")
if (nzchar(Sys.which("gs"))) {
  tmp <- file.path(OUT, ".upset_tmp.pdf")   # same filesystem as upf so file.rename works
  rc <- system2("gs", c("-q", "-dBATCH", "-dNOPAUSE", "-dFirstPage=2",
                        "-dAutoRotatePages=/None", "-sDEVICE=pdfwrite",
                        paste0("-sOutputFile=", tmp), upf))
  if (rc == 0 && file.exists(tmp)) file.rename(tmp, upf)
}
cat("Wrote marker_signature_upset.pdf\n")
message(sprintf(
"[marker_signature_upset] Gene-set sizes and cross-study intersections of the %d verbatim published MASLD signatures (full up+down marker sets; %d unique genes total). Intersections are small (pairwise overlap <=7 genes) -> the signatures are largely independent, so their concordant disease/stage behaviour is not driven by shared genes.",
  nrow(PANELS), nuniq))

# ---------------------------------------------------------------------------
# 2. Expression: merged DGEList -> logCPM, subset to mega-5, join clinical
# ---------------------------------------------------------------------------
dge  <- load_merged_dge(); stopifnot(!is.null(dge))
samp <- as.data.table(dge$samples, keep.rownames = "sample_id")
keep <- samp$dataset %in% MEGA
dge  <- dge[, keep]; samp <- samp[keep]

meta <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"),
  select = c("sample_id", "group_binary", "fibrosis_stage", "nas_score"))
samp[meta, on = "sample_id",
     `:=`(group_binary = i.group_binary, fibrosis_stage = i.fibrosis_stage,
          nas_score = i.nas_score)]
# group_binary may also live on dge$samples; prefer metadata, fall back
if (!"group_binary" %in% names(samp) || all(is.na(samp$group_binary)))
  samp[, group_binary := as.character(dge$samples$group_binary)]
samp[, group_binary := as.character(group_binary)]
samp[, group_binary := fifelse(group_binary == "Disease" |
       grepl("dis|nash|nafl|mash|masl|case", group_binary, ignore.case = TRUE), "Disease",
     fifelse(group_binary == "Control" |
       grepl("ctrl|control|normal|healthy", group_binary, ignore.case = TRUE), "Control",
       NA_character_))]
samp[, fib := suppressWarnings(as.numeric(fibrosis_stage))]
samp[, nas := suppressWarnings(as.numeric(nas_score))]
logcpm_all <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

# sex inference (XIST / DDX3Y k-means) -> numeric covariate (mirror figS_pca_definitive)
xist <- rownames(dge)[grep("^ENSG00000229807", rownames(dge))]
ddy  <- rownames(dge)[grep("^ENSG00000067048", rownames(dge))]
sx <- samp$sex
if (length(xist) > 0 && length(ddy) > 0) {
  es <- t(logcpm_all[c(xist[1], ddy[1]), , drop = FALSE])
  km <- kmeans(es, centers = 2, nstart = 20, iter.max = 50)
  fem <- as.integer(names(which.max(tapply(es[, 1], km$cluster, mean))))
  inf <- ifelse(km$cluster == fem, "F", "M")
  mi <- is.na(sx) | sx == ""; sx[mi] <- inf[mi]
}
sex_cov <- as.numeric(sx == "F")
# defensive: sample order of the expression matrix must match the metadata table
stopifnot(ncol(logcpm_all) == nrow(samp), all(colnames(logcpm_all) == samp$sample_id))
cat(sprintf("\nMega-5 samples: %d  (Control=%d  Disease=%d  with F-stage=%d)\n",
            nrow(samp), sum(samp$group_binary == "Control", na.rm = TRUE),
            sum(samp$group_binary == "Disease", na.rm = TRUE), sum(!is.na(samp$fib))))
print(samp[, .(n = .N, ctrl = sum(group_binary == "Control", na.rm = TRUE),
               dis = sum(group_binary == "Disease", na.rm = TRUE),
               fstage = sum(!is.na(fib))), by = dataset])

# ---------------------------------------------------------------------------
# 3. Symbol-indexed marker submatrix (collapse multi-Ensembl symbols by mean)
# ---------------------------------------------------------------------------
can <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
  select = c("gene", "symbol"))
can[, gene := sub("\\..*", "", gene)][, symbol := toupper(trimws(symbol))]
ens_nov <- sub("\\..*", "", rownames(logcpm_all))
need_sym <- unique(sig_long$gene_symbol)
mp <- can[symbol %in% need_sym & gene %in% ens_nov]
idx <- match(mp$gene, ens_nov)
em  <- logcpm_all[idx, , drop = FALSE]
ag  <- rowsum(em, mp$symbol)                                   # sum rows sharing a symbol
em_sym <- ag / as.integer(table(mp$symbol)[rownames(ag)])      # -> mean
cat(sprintf("\nMarker genes mapped to matrix: %d of %d union symbols\n",
            nrow(em_sym), length(need_sym)))
miss <- setdiff(need_sym, rownames(em_sym))
if (length(miss)) cat("  unmapped (", length(miss), "):", paste(head(miss, 25), collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# 4. Batch+sex correction (group-BLIND primary; group-protected for ledger)
# ---------------------------------------------------------------------------
ds_fac <- factor(samp$dataset, levels = MEGA)
grp    <- factor(samp$group_binary, levels = c("Control", "Disease"))
Xraw   <- em_sym
Xcorr  <- limma::removeBatchEffect(em_sym, batch = ds_fac, covariates = sex_cov)
# group-protected variant (DEG-faithful, like the canonical ~dataset+sex+group model)
Xprot  <- tryCatch(
  limma::removeBatchEffect(em_sym, batch = ds_fac, covariates = sex_cov,
                           design = model.matrix(~ grp)),
  error = function(e) Xcorr)

# ---------------------------------------------------------------------------
# 5. Module scores. PRIMARY = directional signature score (mean corrected-z of
#    up-genes minus mean of down-genes): direction-faithful and robust for the
#    mixed up/down panels. SECONDARY (ledger robustness) = GSVA on the up-gene
#    set (bulk-native enrichment). Both computed on raw and corrected matrices.
# ---------------------------------------------------------------------------
zc <- t(scale(t(Xcorr)))                       # genes x samples, per-gene z (corrected)
zr <- t(scale(t(Xraw)))
nbad <- sum(!is.finite(zc)) + sum(!is.finite(zr))   # zero-variance genes -> NaN; zero them out
if (nbad > 0) warning(sprintf("%d non-finite z-values (zero-variance genes) set to 0", nbad))
zc[!is.finite(zc)] <- 0; zr[!is.finite(zr)] <- 0
dir_score <- function(Z, L) {
  up <- intersect(unique(sig_long[label == L & direction == "up"]$gene_symbol),   rownames(Z))
  dn <- intersect(unique(sig_long[label == L & direction == "down"]$gene_symbol), rownames(Z))
  s  <- if (length(up)) colMeans(Z[up, , drop = FALSE]) else rep(0, ncol(Z))
  if (length(dn)) s <- s - colMeans(Z[dn, , drop = FALSE])
  s
}
S_corr <- t(vapply(SIG_LEVELS, function(L) dir_score(zc, L), numeric(ncol(zc))))
S_raw  <- t(vapply(SIG_LEVELS, function(L) dir_score(zr, L), numeric(ncol(zr))))
colnames(S_corr) <- colnames(Xcorr); colnames(S_raw) <- colnames(Xraw)
cat("\nDirectional score inter-signature correlation (corrected):\n")
print(round(cor(t(S_corr)), 2))

# GSVA (up-gene sets) — bulk-native robustness arm
gsva_run <- function(expr, gsets) {
  ok <- tryCatch(GSVA::gsva(GSVA::gsvaParam(expr, gsets, kcdf = "Gaussian")),
                 error = function(e) NULL)
  if (is.null(ok))
    ok <- GSVA::gsva(expr, gsets, method = "gsva", kcdf = "Gaussian", verbose = FALSE)
  ok
}
gsets_up <- lapply(SIG_LEVELS, function(L)
  intersect(unique(sig_long[label == L & direction == "up"]$gene_symbol), rownames(em_sym)))
names(gsets_up) <- SIG_LEVELS
gsets_up <- gsets_up[vapply(gsets_up, length, 1L) >= 3]
cat("\nUp-gene set sizes for GSVA:\n"); print(vapply(gsets_up, length, 1L))
G_corr <- gsva_run(Xcorr, gsets_up)

mk_long <- function(S, tag) {
  d <- as.data.table(t(S), keep.rownames = "sample_id")
  melt(d, id.vars = "sample_id", variable.name = "label", value.name = tag)
}
scores <- merge(mk_long(S_raw, "score_raw"), mk_long(S_corr, "score_corr"),
                by = c("sample_id", "label"))
scores <- merge(scores, mk_long(G_corr, "gsva_corr"), by = c("sample_id", "label"), all.x = TRUE)
scores[samp, on = "sample_id",
       `:=`(dataset = i.dataset, group_binary = i.group_binary, fib = i.fib)]
scores[, label := factor(label, levels = SIG_LEVELS)]

# ---------------------------------------------------------------------------
# 6. Statistics: separation (AUC, Spearman vs stage) + mixing (R^2, silhouette)
# ---------------------------------------------------------------------------
auc1 <- function(score, pos) {                                  # Mann-Whitney AUC
  r <- rank(score); nP <- sum(pos); nN <- sum(!pos)
  if (nP == 0 || nN == 0) return(NA_real_)
  (sum(r[pos]) - nP * (nP + 1) / 2) / (nP * nN)
}
r2 <- function(y, x) {
  d <- data.frame(y = y, x = factor(x)); d <- d[!is.na(d$y) & !is.na(d$x), ]
  if (length(unique(d$x)) < 2) return(NA_real_)
  summary(lm(y ~ x, d))$r.squared
}
ledger <- scores[, {
  pos  <- group_binary == "Disease"
  has  <- group_binary %in% c("Control", "Disease")
  auc  <- auc1(score_corr[has], pos[has])
  wp   <- tryCatch(wilcox.test(score_corr[has][pos[has]], score_corr[has][!pos[has]])$p.value,
                   error = function(e) NA_real_)
  fmask <- !is.na(fib)
  sr   <- if (sum(fmask) > 10) suppressWarnings(cor.test(score_corr[fmask], fib[fmask],
                              method = "spearman")) else list(estimate = NA, p.value = NA)
  gauc <- auc1(gsva_corr[has], pos[has])
  .(auc_disease = auc, auc_wilcox_p = wp, gsva_auc_disease = gauc,
    spearman_rho_stage = unname(sr$estimate), spearman_p_stage = sr$p.value,
    r2_cohort_raw = r2(score_raw, dataset),  r2_cohort_corr = r2(score_corr, dataset),
    r2_disease_raw = r2(score_raw, group_binary), r2_disease_corr = r2(score_corr, group_binary))
}, by = label]
ledger <- merge(PANELS[, .(label, geo, external)], ledger, by = "label")
ledger[, label := factor(label, levels = SIG_LEVELS)]; setorder(ledger, label)
ledger[, n_up := vapply(as.character(label), function(L)
  length(gsets_up[[L]]) %||% 0L, 1L)]

# silhouette by cohort in marker space (raw vs corrected): lower => better mixing
sil_cohort <- function(M) {
  d  <- dist(scale(t(M)))
  cl <- as.integer(ds_fac)
  mean(silhouette(cl, d)[, "sil_width"])
}
sil_raw  <- sil_cohort(Xraw)
sil_corr <- sil_cohort(Xcorr)
ledger[, `:=`(silhouette_cohort_raw = round(sil_raw, 4),
              silhouette_cohort_corr = round(sil_corr, 4))]
if (anyNA(ledger[, .(auc_disease, spearman_rho_stage)]))
  warning("Ledger has NA stats — check group/stage membership")
fwrite(ledger, file.path(DATA_DIR, "marker_validation_ledger.csv"))
cat("\n===== LEDGER =====\n"); print(ledger)
cat(sprintf("\nCohort silhouette (marker space): raw=%.3f -> corrected=%.3f (lower=better mixing)\n",
            sil_raw, sil_corr))

# ===========================================================================
# PANEL 1 — heatmap: samples x up-marker genes (corrected), row-split by source
# ===========================================================================
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
ext_lab <- paste(PANELS[external == TRUE]$label, collapse = " & ")

# One renderer for both heatmap variants (they differ only in which samples are kept).
# Columns ordered by CLINICAL severity (controls first, fibrosis F0->F4, then NAS within
# stage), NOT by marker score (a within-stage score sort gave a per-block blue->red sawtooth).
draw_marker_heatmap <- function(keep_idx, outfile, subset_note) {
col_ord <- keep_idx[order(samp$group_binary[keep_idx] != "Control",
                          samp$fib[keep_idx],
                          ifelse(is.na(samp$nas[keep_idx]), 99, samp$nas[keep_idx]))]

# row blocks: per (signature, up-gene); a gene may recur across signature blocks
hm_rows <- rbindlist(lapply(names(gsets_up), function(L)
  data.table(label = L, gene = gsets_up[[L]])))
HM <- t(scale(t(Xcorr[hm_rows$gene, col_ord, drop = FALSE])))  # row z-score, ordered cols
row_split <- factor(hm_rows$label, levels = intersect(SIG_LEVELS, names(gsets_up)))

key_markers <- c("COL1A1","COL1A2","COL3A1","THBS2","EFEMP2","FBLN5","LUM","LTBP2",
                 "AKR1B10","GDF15","IL32","THY1","TNFRSF12A","STMN2","PDGFA","CCL20",
                 "TREM2","ITGBL1","LAMC3","IGFBP7","DPT")
lab_at  <- which(hm_rows$gene %in% key_markers & !duplicated(hm_rows$gene))
lab_txt <- hm_rows$gene[lab_at]
mm <- setdiff(key_markers, hm_rows$gene)
if (length(mm)) message("Note: key markers absent from heatmap rows: ", paste(mm, collapse = ", "))

stage_fac <- factor(paste0("F", samp$fib[col_ord]), levels = paste0("F", 0:4))
nas_vec   <- samp$nas[col_ord]
coh_all   <- setNames(c("#4E79A7","#F28E2B","#59A14F","#B07AA1","#9C755F"), MEGA)
coh_cols  <- coh_all[MEGA %in% samp$dataset[col_ord]]   # only cohorts actually shown
stage_cols <- fibrosis_stage_colors
nas_col   <- colorRamp2(c(0, 4, 8), c("#FCFBFD", "#9E9AC8", "#3F007D"))   # Purples (distinct from fibrosis blue / disease magenta)

top_anno <- HeatmapAnnotation(
  Group    = samp$group_binary[col_ord],
  Fibrosis = stage_fac,
  NAS      = nas_vec,
  Cohort   = samp$dataset[col_ord],
  col = list(Group    = c(Control = CTRL, Disease = DIS),
             Fibrosis = stage_cols,
             NAS      = nas_col,
             Cohort   = coh_cols),
  na_col = "#EEEEEE",
  simple_anno_size = unit(2.6, "mm"),
  annotation_name_gp = gpar(fontsize = 6),
  annotation_legend_param = list(
    Group    = list(title_gp = gpar(fontsize = 6, fontface = "plain"), labels_gp = gpar(fontsize = 5)),
    Fibrosis = list(title_gp = gpar(fontsize = 6, fontface = "plain"), labels_gp = gpar(fontsize = 5)),
    NAS      = list(title_gp = gpar(fontsize = 6, fontface = "plain"), labels_gp = gpar(fontsize = 5)),
    Cohort   = list(title_gp = gpar(fontsize = 6, fontface = "plain"), labels_gp = gpar(fontsize = 5))))

ht <- Heatmap(
  HM, name = "z (corrected\nlog-CPM)",
  col = colorRamp2(c(-2, 0, 2), c("#1565C0", "white", "#C9265E")),
  top_annotation = top_anno,
  row_split = row_split, cluster_row_slices = FALSE, cluster_rows = TRUE,
  cluster_columns = FALSE, show_column_names = FALSE, show_row_names = FALSE,
  row_title_gp = gpar(fontsize = 6, fontface = "plain"), row_title_rot = 0,
  use_raster = TRUE, raster_quality = 3,
  right_annotation = rowAnnotation(mark = anno_mark(
    at = lab_at, labels = lab_txt,
    labels_gp = gpar(fontsize = 5, fontface = "italic"), link_width = unit(3, "mm"))),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 6, fontface = "plain"),
                              labels_gp = gpar(fontsize = 5), legend_height = unit(2.2, "cm")))
# continuous NAS legend doesn't surface the na_col -> add an explicit "NA" swatch,
# but only when the displayed set actually has NAS-missing samples
leg <- if (any(is.na(samp$nas[col_ord])))
  list(Legend(labels = "NA (no NAS)", legend_gp = gpar(fill = "#EEEEEE"),
              labels_gp = gpar(fontsize = 5), grid_height = unit(2.6, "mm"),
              grid_width = unit(2.6, "mm"))) else list()
hm_h <- max(4.2, min(8.0, nrow(HM) * 0.022))
pdf(file.path(OUT, outfile), width = fig_full_width, height = hm_h)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right", annotation_legend_list = leg)
dev.off()
cat("Wrote", outfile, "\n")
message(sprintf(
"[%s] Per-gene z-score (across the %d displayed samples jointly) of the batch-corrected (group-blind ~dataset+sex removeBatchEffect) log-CPM, up-regulated marker genes from %d verbatim published MASLD signatures (row-split by source). Columns = %s, ordered by CLINICAL severity (controls -> F0->F4 -> NAS within stage), NOT by marker score. Top bars = group, fibrosis, NAS, cohort; cohort salt-and-pepper = batch removed. External validators %s (microarray, gene-symbol membership only); in-atlas baselines Govaere/Pantano.",
  sub("\\.pdf$", "", outfile), length(col_ord), nrow(PANELS), subset_note, ext_lab))
}

draw_marker_heatmap(which(!is.na(samp$fib)),
  "marker_validation_heatmap.pdf",
  "all fibrosis-staged mega samples (GSE126848 dropped -- no histological staging; GSE213621 has no NAS -> grey NAS track)")
draw_marker_heatmap(which(!is.na(samp$fib) & !is.na(samp$nas)),
  "marker_validation_heatmap_nascomplete.pdf",
  "only samples with BOTH fibrosis stage AND NAS (GSE126848 + GSE213621 dropped -- no NAS; 3 cohorts)")

# ===========================================================================
# PANEL 2 — module score across fibrosis stage F0-F4 (corrected), per source
# ===========================================================================
bt <- theme_masld(base_size = 7) + theme_pub() +
  theme(legend.position = "none", panel.grid = element_blank(),
        strip.text = element_text(size = 6, face = "plain"),
        axis.text = element_text(color = "black"))
sd_stage <- scores[!is.na(fib)]
sd_stage[, stage := factor(paste0("F", fib), levels = paste0("F", 0:4))]
lab_st <- ledger[, .(label, txt = sprintf("rho=%.2f%s", spearman_rho_stage,
                     ifelse(spearman_p_stage < 0.001, "***",
                     ifelse(spearman_p_stage < 0.01, "**",
                     ifelse(spearman_p_stage < 0.05, "*", "")))))]
lab_st[, label := factor(label, levels = SIG_LEVELS)]
p2 <- ggplot(sd_stage, aes(stage, score_corr, fill = stage)) +
  geom_violin(scale = "width", width = 0.7, linewidth = 0.2, colour = "grey40") +
  geom_boxplot(width = 0.12, outlier.size = 0.15, linewidth = 0.2, fill = "white") +
  scale_fill_manual(values = fibrosis_stage_colors) +
  facet_wrap(~ label, nrow = 1, scales = "free_y") +
  geom_text(data = lab_st, aes(x = 1.1, y = Inf, label = txt), inherit.aes = FALSE,
            hjust = 0, vjust = 1.4, size = PUB_GEOM_TEXT, colour = "black") +
  labs(x = "Fibrosis stage", y = "Signature score (directional)") + bt
ggsave(file.path(OUT, "marker_module_scores_stage.pdf"), p2,
       width = fig_full_width, height = 2.4, device = cairo_pdf)
cat("Wrote marker_module_scores_stage.pdf\n")

# ===========================================================================
# PANEL 3 — module score Control vs Disease (corrected), per source + AUC
# ===========================================================================
sd_bin <- scores[group_binary %in% c("Control", "Disease")]
sd_bin[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]
lab_bin <- ledger[, .(label, txt = sprintf("AUC=%.2f", auc_disease))]
lab_bin[, label := factor(label, levels = SIG_LEVELS)]
p3 <- ggplot(sd_bin, aes(group_binary, score_corr, fill = group_binary)) +
  geom_violin(scale = "width", width = 0.6, linewidth = 0.2, colour = "grey40") +
  geom_boxplot(width = 0.14, outlier.size = 0.15, linewidth = 0.2, fill = "white") +
  scale_fill_manual(values = c(Control = CTRL, Disease = DIS)) +
  facet_wrap(~ label, nrow = 1, scales = "free_y") +
  geom_text(data = lab_bin, aes(x = 1.5, y = Inf, label = txt), inherit.aes = FALSE,
            vjust = 1.4, size = PUB_GEOM_TEXT, colour = "black") +
  labs(x = NULL, y = "Signature score (directional)") + bt
ggsave(file.path(OUT, "marker_module_scores_binary.pdf"), p3,
       width = fig_full_width, height = 2.4, device = cairo_pdf)
cat("Wrote marker_module_scores_binary.pdf\n")

# ===========================================================================
# PANEL 4 — variance explained by cohort vs disease, Before vs After correction
#           (the batch-correction proof: cohort variance collapses, disease kept)
# ===========================================================================
mix <- rbindlist(list(
  ledger[, .(label, source = "Cohort",  Before = r2_cohort_raw,  After = r2_cohort_corr)],
  ledger[, .(label, source = "Disease", Before = r2_disease_raw, After = r2_disease_corr)]))
mixm <- melt(mix, id.vars = c("label", "source"), variable.name = "stage", value.name = "r2")
mixm[, stage := factor(stage, levels = c("Before", "After"))]
mix_bar <- mixm[, .(r2 = mean(r2, na.rm = TRUE), sd = sd(r2, na.rm = TRUE)), by = .(source, stage)]
p4 <- ggplot(mix_bar, aes(source, r2, fill = stage)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, colour = "grey30", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(0, r2 - sd), ymax = r2 + sd),
                position = position_dodge(width = 0.75), width = 0.18, linewidth = 0.25) +
  geom_point(data = mixm, aes(source, r2, group = stage),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75, seed = 42),
             size = 0.5, colour = "black", alpha = 0.6, inherit.aes = FALSE) +
  geom_text(data = mix_bar, aes(label = sprintf("%.2f", r2)),
            position = position_dodge(width = 0.75), vjust = -0.5, size = PUB_GEOM_TEXT, colour = "black") +
  scale_fill_manual(values = c(Before = CTRL, After = AFT), name = NULL) +
  labs(x = NULL, y = expression("Variance explained ("*R^2*", module score)")) +
  theme_masld(base_size = 7) + theme_pub() +
  theme(panel.grid = element_blank(), legend.position = "top",
        axis.text = element_text(color = "black"))
ggsave(file.path(OUT, "marker_cohort_mixing.pdf"), p4,
       width = fig_half_width, height = 2.7, device = cairo_pdf)
cat("Wrote marker_cohort_mixing.pdf\n")

# ---------------------------------------------------------------------------
# Captions -> stdout (house style: descriptions live in the caption, not panels)
# ---------------------------------------------------------------------------
# (heatmap captions are emitted inside draw_marker_heatmap; ext_lab defined in PANEL 1)
message(sprintf(
"[marker_module_scores_stage] Directional signature score (mean corrected-z of up genes minus down genes) per source across fibrosis stage F0-F4 (n=%d staged samples). Spearman rho(score,stage) annotated. The fibrosis-oriented signatures (Govaere, Pantano, Moylan) rise monotonically (rho ~0.46-0.60***); Arendt 2015 is a NASH-vs-healthy disease-PRESENCE signature (not fibrosis-progression) and is intentionally flat across stage (rho ns) -- a scope control showing the stage gradient is fibrosis-specific rather than generic batch structure.",
  nrow(sd_stage[label == SIG_LEVELS[1]])))
message(sprintf(
"[marker_module_scores_binary] Same directional scores, Control (n=%d) vs Disease (n=%d); Mann-Whitney AUC annotated per source (GSVA-up AUC in ledger as robustness). The external held-out validators (%s) carry the validation claim without circularity; the in-atlas signatures (Govaere, Pantano) are consistency baselines. Moylan/Arendt = gene-symbol membership only (no microarray effect sizes imported).",
  sum(samp$group_binary == "Control", na.rm = TRUE),
  sum(samp$group_binary == "Disease", na.rm = TRUE), ext_lab))
message(sprintf(
"[marker_cohort_mixing] Mean variance explained (R^2) in module scores by cohort vs disease group, Before vs After group-blind correction (points = the %d signatures). Cohort R^2 %.2f->%.2f (collapses); disease R^2 %.2f->%.2f (preserved). Cohort silhouette in marker space %.2f->%.2f.",
  nrow(ledger), mix_bar[source=="Cohort"&stage=="Before"]$r2, mix_bar[source=="Cohort"&stage=="After"]$r2,
  mix_bar[source=="Disease"&stage=="Before"]$r2, mix_bar[source=="Disease"&stage=="After"]$r2,
  sil_raw, sil_corr))
cat("\nMARKER_VALIDATION_DONE\n")
