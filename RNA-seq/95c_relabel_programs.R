#!/usr/bin/env Rscript
# 95c_relabel_programs.R — Bravo Task 9
#
# Apply the rewritten tag_program() (95_nmf_clean_ksweep.R) to an existing
# cached NMF fit and emit an updated program_labels.csv.  Lineage-marker
# overrides are inspected BEFORE the stage_rho directional rule.
#
# Input:
#   RNA-seq/results/subtypes/nmf_results_cache_clean.rds  (default; override via NMF_CACHE_PATH)
# Output:
#   RNA-seq/results/subtypes/program_labels.csv  (suffix follows cache name)
#
# Env:
#   NMF_CACHE_PATH   — path to cache rds (default canonical clean)
#   NMF_CHOSEN_K     — k to label (default 6)
#
# Idempotent: safe to re-run after bravo-prot-only-refit to regenerate labels
# on the new cache.

suppressPackageStartupMessages({
  library(data.table); library(dplyr, warn.conflicts = FALSE)
  library(NMF); library(fgsea); library(msigdbr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
out_dir    <- file.path(BASE, "RNA-seq/results/subtypes")
cache_path <- Sys.getenv("NMF_CACHE_PATH",
                         file.path(out_dir, "nmf_results_cache_clean.rds"))
K <- as.integer(Sys.getenv("NMF_CHOSEN_K", "6"))

atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
meta_path  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")

# Derive labels output suffix from cache filename so prot_only variants don't
# clobber canonical labels. Also suffix by k when k != 6 (the default) so
# protonly k=6 and protonly k=8 label files don't collide.
.cache_base <- sub("\\.rds$", "", basename(cache_path))
.suffix <- sub("^nmf_results_cache_clean", "", .cache_base)
.k_tag <- if (K == 6L) "" else sprintf("_k%d", K)
labels_out <- file.path(out_dir, sprintf("program_labels%s%s.csv", .suffix, .k_tag))

cat(sprintf("=== 95c relabel: cache=%s, k=%d ===\n", cache_path, K))
stopifnot("Cache not found" = file.exists(cache_path))
cached <- readRDS(cache_path)
fit <- cached$nmf_results[[as.character(K)]]
stopifnot("No fit at requested k in cache" = !is.null(fit))

gene_ids <- rownames(cached$mat_nn)
sample_ids <- colnames(cached$mat_nn)

W <- basis(fit); H <- coef(fit)
rownames(W) <- gene_ids; colnames(H) <- sample_ids

# Symbol lookup
atlas <- fread(atlas_path, select = c("ensembl_id","human_symbol"))
id2sym <- setNames(atlas$human_symbol, atlas$ensembl_id)
strip_ver <- function(x) sub("\\.[0-9]+$", "", x)
sym_W <- id2sym[strip_ver(gene_ids)]; names(sym_W) <- gene_ids

# ================= Lineage signatures (source of truth) =================
FIBROTIC_CANON <- c(
  "COL1A1","COL1A2","COL3A1","COL5A1","COL6A3","COL15A1",
  "LUM","DCN","FBN1","LOX","SPARC","ELN","FN1",
  "ACTA2","DES","MYOCD","CNN1","TAGLN","MYL9","MYH11","FLNC",
  "PDGFRB","THY1","VIM","CDH11","LTBP2","CTSK",
  "TIMP1","PDPN","SMOC2","EDIL3","POSTN","RCN3","BGN",
  "LRRC15","SFRP4","PTPRQ")
NEUTROPHIL_MARKERS <- c("S100A8","S100A9","CXCR1","CXCR2","FCGR3B",
                        "MPO","ELANE","LCN2","LTF","CSF3R")
VSMC_PERICYTE_MARKERS <- c("DES","CNN1","MYOCD","ACTG2","ACTA2","MYH11",
                           "TAGLN","MYL9","RGS5")
HSC_COLLAGEN_MARKERS <- c("COL1A1","COL3A1","CTHRC1","LUM","BGN")

tag_program <- function(top50_syms, fgsea_tbl, stage_rho, switch_ratio, pct_dom) {
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

  if (n_neut >= 3) return("Innate-immune-inflammation")
  if (n_vsmc_top50 >= 3 && n_hsc_col == 0) return("Vascular-smooth-muscle / pericyte")
  if (n_vsmc_top50 >= 2 && n_hsc_col >= 1) return("Activated-HSC / myofibroblast")
  # NOTE: "Fibrogenic" -> "Fibrotic-ECM" per 2026-06-24 pathway curation in 95d.
  # 95d (the FINAL labeling authority in 44->95->95c->95d) overwrites these labels
  # with CURATED_LABELS. Values here kept consistent so 95c standalone is coherent.
  if (n_fib >= 3 && emt_NES >= 2 && !is.na(stage_rho) && stage_rho > 0.5) return("Fibrotic-ECM")
  if (is.na(stage_rho)) return("Unlabeled")
  if (stage_rho > 0.6) {
    if (inf_NES >= 2) return("Inflammatory-EMT")
    return("Progression-Other")
  }
  if (stage_rho < -0.6) return("Hepatocyte-Metabolic")
  if (!is.na(switch_ratio) && switch_ratio > 1.3) return("Transition-Late")
  if (!is.na(switch_ratio) && switch_ratio < 0.77) return("Transition-Early")
  "Stable"
}

# ================= Per-program characterization =================
hallmark  <- msigdbr(species = "Homo sapiens", collection = "H")
hall_list <- split(hallmark$gene_symbol, hallmark$gs_name)

meta <- fread(meta_path)[sample_id %in% sample_ids]; setkey(meta, sample_id)
fib_int <- as.integer(meta$fibrosis_stage[match(sample_ids, meta$sample_id)])

W_centered <- sapply(seq_len(K), function(j) W[, j] - rowMeans(W[, -j, drop = FALSE]))
rownames(W_centered) <- rownames(W)
dom <- apply(H, 2, which.max)

stage_rho <- rep(NA_real_, K)
for (j in seq_len(K)) {
  prop <- vapply(0:4, function(s) {
    n_s <- sum(fib_int == s, na.rm = TRUE)
    if (n_s == 0) NA_real_ else mean(dom[fib_int == s] == j, na.rm = TRUE)
  }, numeric(1))
  if (sum(!is.na(prop)) >= 3)
    stage_rho[j] <- suppressWarnings(cor(prop, 0:4, method = "spearman", use = "pairwise"))
}

rows <- list()
# Hard-coded prior-labels snapshot (pre-Task-9 rewrite, 2026-04-22 canonical
# clean cache). Pinned here so the audit column remains stable across re-runs.
# Only applied when the cache being relabeled IS the canonical clean cache;
# NMF program indices are permuted across fits, so the Px labels from the
# canonical cache do not correspond to Px in protonly/nonprotonly caches.
# NOTE: these are the historical pre-pathway-curation labels for audit; the
# 2026-06-24 pathway-curated names are in 95d's CURATED_LABELS (the authority).
HISTORICAL_PRIOR <- c(
  P1 = "Progression-Inflammatory",
  P2 = "Quiescent-Parenchyma_1",
  P3 = "Quiescent-Parenchyma_2",
  P4 = "Stable",
  P5 = "Progression-Other",
  P6 = "Fibrogenic"
)
IS_CANONICAL_CACHE <- basename(cache_path) == "nmf_results_cache_clean.rds"
if (!IS_CANONICAL_CACHE) {
  cat("Non-canonical cache — prior_label column left blank (program indices\n",
      "  permute across fits; cross-cache comparison requires Hungarian match).\n",
      sep = "")
}

for (j in seq_len(K)) {
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
  ok_f <- !is.na(fib_int)
  m_lo <- mean(h_j[ok_f & fib_int <= 2], na.rm = TRUE)
  m_hi <- mean(h_j[ok_f & fib_int >= 3], na.rm = TRUE)
  sw_ratio <- m_hi / (m_lo + 1e-10)
  pct_dom_j <- mean(dom == j) * 100

  new_label <- tag_program(top50_syms, fg, stage_rho[j], sw_ratio, pct_dom_j)

  sigp <- as.data.table(fg)[padj < 0.05 & NES > 0][order(-NES)][1, .(pathway, NES)]
  top_path <- if (nrow(sigp)) sprintf("%s(NES=%.2f)", sigp$pathway, sigp$NES) else "ns"

  prog_code <- paste0("P", j)
  prior_label <- if (IS_CANONICAL_CACHE) unname(HISTORICAL_PRIOR[prog_code]) else ""
  if (is.na(prior_label)) prior_label <- ""

  rows[[j]] <- data.table(
    program_code = prog_code,
    biological_label = new_label,
    label_category = new_label,
    top_pathway = top_path,
    top_genes = paste(head(top50_syms, 15), collapse = ","),
    stage_rho = round(stage_rho[j], 3),
    switch_ratio = round(sw_ratio, 3),
    pct_dominant = round(pct_dom_j, 2),
    prior_label = prior_label,
    relabel_comment = ifelse(nzchar(prior_label) && prior_label != new_label,
      sprintf("RELABELED: prior '%s' was mis-assigned by stage_rho-only rule; new lineage-aware tag_program() in 95_nmf_clean_ksweep.R", prior_label),
      "")
  )
}

out_dt <- rbindlist(rows)
# Header comment documenting the source of truth
hdr <- c(
  paste0("# program_labels.csv regenerated by 95c_relabel_programs.R on ", Sys.Date()),
  "# Labels from lineage-aware tag_program() (95_nmf_clean_ksweep.R:tag_program, Task 9 rewrite).",
  "# Prior labels (pre-rewrite) preserved in prior_label column for audit."
)
writeLines(hdr, labels_out)
fwrite(out_dt, labels_out, append = TRUE, col.names = TRUE)
cat("Wrote:", labels_out, "\n")
print(out_dt[, .(program_code, biological_label, prior_label, relabel_comment)])
