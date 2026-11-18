# In-memory 22-autosome LOCO InstaPrism across all 5 Resource bulk cohorts + SteatoSITE (515 samples),
# using the 100,000-cell NUCLEUS-FIX SingleCellExperiment reference loaded once into RAM.
# Outputs separated 8-part ILR (Hepatocytes and Cholangiocytes split -> ilr1..ilr7) for standard and chr1..chr22.

suppressPackageStartupMessages({
  library(InstaPrism)
  library(SingleCellExperiment)
  library(Matrix)
  library(data.table)
})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CAND <- file.path(ROOT, "Analysis/Deconvolution/candidates/loco_20261118")
CHR_DIR <- file.path(CAND, "chr")
dir.create(CHR_DIR, recursive = TRUE, showWarnings = FALSE)

NCORE <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
MAX_RSS_GB <- 190.0

read_proc_mem_gb <- function() {
  lines <- readLines("/proc/self/status", warn = FALSE)
  get_val <- function(key) {
    hit <- grep(paste0("^", key, ":"), lines, value = TRUE)
    if (!length(hit)) return(NA_real_)
    kb <- as.numeric(gsub("[^0-9]", "", hit[1]))
    kb / (1024^2)
  }
  c(vmhwm_gb = get_val("VmHWM"), vmrss_gb = get_val("VmRSS"))
}

check_rss_limit <- function(stage) {
  m <- read_proc_mem_gb()
  cat(sprintf("[MEM %s] VmHWM=%.2f GB, VmRSS=%.2f GB\n", stage, m["vmhwm_gb"], m["vmrss_gb"]))
  if (!is.na(m["vmhwm_gb"]) && m["vmhwm_gb"] >= MAX_RSS_GB) {
    stop(sprintf("REPLAN TRIGGER #2: RSS %.2f GB reached limit %.2f GB (< 200 GB)", m["vmhwm_gb"], MAX_RSS_GB))
  }
  m
}

# 1. Load 100k SCE reference once into RAM
t_start <- Sys.time()
sce <- readRDS(file.path(CAND, "reference_nucleus_fix_100k_sce.rds"))
sc_mat <- counts(sce)
sc_ens <- as.character(rowData(sce)$ensembl_base)
sc_sym <- rownames(sce)
ct <- as.character(colData(sce)$cell_type)
rm(sce); gc(FALSE)
check_rss_limit("after_sce_load")

# 2. Map SCE and bulk features to chromosomes via gencode_v49_gene_metadata.tsv.gz
gmeta_path <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
gmeta <- fread(cmd = paste("gzip -dc", shQuote(gmeta_path)))
# Prefer ensembl_base match for SCE, fallback gene_name
sc_chr_ens <- gmeta$chromosome[match(sc_ens, gmeta$ensembl_base)]
sc_chr_sym <- gmeta$chromosome[match(sc_sym, gmeta$gene_name)]
sc_chr <- ifelse(!is.na(sc_chr_ens) & nzchar(sc_chr_ens), sc_chr_ens, sc_chr_sym)
cat("SCE genes mapped to chr1..22:", sum(sc_chr %in% paste0("chr", 1:22)), "/", length(sc_sym), "\n")

# For bulk symbol rows, build conservative set of chromosomes per symbol or ensembl_base
g_sym_chr <- unique(rbind(
  gmeta[!is.na(gene_name) & nzchar(gene_name), .(key = gene_name, chromosome)],
  gmeta[!is.na(ensembl_base) & nzchar(ensembl_base), .(key = ensembl_base, chromosome)],
  gmeta[!is.na(gene_id) & nzchar(gene_id), .(key = gene_id, chromosome)]
))
sym_to_chrs <- split(g_sym_chr$chromosome, g_sym_chr$key)

# 3. Load all 5 Resource cohorts + SteatoSITE bulk matrices
read_bulk_tsv <- function(path) {
  dt <- fread(path)
  sym <- as.character(dt[[1L]])
  mat <- as.matrix(dt[, -1L, with = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- sym
  mat
}

cohorts <- list(
  GSE126848  = read_bulk_tsv(file.path(ROOT, "Analysis/Deconvolution/corrected_774_20261009T141331Z/bulk/GSE126848/GSE126848_counts.tsv")),
  GSE130970  = read_bulk_tsv(file.path(ROOT, "Analysis/Deconvolution/corrected_774_20261009T141331Z/bulk/GSE130970/GSE130970_counts.tsv")),
  GSE135251  = read_bulk_tsv(file.path(ROOT, "Analysis/Deconvolution/corrected_774_20261009T141331Z/bulk/GSE135251/GSE135251_counts.tsv")),
  GSE162694  = read_bulk_tsv(file.path(ROOT, "Analysis/Deconvolution/corrected_774_20261009T141331Z/bulk/GSE162694/GSE162694_counts.tsv")),
  GSE213621  = read_bulk_tsv(file.path(ROOT, "Analysis/Deconvolution/corrected_774_20261009T141331Z/bulk/GSE213621/GSE213621_counts.tsv")),
  SteatoSITE = read_bulk_tsv(file.path(ROOT, "Analysis/Hypotheses/SS-20261025T135237Z/out_s3e-20261026T233706Z/bulk/SteatoSITE_counts.tsv"))
)

blk_rows <- rownames(cohorts[[1L]])
for (nm in names(cohorts)) {
  stopifnot(identical(rownames(cohorts[[nm]]), blk_rows))
  cat(sprintf("Loaded bulk cohort %-11s: %d x %d\n", nm, nrow(cohorts[[nm]]), ncol(cohorts[[nm]])))
}
stopifnot(sum(vapply(cohorts, ncol, integer(1))) == 774L + 515L)

# Map each bulk rowname to chromosomes
blk_chr_idx <- sym_to_chrs[blk_rows]
names(blk_chr_idx) <- blk_rows

# 4. Separated 8-part ILR transformation (Hepatocytes & Cholangiocytes split -> ilr1..ilr7)
LINEAGES_8 <- c("Hepatocytes", "Cholangiocytes", "endothelial", "fibroblasts", "mac_mono", "other_myeloid", "t_nk", "b_plasma")
MAP_8 <- list(
  Hepatocytes    = "Hepatocytes",
  Cholangiocytes = "Cholangiocytes",
  endothelial    = "Endothelial cells",
  fibroblasts    = "Fibroblasts",
  mac_mono       = c("Macrophages", "Mono+mono derived cells"),
  other_myeloid  = c("Basophils", "Neutrophils", "cDC1s", "cDC2s", "pDCs"),
  t_nk           = c("T cells", "Circulating NK/NKT", "Resident NK"),
  b_plasma       = c("B cells", "Plasma cells")
)
H8 <- contr.helmert(8)
H8 <- sweep(H8, 2, sqrt(colSums(H8^2)), "/")

to_separated_8_ilr <- function(frac_mat, dataset_name) {
  P8 <- vapply(LINEAGES_8, function(lin) {
    cols <- intersect(MAP_8[[lin]], colnames(frac_mat))
    if (!length(cols)) rep(0, nrow(frac_mat)) else rowSums(frac_mat[, cols, drop = FALSE])
  }, numeric(nrow(frac_mat)))
  rownames(P8) <- rownames(frac_mat)
  P8_fl <- pmax(P8, 1e-6)
  P8_fl <- P8_fl / rowSums(P8_fl)
  ilr_m <- log(P8_fl) %*% H8
  colnames(ilr_m) <- paste0("ilr", 1:7)
  dt <- data.table(dataset = dataset_name, sample_id = rownames(frac_mat))
  dt <- cbind(dt, as.data.table(P8), as.data.table(ilr_m))
  dt
}

run_one_set <- function(tag, chr_holdout = NULL) {
  t0 <- Sys.time()
  if (is.null(chr_holdout)) {
    keep_sc <- rep(TRUE, nrow(sc_mat))
    keep_blk <- rep(TRUE, length(blk_rows))
  } else {
    keep_sc <- is.na(sc_chr) | (sc_chr != chr_holdout)
    keep_blk <- !vapply(blk_chr_idx, function(ch) !is.null(ch) && (chr_holdout %in% ch), logical(1))
  }
  ref <- InstaPrism::refPrepare(sc_Expr = sc_mat[keep_sc, , drop = FALSE],
                                cell.type.labels = ct, cell.state.labels = ct)
  frac_list <- list()
  ilr_list <- list()
  for (cname in names(cohorts)) {
    blk_sub <- cohorts[[cname]][keep_blk, , drop = FALSE]
    res <- InstaPrism::InstaPrism(bulk_Expr = blk_sub, refPhi_cs = ref,
                                  n.core = NCORE, verbose = FALSE)
    fm <- t(res@Post.ini.ct@theta)
    fdt <- cbind(data.table(dataset = cname, sample_id = rownames(fm)), as.data.table(fm))
    frac_list[[cname]] <- fdt
    ilr_list[[cname]] <- to_separated_8_ilr(fm, cname)
    rm(res, blk_sub)
  }
  rm(ref); gc(FALSE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  mem <- check_rss_limit(tag)
  list(
    frac = rbindlist(frac_list, fill = TRUE),
    ilr  = rbindlist(ilr_list, fill = TRUE),
    audit = data.table(
      run = tag,
      held_out_chr = if (is.null(chr_holdout)) "NONE" else chr_holdout,
      n_sc_genes_kept = sum(keep_sc),
      n_sc_genes_dropped = sum(!keep_sc),
      n_blk_genes_kept = sum(keep_blk),
      n_blk_genes_dropped = sum(!keep_blk),
      n_samples = sum(vapply(cohorts, ncol, integer(1))),
      elapsed_sec = elapsed,
      vmhwm_gb = unname(mem["vmhwm_gb"]),
      vmrss_gb = unname(mem["vmrss_gb"])
    )
  )
}

# Standard non-LOCO separated 8-part ILR
cat("Running standard (non-LOCO) InstaPrism across 1,289 samples...\n")
std <- run_one_set("standard", NULL)
fwrite(std$frac, file.path(CAND, "standard_fractions.tsv.gz"), sep = "\t")
fwrite(std$ilr,  file.path(CAND, "standard_ilr8.tsv.gz"), sep = "\t")
ss_std <- copy(std$ilr[dataset == "SteatoSITE"])
setnames(ss_std, "sample_id", "sample_alias")
ss_std[, dataset := NULL]
fwrite(ss_std, file.path(CAND, "standard_steatosite_ilr8.tsv"), sep = "\t")

audits <- list(std$audit)
all_loco_ilr <- list()

for (k in 1:22) {
  ck <- paste0("chr", k)
  cat(sprintf("Running LOCO %s (%d/22)...\n", ck, k))
  out_k <- run_one_set(ck, ck)
  fwrite(out_k$frac, file.path(CHR_DIR, paste0(ck, "_fractions.tsv.gz")), sep = "\t")
  fwrite(out_k$ilr,  file.path(CHR_DIR, paste0(ck, "_ilr8.tsv.gz")), sep = "\t")
  ss_k <- copy(out_k$ilr[dataset == "SteatoSITE"])
  setnames(ss_k, "sample_id", "sample_alias")
  ss_k[, dataset := NULL]
  fwrite(ss_k, file.path(CHR_DIR, paste0(ck, "_steatosite_ilr8.tsv")), sep = "\t")
  audits[[ck]] <- out_k$audit
  all_loco_ilr[[ck]] <- cbind(data.table(held_out_chr = ck), out_k$ilr)
  fwrite(rbindlist(audits), file.path(CAND, "loco_run_audit.tsv"), sep = "\t")
}

fwrite(rbindlist(all_loco_ilr), file.path(CAND, "loco_ilr8_all_chr.tsv.gz"), sep = "\t")
audit_df <- rbindlist(audits)
fwrite(audit_df, file.path(CAND, "loco_run_audit.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(CAND, "sessionInfo.txt"))
cat(sprintf("LOCO_INSTAPRISM_COMPLETE total_sec=%.1f max_VmHWM_GB=%.2f\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs")),
            max(audit_df$vmhwm_gb, na.rm = TRUE)))
