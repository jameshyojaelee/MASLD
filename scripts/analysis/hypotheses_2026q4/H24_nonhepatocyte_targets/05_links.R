#!/usr/bin/env Rscript
# H24 step 5: donor-level peak-gene links (layers L1 and L2). See PRESPEC_H24.md.
#
# L1, GSE296875 (same-nucleus multiome). Per lineage, donor units with >= 20
#   nuclei; lineage tested if >= 15 donors.
#   RNA: raw UMIs, edgeR filterByExpr (intercept design), TMM, logCPM (prior 2).
#   ATAC: counts in the authors' peaks. A peak is kept if its raw count is >= 5
#   in >= 50% of the lineage's donor units (filterByExpr is calibrated for genes
#   and is not used for peaks); TMM on kept peaks; logCPM (prior 2).
#   Both residualized on well (pooled-donor batch) and log(nuclei); statistic is
#   Pearson r between residuals across donors.
#   Null: RNA residual columns permuted within well. Seeds: 20260923 + offsets
#   (sampling: +k, test permutations: +1000+k, control permutations: +2000+k,
#   k = lineage index in sorted order; L2: +3000).
#   Background: r of the same gene with 200 peaks on other chromosomes matched on
#   mean logCPM (20 bins) and width (3 bins); p_bg = (1 + #|r_bg| >= |r|)/201.
#   BH within each family (H24, REF_PNPLA3, REF_ZHU) across lineages.
#   L1 support: BH q < 0.05, p_bg < 0.05, lineage informative (own-promoter
#   nominal positive rate >= 2 x trans nominal rate).
# L2, GSE244832 (donor-matched, different aliquots), only for L1-supported pairs
#   in lineages with ATAC: hepatocyte, fibroblast->stellate, macrophage,
#   cholangiocyte, t_cell->t_nk. Donors with >= 20 nuclei in both modalities,
#   >= 12 donors, same filters and TMM logCPM, no covariates, 10,000
#   unrestricted permutations, BH over the L2-tested pairs of each family;
#   support: q < 0.05 and the same sign as L1.
# Usage: Rscript 05_links.R <OUT_DIR> <N_PERM>
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer); library(edgeR)
  library(jsonlite); library(Matrix)
})
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[[1]]; N_PERM <- as.integer(args[[2]])
SEED <- 20260923L; MIN_NUC <- 20L; MIN_DONORS <- 15L; N_CTRL <- 2000L; N_PERM_CTRL <- 10000L; N_BG <- 200L
N_PERM_L2 <- 10000L; MIN_DONORS_L2 <- 12L; ATAC_MIN_COUNT <- 5; ATAC_MIN_FRAC <- 0.5
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
GTF <- "/gpfs/commons/home/jameslee/reference_genome/cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
PB <- file.path(ROOT, "Analysis/MASLD_Model_Benchmark/executions/gse296875-donor-lineage-pseudobulk-21107249/pseudobulk")
ATAC244 <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/counts/GSE244832")
L2MAP <- c(hepatocyte = "hepatocyte", fibroblast = "stellate", macrophage = "macrophage",
           cholangiocyte = "cholangiocyte", t_cell = "t_nk")
res_dir <- file.path(OUT, "results"); dir.create(res_dir, showWarnings = FALSE)
atac_keep <- function(m) rowSums(m >= ATAC_MIN_COUNT) >= ATAC_MIN_FRAC * ncol(m)

read_f32 <- function(path, dims) matrix(readBin(path, "numeric", n = prod(dims), size = 4, endian = "little"), nrow = dims[1])
units <- fread(file.path(OUT, "gse296875/units.tsv"))
dims <- fromJSON(file.path(OUT, "gse296875/dims.json"))
rna <- read_f32(file.path(OUT, "gse296875/rna_gene_by_unit.f32"), dims$rna)
atac <- read_f32(file.path(OUT, "gse296875/atac_peak_by_unit.f32"), dims$atac)
genes <- fread(file.path(PB, "genes.tsv"))
peaks <- fread(file.path(OUT, "gse296875/peaks.tsv"))
pp <- tstrsplit(peaks$peak_id, "-")
pgr <- GRanges(pp[[1]], IRanges(as.integer(pp[[2]]), as.integer(pp[[3]])))
uni <- fread(file.path(OUT, "inventory_v2/l1_test_universe.tsv"))
gtf <- import(GTF, feature.type = "gene"); gtf$gene_id <- sub("\\..*$", "", gtf$gene_id)
tss <- resize(gtf, 1L, fix = "start"); tss <- tss[tss$gene_id %in% genes$ensembl_stable_id]
gene_chr <- setNames(as.character(seqnames(tss)), tss$gene_id)

std_resid <- function(mat, X) {
  r <- t(qr.resid(qr(X), t(mat)))
  r / sqrt(rowSums(r^2))
}
make_perm <- function(strata) {
  idx <- split(seq_along(strata), strata)
  function() {
    p <- seq_along(strata)
    for (i in idx) if (length(i) > 1L) p[i] <- i[sample.int(length(i))]
    p
  }
}
# Exceedance counts for paired rows (Re[j,], Ae[j,]) under permutations of the RNA columns.
perm_exceed <- function(Re, Ae, robs, nperm, perm_fun) {
  a <- abs(robs) - 1e-12; cnt <- integer(length(robs))
  for (b in seq_len(nperm)) {
    p <- perm_fun()
    cnt <- cnt + (abs(rowSums(Re[, p, drop = FALSE] * Ae)) >= a)
  }
  cnt
}

all_tests <- list(); diag <- list(); lin_order <- sort(unique(units$lineage_id))
for (k in seq_along(lin_order)) {
  lin <- lin_order[k]
  u <- units[lineage_id == lin & nuclei >= MIN_NUC]
  t_lin <- uni[lineage == lin]
  if (nrow(u) < MIN_DONORS) { diag[[lin]] <- data.table(lineage = lin, donors = nrow(u), status = "untestable_donors"); next }
  if (!nrow(t_lin)) { diag[[lin]] <- data.table(lineage = lin, donors = nrow(u), tests = 0L, status = "no_accessible_variant_peak"); next }
  cols <- u$unit_index + 1L
  X <- model.matrix(~ factor(well_id) + log(nuclei), data = u)
  y <- DGEList(rna[, cols]); kg <- filterByExpr(y, design = matrix(1, length(cols), 1))
  y <- calcNormFactors(y[kg, , keep.lib.sizes = FALSE]); lr <- cpm(y, log = TRUE, prior.count = 2)
  z <- DGEList(atac[, cols]); kp <- atac_keep(atac[, cols, drop = FALSE])
  z <- calcNormFactors(z[kp, , keep.lib.sizes = FALSE]); la <- cpm(z, log = TRUE, prior.count = 2)
  gid <- genes$ensembl_stable_id[kg]; pid <- which(kp)
  R <- std_resid(lr, X); A <- std_resid(la, X)
  g_row <- function(e) match(e, gid); p_row <- function(i) match(i, pid)
  perm_fun <- make_perm(u$well_id)

  # test pairs
  Re <- R[g_row(t_lin$ensembl), , drop = FALSE]; Ae <- A[p_row(t_lin$peak_index), , drop = FALSE]
  robs <- rowSums(Re * Ae)
  set.seed(SEED + 1000L + k)
  ex <- perm_exceed(Re, Ae, robs, N_PERM, perm_fun)
  t_lin[, `:=`(n_donors = length(cols), n_wells = uniqueN(u$well_id), r = robs, p_perm = (1 + ex) / (1 + N_PERM))]

  # accessibility- and width-matched background on other chromosomes
  set.seed(SEED + k)
  mean_la <- rowMeans(la); wid <- width(pgr[pid])
  abin <- cut(mean_la, unique(quantile(mean_la, 0:20 / 20)), include.lowest = TRUE, labels = FALSE)
  wbin <- cut(wid, unique(quantile(wid, 0:3 / 3)), include.lowest = TRUE, labels = FALSE)
  pchr <- as.character(seqnames(pgr[pid]))
  tp <- unique(t_lin$peak_index)
  bg_idx <- lapply(tp, function(p) {
    j <- p_row(p)
    cand <- which(abin == abin[j] & wbin == wbin[j] & pchr != pchr[j])
    cand[sample.int(length(cand), min(N_BG, length(cand)))]
  })
  names(bg_idx) <- tp
  t_lin[, p_bg := NA_real_]
  for (j in seq_len(nrow(t_lin))) {
    rb <- as.vector(A[bg_idx[[as.character(t_lin$peak_index[j])]], , drop = FALSE] %*% Re[j, ])
    t_lin$p_bg[j] <- (1 + sum(abs(rb) >= abs(robs[j]))) / (1 + length(rb))
  }
  all_tests[[lin]] <- t_lin

  # positive control: gene and its own promoter peak (TSS +/- 1 kb), genes outside the test universe
  tss_k <- tss[tss$gene_id %in% setdiff(gid, uni$ensembl)]
  ov <- findOverlaps(resize(tss_k, 2001L, fix = "center"), pgr[pid], ignore.strand = TRUE)
  prom <- unique(data.table(g = tss_k$gene_id[queryHits(ov)], p = pid[subjectHits(ov)]), by = "g")
  prom <- prom[sample.int(nrow(prom), min(N_CTRL, nrow(prom)))]
  # negative control: test peaks x random expressed genes on other chromosomes
  tpk <- rep(tp, length.out = N_CTRL)
  trans <- rbindlist(lapply(tpk, function(p) {
    chr <- as.character(seqnames(pgr[p])); cand <- gid[!is.na(gene_chr[gid]) & gene_chr[gid] != chr]
    data.table(g = cand[sample.int(length(cand), 1L)], p = p)
  }))
  set.seed(SEED + 2000L + k)
  ctrl_stats <- function(d) {
    Rc <- R[g_row(d$g), , drop = FALSE]; Ac <- A[p_row(d$p), , drop = FALSE]
    ro <- rowSums(Rc * Ac)
    e <- perm_exceed(Rc, Ac, ro, N_PERM_CTRL, perm_fun)
    data.table(r = ro, p = (1 + e) / (1 + N_PERM_CTRL))
  }
  ps <- ctrl_stats(prom); ts <- ctrl_stats(trans)
  fwrite(cbind(prom, ps), file.path(res_dir, paste0("ctrl_promoter_", lin, ".tsv")), sep = "\t")
  fwrite(cbind(trans, ts), file.path(res_dir, paste0("ctrl_trans_", lin, ".tsv")), sep = "\t")
  pos_rate <- mean(ps$p < 0.05 & ps$r > 0); trans_rate <- mean(ts$p < 0.05)
  diag[[lin]] <- data.table(lineage = lin, donors = length(cols), wells = uniqueN(u$well_id),
                            expressed_genes = length(gid), accessible_peaks = length(pid),
                            tests = nrow(t_lin), promoter_pairs = nrow(ps),
                            promoter_median_r = median(ps$r), promoter_pos_nominal_rate = pos_rate,
                            trans_pairs = nrow(ts), trans_nominal_rate = trans_rate,
                            trans_median_abs_r = median(abs(ts$r)),
                            l1_informative = pos_rate >= 2 * trans_rate,
                            status = "tested")
  gc()
  cat(sprintf("[%s] donors %d tests %d promoter+ %.3f trans %.3f\n", lin, length(cols), nrow(t_lin), pos_rate, trans_rate))
}
tests <- rbindlist(all_tests, fill = TRUE); diag <- rbindlist(diag, fill = TRUE)
tests[, q_bh := p.adjust(p_perm, "BH"), by = family]
inform <- diag[status == "tested" & l1_informative == TRUE, lineage]
tests[, l1_support := q_bh < 0.05 & p_bg < 0.05 & lineage %in% inform]
fwrite(diag, file.path(res_dir, "l1_lineage_diagnostics.tsv"), sep = "\t")

# ---- L2: GSE244832 confirmation of L1-supported pairs ----
set.seed(SEED + 3000L)
sup <- tests[l1_support == TRUE & lineage %in% names(L2MAP)]
l2 <- list()
if (nrow(sup)) {
  r244u <- fread(file.path(OUT, "gse244832/rna_units.tsv"))
  d244 <- fromJSON(file.path(OUT, "gse244832/dims.json"))
  r244 <- read_f32(file.path(OUT, "gse244832/rna_gene_by_unit.f32"), d244$rna)
  g244 <- fread(file.path(OUT, "gse244832/genes.tsv"))$ensembl_id
  for (lin in sort(unique(sup$lineage))) {
    al <- L2MAP[[lin]]
    dn <- fread(file.path(ATAC244, paste0(al, ".donors.tsv")))
    pk <- fread(file.path(ATAC244, paste0(al, ".peaks.tsv")))
    mm <- t(as(readMM(file.path(ATAC244, paste0(al, ".mtx.gz"))), "CsparseMatrix"))  # peaks x donors
    stopifnot(ncol(mm) == nrow(dn), nrow(mm) == nrow(pk))
    ru <- r244u[lineage_id == al]
    both <- merge(ru[, .(donor_id, rna_col = unit_index + 1L, rna_nuc = nuclei)],
                  data.table(donor_id = dn$donor_id, atac_col = seq_len(nrow(dn)), atac_nuc = dn$n_cells), by = "donor_id")
    both <- both[rna_nuc >= MIN_NUC & atac_nuc >= MIN_NUC]
    s <- sup[lineage == lin]
    if (nrow(both) < MIN_DONORS_L2) {
      l2[[lin]] <- s[, .(lineage, family, peak_id, ensembl, l2_status = "untestable_donors", l2_donors = nrow(both))]; next
    }
    y <- DGEList(r244[, both$rna_col]); kg <- filterByExpr(y, design = matrix(1, nrow(both), 1))
    y <- calcNormFactors(y[kg, , keep.lib.sizes = FALSE]); lr <- cpm(y, log = TRUE, prior.count = 2)
    am <- as.matrix(mm[, both$atac_col]); kp <- atac_keep(am)
    z <- calcNormFactors(DGEList(am[kp, , drop = FALSE]))
    lib <- z$samples$lib.size * z$samples$norm.factors
    pc <- tstrsplit(sub(":", "-", pk$peak_coordinate), "-")
    cgr <- GRanges(pc[[1]], IRanges(as.integer(pc[[2]]), as.integer(pc[[3]])))
    X1 <- matrix(1, nrow(both), 1)
    out <- rbindlist(lapply(seq_len(nrow(s)), function(j) {
      gr <- which(g244[kg] == s$ensembl[j])
      hit <- intersect(subjectHits(findOverlaps(pgr[s$peak_index[j]], cgr)), which(kp))
      if (!length(gr) || !length(hit)) return(data.table(lineage = lin, family = s$family[j], peak_id = s$peak_id[j], ensembl = s$ensembl[j],
                                                          l2_status = if (!length(gr)) "gene_not_expressed" else "no_accessible_peak",
                                                          l2_donors = nrow(both)))
      cnt <- colSums(am[hit, , drop = FALSE])
      av <- log2((cnt + 0.5) / (lib + 1) * 1e6)
      rv <- lr[gr[1], ]
      rr <- std_resid(rbind(rv), X1); aa <- std_resid(rbind(av), X1)
      ro <- sum(rr * aa)
      e <- sum(vapply(seq_len(N_PERM_L2), function(i) abs(sum(rr[sample.int(ncol(rr))] * aa)) >= abs(ro) - 1e-12, logical(1)))
      data.table(lineage = lin, family = s$family[j], peak_id = s$peak_id[j], ensembl = s$ensembl[j], l2_status = "tested",
                 l2_donors = nrow(both), l2_peaks = length(hit), l2_r = ro, l2_p = (1 + e) / (1 + N_PERM_L2))
    }), fill = TRUE)
    l2[[lin]] <- out
  }
}
l2 <- rbindlist(l2, fill = TRUE)
if (nrow(l2)) {
  if (!"l2_p" %in% names(l2)) l2[, `:=`(l2_r = NA_real_, l2_p = NA_real_)]
  l2[l2_status == "tested", l2_q := p.adjust(l2_p, "BH"), by = family]
  tests <- merge(tests, l2, by = c("lineage", "family", "peak_id", "ensembl"), all.x = TRUE)
  tests[, l2_support := !is.na(l2_q) & l2_q < 0.05 & sign(l2_r) == sign(r)]
} else tests[, `:=`(l2_status = NA_character_, l2_support = FALSE)]
fwrite(tests, file.path(res_dir, "l1_l2_link_tests.tsv"), sep = "\t")
write_json(list(seed = SEED, n_perm = N_PERM, tests = tests[, .N, by = family],
                l1_supported = tests[, .(n = sum(l1_support)), by = family],
                l2_tested = sum(tests$l2_status == "tested", na.rm = TRUE), l2_supported = sum(tests$l2_support),
                informative_lineages = inform), file.path(res_dir, "links_summary.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_links.txt"))
cat("links done\n")
