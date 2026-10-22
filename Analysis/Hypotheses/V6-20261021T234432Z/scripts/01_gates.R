# V6 job 1: rebuild GeoMx counts from DCC + PKC, gates K0 (rebuild vs deposited),
# K3 (stain orientation), K4 (LOQ detection and program coverage). No program score.
# Follows PRESPEC_V6.md. Usage: Rscript 01_gates.R
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/V6-20261021T234432Z")
DL <- file.path(RUN, "downloads_untrusted")
OUT <- file.path(RUN, "outputs")
MEMBERSHIP <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv")
PROGRAMS <- c(stromal_ecm = "hotspot_hepatocytes_f05c535ae5bbc0b9", ductular = "hotspot_hepatocytes_48f39dd4d817a10e")
K3_SETS <- list(Albumin = c("ALB", "APOA1"), Vinculin = c("COL1A1", "COL1A2", "COL3A1"), CD45 = "PTPRC")
MIN_GENES <- 8; MIN_WEIGHT <- 0.80; MIN_DET_AOIS <- 5

if (dir.exists(OUT)) stop("outputs/ exists; refusing to overwrite")
sha <- function(f) strsplit(system2("sha256sum", f, stdout = TRUE), " ")[[1]][1]
if (sha(file.path(RUN, "PRESPEC_V6.md")) != strsplit(readLines(file.path(RUN, "PRESPEC_V6.md.sha256")), " ")[[1]][1])
  stop("STOP: prespec hash differs")
if (sha(MEMBERSHIP) != "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b") stop("membership hash")
dir.create(OUT)
wr <- function(x, f) fwrite(as.data.table(x), file.path(OUT, f), sep = "\t", na = "NA")

# ---- AOI metadata from SOFT ------------------------------------------------------------
soft <- readLines(gzfile(file.path(DL, "GSE255312_family.soft.gz")))
starts <- grep("^\\^SAMPLE = ", soft)
ends <- c(starts[-1] - 1, length(soft))
val <- function(block, key) sub(paste0("^", key, " = "), "", grep(paste0("^", key, " = "), block, value = TRUE))
aoi <- rbindlist(lapply(seq_along(starts), function(i) {
  b <- soft[starts[i]:ends[i]]
  ch <- val(b, "!Sample_characteristics_ch1")
  desc <- val(b, "!Sample_description")
  data.table(gsm = sub("^\\^SAMPLE = ", "", b[1]), title = val(b, "!Sample_title"),
             stain = sub("^stained with: ", "", grep("^stained with", ch, value = TRUE)),
             disease = sub("^disease state: ", "", grep("^disease state", ch, value = TRUE)),
             dcc = grep("^DSP-", desc, value = TRUE))
}))
parts <- tstrsplit(aoi$title, " | ", fixed = TRUE)
aoi[, slide := parts[[1]]][, roi_label := parts[[2]]][, roi := paste0(slide, sub("[A-Z]$", "", roi_label))]
if (nrow(aoi) != 72 || any(parts[[3]] != aoi$stain)) stop("STOP: SOFT parse")
if (!all(aoi[, .N, by = roi]$N == 3) || uniqueN(aoi$roi) != 24) stop("STOP: ROIs do not have 3 AOIs each")
if (!all(aoi[, uniqueN(stain), by = roi]$V1 == 3)) stop("STOP: ROI without one AOI per stain")

# ---- PKC and DCC ------------------------------------------------------------------------
pkc <- fromJSON(gzfile(file.path(DL, "GSE255312_Hs_R_NGS_WTA_v1.0.pkc.gz")), simplifyVector = FALSE)
probe_map <- rbindlist(lapply(pkc$Targets, function(t) rbindlist(lapply(t$Probes, function(p)
  data.table(rts = p$RTS_ID, target = t$DisplayName, code_class = t$CodeClass)))))
read_dcc <- function(f) {
  x <- readLines(gzfile(f))
  cs <- x[(grep("<Code_Summary>", x) + 1):(grep("</Code_Summary>", x) - 1)]
  cs <- cs[nzchar(cs)]
  aligned <- as.numeric(sub("^Aligned,", "", grep("^Aligned,", x, value = TRUE)))
  id <- sub("^ID,", "", grep("^ID,", x, value = TRUE))
  list(id = id, aligned = aligned,
       counts = data.table(rts = sub(",.*", "", cs), count = as.numeric(sub(".*,", "", cs))))
}
files <- list.files(file.path(DL, "raw"), pattern = "\\.dcc\\.gz$", full.names = TRUE)
dccs <- lapply(files, read_dcc)
ids <- vapply(dccs, `[[`, "", "id")
if (!setequal(ids, aoi$dcc)) stop("STOP: DCC IDs differ from SOFT descriptions")
aoi[, aligned := vapply(dccs, `[[`, 0, "aligned")[match(dcc, ids)]]
setorder(aoi, slide, roi, stain)
P <- matrix(0, nrow(probe_map), nrow(aoi), dimnames = list(probe_map$rts, aoi$dcc))
for (d in dccs) {
  cnt <- d$counts[rts %in% probe_map$rts]
  P[cnt$rts, d$id] <- cnt$count
}
unknown <- sum(vapply(dccs, function(d) sum(!d$counts$rts %in% probe_map$rts), 0))

geomean <- function(x) exp(mean(log(pmax(x, 1))))
endo <- probe_map[code_class == "Endogenous01"]
if (anyDuplicated(endo$target)) {
  G <- t(sapply(split(seq_len(nrow(endo)), endo$target), function(i) apply(P[endo$rts[i], , drop = FALSE], 2, geomean)))
} else {
  G <- P[endo$rts, , drop = FALSE]; rownames(G) <- endo$target
}
neg <- P[probe_map[code_class == "Negative01", rts], , drop = FALSE]
q3 <- apply(G, 2, quantile, 0.75)
L <- log2(sweep(G, 2, q3, "/") * geomean(q3) + 1)
loq <- apply(neg, 2, function(x) max(2, geomean(x) * exp(sd(log(pmax(x, 1))))^2))
aoi[, `:=`(q3 = q3[dcc], neg_geomean = apply(neg, 2, geomean)[dcc], loq = loq[dcc])]

# ---- K0 ---------------------------------------------------------------------------------
dep <- fread(cmd = paste("zcat", shQuote(file.path(DL, "GSE255312_GeoMx_analysis_log_Q3.csv.gz"))))
dm <- as.matrix(dep[, -1]); rownames(dm) <- dep$Gene; colnames(dm) <- sub("\\.dcc$", "", colnames(dm))
sg <- intersect(rownames(dm), rownames(L))
sa <- intersect(colnames(dm), colnames(L))
R <- cor(L[sg, sa], dm[sg, sa])
best <- colnames(R)[apply(R, 1, which.max)]
k0_identity <- length(sa) == 72 && all(best == rownames(R))
k0_median_r <- median(diag(R))
wr(data.table(dcc = rownames(R), r_own = diag(R), best_match = best, r_best = apply(R, 1, max)), "K0_rebuild_vs_deposited.tsv")

# ---- K3 ---------------------------------------------------------------------------------
Z <- t(scale(t(L)))
k3 <- rbindlist(lapply(names(K3_SETS), function(target) {
  genes <- intersect(K3_SETS[[target]], rownames(Z))
  m <- colMeans(Z[genes, aoi$dcc, drop = FALSE])
  a <- copy(aoi)[, mz := m]
  a[, .(target = target, genes = paste(genes, collapse = ";"),
        ok = mz[stain == target] > max(mz[stain != target])), by = roi]
}))
k3_sum <- k3[, .(n_roi_ok = sum(ok), n_roi = .N), by = .(target, genes)][, pass := n_roi_ok >= 20]
wr(k3, "K3_orientation_by_roi.tsv"); wr(k3_sum, "K3_orientation.tsv")

# ---- K4 ---------------------------------------------------------------------------------
av <- aoi[stain %in% c("Albumin", "Vinculin"), dcc]
det_n <- rowSums(sweep(G[, av], 2, loq[av], ">"))
det <- data.table(gene = rownames(G), n_aoi_above_loq = det_n, frac_above_loq = det_n / length(av),
                  mean_log2q3 = rowMeans(L[, av]), detected = det_n >= MIN_DET_AOIS)
wr(det, "K4_gene_detection_alb_vin.tsv")
mem <- fread(MEMBERSHIP)[mapped_symbol_status == "gencode_v49_unique_symbol_confirmed"]
cov <- rbindlist(lapply(names(PROGRAMS), function(p) {
  w <- mem[program_uid == PROGRAMS[[p]], .(gene = mapped_symbol, weight = original_l1_weight)]
  w <- merge(w, det, by = "gene", all.x = TRUE)
  w[, on_axis := !is.na(detected)]
  wr(w[, program := p], sprintf("K4_members_%s.tsv", p))
  data.table(program = p, uid = PROGRAMS[[p]], n_members = nrow(w),
             axis_genes = sum(w$on_axis), axis_weight = sum(w[on_axis == TRUE, weight]),
             detected_genes = sum(w$detected, na.rm = TRUE), detected_weight = sum(w[detected == TRUE, weight]))
}))
cov[, testable := detected_genes >= MIN_GENES & detected_weight >= MIN_WEIGHT]
wr(cov, "K4_coverage.tsv")

# ---- save derived and gates ------------------------------------------------------------
dir.create(file.path(RUN, "derived"), showWarnings = FALSE)
saveRDS(list(G = G, L = L, neg = neg, aoi = aoi), file.path(RUN, "derived/geomx_rebuilt.rds"))
wr(aoi, "aoi_metadata.tsv")
gates <- data.table(
  gate = c("n_aoi", "n_roi", "unknown_rts_reads_rows", "K0_shared_genes", "K0_shared_aois", "K0_identity_pass",
          "K0_median_r", "K0_min_r", "K0_r_ge_0.999_all", "K3_pass", "K4_stromal_ecm_testable", "K4_ductular_testable"),
  value = c(nrow(aoi), uniqueN(aoi$roi), unknown, length(sg), length(sa), k0_identity, signif(k0_median_r, 6),
            signif(min(diag(R)), 6), all(diag(R) >= 0.999), all(k3_sum$pass),
            cov[program == "stromal_ecm", testable], cov[program == "ductular", testable]))
wr(gates, "gates_K0_K4.tsv")
print(gates); print(k3_sum); print(cov)
if (!k0_identity) { writeLines("K0 identity failed", file.path(OUT, "STOP_K0.txt")); stop("STOP K0") }
if (!all(k3_sum$pass)) { writeLines("K3 orientation failed", file.path(OUT, "STOP_K3.txt")); stop("STOP K3") }
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_job1.txt"))
