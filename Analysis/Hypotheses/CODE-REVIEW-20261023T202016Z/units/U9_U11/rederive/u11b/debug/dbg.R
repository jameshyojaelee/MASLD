#!/usr/bin/env Rscript
# U11b independent re-derivation of Figure 5E protein column (code review 2026-10-23).
# From METHODS.md "For Figure 5E ..." paragraph:
#   single-gene protein groups, <=50% missing, log2, quantile normalized, median per symbol,
#   z per protein over 58 samples; program score = L1-weighted mean of measured member z,
#   members once per primary-assembly gene; scorable: >=8 measured and >=20% L1 weight;
#   burden = PC1 of steatosis, ballooning, inflammation, fibrosis in 46 MASLD participants;
#   rank-residual partial Spearman adjusted for batch, age, BMI, sex; batch-stratified
#   Freedman-Lane permutation (9,999); participant bootstrap within batch (4,999);
#   BH over scorable programs; robust = q<0.05 + sign agreement (equal-weight, leave-top,
#   both batches). Seeds 20260923 / 30260923 + program index (position 1..94 among programs
#   scorable under the earlier confirmed-symbol membership; 1032 for Portal fibroblast).
# Parts: (1) same-seed reproduction incl. bootstrap CI; (2) 20 alternate permutation seed
# sets; (3) one reordering of the program index (reverse registry order).
suppressPackageStartupMessages({ library(data.table); library(limma); library(parallel) })
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out <- Sys.getenv("U11B_OUT"); stopifnot(nzchar(out)); dir.create(out, recursive = TRUE, showWarnings = FALSE)
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
B_PERM <- 9999L; B_BOOT <- 50L
inp <- c(
  raw = file.path(R, "data/PXD051911/liver_protein_quant.txt"),
  meta = file.path(R, "data/PXD051911/meta_data.txt"),
  membership = file.path(R, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv"),
  registry = file.path(R, "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv"),
  gencode = file.path(R, "data/gencode_v49_gene_metadata.tsv.gz"),
  printed = file.path(R, "figures/candidates/panel-deep-fixes-20261017T035633Z/fig5/work/protein_histology_117_primary.tsv"),
  printed_5E = file.path(R, "figures/candidates/panel-deep-fixes-20261017T035633Z/fig5/source_tables/5E.tsv")
)
fwrite(data.table(input = names(inp), path = inp,
                  sha256 = sapply(inp, function(f) sub(" .*", "", system2("sha256sum", f, stdout = TRUE)))),
       file.path(out, "input_sha256.tsv"), sep = "\t")
log <- function(...) cat(format(Sys.time(), "%H:%M:%S "), ..., "\n", sep = "")

# ---- protein matrix (58 samples, raw-file column order)
raw <- fread(inp[["raw"]]); md <- fread(inp[["meta"]])
cols <- setdiff(names(raw), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% cols], by = "liver_proteomics_filename")
md <- md[order(match(liver_proteomics_filename, cols))]
stopifnot(nrow(md) == 58L)
md[, `:=`(sid = liver_proteomics_filename, masld = saf_diagnosis != "No_MASLD",
          batch = ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020"),
          age = as.numeric(alder), bmi = as.numeric(bmi), sex = gender,
          Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)), Steatosis = as.numeric(steatosis_score),
          Ballooning = as.numeric(hepatocellular_ballooning_score), Inflammation = as.numeric(lobular_inflammation_score))]
r1 <- raw[!is.na(Genes) & nzchar(Genes) & !grepl(";", Genes, fixed = TRUE)]
E <- as.matrix(r1[, md$sid, with = FALSE]); storage.mode(E) <- "double"
E[!is.finite(E) | E <= 0] <- NA; E <- log2(E); rownames(E) <- r1$Genes
E <- normalizeBetweenArrays(E[rowMeans(is.na(E)) <= 0.5, ], method = "quantile")
G <- do.call(rbind, lapply(split(seq_len(nrow(E)), rownames(E)), function(i)
  if (length(i) == 1L) E[i, ] else apply(E[i, , drop = FALSE], 2, median, na.rm = TRUE)))
G[!is.finite(G)] <- NA
Zp <- t(apply(G, 1, function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)))
Zp[!is.finite(Zp)] <- NA
log("proteins ", nrow(Zp))

# ---- membership
reg <- fread(inp[["registry"]], select = c("cell_type", "module", "program_uid", "module_name"))
setorder(reg, cell_type, module); reg[, registry_row := .I]
mem <- fread(inp[["membership"]])
gc5 <- fread(inp[["gencode"]], select = c("gene_id", "gene_name", "chromosome"))
nprim <- gc5[chromosome %in% paste0("chr", c(1:22, "X", "Y", "M")), .N, by = gene_name]
mem[, sym_conf := fifelse(mapped_symbol_status == "gencode_v49_unique_symbol_confirmed", mapped_symbol, NA_character_)]
mem[, sym_prim := fifelse(mapped_symbol_status == "gencode_v49_unique_symbol_confirmed", mapped_symbol,
                   fifelse(mapped_symbol_status == "unmapped_or_ambiguous" & canonical_gene %in% nprim[N == 1L, gene_name],
                           canonical_gene, NA_character_))]
ens_named <- gc5$gene_name[match(mem$mapped_symbol, sub("\\..*$", "", gc5$gene_id))]
log("Ensembl-ID members whose GENCODE name is a measured protein: ",
    sum(mem$mapped_symbol_status == "gencode_v49_unambiguous_ensembl_to_symbol" & ens_named %in% rownames(Zp), na.rm = TRUE))

score_programs <- function(symcol) {
  rbindlist(lapply(reg$program_uid, function(p) {
    m <- mem[program_uid == p & !is.na(get(symcol)), .(w = sum(original_l1_weight)), by = .(g = get(symcol))]
    ms <- m[g %in% rownames(Zp)]
    ok <- nrow(ms) >= 8L && sum(ms$w) >= 0.20
    if (!ok) return(data.table(program_uid = p, n_measured = nrow(ms), retained = sum(ms$w), scorable = FALSE))
    w <- ms$w / sum(ms$w); Zm <- Zp[ms$g, md$sid, drop = FALSE]; fin <- is.finite(Zm)
    wsc <- colSums(ifelse(fin, Zm * w, 0)) / colSums(fin * w)
    esc <- colSums(ifelse(fin, Zm, 0)) / colSums(fin)
    top <- ms$g[which.max(ms$w)]; ml <- ms[g != top]; wl <- ml$w / sum(ml$w)
    Zl <- Zp[ml$g, md$sid, drop = FALSE]; fl <- is.finite(Zl)
    lsc <- colSums(ifelse(fl, Zl * wl, 0)) / colSums(fl * wl)
    data.table(program_uid = p, n_measured = nrow(ms), retained = sum(ms$w), scorable = TRUE,
               sid = md$sid, score = wsc, equal = esc, leave = lsc)
  }), fill = TRUE)
}
S_conf <- score_programs("sym_conf"); S_prim <- score_programs("sym_prim")
sc_conf <- unique(S_conf[scorable == TRUE, program_uid]); sc_prim <- unique(S_prim[scorable == TRUE, program_uid])
log("scorable: confirmed-only ", length(sc_conf), "; primary ", length(sc_prim))

# ---- seed index from Methods: position among confirmed-only scorable programs in registry order
pos <- reg[program_uid %in% sc_conf, .(program_uid, pos = seq_len(.N))]
seed_idx <- merge(reg[program_uid %in% sc_prim, .(program_uid, registry_row)], pos, by = "program_uid", all.x = TRUE)
seed_idx[is.na(pos), pos := 1000L + registry_row]
printed <- fread(inp[["printed"]])
chk_seed <- merge(seed_idx, printed[testable == TRUE, .(program_uid, perm_seed, boot_seed)], by = "program_uid")
log("seeds from Methods equal stored seeds: perm ", all(chk_seed$perm_seed == 20260923L + chk_seed$pos),
    " boot ", all(chk_seed$boot_seed == 30260923L + chk_seed$pos), " (n=", nrow(chk_seed), ")")

# ---- analysis set: 46 MASLD participants with complete histology
comp <- c("Steatosis", "Ballooning", "Inflammation", "Fibrosis")
mm <- md[masld == TRUE & complete.cases(md[, ..comp])]
log("MASLD participants with complete histology ", nrow(mm), " (MASLD total ", sum(md$masld), ")")
burden <- function(d) {
  z <- scale(as.matrix(d[, ..comp])); if (any(!is.finite(z))) return(NULL)
  pc <- prcomp(z, center = FALSE, scale. = FALSE); l <- pc$rotation[, 1]; s <- pc$x[, 1]
  if (sum(l) < 0) { l <- -l; s <- -s }
  list(score = as.numeric(s), loading = l)
}
bf <- burden(mm); stopifnot(all(bf$loading > 0)); mm[, hb := bf$score]
log("burden loadings: ", paste(names(bf$loading), round(bf$loading, 4), collapse = " "))
Zn <- function(d, with_batch = TRUE) {
  X <- cbind(1, d$age, d$bmi)
  if (with_batch && length(unique(d$batch)) > 1) X <- cbind(X, d$batch == "2020")
  if (length(unique(d$sex)) > 1) X <- cbind(X, d$sex == unique(d$sex)[1])
  X
}
pspear <- function(x, y, d, with_batch = TRUE) {
  ok <- is.finite(x) & is.finite(y); X <- Zn(d[ok], with_batch)
  if (sum(ok) <= ncol(X) + 2L) return(NA_real_)
  Q <- qr(X); cor(qr.resid(Q, rank(x[ok])), qr.resid(Q, rank(y[ok])))
}
# Freedman-Lane, draws mirror sample(ii, length(ii)) within batch in batch-level order
fl_p <- function(x, d, seed, B = B_PERM) {
  X <- Zn(d); Q <- qr(X)
  xr <- qr.resid(Q, rank(x)); yrank <- rank(d$hb); yr <- qr.resid(Q, yrank); fit0 <- yrank - yr
  obs <- cor(xr, yr)
  bb <- split(seq_along(yr), factor(d$batch))
  set.seed(seed)
  P <- matrix(0L, length(yr), B)
  for (i in seq_len(B)) { idx <- seq_along(yr); for (ii in bb) idx[ii] <- ii[sample.int(length(ii))]; P[, i] <- idx }
  Ys <- qr.resid(Q, fit0 + matrix(yr[P], nrow = length(yr)))
  st <- as.numeric(cor(xr, Ys))
  c(obs = obs, p = (sum(abs(st) >= abs(obs)) + 1) / (B + 1))
}
boot_ci <- function(x, d, seed, B = B_BOOT) {
  set.seed(seed); bb <- split(seq_len(nrow(d)), d$batch)
  r <- vapply(seq_len(B), function(b) {
    idx <- unlist(lapply(bb, function(ii) ii[sample.int(length(ii), length(ii), replace = TRUE)]), use.names = FALSE)
    db <- copy(d[idx]); f <- burden(db)
    if (is.null(f) || any(f$loading <= 0)) return(NA_real_)
    db[, hb := f$score]; pspear(x[idx], db$hb, db)
  }, numeric(1))
  c(lo = quantile(r, 0.025, na.rm = TRUE, names = FALSE), hi = quantile(r, 0.975, na.rm = TRUE, names = FALSE), nfin = sum(is.finite(r)))
}

# ---- per-program data with explicit merge row counts
progs <- sc_prim
sp <- S_prim[scorable == TRUE]
n_before <- sp[program_uid == progs[1] & sid %in% md[masld == TRUE, sid], .N]
an <- merge(sp, mm[, c("sid", "batch", "age", "bmi", "sex", "hb", comp), with = FALSE], by = "sid", all = FALSE)
nper <- an[, .N, by = program_uid]
log("merge: MASLD score rows per program ", n_before, "; analysis rows per program range ", min(nper$N), "-", max(nper$N),
    "; programs ", nrow(nper))
fwrite(nper, file.path(out, "merge_rows_per_program.tsv"), sep = "\t")

one_program <- function(p, perm_seed, boot_seed = NA, do_boot = FALSE) {
  d <- an[program_uid == p][order(match(sid, md$sid))]
  ok <- is.finite(d$score); d <- d[ok]
  f <- fl_p(d$score, d, perm_seed)
  res <- data.table(program_uid = p, n = nrow(d), rho = pspear(d$score, d$hb, d), p = f[["p"]],
                    equal_rho = pspear(d$equal, d$hb, d), leave_rho = pspear(d$leave, d$hb, d),
                    rho_2019 = pspear(d[batch == "2019"]$score, d[batch == "2019"]$hb, d[batch == "2019"], FALSE),
                    rho_2020 = pspear(d[batch == "2020"]$score, d[batch == "2020"]$hb, d[batch == "2020"], FALSE))
  if (do_boot) { b <- boot_ci(d$score, d, boot_seed); res[, `:=`(ci_low = b[["lo"]], ci_high = b[["hi"]], n_boot_finite = b[["nfin"]])] }
  res
}
call_set <- function(r) {
  r[, q := p.adjust(p, "BH")]
  r[, robust := q < 0.05 & sign(rho) == sign(equal_rho) & sign(rho) == sign(leave_rho) &
      sign(rho) == sign(rho_2019) & sign(rho) == sign(rho_2020)]
  r
}
seed_idx <- seed_idx[match(progs, program_uid)]
r <- tryCatch(one_program(progs[1], 20260923L + seed_idx$pos[1], 30260923L + seed_idx$pos[1], TRUE), error = function(e) {print(e); traceback(); NULL})
print(r)
d <- an[program_uid == progs[1]]; str(d)
