#!/usr/bin/env Rscript
# H9: published MASLD locus classes (Jamialahmadi 2024) versus the 117 frozen Hotspot programs
# and their cell types. Implements PRESPEC_H9.md, written before any overlap was computed.
# Provisional on the unadopted COLOC release susie_coloc_r2_20260924T053244Z.
suppressPackageStartupMessages({library(data.table); library(readxl); library(jsonlite)})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
H9   <- file.path(ROOT, "Analysis/Hypotheses/H9-20261009T140001Z")
stopifnot(file.exists(file.path(H9, "PRESPEC_H9.md")))
OUT  <- file.path(H9, "results", paste0("run-", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")))
if (dir.exists(OUT)) stop("Refusing to overwrite ", OUT)
dir.create(OUT, recursive = TRUE)
SEED <- 20260923L
B <- 1000L
cat("H9 start", format(Sys.time(), tz = "UTC", usetz = TRUE), "\nout:", OUT, "\n")

XL  <- file.path(H9, "downloads_untrusted/jamialahmadi2024_nat_med/jamialahmadi2024_MOESM3.xlsx")
HOT <- file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot")
HM  <- file.path(ROOT, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
H10 <- file.path(ROOT, "Analysis/Hypotheses/H10-20261007T210859Z/loci.tsv")
GMT <- file.path(ROOT, "Analysis/downstream_analysis/pathway_analysis/data/msigdb.v2025.1.Hs.symbols.gmt")

checks <- data.table(check = character(), pass = logical(), detail = character())
chk <- function(name, ok, detail = "", fatal = TRUE) {
  checks <<- rbind(checks, data.table(check = name, pass = isTRUE(ok), detail = as.character(detail)))
  fwrite(checks, file.path(OUT, "checks.tsv"), sep = "\t")
  cat(if (isTRUE(ok)) "PASS" else "FAIL", name, detail, "\n")
  if (!isTRUE(ok) && fatal) stop("check failed: ", name, " ", detail)
}
sha <- function(f) sub(" .*", "", system2("sha256sum", shQuote(f), stdout = TRUE))

# ---- 0. input identity --------------------------------------------------------------------
chk("xlsx sha256", startsWith(sha(XL), "56e39d6b"), sha(XL))
chk("H10 loci sha256", startsWith(sha(H10), "7aaa7658"), sha(H10))
rd <- function(s, skip) as.data.table(suppressMessages(read_excel(XL, sheet = s, skip = skip, col_types = "text")))

# ---- 1. class derivation (PRESPEC section 3) ---------------------------------------------
s5 <- as.data.table(suppressMessages(read_excel(XL, sheet = "S5", col_names = FALSE, col_types = "text")))
hdr5 <- as.character(unlist(s5[3, .SD, .SDcols = 1:16]))
chk("S5 header row", identical(hdr5, c("Trait", "CHR", "POS", "SNP", "Locus", "A1", "A2", "A1Freq", "Beta", "P", "Beta", "P", "Beta", "P", "Beta", "P")),
    paste(hdr5, collapse = "|"))
s5 <- s5[-(1:3), .SD, .SDcols = 1:16]
setnames(s5, c("trait", "chr", "pos", "snp", "locus", "a1", "a2", "a1freq",
               "b_bmi", "p_bmi", "b_na", "p_na", "b_vat", "p_vat", "b_wfm", "p_wfm"))
s5 <- s5[!is.na(snp)]
chk("S5 rows: 37 PDFF and 18 cT1", s5[, sum(trait == "PDFF")] == 37L && s5[, sum(trait == "Liver iron corrected T1")] == 18L,
    paste(names(table(s5$trait)), table(s5$trait), collapse = ", "))
s5[, `:=`(chr = as.integer(chr), pos = as.integer(pos))]
numc <- c("a1freq", "b_bmi", "p_bmi", "b_na", "p_na", "b_vat", "p_vat", "b_wfm", "p_wfm")
s5[, (numc) := lapply(.SD, as.numeric), .SDcols = numc]
pm <- as.matrix(s5[, .(p_bmi, p_na, p_vat, p_wfm)]); bm <- as.matrix(s5[, .(b_bmi, b_na, b_vat, b_wfm)])
pm[is.na(pm)] <- Inf
best <- max.col(-pm, ties.method = "first")
s5[, `:=`(adj_best = c("BMI", "NA", "VAT", "WFM")[best],
          pdff_beta = bm[cbind(seq_len(nrow(bm)), best)],
          pdff_p = pm[cbind(seq_len(nrow(pm)), best)],
          pdff_sign_all_agree = apply(bm, 1, function(x) length(unique(sign(x[!is.na(x)]))) == 1L))]

s15 <- rd("S15", 1)
s15 <- s15[!is.na(SNP) & nzchar(SNP), .(snp = SNP, tg_a1 = A1, tg_beta = as.numeric(Beta), tg_p = as.numeric(P))]
chk("S15 has 32 variants", nrow(s15) == 32L, nrow(s15))
s17 <- rd("S17", 1)
s17v <- unique(s17[!is.na(SNP), .(snp = SNP, gene_label = Locus, chr = as.integer(CHR), pos = as.integer(POS))])
chk("S17 has 32 variants", nrow(s17v) == 32L, nrow(s17v))

u10 <- function(s) { g <- trimws(unlist(strsplit(na.omit(rd(s, 1)$Genes), ";"))); unique(g[nzchar(g)]) }
u10C <- u10("S10C"); u10D <- u10("S10D")

cand <- merge(s5[trait == "PDFF"], s15, by = "snp")
chk("S5 and S15 allele A1 agree", all(cand$a1 == cand$tg_a1), nrow(cand))
cand <- cand[tg_p < 0.05]
cand[, `:=`(class = fifelse(sign(pdff_beta) != sign(tg_beta), "discordant", "concordant"), rule = "PDFF row of S5 x S15")]
ct1 <- merge(s5[trait != "PDFF" & snp %in% s17v$snp & !snp %in% cand$snp], s15, by = "snp")
ct1[, class := fifelse(locus %in% u10C, "discordant", fifelse(locus %in% u10D, "concordant", NA_character_))]
ct1 <- ct1[!is.na(class)][, rule := "cT1 row of S5, gene in S10C/S10D union"]
cls <- rbind(cand, ct1, fill = TRUE)
cls <- cls[, .(snp, chr, pos, gene_label = locus, a1, class, rule, s5_trait = trait, adj_best,
               s5_beta = pdff_beta, s5_p = pdff_p, s5_signs_agree_all_adjustments = pdff_sign_all_agree, tg_beta, tg_p)]
chk("every classed variant is in S17", all(cls$snp %in% s17v$snp), paste(setdiff(cls$snp, s17v$snp), collapse = ","))

s9 <- rbind(rd("S9A", 1)[, src := "S9A"], rd("S9B", 1)[, src := "S9B"])
setnames(s9, c("Symbol (GRCh38)", "Locus (GRCh37)"), c("sym", "label"))
s9[, rk := as.numeric(Rank)]
paper_genes <- function(snp) {
  x <- s9[startsWith(label, paste0(snp, " "))]
  if (!nrow(x)) return(character())
  unique(x[, .SD[rk == max(rk)], by = .(src, Trait)]$sym)
}
cls[, genes_paper := vapply(snp, function(s) paste(sort(paper_genes(s)), collapse = ","), "")]
chk("every classed variant has paper-mapped genes", all(nzchar(cls$genes_paper)), paste(cls[!nzchar(genes_paper), snp], collapse = ","))
cls[, unit_genes := strsplit(genes_paper, ",")]
chk("10 discordant loci", sum(cls$class == "discordant") == 10L, sum(cls$class == "discordant"))
chk("13 concordant loci", sum(cls$class == "concordant") == 13L, sum(cls$class == "concordant"))
gD <- unique(unlist(cls[class == "discordant", unit_genes])); gC <- unique(unlist(cls[class == "concordant", unit_genes]))
chk("S10C genes within discordant paper-mapped genes", all(u10C %in% gD), paste(setdiff(u10C, gD), collapse = ","))
chk("S10D genes within concordant paper-mapped genes", all(u10D %in% gC), paste(setdiff(u10D, gC), collapse = ","))
fwrite(cls[order(class, chr, pos), !"unit_genes"], file.path(OUT, "classes.tsv"), sep = "\t")

# ---- 2. COLOC-anchored loci (arm A) -------------------------------------------------------
h10 <- fread(H10)
h10[, mhc := as.character(mhc) %in% c("True", "TRUE", "true", "T", "1")]
chk("3 MHC loci", sum(h10$mhc) == 3L, paste(h10[mhc == TRUE, locus_id], collapse = ","))
chk("289 COLOC loci", nrow(h10) == 289L, nrow(h10))
anch <- unique(unlist(strsplit(h10$genes, ",")))
chk("436 anchored genes", length(anch) == 436L, length(anch))
h10[, gene_list := strsplit(genes, ",")]

map_variants <- function(window) {
  rbindlist(lapply(seq_len(nrow(cls)), function(i) {
    v <- cls[i]; gs <- unique(c(v$unit_genes[[1]], v$gene_label))
    L <- h10[as.character(chr) == as.character(v$chr)]
    if (!nrow(L)) return(NULL)
    d <- pmax(0, pmax(L$start - v$pos, v$pos - L$end))
    byg <- vapply(L$gene_list, function(g) any(g %in% gs), logical(1))
    k <- which(d <= window | byg)
    if (!length(k)) return(NULL)
    data.table(snp = v$snp, gene_label = v$gene_label, class = v$class, locus_id = L$locus_id[k], locus_genes = L$genes[k],
               dist_bp = d[k], via = fifelse(d[k] <= window, "window", "gene"), mhc = L$mhc[k])
  }), fill = TRUE)
}
arm_units <- list(); map_tabs <- list()
for (w in c(100000L, 500000L)) {
  nm <- paste0("A", w / 1000L)
  m <- map_variants(w); m[, arm := nm]; map_tabs[[nm]] <- m
  mm <- m[mhc == FALSE]
  u <- mm[, .(classes = paste(sort(unique(class)), collapse = "+"), variants = paste(unique(snp), collapse = ","), genes = locus_genes[1]), by = locus_id]
  cat(nm, ": mapped COLOC loci", nrow(u), "| ambiguous (both classes):", sum(grepl("\\+", u$classes)),
      "| MHC-excluded loci:", uniqueN(m[mhc == TRUE, locus_id]), "\n")
  u <- u[!grepl("\\+", classes)]
  arm_units[[nm]] <- data.table(unit_id = u$locus_id, class = u$classes, variants = u$variants, genes = strsplit(u$genes, ","))
}
fwrite(rbindlist(map_tabs, fill = TRUE), file.path(OUT, "variant_to_coloc_locus.tsv"), sep = "\t")
arm_units[["B"]] <- data.table(unit_id = cls$snp, class = cls$class, variants = cls$snp, genes = cls$unit_genes)
arm_summary <- rbindlist(lapply(names(arm_units), function(a) {
  u <- arm_units[[a]]
  data.table(arm = a, n_units_discordant = sum(u$class == "discordant"), n_units_concordant = sum(u$class == "concordant"),
             n_genes_discordant = uniqueN(unlist(u[class == "discordant", genes])), n_genes_concordant = uniqueN(unlist(u[class == "concordant", genes])))
}))
arm_summary[, `:=`(tested_discordant = n_units_discordant >= 3L, tested_concordant = n_units_concordant >= 3L)]
fwrite(arm_summary, file.path(OUT, "arm_summary.tsv"), sep = "\t")
print(arm_summary)
fwrite(rbindlist(lapply(names(arm_units), function(a) arm_units[[a]][, .(arm = a, unit_id, class, variants, genes = vapply(genes, paste, "", collapse = ","))])),
       file.path(OUT, "arm_units.tsv"), sep = "\t")

# ---- 3. programs, universes, incidence matrices ------------------------------------------
reg <- fread(file.path(HOT, "program_registry_v2.tsv"), colClasses = "character")
mem <- fread(file.path(HOT, "program_membership_v2.tsv"), colClasses = "character")
cts <- c("hepatocytes", "fibroblasts", "cholangiocytes", "macrophages", "tcells")
chk("117 programs", nrow(reg) == 117L && uniqueN(reg$program_uid) == 117L, nrow(reg))
chk("program cell-type counts", identical(as.integer(table(factor(reg$cell_type, levels = cts))), c(30L, 29L, 28L, 19L, 11L)),
    paste(table(reg$cell_type), collapse = "/"))
reg[, pid := paste0(cell_type, "_", module)]; mem[, pid := paste0(cell_type, "_", module)]
chk("membership programs equal registry programs", setequal(reg$pid, mem$pid) && uniqueN(reg$pid) == 117L)
prog_meta <- reg[, .(pid, cell_type, module, module_name, K = as.integer(n_source_genes))]
Kmem <- mem[, .(Kobs = uniqueN(source_gene)), by = pid]
chk("program sizes equal n_source_genes", all(merge(prog_meta, Kmem, by = "pid")[, K == Kobs]))
ct_of_prog <- match(prog_meta$cell_type, cts)
univ <- lapply(cts, function(ct) unique(fread(file.path(HM, ct, "autocorr.tsv"), select = "gene")$gene)); names(univ) <- cts
N <- vapply(univ, length, 1L)
chk("universe sizes", identical(unname(N), c(12969L, 10675L, 12128L, 10663L, 9364L)), paste(N, collapse = "/"))
memb_in_univ <- mem[, mean(source_gene %in% univ[[cell_type[1]]]), by = cell_type]
chk("all program members are in their cell-type universe", all(memb_in_univ$V1 == 1), paste(round(memb_in_univ$V1, 4), collapse = "/"), fatal = FALSE)
Ukeys <- lapply(cts, function(ct) unique(c(univ[[ct]], mem[cell_type == ct, c(mapped_symbol, canonical_gene)])))
Kct <- vapply(cts, function(ct) uniqueN(mem[cell_type == ct, source_gene]), 1L)
prog_keys <- lapply(prog_meta$pid, function(p) { d <- mem[pid == p]; setdiff(unique(na.omit(c(d$source_gene, d$mapped_symbol, d$canonical_gene))), c("", "NA")) })

pool <- setdiff(unique(unlist(univ)), grep("^ENSG", unique(unlist(univ)), value = TRUE))
all_class_genes <- unique(unlist(lapply(arm_units, function(u) unlist(u$genes))))
G_all <- unique(c(pool, anch, all_class_genes))
inU <- vapply(Ukeys, function(k) as.numeric(G_all %in% k), numeric(length(G_all)))
Hm  <- vapply(prog_keys, function(k) as.numeric(G_all %in% k), numeric(length(G_all)))
colnames(inU) <- cts; colnames(Hm) <- prog_meta$pid
Hct <- vapply(cts, function(ct) as.numeric(rowSums(Hm[, prog_meta$cell_type == ct, drop = FALSE]) > 0), numeric(length(G_all)))
chk("cell-type union matches member counts", all(colSums(Hct) <= Kct + 0 & colSums(Hct) > 0), paste(colSums(Hct), collapse = "/"), fatal = FALSE)
cat("incidence:", nrow(Hm), "genes x", ncol(Hm), "programs; pool", length(pool), "\n")

pvals_from <- function(X_p, X_ct, n_ct) {
  n_p <- n_ct[ct_of_prog]; N_p <- N[ct_of_prog]
  list(prog = phyper(X_p - 1, prog_meta$K, N_p - prog_meta$K, n_p, lower.tail = FALSE),
       ct = phyper(X_ct - 1, Kct, N - Kct, n_ct, lower.tail = FALSE))
}
or_of <- function(a, n, K, Nn) {
  if (n == 0) return(NA_real_)
  tryCatch(unname(fisher.test(matrix(c(a, n - a, K - a, Nn - K - (n - a)), 2))$estimate), error = function(e) NA_real_)
}
# phyper must equal the one-sided Fisher exact test
set.seed(SEED)
tst <- list(c(3, 60, 800, 12969), c(0, 10, 100, 9000), c(2, 15, 200, 10000))
chk("phyper equals one-sided fisher.test", all(vapply(tst, function(z) abs(phyper(z[1] - 1, z[3], z[4] - z[3], z[2], lower.tail = FALSE) -
  fisher.test(matrix(c(z[1], z[2] - z[1], z[3] - z[1], z[4] - z[3] - (z[2] - z[1])), 2), alternative = "greater")$p.value) < 1e-10, logical(1))))

# ---- 4. tests per arm -----------------------------------------------------------------------
design_power <- function(nD, nC) {
  rbindlist(lapply(seq_len(nD + nC), function(Kh) {
    a <- max(0, Kh - nC):min(nD, Kh)
    p <- vapply(a, function(x) fisher.test(matrix(c(x, nD - x, Kh - x, nC - (Kh - x)), 2))$p.value, 1)
    data.table(hit_loci_total = Kh, n_discordant = nD, n_concordant = nC, smallest_exact_p = min(p))
  }))
}
run_arm <- function(arm, units) {
  set.seed(SEED)
  nU <- nrow(units); lab <- units$class
  levs <- c("discordant", "concordant")
  nby <- vapply(levs, function(l) sum(lab == l), 1L)
  tested <- levs[nby >= 3L]
  Gu <- unique(unlist(units$genes)); gi <- match(Gu, G_all)
  Hm_u <- Hm[gi, , drop = FALSE]; Hct_u <- Hct[gi, , drop = FALSE]; inU_u <- inU[gi, , drop = FALSE]
  Mloc <- do.call(rbind, lapply(units$genes, function(g) as.numeric(Gu %in% g)))
  stat <- function(idx) {
    sub <- Mloc[idx, , drop = FALSE]; v <- as.numeric(colSums(sub) > 0)
    list(X_p = as.numeric(v %*% Hm_u), X_ct = as.numeric(v %*% Hct_u), n_ct = as.numeric(v %*% inU_u),
         hit_p = colSums((sub %*% Hm_u) > 0), hit_ct = colSums((sub %*% Hct_u) > 0), n_genes = sum(v))
  }
  both <- length(tested) == 2L
  obs <- lapply(setNames(levs, levs), function(l) stat(which(lab == l)))
  omni <- function(sp, sct) {
    c(prog = sum((sp$D$hit_p / nby[1] - sp$C$hit_p / nby[2])^2), ct = sum((sp$D$hit_ct / nby[1] - sp$C$hit_ct / nby[2])^2))
  }
  O_obs <- if (both) omni(list(D = obs$discordant, C = obs$concordant), NULL) else c(prog = NA_real_, ct = NA_real_)
  # permutations of class labels across the arm's classed loci
  pX_p <- pX_ct <- list(); pnp <- list()
  for (l in levs) { pX_p[[l]] <- matrix(0, B, 117); pX_ct[[l]] <- matrix(0, B, 5); pnp[[l]] <- numeric(B) }
  pO <- matrix(NA_real_, B, 2)
  for (b in seq_len(B)) {
    pl <- sample(lab)
    ps <- lapply(setNames(levs, levs), function(l) stat(which(pl == l)))
    for (l in tested) {
      pX_p[[l]][b, ] <- ps[[l]]$X_p; pX_ct[[l]][b, ] <- ps[[l]]$X_ct
      pnp[[l]][b] <- sum(pvals_from(ps[[l]]$X_p, ps[[l]]$X_ct, ps[[l]]$n_ct)$prog < 0.05)
    }
    if (both) pO[b, ] <- omni(list(D = ps$discordant, C = ps$concordant), NULL)
  }
  tests <- rbindlist(lapply(tested, function(l) {
    o <- obs[[l]]; pv <- pvals_from(o$X_p, o$X_ct, o$n_ct)
    pp <- (1 + colSums(sweep(pX_p[[l]], 2, o$X_p, ">="))) / (B + 1)
    pc <- (1 + colSums(sweep(pX_ct[[l]], 2, o$X_ct, ">="))) / (B + 1)
    nu <- nby[[l]]
    rbind(
      data.table(level = "program", id = prog_meta$pid, cell_type = prog_meta$cell_type, module_name = prog_meta$module_name,
                 K = prog_meta$K, N = N[ct_of_prog], n_class_genes_in_universe = o$n_ct[ct_of_prog], X = o$X_p, hit_loci = o$hit_p,
                 fisher_p = pv$prog, perm_p = pp),
      data.table(level = "celltype", id = cts, cell_type = cts, module_name = "union of programs",
                 K = Kct, N = N, n_class_genes_in_universe = o$n_ct, X = o$X_ct, hit_loci = o$hit_ct,
                 fisher_p = pv$ct, perm_p = pc))[, `:=`(arm = arm, class = l, n_loci_class = nu, n_genes_class = o$n_genes)]
  }), fill = TRUE)
  if (nrow(tests)) {
    tests[, fisher_or := mapply(or_of, X, n_class_genes_in_universe, K, N)]
    tests[, fisher_bh_q := p.adjust(fisher_p, "BH"), by = .(arm, class, level)]
    tests[, supported := fisher_bh_q < 0.05 & perm_p <= 0.05]
  }
  contrast <- data.table(); omnibus <- data.table(); pw <- data.table()
  if (both) {
    D <- obs$discordant; C <- obs$concordant; nD <- nby[[1]]; nC <- nby[[2]]
    cf <- function(hD, hC) vapply(seq_along(hD), function(i) fisher.test(matrix(c(hD[i], nD - hD[i], hC[i], nC - hC[i]), 2))$p.value, 1)
    contrast <- rbind(
      data.table(level = "program", id = prog_meta$pid, cell_type = prog_meta$cell_type, hit_discordant = D$hit_p, n_discordant = nD,
                 hit_concordant = C$hit_p, n_concordant = nC, exact_p = cf(D$hit_p, C$hit_p)),
      data.table(level = "celltype", id = cts, cell_type = cts, hit_discordant = D$hit_ct, n_discordant = nD,
                 hit_concordant = C$hit_ct, n_concordant = nC, exact_p = cf(D$hit_ct, C$hit_ct)))
    contrast[, `:=`(arm = arm, bh_q = p.adjust(exact_p, "BH")), by = level]
    omnibus <- data.table(arm = arm, level = c("program", "celltype"), statistic = O_obs,
                          perm_p = c((1 + sum(pO[, 1] >= O_obs[1])) / (B + 1), (1 + sum(pO[, 2] >= O_obs[2])) / (B + 1)))
    pw <- design_power(nD, nC)[, arm := arm]
  }
  nc1 <- rbindlist(lapply(tested, function(l) {
    o <- obs[[l]]; ob <- sum(pvals_from(o$X_p, o$X_ct, o$n_ct)$prog < 0.05)
    data.table(arm = arm, class = l, observed_programs_p_lt_0.05 = ob, permuted_mean = mean(pnp[[l]]), permuted_p95 = unname(quantile(pnp[[l]], 0.95)),
               permuted_max = max(pnp[[l]]), empirical_p_count = (1 + sum(pnp[[l]] >= ob)) / (B + 1))
  }), fill = TRUE)
  # descriptive: which cell types and programs each classed locus falls in
  unit_hits <- rbindlist(lapply(seq_len(nU), function(i) {
    g <- units$genes[[i]]; gi2 <- match(g, G_all)
    rbindlist(lapply(seq_along(cts), function(k) {
      hit <- which(prog_meta$cell_type == cts[k] & colSums(Hm[gi2, , drop = FALSE]) > 0)
      data.table(arm = arm, unit_id = units$unit_id[i], class = lab[i], genes = paste(g, collapse = ","), cell_type = cts[k],
                 genes_in_universe = sum(inU[gi2, k]), n_programs = length(hit),
                 programs = paste0(prog_meta$module[hit], ":", prog_meta$module_name[hit], collapse = "; "),
                 genes_in_programs = paste(G_all[gi2][rowSums(Hm[gi2, prog_meta$cell_type == cts[k], drop = FALSE]) > 0], collapse = ","))
    }))
  }))
  list(tests = tests, contrast = contrast, omnibus = omnibus, power = pw, nc1 = nc1, unit_hits = unit_hits, tested = tested)
}

res <- lapply(names(arm_units), function(a) { cat("arm", a, "\n"); run_arm(a, arm_units[[a]]) })
names(res) <- names(arm_units)
cmb <- function(nm) rbindlist(lapply(res, `[[`, nm), fill = TRUE)
tests <- cmb("tests"); contrast <- cmb("contrast"); omnibus <- cmb("omnibus"); pwr <- cmb("power"); nc1 <- cmb("nc1"); unit_hits <- cmb("unit_hits")
if (nrow(omnibus)) omnibus[, holm_p := p.adjust(perm_p, "holm")]
fwrite(tests, file.path(OUT, "program_tests.tsv"), sep = "\t")
fwrite(contrast, file.path(OUT, "contrast_tests.tsv"), sep = "\t")
fwrite(omnibus, file.path(OUT, "omnibus.tsv"), sep = "\t")
fwrite(pwr, file.path(OUT, "design_power.tsv"), sep = "\t")
fwrite(nc1, file.path(OUT, "nc1_label_permutation.tsv"), sep = "\t")
fwrite(unit_hits, file.path(OUT, "locus_cell_type_programs.tsv"), sep = "\t")

# ---- 5. controls ------------------------------------------------------------------------------
set.seed(SEED)
ctl <- list()
# PC1: pooled anchored genes against cell-type unions
gi <- match(anch, G_all)
n_ct <- as.numeric(colSums(inU[gi, , drop = FALSE])); X_ct <- as.numeric(colSums(Hct[gi, , drop = FALSE]))
pc1 <- data.table(cell_type = cts, K = Kct, N = N, n_anchored_in_universe = n_ct, X = X_ct,
                  fisher_p = pvals_from(rep(0, 117), X_ct, n_ct)$ct, odds_ratio = mapply(or_of, X_ct, n_ct, Kct, N))
pc1_pass <- pc1[cell_type == "hepatocytes", fisher_p < 0.05] && pc1[cell_type == "hepatocytes", odds_ratio] == max(pc1$odds_ratio, na.rm = TRUE)
fwrite(pc1, file.path(OUT, "pc1_anchored_by_cell_type.tsv"), sep = "\t")
ctl$PC1 <- list(pass = pc1_pass, detail = "436 anchored genes: hepatocyte union P < 0.05 and largest odds ratio of five cell types")
chk("PC1 anchored genes land in hepatocyte programs", pc1_pass, paste(round(pc1$odds_ratio, 2), collapse = "/"), fatal = FALSE)

# PC2: lipid class (Hallmark) against lipid hepatocyte programs
gm <- strsplit(readLines(GMT), "\t")
hall <- unique(unlist(lapply(gm[vapply(gm, function(z) z[1] %in% c("HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_BILE_ACID_METABOLISM"), NA)], `[`, -(1:2))))
chk("Hallmark lipid sets found", length(hall) > 100, length(hall))
lipid_loci <- h10[mhc == FALSE & vapply(gene_list, function(g) any(g %in% hall), logical(1))]
lip_genes <- unique(unlist(lipid_loci$gene_list))
lip_prog <- paste0("hepatocytes_", c(5L, 10L, 16L, 19L, 28L))   # fixed in PRESPEC_H9.md section 7
lip_regex <- prog_meta[cell_type == "hepatocytes" & grepl("lipid|lipoprotein|fatty", module_name, ignore.case = TRUE), pid]
chk("regex lipid|lipoprotein|fatty reproduces modules 5, 10, 16, 19, 28", setequal(lip_prog, lip_regex), paste(lip_regex, collapse = ","), fatal = FALSE)
chk("lipid programs exist in the registry", all(lip_prog %in% prog_meta$pid))
if (nrow(lipid_loci) >= 3L) {
  gi <- match(lip_genes, G_all); n1 <- sum(inU[gi, "hepatocytes"])
  hp <- which(prog_meta$cell_type == "hepatocytes")
  Xh <- as.numeric(colSums(Hm[gi, hp, drop = FALSE]))
  ph <- phyper(Xh - 1, prog_meta$K[hp], N["hepatocytes"] - prog_meta$K[hp], n1, lower.tail = FALSE)
  orh <- mapply(or_of, Xh, n1, prog_meta$K[hp], N["hepatocytes"])
  pc2 <- data.table(pid = prog_meta$pid[hp], module_name = prog_meta$module_name[hp], lipid_program = prog_meta$pid[hp] %in% lip_prog, K = prog_meta$K[hp], X = Xh, fisher_p = ph, odds_ratio = orh)
  fwrite(pc2, file.path(OUT, "pc2_lipid_class_hepatocyte_programs.tsv"), sep = "\t")
  pc2_pass <- any(pc2[lipid_program == TRUE, fisher_p] < 0.05) && median(pc2[lipid_program == TRUE, odds_ratio], na.rm = TRUE) > median(pc2[lipid_program == FALSE, odds_ratio], na.rm = TRUE)
  ctl$PC2 <- list(pass = pc2_pass, n_lipid_loci = nrow(lipid_loci), n_lipid_genes = length(lip_genes))
  chk("PC2 lipid class aligns with lipid hepatocyte programs", pc2_pass, paste0(nrow(lipid_loci), " lipid loci, ", length(lip_genes), " genes"), fatal = FALSE)
} else ctl$PC2 <- list(pass = NA, n_lipid_loci = nrow(lipid_loci))

# NC2: random gene sets of the arm B class sizes
nB <- arm_summary[arm == "B"]
sizes <- c(discordant = nB$n_genes_discordant, concordant = nB$n_genes_concordant)
pool_i <- match(pool, G_all)
nc2 <- rbindlist(lapply(names(sizes), function(l) {
  r <- vapply(seq_len(B), function(b) {
    i <- sample(pool_i, sizes[[l]]); n_ct <- as.numeric(colSums(inU[i, , drop = FALSE]))
    p <- pvals_from(as.numeric(colSums(Hm[i, , drop = FALSE])), as.numeric(colSums(Hct[i, , drop = FALSE])), n_ct)$prog
    c(mean(p < 0.05), min(p.adjust(p, "BH")) < 0.05)
  }, numeric(2))
  data.table(size_of = l, n_genes = sizes[[l]], draws = B, mean_fraction_programs_p_lt_0.05 = mean(r[1, ]), draws_with_any_bh_q_lt_0.05 = sum(r[2, ]))
}))
fwrite(nc2, file.path(OUT, "nc2_random_gene_sets.tsv"), sep = "\t")
chk("NC2 random sets: fraction of program tests with P < 0.05 at most 0.055 (about 5%)", all(nc2$mean_fraction_programs_p_lt_0.05 <= 0.055), paste(round(nc2$mean_fraction_programs_p_lt_0.05, 4), collapse = "/"), fatal = FALSE)

# ---- 6. summary -----------------------------------------------------------------------------------
tests_nom <- tests[, .(cells = .N, nominal_p_lt_0.05 = sum(fisher_p < 0.05), min_fisher_p = min(fisher_p), min_bh_q = min(fisher_bh_q),
                       supported = sum(supported), min_perm_p = min(perm_p)), by = .(arm, class, level)]
fwrite(tests_nom, file.path(OUT, "test_summary.tsv"), sep = "\t")
summ <- list(out = OUT, seed = SEED, permutations = B, arm_summary = arm_summary, controls = ctl,
             omnibus = omnibus, nc1 = nc1, nc2 = nc2, test_summary = tests_nom,
             contrast_min_bh_q = if (nrow(contrast)) contrast[, .(min_exact_p = min(exact_p), min_bh_q = min(bh_q)), by = .(arm, level)] else NULL,
             design_power_min = if (nrow(pwr)) pwr[, .(smallest_p_overall = min(smallest_exact_p)), by = arm] else NULL,
             checks_failed = checks[pass == FALSE, check])
write_json(summ, file.path(OUT, "summary.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE, null = "null", na = "null")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("H9 done", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
