#!/usr/bin/env Rscript
# Independent re-derivation of Fig 2B point positions (claim R02-04).
# Replicates scripts/figures/fig2_resolution_plane.R lines 40-169 (read-only copy of logic),
# then recomputes x under alternative definitions. Writes only to this job dir.
suppressPackageStartupMessages(library(data.table))
set.seed(1)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT  <- Sys.getenv("OUT_DIR", getwd())
DRY  <- nzchar(Sys.getenv("DRYRUN"))
FM   <- file.path(BASE, "GWAS/finemapping/results")
COLOC_INPUTS <- c(
  r2      = file.path(FM, "susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"),
  default = file.path(FM, "susie_coloc/susie_coloc_all_gwas.csv"))
SRC <- fread(file.path(BASE, "figures/main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc_source.csv"))
MAIN <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[placement == "main", study_name]
cat("MAIN studies:", length(MAIN), "\n")

read_coloc <- function(path) {
  cl <- fread(path, select = c("gwas_name", "gene", "PP.H4.susie", "PP.H4.abf", "top_snp", "top_snp_PP"))
  cl <- cl[gwas_name %in% MAIN]
  setnames(cl, c("PP.H4.susie", "PP.H4.abf"), c("pp4_susie", "pp4_abf"))
  cl <- cl[!is.na(gene) & gene != ""]
  cl[, pp4_susie := suppressWarnings(as.numeric(pp4_susie))]
  cl[, pp4_abf := suppressWarnings(as.numeric(pp4_abf))]
  cl <- cl[!is.na(top_snp) & top_snp != ""]
  cl[, c("snp_chr", "snp_pos") := tstrsplit(top_snp, ":", fixed = TRUE)]
  cl[, snp_pos := as.integer(snp_pos)]
  cl[, row_id := .I]
  cl
}
colocs <- lapply(COLOC_INPUTS, read_coloc)
keys <- unique(rbindlist(lapply(colocs, function(x) x[, .(study = gwas_name, chromosome = snp_chr, position = snp_pos)])))
cat("distinct coloc lead keys:", nrow(keys), "\n")

t0 <- Sys.time()
canon_cols <- c("study", "chromosome", "position", "ancestry", "susie_pip", "carma_pip",
                "recommended_pip", "susiex_pip", "mesusie_pip")
canon <- fread(file.path(FM, "combined_finemapping.csv"), select = canon_cols,
               nrows = if (DRY) 3e6 else Inf, colClasses = list(character = "chromosome"))
cat("canonical rows read:", nrow(canon), " secs:", round(as.numeric(Sys.time() - t0, units = "secs")), "\n")
canon <- canon[keys, on = .(study, chromosome, position), nomatch = 0L]
cat("canonical rows at coloc keys:", nrow(canon), "\n")
canon[, x_source := "canonical_finemapping"]
targ <- fread(file.path(FM, "combined_finemapping_coloc_targeted.csv"),
              select = c("study", "chromosome", "position", "ancestry", "susie_pip", "carma_pip", "recommended_pip"),
              colClasses = list(character = "chromosome"))
targ <- targ[keys, on = .(study, chromosome, position), nomatch = 0L]
targ[, `:=`(susiex_pip = NA_real_, mesusie_pip = NA_real_, x_source = "coloc_targeted_finemapping")]
fm <- rbindlist(list(canon, targ), use.names = TRUE)
num <- c("susie_pip", "carma_pip", "recommended_pip", "susiex_pip", "mesusie_pip")
fm[, (num) := lapply(.SD, function(v) suppressWarnings(as.numeric(v))), .SDcols = num]
# Producer x definition (fig2_resolution_plane.R:75-76)
fm[, causal_pip := pmax(recommended_pip, susie_pip, na.rm = TRUE)]
fm[!is.finite(causal_pip), causal_pip := NA_real_]
fm[, source_rank := fifelse(x_source == "canonical_finemapping", 1L, 2L)]
# Targeted table: does recommended_pip exceed max(SuSiE, CARMA)?
cat("targeted rows with recommended_pip > max(susie,carma)+1e-9:",
    fm[x_source != "canonical_finemapping", sum(recommended_pip > pmax(susie_pip, carma_pip, na.rm = TRUE) + 1e-9, na.rm = TRUE)],
    "of", fm[x_source != "canonical_finemapping", .N], "\n")

max_or_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
label_at_max <- function(value, label) if (all(is.na(value))) NA_character_ else label[which.max(fifelse(is.na(value), -Inf, value))]
roster <- c("HSD17B13", "TM6SF2", "MBOAT7", "GCKR", "RORA")
gate <- function(pp4, x, gene) (pp4 > 0.05) | (!is.na(x) & x >= 0.99) | gene %in% roster

run_one <- function(nm, coloc) {
  cat("\n================ COLOC input:", nm, "================\n")
  cat("coloc MAIN rows:", nrow(coloc), " genes:", uniqueN(coloc$gene), "\n")
  # --- producer path
  f1 <- fm[is.finite(causal_pip)][order(study, chromosome, position, -causal_pip, source_rank)]
  fu <- unique(f1, by = c("study", "chromosome", "position"))
  cm <- merge(coloc, fu, by.x = c("gwas_name", "snp_chr", "snp_pos"), by.y = c("study", "chromosome", "position"))
  cat("COLOC top-SNP fine-mapping coverage:", nrow(cm), "/", nrow(coloc), "\n")
  gy <- coloc[, .(pp4_susie = max_or_na(pp4_susie),
                  susie_gwas = label_at_max(pp4_susie, gwas_name),
                  y_top_snp  = label_at_max(pp4_susie, top_snp)), by = gene]
  gx <- cm[order(gene, -causal_pip, x_source), .SD[1], by = gene]
  gx <- gx[, .(gene, causal_pip, x_source, x_top_snp = top_snp, x_gwas = gwas_name, x_anc = ancestry,
               x_susie = susie_pip, x_carma = carma_pip, x_rec = recommended_pip,
               x_susiex = susiex_pip, x_mesusie = mesusie_pip)]
  dt <- merge(gy, gx, by = "gene", all.x = TRUE)
  pl <- dt[!is.na(pp4_susie) & is.finite(causal_pip)]
  pl[, pp4 := pmax(pp4_susie, 0)]
  pl[, plotted := gate(pp4, causal_pip, gene)]
  cat("rows (genes with finite x and SuSiE y):", nrow(pl), "  plotted:", sum(pl$plotted), "\n")
  # --- agreement with the installed source CSV
  chk <- merge(pl[, .(gene, cx = round(causal_pip, 4), cy = round(pp4, 4), xs = x_source, pt = plotted, xt = x_top_snp, yg = susie_gwas)],
               SRC[, .(gene, sx = causal_pip, sy = pp4, sxs = x_source, spt = plotted, sxt = x_top_snp, syg = pp4_gwas)],
               by = "gene", all = TRUE)
  cat("source CSV match: genes in both", chk[!is.na(cx) & !is.na(sx), .N],
      "| only rederived", chk[is.na(sx), .N], "| only source", chk[is.na(cx), .N],
      "| x equal", chk[abs(cx - sx) < 1e-9, .N], "| y equal", chk[abs(cy - sy) < 1e-9, .N],
      "| x_source equal", chk[xs == sxs, .N], "| plotted equal", chk[pt == spt, .N],
      "| x_top_snp equal", chk[xt == sxt, .N], "| y gwas equal", chk[yg == syg, .N], "\n")
  P <- pl[plotted == TRUE]
  tol <- 1e-9
  P[, base_sc := pmax(x_susie, x_carma, na.rm = TRUE)]
  P[, joint_max := pmax(x_susiex, x_mesusie, na.rm = TRUE)]
  P[, origin := fcase(
        !is.na(joint_max) & abs(causal_pip - joint_max) < tol & (is.na(base_sc) | causal_pip > base_sc + tol), "joint_above_susie_carma",
        x_anc != "EUR" & !is.na(x_susie) & abs(causal_pip - x_susie) < tol & (is.na(x_carma) | x_susie > x_carma + tol), "susie_nonEUR_above_carma",
        !is.na(x_susie) & abs(causal_pip - x_susie) < tol, "susie",
        !is.na(x_carma) & abs(causal_pip - x_carma) < tol, "carma",
        default = "other")]
  cat("\nPlotted x origin (n =", nrow(P), "):\n"); print(P[, .N, by = origin][order(-N)])
  cat("joint origin by source:\n"); print(P[origin == "joint_above_susie_carma", .N, by = .(x_source, susiex_wins = !is.na(x_susiex) & abs(causal_pip - x_susiex) < tol)])
  cat("x ancestry of plotted points:\n"); print(P[, .N, by = x_anc][order(-N)])
  cat("\nplotted: x GWAS != y GWAS:", P[x_gwas != susie_gwas, .N],
      "  x lead SNP != y lead SNP:", P[x_top_snp != y_top_snp, .N],
      "  either differs:", P[x_gwas != susie_gwas | x_top_snp != y_top_snp, .N], "\n")
  corner <- function(x, y) sum(x >= 0.9 & y < 0.5, na.rm = TRUE)
  cat("corner (x>=0.9 & PP.H4<0.5) producer:", corner(P$causal_pip, P$pp4), "\n")

  # --- Alternative A: SuSiE-only x, same gene-level max over lead SNPs/GWAS
  fA <- fm[is.finite(susie_pip)][order(study, chromosome, position, -susie_pip, source_rank)]
  fuA <- unique(fA, by = c("study", "chromosome", "position"))
  cmA <- merge(coloc, fuA[, .(study, chromosome, position, susie_pip)], by.x = c("gwas_name", "snp_chr", "snp_pos"), by.y = c("study", "chromosome", "position"))
  gxA <- cmA[, .(xA = max(susie_pip)), by = gene]
  A <- merge(pl[, .(gene, pp4, causal_pip, plotted)], gxA, by = "gene", all.x = TRUE)
  A[, plotA := !is.na(xA) & gate(pp4, xA, gene)]
  PA <- A[plotted == TRUE]
  cat("\nAlt A (SuSiE-only x, gene-level max):\n")
  cat("  of producer-plotted points: x changes (|dx|>5e-5):", PA[!is.na(xA) & abs(xA - causal_pip) > 5e-5, .N],
      " x missing (no SuSiE PIP):", PA[is.na(xA), .N],
      " changed or missing:", PA[is.na(xA) | abs(xA - causal_pip) > 5e-5, .N], "\n")
  cat("  of producer-plotted: drop out of plot under A:", PA[plotA == FALSE, .N],
      " of which had x>=0.99 & pp4<=0.05 & not roster under producer:", PA[plotA == FALSE & causal_pip >= 0.99 & pp4 <= 0.05 & !gene %in% roster, .N], "\n")
  cat("  cross 0.99 (producer x>=0.99 and A x<0.99 or NA):", PA[causal_pip >= 0.99 & (is.na(xA) | xA < 0.99), .N],
      " reverse:", PA[causal_pip < 0.99 & !is.na(xA) & xA >= 0.99, .N], "\n")
  cat("  plotted under A:", A[plotA == TRUE, .N], " corner:", A[plotA == TRUE, corner(xA, pp4)], "\n")

  # --- Alternative B: SuSiE PIP at the SAME GWAS and lead SNP as y
  gyB <- gy[!is.na(pp4_susie)]
  gyB[, c("ychr", "ypos") := tstrsplit(y_top_snp, ":", fixed = TRUE)]
  gyB[, ypos := as.integer(ypos)]
  lkA <- fuA[, .(study, chromosome, position, xB = susie_pip)]
  B <- merge(gyB, lkA, by.x = c("susie_gwas", "ychr", "ypos"), by.y = c("study", "chromosome", "position"), all.x = TRUE)
  # Alternative C: producer causal_pip at the same GWAS and lead SNP as y
  lkC <- fu[, .(study, chromosome, position, xC = causal_pip)]
  B <- merge(B, lkC, by.x = c("susie_gwas", "ychr", "ypos"), by.y = c("study", "chromosome", "position"), all.x = TRUE)
  B[, pp4 := pmax(pp4_susie, 0)]
  B[, plotB := !is.na(xB) & gate(pp4, xB, gene)]
  B[, plotC := !is.na(xC) & gate(pp4, xC, gene)]
  cat("\nAlt B (SuSiE x at y's GWAS + lead SNP): genes with finite x:", B[!is.na(xB), .N],
      " plotted:", B[plotB == TRUE, .N], " corner:", B[plotB == TRUE, corner(xB, pp4)], "\n")
  cat("Alt C (producer causal_pip at y's GWAS + lead SNP): genes with finite x:", B[!is.na(xC), .N],
      " plotted:", B[plotC == TRUE, .N], " corner:", B[plotC == TRUE, corner(xC, pp4)], "\n")
  cat("\nRoster genes:\n")
  print(merge(merge(pl[gene %in% roster, .(gene, pp4 = round(pp4, 4), y_gwas = susie_gwas, y_top_snp, x_prod = round(causal_pip, 4),
                                          x_gwas, x_top_snp, x_anc, x_susie = round(x_susie, 4), x_carma = round(x_carma, 4),
                                          x_susiex = round(x_susiex, 4), x_mesusie = round(x_mesusie, 4), x_source)],
                    A[, .(gene, xA = round(xA, 4))], by = "gene", all.x = TRUE),
              B[, .(gene, xB = round(xB, 4), xC = round(xC, 4))], by = "gene", all.x = TRUE))
  fwrite(P, file.path(OUT, sprintf("plotted_points_%s.tsv", nm)), sep = "\t")
  fwrite(A, file.path(OUT, sprintf("altA_susie_only_%s.tsv", nm)), sep = "\t")
  fwrite(B, file.path(OUT, sprintf("altBC_same_row_%s.tsv", nm)), sep = "\t")
  invisible(NULL)
}
for (nm in names(colocs)) run_one(nm, colocs[[nm]])
cat("\nsessionInfo:\n"); print(sessionInfo())
