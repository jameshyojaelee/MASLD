# RV7 post hoc (written after reading conditional_C1/C2): is there a local enzyme association at GNMT once the
# distant window lead (6:43.74-43.81 Mb) is removed? Not prespecified. Usage: Rscript rv7_posthoc_local_enzyme.R <rv7_out_dir>
suppressPackageStartupMessages({library(data.table); library(coloc)})
OUT <- commandArgs(trailingOnly = TRUE)[1]
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
reg <- fread(file.path(FM, "config/gwas_registry.tsv"))
ANCHOR <- "6:42931261"; APOS <- 42931261L; EQTL_N <- 1183
S <- c("2021_34841290_NAFLD_EUR", "MVP_ALT_EUR", "MVP_AST_EUR", "UKBB_ALT", "UKBB_AST", "UKBB_GGT")
dat <- setNames(lapply(S, function(s) fread(file.path(OUT, sprintf("variants_%s.tsv.gz", s)))), S)
leads <- fread(file.path(OUT, "leads.tsv"))

raw <- fread(file.path(OUT, "ld_window.raw")); G <- as.matrix(raw[, -(1:6)])
bim <- fread(file.path(FM, "data/ld_ref/1kg_eur/chr6_eur.bim"), header = FALSE, col.names = c("chr", "id", "cm", "pos", "b1", "b2"))
cid <- sub("_[ACGT]+$", "", colnames(G)); counted <- toupper(sub("^.*_", "", colnames(G)))
pinfo <- bim[match(cid, id)][, `:=`(counted = counted, col = .I)]
pinfo[, other := toupper(fifelse(counted == toupper(b1), b2, b1))][, merge_key := paste(chr, pos, sep = ":")]
ea_tab <- unique(rbindlist(lapply(dat, function(m) m[, .(merge_key, ea, nea)])))
pinfo <- merge(pinfo, ea_tab, by = "merge_key")
pinfo[, sgn := fifelse(counted == ea & other == nea, 1, fifelse(counted == nea & other == ea, -1, NA_real_))]
pinfo <- pinfo[!is.na(sgn)][!duplicated(merge_key)]
X <- sweep(G[, pinfo$col, drop = FALSE], 2, pinfo$sgn, `*`)
X <- apply(X, 2, function(v) { v[is.na(v)] <- mean(v, na.rm = TRUE); v })
keep <- apply(X, 2, sd) > 0; X <- X[, keep]; colnames(X) <- pinfo$merge_key[keep]

abf <- function(m, s, beta, varbeta) {
  r <- reg[study_name == s]
  d1 <- list(beta = beta, varbeta = varbeta, N = r$N_tot, snp = m$merge_key)
  if (r$trait_type == "binary") { d1$type <- "cc"; d1$s <- r$N_cases / r$N_tot } else { d1$type <- "quant"; d1$sdY <- 1 }
  d2 <- list(beta = m$eqtl_beta, varbeta = m$eqtl_se^2, N = EQTL_N, type = "quant", sdY = 1, snp = m$merge_key)
  a <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6))
  c(PP.H3 = unname(a$summary["PP.H3.abf"]), PP.H4 = unname(a$summary["PP.H4.abf"]),
    top = a$results$snp[which.max(a$results$SNP.PP.H4)], top_snp_pp4 = max(a$results$SNP.PP.H4), n = nrow(m))
}
res <- rbindlist(lapply(S, function(s) {
  m <- dat[[s]]
  full <- abf(m, s, m$gwas_beta, m$se^2)
  loc <- m[abs(pos - APOS) <= 250000]
  local <- abf(loc, s, loc$gwas_beta, loc$se^2)
  out <- data.table(study = s, z_anchor = m[merge_key == ANCHOR, z], p_anchor = m[merge_key == ANCHOR, pval],
                    local_min_p = min(loc$pval), local_min_variant = loc[which.min(pval), merge_key],
                    full_PP.H3 = full[["PP.H3"]], full_PP.H4 = full[["PP.H4"]], full_top = full[["top"]],
                    local250_PP.H3 = local[["PP.H3"]], local250_PP.H4 = local[["PP.H4"]], local250_top = local[["top"]],
                    local250_top_snp_pp4 = local[["top_snp_pp4"]], local250_n = local[["n"]])
  if (s %in% leads$study) {
    L <- leads[study == s, lead]
    mc <- m[merge_key %in% colnames(X)]; zL <- mc[merge_key == L, z]
    r <- drop(cor(X[, mc$merge_key], X[, L])); mc[, r := r]; mc <- mc[r^2 <= 0.95]
    mc[, `:=`(z_c = (z - r * zL) / sqrt(1 - r^2), se_c = se / sqrt(1 - r^2))]
    cc <- abf(mc, s, mc$z_c * mc$se_c, mc$se_c^2)
    out[, `:=`(cond_on_lead = L, cond_PP.H3 = cc[["PP.H3"]], cond_PP.H4 = cc[["PP.H4"]], cond_top = cc[["top"]],
               cond_z_anchor = mc[merge_key == ANCHOR, z_c])]
  }
  out
}), fill = TRUE)
fwrite(res, file.path(OUT, "posthoc_local_enzyme.tsv"), sep = "\t")
sink(file.path(OUT, "posthoc_summary.txt")); print(res); sink()
cat("done\n")
