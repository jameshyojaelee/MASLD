# RV7: GNMT direct-NAFLD association conditioned on liver-enzyme leads. Follows PRESPEC_RV7.md.
# Usage: Rscript rv7_gnmt_conditional.R <out_dir>   (expects <out_dir>/ld_window.raw from plink --recode A)
suppressPackageStartupMessages({library(data.table); library(coloc)})
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
EQ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/broadaway_eqtl/chr6_marginal_summary_results.tsv"
BIM <- file.path(FM, "data/ld_ref/1kg_eur/chr6_eur.bim")
GENE <- "ENSG00000124713"; FROM <- 41928548L; TO <- 43928108L
ANCHOR <- "6:42931261"; EQTL_N <- 1183
NAFLD <- "2021_34841290_NAFLD_EUR"
ENZ <- c("MVP_ALT_EUR", "MVP_AST_EUR", "UKBB_ALT", "UKBB_AST", "UKBB_GGT")
reg <- fread(file.path(FM, "config/gwas_registry.tsv"))
wr <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t", na = "NA")

eq <- fread(cmd = sprintf("awk -F'\\t' 'NR==1 || $14==\"%s\"' %s", GENE, EQ))
stopifnot(nrow(eq) > 4000, all(eq$CHR == 6))
eq[, merge_key := paste(CHR, POS, sep = ":")]

load_study <- function(s) {
  path <- file.path(FM, reg[study_name == s, sumstats_path])
  g <- fread(cmd = sprintf("awk -F'\\t' 'NR==1 || ($1==\"6\" && $2>=%d && $2<=%d)' %s", FROM, TO, path))
  g[, merge_key := paste(chromosome, position, sep = ":")]
  m <- merge(eq[, .(merge_key, pos = POS, ea = toupper(EA), nea = toupper(NEA), eqtl_beta = Beta, eqtl_se = SE, eqtl_p = PVAL)],
             g[, .(merge_key, a1 = toupper(allele1), a2 = toupper(allele2), beta, se, pval)], by = "merge_key")
  m <- m[order(pval)][!duplicated(merge_key)]
  m[, h := fifelse(ea == a1 & nea == a2, 1, fifelse(ea == a2 & nea == a1, -1, NA_real_))]
  m <- m[!is.na(h)]
  m[, `:=`(gwas_beta = h * beta, z = h * beta / se)]  # GWAS effect aligned to the eQTL EA
  m[, study := s][]
}
S <- c(NAFLD, ENZ)
dat <- setNames(lapply(S, load_study), S)
for (s in S) wr(dat[[s]], sprintf("variants_%s.tsv.gz", s))

run_abf <- function(m, s, beta = m$gwas_beta, varbeta = m$se^2) {
  r <- reg[study_name == s]
  d1 <- list(beta = beta, varbeta = varbeta, N = r$N_tot, snp = m$merge_key)
  if (r$trait_type == "binary") { d1$type <- "cc"; d1$s <- r$N_cases / r$N_tot } else { d1$type <- "quant"; d1$sdY <- 1 }
  d2 <- list(beta = m$eqtl_beta, varbeta = m$eqtl_se^2, N = EQTL_N, type = "quant", sdY = 1, snp = m$merge_key)
  a <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6))
  list(pp4 = unname(a$summary["PP.H4.abf"]), pp3 = unname(a$summary["PP.H3.abf"]),
       top = a$results$snp[which.max(a$results$SNP.PP.H4)], n = nrow(m))
}
repro <- run_abf(dat[[NAFLD]], NAFLD)
cat(sprintf("Reproduction: NAFLD n %d (released 3384), PP.H4 %.4f (released 0.7058), top %s\n", repro$n, repro$pp4, repro$top))

# ---- LD from the EUR panel, dosage aligned to the eQTL effect allele ----
raw <- fread(file.path(OUT, "ld_window.raw"))
G <- as.matrix(raw[, -(1:6)])
bim <- fread(BIM, header = FALSE, col.names = c("chr", "id", "cm", "pos", "b1", "b2"))
cid <- sub("_[ACGT]+$", "", colnames(G)); counted <- sub("^.*_", "", colnames(G))
pinfo <- bim[match(cid, id)]
pinfo[, `:=`(counted = toupper(counted), other = toupper(fifelse(toupper(counted) == toupper(b1), b2, b1)))]
pinfo[, merge_key := paste(chr, pos, sep = ":")]
ea_tab <- unique(rbindlist(lapply(dat, function(m) m[, .(merge_key, ea, nea)])))
pinfo <- merge(pinfo[, .(col = seq_len(.N), merge_key, counted, other)], ea_tab, by = "merge_key")
pinfo[, sgn := fifelse(counted == ea & other == nea, 1, fifelse(counted == nea & other == ea, -1, NA_real_))]
pinfo <- pinfo[!is.na(sgn)][!duplicated(merge_key)]
X <- sweep(G[, pinfo$col, drop = FALSE], 2, pinfo$sgn, `*`)
X <- apply(X, 2, function(v) { v[is.na(v)] <- mean(v, na.rm = TRUE); v })
keep <- apply(X, 2, sd) > 0
X <- X[, keep, drop = FALSE]; colnames(X) <- pinfo$merge_key[keep]
cat(sprintf("Panel variants aligned: %d of %d in the window\n", ncol(X), ncol(G)))
rcol <- function(k) { r <- drop(cor(X, X[, k])); names(r) <- colnames(X); r }
wr(data.table(merge_key = colnames(X), maf = pmin(colMeans(X) / 2, 1 - colMeans(X) / 2)), "ld_window.tsv.gz")

cond <- function(m, L) {
  m <- m[merge_key %in% colnames(X)]
  zL <- m[merge_key == L, z]; stopifnot(length(zL) == 1)
  r <- rcol(L)[m$merge_key]
  m <- copy(m)[, r := r][r^2 <= 0.95]
  m[, `:=`(z_cond = (z - r * zL) / sqrt(1 - r^2), se_cond = se / sqrt(1 - r^2))]
  m[, `:=`(beta_cond = z_cond * se_cond, p_cond = 2 * pnorm(-abs(z_cond)))][]
}
# conditioning variant must be in the panel and in the NAFLD set (z_L comes from the GWAS being conditioned)
lead_in_panel <- function(m) m[merge_key %in% intersect(colnames(X), dat[[NAFLD]]$merge_key)][which.min(pval), merge_key]

leads <- rbindlist(lapply(ENZ, function(s) {
  m <- dat[[s]]
  data.table(study = s, lead_any = m[which.min(pval), merge_key], lead_p_any = min(m$pval),
             lead = lead_in_panel(m), lead_p = m[merge_key == lead_in_panel(m), pval])
}))
wr(leads, "leads.tsv")
stopifnot(ANCHOR %in% colnames(X), all(sapply(dat, function(m) ANCHOR %in% m$merge_key)))
eqtl_lead <- dat[[NAFLD]][which.min(eqtl_p), merge_key]

# ---- C1: NAFLD conditioned on each distinct enzyme lead ----
nf <- dat[[NAFLD]]; nfp <- nf[merge_key %in% colnames(X)]
base_panel <- run_abf(nfp, NAFLD)
c1 <- rbindlist(lapply(unique(leads$lead), function(L) {
  m <- cond(nf, L)
  a <- run_abf(m, NAFLD, beta = m$beta_cond, varbeta = m$se_cond^2)
  za <- m[merge_key == ANCHOR, z_cond]
  data.table(enzyme_lead = L, enzyme_studies = paste(leads[lead == L, study], collapse = ","),
             r_anchor_lead = rcol(L)[[ANCHOR]], nafld_z_lead = nf[merge_key == L, z],
             anchor_z = nf[merge_key == ANCHOR, z], anchor_z_cond = if (length(za)) za else NA_real_,
             region_min_p = min(nfp$pval), region_min_p_cond = min(m$p_cond), region_min_cond_variant = m[which.min(p_cond), merge_key],
             n_unconditional_panel = base_panel$n, pp4_unconditional_panel = base_panel$pp4,
             n_cond = a$n, pp4_cond = a$pp4, pp3_cond = a$pp3, top_cond = a$top)
}))
c1[, call := fifelse(abs(anchor_z_cond) >= 3 & pp4_cond >= 0.5, "separate",
             fifelse(abs(anchor_z_cond) < 2 & pp4_cond < 0.5, "not separable", "unresolved"))]
wr(c1, "conditional_C1.tsv")
overall <- if (length(unique(c1$call)) == 1) c1$call[1] else "unresolved"

# ---- C2: each enzyme conditioned on the NAFLD anchor ----
c2 <- rbindlist(lapply(ENZ, function(s) {
  m <- dat[[s]]; mp <- m[merge_key %in% colnames(X)]; mc <- cond(m, ANCHOR)
  data.table(study = s, enzyme_z_anchor = m[merge_key == ANCHOR, z], region_min_p = min(mp$pval),
             region_min_p_cond = min(mc$p_cond), region_min_cond_variant = mc[which.min(p_cond), merge_key],
             lead_z = m[merge_key == leads[study == s, lead], z], lead_z_cond = mc[merge_key == leads[study == s, lead], z_cond])
}))
wr(c2, "conditional_C2.tsv")

# ---- C3: LD among anchor, eQTL lead, enzyme leads ----
pts <- unique(c(ANCHOR, eqtl_lead, leads$lead))
c3 <- CJ(a = pts, b = pts)[a < b]
c3[, r := mapply(function(a, b) cor(X[, a], X[, b]), a, b)][, r2 := r^2]
wr(c3, "ld_C3.tsv")

sink(file.path(OUT, "summary.txt"))
cat(sprintf("Reproduction: NAFLD variants %d (released 3384); PP.H4 %.4f (released 0.7058); top %s\n", repro$n, repro$pp4, repro$top))
cat(sprintf("EUR panel (379) variants aligned in window: %d; NAFLD variants in panel: %d\n", ncol(X), nrow(nfp)))
cat(sprintf("eQTL lead %s; NAFLD anchor %s\n\nEnzyme leads\n", eqtl_lead, ANCHOR)); print(leads)
cat("\nC1: NAFLD conditioned on each enzyme lead\n"); print(c1)
cat(sprintf("\nOverall call: %s\n", overall))
cat("\nC2: enzyme GWAS conditioned on the NAFLD anchor\n"); print(c2)
cat("\nC3: LD (EUR, 379)\n"); print(c3)
sink()
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("done\n")
