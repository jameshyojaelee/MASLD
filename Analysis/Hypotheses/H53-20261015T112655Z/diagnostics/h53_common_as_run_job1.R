# h53_common.R -- shared helpers for H53 (prespec inputs/PRESPEC_H53.md).
# Sourced by h53_job1.R and h53_job2.R. No outcome is computed here.
suppressPackageStartupMessages({
  library(data.table); library(coloc); library(parallel)
})

PROJ   <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN    <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H53-20261015T112655Z"
FM     <- file.path(PROJ, "GWAS/finemapping")
R2DIR  <- file.path(FM, "results/susie_coloc_r2_20260924T053244Z")
CQDIR  <- file.path(PROJ, "data/external/currin_2025_caqtl")
SSDIR  <- file.path(FM, "data/sumstats")
EQTL_SUSIE_DIR <- file.path(FM, "results/eqtl_susie_polyfun")
PRESPEC_SHA <- "f575ee90f2f895b2c00740c780d9ab78261f55737a21af4d9c57329c5b302a26"

# Fixed analysis constants (prespec "Modality and method", "Unit, sets and inclusion rules")
WIN        <- 5e5      # peak-centre window around hit1, locus merge gap
MAF_MIN    <- 0.01
W_PRIOR    <- 0.15     # caQTL prior SD (sdY = 1)
P1 <- 1e-4; P2 <- 1e-4; P12 <- 5e-6; P12_ALT <- 1e-5
MIN_SNPS   <- 100L
OVERLAP_MIN <- 0.5
PP_ENTRY   <- 0.8; PP_SECOND <- 0.5; LD_ENTRY <- 0.5
SEED       <- 20260923L
N_CAQTL    <- 138
PILOT_CHR  <- c("1", "19")

stamp <- function(...) cat(format(Sys.time(), "[%Y-%m-%dT%H:%M:%S]"), ..., "\n")

sha256 <- function(paths, cores = 8) {
  unlist(mclapply(paths, function(p) {
    out <- system2("sha256sum", shQuote(p), stdout = TRUE)
    sub(" .*", "", out[1])
  }, mc.cores = cores))
}

# Refuses to continue on any mismatch (prespec hash + every prespec-listed input).
verify_inputs <- function(cores = 8) {
  ph <- sha256(file.path(RUN, "inputs/PRESPEC_H53.md"), 1)
  if (ph != PRESPEC_SHA) stop("PRESPEC hash mismatch: ", ph)
  stamp("prespec sha256 OK", ph)
  ex <- fread(file.path(RUN, "inputs/expected_sha256.tsv"))
  got <- sha256(file.path(PROJ, ex$path), cores)
  ex[, observed := got]
  bad <- ex[sha256 != observed]
  fwrite(ex, file.path(RUN, "diagnostics", paste0("input_sha256_check_", Sys.getenv("SLURM_JOB_ID", "nojob"), ".tsv")), sep = "\t")
  if (nrow(bad)) { print(bad); stop("INPUT HASH MISMATCH for ", nrow(bad), " file(s); refusing to run") }
  stamp("input sha256 OK for", nrow(ex), "prespec-listed files")
  invisible(ex)
}

split_hit <- function(h) {
  data.table(chr = sub(":.*", "", h), pos = as.integer(sub(".*:", "", h)))
}

cluster_leads <- function(leads, gap = WIN) {
  d <- split_hit(sort(unique(leads)))
  d[, chrn := as.integer(chr)]
  setorder(d, chrn, pos)
  cid <- integer(nrow(d)); c <- 0L; pc <- ""; pp <- -Inf
  for (i in seq_len(nrow(d))) {
    if (d$chr[i] != pc || d$pos[i] - pp > gap) c <- c + 1L
    cid[i] <- c; pc <- d$chr[i]; pp <- d$pos[i]
  }
  d[, locus := cid]
  d[, hit1 := paste0(chr, ":", pos)]
  d[]
}

near_any <- function(h, ref) {  # TRUE if hit h within 500 kb of any ref lead on same chr
  s <- split_hit(h); r <- split_hit(unique(ref))
  vapply(seq_len(nrow(s)), function(i) {
    p <- r[chr == s$chr[i], pos]; length(p) > 0 && any(abs(p - s$pos[i]) <= WIN)
  }, logical(1))
}

# r2 tables -> liver arm, platelet arm (pair level). Mirrors scratch/plan-H53/count_n.py
load_arms <- function() {
  t <- fread(file.path(R2DIR, "gene_level_coloc_tier12.csv"))
  pos_genes <- t[susie_state_t12 == "positive", ensembl]
  s <- fread(file.path(R2DIR, "susie_coloc_signal_pairs.tsv.gz"))
  tier <- fread(file.path(FM, "config/gwas_trait_tier.tsv"))[, .(study_name, tier)]
  reg  <- fread(file.path(FM, "config/gwas_registry.tsv"))[, .(study_name, ancestry, trait_type, N_tot, N_cases)]
  s <- merge(s, tier, by.x = "gwas", by.y = "study_name", all.x = TRUE)
  s <- merge(s, reg,  by.x = "gwas", by.y = "study_name", all.x = TRUE)
  s[, chr := as.character(chr)]
  k <- s[PP.H4 > 0.5 & eligible == TRUE]
  L <- k[ensembl %in% pos_genes & tier %in% c(1, 2) & ancestry == "EUR"]
  P <- k[gwas == "MVP_Platelet_EUR"]
  list(all = s, k = k, L = L, P = P, pos_genes = pos_genes)
}

# ---------------------------------------------------------------------------
# Liftover (hg19 -> hg38) with the prespec chain; pilot-coded (first mapping kept,
# must stay on the same chromosome, otherwise NA).
# ---------------------------------------------------------------------------
.CHAIN <- NULL
get_chain <- function() {
  if (is.null(.CHAIN)) {
    suppressPackageStartupMessages({ library(rtracklayer); library(GenomicRanges) })
    f <- file.path(RUN, "tmp/hg19ToHg38.over.chain")
    if (!file.exists(f)) {
      system2("bash", c("-c", shQuote(paste("zcat", shQuote(file.path(CQDIR, "hg19ToHg38.over.chain.gz")), ">", shQuote(paste0(f, ".part")), "&& mv", shQuote(paste0(f, ".part")), shQuote(f)))))
    }
    .CHAIN <<- import.chain(f)
  }
  .CHAIN
}
lift19to38 <- function(chr, pos) {
  if (!length(pos)) return(integer(0))
  chain <- get_chain()
  gr <- GRanges(paste0("chr", chr), IRanges(pos, width = 1), idx = seq_along(pos))
  l <- unlist(liftOver(gr, chain))
  l <- l[!duplicated(l$idx)]
  out <- rep(NA_integer_, length(pos)); out[l$idx] <- start(l)
  ok <- rep(FALSE, length(pos)); ok[l$idx] <- as.character(seqnames(l)) == paste0("chr", chr[1])
  out[!ok] <- NA_integer_
  out
}

source(file.path(FM, "src/coloc_signal_eligibility.R"))  # signal_shared_prop() (hash-checked input)

# ---------------------------------------------------------------------------
# Primary route (prespec "Modality and method", steps 2-5)
# ---------------------------------------------------------------------------
# GWAS side: keys for every SNP of the refit (hg19 order kept). Shareable = SNV that
# lifts to the same chromosome; key = hg38 pos + unordered alleles. Others keep a key
# that cannot match.
gwas_keys <- function(fit) {
  ss <- copy(fit$ss)
  ss[, snp_id := fit$snp_ids]
  ss[, pos38 := lift19to38(as.character(fit$chr), position)]
  ss[, snv := nchar(allele1) == 1 & nchar(allele2) == 1]
  ss[, key := ifelse(snv & !is.na(pos38),
                     paste(pos38, pmin(allele1, allele2), pmax(allele1, allele2), sep = ":"),
                     paste0("g|", snp_id))]
  d <- duplicated(ss$key)
  if (any(d)) ss[d, key := paste0("gdup|", snp_id, "|", seq_len(sum(d)))]
  ss[]
}
gwas_row <- function(fit, gk, idx1) {
  v <- fit$lbf_variable[idx1, ]
  names(v) <- gk$key   # same order as fit$snp_ids
  v
}
# caQTL side: single-effect log ABF per variant (coloc lABF.df2 form), every variant of
# the peak with MAF >= 0.01; final column `null` = 0.
caqtl_lbf <- function(beta, varbeta, W = W_PRIOR) {
  z <- beta / sqrt(varbeta); r <- W^2 / (varbeta + W^2)
  0.5 * log(1 - r) + 0.5 * z^2 * r
}
caqtl_prepare <- function(q) {  # q: rows of ONE peak (data.table from the nominal file)
  q <- q[MAF >= MAF_MIN & is.finite(varbeta) & varbeta > 0 & is.finite(beta)]
  sp <- tstrsplit(q$variant, ":", fixed = TRUE)
  q[, `:=`(pos38 = as.integer(sp[[2]]), ref = sp[[3]], alt = sp[[4]])]
  q[, snv := nchar(ref) == 1 & nchar(alt) == 1]
  q[, key := ifelse(snv, paste(pos38, pmin(ref, alt), pmax(ref, alt), sep = ":"), paste0("q|", variant))]
  d <- duplicated(q$key)
  if (any(d)) q[d, key := paste0("qdup|", variant, "|", seq_len(sum(d)))]
  q[, lbf := caqtl_lbf(beta, varbeta)]
  setorder(q, pos38)
  q[]
}
caqtl_row <- function(qp, lbf = qp$lbf) {
  v <- c(lbf, 0); names(v) <- c(qp$key, "null"); v
}

# One (GWAS signal, peak) test. gk: gwas_keys(); grow: GWAS lbf row; qp: caqtl_prepare();
# lead_snp: hg19 snp id of the GWAS lead (top-alpha of idx1); ld_rows from the fit.
# qp_p / qp_lbf allow NC1 to pass shifted values without touching anything else.
pair_test <- function(grow, gk, qp, ld_rows, lead_snp, p12s = c(P12, P12_ALT),
                      qp_lbf = qp$lbf, qp_p = qp$pvalue) {
  crow <- caqtl_row(qp, qp_lbf)
  shared <- setdiff(intersect(names(grow), names(crow)), "null")
  out <- list(n_shared = length(shared), prop_gwas = NA_real_, prop_caqtl = NA_real_,
              testable = FALSE, reason = NA_character_)
  b1 <- matrix(grow, nrow = 1, dimnames = list(NULL, names(grow)))
  b2 <- matrix(crow, nrow = 1, dimnames = list(NULL, names(crow)))
  if (length(shared)) {
    out$prop_gwas  <- signal_shared_prop(b1, shared, P1)$prop
    out$prop_caqtl <- signal_shared_prop(b2, shared, P2)$prop
  }
  if (length(shared) < MIN_SNPS) { out$reason <- "lt100_shared"; return(out) }
  for (p12 in p12s) {
    r <- suppressWarnings(coloc.bf_bf(b1, b2, p1 = P1, p2 = P2, p12 = p12))
    sm <- r$summary
    tag <- if (p12 == P12) "" else "_p12_1e5"
    if (is.na(sm$PP.H4.abf[1])) {
      out[[paste0("dropped", tag)]] <- TRUE
    } else {
      out[[paste0("dropped", tag)]] <- FALSE
      for (h in 0:4) out[[paste0("PP.H", h, tag)]] <- sm[[paste0("PP.H", h, ".abf")]][1]
    }
  }
  if (isTRUE(out$dropped)) {
    out$reason <- if (out$prop_gwas < OVERLAP_MIN && out$prop_caqtl < OVERLAP_MIN) "overlap_both" else
                  if (out$prop_gwas < OVERLAP_MIN) "overlap_gwas" else "overlap_caqtl"
    return(out)
  }
  out$testable <- TRUE
  # caQTL lead = shared SNP with the smallest caQTL p (ties: first in hg38 order)
  m <- match(shared, qp$key)
  pv <- qp_p[m]
  lk <- shared[which.min(pv)]
  out$caqtl_lead_key <- lk
  out$caqtl_lead_p <- min(pv)
  out$caqtl_lead_variant <- qp$variant[match(lk, qp$key)]
  out$caqtl_lead_impR2 <- qp$imputation_R2[match(lk, qp$key)]
  out$min_caqtl_p_shared <- min(qp$pvalue[m])   # unshifted, for covariates / NC2
  lsnp <- gk$snp_id[match(lk, gk$key)]
  out$caqtl_lead_snp19 <- lsnp
  out$lead_ld_r2 <- if (!is.na(lsnp) && lead_snp %in% rownames(ld_rows) && lsnp %in% colnames(ld_rows))
                      ld_rows[lead_snp, lsnp]^2 else NA_real_
  out
}
