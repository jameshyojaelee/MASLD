#!/usr/bin/env Rscript
# cross_ancestry_replication_lookup.R
# ---------------------------------------------------------------------------
# Cross-Ancestry variant-level REPLICATION LOOKUP (complement to COLOC-based
# cross-ancestry replication, Script 48).
#
# Question: for each EUR-DISCOVERED genome-wide-significant lead locus, does the
#   effect REPLICATE (same direction, Bonferroni-significant one-sided) in the
#   smaller non-EUR MVP GWAS of the SAME phenotype (AFR / AMR / EAS)?
#
# Method (locked; literature-backed):
#   - Discovery-and-lookup replication as used by Vujkovic et al. 2022 Nat Genet
#     (PMID 35654975) for MVP NAFLD, and the directional-generalization framing
#     of Sofer et al. 2017 (PMID 28090672).
#   1. Discovery set = EUR lead loci genome-wide significant (P < 5e-8) in the
#      EUR arm. Leads are taken from an existing *_leadSNPs.tsv when present,
#      else derived here by distance clumping the EUR sumstats (500 kb window),
#      matching the pipeline's EAS lead-SNP convention (src/01_identify_lead_snps.sh).
#   2. Each EUR lead variant is looked up in each non-EUR arm by chr:pos (hg19).
#      Alleles are oriented to the EUR effect allele (direct / swap->flip sign /
#      strand-complement for unambiguous SNVs). If the exact variant is absent,
#      a positional-window proxy is used ONLY when it carries the same allele
#      pair as the EUR lead (i.e. the same variant re-represented nearby), so the
#      effect stays phase-valid; proxy_used is flagged. (A true LD proxy would
#      need a per-ancestry LD panel and is out of scope for this lookup; untested
#      loci instead report window support columns for transparency.)
#   3. Replication test = ONE-SIDED in the direction of the EUR effect:
#         z_dir     = sign(beta_eur) * beta_neur_oriented / se_neur
#         p_oneside = pnorm(z_dir, lower.tail = FALSE)
#      (folds effect direction into the p-value: concordant -> p<0.5, discordant
#       -> p>0.5). Declared REPLICATED when direction-concordant AND
#      p_oneside < 0.05 / k, where k = number of pre-specified EUR loci with a
#      usable (exact or phase-valid proxy) lookup in that non-EUR arm.
#   4. Aggregate binomial SIGN test per non-EUR arm: of the k loci, how many are
#      directionally concordant (expected 50% under the null), one-sided
#      binom.test(n_concordant, k, 0.5, alternative = "greater") -> portfolio
#      transferability p.
#
# Parameterized by phenotype family (matches each MVP EUR arm to the non-EUR
#   arms of the SAME phenotype), so it runs NAFLD, Cirrhosis, and each enzyme
#   (ALT / AST / Albumin / Platelet) uniformly.
#
# Inputs:
#   - GWAS/MR_Data/MVP/mvp_manifest.tsv               (arm -> sumstats/leadsnps/ancestry)
#   - GWAS/finemapping/data/sumstats/<arm>_reformatted_hg19.tsv
#         cols: chromosome position allele1 allele2 beta se pval  (allele1 = effect allele; hg19)
#   - GWAS/finemapping/data/lead_snps/<arm>_leadSNPs.tsv  (optional; CHR BP locus)
#
# Output: GWAS/MR_Data/MVP/results/replication_lookup/
#   - <PHENO>_replication_lookup.tsv        per-(locus x non-EUR-ancestry) rows
#   - cross_ancestry_replication_all.tsv    all phenotypes stacked
#   - replication_summary.tsv               per (pheno x non-EUR-ancestry) summary
#
# Env vars (all optional):
#   MASLD_PROJECT_ROOT  project root (default below)
#   REPL_PHENOS         comma list of phenotypes, or "ALL" (default:
#                       NAFLD,Cirrhosis,ALT,AST,Albumin,Platelet)
#   REPL_GWSIG          discovery p threshold           (default 5e-8)
#   REPL_CLUMP_KB       distance-clump half-window kb   (default 500)
#   REPL_PROXY_KB       proxy search half-window kb     (default 250; 0 = exact only)
#   REPL_OUTDIR         output dir (default results/replication_lookup under MVP)
#   GENE_BED_HG19       optional hg19 gene bed (chr start end gene) for nearest_gene
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASE     <- Sys.getenv("MASLD_PROJECT_ROOT",
                       "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE, "GWAS/finemapping")
MVP_DIR  <- file.path(BASE, "GWAS/MR_Data/MVP")
MANIFEST <- file.path(MVP_DIR, "mvp_manifest.tsv")

GWSIG    <- as.numeric(Sys.getenv("REPL_GWSIG",   "5e-8"))
CLUMP_KB <- as.numeric(Sys.getenv("REPL_CLUMP_KB", "500"))
PROXY_KB <- as.numeric(Sys.getenv("REPL_PROXY_KB", "250"))
CLUMP_BP <- CLUMP_KB * 1000L
PROXY_BP <- PROXY_KB * 1000L

OUTDIR   <- Sys.getenv("REPL_OUTDIR",
                       file.path(MVP_DIR, "results", "replication_lookup"))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

phenos_env <- Sys.getenv("REPL_PHENOS", "NAFLD,Cirrhosis,ALT,AST,Albumin,Platelet")
GENE_BED   <- Sys.getenv("GENE_BED_HG19", "")

NONEUR <- c("AFR", "AMR", "EAS")   # non-EUR target ancestries, in report order

cat("=== Cross-Ancestry Replication Lookup ===\n")
cat("Start:", format(Sys.time()), "\n")
cat(sprintf("GW-sig p<%.0e | clump=%gkb | proxy=%gkb | out=%s\n\n",
            GWSIG, CLUMP_KB, PROXY_KB, OUTDIR))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
comp_base <- function(a) {
  # single-nucleotide complement; multi-char (indel) -> NA
  out <- rep(NA_character_, length(a))
  m <- c(A = "T", T = "A", C = "G", G = "C")
  single <- nchar(a) == 1L
  out[single] <- m[a[single]]
  out
}
is_ambiguous <- function(a1, a2) {
  # A/T or C/G SNVs are strand-ambiguous -> do NOT use complement matching
  key <- paste0(pmin(a1, a2), pmax(a1, a2))
  key %in% c("AT", "CG")
}

# Resolve a sumstats path from a manifest row (out_sumstats is relative to FM_DIR)
sumstats_path <- function(arm_row) {
  file.path(FM_DIR, arm_row$out_sumstats)
}
leadsnps_path <- function(arm_row) {
  file.path(FM_DIR, arm_row$leadsnps)
}

# awk pre-filter a reformatted sumstats file to a set of chromosomes, then fread.
read_sumstats_chr <- function(path, chrs) {
  if (!file.exists(path)) return(NULL)
  cols <- c("chromosome", "position", "allele1", "allele2", "beta", "se", "pval")
  if (length(chrs) == 0) {
    dt <- fread(path, select = cols, showProgress = FALSE)
  } else {
    chrs <- unique(as.character(chrs))
    cond <- paste(sprintf('$1=="%s"', chrs), collapse = " || ")
    cmd  <- sprintf("awk -F'\\t' 'NR==1 || (%s)' %s", cond, shQuote(path))
    dt   <- fread(cmd = cmd, select = cols, showProgress = FALSE)
  }
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  dt[, chromosome := as.character(chromosome)]
  dt[, allele1 := toupper(allele1)]
  dt[, allele2 := toupper(allele2)]
  dt
}

# Distance-clump a GW-sig table into independent lead rows (greedy, min-p first).
clump_leads <- function(sig, half_bp) {
  if (nrow(sig) == 0) return(sig[0])
  setorder(sig, chromosome, pval)
  keep <- vector("list", 0L)
  for (chr in unique(sig$chromosome)) {
    s <- sig[chromosome == chr]            # already p-sorted within chr
    taken_pos <- numeric(0)
    for (i in seq_len(nrow(s))) {
      p <- s$position[i]
      if (length(taken_pos) == 0 || all(abs(taken_pos - p) > half_bp)) {
        keep[[length(keep) + 1L]] <- s[i]
        taken_pos <- c(taken_pos, p)
      }
    }
  }
  rbindlist(keep)
}

# Optional nearest-gene annotation (hg19). Returns character vector.
load_gene_bed <- function(path) {
  if (path == "" || !file.exists(path)) return(NULL)
  g <- fread(path, header = FALSE)
  setnames(g, 1:4, c("chr", "start", "end", "gene"))
  g[, chr := sub("^chr", "", as.character(chr))]
  g[, mid := (start + end) / 2]
  g
}
nearest_gene <- function(gene_bed, chr, pos) {
  if (is.null(gene_bed)) return(NA_character_)
  cand <- gene_bed[chr == as.character(chr)]
  if (nrow(cand) == 0) return(NA_character_)
  # overlap first, else nearest midpoint
  ov <- cand[start <= pos & end >= pos]
  if (nrow(ov) > 0) return(ov$gene[which.min(abs(ov$mid - pos))])
  cand$gene[which.min(abs(cand$mid - pos))]
}

# ---------------------------------------------------------------------------
# Read manifest, build phenotype families (MVP arms: MVP_<PHENO>_<ANC>)
# ---------------------------------------------------------------------------
if (!file.exists(MANIFEST)) stop("Manifest not found: ", MANIFEST)
man <- fread(MANIFEST)
man <- man[grepl("^MVP_", gwas_name)]
man[, pheno    := sub("^MVP_(.*)_[A-Z]+$", "\\1", gwas_name)]
man[, ancestry := sub("^MVP_.*_([A-Z]+)$", "\\1", gwas_name)]

fams_avail <- man[ancestry == "EUR", unique(pheno)]
fams_avail <- intersect(fams_avail, man[ancestry %in% NONEUR, unique(pheno)])

if (toupper(phenos_env) == "ALL") {
  phenos <- fams_avail
} else {
  phenos <- trimws(strsplit(phenos_env, ",")[[1]])
  miss   <- setdiff(phenos, fams_avail)
  if (length(miss) > 0)
    cat("  NOTE: requested phenotypes with no EUR+non-EUR pair skipped:",
        paste(miss, collapse = ", "), "\n")
  phenos <- intersect(phenos, fams_avail)
}
cat("Phenotype families to run:", paste(phenos, collapse = ", "), "\n\n")

gene_bed <- load_gene_bed(GENE_BED)
if (!is.null(gene_bed)) cat("Nearest-gene annotation from:", GENE_BED, "\n")

# ---------------------------------------------------------------------------
# Per-phenotype replication
# ---------------------------------------------------------------------------
all_rows    <- list()
all_summary <- list()

for (ph in phenos) {
  cat("========================================================\n")
  cat("Phenotype:", ph, "\n")
  eur_row <- man[pheno == ph & ancestry == "EUR"]
  if (nrow(eur_row) == 0) { cat("  no EUR arm; skip\n"); next }
  eur_row  <- eur_row[1]
  eur_arm  <- eur_row$gwas_name
  eur_path <- sumstats_path(eur_row)
  if (!file.exists(eur_path)) { cat("  EUR sumstats missing:", eur_path, "; skip\n"); next }

  # --- Discovery lead loci (EUR) ---------------------------------------------
  lp <- leadsnps_path(eur_row)
  if (file.exists(lp)) {
    ld <- fread(lp)
    setnames(ld, tolower(names(ld)))
    lead_chr <- as.character(ld$chr)
    cat(sprintf("  discovery: %d leads from lead-SNP file %s\n", nrow(ld), basename(lp)))
    eur <- read_sumstats_chr(eur_path, unique(lead_chr))
    if (is.null(eur)) { cat("  cannot read EUR sumstats; skip\n"); next }
    setkey(eur, chromosome, position)
    key_dt <- data.table(chromosome = as.character(ld$chr), position = as.integer(ld$bp))
    leads  <- eur[key_dt, on = .(chromosome, position), nomatch = 0L]
    leads  <- leads[pval < GWSIG]                    # enforce GW-sig discovery
  } else {
    cat("  discovery: no lead-SNP file; clumping EUR sumstats (this reads the full arm)\n")
    eur_full <- fread(eur_path,
                      select = c("chromosome","position","allele1","allele2","beta","se","pval"),
                      showProgress = FALSE)
    eur_full[, chromosome := as.character(chromosome)]
    sig <- eur_full[pval < GWSIG]
    sig[, allele1 := toupper(allele1)][, allele2 := toupper(allele2)]
    rm(eur_full); gc()
    leads <- clump_leads(sig, CLUMP_BP)
  }

  if (nrow(leads) == 0) { cat("  no genome-wide-significant EUR leads; skip\n"); next }
  # one row per locus; canonical fields
  leads <- leads[, .(chromosome = as.character(chromosome), position = as.integer(position),
                     a1_eur = toupper(allele1), a2_eur = toupper(allele2),
                     beta_eur = beta, se_eur = se, p_eur = pval)]
  setorder(leads, chromosome, position)
  leads[, locus_id  := paste0(chromosome, ":", position)]
  leads[, variant_id := paste0(chromosome, ":", position, ":", a1_eur, ":", a2_eur)]
  leads[, nearest_gene := vapply(seq_len(.N),
                                 function(i) nearest_gene(gene_bed, chromosome[i], position[i]),
                                 character(1))]
  cat(sprintf("  %d GW-sig EUR lead loci (chr: %s)\n",
              nrow(leads), paste(sort(unique(as.integer(leads$chromosome))), collapse = ",")))

  lead_chr <- unique(leads$chromosome)

  # --- Lookup in each non-EUR arm --------------------------------------------
  pheno_rows <- list()
  for (anc in NONEUR) {
    trow <- man[pheno == ph & ancestry == anc]
    if (nrow(trow) == 0) next
    trow  <- trow[1]
    tarm  <- trow$gwas_name
    tpath <- sumstats_path(trow)
    if (!file.exists(tpath)) { cat("   ", tarm, ": sumstats missing; skip\n"); next }

    tgt <- read_sumstats_chr(tpath, lead_chr)
    if (is.null(tgt) || nrow(tgt) == 0) { cat("   ", tarm, ": no rows on lead chr; skip\n"); next }
    setkey(tgt, chromosome, position)

    rows <- vector("list", nrow(leads))
    for (i in seq_len(nrow(leads))) {
      L <- leads[i]
      chr <- L$chromosome; pos <- L$position
      a1e <- L$a1_eur; a2e <- L$a2_eur
      amb <- is_ambiguous(a1e, a2e)

      # candidate exact-position rows
      exact <- tgt[.(chr, pos), nomatch = 0L]
      matched <- "none"; proxy_used <- FALSE; proxy_dist <- NA_integer_
      a1n <- NA_character_; a2n <- NA_character_
      beta_or <- NA_real_; se_n <- NA_real_; p_two <- NA_real_

      orient <- function(cand) {
        # return oriented beta (to EUR effect allele) or NA if alleles don't match
        for (j in seq_len(nrow(cand))) {
          x1 <- cand$allele1[j]; x2 <- cand$allele2[j]; b <- cand$beta[j]
          if (x1 == a1e && x2 == a2e) return(list(b =  b, a1 = x1, a2 = x2, se = cand$se[j], p = cand$pval[j]))
          if (x1 == a2e && x2 == a1e) return(list(b = -b, a1 = x1, a2 = x2, se = cand$se[j], p = cand$pval[j]))
          if (!amb) {
            c1 <- comp_base(x1); c2 <- comp_base(x2)
            if (!is.na(c1) && !is.na(c2)) {
              if (c1 == a1e && c2 == a2e) return(list(b =  b, a1 = x1, a2 = x2, se = cand$se[j], p = cand$pval[j]))
              if (c1 == a2e && c2 == a1e) return(list(b = -b, a1 = x1, a2 = x2, se = cand$se[j], p = cand$pval[j]))
            }
          }
        }
        NULL
      }

      window_min_p <- NA_real_; n_in_window <- 0L
      if (nrow(exact) > 0) {
        o <- orient(exact)
        if (!is.null(o)) {
          matched <- "exact"; beta_or <- o$b; se_n <- o$se; p_two <- o$p; a1n <- o$a1; a2n <- o$a2
        }
      }
      if (matched == "none" && PROXY_BP > 0) {
        win <- tgt[chromosome == chr & abs(position - pos) <= PROXY_BP]
        n_in_window <- nrow(win)
        if (n_in_window > 0) window_min_p <- min(win$pval, na.rm = TRUE)
        # phase-valid proxy = same allele pair as EUR lead, nearest to lead
        if (n_in_window > 0) {
          setorder(win, position)
          win[, dist := abs(position - pos)]
          setorder(win, dist, pval)
          o <- orient(win)
          if (!is.null(o)) {
            matched <- "proxy"; proxy_used <- TRUE
            beta_or <- o$b; se_n <- o$se; p_two <- o$p; a1n <- o$a1; a2n <- o$a2
            # distance of the chosen proxy
            hit <- win[allele1 == a1n & allele2 == a2n]
            proxy_dist <- as.integer(hit$dist[1])
          }
        }
      } else if (matched == "none") {
        win <- tgt[chromosome == chr & abs(position - pos) <= PROXY_BP]
        n_in_window <- nrow(win)
        if (n_in_window > 0) window_min_p <- min(win$pval, na.rm = TRUE)
      }

      # directional test
      dir_conc <- NA; p_one <- NA_real_
      if (matched %in% c("exact", "proxy") && !is.na(beta_or) && !is.na(se_n) && se_n > 0) {
        z_dir <- sign(L$beta_eur) * beta_or / se_n
        p_one <- pnorm(z_dir, lower.tail = FALSE)
        dir_conc <- (sign(beta_or) == sign(L$beta_eur)) && beta_or != 0
      }

      rows[[i]] <- data.table(
        pheno = ph, eur_arm = eur_arm, neur_arm = tarm, ancestry = anc,
        chr = chr, pos = pos, locus_id = L$locus_id, variant_id = L$variant_id,
        nearest_gene = L$nearest_gene,
        a1_eur = a1e, a2_eur = a2e, beta_eur = L$beta_eur, se_eur = L$se_eur, p_eur = L$p_eur,
        matched = matched, proxy_used = proxy_used, proxy_dist_bp = proxy_dist,
        a1_neur = a1n, a2_neur = a2n, beta_neur_oriented = beta_or, se_neur = se_n,
        p_neur_twosided = p_two, p_oneside = p_one, direction_concordant = dir_conc,
        window_min_p = window_min_p, n_snps_in_window = n_in_window
      )
    }
    rt <- rbindlist(rows)

    # k, Bonferroni threshold, replication flag (per non-EUR arm)
    usable <- rt[matched %in% c("exact", "proxy")]
    k      <- nrow(usable)
    bonf   <- if (k > 0) 0.05 / k else NA_real_
    rt[, k_tested := k]
    rt[, bonf_threshold := bonf]
    rt[, replicated := (matched %in% c("exact","proxy")) &
                       !is.na(direction_concordant) & direction_concordant &
                       !is.na(p_oneside) & p_oneside < bonf]
    rt[is.na(replicated), replicated := FALSE]

    pheno_rows[[anc]] <- rt

    # per-arm summary + binomial sign test
    n_conc <- sum(usable$direction_concordant, na.rm = TRUE)
    sign_p <- if (k > 0) binom.test(n_conc, k, 0.5, alternative = "greater")$p.value else NA_real_
    all_summary[[paste(ph, anc)]] <- data.table(
      pheno = ph, eur_arm = eur_arm, neur_arm = tarm, ancestry = anc,
      n_eur_leads = nrow(leads), k_tested = k,
      n_exact = sum(rt$matched == "exact"), n_proxy = sum(rt$matched == "proxy"),
      n_untested = sum(rt$matched == "none"),
      n_concordant = n_conc,
      prop_concordant = if (k > 0) round(n_conc / k, 4) else NA_real_,
      n_replicated = sum(rt$replicated),
      bonf_threshold = bonf, sign_test_p = sign_p
    )
    cat(sprintf("    %-16s N=%d(cases=%s): k=%d exact=%d proxy=%d untested=%d | conc=%d/%d (p_sign=%s) | replicated=%d\n",
                tarm, trow$N_tot, trow$N_cases, k,
                sum(rt$matched=="exact"), sum(rt$matched=="proxy"), sum(rt$matched=="none"),
                n_conc, k, formatC(sign_p, format = "g", digits = 3), sum(rt$replicated)))
    rm(tgt); gc()
  }

  if (length(pheno_rows) > 0) {
    ptab <- rbindlist(pheno_rows)
    fwrite(ptab, file.path(OUTDIR, paste0(ph, "_replication_lookup.tsv")), sep = "\t")
    all_rows[[ph]] <- ptab
  }
  cat("\n")
}

# ---------------------------------------------------------------------------
# Combined outputs
# ---------------------------------------------------------------------------
if (length(all_rows) > 0) {
  allt <- rbindlist(all_rows, fill = TRUE)
  fwrite(allt, file.path(OUTDIR, "cross_ancestry_replication_all.tsv"), sep = "\t")
  cat("Wrote:", file.path(OUTDIR, "cross_ancestry_replication_all.tsv"),
      sprintf("(%d locus x ancestry rows)\n", nrow(allt)))
}
if (length(all_summary) > 0) {
  summ <- rbindlist(all_summary, fill = TRUE)
  fwrite(summ, file.path(OUTDIR, "replication_summary.tsv"), sep = "\t")
  cat("Wrote:", file.path(OUTDIR, "replication_summary.tsv"), "\n\n")
  cat("=== Per-ancestry replication summary ===\n")
  print(summ)
}

cat("\nDone:", format(Sys.time()), "\n")
