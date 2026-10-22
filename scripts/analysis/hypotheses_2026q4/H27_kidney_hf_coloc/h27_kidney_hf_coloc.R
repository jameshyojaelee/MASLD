#!/usr/bin/env Rscript
# H27: kidney function (eGFR), CKD and heart failure colocalization at anchored
# MASLD loci. Implements PRESPEC_H27.md in the output directory. Modality: GWAS
# summary statistics on both sides; unit: independent locus. All results are
# provisional on the COLOC release susie_coloc_r2_20260924T053244Z.
#
# Usage:
#   Rscript h27_kidney_hf_coloc.R prep    <repo> <out>
#   Rscript h27_kidney_hf_coloc.R coloc   <repo> <out> <worker_id> <n_workers>
#   Rscript h27_kidney_hf_coloc.R combine <repo> <out>

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]; REPO <- args[2]; OUT <- args[3]

suppressMessages({
  library(data.table)
  library(coloc)
  library(susieR)
})
stopifnot(as.character(packageVersion("coloc")) == "5.2.3",
          as.character(packageVersion("susieR")) == "0.14.2")

FM   <- file.path(REPO, "GWAS/finemapping")
REL  <- file.path(FM, "results/susie_coloc_r2_20260924T053244Z")
H10_LOCI <- file.path(REPO, "Analysis/Hypotheses/H10-20261007T210859Z/loci.tsv")
WORK <- file.path(OUT, "work"); RES <- file.path(OUT, "results")
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
dir.create(RES,  showWarnings = FALSE, recursive = TRUE)

# --- release parameters (06_susie_coloc.R) -----------------------------------
COLOC_P1 <- 1e-4; COLOC_P2 <- 1e-4; COLOC_P12 <- 5e-6
MIN_SNPS <- 100L; MIN_TRIPLE_SNPS <- 50L
LD_REGULARIZE <- 1e-3; SUSIE_L <- 10L
LAMBDA_S_HIGH <- 0.20
PP_THR <- 0.5
WINDOW <- 500000L
MHC <- c(6, 25e6, 34e6)
SEED_PERM <- 20260923L
N_PERM <- 200L
# release palindrome AF-guard constants (06_susie_coloc.R lines 84-86)
AF_PAL_MAF_MAX <- 0.40; AF_PAL_MAX_DIST <- 0.15; AF_PAL_MARGIN <- 0.10

# --- outcomes (PRESPEC section 2) ---------------------------------------------
OUTCOMES <- list(
  eGFR = list(acc = "GCST90103634", type = "quant", N = 1004040, Ncase = NA,
              file = "GCST90103634_buildGRCh37.tsv"),
  CKD  = list(acc = "GCST008065", type = "cc", N = 64164 + 561055, Ncase = 64164,
              file = "CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz"),
  HF   = list(acc = "GCST90728695", type = "cc", N = 139533 + 1568809, Ncase = 139533,
              file = "GCST90728695.tsv.gz")
)
n_eff <- function(N, Ncase) if (is.na(Ncase) || Ncase <= 0) N else 4 / (1 / Ncase + 1 / (N - Ncase))

registry <- fread(file.path(FM, "config/gwas_registry.tsv"))
af_src   <- fread(file.path(FM, "config/gwas_af_sources.tsv"))
tier     <- fread(file.path(FM, "config/gwas_trait_tier.tsv"))

liver_path <- function(study) {
  hit <- af_src[study_name == study & file.exists(af_sumstats_path)]
  if (nrow(hit) == 1L) return(hit$af_sumstats_path)
  file.path(FM, registry[study_name == study]$sumstats_path)
}
liver_spec <- function(study) {
  r <- registry[study_name == study]
  stopifnot(nrow(r) == 1L)
  list(type = ifelse(r$trait_type == "binary", "cc", "quant"), N = r$N_tot,
       Ncase = r$N_cases, trait = tier[study_name == study]$trait)
}

# ==============================================================================
if (stage == "prep") {
  cat("== prep\n")
  # ---- 1. loci (H10 rule) -----------------------------------------------------
  genes <- fread(file.path(REL, "gene_level_coloc_tier12.csv"))
  pos_genes <- genes[susie_state_t12 == "positive"]
  stopifnot(nrow(pos_genes) == 436L)
  t12 <- tier[tier <= 2]$study_name
  pairs <- fread(file.path(REL, "susie_coloc_signal_pairs.tsv.gz"))
  pairs <- pairs[eligible %in% TRUE & PP.H4 > PP_THR & gwas %in% t12 & ensembl %in% pos_genes$ensembl]
  pairs[, c("hchr", "hpos") := tstrsplit(hit1, ":", fixed = TRUE, keep = 1:2)]
  pairs[, `:=`(hchr = as.integer(hchr), hpos = as.integer(hpos))]
  cand <- unique(pairs[, .(hchr, hpos)])[order(hchr, hpos)]
  parent <- seq_len(nrow(cand))
  find <- function(i) { while (parent[i] != i) { parent[i] <<- parent[parent[i]]; i <- parent[i] }; as.integer(i) }
  unite <- function(i, j) parent[find(i)] <<- find(j)
  for (i in seq_len(nrow(cand))[-1]) {
    if (cand$hchr[i] == cand$hchr[i - 1L] && cand$hpos[i] - cand$hpos[i - 1L] <= WINDOW) unite(i, i - 1L)
  }
  cand[, key := paste(hchr, hpos)]
  idx <- setNames(seq_len(nrow(cand)), cand$key)
  for (g in split(pairs, pairs$ensembl)) {
    ids <- idx[paste(g$hchr, g$hpos)]
    if (length(ids) > 1) for (j in ids[-1]) unite(ids[1], j)
  }
  cand[, root := vapply(seq_len(.N), find, integer(1))]
  pairs[, root := cand$root[idx[paste(hchr, hpos)]]]
  sym <- pairs[, .(gene = gene[1]), by = ensembl]
  loci <- cand[, .(chr = hchr[1], start = min(hpos), end = max(hpos)), by = root]
  # radix sort = C-locale byte order, as Python sorted() in the H10 script
  loci[, genes := vapply(root, function(r) {
    e <- unique(pairs[root == r]$ensembl)
    paste(sort(as.character(sym[ensembl %in% e]$gene), method = "radix"), collapse = ",")
  }, character(1))]
  # H10 numbers loci in the order clusters first appear in the sorted positions
  loci[, first := vapply(root, function(r) min(which(cand$root == r)), integer(1))]
  setorder(loci, first)
  loci[, locus_id := sprintf("L%03d", seq_len(.N))]
  h10 <- fread(H10_LOCI)
  chk <- merge(loci[, .(locus_id, chr, start, end, genes)],
               h10[, .(locus_id, chr, start, end, genes)], by = "locus_id", all = TRUE)
  ok <- nrow(loci) == 289L && nrow(chk) == 289L &&
        all(chk$chr.x == chk$chr.y & chk$start.x == chk$start.y &
            chk$end.x == chk$end.y & chk$genes.x == chk$genes.y)
  if (!isTRUE(ok)) stop("locus rebuild does not match H10 loci.tsv")
  cat("loci rebuilt:", nrow(loci), "(matches H10)\n")
  loci[, mhc := chr == MHC[1] & end >= MHC[2] & start <= MHC[3]]

  # ---- 2. EUR liver GWAS and lead per locus ------------------------------------
  eur <- registry[ancestry == "EUR"]$study_name
  pe <- pairs[gwas %in% eur]
  setorder(pe, -PP.H4, gwas, hpos)
  lead <- pe[, .SD[1], by = root][, .(root, liver_study = gwas, lead_chr = hchr,
                                       lead_pos = hpos, lead_pp4 = PP.H4, lead_gene = gene)]
  loci <- merge(loci, lead, by = "root", all.x = TRUE)
  loci[, status := fifelse(mhc, "not_tested_mhc",
                    fifelse(is.na(liver_study), "not_testable_no_eur_liver_signal", "test"))]
  setorder(loci, locus_id)
  # cohort of the liver GWAS as named in the registry; PubMed-ID-named studies
  # have no cohort list on disk and are recorded as such (not guessed)
  loci[, liver_cohort := fifelse(is.na(liver_study), NA_character_,
                          fifelse(grepl("UKBB", liver_study), "UKB",
                          fifelse(grepl("^MVP_", liver_study), "MVP",
                          fifelse(grepl("^FinnGen_", liver_study), "FinnGen",
                          fifelse(grepl("deCode", liver_study), "deCODE",
                          fifelse(grepl("Intermountain", liver_study), "Intermountain",
                                  "not_recorded_on_disk"))))))]
  fwrite(loci[, !c("root", "first")], file.path(RES, "loci_h27.tsv"), sep = "\t")
  cat("status:", paste(names(table(loci$status)), table(loci$status), collapse = "; "), "\n")

  # ---- 3. tasks: test loci and controls -----------------------------------------
  tasks <- loci[status == "test", .(task_id = locus_id, role = "test", chr = lead_chr,
                                    lead = lead_pos, liver_study,
                                    traits = "eGFR,CKD,HF", genes)]
  shroom <- loci[grepl("(^|,)SHROOM3(,|$)", genes)]
  stopifnot(nrow(shroom) == 1L)
  # GPT negative-control lead: smallest MVP_ALT_EUR p in gene +/- 100 kb (hg19)
  mvp <- fread(liver_path("MVP_ALT_EUR"), select = c("chromosome", "position", "pval"))
  g <- mvp[chromosome == 8 & position >= 145728356 - 100000 & position <= 145732557 + 100000]
  gpt_lead <- g[which.min(pval)]$position
  cat("GPT lead", gpt_lead, "p", min(g$pval), "\n")
  rm(mvp, g); gc(FALSE)
  ctrl <- data.table(
    task_id = c("PC2_TM6SF2", "PC3_ABO", "NC1_GPT"),
    role    = c("positive_control", "positive_control", "negative_control"),
    chr     = c(19L, 9L, 8L),
    lead    = c(19379549L, 136137065L, as.integer(gpt_lead)),
    liver_study = "MVP_ALT_EUR",
    traits  = c("eGFR", "HF", "eGFR,CKD,HF"),
    genes   = c("TM6SF2", "ABO", "GPT"))
  tasks <- rbind(tasks, ctrl)
  tasks[, `:=`(start = pmax(1L, lead - WINDOW), end = lead + WINDOW)]
  setorder(tasks, chr, lead)
  fwrite(tasks, file.path(WORK, "tasks.tsv"), sep = "\t")
  fwrite(data.table(control = "PC1_SHROOM3", locus_id = shroom$locus_id, trait = "eGFR"),
         file.path(WORK, "pc1.tsv"), sep = "\t")

  # rows of d inside any window of w (non-equi join; row indices, no duplicates)
  in_windows <- function(d, w) {
    d[, row_id__ := .I]
    ids <- d[w, on = .(chromosome = chr, position >= start, position <= end),
             nomatch = NULL, x.row_id__]
    d[, row_id__ := NULL]
    sort(unique(ids))
  }

  # ---- 4. extract windows: liver side --------------------------------------------
  for (st in unique(tasks$liver_study)) {
    p <- liver_path(st)
    cat("liver", st, p, "\n")
    cols <- c("chromosome", "position", "allele1", "allele2", "beta", "se", "pval")
    hdr <- names(fread(p, nrows = 0))
    d <- fread(p, select = intersect(c(cols, "af"), hdr))
    if (!"af" %in% names(d)) d[, af := NA_real_]
    d[, chromosome := suppressWarnings(as.integer(chromosome))]
    d <- d[!is.na(chromosome)]
    d <- d[in_windows(d, tasks[liver_study == st, .(chr, start, end)])]
    d[, `:=`(allele1 = toupper(allele1), allele2 = toupper(allele2))]
    fwrite(d, file.path(WORK, paste0("liver_", st, ".tsv.gz")), sep = "\t")
    rm(d); gc(FALSE)
  }

  # ---- 5. extract windows: outcomes ------------------------------------------------
  win <- tasks[, .(chr, start, end)]
  read_outcome <- function(tr) {
    f <- file.path(OUT, "inputs/outcomes", OUTCOMES[[tr]]$acc, OUTCOMES[[tr]]$file)
    if (tr == "eGFR") {
      d <- fread(f, select = c("chromosome", "base_pair_location", "effect_allele", "other_allele",
                               "beta", "standard_error", "effect_allele_frequency", "p_value", "n"))
      setnames(d, c("chromosome", "position", "ea", "oa", "beta", "se", "eaf", "pval", "n"))
    } else if (tr == "CKD") {
      d <- fread(f, select = c("Chr", "Pos_b37", "Allele1", "Allele2", "Effect", "StdErr", "Freq1", "P-value", "n_total_sum"))
      setnames(d, c("chromosome", "position", "ea", "oa", "beta", "se", "eaf", "pval", "n"))
    } else {
      d <- fread(f, select = c("chromosome", "base_pair_location", "effect_allele", "other_allele",
                               "beta", "standard_error", "effect_allele_frequency", "p_value", "N_total"))
      setnames(d, c("chromosome", "position", "ea", "oa", "beta", "se", "eaf", "pval", "n"))
    }
    d[, chromosome := suppressWarnings(as.integer(chromosome))]
    d <- d[!is.na(chromosome)]
    d <- d[in_windows(d, win)]
    d[, `:=`(ea = toupper(ea), oa = toupper(oa))]
    d
  }
  for (tr in names(OUTCOMES)) {
    d <- read_outcome(tr)
    cat("outcome", tr, "rows in windows", nrow(d), "\n")
    fwrite(d, file.path(WORK, paste0("outcome_", tr, ".tsv.gz")), sep = "\t")
    rm(d); gc(FALSE)
  }

  # ---- 6. published lists ----------------------------------------------------------
  suppressMessages(library(readxl))
  xl_j <- file.path(OUT, "inputs/published_lists/jamialahmadi2024_MOESM3.xlsx")
  xl_v <- file.path(OUT, "inputs/published_lists/vujkovic2022_MOESM4.xlsx")
  rd <- function(f, s, skip) as.data.table(suppressMessages(read_excel(f, sheet = s, skip = skip, col_types = "text")))
  j17 <- rd(xl_j, "S17", 1)
  j17 <- unique(j17[!is.na(SNP) & !is.na(POS), .(rsid = SNP, genes = Locus, chr = as.integer(CHR), pos = as.integer(POS))])
  jl <- rbind(j17[, .(source = "Jamialahmadi2024_S17", list = "kidney", rsid, genes, chr, pos, trait = "pPRS variant tested vs CKF/HF")],
              j17[, .(source = "Jamialahmadi2024_S17", list = "HF", rsid, genes, chr, pos, trait = "pPRS variant tested vs CKF/HF")])
  # rsid -> hg19 map from Vujkovic S2 (scan), S29, S32
  v2 <- as.matrix(suppressMessages(read_excel(xl_v, sheet = "2", col_names = FALSE, col_types = "text")))
  m2 <- rbindlist(lapply(seq_len(nrow(v2)), function(i) {
    r <- v2[i, ]; j <- which(grepl("^rs[0-9]+$", r))[1]
    if (is.na(j) || j + 2 > length(r)) return(NULL)
    data.table(rsid = r[j], chr = suppressWarnings(as.integer(r[j + 1])), pos = suppressWarnings(as.integer(r[j + 2])))
  }))
  v29 <- rd(xl_v, "29", 2); v32 <- rd(xl_v, "32", 2); v30 <- rd(xl_v, "30", 2); v33 <- rd(xl_v, "33", 2)
  map19 <- unique(rbind(m2, v29[, .(rsid = RSID, chr = as.integer(Chr), pos = as.integer(Position))],
                        v32[, .(rsid = SNP, chr = as.integer(CHR), pos = as.integer(Position))])[!is.na(pos)], by = "rsid")
  # PRESPEC section 5 terms; \\brenal excludes "adrenal"
  kid_re <- "kidney|\\brenal|glomerul|nephr|creatinine|urea|cystatin|glomerular filtration|urinary albumin"
  hf_re  <- "heart failure|cardiomyopathy"
  is_kid <- function(x) grepl(kid_re, x, ignore.case = TRUE)
  v29[, P := suppressWarnings(as.numeric(P))]
  a <- v29[Pheno %in% c("Creatinine", "Estimated Glomerular Filtration Rate (eGFR)", "Blood Urea Nitrogen (BUN)") & P < 3.09e-5,
           .(source = "Vujkovic2022_S29", list = "kidney", rsid = RSID,
             genes = paste(Nominated.Gene, Biological.Prior, sep = ";"), chr = as.integer(Chr), pos = as.integer(Position), trait = Pheno)]
  v30 <- v30[!is.na(SNP) & grepl("^rs", SNP)]
  b <- rbind(v30[is_kid(Phenotype), .(source = "Vujkovic2022_S30", list = "kidney", rsid = SNP,
                                      genes = paste(`Nominated Gene`, `Biological Prior`, sep = ";"), trait = Phenotype)],
             v30[grepl(hf_re, Phenotype, ignore.case = TRUE), .(source = "Vujkovic2022_S30", list = "HF", rsid = SNP,
                                      genes = paste(`Nominated Gene`, `Biological Prior`, sep = ";"), trait = Phenotype)])
  b <- merge(b, map19, by = "rsid", all.x = TRUE)
  v32 <- v32[!is.na(SNP) & grepl("^rs", SNP)]
  c32 <- rbind(v32[is_kid(Trait), .(source = "Vujkovic2022_S32", list = "kidney", rsid = SNP,
                                    genes = paste(Nominated.Gene, Biological.Prior, sep = ";"),
                                    chr = as.integer(CHR), pos = as.integer(Position), trait = Trait)],
               v32[grepl(hf_re, Trait, ignore.case = TRUE), .(source = "Vujkovic2022_S32", list = "HF", rsid = SNP,
                                    genes = paste(Nominated.Gene, Biological.Prior, sep = ";"),
                                    chr = as.integer(CHR), pos = as.integer(Position), trait = Trait)])
  v33 <- v33[!is.na(gwasSNP) & is_kid(Trait)]
  d33 <- v33[, .(source = "Vujkovic2022_S33", list = "kidney", rsid = gwasSNP,
                 genes = paste(Nominated.Gene, Biological.Prior, sep = ";"),
                 chr = as.integer(sub("^chr", "", tstrsplit(gwasLocus, ":")[[1]])),
                 pos = as.integer(tstrsplit(gwasLocus, ":")[[2]]), trait = Trait)]
  pub <- rbind(jl, a, b[, .(source, list, rsid, genes, chr, pos, trait)], c32, d33)
  fwrite(pub, file.path(RES, "published_lists.tsv"), sep = "\t")
  cat("published list rows:", nrow(pub), "| unplaced (no hg19):", sum(is.na(pub$pos)), "\n")
  cat("prep done\n")
}

# ==============================================================================
# helpers for coloc stages
# panel_freq_of() and strand_verdict_vs_panel() are copied verbatim from the
# release GWAS/finemapping/src/06_susie_coloc.R (lines 76-98); that file cannot be
# sourced because it runs its main loop on source.
panel_freq_of <- function(allele, other, bim_a1, bim_a2, af_a1) {
  pair_ok <- (allele == bim_a1 & other == bim_a2) | (allele == bim_a2 & other == bim_a1)
  fifelse(!pair_ok | is.na(pair_ok), NA_real_, fifelse(allele == bim_a1, af_a1, 1 - af_a1))
}
strand_verdict_vs_panel <- function(obs_af, panel_af) {
  d_same <- abs(obs_af - panel_af); d_opp <- abs(obs_af - (1 - panel_af))
  testable <- is.finite(obs_af) & is.finite(panel_af) &
    pmin(obs_af, 1 - obs_af) <= AF_PAL_MAF_MAX & pmin(panel_af, 1 - panel_af) <= AF_PAL_MAF_MAX &
    pmin(d_same, d_opp) <= AF_PAL_MAX_DIST
  fifelse(!testable, NA_character_,
    fifelse((d_opp - d_same) >= AF_PAL_MARGIN, "same",
      fifelse((d_same - d_opp) >= AF_PAL_MARGIN, "opposite", NA_character_)))
}
load_panel_af <- function(chr) {
  p <- file.path(FM, "data/ld_ref/panel_af/eur", sprintf("chr%s.af.tsv.gz", chr))
  a <- fread(p, select = c("position", "bim_a1", "bim_a2", "af_a1"), showProgress = FALSE)
  a[, `:=`(bim_a1 = toupper(bim_a1), bim_a2 = toupper(bim_a2))]
  unique(a, by = "position")
}

harmonise <- function(liv, out, panel) {
  m <- merge(liv[, .(position, a1 = allele1, a2 = allele2, lb = beta, lse = se, lp = pval, laf = af)],
             out[, .(position, ea, oa, ob = beta, ose = se, oaf = eaf, op = pval, on = n)],
             by = "position", allow.cartesian = TRUE)
  n_pre <- nrow(m)
  m[, match := ea == a1 & oa == a2]
  m[, flip := ea == a2 & oa == a1]
  m <- m[match | flip]
  m[flip == TRUE, `:=`(ob = -ob, oaf = 1 - oaf)]
  m <- m[is.finite(lb) & is.finite(lse) & lse > 0 & is.finite(ob) & is.finite(ose) & ose > 0]
  setorder(m, lp)
  m <- m[!duplicated(position)]
  m[, pal := (a1 %in% c("A", "T") & a2 %in% c("A", "T")) | (a1 %in% c("C", "G") & a2 %in% c("C", "G"))]
  m[panel, on = .(position), `:=`(pa1 = i.bim_a1, pa2 = i.bim_a2, paf = i.af_a1)]
  pf <- panel_freq_of(m$a1, m$a2, m$pa1, m$pa2, m$paf)
  m[, l_strand := strand_verdict_vs_panel(laf, pf)]
  m[, o_strand := strand_verdict_vs_panel(oaf, pf)]
  n_pal_drop <- sum(m$pal & (m$l_strand %in% "opposite" | m$o_strand %in% "opposite"))
  m <- m[!(pal & (l_strand %in% "opposite" | o_strand %in% "opposite"))]
  setorder(m, position)
  attr(m, "n_pre") <- n_pre; attr(m, "n_pal_drop") <- n_pal_drop
  m
}

abf_datasets <- function(m, lspec, tr) {
  d1 <- list(beta = m$lb, varbeta = m$lse^2, N = lspec$N, type = lspec$type, snp = as.character(m$position))
  if (lspec$type == "cc") d1$s <- lspec$Ncase / lspec$N else d1$sdY <- 1
  o <- OUTCOMES[[tr]]
  d2 <- list(beta = m$ob, varbeta = m$ose^2, N = o$N, type = o$type, snp = as.character(m$position))
  if (o$type == "cc") {
    d2$s <- o$Ncase / o$N
  } else {
    maf <- pmin(m$oaf, 1 - m$oaf)
    ok <- is.finite(maf) & maf > 0 & is.finite(m$on)
    d2$sdY <- suppressWarnings(coloc:::sdY.est(m$ose[ok]^2, maf[ok], m$on[ok]))
  }
  list(d1 = d1, d2 = d2)
}
run_abf <- function(m, lspec, tr) {
  ds <- abf_datasets(m, lspec, tr)
  r <- suppressMessages(suppressWarnings(coloc.abf(ds$d1, ds$d2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12)))
  pp <- r$summary[paste0("PP.H", 0:4, ".abf")]
  top <- r$results$snp[which.max(r$results$SNP.PP.H4)]
  list(pp = as.numeric(pp), top = as.integer(top), top_pp = max(r$results$SNP.PP.H4),
       sdY_out = if (is.null(ds$d2$sdY)) NA_real_ else ds$d2$sdY)
}

# ==============================================================================
if (stage == "coloc") {
  WID <- as.integer(args[4]); NW <- as.integer(args[5])
  Sys.setenv(MASLD_PROJECT_ROOT = REPO)
  suppressMessages(source(file.path(FM, "src/finemapping_functions.R")))
  source(file.path(FM, "src/perf/get_ld_per_locus_cached.R"))
  source(file.path(FM, "src/coloc_signal_eligibility.R"))
  tasks <- fread(file.path(WORK, "tasks.tsv"))
  chunk <- split(seq_len(nrow(tasks)), cut(seq_len(nrow(tasks)), NW, labels = FALSE))[[WID + 1L]]
  cat("worker", WID, "of", NW, ":", length(chunk), "tasks\n")
  outc <- lapply(setNames(names(OUTCOMES), names(OUTCOMES)),
                 function(tr) fread(file.path(WORK, paste0("outcome_", tr, ".tsv.gz"))))
  liv_cache <- list(); panel_chr <- NA; panel <- NULL
  dir.create(file.path(WORK, "merged"), showWarnings = FALSE)
  rows <- list(); prs <- list()
  for (k in chunk) {
    tk <- tasks[k]
    st <- tk$liver_study
    if (is.null(liv_cache[[st]])) liv_cache[[st]] <- fread(file.path(WORK, paste0("liver_", st, ".tsv.gz")))
    lv <- liv_cache[[st]][chromosome == tk$chr & position >= tk$start & position <= tk$end]
    lspec <- liver_spec(st)
    if (!identical(panel_chr, tk$chr)) { panel <- load_panel_af(tk$chr); panel_chr <- tk$chr }
    for (tr in strsplit(tk$traits, ",")[[1]]) {
      t0 <- Sys.time()
      o <- outc[[tr]][chromosome == tk$chr & position >= tk$start & position <= tk$end]
      m <- harmonise(lv, o, panel)
      row <- data.table(task_id = tk$task_id, role = tk$role, trait = tr, liver_study = st,
                        liver_trait = lspec$trait, chr = tk$chr, lead = tk$lead, genes = tk$genes,
                        n_liver = nrow(lv), n_outcome = nrow(o), n_merged_pre = attr(m, "n_pre"),
                        n_pal_dropped = attr(m, "n_pal_drop"), n_snps = nrow(m),
                        PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
                        PP.H3.abf = NA_real_, PP.H4.abf = NA_real_, abf_top = NA_integer_,
                        abf_top_pp = NA_real_, sdY_outcome = NA_real_,
                        susie_method = "not_run", PP.H4.susie = NA_real_, PP.H3.susie = NA_real_,
                        n_ld = NA_integer_, n_cs_liver = NA_integer_, n_cs_outcome = NA_integer_,
                        n_pairs = NA_integer_, n_pairs_eligible = NA_integer_,
                        susie_hit1 = NA_character_, susie_hit2 = NA_character_,
                        lambda_s_liver = NA_real_, lambda_s_outcome = NA_real_,
                        lead_liver_p = if (nrow(lv)) min(lv$pval, na.rm = TRUE) else NA_real_,
                        min_outcome_p = if (nrow(o)) min(o$pval, na.rm = TRUE) else NA_real_)
      if (nrow(m) < MIN_SNPS) {
        row$susie_method <- "too_few_snps"; rows[[length(rows) + 1L]] <- row; next
      }
      if (tk$role == "test") fwrite(m[, .(position, a1, a2, lb, lse, ob, ose, oaf, on)],
                                     file.path(WORK, "merged", paste0(tk$task_id, "_", tr, ".tsv.gz")), sep = "\t")
      ab <- tryCatch(run_abf(m, lspec, tr), error = function(e) { cat("ABF error", tk$task_id, tr, conditionMessage(e), "\n"); NULL })
      if (!is.null(ab)) {
        row[, c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf") := as.list(ab$pp)]
        row[, `:=`(abf_top = ab$top, abf_top_pp = ab$top_pp, sdY_outcome = ab$sdY_out)]
      }
      # ---- SuSiE arm (release settings) ----
      res <- tryCatch({
        ss <- data.frame(chromosome = tk$chr, position = m$position, allele1 = m$a1, allele2 = m$a2,
                         beta = m$lb, standard_error = m$lse)
        ld <- get_ld_per_locus(ss_per_locus = ss, LOCUS = tk$task_id, CHR = tk$chr,
                               START = min(m$position), END = max(m$position), ancestry = "EUR")
        if (is.null(ld) || nrow(ld[[1]]) < MIN_TRIPLE_SNPS) stop("ld_unavailable")
        sl <- as.data.table(ld[[1]])
        R <- as.matrix(ld[[2]]) + LD_REGULARIZE * diag(nrow(ld[[2]]))
        mm <- m[match(sl$position, m$position)]
        # get_ld_per_locus flips beta (and swaps alleles) to the panel orientation
        flipped <- sl$allele1 != mm$a1
        z1 <- sl$beta / sl$standard_error
        z2 <- fifelse(flipped, -mm$ob, mm$ob) / mm$ose
        ids <- paste0(tk$chr, ":", sl$position)
        names(z1) <- names(z2) <- ids; colnames(R) <- rownames(R) <- ids
        n1 <- n_eff(lspec$N, lspec$Ncase); o <- OUTCOMES[[tr]]; n2 <- n_eff(o$N, o$Ncase)
        ls1 <- tryCatch(as.numeric(susieR::estimate_s_rss(z = z1, R = R, n = n1)), error = function(e) NA_real_)
        ls2 <- tryCatch(as.numeric(susieR::estimate_s_rss(z = z2, R = R, n = n2)), error = function(e) NA_real_)
        fit <- function(z, n) { set.seed(42); tryCatch(susie_rss(z = z, R = R, n = n, L = SUSIE_L,
                     estimate_residual_variance = FALSE, residual_variance = 1, check_R = FALSE,
                     min_abs_corr = 0.5, max_iter = 500), error = function(e) NULL) }
        s1 <- fit(z1, n1); s2 <- fit(z2, n2)
        out <- list(n_ld = length(ids), ls1 = ls1, ls2 = ls2, method = "susie_not_converged",
                    pp4 = NA_real_, pp3 = NA_real_, ncs1 = NA_integer_, ncs2 = NA_integer_,
                    np = NA_integer_, npe = NA_integer_, h1 = NA_character_, h2 = NA_character_, pairs = NULL)
        if (!is.null(s1) && !is.null(s2) && isTRUE(s1$converged) && isTRUE(s2$converged)) {
          out$ncs1 <- length(s1$sets$cs); out$ncs2 <- length(s2$sets$cs)
          if (out$ncs1 == 0L || out$ncs2 == 0L) {
            out$method <- "susie_no_credible_set"
          } else {
            set.seed(42)
            sr <- tryCatch(coloc.susie(s1, s2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12), error = function(e) NULL)
            if (!is.null(sr) && !is.null(sr$summary) && nrow(sr$summary) > 0 &&
                "PP.H4.abf" %in% names(sr$summary)) {
              pt <- coloc_susie_pair_table(s1, s2, sr)
              kept <- pt$kept
              out$np <- nrow(pt$pairs); out$npe <- sum(pt$pairs$eligible %in% TRUE)
              out$pairs <- pt$pairs
              out$method <- coloc_susie_method(kept)
              if (out$method == "susie") {
                b <- kept[which.max(kept$PP.H4.abf)]
                out$pp4 <- b$PP.H4.abf; out$pp3 <- b$PP.H3.abf; out$h1 <- b$hit1; out$h2 <- b$hit2
              }
            } else out$method <- "coloc_susie_error"
          }
        }
        out
      }, error = function(e) list(method = paste0("susie_failed:", conditionMessage(e))))
      row$susie_method <- res$method
      if (!is.null(res$n_ld)) {
        row[, `:=`(n_ld = res$n_ld, lambda_s_liver = res$ls1, lambda_s_outcome = res$ls2,
                   PP.H4.susie = res$pp4, PP.H3.susie = res$pp3, n_cs_liver = res$ncs1,
                   n_cs_outcome = res$ncs2, n_pairs = res$np, n_pairs_eligible = res$npe,
                   susie_hit1 = res$h1, susie_hit2 = res$h2)]
        if (!is.null(res$pairs)) prs[[length(prs) + 1L]] <- cbind(data.table(task_id = tk$task_id, trait = tr), res$pairs)
      }
      # direction at the top shared variant (descriptive)
      tv <- if (!is.na(row$PP.H4.susie) && row$PP.H4.susie > PP_THR) as.integer(sub("^[0-9]+:", "", row$susie_hit1)) else row$abf_top
      if (!is.na(tv) && tv %in% m$position) {
        r1 <- m[position == tv]
        row[, `:=`(top_variant = tv, top_a1 = r1$a1, top_a2 = r1$a2, top_liver_beta = r1$lb,
                   top_liver_p = r1$lp, top_outcome_beta = r1$ob, top_outcome_se = r1$ose, top_outcome_p = r1$op)]
      }
      rows[[length(rows) + 1L]] <- row
      cat(sprintf("%s %s n=%d abf4=%.3f susie=%s/%.3f (%.0fs)\n", tk$task_id, tr, nrow(m),
                  row$PP.H4.abf, row$susie_method, row$PP.H4.susie, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    gc(FALSE)
  }
  fwrite(rbindlist(rows, fill = TRUE), file.path(WORK, sprintf("coloc_worker%02d.tsv", WID)), sep = "\t")
  if (length(prs)) fwrite(rbindlist(prs, fill = TRUE), file.path(WORK, sprintf("pairs_worker%02d.tsv", WID)), sep = "\t")
  cat("worker done\n")
}

# ==============================================================================
if (stage == "combine") {
  suppressMessages(library(jsonlite))
  tasks <- fread(file.path(WORK, "tasks.tsv"))
  wf <- list.files(WORK, "^coloc_worker[0-9]+\\.tsv$", full.names = TRUE)
  res <- rbindlist(lapply(wf, fread), fill = TRUE)
  expected <- sum(lengths(strsplit(tasks$traits, ",")))
  if (nrow(res) != expected) stop("result rows ", nrow(res), " != expected ", expected)
  pf <- list.files(WORK, "^pairs_worker[0-9]+\\.tsv$", full.names = TRUE)
  if (length(pf)) fwrite(rbindlist(lapply(pf, fread), fill = TRUE), file.path(RES, "susie_signal_pairs.tsv"), sep = "\t")

  res[, call_susie := !is.na(PP.H4.susie) & PP.H4.susie > PP_THR]
  res[, call_abf := !is.na(PP.H4.abf) & PP.H4.abf > PP_THR]
  res[, call_union := call_susie | call_abf]
  res[, lambda_flag := (lambda_s_liver > LAMBDA_S_HIGH) %in% TRUE | (lambda_s_outcome > LAMBDA_S_HIGH) %in% TRUE]
  # direction: liver-risk-raising allele = a1 if liver beta > 0
  res[, outcome_beta_liver_risk := fifelse(top_liver_beta > 0, top_outcome_beta, -top_outcome_beta)]
  res[, direction := fifelse(is.na(outcome_beta_liver_risk), NA_character_,
                       fifelse(trait == "eGFR",
                               fifelse(outcome_beta_liver_risk > 0, "trade_off", "concordant_harm"),
                               fifelse(outcome_beta_liver_risk < 0, "trade_off", "concordant_harm")))]

  # published-list overlap
  loci <- fread(file.path(RES, "loci_h27.tsv"))
  pub <- fread(file.path(RES, "published_lists.tsv"))
  split_genes <- function(x) unique(trimws(unlist(strsplit(x, "[;,()]"))))
  overlap <- function(lid, tr) {
    L <- loci[locus_id == lid]
    lst <- if (tr == "HF") "HF" else "kidney"
    P <- pub[list == lst]
    pos_hit <- P[!is.na(pos) & chr == L$chr & pos >= L$start - WINDOW & pos <= L$end + WINDOW]
    lg <- split_genes(L$genes)
    gene_hit <- P[vapply(genes, function(g) any(split_genes(g) %in% lg), logical(1))]
    h <- unique(rbind(pos_hit, gene_hit))
    if (nrow(h) == 0) return(list(on = FALSE, src = ""))
    list(on = TRUE, src = paste(unique(paste0(h$source, ":", h$rsid)), collapse = ";"))
  }
  tst <- res[role == "test"]
  ov <- tst[, overlap(task_id, trait), by = .(task_id, trait)]
  setnames(ov, c("on", "src"), c("on_published_list", "published_hits"))
  res <- merge(res, ov, by = c("task_id", "trait"), all.x = TRUE)
  res <- merge(res, loci[, .(task_id = locus_id, liver_cohort)], by = "task_id", all.x = TRUE)
  fwrite(res[role == "test"], file.path(RES, "locus_trait_coloc.tsv"), sep = "\t")

  # controls
  pc1 <- fread(file.path(WORK, "pc1.tsv"))
  ctl <- rbind(res[role != "test"],
               res[task_id == pc1$locus_id & trait == pc1$trait][, role := "positive_control_PC1"], fill = TRUE)
  ctl[, pass := fifelse(grepl("positive", role), call_union, !call_union)]
  fwrite(ctl, file.path(RES, "controls.tsv"), sep = "\t")

  # ---- permutation null (ABF arm), seed 20260923 ----
  set.seed(SEED_PERM)
  perm <- list()
  for (tr in names(OUTCOMES)) {
    ids <- tst[trait == tr & n_snps >= MIN_SNPS & !is.na(PP.H4.abf)]$task_id
    mf <- lapply(ids, function(i) fread(file.path(WORK, "merged", paste0(i, "_", tr, ".tsv.gz"))))
    names(mf) <- ids
    sp <- lapply(tasks[match(ids, task_id)]$liver_study, liver_spec)
    n <- length(ids)
    stopifnot(n >= 2L)
    for (b in seq_len(N_PERM)) {
      repeat { pi <- sample.int(n); if (all(pi != seq_len(n))) break }
      cnt <- 0L
      for (i in seq_len(n)) {
        a <- mf[[i]]; o <- mf[[pi[i]]]; mlen <- min(nrow(a), nrow(o))
        mm <- data.table(position = seq_len(mlen), lb = a$lb[1:mlen], lse = a$lse[1:mlen],
                         ob = o$ob[1:mlen], ose = o$ose[1:mlen], oaf = o$oaf[1:mlen], on = o$on[1:mlen])
        r <- tryCatch(run_abf(mm, sp[[i]], tr), error = function(e) NULL)
        if (!is.null(r) && r$pp[5] > PP_THR) cnt <- cnt + 1L
      }
      perm[[length(perm) + 1L]] <- data.table(trait = tr, perm = b, n_loci = n, n_calls_abf = cnt)
    }
    cat("perm", tr, "done\n")
  }
  perm <- rbindlist(perm)
  fwrite(perm, file.path(RES, "permutation_null.tsv"), sep = "\t")
  cat("combine tables done\n")
}

# ==============================================================================
# summary: reads only the tables written by combine (no recomputation)
if (stage %in% c("combine", "summary")) {
  suppressMessages(library(jsonlite))
  loci <- fread(file.path(RES, "loci_h27.tsv"))
  res  <- fread(file.path(RES, "locus_trait_coloc.tsv"))
  ctl  <- fread(file.path(RES, "controls.tsv"))
  perm <- fread(file.path(RES, "permutation_null.tsv"))
  summ <- list(label = "provisional on the current COLOC release (not adopted)",
               coloc_release = REL, threshold_PP.H4 = PP_THR,
               loci_status = as.list(table(loci$status)))
  for (tr in names(OUTCOMES)) {
    x <- res[role == "test" & trait == tr]
    tested_abf <- x[!is.na(PP.H4.abf)]
    pn <- perm[trait == tr]
    obs_abf <- sum(x$call_abf)
    summ[[tr]] <- list(
      n_loci_in_family = nrow(x), n_tested_abf = nrow(tested_abf),
      n_tested_susie_pair = nrow(x[susie_method == "susie"]),
      susie_method_counts = as.list(table(x$susie_method)),
      n_call_susie = sum(x$call_susie), n_call_abf = obs_abf, n_call_union = sum(x$call_union),
      n_call_susie_lambda_flagged = sum(x$call_susie & x$lambda_flag),
      expected_false_calls_susie = sum(1 - x[call_susie == TRUE]$PP.H4.susie),
      expected_false_calls_abf = sum(1 - x[call_abf == TRUE]$PP.H4.abf),
      perm_mean_calls_abf = mean(pn$n_calls_abf), perm_q95_calls_abf = unname(quantile(pn$n_calls_abf, 0.95)),
      perm_empirical_fdr_abf = if (obs_abf > 0) mean(pn$n_calls_abf) / obs_abf else NA,
      perm_frac_ge_observed = mean(pn$n_calls_abf >= obs_abf),
      union_calls_on_published_list = sum(x$call_union & x$on_published_list %in% TRUE),
      union_calls_new = sum(x$call_union & !(x$on_published_list %in% TRUE)),
      new_union_calls = x[call_union == TRUE & !(on_published_list %in% TRUE), paste0(task_id, "(", genes, ")")],
      published_union_calls = x[call_union == TRUE & on_published_list %in% TRUE, paste0(task_id, "(", genes, ")")],
      direction_of_union_calls = as.list(table(x[call_union == TRUE]$direction, useNA = "ifany")))
  }
  summ$controls <- lapply(split(ctl, paste(ctl$task_id, ctl$trait)), function(z)
    list(role = z$role, PP.H4.abf = z$PP.H4.abf, PP.H4.susie = z$PP.H4.susie,
         susie_method = z$susie_method, pass = z$pass))
  write_json(summ, file.path(OUT, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  writeLines(c(paste("H27", stage, format(Sys.time(), tz = "UTC", usetz = TRUE)),
               paste("permutation seed", SEED_PERM, "; SuSiE seed 42 before each fit"),
               capture.output(sessionInfo())), file.path(OUT, "sessionInfo.txt"))
  cat("summary done\n")
}
