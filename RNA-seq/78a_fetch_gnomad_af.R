#!/usr/bin/env Rscript
# ============================================================================
# 78a_fetch_gnomad_af.R
# One-time network fetcher for Script 78 [2] gnomAD-AF component.
#
# Takes MASLD-GWAS credible-set variants from combined_finemapping.csv, lifts
# hg19 -> hg38 (combined_finemapping is hg19), assigns each variant to its
# nearest-TSS gene (GENCODE v49, the same convention as Script 55), then pulls
# ancestry-stratified allele frequencies (NFE/EAS/global) from the gnomAD v4
# GraphQL API and aggregates per gene.
#
# Restriction: queries the higher-confidence credible-set core (recommended_pip
# > 0.1; ~1,600 unique hg38 variants) to keep the per-variant API call count
# tractable and rate-limit-safe. This is the scientifically defensible subset
# (the fine-mapped signal), and is recorded in the output provenance.
#
# Output (cache read by Script 78): RNA-seq/results/multi_evidence/gnomad_af_credset.tsv
#   one row per gene: gnomad_af_eur, gnomad_af_eas, gnomad_af_global,
#                     gnomad_ancestry_divergence, gnomad_af_n_variants, gnomad_af_lead_variant
# Lightweight (network-bound); safe on the login node.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(httr); library(jsonlite)
  library(rtracklayer); library(GenomicRanges)
})

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM    <- file.path(BASE, "GWAS/finemapping/results/combined_finemapping.csv")
CHAIN <- file.path(BASE, "data/broadaway_eqtl/hg19ToHg38.over.chain")
GTF   <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
OUT   <- file.path(BASE, "RNA-seq/results/multi_evidence/gnomad_af_credset.tsv")
PIP_MIN <- as.numeric(Sys.getenv("GNOMAD_PIP_MIN", "0.1"))

stopifnot(file.exists(FM), file.exists(CHAIN), file.exists(GTF))
cat("== 78a_fetch_gnomad_af ==\n")

# ---- 1. credible-set core variants ----------------------------------------
fm <- fread(FM)
cs <- fm[(either_in_cs %in% TRUE | recommended_pip > PIP_MIN) & recommended_pip > PIP_MIN]
cs <- unique(cs[, .(chromosome, position, allele1, allele2, recommended_pip)])
cat(sprintf("  credible-set core (recommended_pip > %.2f): %d unique hg19 variants\n",
            PIP_MIN, nrow(cs)))

# ---- 2. liftOver hg19 -> hg38 ----------------------------------------------
chain <- import.chain(CHAIN)
gr19 <- GRanges(paste0("chr", cs$chromosome), IRanges(cs$position, width = 1))
mcols(gr19)$idx <- seq_len(nrow(cs))
lo <- liftOver(gr19, chain)
ok <- lengths(lo) == 1L
gr38 <- unlist(lo[ok])
cs38 <- cs[mcols(gr38)$idx]
cs38[, `:=`(chr_hg38 = as.character(seqnames(gr38)), pos_hg38 = start(gr38))]
std <- paste0("chr", c(1:22, "X", "Y"))
cs38 <- cs38[chr_hg38 %in% std]
cat(sprintf("  liftOver hg19->hg38: %d / %d on standard chromosomes\n", nrow(cs38), nrow(cs)))

# ---- 3. nearest-TSS gene (GENCODE v49; same convention as Script 55) -------
genes <- rtracklayer::import(GTF, feature.type = "gene")
genes <- genes[genes$gene_type %in% c("protein_coding", "lncRNA")]
tss   <- resize(genes, width = 1L, fix = "start")
gv    <- GRanges(cs38$chr_hg38, IRanges(cs38$pos_hg38, width = 1))
ni    <- nearest(gv, tss, ignore.strand = TRUE)
cs38[, gene := genes$gene_name[ni]]
cs38 <- cs38[!is.na(gene)]
# gnomAD variant id (chrom-pos-ref-alt, no "chr" prefix); allele1 = ref, allele2 = alt
cs38[, vid := sprintf("%s-%d-%s-%s", gsub("chr", "", chr_hg38), pos_hg38, allele1, allele2)]
cs38 <- unique(cs38, by = "vid")
cat(sprintf("  nearest-TSS gene assigned: %d variants -> %d genes\n",
            nrow(cs38), uniqueN(cs38$gene)))

# ---- 4. gnomAD v4 GraphQL: batched per-variant AF --------------------------
API  <- "https://gnomad.broadinstitute.org/api"
frag <- "genome{ac an populations{id ac an}} exome{ac an populations{id ac an}}"

query_batch <- function(vids) {
  aliases <- vapply(seq_along(vids), function(i)
    sprintf('v%d: variant(variantId:"%s",dataset:gnomad_r4){variant_id %s}',
            i, vids[i], frag), character(1))
  q <- paste0("query{", paste(aliases, collapse = " "), "}")
  for (attempt in 1:4) {
    r <- tryCatch(
      POST(API, add_headers("Content-Type" = "application/json"),
           body = toJSON(list(query = q), auto_unbox = TRUE), timeout(90)),
      error = function(e) NULL)
    if (!is.null(r) && status_code(r) == 200) {
      ct <- content(r, as = "parsed", type = "application/json")
      if (!is.null(ct$data)) return(ct$data)
    }
    Sys.sleep(2 * attempt)   # backoff on transient error / rate-limit
  }
  NULL
}

# pop AF helper: gnomad pop ids -> nfe (EUR), eas (EAS); combine genome+exome AC/AN
pop_af <- function(v, pid) {
  ac <- 0; an <- 0
  for (src in c("genome", "exome")) {
    s <- v[[src]]; if (is.null(s)) next
    for (p in s$populations) if (!is.null(p$id) && p$id == pid) { ac <- ac + p$ac; an <- an + p$an }
  }
  if (an > 0) ac / an else NA_real_
}
global_af <- function(v) {
  ac <- 0; an <- 0
  for (src in c("genome", "exome")) {
    s <- v[[src]]; if (is.null(s)) next
    if (!is.null(s$ac) && !is.null(s$an)) { ac <- ac + s$ac; an <- an + s$an }
  }
  if (an > 0) ac / an else NA_real_
}

vids  <- cs38$vid
BATCH <- 12L
recs  <- vector("list", length(vids)); names(recs) <- vids
n_batches <- ceiling(length(vids) / BATCH)
cat(sprintf("  querying gnomAD v4 GraphQL: %d variants in %d batches of %d\n",
            length(vids), n_batches, BATCH))
for (b in seq_len(n_batches)) {
  idx  <- ((b - 1) * BATCH + 1):min(b * BATCH, length(vids))
  vb   <- vids[idx]
  data <- query_batch(vb)
  if (!is.null(data)) {
    for (i in seq_along(vb)) {
      v <- data[[sprintf("v%d", i)]]
      if (is.null(v)) next
      recs[[vb[i]]] <- data.table(
        vid    = vb[i],
        af_eur = pop_af(v, "nfe"),
        af_eas = pop_af(v, "eas"),
        af_glb = global_af(v))
    }
  }
  if (b %% 20 == 0) cat(sprintf("    ...batch %d/%d (%d variants found so far)\n",
                                b, n_batches, sum(!vapply(recs, is.null, logical(1)))))
  Sys.sleep(0.4)   # be polite to the public API
}

af <- rbindlist(recs, fill = TRUE)
cat(sprintf("  gnomAD returned AF for %d / %d variants\n", nrow(af), length(vids)))
if (nrow(af) == 0) stop("gnomAD API returned nothing — aborting (cache not written)")

# ---- 5. aggregate per gene -------------------------------------------------
vg <- merge(cs38[, .(vid, gene, recommended_pip)], af, by = "vid")
setorder(vg, gene, -recommended_pip)
gnomad <- vg[, .(
  gnomad_af_eur          = max(af_eur, na.rm = TRUE),
  gnomad_af_eas          = max(af_eas, na.rm = TRUE),
  gnomad_af_global       = max(af_glb, na.rm = TRUE),
  gnomad_af_n_variants   = .N,
  gnomad_af_lead_variant = vid[which.max(recommended_pip)]
), by = .(human_symbol = gene)]
# ancestry divergence: EUR vs EAS differ >2x (and both observed > 0)
gnomad[, gnomad_ancestry_divergence :=
         is.finite(gnomad_af_eur) & is.finite(gnomad_af_eas) &
         gnomad_af_eur > 0 & gnomad_af_eas > 0 &
         (pmax(gnomad_af_eur, gnomad_af_eas) / pmax(pmin(gnomad_af_eur, gnomad_af_eas), 1e-9) > 2)]
for (c in c("gnomad_af_eur","gnomad_af_eas","gnomad_af_global"))
  gnomad[!is.finite(get(c)), (c) := NA_real_]

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(gnomad, OUT, sep = "\t")
cat(sprintf("WROTE %s : %d genes\n", OUT, nrow(gnomad)))
for (g in c("PNPLA3","TM6SF2","HSD17B13","GCKR","APOB","LIPA","ATP7B","SAMM50","PARVB"))
  if (g %in% gnomad$human_symbol) {
    r <- gnomad[human_symbol == g]
    cat(sprintf("  [%s] eur=%.4g eas=%.4g global=%.4g divergent=%s n=%d (%s)\n",
                g, r$gnomad_af_eur, r$gnomad_af_eas, r$gnomad_af_global,
                r$gnomad_ancestry_divergence, r$gnomad_af_n_variants, r$gnomad_af_lead_variant))
  }
