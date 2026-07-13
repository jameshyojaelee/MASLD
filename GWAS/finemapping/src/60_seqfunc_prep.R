#!/usr/bin/env Rscript
# =============================================================================
# 60_seqfunc_prep.R  —  Stage 0: model-ready hg38 VARIANT SUBSTRATE
# =============================================================================
# Builds the single hg38 variant table that every downstream sequence-function
# stage consumes. ADDITIVE new file (deliverable 1); does NOT edit any existing
# script.
#
# Pipeline:
#   1. Load credible-set variants (hg19) + their VEP consequence/class.
#   2. Classify var_class in {coding, splice, regulatory} from Consequence/class.
#   3. Attach coloc effector gene + coloc PP.H4 (cs_member_annotation.csv).
#   4. Attach is_DEG (canonical limma-voom-qw C2 TREAT gate, treat_fdr<0.05) and
#      ensembl (GENCODE v49 metadata + gene_level_coloc).
#   5. LiftOver hg19 -> hg38 REUSING the repo pattern
#      (rtracklayer::import.chain + liftOver, chain
#       data/broadaway_eqtl/hg19ToHg38.over.chain; same as src/55, 55c, 55d).
#   6. REF-allele consistency check of the lifted position vs the hg38 FASTA
#      (Rsamtools::scanFa); flag ref_match / alt_match / mismatch; drop true
#      SNV mismatches to a rejected sidecar.
#   7. Emit results/seqfunc/variant_substrate_hg38.tsv.
#
# NOTE on is_DEG source: the task named cs_member_annotation.csv /
# lead_causal_annotation.csv for is_DEG, but those columns are unusable
# (lead_causal.is_DEG is uniformly FALSE and its deg_logFC/pp4_best/ensembl are
# empty; cs_member has no is_DEG column). is_DEG is therefore derived per
# effector gene from the CANONICAL DEG table (treat_fdr<0.05), which is the
# authoritative gene-level DEG call. Documented here + in the run log.
# =============================================================================

suppressMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(Rsamtools)
  library(Biostrings)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(BASE, "GWAS/finemapping/results/seqfunc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FA    <- "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
CHAIN <- file.path(BASE, "data/broadaway_eqtl/hg19ToHg38.over.chain")

CS_FILE     <- file.path(BASE, "GWAS/finemapping/results/credible_set_variant_consequences.csv")
MEMBER_FILE <- file.path(BASE, "RNA-seq/results/coloc_variant_classes/cs_member_annotation.csv")
DEG_FILE    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
META_FILE   <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
GLC_FILE    <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")

stopifnot(file.exists(FA), file.exists(paste0(FA, ".fai")), file.exists(CHAIN),
          file.exists(CS_FILE), file.exists(MEMBER_FILE), file.exists(DEG_FILE))

log <- function(...) cat(sprintf(...), "\n")

# -----------------------------------------------------------------------------
# 1. Credible-set variants (hg19)
# -----------------------------------------------------------------------------
cs <- fread(CS_FILE)
# expected columns: variant_id, chr, pos, ref, alt, max_pip, SYMBOL, Consequence, class, is_indel
setnames(cs, "class", "vep_class", skip_absent = TRUE)
N_IN <- nrow(cs)
log("[1] credible-set variants in (hg19): %d  (unique variant_id: %d)", N_IN, uniqueN(cs$variant_id))

# -----------------------------------------------------------------------------
# 2. var_class in {coding, splice, regulatory}
#    Priority splice > coding > regulatory. "coding" = coding-sequence variant
#    (protein-altering + synonymous, both routed to the coding/protein model);
#    the finer 4-level VEP call is preserved in var_class_detail.
# -----------------------------------------------------------------------------
cs[, var_class := fifelse(
      grepl("splice", Consequence, ignore.case = TRUE) | vep_class == "splice_region", "splice",
   fifelse(grepl("^coding", vep_class), "coding", "regulatory"))]
cs[, var_class_detail := vep_class]
log("[2] var_class: %s", paste(sprintf("%s=%d", names(table(cs$var_class)), table(cs$var_class)), collapse = "  "))

# -----------------------------------------------------------------------------
# 3. Coloc effector gene + coloc PP.H4 (cs_member_annotation.csv)
#    cs_member is pre-filtered to colocalizing CS members (all PP.H4>0.5).
#    Match on full variant_id (chr:pos:ref:alt); fall back to chr:pos.
# -----------------------------------------------------------------------------
member <- fread(MEMBER_FILE)
member[, coloc_pp4_best := as.numeric(coloc_pp4_best)]
member[, chrpos := paste0(chromosome, ":", position)]
setorder(member, -coloc_pp4_best)
m_vid <- member[!duplicated(variant_id), .(variant_id, vid_gene = coloc_genes, vid_pp4 = coloc_pp4_best)]
m_cp  <- member[!duplicated(chrpos),     .(chrpos,     cp_gene  = coloc_genes, cp_pp4  = coloc_pp4_best)]

cs[, chrpos := paste0(chr, ":", pos)]
cs <- merge(cs, m_vid, by = "variant_id", all.x = TRUE, sort = FALSE)
cs <- merge(cs, m_cp,  by = "chrpos",     all.x = TRUE, sort = FALSE)

cs[, coloc_gene     := fifelse(!is.na(vid_gene), vid_gene, cp_gene)]
cs[, coloc_pp4_best := fifelse(!is.na(vid_pp4),  vid_pp4,  cp_pp4)]
cs[, coloc_match    := fifelse(!is.na(vid_pp4), "variant_id",
                        fifelse(!is.na(cp_pp4), "chr_pos", NA_character_))]
cs[, colocalizes    := !is.na(coloc_pp4_best) & coloc_pp4_best > 0.5]
log("[3] colocalizing (PP.H4>0.5): %d  |  eQTL-absent (no coloc): %d  [match: vid=%d chrpos-only=%d]",
    sum(cs$colocalizes), sum(!cs$colocalizes),
    sum(cs$coloc_match == "variant_id", na.rm = TRUE),
    sum(cs$coloc_match == "chr_pos", na.rm = TRUE))

# gene = coloc effector when colocalizing, else the VEP-assigned gene symbol
cs[, gene := fifelse(colocalizes & !is.na(coloc_gene) & coloc_gene != "", coloc_gene, SYMBOL)]
cs[gene == "" | is.na(gene), gene := NA_character_]

# -----------------------------------------------------------------------------
# 4a. is_DEG from canonical DEG table (treat_fdr < 0.05); any effector token
# -----------------------------------------------------------------------------
deg <- fread(DEG_FILE)
deg_set <- unique(deg[!is.na(treat_fdr) & treat_fdr < 0.05, symbol])
log("[4a] canonical DEG set (treat_fdr<0.05): %d genes", length(deg_set))
first_tok <- function(x) trimws(vapply(strsplit(x, ";", fixed = TRUE), `[`, character(1), 1))
any_deg   <- function(x) vapply(strsplit(x, ";", fixed = TRUE),
                                function(g) any(trimws(g) %in% deg_set), logical(1))
cs[, is_DEG := fifelse(is.na(gene), FALSE, any_deg(gene))]

# -----------------------------------------------------------------------------
# 4b. ensembl for the primary effector gene (GENCODE v49 metadata; glc fallback)
# -----------------------------------------------------------------------------
meta <- fread(META_FILE)                       # gene_id, gene_name, chromosome, gene_biotype, ensembl_base
g2e  <- meta[!duplicated(gene_name), setNames(ensembl_base, gene_name)]
glc  <- tryCatch(fread(GLC_FILE), error = function(e) NULL)
if (!is.null(glc)) {
  glc2 <- glc[gene != "" & !is.na(ensembl) & ensembl != ""]
  g2e_glc <- glc2[!duplicated(gene), setNames(ensembl, gene)]
  add <- setdiff(names(g2e_glc), names(g2e))
  g2e <- c(g2e, g2e_glc[add])
}
cs[, primary_gene := fifelse(is.na(gene), NA_character_, first_tok(gene))]
cs[, ensembl := unname(g2e[primary_gene])]
log("[4b] ensembl mapped for %d / %d rows with a gene", sum(!is.na(cs$ensembl)), sum(!is.na(cs$gene)))

# -----------------------------------------------------------------------------
# 5. LiftOver hg19 -> hg38 (repo pattern: import.chain + liftOver, keep 1:1)
# -----------------------------------------------------------------------------
cs[, rowid := seq_len(.N)]
chain <- import.chain(CHAIN)
gr19 <- GRanges(paste0("chr", cs$chr), IRanges(cs$pos, width = 1), strand = "*")
mcols(gr19)$rowid <- cs$rowid
lifted_gr <- liftOver(gr19, chain)
nmap   <- lengths(lifted_gr)
ok_idx <- which(nmap == 1L)
gr38   <- unlist(lifted_gr[ok_idx])
lift_dt <- data.table(rowid       = mcols(gr38)$rowid,
                      chr_hg38     = as.character(seqnames(gr38)),
                      pos_hg38     = start(gr38),
                      strand_hg38  = as.character(strand(gr38)))
cs <- merge(cs, lift_dt, by = "rowid", all.x = TRUE, sort = FALSE)
cs[, lifted := !is.na(pos_hg38)]
n_lift <- sum(cs$lifted)
log("[5] liftOver hg19->hg38: %d / %d mapped 1:1 (%.1f%%); %d failed/multi-mapped",
    n_lift, N_IN, 100 * n_lift / N_IN, N_IN - n_lift)

# Restrict REF-check to contigs present in the FASTA index (drop alt/hap contigs)
fai_seqs <- fread(paste0(FA, ".fai"), header = FALSE)$V1
cs[, on_contig := lifted & (chr_hg38 %in% fai_seqs)]
n_offcontig <- sum(cs$lifted & !cs$on_contig)
if (n_offcontig > 0) log("[5] %d lifted to contigs absent from FASTA (excluded from REF-check)", n_offcontig)

# -----------------------------------------------------------------------------
# 6. REF-allele consistency check vs hg38 FASTA
# -----------------------------------------------------------------------------
COMP <- c(A = "T", T = "A", C = "G", G = "C", N = "N")
comp1 <- function(b) { b <- toupper(b); ifelse(b %in% names(COMP), COMP[b], "N") }

chk <- cs[on_contig == TRUE]
gr_q <- GRanges(chk$chr_hg38, IRanges(chk$pos_hg38, width = 1), strand = "+")
hg38_base <- toupper(as.character(scanFa(FaFile(FA), gr_q)))
chk[, hg38_base := hg38_base]

# strand-aware first-base of ref/alt (minus-strand lift => compare complement)
chk[, ref1 := toupper(substr(ref, 1, 1))]
chk[, alt1 := toupper(substr(alt, 1, 1))]
chk[strand_hg38 == "-", ref1 := unname(comp1(ref1))]
chk[strand_hg38 == "-", alt1 := unname(comp1(alt1))]
chk[, ref_match := hg38_base == ref1]
chk[, alt_match := hg38_base == alt1]
chk[, ref_status := fifelse(ref_match, "ref_match",
                     fifelse(alt_match, "alt_match", "mismatch"))]

# hg38-reference-oriented alleles (+ strand). SNVs only: ref_hg38 = the hg38
# reference base; alt_hg38 = the other file allele. Indels are left NA (need
# explicit VCF normalization downstream) but flagged via is_indel.
is_snv <- (nchar(chk$ref) == 1L & nchar(chk$alt) == 1L)
chk[, ref_hg38 := NA_character_]
chk[, alt_hg38 := NA_character_]
chk[is_snv & ref_status == "ref_match", `:=`(ref_hg38 = ref1, alt_hg38 = alt1)]
chk[is_snv & ref_status == "alt_match", `:=`(ref_hg38 = alt1, alt_hg38 = ref1)]

cs <- merge(cs, chk[, .(rowid, hg38_base, ref_match, ref_status, ref_hg38, alt_hg38)],
            by = "rowid", all.x = TRUE, sort = FALSE)
cs[lifted == TRUE & on_contig == FALSE, ref_status := "off_contig"]
cs[lifted == FALSE,                     ref_status := "liftover_failed"]
cs[is.na(ref_match), ref_match := FALSE]

# is_indel present? indels never dropped for REF reasons (multi-base ref, lenient first-base check)
if (!"is_indel" %in% names(cs)) cs[, is_indel := (nchar(ref) != 1 | nchar(alt) != 1)]

pct <- function(n, d) sprintf("%.1f%%", 100 * n / max(d, 1))
n_chk <- nrow(chk)
log("[6] REF-check (on-contig lifted, n=%d): ref_match=%d (%s)  alt_match=%d (%s)  mismatch=%d (%s)",
    n_chk,
    sum(chk$ref_status == "ref_match"), pct(sum(chk$ref_status == "ref_match"), n_chk),
    sum(chk$ref_status == "alt_match"), pct(sum(chk$ref_status == "alt_match"), n_chk),
    sum(chk$ref_status == "mismatch"),  pct(sum(chk$ref_status == "mismatch"),  n_chk))

# -----------------------------------------------------------------------------
# 7. Split substrate vs rejected and write
#    Keep: lifted, on-contig, and (ref_match | alt_match | is_indel).
#    Reject: liftover_failed, off_contig, or true SNV mismatch.
# -----------------------------------------------------------------------------
cs[, keep := lifted & on_contig &
             (ref_status %in% c("ref_match", "alt_match") | is_indel == TRUE)]

out_cols <- c("variant_id_hg19", "chr", "pos_hg38", "ref", "alt", "gene", "ensembl",
              "max_pip", "var_class", "coloc_pp4_best", "is_DEG", "colocalizes")
cs[, variant_id_hg19 := variant_id]

sub <- cs[keep == TRUE]
# assemble in required schema order, then provenance extras
setnames(sub, "chr", "chr_hg19_num", skip_absent = TRUE)  # avoid clobber; output 'chr' = hg38 contig
sub[, chr := chr_hg38]
out <- sub[, .(variant_id_hg19, chr, pos_hg38, ref, alt, gene, ensembl, max_pip,
               var_class, coloc_pp4_best, is_DEG, colocalizes,
               # provenance / extras (additive, downstream may ignore)
               pos_hg19 = variant_id_hg19,      # placeholder, replaced below
               var_class_detail, consequence = Consequence, primary_gene,
               coloc_match, hg38_base, ref_hg38, alt_hg38, ref_status, strand_hg38, is_indel)]
out[, pos_hg19 := as.integer(tstrsplit(variant_id_hg19, ":", fixed = TRUE)[[2]])]

setorder(out, chr, pos_hg38)
OUT_TSV <- file.path(OUT_DIR, "variant_substrate_hg38.tsv")
fwrite(out, OUT_TSV, sep = "\t")
log("[7] wrote substrate: %s  (%d variants, %d columns)", OUT_TSV, nrow(out), ncol(out))

# rejected sidecar (audit)
rej <- cs[keep == FALSE, .(variant_id_hg19 = variant_id, chr_hg19 = chr, pos_hg19 = pos,
                           ref, alt, chr_hg38, pos_hg38, hg38_base, ref_status,
                           var_class, gene, colocalizes)]
REJ_TSV <- file.path(OUT_DIR, "variant_substrate_hg38.rejected.tsv")
fwrite(rej, REJ_TSV, sep = "\t")
log("[7] wrote rejected sidecar: %s  (%d variants)", REJ_TSV, nrow(rej))

# -----------------------------------------------------------------------------
# 8. Summary + PNPLA3 rs738409 spot-check
# -----------------------------------------------------------------------------
log("\n================ SUMMARY ================")
log("total variants in            : %d", N_IN)
log("liftover 1:1 success         : %d (%s)", n_lift, pct(n_lift, N_IN))
log("REF-match (on-contig lifted) : %d (%s of %d)", sum(chk$ref_status=="ref_match"),
    pct(sum(chk$ref_status=="ref_match"), n_chk), n_chk)
log("  (ref|alt)-consistent        : %d (%s)", sum(chk$ref_status %in% c("ref_match","alt_match")),
    pct(sum(chk$ref_status %in% c("ref_match","alt_match")), n_chk))
log("substrate written            : %d", nrow(out))
log("rejected                     : %d", nrow(rej))
log("--- substrate var_class ---")
print(table(out$var_class))
log("--- substrate coloc ---")
log("colocalizing : %d   eQTL-absent : %d", sum(out$colocalizes), sum(!out$colocalizes))
log("is_DEG TRUE   : %d", sum(out$is_DEG))

log("\n--- PNPLA3 rs738409 spot-check (hg19 22:44324727 C>G, missense) ---")
pn <- out[variant_id_hg19 == "22:44324727:C:G"]
if (nrow(pn) == 1) {
  log("lifted -> %s:%d  | var_class=%s (detail=%s) | ref=%s alt=%s hg38_base=%s ref_status=%s",
      pn$chr, pn$pos_hg38, pn$var_class, pn$var_class_detail, pn$ref, pn$alt, pn$hg38_base, pn$ref_status)
  log("gene=%s ensembl=%s max_pip=%.4f colocalizes=%s coloc_pp4=%s is_DEG=%s",
      pn$gene, pn$ensembl, pn$max_pip, pn$colocalizes,
      ifelse(is.na(pn$coloc_pp4_best), "NA", sprintf("%.3f", pn$coloc_pp4_best)), pn$is_DEG)
  ok <- pn$chr == "chr22" && pn$pos_hg38 == 43928847L && pn$var_class == "coding"
  log("SPOT-CHECK %s (expected chr22:43928847, coding)", ifelse(ok, "PASS", "FAIL"))
} else {
  log("SPOT-CHECK FAIL: PNPLA3 22:44324727:C:G not found in substrate (rows=%d)", nrow(pn))
}
log("=========================================")
