#!/usr/bin/env Rscript
# Source-faithful Hu MPRA reconstruction. Official DAV calls remain authoritative.
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(readxl)
  library(jsonlite)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CELL <- Sys.getenv("MPRA_CELL_LINE", "")
if (!CELL %in% c("HepG2", "LX2")) stop("MPRA_CELL_LINE must be HepG2 or LX2")
SF <- file.path(ROOT, "GWAS/finemapping/results/seqfunc/mpra_benchmark/v2")
COUNT_DIR <- file.path(SF, "counts")
OUT <- file.path(SF, "source_reproduction")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

counts_path <- file.path(COUNT_DIR, paste0(CELL, ".barcode_counts.tsv.gz"))
manifest_path <- file.path(COUNT_DIR, paste0(CELL, ".sample_manifest.tsv"))
stopifnot(file.exists(counts_path), file.exists(manifest_path))
message("Reading ", counts_path)
d <- fread(counts_path, nThread=max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))))
m <- fread(manifest_path)
stopifnot(nrow(m) == 8L, all(c("construct","element_id","allele","barcode") %in% names(d)))

# Source methods exclude a barcode in a replicate when its DNA count is zero,
# then sum retained barcodes by oligo. Preserve that rule independently per sample.
agg <- rbindlist(lapply(seq_len(nrow(m)), function(i) {
  sid <- m$sample_id[i]
  dc <- paste0("dna__", sid); rc <- paste0("rna__", sid)
  stopifnot(dc %in% names(d), rc %in% names(d))
  d[get(dc) > 0, .(
    DNA=sum(as.numeric(get(dc))), RNA=sum(as.numeric(get(rc))),
    n_barcodes=.N
  ), by=.(element_id, allele)][, `:=`(
    sample_id=sid, condition=m$condition[i], replicate=as.integer(m$replicate[i])
  )]
}))
fwrite(agg, file.path(OUT, paste0(CELL, ".allele_replicate_counts.tsv.gz")), sep="\t")

run_context <- function(ctx) {
  a <- agg[condition == ctx]
  if (uniqueN(a$replicate) != 4L) stop("Expected four replicates for ", CELL, " ", ctx)
  constructs <- unique(a[, .(element_id, allele)])
  constructs[, key := paste(element_id, allele, sep="__")]
  sample_rows <- CJ(replicate=1:4, assay=c("DNA","RNA"))
  sample_rows[, sample := paste0(assay, "_r", replicate)]
  mat <- matrix(0L, nrow=nrow(constructs), ncol=nrow(sample_rows),
                dimnames=list(constructs$key, sample_rows$sample))
  for (j in seq_len(nrow(sample_rows))) {
    z <- a[replicate == sample_rows$replicate[j], .(key=paste(element_id,allele,sep="__"),
          value=as.numeric(get(sample_rows$assay[j])))]
    mat[z$key, j] <- round(z$value)
  }
  coldata <- data.frame(replicate=factor(sample_rows$replicate), assay=factor(sample_rows$assay, levels=c("DNA","RNA")),
                        row.names=sample_rows$sample)
  dds <- DESeqDataSetFromMatrix(mat, coldata, design=~replicate + assay)
  dds <- DESeq(dds, quiet=TRUE, parallel=FALSE)
  active <- as.data.table(as.data.frame(results(dds, contrast=c("assay","RNA","DNA"))), keep.rownames="key")
  active <- merge(active, constructs, by="key", all.x=TRUE)
  active[, active_oligo := !is.na(padj) & padj < 0.05 & log2FoldChange > 0]
  fwrite(active, file.path(OUT, paste0(CELL, ".", ctx, ".active_oligos.tsv.gz")), sep="\t")

  active_variant <- active[, .(active_variant=any(active_oligo)), by=element_id]
  # Retain every assayed ref/alt pair in the comparison table. The active-element
  # rule gates DAV calls, but must not redefine published-variant coverage.
  wide <- dcast(a, element_id + allele ~ replicate,
                value.var=c("DNA","RNA","n_barcodes"), fill=0)
  ref <- wide[allele == "ref"]; alt <- wide[allele == "alt"]
  setnames(ref, setdiff(names(ref), c("element_id","allele")), paste0(setdiff(names(ref), c("element_id","allele")), "_ref"))
  setnames(alt, setdiff(names(alt), c("element_id","allele")), paste0(setdiff(names(alt), c("element_id","allele")), "_alt"))
  paired <- merge(ref, alt, by="element_id", all=FALSE)
  dav <- paired[, {
    # Hu Table S2-S5 reports alternative-over-reference activity.
    dna_ratio <- vapply(1:4, function(r) get(paste0("DNA_",r,"_alt"))/get(paste0("DNA_",r,"_ref")), numeric(1))
    rna_ratio <- vapply(1:4, function(r) get(paste0("RNA_",r,"_alt"))/get(paste0("RNA_",r,"_ref")), numeric(1))
    dna_control <- median(dna_ratio[is.finite(dna_ratio) & dna_ratio > 0], na.rm=TRUE)
    norm <- rna_ratio/dna_control
    l2 <- log2(norm[is.finite(norm) & norm > 0])
    tt <- if (length(l2) >= 2L && sd(l2) > 0) t.test(l2, mu=0) else NULL
    .(log2FC_alt_over_ref=if(length(l2)) mean(l2) else NA_real_,
      p_value=if(is.null(tt)) NA_real_ else tt$p.value,
      n_ratio_replicates=length(l2), median_DNA_alt_over_ref=dna_control,
      min_barcodes_ref=min(vapply(1:4, function(r) get(paste0("n_barcodes_",r,"_ref")), numeric(1))),
      min_barcodes_alt=min(vapply(1:4, function(r) get(paste0("n_barcodes_",r,"_alt")), numeric(1))))
  }, by=element_id]
  dav <- merge(dav, active_variant, by="element_id", all.x=TRUE)
  dav[is.na(active_variant), active_variant := FALSE]
  dav[, fdr := NA_real_]
  dav[active_variant == TRUE, fdr := p.adjust(p_value, method="BH")]
  dav[, reproduced_dav := active_variant == TRUE & !is.na(fdr) & fdr < 0.01]
  fwrite(dav, file.path(OUT, paste0(CELL, ".", ctx, ".source_DAV_reproduction.tsv.gz")), sep="\t")

  table_no <- if (CELL == "HepG2" && ctx == "control") 2L else if (CELL == "LX2" && ctx == "control") 3L else if (CELL == "HepG2") 4L else 5L
  off <- as.data.table(read_excel(file.path(ROOT,"GWAS/finemapping/data/seqfunc_external/hu2025_mpra/raw",paste0("TableS",table_no,".xlsx")), skip=2))
  setnames(off, names(off)[1:5], c("element_id","official_log2FC","official_p","official_fdr","rsid"))
  off <- off[!is.na(element_id)]
  cmp <- merge(off, dav, by="element_id", all.x=TRUE)
  cmp[, direction_match := sign(as.numeric(official_log2FC)) == sign(log2FC_alt_over_ref)]
  cmp[, reproduced_positive := reproduced_dav %in% TRUE]
  fwrite(cmp, file.path(OUT, paste0(CELL, ".",ctx,".official_comparison.tsv.gz")), sep="\t")
  coverage <- mean(!is.na(cmp$log2FC_alt_over_ref))
  direction <- mean(cmp$direction_match,na.rm=TRUE)
  recall <- mean(cmp$reproduced_positive,na.rm=TRUE)
  list(
    cell_line=CELL, context=ctx, n_tested=uniqueN(dav$element_id), n_reproduced=sum(dav$reproduced_dav,na.rm=TRUE),
    n_official=nrow(off), official_coverage=coverage,
    direction_concordance=direction,
    official_positive_recall=recall,
    coverage_threshold=0.98, direction_threshold=0.99, call_recall_threshold=0.95,
    pass_coverage=coverage >= 0.98,
    pass_direction=direction >= 0.99,
    pass_call_recall=recall >= 0.95,
    effect_direction="alternative_over_reference",
    official_calls_authoritative=TRUE
  )
}

contexts <- if (CELL == "HepG2") c("control","PAOA") else c("control","TGFb")
qc <- lapply(contexts, run_context)
if (!all(vapply(qc, function(x) x$pass_coverage && x$pass_direction && x$pass_call_recall, logical(1)))) {
  warning("One or more source-reproduction gates failed; see QC JSON")
}
write_json(qc, file.path(OUT, paste0(CELL,".source_reproduction_qc.json")), pretty=TRUE, auto_unbox=TRUE)
message(toJSON(qc, pretty=TRUE, auto_unbox=TRUE))
