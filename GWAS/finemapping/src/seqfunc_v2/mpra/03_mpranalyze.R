#!/usr/bin/env Rscript
# MPRAnalyze sensitivity for Hu MPRA. Source calls remain authoritative.
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(MPRAnalyze)
  library(BiocParallel)
  library(jsonlite)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CELL <- Sys.getenv("MPRA_CELL_LINE", "")
if (!CELL %in% c("HepG2", "LX2")) stop("MPRA_CELL_LINE must be HepG2 or LX2")
NCPU <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "32")))
MODE <- Sys.getenv("MPRA_MODE", "full")
if (!MODE %in% c("smoke", "full")) stop("MPRA_MODE must be smoke or full")
MAX_BARCODES <- if (MODE == "smoke") 10L else 100L
SF <- file.path(ROOT, "GWAS/finemapping/results/seqfunc/mpra_benchmark/v2")
COUNT_DIR <- file.path(SF, "counts")
OUT <- file.path(SF, if (MODE == "smoke") "mpranalyze_smoke" else "mpranalyze")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

counts_path <- file.path(COUNT_DIR, paste0(CELL,".barcode_counts.tsv.gz"))
# The isolated data.table build delegates direct gzip reads to the optional
# R.utils package. Stream through system gzip instead of adding an otherwise
# unnecessary R dependency to this analysis environment.
d <- fread(cmd=paste("gzip -dc", shQuote(counts_path)), nThread=NCPU)
m <- fread(file.path(COUNT_DIR, paste0(CELL,".sample_manifest.tsv")))
stopifnot(nrow(m)==8L, all(c("ref","alt") %in% d$allele))
dna_cols <- paste0("dna__",m$sample_id); rna_cols <- paste0("rna__",m$sample_id)
stopifnot(all(c(dna_cols,rna_cols) %in% names(d)))

# Frozen support rule: at least ten barcodes per allele and nonzero aggregate
# DNA in at least three of four replicates for each allele/context.
support <- d[,.(n_barcodes=uniqueN(barcode)),by=.(element_id,allele)]
support_w <- dcast(support,element_id~allele,value.var="n_barcodes",fill=0)
eligible <- support_w[ref>=10 & alt>=10,element_id]
for (i in seq_len(nrow(m))) {
  dc <- dna_cols[i]
  z <- d[element_id %in% eligible,.(DNA=sum(as.numeric(get(dc)))),by=.(element_id,allele)]
  z[, nonzero := DNA>0]
  z[, `:=`(condition=m$condition[i], replicate=as.integer(m$replicate[i]))]
  if (i==1L) dna_support <- z else dna_support <- rbind(dna_support,z)
}
dna_ok <- dna_support[,.(n_nonzero=sum(nonzero)),by=.(element_id,allele,condition)]
dna_ok <- dna_ok[,.(all_contexts_pass=uniqueN(condition)==2L & all(n_nonzero>=3L)),by=.(element_id,allele)]
dna_ok <- dna_ok[,.(both_alleles_pass=uniqueN(allele)==2L & all(all_contexts_pass)),by=element_id]
eligible <- intersect(eligible,dna_ok[both_alleles_pass==TRUE,element_id])
if (length(eligible)<100L) stop("Fewer than 100 elements satisfy frozen MPRAnalyze support rules")
eligible <- sort(eligible)
if (MODE == "smoke") eligible <- head(eligible, 100L)
d <- d[element_id %in% eligible]

# MPRAnalyze recommends limiting very large barcode designs. Rank within each
# allele by total DNA support, retain the deterministic top 100, and preserve the
# same barcode slots across control/stimulus samples.
d[, total_DNA := rowSums(.SD), .SDcols=dna_cols]
setorder(d,element_id,allele,-total_DNA,barcode)
d[, barcode_slot:=seq_len(.N),by=.(element_id,allele)]
d <- d[barcode_slot<=MAX_BARCODES]
elements <- sort(unique(d$element_id))
element_index <- setNames(seq_along(elements),elements)

make_block <- function(sid, allele_value) {
  z <- d[d$allele == allele_value]
  i <- unname(element_index[z$element_id]); j <- z$barcode_slot
  dna <- matrix(0,nrow=length(elements),ncol=MAX_BARCODES)
  rna <- matrix(0,nrow=length(elements),ncol=MAX_BARCODES)
  dna[cbind(i,j)] <- as.numeric(z[[paste0("dna__",sid)]])
  rna[cbind(i,j)] <- as.numeric(z[[paste0("rna__",sid)]])
  list(dna=dna,rna=rna)
}

dna_blocks <- list(); rna_blocks <- list(); annot <- list(); k <- 0L
for (s in seq_len(nrow(m))) for (al in c("ref","alt")) {
  k <- k+1L; b <- make_block(m$sample_id[s],al)
  dna_blocks[[k]] <- b$dna; rna_blocks[[k]] <- b$rna
  annot[[k]] <- data.table(
    sample_id=m$sample_id[s], cell_line=CELL,
    stimulus=m$condition[s], replicate=as.integer(m$replicate[s]),
    allele=al, barcode_slot=seq_len(MAX_BARCODES)
  )
}
dna <- do.call(cbind,dna_blocks); rna <- do.call(cbind,rna_blocks)
rownames(dna)<-rownames(rna)<-elements
ann <- rbindlist(annot)
ann[, `:=`(
  sample_id=factor(sample_id),
  stimulus=factor(stimulus,levels=c("control",setdiff(unique(stimulus),"control"))),
  replicate=factor(replicate), allele=factor(allele,levels=c("ref","alt")),
  barcode_slot=factor(barcode_slot)
)]
ann[, barcode_allelic:=interaction(allele,barcode_slot,drop=TRUE)]
stopifnot(ncol(dna)==nrow(ann),identical(dim(dna),dim(rna)))

bp <- MulticoreParam(workers=NCPU, progressbar=TRUE)
depth_qc <- list()
fit_contrast <- function(label, keep, rna_design, reduced_design, dna_design) {
  subset_annot <- droplevels(as.data.frame(ann[keep]))
  obj <- MpraObject(
    dnaCounts=dna[,keep,drop=FALSE], rnaCounts=rna[,keep,drop=FALSE],
    dnaAnnot=subset_annot, rnaAnnot=subset_annot
  )
  # DNA and RNA are distinct sequencing libraries even when they share a sample
  # label. Estimate their depths separately, grouping the barcode columns by the
  # exact physical library ID. This avoids accidentally recycling one assay's
  # depth factors into the other and is robust to subsetted contrast designs.
  obj <- estimateDepthFactors(obj,lib.factor="sample_id",which.lib="dna")
  obj <- estimateDepthFactors(obj,lib.factor="sample_id",which.lib="rna")
  dd <- dnaDepth(obj); rd <- rnaDepth(obj)
  if (length(dd) != sum(keep) || length(rd) != sum(keep) ||
      any(!is.finite(dd) | dd <= 0) || any(!is.finite(rd) | rd <= 0)) {
    stop("Invalid MPRAnalyze DNA/RNA depth factors for ", label)
  }
  # Use the exact annotation frame supplied to MpraObject. Besides keeping the
  # QC aligned with the fitted columns, this guarantees that contrast-specific
  # subsets cannot reintroduce unused sample levels after depth estimation.
  aq <- as.data.table(subset_annot)[, .(sample_id)]
  aq[, `:=`(dna_depth=as.numeric(dd),rna_depth=as.numeric(rd))]
  within_sample <- aq[, .(
    dna_unique=uniqueN(signif(dna_depth,12)), rna_unique=uniqueN(signif(rna_depth,12)),
    dna_depth=first(dna_depth), rna_depth=first(rna_depth)
  ),by=sample_id]
  if (any(within_sample$dna_unique != 1L) || any(within_sample$rna_unique != 1L)) {
    stop("Depth factors are not constant within physical library for ", label)
  }
  depth_qc[[label]] <<- within_sample[, `:=`(contrast=label,mode=MODE)]
  obj <- analyzeComparative(
    obj=obj,dnaDesign=dna_design,rnaDesign=rna_design,reducedDesign=reduced_design,
    fit.se=TRUE,mode="classic",BPPARAM=bp
  )
  ans <- as.data.table(testLrt(obj),keep.rownames="element_id")
  if (!"fdr" %in% names(ans) && "pval" %in% names(ans)) ans[,fdr:=p.adjust(pval,"BH")]
  ans[,`:=`(cell_line=CELL,contrast=label,max_barcodes_per_allele=MAX_BARCODES,
            support_rule="n_barcodes>=10; aggregate_DNA_nonzero_in>=3/4_replicates_each_context")]
  out_path <- file.path(OUT,paste0(CELL,".",label,".tsv.gz"))
  tmp <- paste0(out_path,".tmp.",Sys.getpid())
  fwrite(ans,tmp,sep="\t",compress="gzip")
  if (!file.rename(tmp,out_path)) stop("Atomic output rename failed: ",out_path)
  ans
}

control <- ann$stimulus=="control"
stimulus_name <- setdiff(levels(ann$stimulus),"control")
stim <- ann$stimulus==stimulus_name
res_control <- fit_contrast("alt_vs_ref_control",control,
  ~replicate+allele,~replicate,~barcode_allelic+replicate)
res_stim <- fit_contrast("alt_vs_ref_stimulus",stim,
  ~replicate+allele,~replicate,~barcode_allelic+replicate)
res_interaction <- fit_contrast("allele_by_stimulus",rep(TRUE,nrow(ann)),
  ~replicate+stimulus*allele,~replicate+stimulus+allele,
  ~barcode_allelic+replicate+stimulus)

contract <- list(
  cell_line=CELL,stimulus=stimulus_name,mode=MODE,completed=TRUE,n_elements=length(elements),
  max_barcodes_per_allele=MAX_BARCODES,n_observations=ncol(dna),
  support_rule="at least 10 barcodes per allele and aggregate DNA nonzero in at least 3/4 replicates for each allele/context",
  barcode_selection=paste0("top ",MAX_BARCODES," per element/allele by total DNA across all eight samples; deterministic barcode tie-break"),
  effect_direction="alternative_over_reference; positive logFC means alternative allele has higher activity",
  depth_factors="DNA and RNA estimated separately by exact physical sample_id; positive finite and constant within library",
  dna_model="barcode nested within allele + replicate (+ stimulus for interaction)",
  rna_models=c(control="replicate + allele",stimulus="replicate + allele",interaction="replicate + stimulus * allele"),
  multiple_testing="BH within cell-line/contrast family",
  official_source_calls_authoritative=TRUE,MPRAnalyze_sensitivity_only=TRUE,
  MPRAnalyze_source_authoritative=FALSE
)
fwrite(rbindlist(depth_qc),file.path(OUT,paste0(CELL,".depth_factor_qc.tsv")),sep="\t")
write_json(contract,file.path(OUT,paste0(CELL,".contract.json")),pretty=TRUE,auto_unbox=TRUE)
message(toJSON(contract,pretty=TRUE,auto_unbox=TRUE))
