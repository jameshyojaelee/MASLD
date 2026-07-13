#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
V2 <- file.path(ROOT,"GWAS/finemapping/results/seqfunc/mpra_benchmark/v2")
dir.create(V2,recursive=TRUE,showWarnings=FALSE)

cells <- c("HepG2","LX2")
contexts <- list(HepG2=c("control","PAOA"),LX2=c("control","TGFb"))
source_rows <- list()
source_complete <- TRUE
for (cell in cells) {
  p <- file.path(V2,"source_reproduction",paste0(cell,".source_reproduction_qc.json"))
  if (!file.exists(p)) {
    source_complete <- FALSE
    next
  }
  x <- as.data.table(fromJSON(p,simplifyDataFrame=TRUE))
  x[,qc_path:=p]
  source_rows[[cell]] <- x
}
source_qc <- rbindlist(source_rows,fill=TRUE)
for (nm in c("cell_line","context","n_tested","n_reproduced","n_official","official_coverage",
             "direction_concordance","official_positive_recall","pass_coverage","pass_direction",
             "pass_call_recall","effect_direction")) {
  if (!nm %in% names(source_qc)) source_qc[, (nm) := NA]
}
source_pass <- source_complete && nrow(source_qc)==4L && all(
  source_qc$pass_coverage & source_qc$pass_direction & source_qc$pass_call_recall
)

env_path <- file.path(V2,"environment_qc.json")
environment_complete <- file.exists(env_path) && isTRUE(fromJSON(env_path)$environment_complete)

check_mpr <- function(mode) {
  subdir <- if (mode=="smoke") "mpranalyze_smoke" else "mpranalyze"
  out <- file.path(V2,subdir)
  rows <- list()
  complete <- TRUE
  for (cell in cells) {
    contract_path <- file.path(out,paste0(cell,".contract.json"))
    depth_path <- file.path(out,paste0(cell,".depth_factor_qc.tsv"))
    result_paths <- file.path(out,paste0(cell,".",c("alt_vs_ref_control","alt_vs_ref_stimulus","allele_by_stimulus"),".tsv.gz"))
    ok <- file.exists(contract_path) && file.exists(depth_path) && all(file.exists(result_paths))
    contract <- if (file.exists(contract_path)) fromJSON(contract_path) else list()
    depth <- if (file.exists(depth_path)) fread(depth_path) else data.table()
    depth_ok <- nrow(depth)>0L && all(is.finite(depth$dna_depth) & depth$dna_depth>0 & is.finite(depth$rna_depth) & depth$rna_depth>0 & depth$dna_unique==1L & depth$rna_unique==1L)
    result_nonempty <- ok && all(file.info(result_paths)$size>0)
    cell_ok <- ok && isTRUE(contract$completed) && identical(contract$mode,mode) && depth_ok && result_nonempty
    complete <- complete && cell_ok
    rows[[cell]] <- data.table(
      cell_line=cell, mode=mode, complete=cell_ok,
      n_elements=if(is.null(contract$n_elements)) NA_integer_ else as.integer(contract$n_elements),
      max_barcodes=if(is.null(contract$max_barcodes_per_allele)) NA_integer_ else as.integer(contract$max_barcodes_per_allele),
      n_depth_libraries=nrow(depth), depth_factors_pass=depth_ok,
      n_result_files=sum(file.exists(result_paths))
    )
  }
  list(complete=complete,qc=rbindlist(rows))
}

smoke <- check_mpr("smoke")
full <- check_mpr("full")
full_submission_allowed <- source_pass && environment_complete && smoke$complete

verdict <- list(
  generated_at_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
  source_reproduction_complete=source_complete,
  source_reproduction_pass=source_pass,
  source_reproduction_qc=if(nrow(source_qc)) source_qc[,.(cell_line,context,n_tested,n_reproduced,n_official,official_coverage,direction_concordance,official_positive_recall,pass_coverage,pass_direction,pass_call_recall,effect_direction)] else list(),
  mpranalyze_environment_complete=environment_complete,
  mpranalyze_smoke_complete=smoke$complete,
  mpranalyze_smoke_qc=smoke$qc,
  mpranalyze_full_submission_allowed=full_submission_allowed,
  mpranalyze_completion=full$complete,
  mpranalyze_full_qc=full$qc,
  official_source_calls_authoritative=TRUE,
  mpranalyze_source_authoritative=FALSE,
  canonical_outputs_mutated=FALSE
)
path <- file.path(V2,"gate_verdict.json")
tmp <- paste0(path,".tmp.",Sys.getpid())
write_json(verdict,tmp,pretty=TRUE,auto_unbox=TRUE,null="null",digits=NA)
if (!file.rename(tmp,path)) stop("Could not atomically write ",path)
message(toJSON(verdict,pretty=TRUE,auto_unbox=TRUE,null="null"))
