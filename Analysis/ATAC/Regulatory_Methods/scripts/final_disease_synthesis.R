#!/usr/bin/env Rscript
# Same-donor/peak disease sensitivities; candidate scientific context only.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, nzchar(Sys.getenv('SLURM_JOB_ID')))
script_root <- normalizePath(args[1L])
output <- args[2L]
stopifnot(!file.exists(output))
dir.create(output, recursive = TRUE)
source(file.path(script_root, 'scripts/compare_da.R'))
data.table::setDTthreads(1L)
cfg <- fromJSON(file.path(script_root, 'config/analysis.json'), simplifyVector = FALSE)
baseline <- file.path(script_root, 'results/20261023/da_full/full')
existing <- fread(file.path(baseline, 'peak_results.tsv.gz'))
existing <- existing[representation == 'fragment' & filter == 'primary' & method == 'voom_qw']
results <- list(); models <- list(); availability <- list(); overlaps <- list()
influence <- list(); stability <- list(); roster <- list(); checks <- list()
for (lineage in c('hepatocyte', 'stellate', 'macrophage')) {
  inputs <- lapply(names(cfg$cohorts), function(cohort) {
    x <- load_da_counts(cfg$counts_root, cohort, lineage)
    pairing <- fread(cfg$cohorts[[cohort]]$metadata)
    for (column in intersect(c('batch_mm','batch','sex','age','BMI','bmi'), names(pairing))) {
      x$donors[[column]] <- pairing[[column]][match(x$donors$donor_id, pairing$donor_id)]
    }
    path <- cfg$cohort_covariate_metadata[[cohort]]
    x$donors <- join_donor_metadata(x$donors, path)
    # Expression-derived stage is deliberately excluded from analytical metadata.
    x$donors <- x$donors[, intersect(names(x$donors), c('donor_id','condition','n_cells','contrast_eligible',
      'batch_mm','batch','sex','age','BMI','bmi')), with = FALSE]
    x
  })
  names(inputs) <- names(cfg$cohorts)
  keep <- lapply(inputs, function(x) x$peaks$peak_coordinate[filter_da_peaks(x$counts, 'primary')])
  common <- Reduce(intersect, keep)
  for (cohort in names(inputs)) {
    x <- inputs[[cohort]]; d <- x$donors
    roster[[length(roster)+1L]] <- cbind(data.table(cohort, lineage), d)
    if (!length(common)) {
      models[[length(models)+1L]] <- data.table(cohort, lineage, adjustment='condition_only',
        status='untestable', reason='empty_joint_primary_peak_universe', n_peaks=0L)
      next
    }
    counts <- x$counts[, common, drop=FALSE]
    base <- fit_da(counts, d, 'voom_qw')
    expected_index <- which(existing$cohort == cohort & existing$lineage == lineage)
    expected <- existing[expected_index]
    anchor <- merge(base$results, expected, by='peak_coordinate', suffixes=c('_new','_saved'))
    stopifnot(nrow(anchor) == length(common), nrow(expected) == length(common))
    delta <- max(abs(anchor$logFC_new-anchor$logFC_saved), abs(anchor$pvalue_new-anchor$pvalue_saved),
      abs(anchor$standard_error_new-anchor$standard_error_saved))
    stopifnot(delta < 1e-7)
    checks[[length(checks)+1L]] <- data.table(cohort,lineage,status='passed',n_peaks=length(common),maximum_difference=delta)
    add_fit <- function(fit, adjustment, n) {
      value <- copy(fit$results)
      value[, `:=`(cohort=cohort,lineage=lineage,adjustment=adjustment,n_donors=n,
        effect_unit='log2_normalized_fragment_count',bh_family='all_fixed_joint_primary_peaks_within_cohort_lineage_model',
        evidence_status='candidate_not_adopted')]
      results[[length(results)+1L]] <<- value
      q <- fit$qc[1L, .(design_rank,n_model_columns,residual_df,n_tested_peaks)]
      models[[length(models)+1L]] <<- cbind(data.table(cohort,lineage,adjustment,status='testable',reason='',n_donors=n),q)
    }
    add_fit(base,'condition_only',nrow(d))
    for (column in c('batch_mm','batch','sex','age','BMI','bmi')) {
      if (!column %in% names(d)) {
        availability[[length(availability)+1L]] <- data.table(cohort,lineage,covariate=column,status='not_recorded',n_complete=0L)
        next
      }
      ok <- !is.na(d[[column]]) & nzchar(trimws(as.character(d[[column]])))
      dc <- copy(d[ok]); cc <- counts[ok,,drop=FALSE]
      if (grepl('batch|sex',column)) dc[[column]] <- as.character(dc[[column]])
      if (!grepl('batch|sex',column)) dc[[column]] <- suppressWarnings(as.numeric(dc[[column]]))
      if (grepl('batch',column) && nrow(dc)) {
        cross <- as.data.table(table(batch=as.character(dc[[column]]),condition=dc$condition))
        cross[, `:=`(cohort=cohort,lineage=lineage,covariate=column)]
        overlaps[[length(overlaps)+1L]] <- cross
      }
      eligible <- nrow(dc) > 0 && !anyNA(dc[[column]]) && all(table(factor(dc$condition,levels=c('NORMAL','MASH'))) >= 4L)
      design <- if (eligible) tryCatch(make_da_design(dc,column),error=function(e)e) else simpleError('fewer_than_four_complete_donors_per_condition_or_invalid_value')
      status <- if (inherits(design,'error')) conditionMessage(design) else 'estimable'
      availability[[length(availability)+1L]] <- data.table(cohort,lineage,covariate=column,status=status,n_complete=nrow(dc))
      if (inherits(design,'error')) next
      if (nrow(dc) != nrow(d)) add_fit(fit_da(cc,dc,'voom_qw'),paste0(column,'_matched_condition_only'),nrow(dc))
      add_fit(fit_da(cc,dc,'voom_qw',column),column,nrow(dc))
    }
    delta_min <- rep(Inf,length(common)); delta_max <- rep(-Inf,length(common))
    sign_stable <- rep(TRUE,length(common)); max_delta <- rep(0,length(common))
    for (donor in d$donor_id) {
      idx <- d$donor_id != donor
      reduced <- fit_da(counts[idx,,drop=FALSE],d[idx],'voom_qw')$results
      stopifnot(identical(reduced$peak_coordinate,base$results$peak_coordinate))
      change <- reduced$logFC-base$results$logFC
      delta_min <- pmin(delta_min,reduced$logFC); delta_max <- pmax(delta_max,reduced$logFC)
      sign_stable <- sign_stable & reduced$logFC*base$results$logFC > 0
      max_delta <- pmax(max_delta,abs(change))
      influence[[length(influence)+1L]] <- data.table(cohort,lineage,omitted_donor=donor,
        median_absolute_effect_change=median(abs(change)),maximum_absolute_effect_change=max(abs(change)),
        sign_agreement=mean(reduced$logFC*base$results$logFC > 0),n_peaks=length(common))
    }
    stability[[length(stability)+1L]] <- data.table(cohort,lineage,peak_coordinate=common,
      lodo_min_logFC=delta_min,lodo_max_logFC=delta_max,lodo_all_same_sign=sign_stable,lodo_max_absolute_change=max_delta)
    message('Completed ',cohort,' / ',lineage)
  }
}
for (name in c('results','models','availability','overlaps','influence','stability','roster','checks')) {
  value <- rbindlist(get(name),fill=TRUE)
  fwrite(value,file.path(output,paste0(name,if(name %in% c('results','stability')) '.tsv.gz' else '.tsv')),sep='\t')
}
writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))
write_json(list(status='completed',candidate_only=TRUE,models=length(results),
  no_new_primary_multiplicity_family=TRUE),file.path(output,'completion.json'),auto_unbox=TRUE,pretty=TRUE)
