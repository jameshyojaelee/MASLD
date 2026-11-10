#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
source('GWAS/finemapping/src/coloc_variant_identity.R')
stopifnot(identical(coloc_variant_identity_collisions(c('1:1','1:2')),integer()),
          identical(coloc_variant_identity_collisions(c('1:1','1:2','1:1')),c(1L,3L)))
stopifnot(coloc_reference_ld_valid(diag(2)),
          !coloc_reference_ld_valid(matrix(c(1,2,2,1),2)),
          !coloc_reference_ld_valid(matrix(c(1,NA,NA,1),2)),
          !coloc_reference_ld_valid(matrix(c(1,.2,.4,1),2)))
tmp <- tempfile('eqtl-identity-');dir.create(tmp)
# Execute the injected eQTL control flow with the actual MPRIP allele collision.
# A unique-ID control must reach the fit branch unchanged; MPRIP must not.
fixture <- c(
  'summary_file <- file.path(tmp,"summary.tsv")',
  'fwrite(data.table(gene=character(),ensg=character(),chr=integer(),n_snps=integer(),converged=logical(),n_cs=integer(),max_pip=numeric(),status=character()),summary_file,sep="\t")',
  'n_failed <- 0L; reached_fit <- character(); chr_out_dir <- tmp; CHR <- 17L',
  'for (kind in c("MPRIP","invalid_LD","unique_control")) {',
  '  gene_name <- ensg_id <- kind',
  '  ss_filtered <- data.frame(position=if(kind=="MPRIP") c(17164840L,17164840L) else c(1L,2L),allele1=c("CTAAA","C"),allele2=c("C","CTAAA"))',
  '  n_snps <- nrow(ss_filtered)',
  '  R_mat <- if(kind=="invalid_LD") matrix(c(1,2,2,1),2) else diag(2); LD_REGULARIZE <- .001',
  '  snp_ids <- paste0(CHR, ":", ss_filtered$position)',
  '  reached_fit <- c(reached_fit,kind)',
  '}'
)
eval(parse(text=coloc_identity_guarded_eqtl_runner(fixture)))
stopifnot(identical(reached_fit,'unique_control'),n_failed==2L,
          identical(fread(summary_file)$status,c('ambiguous_variant_identity','invalid_reference_LD')),
          nrow(fread(file.path(tmp,'MPRIP_identity_collisions.tsv')))==2L)
invisible(parse(text=coloc_identity_guarded_eqtl_runner(readLines('GWAS/finemapping/src/08_run_eqtl_susie.R'))))
invisible(parse(text=coloc_identity_guarded_runner(readLines('GWAS/finemapping/src/06_susie_coloc.R'))))
cat('PASS: duplicate identities and invalid LD refused before fitting; valid unique control preserved; both transformed producers parse\n')
