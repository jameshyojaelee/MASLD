#!/usr/bin/env Rscript
# Machine-readable validation companion to the detailed methods audit.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
ROOT<-Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SF<-file.path(ROOT,"GWAS/finemapping/results/seqfunc")
OUT<-file.path(SF,"methods_audit");dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

classification<-data.table(
 component=c("Script 93 Nasser ABC","Script 93 MASLD SCENIC+","Script 95 coding prediction",
   "Script 95 bulk disease proteomics","Script 95 external liver pQTL colocalization",
   "Script 100 substrate annotation overlaps","Script 107 thresholded-PIP lasso",
   "Script 108 bespoke-prior SuSiE pilot","Canonical uniform-prior SuSiE"),
 classification=c("SENSITIVITY","SENSITIVITY","PRIMARY_WITHIN_CODING_ARM","CONTEXT_ONLY","SENSITIVITY_CONTEXT_ONLY",
   "DESCRIPTIVE_ONLY","RETIRED","RETIRED_PROVENANCE","PRIMARY"),
 legitimate_use=c("candidate enhancer-gene prediction","candidate disease-context eGRN link",
   "protein consequence hypothesis","gene-level disease abundance context","gene-level pQTL context; not same-variant validation",
   "overlap description within selected substrate","none for inference","reproduce two pilot loci only","canonical statistical fine-mapping"),
 reason=c("published ABC prediction, not direct link measurement","inferred region-to-gene relation, not physical contact",
   "sequence/protein model appropriate to coding alleles","disease abundance cannot establish germline allele mechanism",
   "gene-level join does not connect nominated coding allele to protein level",
   "PIP-enriched substrate is not a genome-wide negative universe",
   "thresholded posterior used as truth; phenotype/locus duplication; arbitrary balancing/cap",
   "inherits invalid prior weights despite technically valid susieR prior_weights API",
   "independent of seqfunc annotations and protected by firewall"))
fwrite(classification,file.path(OUT,"component_classification.tsv"),sep="\t")

links<-fread(file.path(SF,"activity_by_contact/measured_variant_gene_links.tsv"))
coding<-fread(file.path(SF,"coding_protein_bridge/coding_variant_protein_evidence.tsv"))
ann<-fread(file.path(SF,"functional_prior_sensitivity/measured_annotation_matrix.tsv"))
checks<-data.table(
 check=c("all links explicitly unvalidated","ABC and SCENIC scores not summed","coding rows explicitly not same-variant validated",
   "retired prior guards present","canonical consumers do not read retired weights","selected annotation substrate size"),
 pass=c(all(links$effector_gene_validated==FALSE),
   !any(c("combined_link_score","summed_score") %in% names(links)),
   all(coding$same_variant_protein_mechanism_validated==FALSE),
   all(vapply(c("107_functional_prior_crossfit.R","108_functional_prior_susie.R"),function(f)
     any(grepl("ALLOW_RETIRED_BESPOKE_PRIOR",readLines(file.path(ROOT,"GWAS/finemapping/src",f),warn=FALSE))),logical(1))),
   TRUE,nrow(ann)==6579L),
 detail=c(sprintf("%d predicted links",nrow(links)),"source-specific link_strength retained",
   sprintf("%d coding variants",nrow(coding)),"explicit opt-in required","validated by source scan below",
   sprintf("n=%d; not genome-wide",nrow(ann))))

# Static source scan: the canonical atlas/convergence scripts must not consume
# retired prior artifacts or sensitivity PIPs.
consumers<-c("RNA-seq/27a_assemble_evidence_atlas.R","RNA-seq/46d_convergence_evidence.R",
             "RNA-seq/75_integrate_causal_overhaul.R","RNA-seq/217_stratified_causal_atlas.R")
consumers<-consumers[file.exists(file.path(ROOT,consumers))]
bad<-character()
for(f in consumers){z<-readLines(file.path(ROOT,f),warn=FALSE);if(any(grepl("crossfit_prior|functional_prior_pip|susie_rerun",z)))bad<-c(bad,f)}
checks[check=="canonical consumers do not read retired weights",`:=`(pass=length(bad)==0,
  detail=if(length(bad)) paste(bad,collapse=",") else paste("PASS; scanned",length(consumers),"consumers"))]
fwrite(checks,file.path(OUT,"validation_checks.tsv"),sep="\t")
write_json(list(status=if(all(checks$pass))"PASS" else "FAIL",canonical_uniform_prior_authority=TRUE,
 retired=c("107 thresholded-PIP lasso","108 bespoke-prior pilot inference"),
 n_predicted_variant_gene_links=nrow(links),n_coding_context_rows=nrow(coding),
 n_selected_substrate_annotations=nrow(ann),checks=checks),file.path(OUT,"validation_summary.json"),
 pretty=TRUE,auto_unbox=TRUE)
if(!all(checks$pass))stop("methods-audit validation failed")
message("[109] methods audit validation PASS")
