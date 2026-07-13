#!/usr/bin/env Rscript
# Select LD-aware multi-variant and saturation targets. This does NOT label
# correlated variants as observed haplotypes because phased genotypes are absent.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
ROOT <- Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SF <- file.path(ROOT,"GWAS/finemapping/results/seqfunc")
OUT <- file.path(SF,"haplotype_saturation"); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
fm <- fread(file.path(ROOT,"GWAS/finemapping/results/combined_finemapping.csv"),
            select=c("chromosome","position","allele1","allele2","trait","locus",
                     "study","ancestry","susie_pip_clean","max_pip","either_in_cs",
                     "variant_id","recommended_pip"))
fm[, pip := fcoalesce(suppressWarnings(as.numeric(recommended_pip)),
                      suppressWarnings(as.numeric(max_pip)),
                      suppressWarnings(as.numeric(susie_pip_clean)),0)]
# Use the version-controlled portfolio authority, not a phenotype-name regex.
tier <- fread(file.path(ROOT,"GWAS/finemapping/config/gwas_trait_tier.tsv"))
primary_studies <- sort(tier[placement=="main",study_name])
stopifnot(length(primary_studies) == 35L)
fm[, study_canonical := sub("_coloc_targeted$","",study,ignore.case=TRUE)]
fm <- fm[toupper(study_canonical) %in% toupper(primary_studies)]
available_primary_studies <- sort(unique(fm$study_canonical))
missing_primary_studies <- primary_studies[!toupper(primary_studies) %in%
                                            toupper(available_primary_studies)]
# The current combined variant substrate lacks the three MVP-EUR fine-mapping
# arms. Record that coverage loss explicitly; never fill it with non-primary studies.
stopifnot(length(available_primary_studies) == 32L,
          setequal(toupper(missing_primary_studies),
                   c("MVP_ALT_EUR","MVP_AST_EUR","MVP_NAFLD_EUR")))
fm <- fm[pip >= 0.05]
fm[, variant_key_hg19 := paste(chromosome,position,pmin(allele1,allele2),
                               pmax(allele1,allele2),sep=":")]
fm <- fm[pip >= 0.05 & (either_in_cs == TRUE | pip >= 0.1)]
fm[, locus_key := paste(study,locus,sep="|")]
loci <- fm[, .(n_variants=.N, pip_sum=sum(pip), lead_pip=max(pip),
               chromosome=first(chromosome), ancestry=first(ancestry),
               span_bp=max(position)-min(position),
               variants=paste(variant_id[order(-pip)],collapse=";")), by=locus_key]
loci[, `:=`(design_class=fifelse(n_variants>=2,"credible_set_multi_variant_candidate","single_variant"),
            upstream_finemapping_used_ancestry_ld=TRUE,phase_used=FALSE,
            joint_haplotype_effect_estimable=FALSE)]
fwrite(loci[order(-n_variants,-pip_sum)],file.path(OUT,"locus_design_inventory.tsv"),sep="\t")
sub <- fread(file.path(SF,"variant_substrate_hg38.tsv")); sub[, max_pip:=as.numeric(max_pip)]
# Recompute the saturation PIP from the 35 paper-primary studies.  The substrate
# max_pip spans the full portfolio and must not determine this experiment.
primary_pip <- fm[, .(primary_max_pip=max(pip,na.rm=TRUE),
                      primary_n_studies=uniqueN(study),
                      primary_studies=paste(sort(unique(study)),collapse=";")),
                  by=variant_key_hg19]
sub[, variant_key_hg19 := paste(sub("^chr","",chr),pos_hg19,pmin(ref,alt),
                                pmax(ref,alt),sep=":")]
sub <- primary_pip[sub,on="variant_key_hg19",nomatch=0]
mp <- file.path(SF,"mpra_benchmark/mpra_nomination_audit.tsv")
if(file.exists(mp)){m<-fread(mp); dav_ids<-unique(m$variant_id_hg19)} else dav_ids<-character()
sub[, official_hu_dav := variant_id_hg19 %in% dav_ids]
# Avoid evaluation leakage: Hu MPRA truth and eQTL status remain control/
# stratification annotations and never determine which windows are selected.
sub[, saturation_priority := primary_max_pip]
sat <- sub[var_class=="regulatory" & !is_indel & primary_max_pip >= 0.05]
sat <- sat[order(-saturation_priority,variant_key_hg19)]
sat[, window_bp := 501L]
sat[, design := "all_SNV_substitutions_plus_observed_credible_alleles"]
fwrite(sat[1:min(.N,500)],file.path(OUT,"saturation_candidate_manifest.tsv"),sep="\t")
write_json(list(status="candidate_design_complete",apply_only_firewall=TRUE,
  n_primary_studies_authority=length(primary_studies), primary_studies=primary_studies,
  n_primary_studies_available=length(available_primary_studies),
  missing_primary_studies=missing_primary_studies,
  haplotype_status="BLOCKED_no_phased_reference_genotypes",
  multi_variant_status="CANDIDATE_ONLY_upstream_finemapping_used_LD_but_no_phase_or_joint_sequence_model",
  saturation_selection="top primary-trait PIP; MPRA truth excluded from ranking",
  saturation_interpretation="mechanistic_nomination_not_causal_evidence",
  saturation_status="500 regulatory SNV windows selected for in-silico mechanistic nomination",
  no_alphagenome_training=TRUE),file.path(OUT,"contract.json"),pretty=TRUE,auto_unbox=TRUE)
message("[99] wrote ",nrow(loci)," locus designs and ",min(nrow(sat),500)," saturation targets")
