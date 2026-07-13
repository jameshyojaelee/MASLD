#!/usr/bin/env Rscript
# Strict post-scoring QC and per-anchor summary for adult-liver ChromBPNet
# saturation mutagenesis. Apply-only: never updates fine-mapping or convergence.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BASE <- file.path(ROOT,"GWAS/finemapping/results/seqfunc/haplotype_saturation")
SEQ <- Sys.getenv("SATURATION_SEQUENCE_DIR",file.path(BASE,"sequences"))
SCORES <- Sys.getenv("SATURATION_SCORE_DIR",file.path(SEQ,"scores"))
OUT <- Sys.getenv("SATURATION_AGG_DIR",SCORES)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

stopf <- function(...) stop(sprintf(...),call.=FALSE)
expected_n <- 671187L; target_map_n <- 751500L; expected_shards <- sprintf("%03d",0:26)
expected_files <- file.path(SEQ,paste0("genomic_variants.shard",expected_shards,".tsv"))
score_files <- file.path(SCORES,paste0("genomic_shard",expected_shards,".variant_scores.tsv"))
missing_expected <- expected_files[!file.exists(expected_files) | file.info(expected_files)$size==0]
missing_scores <- score_files[!file.exists(score_files) | file.info(score_files)$size==0]
if(length(missing_expected)) stopf("missing expected-manifest shards: %s",paste(basename(missing_expected),collapse=","))
if(length(missing_scores)) stopf("missing score shards: %s",paste(basename(missing_scores),collapse=","))
extra <- list.files(SCORES,pattern="^genomic_shard[0-9]+\\.variant_scores\\.tsv$",full.names=TRUE)
extra <- setdiff(normalizePath(extra),normalizePath(score_files))
if(length(extra)) stopf("unexpected score shards: %s",paste(basename(extra),collapse=","))

required_expected <- c("chrom","pos_hg38","ref","alt","variant_id")
required_score <- c("chr","pos","allele1","allele2","variant_id","logfc","abs_logfc",
  "jsd","abs_logfc_x_jsd","active_allele_quantile",
  "abs_logfc_x_jsd_x_active_allele_quantile","abs_logfc.pval","jsd.pval",
  "abs_logfc_x_jsd.pval","abs_logfc_x_jsd_x_active_allele_quantile.pval")

expected_list <- vector("list",27); score_list <- vector("list",27); shard_qc <- vector("list",27)
for(i in seq_along(expected_files)) {
  e <- fread(expected_files[i]); s <- fread(score_files[i])
  me <- setdiff(required_expected,names(e)); ms <- setdiff(required_score,names(s))
  if(length(me)) stopf("%s missing expected columns: %s",basename(expected_files[i]),paste(me,collapse=","))
  if(length(ms)) stopf("%s missing score columns: %s",basename(score_files[i]),paste(ms,collapse=","))
  if(anyDuplicated(e$variant_id)) stopf("duplicate expected variant_id in shard %s",expected_shards[i])
  if(anyDuplicated(s$variant_id)) stopf("duplicate scored variant_id in shard %s",expected_shards[i])
  if(nrow(e)!=nrow(s)) stopf("row mismatch shard %s: expected=%d scored=%d",expected_shards[i],nrow(e),nrow(s))
  # Set equality first; ordering is deliberately not assumed.
  if(!setequal(e$variant_id,s$variant_id)) stopf("variant set mismatch in shard %s",expected_shards[i])
  chk <- e[s,on="variant_id"]
  bad <- chk[chrom!=chr | pos_hg38!=pos | ref!=allele1 | alt!=allele2]
  if(nrow(bad)) stopf("coordinate/allele mismatch in shard %s (%d rows)",expected_shards[i],nrow(bad))
  prov_file <- file.path(SCORES,paste0("genomic_shard",expected_shards[i],".provenance.json"))
  if(!file.exists(prov_file) || file.info(prov_file)$size==0) stopf("missing provenance for shard %s",expected_shards[i])
  prov <- fromJSON(prov_file)
  if(!isTRUE(prov$apply_only_firewall) || prov$interpretation!="mechanistic_nomination_not_causal_evidence")
    stopf("invalid provenance/firewall for shard %s",expected_shards[i])
  if(any(!unique(e$chrom) %in% prov$heldout_test_chromosomes))
    stopf("shard %s contains chromosomes not held out by its model",expected_shards[i])
  e[, shard:=expected_shards[i]]; s[, shard:=expected_shards[i]]
  expected_list[[i]] <- e; score_list[[i]] <- s
  shard_qc[[i]] <- data.table(shard=expected_shards[i],n_expected=nrow(e),n_scored=nrow(s),
    n_unique_expected=uniqueN(e$variant_id),n_unique_scored=uniqueN(s$variant_id),
    coordinate_allele_mismatches=0L,status="PASS")
}
expected <- rbindlist(expected_list); scored <- rbindlist(score_list,fill=TRUE)
if(nrow(expected)!=expected_n || uniqueN(expected$variant_id)!=expected_n)
  stopf("global expected-set failure: rows=%d unique=%d required=%d",nrow(expected),uniqueN(expected$variant_id),expected_n)
if(nrow(scored)!=expected_n || uniqueN(scored$variant_id)!=expected_n)
  stopf("global score-set failure: rows=%d unique=%d required=%d",nrow(scored),uniqueN(scored$variant_id),expected_n)
if(!setequal(expected$variant_id,scored$variant_id)) stopf("global expected/scored set mismatch")

keep_score <- unique(c(required_score,setdiff(names(scored),names(expected))))
genomic_scores <- expected[scored[,..keep_score],on="variant_id"]
setnames(genomic_scores,"variant_id","genomic_variant_id")
# Empirical shuffled-sequence p-values are scorer-standard. Correct each metric
# family across unique genomic SNVs before expansion; target-wise BH is added
# after expansion for local map annotation. These are not causal p-values.
pmap <- c(abs_logfc="abs_logfc.pval",jsd="jsd.pval",ies="abs_logfc_x_jsd.pval",
          ips="abs_logfc_x_jsd_x_active_allele_quantile.pval")
for(nm in names(pmap)) {
  pcol <- pmap[[nm]]; genomic_scores[,(pcol):=as.numeric(get(pcol))]
  if(any(!is.finite(genomic_scores[[pcol]]) | genomic_scores[[pcol]]<0 | genomic_scores[[pcol]]>1))
    stopf("invalid empirical p-values in %s",pcol)
  genomic_scores[,paste0("q_global_",nm):=p.adjust(get(pcol),method="BH")]
}
target_map_file <- file.path(SEQ,"target_variant_map.tsv.gz")
if(!file.exists(target_map_file) || file.info(target_map_file)$size==0) stopf("missing target_variant_map.tsv.gz")
target_map <- fread(target_map_file)
required_map <- c("target_perturbation_id","genomic_variant_id","target_id","target_rank",
  "offset_from_anchor","anchor_pos_hg38","primary_max_pip","primary_n_studies",
  "source_variant_id_hg19","is_credible_allele","apply_only","alphagenome_training_prohibited")
mm <- setdiff(required_map,names(target_map)); if(length(mm)) stopf("target map missing columns: %s",paste(mm,collapse=","))
if(nrow(target_map)!=target_map_n || uniqueN(target_map$target_perturbation_id)!=target_map_n)
  stopf("target-map identity failure: rows=%d unique=%d",nrow(target_map),uniqueN(target_map$target_perturbation_id))
if(uniqueN(target_map$genomic_variant_id)!=expected_n) stopf("target map must reference exactly %d genomic SNVs",expected_n)
if(any(toupper(as.character(target_map$apply_only))!="TRUE") ||
   any(toupper(as.character(target_map$alphagenome_training_prohibited))!="TRUE")) stopf("target-map firewall violation")
joined <- genomic_scores[target_map,on="genomic_variant_id"]
if(nrow(joined)!=target_map_n || any(!is.finite(joined$abs_logfc))) stopf("genomic-score expansion failed")
joined[, variant_id:=target_perturbation_id]
joined[, `:=`(ies=abs_logfc_x_jsd,
  ips=abs_logfc_x_jsd_x_active_allele_quantile,
  apply_only=TRUE, alphagenome_training_prohibited=TRUE)]
numeric_required <- c("logfc","abs_logfc","jsd","ies","ips")
for(x in numeric_required) {
  joined[,(x):=as.numeric(get(x))]
  if(any(!is.finite(joined[[x]]))) stopf("non-finite values in required score %s",x)
}
for(nm in names(pmap)) joined[,paste0("q_target_",nm):=p.adjust(get(pmap[[nm]]),method="BH"),by=target_id]

rank_desc <- function(x) frank(-x,ties.method="min",na.last="keep")
joined[, `:=`(rank_abs_logfc=rank_desc(abs_logfc),rank_jsd=rank_desc(jsd),
              rank_ies=rank_desc(ies),rank_ips=rank_desc(ips)),by=target_id]
if(uniqueN(joined$target_id)!=500L) stopf("expected 500 targets; saw %d",uniqueN(joined$target_id))
target_sizes <- joined[, .N,by=target_id]
if(any(target_sizes$N!=1503L)) stopf("target size failure: %d targets do not contain 1503 SNVs",sum(target_sizes$N!=1503L))
anchors <- joined[toupper(as.character(is_credible_allele))=="TRUE"]
if(nrow(anchors)!=500L || any(anchors[, .N,by=target_id]$N!=1L)) stopf("credible-anchor uniqueness failure")

toprow <- function(d,metric,prefix,qname) {
  z <- d[which.max(get(metric))]
  z <- z[,.(variant_id,offset_from_anchor,value=get(metric),
            q_global=get(paste0("q_global_",qname)),q_target=get(paste0("q_target_",qname)))]
  setnames(z,c("variant_id","offset_from_anchor","value","q_global","q_target"),
           paste0(prefix,c("_variant","_offset","_value","_q_global","_q_target")))
  z
}
summary <- joined[, c(list(n_saturation_snvs=.N),
  toprow(.SD,"abs_logfc","top_abs_logfc","abs_logfc"),toprow(.SD,"jsd","top_jsd","jsd"),
  toprow(.SD,"ies","top_ies","ies"),toprow(.SD,"ips","top_ips","ips")),by=.(target_id,target_rank)]
anchor_keep <- anchors[,.(target_id,anchor_variant_id=variant_id,anchor_offset=offset_from_anchor,
  anchor_logfc=logfc,anchor_abs_logfc=abs_logfc,anchor_jsd=jsd,anchor_ies=ies,anchor_ips=ips,
  anchor_rank_abs_logfc=rank_abs_logfc,anchor_rank_jsd=rank_jsd,
  anchor_rank_ies=rank_ies,anchor_rank_ips=rank_ips,
  anchor_q_global_abs_logfc=q_global_abs_logfc,anchor_q_global_jsd=q_global_jsd,
  anchor_q_global_ies=q_global_ies,anchor_q_global_ips=q_global_ips,
  anchor_q_target_abs_logfc=q_target_abs_logfc,anchor_q_target_jsd=q_target_jsd,
  anchor_q_target_ies=q_target_ies,anchor_q_target_ips=q_target_ips)]
summary <- anchor_keep[summary,on="target_id"]
sig_counts <- joined[,.(n_qtarget_abs_logfc_005=sum(q_target_abs_logfc<0.05),
  n_qtarget_jsd_005=sum(q_target_jsd<0.05),n_qtarget_ies_005=sum(q_target_ies<0.05),
  n_qtarget_ips_005=sum(q_target_ips<0.05),n_qglobal_any_005=sum(pmin(q_global_abs_logfc,
  q_global_jsd,q_global_ies,q_global_ips)<0.05)),by=target_id]
summary <- sig_counts[summary,on="target_id"]

# Add the complete anchor-design provenance by the immutable target rank.
design <- fread(file.path(BASE,"saturation_candidate_manifest.tsv"))
design[, target_rank:=.I]
design_cols <- setdiff(names(design),c("target_id","apply_only","alphagenome_training_prohibited"))
summary <- design[,..design_cols][summary,on="target_rank"]
summary[, `:=`(apply_only=TRUE,alphagenome_training_prohibited=TRUE,
  canonical_pip_update_allowed=FALSE,convergence_update_allowed=FALSE,
  evidence_class="mechanistic_nomination",causal_evidence=FALSE)]
setorder(summary,target_rank)

fwrite(rbindlist(shard_qc),file.path(OUT,"saturation_score_shard_qc.tsv"),sep="\t")
fwrite(summary,file.path(OUT,"saturation_target_summary.tsv"),sep="\t")
fwrite(joined,file.path(OUT,"saturation_scores_joined.tsv.gz"),sep="\t",compress="gzip")
contract <- list(status="PASS",n_shards=27,n_unique_genomic_expected=expected_n,n_unique_genomic_scored=nrow(scored),
  n_target_variant_rows=nrow(joined),unique_target_perturbation_ids=uniqueN(joined$target_perturbation_id),
  n_targets=nrow(summary),snvs_per_target=1503,overlapping_genomic_snvs_scored_once=TRUE,
  metrics=list(abs_logfc="abs_logfc",jsd="jsd",IES="abs_logfc_x_jsd",
    IPS="abs_logfc_x_jsd_x_active_allele_quantile"),
  multiple_testing="BH separately per metric over 671187 unique genomic SNVs; descriptive target-wise BH also reported",
  uncertainty="empirical shuffled-sequence null only; no model-training uncertainty from a single held-out model per chromosome",
  evidence_class="mechanistic_nomination",causal_evidence=FALSE,
  apply_only_firewall=TRUE,canonical_pip_update_allowed=FALSE,
  convergence_update_allowed=FALSE,
  alphagenome_mode="zero_shot_only_no_training_or_calibration")
write_json(contract,file.path(OUT,"saturation_score_qc_contract.json"),pretty=TRUE,auto_unbox=TRUE)
message("[102] PASS: ",nrow(scored)," unique scores; ",nrow(summary)," targets")
