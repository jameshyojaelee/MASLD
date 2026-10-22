suppressPackageStartupMessages({library(data.table);library(lme4);library(lmerTest)})
out <- Sys.getenv('CCC_QUALIFIED_ROOT'); stopifnot(nzchar(out), !dir.exists(out));dir.create(out,recursive=TRUE)
root <- Sys.getenv('MASLD_PROJECT_ROOT');stopifnot(nzchar(root))
old <- file.path(root,'figures/candidates/fig4cd-ccc-bridge-2026-08-17-v13/source_tables')
roster_path <- file.path(root,'Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/current_source_diagnosis_donor_roster.tsv')
ccc_path <- file.path(root,'Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2/all_donor_lr_scores_v2_dc.tsv.gz')
source_full <- file.path(root,'scripts/figures/ccc_v3_dc_rebuild.R')
source_bridge <- file.path(root,'scripts/figures/render_fig4cd_ccc_bridge.R')
sha <- function(p) strsplit(system2('sha256sum',normalizePath(p),stdout=TRUE),' +')[[1]][1]
inputs <- c(roster_path,ccc_path,source_full,source_bridge,file.path(old,'fig4d_program_ccc_matched_donors.tsv'),file.path(old,'fig4d_program_ccc_associations.tsv'))
checks <- data.table(path=inputs,sha256=vapply(inputs,sha,character(1)));fwrite(checks,file.path(out,'input_sha256_before.tsv'),sep='\t')
extract <- function(path,start,end) {s<-paste(readLines(path,warn=FALSE),collapse='\n');a<-regexpr(start,s,fixed=TRUE)[1];stopifnot(a>0);tail<-substring(s,a);b<-regexpr(end,tail,fixed=TRUE)[1];stopifnot(b>0);eval(parse(text=substring(tail,1,b-1)),envir=.GlobalEnv)}
extract(source_bridge,'fit_bridge <- function(d) {','bridge <- matched[')
extract(source_bridge,'fit_lodo_beta <- function(d) {','lodo <- matched[')
extract(source_full,'fit_one <- function(d, axis_term) {','# ============================================================================\n# 1. Fit coarse LMM')
STAGE_LEVELS <- c('Healthy','Steatosis','Steatohepatitis','Cirrhosis');MIN_DONORS_PER_LR<-30;MIN_DONORS_PER_STRATUM<-10;N_CORES<-1
qual<-fread(roster_path);stopifnot(!anyDuplicated(qual$donor));qual[,exclude:=tolower(as.character(exclude))=='true']
original<-fread(file.path(old,'fig4d_program_ccc_matched_donors.tsv'));original[,`:=`(stage=factor(stage,levels=STAGE_LEVELS[1:3]),dataset=factor(dataset))]
prior<-fread(file.path(old,'fig4d_program_ccc_associations.tsv'))
reproduced<-original[,fit_bridge(.SD),by=.(program_uid,ct_pair,lr_pair)]
c<-merge(reproduced,prior,by=c('program_uid','ct_pair','lr_pair'),suffixes=c('_reproduced','_old'))
stopifnot(nrow(c)==26)
for (n in c('beta','se','pvalue','hc3_se','hc3_pvalue')) stopifnot(isTRUE(all.equal(c[[paste0(n,'_reproduced')]],c[[paste0(n,'_old')]],tolerance=1e-8)))
fwrite(c,file.path(out,'original_26_fit_reproduction.tsv'),sep='\t')
matched<-copy(original);matched[,c('disease_stage_coarse','stage'):=NULL]
matched<-merge(matched,qual[,.(donor,disease_stage_coarse,source_exclude=exclude)],by='donor',all.x=TRUE,sort=FALSE)
stopifnot(!anyNA(matched$source_exclude));dropped<-matched[source_exclude | !disease_stage_coarse%in%STAGE_LEVELS[1:3]]
fwrite(dropped,file.path(out,'excluded_26_input_rows.tsv'),sep='\t')
matched<-matched[!source_exclude & disease_stage_coarse%in%STAGE_LEVELS[1:3]]
matched[,`:=`(stage=factor(disease_stage_coarse,levels=STAGE_LEVELS[1:3]),dataset=factor(as.character(dataset)))]
bridge<-matched[,fit_bridge(.SD),by=.(program_uid,ct_pair,lr_pair)]
roster<-unique(prior[,.(program_uid,ct_pair,lr_pair,program_label,program_short,headline_label,sender,receiver,ligand_complex,receptor_complex,overlap_gene)])
bridge<-merge(bridge,roster,by=c('program_uid','ct_pair','lr_pair'),all=TRUE)
stopifnot(nrow(bridge)==26)
bridge[,`:=`(qvalue=p.adjust(pvalue,method='BH',n=26),hc3_qvalue=p.adjust(hc3_pvalue,method='BH',n=26))]
bridge[,evidence_state:=fifelse(!testable,'untestable',fifelse(qvalue<.05&hc3_qvalue<.05,'supported','indeterminate'))]
lodo<-matched[,rbindlist(lapply(sort(unique(as.character(dataset))),function(left_out){r<-fit_lodo_beta(.SD[as.character(dataset)!=left_out]);r[,excluded_dataset:=left_out];r})),by=.(program_uid,ct_pair,lr_pair)]
lodo<-merge(lodo,bridge[,.(program_uid,ct_pair,lr_pair,primary_beta=beta)],by=c('program_uid','ct_pair','lr_pair'))
lodo[,same_direction:=testable_lodo&sign(beta_lodo)==sign(primary_beta)]
lo<-lodo[,.(n_lodo_estimable=sum(testable_lodo),n_lodo_same_direction=sum(same_direction),lodo_beta_min=if(any(testable_lodo))min(beta_lodo[testable_lodo])else NA_real_,lodo_beta_max=if(any(testable_lodo))max(beta_lodo[testable_lodo])else NA_real_),by=.(program_uid,ct_pair,lr_pair)]
bridge<-merge(bridge,lo,by=c('program_uid','ct_pair','lr_pair'));bridge[,`:=`(biological_unit='biological donor',multiple_testing_family='2 fixed programs x13 fixed communication pairs;26 planned tests',claim_boundary='Candidate contextual covariation; no ligand secretion or receptor activation established')]
fwrite(bridge,file.path(out,'fig4d_program_ccc_associations.tsv'),sep='\t');fwrite(matched,file.path(out,'fig4d_program_ccc_matched_donors.tsv'),sep='\t');fwrite(lodo,file.path(out,'fig4d_program_ccc_lodo.tsv'),sep='\t')
# Full donor-level original LR family, not the selected thirteen pairs.
lr<-fread(ccc_path);stopifnot(all(lr$sample==lr$biological_donor));all_groups<-unique(lr[,.(ct_pair=paste(source,target,sep='->'),lr_pair=paste(ligand_complex,receptor_complex,sep='__'))]);pre<-unique(lr[,.(sample,dataset)])
pre<-merge(pre,qual[,.(sample=donor,source_stage=disease_stage_coarse,source_exclude=exclude)],by='sample',all.x=TRUE);stopifnot(all(qual[dataset=='Liver_Atlas']$exclude),all(pre[is.na(source_exclude)]$dataset=='Liver_Atlas'));pre[is.na(source_exclude),`:=`(source_exclude=TRUE,source_stage=NA_character_)];pre[,source_unit_qualification:=fifelse(grepl('^SRR',sample)&dataset=='Liver_Atlas','Original uncollapsed run ID; source-study excluded','Qualified donor key')];fwrite(pre,file.path(out,'full_ccc_source_donor_membership.tsv'),sep='\t')
lr<-merge(lr,qual[,.(sample=donor,disease_stage_coarse,source_exclude=exclude)],by='sample',all.x=TRUE);stopifnot(all(lr[is.na(source_exclude)]$dataset=='Liver_Atlas'));lr[is.na(source_exclude),source_exclude:=TRUE];lr<-lr[!source_exclude & disease_stage_coarse%in%STAGE_LEVELS]
lr[,`:=`(dataset_pooled=ifelse(dataset%in%c('GSE174748','GSE189600'),'Other_small',as.character(dataset)),score=-log10(pmax(magnitude_rank,1e-4)),disease_stage_coarse=factor(disease_stage_coarse,levels=STAGE_LEVELS),lr_pair=paste(ligand_complex,receptor_complex,sep='__'),ct_pair=paste(source,target,sep='->'))]
gate<-lr[,.(n_donors=.N,n_Healthy=sum(disease_stage_coarse=='Healthy'),n_Steatosis=sum(disease_stage_coarse=='Steatosis'),n_Steatohepatitis=sum(disease_stage_coarse=='Steatohepatitis'),n_Cirrhosis=sum(disease_stage_coarse=='Cirrhosis')),by=.(ct_pair,lr_pair)]
gate<-merge(all_groups,gate,by=c('ct_pair','lr_pair'),all.x=TRUE);for(n in c('n_donors','n_Healthy','n_Steatosis','n_Steatohepatitis','n_Cirrhosis'))set(gate,which(is.na(gate[[n]])),n,0L)
gate[,eligible_stage_counts:=n_donors>=30&n_Healthy>=10&n_Steatosis>=10&n_Steatohepatitis>=10&n_Cirrhosis>=10]
cat('Full source-qualified groups:',nrow(gate),'eligible stage counts:',sum(gate$eligible_stage_counts),'donors:',uniqueN(lr$sample),'\n')
full<-fit_coarse(lr,'source_qualified_dc');if(!nrow(full))full<-data.table(term=character(),ct_pair=character(),lr_pair=character(),Estimate=numeric(),pval=numeric(),q_bonferroni_family=numeric())
fwrite(full,file.path(out,'full_donor_lr_coarse_family.tsv'),sep='\t')
gate[,disposition:=fifelse(!eligible_stage_counts,'untestable_insufficient_stage_counts','eligible')]
fitted<-unique(full[,.(ct_pair,lr_pair)]);gate[fitted,on=.(ct_pair,lr_pair),disposition:='fitted'];gate[eligible_stage_counts&disposition!='fitted',disposition:='untestable_native_fit_failure']
fwrite(gate,file.path(out,'full_donor_lr_family_dispositions.tsv'),sep='\t')
fixed<-unique(prior[,.(ct_pair,lr_pair,headline_label)]);stopifnot(nrow(fixed)==13)
fixed<-merge(fixed,gate,by=c('ct_pair','lr_pair'),all.x=TRUE);fixed[is.na(disposition),disposition:='untestable_no_qualified_input_rows']
fixed<-merge(fixed,full[term=='disease_stage_coarseSteatohepatitis',.(ct_pair,lr_pair,dc_Estimate=Estimate,dc_pval=pval,dc_q=q_bonferroni_family)],by=c('ct_pair','lr_pair'),all.x=TRUE)
fwrite(fixed,file.path(out,'fixed_13_current_stage_context.tsv'),sep='\t')
stopifnot(identical(checks$sha256,unname(vapply(inputs,sha,character(1)))));fwrite(checks,file.path(out,'input_sha256_after.tsv'),sep='\t');writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(c('Source exclusions follow344_build_donor_metadata_extended.R and346_stage_ccc_mixedmodel.R; protocol/source-stage flags, not an arbitrary three-stage filter. Four-stage cirrhotic donors in qualifiedGSE202379 andGSE174748 retained. OriginalLiverAtlas42runrecords excluded; not counted as42donors.'),file.path(out,'source_eligibility_decision.txt'))
writeLines(sprintf('Original26fit reproduced;26planned models rerun, %d estimable,%d supported. Full source-qualified LR groups%d: %d stage-count eligible,%d fitted. Native within(ct_pair,term)BH then acrossctfamilies Bonferroni retained. Candidate only.',sum(bridge$testable),sum(bridge$evidence_state=='supported'),nrow(gate),sum(gate$eligible_stage_counts),nrow(fitted)),file.path(out,'CHECKED.txt'))
print(bridge[evidence_state=='supported',.(program_uid,ct_pair,lr_pair,n_donors,beta,qvalue,hc3_qvalue,n_lodo_same_direction,n_lodo_estimable,lodo_beta_min,lodo_beta_max)])
