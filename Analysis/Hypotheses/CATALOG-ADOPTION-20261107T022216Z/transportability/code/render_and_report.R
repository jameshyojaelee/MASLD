# KEY MESSAGE: Fixed disease-state effects have cohort-dependent transfer; stage replication is a separate contrast.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(jsonlite)})
OUT <- '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/transportability'
R<-file.path(OUT,'results');stopifnot(file.exists(file.path(R,'COMPLETE.json')))
T<-fread(file.path(R,'gene_transportability.tsv.gz'));C<-fread(file.path(R,'fixed_signature_cohort_summary.tsv'))
S<-fread(file.path(R,'existing_steatosite_grouped_stage_family.tsv'))
FS<-6;pdf.options(useDingbats=FALSE)
th<-theme_classic(base_size=FS,base_family='Helvetica')+theme(text=element_text(size=FS,colour='black'),
 axis.text=element_text(size=FS,colour='black'),axis.title=element_text(size=FS),plot.title=element_text(size=FS),
 legend.title=element_blank(),legend.text=element_text(size=FS),legend.position='bottom',plot.margin=margin(4,4,4,4))
P<-melt(C,id.vars=c('dataset','fixed_calls','testable'),measure.vars=c('same_direction','called_same_direction'),
 variable.name='measure',value.name='n')
P[,percent:=100*n/testable]
P[,measure:=factor(measure,levels=c('same_direction','called_same_direction'),labels=c('Same effect sign','Called with same sign'))]
fwrite(P,file.path(OUT,'figures/cohort_transfer_source.tsv'),sep='\t')
g<-ggplot(P,aes(percent,dataset,shape=measure,colour=measure))+geom_point(size=1.5,stroke=.4)+
 scale_shape_manual(values=c(16,1))+scale_colour_manual(values=c('#1565C0','#C9265E'))+
 scale_x_continuous(limits=c(0,100),breaks=c(0,50,100))+labs(x='Percent of testable fixed disease-state calls',y=NULL,
 title='Five discovery cohorts: retrospective agreement')+th
cairo_pdf(file.path(OUT,'figures/fixed_signature_cohort_transfer.pdf'),width=3.2,height=2,family='Helvetica');print(g);dev.off()
S[,contrast_label:=factor(contrast,levels=c('F2','F3F4'),labels=c('F2 vs F0-F1','F3-F4 vs F0-F1'))]
fwrite(S,file.path(OUT,'figures/existing_stage_replication_source.tsv'),sep='\t')
g<-ggplot(S,aes(slope,contrast_label,shape=testable))+geom_vline(xintercept=0,linewidth=.25,colour='#9E9E9E')+
 geom_errorbar(aes(xmin=ci_norm_lo,xmax=ci_norm_hi),width=.08,linewidth=.3,orientation='y')+
 geom_point(size=1.5,colour='#C9265E')+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),
 labels=c('FALSE'='Underpowered; descriptive','TRUE'='Passed fixed power rule'))+
 labs(x='SteatoSITE / discovery stage-effect slope (95% CI)',y=NULL,title='Previously observed independent-cohort stage comparison')+th
cairo_pdf(file.path(OUT,'figures/existing_stage_replication.pdf'),width=3.2,height=1.8,family='Helvetica');print(g);dev.off()
# These three examples were fixed in the manuscript before this sensitivity.
L<-fread(file.path(OUT,'cohort_effects_long.tsv.gz'))
examples<-T[symbol %in% c('GNMT','MAT1A','CYP2C19')]
F<-merge(L[native_testable==TRUE],examples[,.(gene_id_versioned=gene,symbol)],by='gene_id_versioned')
F<-F[,.(gene_id_versioned,symbol,source=dataset,effect=logFC,low=logFC-1.96*moderated_SE,
        high=logFC+1.96*moderated_SE,interval='cohort_normal_approximation_from_moderated_SE')]
F<-rbind(F,examples[,.(gene_id_versioned=gene,symbol,source='REML sensitivity',effect=beta,low=ci_low,high=ci_high,interval='modified_KH_95_percent')])
F[,source:=factor(source,levels=c(sort(unique(L$dataset)),'REML sensitivity'))]
fwrite(F,file.path(OUT,'figures/prespecified_examples_forest_source.tsv'),sep='\t')
g<-ggplot(F,aes(effect,source))+geom_vline(xintercept=0,linewidth=.25,colour='#9E9E9E')+
 geom_errorbar(aes(xmin=low,xmax=high),orientation='y',width=.08,linewidth=.3)+
 geom_point(size=1,colour='#1565C0')+facet_wrap(~symbol,nrow=1)+
 labs(x='Disease-minus-control log2FC',y=NULL,title='Prespecified tissue examples: cohort effects and REML sensitivity')+
 th+theme(strip.background=element_blank(),strip.text=element_text(size=FS))
cairo_pdf(file.path(OUT,'figures/prespecified_examples_forest.pdf'),width=5.1,height=2.2,family='Helvetica');print(g);dev.off()
fixed<-T[primary_call==TRUE];states<-fixed[,.N,by=sensitivity_state]
lines<-c('# Retrospective tissue-effect transportability',
 '', 'This analysis preserves the 1,854 primary disease-state calls. It adds a summary-data sensitivity using moderated input SEs; it is not a new blind confirmation.',
 '',sprintf('All 24,115 genes were retained in the BH family. Fits estimated %s genes; unavailable fits %s.',
 format(sum(T$status=='estimated'),big.mark=','),sum(T$status!='estimated')),
 sprintf('Cohort coverage: %s genes have five testable cohorts; %s have four because GSE126848 native expression filtering makes its other slots untestable.',sum(T$k==5),sum(T$k==4)),
 sprintf('Among the fixed 1,854 calls, %s retain same-direction random-effects BH support; %s have opposite-direction support; %s remain uncertain in this sensitivity.',
 sum(fixed$sensitivity_state=='same_direction_RE_support'),sum(fixed$sensitivity_state=='opposite_direction_RE_support'),sum(fixed$sensitivity_state=='RE_uncertain')),
 sprintf('This does not show absent effects: %s of 1,854 random-effects means keep the pooled sign; %s have one common sign in every testable cohort; %s have nominal random-effects p<0.05. These nominal values are not gene discoveries under the complete-family BH rule.',
 sum(fixed$same_sign_as_pooled),sum(fixed$n_positive==fixed$k | fixed$n_negative==fixed$k),sum(fixed$p<.05)),
 sprintf('The median I-squared among fixed calls is %.1f%%. Cochran Q is nominal p<0.05 for %s of them; these Q calls are descriptive and are not a corrected discovery family.',median(fixed$I2,na.rm=TRUE),sum(fixed$Q_p<.05,na.rm=TRUE)),
 '', 'The random-effects mean averages the effects of the observed study contrasts under an exchangeability approximation. It does not identify cell-autonomous effects or remove ancestry, control-definition, protocol, specimen, treatment, or clinical confounding.',
 'The small-study interval is modified Knapp-Hartung with variance inflation floored at one. An exploratory prediction interval is reported only with five cohorts and a positive heterogeneity estimate; it assumes a hypothetical exchangeable cohort and is not calibrated external prediction.',
 'The prespecified three-gene forest uses normal-approximation intervals (estimate +/-1.96 moderated SE) for the cohort rows and the modified Knapp-Hartung intervals for the summary rows; their inferential assumptions differ.',
 '', 'SteatoSITE remains a separate stage comparison. Its existing F3-F4 vs F0-F1 slope is 1.462 (95% donor-bootstrap interval 1.323-1.601), Holm p=0.001998 over two grouped contrasts, and 93.5% sign agreement. The F2 comparison fails its fixed power rule and is descriptive. Existing results were not rebranded as new confirmation.',
 'SteatoSITE source findings belong to Kendall et al. (2023); this output summarizes the already completed comparison of the fixed Resource effects with that public cohort.',
 '', 'The exposure inventory establishes no qualified unexposed disease/control cohort for this bounded task. The new blind-confirmation arm stops without opening additional expression outcomes. This does not establish that no suitable public cohort exists.',
 '', 'Outputs: gene_transportability.tsv.gz, fixed_1854_transportability.tsv.gz, proposed Catalog sidecar, complete long cohort-effect table, immutable unit metadata, compact individual PDFs and sources.',
 'Checks: five exact cohort refits verify effect/t reconstruction against eBayes moderated SEs; planted controls verify the variance floor and sign behavior; independent restricted-likelihood optimization checks the estimator; every fitted gene is checked against the explicit final-weight mean, variance and t-test formula.')
writeLines(lines,file.path(OUT,'REPORT.md'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'figures/sessionInfo.txt'))
