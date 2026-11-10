suppressPackageStartupMessages({library(metafor);library(jsonlite)})
OUT <- '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/transportability'
ex<-parse(file=file.path(OUT,'code/fit_transportability.R'))
for(e in ex)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('fit_one')))eval(e)
equal<-fit_one(rep(.6,5),rep(.1,5))
y<-c(-1,-.5,0,.5,1);se<-c(.2,.3,.2,.1,.2)
hetero<-fit_one(y,se)
objective<-function(tau){v<-se^2+tau;w<-1/v;b<-sum(w*y)/sum(w);.5*(sum(log(v))+log(sum(w))+sum(w*(y-b)^2))}
tau<-optimize(objective,c(0,100),tol=1e-10)$minimum
stopifnot(abs(hetero$tau2-tau)<1e-5,equal$kh_scale>=1,equal$SE>=.1/sqrt(5),
 abs(equal$beta-.6)<1e-10,is.na(fit_one(y[-1],se[-1])$pi_low))
write_json(list(equal_effect=TRUE,KH_floor=TRUE,k4_PI_suppressed=TRUE,independent_REML_tau2=tau,
 metafor_REML_tau2=hetero$tau2,abs_difference=abs(tau-hetero$tau2)),
 file.path(OUT,'summary_control_check.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
cat('PASS: known effect, KH floor, k4 PI suppression, independent restricted-likelihood optimum\n')
