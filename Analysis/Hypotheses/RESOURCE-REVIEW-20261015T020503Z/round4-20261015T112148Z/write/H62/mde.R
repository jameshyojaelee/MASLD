# H62 prespec power: two-sided alpha 0.05, 34 vs 33, noncentral t; df 65 (unadjusted) and 62 (sex, age, BMI).
f<-function(d,df){ncp<-d/sqrt(1/34+1/33);q<-qt(.975,df);1-pt(q,df,ncp)+pt(-q,df,ncp)}
for(df in c(65,62)) cat(df, uniroot(function(d) f(d,df)-.8,c(.1,2))$root, "power at 0.67/0.97/1.32:", f(0.67,df), f(0.97,df), f(1.32,df), "power at 0.5:", f(0.5,df),"\n")
cat(power.t.test(n=33.5,power=.8)$delta,"\n")
