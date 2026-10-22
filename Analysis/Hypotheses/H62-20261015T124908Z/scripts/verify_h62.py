# Independent re-derivation of the H62 GSE66676 score and T1 from raw files (numpy/pandas, no R code reused).
import gzip, sys, numpy as np, pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; O=sys.argv[1]
L=gzip.open(R+"/data/public_r4/GSE66676/GSE66676_series_matrix.txt.gz","rt").read().split("\n")
b=[i for i,l in enumerate(L) if l.startswith("!series_matrix_table_begin")][0]; e=[i for i,l in enumerate(L) if l.startswith("!series_matrix_table_end")][0]
meta={}
for l in L[:b]:
    if l.startswith("!Sample_characteristics_ch1"):
        v=[x.strip('"') for x in l.split("\t")[1:]]; meta[v[0].split(":")[0]]=[x.split(": ",1)[1] for x in v]
gsm=[x.strip('"') for x in L[b+1].split("\t")[1:]]
X=pd.DataFrame([r.split("\t") for r in L[b+2:e]]); X.index=X[0].str.strip('"'); X=X.drop(columns=0).astype(float); X.columns=gsm
hist=np.array(meta["histology"]); sex=np.array(meta["gender"]); age=np.array(meta["age (yrs)"],float); bmi=np.array(meta["bmi"],float)
drop=["GSM1627766","GSM1627794"]; keep=np.array([g not in drop for g in gsm])
X=X.loc[:,keep]; hist=hist[keep]; sex=sex[keep]; age=age[keep]; bmi=bmi[keep]
naf=(hist!="No NAFLD").astype(float); ref=naf==0
m=pd.read_csv(O+"/signature_mapping.tsv",sep="\t",dtype={"probe":str})
Xs=X.loc[m.probe]; mu=Xs.loc[:,ref].mean(1); sd=Xs.loc[:,ref].std(1,ddof=1)
z=Xs.sub(mu,axis=0).div(sd,axis=0); z.index=pd.MultiIndex.from_arrays([m.direction,m.symbol])
A=z.groupby(level=[0,1]).mean()
score=A.loc["up"].mean(0)-A.loc["down"].mean(0)
s_run=pd.read_csv(O+"/scores_GSE66676.tsv",sep="\t").set_index("gsm")["score"]
print("n",len(score),"NoNAFLD",int(ref.sum()),"NAFLD",int(naf.sum()),"max |score diff| vs run",np.abs(score.values-s_run.loc[score.index].values).max())
D=np.column_stack([np.ones(len(naf)),naf,(sex=="male").astype(float),age,bmi]); y=score.values
beta,res,rk,_=np.linalg.lstsq(D,y,rcond=None); r=y-D@beta; df=len(y)-D.shape[1]; s2=r@r/df
se=np.sqrt(s2*np.linalg.inv(D.T@D)[1,1]); t=beta[1]/se
from scipy import stats
print("coef %.6f se %.6f t %.4f df %d p %.5f"%(beta[1],se,t,df,2*stats.t.sf(abs(t),df)))
a=y[naf==1]; c=y[naf==0]; sp=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(c)-1)*c.var(ddof=1))/(len(y)-2))
print("d_adj %.4f ci %.4f to %.4f"%(beta[1]/sp,(beta[1]-stats.t.ppf(.975,df)*se)/sp,(beta[1]+stats.t.ppf(.975,df)*se)/sp))
