exec(open('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify2-R15-02/check_r1502.py').read().split("print(pd.crosstab")[0])
ind=m[m.call_state=='indeterminate']
x=ind[ind.susie_state_t12=='tested_le_0.5']
print('tested_le_0.5 & susie shown: max abf', x[x.src=='susie'].coloc_best_abf_pp4.max())
print('tested_le_0.5 & abf shown: min abf', x[x.src=='abf'].coloc_best_abf_pp4.min())
print(ind[ind.src=='other'][['symbol','dataset_id','est','coloc_best_susie_pp4','coloc_best_abf_pp4','susie_state_t12','driving_gwas']].to_string())
for g in ['GLP1R','THRB','HKDC1','MTARC1','PNPLA3','TM6SF2']:
    print(m[m.symbol==g][['symbol','dataset_id','est','coloc_best_susie_pp4','coloc_best_abf_pp4','susie_state_t12','call_state','src']].to_string(header=False))
