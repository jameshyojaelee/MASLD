import pandas as pd, numpy as np
g=pd.read_csv('panel_7f_exposure_grid.tsv',sep='\t')
clean=['GSE136103','GSE174748','GSE189600','GSE202379','GSE244832']
for m in ['transcriptformer_tf_sapiens::linear','transcriptformer_tf_sapiens::two_layer_mlp','uce_4l::two_layer_mlp','uce_33l::two_layer_mlp','uce_4l::linear','uce_33l::linear','scimilarity_v1_1::linear','scimilarity_v1_1::two_layer_mlp']:
    d=g[(g.model_id==m)&(g.held_out_study.isin(clean))].set_index('held_out_study')
    eq=d.per_study_delta.mean(); dw=np.average(d.per_study_delta,weights=d.donors)
    drop=d.drop('GSE189600').per_study_delta.mean()
    print(f'{m}: clean5 equal-study {eq:+.6f}  donor-weighted {dw:+.6f}  drop GSE189600 {drop:+.6f}  GSE189600 delta {d.loc["GSE189600","per_study_delta"]:+.4f}')
