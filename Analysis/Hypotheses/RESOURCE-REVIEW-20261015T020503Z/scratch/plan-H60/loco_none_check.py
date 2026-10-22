import pandas as pd, numpy as np
C="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/"
a=pd.read_csv(C+"bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
b=pd.read_csv(C+"bulk-loco-NONE-20261009T141221Z/deg_results.csv")
m=a.merge(b,on="gene")
print("primary vs LOCO-NONE genes",len(a),len(b),"joined",len(m),"max|dt|",np.abs(m.t_x-m.t_y).max())
