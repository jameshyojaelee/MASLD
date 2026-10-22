# DEBUG ONLY: stand-in refit outputs so the python stage can be exercised while the R refit runs.
import sys, pandas as pd, numpy as np
O=sys.argv[1]
d=pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
d=d.rename(columns={"padj":"adj.P.Val"})[["gene","logFC","AveExpr","t","P.Value","adj.P.Val"]]
d["t"]=d.t+np.random.default_rng(5).normal(0,0.5,len(d))
d.to_csv(O+"/refit_R1_hepadj_toptable.tsv.gz",sep="\t",index=False)
pd.DataFrame({"item":["R0_design_equals_model_design","R0_same_gene_order","R0_max_abs_t_diff","R0_pass","R1_VIF_hepatocyte_fraction","R1_canonical_calls_padj005_abslogFC05","R1_n_genes","hep_median_control","hep_median_disease"],"value":["TRUE","TRUE","0","TRUE","1.1","0","24115","0.9","0.88"]}).to_csv(O+"/refit_checks.tsv",sep="\t",index=False)
