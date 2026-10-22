"""Fig 4D recomputation: discovery lineage x content-associated lineage counts for 117 programs.
Upstream: Analysis/SingleCell/results_gpu_v2/program_label_vs_content/plc-20260907T105813Z/program_label_vs_content.tsv
(per-program lineage content z-scores). Content lineage recomputed as argmax of the six z columns.
Bootstrap stability per program from R04-09 joint bootstrap (seed 20260923, B=1000)."""
import pandas as pd, numpy as np
root = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out = f"{root}/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
T = pd.read_csv(f"{root}/Analysis/SingleCell/results_gpu_v2/program_label_vs_content/plc-20260907T105813Z/program_label_vs_content.tsv", sep="\t")
zc = [c for c in T.columns if c.startswith("z_")]
T["content_re"] = T[zc].idxmax(axis=1).str.replace("z_", "").str.replace("_", " ")
T["content_re"] = T.content_re.replace({"T cells": "T_cells", "Endothelial cells": "Endothelial_cells"})
short = {"Hepatocytes": "Hep", "Fibroblasts": "Fib", "Cholangiocytes": "Chol", "Macrophages": "Mac", "T_cells": "T", "T cells": "T",
         "Endothelial_cells": "Endo", "Endothelial cells": "Endo"}
T["disc"] = T.declared_lineage.map(short); T["cont"] = T.content_re.map(short)
T["cont_file"] = T.content_lineage.map(short)
print("n programs", len(T), "; recomputed content == file content:", int((T.cont == T.cont_file).sum()))
order_r = ["Hep", "Fib", "Chol", "Mac", "T"]; order_c = ["Hep", "Fib", "Chol", "Mac", "T", "Endo"]
X = pd.crosstab(T.disc, T.cont).reindex(index=order_r, columns=order_c, fill_value=0)
print(X); print("diagonal", int(np.trace(X.values[:, :5])), "off", int(X.values.sum() - np.trace(X.values[:, :5])))
X.to_csv(f"{out}/r4d_counts.tsv", sep="\t")
B = pd.read_csv(f"{root}/figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/source_tables/R04-09_content_lineage_bootstrap_joint_seed20260923_B1000.tsv", sep="\t")
T = T.merge(B[["program_uid", "boot_same_call"]], on="program_uid", how="left")
T["unstable"] = T.boot_same_call < 0.95
U = pd.crosstab(T[T.unstable].disc, T[T.unstable].cont).reindex(index=order_r, columns=order_c, fill_value=0)
print("unstable calls by cell (boot_same_call < 0.95):"); print(U)
print(T[T.unstable][["disc", "cont", "module", "module_name", "boot_same_call", "margin"]].to_string())
# margin: how close the call is
print("programs with margin < 0.1:", int((T.margin < 0.1).sum()))
T.to_csv(f"{out}/r4d_programs.tsv", sep="\t", index=False)
