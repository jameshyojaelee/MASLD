"""Re-check DR-main design numbers with a full-rank design (drop one cohort dummy)."""
import numpy as np, pandas as pd
exec(open(__file__.replace("dr_main_check.py", "n_and_design.py")).read())
keep = main.to_numpy()
Xk = X0[keep][:, np.abs(X0[keep]).sum(0) > 0]
print("DR-main design columns before fix", Xk.shape[1], "rank", np.linalg.matrix_rank(Xk))
cols = [c for c in dz.columns if c != "group_binaryDisease"]
cols_k = [c for c, j in zip(cols, np.abs(X0[keep]).sum(0) > 0) if j]
Xf = Xk[:, [i for i, c in enumerate(cols_k) if c != "datasetGSE130970"]]
print("fixed columns", [c for c in cols_k if c != "datasetGSE130970"], "rank", np.linalg.matrix_rank(Xf))
cse = calls.SE.median()
for nm, fr in [("rebuilt", new), ("earlier", old)]:
    Z = ilr(amalg(fr, DR5), 1e-6)[keep]
    k = infl(Xf, Z, y[keep], wt[keep])
    su = unscaled(np.column_stack([Xf, y[keep]]), wt[keep]) / unscaled(np.column_stack([X0, y]), wt)
    se_u = cse * su
    print(f"{nm} DRmain: inflation {k:.3f}; unadjusted SE ratio vs 774 {su:.3f}; median adjusted SE {se_u*k:.4f}; "
          f"MDE attenuation nominal {2.80*se_u*np.sqrt(k**2-1):.3f}; SE keep contrast {se_u*np.sqrt(k**2-0.75):.4f}")
