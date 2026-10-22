"""Fig 4H: collection-balanced feature-first aggregation of my independent node fits (r4g_nodes_367_grid.tsv.gz,
job 9841) to test whether the panel's system coordinates are reproduced, and sensitivity of rho / same-direction."""
import pandas as pd, numpy as np
from scipy import stats
out = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
S = pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/fig4h_disease_only_source.tsv", sep="\t")
N = pd.read_csv(f"{out}/r4g_nodes_367_grid.tsv.gz", sep="\t")
N = N[N.component_role == "main_system_component"]
def cb(df, col):
    return df.dropna(subset=[col]).groupby(["community_id", "collection"])[col].median().groupby("community_id").mean()
rows = []
for xm in ["zall", "zdis"]:
    for ym in ["zall", "zdis", "panel_mixed"]:
        for ax in ["axis_fp", "axis_sigpc1"]:
            gx = N[(N.ymode == xm) & (N.axis == ax)]
            if ym == "panel_mixed":   # pathways re-standardized over kept units; hotspot outcome_z unchanged
                gy = pd.concat([N[(N.ymode == "zdis") & (N.axis == ax) & (N.collection != "hotspot")],
                                N[(N.ymode == "zall") & (N.axis == ax) & (N.collection == "hotspot")]])
            else:
                gy = N[(N.ymode == ym) & (N.axis == ax)]
            x = cb(gx, "stage_beta"); y = cb(gy, "cont_beta")
            d = pd.DataFrame(dict(x=x, y=y)).reset_index().merge(S, on="community_id")
            rows.append(dict(x_ymode=xm, y_ymode=ym, axis=ax, max_dx=np.abs(d.x - d.stage_beta).max(), max_dy=np.abs(d.y - d.continuum_beta).max(),
                             rho=stats.spearmanr(d.x, d.y).correlation, same=int((np.sign(d.x) == np.sign(d.y)).sum()),
                             both_pos=int(((d.x > 0) & (d.y > 0)).sum()), both_neg=int(((d.x < 0) & (d.y < 0)).sum()),
                             n_absx_lt_0p02=int((d.x.abs() < 0.02).sum()), median_x=d.x.median(), median_y=d.y.median()))
            if xm == "zall" and ym == "panel_mixed" and ax == "axis_fp":
                d[["community_id", "figure_label", "x", "y", "stage_beta", "continuum_beta"]].to_csv(f"{out}/r4h_featurefirst_panel_estimand.tsv", sep="\t", index=False)
R = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(R.round(4).to_string())
R.to_csv(f"{out}/r4h_featurefirst_grid.tsv", sep="\t", index=False)
