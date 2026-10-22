"""Fig 4G geometry/label check. Node centres (filled circles = pathway nodes, diamonds = Hotspot programs) are inverted
through the UMAP ticks and matched to the frozen latent coordinates (nearest neighbour); node fill colours are inverted
through the drawn colour bar (nearest RGB among its segments) and compared to my recomputed stage-adjusted continuum beta
(r4g_nodes_367_grid.tsv.gz; panel estimand = fixed projection, pathways re-standardized over disease units, Hotspot
outcome_z unchanged). Saturated nodes are compared after clipping. Also: S-label positions vs community medoids/centroids,
nodes hidden by label boxes. Usage: python g4g_geometry.py <pdf> <tag>"""
import sys
import numpy as np
from scipy.spatial import cKDTree
from pdfmarks import *
pdfp, tag = sys.argv[1], sys.argv[2]
out = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
N = pd.read_csv(f"{out}/r4g_nodes_367_grid.tsv.gz", sep="\t")
N = pd.concat([N[(N.ymode == "zdis") & (N.axis == "axis_fp") & (N.collection != "hotspot")],
               N[(N.ymode == "zall") & (N.axis == "axis_fp") & (N.collection == "hotspot")]])
G = pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_atlas_20260824/geometry/latent_coordinates.tsv.gz", sep="\t")
print("geometry nodes", len(G), G.component_role.value_counts().to_dict())
G = G.merge(N[["feature_id", "collection", "cont_beta"]], on=["feature_id", "collection"], how="left")
print("geometry nodes without a recomputed beta (untestable):", int(G.cont_beta.isna().sum()), G[G.cont_beta.isna()][["feature_id", "component_role"]].values.tolist())
pdf, pg, ch, w, r, l, c = load(pdfp)
num = w[w.text.str.fullmatch(r"-?\d+(\.\d+)?")].copy(); num["v"] = num.text.astype(float)
# UMAP axes
yaxis_x = l[(l.x0 == l.x1) & ((l.bottom - l.top) > 50)].x0.iloc[0]
xaxis_y = l[(l.top == l.bottom) & ((l.x1 - l.x0) > 50)].top.iloc[0]
tx = l[(l.x0 == l.x1) & ((l.bottom - l.top) < 2) & (l.top >= xaxis_y - 0.01) & (l.stroke == "#000000")]
ty = l[(l.top == l.bottom) & ((l.x1 - l.x0) < 2) & ((l.x1 - yaxis_x).abs() < 0.05)]
xl = num[(num.top > xaxis_y) & (num.top < xaxis_y + 10)]; yl = num[num.x1 < yaxis_x]
mx = [(t.x0, xl.loc[(xl.xc - t.x0).abs().idxmin(), "v"]) for _, t in tx.iterrows()]
my = [(t.top, yl.loc[(yl.yc - t.top).abs().idxmin(), "v"]) for _, t in ty.iterrows()]
ax, bx, rx = axis_fit(*zip(*mx)); ay, by, ry = axis_fit(*zip(*my))
print("UMAP x ticks", mx, "y ticks", my)
# colour bar
cb = r[(r.w < 0.5) & (r.h > 4)].sort_values("x0")
cbt = l[(l.stroke == "#FFFFFF")].drop_duplicates("x0")
cbl = num[num.top > cb.bottom.max()]
mc = [(t.x0, cbl.loc[(cbl.xc - t.x0).abs().idxmin(), "v"]) for _, t in cbt.iterrows()]
ac, bc, rc = axis_fit(*zip(*mc))
cb_vals = ac * ((cb.x0 + cb.x1) / 2).values + bc
cb_rgb = np.vstack([hex2rgb(h) for h in cb.fill])
lim_lo, lim_hi = ac * cb.x0.min() + bc, ac * cb.x1.max() + bc
print("colour bar ticks", mc, "| segments", len(cb), "| range", round(lim_lo, 4), round(lim_hi, 4), "| step", round(np.diff(cb_vals).mean(), 4))
nodes = c[(c.filled == True) & (c.w < 3)].copy()
nodes["shape"] = np.where(nodes.w > 2, "diamond", "circle")
nodes["u1"], nodes["u2"] = ax * nodes.xc + bx, ay * nodes.yc + by
tree = cKDTree(G[["umap_1", "umap_2"]].values)
d, idx = tree.query(nodes[["u1", "u2"]].values)
nodes["feature_id"] = G.feature_id.values[idx]; nodes["collection"] = G.collection.values[idx]; nodes["community_id"] = G.community_id.values[idx]
nodes["match_dist"] = d; nodes["beta_exp"] = G.cont_beta.values[idx]
rgb = np.vstack([hex2rgb(h) for h in nodes.fill])
k = np.argmin(((rgb[:, None, :] - cb_rgb[None]) ** 2).sum(-1), axis=1)
nodes["beta_drawn"] = cb_vals[k]; nodes["rgb_err"] = np.sqrt(((rgb - cb_rgb[k]) ** 2).sum(-1))
nodes["beta_exp_clip"] = nodes.beta_exp.clip(lim_lo, lim_hi)
step = np.diff(cb_vals).mean(); tol = max(0.005 * (lim_hi - lim_lo), 1.5 * step)
nodes["off"] = (nodes.beta_drawn - nodes.beta_exp_clip).abs() > tol
print("nodes drawn", len(nodes), nodes["shape"].value_counts().to_dict(), "| unique matched", nodes.feature_id.nunique(), "| max match dist (UMAP units)", round(nodes.match_dist.max(), 4))
print("diamonds are Hotspot:", (nodes[nodes["shape"] == "diamond"].collection == "hotspot").mean(), "| circles are pathways:", (nodes[nodes["shape"] == "circle"].collection != "hotspot").mean())
print("tolerance", round(tol, 4), "| checked", len(nodes), "| off", int(nodes.off.sum()), "| max|d|", round((nodes.beta_drawn - nodes.beta_exp_clip).abs().max(), 4), "| max rgb err", round(nodes.rgb_err.max(), 4))
print(nodes[nodes.off][["feature_id", "collection", "beta_drawn", "beta_exp", "rgb_err"]].head(15).to_string())
sat = nodes[(nodes.beta_exp < lim_lo) | (nodes.beta_exp > lim_hi)]
print("saturated nodes", len(sat), "| above", int((nodes.beta_exp > lim_hi).sum()), "below", int((nodes.beta_exp < lim_lo).sum()), "| max |beta|", round(nodes.beta_exp.abs().max(), 3))
print("nodes with NaN expected beta drawn:", int(nodes.beta_exp.isna().sum()))
# label boxes and hidden nodes
boxes = c[(c.w > 10) & (c.fill == "#FFFFFF")]
hid = 0
for _, b in boxes.iterrows():
    hid += int(((nodes.xc > b.x0) & (nodes.xc < b.x1) & (nodes.yc > b.top) & (nodes.yc < b.bottom)).sum())
print("label boxes", len(boxes), "| node centres under label boxes", hid)
# labels vs community centroid: nearest community centroid to the leader-line end that is away from the label
cent = G.groupby("community_id")[["umap_1", "umap_2"]].median()
lead = l[(l.lw.round(2) == 0.3) & (l.stroke == "#000000")]
lab = w[w.text.str.fullmatch(r"S\d\d")]
lr = []
for _, t in lab.iterrows():
    cid = "system_%s" % t.text[1:]
    # leader end farthest from the label centre among leaders touching the label box
    best = None
    for _, ln in lead.iterrows():
        ends = [(ln.x0, ln.top), (ln.x1, ln.bottom), (ln.x0, ln.bottom), (ln.x1, ln.top)]
        dlab = min(np.hypot(e[0] - t.xc, e[1] - t.yc) for e in ends)
        if best is None or dlab < best[0]:
            best = (dlab, ln)
    ln = best[1]
    ends = [(ln.x0, ln.top), (ln.x1, ln.bottom)] if True else None
    far = max([(ln.x0, ln.top), (ln.x1, ln.bottom), (ln.x0, ln.bottom), (ln.x1, ln.top)], key=lambda e: np.hypot(e[0] - t.xc, e[1] - t.yc))
    fu = (ax * far[0] + bx, ay * far[1] + by)
    dist_c = np.hypot(cent.umap_1 - fu[0], cent.umap_2 - fu[1])
    own = G[G.community_id == cid]
    dn = np.hypot(own.umap_1 - fu[0], own.umap_2 - fu[1]).min()
    near = nodes.assign(dd=np.hypot(nodes.u1 - fu[0], nodes.u2 - fu[1])).nsmallest(5, "dd").community_id.tolist()
    lr.append(dict(label=t.text, leader_end_u1=fu[0], leader_end_u2=fu[1], nearest_centroid=dist_c.idxmin(), dist_to_own_centroid=dist_c[cid],
                   dist_to_nearest_own_node=dn, five_nearest_nodes_own=sum(x == cid for x in near)))
L = pd.DataFrame(lr); print(L.round(3).to_string())
txt = " ".join(w.text)
print("x title", "UMAP 1" in txt, "; y title", "".join(ch[ch.x1 < 53].sort_values("top", ascending=False).text), "; colour-bar title", "continuum beta" in txt, "; n printed", any(k in txt for k in ["367", "n ="]))
nodes.drop(columns=["pts"], errors="ignore").to_csv(f"{out}/g4g_{tag}.tsv.gz", sep="\t", index=False); L.to_csv(f"{out}/g4g_{tag}_labels.tsv", sep="\t", index=False)
