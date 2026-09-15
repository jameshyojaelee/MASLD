"""Hand re-derivation of step 43 features for one pilot region, by a separate code path (pandas groupby,
explicit ontology lists instead of lib_atlas.track_class)."""
import sys, gzip, csv, pathlib
import numpy as np, pandas as pd, anndata
pilot = pathlib.Path(sys.argv[1]); sat = pathlib.Path(sys.argv[2]); which = int(sys.argv[3])
rows = list(csv.DictReader(gzip.open(pilot / "tables/saturation_region_features.tsv.gz", "rt"), delimiter="\t"))
r = rows[which]
d = sat / "raw/saturation" / r["universe"] / f"{r['chrom']}_{r['start0']}-{r['end']}"
LIVER = {"UBERON:0002107", "UBERON:0001114", "UBERON:0001115"}
def feats(scorer, select, scale):
    ad = anndata.read_h5ad(d / f"{scorer}.h5ad")
    cols = [i for i in range(ad.n_vars) if select(ad.var.iloc[i])]
    if scale == "q":
        m = np.abs(np.asarray(ad.layers["quantiles"])[:, cols])
    elif scorer == "AVI_SCORE":
        m = np.abs(np.asarray(ad.X)[:, cols])
    else:
        m = np.abs(np.asarray(ad.X)[:, cols]) / ad.var["nonzero_mean"].astype(float).to_numpy()[cols]
    q = pd.DataFrame(m)
    q["pos"] = [int(v.split(":")[1]) for v in ad.obs["variant"]]
    per = q.groupby("pos").max()               # max over substitutions, per track
    prof = per.mean(axis=1)                    # mean over tracks
    p = prof.index.to_numpy(); x = prof.to_numpy()
    best = max(x[(p >= s) & (p < s + 50)].sum() for s in p)
    return len(cols), x.mean(), x.max(), x.sum(), best / x.sum(), (x > (0.25 if scale == 'q' else 0.10)).mean()
checks = {
 "liver_atac": ("ATAC", lambda v: v["ontology_curie"] in LIVER),
 "liver_dnase": ("DNASE", lambda v: v["ontology_curie"] in LIVER),
 "liver_h3k27ac": ("CHIP_HISTONE", lambda v: v["ontology_curie"] in LIVER and str(v["histone_mark"]).upper() == "H3K27AC"),
 "cage_liver": ("CAGE", lambda v: v["ontology_curie"] in LIVER),
 "chip_tf_liver": ("CHIP_TF", lambda v: v["ontology_curie"] in LIVER),
 "chip_tf_hepg2": ("CHIP_TF", lambda v: v["ontology_curie"] == "EFO:0001187"),
 "avi": ("AVI_SCORE", lambda v: True),
}
worst = 0.0
for (g, (scorer, sel)), scale in [(c, sc) for c in checks.items() for sc in ("rel", "q")]:
    n, *vals = feats(scorer, sel, scale)
    mine = [float(r[f"{g}_{scale}_{k}"]) for k in ("mean", "max", "sum", "share_best50", "frac_above_high")]
    diff = max(abs(a - b) for a, b in zip(vals, mine))
    worst = max(worst, diff)
    print(g, scale, "n_tracks", n, r[f"{g}_n_tracks"], "max_abs_diff", f"{diff:.2e}", [round(v, 4) for v in vals])
print("WORST", worst, r["region_key"])
