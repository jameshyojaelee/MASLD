import sys, pathlib, numpy as np, anndata, csv, gzip
sat = pathlib.Path(sys.argv[1])
univ = list(csv.DictReader(gzip.open(sat / "tables/saturation_universe.tsv.gz", "rt"), delimiter="\t"))
rng = np.random.default_rng(3)
pick = [univ[i] for i in rng.choice(len(univ), 400, replace=False)]
pick = [u for u in pick if (sat / "raw/saturation" / u["universe"] / f"{u['chrom']}_{u['start0']}-{u['end']}" / "ATAC.h5ad").exists()][:6]
for u in pick:
    d = sat / "raw/saturation" / u["universe"] / f"{u['chrom']}_{u['start0']}-{u['end']}"
    for scorer in ("ATAC", "CHIP_HISTONE"):
        ad = anndata.read_h5ad(d / f"{scorer}.h5ad")
        liver = [i for i in range(ad.n_vars) if ad.var["ontology_curie"].iloc[i] in
                 ("UBERON:0002107", "UBERON:0001114", "UBERON:0001115")]
        q = np.abs(np.asarray(ad.layers["quantiles"])[:, liver]).ravel()
        x = np.abs(np.asarray(ad.X)[:, liver]).ravel()
        pc = lambda v: [round(float(np.percentile(v, p)), 4) for p in (5, 25, 50, 75, 95)]
        print(u["region_key"], scorer, "n_tracks", len(liver), "|quantile| pct", pc(q), "|raw| pct", pc(x))
