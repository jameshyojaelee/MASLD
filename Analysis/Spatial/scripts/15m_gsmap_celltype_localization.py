#!/usr/bin/env python
"""15m — Depth-corrected cell-type localization of MASLD spatial GWAS risk.

Which cell types carry each MASLD trait-family's genetic risk, in tissue? gsMap's
per-spot -log10P heritability is strongly confounded by sequencing depth
(r(-log10P, total_counts) ~= 0.7-0.76), so a raw correlation with cell-type
abundance is an artifact. This script computes the DEPTH-CORRECTED partial
Spearman correlation between per-spot -log10P and each cell type's cell2location
proportion, controlling total_counts, per section, then aggregates within each
Visium cohort and reports cross-cohort concordance.

Verified finding this reproduces (both cohorts, after depth correction):
  NAFLD-diagnosis risk -> Hepatocytes / anti-Fibroblast
  PDFF (liver-fat)     risk -> Fibroblasts / anti-Hepatocyte
  (enzyme + Endothelial do NOT replicate across cohorts -> QC/negative contrast)

Inputs:
  results/gsmap/{gse192741,vu}/<section>/spatial_ldsc/<section>_<trait>.csv.gz  (per-spot p)
  results/cell2location/spatial_model/spatial_deconvolved.h5ad                  (GSE: obs c2l cols)
  results/cell2location/spatial_model_vu/spatial_deconvolved_vu.h5ad            (Vu: obsm q05)
Output:
  results/gsmap/gsmap_celltype_localization.csv  (cohort x trait_family x cell_type: mean/CI/sign)
Env: spatial (scanpy/anndata). Compute node, ~48-64 G.
"""
import os
import re

import numpy as np
import pandas as pd
import anndata as ad

BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
GS = os.path.join(BASE, "Analysis/Spatial/results/gsmap")
C2L = os.path.join(BASE, "Analysis/Spatial/results/cell2location")
PREF = "c2l_q05cell_abundance_w_sf_means_per_cluster_mu_fg_"

FAMILIES = {
    "enzyme": ["ukbb_alt", "mvp_alt", "ukbb_ast", "mvp_ast", "ukbb_ggt"],
    "NAFLD": ["finngen_nafld", "ghodsian_nafld", "mvp_nafld", "ukbb2023_nafld", "anstee2020_nafld"],
    "PDFF": ["pdff", "pdff_2021a", "pdff_2021b", "pdff_2022"],
    "NASH": ["finngen_nash"],
}
# all cell types are tested (derived from the cell2location model columns at runtime);
# this list only controls the compact console summary at the end.
CELLTYPES = ["Hepatocytes", "Fibroblasts", "Cholangiocytes", "Macrophages",
             "Endothelial cells", "Mono+mono derived cells", "B cells"]


def strip_bc(bcs):
    out = []
    for b in bcs:
        m = re.match(r"^([ACGTN]{16}-\d+)", str(b))
        out.append(m.group(1) if m else str(b))
    return out


def rk(x):
    return pd.Series(x).rank().values


def _resid(a, c):
    a, c = rk(a), rk(c)
    b = np.polyfit(c, a, 1)
    return a - (b[0] * c + b[1])


def partial(y, x, z):
    """partial Spearman r(y, x | z)."""
    return np.corrcoef(_resid(y, z), _resid(x, z))[0, 1]


def load_cohort(cohort):
    """Return per-spot DataFrame [bc, sample_id, total_counts, <celltype proportions>]."""
    if cohort == "gse192741":
        A = ad.read_h5ad(os.path.join(C2L, "spatial_model/spatial_deconvolved.h5ad"))
        ctcols = [c for c in A.obs.columns if c.startswith(PREF)]
        names = [c[len(PREF):] for c in ctcols]
        ab = A.obs[ctcols].copy(); ab.columns = names
    else:
        # Vu: proportions in obsm; cell-type name order == GSE model order (verified corr>0.99)
        G = ad.read_h5ad(os.path.join(C2L, "spatial_model/spatial_deconvolved.h5ad"), backed="r")
        names = [c[len(PREF):] for c in G.obs.columns if c.startswith(PREF)]
        A = ad.read_h5ad(os.path.join(C2L, "spatial_model_vu/spatial_deconvolved_vu.h5ad"))
        # obsm may be a pandas DataFrame with its OWN column names; take values positionally
        # (column order == GSE model order, verified corr>0.99) and rename to `names`.
        ab = pd.DataFrame(np.asarray(A.obsm["q05_cell_abundance_w_sf"]),
                          index=A.obs_names, columns=names)
    prop = ab.div(ab.sum(1), axis=0)
    prop["bc"] = strip_bc(A.obs_names)
    prop["sample_id"] = A.obs["sample_id"].astype(str).values
    prop["tc"] = A.obs["total_counts"].values
    return prop


def consensus_logp(section_dir, traits, bc):
    mats = []
    for t in traits:
        f = os.path.join(section_dir, "spatial_ldsc", f"{os.path.basename(section_dir)}_{t}.csv.gz")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        d["spot"] = strip_bc(d["spot"])
        s = d.set_index("spot")["p"].apply(lambda x: -np.log10(max(x, 1e-300)))
        mats.append(s.reindex(bc).values)
    if not mats:
        return None
    return np.nanmean(np.vstack(mats), axis=0)


def main():
    rows = []
    depth_r = []
    for cohort in ["gse192741", "vu"]:
        prop = load_cohort(cohort)
        sections = sorted(os.listdir(os.path.join(GS, cohort)))
        sections = [s for s in sections if os.path.isdir(os.path.join(GS, cohort, s))]
        print(f"\n[{cohort}] {len(sections)} sections, {prop['sample_id'].nunique()} c2l samples")
        all_celltypes = [c for c in prop.columns if c not in ("bc", "sample_id", "tc")]
        for fam, traits in FAMILIES.items():
            for ct in all_celltypes:
                vals = []
                for sec in sections:
                    sid = sec.replace(f"{cohort}_", "").replace("gse192741_", "") if cohort == "gse192741" \
                        else sec.replace("vu_", "")
                    sp = prop[prop["sample_id"] == sid]
                    if len(sp) < 50:
                        continue
                    sp = sp.set_index("bc")
                    p = consensus_logp(os.path.join(GS, cohort, sec), traits, list(sp.index))
                    if p is None:
                        continue
                    j = sp.assign(p=p).dropna(subset=["tc", "p", ct])
                    if len(j) < 50:
                        continue
                    vals.append(partial(j["p"].values, j[ct].values, j["tc"].values))
                    if fam == "NAFLD" and ct == "Hepatocytes":
                        depth_r.append((cohort, sec,
                                        np.corrcoef(rk(j["p"]), rk(j["tc"]))[0, 1]))
                v = np.array(vals)
                if len(v) == 0:
                    continue
                rows.append(dict(cohort=cohort, trait_family=fam, cell_type=ct,
                                 mean_partial_r=v.mean(), sd=v.std(), n_sections=len(v),
                                 n_pos=int((v > 0).sum()), n_neg=int((v < 0).sum()),
                                 n_same_sign=int(max((v > 0).sum(), (v < 0).sum()))))
    df = pd.DataFrame(rows)
    print("rows per cohort:", df["cohort"].value_counts().to_dict())

    # cross-cohort concordance: same sign of mean_partial_r in both cohorts
    piv = df.pivot_table(index=["trait_family", "cell_type"], columns="cohort",
                         values="mean_partial_r")
    if "gse192741" in piv.columns and "vu" in piv.columns:
        conc = (np.sign(piv["gse192741"]) == np.sign(piv["vu"]))
    else:
        conc = pd.Series(False, index=piv.index)
    df = df.merge(conc.rename("cross_cohort_concordant").reset_index(),
                  on=["trait_family", "cell_type"], how="left")

    out = os.path.join(GS, "gsmap_celltype_localization.csv")
    df.to_csv(out, index=False)
    print(f"\nwrote {out} ({len(df)} rows)")
    print(f"depth confound r(-log10P NAFLD, total_counts) mean: "
          f"{np.mean([r for _,_,r in depth_r]):.2f}")

    print("\n=== depth-corrected partial r (mean | n_same_sign/n_sections), key contrast ===")
    print(f"{'family':8}{'cell_type':16}{'GSE':>16}{'Vu':>16}{'concord':>9}")
    for fam in ["NAFLD", "PDFF", "enzyme"]:
        for ct in ["Hepatocytes", "Fibroblasts", "Macrophages"]:
            g = df[(df.cohort == "gse192741") & (df.trait_family == fam) & (df.cell_type == ct)]
            v = df[(df.cohort == "vu") & (df.trait_family == fam) & (df.cell_type == ct)]
            if g.empty or v.empty:
                continue
            gs = f"{g.mean_partial_r.iloc[0]:+.2f}({g.n_same_sign.iloc[0]}/{g.n_sections.iloc[0]})"
            vs = f"{v.mean_partial_r.iloc[0]:+.2f}({v.n_same_sign.iloc[0]}/{v.n_sections.iloc[0]})"
            cc = "yes" if bool(g.cross_cohort_concordant.iloc[0]) else "NO"
            print(f"{fam:8}{ct:16}{gs:>16}{vs:>16}{cc:>9}")


if __name__ == "__main__":
    main()
