#!/usr/bin/env python3
"""R14-02 check: split the Fig 2D direct class by trait (PDFF vs NAFLD/NASH).

Reuses the adopted producer functions unchanged (copied from the
fig2-coordinate-rerun-20261013T213349Z_22218816 snapshot). Input is the gated,
coordinate-restored table that the adopted robustness run wrote. The reference
fit is reproduced first with the producer's own fit id, so the loader is shown
to regenerate the adopted ratio and interval before any subset is fitted.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from fig2_multiplicity_correction import (SCOPES, CORRECTED_CODING, SEED, N_BOOT,  # noqa: E402
                                          fit_rng, assign_clusters, cluster_table, analyse)
from fig2_coding_mass_robustness import coding_frame  # noqa: E402

RUN = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/"
           "fig2-coordinate-rerun-20261013T213349Z_22218816")
arch = pd.read_csv(RUN / "robustness/corrected_locus_coordinates.tsv", sep="\t",
                   dtype={"locus": str})
arch = arch[arch.trait_scope.isin(SCOPES)].copy()
print("rows by scope:", arch.groupby("trait_scope").size().to_dict(), flush=True)
assert len(arch) == 6415, len(arch)
arch = coding_frame(arch)

direct = arch.trait_scope == SCOPES[0]
print("direct sets by trait:", arch[direct].groupby("trait").size().to_dict())
print("direct studies by trait:",
      arch[direct].groupby("trait")["study"].nunique().to_dict(), flush=True)

rows = []


def fit(frame, arm, dropped, recluster=True):
    sub = assign_clusters(frame.copy()) if recluster else frame.copy()
    fit_id = f"{arm}|{dropped}"
    r = analyse(cluster_table(sub, "coding_mass", "arch_total"), fit_rng(fit_id))
    d = sub[sub.trait_scope == SCOPES[0]]
    r.update(arm=arm, dropped=dropped, fit_id=fit_id, recluster=recluster,
             n_direct_sets=int(len(d)), n_direct_studies=int(d.study.nunique()),
             direct_traits=",".join(sorted(d.trait.unique())))
    rows.append(r)
    print(f"[{arm}] {dropped:<28} sets {r['n_direct_sets']:>4} studies "
          f"{r['n_direct_studies']:>2} clusters d/e/phys/shared "
          f"{r.get('n_clusters_direct')}/{r.get('n_clusters_enzyme')}/"
          f"{r.get('n_physical_clusters')}/{r.get('n_shared_clusters')} "
          f"frac_d {r.get('fraction_direct', np.nan):.4f} frac_e {r.get('fraction_enzyme', np.nan):.4f} "
          f"ratio {r.get('ratio', np.nan):.3f} [{r.get('ratio_boot_lo', np.nan):.3f}, "
          f"{r.get('ratio_boot_hi', np.nan):.3f}] diff {r.get('difference', np.nan):.4f} "
          f"[{r.get('diff_boot_lo', np.nan):.4f}, {r.get('diff_boot_hi', np.nan):.4f}] "
          f"nominal_p {r.get('nominal_p_bootstrap_diff', np.nan):.4g}", flush=True)
    return r


# 0. reference, same fit id as the adopted robustness producer
r0 = fit(arch, "reference_1Mb", "none")
for key, expected in CORRECTED_CODING.items():
    tol = 1e-9 if key.startswith("fraction") else 1e-7 * expected
    assert abs(r0[key] - expected) <= tol, (key, r0[key], expected)
adopted = pd.read_csv(RUN / "robustness/coding_mass_robustness.tsv", sep="\t")
a0 = adopted[adopted.fit_id == "reference_1Mb|none"].iloc[0]
assert abs(a0.ratio_boot_lo - r0["ratio_boot_lo"]) < 1e-9, (a0.ratio_boot_lo, r0["ratio_boot_lo"])
assert abs(a0.ratio_boot_hi - r0["ratio_boot_hi"]) < 1e-9
print("[check] adopted reference ratio and interval reproduced exactly", flush=True)
# also reproduce the adopted PDFF-drop fit
st = "2022_36402844_PDFF_EUR"
rp = fit(arch[~(direct & (arch.study == st))], "leave_one_study_out_direct", st)
ap = adopted[adopted.fit_id == f"leave_one_study_out_direct|{st}"].iloc[0]
assert abs(ap.ratio - rp["ratio"]) < 1e-9 and abs(ap.ratio_boot_lo - rp["ratio_boot_lo"]) < 1e-9
print("[check] adopted leave-out-36402844 fit reproduced exactly", flush=True)

pdff = direct & (arch.trait == "PDFF")
nafld = direct & (arch.trait == "NAFLD")
nash = direct & (arch.trait == "NASH")

# 1. trait split, reclustered on pooled positions (as the LOSO arm does)
fit(arch[~(nafld | nash)], "trait_split", "PDFF_only")
fit(arch[~pdff], "trait_split", "NAFLD_NASH_only")
fit(arch[~(pdff | nash)], "trait_split", "NAFLD_only")
fit(arch[~(pdff | nafld)], "trait_split", "NASH_only")

# 2. same split on the reference clustering (no reclustering)
ref = assign_clusters(arch.copy())
rd = ref.trait_scope == SCOPES[0]
fit(ref[~(rd & ref.trait.isin(["NAFLD", "NASH"]))], "trait_split_refclusters", "PDFF_only", recluster=False)
fit(ref[~(rd & (ref.trait == "PDFF"))], "trait_split_refclusters", "NAFLD_NASH_only", recluster=False)

# 3. seed sensitivity of the NAFLD/NASH interval: five alternative fit ids
for k in range(5):
    fit(arch[~pdff], f"trait_split_seed{k}", "NAFLD_NASH_only")
for k in range(3):
    fit(arch[~(nafld | nash)], f"trait_split_seed{k}", "PDFF_only")

# 4. cluster overlap between PDFF and NAFLD/NASH direct clusters on reference clusters
cp = set(ref.loc[rd & (ref.trait == "PDFF"), "cluster_id"])
cn = set(ref.loc[rd & ref.trait.isin(["NAFLD", "NASH"]), "cluster_id"])
print(f"[clusters] reference direct clusters PDFF {len(cp)} NAFLD/NASH {len(cn)} "
      f"both {len(cp & cn)} union {len(cp | cn)}", flush=True)
# coding mass per cluster by trait group (cluster means), reference clustering
for name, mask in (("PDFF", rd & (ref.trait == "PDFF")),
                   ("NAFLD_NASH", rd & ref.trait.isin(["NAFLD", "NASH"]))):
    g = ref[mask].groupby("cluster_id")[["coding_mass", "arch_total"]].mean()
    g["frac"] = g.coding_mass / g.arch_total
    print(f"[percluster] {name}: clusters {len(g)}, sum coding {g.coding_mass.sum():.4f}, "
          f"sum total {g.arch_total.sum():.4f}, pooled frac {g.coding_mass.sum()/g.arch_total.sum():.4f}")
    print(g.sort_values("coding_mass", ascending=False).round(4).to_string(), flush=True)

out = pd.DataFrame(rows)
out["seed"] = SEED
out["n_boot"] = N_BOOT
out.to_csv(HERE / "trait_split_results.tsv", sep="\t", index=False)
print("[done] wrote", HERE / "trait_split_results.tsv")
