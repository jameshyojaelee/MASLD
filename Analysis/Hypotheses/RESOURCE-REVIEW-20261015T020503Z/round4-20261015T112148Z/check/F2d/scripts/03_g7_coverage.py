"""Independent fine-mapping coverage of lead windows per ancestry and class.
Rule: lead windows = 1-Mb single linkage on lead BP; covered if a lead locus
string (trailing zeros then trailing dot stripped) equals a credible-set locus
string of the same study. Also a position rule: any credible-set member within
the window span +-500 kb."""
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = ROOT + "/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2d/out"
reg = pd.read_csv(ROOT + "/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/snapshot/GWAS/finemapping/config/gwas_registry.tsv", sep="\t")
tier = pd.read_csv(ROOT + "/GWAS/finemapping/config/gwas_trait_tier.tsv", sep="\t")
coords = pd.read_csv(ROOT + "/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv", sep="\t")
studies = sorted(coords.study.unique())
cs = pd.read_csv(ROOT + "/GWAS/finemapping/results/credible_sets.csv",
                 usecols=["chromosome", "position", "locus", "study"], dtype={"locus": str})


def norm(s):
    s = str(s)
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s

cs["lk"] = cs.locus.map(norm)
rows = []
for st in studies:
    r = reg[reg.study_name == st].iloc[0]
    sc = coords[coords.study == st].trait_scope.iloc[0]
    anc = coords[coords.study == st].ancestry.iloc[0]
    path = ROOT + "/GWAS/finemapping/" + r.leadsnps_path
    try:
        lead = pd.read_csv(path, sep="\t", dtype={"locus": str})
    except Exception as ex:
        lead = pd.DataFrame(columns=["CHR", "BP", "locus"])
    csst = cs[cs.study == st]
    keys = set(csst.lk)
    nwin = ncov = ncov_pos = 0
    if len(lead):
        lead = lead.sort_values(["CHR", "BP"])
        for c, b in lead.groupby("CHR"):
            bp = b.BP.to_numpy()
            w = np.cumsum(np.concatenate([[True], np.diff(bp) > 1_000_000]))
            mem = csst[csst.chromosome.astype(str) == str(c)].position.to_numpy()
            for wi in np.unique(w):
                bb = b[w == wi]
                nwin += 1
                if any(norm(x) in keys for x in bb.locus):
                    ncov += 1
                lo, hi = bb.BP.min() - 500_000, bb.BP.max() + 500_000
                if ((mem >= lo) & (mem <= hi)).any():
                    ncov_pos += 1
    nfm = csst.groupby(["chromosome", "lk"]).ngroups
    rows.append(dict(study=st, scope=sc, ancestry=anc, lead_rows=len(lead), lead_windows=nwin,
                     covered=ncov, covered_pos500=ncov_pos, finemapped_windows_in_cs=nfm,
                     n_cs_gated=int((coords.study == st).sum())))
t = pd.DataFrame(rows)
print(t.to_string())
s = t[t.lead_windows > 0].groupby(["ancestry", "scope"])[["lead_windows", "covered", "covered_pos500"]].sum()
s["cov"] = s.covered / s.lead_windows
print(s.to_string())
t.to_csv(OUT + "/g7_per_study_check.tsv", sep="\t", index=False)
s.to_csv(OUT + "/g7_coverage_check.tsv", sep="\t")
