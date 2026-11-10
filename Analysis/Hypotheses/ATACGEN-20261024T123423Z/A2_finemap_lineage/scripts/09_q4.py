#!/usr/bin/env python3
"""A2 Q4: worked loci with the largest share of fine-mapped posterior in lineage-specific
peaks (liver classes), with the COLOC r2 genes colocalized at the same variant. Descriptive."""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from a2common import LIN, RES, ROOT, W

TOPN = 15
cs = pd.read_csv(RES / "q1_cs_level.tsv.gz", sep="\t")
cs = cs[(cs.method == "matched_tss_maf_ld") &
        cs.trait_class.isin(["direct", "liver_enzyme", "sequela"])]
mem = pd.read_parquet(W / "q1_members_annot.parquet")
pk = pd.read_parquet(W / "peaks_index.parquet").set_index("pk_index")
cv = pd.read_parquet(W / "q2_coloc_annot.parquet",
                     columns=["gene", "gwas_name", "pp_h4", "chromosome", "position", "w"])
cvh = cv[cv.w >= 0.05]
REF = ROOT / "GWAS/finemapping/data/ld_ref"
POPS = ["eur", "eas", "afr", "amr", "sas"]
bims = {}


def rsid(c, p, anc):
    """rsID from the 1000G panel bims (study ancestry first); "" if only a positional ID."""
    for pop in [anc.lower()] + [q for q in POPS if q != anc.lower()]:
        if (pop, c) not in bims:
            b = pd.read_csv(REF / f"1kg_{pop}/chr{c}_{pop}.bim", sep="\t", header=None,
                            usecols=[1, 3], names=["rs", "bp"])
            b = b[b.rs.str.startswith("rs")]
            bims[(pop, c)] = b.drop_duplicates("bp").set_index("bp").rs
        r = bims[(pop, c)].get(p, "")
        if r:
            return r
    return ""


mem = mem[mem.trait_class.isin(["direct", "liver_enzyme", "sequela"])]
rows = []
for lin in LIN:
    col = f"spec_{lin}"
    # cross-study support: distinct liver-class studies with >= 0.1 CS posterior on a variant
    # in the same lineage-specific peak (tie-break among CSs with equal in-peak mass)
    sup = mem[(mem[col] == 1) & (mem.w >= 0.1)].groupby(f"pk_{lin}").study.nunique()
    c = cs[cs[f"obs_{col}"] > 0].copy()
    c["oe"] = c[f"obs_{col}"] / c[f"exp_{col}"]
    top = (mem[mem[col] == 1].sort_values("w", ascending=False).drop_duplicates("cs_uid")
           .set_index("cs_uid"))
    c = c[c.cs_uid.isin(top.index)]
    c["top_pk"] = c.cs_uid.map(top[f"pk_{lin}"]).astype(int)
    c["n_studies_same_peak"] = c.top_pk.map(sup).fillna(0).astype(int)
    c["mass_r"] = c[f"obs_{col}"].round(2)
    c = c.sort_values(["mass_r", "n_studies_same_peak", "lead_pip"], ascending=False)
    seen, n = set(), 0
    for r in c.itertuples():
        mm = mem[(mem.cs_uid == r.cs_uid) & (mem[col] == 1)].sort_values("w", ascending=False)
        v = mm.iloc[0]
        vk = (int(v.chromosome), int(v.position))
        if r.top_pk in seen:
            continue
        seen.add(r.top_pk)
        peak = pk.loc[int(v[f"pk_{lin}"])]
        genes = cvh[(cvh.chromosome == vk[0]) & (cvh.position == vk[1])]
        gl = ";".join(f"{g.gene}({g.gwas_name},PP.H4={g.pp_h4:.2f},w={g.w:.2f})"
                      for g in genes.sort_values("w", ascending=False).itertuples())
        rows.append({"lineage": lin, "rank": n + 1, "study": r.study, "trait": r.trait,
                     "trait_class": r.trait_class, "ancestry": r.ancestry, "cs_uid": r.cs_uid,
                     "cs_mass_in_lineage_specific_peaks": getattr(r, f"obs_{col}"),
                     "cs_expected": getattr(r, f"exp_{col}"), "oe": r.oe,
                     "n_cs_variants_in_specific_peaks": len(mm),
                     "n_liver_studies_same_peak": r.n_studies_same_peak,
                     "cs_lead_pip": r.lead_pip, "cs_n_members": r.n_members,
                     "cs_coding_splice_mass": r.coding_splice_mass,
                     "variant_hg19": f"{v.chromosome}:{v.position}:{v.allele1}:{v.allele2}",
                     "rsid_1kg": rsid(int(v.chromosome), int(v.position), r.ancestry),
                     "variant_hg38": f"chr{v.chromosome}:{v.pos38}",
                     "variant_pip_raw": v.susie_pip, "variant_pip_norm": v.pip_norm,
                     "peak_id": peak.peak_id, "peak_hg38": f"{peak.chrom}:{peak.start0}-{peak.end}",
                     "peak_both_cohorts": bool(peak.both), "tss_dist_bp": int(v.tss_dist),
                     "nearest_pc_gene": v.nearest_gene,
                     "hep_da_class": v.da_hepatocyte, "stellate_da_class": v.da_stellate,
                     "coloc_r2_genes_at_variant": gl})
        n += 1
        if n >= TOPN:
            break
out = pd.DataFrame(rows)
out.to_csv(RES / "q4_worked_loci.tsv", sep="\t", index=False)
pd.set_option("display.width", 300)
print(out[["lineage", "rank", "study", "cs_mass_in_lineage_specific_peaks", "oe", "variant_hg19",
           "rsid_1kg", "variant_pip_norm", "n_liver_studies_same_peak", "nearest_pc_gene",
           "coloc_r2_genes_at_variant"]]
      .round(3).to_string())
