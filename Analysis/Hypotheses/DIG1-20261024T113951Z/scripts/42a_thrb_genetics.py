"""DIG1 job 42a (exploratory; follows job 41). THRB colocalizes with UK Biobank GGT (SuSiE PP.H4 > 0.999;
signal SNP 3:24520283, hg19). At that SNP the G allele raises liver THRB expression (eQTL z +5.9) and
lowers GGT (z -5.7). Questions:
  1. Across all 50 GWAS in the registry, what does the THRB-raising allele do to each liver trait?
  2. Is the SNP also an eQTL for other genes on chromosome 3 (alternative gene at the locus)?
  3. Across THRB eQTL SNPs (p < 1e-4) within +-1 Mb, does eQTL z correlate with GWAS z in each study?
Data: Broadaway liver eQTL marginal results (EUR, n 1,183); registry sumstats (allele1 = effect allele).
Alleles aligned to eQTL EA/NEA (match or swap; others dropped; palindromic SNPs dropped).
Descriptive only; no test family.
Usage: python -I 42a_thrb_genetics.py <new_out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
FM = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping")
EQ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/broadaway_eqtl")
CO = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
ENSG, CHR, SNP = "ENSG00000151090", 3, 24520283
WIN = 1_000_000
PAL = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

reg = pd.read_csv(FM / "config/gwas_registry.tsv", sep="\t").set_index("study_name")
tier = pd.read_csv(FM / "config/gwas_trait_tier.tsv", sep="\t").set_index("study_name")
ab = pd.read_csv(CO / "susie_coloc_all_gwas.csv")
ab = ab[ab.gene == "THRB"].set_index("gwas_name")[["PP.H4.abf", "PP.H4.susie", "method", "top_snp", "top_snp_PP"]]

# eQTL: every gene at the SNP, and THRB across the window
at_snp, thrb = [], []
for ch in pd.read_csv(EQ / f"chr{CHR}_marginal_summary_results.tsv", sep="\t", chunksize=3_000_000,
                      usecols=["POS", "NEA", "EA", "EAF", "Beta", "SE", "PVAL", "GeneSymbol", "ENSG"]):
    at_snp.append(ch[ch.POS == SNP])
    thrb.append(ch[(ch.ENSG == ENSG) & ch.POS.between(SNP - WIN, SNP + WIN)])
at_snp = pd.concat(at_snp); at_snp["z_eqtl"] = at_snp.Beta / at_snp.SE
at_snp.sort_values("PVAL").to_csv(OUT / "eqtl_all_genes_at_snp.tsv", sep="\t", index=False)
thrb = pd.concat(thrb); thrb["z_eqtl"] = thrb.Beta / thrb.SE
thrb["EA"] = thrb.EA.str.upper(); thrb["NEA"] = thrb.NEA.str.upper()
ref = thrb[thrb.POS == SNP].iloc[0]
print("THRB eQTL at SNP", ref.EA, ref.NEA, round(ref.z_eqtl, 2), "EAF", ref.EAF, flush=True)
print("top THRB eQTL", thrb.loc[thrb.PVAL.idxmin(), ["POS", "z_eqtl", "PVAL"]].to_dict(), flush=True)

rows, loc = [], []
for s in reg.index:
    gw = []
    for ch in pd.read_csv(FM / reg.loc[s, "sumstats_path"], sep="\t", chunksize=4_000_000, dtype={"chromosome": str}):
        gw.append(ch[(ch.chromosome == str(CHR)) & ch.position.between(SNP - WIN, SNP + WIN)])
    gw = pd.concat(gw)
    m = thrb.merge(gw, left_on="POS", right_on="position")
    a1, a2 = m.allele1.str.upper(), m.allele2.str.upper()
    same = (m.EA == a1) & (m.NEA == a2); swap = (m.EA == a2) & (m.NEA == a1)
    m = m[(same | swap) & [(a, b) not in PAL for a, b in zip(m.EA, m.NEA)]].copy()
    m["z_gwas"] = np.where(((m.EA == m.allele1.str.upper()) & (m.NEA == m.allele2.str.upper())), 1, -1) * m.beta / m.se
    r = m[m.POS == SNP]
    s5 = m[m.PVAL < 1e-4]
    rows.append(dict(gwas=s, trait=tier.trait.get(s, np.nan), tier=tier.tier.get(s, np.nan), ancestry=reg.loc[s, "ancestry"],
                     z_eqtl_at_snp=ref.z_eqtl, z_gwas_thrb_up_allele=r.z_gwas.iloc[0] if len(r) else np.nan,
                     n_snps=len(m), n_eqtl_p1e4=len(s5), cor_z_p1e4=np.corrcoef(s5.z_eqtl, s5.z_gwas)[0, 1] if len(s5) > 5 else np.nan,
                     min_gwas_p_window=m.pval.min() if len(m) else np.nan,
                     PP4_abf=ab["PP.H4.abf"].get(s, np.nan), PP4_susie=ab["PP.H4.susie"].get(s, np.nan), coloc_method=ab["method"].get(s, np.nan),
                     top_snp=ab["top_snp"].get(s, np.nan)))
    print(rows[-1], flush=True)
R = pd.DataFrame(rows)
R.to_csv(OUT / "thrb_allele_across_gwas.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 100)
print(R.sort_values(["tier", "trait"]).round(3).to_string())
print(at_snp.sort_values("PVAL")[["GeneSymbol", "ENSG", "EA", "NEA", "z_eqtl", "PVAL"]].head(15).to_string())
