"""DIG1 job 40 (exploratory; follows job 37). COLOC r2 places HSPA4 (ALT, GGT, albumin), KRT8 (AST, GGT,
platelets) and ACTG1 (ALT, AST, GGT in EUR, AFR, EAS) at liver-trait loci with liver eQTLs (PP4 > 0.5).
COLOC gives no direction. Here: does the allele that raises liver expression raise or lower the trait?
Data: Broadaway liver eQTL marginal results (EUR, n 1183) and the reformatted GWAS sumstats used by
COLOC (hg19; allele1 = effect allele). Alleles aligned as in 06_susie_coloc.R (match or swap; others
dropped; palindromic A/T and C/G SNPs dropped here).
Per pair, within +-1 Mb of the SuSiE eQTL signal SNP: GWAS z and eQTL z at the eQTL signal SNP and at
the GWAS signal SNP (hit2 and hit1 from susie_coloc_signal_pairs), and the correlation of GWAS z with
eQTL z across SNPs with eQTL p < 1e-4. Descriptive only.
Usage: python -I 40_coloc_direction.py <new_out_dir>
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
FM = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping")
EQ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/broadaway_eqtl")
PAIRS = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
GENES = {"HSPA4": "ENSG00000170606", "KRT8": "ENSG00000170421", "ACTG1": "ENSG00000184009"}
WIN = 1_000_000
PAL = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

pairs = pd.read_csv(PAIRS, sep="\t")
pairs = pairs[pairs.gene.isin(GENES) & (pairs["PP.H4"] > 0.5)].copy()
reg = pd.read_csv(FM / "config/gwas_registry.tsv", sep="\t")
reg = reg.set_index("study_name")
rows = []
for gene, ensg in GENES.items():
    sub = pairs[pairs.gene == gene]
    if sub.empty:
        continue
    chrom = int(sub.chr.iloc[0])
    eq = []
    for ch in pd.read_csv(EQ / f"chr{chrom}_marginal_summary_results.tsv", sep="\t", chunksize=2_000_000,
                          usecols=["CHR", "POS", "NEA", "EA", "EAF", "Beta", "SE", "PVAL", "ENSG"]):
        eq.append(ch[ch.ENSG == ensg])
    eq = pd.concat(eq)
    eq["z_eqtl"] = eq.Beta / eq.SE
    for _, p in sub.iterrows():
        gw_path = FM / reg.loc[p.gwas, "sumstats_path"]
        h1, h2 = int(str(p.hit1).split(":")[1]), int(str(p.hit2).split(":")[1])
        lo, hi = h2 - WIN, h2 + WIN
        gw = []
        for ch in pd.read_csv(gw_path, sep="\t", chunksize=4_000_000):
            ch = ch[(ch.chromosome.astype(str) == str(chrom)) & ch.position.between(lo, hi)]
            gw.append(ch)
        gw = pd.concat(gw)
        m = eq[eq.POS.between(lo, hi)].merge(gw, left_on="POS", right_on="position")
        m["EA"] = m.EA.str.upper(); m["NEA"] = m.NEA.str.upper()
        match = (m.EA == m.allele1) & (m.NEA == m.allele2)
        swap = (m.EA == m.allele2) & (m.NEA == m.allele1)
        m = m[match | swap].copy()
        m["z_gwas"] = np.where(((m.EA == m.allele1) & (m.NEA == m.allele2)), 1, -1) * m.beta / m.se
        m = m[[(a, b) not in PAL for a, b in zip(m.EA, m.NEA)]]

        def at(pos):
            r = m[m.POS == pos]
            return (np.nan, np.nan, np.nan) if r.empty else (r.z_eqtl.iloc[0], r.z_gwas.iloc[0], r.EAF.iloc[0])
        ze2, zg2, f2 = at(h2); ze1, zg1, f1 = at(h1)
        s = m[m.PVAL < 1e-4]
        rows.append(dict(gene=gene, gwas=p.gwas, ancestry=reg.loc[p.gwas, "ancestry"], PP4=p["PP.H4"], eqtl_hit=p.hit2, gwas_hit=p.hit1,
                         n_snps=len(m), z_eqtl_at_eqtl_hit=ze2, z_gwas_at_eqtl_hit=zg2, eaf_eqtl_hit=f2,
                         z_eqtl_at_gwas_hit=ze1, z_gwas_at_gwas_hit=zg1,
                         n_eqtl_p1e4=len(s), cor_z_eqtl_p1e4=np.corrcoef(s.z_eqtl, s.z_gwas)[0, 1] if len(s) > 5 else np.nan,
                         top_eqtl_snp=m.loc[m.PVAL.idxmin(), "POS"] if len(m) else np.nan,
                         z_gwas_at_top_eqtl=m.loc[m.PVAL.idxmin(), "z_gwas"] if len(m) else np.nan,
                         z_eqtl_top=m.loc[m.PVAL.idxmin(), "z_eqtl"] if len(m) else np.nan))
        print(rows[-1], flush=True)
R = pd.DataFrame(rows)
R["direction_at_eqtl_hit"] = np.sign(R.z_eqtl_at_eqtl_hit * R.z_gwas_at_eqtl_hit)
R.to_csv(OUT / "coloc_direction.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(R.round(3).to_string())
