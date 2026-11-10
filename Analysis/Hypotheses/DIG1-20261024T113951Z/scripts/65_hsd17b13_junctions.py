"""DIG1 job 65. HSD17B13 splice junctions per release sample from STAR SJ.out.tab (unique reads, column 7).
rs72613567 (A -> AA, chr4:87310241, GRCh38) sits at the donor site of exon 6 (GENCODE v49 canonical exon
87310243-87310359, minus strand). Readouts per sample, written before genotypes were read:
  skip6 = exon 5 -> exon 7 junction (intron 87305309-87313822);
  canon67 = exon 6 -> exon 7 junction (87305309-87310242);
  skip_ratio = skip6 / (skip6 + canon67); every other junction with a start or end in 87305309-87313822 is
  kept in the long table.
Usage: python -I 65_hsd17b13_junctions.py <new_out_dir>
"""
import sys
from pathlib import Path

import pandas as pd

assert len(sys.argv) == 2
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=False)
B = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results")
S = pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/out_rna_genotypes/samples.tsv", sep="\t")
rows, summ = [], []
for s in S.itertuples(index=False):
    f = B / s.dataset / "alignments/star" / s.sample_id / f"{s.sample_id}.SJ.out.tab"
    skip6 = canon67 = 0
    with open(f) as fh:
        for line in fh:
            if not line.startswith("chr4\t"):
                continue
            c = line.split("\t"); a, b, n = int(c[1]), int(c[2]), int(c[6])
            if 87300000 <= a <= 87335000 and 87300000 <= b <= 87335000:
                rows.append(dict(sample_id=s.sample_id, dataset=s.dataset, start=a, end=b, unique=n, annotated=int(c[5])))
                if (a, b) == (87305309, 87313822):
                    skip6 = n
                if (a, b) == (87305309, 87310242):
                    canon67 = n
    summ.append(dict(sample_id=s.sample_id, dataset=s.dataset, skip6=skip6, canon67=canon67,
                     skip_ratio=skip6 / (skip6 + canon67) if skip6 + canon67 else float("nan")))
pd.DataFrame(rows).to_csv(OUT / "hsd17b13_junctions_long.tsv.gz", sep="\t", index=False)
P = pd.DataFrame(summ); P.to_csv(OUT / "hsd17b13_skip6.tsv", sep="\t", index=False)
print(P.describe().to_string())
print(pd.cut(P.skip_ratio, [-0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1]).value_counts().sort_index().to_string())
