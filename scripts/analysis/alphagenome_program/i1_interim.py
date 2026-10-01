#!/usr/bin/env python3
"""Interim read of a partial gate score table against the archive, while the GPU job still runs.

Diagnostic only. It prints the running Spearman, sign agreement and rejected-window count from the
rows flushed so far so that a systematic problem surfaces in minutes rather than after the full
1.7-hour run. Nothing it prints enters RESULTS.md; the reported numbers all come from
`i1_04_gate_stats.py` on the complete table.

    python i1_interim.py <run_dir> <length_bp>
"""
import sys, numpy as np, pandas as pd, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import i1_common as C
R = pathlib.Path(sys.argv[1]); L = sys.argv[2]
g = pd.read_csv(R/"inputs/gate_variants.tsv", sep="\t").set_index("key")
d = pd.read_csv(R/f"raw/gate_{L}.tsv", sep="\t")
m = d.join(g[["archived_atac_liver","archived_dnase_liver"]], on="key")
a = m.archived_atac_liver.values; l = m.local_atac_liver.values
print("INTERIM n=%d blocks=%d" % (len(m), m.block_1mb.nunique()))
print("  ATAC  rho %.5f  sign %.4f  medAbsDiff %.5f  ratio %.4f" % (
    C.fast_spearman(a,l), np.mean(np.sign(a)==np.sign(l)), np.median(np.abs(l-a)), np.median(l/a)))
print("  DNASE rho %.5f" % C.fast_spearman(m.archived_dnase_liver.values, m.local_dnase_liver.values))
print("  archived vs beta_alt %.4f | local vs beta_alt %.4f" % (
    C.fast_spearman(m.beta_alt.values, a), C.fast_spearman(m.beta_alt.values, l)))
rp = R/f"raw/gate_{L}_rejected.tsv"
if rp.exists():
    rd = pd.read_csv(rp, sep="\t")
    print("  rejected so far %d:" % len(rd), dict(rd.reason.value_counts()))
