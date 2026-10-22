"""H41: draw the negative-control knockdown genes from LINCS metadata only.

Usage: python3 -I h41_draw_negative_controls.py <downloads_untrusted_dir> <out_tsv>

Reads GSE92742 sig_info and gene_info (metadata). No expression value is read.
Matched on hairpin count to the candidate TFs present in HEPG2 (counted here),
excluding the candidate list, an extended list of hepatocyte transcription
factors and nuclear receptors, and genes that are themselves L1000 landmarks.
Seed 20260923 (numpy default_rng). Pools are sorted by symbol before drawing.
"""
import os
import sys

import numpy as np
import pandas as pd

SEED = 20260923
CANDIDATES = ["HNF4A", "HNF1A", "HNF1B", "FOXA1", "FOXA2", "FOXA3", "CEBPA", "CEBPB",
              "ONECUT1", "ONECUT2", "NR1H4", "PPARA", "RXRA", "ESRRA", "ESRRG", "HNF4G", "NR5A2"]
EXTENDED_EXCLUDE = ["ONECUT3", "PROX1", "HHEX", "GATA4", "GATA6", "TBX3", "CEBPD", "CEBPG",
                    "CEBPE", "NR1H2", "NR1H3", "NR1I2", "NR1I3", "NR0B2", "PPARD", "PPARG",
                    "PPARGC1A", "PPARGC1B", "RXRB", "RXRG", "ESRRB", "NR2F1", "NR2F2", "NR3C1",
                    "XBP1", "ATF4", "ATF5", "ATF6", "CREB3L3", "MLXIPL", "MLXIP", "MLX", "SREBF1",
                    "SREBF2", "KLF15", "TCF7L2", "NFE2L2", "AHR", "ARNT", "ZHX2", "HLF", "DBP",
                    "TEF", "NFIL3", "THRA", "THRB", "VDR", "NR1D1", "NR1D2", "RORA", "RORC"]


def main(dl_dir, out_tsv):
    if os.path.exists(out_tsv):
        sys.exit(f"refusing to overwrite {out_tsv}")
    si = pd.read_csv(os.path.join(dl_dir, "GSE92742_Broad_LINCS_sig_info.txt.gz"), sep="\t", dtype=str)
    gi = pd.read_csv(os.path.join(dl_dir, "GSE92742_Broad_LINCS_gene_info.txt.gz"), sep="\t", dtype=str)
    landmarks = set(gi.loc[gi.pr_is_lm == "1", "pr_gene_symbol"])
    sh = si[(si.cell_id == "HEPG2") & (si.pert_type == "trt_sh")]
    k = sh.groupby("pert_iname").pert_id.nunique()
    present = [g for g in CANDIDATES if g in k.index]
    need = k.loc[present].value_counts().sort_index()
    exclude = set(CANDIDATES) | set(EXTENDED_EXCLUDE) | landmarks
    rng = np.random.default_rng(SEED)
    rows = []
    for kk, n in need.items():
        pool = sorted(g for g, v in k.items() if v == kk and g not in exclude)
        pick = rng.choice(len(pool), size=int(n), replace=False)
        for i in sorted(pick):
            rows.append({"gene": pool[i], "n_hairpins": int(kk), "pool_size": len(pool)})
    out = pd.DataFrame(rows)
    out.to_csv(out_tsv, sep="\t", index=False)
    print("candidates present:", present)
    print("hairpin counts needed:", need.to_dict())
    print(out.to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
