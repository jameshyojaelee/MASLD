"""H58 step 04: Figure 2A support count (descriptive, no test). PRESPEC_H58.md
"Figure 2A support count". Appends a section to RESULT_H58.md."""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import h58_common as h  # noqa: E402


def windows(cs: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """[min, max] of member positions per (study, chromosome, locus string); a
    window spanning more than 1 Mb is split at its largest gap (repeatedly)."""
    rows, n_split = [], 0
    for (study, ch, loc), b in cs.groupby(["study", "chromosome", "locus"]):
        stack = [np.sort(b.position.to_numpy())]
        while stack:
            p = stack.pop()
            if p[-1] - p[0] > 1_000_000:
                n_split += 1
                j = int(np.argmax(np.diff(p)))
                stack += [p[:j + 1], p[j + 1:]]
            else:
                rows.append(dict(study=study, chrom=str(ch), locus=loc, mn=int(p[0]), mx=int(p[-1])))
    return pd.DataFrame(rows), n_split


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    out = Path(ap.parse_args().out)
    h.check_prespec(out)
    tier = pd.read_csv(h.TIER, sep="\t")
    main_t = tier[tier.placement == "main"]
    MAIN = set(main_t.study_name)
    reg = pd.read_csv(h.REGISTRY, sep="\t").set_index("study_name")
    cs = pd.read_csv(h.CREDSETS, usecols=["chromosome", "position", "study", "locus"],
                     dtype={"locus": str, "chromosome": str})
    cs = cs[cs.study.isin(MAIN)].copy()
    cs["position"] = cs.position.astype(float).astype(int)
    win, n_split = windows(cs)
    strata_cs = set(win.study)
    co = pd.read_csv(h.COLOC, usecols=["gwas_name", "gene", "PP.H4.susie", "top_snp"])
    co = co[co.gwas_name.isin(MAIN) & (co.gene.fillna("").str.strip() != "")]
    gmax = co.groupby("gene")["PP.H4.susie"].max()
    multi = set(gmax[gmax > 0.5].index)
    pos = co[co["PP.H4.susie"] > 0.5].copy()
    pos["chrom"] = pos.top_snp.str.split(":").str[0].astype(str)
    pos["bp"] = pos.top_snp.str.split(":").str[1].astype(float)
    wd, wa = defaultdict(list), defaultdict(list)
    for r in win.itertuples():
        wd[(r.study, r.chrom)].append((r.mn, r.mx))
        wa[r.chrom].append((r.mn, r.mx))

    def inside(lst, bp, pad):
        return any(lo - pad <= bp <= hi + pad for lo, hi in lst)

    res = dict(main_strata=len(MAIN), main_strata_with_cs=len(strata_cs),
               main_strata_without_cs=sorted(MAIN - strata_cs), coloc_strata=co.gwas_name.nunique(),
               genes_tested=co.gene.nunique(), multi_signal_genes=len(multi),
               positive_pairs=len(pos), windows=len(win), windows_split=n_split)
    for pad in (0, 250_000):
        pos[f"same_{pad}"] = [inside(wd[(s, c)], b, pad) for s, c, b in zip(pos.gwas_name, pos.chrom, pos.bp)]
        pos[f"any_{pad}"] = [inside(wa[c], b, pad) for c, b in zip(pos.chrom, pos.bp)]
        res[f"genes_same_study_window_pad{pad}"] = int(pos.groupby("gene")[f"same_{pad}"].any().sum())
        res[f"genes_any_window_pad{pad}"] = int(pos.groupby("gene")[f"any_{pad}"].any().sum())
    gs = pos.groupby("gene").gwas_name.apply(set)
    res["genes_all_positive_pairs_in_strata_without_cs"] = int(sum(1 for s in gs if not (s & strata_cs)))
    for st in ["UKBB_ALT", "UKBB_AST", "UKBB_GGT"]:
        res[f"windows_{st}"] = int((win.study == st).sum())
    # per ancestry-by-trait cell
    cell = main_t.set_index("study_name").trait.astype(str) + "|" + reg.reindex(main_t.study_name).ancestry.astype(str).to_numpy()
    win["cell"] = win.study.map(cell)
    pos["cell"] = pos.gwas_name.map(cell)
    cells = pd.DataFrame({
        "finemapped_windows": win.groupby("cell").size(),
        "strata": pd.Series(cell).groupby(cell.values).size(),
        "multi_signal_genes_positive": pos.groupby("cell").gene.nunique(),
        "genes_in_same_study_window_pad250kb": pos[pos.same_250000].groupby("cell").gene.nunique(),
    }).fillna(0).astype(int)
    cells.index.name = "trait|ancestry"
    cells.to_csv(out / "FIG2A_cells.tsv", sep="\t")
    pos.to_csv(out / "FIG2A_positive_pairs.tsv", sep="\t", index=False)
    pd.Series(res).to_csv(out / "FIG2A_support.tsv", sep="\t", header=False)
    txt = ["", "## Figure 2A support count (descriptive, no test)", ""]
    txt += [f"- {k}: {v}" for k, v in res.items()]
    txt += ["", "Per ancestry-by-trait cell:", "", "```", cells.to_string(), "```"]
    with open(out / "RESULT_H58.md", "a") as fh:
        fh.write("\n".join(txt) + "\n")
    print("\n".join(txt))
    return 0


if __name__ == "__main__":
    sys.exit(main())
